/*   sequin2.c
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
* File Name:  sequin2.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.719 $
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

#include "sequin.h"
#include <document.h>
#include <sequtil.h>
#include <biosrc.h>
#include <cdrgn.h>
#include <seqsub.h>
#include <tofasta.h>
#include <gather.h>
#include <subutil.h>
#include <suggslp.h>
#include <toasn3.h>
#include <toporg.h>
#include <salfiles.h>
#include <salsap.h>
#include <salign.h>
#include <edutil.h>
#include <vsm.h>
//#include <accentr.h>
//#include <accutils.h>
#include <pmfapi.h>
#include <explore.h>
#include <aliparse.h>
#include <algo/blast/api/twoseq_api.h>
#ifdef WIN_MOTIF
#include <netscape.h>
#endif
#include <actutils.h>
#include <salpanel.h>
#include <findrepl.h>
#include <macrodlg.h>
#include <macroapi.h>

extern EnumFieldAssoc  biosource_genome_simple_alist [];
extern EnumFieldAssoc  biosource_origin_alist [];

static ENUM_ALIST(biomol_nucX_alist)
  {"Genomic DNA",            253},
  {"Genomic RNA",            254},
  {"Precursor RNA",            2},
  {"mRNA [cDNA]",              3},
  {"Ribosomal RNA",            4},
  {"Transfer RNA",             5},
  {"Other-Genetic",            9},
  {"cRNA",                    11},
  {"Transcribed RNA",         13}, 
  {"Transfer-messenger RNA", MOLECULE_TYPE_TMRNA },
  {"ncRNA",                  MOLECULE_TYPE_NCRNA},
END_ENUM_ALIST

static ENUM_ALIST(biomol_nucGen_alist)
  {"Genomic DNA",            253},
  {"Genomic RNA",            254},
END_ENUM_ALIST

static ENUM_ALIST(topology_nuc_alist)
{"Linear",          TOPOLOGY_LINEAR},
{"Circular",        TOPOLOGY_CIRCULAR},
END_ENUM_ALIST

static ENUM_ALIST(molecule_alist)
{"DNA",             Seq_mol_dna },
{"RNA",             Seq_mol_rna },
END_ENUM_ALIST

#define PRINTED_INT_MAX_LEN 15

#define CREATE_FASTA_REQUIRED 0
#define CREATE_FASTA_WARNING  1

/* These functions are for creating, copying, and freeing lists
 * of titles and IDs.
 */
static IDAndTitleEditPtr IDAndTitleEditNew (void)
{
  IDAndTitleEditPtr iatep;
  
  iatep = (IDAndTitleEditPtr) MemNew (sizeof (IDAndTitleEditData));
  if (iatep != NULL)
  {
    iatep->id_list = NULL;
    iatep->title_list = NULL;
    iatep->is_seg = NULL;
    iatep->num_sequences = 0;
    iatep->nuc_only = FALSE;
  }
  return iatep;
}

static void IDAndTitleEditInit (IDAndTitleEditPtr iatep, Int4 new_num_sequences)
{
  Int4 seq_num;
  if (iatep == NULL)
  {
    return;
  }
  
  /* free old lists, if any */
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    iatep->id_list [seq_num] = MemFree (iatep->id_list [seq_num]);
    iatep->title_list [seq_num] = MemFree (iatep->title_list [seq_num]);
  }
  iatep->id_list = MemFree (iatep->id_list);
  iatep->title_list = MemFree (iatep->title_list);
  iatep->is_seg = MemFree (iatep->is_seg);
  
  /* now create blanks for num_sequences entries */
  iatep->num_sequences = MAX (0, new_num_sequences);
  if (iatep->num_sequences > 0)
  {
    iatep->id_list = (CharPtr PNTR) MemNew (iatep->num_sequences * sizeof (CharPtr));
    iatep->title_list = (CharPtr PNTR) MemNew (iatep->num_sequences * sizeof (CharPtr));
    iatep->is_seg = (BoolPtr) MemNew (iatep->num_sequences * sizeof (Boolean));
    for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
    {
      iatep->id_list [seq_num] = NULL;
      iatep->title_list  [seq_num] = NULL;
      iatep->is_seg [seq_num] = FALSE;
    }
  } 
}

static IDAndTitleEditPtr IDAndTitleEditCopy (IDAndTitleEditPtr iatep_orig)
{
  IDAndTitleEditPtr iatep_copy;
  Int4              seq_num;
  
  if (iatep_orig == NULL)
  {
    return NULL;
  }
  
  iatep_copy = IDAndTitleEditNew ();
  if (iatep_copy == NULL)
  {
    return NULL;
  }
  
  IDAndTitleEditInit (iatep_copy, iatep_orig->num_sequences);
  for (seq_num = 0; seq_num < iatep_copy->num_sequences; seq_num++)
  {
    iatep_copy->id_list [seq_num] = StringSave (iatep_orig->id_list [seq_num]);
    iatep_copy->title_list [seq_num] = StringSave (iatep_orig->title_list [seq_num]);
    if (iatep_orig->is_seg != NULL)
    {
      iatep_copy->is_seg [seq_num] = iatep_orig->is_seg [seq_num];
    }
  }
  
  return iatep_copy;
}

NLM_EXTERN IDAndTitleEditPtr IDAndTitleEditFree (IDAndTitleEditPtr iatep)
{
  Int4 i;
  
  if (iatep != NULL)
  {
    for (i = 0; i < iatep->num_sequences; i++)
    {
      iatep->id_list [i] = MemFree (iatep->id_list [i]);
      iatep->title_list [i] = MemFree (iatep->title_list [i]);
    }
    iatep->id_list = MemFree (iatep->id_list);
    iatep->title_list = MemFree (iatep->title_list);
    iatep->is_seg = MemFree (iatep->is_seg);
    iatep = MemFree (iatep);
  }
  return iatep;
}

/* These functions are for applying lists of titles and IDs
 * to a SeqEntry list.
 */
NLM_EXTERN Int4 CountSequencesAndSegments (SeqEntryPtr list, Boolean nuc_only)
{
  Int4         num_seqs = 0;
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
  
  while (list != NULL)
  {
    if (list->data.ptrvalue != NULL)
    {
      if (IS_Bioseq (list))
      {
        bsp = (BioseqPtr) list->data.ptrvalue;
        if (!nuc_only || ISA_na (bsp->mol)) {
          num_seqs ++;
        }
      }
      else if (IS_Bioseq_set (list))
      {
        bssp = (BioseqSetPtr) list->data.ptrvalue;
        num_seqs += CountSequencesAndSegments (bssp->seq_set, nuc_only);
      }
    }
    list = list->next;
  }
  return num_seqs;
}

static BioseqPtr FindNthSequenceInSet (SeqEntryPtr seq_list, Int4 nth, BoolPtr is_seg, Boolean nuc_only)
{
  Int4         pos = 0;
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp;
  SeqEntryPtr  sep;
  
  while (seq_list != NULL && bsp == NULL)
  {
    if (seq_list->data.ptrvalue != NULL)
    {
      if (IS_Bioseq (seq_list) && seq_list->data.ptrvalue != NULL 
          && (!nuc_only || ISA_na(((BioseqPtr)seq_list->data.ptrvalue)->mol)))
      {
        if (nth == pos)
        {
          bsp = seq_list->data.ptrvalue;
        }
        else
        {
          pos ++;
        }
      }
      else if (IS_Bioseq_set (seq_list))
      {
        bssp = (BioseqSetPtr) seq_list->data.ptrvalue;
        if (bssp->_class == BioseqseqSet_class_parts && is_seg != NULL)
        {
          *is_seg = TRUE;
        }
        sep = bssp->seq_set;
        while (sep != NULL && bsp == NULL)
        {
          bsp = FindNthSequenceInSet (sep, nth - pos, is_seg, nuc_only);
          if (bsp == NULL)
          {
            if (IS_Bioseq_set (sep))
            {
              bssp = (BioseqSetPtr) sep->data.ptrvalue;
              pos += CountSequencesAndSegments (bssp->seq_set, nuc_only);
            }
            else if (IS_Bioseq (sep) && (!nuc_only || ISA_na (((BioseqPtr)(sep->data.ptrvalue))->mol)))
            {
              pos ++;
            }
          }
          sep = sep->next;
        }
        if (bsp == NULL && is_seg != NULL)
        {
          *is_seg = FALSE;
        }
      }      
    }
    seq_list = seq_list->next;
  }
  return bsp;
}


static void AddOneSequenceToIatep (BioseqPtr bsp, IDAndTitleEditPtr iatep, Int4 pos)
{
  SeqIdPtr sip;
  SeqDescPtr sdp;

  sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  if (sip != NULL)
  {
    if (sip->choice == SEQID_LOCAL) {
      iatep->id_list [pos] = SeqIdWholeLabel (sip, PRINTID_REPORT);
    } else {
      iatep->id_list [pos] = SeqIdWholeLabel (sip, PRINTID_FASTA_SHORT);
    }
  }
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp != NULL && !StringHasNoText (sdp->data.ptrvalue))
  {
    iatep->title_list [pos] = StringSave (sdp->data.ptrvalue);
  }
}


static void AddSequencesInSetToIatep (SeqEntryPtr seq_list, Int4Ptr nth, Boolean nuc_only, IDAndTitleEditPtr iatep)
{
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
  
  while (seq_list != NULL)
  {
    if (seq_list->data.ptrvalue != NULL)
    {
      if (IS_Bioseq (seq_list) 
          && (bsp = (BioseqPtr) seq_list->data.ptrvalue) != NULL 
          && (!nuc_only || ISA_na(bsp->mol)))
      {
        AddOneSequenceToIatep (bsp, iatep, *nth);
        (*nth)++;
      }
      else if (IS_Bioseq_set (seq_list))
      {
        bssp = (BioseqSetPtr) seq_list->data.ptrvalue;
        AddSequencesInSetToIatep (bssp->seq_set, nth, nuc_only, iatep);
      }      
    }
    seq_list = seq_list->next;
  }
}


NLM_EXTERN IDAndTitleEditPtr SeqEntryListToIDAndTitleEditEx (SeqEntryPtr list, Boolean nuc_only)
{
  IDAndTitleEditPtr iatep;
  Int4              num_sequences, i = 0;
  
  num_sequences = CountSequencesAndSegments (list, nuc_only);
  if (num_sequences == 0)
  {
    return NULL;
  }
  
  iatep = IDAndTitleEditNew ();
  if (iatep == NULL)
  {
    return NULL;
  }

  iatep->nuc_only = nuc_only;
  IDAndTitleEditInit (iatep, num_sequences);
  AddSequencesInSetToIatep (list, &i, nuc_only, iatep);

  return iatep;
}


static IDAndTitleEditPtr SeqEntryListToIDAndTitleEdit (SeqEntryPtr list)
{
  return SeqEntryListToIDAndTitleEditEx (list, FALSE);
}

static void ReplaceIDAndTitleForBioseq (BioseqPtr bsp, SeqIdPtr new_sip, CharPtr title)
{
  SeqDescrPtr sdp;
  SeqEntryPtr sep;
  
  if (bsp == NULL)
  {
    return;
  }
  
  /* replace ID */
  
  if (new_sip != NULL)
  {
    if (SeqIdComp (new_sip, bsp->id) == SIC_YES) 
    {
      /* do nothing, no change */
    }
    else
    {
      if (bsp->id != NULL)
      {
        new_sip->next = bsp->id->next;
        bsp->id->next = NULL;
        bsp->id = SeqIdFree (bsp->id);
      }
      bsp->id = new_sip;
      SeqMgrReplaceInBioseqIndex(bsp);
    }
  }
  else
  {
    bsp->id = SeqIdFree (bsp->id);
  }
    
  /* replace title */
  if (title == NULL)
  {
    title = StringSave ("");
  }
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp == NULL)
  {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_title);
    sdp->data.ptrvalue = title;
  }
  else
  {
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = title;
  }  
}

static void ResetSegSetIDLists (SeqEntryPtr list)
{
  BioseqSetPtr bssp, parts;
  BioseqPtr    seg_bsp;
  SeqEntryPtr  sep;
  SeqLocPtr    loc, next_loc, last_loc;
  
  while (list != NULL)
  {  
    if (IS_Bioseq_set (list) && (bssp = (BioseqSetPtr) list->data.ptrvalue) != NULL)
    {
      if (bssp->_class == BioseqseqSet_class_segset)
      {
        sep = bssp->seq_set;
        seg_bsp = NULL;
        parts = NULL;
        while (sep != NULL && (seg_bsp == NULL || parts == NULL))
        {
          if (IS_Bioseq (sep))
          {
            seg_bsp = sep->data.ptrvalue;
          }
          else if (IS_Bioseq_set (sep))
          {
            parts = sep->data.ptrvalue;
            if (parts != NULL && parts->_class != BioseqseqSet_class_parts)
            {
              parts = NULL;
            }
          }
          sep = sep->next;
        }
        if (seg_bsp != NULL)
        {
          /* remove old location */
          loc = (SeqLocPtr) seg_bsp->seq_ext;
          while (loc != NULL)
          {
            next_loc = loc->next;
            loc->next = NULL;
            loc = SeqLocFree (loc);
            loc = next_loc;
          }
          seg_bsp->seq_ext = NULL;
          /* put in new locations */
          sep = parts->seq_set;
          last_loc = NULL;
          while (sep != NULL)
          {
            if (IS_Bioseq (sep) && sep->data.ptrvalue != NULL)
            {
              loc = SeqLocWholeNew (sep->data.ptrvalue);
              if (loc != NULL)
              {
                if (last_loc == NULL)
                {
                  seg_bsp->seq_ext = loc;
                }
                else
                {
                  last_loc->next = loc;
                }
                last_loc = loc;
              }
            }
            sep = sep->next;
          }
        }
      }
      else
      {
        ResetSegSetIDLists (bssp->seq_set);
      }
    }
    list = list->next;
  }
}


static void ApplyIatepToOneSequence (BioseqPtr bsp, IDAndTitleEditPtr iatep, Int4 pos)
{
  SeqIdPtr new_sip;

  new_sip = NULL;
  if (StringChr (iatep->id_list[pos], '|') != NULL) {
    new_sip = SeqIdParse (iatep->id_list[pos]);
  }
  if (new_sip == NULL) {
    new_sip = MakeSeqID (iatep->id_list [pos]);
  }
  ReplaceIDAndTitleForBioseq (bsp, new_sip, StringSave (iatep->title_list [pos]));

}


static void ApplyIatepToSequencesInSet (SeqEntryPtr seq_list, Int4Ptr nth, Boolean nuc_only, IDAndTitleEditPtr iatep)
{
  BioseqSetPtr bssp;
  
  while (seq_list != NULL)
  {
    if (seq_list->data.ptrvalue != NULL)
    {
      if (IS_Bioseq (seq_list) && seq_list->data.ptrvalue != NULL 
          && (!nuc_only || ISA_na(((BioseqPtr)seq_list->data.ptrvalue)->mol)))
      {
        ApplyIatepToOneSequence (seq_list->data.ptrvalue, iatep, *nth);
        (*nth)++;
      }
      else if (IS_Bioseq_set (seq_list))
      {
        bssp = (BioseqSetPtr) seq_list->data.ptrvalue;
        ApplyIatepToSequencesInSet (bssp->seq_set, nth, nuc_only, iatep);
      }      
    }
    seq_list = seq_list->next;
  }
}


NLM_EXTERN Boolean ApplyIDAndTitleEditToSeqEntryList (SeqEntryPtr list, IDAndTitleEditPtr iatep)
{
  Int4      i = 0;
  
  if (list == NULL || iatep == NULL)
  {
    return FALSE;
  }
  
  if (CountSequencesAndSegments (list, iatep->nuc_only) != iatep->num_sequences)
  {
    return FALSE;
  }
  
  ApplyIatepToSequencesInSet (list, &i, iatep->nuc_only, iatep); 
  ResetSegSetIDLists (list);
  return TRUE;
}


static ValNodePtr BuildModifierTypeList (ValNodePtr type_list, CharPtr new_title, Boolean allow_prot);

NLM_EXTERN void RemoveSourceModifiersFromIdAndTitleEdit (IDAndTitleEditPtr iatep)
{
  ValNodePtr  found_modifiers = NULL, vnp;
  Int4        seq_num;

  if (iatep != NULL)
  {
    /* get list of modifiers from titles */    
    for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
    {
      found_modifiers = BuildModifierTypeList (found_modifiers, 
                                               iatep->title_list [seq_num],
                                               FALSE);
      for (vnp = found_modifiers; vnp != NULL; vnp = vnp->next)
      { 
        if (StringICmp (vnp->data.ptrvalue, "genetic_code") == 0
            || StringICmp (vnp->data.ptrvalue, "organism") == 0
            || StringICmp (vnp->data.ptrvalue, "location") == 0
            || StringICmp (vnp->data.ptrvalue, "gencode_comment") == 0
            || StringICmp (vnp->data.ptrvalue, "moltype") == 0
            || StringICmp (vnp->data.ptrvalue, "topology") == 0)
        {
          continue;
        }
        RemoveValueFromDefline (vnp->data.ptrvalue, iatep->title_list [seq_num]);
      }
      found_modifiers = ValNodeFreeData (found_modifiers);
    }
  }
}


/* this section of code is used to read and parse the taxlist.txt 
 * and lineages.txt files */
static ValNodePtr orglist = NULL;

typedef struct orginfo 
{
  CharPtr taxname;
  CharPtr common;
  Int4    ngcode;
  Int4    mgcode;
  CharPtr div;
  Int4    taxnum;
  CharPtr lineage;
} OrgInfoData, PNTR OrgInfoPtr;

static FILE *OpenSequinDataFile (CharPtr filename)
{
  Char              str [PATH_MAX];
  CharPtr           ptr;
  FILE              *f = NULL;

  if (StringHasNoText (filename))
  {
    return NULL;
  }

  ProgramPath (str, sizeof (str));
  ptr = StringRChr (str, DIRDELIMCHR);
  if (ptr == NULL)
  {
    return NULL;
  }
  
  *ptr = '\0';
  FileBuildPath (str, NULL, filename);
  f = FileOpen (str, "r");
  if (f == NULL) {
    if (GetAppParam ("NCBI", "NCBI", "DATA", "", str, sizeof (str))) {
      FileBuildPath (str, NULL, filename);
      f = FileOpen (str, "r");
    }
  }
  return f;
}

static OrgInfoPtr FindByTaxNum (Int4 taxnum)
{
  ValNodePtr vnp;
  OrgInfoPtr oip;
  
  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
    oip = (OrgInfoPtr) vnp->data.ptrvalue;
    if (oip != NULL && oip->taxnum == taxnum)
    {
      return oip;
    }
  }
  return NULL;
}

static OrgInfoPtr FindByTaxName (CharPtr taxname)
{
  ValNodePtr vnp;
  OrgInfoPtr oip;
  
  if (StringHasNoText (taxname))
  {
    return NULL;
  }
  
  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
    oip = (OrgInfoPtr) vnp->data.ptrvalue;
    if (oip != NULL && StringICmp (oip->taxname, taxname) == 0)
    {
      return oip;
    }
  }
  return NULL;
}

static void AddLineagesToOrganismList (void)
{
  ReadBufferData    rbd;
  CharPtr           line;
  CharPtr           ptr;
  FILE              *f;
  OrgInfoPtr        oip;
  Int4              taxnum;

  /* can only add lineages to existing list */
  if (orglist == NULL) return;

  /* now read in lineages */
  f = OpenSequinDataFile ("lineages.txt");

  if (f != NULL) 
  {
    rbd.fp = f;
    rbd.current_data = NULL;
    line = AbstractReadFunction (&rbd);
    line = AbstractReadFunction (&rbd);
    while (line != NULL)
    {
      ptr = StringChr (line, '\t');
      if (ptr != NULL) 
      {
        *ptr = '\0';
        if (StrToLong (line, &taxnum)) 
        {
          oip = FindByTaxNum (taxnum);
          if (oip != NULL)
          {
            oip->lineage = StringSave (ptr + 1);
          }
        }
      }
    	line = AbstractReadFunction (&rbd);
    }
    FileClose (f);
  }
}

static CharPtr GetNextToken (CharPtr PNTR pstart)
{
  CharPtr pend;
  CharPtr newval = NULL;
  
  if (pstart == NULL || *pstart == NULL)
  {
    return NULL;
  }
  
  pend = StringChr (*pstart, '\t');
  if (pend != NULL)
  {
    *pend = 0;
  }
  newval = StringSave (*pstart);
  if (pend == NULL)
  {
    *pstart = NULL;
  }
  else
  {
    *pstart = pend + 1;
  }
  return newval;
}

static void LoadOrganismList (void)
{
  ReadBufferData    rbd;
  CharPtr           line;
  CharPtr           p_start, numval;
  FILE              *f;
  OrgInfoPtr        oip;

  if (orglist != NULL) return;
  
  f = OpenSequinDataFile ("taxlist.txt");

  if (f != NULL) {
    rbd.fp = f;
    rbd.current_data = NULL;
    line = AbstractReadFunction (&rbd);
    line = AbstractReadFunction (&rbd);
    while (line != NULL)
    {
      oip = (OrgInfoPtr) MemNew (sizeof (OrgInfoData));
      if (oip != NULL)
      {
        p_start = line;
        /* read in tax name */
        oip->taxname = GetNextToken (&p_start);
        
        /* read in common name */
        oip->common = GetNextToken (&p_start);
         
        /* read in nuclear genetic code */
        numval = GetNextToken (&p_start);
        if (numval != NULL)
        {
          StrToLong (numval, &(oip->ngcode));
          numval = MemFree (numval);
        }
        /* read in mitochondrial genetic code */
        numval = GetNextToken (&p_start);
        if (numval != NULL)
        {
          StrToLong (numval, &(oip->mgcode));
          numval = MemFree (numval);
        }
        
        /* read in div */
        oip->div = GetNextToken (&p_start);
        
        /* read in taxnum */
        numval = GetNextToken (&p_start);
        if (numval != NULL)
        {
          StrToLong (numval, &(oip->taxnum));
          numval = MemFree (numval);
        }
                
        ValNodeAddPointer (&orglist, 0, oip);
      }
      line = MemFree (line);
    	line = AbstractReadFunction (&rbd);
    }
    FileClose (f);
  }
  AddLineagesToOrganismList ();
}

/* This section of code is used for determining genetic codes based on
 * FASTA-defline values.
 */
#define USE_NUCLEAR_GENETIC_CODE       1
#define USE_MITOCHONDRIAL_GENETIC_CODE 2
#define USE_OTHER_GENETIC_CODE         3

static Int4 UseGeneticCodeForLocation (CharPtr location)
{
  if (StringHasNoText (location))
  {
    return USE_NUCLEAR_GENETIC_CODE;
  }
  else if (StringICmp (location, "Mitochondrion") == 0
           || StringICmp (location, "Kinetoplast") == 0
           || StringICmp (location, "Hydrogenosome") == 0)
  {
    return USE_MITOCHONDRIAL_GENETIC_CODE;
  }
  else if (StringICmp (location, "Chloroplast") == 0
           || StringICmp (location, "Chromoplast") == 0
           || StringICmp (location, "plastid") == 0
           || StringICmp (location, "cyanelle") == 0
           || StringICmp (location, "apicoplast") == 0
           || StringICmp (location, "leucoplast") == 0
           || StringICmp (location, "proplastid") == 0)
  {
    return USE_OTHER_GENETIC_CODE;
  }
  else
  {
    return USE_NUCLEAR_GENETIC_CODE;
  }
}


static Int4 GetGeneticCodeForTaxNameAndLocation (CharPtr taxname, CharPtr location)
{
  ValNodePtr vnp;
  OrgInfoPtr oip;
  Int4       use_code;
  
  use_code = UseGeneticCodeForLocation (location);
  if (use_code == USE_OTHER_GENETIC_CODE)
  {
    return 11;
  }
  else if (StringHasNoText (taxname))
  {
    return -1;
  }
  
  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    oip = (OrgInfoPtr) vnp->data.ptrvalue;
    if (StringICmp (oip->taxname, taxname) == 0)
    {
      if (use_code == USE_NUCLEAR_GENETIC_CODE)
      {
        return oip->ngcode;
      }
      else
      {
        return oip->mgcode;
      }
    }
  }
  
  return -1;
}

static CharPtr GeneticCodeStringFromIntAndList (Int4 num, ValNodePtr list)
{
  while (list != NULL)
  {
    if (list->choice == num)
    {
      return list->data.ptrvalue;
    }
    list = list->next;
  }
  return NULL;
}


/* these functions deal with commonly asked questions about package types - 
 * which ones are sets, which ones are single sequences, which ones have
 * which default molecule types.
 */
static Boolean PackageTypeIsSet (Int2 seqPackage)
{
  if (seqPackage == SEQ_PKG_POPULATION
      || seqPackage == SEQ_PKG_PHYLOGENETIC 
      || seqPackage == SEQ_PKG_MUTATION
      || seqPackage == SEQ_PKG_ENVIRONMENT
      || seqPackage == SEQ_PKG_GENBANK
      || seqPackage == SEQ_PKG_TSA)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
  
}

static Boolean PackageTypeIsSingle (Int2 seqPackage)
{
  if (seqPackage == SEQ_PKG_SINGLE
      || seqPackage == SEQ_PKG_GAPPED)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/* These functions are used to find titles in SeqEntries */
static void FindFirstTitle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  CharPtr PNTR  ttlptr;

  if (mydata == NULL) return;
  ttlptr = (CharPtr PNTR) mydata;
  if (*ttlptr != NULL) return;
  *ttlptr = SeqEntryGetTitle (sep);
}

static void FindFirstSeqEntryTitle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  SeqEntryPtr PNTR  sepptr;

  if (mydata == NULL) return;
  sepptr = (SeqEntryPtr PNTR) mydata;
  if (*sepptr != NULL) return;
  if (SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL) != NULL) {
   *sepptr = sep;
  }
}

/* These functions are used to change the values of modifiers in definition lines */

extern void MakeSearchStringFromAlist (CharPtr str, CharPtr name)

{
  Char     ch;
  CharPtr  ptr;

  StringCpy (str, "[");
  StringCat (str, name);
  StringCat (str, "=");
  ptr = str;
  ch = *ptr;
  while (*ptr != '\0') {
    *ptr = TO_LOWER (ch);
    ptr++;
    ch = *ptr;
  }
}

/* This section of code is used for parsing well-formatted definition lines.
 */
typedef struct modifieralias 
{
  CharPtr alias;
  CharPtr modifier;
} ModifierAlias, PNTR ModifierAliasPtr;

static ModifierAlias alias_list [] =
{
  { "org", "organism" },
  { "mol-type", "moltype" },
  { "mol_type", "moltype" },
  { "note", "note-orgmod" },
  { "comment", "note-orgmod" },
  { "common-name", "common name"},
  { "subsource", "note-subsrc" },
  { "technique", "tech" },
  { "prot", "protein" },
  { "prot_desc", "protein_desc" }
};

static Int4 num_aliases = sizeof (alias_list) / sizeof (ModifierAlias);

static CharPtr protein_modifier_names [] = 
{
  "gene",
  "gene_syn",
  "protein",
  "protein_desc",
  "note",
  "comment",
  "orf",
  "function",
  "EC_number"
};

static Int4 num_protein_modifier_names = sizeof (protein_modifier_names) / sizeof (CharPtr);

static CharPtr GetCanonicalName (CharPtr mod_name)
{
  Int4 j;
  Uint1   subtype;
  
  if (StringHasNoText (mod_name))
  {
    return StringSave ("");
  }
  
  for (j = 0; j < num_protein_modifier_names; j++) {
    if (StringsAreEquivalent (mod_name, protein_modifier_names[j])) {
      return StringSave (protein_modifier_names[j]);
    }
  }

  subtype = EquivalentOrgMod (mod_name);
  if (subtype != 0) {
    return StringSave (GetOrgModQualName (subtype));
  }

  if (StringISearch (mod_name, "primer") == NULL) {
    subtype = EquivalentSubSource (mod_name);
  } else {
    subtype = EquivalentSubSourceEx (mod_name, TRUE);
  }
  if (subtype != 0) {
    return StringSave (GetSubsourceQualName (subtype));
  }
 
  for (j = 0; j < num_aliases; j++)
  {
    if (StringsAreEquivalent (alias_list [j].alias, mod_name))
    {
      return StringSave (alias_list [j].modifier);
    }
  }
  return StringSave (mod_name);
}

typedef enum {
    eModifierType_SourceQual = 0,
    eModifierType_Organism,
    eModifierType_Location,
    eModifierType_Lineage,
    eModifierType_GeneticCode,
    eModifierType_GeneticCodeComment,
    eModifierType_NucGeneticCode,
    eModifierType_MitoGeneticCode,
    eModifierType_MolType,
    eModifierType_Molecule,
    eModifierType_Origin,
    eModifierType_Topology,
    eModifierType_CommonName,
    eModifierType_Technique,
    eModifierType_Protein
} EModifierType;

typedef struct modifierinfo 
{
  CharPtr       name;
  Uint1         subtype;
  CharPtr       value;
  EModifierType modtype;
} ModifierInfoData, PNTR ModifierInfoPtr;

static ModifierInfoPtr ModifierInfoNew (void)
{
  ModifierInfoPtr mip;
  mip = (ModifierInfoPtr) MemNew (sizeof (ModifierInfoData));
  if (mip == NULL) return NULL;
  mip->name = NULL;
  mip->value = NULL;
  mip->modtype = eModifierType_SourceQual;
  return mip;
}

static ModifierInfoPtr ModifierInfoFree (ModifierInfoPtr mip)
{
  if (mip == NULL) return NULL;
  mip->name = MemFree (mip->name);
  mip->value = MemFree (mip->value);
  mip = MemFree (mip);
  return mip;
}

static ValNodePtr ModifierInfoListFree (ValNodePtr list)
{
  if (list == NULL) return NULL;
  ModifierInfoListFree (list->next);
  list->next = NULL;
  list->data.ptrvalue = ModifierInfoFree (list->data.ptrvalue);
  ValNodeFree (list);
  return NULL;
}

static EModifierType GetModifierType (CharPtr mod_name)
{
  Int4 i;
  CharPtr canonical_name;
  EModifierType returntype;
  
  canonical_name = GetCanonicalName (mod_name);
  
  if (StringHasNoText (canonical_name))
  {
    returntype = eModifierType_SourceQual;
  }
  else if (StringICmp (canonical_name, "organism") == 0
           || StringICmp (canonical_name, "org") == 0)
  {
	returntype = eModifierType_Organism;
  }
  else if (StringICmp (canonical_name, "location") == 0)
  {
    returntype = eModifierType_Location;
  }
  else if (StringICmp (canonical_name, "lineage") == 0)
  {
    returntype = eModifierType_Lineage;
  }
  else if (StringICmp (canonical_name, "gcode") == 0)
  {
    returntype = eModifierType_NucGeneticCode;
  }
  else if (StringICmp (canonical_name, "mgcode") == 0)
  {
    returntype = eModifierType_MitoGeneticCode;
  }
  else if (StringICmp (canonical_name, "genetic_code") == 0)
  {
    returntype = eModifierType_GeneticCode;
  }
  else if (StringICmp (canonical_name, "gencode_comment") == 0)
  {
    returntype = eModifierType_GeneticCodeComment;
  }
  else if (StringICmp (canonical_name, "moltype") == 0)
  {
    returntype = eModifierType_MolType;
  }
  else if (StringICmp (canonical_name, "molecule") == 0)
  {
    returntype = eModifierType_Molecule;
  }
  else if (StringICmp (canonical_name, "origin") == 0)
  {
    returntype = eModifierType_Origin;
  }
  else if (StringICmp (canonical_name, "topology") == 0)
  {
    returntype = eModifierType_Topology;
  }
  else if (StringICmp (canonical_name, "common name") == 0)
  {
    returntype = eModifierType_CommonName;
  }
  else if (StringICmp (canonical_name, "tech") == 0)
  {
    returntype = eModifierType_Technique;
  }
  else
  {
    for (i = 0; i < num_protein_modifier_names; i++)
    {
      if (StringICmp (canonical_name, protein_modifier_names[i]) == 0)
      {
        returntype = eModifierType_Protein;
        canonical_name = MemFree (canonical_name);
        return returntype;
      }
    }
    returntype = eModifierType_SourceQual;   
  }
  
  canonical_name = MemFree (canonical_name);
  return returntype;
}

static Boolean AllowMultipleValues (CharPtr mod_name)
{
  EModifierType mod_type;
  Boolean       rval = FALSE;
  
  mod_type = GetModifierType (mod_name);
  switch (mod_type)
  {
    case eModifierType_SourceQual:
      if (! IsNonTextModifier (mod_name))
      {
        rval = TRUE;
      }
      break;
    case eModifierType_CommonName:
      rval = TRUE;
      break;
    case eModifierType_Organism:
      rval = TRUE;
      break;
    default:
      rval = FALSE;
      break;
  }
  return rval;
}

typedef enum 
{
  BRACKET_ERR_NO_ERR = 0,
  BRACKET_ERR_MISMATCHED_BRACKETS,
  BRACKET_ERR_MISSING_EQUALS,
  BRACKET_ERR_MULT_EQUALS,
  BRACKET_ERR_NO_MOD_NAME,
  BRACKET_ERR_MISMATCHED_QUOTES
} bracketing_err_num;

static Char ExpectToken (CharPtr cp)
{
  CharPtr valstart;
  
  if (cp == NULL)
  {
    return 0;
  }
  else if (*cp == '[')
  {
    valstart = cp + 1 + StringSpn (cp + 1, " \t");
    if (StringLen (valstart) > 3 
        && (StringNICmp (valstart, "dna", 3) == 0 
            || StringNICmp (valstart, "rna", 3) == 0
            || StringNICmp (valstart, "orf", 3) == 0)
        && *(valstart + 3 + StringSpn (valstart + 3, " \t")) == ']')
    {
      return ']';    
    }
    else
    {
      return '=';
    }
  }
  else if (*cp == '=')
  {
    return ']';
  }
  else if (*cp == ']')
  {
    return '[';
  }
  else
  {
    return 0;
  }
}

/* When we are looking for double-quotation marks to use for delimiting
 * sections of a title that should not be parsed or values that may contain
 * brackets, equals signs, or other reserved characters, skip over
 * quotation marks that are preceded by the escape character (backslash).
 * This allows quotation marks to be included in a quoted string.
 */
static CharPtr NextUnescapedQuote (CharPtr str)
{
  CharPtr cp;
  
  if (StringHasNoText (str))
  {
    return NULL;
  }
  cp = StringChr (str, '"');
  if (cp != NULL && cp != str)
  {
    while (cp != NULL && *(cp - 1) == '\\')
    {
      cp = StringChr (cp + 1, '"');
    }
  }
  return cp;  
}

/* This function steps backward from str_end until it has located 
 * an unescaped double-quotation mark or it has reached the 
 * start of the string (str_start).
 */
static CharPtr FindPreviousUnescapedQuote (CharPtr str_start, CharPtr str_end)
{
  CharPtr cp;
  if (str_start == NULL || str_end == NULL || str_end < str_start)
  {
    return NULL;
  }
  
  cp = str_end;
  while (cp > str_start && (*cp != '"' || *(cp - 1) == '\\'))
  {
    cp--;
  }
  if (*cp != '"')
  {
    cp = NULL;
  }
  return cp;
}


/* This function finds the next bracketing token ([, =, or ]) in
 * the string that is not enclosed by unescaped quotation marks.
 */
static CharPtr NextBracketToken (CharPtr str)
{
  CharPtr next_quote;
  CharPtr cp;
  
  if (StringHasNoText (str))
  {
    return NULL;
  }
  
  cp = str;
  while (*cp != 0)
  {
    switch (*cp)
    {
      case '"':
        if (cp == str || (*(cp - 1) != '\\'))
        {
          next_quote = NextUnescapedQuote (cp + 1);
          if (next_quote == NULL)
          {
            return cp;
          }
          else
          {
            cp = next_quote + 1;;
          }
        }
        else
        {
          cp++;
        }
        break;
      case '[':
      case ']':
      case '=':
        return cp;
      default:
        cp++;
    }
  }
    
  return NULL;
}

static Int4 DetectBadBracketing (CharPtr str)
{
  CharPtr cp;
  Char    expected_token;
  CharPtr last_token = NULL, namestart;
  
  if (StringHasNoText (str))
  {
    return BRACKET_ERR_NO_ERR;
  }
  
  expected_token = '[';
  cp = NextBracketToken (str);
  while (cp != NULL)
  {
    switch (*cp)
    {
      case '"':
        return BRACKET_ERR_MISMATCHED_QUOTES;
        break;
      case '[':
      case ']':
      case '=':
        if (expected_token == *cp)
        {
          if (expected_token == '=' && last_token != NULL)
          {
            namestart = last_token + 1 + StringSpn (last_token + 1, " \t");
            if (namestart == cp)
            {
              return BRACKET_ERR_NO_MOD_NAME;
            }
          }
          expected_token = ExpectToken (cp);
          last_token = cp;
        }
        else if (expected_token == '=')
        {
          if (cp - last_token - 1 == StringSpn (last_token + 1, " \t"))
          {
            return BRACKET_ERR_MISMATCHED_BRACKETS;
          }
          else
          {
            return BRACKET_ERR_MISSING_EQUALS;
          }
        }
        else if (*cp == '=')
        {
          if (expected_token == ']')
          {
            return BRACKET_ERR_MULT_EQUALS;
          }
          else
          {
            return BRACKET_ERR_MISMATCHED_BRACKETS;
          }
        }
        else
        {
          return BRACKET_ERR_MISMATCHED_BRACKETS;
        }
        break;
    }
    cp = NextBracketToken (cp + 1);
  }
  
  if (cp == NULL && expected_token != '[')
  {
    return BRACKET_ERR_MISMATCHED_BRACKETS;
  }
  
  return BRACKET_ERR_NO_ERR;
}

static ModifierInfoPtr 
ParseOneBracketedModifier 
(CharPtr      str, 
 CharPtr PNTR bracket_start,
 CharPtr PNTR bracket_stop)
{
  CharPtr         start, stop, eq_loc;
  ModifierInfoPtr mip;
  Int4            value_len, name_len;
  CharPtr         canonical_name;
  
  start = NextBracketToken (str);
  while (start != NULL && *start != '[')
  {
    start = NextBracketToken (start + 1);
  }
  if (start == NULL) return NULL;
  eq_loc = NextBracketToken (start + 1);
  if (eq_loc == NULL) return NULL;
  if (*eq_loc == ']')
  {
    stop = eq_loc;
  }
  else if (*eq_loc == '=')
  {
    stop = NextBracketToken (eq_loc + 1);
  }
  else
  {
    return NULL;
  }
  
  if (stop == NULL || *stop != ']') return NULL;
      
  mip = ModifierInfoNew();
  if (mip == NULL) return NULL;
  
  /* copy in modifier name */
  name_len = eq_loc - start + 1;
  mip->name = (CharPtr) MemNew (name_len * sizeof (Char));
  if (mip->name == NULL)
  {
    mip = ModifierInfoFree (mip);
    return NULL;
  }
  StringNCpy (mip->name, start + 1, name_len - 2);
  mip->name [name_len - 1] = 0;
  TrimSpacesAroundString (mip->name);
  canonical_name = GetCanonicalName (mip->name);
  mip->name = MemFree (mip->name);
  mip->name = canonical_name;
  if (StringICmp (mip->name, "note") == 0)
  {
    mip->name = MemFree (mip->name);
    mip->name = StringSave ("Note-SubSrc");
  }

  /* [orf], [rna], and [dna] don't have values */  
  if (stop > eq_loc)
  {
    value_len = stop - eq_loc + 1;
    mip->value = (CharPtr) MemNew (value_len * sizeof (Char));
    if (mip->value == NULL)
    {
      mip = ModifierInfoFree (mip);
      return NULL;
    }
  
    StringNCpy (mip->value, eq_loc + 1, value_len - 2);
    mip->value [value_len - 1] = 0;
    TrimSpacesAroundString (mip->value);
  }
  
  mip->modtype = GetModifierType (mip->name);
  if (mip->modtype == eModifierType_SourceQual)
  {
    mip->subtype = FindTypeForModNameText (mip->name);
  }
  else
  {
    mip->subtype = 0;
  }
  
  if (bracket_start != NULL)
  {
    *bracket_start = start;
  }
  
  if (bracket_stop != NULL)
  {
    *bracket_stop = stop;
  }
  
  return mip;
}

static ValNodePtr ParseAllBracketedModifiers (CharPtr str)
{
  CharPtr         stop, cp;
  ValNodePtr      list = NULL;
  ModifierInfoPtr mip;
  
  cp = str;
  mip = ParseOneBracketedModifier (cp, NULL, &stop);
  while (mip != NULL && stop != NULL)
  {
    ValNodeAddPointer (&list, 0, mip);
    cp = stop + 1;
    mip = ParseOneBracketedModifier (cp, NULL, &stop);  
  }
  return list;
}


NLM_EXTERN CharPtr GetPresentModifierNames (CharPtr str)
{
  ValNodePtr      list = NULL, vnp;
  Int4            len = 1;
  CharPtr         text = NULL;
  ModifierInfoPtr mip;
  
  list = ParseAllBracketedModifiers(str);
  if (list == NULL) {
    return StringSave ("");
  }
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    mip = (ModifierInfoPtr)vnp->data.ptrvalue;
    len += StringLen (mip->name) + 2;
  }

  text = (CharPtr) MemNew (sizeof (Char) * len);
  mip = (ModifierInfoPtr)list->data.ptrvalue;
  StringCpy (text, mip->name);
  for (vnp = list->next; vnp != NULL; vnp = vnp->next) {
    mip = (ModifierInfoPtr)vnp->data.ptrvalue;
    StringCat (text, ", ");
    StringCat (text, mip->name);
  }

  list = ModifierInfoListFree(list);

  return text;
}


static Boolean IsValueInEnumAssoc (CharPtr value, EnumFieldAssocPtr eap)
{
  while (eap != NULL && eap->name != NULL) 
  {
    if (StringICmp (eap->name, value) == 0)
    {
      return TRUE;
    }
    eap++;
  }
  return FALSE;
}

static Int4 GeneticCodeFromStringAndList (CharPtr str, ValNodePtr list)
{
  while (list != NULL)
  {
    if (StringICmp (str, list->data.ptrvalue) == 0)
    {
      return list->choice;
    }
    list = list->next;
  }
  return 0;
}

static Int4 GeneticCodeFromString (CharPtr str)
{
  ValNodePtr gencodelist;
  Int4       gcode = 0;
  
  if (StringHasNoText (str))
  {
    gcode = 0;
  }
  else if (isdigit (str[0]))
  {
    gcode = atoi (str);
  }
  else
  {
    gencodelist = GetGeneticCodeValNodeList ();
    gcode = GeneticCodeFromStringAndList (str, gencodelist);
    gencodelist = ValNodeFreeData (gencodelist);
  }
  return gcode;
}

static Int4 MolTypeFromString (CharPtr str)
{
  EnumFieldAssocPtr  eap;

  if (StringICmp (str, "dna") == 0)
  {
    return 253;
  }
  else if (StringICmp (str, "rna") == 0)
  {
    return 254;
  }
  else if (StringICmp (str, "genomic") == 0)
  {
    return 253;
  }
  for (eap = biomol_nucGen_alist; eap != NULL && eap->name != NULL; eap++)
  {
    if (StringsAreEquivalent (eap->name, str))
    {
      return eap->value;
    }
  }
  for (eap = biomol_nucX_alist; eap != NULL && eap->name != NULL; eap++)
  {
    if (StringsAreEquivalent (eap->name, str))
    {
      return eap->value;
    }
    else if (eap->name [0] == 'm'
             && StringICmp (eap->name, "mRNA [cDNA]") == 0
             && StringICmp (str, "mRNA") == 0)
    {
      return eap->value;
    }
  }
  return 0;
}


/* This function looks at a parsed modifier structure to determine whether the
 * value is acceptable for this modifier type.
 */
static Boolean ModifierHasInvalidValue (ModifierInfoPtr mip)
{
  Boolean rval = FALSE;
  
  if (mip != NULL
      && ((mip->modtype == eModifierType_Location
  	          && !IsValueInEnumAssoc (mip->value, biosource_genome_simple_alist))
  	    || (mip->modtype == eModifierType_Origin
  	          && !IsValueInEnumAssoc (mip->value, biosource_origin_alist))      
  	    || (mip->modtype == eModifierType_Topology
  	          && !IsValueInEnumAssoc (mip->value, topology_nuc_alist))
  	    || (mip->modtype == eModifierType_Molecule
  	          && !IsValueInEnumAssoc (mip->value, molecule_alist))
  	    || ((mip->modtype == eModifierType_GeneticCode
  	              || mip->modtype == eModifierType_NucGeneticCode
  	              || mip->modtype == eModifierType_MitoGeneticCode)
  	             && GeneticCodeFromString (mip->value) == 0)
  	    || (mip->modtype == eModifierType_MolType
  	             && MolTypeFromString (mip->value) == 0)
  	    || (mip->modtype == eModifierType_SourceQual
  	             && IsNonTextModifier (mip->name)
  	             && !StringHasNoText (mip->value)
                 && StringICmp (mip->value, "TRUE") != 0
                 && StringICmp (mip->value, "FALSE") != 0)))
  {
    rval = TRUE;
  }

  return rval;
}

/* This section contains functions for finding, changing, and removing 
 * bracketed value pairs in definition lines.
 * These functions include:
 *
 * FindValuePairInDefLine - returns pointer to position in title where 
 *                          the first bracketed pair with the specified 
 *                          modifier name (or one of its aliases) occurs.
 *                          Useful for non-text modifiers, which do not
 *                          have values.
 *
 * FindValueFromPairInDefline - returns value from the first bracketed
 *                              pair in the title with the specified
 *                              modifier name (or one of its aliases).
 *
 * RemoveValueFromDefline - removes the first bracketed pair in the title
 *                          with the specified modifier name (or one of its aliases)
 *
 * ReplaceValueInThisValuePair - replaces the value in the specified value pair.
 *                               if new value is empty, pair is removed.
 *
 * ReplaceValueInOneDefLine - finds the first bracketed pair in the title
 *                            with the specified modifier name (or one of its aliases).
 *                            If a pair is found, the value in that pair is replaced
 *                            with the new value; otherwise a new pair is added to
 *                            the title.  
 *
 * ReplaceOneModifierValue - finds all bracketed pairs in a title with the specified
 *                           modifier name or one of its aliases and the specified value
 *                           and replaces that value with the new value (or removes the
 *                           pair, if the new value is empty.
 *
 * RemoveAllDuplicatePairsFromOneTitle - removes all bracketed pairs that are duplicates
 *                                       in name and value of another pair already in
 *                                       the title.
 *
 * RemoveMeaninglessEmptyPairsFromOneTitle - removes bracketed pairs without values
 *                                           that are not non-text modifiers
 *
 * StripAllInstancesOfModNameFromTitle - removes all mentions of specified modifier
 *                                       name from title
 *
 */

static CharPtr FindValuePairInDefLine (CharPtr mod_name, CharPtr def_line, CharPtr PNTR valstop)
{
  CharPtr         cp, start, stop;
  ModifierInfoPtr mip;
  CharPtr         canonical_name;
  
  if (mod_name == NULL || def_line == NULL)
  {
    return NULL;
  }
  
  cp = NextBracketToken (def_line);
  if (cp == NULL)
  {
    return NULL;
  }
  
  canonical_name = GetCanonicalName (mod_name);
  
  mip = ParseOneBracketedModifier (cp, &start, &stop);
  while (mip != NULL && start != NULL && stop != NULL 
         && StringICmp (mip->name, canonical_name) != 0)
  {
    cp = NextBracketToken (stop + 1);
    mip = ModifierInfoFree (mip);
    mip = ParseOneBracketedModifier (cp, &start, &stop);
  }
  
  if (mip != NULL && StringICmp (mip->name, canonical_name) == 0)
  {
    mip = ModifierInfoFree (mip);
    if (valstop != NULL)
    {
      *valstop = stop;
    }
    canonical_name = MemFree (canonical_name);
    return start;
  }
  else
  {
    mip = ModifierInfoFree (mip);
    canonical_name = MemFree (canonical_name);
    return NULL;
  }
}

static CharPtr FindNthValuePairInDefLine (CharPtr title, CharPtr val_name, Int4 val_num, CharPtr PNTR p_val_end)
{
  CharPtr val_loc, val_end = NULL;
  Int4    title_val_num;
  
  if (StringHasNoText (val_name))
  {
    return NULL;
  }
  
  val_loc = FindValuePairInDefLine (val_name, title, &val_end);
  title_val_num = 0;
  while (val_loc != NULL && val_end != NULL && title_val_num != val_num)
  {
    val_loc = FindValuePairInDefLine (val_name, val_end + 1, &val_end);
    title_val_num++;
  }
  if (p_val_end != NULL)
  {
    *p_val_end = val_end;
  }
  return val_loc;
}

NLM_EXTERN CharPtr FindValueFromPairInDefline (CharPtr mod_name, CharPtr def_line)
{
  CharPtr bracket_start, eq_loc, bracket_end;
  CharPtr new_val = NULL;
  Int4 new_val_len;
  
  bracket_start = FindValuePairInDefLine (mod_name, def_line, &bracket_end);
  if (bracket_start == NULL || bracket_end == NULL)
  {
    return NULL;
  }
  
  eq_loc = NextBracketToken (bracket_start + 1);
  if (eq_loc == NULL || *eq_loc != '=')
  {
    return NULL;
  }
    
  new_val_len = bracket_end - eq_loc;
  new_val = (CharPtr) MemNew (new_val_len * sizeof (Char));
  if (new_val != NULL)
  {
    StringNCpy (new_val, eq_loc + 1, new_val_len - 1);
    new_val [new_val_len - 1] = 0;
  }
  TrimSpacesAroundString (new_val);
  return new_val;
}

static CharPtr FindValueFromPairInDeflineBeforeCharPtr (CharPtr mod_name, CharPtr def_line, CharPtr cp)
{
  CharPtr bracket_start, bracket_end;
  
  bracket_start = FindValuePairInDefLine (mod_name, def_line, &bracket_end);
  if (bracket_start == NULL || (cp != NULL && bracket_start > cp))
  {
    return NULL;
  }
  else
  {
    return FindValueFromPairInDefline (mod_name, bracket_start);
  }
}

static void RemoveValuePairFromDefline (CharPtr pair_start, CharPtr pair_end, CharPtr defline)
{
  CharPtr src, dst;

  if (pair_start == NULL || pair_end == NULL || defline == NULL
      || pair_end <= pair_start)
  {
    return;
  }
  
  dst = pair_start;
  src = pair_end;
  while (isspace (*src))
  {
    src++;
  }
  
  while (*src != 0)
  {
    *dst = *src;
    dst++;
    src++;
  }
  *dst = 0;
}

NLM_EXTERN void RemoveValueFromDefline (CharPtr mod_name, CharPtr def_line)
{
  CharPtr bracket_start, bracket_end;
  
  bracket_start = FindValuePairInDefLine (mod_name, def_line, &bracket_end);
  if (bracket_start == NULL || bracket_end == NULL)
  {
    return;
  }
  
  RemoveValuePairFromDefline (bracket_start, bracket_end + 1, def_line);
}

static CharPtr AddQuotesToValueWithBrackets (CharPtr orig_value)
{
  CharPtr first_bracket, first_quote;
  CharPtr cp, new_value = NULL, tmp_value;
  Char    bracket_buf [2];
  Int4    offset;
  
  if (orig_value == NULL)
  {
    return NULL;
  }
  else if (StringHasNoText (orig_value))
  {
    return StringSave (orig_value);
  }
  
  new_value = StringSave (orig_value);
  
  first_bracket = StringChr (new_value, '[');
  if (first_bracket == NULL)
  {
    first_bracket = StringChr (new_value, ']');
  }
  
  first_quote = NextUnescapedQuote (new_value);
  
  if (first_bracket == NULL && first_quote == NULL)
  {
    return new_value;
  }
  else if (first_bracket != NULL && first_quote == NULL)
  {
    tmp_value = (CharPtr) MemNew ((StringLen (new_value) + 3) * sizeof (Char));
    if (tmp_value == NULL)
    {
      new_value = MemFree (new_value);
      return NULL;
    }
    StringCat (tmp_value, "\"");
    StringCat (tmp_value, new_value);
    StringCat (tmp_value, "\"");
    new_value = MemFree (new_value);
    new_value = tmp_value;
    return new_value;
  }
  
  cp = orig_value;
  
  bracket_buf [0] = 0;
  bracket_buf [1] = 0;
  
  while (*cp != 0)
  {
    if (*cp == '"' && (cp == orig_value || *(cp - 1) != '\\'))
    {
      cp = NextUnescapedQuote (cp + 1);
      if (cp == NULL)
      {
        tmp_value = (CharPtr) MemNew ((StringLen (new_value) + 3) * sizeof (Char));
        if (tmp_value == NULL)
        {
          new_value = MemFree (new_value);
          return NULL;
        }
        StringCpy (tmp_value, new_value);
        if (new_value [StringLen (new_value) - 1] == '\\')
        {
          StringCat (tmp_value, " ");
        }
        StringCat (tmp_value, "\"");
        return tmp_value;
      }
      else
      {
        cp++;
      }
    }
    else if (*cp == '[' || *cp == ']')
    {
      tmp_value = (CharPtr) MemNew ((StringLen (new_value) + 3) * sizeof (Char));
      if (tmp_value == NULL)
      {
        new_value = MemFree (new_value);
        return new_value;
      }
      offset = cp - new_value;
      StringNCpy (tmp_value, new_value, offset);
      StringCat (tmp_value, "\"");
      bracket_buf [0] = *cp;
      StringCat (tmp_value, bracket_buf);
      StringCat (tmp_value, "\"");
      StringCat (tmp_value, cp + 1);
      new_value = MemFree (new_value);
      new_value = tmp_value;
      cp = new_value + offset + 3;      
    }
    else
    {
      cp++;
    }
  }
  
  return new_value;
}

static CharPtr
ReplaceValueInThisValuePair 
(CharPtr orig_defline,
 CharPtr value_loc, 
 CharPtr value_name,
 CharPtr end_loc, 
 CharPtr new_value)
{
  CharPtr new_title;
  Int4    new_title_len = 0;
  Boolean is_nontext;
  CharPtr tmp_name;
  CharPtr fixed_value;

  if (StringHasNoText (orig_defline) || value_loc == NULL || end_loc == NULL
      || *value_loc != '[' || *end_loc != ']')
  {
    return orig_defline;
  }

  fixed_value = AddQuotesToValueWithBrackets (new_value);
  
  if (StringHasNoText (fixed_value))
  {
    RemoveValuePairFromDefline (value_loc, end_loc, orig_defline);
  }
  else
  {
    /* keep part before pair and after pair, insert new value in position */
    new_title_len = StringLen (orig_defline)
                               + StringLen (value_name)
                               + StringLen (fixed_value)
                               + 5;
    new_title = MemNew (new_title_len * sizeof (Char));
    if (new_title != NULL)
    {
      if (value_loc > orig_defline)
      {
        StringNCpy (new_title, orig_defline, value_loc - orig_defline);
      }
      StringCat (new_title, "[");
      tmp_name = StringSave (value_name);
      tmp_name [0] = TO_LOWER (tmp_name [0]);
      StringCat (new_title, tmp_name);
      is_nontext = IsNonTextModifier (tmp_name);
      tmp_name = MemFree (tmp_name);
      StringCat (new_title, "=");
      if (!is_nontext)
      {
        StringCat (new_title, fixed_value);
      }
      StringCat (new_title, "]");
      if (end_loc != NULL && *end_loc != 0)
      {
        if (*end_loc == ']')
        {
          StringCat (new_title, end_loc + 1);
        }
        else
        {
          StringCat (new_title, end_loc);
        }
      }
      orig_defline = MemFree (orig_defline);
      orig_defline = new_title;
    }
  }  
  TrimSpacesAroundString (orig_defline);
  
  fixed_value = MemFree (fixed_value);
  
  return orig_defline;
}

static CharPtr InsertStringAtOffset (CharPtr old_string, CharPtr new_string, Int4 offset)
{
  Int4    new_len;
  CharPtr new_str = NULL;
  
  if (old_string == NULL)
  {
    new_str = StringSave (new_string);
  }
  else if (new_string == NULL)
  {
    new_str =  StringSave (old_string);
  }
  else
  {
    new_len = StringLen (old_string) + StringLen (new_string) + 1;
    new_str = (CharPtr) MemNew (new_len * sizeof (Char));
    if (new_str != NULL)
    {
      StringNCpy (new_str, old_string, offset);
      StringCat (new_str, new_string);
      if ((Uint4)offset < StringLen (old_string))
      {
        StringCat (new_str, old_string + offset);
      }
    }
  }
  return new_str;
}

static CharPtr 
InsertValuePairAtOffset 
(CharPtr orig_defline, 
 CharPtr value_name, 
 CharPtr value_str,
 Int4    offset)
{
  CharPtr pair_string, fixed_value;
  
  if (StringHasNoText (value_name) || offset < 0)
  {
    return orig_defline;
  }
  
  fixed_value = AddQuotesToValueWithBrackets (value_str);
  
  pair_string = (CharPtr) MemNew ((StringLen (value_name) + StringLen (fixed_value) + 6) * sizeof (Char));
  if (pair_string != NULL)
  {
    if (IsNonTextModifier (value_name))
    {
      sprintf (pair_string, "[%s=]", value_name);
    }
    else
    {
      sprintf (pair_string, "[%s=%s]", value_name, fixed_value);
    }
    orig_defline = InsertStringAtOffset (orig_defline, pair_string, offset);
    pair_string = MemFree (pair_string);
  }
  fixed_value = MemFree (fixed_value);
  return orig_defline;
}


static CharPtr
ReplaceValueInOneDefLineForOrganism
(CharPtr orig_defline,
 CharPtr value_name, 
 CharPtr new_value,
 CharPtr organism)
{ 
  CharPtr value_loc = NULL, end_loc = NULL;
  CharPtr fixed_value;
  CharPtr next_org_loc = NULL, org_stop = NULL, first_org_stop = NULL;
  CharPtr first_organism;
  
  if (StringHasNoText (value_name))
  {
    return orig_defline;
  }
  
  /* if we want to add a value to a specific organism, we need to make sure
   * that we insert or replace a value after that organism name but before
   * the next organism name.
   */
   
  if (organism != NULL)
  {
    if (organism < orig_defline || organism - orig_defline > (Int4) StringLen (orig_defline))
    {
      organism = NULL;
    }
  }
  
  if (organism != NULL)
  {
    if (organism != FindValuePairInDefLine ("organism", organism, &org_stop))
    {
      return orig_defline;
    }
  }
  
  first_organism = FindValuePairInDefLine ("organism", orig_defline, &first_org_stop);

  
  if (organism == NULL)
  {
    organism = first_organism;
    org_stop = first_org_stop;
  }
  
  if (org_stop != NULL)
  {
    next_org_loc = FindValuePairInDefLine ("organism", org_stop + 1, NULL);
  }
  
  fixed_value = AddQuotesToValueWithBrackets (new_value);
  
  /* if this is the first organism, or if we have no organism, start looking for
   * a value to replace at the beginning of the line.
   */
  if (organism == NULL || organism == first_organism)
  {
    value_loc = FindValuePairInDefLine (value_name, orig_defline, &end_loc);
  }
  else
  {     
    value_loc = FindValuePairInDefLine (value_name, organism, &end_loc);
  }
  
  if (next_org_loc != NULL && value_loc > next_org_loc)
  {
    value_loc = NULL;
  }
  
  if (StringHasNoText (fixed_value))
  {
    if (value_loc == NULL)
    {
      /* old line had no value, no new value provided, no change */
    }
    else
    {
      RemoveValuePairFromDefline (value_loc, end_loc, orig_defline);
    }
  }
  else
  {
    if (value_loc == NULL)
    {
      /* add new value just before next organism */
      if (next_org_loc == NULL)
      {
        orig_defline = InsertValuePairAtOffset (orig_defline, value_name, new_value,
                                                StringLen (orig_defline));
      }
      else
      {
        orig_defline = InsertValuePairAtOffset (orig_defline, value_name, new_value,
                                                next_org_loc - orig_defline);
      }
    }
    else
    {
      /* replace this value */
      orig_defline = ReplaceValueInThisValuePair (orig_defline, value_loc, value_name,
                                                  end_loc, new_value);
    }
  }  
  TrimSpacesAroundString (orig_defline);
  
  fixed_value = MemFree (fixed_value);
  
  return orig_defline;
}

NLM_EXTERN CharPtr 
ReplaceValueInOneDefLine 
(CharPtr orig_defline,
 CharPtr value_name, 
 CharPtr new_value)
{
  CharPtr value_loc = NULL, end_loc = NULL;
  
  if (StringHasNoText (value_name))
  {
    return orig_defline;
  }
  
  value_loc = FindValuePairInDefLine (value_name, orig_defline, &end_loc);

  if (value_loc == NULL)
  {
    if (StringHasNoText (new_value))
    {
      /* old line had no value, no new value provided, no change */    
      return orig_defline;
    }
    else
    {
      /* make sure value is added for first organism */
      orig_defline = ReplaceValueInOneDefLineForOrganism (orig_defline, value_name,
                                                          new_value, NULL);        
    }
  }
  else
  {
    orig_defline = ReplaceValueInThisValuePair (orig_defline, value_loc, value_name, end_loc, new_value);
  }  
  
  return orig_defline;
}

static CharPtr 
ReplaceOneModifierValue 
(CharPtr title,
 CharPtr orig_name, 
 CharPtr orig_value,
 CharPtr repl_value,
 Boolean is_nontext,
 Boolean copy_to_note)
{
  CharPtr bracket_loc, eq_loc, end_bracket_loc, new_title;
  Int4    new_title_len;
  CharPtr orig_note, new_note;
  Boolean any_replaced = FALSE;
  
  if (StringHasNoText (title)
      || StringHasNoText (orig_name))
  {
    return title;
  }
  
  bracket_loc = FindValuePairInDefLine (orig_name, title, &end_bracket_loc);
  while (bracket_loc != NULL && end_bracket_loc != NULL)
  {  
    eq_loc = NextBracketToken (bracket_loc + 1);
    if (eq_loc == NULL || *eq_loc != '=')
    {
      return title;
    }
    if ((StringNCmp (orig_value, eq_loc + 1, StringLen (orig_value)) == 0
        && StringLen (orig_value) == end_bracket_loc - eq_loc - 1)
        || (StringHasNoText (orig_value) 
            && StringSpn (eq_loc + 1, " \t") == end_bracket_loc - eq_loc - 1))
    {
      new_title_len = StringLen (title) + StringLen (repl_value) - StringLen (orig_value) + 1;
      new_title = (CharPtr) MemNew (new_title_len * sizeof (Char));
      if (new_title == NULL)
      {
        return title;
      }
      if (is_nontext)
      {
        if (StringHasNoText (repl_value))
        {
          StringNCpy (new_title, title, bracket_loc - title);
          StringCat (new_title, end_bracket_loc + 1 + StringSpn (end_bracket_loc, " "));
        }
        else
        {
          StringNCpy (new_title, title, eq_loc - title + 1);
          StringCat (new_title, end_bracket_loc);
        }
      }
      else if (StringHasNoText (repl_value))
      {
        /* remove pair completely */
        StringNCpy (new_title, title, bracket_loc - title);
        StringCat (new_title, end_bracket_loc + 1);
      }
      else
      {
        StringNCpy (new_title, title, eq_loc - title + 1);
        StringCat (new_title, repl_value);
        StringCat (new_title, end_bracket_loc);        
      }

      title = MemFree (title);
      title = new_title;
      any_replaced = TRUE;
      bracket_loc = FindValuePairInDefLine (orig_name, title, &end_bracket_loc);
    }
    else
    {
      bracket_loc = FindValuePairInDefLine (orig_name, end_bracket_loc, &end_bracket_loc);
    }
  }
  
  if (any_replaced && copy_to_note && !StringHasNoText (repl_value) && !StringHasNoText (orig_value))
  {
    orig_note = FindValueFromPairInDefline ("note", title);
    if (StringHasNoText (orig_note))
    {
      new_note = (CharPtr) MemNew ((StringLen (orig_name) 
                                    + StringLen (orig_value) + 8) * sizeof (Char));
      if (new_note != NULL)
      {
        sprintf (new_note, "%s was %s", orig_name, orig_value);
      }
    }
    else
    {
      new_note = (CharPtr) MemNew ((StringLen (orig_note)
                                    + StringLen (orig_name) 
                                    + StringLen (orig_value) + 8) * sizeof (Char));
      if (new_note != NULL)
      {
        sprintf (new_note, "%s; %s was %s", orig_note, orig_name, orig_value);
      }
    }

    if (new_note != NULL)
    {
      title = ReplaceValueInOneDefLine (title, "note", new_note); 
    }
    
    orig_note = MemFree (orig_note);
    new_note = MemFree (new_note);
  }
  
  return title;
}

static Boolean IsUnrecognizedModifierName (ModifierInfoPtr mip, Boolean is_nuc);

static void RemoveRecognizedModifiersFromTitle (CharPtr title, ValNodePtr modifier_info_list, Boolean is_nuc)
{
  ValNodePtr      vnp;
  ModifierInfoPtr mip;

  for (vnp = modifier_info_list; vnp != NULL; vnp = vnp->next) {
    mip = (ModifierInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL && ! IsUnrecognizedModifierName (mip, is_nuc)
        && (!is_nuc || mip->modtype != eModifierType_Protein)) {
      RemoveValueFromDefline (mip->name, title);
    }
  }
}

static void StripAllInstancesOfModNameFromTitle (CharPtr mod_name, CharPtr title)
{
  CharPtr         valstr;

  valstr = FindValueFromPairInDefline (mod_name, title);
  while (valstr != NULL)
  {
    RemoveValueFromDefline (mod_name, title);
    valstr = MemFree (valstr);
    valstr = FindValueFromPairInDefline (mod_name, title);
  }	
}

static CharPtr RemoveAllDuplicatePairsFromOneTitle (CharPtr title)
{
  CharPtr         start_bracket, end_bracket, tmp_title, new_title;
  ModifierInfoPtr mip;
  Int4            offset;

  mip = ParseOneBracketedModifier (title, &start_bracket, &end_bracket);
  while (mip != NULL && start_bracket != NULL && end_bracket != NULL)
  {
    offset = end_bracket - title + 1;
    tmp_title = StringSave (title + offset);
    tmp_title = ReplaceOneModifierValue (tmp_title, mip->name, mip->value, NULL,
                                     IsNonTextModifier (mip->name), FALSE);
    new_title = (CharPtr) MemNew ((StringLen (tmp_title) + offset + 1)* sizeof (Char));
    if (new_title != NULL)
    {
      StringNCpy (new_title, title, offset);
      StringCat (new_title, tmp_title);
    }
    tmp_title = MemFree (tmp_title);
    title = MemFree (title);
    title = new_title;
    mip = ModifierInfoFree (mip);
    mip = ParseOneBracketedModifier (title + offset, &start_bracket, &end_bracket);
  }
  mip = ModifierInfoFree (mip);
  return title;
}

static void ShiftString (CharPtr str, Int4 shift_size)
{
  CharPtr src, dst;
  
  if (str == NULL)
  {
    return;
  }
  
  if (shift_size > (Int4) StringLen (str))
  {
    *str = 0;
  }
  else
  {
    src = str + shift_size;
    dst = str;
    while (*src != 0)
    {
      *dst = *src;
      dst++;
      src++;
    }
    *dst = 0;
  }  
}

static void RemoveMeaninglessEmptyPairsFromOneTitle (CharPtr title)
{
  CharPtr         start_bracket, end_bracket;
  ModifierInfoPtr mip;

  mip = ParseOneBracketedModifier (title, &start_bracket, &end_bracket);
  while (mip != NULL && start_bracket != NULL && end_bracket != NULL)
  {
    if (StringHasNoText (mip->value) && ! IsNonTextModifier (mip->name))
    {
      ShiftString (start_bracket, end_bracket - start_bracket + 1);
      mip = ModifierInfoFree (mip);
      mip = ParseOneBracketedModifier (start_bracket, &start_bracket, &end_bracket);
    }
    else
    {
      mip = ModifierInfoFree (mip);
      mip = ParseOneBracketedModifier (end_bracket + 1, &start_bracket, &end_bracket);
    }
  }
  mip = ModifierInfoFree (mip);
}

static void ApplyOneModToSeqEntry (SeqEntryPtr sep, CharPtr mod_name, CharPtr mod_value)
{
  BioseqPtr    bsp = NULL;
  SeqDescrPtr  sdp = NULL;
  
  if (sep == NULL || StringHasNoText (mod_name))
  {
    return;
  }
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  }
  else if (IS_Bioseq_set (sep))
  {
    sep = FindNucSeqEntry (sep);
    if (sep != NULL && IS_Bioseq (sep))
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
    }
  }
  
  if (bsp == NULL)
  {
    return;
  }
  
  for (sdp = bsp->descr; sdp != NULL && sdp->choice != Seq_descr_title; sdp = sdp->next)
  {
  }

  if (sdp == NULL)
  {
    sdp = SeqDescrNew (NULL);
    sdp->choice = Seq_descr_title;
    if (bsp->descr == NULL)
    {
      bsp->descr = sdp;
    }
  }
  if (sdp != NULL)
  {
    sdp->data.ptrvalue = ReplaceValueInOneDefLine (sdp->data.ptrvalue,
                                                   mod_name, mod_value);
  }
  
  
}

static ModifierInfoPtr MakeModifierInfoFromNameAndValue (CharPtr value_name, CharPtr value_string)
{
  ModifierInfoPtr mip;
  CharPtr tmp_pair;
  CharPtr fixed_value;

  fixed_value = AddQuotesToValueWithBrackets (value_string);  
  tmp_pair = (CharPtr) MemNew ((StringLen (value_name) + StringLen (fixed_value) + 4));
  if (tmp_pair == NULL)
  {
    return NULL;
  }
  sprintf (tmp_pair, "[%s=%s]", value_name == NULL ? "" : value_name,
                             fixed_value == NULL ? "" : fixed_value);
  mip = ParseOneBracketedModifier (tmp_pair, NULL, NULL);
  tmp_pair = MemFree (tmp_pair);
  fixed_value = MemFree (fixed_value);
  return mip;
}

/* This section is used to import tables of modifiers. */
static CharPtr 
ApplyImportModToTitle 
(CharPtr title, 
 CharPtr value_name, 
 CharPtr value_string,
 Boolean erase_where_blank,
 Boolean parse_multiple)
{
  ModifierInfoPtr mip;
  CharPtr next_semi, val_start, title_loc, title_end;
  CharPtr insert_point;
  Int4    insert_offset, title_val_num;
  Char    val_save_ch;

  if (StringHasNoText (value_name))
  {
    return title;
  }
  
  if (!erase_where_blank && StringHasNoText (value_string))
  {
    return title;
  }
  
  mip = MakeModifierInfoFromNameAndValue (value_name, value_string);

  if (mip == NULL 
      || (mip->modtype == eModifierType_SourceQual
        	&& mip->subtype == 255
  	      && StringICmp (mip->name, "note-subsrc") != 0 
  	      && StringICmp (mip->name, "note-orgmod") != 0))
  {
    mip = ModifierInfoFree (mip);
    return title;
  }
  
  if (erase_where_blank && StringHasNoText (value_string))
  {
    RemoveValueFromDefline (value_name, title);
  }
  else if (parse_multiple
           && value_string [0] == '(' 
           && value_string [StringLen (value_string) - 1] == ')'
           && (next_semi = StringChr (value_string, ';')) != NULL)
  {
    val_start = value_string + 1;
    title_val_num = 0;
    while (next_semi != NULL)
    {
      /* temporarily truncate at end of value */
      val_save_ch = *next_semi;
      *next_semi = 0;
      
      title_loc = FindNthValuePairInDefLine (title, value_name, title_val_num, &title_end);
      if (StringHasNoText (val_start))
      {
        if (title_loc != NULL)
        {
          RemoveValuePairFromDefline (title_loc, title_end, title);
        }
        else
        {
          /* if text is empty and there is no value pair, nothing to do */
        }
        /* note - we do not increment title_val_num here because either we've
         * removed a value or there are no values left.
         */
      }
      else
      {
        if (title_loc == NULL)
        {
          /* need to insert a new value - if organism name, put at end of title,
           * otherwise insert before second organism name if any
           */
          if (StringICmp (value_name, "organism") == 0)
          {
            insert_offset = StringLen (title);
          }
          else
          {
            insert_point = FindNthValuePairInDefLine (title, "organism", 1, NULL);
            if (insert_point == NULL) 
            {
              insert_offset = StringLen (title);
            }
            else
            {
              insert_offset = insert_point - title;
            }
          }
          title = InsertValuePairAtOffset (title, value_name, val_start, insert_offset);
        }
        else
        {
          /* replace values in order */
          title = ReplaceValueInThisValuePair (title, title_loc, value_name, 
                                               title_end, val_start);
        }
        
        title_val_num++;
      }
      
      /* replace character */
      *next_semi = val_save_ch;
      /* advance to next value in list */
      val_start = next_semi + 1;      
      if (*next_semi == ';')
      {
        next_semi = StringChr (next_semi + 1, ';');
        if (next_semi == NULL)
        {
          next_semi = value_string + StringLen (value_string) - 1;
        }
      }
      else
      {
        next_semi = NULL;
      }
    }
  }
  else if (StringCmp (value_name, "organism") == 0)
  {
    title = ReplaceValueInOneDefLine (title, value_name, value_string);
  }
  else
  {
    title = ReplaceValueInOneDefLineForOrganism (title, value_name, value_string, NULL);
  }
 
  mip = ModifierInfoFree (mip);
  return title;
}

static ValNodePtr ReadRowListFromFile (void)
{
  Char          path [PATH_MAX];
  ValNodePtr    header_line = NULL;
  FILE           *fp;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return NULL;
  fp = FileOpen (path, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", path);
  } else {
    header_line = ReadTabTableFromFile (fp);
    FileClose (fp);
  }
  return header_line;
}

/* This function will find the sequence number in the IDAndTitleEdit
 * to use for each row and put that value in the sequence_numbers array.
 */
static Boolean 
ValidateModifierTableSequenceIDs 
(ValNodePtr        header_line,
 IDAndTitleEditPtr iatep,
 Int4Ptr           sequence_numbers,
 Int4Ptr           num_rows)
{
  ValNodePtr   not_found = NULL;
  ValNodePtr   found_more_than_once = NULL;
  CharPtr      too_many_msg = NULL, not_found_msg = NULL;
  Boolean      rval = TRUE;
  Int4         msg_len = 0;
  CharPtr      too_many_fmt = " found more than once\n";
  CharPtr      not_found_fmt = " not found\n";
  CharPtr      err_msg = NULL;
  ValNodePtr   row_vnp, col_vnp, prev_row, next_row;
  Int4         i, seq_num, other_instances;
  Boolean      found;
  Int4         row_number;
  Int4         deleted_rows;
  
  if (header_line == NULL || header_line->next == NULL || iatep == NULL
      || sequence_numbers == NULL || num_rows == NULL || *num_rows < ValNodeLen (header_line->next))
  {
    return FALSE;
  }
  
  for (row_vnp = header_line->next, row_number = 0; 
       row_vnp != NULL && row_number < *num_rows;
       row_vnp = row_vnp->next, row_number++)
  {
    col_vnp = row_vnp->data.ptrvalue;
    if (col_vnp == NULL || col_vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    
    /* find correct sequence number */
    seq_num = -1;
    for (i = 0, found = FALSE; i < iatep->num_sequences && !found; i++)
    {
      if (StringCmp (iatep->id_list [i], col_vnp->data.ptrvalue) == 0)
      {
        seq_num = i;
        found = TRUE;
      }
    }
    sequence_numbers[row_number] = seq_num;
    
    if (!found)
    {
      ValNodeAddPointer (&not_found, 0, StringSave (col_vnp->data.ptrvalue));
    }
    else 
    {
      /* count the number of times this seq_num has already appeared in the list.*/
      other_instances = 0;
      for (prev_row = header_line->next; prev_row != row_vnp; prev_row = prev_row->next)
      {
        if (prev_row->choice == seq_num)
        {
          other_instances++;
        }
      }
      /* if the value was found exactly once, add this to the list of duplicates.
       * if the value was found more than once, it will already have been reported.
       */
      if (other_instances == 1)
      {
        ValNodeAddPointer (&found_more_than_once, 0, StringSave (col_vnp->data.ptrvalue));
      }
    }
  }
  
  if (found_more_than_once != NULL || not_found != NULL)
  {
    if (found_more_than_once != NULL)
    {
      too_many_msg = CreateListMessage ("Sequence ID", NULL, found_more_than_once);
      rval = FALSE;
      msg_len += StringLen (too_many_msg) + StringLen (too_many_fmt) + 5;
    }
    if (not_found != NULL)
    {
      not_found_msg = CreateListMessage ("Sequence ID", NULL, not_found);
      msg_len += StringLen (not_found_msg) + StringLen (not_found_fmt) + 5;
    }
    
    err_msg = (CharPtr) MemNew ((msg_len + 1) * sizeof (Char));
    if (err_msg != NULL)
    {
      if (too_many_msg != NULL)
      {
        StringCat (err_msg, too_many_msg);
        if (found_more_than_once->next != NULL)
        {
          StringCat (err_msg, " were");
        }
        else
        {
          StringCat (err_msg, " was");
        }
        StringCat (err_msg, too_many_fmt);
      }
      if (not_found_msg != NULL)
      {
        StringCat (err_msg, not_found_msg);
        if (not_found->next != NULL)
        {
          StringCat (err_msg, " were");
        }
        else
        {
          StringCat (err_msg, " was");
        }
        StringCat (err_msg, not_found_fmt);
      }
      if (rval)
      {
        if (ANS_NO == Message (MSG_YN, "%sContinue anyway?", err_msg))
        {
          rval = FALSE;
        }
      }
      else
      {
        Message (MSG_ERROR, "%sPlease correct your file.", err_msg);
      }
    }
    too_many_msg = MemFree (too_many_msg);
    not_found_msg = MemFree (not_found_msg);
    err_msg = MemFree (err_msg);
  }
  
  /* remove rows for sequence IDs that are not found */
  for (row_vnp = header_line->next, row_number = 0, prev_row = header_line; 
       row_vnp != NULL && row_number < *num_rows;
       row_vnp = next_row, row_number++)
  {
    next_row = row_vnp->next;
    if (sequence_numbers[row_number] < 0) {
      prev_row->next = next_row;
      row_vnp->next = NULL;
      row_vnp = FreeTableDisplayRowList (row_vnp);
    } else {
      prev_row = row_vnp;
    } 
  }
  
  /* now remove sequence_numbers entries */
  deleted_rows = 0;
  row_number = 0;
  while (row_number < *num_rows)
  {
    if (sequence_numbers[row_number] < 0) {
      for (i = row_number + 1; i < *num_rows; i++) {
        sequence_numbers[i - 1] = sequence_numbers[i];
      }
      (*num_rows)--;
    } else {
      row_number++;
    }
  }   
  
  return rval;
}

/* This checks the column names and puts the modifier type in the choice for each column */
static Boolean ValidateImportModifierColumnNames (ValNodePtr header_line, ValNodePtr preferred_list)
{
  ValNodePtr      header_vnp;
  Boolean         rval = TRUE;
  ModifierInfoPtr mip;
  CharPtr         orig_name;
  Int4            col_num;
  
  if (header_line == NULL)
  {
    return FALSE;
  }
  
  header_vnp = header_line->data.ptrvalue;
  if (header_vnp == NULL || header_vnp->next == NULL)
  {
    return FALSE;
  }
  
  /* check ID column */
  if (StringICmp (header_vnp->data.ptrvalue, "local_id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "local id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "local-id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "seq_id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "seq id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "seq-id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "seqid") != 0
      && StringICmp (header_vnp->data.ptrvalue, "sequence_id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "sequence id") != 0
      && StringICmp (header_vnp->data.ptrvalue, "sequence-id") != 0
      )
  {
    Message (MSG_ERROR, "Table file is missing header line!  Make sure first column header is seq_id");
    return FALSE;      
  }
  header_vnp = header_vnp->next;
  col_num = 1;
  while (header_vnp != NULL && rval)
  {
    mip = MakeModifierInfoFromNameAndValue (header_vnp->data.ptrvalue, NULL);
    if (mip == NULL 
      || (mip->modtype == eModifierType_SourceQual
        	&& mip->subtype == 255
  	      && StringICmp (mip->name, "note-subsrc") != 0 
  	      && StringICmp (mip->name, "note-orgmod") != 0))
    {
      orig_name = (CharPtr) header_vnp->data.ptrvalue;
      rval = ReplaceImportModifierName (&orig_name, col_num, preferred_list);
      header_vnp->data.ptrvalue = orig_name;
    }
    else
    {
      header_vnp->data.ptrvalue = MemFree (header_vnp->data.ptrvalue);
      header_vnp->data.ptrvalue = StringSave (mip->name);
      header_vnp->choice = mip->modtype;
    }
    mip = ModifierInfoFree (mip);
    header_vnp = header_vnp->next;
    col_num++;
  }
  return rval;
}

static Boolean StringAlreadyInList (ValNodePtr list, CharPtr str)
{
  while (list != NULL)
  {
    if (StringICmp (list->data.ptrvalue, str) == 0)
    {
      return TRUE;
    }
    list = list->next;
  }
  return FALSE;
}

static Boolean ValidateTableValues (ValNodePtr header_line)
{
  ValNodePtr      header_vnp, row_vnp, col_vnp;
  Boolean         rval = TRUE;
  ModifierInfoPtr mip;
  Int4            col_num;
  ValNodePtr      bad_value_columns = NULL;
  ValNodePtr      bad_nontext_columns = NULL;
  CharPtr         err_msg;
  
  if (header_line == NULL || header_line->next == NULL 
      || header_line->data.ptrvalue == NULL)
  {
    return FALSE;
  }
    
  for (row_vnp = header_line->next; row_vnp != NULL; row_vnp = row_vnp->next)
  {
    /* skip rows with bad sequence IDs */
    if (row_vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    
    header_vnp = header_line->data.ptrvalue;
    col_vnp = row_vnp->data.ptrvalue;
    /* skip ID column */
    header_vnp = header_vnp->next;
    col_vnp = col_vnp->next;
    for (col_num = 1; 
         header_vnp != NULL && col_vnp != NULL; 
         header_vnp = header_vnp->next, col_vnp = col_vnp->next, col_num++)
    {
      mip = MakeModifierInfoFromNameAndValue (header_vnp->data.ptrvalue, 
                                              col_vnp->data.ptrvalue);
      if (mip->modtype == eModifierType_SourceQual
  	             && IsNonTextModifier (mip->name))
      {
        if (StringICmp (mip->value, "TRUE") != 0
            && StringICmp (mip->value, "FALSE") != 0)
        {
          if (!StringAlreadyInList (bad_nontext_columns, header_vnp->data.ptrvalue))
          {
            ValNodeAddPointer (&bad_nontext_columns, col_num, StringSave (header_vnp->data.ptrvalue));
          }
        }
      }
      else if (ModifierHasInvalidValue (mip))
      {
        if (!StringAlreadyInList (bad_value_columns, header_vnp->data.ptrvalue))
        {
          ValNodeAddPointer (&bad_value_columns, col_num, StringSave (header_vnp->data.ptrvalue));
        }
      }
      mip = ModifierInfoFree (mip);
    }
  }
  
  if (bad_value_columns != NULL)
  {
    err_msg = CreateListMessage ("Your file contains invalid values for column",
                                 ". Please edit your file to list valid values.",
                                 bad_value_columns);
    Message (MSG_ERROR, err_msg);                                 
    rval = FALSE;
  }
  if (bad_nontext_columns != NULL && rval)
  {
    err_msg = CreateListMessage ("Your file contains values other than TRUE or FALSE for column",
                                 ". These modifiers do not allow other text.  Click OK to "
                                 "discard this text and mark the values as TRUE.  If you "
                                 "wish to preserve this text under another modifier, click "
                                 "Cancel and change the column header in your file.",
                                 bad_nontext_columns);
    if (ANS_CANCEL == Message (MSG_OKC, err_msg))
    {
      rval = FALSE;
    }
  }
  
  bad_value_columns = ValNodeFreeData (bad_value_columns);
  bad_nontext_columns = ValNodeFreeData (bad_nontext_columns);
  return rval;
}

static Boolean 
CheckModifiersForOverwrite 
(ValNodePtr        header_line,
 IDAndTitleEditPtr iatep,
 Int4Ptr           sequence_numbers,
 Int4              num_rows,
 BoolPtr           erase_where_blank,
 BoolPtr           parse_multiple)
{
  ValNodePtr row_vnp, header_vnp, col_vnp;
  CharPtr    title_val, data_val;
  ValNodePtr blank_column_list = NULL;
  ValNodePtr replace_column_list = NULL;
  ValNodePtr parse_multi_list = NULL;
  Int4       col_num, row_num;
  Boolean    rval = TRUE;
  CharPtr    err_msg;
  MsgAnswer  ans;
  
  if (header_line == NULL || header_line->next == NULL || iatep == NULL
      || sequence_numbers == NULL || num_rows < ValNodeLen (header_line->next)
      || erase_where_blank == NULL || parse_multiple == NULL)
  {
    return FALSE;
  }
  
  *erase_where_blank = FALSE;
  *parse_multiple = FALSE;
  
  for (row_vnp = header_line->next, row_num = 0;
       row_vnp != NULL && row_num < num_rows;
       row_vnp = row_vnp->next, row_num++)
  {
    if (row_vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    header_vnp = header_line->data.ptrvalue;
    col_vnp = row_vnp->data.ptrvalue;
    
    /* skip ID column */
    header_vnp = header_vnp->next;
    col_vnp = col_vnp->next;
    
    col_num = 1;
    while (header_vnp != NULL && col_vnp != NULL)
    {
      /* if column name is blank, skip */
      if (header_vnp->data.ptrvalue != NULL)
      {
        title_val = FindValueFromPairInDefline (header_vnp->data.ptrvalue,
                                                iatep->title_list [sequence_numbers[row_num]]);
        data_val = col_vnp->data.ptrvalue;
        if (!StringHasNoText (title_val))
        {
          if (StringHasNoText (data_val))
          {
            /* add to list of possible erasures */
            if (!StringAlreadyInList (blank_column_list, header_vnp->data.ptrvalue))
            {
              ValNodeAddPointer (&blank_column_list, col_num, StringSave (header_vnp->data.ptrvalue));
            }
          }
          else if (StringCmp (data_val, title_val) != 0)
          {
            /* add to list of possible replacements */
            if (!StringAlreadyInList (replace_column_list, header_vnp->data.ptrvalue))
            {
              ValNodeAddPointer (&replace_column_list, col_num, StringSave (header_vnp->data.ptrvalue));
            }
          }
        }
        title_val = MemFree (title_val);
        /* check for multival parsing */
        if (data_val != NULL 
            && data_val [0] == '(' && data_val [StringLen (data_val) - 1] == ')'
            && StringChr (data_val, ';') != NULL
            && !StringAlreadyInList (parse_multi_list, header_vnp->data.ptrvalue))
        {
          ValNodeAddPointer (&parse_multi_list, col_num, StringSave (header_vnp->data.ptrvalue));
        }
      }
      header_vnp = header_vnp->next;
      col_vnp = col_vnp->next;
      col_num++;
    }
  }
    
  if (replace_column_list != NULL)
  {
    err_msg = CreateListMessage ("Record already contains values for column",
                                 " also found in the import table.\n"
                                 "Do you wish to overwrite these values?",
                                 replace_column_list);
    if (ANS_NO == Message (MSG_YN, err_msg))
    {
      rval = FALSE;
    }
    err_msg = MemFree (err_msg);
  }
  
  if (blank_column_list != NULL && rval)
  {
    err_msg = CreateListMessage ("Your import table contains blanks in column",
                                 " where data already exists in the sequences.\n"
                                 "Do you wish to erase these values in the sequences?\n"
                                 "If you say no, the old values will remain.",
                                 blank_column_list);
    ans = Message (MSG_YNC, err_msg);
    err_msg = MemFree (err_msg);
    if (ans == ANS_CANCEL)
    {
      rval = FALSE;
    }
    else if (ans == ANS_YES)
    {
      *erase_where_blank = TRUE;
    }
  }

#if 0  
  /* ability to parse multiple entry format removed (for now) */
  if (parse_multi_list != NULL && rval)
  {
    err_msg = CreateListMessage ("Your import table contains values in column",
                                 " where the values are in form '(value1;value2)'.\n"
                                 "Do you wish to parse these values into multiple modifiers?\n"
                                 "If you say no, the values will be applied to a single modifier.",
                                 parse_multi_list);
    ans = Message (MSG_YNC, err_msg);
    err_msg = MemFree (err_msg);
    if (ans == ANS_CANCEL)
    {
      rval = FALSE;
    }
    else if (ans == ANS_YES)
    {
      *parse_multiple = TRUE;
    }
  }
#endif  
  
  blank_column_list = ValNodeFree (blank_column_list);  
  replace_column_list = ValNodeFreeData (replace_column_list);
  parse_multi_list = ValNodeFreeData (parse_multi_list);
  
  return rval;
}

NLM_EXTERN Boolean ImportModifiersToIDAndTitleEditEx (IDAndTitleEditPtr iatep, ValNodePtr preferred_list)
{
  ValNodePtr   header_line, row_vnp, col_vnp, header_vnp, special_list;
  Boolean      erase_where_blank = FALSE, parse_multi = FALSE;
  Int4Ptr      sequence_numbers;
  Int4         num_rows, row_number;
  
  if (iatep == NULL)
  {
    return FALSE;
  }

  SendHelpScrollMessage (helpForm, "Organism Page", "Import Source Modifiers");  
  
  header_line = ReadRowListFromFile ();
  if (header_line == NULL || header_line->next == NULL)
  {
    header_line = FreeTableDisplayRowList (header_line);
    return FALSE;
  }
  
  header_vnp = header_line->data.ptrvalue;
  if (header_vnp == NULL || header_vnp->next == NULL)
  {
    header_line = FreeTableDisplayRowList (header_line);
    return FALSE;
  }

  /* remove special characters */
  special_list = ScanTabTableForSpecialCharacters (header_line);
  if (special_list != NULL 
      && !FixSpecialCharactersForStringsInList (special_list, 
                                                "The table contains special characters\nand cannot be used until they are replaced.",
                                                FALSE)) {
    special_list = FreeContextList (special_list);
    header_line = FreeTableDisplayRowList (header_line);
    return FALSE;
  }
  special_list = FreeContextList (special_list);
  
  num_rows = ValNodeLen (header_line->next);
  sequence_numbers = (Int4Ptr) MemNew (num_rows * sizeof (Int4));
  
  if (!ValidateModifierTableSequenceIDs (header_line, iatep, sequence_numbers, &num_rows))
  {
    header_line = FreeTableDisplayRowList (header_line);
    sequence_numbers = MemFree (sequence_numbers);
    return FALSE;
  }
  
  /* first, validate all column names and values */
  if (!ValidateImportModifierColumnNames (header_line, preferred_list))
  {
    header_line = FreeTableDisplayRowList (header_line);
    sequence_numbers = MemFree (sequence_numbers);
    return FALSE;
  }
  
  if (!ValidateTableValues (header_line))
  {
    header_line = FreeTableDisplayRowList (header_line);
    sequence_numbers = MemFree (sequence_numbers);
    return FALSE;
  }
  
  if (!CheckModifiersForOverwrite (header_line, iatep, 
                                   sequence_numbers, num_rows, 
                                   &erase_where_blank, &parse_multi))
  {
    header_line = FreeTableDisplayRowList (header_line);
    sequence_numbers = MemFree (sequence_numbers);
    return FALSE;
  }
  
  /* now apply */
  for (row_vnp = header_line->next, row_number = 0;
       row_vnp != NULL && row_number < num_rows; 
       row_vnp = row_vnp->next, row_number++)
  {
    if (row_vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    header_vnp = header_line->data.ptrvalue;
    col_vnp = row_vnp->data.ptrvalue;
    
    /* skip the ID column */
    header_vnp = header_vnp->next;
    col_vnp = col_vnp->next;
    
    for (;
         header_vnp != NULL && col_vnp != NULL;
         header_vnp = header_vnp->next, col_vnp = col_vnp->next)
    {
      iatep->title_list [sequence_numbers [row_number]] = ApplyImportModToTitle (iatep->title_list [sequence_numbers[row_number]],
                                                                   header_vnp->data.ptrvalue,
                                                                   col_vnp->data.ptrvalue,
                                                                   erase_where_blank,
                                                                   parse_multi);
    }
  }
  sequence_numbers = MemFree (sequence_numbers);  
  return TRUE;
}


NLM_EXTERN Boolean ImportModifiersToIDAndTitleEdit (IDAndTitleEditPtr iatep)
{
  return ImportModifiersToIDAndTitleEditEx (iatep, NULL);
}


typedef struct fastapage {
  DIALOG_MESSAGE_BLOCK
  Char         path [PATH_MAX];
  SeqEntryPtr  list;
  ValNodePtr   errmsgs;
  DoC          doc;
  GrouP        instructions;
  GrouP        have_seq_instr_grp;
  GrouP        singleIdGrp;
  TexT         singleSeqID;  
  Boolean      is_na;
  Boolean      is_mrna;
  Boolean      is_delta;
  Boolean      parseSeqId;
  Boolean      single;
  Int2Ptr      seqPackagePtr;
  ButtoN       import_btn;
  ButtoN       clear_btn;
} FastaPage, PNTR FastaPagePtr;

static ParData faParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData faColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static void ResetFastaPage (FastaPagePtr fpp)

{
  SeqEntryPtr  next;
  SeqEntryPtr  sep;

  if (fpp != NULL) {
    sep = fpp->list;
    while (sep != NULL) {
      next = sep->next;
      sep->next = NULL;
      SeqEntryFree (sep);
      sep = next;
    }
    fpp->list = NULL;
    fpp->errmsgs = ValNodeFreeData (fpp->errmsgs);
  }
}

static CharPtr GetModValueFromSeqEntry (SeqEntryPtr sep, CharPtr mod_name)
{
  CharPtr ttl = NULL;
  CharPtr value = NULL;
  
  if (sep == NULL || StringHasNoText (mod_name))
  {
    return NULL;
  }

  SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
  if (StringHasNoText (ttl))
  {
    return NULL;
  }
  
  value =  FindValueFromPairInDefline (mod_name, ttl);
  
  return value;  
}

static void AddReportLine (CharPtr str, CharPtr name, CharPtr tmp)

{
  StringCat (str, name);
  StringCat (str, ": ");
  StringCat (str, tmp);
  StringCat (str, "\n");
}

static CharPtr GetDisplayValue (CharPtr mod_name, CharPtr title, BoolPtr multi_found);
static CharPtr GetDisplayValueFromModifierInfoList (CharPtr mod_name, ValNodePtr modifier_info_list, BoolPtr multi_found);

static void ReportModifiers (CharPtr str, CharPtr report_name,
                             ValNodePtr modifier_info_list, CharPtr mod_name, CharPtr not_found_msg)
{
  CharPtr valstr;
  Boolean multi_found = TRUE;

  valstr = GetDisplayValueFromModifierInfoList (mod_name, modifier_info_list, &multi_found);
  if (IsNonTextModifier (mod_name) && StringICmp (valstr, "FALSE") == 0)
  {
  	valstr = MemFree (valstr);
  }

  if (!StringHasNoText (valstr)) {
    AddReportLine (str, report_name, valstr);
  } else if (!StringHasNoText (not_found_msg)) {
    StringCat (str, not_found_msg);
  }
  valstr = MemFree (valstr);
}

static void LookupAndAddReportLine (CharPtr str, CharPtr report_name, 
                                    CharPtr title, CharPtr mod_name, CharPtr not_found_msg)
{
  CharPtr valstr;
  Boolean multi_found = TRUE;
  
  valstr = GetDisplayValue (mod_name, title, &multi_found);
  if (IsNonTextModifier (mod_name) && StringICmp (valstr, "FALSE") == 0)
  {
  	valstr = MemFree (valstr);
  }

  if (!StringHasNoText (valstr)) {
    AddReportLine (str, report_name, valstr);
  } else if (!StringHasNoText (not_found_msg)) {
    StringCat (str, not_found_msg);
  }
  valstr = MemFree (valstr);
}

static void LookupAndAddLocationReportLine (CharPtr str, CharPtr title)
{
  CharPtr valstr;
  
  valstr = FindValueFromPairInDefline ("location", title);
  if (!StringHasNoText (valstr) && StringICmp (valstr, "genomic") != 0) {
    AddReportLine (str, "Location", valstr);
  }
  valstr = MemFree (valstr);
}

static CharPtr singlewarn = "\
ERROR - You may not enter multiple segments for a single sequence submission.\
You should either clear the nucleotide and import a single FASTA record, or\
return to the Sequence Format form and choose the proper submission type.\n\n";

#define FastaFormatBufLen 2000

static Int4 CountSegSetSegments (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL || ! IS_Bioseq_set (sep))
  {
    return 0;
  }
  
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_segset)
  {
    return 0;
  }
  sep = bssp->seq_set;
  
  while (sep != NULL)
  {
    if (IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL)
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp->_class == BioseqseqSet_class_parts)
      {
        return ValNodeLen (bssp->seq_set);
      }
    }
    sep = sep->next;
  }
  return 0;
}

static void FormatFastaDoc (FastaPagePtr fpp)

{
  Nlm_QualNameAssocPtr ap;
  BioseqPtr          bsp;
  Boolean            hasErrors;
  CharPtr            label;
  Int4               len;
  CharPtr            measure;
  SeqEntryPtr        nsep = NULL;
  Int2               num;
  CharPtr            plural;
  CharPtr            ptr;
  SeqIdPtr           sip;
  SeqEntryPtr        sep;
  CharPtr            str;
  CharPtr            title;
  CharPtr            ttl;
  CharPtr            tmp;
  ValNodePtr         vnp;
  Int4               num_seg;
  CharPtr            valstr;
  ValNodePtr         modifier_info_list = NULL;

  if (fpp != NULL) {
    str = MemNew (sizeof (char) * FastaFormatBufLen);
    tmp = MemNew (sizeof (char) * FastaFormatBufLen);
    if (str == NULL || tmp == NULL) return;
    num = 0;
    len = 0;
    hasErrors = FALSE;
    for (sep = fpp->list; sep != NULL; sep = sep->next) {
      num++;
      if (IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          len += bsp->length;
        }
      } else if (IS_Bioseq_set (sep)) {
        nsep = FindNucSeqEntry (sep);
        if (nsep != NULL && IS_Bioseq (nsep)) {
          bsp = (BioseqPtr) nsep->data.ptrvalue;
          if (bsp != NULL) {
            len += bsp->length;
          }
        }
      }
    }
    if (num > 1) {
      plural = "s";
    } else {
      plural = "";
    }
    if (fpp->single && num > 1) {
      AppendText (fpp->doc, singlewarn, &faParFmt, &faColFmt, programFont);
      hasErrors = TRUE;
    }
    if (fpp->is_mrna) {
      label = "Message";
      measure = "nucleotides";
    } else if (fpp->is_na) {
      label = "Sequence";
      measure = "bases";
    } else {
      label = "Sequence";
      measure = "amino acids";
    }
    if (fpp->is_mrna) {
      sprintf (str, "%ld transcript sequence%s, total length %ld %s\n",
               (long) num, plural, (long) len, measure);
    } else if (fpp->is_na) {
      sprintf (str, "%ld nucleotide sequence%s, total length %ld %s\n",
               (long) num, plural, (long) len, measure);
    } else {
      sprintf (str, "%ld protein sequence%s, total length %ld %s\n",
               (long) num, plural, (long) len, measure);
    }
    AppendText (fpp->doc, str, &faParFmt, &faColFmt, programFont);
    vnp = fpp->errmsgs;
    num = 0;
    for (sep = fpp->list; sep != NULL; sep = sep->next) {
      num++;
      len = 0;
      num_seg = CountSegSetSegments (sep);
      sip = NULL;
      tmp [0] = '\0';
      if (IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          len = bsp->length;
          sip = SeqIdFindWorst (bsp->id);
          SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
        }
        nsep = sep;
      } else if (IS_Bioseq_set (sep)) {
        nsep = FindNucSeqEntry (sep);
        if (nsep != NULL && IS_Bioseq (nsep)) {
          bsp = (BioseqPtr) nsep->data.ptrvalue;
          if (bsp != NULL) {
            len = bsp->length;
            sip = SeqIdFindWorst (bsp->id);
            SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
          }
        }
      }
      
      /* if segmented set, show number of segments */
      if (num_seg > 0)
      {
        sprintf (str, "\nSegset %d Sequence ID: %s\nLength: %ld %s (%d segments)\n",
                 (int) num, tmp, (long) len, measure, num_seg);
      }
      else
      {
        sprintf (str, "\n%s %d Sequence ID: %s\nLength: %ld %s\n", label,
                 (int) num, tmp, (long) len, measure);
      }
      ttl = NULL;
      SeqEntryExplore (nsep, (Pointer) (&ttl), FindFirstTitle);
      title = StringSaveNoNull (ttl);
      modifier_info_list = ParseAllBracketedModifiers (title);
      if (title != NULL && (! fpp->is_na)) {

        ReportModifiers (str, "Gene", modifier_info_list, "gene", "No gene name detected\n");
        ReportModifiers (str, "Protein", modifier_info_list, "protein", "No protein name detected\n");
        ReportModifiers (str, "Gene Syn", modifier_info_list, "gene_syn", NULL);
        ReportModifiers (str, "Protein Desc", modifier_info_list, "protein_desc", NULL);

        ptr = StringISearch (title, "[orf]");
        if (ptr != NULL) {
        StringCat (str, "ORF indicated\n");
        }
        ReportModifiers (str, "Protein Comment", modifier_info_list, "comment", NULL);
      }

      if (title != NULL && fpp->is_na && (! fpp->is_mrna)) {
        ReportModifiers (str, "Organism", modifier_info_list, "organism", NULL);
        ReportModifiers (str, "Lineage", modifier_info_list, "lineage", NULL);
        for (ap = current_orgmod_subtype_alist; ap->name != NULL; ap++) {
          ReportModifiers (str, ap->name, modifier_info_list, ap->name, NULL); 
        }
        for (ap = current_subsource_subtype_alist; ap->name != NULL; ap++) {
          ReportModifiers (str, ap->name, modifier_info_list, ap->name, NULL); 
        }
        LookupAndAddReportLine (str, "Note", title, "note", NULL); 
        LookupAndAddReportLine (str, "Note", title, "subsource", NULL); 
        LookupAndAddReportLine (str, "Molecule", title, "molecule", NULL); 
        LookupAndAddReportLine (str, "MolType", title, "moltype", NULL); 
        LookupAndAddLocationReportLine (str, title); 
        LookupAndAddReportLine (str, "Genetic Code", title, "genetic_code", NULL);
      }

      if (title != NULL && fpp->is_na && fpp->is_mrna) {
        LookupAndAddReportLine (str, "Gene", title, "gene", "No gene name detected\n"); 
        valstr = FindValueFromPairInDefline ("mrna", title);
        if (!StringHasNoText (valstr)) {
          AddReportLine (str, "mRNA", valstr);
          valstr = MemFree (valstr);
        } else {
          valstr = MemFree (valstr);
          valstr = FindValueFromPairInDefline ("cdna", title);
          if (!StringHasNoText (valstr)) {
            AddReportLine (str, "cDNA", valstr);
          } else {
            StringCat (str, "No mRNA name detected\n");
          }
          valstr = MemFree (valstr);
        }
        LookupAndAddReportLine (str, "Comment", title, "comment", NULL); 
      }
      MemFree (title);
      ttl = NULL;
      SeqEntryExplore (nsep, (Pointer) (&ttl), FindFirstTitle);
      title = StringSaveNoNull (ttl);
      if (title != NULL) {
        RemoveRecognizedModifiersFromTitle (title, modifier_info_list, fpp->is_na);
        if (fpp->is_mrna) {
          StripAllInstancesOfModNameFromTitle ("gene", title);
          StripAllInstancesOfModNameFromTitle ("mrna", title);
          StripAllInstancesOfModNameFromTitle ("cdna", title);
          StripAllInstancesOfModNameFromTitle ("comment", title);
        } 
        TrimSpacesAroundString (title);
        if (! StringHasNoText (title)) {
          StringCat (str, "Title: ");
          StringNCat (str, title, 128);
          StringCat (str, "\n");
        } else {
          StringCat (str, "No title detected\n");
        }
      }
      MemFree (title);
      ModifierInfoListFree (modifier_info_list);
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        hasErrors = TRUE;
        StringCat (str, (CharPtr) vnp->data.ptrvalue);
        StringCat (str, "\n");
      }
      AppendText (fpp->doc, str, &faParFmt, &faColFmt, programFont);
      if (vnp != NULL) {
        vnp = vnp->next;
      }
    }
    MemFree (str);
    MemFree (tmp);
    UpdateDocument (fpp->doc, 0, 0);
    if (hasErrors) {
      Beep ();
      Beep ();
      Beep ();
    }
  }
}

extern SeqEntryPtr ImportOneGappedSequence (FILE *fp)
{
  BioseqPtr      bsp;
  Pointer        dataptr;
  Uint2          datatype;
  SeqEntryPtr    topsep;
  SeqSubmitPtr   ssp;
  ErrSev         oldsev;
  
  if (fp == NULL) return NULL;
  
  oldsev = ErrSetMessageLevel (SEV_MAX);
  bsp = ReadDeltaFasta (fp, NULL);
  ErrSetMessageLevel (oldsev);
  if (bsp == NULL)
  {
    topsep = NULL;
    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE,
		  		    TRUE, FALSE);
    if (dataptr != NULL)
    {
      /* Get a pointer to the new SeqEntry */
      if (datatype == OBJ_SEQENTRY)
      {
        topsep = (SeqEntryPtr) dataptr;
      }
      else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
      {
        topsep = SeqMgrGetSeqEntryForData (dataptr);
      }
      else if (datatype == OBJ_SEQSUB) 
      {
        ssp = (SeqSubmitPtr) dataptr;
        if (ssp != NULL && ssp->datatype == 1)
        {
          topsep = (SeqEntryPtr) ssp->data;
        }
      }
    }
  }
  else
  {
    topsep = SeqMgrGetSeqEntryForData (bsp);
  }

  return topsep;
}

static SeqEntryPtr SegsetFromSeqEntryList (SeqEntryPtr list)
{
  SeqEntryPtr  first_sep, tmp_sep, next_sep;
  BioseqPtr    bsp;
  SeqDescrPtr  sdp = NULL, set_sdp;
  
  if (list == NULL)
  {
    return NULL;
  }
  
  first_sep = list;
  next_sep = first_sep->next;
  first_sep->next = NULL;
  
  /* grab title on first sequence to put on segmented bioseq */
  if (IS_Bioseq (first_sep) && first_sep->data.ptrvalue != NULL)
  {
    bsp = (BioseqPtr) first_sep->data.ptrvalue;
    sdp = bsp->descr;
    while (sdp != NULL && sdp->choice != Seq_descr_title)
    {
      sdp = sdp->next;
    }
  }

  while (next_sep != NULL)
  {
    tmp_sep = next_sep;
    next_sep = tmp_sep->next;
    tmp_sep->next = NULL;
    AddSeqEntryToSeqEntry (first_sep, tmp_sep, TRUE);
  }
  
  if (sdp != NULL && IS_Bioseq_set (first_sep))
  {
    tmp_sep = FindNucSeqEntry (first_sep);
    if (tmp_sep != NULL && IS_Bioseq (tmp_sep) && tmp_sep->data.ptrvalue != NULL)
    {
      bsp = tmp_sep->data.ptrvalue;
      set_sdp = bsp->descr;
      while (set_sdp != NULL && set_sdp->choice != Seq_descr_title)
      {
        set_sdp = set_sdp->next;
      }
      if (set_sdp == NULL)
      {
        set_sdp = CreateNewDescriptor (tmp_sep, Seq_descr_title);
      }
      if (set_sdp != NULL && StringHasNoText (set_sdp->data.ptrvalue))
      {
        /* make a copy, rather than removing the segment title */
        set_sdp->data.ptrvalue = MemFree (set_sdp->data.ptrvalue);
        set_sdp->data.ptrvalue = StringSave (sdp->data.ptrvalue);
      }
    }
  }  
  
  return first_sep;
}


static SeqEntryPtr 
ReadOneSegSet 
(FILE            *fp,
 Boolean         parse_id,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped)
{
  SeqEntryPtr nextsep;
  CharPtr     errormsg = NULL;
  Char        lastchar;
  SeqEntryPtr seg_list = NULL, seg_list_last = NULL;
  BioseqPtr   bsp;
  
  if (fp == NULL)
  {
    return NULL;
  }
  
  /* note - we pass in FALSE for parse_id in SequinFastaToSeqEntryEx
   * because we do not want to use Sequin's auto-generated sequence IDs.
   * We then parse the sequence ID from the title ourselves using
   * ReplaceFakeIDWithIDFromTitle if parse_id is TRUE, or leave the ID
   * as blank to force the user to select a real ID later.
   */
  nextsep = SequinFastaToSeqEntryExEx (fp, TRUE, &errormsg, FALSE, &lastchar, chars_stripped);
  while (nextsep != NULL ||
         (lastchar != (Char) EOF && lastchar != NULLB && lastchar != (Char) 255
          && lastchar != ']')) 
  {
    if (nextsep != NULL) 
    {
      /* replace fake ID with ID from title */
      if (IS_Bioseq (nextsep) && nextsep->data.ptrvalue != NULL)
      {
        bsp = (BioseqPtr) nextsep->data.ptrvalue;
        if (parse_id)
        {
          ReplaceFakeIDWithIDFromTitle ((BioseqPtr) nextsep->data.ptrvalue);
        }
        else
        {
          bsp->id = SeqIdFree (bsp->id);
        }
      }
      SeqEntryPack (nextsep); 
      if (seg_list_last == NULL)
      {
        seg_list = nextsep;
      }
      else
      {
        seg_list_last->next = nextsep;
      }
      seg_list_last = nextsep;
      
      ValNodeAddPointer (err_msg_list, 0, errormsg);
      errormsg = NULL;
    }
    nextsep = SequinFastaToSeqEntryExEx (fp, TRUE, &errormsg, FALSE, &lastchar, chars_stripped);
  }
  nextsep = SegsetFromSeqEntryList (seg_list);
  return nextsep;
}

static void AddDefaultMoleculeTypeToIDAndTitleEdit (IDAndTitleEditPtr iatep)
{
  Int4    seq_num;
  CharPtr old_value;
  
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    old_value = FindValueFromPairInDefline("moltype", 
                                           iatep->title_list [seq_num]);
    if (StringHasNoText (old_value) || StringICmp (old_value, "dna") == 0)
    {
      iatep->title_list [seq_num] = ReplaceValueInOneDefLine(iatep->title_list [seq_num],
                                                             "moltype",
                                                             "Genomic DNA");
    }
    old_value = MemFree (old_value);
  }
}

static void AddDefaultLocationToIDAndTitleEdit (IDAndTitleEditPtr iatep)
{
  Int4    seq_num;
  CharPtr old_value, first_organism, next_org_loc = NULL, org_stop;
  
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    first_organism = FindValuePairInDefLine ("organism", iatep->title_list [seq_num], &org_stop);
    if (first_organism != NULL)
    {
      next_org_loc = FindValuePairInDefLine ("organism", org_stop + 1, NULL);
    }
    else
    {
      next_org_loc = NULL;
    }
    old_value = FindValueFromPairInDeflineBeforeCharPtr ("location", 
                                                         iatep->title_list [seq_num],
                                                         next_org_loc);
    if (StringHasNoText (old_value))
    {
      iatep->title_list [seq_num] = ReplaceValueInOneDefLineForOrganism (iatep->title_list [seq_num],
                                                                         "location",
                                                                         "genomic",
                                                                         first_organism);
    }
    old_value = MemFree (old_value);
  }
}

static void AddDefaultTopologyToIDAndTitleEdit (IDAndTitleEditPtr iatep)
{
  Int4    seq_num;
  CharPtr old_value;
  
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    old_value = FindValueFromPairInDefline ("topology", 
                                            iatep->title_list [seq_num]);
    if (StringHasNoText (old_value))
    {
      iatep->title_list [seq_num] = ReplaceValueInOneDefLine(iatep->title_list [seq_num],
                                                             "topology",
                                                             "Linear");
    }
    old_value = MemFree (old_value);
  }
}

static void AddDefaultGeneticCodesToIDAndTitleEdit (IDAndTitleEditPtr iatep)
{
  CharPtr     taxname, location, gcode_name;
  Int4        gcode;
  ValNodePtr  gencodelist;
  Int4        seq_num;
  CharPtr     first_organism, next_org_loc = NULL, org_stop;

  if (iatep == NULL)
  {
    return;
  }
  
  gencodelist = GetGeneticCodeValNodeList ();
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    first_organism = FindValuePairInDefLine ("organism", iatep->title_list [seq_num], &org_stop);
    if (first_organism != NULL)
    {
      next_org_loc = FindValuePairInDefLine ("organism", org_stop + 1, NULL);
    }
    else
    {
      next_org_loc = NULL;
    }
    
    taxname = FindValueFromPairInDefline ("organism", first_organism);
    location = FindValueFromPairInDeflineBeforeCharPtr ("location", 
                                                        iatep->title_list [seq_num],
                                                        next_org_loc);
    if (StringHasNoText (location))
    {
      location = StringSave ("genomic");
    }
    
    gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
    taxname = MemFree (taxname);
    location = MemFree (location);
    
    if (gcode < 0)
    {
      gcode_name = FindValueFromPairInDeflineBeforeCharPtr ("genetic_code",
                                                            iatep->title_list [seq_num],
                                                            next_org_loc);
      if (StringHasNoText (gcode_name))
      {
        gcode_name = MemFree (gcode_name);
        gcode_name = GeneticCodeStringFromIntAndList (1, gencodelist);
        iatep->title_list [seq_num] = ReplaceValueInOneDefLineForOrganism (iatep->title_list [seq_num],
                                                                         "genetic_code",
                                                                         gcode_name,
                                                                         first_organism);
      }
      else
      {
        gcode_name = MemFree (gcode_name);
      }
    }
    else
    {
      gcode_name = GeneticCodeStringFromIntAndList (gcode, gencodelist);
      iatep->title_list [seq_num] = ReplaceValueInOneDefLineForOrganism (iatep->title_list [seq_num],
                                                                         "genetic_code",
                                                                         gcode_name,
                                                                         first_organism);
    }
  }
  ValNodeFreeData (gencodelist);
}

static void AddDefaultModifierValues (SeqEntryPtr seq_list)
{
  IDAndTitleEditPtr iatep;
  
  iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  AddDefaultMoleculeTypeToIDAndTitleEdit (iatep);
  AddDefaultLocationToIDAndTitleEdit (iatep);
  AddDefaultTopologyToIDAndTitleEdit (iatep);
  AddDefaultGeneticCodesToIDAndTitleEdit (iatep);
  ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
  iatep = IDAndTitleEditFree (iatep);
}


static void CleanTitles (SeqEntryPtr sep, ValNodePtr PNTR special_list)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp = NULL;

  while (sep != NULL)
  {
    sdp = NULL;
    if (sep->choice == 1)
    {
      bsp = sep->data.ptrvalue;
      if (bsp != NULL) 
      {
        sdp = bsp->descr;
      }
    }
    else if (sep->choice == 2)
    {
      bssp = sep->data.ptrvalue;
      if (bssp != NULL)
      {
        CleanTitles (bssp->seq_set, special_list);
        sdp = bssp->descr;
      }
    }
    while (sdp != NULL)
    {
      if (sdp->choice == Seq_descr_title)
      {
        SpecialCharFindWithContext ((CharPtr PNTR) &(sdp->data.ptrvalue), special_list, NULL, NULL);
      }
      sdp = sdp->next;
    }
    sep = sep->next;
  }
}


extern SeqEntryPtr 
ImportSequencesFromFileEx
(FILE           *fp, 
 SeqEntryPtr     sep_list,
 Boolean         is_na, 
 Boolean         parse_id,
 CharPtr         supplied_id_txt,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped,
 Boolean         allow_char_stripping)
{  
  SeqEntryPtr new_sep_list, last, oldscope;
  ValNodePtr    special_list = NULL;

  if (chars_stripped != NULL)
  {
    *chars_stripped = FALSE;
  }

  oldscope = SeqEntrySetScope (NULL);
  
  if (is_na) 
  {
    new_sep_list = ImportNucleotideFASTASequencesFromFile (fp, parse_id, supplied_id_txt,
                                                   err_msg_list, chars_stripped, allow_char_stripping);
  }
  else
  {
    new_sep_list = ImportProteinFASTASequences (fp, parse_id, supplied_id_txt, err_msg_list, chars_stripped);
  }

  CleanTitles (new_sep_list, &special_list);
  if (!FixSpecialCharactersForStringsInList (special_list, "Definition lines contain special characters.\nThe sequences cannot be imported unless the characters are replaced.", FALSE))
  {
    new_sep_list = SeqEntryFree (new_sep_list);
  }
  special_list = FreeContextList (special_list);

  last = sep_list;
  while (last != NULL && last->next != NULL) 
  {
    last = last->next;
  }
  if (last == NULL)
  {
    sep_list = new_sep_list;
  }
  else
  {
    last->next = new_sep_list;
  }

  SeqEntrySetScope (oldscope);

  return sep_list;
}


extern SeqEntryPtr 
ImportSequencesFromFile
(FILE           *fp, 
 SeqEntryPtr     sep_list,
 Boolean         is_na, 
 Boolean         parse_id,
 CharPtr         supplied_id_txt,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped)
{
  return ImportSequencesFromFileEx (fp, sep_list, is_na, parse_id, supplied_id_txt, err_msg_list, chars_stripped, FALSE);
}


static Boolean CollectIDsAndTitles (SeqEntryPtr new_list, SeqEntryPtr current_list, Boolean is_nuc);

static SeqEntryPtr RemoveZeroLengthSequences (SeqEntryPtr list, Int4Ptr pnum_seqs, Int4Ptr pnum_zero)
{
  SeqEntryPtr  prev_sep, next_sep, this_sep;
  Int4         num_seqs = 0, num_zero = 0;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  
  if (list == NULL)
  {
    return NULL;
  }

  prev_sep = NULL;
  this_sep = list;
  while (this_sep != NULL)
  {
    num_seqs++;
    next_sep = this_sep->next;
    if (this_sep->data.ptrvalue == NULL)
    {
      num_zero++;
      if (prev_sep == NULL)
      {
        list = next_sep;
      }
      else
      {
        prev_sep->next = next_sep;
      }
      this_sep->next = NULL;
      SeqEntryFree (this_sep);
    }
    else if (IS_Bioseq (this_sep))
    {
      bsp = (BioseqPtr) this_sep->data.ptrvalue;
      if (bsp->length == 0)
      {
        num_zero++;
        
        if (prev_sep == NULL)
        {
          list = next_sep;
        }
        else
        {
          prev_sep->next = next_sep;
        }
        this_sep->next = NULL;
        SeqEntryFree (this_sep);
      }
      else
      {
        prev_sep = this_sep;
      }
    }
    else if (IS_Bioseq_set (this_sep))
    {
      bssp = (BioseqSetPtr) this_sep->data.ptrvalue;
      bssp->seq_set = RemoveZeroLengthSequences (bssp->seq_set, pnum_seqs, pnum_zero);
      if (bssp->seq_set == NULL)
      {
        num_zero++;
        if (prev_sep == NULL)
        {
          list = next_sep;
        }
        else
        {
          prev_sep->next = next_sep;
        }
        this_sep->next = NULL;
        SeqEntryFree (this_sep);
      }
      else
      {
        prev_sep = this_sep;
      }
    }
    else
    {
      prev_sep = this_sep;
    }
    this_sep = next_sep;
  }
  
  if (pnum_seqs != NULL)
  {
    *pnum_seqs += num_seqs;
  }
  if (pnum_zero != NULL)
  {
    *pnum_zero += num_zero;
  }
  return list;
}

static Boolean RejectZeroLengthSequences (SeqEntryPtr PNTR new_list)
{
  SeqEntryPtr next_sep;
  Int4        num_zero = 0, num_seq = 0;
  Boolean     rval = TRUE;
  Boolean     delete_all = FALSE;
  
  if (new_list == NULL)
  {
    return FALSE;
  }
  
  *new_list = RemoveZeroLengthSequences (*new_list, &num_seq, &num_zero);

  if (num_zero > 0)
  {
    ResetSegSetIDLists (*new_list);
    if (num_zero == num_seq)
    {
      Message (MSG_ERROR, "The sequences in your file are empty - you cannot import them.");
      delete_all = TRUE;
      rval = FALSE;
    }
    else if (ANS_CANCEL == Message (MSG_OKC, "%d sequences in your file are empty and cannot be imported.  "
                                    "Would you like to import the remaining sequences?", num_zero))
    {
      delete_all = TRUE;
      rval = FALSE;
    }
    if (delete_all)
    {
      
      while ((*new_list) != NULL)
      {
        next_sep = (*new_list)->next;
        (*new_list)->next = NULL;
        SeqEntryFree (*new_list);
        *new_list = next_sep;
      }
    }
  }
  return rval;
}

static Boolean RejectExtraSequences (SeqEntryPtr new_list, FastaPagePtr fpp)
{
  SeqEntryPtr sep, next_sep;
  
  if (new_list == NULL || fpp == NULL)
  {
    return FALSE;
  }
  else if (!fpp->single || new_list->next == NULL)
  {
    return TRUE;
  }

  if (fpp->is_na 
           && fpp->seqPackagePtr != NULL 
           && *(fpp->seqPackagePtr) != SEQ_PKG_GENOMICCDNA)
  {
    if (Message (MSG_YN, "You are importing multiple sequences - did you intend to create a batch submission?") == ANS_YES)
    {
      *(fpp->seqPackagePtr) = SEQ_PKG_GENBANK;
      fpp->single = FALSE;
      SafeHide (fpp->singleIdGrp);
      return TRUE;
    }
  }
  if (Message (MSG_YN, "You cannot import multiple sequences - import the first one and ignore the rest?") == ANS_YES)
  {
    sep = new_list->next;
    new_list->next = NULL;
    while (sep != NULL)
    {
      next_sep = sep->next;
      sep->next = NULL;
      sep = SeqEntryFree (sep);
      sep = next_sep;
    }
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static void ShowImportHelp (ButtoN b)
{
  CharPtr help_msg;
  
  help_msg = (CharPtr) GetObjectExtra (b);
  if (help_msg == NULL)
  {
    return;
  }
  
  Message (MSG_OK, help_msg);
}

static Boolean OkToImport (CharPtr msg, CharPtr help_msg)
{
  WindoW w;
  GrouP  h, c;
  PrompT p;
  ButtoN b;
  ModalAcceptCancelData acd;
  
  if (msg == NULL)
  {
    return TRUE;
  }
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  
  p = StaticPrompt (h, msg, 0, 0, programFont, 'l');
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Yes", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "No", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  if (help_msg != NULL)
  {
    b = PushButton (c, "Help", ShowImportHelp);
    SetObjectExtra (b, help_msg, NULL);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  if (acd.accepted)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static CharPtr segset_import_help_str = "Segmented sequence: a collection of non-overlapping, non-contiguous sequences that cover a specified genetic region. A standard example is a set of genomic DNA sequences that encode exons from a gene along with fragments of their flanking introns.";
static CharPtr gapped_import_help_str = "Gapped sequence: a sequence with one or more gaps of known or unknown length.";


static Boolean ImportedSequenceTypeOk (SeqEntryPtr list, Int2 seqPackage)
{
  BioseqPtr bsp;
  Boolean   rval = TRUE;
  
  if (list == NULL || seqPackage != SEQ_PKG_SINGLE)
  {
    return TRUE;
  }
  if (list->choice == 1)
  {
    bsp = (BioseqPtr) list->data.ptrvalue;
    if (bsp != NULL && bsp->repr == Seq_repr_delta)
    {
      SendHelpScrollMessage (helpForm, "Sequence Format Form", NULL);
      rval = OkToImport ("You have imported a gapped sequence.  Did you mean to do that?",
                         gapped_import_help_str);
    }
  }
  else if (list->choice == 2)
  {
    SendHelpScrollMessage (helpForm, "Sequence Format Form", NULL);
    rval = OkToImport ("You have imported a segmented sequence.  Did you mean to do that?",
                       segset_import_help_str);
  }
  return rval;
}


static Boolean ImportFastaDialog (DialoG d, CharPtr filename)

{
  CharPtr       extension;
  FILE          *f;
  FastaPagePtr  fpp;
  ValNodePtr    head;
  Char          path [PATH_MAX];
  RecT          r;
  SeqEntryPtr   sep, new_sep_list, new_sep, test_sep;
  Boolean       rval = FALSE;
  BioseqPtr     bsp;
  CharPtr       supplied_id_txt = NULL;
  Boolean       chars_stripped = FALSE;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  fpp = (FastaPagePtr) GetObjectExtra (d);
  if (fpp != NULL) {
    if (fpp->list != NULL && fpp->single)
    {
      if (!fpp->is_na
          || fpp->seqPackagePtr == NULL
          || *fpp->seqPackagePtr == SEQ_PKG_GENOMICCDNA)
      {
        Message (MSG_ERROR, "Can't import additional sequences!");
        return FALSE;
      }
      else
      {
        if (Message (MSG_YN, "You are importing multiple sequences - did you intend to create a batch submission?") == ANS_NO)
        {
          Message (MSG_ERROR, "Can't import additional sequences!");
          return FALSE;
        }
        else
        {
          *(fpp->seqPackagePtr) = SEQ_PKG_GENBANK;
          fpp->single = FALSE;
          SafeHide (fpp->singleIdGrp);
        }
      }
    }
    extension = NULL;
    if (fpp->is_mrna) {
      extension = GetAppProperty ("FastaNucExtension");
    } else if (fpp->is_na) {
      extension = GetAppProperty ("FastaNucExtension");
    } else {
      extension = GetAppProperty ("FastaProtExtension");
    }
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
      WatchCursor ();
      StringCpy (fpp->path, path);
      ObjectRect (fpp->doc, &r);
      InsetRect (&r, 4, 4);
      faColFmt.pixWidth = r.right - r.left;
      /*
      ResetFastaPage (fpp);
      */
      Reset (fpp->doc);
      Update ();
      sep = fpp->list;
      head = fpp->errmsgs;
      f = FileOpen (fpp->path, "r");
      if (f == NULL)
      {
        Message (MSG_ERROR, "Unable to open %s", fpp->path);
        fpp->path[0] = 0;
      }
      else
      {
        if (fpp->singleSeqID != NULL)
        {
          supplied_id_txt = SaveStringFromText (fpp->singleSeqID);
        }
        new_sep_list = ImportSequencesFromFile (f, NULL, fpp->is_na, 
                                                fpp->parseSeqId, 
                                                supplied_id_txt,
                                                &head, &chars_stripped);
        if (chars_stripped && new_sep_list != NULL)
        {
          if (ANS_CANCEL == Message (MSG_OKC, "Illegal characters will be stripped from your sequence data.  Do you want to continue?"))
          {
            new_sep_list = SeqEntryFree (new_sep_list);
            FileClose (f);
            fpp->path [0] = 0;
            ArrowCursor ();
            Update ();
            return FALSE;
          }
        }
        supplied_id_txt = MemFree (supplied_id_txt);                                              
        FileClose (f);
        
        if (new_sep_list != NULL
            && new_sep_list->next == NULL
            && fpp->single
            && fpp->list == NULL
            && fpp->is_na
            && new_sep_list->choice == 1
            && new_sep_list->data.ptrvalue != NULL)
        {
          bsp = (BioseqPtr) new_sep_list->data.ptrvalue;
          
          /* assign a fake ID if there is only one sequence being imported, 
           * the package type is single, and there are no other sequences
           * from previous imports.
           */
          
          if (bsp->id == NULL)
          {
            bsp->id = MakeSeqID ("nuc_1");
          }
        }
      
        if (new_sep_list == NULL)
        {
          Message (MSG_ERROR, "Unable to read sequences from %s", fpp->path);
          fpp->path [0] = 0;
        }
        else if (! RejectZeroLengthSequences (&new_sep_list))
        {
          fpp->path [0] = 0;
        }
        else if (! RejectExtraSequences (new_sep_list, fpp))
        {
          /* if unsuccessful, delete new list */ 
          new_sep = new_sep_list;   
          while (new_sep != NULL)
          {
            test_sep = new_sep->next;
            SeqEntryFree (new_sep);
            new_sep = test_sep;
          }
          fpp->path [0] = 0;
        }
        else if (fpp->seqPackagePtr != NULL
                 && ! ImportedSequenceTypeOk (new_sep_list, *(fpp->seqPackagePtr)))
        {
          /* if unsuccessful, delete new list */ 
          new_sep = new_sep_list;   
          while (new_sep != NULL)
          {
            test_sep = new_sep->next;
            SeqEntryFree (new_sep);
            new_sep = test_sep;
          }
          fpp->path [0] = 0;
        }
        else if (CollectIDsAndTitles (new_sep_list, fpp->list, (fpp->is_na && ! fpp->is_mrna)))
        {
          if (fpp->is_na)
          {
            /* add default molecule type, topology, location, and genetic codes */
            AddDefaultModifierValues (new_sep_list);
          }
        
          /* if successful, link old and new lists */
          ValNodeLink (&(fpp->list), new_sep_list);
          rval = TRUE;
        }
        else
        {
          /* if unsuccessful, delete new list */ 
          new_sep = new_sep_list;   
          while (new_sep != NULL)
          {
            test_sep = new_sep->next;
            SeqEntryFree (new_sep);
            new_sep = test_sep;
          }
          fpp->path [0] = 0;
        }
      }
      
      if (fpp->list == NULL)
      {
        SafeHide (fpp->have_seq_instr_grp);
        Reset (fpp->doc);
        SafeShow (fpp->instructions);
        Update ();
        SetTitle (fpp->import_btn, "Import Nucleotide FASTA");
        Enable (fpp->import_btn);
        Disable (fpp->clear_btn);
      }
      else
      {        
        SafeHide (fpp->instructions);
        Update ();
        if (! fpp->is_na || fpp->single 
            || fpp->seqPackagePtr == NULL 
            || *fpp->seqPackagePtr == SEQ_PKG_GENOMICCDNA)
        {
          Disable (fpp->import_btn);
        }
        else
        {
          Enable (fpp->import_btn);
          SetTitle (fpp->import_btn, "Import Additional Nucleotide FASTA");
        }
        Enable (fpp->clear_btn);
        FormatFastaDoc (fpp);
        SafeShow (fpp->have_seq_instr_grp);
      }
      ArrowCursor ();
      Update ();
      return rval;
    }
  }
  return FALSE;
}

#define EXPORT_PAGE_WIDTH 80

static void ExportSeqIdAndTitle (SeqIdPtr sip, CharPtr title, FILE *fp)
{
  CharPtr id_str = NULL;

  if (fp == NULL)
  {
    return;
  }
  
  id_str [0] = 0;
  if (sip == NULL)
  {
    id_str = StringSave ("unknown_id");
  }
  else
  {
    id_str = SeqIdWholeLabel (sip, PRINTID_REPORT);
  }
  
  if (StringCSpn (id_str, " \t") == StringLen (id_str))
  {
    fprintf (fp, ">%s %s\n", id_str, title == NULL ? "" : title);
  }
  else
  {
    fprintf (fp, ">'%s' %s\n", id_str, title == NULL ? "" : title);
  }
  id_str = MemFree (id_str);
}

static void ExportSeqPort (Int4 from, Int4 to, SeqPortPtr spp, FILE *fp)
{
  Char        buffer [EXPORT_PAGE_WIDTH + 1];
  Int4        seq_offset, txt_out;

  if (spp == NULL || fp == NULL || from < 0 || to <= from)
  {
    return;
  }
  
  seq_offset = from;
  while (seq_offset < to)
  {
    txt_out = ReadBufferFromSep (spp, buffer, seq_offset, 
                                 MIN (seq_offset + EXPORT_PAGE_WIDTH, to), 0);
    if (txt_out == 0) break;
    seq_offset += txt_out;
    fprintf(fp, "%s\n", buffer);
  }
  
}

static void ExportOneRawSequence (BioseqPtr bsp, CharPtr title_master, FILE *fp)
{
  SeqDescrPtr sdp;
  Char        buffer [EXPORT_PAGE_WIDTH + 1];
  SeqPortPtr  spp;
  CharPtr     title = NULL;
  CharPtr     combined_title = NULL;
  Boolean     free_combined = FALSE;
  
  if (bsp == NULL || fp == NULL || bsp->repr != Seq_repr_raw)
  {
    return;
  }  
  
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp != NULL)
  {
    title = sdp->data.ptrvalue;
  }
  
  if (StringHasNoText (title_master))
  {
    combined_title = title;
  }
  else if (StringHasNoText (title))
  {
    combined_title = title_master;
  }
  else
  {
    combined_title = (CharPtr) MemNew ((StringLen (title_master) + StringLen (title) + 2) * sizeof (Char));
    if (combined_title != NULL)
    {
      StringCpy (combined_title, title_master);
      StringCat (combined_title, " ");
      StringCat (combined_title, title);
      free_combined = TRUE;
    }
  }
  
  ExportSeqIdAndTitle (bsp->id, combined_title, fp);
  if (free_combined)
  {
    combined_title = MemFree (combined_title);
  }

  buffer [EXPORT_PAGE_WIDTH] = 0;
  
  spp = SeqPortNew (bsp, 0, bsp->length-1, Seq_strand_plus, Seq_code_iupacna);
  
  ExportSeqPort (0, bsp->length, spp, fp);

  SeqPortFree (spp);
  fprintf (fp, "\n");  
}

static void ExportOneSegmentedBioseq (BioseqPtr bsp, FILE *fp)
{
  SeqLocPtr   slp;
  BioseqPtr   bsp_seg;
  SeqDescrPtr sdp;
  CharPtr     title = NULL;
  
  if (bsp == NULL || fp == NULL || bsp->repr != Seq_repr_seg)
  {
    return;
  }
  
  fprintf (fp, "[\n");
  
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp != NULL)
  {
    title = sdp->data.ptrvalue;
  }

  slp = (SeqLocPtr) bsp->seq_ext;
  while (slp != NULL)
  {
    bsp_seg = BioseqFind (SeqLocId (slp));
    ExportOneRawSequence (bsp_seg, title, fp);
    title = NULL;
    slp = slp->next;
  }
  fprintf (fp, "]\n\n");
}

static Boolean ExportOneDeltaBioseq (BioseqPtr bsp, FILE *fp)
{
  SeqDescrPtr sdp;
  CharPtr     title = NULL;
  DeltaSeqPtr dsp;
  SeqLitPtr   slip;
  SeqPortPtr  spp;
  Char        buffer [EXPORT_PAGE_WIDTH + 1];
  Int4        seq_offset;
  
  if (bsp == NULL || fp == NULL || bsp->repr != Seq_repr_delta
      || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL)
  {
    return FALSE;
  }

  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
    if (dsp->data.ptrvalue == NULL || dsp->choice != 2)
    {
      Message (MSG_ERROR, "Can't export badly formed delta sequence!");
      return FALSE;
    }
    dsp = dsp->next;
  }
  
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp != NULL)
  {
    title = sdp->data.ptrvalue;
  }

  ExportSeqIdAndTitle (bsp->id, title, fp);
  
  buffer [EXPORT_PAGE_WIDTH] = 0;
  
  spp = SeqPortNew (bsp, 0, bsp->length-1, Seq_strand_plus, Seq_code_iupacna);
  
  seq_offset = 0;
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL)
  {
		slip = (SeqLitPtr) (dsp->data.ptrvalue);
    if (IsDeltaSeqGap(dsp)) 
    {
      if (IsDeltaSeqUnknownGap (dsp))
      {
        fprintf (fp, ">?unk100\n");
      }
      else
      {
        fprintf (fp, ">?%d\n", slip->length);
      }
    }
    else
    {
      ExportSeqPort (seq_offset, seq_offset + slip->length, spp, fp);
		}
    seq_offset += slip->length;
    dsp = dsp->next;
  }
  fprintf (fp, "\n");
  return TRUE;
}

static void ExportFASTASeqEntryList (SeqEntryPtr sep, FILE *fp)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL || fp == NULL)
  {
    return;
  }
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_na (bsp->mol))
    {
      if (bsp->repr == Seq_repr_raw)
      {
        if (SeqMgrGetParentOfPart (bsp, NULL) == NULL)
        {
          ExportOneRawSequence (bsp, NULL, fp);
        }
      }
      else if (bsp->repr == Seq_repr_seg)
      {
        ExportOneSegmentedBioseq (bsp, fp);
      }
      else if (bsp->repr == Seq_repr_delta)
      {
        ExportOneDeltaBioseq (bsp, fp);    
      }
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    /* we don't export the parts set because we export them
     * when we do the master segment 
     */
    if (bssp->_class != BioseqseqSet_class_parts)
    {
      ExportFASTASeqEntryList (bssp->seq_set, fp);
    }
  }
  ExportFASTASeqEntryList (sep->next, fp);
}

static Boolean ExportNucleotideFASTADialog (DialoG d, CharPtr filename)
{
  CharPtr       extension;
  FILE          *f;
  FastaPagePtr  fpp;
  Char          path [PATH_MAX];
  Boolean       rval = FALSE;

  fpp = (FastaPagePtr) GetObjectExtra (d);
  if (fpp == NULL) {
    return FALSE;
  }
  
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));

  extension = NULL;
  if (fpp->is_mrna) {
    extension = GetAppProperty ("FastaNucExtension");
  } else if (fpp->is_na) {
    extension = GetAppProperty ("FastaNucExtension");
  } else {
    extension = GetAppProperty ("FastaProtExtension");
  }
  if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), extension)) {
    f = FileOpen (path, "w");
    if (f == NULL)
    {
      Message (MSG_ERROR, "Unable to open %s", path);
    }
    else
    {
      WatchCursor ();
      ExportFASTASeqEntryList (fpp->list, f);    
      FileClose (f);
      
      ArrowCursor ();
      Update ();
      rval = TRUE;
    }
  }
  return rval;
}

static void CleanupFastaDialog (GraphiC g, VoidPtr data)

{
  FastaPagePtr  fpp;

  fpp = (FastaPagePtr) data;
  if (fpp != NULL) {
    ResetFastaPage (fpp);
  }
  MemFree (data);
}

static CharPtr  fastaNucMsg = "\
\nClick on 'Import Nucleotide FASTA' to read a formatted FASTA file \
or 'Add/Modify Sequences' to create the file here.  The FASTA definition \
line must be in the following form:\n\n\
>SeqID [organism=scientific name]\n\n\
where the [ and ] brackets are actually in the text.\n\
Properly formatted modifiers and a title can also be included in the FASTA definition line.";


static CharPtr  fastaGenMsg = "\
\nPlease enter information about the genomic \
sequence in the spaces above.  Then click on either \
'Add/Modify Sequences' to create your sequences with the editor or \
'Import Genomic FASTA' to read a previously generated FASTA file that \
contains the sequence (which can be in segments).  The \
FASTA definition lines may be of the following form:\n\n\
>ID [organism=scientific name] [strain=name] [clone=name] title\n\n\
where the [ and ] brackets are actually in the text.";

static CharPtr  fastaMrnaMsg  = "\
\nPlease enter information about the transcript \
sequences in the spaces above.  Then click on \
'Import Transcript FASTA' to read a FASTA file that \
contains the sequence (which can be in segments).  The \
FASTA definition lines may be of the following form:\n\n\
>ID [gene=symbol] [mrna=name] title\n\n\
where the [ and ] brackets are actually in the text.";

static CharPtr  fastaProtMsg = "\
\nPlease enter information about the protein \
sequences in the spaces above.  Then click on \
'Import Protein FASTA' to read a FASTA file that \
contains the sequences.  The FASTA definition lines should \
be of the following form:\n\n\
>ID [gene=symbol] [protein=name] title\n\n\
where the [ and ] brackets are actually in the text.";

static CharPtr GetFastaSettingName (FastaPagePtr fpp)
{
  if (fpp == NULL)
  {
  	return NULL;
  }
  else if (fpp->is_mrna)
  {
    return "PARSEMRNASEQID";	
  }
  else if (fpp->is_na)
  {
    return "PARSENUCSEQID";
  }
  else
  {
    return "PARSEPROTSEQID";
  }
}

static void ChangeIDParse (ButtoN b)
{
  FastaPagePtr      fpp;
  CharPtr           setting_name;

  fpp = (FastaPagePtr) GetObjectExtra (b);
  if (fpp != NULL) {
    fpp->parseSeqId = GetStatus (b);
  
    setting_name = GetFastaSettingName (fpp);
  
    if (fpp->parseSeqId) {
      SetAppParam ("SEQUINCUSTOM", "PREFERENCES", setting_name, "TRUE");
      SafeHide (fpp->singleIdGrp);
    } else {
      SetAppParam ("SEQUINCUSTOM", "PREFERENCES", setting_name, "FALSE");
      if (fpp->single)
      {
        SafeShow (fpp->singleIdGrp);
      }
      else
      {
      	SafeHide (fpp->singleIdGrp);
      }
    }
  }
}

extern DialoG CreateFastaDialog (GrouP h, CharPtr title,
                                 Boolean is_na, Boolean is_mrna, CharPtr text,
                                 Boolean single, Int2Ptr seqPackagePtr)

{
  FastaPagePtr  fpp;
  GrouP         g;
  GrouP         m;
  GrouP         p;
  GrouP         s;
  PrompT        pr;
  CharPtr       setting_name;
  ButtoN        prs = NULL;
  Char          str [32];
  Boolean       parseSeqId;
#ifdef WIN_MAC
  Int2          wid = 25;
#else
  Int2          wid = 33;
#endif

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  fpp = (FastaPagePtr) MemNew (sizeof (FastaPage));
  if (fpp != NULL) {

    SetObjectExtra (p, fpp, CleanupFastaDialog);
    fpp->dialog = (DialoG) p;
    fpp->todialog = NULL;
    fpp->fromdialog = NULL;
    fpp->importdialog = ImportFastaDialog;
    if (is_na)
    {
      fpp->exportdialog = ExportNucleotideFASTADialog;
    }
    else
    {
      fpp->exportdialog = NULL;
    }

    fpp->seqPackagePtr = seqPackagePtr;
    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);

    fpp->path [0] = '\0';
    fpp->is_na = is_na;
    fpp->is_mrna = is_mrna;
    fpp->single = single;

    setting_name = GetFastaSettingName (fpp);
    
    if (GetAppParam ("SEQUINCUSTOM", "SETTINGS", "ALLOWNOSEQID", NULL, str, sizeof (str))
        && StringICmp (str, "TRUE") == 0)
    {
      prs = CheckBox (m, "Fasta definition line starts with sequence ID", ChangeIDParse);
      SetObjectExtra (prs, fpp, NULL);
    }
    parseSeqId = FALSE;
    if (GetAppParam ("SEQUINCUSTOM", "PREFERENCES", setting_name, NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        parseSeqId = TRUE;
      }
    }
    else
    {
      parseSeqId = TRUE;
    }
    SetStatus (prs, parseSeqId);
    
    fpp->parseSeqId = parseSeqId;
    if (fpp->single) {
      fpp->singleIdGrp = HiddenGroup (m, 2, 0, NULL);
      StaticPrompt (fpp->singleIdGrp, "Enter unique identifier for this sequence", 0, dialogTextHeight, programFont, 'l');
      fpp->singleSeqID = DialogText (fpp->singleIdGrp, "", 6, NULL);
      if (parseSeqId) {
        Hide (fpp->singleIdGrp);
      }
    }

    g = HiddenGroup (m, 0, 0, NULL);
    fpp->instructions = MultiLinePrompt (g, text, 27 * stdCharWidth, programFont);
    fpp->have_seq_instr_grp = HiddenGroup (g, -1, 0, NULL);
    SetGroupSpacing (fpp->have_seq_instr_grp, 10, 10);
    fpp->doc = DocumentPanel (fpp->have_seq_instr_grp, stdCharWidth * wid, stdLineHeight * 12);
    SetDocAutoAdjust (fpp->doc, FALSE);
    pr = StaticPrompt (fpp->have_seq_instr_grp, "Choose Clear from the Edit menu to clear these sequences", 0, dialogTextHeight, systemFont, 'c');
    AlignObjects (ALIGN_CENTER, (HANDLE) fpp->doc, (HANDLE) pr, NULL);
    Hide (fpp->have_seq_instr_grp);
    AlignObjects (ALIGN_CENTER, (HANDLE) fpp->instructions,
                  (HANDLE) fpp->have_seq_instr_grp, NULL);
                  
    AlignObjects (ALIGN_CENTER, (HANDLE) g,
                                (HANDLE) prs,
                                (HANDLE) fpp->singleIdGrp,
                                NULL);                  
  }

  return (DialoG) p;
}

typedef struct phylippage {
  DIALOG_MESSAGE_BLOCK
  Uint1        format;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  ValNodePtr   errmsgs;
  DoC          doc;
  GrouP        instructions;
  Char         extension [10];
  Int4         type;
  TSequenceInfoPtr aln_settings;

} PhylipPage, PNTR PhylipPagePtr;


#define PhylipFormatBufLen 1000

static void FormatPhylipDoc (PhylipPagePtr ppp)

{
  Nlm_QualNameAssocPtr ap;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  CharPtr            label;
  Int4               len;
  CharPtr            measure;
  SeqEntryPtr        nsep;
  Int2               num;
  CharPtr            plural;
  SeqIdPtr           sip;
  SeqEntryPtr        sep;
  CharPtr            str;
  CharPtr            title;
  CharPtr            ttl;
  CharPtr            tmp;
  CharPtr            valstr;
  ValNodePtr         vnp;

  if (ppp != NULL) {
    str = MemNew (sizeof (char) * PhylipFormatBufLen);
    tmp = MemNew (sizeof (char) * PhylipFormatBufLen);
    if (str == NULL || tmp == NULL) return;
    num = 0;
    len = 0;
    sep = ppp->sep;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && (bssp->_class == 7 ||
                           (IsPopPhyEtcSet (bssp->_class)))) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          num++;
          if (IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              len += bsp->length;
            }
          } else if (IS_Bioseq_set (sep)) {
            nsep = FindNucSeqEntry (sep);
            if (nsep != NULL && IS_Bioseq (nsep)) {
              bsp = (BioseqPtr) nsep->data.ptrvalue;
              if (bsp != NULL) {
                len += bsp->length;
              }
            }
          }
        }
      }
    }
    if (num > 1) {
      plural = "s";
    } else {
      plural = "";
    }
    label = "Sequence";
    measure = "nucleotides";
    sprintf (str, "%d nucleotide sequence%s, total length %ld %s\n",
             (int) num, plural, (long) len, measure);
    AppendText (ppp->doc, str, &faParFmt, &faColFmt, programFont);
    vnp = ppp->errmsgs;
    num = 0;
    sep = ppp->sep;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && (bssp->_class == 7 ||
                           (IsPopPhyEtcSet (bssp->_class)))) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          nsep = NULL;
          num++;
          len = 0;
          sip = NULL;
          tmp [0] = '\0';
          if (IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              len = bsp->length;
              sip = SeqIdFindWorst (bsp->id);
              SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
            }
          } else if (IS_Bioseq_set (sep)) {
            nsep = FindNucSeqEntry (sep);
            if (nsep != NULL && IS_Bioseq (nsep)) {
              bsp = (BioseqPtr) nsep->data.ptrvalue;
              if (bsp != NULL) {
                len = bsp->length;
                sip = SeqIdFindWorst (bsp->id);
                SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
              }
            }
          }
          sprintf (str, "\n%s %d\nLength: %ld %s\nSequence ID: %s\n", label,
                   (int) num, (long) len, measure, tmp);
          ttl = NULL;
          SeqEntryExplore (nsep, (Pointer) (&ttl), FindFirstTitle);
          title = StringSaveNoNull (ttl);
          if (title != NULL) {
            valstr = FindValueFromPairInDefline ("organism", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "Organism", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("organism", title);
            
            valstr = FindValueFromPairInDefline ("lineage", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "Lineage", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("lineage", title);

            for (ap = current_orgmod_subtype_alist; ap->name != NULL; ap++) {
              if (IsNonTextModifier (ap->name))
              {
                if (FindValuePairInDefLine (ap->name, title, NULL) != NULL)
                {
                  AddReportLine (str, ap->name, "TRUE");
                  RemoveValueFromDefline (ap->name, title);
                }
              }
              else
              {
                valstr = FindValueFromPairInDefline (ap->name, title);
                if (!StringHasNoText (valstr)) {
                  AddReportLine (str, ap->name, title);  
                }
                valstr = MemFree (valstr);
                RemoveValueFromDefline (ap->name, title);
              }
            }
            for (ap = current_subsource_subtype_alist; ap->name != NULL; ap++) {
              if (IsNonTextModifier (ap->name))
              {
                if (FindValuePairInDefLine (ap->name, title, NULL) != NULL)
                {
                  AddReportLine (str, ap->name, "TRUE");
                  RemoveValueFromDefline (ap->name, title);
                }
              }
              else
              {
                valstr = FindValueFromPairInDefline (ap->name, title);
                if (!StringHasNoText (valstr)) {
                  AddReportLine (str, ap->name, title);  
                }
                valstr = MemFree (valstr);
                RemoveValueFromDefline (ap->name, title);
              }
            }
            
            valstr = FindValueFromPairInDefline ("note-orgmod", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "Note", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("note-orgmod", title);
            
            valstr = FindValueFromPairInDefline ("note-subsrc", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "Note", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("note-subsrc", title);
            
            valstr = FindValueFromPairInDefline ("molecule", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "Molecule", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("molecule", title);
            
            valstr = FindValueFromPairInDefline ("moltype", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "MolType", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("moltype", title);
            
            valstr = FindValueFromPairInDefline ("location", title);
            if (!StringHasNoText (valstr)) {
              AddReportLine (str, "Location", valstr);
            }
            valstr = MemFree (valstr);
            RemoveValueFromDefline ("location", valstr);

            TrimSpacesAroundString (title);
            if (! StringHasNoText (title)) {
              StringCat (str, "Title: ");
              StringNCat (str, title, 128);
              StringCat (str, "\n");
            } else {
              StringCat (str, "No title detected\n");
            }
          }
          MemFree (title);
          if (vnp != NULL && vnp->data.ptrvalue != NULL) {
            StringCat (str, (CharPtr) vnp->data.ptrvalue);
            StringCat (str, "\n");
          }
          AppendText (ppp->doc, str, &faParFmt, &faColFmt, programFont);
          if (vnp != NULL) {
            vnp = vnp->next;
          }
        }
      }
    }
    MemFree (str);
    MemFree (tmp);
    UpdateDocument (ppp->doc, 0, 0);
  }
}

static void ResetPhylipPage (PhylipPagePtr ppp)

{
  if (ppp != NULL) {
    ppp->sep = SeqEntryFree (ppp->sep);
    ppp->errmsgs = ValNodeFreeData (ppp->errmsgs);
  }
}

static CharPtr noOrgInTitleWarning =
"sequences have organism information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Please quit Sequin and read the Sequin Quick Guide section on preparing the data files before proceeding.";

static void CountTitlesWithoutOrganisms (SeqEntryPtr sep)
{
  IDAndTitleEditPtr iatep;
  Int4              seq_num;
  CharPtr           org_name;
  Int4              num_sequences = 0, num_with_orgs = 0;
  
  iatep = SeqEntryListToIDAndTitleEdit (sep);
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    num_sequences ++;
    org_name = FindValueFromPairInDefline ("organism", iatep->title_list [seq_num]);
    if (!StringHasNoText (org_name))
    {
      num_with_orgs ++;
    }
    org_name = MemFree (org_name);
  }
  iatep = IDAndTitleEditFree (iatep);
  if (num_sequences != num_with_orgs && num_with_orgs != 0)
  {
    Message (MSG_OK, "%d of %d %s", num_sequences - num_with_orgs, (int) num_sequences, noOrgInTitleWarning);
  }
  
}

static CharPtr  phylipNucMsg = "\
\nClick 'Import Nucleotide Alignment' to load your \
nucleotide alignment file.\n\nClick on 'Custom Alignment Settings' \
if Sequin has trouble reading your alignment file.";

static void SetPhylipDocInstructions (PhylipPagePtr ppp)
{
  if (ppp == NULL || ppp->doc == NULL) return;
  Reset (ppp->doc);
  AppendText (ppp->doc, phylipNucMsg, &faParFmt, &faColFmt, programFont);
  UpdateDocument (ppp->doc, 0, 0);
  Update ();
}

static Boolean ImportPhylipDialog (DialoG d, CharPtr filename)
{
  Char           path [PATH_MAX];
  PhylipPagePtr  ppp;
  SeqEntryPtr    sep;
  RecT           r;
  FILE           *fp;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  Char           errStr [PATH_MAX + 64];
  CharPtr        no_org_err_msg = NULL;

  if (d == NULL || filename == NULL) return FALSE;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  ppp = (PhylipPagePtr) GetObjectExtra (d);
  if (ppp == NULL) {
    return FALSE;
  }

  if (path [0] != '\0' || GetInputFileName (path, sizeof (path), ppp->extension, "TEXT")) {
    WatchCursor ();
    StringCpy (ppp->path, path);
    ObjectRect (ppp->doc, &r);
    InsetRect (&r, 4, 4);
    faColFmt.pixWidth = r.right - r.left;
    Reset (ppp->doc);
    Update ();
    ppp->sep = SeqEntryFree (ppp->sep);
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      ppp->sep = SeqEntryFromAlignmentFile (fp, ppp->aln_settings,
                                            Seq_mol_na, no_org_err_msg);
                                            
      /* check for bracketing issues here */
      if (CollectIDsAndTitles (ppp->sep, NULL, TRUE))
      {
        /* add default molecule type, topology, location, and genetic codes */
        AddDefaultModifierValues (ppp->sep);
      }        
      else
      {
        ppp->sep = SeqEntryFree (ppp->sep);
      }
                                                  
      sep = ppp->sep;
      if (sep != NULL) {
        SaveSeqEntryObjMgrData (ppp->sep, &omdptop, &omdata);
        GetSeqEntryParent (ppp->sep, &parentptr, &parenttype);
        SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
        RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);

        FormatPhylipDoc (ppp);
        SafeShow (ppp->doc);

        CountTitlesWithoutOrganisms (sep);
      } else {
        SendHelpScrollMessage (helpForm, "Nucleotide Page", "Nucleotide Page for Aligned Data Formats");
        SetPhylipDocInstructions (ppp);
      }
    } else {
      SetPhylipDocInstructions (ppp);
    }
  } else {
	sprintf (errStr, "ERROR: Unable to open file %s\n\n", path);
	AppendText (ppp->doc, errStr, &faParFmt, &faColFmt, programFont);
	AppendText (ppp->doc, strerror(errno), &faParFmt, &faColFmt, programFont);
	SafeShow (ppp->doc);
    Update ();
  }
  ArrowCursor ();
  Update ();
  return TRUE;
}

static void CleanupPhylipDialog (GraphiC g, VoidPtr data)

{
  PhylipPagePtr  ppp;

  ppp = (PhylipPagePtr) data;
  if (ppp != NULL) {
    ResetPhylipPage (ppp);
    SequenceInfoFree (ppp->aln_settings);
    ppp->aln_settings = NULL;
  }
  MemFree (data);
}


static DialoG CreatePhylipDialog (GrouP h, CharPtr title, CharPtr text,
                                  Int2 format, CharPtr extension,
                                  Int4 type)

{
  PhylipPagePtr  ppp;
  GrouP          g;
  GrouP          m;
  GrouP          p;
  GrouP          s;
  RecT          r;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ppp = (PhylipPagePtr) MemNew (sizeof (PhylipPage));
  if (ppp != NULL) {

    SetObjectExtra (p, ppp, CleanupPhylipDialog);
    ppp->dialog = (DialoG) p;
    ppp->todialog = NULL;
    ppp->fromdialog = NULL;
    ppp->importdialog = ImportPhylipDialog;
    ppp->type = type;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);

    ppp->format = (Uint1) format;
    ppp->path [0] = '\0';
    StringNCpy_0 (ppp->extension, extension, sizeof (ppp->extension));
  
    g = HiddenGroup (m, 0, 0, NULL);
    ppp->doc = DocumentPanel (g, stdCharWidth * 27, stdLineHeight * 8);
    ObjectRect (ppp->doc, &r);
    InsetRect (&r, 4, 4);
    faColFmt.pixWidth = r.right - r.left;

    ppp->aln_settings = GetDefaultSequenceInfo();

    SetPhylipDocInstructions (ppp);
  }

  return (DialoG) p;
}

#define NUCLEOTIDE_PAGE        0
#define SEQUENCING_METHOD_PAGE 1
#define ORGANISM_PAGE          2
#define MRNA_PAGE              3
#define PROTEIN_PAGE           4
#define ANNOTATE_PAGE          5

/*---------------------------------------------------------------------*/
/*                                                                     */
/* HasZeroLengthSequence () -- Checks to see if any of a submission's  */
/*                             sequences are missing (ie -- zero       */
/*                             length).                                */
/*                                                                     */
/*---------------------------------------------------------------------*/

extern Boolean HasZeroLengthSequence (ForM newForm)
{
  SequencesFormPtr  sqfp;
  FastaPagePtr      fpp;
  SeqEntryPtr       sep;
  BioseqPtr         bsp;

  /* Get the list of Bioseqs to check */

  sqfp = (SequencesFormPtr) GetObjectExtra (newForm);
  if (NULL == sqfp)
    return TRUE;

  fpp = GetObjectExtra (sqfp->dnaseq);
  sep = fpp->list;

  /* Check the list */

  while (NULL != sep) {
    if (sep->choice == 1) { 
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp->length <= 0)
	return TRUE;
    }
    sep = sep->next;
  }

  /* If we made it to here, then */
  /* there were none found.      */

  return FALSE;
}

extern Boolean SequencesFormHasProteins (ForM f)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    if (PackageTypeIsSet (sqfp->seqPackage)) return TRUE;
    fpp = GetObjectExtra (sqfp->protseq);
    if (fpp != NULL) {
      if (fpp->path [0] != '\0') {
        return TRUE;
      }
    }
  }
  return FALSE;
}

extern SeqEntryPtr GetSequencesFormProteinList (ForM f)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    fpp = GetObjectExtra (sqfp->protseq);
    if (fpp != NULL) {
      return fpp->list;
    }
  }
  return NULL;
}

static SeqEntryPtr GetSeqEntryFromSequencesForm (SequencesFormPtr sqfp)
{
  SeqEntryPtr list = NULL;
  FastaPagePtr       fpp;
  PhylipPagePtr      ppp;
  SeqEntryPtr        sep;
  BioseqSetPtr       bssp;
  
  if (sqfp == NULL) return NULL;

  if (sqfp->seqFormat == SEQ_FMT_FASTA) {
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) 
    {
      list = fpp->list;
    }
  } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
    ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (ppp != NULL) {
      sep = ppp->sep;
      if (sep != NULL && IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL) {
          list = bssp->seq_set;
        }
      }
    }
  }
  return list;
}

extern SeqEntryPtr GetSequencesFormNucleotideList (ForM f)
{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    return GetSeqEntryFromSequencesForm (sqfp);
  }
  return NULL;
}

extern Boolean SequencesFormHasTooManyNucleotides (ForM f)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL && PackageTypeIsSingle (sqfp->seqPackage))
  {
    fpp = GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      if (fpp->list != NULL && fpp->list->next != NULL) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

static ValNodePtr 
BuildModifierTypeList 
(ValNodePtr type_list,
 CharPtr    new_title, 
 Boolean    allow_prot)
{
  ValNodePtr      modifier_info_list;
  ValNodePtr      info_vnp, type_vnp;
  ModifierInfoPtr mip;
  
  modifier_info_list = ParseAllBracketedModifiers (new_title);
  for (info_vnp = modifier_info_list; info_vnp != NULL; info_vnp = info_vnp->next)
  {
    mip = (ModifierInfoPtr)info_vnp->data.ptrvalue;
    if (mip == NULL 
        || mip->modtype == eModifierType_Protein
        || mip->modtype == eModifierType_Organism)
    {
      continue;
    }
    if (mip->modtype == eModifierType_SourceQual)
    {
  	  for (type_vnp = type_list;
  	       type_vnp != NULL 
  	         && (type_vnp->choice != mip->subtype 
  	             || StringICmp (type_vnp->data.ptrvalue, mip->name) != 0); 
  	       type_vnp = type_vnp->next)
  	  {
  	  }
    }
    else
    {
  	  for (type_vnp = type_list;
  	       type_vnp != NULL && StringICmp (type_vnp->data.ptrvalue, mip->name) != 0;
  	       type_vnp = type_vnp->next)
  	  {
  	  }
    }
  	if (type_vnp == NULL)
  	{
  	  type_vnp = ValNodeNew (type_list);
  	  if (type_list == NULL) type_list = type_vnp;
  	  if (type_vnp != NULL)
  	  {
  	  	type_vnp->choice = mip->subtype;
  	  	type_vnp->data.ptrvalue = StringSave (mip->name);
  	  }
  	}
  }
  ModifierInfoListFree (modifier_info_list);
  return type_list;
}


static Uint2 modedit_widths [] = {
  0, 0,
};

ENUM_ALIST(nontextmodedit_alist)
  {"FALSE",             0},
  {"TRUE",              1},
END_ENUM_ALIST

extern void ConfirmSequencesFormParsing (ForM f, FormActnFunc putItAllTogether)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL && putItAllTogether != NULL) {
    putItAllTogether (sqfp->form);
  }
}

extern void AddToSubSource (BioSourcePtr biop, CharPtr title, CharPtr label, Uint1 subtype)

{
  CharPtr       ptr;
  SubSourcePtr  ssp;
  CharPtr       str;
  SubSourcePtr  tmpssp;

  if (biop == NULL || title == NULL || label == NULL) return;
  str = MemNew (StringLen (title));
  if (str == NULL) return;
  ptr = StringISearch (title, label);
  if (ptr != NULL) {
    StringCpy (str, ptr + StringLen (label));
    ptr = StringChr (str, ']');
    if (ptr != NULL) {
      *ptr = '\0';
      TrimSpacesAroundString (str);
      ssp = SubSourceNew ();
      if (biop->subtype == NULL) {
        biop->subtype = ssp;
      } else {
        tmpssp = biop->subtype;
        while (tmpssp->next != NULL) {
          tmpssp = tmpssp->next;
        }
        tmpssp->next = ssp;
      }
      if (ssp != NULL) {
        ssp->subtype = subtype;
        ssp->name = StringSave (str);
      }
    }
  }
  MemFree (str);
}

extern void AddToOrgMod (BioSourcePtr biop, CharPtr title, CharPtr label, Uint1 subtype)

{
  OrgModPtr   mod;
  OrgNamePtr  onp;
  OrgRefPtr   orp;
  CharPtr     ptr;
  CharPtr     str;
  OrgModPtr   tmpmod;

  if (biop == NULL || title == NULL || label == NULL) return;
  str = MemNew (StringLen (title));
  if (str == NULL) return;
  ptr = StringISearch (title, label);
  if (ptr != NULL) {
    StringCpy (str, ptr + StringLen (label));
    ptr = StringChr (str, ']');
    if (ptr != NULL) {
      *ptr = '\0';
      TrimSpacesAroundString (str);
      orp = biop->org;
      if (orp == NULL) {
        orp = OrgRefNew ();
        biop->org = orp;
      }
      if (orp != NULL) {
        onp = orp->orgname;
        if (onp == NULL) {
          onp = OrgNameNew ();
          orp->orgname = onp;
        }
        if (onp != NULL) {
          mod = OrgModNew ();
          if (onp->mod == NULL) {
            onp->mod = mod;
          } else {
            tmpmod = onp->mod;
            while (tmpmod->next != NULL) {
              tmpmod = tmpmod->next;
            }
            tmpmod->next = mod;
          }
          if (mod != NULL) {
            mod->subtype = subtype;
            mod->subname = StringSave (str);
          }
        }
      }
    }
  }
  MemFree (str);
}

#define PROC_NUC_STR_SIZE 4096

static Int4 TopologyFromString (CharPtr str)
{
  EnumFieldAssocPtr  eap;

  for (eap = topology_nuc_alist; eap != NULL && eap->name != NULL; eap++)
  {
    if (StringICmp (eap->name, str) == 0)
    {
      return eap->value;
    }
  }
  return 1; 
}

static BioSourcePtr AddOrgRef (BioSourcePtr biop)
{
  if (biop == NULL)
  {
    biop = BioSourceNew ();
  }
  if (biop == NULL)
  {
    return NULL;
  }
  if (biop->org == NULL)
  {
    biop->org = OrgRefNew ();
  }
  if (biop->org == NULL)
  {
    biop = BioSourceFree (biop);
    return NULL;
  }
  return biop;
}

static BioSourcePtr AddOrgName (BioSourcePtr biop)
{
  biop = AddOrgRef (biop);
  if (biop == NULL || biop->org == NULL)
  {
    biop = BioSourceFree (biop);
    return NULL; 
  }
  if (biop->org->orgname == NULL)
  {
    biop->org->orgname = OrgNameNew ();
    if (biop->org->orgname == NULL)
    {
      biop = BioSourceFree (biop);
      return NULL;
    }
  }
  return biop;
}

static BioSourcePtr SetGeneticCodeForBioSource (BioSourcePtr biop, Int4 gcode, Boolean is_nuc)
{
  OrgRefPtr  orp;
  OrgNamePtr onp;

  if (gcode < 0)
  {
    return biop;
  }

  biop = AddOrgName (biop);
  if (biop == NULL)
  {
    return biop;
  }
  
  orp = biop->org;
  if (biop->org == NULL)
  {
    biop->org = OrgRefNew ();
    orp = biop->org;
  }
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp == NULL) {
      onp = OrgNameNew ();
      orp->orgname = onp;
    }
    if (onp != NULL) {
      if (is_nuc)
      {
        onp->gcode = gcode;
      }
      else
      {
        onp->mgcode = gcode;
      }
    }
  }
  return biop;
}

static BioSourcePtr
SetGeneticCodeFromTitle 
(BioSourcePtr biop,
 CharPtr      title, 
 CharPtr      mod_name, 
 Boolean      is_nuc)
{
  CharPtr    gcode_str;
  Int4       gcode;
  CharPtr    next_org_loc;
  
  if (StringHasNoText (title))
  {
    return biop;
  }
  
  next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  gcode_str = FindValueFromPairInDeflineBeforeCharPtr (mod_name, title, next_org_loc);
  if (!StringHasNoText (gcode_str))
  {
    gcode = GeneticCodeFromString (gcode_str);  
    biop = SetGeneticCodeForBioSource (biop, gcode, is_nuc);     
  }
  if (gcode_str != NULL)
  {
    RemoveValueFromDefline (mod_name, title);
  }
  gcode_str = MemFree (gcode_str);
  return biop;
}

static BioSourcePtr 
SetAllGeneticCodesFromTitle 
(BioSourcePtr biop,
 CharPtr      title)
{
  Int4    code_to_use;
  CharPtr location;
  CharPtr next_org_loc;
  
  if (StringHasNoText (title))
  {
    return biop;
  }
  
  next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  location = FindValueFromPairInDeflineBeforeCharPtr ("location", title, next_org_loc);
  if (!StringHasNoText (location)) 
  {      
    code_to_use = UseGeneticCodeForLocation (location);
    if (code_to_use == USE_OTHER_GENETIC_CODE)
    {
      biop = SetGeneticCodeForBioSource (biop, 11, TRUE); 
      RemoveValueFromDefline ("genetic_code", title);    
    }
    else if (code_to_use == USE_NUCLEAR_GENETIC_CODE)
    {
      biop = SetGeneticCodeFromTitle (biop, title, "genetic_code", TRUE);
    }
    else if (code_to_use == USE_MITOCHONDRIAL_GENETIC_CODE)
    {
      biop = SetGeneticCodeFromTitle (biop, title, "genetic_code", FALSE);
    }
  }
  location = MemFree (location);
  
  biop = SetGeneticCodeFromTitle (biop, title, "gcode", TRUE);
  biop = SetGeneticCodeFromTitle (biop, title, "mgcode", FALSE);

  return biop;
}

static void 
SetMoleculeAndMolTypeFromTitle 
(BioseqPtr   bsp, 
 CharPtr     title,
 Int2        seqPackage)
{
  SeqEntryPtr sep;
  ValNodePtr vnp;
  MolInfoPtr mip = NULL;
  Uint1      biomol;
  Int4       molecule;
  CharPtr    valstr;
  CharPtr    ptr;
  SeqLocPtr  slp;
  BioseqPtr  bsp_seg;
  
  if (bsp == NULL)
  {
    return;
  }
  
  sep = SeqMgrGetSeqEntryForData (bsp); 
  if (sep == NULL)
  {
    return;
  }
    
  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  if (vnp == NULL)
  {
    if (seqPackage == SEQ_PKG_SINGLE)
    {
      biomol = 3;
      molecule = Seq_mol_rna;
    }
    else 
    {
      biomol = 1;
      molecule = Seq_mol_dna;
    }
  }
  else
  {
    mip = (MolInfoPtr) vnp->data.ptrvalue;
    biomol = mip->biomol;
    molecule = bsp->mol;
  }
  
  /* get moltype from defline */
  valstr = FindValueFromPairInDefline ("moltype", title);
  if (!StringHasNoText (valstr))
  {
    biomol = MolTypeFromString (valstr);
    if (biomol == 1)
    {
      molecule = Seq_mol_na;
    }
    else if (biomol >= 2 && biomol <= 7)
    {
      molecule = Seq_mol_rna;
    }
    else if (biomol == 9)
    {
      molecule = Seq_mol_dna;
    }
    else if (biomol == MOLECULE_TYPE_NCRNA)
    {
      molecule = Seq_mol_rna;
    }
    else if (biomol == 253)
    {
      molecule = Seq_mol_dna;
      biomol = 1;
    }
    else if (biomol == 254)
    {
      molecule = Seq_mol_rna;
      biomol = 1;
    }
    else if (biomol == 255)
    {
      molecule = Seq_mol_other;
    }
  }
  valstr = MemFree (valstr);
  
  RemoveValueFromDefline ("moltype", title);

  /* get molecule from defline */ 
  valstr = FindValueFromPairInDefline ("molecule", title);
  if (!StringHasNoText (valstr))
  {
    if (StringICmp (valstr, "dna") == 0) {
      molecule = Seq_mol_dna;
    } else if (StringICmp (valstr, "rna") == 0) {
      molecule = Seq_mol_rna;
    }
  }
  valstr = MemFree (valstr);
  RemoveValueFromDefline ("molecule", title);
  
  ptr = StringISearch (title, "[dna]");
  if (ptr != NULL)
  {
    molecule = Seq_mol_dna;
    ExciseString (title, "[dna", "]");
  }
  
  ptr = StringISearch (title, "[rna]");
  if (ptr != NULL)
  {
    molecule = Seq_mol_rna;
    ExciseString (title, "[rna", "]");
  }
  
  if (mip == NULL)
  {
    vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    mip = MolInfoNew ();
    vnp->data.ptrvalue = mip;
  }
     
  mip->biomol = biomol;
  bsp->mol = molecule;  

  valstr = FindValueFromPairInDefline ("tech", title);
  if (!StringHasNoText (valstr))
  {
    ReadTechFromString (valstr, mip);
  }
  valstr = MemFree (valstr);
  RemoveValueFromDefline ("tech", title);
  
  if (bsp->repr == Seq_repr_seg)
  {
    slp = (SeqLocPtr) bsp->seq_ext;
    while (slp != NULL)
    {
      bsp_seg = BioseqFind (SeqLocId (slp));
      sep = SeqMgrGetSeqEntryForData (bsp_seg);
      if (bsp_seg != NULL)
      {
        bsp_seg->mol = bsp->mol;
      }
      vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
      if (vnp == NULL)
      {
        vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      }
      if (vnp != NULL)
      {
        vnp->data.ptrvalue = MolInfoFree (vnp->data.ptrvalue);
        vnp->data.ptrvalue = (MolInfoPtr) AsnIoMemCopy (mip, (AsnReadFunc) MolInfoAsnRead,
                                                            (AsnWriteFunc) MolInfoAsnWrite);
      }
      slp = slp->next;
    }
  }
}

static void AddGeneticCodeComment (BioseqPtr bsp, CharPtr comment)
{
  SeqDescPtr         sdp;
  UserObjectPtr      uop = NULL;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp, last_ufp = NULL;
  CharPtr            comment_fmt = "Submitter genetic code: %s";
  CharPtr            new_comment;
  Int4               new_comment_len;

  if (bsp == NULL || StringHasNoText (comment))
  {
    return;
  }
  
  sdp = bsp->descr;
  while (sdp != NULL && uop == NULL)
  {
    if (sdp->choice == Seq_descr_user && sdp->data.ptrvalue != NULL)
    {
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      oip = uop->type;
      if (oip == NULL || StringCmp (oip->str, "Submission") != 0)
      {
        uop = NULL;
      }
    }
    sdp = sdp->next;
  }
  
  
  if (uop == NULL)
  {
    uop = UserObjectNew ();
    if (uop == NULL)
    {
      return;
    }
    uop->type = ObjectIdNew ();
    uop->type->str = StringSave ("Submission");
    ValNodeAddPointer (&bsp->descr, Seq_descr_user, uop);  
  }
  
  ufp = uop->data;
  while (ufp != NULL 
         && (ufp->label == NULL 
           || StringCmp (ufp->label->str, "AdditionalComment") != 0))
  {
    last_ufp = ufp;
    ufp = ufp->next;
  }
  
  if (ufp == NULL)
  {
    ufp = UserFieldNew ();
    ufp->label = ObjectIdNew ();
    ufp->label->str = StringSave ("AdditionalComment");
    if (last_ufp == NULL)
    {
      uop->data = ufp;
    }
    else
    {
      last_ufp->next = ufp;
    }
  }
  
  new_comment_len = StringLen (comment) + StringLen (comment_fmt);
  if (!StringHasNoText (ufp->data.ptrvalue))
  {
    new_comment_len += StringLen (ufp->data.ptrvalue);
  }
  new_comment = (CharPtr) MemNew (new_comment_len * sizeof (Char));
  sprintf (new_comment, comment_fmt, comment);
  
  if (!StringHasNoText (ufp->data.ptrvalue))
  {
    StringCat (new_comment, ufp->data.ptrvalue);
  }
  
  ufp->data.ptrvalue = MemFree (ufp->data.ptrvalue);
  ufp->data.ptrvalue = new_comment;
}

static BioSourcePtr AddOrgModValue (BioSourcePtr biop, Uint1 subtype, CharPtr subname)
{
  OrgModPtr    mod;
  
  if (subname == NULL)
  {
    return biop;
  }

  biop = AddOrgName (biop);
  if (biop != NULL)
  {
    mod = OrgModNew ();
    if (mod != NULL)
    {
      mod->subtype = subtype;
      mod->subname = subname;
      subname = NULL;
      mod->next = biop->org->orgname->mod;
      biop->org->orgname->mod = mod;
    }
  }
  subname = MemFree (subname);
  return biop;
}

static BioSourcePtr AddSubSourceValue (BioSourcePtr biop, Uint1 subtype, CharPtr subname)
{
  SubSourcePtr ssp;
  
  if (subname == NULL)
  {
    return biop;
  }

  if (biop == NULL)
  {
    biop = BioSourceNew ();
  }
  if (biop != NULL)
  {
    ssp = SubSourceNew ();
    if (ssp != NULL)
    {
      ssp->subtype = subtype;
      ssp->name = subname;
      subname = NULL;
      ssp->next = biop->subtype;
      biop->subtype = ssp;
    }
  }
  subname = MemFree (subname);
  return biop;
}

extern BioSourcePtr 
ExtractFromTitleToBioSourceOrgMod 
(CharPtr      title,
 BioSourcePtr biop, 
 CharPtr      mod_name,
 Int4         subtype)
{
  CharPtr valstr;
  CharPtr next_org_loc;
  
  next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  while ((valstr = FindValueFromPairInDeflineBeforeCharPtr (mod_name, title, next_org_loc)) != NULL)
  {
    biop = AddOrgModValue (biop, subtype, valstr);
    RemoveValueFromDefline (mod_name, title);
    next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  }  
  return biop;
}

extern BioSourcePtr 
ExtractFromTitleToBioSourceSubSource 
(CharPtr      title,
 BioSourcePtr biop, 
 CharPtr      mod_name,
 Int4         subtype)
{
  CharPtr valstr;
  CharPtr next_org_loc;
  
  next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  while ((valstr = FindValueFromPairInDeflineBeforeCharPtr (mod_name, title, next_org_loc)) != NULL)
  {
    if (IsNonTextModifier (mod_name)) {
      if (StringICmp (valstr, "FALSE") == 0) {
        valstr = MemFree (valstr);
      } else if (StringICmp (valstr, "TRUE") == 0) {
        biop = AddSubSourceValue (biop, subtype, StringSave (""));
        valstr = MemFree (valstr);
      } else {
        biop = AddSubSourceValue (biop, subtype, valstr);
      }
    } else {
      biop = AddSubSourceValue (biop, subtype, valstr);
    }
    RemoveValueFromDefline (mod_name, title);
    next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  }  
  return biop;
}


static BioSourcePtr ExtractFromTitleToBioSourcePrimers (CharPtr title, BioSourcePtr biop)
{
  CharPtr valstr;
  CharPtr next_org_loc;
  ValNode qual;
  CharPtr mod_name;
  Int4    qual_defs[] = { Source_qual_fwd_primer_name, Source_qual_fwd_primer_seq, Source_qual_rev_primer_name, Source_qual_rev_primer_seq};
  Int4    num_quals = 4;
  Int4    i;
  
  if (StringHasNoText (title)) {
    return biop;
  }
  if (biop == NULL) {
    biop = BioSourceNew ();
  }

  next_org_loc = FindValuePairInDefLine ("organism", title, NULL);

  for (i = 0; i < num_quals; i++) {
    MemSet (&qual, 0, sizeof (ValNode));
    qual.choice = SourceQualChoice_textqual;
    qual.data.intvalue = qual_defs[i];
    mod_name = GetSourceQualName (qual_defs[i]);
    while ((valstr = FindValueFromPairInDeflineBeforeCharPtr (mod_name, title, next_org_loc)) != NULL) {
      SetSourceQualInBioSource (biop, &qual, NULL, valstr, ExistingTextOption_add_qual);
      RemoveValueFromDefline (mod_name, title);
    }      
  }

  return biop;
}


/* this function collects all of the common names prior to the next organism name
 * and assembles a semicolon-delimited list.
 */
extern BioSourcePtr 
ExtractFromTitleToBioSourceCommonName 
(CharPtr      title,
 BioSourcePtr biop)
{
  CharPtr valstr, new_val;
  Int4    new_len;
  CharPtr next_org_loc;
  
  next_org_loc = FindValuePairInDefLine ("organism", title, NULL);  
  while ((valstr = FindValueFromPairInDeflineBeforeCharPtr ("common name", title, next_org_loc)) != NULL)
  {
    if (!StringHasNoText (valstr))
    {
      biop = AddOrgRef (biop);
      if (StringHasNoText (biop->org->common))
      {
        biop->org->common = MemFree (biop->org->common);
        biop->org->common = valstr;
        valstr = NULL;
      }
      else
      {
        new_len = StringLen (biop->org->common) + StringLen (valstr) + 3;
        new_val = (CharPtr) MemNew (new_len * sizeof (Char));
        if (new_val != NULL)
        {
          sprintf (new_val, "%s; %s", biop->org->common, valstr);
          biop->org->common = MemFree (biop->org->common);
          biop->org->common = new_val;
        }
      }
    }
    valstr = MemFree (valstr);
    RemoveValueFromDefline ("common name", title);
    next_org_loc = FindValuePairInDefLine ("organism", title, NULL);
  }  
  return biop;
}

/* When the user specifies multiple organisms on the definition line, modifiers after the
 * second organism go with the second organism, after the third organism go with the third
 * organism, etc.
 */
extern BioSourcePtr ExtractFromDeflineToBioSource (CharPtr defline, BioSourcePtr biop)
{
  CharPtr      taxname = NULL;
  OrgInfoPtr   oip = NULL;
  CharPtr      valstr;
  Nlm_EnumFieldAssocPtr ap;
  Nlm_QualNameAssocPtr  qp;
  CharPtr            next_org_loc;
  
  if (StringHasNoText (defline))
  {
    return NULL;
  }
  
  taxname = FindValueFromPairInDefline ("organism", defline);
  RemoveValueFromDefline ("organism", defline);
  if (StringHasNoText (taxname))
  {
    taxname = MemFree (taxname);
    return NULL;
  }
  else
  {
    biop = AddOrgRef (biop);
    if (biop == NULL)
    {
      return biop;
    }
    LoadOrganismList ();
    oip = FindByTaxName (taxname);
    SetTaxNameAndRemoveTaxRef (biop->org, taxname);
  }
  
  /* add division */
  if (oip != NULL && !StringHasNoText (oip->div))
  {
    biop = AddOrgName (biop);
    if (biop == NULL)
    {
      return biop;
    }    
    biop->org->orgname->div = StringSave (oip->div);
  }
  
  /* add common name (s) - if there are multiple entries, separate with semicolon */
  biop = ExtractFromTitleToBioSourceCommonName (defline, biop);
  /* if common name was not supplied in defline, use common name from organism list */
  if (biop->org == NULL || StringHasNoText (biop->org->common))
  {
    if (oip != NULL && !StringHasNoText (oip->common))
    {
      biop = AddOrgRef (biop);
      if (biop == NULL)
      {
        return biop;
      }
      biop->org->common = StringSave (oip->common);
    }
  }
  
  /* add lineage */
  if (oip != NULL && !StringHasNoText (oip->lineage))
  {
    biop = AddOrgName (biop);
    if (biop == NULL)
    {
      return biop;
    }
    biop->org->orgname->lineage = StringSave (oip->lineage);
  }
  
  /* add origin */
  next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);
  valstr = FindValueFromPairInDeflineBeforeCharPtr ("origin", defline, next_org_loc);
  if (!StringHasNoText (valstr))
  {
    for (ap = biosource_origin_alist; ap->name != NULL; ap++) {
      if (StringICmp (valstr, ap->name) == 0) {
        if (biop == NULL)
        {
          biop = BioSourceNew ();
        }
        if (biop == NULL)
        {
          return biop;
        }
        biop->origin = (Uint1) ap->value;
      }
    }
  }
  if (valstr != NULL)
  {
    RemoveValueFromDefline ("origin", defline);
    next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);
  }
  valstr = MemFree (valstr);
  
  valstr = FindValueFromPairInDeflineBeforeCharPtr ("lineage", defline, next_org_loc);
  if (!StringHasNoText (valstr))
  {
    biop = AddOrgName (biop);
  }
  if (!StringHasNoText (valstr) && StringCmp (valstr, biop->org->orgname->lineage) != 0)
  {
    biop = AddOrgModValue (biop, ORGMOD_old_lineage, valstr);
    valstr = NULL;
  }
  if (valstr != NULL)
  {
    RemoveValueFromDefline ("lineage", defline);
  }
  valstr = MemFree (valstr);
  
  biop = SetAllGeneticCodesFromTitle (biop, defline);
  next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);
        
  for (qp = current_orgmod_subtype_alist; qp->name != NULL; qp++) {
    biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, qp->name, qp->value);
  }
  for (qp = current_subsource_subtype_alist; qp->name != NULL; qp++) {
    biop = ExtractFromTitleToBioSourceSubSource (defline, biop, qp->name, qp->value);
  }

  /* parse primers */
  biop = ExtractFromTitleToBioSourcePrimers (defline, biop);
  
  /* parse notes */
  biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, "note-orgmod", 255);
  biop = ExtractFromTitleToBioSourceSubSource (defline, biop, "note-subsrc", 255);
  biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, "note", 255);
  biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, "comment", 255);
  biop = ExtractFromTitleToBioSourceSubSource (defline, biop, "subsource", 255);


  next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);

  /* set location */
  valstr = FindValueFromPairInDeflineBeforeCharPtr ("location", defline, next_org_loc);
  if (StringHasNoText (valstr))
  {
    if (biop == NULL)
    {
      biop = BioSourceNew ();
    }
    if (biop == NULL)
    {
      return biop;
    }
    biop->genome = 1;
  }
  else if (StringICmp (valstr, "Mitochondrial") == 0)
  {
    if (biop == NULL)
    {
      biop = BioSourceNew ();
    }
    if (biop == NULL)
    {
      return biop;
    }
    biop->genome = 5;
  }
  else
  {
    for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
      if (StringICmp (valstr, ap->name) == 0) {
        if (biop == NULL)
        {
          biop = BioSourceNew ();
        }
        if (biop == NULL)
        {
          return biop;
        }
        biop->genome = (Uint1) ap->value;
      }
    }
  }
  if (valstr != NULL)
  {
    RemoveValueFromDefline ("location", defline);
  }
  valstr = MemFree (valstr);
 
  TrimSpacesAroundString (defline);
  
  return biop;
  
}

extern Boolean ProcessOneNucleotideTitle (Int2 seqPackage,
                                          SeqEntryPtr nsep, SeqEntryPtr top);
                                          

static void ParseDeflineToBiop(CharPtr defline, BioSourcePtr biop) 
{
  CharPtr      taxname = NULL;
  OrgInfoPtr   oip = NULL;
  CharPtr      valstr;
  EnumFieldAssocPtr  ap;
  Nlm_QualNameAssocPtr qp;
  CharPtr            next_org_loc;
  
  if (StringHasNoText (defline) || biop == NULL)
  {
    return;
  }
  
  taxname = FindValueFromPairInDefline ("organism", defline);
  RemoveValueFromDefline ("organism", defline);
  if (StringHasNoText (taxname))
  {
    taxname = MemFree (taxname);
  }
  else
  {
    biop = AddOrgRef (biop);
    biop->org->taxname = taxname;
    LoadOrganismList ();
    oip = FindByTaxName (taxname);
  }
  
  /* add division */
  if (oip != NULL && !StringHasNoText (oip->div))
  {
    biop = AddOrgName (biop);
    biop->org->orgname->div = StringSave (oip->div);
  }
  
  /* add common name (s) - if there are multiple entries, separate with semicolon */
  biop = ExtractFromTitleToBioSourceCommonName (defline, biop);
  /* if common name was not supplied in defline, use common name from organism list */
  if (biop->org == NULL || StringHasNoText (biop->org->common))
  {
    if (oip != NULL && !StringHasNoText (oip->common))
    {
      biop = AddOrgRef (biop);
      biop->org->common = StringSave (oip->common);
    }
  }
  
  /* add lineage */
  if (oip != NULL && !StringHasNoText (oip->lineage))
  {
    biop = AddOrgName (biop);
    biop->org->orgname->lineage = StringSave (oip->lineage);
  }
  
  /* add origin */
  next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);
  valstr = FindValueFromPairInDeflineBeforeCharPtr ("origin", defline, next_org_loc);
  if (!StringHasNoText (valstr))
  {
    for (ap = biosource_origin_alist; ap->name != NULL; ap++) {
      if (StringICmp (valstr, ap->name) == 0) {
        biop->origin = (Uint1) ap->value;
      }
    }
  }
  if (valstr != NULL)
  {
    RemoveValueFromDefline ("origin", defline);
    next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);
  }
  valstr = MemFree (valstr);
  
  valstr = FindValueFromPairInDeflineBeforeCharPtr ("lineage", defline, next_org_loc);
  if (!StringHasNoText (valstr))
  {
    biop = AddOrgName (biop);
  }
  if (!StringHasNoText (valstr) && StringCmp (valstr, biop->org->orgname->lineage) != 0)
  {
    biop = AddOrgModValue (biop, ORGMOD_old_lineage, valstr);
    valstr = NULL;
  }
  if (valstr != NULL)
  {
    RemoveValueFromDefline ("lineage", defline);
  }
  valstr = MemFree (valstr);
  
  biop = SetAllGeneticCodesFromTitle (biop, defline);
  next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);
        
  for (qp = current_orgmod_subtype_alist; qp->name != NULL; qp++) {
    biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, qp->name, qp->value);
  }
  for (qp = current_subsource_subtype_alist; qp->name != NULL; qp++) {
    biop = ExtractFromTitleToBioSourceSubSource (defline, biop, qp->name, qp->value);
  }
  
  biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, "note-orgmod", 255);
  biop = ExtractFromTitleToBioSourceSubSource (defline, biop, "note-subsrc", 255);
  biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, "note", 255);
  biop = ExtractFromTitleToBioSourceOrgMod (defline, biop, "comment", 255);
  biop = ExtractFromTitleToBioSourceSubSource (defline, biop, "subsource", 255);


  next_org_loc = FindValuePairInDefLine ("organism", defline, NULL);

  /* set location */
  valstr = FindValueFromPairInDeflineBeforeCharPtr ("location", defline, next_org_loc);
  if (StringHasNoText (valstr))
  {
    /* don't set defaults */
  }
  else if (StringICmp (valstr, "Mitochondrial") == 0)
  {
    biop->genome = 5;
  }
  else
  {
    for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
      if (StringICmp (valstr, ap->name) == 0) {
        biop->genome = (Uint1) ap->value;
      }
    }
  }
  if (valstr != NULL)
  {
    RemoveValueFromDefline ("location", defline);
  }
  valstr = MemFree (valstr);
 
  TrimSpacesAroundString (defline);
}


static void ParseModifiersFromDeflineCallback (BioseqPtr bsp, Pointer userdata)
{
  CharPtr           title;
  SeqDescrPtr       sdp, sdp_biop, prev_sdp = NULL;
  BioSourcePtr      biop = NULL;
  CharPtr           valstr;
  SeqMgrDescContext context;
  
  if (bsp == NULL) return;
  
  if (ISA_aa(bsp->mol)) {
    return;
  }
  
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title) {
    prev_sdp = sdp;
    sdp = sdp->next;
  }
  if (sdp == NULL || sdp->data.ptrvalue == NULL) {
    return;
  }
  
  title = sdp->data.ptrvalue;
  
  if (StringChr(title, '[') == NULL || StringChr(title, ']') == NULL) {
    return;
  }
  
  /* parse moltype values */
  SetMoleculeAndMolTypeFromTitle (bsp, title, SEQ_PKG_GENBANK);
  
  /* get topology from defline */
  valstr = FindValueFromPairInDefline ("topology", title);
  if (valstr != NULL)
  {
    if (!StringHasNoText (valstr))
    {
      bsp->topology = TopologyFromString (valstr);
    }
    RemoveValueFromDefline ("topology", title);
    valstr = MemFree (valstr);
  }
  
  /* add bankit comment for genetic code */
  valstr = FindValueFromPairInDefline ("gencode_comment", title);
  if (valstr != NULL)
  {
    AddGeneticCodeComment (bsp, valstr);
    RemoveValueFromDefline ("gencode_comment", title);
    valstr = MemFree (valstr);
  }
  
  sdp_biop = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);

  if(sdp_biop == NULL) {
    if (bsp->idx.parenttype == OBJ_BIOSEQSET 
        && bsp->idx.parentptr != NULL 
        && ((BioseqSetPtr)bsp->idx.parentptr)->_class == BioseqseqSet_class_parts) {
        /* don't put sources on parts */
    } else {
        biop = BioSourceNew();
        sdp_biop = SeqDescrNew(bsp->descr);
        sdp_biop->choice = Seq_descr_source;
        sdp_biop->data.ptrvalue = biop;
    }
  } else {
    biop = sdp_biop->data.ptrvalue;
  }          
  
  ParseDeflineToBiop (title, biop);
  
  if (StringHasNoText (title)) {
    /* remove empty defline */
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    if (prev_sdp == NULL) {
      bsp->descr = sdp->next;
    } else {
      prev_sdp->next = sdp->next;
    }
    sdp->next = NULL;
    sdp = SeqDescrFree (sdp);
  }
}

                                         
extern void ParseModifiersFromDefline (IteM i)
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
  
  VisitBioseqsInSep (sep, NULL, ParseModifiersFromDeflineCallback);
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}
                                    
                                          
extern Boolean ProcessOneNucleotideTitle (Int2 seqPackage, 
                                          SeqEntryPtr nsep, SeqEntryPtr top)

{
  BioSourcePtr       biop = NULL;
  BioseqSetPtr       bssp;
  BioseqPtr          nbsp;
  Boolean            needbiop;
  SeqEntryPtr        sep;
  CharPtr            str;
  CharPtr            valstr;
  CharPtr            title;
  ValNodePtr         vnp;
  Int4               topology;
#if 0
  SeqFeatPtr         sfp;
#endif  

  if (nsep == NULL || top == NULL) return FALSE;
  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  if (nbsp == NULL) return FALSE;
  if (! ISA_na (nbsp->mol)) return FALSE;
  str = MemNew (PROC_NUC_STR_SIZE * sizeof (Char));
  if (str == NULL) return FALSE;
  sep = NULL;
 
  SeqEntryExplore (top, (Pointer) &sep, FindFirstSeqEntryTitle);
  sep = FindNucSeqEntry (sep);
  if (sep != NULL) {
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
    if (vnp != NULL && vnp->data.ptrvalue != NULL) {
      title = (CharPtr) vnp->data.ptrvalue;
      
      SetMoleculeAndMolTypeFromTitle (nbsp, title, seqPackage);

      if (nbsp->topology == 0)
      {
        topology = TOPOLOGY_LINEAR;
      }
      else
      {
        topology = nbsp->topology;
      }
  
      /* get topology from defline */
      valstr = FindValueFromPairInDefline ("topology", title);
      if (valstr != NULL)
      {
        if (!StringHasNoText (valstr))
        {
          topology = TopologyFromString (valstr);
        }
        RemoveValueFromDefline ("topology", title);
        valstr = MemFree (valstr);
      }
      nbsp->topology = topology;
      
      /* add bankit comment for genetic code */
      valstr = FindValueFromPairInDefline ("gencode_comment", title);
      if (valstr != NULL)
      {
        AddGeneticCodeComment (nbsp, valstr);
        RemoveValueFromDefline ("gencode_comment", title);
        valstr = MemFree (valstr);
      }

      needbiop = FALSE;
      
      if (PackageTypeIsSet (seqPackage))
      {
        needbiop = TRUE;
        if (GetAppParam ("SEQUIN", "PREFERENCES", "BIOSRCONALL", NULL, str, PROC_NUC_STR_SIZE)) {
          if (StringICmp (str, "FALSE") == 0) {
            needbiop = FALSE;
          }
        }
      }
      
      vnp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
      if (vnp == NULL)
      {
        biop = ExtractFromDeflineToBioSource (title, NULL);
        if (biop == NULL && needbiop)
        {
          biop = BioSourceNew ();
        }

        if (biop != NULL)
        {
          vnp = CreateNewDescriptor (top, Seq_descr_source);
          if (vnp != NULL) {
            vnp->data.ptrvalue = (Pointer) biop;
          }
        }
#if 0        
        biop = BioSourceFromDefline (title);
        while (biop != NULL)
        {
          sfp = CreateNewFeature (sep, NULL, SEQFEAT_BIOSRC, NULL);
          if (sfp != NULL)
          {
            sfp->data.value.ptrvalue = biop;
          }
          biop = BioSourceFromDefline (title);
        }
#endif        
      }

      if (StringHasNoText (title) || sep != top) {
        vnp = NULL;
        if (IS_Bioseq (sep)) {
          nbsp = (BioseqPtr) sep->data.ptrvalue;
          vnp = ValNodeExtract (&(nbsp->descr), Seq_descr_title);
        } else if (IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          vnp = ValNodeExtract (&(bssp->descr), Seq_descr_title);
        }
        if (vnp != NULL && StringHasNoText ((CharPtr) vnp->data.ptrvalue)) {
          vnp = ValNodeFreeData (vnp);
        }
        if (sep != top && vnp != NULL) {
          if (IS_Bioseq (top)) {
            nbsp = (BioseqPtr) top->data.ptrvalue;
            ValNodeLink (&(nbsp->descr), vnp);
          } else if (IS_Bioseq_set (top)) {
            bssp = (BioseqSetPtr) top->data.ptrvalue;
            ValNodeLink (&(bssp->descr), vnp);
          }
        }
      }
    }
  } else {
    needbiop = FALSE;
    if (PackageTypeIsSet (seqPackage)
        || seqPackage == SEQ_PKG_GENOMICCDNA)
    {
      needbiop = TRUE;
      if (GetAppParam ("SEQUIN", "PREFERENCES", "BIOSRCONALL", NULL, str, PROC_NUC_STR_SIZE)) {
        if (StringICmp (str, "FALSE") == 0) {
          needbiop = FALSE;
        }
      }
    }
  }
  MemFree (str);
  
  return TRUE;
}

static Boolean AutomaticNucleotideProcess (SequencesFormPtr sqfp, SeqEntryPtr nsep,
                                           SeqEntryPtr top)

{
  BioseqSetPtr  bssp;
  Boolean       rsult;
  SeqEntryPtr   tmp;

  if (sqfp == NULL || nsep == NULL || top == NULL) return FALSE;
  if (IS_Bioseq_set (nsep)) {
    bssp = (BioseqSetPtr) nsep->data.ptrvalue;
    rsult = FALSE;
    if (bssp != NULL) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        if (AutomaticNucleotideProcess (sqfp, tmp, top)) {
          rsult = TRUE;
        }
      }
    }
    return rsult;
  }
  return ProcessOneNucleotideTitle (sqfp->seqPackage, 
                                    nsep, top);
}

typedef struct idlist {
  BioseqPtr  bsp;
  CharPtr    key;
  struct idlist PNTR left;
  struct idlist PNTR right;
} IdList, PNTR IdListPtr;

static void BuildTree (IdListPtr PNTR head, BioseqPtr bsp, CharPtr x)

{
  Int2       comp;
  IdListPtr  idlist;
  SeqIdPtr   sip;
  CharPtr    str;

  if (*head != NULL) {
    idlist = *head;
    comp = StringICmp (idlist->key, x);
    if (comp < 0) {
      BuildTree (&(idlist->right), bsp, x);
    } else if (comp > 0) {
      BuildTree (&(idlist->left), bsp, x);
    } else {
      sip = MakeNewProteinSeqId (NULL, NULL);
      if (sip != NULL) {
        bsp->id = SeqIdFree (bsp->id);
        bsp->id = sip;
        SeqMgrReplaceInBioseqIndex (bsp);
        str = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);
        BuildTree (head, bsp, str);
        str = MemFree (str);
      }
    }
  } else {
    idlist = MemNew (sizeof (IdList));
    if (idlist != NULL) {
      *head = idlist;
      idlist->bsp = bsp;
      idlist->key = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);
      idlist->left = NULL;
      idlist->right = NULL;
    }
  }
}

static void FreeTree (IdListPtr PNTR head)

{
  IdListPtr  idlist;

  if (head != NULL && *head != NULL) {
    idlist = *head;
    FreeTree (&(idlist->left));
    FreeTree (&(idlist->right));
    MemFree (idlist->key);
    MemFree (idlist);
  }
}

static void ResolveCollidingIDs (IdListPtr PNTR head, SeqEntryPtr list)

{
  BioseqPtr  bsp;
  CharPtr    str;

  if (head == NULL) return;
  while (list != NULL) {
    if (IS_Bioseq (list)) {
      bsp = (BioseqPtr) list->data.ptrvalue;
      if (bsp != NULL) {
        str = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);
        BuildTree (head, bsp, str);
        str = MemFree (str);
      }
    }
    list = list->next;
  }
}


static void PutMolInfoOnSeqEntry (SequencesFormPtr sqfp, SeqEntryPtr sep)

{
  BioseqSetPtr bssp;
  MolInfoPtr   mip;
  ValNodePtr   vnp;

  if (sqfp != NULL && sep != NULL) {
    if (IS_Bioseq_set (sep))
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
      {
      	PutMolInfoOnSeqEntry (sqfp, sep);
      }
      return;
    }

    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL)
    {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    }
    if (vnp != NULL)
    {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip == NULL)
      {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = mip;
      }
    }
  }
}

static void PrefixOrgToDefline (SeqEntryPtr sep)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  CharPtr       def;
  OrgRefPtr     orp;
  CharPtr       ptr;
  CharPtr       str;
  Char          taxname [64];
  ValNodePtr    ttl;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        PrefixOrgToDefline (sep);
      }
      return;
    }
  }

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  taxname [0] = '\0';
  orp = NULL;
  biop = NULL;
  ttl = NULL;
  vnp = bsp->descr;
  for (vnp = bsp->descr; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) vnp->data.ptrvalue;
    } else if (vnp->choice == Seq_descr_org) {
      orp = (OrgRefPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == Seq_descr_title) {
      ttl = vnp;
    }
  }
  if (orp == NULL && biop != NULL) {
    orp = biop->org;
  }
  if (orp == NULL) return;
  if (ttl == NULL) return;
  StringNCpy_0 (taxname, orp->taxname, sizeof (taxname));
  ptr = StringSearch (taxname, "(");
  if (ptr != NULL) {
    *ptr = '\0';
  }
  TrimSpacesAroundString (taxname);
  if ((StringICmp (taxname, "Human immunodeficiency virus type 1") == 0) ||
      (StringICmp (taxname, "Human immunodeficiency virus 1") == 0)) {
    StringCpy (taxname, "HIV-1");
  } else if ((StringICmp (taxname,"Human immunodeficiency virus type 2")==0) ||
	     (StringICmp (taxname,"Human immunodeficiency virus 2")==0)) {
    StringCpy (taxname, "HIV-2");
  }

  def = (CharPtr) ttl->data.ptrvalue;
  if (StringHasNoText (def)) return;

  ptr = StringISearch (def, taxname);
  if (ptr != NULL && ptr == def) return;
  str = MemNew ((StringLen (taxname) + StringLen (def) + 4) * sizeof (Char));
  if (str == NULL) return;
  StringCpy (str, taxname);
  StringCat (str, " ");
  StringCat (str, def);
  ttl->data.ptrvalue = MemFree (ttl->data.ptrvalue);
  ttl->data.ptrvalue = str;
}

static CharPtr onecomponent = "\
Multiple sequence components are expected in this submission.\n\
They should all be read in at the same time from the same file.";

static void OnlyOneComponentWarning (SequencesFormPtr sqfp)

{
  CharPtr  type;

  if (sqfp != NULL) {
    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA
        || PackageTypeIsSingle (sqfp->seqPackage))
    {
      return;
    }
    switch (sqfp->seqPackage) {
      case SEQ_PKG_POPULATION :
        type = "population set";
        break;
      case SEQ_PKG_PHYLOGENETIC :
        type = "phylogenetic set";
        break;
      case SEQ_PKG_MUTATION :
        type = "mutation set";
        break;
      case SEQ_PKG_ENVIRONMENT :
        type = "environmental samples";
        break;
      case SEQ_PKG_GENBANK :
        type = "batch submission";
        break;
      case SEQ_PKG_TSA:
        type = "Transcriptome Shotgun Assembly";
        break;
      default :
        type = "unknown set";
        break;
    }
    Message (MSG_OK, "WARNING - There is only one component in this %s.\n%s",
             type, onecomponent);
  }
}

/*---------------------------------*/
/* Parse the gene and gene-related */
/* fields from the title.          */
/*---------------------------------*/
extern void 
AddGeneFeatureFromTitle 
(SeqEntryPtr nucsep,
 CharPtr ttl, 
 SeqLocPtr slp)
{
  CharPtr    gene = NULL;
  CharPtr    gene_desc = NULL;
  CharPtr    allele = NULL;
  CharPtr    gene_syn = NULL;
  GeneRefPtr grp = NULL;
  SeqFeatPtr sfp;
  SeqIdPtr   sip;
  BioseqPtr  nbsp, bsp;
  SeqLocPtr  gslp;
  Boolean    hasNulls;
  
  if (nucsep == NULL || !IS_Bioseq (nucsep) 
      || (nbsp = (BioseqPtr) nucsep->data.ptrvalue) == NULL
      || StringHasNoText (ttl) || slp == NULL)
  {
    return;
  }
  
  gene = FindValueFromPairInDefline ("gene", ttl);
  if (!StringHasNoText (gene))
  {
    gene_desc = StringChr (gene, ';');
    if (gene_desc != NULL) {
      *gene_desc = '\0';
      gene_desc++;
      allele = StringChr (gene_desc, ';');
      if (allele != NULL) {
        *allele = '\0';
        allele++;
      }
    }
    grp = CreateNewGeneRef (gene, allele, gene_desc, FALSE);
  }
  gene = MemFree (gene);
  
  /*-----------------------------------------*/
  /* Parse the gene_syn field from the title */
  /*-----------------------------------------*/
  
  gene_syn = FindValueFromPairInDefline ("gene_syn", ttl);
  if (!StringHasNoText (gene_syn))
  {
    if (grp == NULL) {
      grp = GeneRefNew ();
    }
    ValNodeCopyStr(&(grp->syn),0,gene_syn);
  }
  gene_syn = MemFree (gene_syn);

  /* Create the gene feature */
  if (grp != NULL) {
    if (ExtendGene (grp, nucsep, slp)) {
      grp = GeneRefFree (grp);
    } else {
      sfp = CreateNewFeature (nucsep, NULL, SEQFEAT_GENE, NULL);
      if (sfp != NULL) {
        sfp->data.value.ptrvalue = (Pointer) grp;
        sfp->location = SeqLocFree (sfp->location);
        sfp->location = AsnIoMemCopy ((Pointer) slp,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);
        sip = SeqLocId (sfp->location);
        if (sip != NULL) {
          bsp = BioseqFind (sip);
        } else {
          bsp = nbsp;
        }
        if (bsp != NULL) {
          gslp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
          if (gslp != NULL) {
            sfp->location = SeqLocFree (sfp->location);
            sfp->location = gslp;
            if (bsp->repr == Seq_repr_seg) {
              gslp = SegLocToPartsEx (bsp, sfp->location, TRUE);
              sfp->location = SeqLocFree (sfp->location);
              sfp->location = gslp;
              hasNulls = LocationHasNullsBetween (sfp->location);
              sfp->partial = (sfp->partial || hasNulls);
            }
            FreeAllFuzz (gslp);
          }
        }
      }
    }
    RemoveValueFromDefline ("gene", ttl);
    RemoveValueFromDefline ("gene_syn", ttl);
  }
}

extern SeqFeatPtr AddProteinFeatureFromDefline (SeqEntryPtr psep, CharPtr title)
{
  CharPtr    activity = NULL;
  CharPtr    ec = NULL;
  CharPtr    prot_name = NULL;
  CharPtr    prot_desc = NULL;
  CharPtr    other_prot_desc = NULL, tmp_desc;
  ProtRefPtr prp;
  SeqFeatPtr sfp = NULL;
  
  if (psep == NULL)
  {
    return NULL;
  }
  
	/*-----------------------------------------*/
	/* Parse the function field from the title */
	/*-----------------------------------------*/

  activity = FindValueFromPairInDefline ("function", title);

	/*------------------------------------------*/
	/* Parse the EC_number field from the title */
	/*------------------------------------------*/

  ec = FindValueFromPairInDefline ("EC_number", title);

	/*---------------------------------*/
	/* Parse the protein and prot_desc */
	/* fields from the title.          */
	/*---------------------------------*/

  prot_name = FindValueFromPairInDefline ("protein", title);

	/*---------------------------------*/
	/* If we found a protein value ... */
	/*---------------------------------*/
  if (!StringHasNoText (prot_name))
  {
	  /*----------------------------------------------*/
	  /* ... search for a protein description, either */
	  /*     in the prot field (seperated by a ';')   */
	  /*     or in its own 'prot_desc' field.         */
	  /*----------------------------------------------*/

    prot_desc = StringChr (prot_name, ';');
	  if (prot_desc != NULL)
	  {
		  *prot_desc = '\0';
		  prot_desc++;
		  /* ignore this description if empty */
		  if (StringHasNoText (prot_desc))
		  {
		    prot_desc = NULL;
		  }
		  else
		  {
		    prot_desc = StringSave (prot_desc);
		  }
	  }
  }
	other_prot_desc = FindValueFromPairInDefline ("prot_desc", title);
	if (StringHasNoText (other_prot_desc))
	{
	  other_prot_desc = MemFree (other_prot_desc);
	}
	else
	{
   if (prot_desc == NULL)
	  {
	    prot_desc = other_prot_desc;
	    other_prot_desc = NULL;
	  }
	  else 
	  {
      tmp_desc = (CharPtr) MemNew ((StringLen (prot_desc) + StringLen (other_prot_desc) + 3)
	                               * sizeof (Char));
	    if (tmp_desc != NULL)
	    {
	      StringCpy (tmp_desc, prot_desc);
	      StringCat (tmp_desc, ";");
  	    StringCat (tmp_desc, other_prot_desc);
	      prot_desc = MemFree (prot_desc);
	      other_prot_desc = MemFree (other_prot_desc);
	      prot_desc = tmp_desc;
	    }
	  }
	}
	
	/*--------------------------------*/
	/* ... add the prot and prot_desc */
	/*     to the Seq Features.       */
	/*--------------------------------*/

	prp = CreateNewProtRef (prot_name, prot_desc, ec, activity);
	if (prp != NULL)
	{
		sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
		if (sfp != NULL)
		{
		  sfp->data.value.ptrvalue = (Pointer) prp;
		  RemoveValueFromDefline ("protein", title);
		  RemoveValueFromDefline ("prot_desc", title);
		  RemoveValueFromDefline ("function", title);
		  RemoveValueFromDefline ("EC_number", title);
		}
	}
  return sfp;
}

extern void 
AddCodingRegionFieldsFromProteinTitle 
(CdRegionPtr  crp,
 CharPtr      title, 
 CharPtr PNTR pcomment)
{
  CharPtr comment, comment_loc, total_comment = NULL, tmp_comment;
  
  if (crp == NULL || StringHasNoText (title))
  {
    return;
  }
  
	/*---------------------*/
	/* Parse the ORF field */
	/*---------------------*/
  if (FindValuePairInDefLine ("orf", title, NULL) != NULL)
  {
    crp->orf = TRUE;
    RemoveValueFromDefline ("orf", title);
  }

  if (pcomment == NULL)
  {
    return;
  }

	/*-------------------------------*/
	/* Parse the comment/note fields */
	/*-------------------------------*/
  comment_loc = FindValuePairInDefLine ("comment", title, NULL);
  while (comment_loc != NULL)
  {
    comment = FindValueFromPairInDefline ("comment", comment_loc);
    if (!StringHasNoText (comment))
    {
      if (total_comment == NULL)
      {
        total_comment = comment;
        comment = NULL;
      }
      else
      {
        tmp_comment = (CharPtr) MemNew ((StringLen (total_comment) + StringLen (comment) + 3) * sizeof (Char));
        if (tmp_comment != NULL)
        {
          StringCpy (tmp_comment, total_comment);
          StringCat (tmp_comment, ";");
          StringCat (tmp_comment, comment);
          total_comment = MemFree (total_comment);
          total_comment = tmp_comment;
        }
      }
    }
    comment = MemFree (comment);
    RemoveValueFromDefline ("comment", title);
    comment_loc = FindValuePairInDefLine ("comment", title, NULL);
  }
  
  *pcomment = total_comment;
}

static void AutomaticMrnaProcess (SeqEntryPtr nucsep, SeqEntryPtr mrnasep, Boolean partial5, Boolean partial3)

{
  CharPtr     mrna = NULL;
  CharPtr     comment = NULL;
  BioseqPtr   bsp;
  MolInfoPtr  mip;
  BioseqPtr   mrnabsp;
  BioseqPtr   nucbsp;
  SeqLocPtr   oldslp;
  RnaRefPtr   rrp;
  SeqFeatPtr  sfp;
  SeqIdPtr    sip;
  SeqLocPtr   slp;
  CharPtr     ttl;
  ValNodePtr  vnp;

  if (nucsep == NULL || mrnasep == NULL) return;
  if (IS_Bioseq (nucsep) && IS_Bioseq (mrnasep)) {
    nucbsp = (BioseqPtr) nucsep->data.ptrvalue;
    mrnabsp = (BioseqPtr) mrnasep->data.ptrvalue;
    if (nucbsp == NULL || mrnabsp == NULL) return;
    slp = AlignmRNA2genomic (nucbsp, mrnabsp);
    if (slp == NULL) return;
    sip = SeqLocId (slp);
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        if (bsp->repr == Seq_repr_seg) {
          oldslp = slp;
          slp = SegLocToParts (bsp, oldslp);
          FreeAllFuzz (slp);
          SeqLocFree (oldslp);
        }
      }
    }
    StripLocusFromSeqLoc (slp);
    ttl = NULL;
    vnp = ValNodeFindNext (mrnabsp->descr, NULL, Seq_descr_title);
    if (vnp != NULL) {
      ttl = (CharPtr) vnp->data.ptrvalue;
    }
    if (ttl != NULL) {
      AddGeneFeatureFromTitle (nucsep, ttl, slp);
    
      /* get mRNA name */
      mrna = FindValueFromPairInDefline ("mrna", ttl);
      RemoveValueFromDefline ("mrna", ttl);
      if (StringHasNoText (mrna))
      {
        mrna = MemFree (mrna);
        mrna = FindValueFromPairInDefline ("cdna", ttl);
        RemoveValueFromDefline ("cdna", ttl);
      }
    }
    rrp = RnaRefNew ();
    if (rrp != NULL) {
      rrp->type = 2;
      if (! StringHasNoText (mrna)) {
        rrp->ext.choice = 1;
        rrp->ext.value.ptrvalue = mrna;
        mrna = NULL;
      }
      sfp = CreateNewFeature (nucsep, NULL, SEQFEAT_RNA, NULL);
      if (sfp != NULL) {
        sfp->data.value.ptrvalue = (Pointer) rrp;
        sfp->location = SeqLocFree (sfp->location);
        sfp->location = AsnIoMemCopy ((Pointer) slp,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);
        SetSeqFeatProduct (sfp, mrnabsp);
        SetSeqLocPartial (sfp->location, partial5, partial3);
        sfp->partial = (sfp->partial || partial5 || partial3);
        if (ttl != NULL) {
          comment = FindValueFromPairInDefline ("comment", ttl);
          if (!StringHasNoText (comment)) {
            sfp->comment = comment;
          }
          else
          {
            comment = MemFree (comment);
          }
          RemoveValueFromDefline ("comment", ttl);
        }
      }
    }
    mrna = MemFree (mrna);
    SeqLocFree (slp);
    if (StringHasNoText (ttl)) {
      ValNodeExtract (&(mrnabsp->descr), Seq_descr_title);
    }
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = 3;
      if (partial5 && partial3) {
        mip->completeness = 5;
      } else if (partial5) {
        mip->completeness = 3;
      } else if (partial3) {
        mip->completeness = 4;
      }
      vnp = CreateNewDescriptor (mrnasep, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
      }
    }
    mrnabsp->mol = Seq_mol_rna;
  }
}

static CharPtr LookForValueInBioseq (SeqEntryPtr sep, Uint1 mol, CharPtr valname)
{
  BioseqPtr   bsp;
  CharPtr     title;
  ValNodePtr  vnp;

  if (sep == NULL || StringHasNoText (valname)) return FALSE;
  if (! IS_Bioseq (sep)) return FALSE;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || bsp->mol != mol || bsp->descr == NULL) return FALSE;
  vnp = ValNodeFindNext (bsp->descr, NULL, Seq_descr_title);
  if (vnp == NULL || vnp->data.ptrvalue == NULL) return FALSE;
  title = (CharPtr) vnp->data.ptrvalue;
  return FindValueFromPairInDefline (valname, title);
}

static void FindBioseqWithValue (SeqEntryPtr sep, Uint1 mol, CharPtr valname, CharPtr value, SeqEntryPtr PNTR rsult)
{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  CharPtr       match_value;

  if (sep == NULL || sep->data.ptrvalue == NULL || rsult == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    match_value = LookForValueInBioseq (sep, mol, valname);
    if (StringICmp (match_value, value))
    {
      *rsult = sep;
    }
    match_value = MemFree (match_value);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindBioseqWithValue (sep, mol, valname, value, rsult);
    }
  }
}

static void RemoveValueFromBioseq (SeqEntryPtr sep, CharPtr valname)
{
  BioseqPtr   bsp;
  ValNodePtr  vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || bsp->descr == NULL) return;
  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  RemoveValueFromDefline (valname, vnp->data.ptrvalue);
  if (StringHasNoText (vnp->data.ptrvalue)) {
    ValNodeExtract (&(bsp->descr), Seq_descr_title);
  }  
}

static SeqEntryPtr FindRnaByRefOnRna (SeqEntryPtr sep, SeqEntryPtr psep)

{
  SeqEntryPtr  msep;
  CharPtr      prot_name;

  msep = NULL;
  if (sep == NULL || psep == NULL) return NULL;
  prot_name = LookForValueInBioseq (psep, Seq_mol_aa, "prot");
  if (!StringHasNoText (prot_name))
  {
    FindBioseqWithValue (sep, Seq_mol_rna, "prot", prot_name, &msep);
    RemoveValueFromBioseq (msep, "prot");
  }
  prot_name = MemFree (prot_name);
  return msep;
}

static void FindRnaByName (SeqEntryPtr sep, CharPtr str, SeqEntryPtr PNTR msep)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  RnaRefPtr     rrp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (str == NULL || msep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_RNA) {
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->type == 2 && rrp->ext.choice == 1 && sfp->product != NULL) {
            if (StringICmp (rrp->ext.value.ptrvalue, str) == 0) {
              bsp = BioseqFind (SeqLocId (sfp->product));
              if (bsp != NULL) {
                *msep = SeqMgrGetSeqEntryForData (bsp);
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  if (bssp != NULL) {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindRnaByName (sep, str, msep);
    }
  }
}

static SeqEntryPtr FindRnaByRefOnProtein (SeqEntryPtr sep, SeqEntryPtr psep)

{
  SeqEntryPtr  msep;
  CharPtr      mrna_name;

  msep = NULL;
  if (sep == NULL || psep == NULL) return NULL;
  mrna_name = LookForValueInBioseq (psep, Seq_mol_aa, "mrna");
  if (!StringHasNoText (mrna_name))
  {
    FindRnaByName (sep, mrna_name, &msep);
    RemoveValueFromBioseq (msep, "mrna");
  }
  mrna_name = MemFree (mrna_name);
  return msep;
}

static void FindRnaByLocationOverlap (SeqEntryPtr sep, SeqLocPtr slp,
                                      Int4Ptr mindiff, SeqEntryPtr PNTR msep)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  Int4          diff;
  RnaRefPtr     rrp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (slp == NULL || mindiff == NULL || msep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_RNA) {
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->type == 2 && sfp->product != NULL) {
            diff = SeqLocAinB (slp, sfp->location);
            if (diff >= 0) {
              if (diff < *mindiff) {
                bsp = BioseqFind (SeqLocId (sfp->product));
                if (bsp != NULL) {
                  *mindiff = diff;
                  *msep = SeqMgrGetSeqEntryForData (bsp);
                }
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  if (bssp != NULL) {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindRnaByLocationOverlap (sep, slp, mindiff, msep);
    }
  }
}

static void FuseNucProtBiosources (SeqEntryPtr sep)

{
  BioSourcePtr  biop1, biop2;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    PNTR prev;
  ValNodePtr    sdp1, sdp2;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_nuc_prot) return;
  tmp = FindNucSeqEntry (sep);
  if (tmp == NULL) return;
  if (! IS_Bioseq (tmp)) return;
  bsp = (BioseqPtr) tmp->data.ptrvalue;
  if (bsp == NULL) return;
  prev = &(bssp->descr);
  sdp1 = bssp->descr;
  while (sdp1 != NULL && sdp1->choice != Seq_descr_source) {
    prev = &(sdp1->next);
    sdp1 = sdp1->next;
  }
  if (sdp1 == NULL) return;
  sdp2 = SeqEntryGetSeqDescr (tmp, Seq_descr_source, NULL);
  if (sdp2 == NULL) return;
  biop1 = (BioSourcePtr) sdp1->data.ptrvalue;
  biop2 = (BioSourcePtr) sdp2->data.ptrvalue;
  if (CmpOrgById (biop1, biop2)) {
    *prev = sdp1->next;
    sdp1->next = NULL;
    SeqDescrFree (sdp1);
  }
}

static void AssignOneProtein 
(SeqEntryPtr      prot_sep, 
 SequencesFormPtr sqfp,
 SeqEntryPtr      assign_sep,
 SeqLocPtr        use_this,
 BioseqPtr        nucbsp,
 Int2             code,
 Boolean          makeMRNA)
{
  MolInfoPtr        mip;
  SeqEntryPtr       msep = NULL;
  BioseqPtr         protbsp;
  SeqLocPtr         slp;
  Int4              mindiff;
  Boolean           partialN;
  Boolean           partialC;
  ValNodePtr        vnp;
  
  if (prot_sep == NULL)
  {
    return;
  }
  
  mip = MolInfoNew ();
  if (mip != NULL) {
    mip->biomol = 8;
    if (sqfp == NULL) {
      /* no technique */
    } else if (GetStatus (sqfp->protTechBoth)) {
      mip->tech = 10;
    } else {
      mip->tech = 13;
    }
    if (sqfp == NULL) {
      if (use_this == NULL) {
        partialN = FALSE;
        partialC = FALSE;
      } else {
        CheckSeqLocForPartial (use_this, &partialN, &partialC);
      }
    } else {
      partialN = GetStatus (sqfp->partialN);
      partialC = GetStatus (sqfp->partialC);
    }
    if (partialN && partialC) {
      mip->completeness = 5;
    } else if (partialN) {
      mip->completeness = 3;
    } else if (partialC) {
      mip->completeness = 4;
    }
    vnp = CreateNewDescriptor (prot_sep, Seq_descr_molinfo);
    if (vnp != NULL) {
      vnp->data.ptrvalue = (Pointer) mip;
    }
  }
  if (assign_sep != NULL) {
    if (sqfp != NULL && sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      ClearBatchSuggestNucleotide ();
      msep = FindRnaByRefOnProtein (assign_sep, prot_sep);
      if (msep == NULL) {
        msep = FindRnaByRefOnRna (assign_sep, prot_sep);
      }
      if (msep == NULL && nucbsp != NULL && IS_Bioseq (prot_sep)) {
        protbsp = (BioseqPtr) prot_sep->data.ptrvalue;
        if (protbsp != NULL) {
          slp = PredictCodingRegion (nucbsp, protbsp, code);
          if (slp != NULL) {
            mindiff = INT4_MAX;
            FindRnaByLocationOverlap (assign_sep, slp, &mindiff, &msep);
          }
          SeqLocFree (slp);
        }
      }
    }
    if (msep != NULL) {
      msep = GetBestTopParentForDataEx (ObjMgrGetEntityIDForChoice (msep),
                                        (BioseqPtr) msep->data.ptrvalue, TRUE);
    }
    if (msep == NULL) {
      msep = assign_sep;
      if (IS_Bioseq (msep))
      {
        msep = GetBestTopParentForDataEx (ObjMgrGetEntityIDForChoice (msep),
                                          (BioseqPtr) msep->data.ptrvalue, TRUE);
      }
    }
    AddSeqEntryToSeqEntry (msep, prot_sep, TRUE);
    AutomaticProteinProcess (msep, prot_sep, code, makeMRNA, use_this);
  } else {
    AutomaticProteinProcess (assign_sep, prot_sep, code, makeMRNA, use_this);
  }  
}

static SeqEntryPtr FindSeqEntryWithTranscriptID (SeqEntryPtr sep, CharPtr transcript_id)
{
  SeqEntryPtr  found_sep = NULL;
  BioseqPtr    nbsp;
  SeqIdPtr     sip, sip_next;
  CharPtr      tmp;
  BioseqSetPtr bssp;
  
  if (IS_Bioseq (sep))
  {
    nbsp = sep->data.ptrvalue;
    for (sip = nbsp->id; sip != NULL && found_sep == NULL; sip = sip_next)
    {
      sip_next = sip->next;
      sip->next = NULL;
      tmp = SeqIdWholeLabel (sip, PRINTID_REPORT);
      sip->next = sip_next;
      if (StringCmp (tmp, transcript_id) == 0)
      {
        found_sep = sep;
      }
      tmp = MemFree (tmp);
    }
  }
  else
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL && found_sep == NULL; sep = sep->next)
    {
      found_sep = FindSeqEntryWithTranscriptID (sep, transcript_id);
    }
  }
  return found_sep; 
}

/* This section of code is used for matching up proteins to coding region locations
 * on the nucleotide sequences.
 */

/* A ValNode list will be used to hold the list of pairings between protein and nucleotide
 * sequences.  There will be one ValNode per protein sequence.  The choice for the ValNode
 * indicates the position of the nucleotide sequence in the set plus one - a zero indicates
 * that there is no nucleotide for this protein.  The data.ptrvalue will be used to hold the
 * location of the coding region on the nucleotide.
 */

/* This function frees the AssociationList. */ 
extern NucProtAssocPtr FreeAssociationList (NucProtAssocPtr assoc_list)
{
  if (assoc_list == NULL)
  {
    return NULL;
  }
  assoc_list->next = FreeAssociationList (assoc_list->next);
  assoc_list->loc = SeqLocFree (assoc_list->loc);
  assoc_list = MemFree (assoc_list);
  return assoc_list;
}

static NucProtAssocPtr NewAssociationList (NucProtAssocPtr PNTR assoc_list, Int4 position, SeqLocPtr loc)
{
  NucProtAssocPtr last = NULL;
  NucProtAssocPtr new_assoc = (NucProtAssocPtr) MemNew (sizeof (NucProtAssocData));
  
  if (assoc_list == NULL) {
    return NULL;
  }
  if (new_assoc != NULL) {
    new_assoc->position = position;
    new_assoc->loc = loc;
    new_assoc->next = NULL;
    if (*assoc_list == NULL) {
      *assoc_list = new_assoc;
    } else {
      last = *assoc_list;
      while (last->next != NULL) {
        last = last->next;
      }
      last->next = new_assoc;
    }
  }
  return *assoc_list; 
}

/* This function copies the AssociationList */
static NucProtAssocPtr CopyAssociationList (NucProtAssocPtr orig_assoc_list)
{
  NucProtAssocPtr copy_assoc_list = NULL;
  
  if (orig_assoc_list == NULL)
  {
    return NULL;
  }
  copy_assoc_list = (NucProtAssocPtr) MemNew (sizeof (NucProtAssocData));
  if (copy_assoc_list != NULL)
  {
    copy_assoc_list->position = orig_assoc_list->position;
    copy_assoc_list->loc = SeqLocCopy (orig_assoc_list->loc);
    copy_assoc_list->next = CopyAssociationList (orig_assoc_list->next);
  }
  
  return copy_assoc_list;
}


/* This function determines whether all proteins have been assigned to 
 * nucleotide sequences. 
 */
static Boolean AllLocationsProvided (NucProtAssocPtr vnp)
{
  if (vnp == NULL)
  {
    return FALSE;
  }
  while (vnp != NULL)
  {
    if (vnp->position == 0)
    {
      return FALSE;
    }
    vnp = vnp->next;
  }
  return TRUE;
}

/* This function determines whether any proteins have been assigned to 
 * nucleotide sequences. 
 */
static Boolean AnyLocationsProvided (NucProtAssocPtr vnp)
{
  if (vnp == NULL)
  {
    return FALSE;
  }
  while (vnp != NULL)
  {
    if (vnp->position != 0)
    {
      return TRUE;
    }
    vnp = vnp->next;
  }
  return FALSE;
}

/* Given a nucleotide-protein pair, this function calculates a coding region location
 * using Suggest Intervals.  If no location is found, a location that includes the
 * entire sequence is returned instead.
 */
static SeqLocPtr DefaultPairInterval (BioseqPtr nbsp, BioseqPtr pbsp, Int2 code)
{
  SeqLocPtr slp;
  ErrSev    oldsev;
  Char      prot_str[3];
  Boolean   partial5 = FALSE, partial3 = FALSE;
  
  if (nbsp == NULL || pbsp == NULL)
  {
    return NULL;
  }
    
  /* need to suppress errors */  
  oldsev = ErrSetMessageLevel (SEV_MAX);

  /* try to get location using SuggestIntervals */
  SetBatchSuggestNucleotide (nbsp, code);
  slp = PredictCodingRegion (nbsp, pbsp, code);
  ClearBatchSuggestNucleotide ();  
  
  ErrSetMessageLevel (oldsev);
  
  /* if no location, use entire sequence */
  if (slp == NULL)
  {
    return slp; 
  }

  /* check for start and stop codons */  
  SeqPortStreamInt (pbsp, 0, 1, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) (prot_str), NULL);
  if (prot_str[0] != 'M') {
    partial5 = TRUE;
  }
  
  if (SeqLocLen (slp) / 3 != pbsp->length + 1) {
    partial3 = TRUE;
  }
  
  SetSeqLocPartial (slp, partial5, partial3);
  
  return slp;
}


static Boolean FindFeaturesInIdenticalRegions (NucProtAssocPtr assoc_list)
{
  Char       path [PATH_MAX];
  FILE       *fp;
  NucProtAssocPtr   vnp;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;  
  Char        id_txt [128];
  Boolean     found_any = FALSE;

  if (assoc_list == NULL)
  {
    return FALSE;
  }
  
  TmpNam (path);
  fp = FileOpen (path, "wb");
  
  for (vnp = assoc_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->loc != NULL) {
      sfp = SeqMgrGetOverlappingCDS (vnp->loc, &fcontext);
      if (sfp != NULL && SeqLocCompare (vnp->loc, sfp->location) == SLC_A_EQ_B) {
        if (fp == NULL) {
          return TRUE;
        } else {
          found_any = TRUE;
          SeqIdWrite (SeqLocId (vnp->loc), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);          
          fprintf (fp, "%s\n", id_txt);
        }
      }
    }
  }
  FileClose (fp);

  if (found_any) {
    LaunchGeneralTextViewer (path, "Sequences with pre-existing Coding Regions");
  }
  FileRemove (path);   
  return found_any;
}


static Int2 GetGeneticCodeFromBioseq (BioseqPtr bsp)
{
  Int2 code = 1;
  Int2 which_code, read_code;
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  CharPtr location, genetic_code;
  BioSourcePtr biop;

  if (bsp != NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
    if (sdp != NULL) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop->genome == GENOME_mitochondrion) {
        if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->mgcode != 0) {
          code = biop->org->orgname->mgcode;
        }
      } else if (biop->genome == GENOME_plastid) {
        if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->pgcode != 0) {
          code = biop->org->orgname->pgcode;
        } else {
          code = 11;
        }
      } else if (biop->org != NULL && biop->org->orgname != NULL && biop->org->orgname->gcode != 0) {
        code = biop->org->orgname->gcode;
      }
    } else {        
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);
      if (sdp != NULL) {
        location = FindValueFromPairInDeflineBeforeCharPtr ("location", sdp->data.ptrvalue, NULL);
        if (!StringHasNoText (location)) {      
          which_code = UseGeneticCodeForLocation (location);
          if (which_code == USE_MITOCHONDRIAL_GENETIC_CODE) {
            genetic_code = FindValueFromPairInDeflineBeforeCharPtr ("mgcode", sdp->data.ptrvalue, NULL);
            if (genetic_code == NULL) {
              genetic_code = FindValueFromPairInDeflineBeforeCharPtr ("gcode", sdp->data.ptrvalue, NULL);
            }
            if (genetic_code == NULL) {
              genetic_code = FindValueFromPairInDeflineBeforeCharPtr ("genetic_code", sdp->data.ptrvalue, NULL);
            }
            if (genetic_code != NULL) {
              read_code = GeneticCodeFromString (genetic_code);
              genetic_code = MemFree (genetic_code);
              if (read_code > 0) {
                code = read_code;
              }
            }
          } else if (which_code == USE_OTHER_GENETIC_CODE) {
            /* plastid */
            code = 11;
          } else {
            genetic_code = FindValueFromPairInDeflineBeforeCharPtr ("genetic_code", sdp->data.ptrvalue, NULL);
            read_code = GeneticCodeFromString (genetic_code);
            genetic_code = MemFree (genetic_code);
            if (read_code > 0) {
              code = read_code;
            }
          }
        }
        location = MemFree (location);
      }
    }
  }
  return code;
}


/* This function takes a ValNode list where each ValNode represents
 * a protein in prot_list (in order).  The choice for each ValNode
 * represents the position of the chosen nucleotide in the nuc_list
 * (position includes segments in segmented sets, which is why 
 * FindNthSequenceInSet is used) plus one - zero indicates that there
 * is no nucleotide sequence for this protein.
 * The data.ptrvalue for the ValNode is to be populated with a
 * coding region SeqLoc, or NULL if there is no nucleotide for the protein.
 */
static Boolean 
PickCodingRegionLocationsForProteinNucleotidePairs 
(NucProtAssocPtr  assoc_list,
 SeqEntryPtr nuc_list,
 SeqEntryPtr prot_list)
{
  NucProtAssocPtr vnp_assoc;
  Int4       data_row;
  BioseqPtr  nbsp, pbsp;
  SeqLocPtr  slp;
  Char       path [PATH_MAX];
  FILE       *fp;
  Boolean    errors_found = FALSE;
  Char       n_idstr[128];
  Char       p_idstr[128];
  Int2       code;
  
  if (assoc_list == NULL || nuc_list == NULL || prot_list == NULL)
  {
    return FALSE;
  }

  TmpNam (path);
  fp = FileOpen (path, "wb");

  vnp_assoc = assoc_list;
  for (data_row = 0, vnp_assoc = assoc_list;
       vnp_assoc != NULL; 
       data_row++, vnp_assoc = vnp_assoc->next)
  {
    if (vnp_assoc->position > 0)
    {
      nbsp = FindNthSequenceInSet (nuc_list, vnp_assoc->position - 1, NULL, TRUE);
      pbsp = FindNthSequenceInSet (prot_list, data_row, NULL, FALSE);
      if (nbsp == NULL || pbsp == NULL) {
        slp = NULL;
      } else if ((nbsp->length +1) / 3 < pbsp->length) {
        if (fp != NULL) {
          SeqIdWrite (SeqIdFindWorst (nbsp->id), n_idstr, PRINTID_REPORT,
                      sizeof (n_idstr));
          SeqIdWrite (SeqIdFindWorst (pbsp->id), p_idstr, PRINTID_REPORT,
                      sizeof (p_idstr));
          fprintf (fp, "%s is too short to encode %s\n", n_idstr, p_idstr);
        }
        vnp_assoc->position = 0;
        errors_found = TRUE;
        slp = NULL;
      } else {
        code = GetGeneticCodeFromBioseq (nbsp);
        slp = DefaultPairInterval (nbsp, pbsp, code);
        if (slp == NULL) {
          errors_found = TRUE;
          if (fp != NULL) {
            SeqIdWrite (SeqIdFindWorst (nbsp->id), n_idstr, PRINTID_REPORT,
                        sizeof (n_idstr));
            SeqIdWrite (SeqIdFindWorst (pbsp->id), p_idstr, PRINTID_REPORT,
                        sizeof (p_idstr));
            fprintf (fp, "Unable to determine coding region location on %s for %s\n", n_idstr, p_idstr); 
          }
        }
      }
    }
    else
    {
      slp = NULL;
    }
    vnp_assoc->loc = SeqLocFree (vnp_assoc->loc);
    vnp_assoc->loc = slp;
  }
  
  FileClose (fp);
  if (errors_found) {
    LaunchGeneralTextViewer (path, "Nucleotide-Protein Mismatches");
  }
  FileRemove (path);
  return !errors_found;
}


static Int2 FindGeneticCodeForBioseq (BioseqPtr bsp, Int2 default_code)
{
  Int2         code = default_code;
  BioSourcePtr biop;
  SeqEntryPtr  nsep;
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp = NULL;

  if (bsp == NULL) return default_code;
  nsep = GetBestTopParentForData (ObjMgrGetEntityIDForPointer (bsp), bsp);
  if (nsep == NULL || nsep->data.ptrvalue == NULL) return default_code;
  if (nsep->choice == 1)
  {
    bsp = nsep->data.ptrvalue;
    sdp = bsp->descr;
  }
  else if (nsep->choice == 2)
  {
    bssp = nsep->data.ptrvalue;
    sdp = bssp->descr;
  }
  while (sdp != NULL)
  {
    if (sdp->choice == Seq_descr_source && sdp->data.ptrvalue != NULL)
    {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop->org != NULL && biop->org->orgname != NULL)
      {
        code = BioSourceToGeneticCode (biop);
      }
    }
    sdp = sdp->next;
  }
  return code;
}


/* This function takes a ValNode list of coding region SeqLocs,
 * the list of nucleotide sequences, and the list of protein sequences
 * and creates the nuc-prot sets.
 */
static void 
AssignProteinsToSelectedNucleotides 
(NucProtAssocPtr  assoc_list,
 SeqEntryPtr      nuc_list,
 SeqEntryPtr      prot_list,
 SequencesFormPtr sqfp, 
 Int2             code,
 Boolean          makeMRNA)
{
  SeqEntryPtr prot_sep, nsep, prot_next;
  NucProtAssocPtr  vnp_assoc;
  BioseqPtr   nbsp;
  BioseqPtr PNTR bsp_array;
  Int4           prot_num;
  ValNodePtr     descr = NULL;
  Int2           genCode;

  if (assoc_list == NULL || nuc_list == NULL || prot_list == NULL)
  {
    return;
  }

  /* need to collect bioseqs before we start adding, otherwise the position in
   * the set changes */
  
  bsp_array = (BioseqPtr PNTR) MemNew (ValNodeLen (prot_list) * sizeof (BioseqPtr));
  if (bsp_array == NULL)
  {
    return;
  }
  
  for (prot_num = 0, vnp_assoc = assoc_list;
       vnp_assoc != NULL;
       prot_num++, vnp_assoc = vnp_assoc->next)
  {
    if (vnp_assoc->loc == NULL)
    {
      bsp_array [prot_num] = NULL;
    }
    else
    {
      bsp_array [prot_num] = FindNthSequenceInSet (nuc_list, vnp_assoc->position - 1, NULL, TRUE);
    }
  }
  
  for (prot_sep = prot_list, vnp_assoc = assoc_list, prot_num = 0;
       prot_sep != NULL && vnp_assoc != NULL;
       prot_sep = prot_next, vnp_assoc = vnp_assoc->next, prot_num++)
  {
    prot_next = prot_sep->next;
    prot_sep->next = NULL;
    
    if (vnp_assoc->loc == NULL)
    {
      /* discard protein */
      if (IS_Bioseq (prot_sep))
      {
        SeqMgrDeleteFromBioseqIndex (prot_sep->data.ptrvalue);
      }
      prot_sep = SeqEntryFree (prot_sep);
    }
    else
    {
      nbsp = bsp_array [prot_num];
      nsep = SeqMgrGetSeqEntryForData (nbsp);
      if (nbsp != NULL && nbsp->repr == Seq_repr_seg)
      {
        nsep = GetBestTopParentForData (ObjMgrGetEntityIDForPointer (nbsp), nbsp); 
      }
      genCode = FindGeneticCodeForBioseq (nbsp, code);
      if (nsep != NULL && nsep->data.ptrvalue == nbsp) {
        descr = ExtractBioSourceAndPubs (nsep);
      }
      AssignOneProtein (prot_sep, sqfp, nsep, vnp_assoc->loc, nbsp, 
                        genCode, makeMRNA);
      if (descr != NULL) {
        ReplaceBioSourceAndPubs (nsep, descr);
      }
      vnp_assoc->loc = NULL; /*SeqLoc was freed in AssignOneProtein */
    }
  }
  
  bsp_array = MemFree (bsp_array);
}

/* This function creates a new protein ID based on the nucleotide ID that will be
 * unique within the record - nucleotide and protein sequence IDs are checked
 * for matches.
 */
static CharPtr 
BuildProteinIDUniqueInIDAndTitleEdit 
(CharPtr nuc_id,
 IDAndTitleEditPtr iatep_nuc,
 IDAndTitleEditPtr iatep_prot)
{
  CharPtr new_id, cp;
  Int4    offset, seq_num;
  Boolean unique_found = FALSE;
  
  if (iatep_nuc == NULL || iatep_prot == NULL || StringHasNoText (nuc_id))
  {
    return NULL;
  }
  
  new_id = (CharPtr) MemNew ((StringLen (nuc_id) + 20) * sizeof (Char));
  if (new_id != NULL)
  {
    StringCpy (new_id, nuc_id);
    StringCat (new_id, "_");
    cp = new_id + StringLen (new_id);
    for (offset = 1; offset < INT4_MAX && ! unique_found; offset ++)
    {
      sprintf (cp, "%d", offset);
      unique_found = TRUE;
      for (seq_num = 0; seq_num < iatep_nuc->num_sequences && unique_found; seq_num++)
      {
        if (StringCmp (iatep_nuc->id_list [seq_num], new_id) == 0)
        {
          unique_found = FALSE;
        }
      }
      for (seq_num = 0; seq_num < iatep_prot->num_sequences && unique_found; seq_num++)
      {
        if (StringCmp (iatep_prot->id_list [seq_num], new_id) == 0)
        {
          unique_found = FALSE;
        }
      }
    }
  }
  if (unique_found)
  {
    return new_id;
  }
  else
  {
    new_id = MemFree (new_id);
    return StringSave ("too_many");
  }
}

static Boolean DoIdsMatch (CharPtr id1, CharPtr id2)
{
  CharPtr tmp1, cp1 = NULL;
  CharPtr tmp2, cp2 = NULL;
  Boolean match = FALSE;
  
  tmp1 = StringChr (id1, '|');
  if (tmp1 == NULL) {
    tmp1 = id1;
  } else if (tmp1 == id1 + 2) {
    tmp1++;
    cp1 = StringChr (tmp1, '|');
    if (cp1 != NULL) {
      *cp1 = 0;
    }
  }
  
  tmp2 = StringChr (id2, '|');
  if (tmp2 == NULL) {
    tmp2 = id2;
  } else if (tmp2 == id2 + 2) {
    tmp2++;
    cp2 = StringChr (tmp2, '|');
    if (cp2 != NULL) {
      *cp2 = 0;
    }
  }
 
  if (StringCmp (tmp1, tmp2) == 0)
  {
    match = TRUE;
  }
  if (cp1 != NULL) {
    *cp1 = '|';
  }
  if (cp2 != NULL) {
    *cp2 = '|';
  }
  return match;
}


/* if the user gave the protein sequences the same IDs as the nucleotide sequences,
 * we need to create new sequence IDs for the proteins so that they will be unique.
 * We should also make sure that sequence IDs that don't match nucleotide sequence
 * IDs are unique.
 */
static void ReplaceDuplicateProteinIDs (SeqEntryPtr nuc_list, SeqEntryPtr prot_list)
{
  Int4              nuc_seq_num, prot_seq_num, prot_seq_num_check;
  IDAndTitleEditPtr iatep_nuc, iatep_prot;
  Boolean           found_nuc_match;
  CharPtr           tmp_str, cp;
  BioseqPtr         prot_bsp, nuc_bsp;
  
  if (nuc_list == NULL || prot_list == NULL)
  {
    return;
  }
  
  iatep_nuc = SeqEntryListToIDAndTitleEditEx (nuc_list, TRUE);
  iatep_prot = SeqEntryListToIDAndTitleEdit (prot_list);
  if (iatep_nuc != NULL && iatep_prot != NULL)
  {
    for (prot_seq_num = 0; prot_seq_num < iatep_prot->num_sequences; prot_seq_num++)
    {
      /* This part replaces any protein sequence IDs that match a nucleotide ID with
       * the nucleotide ID plus an underscore plus a number that makes the ID
       * unique.
       */
      found_nuc_match = FALSE;
      prot_bsp = FindNthSequenceInSet (prot_list, prot_seq_num, &(iatep_prot->is_seg[prot_seq_num]), FALSE);
      if (prot_bsp == NULL) continue;
      for (nuc_seq_num = 0;
           nuc_seq_num < iatep_nuc->num_sequences && ! found_nuc_match; 
           nuc_seq_num++)
      {
        nuc_bsp = FindNthSequenceInSet (nuc_list, nuc_seq_num, &(iatep_nuc->is_seg[prot_seq_num]), TRUE);
        if (nuc_bsp == NULL) continue;
      
        if (SeqIdIn (prot_bsp->id, nuc_bsp->id) || RelaxedSeqIdIn (prot_bsp->id, nuc_bsp->id)
            || DoIdsMatch (iatep_nuc->id_list [nuc_seq_num],
                           iatep_prot->id_list [prot_seq_num])) {
          tmp_str = iatep_nuc->id_list [nuc_seq_num];
          cp = StringChr (tmp_str, '|');
          if (cp == tmp_str + 2) {
            tmp_str += 3;
            cp = StringChr (tmp_str, '|');
          }
          if (cp != NULL) {
            *cp = 0;
          }
                        
          iatep_prot->id_list [prot_seq_num] = MemFree (iatep_prot->id_list [prot_seq_num]);
          iatep_prot->id_list [prot_seq_num] = BuildProteinIDUniqueInIDAndTitleEdit (tmp_str,
                                                                                     iatep_nuc,
                                                                                     iatep_prot);
          if (cp != NULL) {
            *cp = '|';
          }
          found_nuc_match = TRUE;
        }
      }
      /* This part replaces a protein sequence ID that matches a previous protein
       * sequence ID with the original protein sequence ID plus an underscore plus
       * a number that makes the ID unique.
       */
      if (!found_nuc_match)
      {
        for (prot_seq_num_check = prot_seq_num + 1; 
             prot_seq_num_check < iatep_prot->num_sequences;
             prot_seq_num_check ++)
        {
          if (StringCmp (iatep_prot->id_list [prot_seq_num],
                         iatep_prot->id_list [prot_seq_num_check]) == 0)
          {
            tmp_str = iatep_prot->id_list [prot_seq_num_check];
            cp = StringChr (tmp_str, '|');
            if (cp == tmp_str + 2) {
              tmp_str += 3;
              cp = StringChr (tmp_str, '|');
            }
            if (cp != NULL) {
              *cp = 0;
            }
            tmp_str = StringSave (tmp_str);
          
            iatep_prot->id_list [prot_seq_num_check] = MemFree (iatep_prot->id_list [prot_seq_num_check]);
            iatep_prot->id_list [prot_seq_num_check] = BuildProteinIDUniqueInIDAndTitleEdit (tmp_str,
                                                                                             iatep_nuc,
                                                                                             iatep_prot);
            tmp_str = MemFree (tmp_str);
          }
        }
      }
    }
  }
  ApplyIDAndTitleEditToSeqEntryList (prot_list, iatep_prot);
  iatep_prot = IDAndTitleEditFree (iatep_prot);
  iatep_nuc = IDAndTitleEditFree (iatep_nuc);
}

static Uint2 nucprotedit_types [] = {
  TAGLIST_PROMPT, TAGLIST_PROMPT, TAGLIST_POPUP, TAGLIST_TEXT, TAGLIST_TEXT
};

static Uint2 nucprotedit_widths [] = {
  5, 20, 10, 15, 15
};

#define NUCPROTEDIT_NUCID_COLUMN 2
#define NUCPROTEDIT_GENE_COLUMN  3
#define NUCPROTEDIT_PROT_COLUMN  4

typedef struct nucprotedit
{
  SeqEntryPtr nuc_list;
  SeqEntryPtr prot_list;
  DialoG      dlg;
  ButtoN      accept_btn;
  NucProtAssocPtr  assoc_list;
  TexT        all_gene_txt;
  TexT        all_prot_txt;
} NucProtEditData, PNTR NucProtEditPtr;

static void PopulateNucProtEdit (NucProtEditPtr npep)
{
  IDAndTitleEditPtr     iatep_nuc, iatep_prot;
  ValNodePtr            row_list = NULL;
  NucProtAssocPtr       vnp_assoc;
  TagListPtr            tlp;
  CharPtr               data_string, gene_locus, prot_name;
  Int4                  data_len;
  Int4                  prot_num;
  Int4                  old_scroll_pos = 0;
  
  if (npep == NULL)
  {
    return;
  }
  
  tlp = (TagListPtr) GetObjectExtra (npep->dlg);
  if (tlp == NULL)
  {
    return;
  }
  
  /* need to get bar value and reset after populating */
  if (tlp->bar != NULL)  
  {
    old_scroll_pos = GetBarValue (tlp->bar);
  }
  
  iatep_nuc = SeqEntryListToIDAndTitleEditEx (npep->nuc_list, TRUE);
  iatep_prot = SeqEntryListToIDAndTitleEdit (npep->prot_list);
  if (iatep_nuc != NULL && iatep_prot != NULL)
  {
    vnp_assoc = npep->assoc_list;
    for (prot_num = 0; prot_num < iatep_prot->num_sequences; prot_num++)
    {
      /* first column is protein ID */
      /* second column is choice for nucleotide ID */
      /* third column is gene locus tag */
      /* fourth column is protein name */
      /* fifth column indicates presence of suggested interval */
      gene_locus = FindValueFromPairInDefline ("gene", iatep_prot->title_list [prot_num]);
      prot_name = FindValueFromPairInDefline ("protein", iatep_prot->title_list [prot_num]);
      
      data_len = StringLen (iatep_prot->id_list [prot_num])
                  + 20
                  + StringLen (gene_locus)
                  + StringLen (prot_name);
      data_string = (CharPtr) MemNew (data_len * sizeof (Char));                  
      if (data_string != NULL)
      {
        sprintf (data_string, "%d\t%s\t%d\t%s\t%s\n",
                               prot_num + 1,
                               iatep_prot->id_list [prot_num],
                               vnp_assoc == NULL ? 0 : vnp_assoc->position,
                               gene_locus == NULL ? "" : gene_locus,
                               prot_name == NULL ? "" : prot_name);
        ValNodeAddPointer (&row_list, 0, data_string);                               
      }
      gene_locus = MemFree (gene_locus);
      prot_name = MemFree (prot_name);
      if (vnp_assoc != NULL)
      {
        vnp_assoc = vnp_assoc->next;
      }
    }
    SendMessageToDialog (npep->dlg, VIB_MSG_RESET);
    tlp->vnp = row_list;

    if (iatep_prot->num_sequences > tlp->rows)
    {
      tlp->max = MAX ((Int2) 0, (Int2) (iatep_prot->num_sequences - tlp->rows));  
      CorrectBarMax (tlp->bar, tlp->max);
      CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1); 
      Enable (tlp->bar);
      SetBarValue (tlp->bar, old_scroll_pos);
    }
    else
    {
      Hide (tlp->bar);
    }
    SendMessageToDialog (npep->dlg, VIB_MSG_REDRAW);    
  }
  
  iatep_nuc = IDAndTitleEditFree (iatep_nuc);
  iatep_prot = IDAndTitleEditFree (iatep_prot);
}

static CharPtr 
GetTagListValueEx (TagListPtr tlp, Int4 seq_num, Int4 col_num);

static void ApplyGeneNameToAllSequences (ButtoN b)
{
  NucProtEditPtr npep;
  CharPtr        all_gene_name, new_val;
  TagListPtr     tlp;
  Int4           seq_num;
  ValNodePtr     vnp;
  
  npep = (NucProtEditPtr) GetObjectExtra (b);
  if (npep == NULL)
  {
    return;
  }
  
  tlp = (TagListPtr) GetObjectExtra (npep->dlg);
  if (tlp == NULL)
  {
    return;
  }
  all_gene_name = SaveStringFromText (npep->all_gene_txt);
  if (ANS_YES == Message (MSG_YN, "Are you sure you want to set all of the gene locus values to %s?",
                          all_gene_name))
  {
    for (vnp = tlp->vnp, seq_num = 0;
         vnp != NULL;
         vnp = vnp->next, seq_num++)
    {
      new_val = ReplaceTagListColumn (vnp->data.ptrvalue, all_gene_name, NUCPROTEDIT_GENE_COLUMN);
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = new_val;
    }
    SendMessageToDialog (npep->dlg, VIB_MSG_REDRAW);
  }
  all_gene_name = MemFree (all_gene_name);
}

static void ApplyProteinNameToAllSequences (ButtoN b)
{
  NucProtEditPtr npep;
  CharPtr        all_prot_name, new_val;
  TagListPtr     tlp;
  Int4           seq_num;
  ValNodePtr     vnp;
  
  npep = (NucProtEditPtr) GetObjectExtra (b);
  if (npep == NULL)
  {
    return;
  }
  
  tlp = (TagListPtr) GetObjectExtra (npep->dlg);
  if (tlp == NULL)
  {
    return;
  }
  all_prot_name = SaveStringFromText (npep->all_prot_txt);
  if (ANS_YES == Message (MSG_YN, "Are you sure you want to set all of the protein names to %s?",
                          all_prot_name))
  {
    for (vnp = tlp->vnp, seq_num = 0;
         vnp != NULL;
         vnp = vnp->next, seq_num++)
    {
      new_val = ReplaceTagListColumn (vnp->data.ptrvalue, all_prot_name, NUCPROTEDIT_PROT_COLUMN);
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = new_val;
    }
    SendMessageToDialog (npep->dlg, VIB_MSG_REDRAW);   
  }
  all_prot_name = MemFree (all_prot_name);
}

static void ApplyNucProtEditGeneAndProt (NucProtEditPtr npep)
{
  TagListPtr tlp;
  Int4       seq_num;
  CharPtr    gene, prot;
  IDAndTitleEditPtr iatep;
  ValNodePtr        vnp;
  
  if (npep == NULL)
  {
    return;
  }
  
  tlp = (TagListPtr) GetObjectExtra (npep->dlg);
  if (tlp == NULL)
  {
    return;
  }
  
  iatep = SeqEntryListToIDAndTitleEdit (npep->prot_list);
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0, vnp = tlp->vnp; 
       seq_num < iatep->num_sequences && vnp != NULL; 
       seq_num++, vnp = vnp->next)
  {
    gene = GetTagListValueEx (tlp, seq_num, NUCPROTEDIT_GENE_COLUMN);
    iatep->title_list [seq_num] = ReplaceValueInOneDefLine (iatep->title_list [seq_num],
                                                            "gene",
                                                            gene);
    gene = MemFree (gene);
    prot = GetTagListValueEx (tlp, seq_num, NUCPROTEDIT_PROT_COLUMN);
    iatep->title_list [seq_num] = ReplaceValueInOneDefLine (iatep->title_list [seq_num],
                                                            "protein",
                                                            prot);
  }
  ApplyIDAndTitleEditToSeqEntryList (npep->prot_list, iatep);
  iatep = IDAndTitleEditFree (iatep);
}

/* This function collects the pairings of nucleotides and proteins from
 * the NucProtEdit dialog.
 * The ValNode list stored in npep->assoc_list has one ValNode for each
 * protein in npep->prot_list.  The choice for each ValNode is the position
 * of the nucleotide in npep->nuc_list plus one - zero indicates that there
 * is no nucleotide for this protein sequence.
 */
static void CollectSequenceAssociationsFromNucProtEdit (NucProtEditPtr npep)
{
  TagListPtr tlp;
  CharPtr    str;
  Int4       num_data_rows, assoc_num, data_row;
  NucProtAssocPtr assoc_list = NULL;
  
  if (npep == NULL)
  {
    return;
  }
  tlp = (TagListPtr) GetObjectExtra (npep->dlg);
  if (tlp == NULL)
  {
    return;
  }
  
  num_data_rows = ValNodeLen (tlp->vnp);
  for (data_row = 0; data_row < num_data_rows; data_row++)
  {
    str  = GetTagListValueEx (tlp, data_row, NUCPROTEDIT_NUCID_COLUMN);
    if (!StringHasNoText (str))
    {
      assoc_num = atoi (str);
    }
    else
    {
      assoc_num = 0;
    }
    str = MemFree (str);
        
    NewAssociationList (&assoc_list, assoc_num, NULL);
  }
  
  npep->assoc_list = FreeAssociationList (npep->assoc_list);
  npep->assoc_list = assoc_list;    
}

static CharPtr GetNumSuffix (Int4 i) 
{
  if (i == 1) {
    return "st";
  } else if (i == 2) {
    return "nd";
  } else if (i == 3) {
    return "rd";
  } else {
    return "th";
  }
}

static CharPtr GetNucProtAssocMessage (NucProtAssocPtr assoc_list)
{
  NucProtAssocPtr   vnp_assoc;
  Int4              num = 0, pos;
  CharPtr           msg = NULL;
  CharPtr           msg_fmt = "The %s protein sequence%s could not be associated with %snucleotide sequence%s.  You must scroll to %s and select the correct nucleotide sequence ID%s for %s.  Clicking Accept at this point will discard %s.";
  Int4              extra_space = 2 * StringLen ("these protein sequences") + StringLen ("these positions");
  Int4              msg_len = 0;
  ValNodePtr        desc_list = NULL, vnp;
  Char              desc[18];
  CharPtr           num_list = NULL;
  Int4              num_list_len = 0;

  if (assoc_list == NULL) {
    return StringSave ("You must select nucleotide sequences for your protein sequences.");
  } else {
    for (vnp_assoc = assoc_list, pos = 1; vnp_assoc != NULL; vnp_assoc = vnp_assoc->next, pos++) {
      if (vnp_assoc->position == 0) {
        sprintf (desc, "%d%s", pos, GetNumSuffix(pos));
        ValNodeAddPointer (&desc_list, 0, StringSave (desc));
        num_list_len += StringLen (desc) + 2;
        num++;
      }
    }
    if (num == 0) {
      return NULL;
    } 
    if (num > 1) {
      num_list_len += 6;
    }
    num_list = (CharPtr) MemNew (sizeof (Char) * num_list_len);
    for (vnp = desc_list; vnp != NULL; vnp = vnp->next) {
      StringCat (num_list, vnp->data.ptrvalue);
      if (vnp->next != NULL) {
        if (num == 2) {
          StringCat (num_list, " and ");
        } else if (vnp->next->next == NULL) {
          StringCat (num_list, ", and ");
        } else {
          StringCat (num_list, ", ");
        }
      }
    }
    desc_list = ValNodeFreeData (desc_list);
    
    msg_len = StringLen (msg_fmt) + extra_space + num_list_len;

    msg = (CharPtr) MemNew (sizeof (Char) * (msg_len));

    if (num == 1) {
      sprintf (msg, msg_fmt, num_list, "",  "a ", "",  "this position",   "",  "this protein sequence",  "this protein sequence");
    } else {
      sprintf (msg, msg_fmt, num_list, "s", "",   "s", "these positions", "s", "these protein sequences", "these protein sequences"); 
    }
  }
  return msg;  
}

/* This function produces a dialog that allows the user to edit the gene and protein names
 * and to select the locations for the coding regions for each protein sequence.
 */
static NucProtAssocPtr 
CollectNucleotideProteinAssociations 
(SeqEntryPtr nuc_list,
 SeqEntryPtr prot_list,
 NucProtAssocPtr default_assoc_list)
{
  WindoW                w;
  GrouP                 h, title_grp, all_gene_grp, all_prot_grp, k, c, p_msg;
  PrompT                p_prot, p_nuc, p_locus, p_name;
  ButtoN                b;
  Int4                  num_prots;
  Int4                  rows_shown = 0;
  TagListPtr            tlp;
  ModalAcceptCancelData acd;
  NucProtEditData       nped;
  IDAndTitleEditPtr     iatep_nuc;
  EnumFieldAssocPtr     nuc_alist;
  EnumFieldAssocPtr     nucprotedit_alists [] = { NULL, NULL, NULL, NULL, NULL};
  Int4                  nuc_num;
  CharPtr               msg;
  
  if (nuc_list == NULL || prot_list == NULL)
  {
    return NULL;
  }
  
  nped.nuc_list = nuc_list;
  nped.prot_list = prot_list;
  nped.assoc_list = CopyAssociationList (default_assoc_list);
  
  num_prots = ValNodeLen (prot_list);
  rows_shown = MIN (num_prots, 5);
  
  /* set up ALIST for nucleotide list.
   * cannot free IDAndTitleEdit until done with ALIST.
   */
  iatep_nuc = SeqEntryListToIDAndTitleEditEx (nped.nuc_list, TRUE);
  nuc_alist = (EnumFieldAssocPtr) MemNew ((iatep_nuc->num_sequences + 2) * sizeof (EnumFieldAssoc));
  nuc_alist [0].name = "";
  nuc_alist [0].value = 0;
  for (nuc_num = 0; nuc_num < iatep_nuc->num_sequences; nuc_num++)
  {
    nuc_alist [nuc_num + 1].name = iatep_nuc->id_list [nuc_num];
    nuc_alist [nuc_num + 1].value = nuc_num + 1;
  }
  nuc_alist [nuc_num + 1].name = NULL;
  nucprotedit_alists [NUCPROTEDIT_NUCID_COLUMN] = nuc_alist;
  
  w = MovableModalWindow (-20, -13, -10, -10, "Map Proteins to Nucleotides", NULL);
  
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  msg = GetNucProtAssocMessage (nped.assoc_list);
  p_msg = MultiLinePrompt (h, msg, 750, systemFont);

  k = HiddenGroup (h, 2, 0, NULL);
  /* text and button for setting all gene locus values */
  all_gene_grp = HiddenGroup (k, -1, 0, NULL);  
  b = PushButton (all_gene_grp, "Set All Gene Locus Values to Value Below", ApplyGeneNameToAllSequences);
  SetObjectExtra (b, &nped, NULL);
  nped.all_gene_txt = DialogText (all_gene_grp, "", 15, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) nped.all_gene_txt, (HANDLE) b, NULL);

  /* text and button for setting all protein names */
  all_prot_grp = HiddenGroup (k, -1, 0, NULL);  
  b = PushButton (all_prot_grp, "Set All Protein Names to Value Below", ApplyProteinNameToAllSequences);
  SetObjectExtra (b, &nped, NULL);
  nped.all_prot_txt = DialogText (all_prot_grp, "", 15, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) nped.all_prot_txt, (HANDLE) b, NULL);  
  
  title_grp = HiddenGroup (h, 5, 0, NULL);
  SetGroupSpacing (title_grp, 10, 10);
     
  p_prot = StaticPrompt (title_grp, "Prot ID", 0, 0, programFont, 'l');    
  p_nuc = StaticPrompt (title_grp, "Nuc ID", 0, 0, programFont, 'l');    
  p_locus = StaticPrompt (title_grp, "Gene Locus", 0, 0, programFont, 'l');    
  p_name = StaticPrompt (title_grp, "Protein Name", 0, 0, programFont, 'l');    
 
  nped.dlg = CreateTagListDialogEx (h, rows_shown, 5, 2,
                                           nucprotedit_types, nucprotedit_widths,
                                           nucprotedit_alists, TRUE, TRUE, 
                                           NULL, NULL);

  
  tlp = (TagListPtr) GetObjectExtra (nped.dlg);  
  if (tlp == NULL) return NULL;
  
  if (num_prots > rows_shown)
  {
    tlp->max = MAX ((Int2) 0, (Int2) (num_prots - tlp->rows));  
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1); 
    Enable (tlp->bar);
  }
  else
  {
    Hide (tlp->bar);
  }

  c = HiddenGroup (h, 2, 0, NULL);
  nped.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (nped.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) nped.dlg, (HANDLE) c, (HANDLE) NULL);

  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [1], (HANDLE) p_prot, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [NUCPROTEDIT_NUCID_COLUMN], (HANDLE) p_nuc, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [NUCPROTEDIT_GENE_COLUMN], (HANDLE) p_locus, 
                               (HANDLE) all_gene_grp, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [NUCPROTEDIT_PROT_COLUMN], (HANDLE) p_name, 
                               (HANDLE) all_prot_grp, NULL);
  
  PopulateNucProtEdit (&nped);

  Show (w);
  Select (w);
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
    if (acd.accepted)
    {
      CollectSequenceAssociationsFromNucProtEdit (&nped);
      if (!AllLocationsProvided (nped.assoc_list))
      {
        if (ANS_NO == Message (MSG_YN, "You have not provided coding region locations for all of your proteins - these proteins will be discarded.  Are you sure you want to continue?"))
        {
          acd.accepted = FALSE;
        }
      }
      if (acd.accepted) 
      {
        if (!PickCodingRegionLocationsForProteinNucleotidePairs (nped.assoc_list,
                                                                 nuc_list,
                                                                 prot_list)) 
        {
          if (ANS_NO == Message (MSG_YN, "Unable to determine coding region location for some proteins.  These proteins will be discarded.  Do you wish to continue?")) 
          {
            acd.accepted = FALSE;
            PopulateNucProtEdit (&nped);
          }
        }  
        if (acd.accepted && FindFeaturesInIdenticalRegions (nped.assoc_list)) 
        {
          if (ANS_CANCEL == Message (MSG_OKC, "Duplicate coding regions were found on one or more sequences - do you wish to continue?")) 
          {
            acd.accepted = FALSE;
          }
        }
      }      
    }
  }
    
  Remove (w);
  
  if (acd.accepted)
  {
    /* apply any gene protein data the user may have entered to the titles */
    ApplyNucProtEditGeneAndProt (&nped);    
  }
  else
  {
    nped.assoc_list = FreeAssociationList (nped.assoc_list);
  }
  
  nuc_alist = MemFree (nuc_alist);
  iatep_nuc = IDAndTitleEditFree (iatep_nuc);
  
  Update ();
  return nped.assoc_list;
}

/* This function tries to build an association list by matching nucleotide sequence IDs 
 * to protein sequence IDs.
 */
static NucProtAssocPtr 
BuildAssociationListByMatch
(SeqEntryPtr       nuc_list,
 SeqEntryPtr       prot_list)
{
  IDAndTitleEditPtr iatep_nuc, iatep_prot;
  Int4       nuc_seq_num, prot_seq_num, found_num;
  NucProtAssocPtr assoc_list = NULL;
  BioseqPtr       nuc_bsp, prot_bsp;
  
  if (nuc_list == NULL || prot_list == NULL)
  {
    return NULL;
  }
  
  iatep_nuc = SeqEntryListToIDAndTitleEditEx (nuc_list, TRUE);
  iatep_prot = SeqEntryListToIDAndTitleEdit (prot_list);
  if (iatep_nuc == NULL || iatep_prot == NULL)
  {
    iatep_nuc = IDAndTitleEditFree (iatep_nuc);
    iatep_prot = IDAndTitleEditFree (iatep_prot);
    return NULL;
  }
  
  for (prot_seq_num = 0; prot_seq_num < iatep_prot->num_sequences; prot_seq_num++)
  {
    found_num = 0;
    if (iatep_nuc->num_sequences == 1)
    {
      found_num = 1;
    }
    prot_bsp = FindNthSequenceInSet (prot_list, prot_seq_num, &(iatep_prot->is_seg[prot_seq_num]), FALSE);
    if (prot_bsp == NULL) continue;
    for (nuc_seq_num = 0;
         nuc_seq_num < iatep_nuc->num_sequences && found_num == 0;
         nuc_seq_num++)
    {
      nuc_bsp = FindNthSequenceInSet (nuc_list, nuc_seq_num, &(iatep_nuc->is_seg [nuc_seq_num]), TRUE);
      if (SeqIdIn (prot_bsp->id, nuc_bsp->id) || RelaxedSeqIdIn (prot_bsp->id, nuc_bsp->id)
          || DoIdsMatch(iatep_nuc->id_list[nuc_seq_num], iatep_prot->id_list [prot_seq_num])) {
        found_num = nuc_seq_num + 1;
      }
    }
 
    NewAssociationList (&assoc_list, found_num, NULL);
  }
  iatep_nuc = IDAndTitleEditFree (iatep_nuc);
  iatep_prot = IDAndTitleEditFree (iatep_prot);
  return assoc_list;
}

/* This function builds an association list by matching nucleotide sequences 
 * to protein sequences by position.
 */
static NucProtAssocPtr 
BuildAssociationListByPosition
(SeqEntryPtr       nuc_list,
 SeqEntryPtr       prot_list)
{
  IDAndTitleEditPtr iatep_nuc, iatep_prot;
  Int4       nuc_seq_num, prot_seq_num;
  NucProtAssocPtr assoc_list = NULL;
  Int4       num_masters = 0, num_segs = 0;
  
  if (nuc_list == NULL || prot_list == NULL)
  {
    return NULL;
  }
  
  iatep_nuc = SeqEntryListToIDAndTitleEditEx (nuc_list, TRUE);
  iatep_prot = SeqEntryListToIDAndTitleEdit (prot_list);
  if (iatep_nuc == NULL || iatep_prot == NULL)
  {
    iatep_nuc = IDAndTitleEditFree (iatep_nuc);
    iatep_prot = IDAndTitleEditFree (iatep_prot);
    return NULL;
  }
  
  for (nuc_seq_num = 0; nuc_seq_num < iatep_nuc->num_sequences; nuc_seq_num++)
  {
    if (iatep_nuc->is_seg != NULL && iatep_nuc->is_seg [nuc_seq_num])
    {
      num_segs ++;
    }
    else
    {
      num_masters ++;
    }
  }
  
  if (num_segs == iatep_prot->num_sequences && iatep_nuc->is_seg != NULL)
  {
    /* assign proteins to segments */
    nuc_seq_num = 0;
    for (prot_seq_num = 0; prot_seq_num < iatep_prot->num_sequences; prot_seq_num++)
    {
      while (! iatep_nuc->is_seg [nuc_seq_num] && nuc_seq_num < iatep_nuc->num_sequences)
      {
        nuc_seq_num++;
      }
      if (nuc_seq_num < iatep_nuc->num_sequences)
      {
        NewAssociationList (&assoc_list, nuc_seq_num + 1, NULL);
        nuc_seq_num++;
      }
      else
      {
        NewAssociationList (&assoc_list, 0, NULL);
      }
    }
  }
  else if (num_masters == iatep_prot->num_sequences)
  {
    /* assign proteins to master sequences */
    nuc_seq_num = 0;
    for (prot_seq_num = 0; prot_seq_num < iatep_prot->num_sequences; prot_seq_num++)
    {
      if (iatep_nuc->is_seg != NULL)
      {
        while (iatep_nuc->is_seg [nuc_seq_num] && nuc_seq_num < iatep_nuc->num_sequences)
        {
          nuc_seq_num++;
        }
      }
      if (nuc_seq_num < iatep_nuc->num_sequences)
      {
        NewAssociationList (&assoc_list, nuc_seq_num + 1, NULL);
        nuc_seq_num ++;
      }
      else
      {
        NewAssociationList (&assoc_list, 0, NULL);
      }
    }
  }
  else if (num_masters == 1)
  {
    /* assign all proteins to one sequence */
    for (prot_seq_num = 0; prot_seq_num < iatep_prot->num_sequences; prot_seq_num++)
    {
      NewAssociationList (&assoc_list, 1, NULL);
    }
  }
  else
  {
    /* can't get a match.  Null will be returned. */
  }
  
  iatep_nuc = IDAndTitleEditFree (iatep_nuc);
  iatep_prot = IDAndTitleEditFree (iatep_prot);
  return assoc_list;
}

/* This function will attempt to make a default assignation of nucleotides to proteins.
 * If it is unable to produce a default mapping, it will prompt the user for the mapping.
 * It will then build the nuc-prot sets and discard any proteins for which no nucleotide
 * was assigned.
 */
 
extern NucProtAssocPtr 
AssignProteinsForSequenceSet 
(SeqEntryPtr nuc_list,
 SeqEntryPtr prot_list,
 Boolean     always_review)
{
  NucProtAssocPtr            assoc_list, tmp_list;
  Boolean          all_provided;
  
  if (nuc_list == NULL || prot_list == NULL)
  {
    return NULL;
  }
  
  assoc_list = BuildAssociationListByMatch (nuc_list, prot_list);
  
  if (! AnyLocationsProvided (assoc_list))
  {
    tmp_list = BuildAssociationListByPosition (nuc_list, prot_list);
    if (tmp_list != NULL)
    {
      assoc_list = FreeAssociationList (assoc_list);
      assoc_list = tmp_list;
    }
  }
  
  all_provided = AllLocationsProvided (assoc_list);
  if (!always_review && all_provided) {
    if (!PickCodingRegionLocationsForProteinNucleotidePairs (assoc_list,
                                                             nuc_list,
                                                             prot_list)) {
      always_review = TRUE;
    } else if (FindFeaturesInIdenticalRegions (assoc_list)) {
      always_review = TRUE;
    }
  }  
  
  if (always_review || !AllLocationsProvided (assoc_list))
  {
    tmp_list = CollectNucleotideProteinAssociations (nuc_list, prot_list, assoc_list);
    assoc_list = FreeAssociationList (assoc_list);
    assoc_list = tmp_list;
  }
  return assoc_list;
}

static void BuildNucProtSets 
(SeqEntryPtr      nuc_list,
 SeqEntryPtr      prot_list,
 SequencesFormPtr sqfp,
 Int2             code)
{
  if (nuc_list == NULL || prot_list == NULL ||  sqfp == NULL )
  {
    return;
  }
  ReplaceDuplicateProteinIDs (nuc_list, prot_list);
  AssignProteinsToSelectedNucleotides (sqfp->nuc_prot_assoc_list,
                                       nuc_list,
                                       prot_list,
                                       sqfp, code,
                                       sqfp->makeMRNA);
  sqfp->nuc_prot_assoc_list = FreeAssociationList (sqfp->nuc_prot_assoc_list);
}

typedef enum {
  eSubmitAnnotType_CDS = 1,
  eSubmitAnnotType_rRNA,
  eSubmitAnnotType_Gene,
  eSubmitAnnotType_None 
} ESubmitAnnotTypen;


static SubmissionFeatureInfoPtr SubmissionFeatureInfoNew(void)
{
  SubmissionFeatureInfoPtr info;

  info = (SubmissionFeatureInfoPtr) MemNew (sizeof (SubmissionFeatureInfoData));
  info->strand = Seq_strand_unknown;
  return info;
}


NLM_EXTERN SubmissionFeatureInfoPtr SubmissionFeatureInfoFree (SubmissionFeatureInfoPtr info)
{
  if (info != NULL) {
    info->gene_name = MemFree (info->gene_name);
    info->product = MemFree (info->product);
    info->prot_desc = MemFree (info->prot_desc);
    info->feat_comment = MemFree (info->feat_comment);
    info = MemFree (info);
  }
  return info;
}


static Boolean IsSubmissionFeatureInfoEmpty(SubmissionFeatureInfoPtr info)
{
  Boolean rval = TRUE;

  if (info == NULL) {
    return TRUE;
  }
  switch (info->feature_type) {
    case FEATDEF_CDS:
      if (!StringHasNoText (info->product)
          || !StringHasNoText (info->prot_desc)
          || !StringHasNoText (info->feat_comment)) {
        rval = FALSE;
      }
      break;
    case FEATDEF_rRNA:
      if (!StringHasNoText (info->product)
          || !StringHasNoText (info->feat_comment)) {
        rval = FALSE;
      }
      break;
    case FEATDEF_GENE:
      if (!StringHasNoText (info->gene_name)) {
        rval = FALSE;
      }
      break;
  }
  return rval;
}


static ProtRefPtr ProtRefFromSubmissionFeatureInfo (SubmissionFeatureInfoPtr info)
{
  ProtRefPtr prp;

  if (info == NULL || info->feature_type != FEATDEF_CDS) {
    return NULL;
  }
  prp = ProtRefNew ();
  if (!StringHasNoText (info->product)) {
    ValNodeAddPointer (&(prp->name), 0, StringSave (info->product));
  }
  if (!StringHasNoText (info->prot_desc)) {
    prp->desc = StringSave (info->prot_desc);
  }
  return prp;
}


static GeneRefPtr GeneRefFromSubmissionFeatureInfo (SubmissionFeatureInfoPtr info)
{
  GeneRefPtr grp = NULL;

  if (info == NULL) {
    return NULL;
  }
  if (!StringHasNoText (info->gene_name)) {
    grp = GeneRefNew ();
    grp->locus = StringSave (info->gene_name);
  }
  return grp;
}


static RnaRefPtr RnaRefFromSubmissionFeatureInfo (SubmissionFeatureInfoPtr info)
{
  RnaRefPtr rrp;

  if (info == NULL) {
    return NULL;
  }
  rrp = RnaRefNew ();
  rrp->type = RNA_TYPE_rRNA;
  if (!StringHasNoText (info->product)) {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave (info->product);
  }
  return rrp;
}


static void ApplySubmissionFeatureInfoToLocation (SeqFeatPtr sfp, SubmissionFeatureInfoPtr info)
{
  if (sfp == NULL || info == NULL) {
    return;
  }

  if (info->strand == Seq_strand_minus) 
  {
    /* reverse strand direction - strand direction is plus by default */
    SetSeqLocStrand (sfp->location, Seq_strand_minus);
  }

  SetSeqLocPartial (sfp->location, info->partial5, info->partial3);
    
  sfp->partial = (info->partial5 || info->partial3);
}


NLM_EXTERN SeqFeatPtr AddCDSToSeqEntry (SeqEntryPtr sep, SubmissionFeatureInfoPtr info, BoolPtr ambig)
{
  Int2               genCode;
  CdRegionPtr        crp;
  SeqFeatPtr         sfp;
  Uint2              parenttype;
  Pointer            parentptr;
  BioseqSetPtr       bssp;
  SeqEntryPtr        setsep = NULL;

  GetSeqEntryParent (sep, &parentptr, &parenttype);
  if (parenttype == OBJ_BIOSEQSET
      && (bssp = (BioseqSetPtr) parentptr) != NULL
      && bssp->_class == BioseqseqSet_class_nuc_prot) {
    setsep = SeqMgrGetSeqEntryForData (bssp);
  }

  /*Create a new CDS feature */

  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (NULL == crp) 
  {
    return NULL;
  }

  sfp = CreateNewFeature (sep, setsep, SEQFEAT_CDREGION, NULL); 	

  if (NULL == sfp)
  {
    crp = CdRegionFree (crp);
    return NULL;
  }
  
  sfp->data.value.ptrvalue = (Pointer) crp;

  ApplySubmissionFeatureInfoToLocation (sfp, info);
  
  /* Choose frame for new CDS feature */
  if (!SetBestFrameByLocation (sfp)) {
    if (ambig != NULL) {
      *ambig = TRUE;
    }
  }

  return sfp;
}


static SeqEntryPtr CreateProteinSeqEntryForCDS (SeqFeatPtr cds, ProtRefPtr prp, SeqEntryPtr sep)
{
  ByteStorePtr       bs;
  CharPtr            prot, ptr;
  Char               ch;
  Int4               i;
  BioseqPtr          bsp;
  SeqEntryPtr        old, psep;
  MolInfoPtr         mip;
  Boolean            partial5, partial3;
  ValNodePtr         vnp;
  SeqFeatPtr         prot_sfp;

  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
  if (NULL == bs)
    return NULL;

  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (NULL == prot)
    return NULL;

  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    BSWrite (bs, (VoidPtr) ptr, StringLen (ptr));
  }
  prot = MemFree (prot);

  /* Create the product protein Bioseq */
  
  bsp = BioseqNew ();
  if (NULL == bsp) {
    bs = BSFree (bs);
    return NULL;
  }

  /* Create a new SeqEntry for the Prot Bioseq */
  psep = SeqEntryNew ();
  if (NULL == psep) {
    bs = BSFree (bs);
    return NULL;
  }

  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  old = SeqEntrySetScope (sep);
  bsp->id = MakeNewProteinSeqId (cds->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  SeqEntrySetScope (old);  
  
  psep->choice = 1;
  psep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
  
  /* Add a descriptor to the protein Bioseq */
  
  mip = MolInfoNew ();
  if (mip != NULL) {  
    mip->biomol = 8;
    mip->tech = 8;
    CheckSeqLocForPartial (cds->location, &partial5, &partial3);
    if (partial5 && partial3) {
      mip->completeness = 5;
    } else if (partial5) {
      mip->completeness = 3;
    } else if (partial3) {
      mip->completeness = 4;
    }
    vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
    if (vnp == NULL) {
      mip = MolInfoFree (mip);
    } else {
      vnp->data.ptrvalue = (Pointer) mip;
    }
  }

  if (prp != NULL) {
    prot_sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (prot_sfp != NULL) {
      prot_sfp->data.value.ptrvalue = (Pointer) prp;
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    }
  }
  return psep;
}


static void AddProteinSequenceToNucSeq (SeqEntryPtr sep, SeqEntryPtr psep, SeqFeatPtr cds)
{
  SeqDescrPtr  descr;
  Uint2        parenttype;
  Pointer      parentptr;
  BioseqSetPtr bssp;
  SeqEntryPtr  setsep = NULL;

  GetSeqEntryParent (sep, &parentptr, &parenttype);
  if (parenttype == OBJ_BIOSEQSET
      && (bssp = (BioseqSetPtr) parentptr) != NULL
      && bssp->_class == BioseqseqSet_class_nuc_prot) {
    setsep = SeqMgrGetSeqEntryForData (bssp);
  }

  if (setsep == NULL) {
    descr = ExtractBioSourceAndPubs (sep);
    AddSeqEntryToSeqEntry (sep, psep, TRUE);
    ReplaceBioSourceAndPubs (sep, descr);
  } else {
    AddSeqEntryToSeqEntry (setsep, psep, TRUE);
  }

  SetSeqFeatProduct (cds, psep->data.ptrvalue);
    
  /* after the feature has been created, then adjust it for gaps */
  /* Note - this step may result in multiple coding regions being created. */
  AdjustCDSLocationsForUnknownGapsCallback (cds, NULL);
}


typedef struct submissionfeatureinfoplus {
  SubmissionFeatureInfoPtr info;
  Int4 ambig_count;
} SubmissionFeatureInfoPlusData, PNTR SubmissionFeatureInfoPlusPtr;

static void AddFeaturesFromSubmissionDialogCallback (BioseqPtr bsp, Pointer data)
{
  SubmissionFeatureInfoPlusPtr plus;
  SeqFeatPtr sfp = NULL, gene = NULL;
  Boolean    ambig = FALSE;
  ProtRefPtr prp;
  RnaRefPtr  rrp;
  GeneRefPtr grp;
  SeqEntryPtr psep, sep;

  if (bsp == NULL || ISA_aa (bsp->mol) 
      || (plus = (SubmissionFeatureInfoPlusPtr) data) == NULL
      || plus->info == NULL) {
    return;
  }

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    return;
  }

  switch (plus->info->feature_type) {
    case FEATDEF_CDS:
      sfp = AddCDSToSeqEntry (sep, plus->info, &ambig);
      if (ambig) {
        plus->ambig_count++;
      }
      if (sfp != NULL) {
        prp = ProtRefFromSubmissionFeatureInfo (plus->info);
        psep = CreateProteinSeqEntryForCDS (sfp, prp, sep);
        AddProteinSequenceToNucSeq (sep, psep, sfp);
      }
      break;
    case FEATDEF_rRNA:
      rrp = RnaRefFromSubmissionFeatureInfo (plus->info);
      if (rrp != NULL) {
        sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_RNA, NULL);
        sfp->data.value.ptrvalue = rrp;
        ApplySubmissionFeatureInfoToLocation (sfp, plus->info);
      }        
      break;
    case FEATDEF_GENE:
      grp = GeneRefFromSubmissionFeatureInfo (plus->info);
      if (grp != NULL) {
        sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
        sfp->data.value.ptrvalue = grp;
        ApplySubmissionFeatureInfoToLocation (sfp, plus->info);
      }        
      break;
  }
  if (sfp != NULL) {
    if (plus->info->feature_type != FEATDEF_GENE) {
      grp = GeneRefFromSubmissionFeatureInfo (plus->info);
      if (grp != NULL) {
        gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
        gene->data.value.ptrvalue = grp;
        ApplySubmissionFeatureInfoToLocation (gene, plus->info);
      }        
    }
    if (!StringHasNoText (plus->info->feat_comment)) {
      sfp->comment = StringSave (plus->info->feat_comment);
    }
  }
}


static void AddFeaturesFromSubmissionDialog (SeqEntryPtr sep, SubmissionFeatureInfoPtr info)
{
  SubmissionFeatureInfoPlusData plus;

  if (sep == NULL || info == NULL) {
    return;
  }
  MemSet (&plus, 0, sizeof (SubmissionFeatureInfoPlusData));
  plus.info = info;

  VisitBioseqsInSep (sep, &plus, AddFeaturesFromSubmissionDialogCallback);
  if (plus.ambig_count > 0) {
     Message (MSG_OK, "Possible ambiguous frames detected in %d record%s",
              plus.ambig_count, plus.ambig_count > 1 ? "s": "");
  }
}


typedef struct submissionfeatureinfodialog
{
  DIALOG_MESSAGE_BLOCK

  GrouP           feature_type;
  GrouP           details;
  GrouP           strand;
  ButtoN          partial5;
  ButtoN          partial3;
  TexT            geneName;
  PrompT          protOrRnaPpt;
  TexT            protOrRnaName;
  PrompT          protDescPpt;
  TexT            protDesc;
  TexT            featcomment;
} SubmissionFeatureInfoDialogData, PNTR SubmissionFeatureInfoDialogPtr;


static Pointer SubmissionFeatureInfoDialogToPointer (DialoG d)
{
  SubmissionFeatureInfoDialogPtr dlg;
  SubmissionFeatureInfoPtr info;
  Int2 val;

  dlg = (SubmissionFeatureInfoDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  info = SubmissionFeatureInfoNew();
  val = GetValue (dlg->strand);
  if (val == 2) {
    info->strand = Seq_strand_minus;
  }
  info->partial5 = GetStatus (dlg->partial5);
  info->partial3 = GetStatus (dlg->partial3);
  val = GetValue (dlg->feature_type);
  switch (val) {
    case 1:
      info->feature_type = FEATDEF_CDS;
      info->product = SaveStringFromText (dlg->protOrRnaName);
      info->prot_desc = SaveStringFromText (dlg->protDesc);
      info->gene_name = SaveStringFromText (dlg->geneName);
      info->feat_comment = SaveStringFromText (dlg->featcomment);
      break;
    case 2:
      info->feature_type = FEATDEF_rRNA;
      info->product = SaveStringFromText (dlg->protOrRnaName);
      info->gene_name = SaveStringFromText (dlg->geneName);
      info->feat_comment = SaveStringFromText (dlg->featcomment);
      break;
    case 3:
      info->feature_type = FEATDEF_GENE;
      info->gene_name = SaveStringFromText (dlg->geneName);
      info->feat_comment = SaveStringFromText (dlg->featcomment);
      break;
    default:
      info = SubmissionFeatureInfoFree(info);
      break;
  }
  return info;
}


static void ChangeAnnotType (GrouP g)

{
  SubmissionFeatureInfoDialogPtr dlg;
  Int2              val;

  dlg = (SubmissionFeatureInfoDialogPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  val = GetValue (dlg->feature_type);
  switch (val) {
    case eSubmitAnnotType_Gene :
      SafeHide (dlg->protOrRnaPpt);
      SafeHide (dlg->protOrRnaName);
      SafeHide (dlg->protDescPpt);
      SafeHide (dlg->protDesc);
      SafeShow (dlg->details);
      Select (dlg->geneName);
      break;
    case eSubmitAnnotType_rRNA :
      SafeSetTitle (dlg->protOrRnaPpt, "rRNA Name");
      SafeShow (dlg->protOrRnaPpt);
      SafeShow (dlg->protOrRnaName);
      SafeHide (dlg->protDescPpt);
      SafeHide (dlg->protDesc);
      SafeShow (dlg->details);
      Select (dlg->protOrRnaName);
      break;
    case eSubmitAnnotType_CDS :
      SafeSetTitle (dlg->protOrRnaPpt, "Protein Name");
      SafeShow (dlg->protOrRnaPpt);
      SafeShow (dlg->protOrRnaName);
      SafeShow (dlg->protDescPpt);
      SafeShow (dlg->protDesc);
      SafeShow (dlg->details);
      Select (dlg->protOrRnaName);
      break;
    default :
      SafeHide (dlg->details);
      break;
  }
  Update ();
}


static void ChooseFeatureTypeForSubmissionFeatureInfoDialog(DialoG d, Int2 feature_type)
{
  SubmissionFeatureInfoDialogPtr dlg;

  dlg = (SubmissionFeatureInfoDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  SetValue (dlg->feature_type, feature_type);
  ChangeAnnotType (dlg->feature_type);
}


static DialoG SubmissionFeatureInfoDialog (GrouP h)
{
  SubmissionFeatureInfoDialogPtr dlg;
  GrouP                          p, x, y;
  
  dlg = (SubmissionFeatureInfoDialogPtr) MemNew (sizeof (SubmissionFeatureInfoDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = SubmissionFeatureInfoDialogToPointer;

  dlg->feature_type = HiddenGroup (p, 5, 0, ChangeAnnotType);
  SetObjectExtra (dlg->feature_type, dlg, NULL);
  RadioButton (dlg->feature_type, "CDS");
  RadioButton (dlg->feature_type, "rRNA");
  RadioButton (dlg->feature_type, "Gene");
  RadioButton (dlg->feature_type, "None");
  SetValue (dlg->feature_type, eSubmitAnnotType_CDS);
  dlg->details = HiddenGroup (p, -1, 0, NULL);
  SetGroupSpacing (dlg->details, 10, 10);
  x = HiddenGroup (dlg->details, 2, 0, NULL);
  dlg->partial5 = CheckBox (x, "Incomplete at 5' end", NULL);
  dlg->partial3 = CheckBox (x, "Incomplete at 3' end", NULL);
  dlg->strand = HiddenGroup (dlg->details, 2, 0, NULL);
  RadioButton (dlg->strand, "Plus strand");
  RadioButton (dlg->strand, "Minus strand");
  SetValue (dlg->strand, 1);
  y = HiddenGroup (dlg->details, 2, 0, NULL);
  dlg->protOrRnaPpt = StaticPrompt (y, "Protein Name", 0, dialogTextHeight, programFont, 'l');
  dlg->protOrRnaName = DialogText (y, "", 20, NULL);
  dlg->protDescPpt = StaticPrompt (y, "Protein Description", 0, dialogTextHeight, programFont, 'l');
  dlg->protDesc = DialogText (y, "", 20, NULL);
  StaticPrompt (y, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
  dlg->geneName = DialogText (y, "", 20, NULL);
  StaticPrompt (y, "Comment", 0, 3 * Nlm_stdLineHeight, programFont, 'l');
  dlg->featcomment = ScrollText (y, 20, 3, programFont, TRUE, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) dlg->strand, (HANDLE) y, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_type, (HANDLE) dlg->details, NULL);

  return (DialoG) p;
}


static void GetSubmissionFeatureInfoFeatCommentRect (DialoG d, RectPtr rect)
{
  SubmissionFeatureInfoDialogPtr dlg;

  if (rect == NULL) {
    return;
  }

  dlg = (SubmissionFeatureInfoDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  ObjectRect (dlg->featcomment, rect);
}


static void ApplySubmissionFeatureInfo (SeqEntryPtr sep, DialoG d)
{
  SubmissionFeatureInfoPtr feature_info;

  feature_info = DialogToPointer (d);
  if (feature_info != NULL 
      && feature_info->feature_type == FEATDEF_CDS
      && IsSubmissionFeatureInfoEmpty (feature_info)) {
    ChooseFeatureTypeForSubmissionFeatureInfoDialog (d, eSubmitAnnotType_Gene);
    feature_info = SubmissionFeatureInfoFree (feature_info);
    feature_info = DialogToPointer (d);
  }

  if (feature_info != NULL) {
    AddFeaturesFromSubmissionDialog (sep, feature_info);
    feature_info = SubmissionFeatureInfoFree(feature_info);
  }
}


NLM_EXTERN void SetTsaCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp, sdp_last = NULL;
  MolInfoPtr mip;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  /* note - can't use indexed functions yet */

  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
    sdp = sdp->next;
  }
  if (sdp == NULL) {
    sdp = SeqDescrNew (NULL);
    if (sdp_last == NULL) {
      bsp->descr = sdp;
    } else {
      sdp_last->next = sdp;
    }
    sdp->choice = Seq_descr_molinfo;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) {
    mip = MolInfoNew ();
    sdp->data.ptrvalue = mip;
  }
  mip->tech = MI_TECH_tsa;
  mip->biomol = MOLECULE_TYPE_MRNA;
}


static void ApplySeqPkgToSet (BioseqSetPtr bssp, Int2 seqPackage)
{
  if (bssp == NULL) {
    return;
  }
  switch (seqPackage) {
    case SEQ_PKG_GENOMICCDNA :
      bssp->_class = BioseqseqSet_class_gen_prod_set;
      break;
    case SEQ_PKG_POPULATION :
      bssp->_class = 14;
      break;
    case SEQ_PKG_PHYLOGENETIC :
      bssp->_class = 15;
      break;
    case SEQ_PKG_MUTATION :
      bssp->_class = 13;
      break;
    case SEQ_PKG_ENVIRONMENT :
      bssp->_class = 16;
      break;
    case SEQ_PKG_GENBANK :
      bssp->_class = 7;
      break;
    case SEQ_PKG_TSA :
      bssp->_class = 7;
      break;
    default :
      bssp->_class = 7;
      break;
  }

}


static Pointer FastaSequencesFormToSeqEntryPtr (ForM f)

{
  BioSourcePtr      biop = NULL;
  BioseqSetPtr      bssp;
  Int2              code;
  DatePtr           dp;
  FastaPagePtr      fpp;
  IdListPtr         head;
  SeqEntryPtr       list;
  SeqEntryPtr       next;
  SeqEntryPtr       nucsep;
  Boolean           partialmRNA5;
  Boolean           partialmRNA3;
  SeqEntryPtr       sep;
  SequencesFormPtr  sqfp;
  ValNodePtr        vnp;

  sep = NULL;
  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp == NULL) {
    return (Pointer) sep;
  }

  WatchCursor ();
  Update ();
  head = NULL;
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp != NULL) {
    ResolveCollidingIDs (&head, fpp->list);
  }
  /* NOTE - we do not resolve colliding IDs for proteins here.
   * Duplicate protein IDs are resolved when they are assigned
   * to nucleotide sequences.
   */
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
  if (fpp != NULL) {
    ResolveCollidingIDs (&head, fpp->list);
  }
  FreeTree (&head);
  code = 1;
  list = NULL;
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp != NULL) {
    list = fpp->list;
    fpp->list = NULL;
  }
  if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA
      || PackageTypeIsSet (sqfp->seqPackage))
  {
    bssp = BioseqSetNew ();
    if (bssp != NULL) {
      ApplySeqPkgToSet (bssp, sqfp->seqPackage);
      sep = SeqEntryNew ();
      if (sep != NULL) {
        sep->choice = 2;
        sep->data.ptrvalue = (Pointer) bssp;
      }
    }
  }
  if (list != NULL) {
    if (list->next == NULL) {
      OnlyOneComponentWarning (sqfp);
    }
    while (list != NULL) {
      next = list->next;
      list->next = NULL;
      if (sep != NULL) {
        AddSeqEntryToSeqEntry (sep, list, TRUE);
        AutomaticNucleotideProcess (sqfp, list, list);
      } else {
        sep = list;
        AutomaticNucleotideProcess (sqfp, list, list);
      }
      if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA
          || PackageTypeIsSet (sqfp->seqPackage)) 
      {
        PutMolInfoOnSeqEntry (sqfp, list);
      }
      list = next;
    }
  }
  if (sep != NULL) {
    sqfp->dnamolfrommolinfo = 0;
    dp = DateCurr ();
    if (dp != NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_create_date);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) dp;
      }
    }
  }
  if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA &&
      sqfp->seqPackage != SEQ_PKG_GAPPED) {
    if (! TextHasNoText (sqfp->defline)) {
      ApplyAnnotationToAll (ADD_TITLE, sep, NULL, NULL,
                            NULL, NULL, NULL, NULL, NULL, sqfp->defline);
    }
    if (GetStatus (sqfp->orgPrefix)) {
      PrefixOrgToDefline (sep);
    }
  }

  if (sep != NULL && sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
    list = NULL;
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
    if (fpp != NULL) {
      list = fpp->list;
      /* now we will keep instantiated mrna bioseqs */
      fpp->list = NULL;
    }

    if (list != NULL) {
      nucsep = FindNucSeqEntry (sep);
      if (nucsep != NULL) {
        partialmRNA5 = GetStatus (sqfp->partialmRNA5);
        partialmRNA3 = GetStatus (sqfp->partialmRNA3);
        while (list != NULL) {
          next = list->next;
          list->next = NULL;
          AddSeqEntryToSeqEntry (sep, list, TRUE);
          AutomaticMrnaProcess (nucsep, list, partialmRNA5, partialmRNA3);
          list = next;
        }
      }
    }
  }

  list = NULL;
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
  if (fpp != NULL) {
    list = fpp->list;
    fpp->list = NULL;
  }
  if (list != NULL) {
    BuildNucProtSets (sep, list, sqfp, code);
  }
  if (biop != NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_source);
    if (vnp != NULL) {
      vnp->data.ptrvalue = (Pointer) biop;
    }
  }
  if (PackageTypeIsSet (sqfp->seqPackage)) {
    if (GetStatus (sqfp->makeAlign)) {
      AddSeqAlignForSeqEntry (sep, ObjMgrGetEntityIDForChoice (sep), FALSE, TRUE);
    }
  }
  if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA &&
      sqfp->seqPackage != SEQ_PKG_GAPPED) {
      ApplySubmissionFeatureInfo (sep, sqfp->feature_info);
  }
  if (sqfp->seqPackage == SEQ_PKG_TSA) {
    VisitBioseqsInSep (sep, NULL, SetTsaCallback);
  }
  ArrowCursor ();
  Update ();

  FuseNucProtBiosources (sep);
  return (Pointer) sep;
}

static void LaunchSequinQuickGuide (void)

{
  Char       str [256];
#ifdef WIN_MOTIF
  NS_Window  window = NULL;
#endif

  sprintf (str,
           "http://www.ncbi.nlm.nih.gov/Sequin/QuickGuide/sequin.htm#before");
#ifdef WIN_MAC
  Nlm_SendURLAppleEvent (str, "MOSS", NULL);
#endif
#ifdef WIN_MSWIN
  Nlm_MSWin_OpenDocument (str);
#endif
#ifdef WIN_MOTIF
  NS_OpenURL (&window, str, NULL, TRUE);
  NS_WindowFree (window);
#endif
}

extern Boolean allowUnableToProcessMessage;

static CharPtr noOrgInTitleAbort =
"sequences have organism information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Sequin will not continue processing this submission. " \
"Please read the Sequin Quick Guide section on preparing the data files before proceeding. " \
"Do you wish to launch your browser on the Sequin Quick Guide automatically?";

static CharPtr pleaseReadLocalGuide =
"Please read your local copy of the Sequin Quick Guide before annotating your data file.";

static CharPtr noSrcInTitleAbort =
"sequences have source information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Sequin will continue processing this submission. " \
"However, please consider reading the Sequin Quick Guide section on preparing the data files before proceeding.";

static Pointer PhylipSequencesFormToSeqEntryPtr (ForM f)

{
  MsgAnswer         ans;
  BioseqSetPtr      bssp;
  Int2              code;
  DatePtr           dp;
  PhylipPagePtr     ppp;
  SeqEntryPtr       sep;
  Int2              seqtitles;
  Int2              seqtotals;
  Char              str [256];
  SequencesFormPtr  sqfp;
  SeqEntryPtr       tmp;
  CharPtr           ttl;
  ValNodePtr        vnp;
  SeqEntryPtr       list;
  FastaPagePtr      fpp;

  sep = NULL;
  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    code = 1;
    ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (ppp != NULL) {
      sep = ppp->sep;
      ppp->sep = NULL;
    }
    if (sep != NULL) {

      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && (bssp->_class == 7 ||
                             (IsPopPhyEtcSet (bssp->_class)))) {
          seqtitles = 0;
          seqtotals = 0;
          for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
            /*
            ttl = SeqEntryGetTitle (tmp);
            */
            ttl = NULL;
            SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
            if (ttl != NULL) {
              if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
                if (StringISearch (ttl, "[org=") != NULL ||
                    StringISearch (ttl, "[organism=") != NULL) {
                  seqtitles++;
                }
              } else if (StringISearch (ttl, "[") != NULL) {
                seqtitles++;
              }
            }
            seqtotals++;
          }
          if (seqtotals != seqtitles) {
            sprintf (str, "None");
            if (seqtitles > 0) {
              sprintf (str, "Only %d", (int) seqtitles);
            }
            ArrowCursor ();
            Update ();
            Beep ();
            if (! indexerVersion) {
              if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
                ans = Message (MSG_YN, "%s of %d %s", str, (int) seqtotals, noOrgInTitleAbort);
                if (ans == ANS_YES) {
                  LaunchSequinQuickGuide ();
                } else {
                  Message (MSG_OK, "%s", pleaseReadLocalGuide);
                }
                allowUnableToProcessMessage = FALSE;
                QuitProgram ();
                return NULL; /* aborting */
              } else {
                Message (MSG_OK, "%s of %d %s", str, (int) seqtotals, noSrcInTitleAbort);
              }
            } else {
              if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
                Message (MSG_OK, "%s of %d %s (Regular version will abort here.)", str, (int) seqtotals, noOrgInTitleAbort);
              } else {
                Message (MSG_OK, "%s of %d %s (Regular version will continue here.)", str, (int) seqtotals, noSrcInTitleAbort);
              }
            }
          }
        }
        if (bssp != NULL) {
          ApplySeqPkgToSet (bssp, sqfp->seqPackage);
          tmp = bssp->seq_set;
          if (tmp == NULL || tmp->next == NULL) {
            OnlyOneComponentWarning (sqfp);
          }
          for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
            AutomaticNucleotideProcess (sqfp, tmp, tmp);
            PutMolInfoOnSeqEntry (sqfp, tmp);
          }
        }
      } else {
        OnlyOneComponentWarning (sqfp);
        PutMolInfoOnSeqEntry (sqfp, sep);
      }
      dp = DateCurr ();
      if (dp != NULL) {
        vnp = CreateNewDescriptor (sep, Seq_descr_create_date);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (Pointer) dp;
        }
      }
    }
    if (PackageTypeIsSet (sqfp->seqPackage)) {
      if (! TextHasNoText (sqfp->defline)) {
        ApplyAnnotationToAll (ADD_TITLE, sep, NULL, NULL,
                              NULL, NULL, NULL, NULL, NULL, sqfp->defline);
      }
      if (GetStatus (sqfp->orgPrefix)) {
        PrefixOrgToDefline (sep);
      }
    }
    ApplySubmissionFeatureInfo (sep, sqfp->feature_info);
  }
  FuseNucProtBiosources (sep);
  
  list = NULL;
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
  if (fpp != NULL) {
    list = fpp->list;
    fpp->list = NULL;
  }
  if (list != NULL) {
    BuildNucProtSets (sep, list, sqfp, code);
  }

  return (Pointer) sep;
}

static void SeqEntryPtrToSourceTab (SequencesFormPtr sqfp);
static SeqEntryPtr GetSeqEntryFromSequencesForm (SequencesFormPtr sqfp);
static void ReplaceAllAliases (SeqEntryPtr sep);
static void ReplaceMolNamesWithMolBracketsInDefinitionLines (SeqEntryPtr sep);
static Boolean CheckSequencesForOrganisms (SequencesFormPtr sqfp);
static void SequencesFormDeleteProc (Pointer formDataPtr);

static void NucleotideImportFinish (SequencesFormPtr sqfp)
{
  FastaPagePtr  fpp = NULL;
  PhylipPagePtr ppp = NULL;
  SeqEntryPtr   sep = NULL;
  Boolean       cancelled = FALSE;
  BioseqSetPtr  bssp;
  
  if (sqfp == NULL) return;
  
  if (sqfp->seqFormat == SEQ_FMT_FASTA) {
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) 
    {
      sep = fpp->list;
    }
  } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
    ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (ppp != NULL) {
      sep = ppp->sep;
      if (sep != NULL && IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL)
      {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        sep = bssp->seq_set;
      }
    }
  }

  if (sep == NULL)
  {
    Disable (sqfp->molecule_btn);
    Disable (sqfp->topology_btn);
    Disable (sqfp->import_mod_btn);
    Disable (sqfp->source_assist_btn);
    Disable (sqfp->specify_orgs_btn);
    Disable (sqfp->specify_locs_btn);
    Disable (sqfp->specify_gcode_btn);
    Disable (sqfp->specify_mgcode_btn);
    Disable (sqfp->clear_mods_btn);
  }
  else
  {
    Enable (sqfp->molecule_btn);
    Enable (sqfp->topology_btn);
    Enable (sqfp->import_mod_btn);
    Enable (sqfp->source_assist_btn);
    Enable (sqfp->specify_orgs_btn);
    Enable (sqfp->specify_locs_btn);
    Enable (sqfp->specify_gcode_btn);
    Enable (sqfp->specify_mgcode_btn);
    Enable (sqfp->clear_mods_btn);
  }

  ReplaceAllAliases (sep);
  ReplaceMolNamesWithMolBracketsInDefinitionLines (sep);
  
  AddDefaultModifierValues (sep);
  
  if (fpp != NULL)
  {
    Reset (fpp->doc);
    FormatFastaDoc (fpp);
  }

  if (cancelled)
  {
    SequencesFormDeleteProc (sqfp);    
  }

  SeqEntryPtrToSourceTab (sqfp);
}

static CharPtr GetFirstModValueFromSeqEntryTitles (SeqEntryPtr sep, CharPtr mod_name)
{
  CharPtr value = NULL;
  
  while (sep != NULL && value == NULL)
  {
    value = GetModValueFromSeqEntry (sep, mod_name);
    sep = sep->next;
  }
  return value;
}

static SeqIdPtr GetSeqIdFromSeqEntryPtr (SeqEntryPtr sep) 
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL) {
  	return NULL;
  }
  if (IS_Bioseq (sep)) {
  	bsp = (BioseqPtr) sep->data.ptrvalue;
    return bsp->id;
  } else if (IS_Bioseq_set(sep)) {
  	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    return GetSeqIdFromSeqEntryPtr (bssp->seq_set);
  }
  return NULL;
}

static void 
ReportMissingOrganismNames 
(ValNodePtr no_org_list)
{
  Char         path [PATH_MAX];
  FILE         *fp;
  ValNodePtr   vnp;

  if (no_org_list == NULL)
  {
    return;
  }
  
  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp == NULL) return;

  if (no_org_list != NULL)
  {
    fprintf (fp, "The following sequences have no organism names.  You must supply one for each sequence listed.\n");
    for (vnp = no_org_list; vnp != NULL; vnp = vnp->next)
    {
      fprintf (fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }
  }

  FileClose (fp);
  LaunchGeneralTextViewer (path, "Organism Errors");
  FileRemove (path);
}


static Boolean CheckSequencesForOrganisms (SequencesFormPtr sqfp)
{
  SeqEntryPtr       sep_list;
  ValNodePtr        no_org_list = NULL;
  Boolean           rval = TRUE;
  IDAndTitleEditPtr iatep;
  Int4              seq_num;
  CharPtr           org_name_from_title;
  
  if (sqfp == NULL) return FALSE;
  sep_list = GetSeqEntryFromSequencesForm (sqfp);
  if (sep_list == NULL) return FALSE;
  
  iatep = SeqEntryListToIDAndTitleEdit (sep_list);
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    org_name_from_title = FindValueFromPairInDefline ("organism", iatep->title_list [seq_num]);
    if (StringHasNoText (org_name_from_title))
 	  {
 	    rval = FALSE;
 	    ValNodeAddPointer (&no_org_list, 0, StringSave (iatep->id_list [seq_num]));
 	  }
    org_name_from_title = MemFree (org_name_from_title);
  }
  
  ReportMissingOrganismNames (no_org_list);
  no_org_list = ValNodeFreeData (no_org_list);
  iatep = IDAndTitleEditFree (iatep);
  return rval;
}

static Boolean ExportSequencesForm (ForM f, CharPtr filename)

{
  SequencesFormPtr  sqfp;
  Boolean           rval = FALSE;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case NUCLEOTIDE_PAGE :
        rval = ExportDialog (sqfp->dnaseq, "");
        break;
      case MRNA_PAGE :
        break;
      case PROTEIN_PAGE :
        break;
      case ORGANISM_PAGE:
        break;
      default :
        break;
    }
  }
  return rval;
}

NLM_EXTERN Boolean IsAnnotTabEmpty (SequencesFormPtr sqfp)

{
  SubmissionFeatureInfoPtr feature_info;
  Boolean rval;
  
  if (sqfp == NULL)
  {
    return TRUE;
  }
  feature_info = (SubmissionFeatureInfoPtr) DialogToPointer (sqfp->feature_info);
  rval = IsSubmissionFeatureInfoEmpty (feature_info);
  feature_info = SubmissionFeatureInfoFree (feature_info);
  return rval;  
}


NLM_EXTERN SubmissionFeatureInfoPtr GetSubmissionFeatureInfo (SequencesFormPtr sqfp)
{
  SubmissionFeatureInfoPtr feature_info;
  
  if (sqfp == NULL)
  {
    return NULL;
  }
  feature_info = (SubmissionFeatureInfoPtr) DialogToPointer (sqfp->feature_info);
  return feature_info;  
}


static void ChooseFeatureTypeForSubmissionFeatureInfoDialog(DialoG d, Int2 feature_type);

static Boolean ImportSequencesForm (ForM f, CharPtr filename)

{
  SequencesFormPtr  sqfp;
  FastaPagePtr      fpp;
  Boolean           rval = FALSE;
  SeqEntryPtr       seq_list;
  IDAndTitleEditPtr iatep;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case NUCLEOTIDE_PAGE :
        rval = ImportDialog (sqfp->dnaseq, "");
        if (rval)
        {
          NucleotideImportFinish (sqfp);
          seq_list = GetSeqEntryFromSequencesForm (sqfp);
          if (sqfp->seqFormat == SEQ_FMT_FASTA && (fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq)) != NULL
              && SuggestJumpingToWizard (seq_list)) {
            fpp->list = NULL;
            Remove (sqfp->form);
          }
        }
        break;
      case MRNA_PAGE :
        rval = ImportDialog (sqfp->mrnaseq, "");
        break;
      case PROTEIN_PAGE :
        rval = ImportDialog (sqfp->protseq, "");
        if (rval && IsAnnotTabEmpty (sqfp))
        {
          ChooseFeatureTypeForSubmissionFeatureInfoDialog (sqfp->feature_info, eSubmitAnnotType_None);
        }
        break;
      case ORGANISM_PAGE:
        seq_list = GetSeqEntryFromSequencesForm (sqfp);
        iatep = SeqEntryListToIDAndTitleEdit (seq_list);
        rval = ImportModifiersToIDAndTitleEdit (iatep);
        if (rval)
        {
          ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
          SeqEntryPtrToSourceTab (sqfp);    
        }
        iatep = IDAndTitleEditFree (iatep);
        break;
      default :
        break;
    }
  }
  return rval;
}

static void ImportBtnProc (ButtoN b)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp) {
    ImportSequencesForm (sqfp->form, "");
  }
}

static void SetOrgNucProtImportExportItems (SequencesFormPtr sqfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (sqfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) sqfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) sqfp, VIB_MSG_EXPORT);
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case ORGANISM_PAGE :
        SafeSetTitle (importItm, "Import Organism Modifiers From File");
        SafeSetTitle (exportItm, "Export...");
        SafeEnable (importItm);
        SafeDisable (exportItm);
        break;
      case NUCLEOTIDE_PAGE :
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (exportItm);
        switch (sqfp->seqFormat) {
          case SEQ_FMT_FASTA :
            if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
              SafeSetTitle (importItm, "Import Genomic FASTA...");
            } else {
              SafeSetTitle (importItm, "Import Nucleotide FASTA...");
              SafeSetTitle (exportItm, "Export Nucleotide FASTA...");
              SafeEnable (exportItm);
            }
            break;
          case SEQ_FMT_ALIGNMENT :
            SafeSetTitle (importItm, "Import Nucleotide Alignment...");
            break;
          default :
            SafeSetTitle (importItm, "Import Nucleotide FASTA...");
            break;
        }
        SafeEnable (importItm);
        break;
      case MRNA_PAGE :
        SafeSetTitle (importItm, "Import Transcript FASTA...");
        SafeSetTitle (exportItm, "Export...");
        SafeEnable (importItm);
        SafeDisable (exportItm);
        break;
      case PROTEIN_PAGE :
        SafeSetTitle (importItm, "Import Protein FASTA...");
        SafeSetTitle (exportItm, "Export...");
        SafeEnable (importItm);
        SafeDisable (exportItm);
        break;
      case ANNOTATE_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeSequencesPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) data;
  if (sqfp == NULL) {
    return;
  }

  if (sqfp->show_annot) {
    newval = newval + 2;
    oldval = oldval + 2;
  }

  sqfp->currentPage = newval;
  SafeHide (sqfp->pages [oldval]);
  Update ();
  switch (sqfp->tagFromPage [newval]) {
    case ORGANISM_PAGE :
      break;
    case NUCLEOTIDE_PAGE :
      SendMessageToDialog (sqfp->dnaseq, VIB_MSG_ENTER);
      break;
    case MRNA_PAGE :
      SendMessageToDialog (sqfp->mrnaseq, VIB_MSG_ENTER);
      break;
    case PROTEIN_PAGE :
      SendMessageToDialog (sqfp->protseq, VIB_MSG_ENTER);
      break;
    default :
      break;
  }
  if (newval == 0) {
    SafeSetTitle (sqfp->prevBtn, "<< Prev Form");
  } else {
    SafeSetTitle (sqfp->prevBtn, "<< Prev Page");
  }
  if (newval == sqfp->numPages - 1) {
    SafeSetTitle (sqfp->nextBtn, "Next Form >>");
  } else {
    SafeSetTitle (sqfp->nextBtn, "Next Page >>");
  }
  SetOrgNucProtImportExportItems (sqfp);
  SafeShow (sqfp->pages [newval]);
  Update ();
  switch (sqfp->tagFromPage [newval]) {
    case ORGANISM_PAGE :
      SendHelpScrollMessage (helpForm, "Organism Page", NULL);
      break;
    case NUCLEOTIDE_PAGE :
      if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
        SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Genomic Page");
      } else {
        if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT)
        {
          SendHelpScrollMessage (helpForm, "Nucleotide Page", "Nucleotide Page for Aligned Data Formats");
        }
        else
        {
          SendHelpScrollMessage (helpForm, "Nucleotide Page", "Nucleotide Page for FASTA Data Format");
        }
      }
      break;
    case MRNA_PAGE :
      SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Transcripts Page");
      break;
    case PROTEIN_PAGE :
      SendHelpScrollMessage (helpForm, "Protein Page", NULL);
      break;
    case ANNOTATE_PAGE :
      SendHelpScrollMessage (helpForm, "Annotation Page", NULL);
      break;
    default :
      break;
  }
}

static Boolean SequenceAssistantValidate (SeqEntryPtr seq_list);
static Boolean FinalSequenceValidation (SequencesFormPtr  sqfp)
{
  FastaPagePtr fpp;
  
  if (sqfp == NULL)
  {
    return FALSE;
  }
  if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT || sqfp->seqPackage == SEQ_PKG_GAPPED)
  {
    /* we can't edit these, so there's no sense pestering the user...*/
    return TRUE;
  }
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp == NULL || fpp->list == NULL)
  {
    return FALSE;
  }
  
  return SequenceAssistantValidate (fpp->list);
}


static void NextSequencesFormBtn (ButtoN b)

{
  SequencesFormPtr  sqfp;
  Boolean force_change = FALSE;
  SequencingMethodInfoPtr info;
  ValNodePtr seq_list;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    if (sqfp->currentPage == 1) {
      /* validate before leaving sequence area */
      seq_list = GetSeqEntryFromSequencesForm (sqfp);
      if (seq_list == NULL) {
        Message (MSG_ERROR, "You must import sequences before continuing.");
        return;
      }
      info = DialogToPointer (sqfp->sequencing_method_dlg);
      if (!IsSequencingMethodInfoValid(info, CountSequencesAndSegments (seq_list, TRUE), TRUE)) {
        if (info->quit_now) {          
          QuitFromWizard (sqfp->form);
        }
        info = SequencingMethodInfoFree (info);
        return;
      }
      info = SequencingMethodInfoFree (info);

      /* now switch to other tab */
      Hide (sqfp->seq_tbs);
      Show (sqfp->annot_tbs);
      sqfp->show_annot = TRUE;
      force_change = TRUE;
    }

    if (sqfp->currentPage < 1) {
      SetValue (sqfp->seq_tbs, sqfp->currentPage + 1);
    } else if (sqfp->currentPage < sqfp->numPages - 1) {
      SetValue (sqfp->annot_tbs, sqfp->currentPage - 1);
      if (force_change) {
        ChangeSequencesPage (sqfp, 0, -1);
      }
    } else if (sqfp->goToNext != NULL) {
      if (!CheckSequencesForOrganisms (sqfp) || !FinalSequenceValidation (sqfp)) return;
      (sqfp->goToNext) (b);
    }
  }
}

static void PrevSequencesFormBtn (ButtoN b)

{
  SequencesFormPtr  sqfp;
  Boolean force_change = FALSE;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    if (sqfp->currentPage == 2) {
      /* switch between tabs */
      Hide (sqfp->annot_tbs);
      Show (sqfp->seq_tbs);
      sqfp->show_annot = FALSE;
      force_change = TRUE;
    }
    if (sqfp->currentPage > 2) {
      SetValue (sqfp->annot_tbs, sqfp->currentPage - 3);
    } else if (sqfp->currentPage > 0) {
      SetValue (sqfp->seq_tbs, sqfp->currentPage - 1);
      if (force_change) {
        ChangeSequencesPage (sqfp, 1, 2);
      }
    } else if (sqfp->goToPrev != NULL) {
      (sqfp->goToPrev) (b);
    }
  }
}

static void SetModifierList (DoC doc, ValNodePtr mod_list);
static void SeqEntryPtrToOrgDoc (SequencesFormPtr sqfp);

static void ClearOrganismModifiers (SequencesFormPtr sqfp)
{
  if (sqfp == NULL) return;
  Disable (sqfp->import_mod_btn);
  Disable (sqfp->source_assist_btn);
  Disable (sqfp->specify_orgs_btn);
  Disable (sqfp->specify_locs_btn);
  Disable (sqfp->specify_gcode_btn);
  Disable (sqfp->specify_mgcode_btn);
  Disable (sqfp->clear_mods_btn);

  SeqEntryPtrToOrgDoc (sqfp);
}


static void DeleteAllSequencesFromForm (SequencesFormPtr sqfp)
{
  FastaPagePtr fpp;
  PhylipPagePtr ppp;

  if (sqfp == NULL) {
    return;
  }
  if (sqfp->seqFormat == SEQ_FMT_FASTA) {
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      ResetFastaPage (fpp);
      fpp->path [0] = '\0';
      SafeHide (fpp->have_seq_instr_grp);
      Reset (fpp->doc);
      SafeShow (fpp->instructions);
      Update ();
      if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA)
      {
        SetTitle (fpp->import_btn, "Import Nucleotide FASTA");
      }
      Enable (fpp->import_btn);
      Disable (fpp->clear_btn);
      Disable (sqfp->molecule_btn);
      Disable (sqfp->topology_btn);
    }
  } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
    ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (ppp != NULL) {
      ResetPhylipPage (ppp);
      ppp->path [0] = '\0';
      SetPhylipDocInstructions (ppp);
    }
  }
  ClearOrganismModifiers (sqfp);
}


static void SequencesFormDeleteProc (Pointer formDataPtr)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) formDataPtr;
  if (sqfp != NULL) {
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case ORGANISM_PAGE :
        ClearText (CurrentVisibleText ());
        break;
      case NUCLEOTIDE_PAGE :
        DeleteAllSequencesFromForm (sqfp);
        break;
      case MRNA_PAGE :
        if (sqfp->seqFormat == SEQ_FMT_FASTA) {
          fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
          if (fpp != NULL) {
            ResetFastaPage (fpp);
            fpp->path [0] = '\0';
            SafeHide (fpp->have_seq_instr_grp);
            Reset (fpp->doc);
            SafeShow (fpp->instructions);
            Update ();
          }
        }
        break;
      case PROTEIN_PAGE :
        if (ANS_YES == Message (MSG_YN, "Are you sure you want to remove all of the protein sequences?"))
        {
          if (sqfp->seqFormat == SEQ_FMT_FASTA) {
            fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
            if (fpp != NULL) {
              ResetFastaPage (fpp);
              fpp->path [0] = '\0';
              SafeHide (fpp->have_seq_instr_grp);
              Reset (fpp->doc);
              Disable (fpp->clear_btn);
              SafeShow (fpp->instructions);
              Update ();
            }
          } else {
            ClearText (CurrentVisibleText ());
          }
        }
        break;
      default :
        break;
    }
  }
}

static CharPtr  seqSegSeqFormTabs [] = {
  "Nucleotide", "Sequencing Method", NULL
};

static CharPtr  cdnaGenSeqFormTabs [] = {
  "Genomic",  "Sequencing Method", NULL
};

static CharPtr  popPhyMutSeqFormTabs [] = {
  "Nucleotide",  "Sequencing Method", NULL
};

static CharPtr  seqSegAnnotFormTabs [] = {
  "Organism", "Proteins", NULL
};

static CharPtr  cdnaGenAnnotFormTabs [] = {
  "Organism", "Transcripts", "Proteins", NULL
};

static CharPtr  popPhyMutAnnotFormTabs [] = {
  "Organism", "Proteins", "Annotation", NULL
};

static void PasteIntoDialog (DialoG seq)

{
  Char     ch;
  FILE     *fp;
  Char     path [PATH_MAX];
  CharPtr  ptr;
  CharPtr  str;

  if (Nlm_ClipboardHasString ()) {
    TmpNam (path);
    fp = FileOpen (path, "w");
    if (fp == NULL) return;
    str = ClipboardToString ();
    if (str != NULL) {
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == '\r') {
          *ptr = '\n';
        }
        ptr++;
        ch = *ptr;
      }
      FilePuts (str, fp);
      MemFree (str);
    }
    FileClose (fp);
    ImportDialog (seq, path);
    FileRemove (path);
  }
}

static void SequencesFormMessage (ForM f, Int2 mssg)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    switch (mssg) {
      case VIB_MSG_IMPORT :
        ImportSequencesForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportSequencesForm (f, NULL);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        switch (sqfp->tagFromPage [sqfp->currentPage]) {
          case ORGANISM_PAGE :
            StdPasteTextProc (NULL);
            break;
          case NUCLEOTIDE_PAGE :
            PasteIntoDialog (sqfp->dnaseq);
            break;
          case MRNA_PAGE :
            PasteIntoDialog (sqfp->mrnaseq);
            break;
          case PROTEIN_PAGE :
            PasteIntoDialog (sqfp->protseq);
            break;
          default :
            StdPasteTextProc (NULL);
            break;
        }
        break;
      case VIB_MSG_DELETE :
        SequencesFormDeleteProc (sqfp);
        break;
      default :
        if (sqfp->appmessage != NULL) {
          sqfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void InitOrgNucProtFormActivate (WindoW w)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (w);
  if (sqfp != NULL) {
    if (sqfp->activate != NULL) {
      sqfp->activate (w);
    }
    SetOrgNucProtImportExportItems (sqfp);
  }
}

static void ChangeMrnaFlag (ButtoN b)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    sqfp->makeMRNA = GetStatus (b);
    if (sqfp->makeMRNA) {
      SetAppParam ("SEQUINCUSTOM", "PREFERENCES", "CREATEMRNA", "TRUE");
    } else {
      SetAppParam ("SEQUINCUSTOM", "PREFERENCES", "CREATEMRNA", "FALSE");
    }
  }
}

extern CharPtr CreateListMessage (CharPtr msg_before, CharPtr msg_after, ValNodePtr id_list)
{
  Int4       num_pos;
  ValNodePtr vnp;
  CharPtr    msg_txt;
  Int4       txt_len = StringLen (msg_before) + StringLen (msg_after) + 3;
  Char       num_buf [14];
  
  if (id_list == NULL) return NULL;
  num_pos = ValNodeLen (id_list);
  for (vnp = id_list; vnp != NULL; vnp = vnp->next)
  {
    if (StringHasNoText (vnp->data.ptrvalue))
    {
      txt_len += PRINTED_INT_MAX_LEN;
    }
    else
    {
      txt_len += StringLen (vnp->data.ptrvalue) + 5;
    }
  }
  msg_txt = (CharPtr) MemNew (txt_len * sizeof (Char));
  if (msg_txt != NULL)
  {
    msg_txt [0] = 0;
    if (msg_before != NULL)
    {
      StringCat (msg_txt, msg_before);
      if (num_pos > 1)
      {
        StringCat (msg_txt, "s ");
      }
      else
      {
        StringCat (msg_txt, " ");
      }
    }
    
    for (vnp = id_list; vnp != NULL; vnp = vnp->next)
    {
      if (num_pos > 1 && vnp->next == NULL)
      {
        StringCat (msg_txt, "and ");
      }
      if (StringHasNoText (vnp->data.ptrvalue))
      {
        sprintf (num_buf, "%d", vnp->choice);
        StringCat (msg_txt, num_buf);
      }
      else
      {
        StringCat (msg_txt, vnp->data.ptrvalue);
      }
      if (vnp->next != NULL)
      {
  	    if (num_pos > 2)
  	    {
	        StringCat (msg_txt, ", ");
  	    }
  	    else
  	    {
  	  	  StringCat (msg_txt, " ");
  	    }
      }
    }
    StringCat (msg_txt, msg_after);
  }
  return msg_txt;
}


static Boolean ContinueWithErrorList (ValNodePtr err_list, Boolean ask_for_continue)
{
  ValNodePtr           vnp;
  GrouP                required_grp = NULL;
  GrouP                warning_grp = NULL;
  GrouP                h, g, c;
  PrompT               p;
  ButtoN               b;
  WindoW                w;
  ModalAcceptCancelData acd;
  Boolean               ok_to_continue = TRUE;

  
  if (err_list == NULL) return TRUE;
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;

  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  g = HiddenGroup (h, 1, 0, NULL);
  /* create required list */
  for (vnp = err_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == CREATE_FASTA_REQUIRED)
    {
      ok_to_continue = FALSE;
      if (required_grp == NULL)
      {
        required_grp = NormalGroup (g, 1, 0, "Required", systemFont, NULL);
      }
      MultiLinePrompt (required_grp, vnp->data.ptrvalue, 600, systemFont);
    }
  }
  /* create warning list */
  for (vnp = err_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == CREATE_FASTA_WARNING)
    {
      if (warning_grp == NULL)
      {
        warning_grp = NormalGroup (g, 1, 0, "Warning", systemFont, NULL);
      }
      MultiLinePrompt (warning_grp, vnp->data.ptrvalue, 600, systemFont);
    }
  }
  
  if (! ask_for_continue)
  {
    p = NULL;
    c = HiddenGroup (h, 1, 0, NULL);
    b = PushButton (c, "OK", ModalCancelButton);
    SetObjectExtra (b, &acd, NULL);
  }
  else if (ok_to_continue)
  {
    p = StaticPrompt (h, "Continue anyway?",
                      0, dialogTextHeight, systemFont, 'c');
    c = HiddenGroup (h, 2, 0, NULL);
    b = PushButton (c, "Yes", ModalAcceptButton);
    SetObjectExtra (b, &acd, NULL);
    b = PushButton (c, "No", ModalCancelButton);
    SetObjectExtra (b, &acd, NULL);
  }
  else
  {
    p = StaticPrompt (h, "You must resolve the required errors before continuing.",
                      0, dialogTextHeight, systemFont, 'c');
    c = HiddenGroup (h, 1, 0, NULL);
    b = PushButton (c, "OK", ModalCancelButton);
    SetObjectExtra (b, &acd, NULL);
  }
  
  if (ask_for_continue)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) p, (HANDLE) c, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  }
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  if (acd.accepted)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
  
}

static CharPtr valid_iupac_characters = "atgcbdhkmnrsuvwy";

static Boolean SeqCharsOk
(CharPtr seq_chars,
 Int4 seq_num,
 CharPtr local_id,
 ValNodePtr PNTR err_list)
{
  CharPtr cp;
  Char    ch;
  Boolean at_least_one = FALSE;
  Int4    len = StringLen (seq_chars);
  CharPtr badchars;
  CharPtr err_msg;
  CharPtr empty_fmt_d = "There are no sequence characters for sequence %d.  Please enter some.";
  CharPtr empty_fmt_s = "There are no sequence characters for sequence %s.  Please enter some.";
  CharPtr bad_char_fmt_d = "There were %d illegal characters were found in sequence %d: %s."
  	         "  Repeated characters are listed only once. "
  	         "  You may only have IUPAC characters in your sequence ";
  CharPtr bad_char_fmt_s = "There were %d illegal characters were found in sequence %s: %s."
  	         "  Repeated characters are listed only once. "
  	         "  You may only have IUPAC characters in your sequence ";
  
  if (StringHasNoText (seq_chars))
  {
    err_msg = (CharPtr) MemNew (sizeof (Char) * 
                (StringLen (empty_fmt_d) + PRINTED_INT_MAX_LEN + StringLen (local_id)));
    if (err_msg != NULL)
    {
      if (StringHasNoText (local_id))
      {
        sprintf (err_msg, empty_fmt_d, seq_num);
      }
      else
      {
        sprintf (err_msg, empty_fmt_s, local_id);
      }
      ValNodeAddPointer (err_list, CREATE_FASTA_REQUIRED, err_msg);
    }
  	return FALSE;
  }
  
  badchars = (CharPtr) MemNew (sizeof (Char) * (len + 1));
  if (badchars == NULL) return FALSE;
  badchars[0] = 0;
  len = 0;
  for (cp = seq_chars; *cp != 0; cp++)
  {
    ch = TO_LOWER (*cp);
  	if (isspace ((Int4)ch))
  	{
  	  /* space allowed */
  	}
  	else if (StringChr (valid_iupac_characters, ch) == NULL)
  	{
  	  if (StringChr (badchars, *cp) == NULL)
  	  {
  	  	badchars [len] = ch;
  	  	len++; 
  	  	badchars [len] = 0;
  	  }
  	}
  	else 
  	{
  	  at_least_one = TRUE;
  	}
  }
  if (len > 0)
  {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_char_fmt_d)
               + (2 * PRINTED_INT_MAX_LEN) + StringLen (local_id) + len
               + StringLen (valid_iupac_characters)));
    if (err_msg != NULL)
    {
      if (StringHasNoText (local_id))
      {
        sprintf (err_msg, bad_char_fmt_d, len, seq_num, badchars, valid_iupac_characters);
      }
      else
      {
        sprintf (err_msg, bad_char_fmt_s, len, local_id, badchars, valid_iupac_characters);
      }
      ValNodeAddPointer (err_list, CREATE_FASTA_REQUIRED, err_msg);
    }
  	return FALSE;
  }
  if (!at_least_one)
  {
    err_msg = (CharPtr) MemNew (sizeof (Char) * 
                (StringLen (empty_fmt_d) + PRINTED_INT_MAX_LEN + StringLen (local_id)));
    if (err_msg != NULL)
    {
      if (StringHasNoText (local_id))
      {
        sprintf (err_msg, empty_fmt_d, seq_num);
      }
      else
      {
        sprintf (err_msg, empty_fmt_s, local_id);
      }
      ValNodeAddPointer (err_list, CREATE_FASTA_REQUIRED, err_msg);
    }
  	return FALSE;
  }
  
  return TRUE;
}

static Boolean IsSequenceAllNs (CharPtr seq_str)
{
  CharPtr cp;
  
  if (StringHasNoText (seq_str)) return FALSE;
  
  for (cp = seq_str; *cp != 0; cp++)
  {
    if (isalpha ((Int4)(*cp)) && *cp != 'n' && *cp != 'N')
    {
      return FALSE;
    }
  }
  return TRUE;
}

static Boolean IsSequenceAllOneCharacter (CharPtr seq_str)
{
  CharPtr cp;
  Char    first_char = 0;
  
  if (StringHasNoText (seq_str)) return FALSE;
  
  for (cp = seq_str; *cp != 0; cp++)
  {
    if (isalpha ((Int4)(*cp)))
    {
      if (first_char == 0)
      {
        first_char = *cp;
      }
      else if (*cp != first_char)
      {
        return FALSE;
      }
    }
  }
  return TRUE;
  
}

static Int4 CountSeqChars (CharPtr seq_str)
{
  CharPtr cp;
  Int4    num_chars = 0;
  
  if (StringHasNoText (seq_str)) return 0;
  for (cp = seq_str; *cp != 0; cp++)
  {
    if (isalpha ((Int4)(*cp)))
    {
      num_chars++;
    }
  }
  return num_chars;
}

static CharPtr ReformatLocalId (CharPtr local_id)
{
  CharPtr cp, new_local_id;
  
  if (local_id == NULL) return NULL;
  
  cp = local_id;
  while (*cp == '>')
  {
    cp ++;
  }
  while (isspace ((Int4)(*cp)))
  {
  	cp++;
  }
  new_local_id = StringSave (cp);
  cp = new_local_id;
  while (*cp != 0)
  {
    if (isspace ((Int4)(*cp)))
    {
      *cp = '_';
    }
    cp++;
  }
  MemFree (local_id);
  return new_local_id;
}

static CharPtr FindPreviousWhitespace (CharPtr str_start, CharPtr str_end)
{
  CharPtr cp;
  if (str_start == NULL || str_end == NULL || str_end < str_start)
  {
    return NULL;
  }

  cp = str_end;
  while (cp > str_start && !isspace (*cp))
  {
    cp--;
  }
  return cp;
}

static CharPtr GetModNameStartFromEqLoc (CharPtr eq_loc, CharPtr prev_eq_loc)
{
  CharPtr cp, prev_quote;
  Char    match_quote;
  
  if (StringHasNoText (eq_loc) || StringHasNoText (prev_eq_loc) || eq_loc < prev_eq_loc)
  {
    return NULL;
  }
  
  cp = eq_loc - 1;
  /* skip over spaces between equals sign and modifier name */
  while (cp > prev_eq_loc && isspace (*cp))
  {
    cp--;
  }
  if (cp != prev_eq_loc)
  {
    /* now backtrack over name */
    if (*cp == '"' && *(cp - 1) != '\\')
    {
      match_quote = *cp;
      /* take everything to previous matching quote */
      cp = FindPreviousUnescapedQuote (prev_eq_loc, cp - 1);
      if (cp == NULL)
      {
        cp = prev_eq_loc;
      }
      if (isspace (*cp))
      {
        cp++;
      }
    }
    else
    {
      /* take everything up to the first space or quote character */
      
      prev_quote = FindPreviousUnescapedQuote (prev_eq_loc, cp);
      cp = FindPreviousWhitespace (prev_eq_loc, cp);
      if (prev_quote != NULL && prev_quote > cp)
      {
        cp = prev_quote + 1;
      }

      if (isspace (*cp) || *cp == '"')
      {
        cp++;
      }
    }
  }
  return cp;  
}

static CharPtr fake_modifier_name = "modifier_name";
static Int4    len_fake_modifier_name = 13;

static CharPtr InsertMissingModifierNames (CharPtr str)
{
  CharPtr start_bracket, next_start_bracket, end_bracket, eq_loc;
  CharPtr new_str, tmp_new;
  Int4    offset;
  
  if (str == NULL) return NULL;

  new_str = StringSave (str);
  start_bracket = StringChr (new_str, '[');
  while (start_bracket != NULL)
  {
    next_start_bracket = StringChr (start_bracket + 1, '[');
    eq_loc = StringChr (start_bracket, '=');
    end_bracket = StringChr (start_bracket, ']');
    
    if (eq_loc == NULL || end_bracket == NULL 
        || eq_loc > end_bracket
        || (next_start_bracket != NULL && end_bracket > next_start_bracket))
    {
      /* can't fix this pair, move along */    
    }
    else if (StringSpn (start_bracket + 1, " \t") == eq_loc - start_bracket - 1)
    {
      offset = start_bracket - new_str;
      tmp_new = InsertStringAtOffset (new_str, fake_modifier_name, offset + 1);
      if (tmp_new != NULL)
      {
        new_str = MemFree (new_str);
        new_str = tmp_new;
      }
      start_bracket = new_str + offset;
      next_start_bracket = StringChr (start_bracket + 1, '[');
    }
    start_bracket = next_start_bracket;
  }
  return new_str;
}

static CharPtr SuggestCorrectBracketing (CharPtr str)
{
  CharPtr cp, next_token;
  CharPtr new_str, tmp_new;
  Int4    offset, name_len, token_offset;
  CharPtr step_back, step_forward, next_next_token, next_next_next_token;
  Char    insert_buf [2];
  
  if (str == NULL) return NULL;
  
  new_str = StringSave (str);
  cp = new_str;
  next_token = NextBracketToken (cp);
  
  while (*cp != 0 && next_token != NULL)
  {
    if (*next_token == '"')
    {
      insert_buf [0] = *next_token;
      insert_buf [1] = 0;
      tmp_new = InsertStringAtOffset (new_str, insert_buf, StringLen (new_str));
      new_str = MemFree (new_str);
      new_str = tmp_new;
      return new_str;
    }
    next_next_token = NextBracketToken (next_token + 1);
    if (next_next_token == NULL)
    {
      next_next_next_token = NULL;
    }
    else
    {
      next_next_next_token = NextBracketToken (next_next_token + 1);
    }
    
    /* skip over correctly formatted bits */
    if (*next_token == '['
        && next_next_token != NULL
        && *next_next_token == '='
        && next_next_next_token != NULL
        && *next_next_next_token == ']')
    {
      cp = next_next_next_token + 1;
    }
    /* remove repeated tokens ([[, ]], or ==) in first pair*/
    else if (next_next_token != NULL 
             && *next_token == *next_next_token
             && next_next_token - next_token - 1== StringSpn (next_token + 1, " \t"))
    {
      ShiftString (next_token, next_next_token - next_token);
    }
    /* remove repeated tokens ([[, ]], or ==) in first pair*/
    else if (next_next_token != NULL 
             && next_next_next_token != NULL
             && *next_next_token == *next_next_next_token
             && next_next_next_token - next_next_token - 1== StringSpn (next_next_token + 1, " \t"))
    {
      ShiftString (next_next_token, next_next_next_token - next_next_token);
    }
    /* no start - either remove end token or insert start */
    else 
    {
      switch (*next_token)
      {
        case '=' :
          /* insert start before equals token */
          if (next_token == cp)
          {
            offset = cp - new_str;
          }
          else
          {
            step_back = GetModNameStartFromEqLoc (next_token, cp);
            offset = step_back - new_str;
          }
          tmp_new = InsertStringAtOffset (new_str, "[", offset);
          if (tmp_new != NULL)
          {
            new_str = MemFree (new_str);
            new_str = tmp_new;
          }
          cp = tmp_new + offset;
          break;
        case ']' :
          /* remove lonely end bracket */
          ShiftString (next_token, 1);
          cp = next_token;
          break;
        case '[' :
          next_next_token = NextBracketToken (next_token + 1);
          if (next_next_token == NULL
              || *next_next_token == '[')
          {
            /* remove unwanted beginning bracket */
            ShiftString (next_token, 1);
            cp = next_token;
          }
          else if (*next_next_token == '=')
          {
            /* find the best place to put a closing bracket */
            next_next_next_token = NextBracketToken (next_next_token + 1);
            if (next_next_next_token == NULL)
            {
              offset = StringLen (new_str);
            }
            else
            {
              token_offset = next_next_token - new_str + 1;
              offset = next_next_next_token - new_str;
              while (offset > token_offset && isspace (new_str[offset - 1]))
              {
                offset--;
              }
              if (*next_next_next_token == '=')
              {
                step_back = GetModNameStartFromEqLoc (next_next_next_token, next_next_token);
                if (step_back == next_next_token) {
                  /* no name before second equals sign, put bracket after value */
                } else {
                  while (step_back > next_next_token + 1
                        && isspace (*(step_back - 1)))
                  {
                    step_back --;
                  }
                  offset = step_back - new_str;
                }
              }
            }
            tmp_new = InsertStringAtOffset (new_str, "]", offset);
            if (tmp_new != NULL)
            {
              new_str = MemFree (new_str);
              new_str = tmp_new;
            }
            cp = tmp_new + offset + 1;            
          }
          else if (*next_next_token == ']')
          {
            /* see if we can insert an equals sign */
            /* skip over empty space after '[', if any */
            step_forward = next_token + 1 + StringSpn (next_token + 1, " \t");
            if (step_forward == next_next_token)
            {
              /* eliminate the empty bracket pair */
              ShiftString (next_token, next_next_token - next_token + 1);
              cp = next_token;
            }
            else
            {
              /* get length of first text token */
              name_len = StringCSpn (step_forward, " \t");
              if (next_next_token - step_forward < name_len)
              {
                name_len = next_next_token - step_forward;
              }
              if (step_forward + name_len < next_next_token && isspace (*(step_forward + name_len)))
              {
                *(step_forward + name_len) = '=';
                cp = next_next_token + 1;
              }
              else if (StringNICmp (step_forward, "DNA", name_len) == 0
                       || StringNICmp (step_forward, "RNA", name_len) == 0
                       || StringNICmp (step_forward, "orf", name_len) == 0)
              {
                cp = next_next_token + 1;
              }
              else
              {
                if ((next_token > new_str && isspace (*(next_token - 1)))
                    || isspace (*(next_token + 1))) 
                {
                  ShiftString (next_token, 1);
                  next_next_token --;
                }
                else
                {
                  *next_token = ' ';
                }
                if (isspace (*(next_next_token - 1)) || isspace (*(next_next_token + 1)))
                {
                  ShiftString (next_next_token, 1);
                }
                else
                {
                  *next_next_token = ' ';
                }
                cp = next_next_token;
              }
            }
          }
          break;
      }
    }
    
    next_token = NextBracketToken (cp);
  }
  
  tmp_new = InsertMissingModifierNames (new_str);
  if (tmp_new != NULL)
  {
    new_str = MemFree (new_str);
    new_str = tmp_new;
  }
  
  return new_str;
}

static void ClearSequencesButton (ButtoN b)
{
  SequencesFormPtr   sqfp;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL) return;
  SequencesFormDeleteProc (sqfp);
}

static void 
SetModifierList
(DoC doc,
 ValNodePtr mod_list)
{
  ValNodePtr vnp;
  Int4       num_modifiers = 0;
  Int4       text_len = 0;
  CharPtr    text;
  CharPtr    text_fmt = "Already have values for:\n";
    
  if (doc == NULL)
  {
    return;
  }
  Reset (doc);
  if (mod_list == NULL) 
  {
    text = StringSave ("No modifiers are present.");
  }
  else
  {
    text_len = StringLen (text_fmt) + 1;
    for (vnp = mod_list; vnp != NULL; vnp = vnp->next)
    {
      num_modifiers ++;
      text_len += StringLen (vnp->data.ptrvalue);
    }
    text_len += num_modifiers * 6;
    text = (CharPtr) MemNew (text_len * sizeof (Char));
    if (text != NULL)
    {
      StringCpy (text, text_fmt);
      for (vnp = mod_list; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->next == NULL && num_modifiers > 1)
        {
          StringCat (text, "and ");
        }
      	StringCat (text, vnp->data.ptrvalue);
      	if (vnp->next != NULL)
      	{
      	  if (num_modifiers > 2)
      	  {
      	  	StringCat (text, ", ");
      	  }
      	  else
      	  {
      	  	StringCat (text, " ");
      	  }
      	}
      }
    }
  }
  AppendText (doc, text, NULL, NULL, programFont);
  InvalDocRows (doc, 0, 0, 0);
  MemFree (text);  
}



static Int4 GetValForEnumName (EnumFieldAssocPtr eap, CharPtr mod_value)
{
  if (StringHasNoText (mod_value) || eap == NULL)
  {
    return 0;
  }
  while (eap != NULL && eap->name != NULL) 
  {
    if (StringICmp (eap->name, mod_value) == 0)
    {
      return eap->value;
    }
    eap++;
  }
  return 0;
}

static CharPtr 
TagListStringFromDefLineValue 
(CharPtr defline_val,
 Boolean is_nontext, 
 Int2    mod_type)
{
  CharPtr taglist_str = NULL;
  Char        text [128];
  
  if (is_nontext)
  {
    if (StringHasNoText (defline_val))
    {
      taglist_str = StringSave ("0");
    }
    else
    {
      taglist_str = StringSave ("1");
    }
  }
  else if (mod_type == eModifierType_Organism)
  {
    taglist_str = StringSave (defline_val);
  }
  else if (mod_type == eModifierType_Location)
  {
    if (StringHasNoText (defline_val))
    {
      taglist_str = StringSave ("1");
    }
    else
    {
      sprintf (text, "%d", GetValForEnumName (biosource_genome_simple_alist,
                                              defline_val));
      taglist_str = StringSave (text);
    }
  }
  else if (mod_type == eModifierType_Origin)
  {
    if (StringHasNoText (defline_val))
    {
      taglist_str = StringSave ("1");
    }
    else
    {
      sprintf (text, "%d", GetValForEnumName (biosource_origin_alist,
                                              defline_val));
      taglist_str = StringSave (text);
    }
  }  
  else if (mod_type == eModifierType_NucGeneticCode
           || mod_type == eModifierType_MitoGeneticCode)
  {
    if (StringHasNoText (defline_val))
    {
      taglist_str = StringSave ("0");
    }
    else
    {
      sprintf (text, "%d", GeneticCodeFromString (defline_val));
      taglist_str = StringSave (text);
    }
  }
  else if (mod_type == eModifierType_MolType)
  {
    if (StringHasNoText (defline_val))
    {
      taglist_str = StringSave ("253");
    }
    else
    {
      sprintf (text, "%d", MolTypeFromString (defline_val));
      taglist_str = StringSave (text);
    }
  }
  else if (mod_type == eModifierType_Molecule)
  {
    if (StringICmp (defline_val, "dna") == 0)
    {
      sprintf (text, "%d", Seq_mol_dna);
      taglist_str = StringSave (text);
    }
    else if (StringICmp (defline_val, "rna") == 0)
    {
      sprintf (text, "%d", Seq_mol_rna);
      taglist_str = StringSave (text);
    }
    else
    {
      sprintf (text, "%d", Seq_mol_dna);
      taglist_str = StringSave (text);
    }
  }
  else if (mod_type == eModifierType_Topology)
  {
    if (StringHasNoText (defline_val))
    {
      sprintf (text, "%d", TopologyFromString (""));
    }
    else
    {
      sprintf (text, "%d", TopologyFromString (defline_val));
    }
    taglist_str = StringSave (text);
  }
  else
  {
    if (StringHasNoText (defline_val))
    {
      taglist_str = StringSave (" ");
    }
    else
    {
      taglist_str = StringSave (defline_val);
    }
  }
  return taglist_str;  
}

static void AddSeqIDAndValueToTagList 
(CharPtr id,
 CharPtr title,
 CharPtr mod_name,
 ValNodePtr PNTR   head)
{
  Char        text [2];
  CharPtr     str;
  Int4        len;
  Int4        mod_type;
  Boolean     is_nontext;
  CharPtr     val_str, taglist_str = NULL;

  if (head == NULL)
  {
    return;
  }

  is_nontext = IsNonTextModifier (mod_name);
  mod_type = GetModifierType (mod_name);

  text [0] = '\0';

  if (is_nontext)
  {
    if (FindValuePairInDefLine (mod_name, title, NULL))
    {
      sprintf (text, "2");
    }
    else
    {
      text [0] = '\0';
    }
    val_str = StringSave (text);
  }
  else
  {
    val_str = FindValueFromPairInDefline (mod_name, title); 
  }
  
  taglist_str = TagListStringFromDefLineValue (val_str, is_nontext, mod_type);
  val_str = MemFree (val_str);
  len = StringLen (id) + StringLen (taglist_str);
  str = MemNew (len + 4);
  if (str != NULL) {
    StringCpy (str, id);
    StringCat (str, "\t");
    StringCat (str, taglist_str);
    StringCat (str, "\n");
  }
  
  taglist_str = MemFree (taglist_str);
  ValNodeAddPointer (head, 0, str);
  
}

NLM_EXTERN CharPtr GetValueFromTitle (CharPtr mod_name, CharPtr title)
{
  Int4        mod_type;
  Boolean     is_nontext;
  CharPtr     valstr;
  
  if (StringHasNoText (mod_name) || StringHasNoText (title))
  {
    return NULL;
  }
  mod_type = GetModifierType (mod_name);
  is_nontext = IsNonTextModifier (mod_name);
  
  if (mod_type == eModifierType_Organism)
  {
    valstr = FindValueFromPairInDefline ("organism", title);
    if (StringHasNoText (valstr))
    {
      valstr = MemFree (valstr);
      valstr = StringSave (" ");
    }
  }
  else if (mod_type == eModifierType_Location)
  {
    valstr = NULL;
    if (FindValuePairInDefLine ("location", title, NULL) != NULL)
    {
      valstr = FindValueFromPairInDefline ("location", title);
    }
    else
    {
      valstr = StringSave ("genomic");
    }
  }
  else if (IsNonTextModifier (mod_name))
  {
    if (FindValuePairInDefLine (mod_name, title, NULL) != NULL)
    {
      valstr = StringSave ("TRUE");
    }
    else
    {
      valstr = StringSave ("FALSE");
    }
  }
  else
  {
    valstr = FindValueFromPairInDefline (mod_name, title);
    if (StringHasNoText (valstr))
    {
      valstr = MemFree (valstr);
      valstr = StringSave (" ");
    }
  }      
  return valstr;  
}

/* This function returns a string suitable for display in a table
 * using the values from the specified modifier name found in the
 * specified table.
 * non-text modifiers are displayed as either TRUE or FALSE,
 * multiple values are listed in parentheses with semicolons between them.
 */
static CharPtr GetDisplayValue (CharPtr mod_name, CharPtr title, BoolPtr multi_found)
{
  CharPtr begin_bracket, end_bracket;
  CharPtr mod_value = NULL, tmp_value;
  Int4    mod_value_len = 0;
  ValNodePtr val_list = NULL, vnp;
  Boolean allow_multi;
  
  allow_multi = AllowMultipleValues (mod_name);
  if (allow_multi)
  {
    begin_bracket = FindValuePairInDefLine (mod_name, title, &end_bracket);
    while (begin_bracket != NULL)
    {
      tmp_value = GetValueFromTitle (mod_name, begin_bracket);
      if (!StringHasNoText (tmp_value))
      {
        mod_value_len += StringLen (tmp_value) + 1;
        ValNodeAddPointer (&val_list, 0, tmp_value);
      }
      else
      {
        tmp_value = MemFree (tmp_value);
      }
      begin_bracket = FindValuePairInDefLine (mod_name, end_bracket + 1, &end_bracket);
    }
    if (val_list == NULL)
    {
      mod_value = StringSave (" ");
    }
    else
    {
      if (val_list->next == NULL)
      {
        mod_value = val_list->data.ptrvalue;
        val_list->data.ptrvalue = NULL;
      }
      else
      {
        mod_value = (CharPtr) MemNew ((mod_value_len + 3) * sizeof (Char));
        if (mod_value != NULL)
        {
          mod_value [0] = '(';
          for (vnp = val_list; vnp != NULL; vnp = vnp->next)
          {
            StringCat (mod_value, vnp->data.ptrvalue);
            if (vnp->next == NULL)
            {
              StringCat (mod_value, ")");
            }
            else
            {
              StringCat (mod_value, ",");
            }
          }
          if (multi_found != NULL)
          {
            *multi_found = TRUE;
          }
        }
      }
      val_list = ValNodeFree (val_list);
    }
  }
  else
  {
    mod_value = GetValueFromTitle (mod_name, title);
  }
  return mod_value;
}

/* This function returns a string suitable for display in a table
 * using the values from the specified modifier name found in the
 * specified table.
 * non-text modifiers are displayed as either TRUE or FALSE,
 * multiple values are listed in parentheses with semicolons between them.
 */
static CharPtr GetDisplayValueFromModifierInfoList (CharPtr mod_name, ValNodePtr modifier_info_list, BoolPtr multi_found)
{
  CharPtr mod_value = NULL;
  Int4    mod_value_len = 0;
  ValNodePtr val_list = NULL, vnp;
  Boolean allow_multi;
  ModifierInfoPtr mip;
  Uint1           subtype = 0;
  
  allow_multi = AllowMultipleValues (mod_name);
  
  subtype = FindTypeForModNameText (mod_name);

  for (vnp = modifier_info_list; vnp != NULL; vnp = vnp->next) {
    mip = (ModifierInfoPtr) vnp->data.ptrvalue;
    if (mip != NULL
        && !StringHasNoText (mip->value)
        && ((mip->subtype != 0 && subtype == mip->subtype)
            || StringsAreEquivalent (mod_name, mip->name))) {
      ValNodeAddPointer (&val_list, 0, mip->value);
      mod_value_len += StringLen (mip->value) + 1;
      if (!allow_multi) {
        break;
      }
    }
  }

  if (val_list == NULL)
  {
    mod_value = StringSave (" ");
  }
  else
  {
    if (val_list->next == NULL)
    {
      mod_value = StringSave (val_list->data.ptrvalue);
    }
    else
    {
      mod_value = (CharPtr) MemNew ((mod_value_len + 3) * sizeof (Char));
      if (mod_value != NULL)
      {
        mod_value [0] = '(';
        for (vnp = val_list; vnp != NULL; vnp = vnp->next)
        {
          StringCat (mod_value, vnp->data.ptrvalue);
          if (vnp->next == NULL)
          {
            StringCat (mod_value, ")");
          }
          else
          {
            StringCat (mod_value, ",");
          }
        }
        if (multi_found != NULL)
        {
          *multi_found = TRUE;
        }
      }
    }
    val_list = ValNodeFree (val_list);
  }
  return mod_value;
}

static Boolean IntValueInValNodeList (Int4 ival, ValNodePtr vnp)
{
  Boolean found_int = FALSE;
  
  while (vnp != NULL && !found_int)
  {
    if (vnp->data.intvalue == ival)
    {
      found_int = TRUE;
    }
    vnp = vnp->next;
  }
  
  return found_int;
}

static Boolean 
DoColumnlistsHaveIdenticalSourceInformation 
(ValNodePtr col1,
 ValNodePtr col2,
 ValNodePtr header)
{
  Boolean are_identical = TRUE;
  Int4    mod_type;
  
  while (col1 != NULL && col2 != NULL && header != NULL && are_identical)
  { 
    mod_type = GetModifierType (header->data.ptrvalue);
    if (mod_type != eModifierType_Protein
        && mod_type != eModifierType_MolType
        && mod_type != eModifierType_Topology
        && mod_type != eModifierType_Molecule
        && StringCmp ((CharPtr) col1->data.ptrvalue, (CharPtr) col2->data.ptrvalue) != 0)
    {
      are_identical = FALSE;
    }
    col1 = col1->next;
    col2 = col2->next;
    header = header->next;
  }
  if ((col1 == NULL && col2 != NULL) || (col1 != NULL && col2 == NULL))
  {
    are_identical = FALSE;
  }
  return are_identical;
}


static CharPtr GetSourceInformationColumns (ValNodePtr header, ValNodePtr columns)
{
  Int4 len = 0;
  CharPtr rval = NULL;
  ValNodePtr col, hcol;
  Int4    mod_type;

  col = columns;
  hcol = header;
  while (col != NULL && hcol != NULL)
  { 
    mod_type = GetModifierType (hcol->data.ptrvalue);
    if (mod_type != eModifierType_Protein
        && mod_type != eModifierType_MolType
        && mod_type != eModifierType_Topology
        && mod_type != eModifierType_Molecule)
    {
      len += StringLen (col->data.ptrvalue) + 1;
    }
    col = col->next;
    hcol = hcol->next;
  }

  rval = (CharPtr) MemNew (sizeof (Char) * (len + 1));

  col = columns;
  hcol = header;
  while (col != NULL && hcol != NULL)
  { 
    mod_type = GetModifierType (hcol->data.ptrvalue);
    if (mod_type != eModifierType_Protein
        && mod_type != eModifierType_MolType
        && mod_type != eModifierType_Topology
        && mod_type != eModifierType_Molecule)
    {
      StringCat (rval, col->data.ptrvalue);
      StringCat (rval, ":");
    }
    col = col->next;
    hcol = hcol->next;
  }
  return rval;
}


static Boolean HasAnySourceInformation (ValNodePtr header_list, ValNodePtr column_list)
{
  Boolean has_any = FALSE;
  Int4    mod_type;
  
  if (header_list == NULL || column_list == NULL)
  {
    return FALSE;
  }
  
  /* skip over SeqID column */
  header_list = header_list->next;
  column_list = column_list->next;
  while (header_list != NULL && column_list != NULL && ! has_any)
  {
    if (!StringHasNoText (column_list->data.ptrvalue))
    {
      mod_type = GetModifierType (header_list->data.ptrvalue);
      if (mod_type != eModifierType_Protein
          && mod_type != eModifierType_MolType
          && mod_type != eModifierType_Molecule
          && mod_type != eModifierType_Topology)
      {
        has_any = TRUE;
      }
    }
    header_list = header_list->next;
    column_list = column_list->next;
  }
  return has_any;
}

static Boolean OrganismMatchesAnotherRow (Int4 row, ValNodePtr row_list, Pointer userdata)
{
  ValNodePtr header_vnp, column_list, check_column_list, row_vnp;
  Int4       row_num;
  
  if (row_list == NULL || row < 1)
  {
    return FALSE;
  }
  
  /* we start with the header of the second column, because the first column
   * is the sequence ID */
  header_vnp = row_list->data.ptrvalue;
  if (header_vnp == NULL) return FALSE;
  header_vnp = header_vnp->next;
  if (header_vnp == NULL) return FALSE;
    
  /* find the row we're interested in */
  for (row_vnp = row_list->next, row_num = 1;
       row_vnp != NULL && row_num != row; 
       row_vnp = row_vnp->next, row_num++)
  {
  }
  if (row_vnp == NULL)
  {
    return FALSE;
  }
  
  column_list = (ValNodePtr) row_vnp->data.ptrvalue;
  if (!HasAnySourceInformation (row_list->data.ptrvalue, column_list))
  {
    return FALSE;
  }
  if (column_list == NULL || column_list->next == NULL)
  {
    return FALSE;
  }
  
  /* don't check when organism name is missing */
  if (StringHasNoText (column_list->next->data.ptrvalue))
  {
    return FALSE;
  }
  
  /* now check it against the other rows */
  for (row_vnp = row_list->next, row_num = 1;
       row_vnp != NULL; 
       row_vnp = row_vnp->next, row_num++)
  {
    if (row_num == row)
    {
      continue;
    }
    
    check_column_list = (ValNodePtr) row_vnp->data.ptrvalue;
    if (check_column_list == NULL || check_column_list->next == NULL)
    {
      continue;
    }

    /* we compare the column lists, starting with the second column
     * because the first column contains the sequence ID
     */
    if (DoColumnlistsHaveIdenticalSourceInformation (column_list->next,
                                                     check_column_list->next,
                                                     header_vnp))
    {
      return TRUE;
    }
  }
  return FALSE;
}

/* Sequence ID is always stored in the first column, organism name
 * is always stored in the second column.
 */
static Boolean AnySequencesHaveMissingOrganisms (ValNodePtr row_list)
{
  Boolean have_missing = FALSE;
  ValNodePtr col_list;
  
  if (row_list == NULL || row_list->next == NULL)
  {
    return FALSE;
  }
  row_list = row_list->next;
  
  while (row_list != NULL && ! have_missing)
  {
    col_list = row_list->data.ptrvalue;
    if (col_list == NULL
        || col_list->next == NULL
        || StringHasNoText (col_list->next->data.ptrvalue))
    {
      have_missing = TRUE;
    }
    row_list = row_list->next;
  }
  return have_missing;
}

static Boolean AnySequencesHaveIdenticalOrganisms (ValNodePtr row_list)
{
  Boolean have_match = FALSE;
  ValNodePtr vnp;
  ValNodeBlock org_list;
  CharPtr prev;
  
  if (row_list == NULL || row_list->next == NULL || row_list->next->next == NULL)
  {
    return FALSE;
  }
  InitValNodeBlock (&org_list, NULL);
  for (vnp = row_list->next; vnp != NULL; vnp = vnp->next) 
  {
    ValNodeAddPointerToEnd (&org_list, 0, GetSourceInformationColumns (row_list->data.ptrvalue, vnp->data.ptrvalue));
  }

  org_list.head = ValNodeSort (org_list.head, SortVnpByString);
  prev = org_list.head->data.ptrvalue;
  for (vnp = org_list.head->next; vnp != NULL && !have_match; vnp = vnp->next)
  {
    if (StringICmp (prev, vnp->data.ptrvalue) == 0) 
    {
      have_match = TRUE;
    }
    prev = vnp->data.ptrvalue;
  }
  org_list.head = ValNodeFreeData (org_list.head);

  return have_match;
}

/* Sequence ID is always stored in the first column, organism name
 * is always stored in the second column.
 */
static void ReportMissingOrganisms (ValNodePtr row_list, DoC doc)
{
  Int4       row_num, num_missing;
  ValNodePtr row_vnp;
  ValNodePtr column_list;
  ValNodeBlock missing_list;
  CharPtr    err_msg;
  CharPtr    missing_fmt = "%d sequences are missing organism names.";
  
  if (row_list == NULL || doc == NULL) return;
      

  InitValNodeBlock (&missing_list, NULL);

  for (row_vnp = row_list->next, row_num = 0;
       row_vnp != NULL; 
       row_vnp = row_vnp->next, row_num++)
  {    
    column_list = (ValNodePtr) row_vnp->data.ptrvalue;
    if (column_list == NULL)
    {
      continue;
    }
    if (column_list->next == NULL || StringHasNoText (column_list->next->data.ptrvalue))
    {
      /* organism is missing */
      ValNodeAddPointerToEnd (&missing_list, row_num, column_list->data.ptrvalue);          
    }
  }
  
  if (missing_list.head != NULL)
  {
    num_missing = ValNodeLen (missing_list.head);
    if (num_missing > 100) 
    {
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + 15));
      sprintf (err_msg, missing_fmt, num_missing);
    } 
    else 
    {
      err_msg = CreateListMessage ("Sequence", 
                                   missing_list.head->next == NULL ? 
                                             " has no organism name."
                                             : " have no organism names.",
                                   missing_list.head);
    }
    AppendText (doc, err_msg, &faParFmt, &faColFmt, programFont);
    err_msg = MemFree (err_msg);                                 
    missing_list.head = ValNodeFree (missing_list.head);
    AppendText (doc, "\n", &faParFmt, &faColFmt, programFont);
  }
}

static void ReportIdenticalOrganisms (ValNodePtr row_list, DoC doc)
{
  ValNodePtr checked_list = NULL;
  Int4       row_num, check_row_num;
  ValNodePtr row_vnp, check_row_vnp;
  ValNodePtr column_list, check_column_list;
  ValNodePtr header_vnp;
  Boolean    skip_this;
  ValNodePtr this_match_list;
  CharPtr    err_msg;
  Boolean    any_data_reported = FALSE;
  
  if (row_list == NULL || doc == NULL) return;
  
  /* we start with the header of the second column, because the first column
   * is the sequence ID */
  header_vnp = row_list->data.ptrvalue;
  if (header_vnp == NULL) return;
  header_vnp = header_vnp->next;
  if (header_vnp == NULL) return;
    
  for (row_vnp = row_list->next, row_num = 0;
       row_vnp != NULL; 
       row_vnp = row_vnp->next, row_num++)
  {
    /* don't need to check rows that have matched a previous row */
    skip_this = IntValueInValNodeList (row_num, checked_list);
    if (skip_this)
    {
      continue;
    }
    
    column_list = (ValNodePtr) row_vnp->data.ptrvalue;
    if (!HasAnySourceInformation (row_list->data.ptrvalue, column_list))
    {
      continue;
    }

    if (column_list == NULL || column_list->next == NULL)
    {
      continue;
    }
    
    if (StringHasNoText (column_list->next->data.ptrvalue))
    {
      /* skip - no organism name, will have already been reported */
      continue;
    }
    
    this_match_list = NULL;
    for (check_row_vnp = row_vnp->next, check_row_num = row_num + 1;
         check_row_vnp != NULL;
         check_row_vnp = check_row_vnp->next, check_row_num++)
    {
      skip_this = IntValueInValNodeList (row_num, checked_list);
      if (skip_this)
      {
        continue;
      }
      check_column_list = (ValNodePtr) check_row_vnp->data.ptrvalue;
      if (check_column_list == NULL || check_column_list->next == NULL)
      {
        continue;
      }

      /* we compare the column lists, starting with the second column
       * because the first column contains the sequence ID
       */
      if (DoColumnlistsHaveIdenticalSourceInformation (column_list->next,
                                                       check_column_list->next,
                                                       header_vnp))
      {
        /* be sure to put the first row to match the other rows in the list */
        if (this_match_list == NULL)
        {
          ValNodeAddPointer (&this_match_list, row_num, column_list->data.ptrvalue);          
        }
        /* add the sequence ID for the check row to the list */
        ValNodeAddPointer (&this_match_list, check_row_num, check_column_list->data.ptrvalue);   
        ValNodeAddInt (&checked_list, 0, check_row_num);       
      }
    }
    
    /* if anything matched this row, put the list in the list of matches */
    if (this_match_list != NULL)
    {
      err_msg = CreateListMessage ("Sequence", 
                     " have identical source information.",
                     this_match_list);
      AppendText (doc, err_msg, &faParFmt, &faColFmt, programFont);
      err_msg = MemFree (err_msg);
      this_match_list = ValNodeFree (this_match_list);
      any_data_reported = TRUE;
    }
  }
  
  checked_list = ValNodeFree (checked_list);
  if (any_data_reported)
  {
    AppendText (doc, "\n", &faParFmt, &faColFmt, programFont);
  }
}

static void SummarizeModifiers (ValNodePtr row_list, DialoG summary_dlg)
{
  ValNodePtr header_vnp, row_vnp, column_vnp;
  Int4       column_offset, col_pos;
  Boolean    any_present;
  Boolean    all_present;
  Boolean    is_unique;
  CharPtr    first_value_seen;
  Boolean    all_unique;
  ValNodePtr values_seen;
  CharPtr    row_status;
  Int4       line_len;
  CharPtr    modifier_line = NULL;
  Int4       num_missing;
  ValNodePtr summary_row_list = NULL;
  ValNodePtr summary_col_list = NULL, summary_header_list = NULL;
  
  if (row_list == NULL || row_list->next == NULL || summary_dlg == NULL)
  {
    return;
  }
  
  summary_col_list = NULL;
  ValNodeAddPointer (&summary_col_list, 8, StringSave ("Modifier"));
  ValNodeAddPointer (&summary_col_list, 6, StringSave ("Status"));
  ValNodeAddPointer (&summary_col_list, 11, StringSave ("First Value"));
  ValNodeAddPointer (&summary_row_list, 0, summary_col_list);
  summary_header_list = summary_col_list;
  
  header_vnp = row_list->data.ptrvalue;
  /* skip over sequence ID column */
  header_vnp = header_vnp->next;
  column_offset = 1;
  while (header_vnp != NULL)
  {
    any_present = FALSE;
    all_present = TRUE;
    is_unique = TRUE;
    all_unique = TRUE;
    first_value_seen = NULL;
    values_seen = NULL;
    num_missing = 0;
    
    /* skip over header line */
    row_vnp = row_list->next;
    while (row_vnp != NULL)
    {
      for (col_pos = 0, column_vnp = row_vnp->data.ptrvalue;
           col_pos < column_offset && column_vnp != NULL;
           col_pos++, column_vnp = column_vnp->next)
      {
      }
      if (column_vnp == NULL)
      {
        continue;
      }
      if (StringHasNoText (column_vnp->data.ptrvalue))
      {
        all_present = FALSE;
        num_missing++;
      }
      else
      {
        any_present = TRUE;
        if (first_value_seen == NULL)
        {
          first_value_seen = StringSave (column_vnp->data.ptrvalue);
          ValNodeAddPointer (&values_seen, 0, first_value_seen);
        }
        else 
        {
          if (StringCmp (first_value_seen, column_vnp->data.ptrvalue) != 0)
          {
            is_unique = FALSE;
          }
          
          if ( FindExactStringListMatch (values_seen, column_vnp->data.ptrvalue)
              == NULL)
          {
            ValNodeAddStr (&values_seen, 0, column_vnp->data.ptrvalue);
          }
          else
          {
            all_unique = FALSE;
          }
        }
      }
      row_vnp = row_vnp->next;
    }
    
    /* add summary line for this modifier */
    if (! any_present)
    {
      row_status = "All missing (%d sequences)";
    }
    else if (all_present && all_unique)
    {
      row_status = "All present, all unique values";
    }
    else if (all_present && is_unique)
    {
      row_status = "All present, one unique value";
    }
    else if (all_present && ! is_unique)
    {
      row_status = "All present, mixed values";
    }
    else if (! all_present && all_unique)
    {
      row_status = "%d missing, all unique values";
    }
    else if (! all_present && is_unique)
    {
      row_status = "%d missing, one unique value present";
    }
    else if (! all_present && ! is_unique)
    {
      row_status = "%d missing, mixed values";
    }
    
    line_len = StringLen (row_status) + 30;
               
    modifier_line = (CharPtr) MemNew (line_len * sizeof (Char));
    if (modifier_line != NULL)
    {
      
      /* add summary row for this modifier */
      summary_col_list = NULL;
      /* add modifier name */
      ValNodeAddPointer (&summary_col_list, 
                         0, 
                         StringSave (header_vnp->data.ptrvalue));
      /* show up to the first fifteen characters of the modifier name */
      summary_header_list->choice = MAX (summary_header_list->choice,
                                         StringLen (header_vnp->data.ptrvalue));
      summary_header_list->choice = MIN (summary_header_list->choice,
                                         15);
      
      /* add status */
      if (all_present)
      {
        ValNodeAddPointer (&summary_col_list, 
                           0, 
                           StringSave (row_status));
        summary_header_list->next->choice = MAX (summary_header_list->next->choice,
                                         StringLen (row_status));
      }
      else
      {
        sprintf (modifier_line, row_status, num_missing);
        ValNodeAddPointer (&summary_col_list,
                           0,
                           StringSave (modifier_line));
        summary_header_list->next->choice = MAX (summary_header_list->next->choice,
                                         StringLen (modifier_line));
      }
      
      /* add sample value */
      if (StringHasNoText (first_value_seen))
      {
        ValNodeAddPointer (&summary_col_list, 0, StringSave (""));
      }
      else
      {
        ValNodeAddPointer (&summary_col_list,
                           0, 
                           StringSave (first_value_seen));
        summary_header_list->next->next->choice = MAX (summary_header_list->next->next->choice,
                                         StringLen (first_value_seen));
      }
      ValNodeAddPointer (&summary_row_list, 0, summary_col_list);
      
      modifier_line = MemFree (modifier_line);
    }
        
    /* free up variables */
    first_value_seen = MemFree (first_value_seen);
    values_seen = ValNodeFree (values_seen);
    
    header_vnp = header_vnp->next;
    column_offset++;    
  }  
  
  PointerToDialog (summary_dlg, summary_row_list);
  summary_row_list = FreeTableDisplayRowList (summary_row_list);
  
}

static ValNodePtr GetListOfCurrentSourceModifiers (IDAndTitleEditPtr iatep)
{
  ValNodePtr  found_modifiers = NULL;
  Int4        seq_num;

  /* we always list organism, and list it first, whether it's present or not */  
  ValNodeAddPointer (&found_modifiers, 0, StringSave ("Organism"));
  if (iatep != NULL)
  {
    /* get list of modifiers from titles */    
    for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
    {
      /* only add modifiers from master sequences */
      if (iatep->is_seg == NULL || !iatep->is_seg [seq_num])
      {
        found_modifiers = BuildModifierTypeList (found_modifiers, 
                                                 iatep->title_list [seq_num],
                                                 FALSE);
      }
    }
  }
  return found_modifiers;
}


static CharPtr multival_explanation = "Note: When there is more than one "
           "modifier of the same type for a single sequence, the value list will "
           "be presented in tables separated by semicolons and enclosed in parentheses.";

static void SeqEntryPtrToOrgDoc (SequencesFormPtr sqfp)
{
  SeqEntryPtr  seq_list;
  ValNodePtr   found_modifiers = NULL, vnp;
  CharPtr      mod_name;
  CharPtr      org_name;
  Int4         seq_num;
  ValNodeBlock column_list;
  ValNodeBlock row_list;
  ValNodePtr   row_vnp;
  ValNodePtr   header_vnp, header_list;
  Int4         column_width;
  CharPtr      mod_value;
  RecT              r;
  IDAndTitleEditPtr iatep;
  Boolean           multi_found = FALSE, have_missing, have_match;
  
  if (sqfp == NULL) return;
  Reset (sqfp->org_doc);
  ObjectRect (sqfp->org_doc, &r);
  InsetRect (&r, 4, 4);
  faColFmt.pixWidth = r.right - r.left;

  InitValNodeBlock (&row_list, NULL);

  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    AppendText (sqfp->org_doc, 
                "You must create sequences before you can add source information.",
                &faParFmt, &faColFmt, programFont);
    Show (sqfp->org_doc);
    Hide (sqfp->ident_org_grp);
    Hide (sqfp->summary_dlg);
  }
  else
  {
    Show (sqfp->summary_dlg);
    /* get list of modifiers */
    iatep = SeqEntryListToIDAndTitleEdit (seq_list);

    found_modifiers = GetListOfCurrentSourceModifiers (iatep);
    
    /* create header line for table */
    /* store max column width in choice */
    InitValNodeBlock (&column_list, NULL);
    ValNodeAddPointerToEnd (&column_list, 6, StringSave ("Seq ID"));
    ValNodeAddPointerToEnd (&column_list, 8, StringSave ("Organism"));
    for (vnp = found_modifiers->next; vnp != NULL; vnp = vnp->next)
    {
      ValNodeAddPointerToEnd (&column_list, StringLen (vnp->data.ptrvalue), StringSave ((CharPtr) vnp->data.ptrvalue));
    } 
    
    ValNodeAddPointerToEnd (&row_list, 0, column_list.head);
    header_list = column_list.head;
    
    /* create data lines for table */
    for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
    {
      /* only add rows for master sequences */
      if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
      {
        continue;
      }

      /* add modifiers from this title */
      InitValNodeBlock (&column_list, NULL);
      header_vnp = header_list;
      
      column_width = MAX (StringLen (iatep->id_list [seq_num]), header_vnp->choice);
      header_vnp->choice = column_width;
      ValNodeAddPointerToEnd (&column_list, 0, StringSave (iatep->id_list [seq_num]));
      
      /* add organism name */
      header_vnp = header_vnp->next;
      org_name = GetDisplayValue ("organism", iatep->title_list [seq_num], &multi_found);
      column_width = MAX (StringLen (org_name), header_vnp->choice);
      header_vnp->choice = column_width;
      ValNodeAddPointerToEnd (&column_list, 0, org_name);
      
      /* get remaining modifiers */
      for (vnp = found_modifiers->next; vnp != NULL; vnp = vnp->next)
      {
        header_vnp = header_vnp->next;
        mod_name = (CharPtr) vnp->data.ptrvalue;
        mod_value = GetDisplayValue (mod_name, iatep->title_list [seq_num], &multi_found);
        column_width = MAX (StringLen (mod_value), header_vnp->choice);
        header_vnp->choice = column_width;
        ValNodeAddPointerToEnd (&column_list, 0, mod_value);
      }
      ValNodeAddPointerToEnd (&row_list, 0, column_list.head);
    }
    have_missing = AnySequencesHaveMissingOrganisms (row_list.head);
    have_match = AnySequencesHaveIdenticalOrganisms (row_list.head);
    if (have_match || have_missing)
    {
      if (have_match)
      {
        Show (sqfp->ident_org_grp);
      }
      Show (sqfp->org_doc);
    }
    else
    {
      Hide (sqfp->ident_org_grp);
      Hide (sqfp->org_doc);
    }
    ReportMissingOrganisms (row_list.head, sqfp->org_doc);
    ReportIdenticalOrganisms (row_list.head, sqfp->org_doc);
    if (multi_found)
    {
      AppendText (sqfp->org_doc, multival_explanation, NULL, NULL, programFont);
      Show (sqfp->org_doc);
    }

    SummarizeModifiers (row_list.head, sqfp->summary_dlg);

    /* free table text */
    for (row_vnp = row_list.head; row_vnp != NULL; row_vnp = row_vnp->next)
    {
      column_list.head = (ValNodePtr) row_vnp->data.ptrvalue;
      row_vnp->data.ptrvalue = ValNodeFreeData (column_list.head);
    }
    row_list.head = ValNodeFree (row_list.head);

    ValNodeFreeData (found_modifiers);
    iatep = IDAndTitleEditFree (iatep);
  }
  /* update document */
  InvalDocRows (sqfp->org_doc, 0, 0, 0);
}

static void SeqEntryPtrToSourceTab (SequencesFormPtr sqfp)

{
  SeqEntryPtrToOrgDoc (sqfp);
}

static ValNodePtr GetFastaModifierList (Boolean allow_nuc, Boolean allow_prot)
{
  ValNodePtr mod_choices = NULL;
  Int4       i;

  if (allow_nuc)
  {
    ValNodeAddPointer (&mod_choices, eModifierType_Organism, StringSave ("Organism"));
  
    ValNodeLink (&mod_choices, GetSourceQualDescListEx (TRUE, TRUE, FALSE, FALSE, FALSE));
    ValNodeAddPointer (&mod_choices, eModifierType_CommonName, StringSave ("Common Name"));
    ValNodeAddPointer (&mod_choices, eModifierType_Location, StringSave ("Location"));
    ValNodeAddPointer (&mod_choices, eModifierType_Origin, StringSave ("Origin"));
    ValNodeAddPointer (&mod_choices, eModifierType_Lineage, StringSave ("Lineage"));
    ValNodeAddPointer (&mod_choices, eModifierType_NucGeneticCode, StringSave ("gcode"));
    ValNodeAddPointer (&mod_choices, eModifierType_MitoGeneticCode, StringSave ("mgcode"));
    ValNodeAddPointer (&mod_choices, eModifierType_Molecule, StringSave ("moltype"));  
    ValNodeAddPointer (&mod_choices, eModifierType_Molecule, StringSave ("molecule"));
    ValNodeAddPointer (&mod_choices, eModifierType_Technique, StringSave ("tech"));
  }
  
  if (allow_prot)
  {
    for (i = 0; i < num_protein_modifier_names; i++)
    {
      if (StringICmp (protein_modifier_names [i], "orf") != 0)
      {
        ValNodeAddPointer (&mod_choices, eModifierType_Protein, StringSave (protein_modifier_names[i]));
      }
    }
  }

  return mod_choices;  
}

static CharPtr ReplaceOneModifierName (CharPtr title, CharPtr orig_name, CharPtr repl_name)
{
  CharPtr bracket_loc, eq_loc, new_title;
  Int4    new_title_len, search_offset;
  
  if (StringHasNoText (title)
      || StringHasNoText (orig_name)
      || StringHasNoText (repl_name)
      || StringICmp (orig_name, repl_name) == 0)
  {
    return title;
  }
  
  bracket_loc = FindValuePairInDefLine (orig_name, title, NULL);
  while (bracket_loc != NULL)
  {  
    eq_loc = NextBracketToken (bracket_loc + 1);
    if (eq_loc == NULL || *eq_loc != '=')
    {
      return title;
    }
    new_title_len = StringLen (title) + StringLen (repl_name) + 1;
    new_title = (CharPtr) MemNew (new_title_len * sizeof (Char));
    if (new_title == NULL)
    {
      return title;
    }
    StringNCpy (new_title, title, bracket_loc - title + 1);
    StringCat (new_title, repl_name);
    search_offset = StringLen (new_title) + 1;
    StringCat (new_title, eq_loc);
    title = MemFree (title);
    title = new_title;
    bracket_loc = FindValuePairInDefLine (orig_name, title + search_offset, NULL);
  }
  return title;
}

static void 
ReplaceAliasInAllDefinitionLines 
(SeqEntryPtr sep, 
 CharPtr alias, 
 CharPtr real_val)
{
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp;

  if (sep == NULL || StringHasNoText (alias) || StringHasNoText (real_val)) return;
    
  while (sep != NULL) {
    if (IS_Bioseq_set (sep))
    {
      sdp = SeqEntryGetSeqDescr(sep, Seq_descr_title, NULL);
      if (sdp != NULL)
      {
        sdp->data.ptrvalue = ReplaceOneModifierName (sdp->data.ptrvalue, 
                                                     alias, real_val);
                                                              
      }

      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL)
      {
        ReplaceAliasInAllDefinitionLines (bssp->seq_set, alias, real_val);
      }
    }
    else
    {
      sdp = SeqEntryGetSeqDescr(sep, Seq_descr_title, NULL);
      if (sdp != NULL)
      {
        sdp->data.ptrvalue = ReplaceOneModifierName (sdp->data.ptrvalue,
                                                     alias, real_val);
      }
    }
    sep = sep->next;
  }
}

static void ReplaceAllAliases (SeqEntryPtr sep)
{
  Int4 j;
  
  for (j = 0; j < num_aliases; j++)
  {
    ReplaceAliasInAllDefinitionLines (sep, alias_list[j].alias, alias_list[j].modifier);
  }
}

static CharPtr ReplaceMolNameWithMolBracketsInOneDefinitionLine (CharPtr title)
{
  CharPtr      ptr;
  
  if (StringHasNoText (title)) 
  {
    return title;
  }
   
  ptr = StringISearch (title, "[dna]"); 
  if (ptr != NULL)
  {
    ExciseString (title, "[dna", "]");
    TrimSpacesAroundString (title);
    title = ReplaceValueInOneDefLine (title, "molecule", "dna");
  }
  
  ptr = StringISearch (title, "[rna]");
  if (ptr != NULL)
  {
    ExciseString (title, "[rna", "]");
    TrimSpacesAroundString (title);
    title = ReplaceValueInOneDefLine (title, "molecule", "rna");
  }
  
  return title;
}

static void ReplaceMolNamesWithMolBracketsInDefinitionLines (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp;

  while (sep != NULL) 
  {
    if (IS_Bioseq_set (sep))
    {
      sdp = SeqEntryGetSeqDescr(sep, Seq_descr_title, NULL);
      if (sdp != NULL)
      {
        sdp->data.ptrvalue = ReplaceMolNameWithMolBracketsInOneDefinitionLine (sdp->data.ptrvalue);
      }

      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL)
      {
        ReplaceMolNamesWithMolBracketsInDefinitionLines (bssp->seq_set);
      }
    }
    else
    {
      sdp = SeqEntryGetSeqDescr(sep, Seq_descr_title, NULL);
      if (sdp != NULL)
      {
        sdp->data.ptrvalue = ReplaceMolNameWithMolBracketsInOneDefinitionLine (sdp->data.ptrvalue);
      }
    }
    sep = sep->next;
  }
}

static void ImportModifiersButtonProc (ButtoN b)
{
  SequencesFormPtr  sqfp;
  Boolean           rval;
  IDAndTitleEditPtr iatep;
  SeqEntryPtr       seq_list;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL) return;
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  rval = ImportModifiersToIDAndTitleEdit (iatep);
  if (rval)
  {
    ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
    SeqEntryPtrToSourceTab (sqfp);    
  }
  iatep = IDAndTitleEditFree (iatep);
}


typedef struct sourceassistant 
{
  CharPtr PNTR    defline_list;
  CharPtr PNTR    id_list;
  Int4            num_deflines;
  Int2            seqPackage;
  DialoG          mod_type_dlg;
  DoC             mod_doc;
  DialoG          orgmod_dlg;
  Boolean         done;
  Boolean         cancelled;
} SourceAssistantData, PNTR SourceAssistantPtr;

/* These functions are used for converting between a SourceAssistant structure and
 * and IDAndTitleEdit structure.
 */
static IDAndTitleEditPtr SourceAssistantToIDAndTitleEdit (SourceAssistantPtr sap)
{
  IDAndTitleEditPtr iatep;
  Int4              j;
  
  if (sap == NULL || sap->num_deflines < 1)
  {
    return NULL;
  }
  
  iatep = IDAndTitleEditNew ();
  if (iatep != NULL)
  {
    iatep->num_sequences = sap->num_deflines;
    iatep->id_list = (CharPtr PNTR) MemNew (iatep->num_sequences * sizeof (CharPtr));
    iatep->title_list = (CharPtr PNTR) MemNew (iatep->num_sequences * sizeof (CharPtr));
    for (j = 0; j < sap->num_deflines; j++)
    {
      iatep->id_list [j] = StringSave (sap->id_list [j]);
      iatep->title_list [j] = StringSave (sap->defline_list [j]);
    }
  }
  return iatep;
}

static void ApplyIDAndTitleEditToSourceAssistant (SourceAssistantPtr sap, IDAndTitleEditPtr iatep)
{
  Int4 seq_num;
  
  if (sap == NULL || iatep == NULL || sap->num_deflines != iatep->num_sequences)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    /* copy sequence IDs */
    sap->id_list [seq_num] = MemFree (sap->id_list [seq_num]);
    sap->id_list [seq_num] = StringSave (iatep->id_list [seq_num]);
    
    /* copy titles */
    sap->defline_list [seq_num] = MemFree (sap->defline_list [seq_num]);
    sap->defline_list [seq_num] = StringSave (iatep->title_list [seq_num]);
  }
}

static ValNodePtr PrepareSourceAssistantTableData (SourceAssistantPtr sap, BoolPtr multi_found)
{
  Int4               i;
  ValNodePtr         found_modifiers = NULL;
  ValNodePtr         vnp;
  ValNodePtr         column_list = NULL, row_list = NULL;
  ValNodePtr         header_list, header_vnp;
  Int4               column_width, num_columns = 0;
  CharPtr            org_name, mod_name, mod_value;
  Int4               max_column_width = 20;
  
  if (sap == NULL)
  {
    return NULL;
  }
  
  /* get list of modifiers */   
  /* location will be listed whether present or not */
  ValNodeAddPointer (&found_modifiers, 0, StringSave ("location")); 
  for (i = 0; i < sap->num_deflines; i++)
  {
    found_modifiers = BuildModifierTypeList (found_modifiers,
                                             sap->defline_list[i],
                                             FALSE);
  }
    
  /* create header line for table */
  /* store max column width in choice */
  ValNodeAddPointer (&column_list, 6, StringSave ("Seq ID"));
  ValNodeAddPointer (&column_list, 8, StringSave ("organism"));
  for (vnp = found_modifiers; vnp != NULL; vnp = vnp->next)
  {
    ValNodeAddPointer (&column_list, StringLen (vnp->data.ptrvalue), StringSave ((CharPtr) vnp->data.ptrvalue));
  } 
    
  ValNodeAddPointer (&row_list, 0, column_list);
  header_list = column_list;
  
  num_columns = ValNodeLen (column_list);

  /* create data lines for table */
  for (i = 0; i < sap->num_deflines; i++)
  {
    column_list = NULL;
    header_vnp = header_list;
    /* add Sequence ID */
    column_width = MAX (StringLen (sap->id_list[i]), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, StringSave (sap->id_list[i]));
      
    /* add organism name */
    header_vnp = header_vnp->next;
    org_name = GetDisplayValue ("organism", sap->defline_list[i], multi_found);
    column_width = MAX (StringLen (org_name), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, org_name);
            
    /* get remaining modifiers */
    for (vnp = found_modifiers; vnp != NULL; vnp = vnp->next)
    {
      header_vnp = header_vnp->next;
      mod_name = (CharPtr) vnp->data.ptrvalue;
      mod_value = GetDisplayValue (mod_name, sap->defline_list[i], multi_found);
      if (StringICmp (mod_name, "location") == 0 && StringHasNoText (mod_value))
      {
        /* display default value for location */
        mod_value = MemFree (mod_value);
        mod_value = StringSave ("genomic");
      }
      column_width = MAX (StringLen (mod_value), header_vnp->choice);
      column_width = MIN (column_width, max_column_width);
      header_vnp->choice = column_width;
      ValNodeAddPointer (&column_list, 0, mod_value);
    }
    ValNodeAddPointer (&row_list, 0, column_list);
  }
  ValNodeFreeData (found_modifiers);
  return row_list;
}


/* code for scientific name selection controls */

typedef struct organismselectiondialog
{
  DIALOG_MESSAGE_BLOCK
  TexT       tax_name_txt;
  DoC        org_list;
  Int4       org_row;
  CharPtr    tax_name_val;
} OrganismSelectionDialogData, PNTR OrganismSelectionDialogPtr;

static void CleanupOrganismSelectionDialog (GraphiC g, VoidPtr data)

{
  OrganismSelectionDialogPtr dlg;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (g);
  if (dlg != NULL)
  {
    dlg->tax_name_val = MemFree (dlg->tax_name_val);
  }

  StdCleanupExtraProc (g, data);
}

static Boolean OrgNameHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  OrganismSelectionDialogPtr dlg;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
  
  if (item == dlg->org_row) return TRUE;
  return FALSE;
}

static CharPtr GetTextForOrgPos (Int4 pos)
{
  ValNodePtr vnp;
  Int4       val;
  OrgInfoPtr oip;
  
  for (vnp = orglist, val = 1; vnp != NULL && val < pos; vnp = vnp->next, val++)
  {
  }
  if (vnp != NULL && vnp->data.ptrvalue != NULL)
  {
    oip = (OrgInfoPtr) vnp->data.ptrvalue;
  	return oip->taxname;;
  }
  else
  {
  	return NULL;
  }
}

static void GetOrgPosForText (CharPtr cp, Int4Ptr pos, Boolean PNTR match)
{
  ValNodePtr vnp;
  Int4       val = 1;
  CharPtr    dat;
  Int4       res;
  OrgInfoPtr oip;
  
  if (cp == NULL || pos == NULL || match == NULL) return;
  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    oip = (OrgInfoPtr) vnp->data.ptrvalue;
  	dat = oip->taxname;
  	res = StringCmp (cp, dat);
  	if (res < 0)
  	{
  	  *pos = val;
  	  *match = FALSE;
  	  return;
  	}
  	else if (res == 0)
  	{
  	  *pos = val;
  	  *match = TRUE;
  	  return;
  	}
  	val++;
  }
  *pos = val - 1;
  *match = FALSE;
}

static void OrgNameOnKey (SlatE s, Char ch)
{
  OrganismSelectionDialogPtr dlg;
  CharPtr                    str;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (s);
  if (dlg == NULL) return;

  if ( (int) ch == 0 ) return;
  
  /* later, handle control key combos */
#ifdef WIN_MSWIN
  if (ch == 3)
  {
    str = SaveStringFromText (dlg->tax_name_txt);
    StringToClipboard (str);
    str = MemFree (str);
  }
#else
  if (ctrlKey && ch == 'c')
  {
    str = SaveStringFromText (dlg->tax_name_txt);
    StringToClipboard (str);
    str = MemFree (str);
  }
#endif
}

static void SetOrganismText (TexT t)
{
  OrganismSelectionDialogPtr dlg;
  Int4                       pos, prevpos;
  Boolean                    match;
  CharPtr                    old_val;
  Boolean                    changed_val = FALSE;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (t);
  if (dlg == NULL) return;
  old_val = dlg->tax_name_val;
  dlg->tax_name_val = SaveStringFromText (dlg->tax_name_txt);
  if (dlg->tax_name_val != NULL)
  {
  	dlg->tax_name_val [0] = TO_UPPER (dlg->tax_name_val [0]);
  }
  if (!StringHasNoText (old_val) && StringCmp (old_val, dlg->tax_name_val) != 0)
  {
    changed_val = TRUE;
  }
  if (old_val != NULL)
  {
  	MemFree (old_val);
  }

  pos = -1;
  match = FALSE;
  GetOrgPosForText (dlg->tax_name_val, &pos, &match);
  SetOffset (dlg->org_list, 0, pos - 1);
  if (pos != dlg->org_row)
  {
    prevpos = dlg->org_row;
    if (match)
    { 
      dlg->org_row = pos;
      SetTitle (dlg->tax_name_txt, dlg->tax_name_val);
    }
    else
    {
      dlg->org_row = -1;
    }
  	if (prevpos != -1)
    {
  	  InvalDocRows (dlg->org_list, prevpos, 1, 1);
    }
    if (match)
    {
      InvalDocRows (dlg->org_list, dlg->org_row, 1, 1);
    }
  }
  else if (!match)
  {
  	dlg->org_row = -1;
    InvalDocRows (dlg->org_list, pos, 1, 1);	
  }
}

static void SetOrganismDoc (DoC d, PoinT pt)
{
  Int2      item, row, prevrow;

  OrganismSelectionDialogPtr dlg;
  CharPtr           old_name;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, NULL, NULL);
  if (item > 0 && row > 0) {
    prevrow = dlg->org_row;
    dlg->org_row = item;
    if (item != prevrow)
    {
      if (prevrow != -1)
      {
        InvalDocRows (d, prevrow, 1, 1);
      }
      InvalDocRows (d, item, 1, 1);
      old_name = SaveStringFromText (dlg->tax_name_txt);
      SetTitle (dlg->tax_name_txt, GetTextForOrgPos (item));
      old_name = MemFree (old_name);
      dlg->tax_name_val = SaveStringFromText (dlg->tax_name_txt);
    }  	
  }
}

static void DataToOrganismSelectionDialog (DialoG d, Pointer data)
{
  OrganismSelectionDialogPtr dlg;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  if (StringHasNoText (data))
  {
    SetTitle (dlg->tax_name_txt, "");
  }
  else
  {
    SetTitle (dlg->tax_name_txt, data);
  }
  SetOrganismText (dlg->tax_name_txt);  
}

static Pointer OrganismSelectionDialogToData (DialoG d)
{
  OrganismSelectionDialogPtr dlg;
  
  dlg = (OrganismSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  else
  {
    return SaveStringFromText (dlg->tax_name_txt);
  }
}

static ParData orgListPar = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData orgListCol = {0, 0, 160, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

extern DialoG OrganismSelectionDialog (GrouP parent, CharPtr org_name)
{
  GrouP           grp;
  Int2            height;
  ValNodePtr      vnp;
  OrganismSelectionDialogPtr dlg;
  RecT                       r;
  OrgInfoPtr                 oip;

  dlg = (OrganismSelectionDialogPtr) MemNew (sizeof (OrganismSelectionDialogData));

  grp = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (grp, dlg, CleanupOrganismSelectionDialog);
  SetGroupSpacing (grp, 10, 10);

  dlg->dialog = (DialoG) grp;
  dlg->todialog = DataToOrganismSelectionDialog;
  dlg->fromdialog = OrganismSelectionDialogToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;

  LoadOrganismList ();	
  
  dlg->tax_name_txt = DialogText (grp, "", 20, SetOrganismText);
  SetObjectExtra (dlg->tax_name_txt, dlg, NULL);
  dlg->org_row = -1;
  if (org_name != NULL)
  {
  	SetTitle (dlg->tax_name_txt, org_name);
  }
  SetOrganismText (dlg->tax_name_txt);

  SelectFont (programFont);
  height = LineHeight ();
  SelectFont (systemFont);
  dlg->org_list = DocumentPanel (grp, stdCharWidth * 25, height * 6);
  SetObjectExtra (dlg->org_list, dlg, NULL);
  
  ObjectRect (dlg->org_list, &r);
  InsetRect (&r, 4, 4);
  orgListCol.pixWidth = r.right - r.left;

  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
    oip = (OrgInfoPtr) vnp->data.ptrvalue;
    if (oip != NULL)
    {
  	  AppendText (dlg->org_list, oip->taxname, &orgListPar, &orgListCol, programFont);
    }
  }
  SetDocAutoAdjust (dlg->org_list, FALSE);
  SetDocProcs (dlg->org_list, SetOrganismDoc, NULL, NULL, NULL);
  SetDocShade (dlg->org_list, NULL, NULL, OrgNameHighlight, NULL);
  SetSlateChar ((SlatE) dlg->org_list, OrgNameOnKey);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->tax_name_txt, (HANDLE) dlg->org_list, NULL);
  InvalDocument (dlg->org_list);
  return (DialoG) grp;
}

#define NUM_ORGS_DISPLAYED 5
typedef struct multiorganismselectiondialog
{
  DIALOG_MESSAGE_BLOCK
  DialoG       org_select_dlg;
  TexT         tax_name_txt [NUM_ORGS_DISPLAYED];
  ButtoN       copy_btn [NUM_ORGS_DISPLAYED];
  PrompT       id_txt [NUM_ORGS_DISPLAYED];
  DialoG       location_dlg [NUM_ORGS_DISPLAYED];
  GrouP        gcode_grp [NUM_ORGS_DISPLAYED];
  ButtoN       gcode_btn [NUM_ORGS_DISPLAYED];
  DialoG       gcode_dlg [NUM_ORGS_DISPLAYED];
  BaR          id_scroll;
  ValNodePtr   row_list;
  Int4         num_vals;
  ValNodePtr   geneticcodelist;
} MultiOrganismSelectionDialogData, PNTR MultiOrganismSelectionDialogPtr;

typedef struct multiorgcopybtn
{
  MultiOrganismSelectionDialogPtr dlg;
  Int4                            pos;
} MultiOrgCopyBtnData, PNTR MultiOrgCopyBtnPtr;

static void CleanupMultiOrganismSelectionDialog (GraphiC g, VoidPtr data)

{
  MultiOrganismSelectionDialogPtr dlg;
  
  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (g);
  if (dlg != NULL)
  {
    dlg->row_list = FreeTableDisplayRowList (dlg->row_list);
    dlg->geneticcodelist = ValNodeFree (dlg->geneticcodelist);
  }

  StdCleanupExtraProc (g, data);
}

static CharPtr GetTableDisplayCellValue (ValNodePtr row_list, Int4 row_num, Int4 col_num)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       j;
  
  if (row_list == NULL)
  {
    return NULL;
  }
  
  for (row_vnp = row_list, j = 0;
       row_vnp != NULL && j < row_num; 
       row_vnp = row_vnp->next, j++)
  {
  }
  
  if (row_vnp == NULL)
  {
    return NULL;
  }
  
  for (col_vnp = row_vnp->data.ptrvalue, j = 0;
       col_vnp != NULL && j < col_num;
       col_vnp = col_vnp->next, j++)
  {
  }
  if (col_vnp == NULL)
  {
    return NULL;
  }
  else
  {
    return col_vnp->data.ptrvalue;
  }
}

static void 
UpdateTableDisplayCellValue 
(ValNodePtr row_list,
 Int4 row_num,
 Int4 col_num,
 CharPtr new_value)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       j;
  
  if (row_list == NULL)
  {
    return;
  }
  
  for (row_vnp = row_list, j = 0;
       row_vnp != NULL && j < row_num; 
       row_vnp = row_vnp->next, j++)
  {
  }
  
  if (row_vnp == NULL)
  {
    return;
  }
  
  for (col_vnp = row_vnp->data.ptrvalue, j = 0;
       col_vnp != NULL && j < col_num;
       col_vnp = col_vnp->next, j++)
  {
  }
  if (col_vnp == NULL)
  {
    return;
  }
  else
  {
    col_vnp->data.ptrvalue = StringSave (new_value);
  }
}

static void 
UpdateGeneticCodePosition 
(MultiOrganismSelectionDialogPtr dlg,
 Int4    row_num,
 CharPtr taxname,
 CharPtr location)
{
  ValNode    vn;
  Int4       gcode = -1;
  Int4       offset;
  CharPtr    gcode_name;
  
  if (dlg == NULL || row_num < 0 || row_num >= NUM_ORGS_DISPLAYED)
  {
    return;
  }

  gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
  
  if (gcode < 0)
  {
    offset = GetBarValue (dlg->id_scroll);
    gcode_name = GetTableDisplayCellValue (dlg->row_list, 
                                           offset + row_num, 3);
              
    vn.choice = GeneticCodeFromStringAndList (gcode_name, dlg->geneticcodelist);
    vn.next = NULL;
    vn.data.ptrvalue = gcode_name;
    PointerToDialog (dlg->gcode_dlg [row_num], &vn);
    Hide (dlg->gcode_btn [row_num]);
    Show (dlg->gcode_dlg [row_num]);
  }
  else
  {
    offset = GetBarValue (dlg->id_scroll);
    gcode_name = GeneticCodeStringFromIntAndList (gcode, dlg->geneticcodelist);
    UpdateTableDisplayCellValue (dlg->row_list, offset + row_num, 3, gcode_name);
    SetTitle (dlg->gcode_btn [row_num], gcode_name);
    Hide (dlg->gcode_dlg [row_num]);
    Show (dlg->gcode_btn [row_num]);
  }
}

static void DisplayPosition (MultiOrganismSelectionDialogPtr dlg, Int4 pos)
{
  Int4       row_num;
  ValNodePtr row_vnp, col_vnp;
  ValNode    vn;
  CharPtr    taxname = NULL;
  CharPtr    location = NULL;
  Int4       gcode;

  if (dlg == NULL)
  {
    return;
  }
  
  for (row_num = 0, row_vnp = dlg->row_list;
       row_num < pos && row_vnp != NULL;
       row_num++, row_vnp = row_vnp->next)
  {
  }
  
  for (row_num = 0;
       row_num < NUM_ORGS_DISPLAYED && row_vnp != NULL; 
       row_num++, row_vnp = row_vnp->next)
  {
    /* set ID */
    col_vnp = row_vnp->data.ptrvalue;
    SetTitle (dlg->id_txt [row_num], col_vnp->data.ptrvalue);
    /* set tax name */
    col_vnp = col_vnp->next;
    taxname = col_vnp->data.ptrvalue;
    SetTitle (dlg->tax_name_txt [row_num], taxname);
    /* set location */
    col_vnp = col_vnp->next;
    location = col_vnp->data.ptrvalue;
    vn.data.ptrvalue = NULL;
    if (StringHasNoText (location))
    {
      vn.choice = Source_location_genomic;
    }
    else
    {
      vn.choice = SrcLocFromGenome (GenomeFromLocName (location));
    }
    vn.next = NULL;
    PointerToDialog (dlg->location_dlg[row_num], &vn);
    
    Show (dlg->copy_btn [row_num]);
    Show (dlg->id_txt [row_num]);
    Show (dlg->tax_name_txt [row_num]);
    Show (dlg->location_dlg [row_num]);

    /* display genetic code */  
    col_vnp = col_vnp->next;   
    Show (dlg->gcode_grp [row_num]);
    
    gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
    if (gcode < 0)
    {
      vn.choice = GeneticCodeFromStringAndList (col_vnp->data.ptrvalue, dlg->geneticcodelist);
      vn.next = NULL;
      vn.data.ptrvalue = col_vnp->data.ptrvalue;
      PointerToDialog (dlg->gcode_dlg [row_num], &vn);
      Hide (dlg->gcode_btn [row_num]);
      Show (dlg->gcode_dlg [row_num]);
    }
    else
    {
      SetTitle (dlg->gcode_btn [row_num], GeneticCodeStringFromIntAndList (gcode, dlg->geneticcodelist));
      Hide (dlg->gcode_dlg [row_num]);
      Show (dlg->gcode_btn [row_num]);
    }
  }
  
  while (row_num < NUM_ORGS_DISPLAYED)
  {
    Hide (dlg->copy_btn [row_num]);
    Hide (dlg->id_txt [row_num]);
    Hide (dlg->tax_name_txt [row_num]);
    Hide (dlg->location_dlg [row_num]);
    Hide (dlg->gcode_grp [row_num]);
    row_num++;
  }
}

static void CollectPositionValues (MultiOrganismSelectionDialogPtr dlg, Int4 pos)
{
  Int4       row_num;
  ValNodePtr row_vnp, col_vnp, val_vnp;
  Int4       gcode;
  CharPtr    taxname, location, gcode_name;

  if (dlg == NULL)
  {
    return;
  }
  for (row_num = 0, row_vnp = dlg->row_list;
       row_num < pos && row_vnp != NULL;
       row_num++, row_vnp = row_vnp->next)
  {
  }
  
  for (row_num = 0;
       row_num < NUM_ORGS_DISPLAYED && row_vnp != NULL; 
       row_num++, row_vnp = row_vnp->next)
  {
    col_vnp = row_vnp->data.ptrvalue;
    /* skip ID - it can't be edited */
    col_vnp = col_vnp->next;
    
    /* get tax name */
    col_vnp->data.ptrvalue = MemFree (col_vnp->data.ptrvalue);
    col_vnp->data.ptrvalue = SaveStringFromText (dlg->tax_name_txt [row_num]);
    taxname = col_vnp->data.ptrvalue;
    col_vnp = col_vnp->next;
    
    /* get location */
    val_vnp = DialogToPointer (dlg->location_dlg [row_num]);
    if (val_vnp == NULL)
    {
      location = NULL;
    }
    else
    {
      location = val_vnp->data.ptrvalue;
    }
    StringToLower (location);
    val_vnp = ValNodeFree (val_vnp);
    col_vnp->data.ptrvalue = MemFree (col_vnp->data.ptrvalue);
    col_vnp->data.ptrvalue = location;
    col_vnp = col_vnp->next;
    
    /* get genetic code */
    gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
    if (gcode < 0)
    {
      val_vnp = DialogToPointer (dlg->gcode_dlg [row_num]);
      if (val_vnp == NULL)
      {
        gcode_name = NULL;
      }
      else
      {
        gcode_name = val_vnp->data.ptrvalue;
      }
      ValNodeFree (val_vnp);
    }
    else
    {
      gcode_name = StringSave (GeneticCodeStringFromIntAndList (gcode, dlg->geneticcodelist));
    }
    col_vnp->data.ptrvalue = MemFree (col_vnp->data.ptrvalue);
    col_vnp->data.ptrvalue = gcode_name;
    col_vnp = col_vnp->next;
    
  }  
}

static void MultiOrgScroll (BaR sb, GraphiC g, Int4 newval, Int4 oldval)
{
  MultiOrganismSelectionDialogPtr dlg;
  
  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (sb);
  if (dlg == NULL)
  {
    return;
  }
  
  /* first, collect old values */
  CollectPositionValues (dlg, oldval);
  
  /* set newly visible values */
  DisplayPosition (dlg, newval);  
}

static void MultiOrgCopy (ButtoN b)
{
  MultiOrgCopyBtnPtr mp;
  CharPtr            tax_name;
  ValNodePtr         val_vnp;
  CharPtr            location;

  mp = (MultiOrgCopyBtnPtr) GetObjectExtra (b);
  if (mp == NULL || mp->dlg == NULL || mp->pos < 0 || mp->pos >= NUM_ORGS_DISPLAYED)
  {
    return;
  }
  
  tax_name = (CharPtr) DialogToPointer (mp->dlg->org_select_dlg);
  SetTitle (mp->dlg->tax_name_txt [mp->pos], tax_name);
  
  /* get location for this row */
  val_vnp = DialogToPointer (mp->dlg->location_dlg [mp->pos]);
  if (val_vnp != NULL)
  {
    location = val_vnp->data.ptrvalue;
    UpdateGeneticCodePosition (mp->dlg, mp->pos, tax_name, location);
      
    val_vnp = ValNodeFreeData (val_vnp);
    /* location is freed when we free val_vnp */
    location = NULL;
  }
  
  tax_name = MemFree (tax_name);
}

static void DataToMultiOrganismSelectionDialog (DialoG d, Pointer userdata)
{
  MultiOrganismSelectionDialogPtr dlg;
  ValNodePtr                      row_list;

  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  row_list = (ValNodePtr) userdata;
  dlg->row_list = FreeTableDisplayRowList (dlg->row_list);
  dlg->row_list = CopyTableDisplayRowList (row_list);
  dlg->num_vals = ValNodeLen (dlg->row_list);
  
  CorrectBarMax (dlg->id_scroll, dlg->num_vals - NUM_ORGS_DISPLAYED);  
  CorrectBarValue (dlg->id_scroll, 0);
  DisplayPosition (dlg, 0);    
}

static Pointer MultiOrganismSelectionDialogToData (DialoG d)
{
  MultiOrganismSelectionDialogPtr dlg;
  Int4                            pos;
  ValNodePtr                      row_list;

  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return NULL;
  }
  
  pos = GetBarValue (dlg->id_scroll);
  CollectPositionValues (dlg, pos);
  row_list = CopyTableDisplayRowList (dlg->row_list);
  return row_list;
}

static void SetRowListColumn (ValNodePtr row_list, Int4 column, CharPtr new_value)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       col_num;
  
  for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
  {
    for (col_vnp = row_vnp->data.ptrvalue, col_num = 0;
         col_vnp != NULL && col_num < column;
         col_vnp = col_vnp->next, col_num++)
    {
    }
    if (col_vnp != NULL)
    {
      col_vnp->data.ptrvalue = MemFree (col_vnp->data.ptrvalue);
      col_vnp->data.ptrvalue = StringSave (new_value);
    }
  }
}

static void 
ApplyOrgModColumnOrCell 
(CharPtr            mod_name,
 CharPtr            suggested_value,
 Int4               row, 
 SourceAssistantPtr sap,
 SeqEntryPtr        seq_list,
 ValNodePtr         row_list,
 Int4               row_list_column,
 Int2               seqPackage);
 
static Boolean 
ContinueWithAutopopulatedGeneticCodes 
(SeqEntryPtr        seq_list,
 SourceAssistantPtr sap,
 ValNodePtr         row_list,
 Int4               affected_row);

static void SetAllGeneticCodes (ButtoN b)
{
  MultiOrganismSelectionDialogPtr dlg;
  Int4                            scroll_pos;

  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  scroll_pos = GetBarValue (dlg->id_scroll);
  CollectPositionValues (dlg, scroll_pos);

  if (ContinueWithAutopopulatedGeneticCodes (NULL, NULL, dlg->row_list, -1))
  {
    ApplyOrgModColumnOrCell ("genetic_code", "Standard", -1, NULL, NULL, dlg->row_list, 3, 0);
    DisplayPosition (dlg, scroll_pos);
  }

}

static void AddGcodeCommentBtn (ButtoN b)
{
  MultiOrgCopyBtnPtr bp;
  Int4               scroll_pos;
  CharPtr            orig_val = NULL;

  bp = (MultiOrgCopyBtnPtr) GetObjectExtra (b);
  if (bp == NULL || bp->dlg == NULL)
  {
    return;
  }

  scroll_pos = GetBarValue (bp->dlg->id_scroll);
  /* first, collect current values from dialog*/
  CollectPositionValues (bp->dlg, scroll_pos);

  orig_val = GetTableDisplayCellValue (bp->dlg->row_list, bp->pos, 4);
  ApplyOrgModColumnOrCell ("gencode_comment", orig_val, bp->pos, NULL, NULL,
                           bp->dlg->row_list, 4, 0);
  /* now repopulate */
  DisplayPosition (bp->dlg, scroll_pos);
}

static void SetAllLocations (ButtoN b)
{
  MultiOrganismSelectionDialogPtr dlg;
  Int4                            scroll_pos;

  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  
  scroll_pos = GetBarValue (dlg->id_scroll);
  /* first, collect current values from dialog*/
  CollectPositionValues (dlg, scroll_pos);
  
  ApplyOrgModColumnOrCell ("location", "genomic", -1, NULL, NULL, dlg->row_list, 2, 0);
  DisplayPosition (dlg, scroll_pos);
}

static void SetAllOrganisms (ButtoN b)
{
  MultiOrganismSelectionDialogPtr dlg;
  Int4                            scroll_pos;

  dlg = (MultiOrganismSelectionDialogPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  scroll_pos = GetBarValue (dlg->id_scroll);
  /* first, collect current values from dialog*/
  CollectPositionValues (dlg, scroll_pos);

  ApplyOrgModColumnOrCell ("organism", NULL, -1, NULL, NULL, dlg->row_list, 1, 0);
  DisplayPosition (dlg, scroll_pos);
}

static void ChangeLocationOrTaxName (MultiOrgCopyBtnPtr bp)
{
  ValNodePtr           val_vnp;
  CharPtr              tax_name = NULL, location = NULL;

  if (bp == NULL)
  {
    return;
  }

  /* get taxname for this row */
  tax_name = SaveStringFromText (bp->dlg->tax_name_txt [bp->pos]);

  /* get location for this row */
  val_vnp = DialogToPointer (bp->dlg->location_dlg [bp->pos]);
  if (val_vnp->data.ptrvalue == NULL)
  {
    UpdateGeneticCodePosition (bp->dlg, bp->pos, tax_name, NULL);
  }
  else
  {
    location = val_vnp->data.ptrvalue;
    UpdateGeneticCodePosition (bp->dlg, bp->pos, tax_name, location);
    
    val_vnp = ValNodeFreeData (val_vnp);
    /* location is freed when we free val_vnp */
    location = NULL;
  }
  
  tax_name = MemFree (tax_name); 
}

static void ChangeLocationPopup (Pointer userdata)
{
  MultiOrgCopyBtnPtr   bp;
  
  bp = (MultiOrgCopyBtnPtr) userdata;
  if (bp == NULL)
  {
    return;
  }

  ChangeLocationOrTaxName (bp);
}

static void MultiOrgText (TexT t)
{
  MultiOrgCopyBtnPtr bp;
  CharPtr            cp;

  bp = (MultiOrgCopyBtnPtr) GetObjectExtra (t);
  if (bp == NULL)
  {
    return;
  }
  cp = SaveStringFromText (t);
  PointerToDialog (bp->dlg->org_select_dlg, cp);
  cp = MemFree (cp);
  ChangeLocationOrTaxName (bp);
}

static void ChangeGeneticCodePopup (Pointer userdata)
{
  MultiOrgCopyBtnPtr   bp;
  ValNodePtr           val_vnp;
  CharPtr              gcode_name;
  Int4                 offset;

  bp = (MultiOrgCopyBtnPtr) userdata;
  if (bp == NULL)
  {
    return;
  }
  
  val_vnp = DialogToPointer (bp->dlg->gcode_dlg [bp->pos]);
  if (val_vnp == NULL)
  {
    gcode_name = NULL;
  }
  else
  {
    gcode_name = val_vnp->data.ptrvalue;
  }
  
  offset = GetBarValue (bp->dlg->id_scroll);

  UpdateTableDisplayCellValue (bp->dlg->row_list, bp->pos + offset, 3, gcode_name);
  
  val_vnp = ValNodeFreeData (val_vnp);
  
}

static DialoG MultiOrganismSelectionDialog (GrouP parent)
{
  MultiOrganismSelectionDialogPtr dlg;
  GrouP                           grp, id_grp, scroll_grp;
  Int4                            k;
  MultiOrgCopyBtnPtr              bp;
  RecT                            r1, r2, r3;
  ValNodePtr                      gencodelist;
  ButtoN                          b;
  PrompT                          p1;
#ifdef WIN_MAC
  Int2                            wid = 12;
#else
  Int2                            wid = 20;
#endif

  dlg = (MultiOrganismSelectionDialogPtr) MemNew (sizeof (MultiOrganismSelectionDialogData));

  grp = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (grp, dlg, CleanupMultiOrganismSelectionDialog);
  SetGroupSpacing (grp, 10, 10);

  dlg->dialog = (DialoG) grp;
  dlg->todialog = DataToMultiOrganismSelectionDialog;
  dlg->fromdialog = MultiOrganismSelectionDialogToData;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  
  dlg->num_vals = 0;
  dlg->row_list = NULL;
  
  dlg->geneticcodelist = GetGeneticCodeValNodeList ();

  dlg->org_select_dlg = OrganismSelectionDialog (grp, "");
  p1 = StaticPrompt (grp, "You can use the Copy buttons to populate the organism field from the selector above.",
                     0, 0, programFont, 'l');
  scroll_grp = NormalGroup (grp, 2, 0, "", NULL, NULL);
  
  id_grp = HiddenGroup (scroll_grp, 5, 0, NULL);
  SetGroupSpacing (id_grp, 10, 10);
  StaticPrompt (id_grp, "SeqID", 7 * stdCharWidth, 0, programFont, 'l');
  StaticPrompt (id_grp, "Copy", 0, 0, programFont, 'l');
  b = PushButton (id_grp, "Organism", SetAllOrganisms);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (id_grp, "Location", SetAllLocations);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (id_grp, "Genetic Code", SetAllGeneticCodes);
  SetObjectExtra (b, dlg, NULL);
  
  for (k = 0; k < NUM_ORGS_DISPLAYED; k++)
  {
    /* prompt with sequence ID */
     dlg->id_txt [k] = StaticPrompt (id_grp, "", 7 * stdCharWidth, 0, programFont, 'l');
    /* button for copying from organism selector */
    dlg->copy_btn [k] = PushButton (id_grp, "->", MultiOrgCopy);
    bp = (MultiOrgCopyBtnPtr) MemNew (sizeof (MultiOrgCopyBtnData));
    if (bp != NULL)
    {
      bp->dlg = dlg;
      bp->pos = k;
    }
    SetObjectExtra (dlg->copy_btn [k], bp, StdCleanupExtraProc);
    dlg->tax_name_txt [k] = DialogText (id_grp, "", wid, MultiOrgText);
    SetTextSelect (dlg->tax_name_txt [k], MultiOrgText, NULL);
    SetObjectExtra (dlg->tax_name_txt [k], bp, NULL);
    
    dlg->location_dlg [k] = ValNodeSelectionDialogExEx (id_grp, 
                                                 GetLocListForBioSource (NULL), 6,
                                                 ValNodeStringName,
                                                 ValNodeSimpleDataFree,
                                                 ValNodeStringCopy,
                                                 ValNodeChoiceMatch,
                                                 "location code",
                                                  ChangeLocationPopup, bp, FALSE,
                                                  FALSE, TRUE, NULL);
                                             
    dlg->gcode_grp [k] = HiddenGroup (id_grp, 0, 0, NULL);
    dlg->gcode_btn [k] = PushButton (dlg->gcode_grp [k], 
                                     "                                                ", 
                                     AddGcodeCommentBtn);
    SetObjectExtra (dlg->gcode_btn [k], bp, NULL);
    Hide (dlg->gcode_btn [k]);
    /* NOTE - need separate list because ValNodeSelectionDialog will free this one */
    gencodelist = GetGeneticCodeValNodeList ();
    dlg->gcode_dlg [k] = ValNodeSelectionDialog (dlg->gcode_grp [k], gencodelist, 6,
                                                 ValNodeStringName,
                                                 ValNodeSimpleDataFree,
                                                 ValNodeStringCopy,
                                                 ValNodeChoiceMatch,
                                                 "genetic code",
                                                  ChangeGeneticCodePopup, bp, FALSE);
    Hide (dlg->gcode_dlg [k]);
    AlignObjects (ALIGN_LEFT, (HANDLE) dlg->gcode_btn [k], (HANDLE) dlg->gcode_dlg [k], NULL);
  }
 
  dlg->id_scroll = ScrollBar4 (scroll_grp, 0, 10, MultiOrgScroll);
  SetObjectExtra (dlg->id_scroll, dlg, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->org_select_dlg, (HANDLE) p1, (HANDLE) scroll_grp, NULL);
  
  ObjectRect (dlg->copy_btn [0], &r1);
  ObjectRect (dlg->copy_btn [4], &r2);
  ObjectRect (dlg->id_scroll, &r3);
  r3.top = r1.top;
  r3.bottom = r2.bottom;
  SetPosition (dlg->id_scroll, &r3);
  
  return (DialoG) grp;
}


static void SourceAssistantExport (ButtoN b)
{
  SourceAssistantPtr sap;
  ValNodePtr         row_list = NULL;
  FILE               *fp;
  Char               path [PATH_MAX];
  
  sap = (SourceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL) return;

  if (! GetOutputFileName (path, sizeof (path), NULL)) return;
  fp = FileOpen (path, "w");
  if (fp == NULL) 
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }
  
  row_list = PrepareSourceAssistantTableData (sap, NULL);

  PrintTableDisplayRowListToFile (row_list, fp);
  row_list = FreeTableDisplayRowList (row_list);
  FileClose (fp);
}

static void SourceAssistantOk (ButtoN b)
{
  SourceAssistantPtr sap;
  
  sap = (SourceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL) return;
  sap->cancelled = FALSE;
  sap->done = TRUE;
}

static void SourceAssistantCancel (ButtoN b)
{
  SourceAssistantPtr sap;
  
  sap = (SourceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL) return;
  if (Message (MSG_YN, "You will lose your changes if you cancel.  Are you sure?")
      == ANS_NO)
  {
    return;
  }
  
  sap->cancelled = TRUE;
  sap->done = TRUE;
}

static CharPtr GetFirstDeflineValue (SourceAssistantPtr sap, CharPtr mod_name)
{
  Int4         i;
  CharPtr      valstr = NULL;

  if (sap == NULL)
  {
    return NULL;
  }
  
  for (i = 0; i < sap->num_deflines && valstr == NULL; i++)
  {
    valstr = FindValueFromPairInDefline (mod_name, sap->defline_list [i]);
    if (StringHasNoText (valstr))
    {
      valstr = MemFree (valstr);
    }
  }
  return valstr;
}

static CharPtr 
GetTagListValueEx (TagListPtr tlp, Int4 seq_num, Int4 col_num)
{
  Int4         seq_pos;
  CharPtr      str = NULL;
  ValNodePtr   vnp;
  
  if (tlp == NULL) return NULL;
  
  for (vnp = tlp->vnp, seq_pos = 0;
       vnp != NULL && seq_pos != seq_num;
       vnp = vnp->next, seq_pos++)
  {
    
  }
  if (vnp != NULL)
  {
    str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, col_num);
  }
  return str;
}

static void SetTagListValue (TagListPtr tlp, Int4 row, Int4 column, CharPtr new_value)
{
  ValNodePtr vnp;
  Int4       row_num;
  CharPtr    new_val;
  
  if (tlp == NULL)
  {
    return;
  }
  
  for (vnp = tlp->vnp, row_num = 0;
       vnp != NULL && row_num != row;
       vnp = vnp->next, row_num++)
  {
  }
  if (vnp == NULL)
  {
    return;
  }
  
  new_val = ReplaceTagListColumn (vnp->data.ptrvalue, new_value, column);
  if (new_val != vnp->data.ptrvalue)
  {
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    vnp->data.ptrvalue = new_val;
  }
}

static void SetTagListColumnValue (TagListPtr tlp, Int4 column, CharPtr new_value)
{
  ValNodePtr vnp;
  CharPtr    new_val;
  
  for (vnp = tlp->vnp;
       vnp != NULL;
       vnp = vnp->next)
  {  
    new_val = ReplaceTagListColumn (vnp->data.ptrvalue, new_value, column);
    if (new_val != vnp->data.ptrvalue)
    {
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = new_val;
    }
  }
}

static void UpdateOrgModDlg (SourceAssistantPtr sap)
{
  Int4               j;
  Boolean            found_organism = FALSE;
  Boolean            multi_found = FALSE;
  ValNodePtr         row_list, header_list = NULL;

  if (sap == NULL)
  {
    return;
  }
  row_list = PrepareSourceAssistantTableData (sap, &multi_found);
  PointerToDialog (sap->orgmod_dlg, row_list);
  
  if (row_list == NULL)
  {
    SetModifierList (sap->mod_doc, NULL);
    return;
  }
  
  for (j = 0; j < sap->num_deflines && !found_organism; j++)
  {
    if (FindValuePairInDefLine ("organism", sap->defline_list[j], NULL) != NULL)
    {
      found_organism = TRUE;
    }
  }
  header_list = row_list->data.ptrvalue;
  header_list = header_list->next;
  if (!found_organism)
  {
    header_list = header_list->next;
  }
  SetModifierList (sap->mod_doc, header_list);
  FreeTableDisplayRowList (row_list);
  
  if (multi_found)
  {
    AppendText (sap->mod_doc, multival_explanation, NULL, NULL, programFont); 
  }
}

static EnumFieldAssocPtr BuildGeneticCodeEnum (void)
{
  ValNodePtr        gencodelist = NULL;
  Int4              num_gencodes = 0, index;
  EnumFieldAssocPtr gencode_alist = NULL;
  ValNodePtr        vnp;

  gencodelist = GetGeneticCodeValNodeList ();
  num_gencodes = ValNodeLen (gencodelist);
  gencode_alist = (EnumFieldAssocPtr) MemNew ((num_gencodes + 2) * sizeof (EnumFieldAssoc));

  gencode_alist [0].name = StringSave (" ");
  gencode_alist [0].value = 0;

  for (index = 1, vnp = gencodelist;
       index <= num_gencodes && vnp != NULL; 
       index++, vnp = vnp->next)
  {
    gencode_alist [index].name = StringSave (vnp->data.ptrvalue);
    gencode_alist [index].value = vnp->choice;
  }

  gencode_alist [index].name = NULL;
  ValNodeFreeData (gencodelist);
  return gencode_alist;
}

static EnumFieldAssocPtr FreeGeneticCodeEnum (EnumFieldAssocPtr gcode_alist)
{
  EnumFieldAssocPtr eap;
  
  for (eap = gcode_alist; eap != NULL && eap->name != NULL; eap++)
  {
    eap->name = MemFree (eap->name);
  }
  gcode_alist = MemFree (gcode_alist);
  return gcode_alist;
}

static GrouP MakeSourceInstructionGroup (GrouP parent)
{
  GrouP instr_grp;
  
  instr_grp = HiddenGroup (parent, 1, 0, NULL);
  StaticPrompt (instr_grp, "Scientific names should not be abbreviated.",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "(use 'Drosophila melanogaster' instead of 'D. melanogaster')",
                0, 0, programFont, 'l');
  return instr_grp;
}

static GrouP MakeGeneticCodeInstructionGroup (GrouP parent)
{
  GrouP instr_grp;
  
  instr_grp = HiddenGroup (parent, 1, 0, NULL);
  StaticPrompt (instr_grp, "Please choose the translation table for your sequence.",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "Examples: Standard, Bacterial, Archaeal and Plant Plastid, Vertebrate Mitochondrial",
                0, 0, programFont, 'l');
  
  return instr_grp;
}

static GrouP MakeGeneticCodeCommentInstructionGroup (GrouP parent)
{
  GrouP instr_grp;
  
  instr_grp = HiddenGroup (parent, 1, 0, NULL);
  StaticPrompt (instr_grp, "When a genetic code is determined automatically from the organism name and location,",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "you cannot edit the genetic code directly.  You may provide an alternate genetic code",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "and the evidence to support it here and it will be reviewed by the GenBank staff.",
                0, 0, programFont, 'l');
  
  return instr_grp;
}

static GrouP MakeNontextInstructionGroup (GrouP parent)
{
  GrouP instr_grp;
  
  instr_grp = HiddenGroup (parent, 1, 0, NULL);
  StaticPrompt (instr_grp, "This modifier allows only TRUE/FALSE values.",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "The modifier will only appear in the file if it is set to TRUE,",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "otherwise it will be absent.",
                0, 0, programFont, 'l');
  
  return instr_grp;
}

static GrouP MakeLocationInstructionGroup (GrouP parent)
{
  GrouP instr_grp;
  
  instr_grp = HiddenGroup (parent, 1, 0, NULL);
  StaticPrompt (instr_grp, "Use this to specify the subcellular location or viral origin of the sequences.",
                0, 0, programFont, 'l');
  StaticPrompt (instr_grp, "Example: Use 'Genomic' for a sequence encoded by a nuclear gene.",
                0, 0, programFont, 'l');
  
  return instr_grp;
}

static GrouP MakeInstructionGroup (GrouP parent, Boolean is_nontext, Int2 mod_type)
{
  GrouP instr_grp = NULL;
    
  if (is_nontext)
  {
    instr_grp = MakeNontextInstructionGroup (parent);
  }
  else if (mod_type == eModifierType_Location)
  {
    instr_grp = MakeLocationInstructionGroup (parent);
  }
  else if (mod_type == eModifierType_Organism)
  {
    instr_grp = MakeSourceInstructionGroup (parent);
  }
  else if (mod_type == eModifierType_NucGeneticCode
           || mod_type == eModifierType_MitoGeneticCode
           || mod_type == eModifierType_GeneticCode)
  {
    instr_grp = MakeGeneticCodeInstructionGroup (parent);
  }
  else if (mod_type == eModifierType_GeneticCodeComment)
  {
    instr_grp = MakeGeneticCodeCommentInstructionGroup (parent);
  }
  return instr_grp;
}


/* This section of code prepares a dialog for editing one value for
 * the specified modifier type.  It can be used for setting the value
 * for a single sequence or the values for all sequences.
 */
typedef struct singlemodvaldlg 
{
  DIALOG_MESSAGE_BLOCK
  Boolean is_nontext;
  Int2    mod_type;
  Int2    seqPackage;
  
  PopuP   nontext_popup;
  DialoG  strvalue_dlg;
  DialoG  org_dlg;
  TexT    text_txt;  
} SingleModValDlgData, PNTR SingleModValDlgPtr;

static void SingleModValToDialog (DialoG d, Pointer userdata)
{
  SingleModValDlgPtr dlg;
  CharPtr            suggested_value;
  ValNode            vn;
  ValNodePtr         gencodelist;
  
  dlg = (SingleModValDlgPtr) GetObjectExtra (d);
  
  if (dlg == NULL)
  {
    return;
  }
  
  suggested_value = (CharPtr) userdata;
  
  if (dlg->is_nontext)
  {
    if (StringICmp (suggested_value, "TRUE") == 0)
    {
      SetValue (dlg->nontext_popup, 2);
    }
    else
    {
      SetValue (dlg->nontext_popup, 1);
    }
  }
  else if (dlg->mod_type == eModifierType_Location)
  {
    if (StringHasNoText (suggested_value))
    {
      vn.choice = 1;
      vn.data.ptrvalue = "genomic";
    }
    else
    {
      vn.choice = GetValForEnumName (biosource_genome_simple_alist, suggested_value);
      vn.data.ptrvalue = suggested_value;
    }
    vn.next = NULL;
    PointerToDialog (dlg->strvalue_dlg, &vn);
  }
  else if (dlg->mod_type == eModifierType_Origin)
  {
    vn.choice = GetValForEnumName (biosource_origin_alist, suggested_value);
    vn.data.ptrvalue = suggested_value;
    vn.next = NULL;
    PointerToDialog (dlg->strvalue_dlg, &vn);
  }
  else if (dlg->mod_type == eModifierType_Organism)
  {
    PointerToDialog (dlg->org_dlg, suggested_value);
  }
  else if (dlg->mod_type == eModifierType_NucGeneticCode
           || dlg->mod_type == eModifierType_MitoGeneticCode
           || dlg->mod_type == eModifierType_GeneticCode)
  {
    gencodelist = GetGeneticCodeValNodeList ();
    if (StringHasNoText (suggested_value))
    {
      vn.choice = 0;
    }
    else if (isdigit (suggested_value[0]))
    {
      vn.choice = atoi (suggested_value);
    }
    else
    {
      vn.choice = GeneticCodeFromStringAndList (suggested_value, gencodelist);
    }
    vn.next = NULL;
    vn.data.ptrvalue = suggested_value;
    PointerToDialog (dlg->strvalue_dlg, &vn);
    gencodelist = ValNodeFreeData (gencodelist);
  }
  else if (dlg->mod_type == eModifierType_MolType)
  {
    if (StringHasNoText (suggested_value))
    {
      vn.choice = 253;
    }
    else if (isdigit (suggested_value[0]))
    {
      vn.choice = atoi (suggested_value);
    }
    else
    {
      vn.choice = MolTypeFromString (suggested_value);
    }
    vn.next = NULL;
    vn.data.ptrvalue = suggested_value;
    PointerToDialog (dlg->strvalue_dlg, &vn);
  }
  else if (dlg->mod_type == eModifierType_Molecule)
  {
    if (StringICmp (suggested_value, "dna") == 0)
    {
      vn.choice = Seq_mol_dna;
    }
    else if (StringICmp (suggested_value, "rna") == 0)
    {
      vn.choice = Seq_mol_rna;
    }
    else
    {
      vn.choice = Seq_mol_dna;
    }
  }
  else if (dlg->mod_type == eModifierType_Topology)
  {
    if (StringHasNoText (suggested_value))
    {
      vn.choice = 1;
    }
    else if (isdigit (suggested_value[0]))
    {
      vn.choice = atoi (suggested_value);
    }
    else
    {
      vn.choice = TopologyFromString (suggested_value);
    }
    vn.next = NULL;
    vn.data.ptrvalue = suggested_value;
    PointerToDialog (dlg->strvalue_dlg, &vn);
  }
  else
  {
    if (StringHasNoText (suggested_value))
    {
      SetTitle (dlg->text_txt, "");
    }
    else
    {
      SetTitle (dlg->text_txt, suggested_value);
    }
  }
}

static Pointer DialogToSingleModVal (DialoG d)
{
  SingleModValDlgPtr dlg;
  CharPtr            new_value = NULL;
  ValNodePtr         value_vnp;
  
  dlg = (SingleModValDlgPtr) GetObjectExtra (d);
  
  if (dlg == NULL)
  {
    return NULL;
  }
  
  /* prepare value */
  if (dlg->is_nontext)
  {
    if (GetValue (dlg->nontext_popup) == 2)
    {
      new_value = StringSave ("2");
    }
  }
  else if (dlg->mod_type == eModifierType_Location
           || dlg->mod_type == eModifierType_Origin
           || dlg->mod_type == eModifierType_NucGeneticCode
           || dlg->mod_type == eModifierType_MitoGeneticCode
           || dlg->mod_type == eModifierType_GeneticCode
           || dlg->mod_type == eModifierType_MolType
           || dlg->mod_type == eModifierType_Molecule
           || dlg->mod_type == eModifierType_Topology)
  {
    value_vnp = DialogToPointer (dlg->strvalue_dlg);
    new_value = value_vnp->data.ptrvalue;
    if (dlg->mod_type == eModifierType_Location
        || dlg->mod_type == eModifierType_Origin)
    {
      StringToLower (new_value);
    }
    value_vnp = ValNodeFree (value_vnp);
    if (dlg->mod_type == eModifierType_MolType
        && StringICmp (new_value, "mRNA [cDNA]") == 0)
    {
      new_value = MemFree (new_value);
      new_value = StringSave ("mRNA");
    }
  }
  else if (dlg->mod_type == eModifierType_Organism)
  {
    new_value = DialogToPointer (dlg->org_dlg);
  }
  else
  {
    new_value = SaveStringFromText (dlg->text_txt);
  }
  return (Pointer) new_value;
}


static DialoG SingleModValDialog (GrouP parent, Boolean is_nontext, Int2 mod_type, Int2 seqPackage)
{
  SingleModValDlgPtr dlg;
  GrouP              grp;
  ValNodePtr         gencodelist;
  
  dlg = (SingleModValDlgPtr) MemNew (sizeof (SingleModValDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  grp = HiddenGroup (parent, 1, 0, NULL);
  SetObjectExtra (grp, dlg, StdCleanupExtraProc);
  SetGroupSpacing (grp, 10, 10);

  dlg->dialog = (DialoG) grp;
  dlg->todialog = SingleModValToDialog;
  dlg->fromdialog = DialogToSingleModVal;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;

  dlg->is_nontext = is_nontext;
  dlg->mod_type = mod_type;
  dlg->nontext_popup = NULL;
  dlg->strvalue_dlg = NULL;
  dlg->org_dlg = NULL;
  dlg->text_txt = NULL;  
    

  if (dlg->is_nontext)
  {
    dlg->nontext_popup = PopupList (grp, TRUE, NULL);
    PopupItem (dlg->nontext_popup, "FALSE");
    PopupItem (dlg->nontext_popup, "TRUE");
  }
  else if (dlg->mod_type == eModifierType_Location)
  {
    dlg->strvalue_dlg = ValNodeSelectionDialogExEx (grp, 
                                                 GetLocListForBioSource (NULL), 6,
                                                 ValNodeStringName,
                                                 ValNodeSimpleDataFree,
                                                 ValNodeStringCopy,
                                                 ValNodeChoiceMatch,
                                                 "location code",
                                                  NULL, NULL, FALSE,
                                                  FALSE, TRUE, NULL);

  }
  else if (dlg->mod_type == eModifierType_Origin)
  {
    dlg->strvalue_dlg = EnumAssocSelectionDialog (grp, biosource_origin_alist,
                                                  "origin", FALSE, NULL, NULL);
  }
  else if (dlg->mod_type == eModifierType_Organism)
  {
    dlg->org_dlg = OrganismSelectionDialog (grp, "");
  }
  else if (dlg->mod_type == eModifierType_NucGeneticCode
           || dlg->mod_type == eModifierType_MitoGeneticCode
           || dlg->mod_type == eModifierType_GeneticCode)
  {
    gencodelist = GetGeneticCodeValNodeList ();
    dlg->strvalue_dlg = ValNodeSelectionDialog (grp, gencodelist, 6,
                                        ValNodeStringName,
                                        ValNodeSimpleDataFree,
                                        ValNodeStringCopy,
                                        ValNodeChoiceMatch,
                                        "genetic code",
                                        NULL, NULL, FALSE);
  }
  else if (dlg->mod_type == eModifierType_MolType)
  {
    if (seqPackage == SEQ_PKG_GENOMICCDNA)
    {
      dlg->strvalue_dlg = EnumAssocSelectionDialog (grp, biomol_nucGen_alist,
                                               "moltype", FALSE, NULL, NULL);
    }
    else
    {
      dlg->strvalue_dlg = EnumAssocSelectionDialog (grp, biomol_nucX_alist,
                                               "moltype", FALSE, NULL, NULL);
    }
  }
  else if (dlg->mod_type == eModifierType_Molecule)
  {
    dlg->strvalue_dlg = EnumAssocSelectionDialog (grp, molecule_alist, 
                                                  "molecule", FALSE, NULL, NULL);
  }
  else if (dlg->mod_type == eModifierType_Topology)
  {
    dlg->strvalue_dlg = EnumAssocSelectionDialog (grp, topology_nuc_alist,
                                               "topology", FALSE, NULL, NULL);
  }
  else
  {
    if (dlg->mod_type == eModifierType_GeneticCodeComment)
    {
      dlg->text_txt = DialogText (grp, "", 40, NULL);
    }
    else
    {
      dlg->text_txt = DialogText (grp, "", 20, NULL);
    }
  }
  
  return (DialoG) grp;
}

static void AddSeqIDAndValueToRowList 
(CharPtr         id,
 CharPtr         title,
 ValNodePtr PNTR row_list)
{
  CharPtr      str = NULL;
  ValNodeBlock column_list;
  ValNodeBlock row_block;
  CharPtr      org_loc, org_end = NULL, next_org = NULL;

  if (row_list == NULL)
  {
    return;
  }
  InitValNodeBlock (&column_list, NULL);
  InitValNodeBlock (&row_block, *row_list);

  /* put ID in first location */
  ValNodeAddPointerToEnd (&column_list, 0, StringSave (id));

  /* get organism */
  org_loc = FindValuePairInDefLine ("organism", title, &org_end);
  str = FindValueFromPairInDefline ("organism", title);
  ValNodeAddPointerToEnd (&column_list, 0, str);
  
  if (org_end != NULL)
  {
    next_org = FindValuePairInDefLine ("organism", org_end + 1, &org_end);
  }
  
  /* get location */
  str = FindValueFromPairInDeflineBeforeCharPtr ("location", title, next_org);
  if (StringHasNoText (str))
  {
    str = MemFree (str);
    str = StringSave ("genomic");
  }
  ValNodeAddPointerToEnd (&column_list, 0, str);
  
  /* get genetic code */
  str = FindValueFromPairInDeflineBeforeCharPtr ("genetic_code", title, next_org);
  if (StringHasNoText (str))
  {
    str = MemFree (str);
    str = StringSave ("Standard");
  }
  ValNodeAddPointerToEnd (&column_list, 0, str);
  
  /* get genetic code comment */
  str = FindValueFromPairInDeflineBeforeCharPtr ("gencode_comment", title, next_org);
  ValNodeAddPointerToEnd (&column_list, 0, str);
  
  ValNodeAddPointerToEnd (&row_block, 0, column_list.head); 
  
  while (next_org != NULL)
  {
    InitValNodeBlock (&column_list, NULL);
    /* put blank ID in first location */
    ValNodeAddPointerToEnd (&column_list, 0, StringSave (""));
    
    /* get organism */
    org_loc = FindValuePairInDefLine ("organism", next_org, &org_end);
    str = FindValueFromPairInDefline ("organism", next_org);
    ValNodeAddPointerToEnd (&column_list, 0, str);
  
    next_org = FindValuePairInDefLine ("organism", org_end + 1, &org_end);
  
    /* get location */
    str = FindValueFromPairInDeflineBeforeCharPtr ("location", org_loc, next_org);
    if (StringHasNoText (str))
    {
      str = MemFree (str);
      str = StringSave ("genomic");
    }
    ValNodeAddPointerToEnd (&column_list, 0, str);
  
    /* get genetic code */
    str = FindValueFromPairInDeflineBeforeCharPtr ("genetic_code", org_loc, next_org);
    if (StringHasNoText (str))
    {
      str = MemFree (str);
      str = StringSave ("Standard");
    }
    ValNodeAddPointerToEnd (&column_list, 0, str);
  
    /* get genetic code comment */
    str = FindValueFromPairInDeflineBeforeCharPtr ("gencode_comment", org_loc, next_org);
    ValNodeAddPointerToEnd (&column_list, 0, str);
  
    ValNodeAddPointerToEnd (&row_block, 0, column_list.head); 
  }
  *row_list = row_block.head;
}

static CharPtr FindNthOrgPair (CharPtr title, Int4 org_num, CharPtr PNTR p_org_end)
{
  return FindNthValuePairInDefLine (title, "organism", org_num, p_org_end);
}

static void ApplyRowListToIDAndTitleEdit (ValNodePtr row_list, IDAndTitleEditPtr iatep)
{
  ValNodePtr row_vnp, col_vnp;
  CharPtr    last_id = NULL, id_txt;
  Int4       j, seq_num;
  Int4       org_num = 0;
  CharPtr    org_loc, org_end, next_org;
  
  if (row_list == NULL || iatep == NULL)
  {
    return;
  }
  
  for (row_vnp = row_list;
       row_vnp != NULL; 
       row_vnp = row_vnp->next)
  {
    col_vnp = row_vnp->data.ptrvalue;
    seq_num = -1;
    /* read sequence ID */
    if (col_vnp != NULL)
    {
      id_txt = col_vnp->data.ptrvalue;
      if (StringHasNoText (id_txt))
      {
        id_txt = last_id;
        org_num++;
      }
      else
      {
        last_id = id_txt;
        org_num = 0;
      }
      for (j = 0; j < iatep->num_sequences && seq_num == -1; j++)
      {
        if (StringCmp (iatep->id_list [j], id_txt) == 0)
        {
          seq_num = j;
        }
      }
      col_vnp = col_vnp->next;
    }
    
    if (seq_num < 0 || seq_num > iatep->num_sequences)
    {
      continue;
    }
    
    /* find organism name # org_num, find next_org and make sure
     * all values are added before it.
     */
    org_loc = FindNthOrgPair (iatep->title_list [seq_num], org_num, &org_end);
    
    if (org_end == NULL)
    {
      next_org = NULL;
    }
    else
    {
      next_org = FindValuePairInDefLine ("organism", org_end + 1, NULL);
    }
    
    /* add tax name */
    if (col_vnp != NULL)
    {
      if (org_loc == NULL)
      {
        iatep->title_list [seq_num] = ReplaceValueInOneDefLine (iatep->title_list [seq_num],
                                                                "organism",
                                                                col_vnp->data.ptrvalue);
      }
      else
      {
        iatep->title_list [seq_num] = ReplaceValueInThisValuePair (iatep->title_list [seq_num],
                                                                   org_loc, "organism", 
                                                                   org_end, 
                                                                   col_vnp->data.ptrvalue);

      }
      col_vnp = col_vnp->next;
    }
       
    /* add location */
    if (col_vnp != NULL)
    {
      org_loc = FindNthOrgPair (iatep->title_list [seq_num], org_num, &org_end);
      iatep->title_list [seq_num] = ReplaceValueInOneDefLineForOrganism (iatep->title_list [seq_num],
                                                                         "location",
                                                                         col_vnp->data.ptrvalue,
                                                                         org_loc);
      col_vnp = col_vnp->next;
    }
        
    /* add genetic code */
    if (col_vnp != NULL)
    {
      org_loc = FindNthOrgPair (iatep->title_list [seq_num], org_num, &org_end);
      iatep->title_list [seq_num] = ReplaceValueInOneDefLineForOrganism (iatep->title_list [seq_num],
                                                                         "genetic_code",
                                                                         col_vnp->data.ptrvalue,
                                                                         org_loc);
      col_vnp = col_vnp->next;
    }
    /* add genetic code comment */
    if (col_vnp != NULL)
    {
      org_loc = FindNthOrgPair (iatep->title_list [seq_num], org_num, &org_end);
      iatep->title_list [seq_num] = ReplaceValueInOneDefLineForOrganism (iatep->title_list [seq_num],
                                                                         "gencode_comment",
                                                                         col_vnp->data.ptrvalue,
                                                                         org_loc);
      col_vnp = col_vnp->next;
    }
  }
}


/* This function allows the user to edit the organisms, locations, and genetic codes for
 * all of the sequences in the set.
 */
static void EditOrganismColumn (SourceAssistantPtr sap, SeqEntryPtr seq_list)
{
  WindoW                 w;
  GrouP                  h, instr_grp, c;
  DialoG                 dlg;
  ValNodePtr             row_list;
  ModalAcceptCancelData  acd;
  ButtoN                 b;
  Int4                   j;
  IDAndTitleEditPtr      iatep;
  
  if (sap == NULL && seq_list == NULL)
  {
    return;
  }
  
  if (sap == NULL)
  {
    iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  }
  else
  {
    iatep = SourceAssistantToIDAndTitleEdit (sap);
  }
  
  if (iatep == NULL || iatep->num_sequences < 1)
  {
    iatep = IDAndTitleEditFree (iatep);
    return;
  }

  SendHelpScrollMessage (helpForm, "Organism Page", "Add Organisms, Locations, and Genetic Codes");
  
  if (iatep->num_sequences == 1)
  {
    w = MovableModalWindow (-20, -13, -10, -10, "Organism Editor", NULL);
  }
  else
  {
    w = MovableModalWindow (-20, -13, -10, -10, "Multiple Organism Editor", NULL);
  }
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  dlg = MultiOrganismSelectionDialog (h);
    
  instr_grp = MakeSourceInstructionGroup (h);
  
  row_list = NULL;
  
  for (j = 0; j < iatep->num_sequences; j++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [j])
    {
      continue;
    }
    AddSeqIDAndValueToRowList (iatep->id_list[j], iatep->title_list[j],
                                 &row_list);
  }
  
  PointerToDialog (dlg, row_list);
  row_list = FreeTableDisplayRowList (row_list);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg, (HANDLE) instr_grp,
                              (HANDLE) c, (HANDLE) NULL);

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
    Remove (w);
    return;
  }
  else
  {
    row_list = DialogToPointer (dlg);
    ApplyRowListToIDAndTitleEdit (row_list, iatep);
    if (sap == NULL)
    {
      ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
    }
    else
    {
      ApplyIDAndTitleEditToSourceAssistant (sap, iatep);
    }
    row_list = FreeTableDisplayRowList (row_list);
    UpdateOrgModDlg (sap);  
    Remove (w);
  }
  iatep = IDAndTitleEditFree (iatep);
}

typedef struct setcolumnvaluesdata
{
  TagListPtr tlp;
  DialoG     all_val_dlg;
  Boolean    is_nontext;
  Int2       mod_type;
} SetColumnValuesData, PNTR SetColumnValuesPtr;

static void SetAllColumnValues (SetColumnValuesPtr scvp, CharPtr new_value)
{
  CharPtr    taglist_str;
  Int4       j;
  ValNodePtr vnp;
  
  if (scvp == NULL || scvp->tlp == NULL)
  {
    return;
  }
  
  taglist_str = TagListStringFromDefLineValue (new_value, scvp->is_nontext, scvp->mod_type);
  SetTagListColumnValue (scvp->tlp, scvp->tlp->cols - 1, taglist_str);
  taglist_str = MemFree (taglist_str);

  SendMessageToDialog (scvp->tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = scvp->tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
  }
  scvp->tlp->max = MAX ((Int2) 0, (Int2) (j - scvp->tlp->rows));
  CorrectBarMax (scvp->tlp->bar, scvp->tlp->max);
  CorrectBarPage (scvp->tlp->bar, scvp->tlp->rows - 1, scvp->tlp->rows - 1);
  if (scvp->tlp->max > 0) {
    SafeShow (scvp->tlp->bar);
  } else {
    SafeHide (scvp->tlp->bar);
  }
  SendMessageToDialog (scvp->tlp->dialog, VIB_MSG_ENTER);
}

static void ClearColumnValues (ButtoN b)
{ 
  SetColumnValuesPtr scvp;
  
  scvp = (SetColumnValuesPtr) GetObjectExtra (b);
  
  if (ANS_NO == Message (MSG_YN, "Are you sure you want to clear all the values?"))
  {
    return;
  }

  SetAllColumnValues (scvp, NULL);  
}

static void EditOrgModApplyAll (ButtoN b)
{
  SetColumnValuesPtr scvp;
  CharPtr            new_value;
  
  scvp = (SetColumnValuesPtr) GetObjectExtra (b);
  if (scvp == NULL || scvp->all_val_dlg == NULL || scvp->tlp == NULL)
  {
    return;
  }
  
  if (ANS_NO == Message (MSG_YN, "Are you sure you want to set all the values?"))
  {
    return;
  }

  new_value = DialogToPointer (scvp->all_val_dlg);
  SetAllColumnValues (scvp, new_value);  
  new_value = MemFree (new_value);
}



/* The following functions are used for editing values for all of the sequences in 
 * a record for a specified modifier name.
 * Some modifiers can have multiple values, so for these we list the original value
 * that will be replaced by the new value.
 */
/* when id_str has no text, this creates one column,
 * otherwise it creates two columns, sequence ID and new value.
 * if there are multiple values, the sequence ID is only listed in the first row.
 */
static void 
AddRowsForModifierValuesForDefline 
(CharPtr         mod_name, 
 Boolean         is_nontext, 
 Int2            mod_type,
 CharPtr         id_str,
 CharPtr         defline,
 ValNodePtr PNTR list)
{
  Boolean added_any = FALSE;
  CharPtr valstr;
  CharPtr bracket_start, bracket_end;
  Int4    len;
  CharPtr taglist_str, tag_str;
  Boolean is_first = TRUE;
  
  if (list == NULL || StringHasNoText (mod_name))
  {
    return;
  }
  
  bracket_start = FindValuePairInDefLine (mod_name, defline, &bracket_end);
  while (bracket_start != NULL)
  {
    valstr = FindValueFromPairInDefline (mod_name, bracket_start);
    if (!StringHasNoText (valstr))
    {
      taglist_str = TagListStringFromDefLineValue (valstr, is_nontext, mod_type);

      len = StringLen (id_str) + StringLen (taglist_str) + 4;
      tag_str = (CharPtr) MemNew (len * sizeof (Char));
      if (tag_str != NULL)
      {
        if (StringHasNoText (id_str))
        {
          sprintf (tag_str, "%s\n", taglist_str);
        }
        else
        {
          if (is_first)
          {
            sprintf (tag_str, "%s\t%s\n", id_str, taglist_str);
            is_first = FALSE;
          }
          else
          {
            sprintf (tag_str, "\t%s\n", taglist_str);
          }
        }
        ValNodeAddPointer (list, 0, tag_str);
        added_any = TRUE;
      }
      taglist_str = MemFree (taglist_str);
    }
    valstr = MemFree (valstr);
    bracket_start = FindValuePairInDefLine (mod_name, bracket_end + 1, &bracket_end);
  }
  
  if (!added_any)
  {
    len = StringLen (id_str) + 4;
    tag_str = (CharPtr) MemNew (len * sizeof (Char));
    if (tag_str != NULL)
    {
      if (StringHasNoText (id_str))
      {
        sprintf (tag_str, "\n");
      }
      else
      {
        sprintf (tag_str, "%s\t\n", id_str);
      }
      ValNodeAddPointer (list, 0, tag_str);
    }
  }  
}

static ValNodePtr 
IDAndTitleEditToModifierColumnTagList 
(IDAndTitleEditPtr iatep,
 CharPtr           mod_name,
 Boolean           is_nontext,
 Int2              mod_type,
 Boolean           allow_multi)
{
  ValNodePtr list = NULL;
  Int4       seq_num;
  
  if (iatep == NULL)
  {
    return NULL;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    if (allow_multi)
    {
      AddRowsForModifierValuesForDefline (mod_name, is_nontext, mod_type, 
                                          iatep->id_list [seq_num],
                                          iatep->title_list [seq_num],
                                          &list);
    }
    else
    {
      AddSeqIDAndValueToTagList (iatep->id_list [seq_num], 
                                 iatep->title_list [seq_num],
                                 mod_name, &list);
    }
  }
  return list;  
}

static CharPtr 
GetValueFromTagListString 
(CharPtr           new_value, 
 Boolean           is_nontext, 
 Int4              mod_type,
 Int2              seqPackage,
 EnumFieldAssocPtr gencode_alist)
{
  CharPtr tmp_value;
  
  if (is_nontext)
  {
    if (StringCmp (new_value, "1") != 0)
    {
      new_value = MemFree (new_value);
    }
  }
  else if (mod_type == eModifierType_Location)
  {
    tmp_value = GetEnumName (atoi(new_value), biosource_genome_simple_alist);
    StringToLower (tmp_value);
    new_value = MemFree (new_value);
    new_value = StringSave (tmp_value);
  }
  else if (mod_type == eModifierType_Origin)
  {
    tmp_value = GetEnumName (atoi (new_value), biosource_origin_alist);
    StringToLower (tmp_value);
    new_value = MemFree (new_value);
    new_value = StringSave (tmp_value);
  }
  else if (mod_type == eModifierType_NucGeneticCode
           || mod_type == eModifierType_MitoGeneticCode)
  {
    tmp_value = GetEnumName (atoi(new_value), gencode_alist);
    new_value = MemFree (new_value);
    new_value = StringSave (tmp_value);
  }
  else if (mod_type == eModifierType_MolType)
  {
    if (seqPackage == SEQ_PKG_GENOMICCDNA)
    {
      tmp_value = GetEnumName (atoi (new_value), biomol_nucGen_alist);
    }
    else
    {
      tmp_value = GetEnumName (atoi (new_value), biomol_nucX_alist);
      if (StringICmp (tmp_value, "mRNA [cDNA]") == 0)
      {
        tmp_value = StringSave ("mRNA");
      }            
    }
    new_value = MemFree (new_value);
    new_value = StringSave (tmp_value);        
  }
  else if (mod_type == eModifierType_Molecule)
  {
    tmp_value = GetEnumName (atoi (new_value), molecule_alist);
    new_value = MemFree (new_value);
    new_value = StringSave (tmp_value);        
  }
  else if (mod_type == eModifierType_Topology)
  {
    tmp_value = GetEnumName (atoi (new_value), topology_nuc_alist);
    new_value = MemFree (new_value);
    new_value = StringSave (tmp_value);        
  }
  return new_value;
}

static void RemoveAllValuePairs (CharPtr value_name, CharPtr title)
{
  CharPtr value_loc;
  CharPtr end_loc;
  
  if (StringHasNoText (value_name) || StringHasNoText (title))
  {
    return;
  }
  
  value_loc = FindValuePairInDefLine (value_name, title, &end_loc);
  while (value_loc != NULL)
  {
    RemoveValuePairFromDefline (value_loc, end_loc, value_loc);
    value_loc = FindValuePairInDefLine (value_name, value_loc, &end_loc);
  }
}

static CharPtr 
ApplyValueListToTitle 
(CharPtr orig_title, 
 CharPtr value_name, 
 ValNodePtr value_list)
{
  Int4       val_num;
  CharPtr    value_loc, end_loc = NULL;
  ValNodePtr vnp;
  
  if (StringHasNoText (value_name))
  {
    return orig_title;
  }
  
  /* if our value list is NULL, remove all values and done. */
  if (value_list == NULL)
  {
    RemoveAllValuePairs (value_name, orig_title);
    return orig_title;
  }
  
  /* if there are no values in the title, make sure the new value is added
   * to the first organism.
   * otherwise, replace values where they are in the title.
   */
  value_loc = FindValuePairInDefLine (value_name, orig_title, &end_loc);
  if (value_loc == NULL)
  {
    orig_title = ReplaceValueInOneDefLineForOrganism (orig_title, value_name, 
                                                      value_list->data.ptrvalue,
                                                      NULL);
  }
  else
  {
    vnp = value_list;
    val_num = 0;
    while (value_loc != NULL && vnp != NULL)
    {
      orig_title = ReplaceValueInThisValuePair (orig_title, value_loc, value_name,
                                                end_loc, vnp->data.ptrvalue);
      /* if the value was empty, it will have been removed */                                                
      if (!StringHasNoText (vnp->data.ptrvalue))
      {
        val_num++;
      }
      vnp = vnp->next;
      value_loc = FindNthValuePairInDefLine (orig_title, value_name, val_num, &end_loc);                                                        
    }
  }
  
  return orig_title;
}

static Int4 GetRowForIDText (CharPtr id_txt, IDAndTitleEditPtr iatep)
{
  Int4 seq_num, row_num = -1;
  
  if (StringHasNoText (id_txt) || iatep == NULL)
  {
    return -1;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences && row_num == -1; seq_num++)
  {
    if (StringCmp (id_txt, iatep->id_list [seq_num]) == 0)
    {
      row_num = seq_num;
    }
  }
  return row_num;
}

static void 
ApplyModifierColumnTagListToIDAndTitleEdit 
(CharPtr           mod_name,
 Boolean           is_nontext,
 Int4              mod_type,
 Int2              seqPackage,
 EnumFieldAssocPtr gencode_alist,
 ValNodePtr        list,
 Int4              num_columns,
 IDAndTitleEditPtr iatep,
 Int4              seq_num)
{
  ValNodePtr vnp, value_list = NULL, row_list;
  CharPtr    id_txt, new_value, last_id = NULL;
  Int4       row_num;
  
  if (list == NULL || iatep == NULL)
  {
    return;
  }

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (seq_num < 0)
    {
      id_txt = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
      if (StringHasNoText (id_txt))
      {
        id_txt = MemFree (id_txt);
        id_txt = StringSave (last_id);
      }
      else
      {
        last_id = MemFree (last_id);
        last_id = StringSave (id_txt);
      }
      
      /* find sequence that corresponds to this id */
      row_num = GetRowForIDText (id_txt, iatep);
      id_txt = MemFree (id_txt);
    }
    else
    {
      row_num = seq_num;
    }
    if (row_num >= iatep->num_sequences || row_num < 0)
    {
      continue;
    }
    
    /* extract new value from column */    
    new_value = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, num_columns - 1);
    /* translate from list value (may be number from popup) to real value */
    new_value = GetValueFromTagListString (new_value, is_nontext, mod_type, 
                                           seqPackage, gencode_alist);
    
    /* add to list */
    /* add NULL if value is blank */
    if (StringHasNoText (new_value))
    {
      new_value = MemFree (new_value);
    }
    ValNodeAddPointer (&value_list, row_num, new_value);
  }

  if (seq_num < 0)
  {
    for (row_num = 0; row_num < iatep->num_sequences; row_num++)
    {
      row_list = ValNodeExtractList (&value_list, row_num);
      iatep->title_list [row_num] = ApplyValueListToTitle (iatep->title_list [row_num],
                                                         mod_name, row_list);
      row_list = ValNodeFreeData (row_list);                                                         
    }
  }
  else
  {
    iatep->title_list [seq_num] = ApplyValueListToTitle (iatep->title_list [seq_num],
                                                         mod_name, value_list);

  }
  value_list = ValNodeFreeData (value_list);

  if (seq_num < 0)
  {
    for (row_num = 0; row_num < iatep->num_sequences; row_num++)
    {
      iatep->title_list [row_num] = RemoveAllDuplicatePairsFromOneTitle (iatep->title_list [row_num]);
    }
  }
  else
  {
    iatep->title_list [seq_num] = RemoveAllDuplicatePairsFromOneTitle (iatep->title_list [seq_num]);
  }
}

static Int4 GetTaglistType (Boolean is_nontext, Int4 mod_type)
{
  if (is_nontext
      || mod_type == eModifierType_Location
      || mod_type == eModifierType_Origin
      || mod_type == eModifierType_NucGeneticCode
      || mod_type == eModifierType_MitoGeneticCode
      || mod_type == eModifierType_MolType
      || mod_type == eModifierType_Molecule
      || mod_type == eModifierType_Topology)
  {
    return TAGLIST_POPUP;
  }
  else
  {
    return TAGLIST_TEXT;
  }
}

static EnumFieldAssocPtr GetTaglistAlist (Boolean is_nontext, Int4 mod_type, Int2 seqPackage)
{
  if (is_nontext)
  {
    return nontextmodedit_alist;
  }
  else if (mod_type == eModifierType_Location)
  {
    return  biosource_genome_simple_alist;
  }
  else if (mod_type == eModifierType_Origin)
  {
    return  biosource_origin_alist;
  }
  else if (mod_type == eModifierType_NucGeneticCode
           || mod_type == eModifierType_MitoGeneticCode)
  {
    return BuildGeneticCodeEnum ();
  }
  else if (mod_type == eModifierType_MolType)
  {
    if (seqPackage == SEQ_PKG_GENOMICCDNA)
    {
      return biomol_nucGen_alist;
    }
    else
    {
      return biomol_nucX_alist;
    }
  }
  else if (mod_type == eModifierType_Molecule)
  {
    return  molecule_alist;    
  }
  else if (mod_type == eModifierType_Topology)
  {
    return  topology_nuc_alist;    
  }
  else
  {
    return NULL;
  }    
}

static Uint2 taglist_types [] = 
{ TAGLIST_PROMPT, TAGLIST_PROMPT};

static Uint2 taglist_textWidths [] =
{ 10, 20};
static EnumFieldAssocPtr taglist_alists [] =
{ NULL, NULL};

static Int4 GetMaxTagListValueWidth (ValNodePtr taglist_list, Int4 col_num)
{
  ValNodePtr vnp;
  Int4       max_len = 0;
  CharPtr    tmp_value;
  
  for (vnp = taglist_list; vnp != NULL; vnp = vnp->next)
  {
    tmp_value = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, col_num);
    max_len = MAX (max_len, (Int4) StringLen (tmp_value));
    tmp_value = MemFree (tmp_value);
  }
  return max_len;
}

static DialoG 
CreateValueListDialog 
(GrouP             parent_grp,
 CharPtr           mod_name, 
 Int2              seqPackage,
 IDAndTitleEditPtr iatep,
 Int4              seq_num)
{
  GrouP                  k, g;
  PrompT                 ppt1 = NULL, ppt2;
  Int4                   num_columns, j;
  ValNodePtr             row_list = NULL;
  Boolean                allow_multi;
  Boolean                is_nontext;
  Int4                   mod_type;
  Int4                   rows_shown;
  DialoG                 dlg = NULL;
  TagListPtr             tlp;
  Int4                   first_colwidth, val_width;
  
  if (iatep == NULL || seq_num >= iatep->num_sequences)
  {
    return NULL;
  }
  
  is_nontext = IsNonTextModifier (mod_name);
  mod_type = GetModifierType (mod_name);

  /* set up row list and number of columns */
  if (seq_num < 0)
  {
    allow_multi = AllowMultipleValues (mod_name);
    row_list = IDAndTitleEditToModifierColumnTagList (iatep, mod_name, 
                                                      is_nontext, mod_type, 
                                                      allow_multi);
                                                    
    num_columns = 2;
  }
  else
  {
    allow_multi = TRUE;
    num_columns = 1;
    AddRowsForModifierValuesForDefline (mod_name, is_nontext, mod_type, 
                                        NULL,
                                        iatep->title_list [seq_num],
                                        &row_list);    
  }
  
  k = HiddenGroup (parent_grp, 1, 0, NULL);
  g = HiddenGroup (k, 2, 0, NULL);
  
  if (seq_num < 0)
  {
    ppt1 = StaticPrompt (g, "SeqID", 0, 0, programFont, 'l');
  }
  if (mod_type == eModifierType_MolType)
  {
    ppt2 = StaticPrompt (g, "molecule type", 0, 0, programFont, 'l');
  }
  else
  {
    ppt2 = StaticPrompt (g, mod_name, 0, 0, programFont, 'l');
  }
  
  rows_shown = ValNodeLen (row_list);
  rows_shown = MIN (rows_shown, 5);
  
  /* calculate appropriate column widths */
  if (seq_num < 0)
  {
    first_colwidth = 5;
    for (j = 0; j < iatep->num_sequences; j++)
    {
      first_colwidth = MAX (first_colwidth, (Int4)StringLen (iatep->id_list [j]));
    }
    
    if (mod_type == eModifierType_MolType)
    {
      val_width = StringLen ("molecule type");
    }
    else
    {
      val_width = StringLen (mod_name);
    }
    val_width = MAX (val_width, GetMaxTagListValueWidth (row_list, num_columns - 1));
  }
  else
  {
    val_width = MAX (14, GetMaxTagListValueWidth (row_list, num_columns - 1));
    first_colwidth = val_width;
  }
  
  taglist_textWidths [0] = first_colwidth;
  taglist_textWidths [1] = val_width;
  
  taglist_types [0] = TAGLIST_PROMPT; /* sequence ID */
  taglist_types [1] = TAGLIST_PROMPT; 
  
  taglist_alists [0] = NULL;
  taglist_alists [1] = NULL;
  
  taglist_types [num_columns - 1] = GetTaglistType (is_nontext, mod_type);
  taglist_alists [num_columns - 1] = GetTaglistAlist (is_nontext, mod_type, seqPackage);
  
  dlg = CreateTagListDialogEx (k, rows_shown, num_columns, 2,
                               taglist_types, taglist_textWidths,
                               taglist_alists, TRUE, TRUE, NULL, NULL);
  
  
  tlp = (TagListPtr) GetObjectExtra (dlg);  
  if (tlp == NULL) return NULL;

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = row_list;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  tlp->max = MAX ((Int2) 0, (Int2) (ValNodeLen (row_list) - tlp->rows));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  if (tlp->max > 0) {
    SafeShow (tlp->bar);
  } else {
    SafeHide (tlp->bar);
  }
  SendMessageToDialog (tlp->dialog, VIB_MSG_ENTER);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control[0], (HANDLE) ppt1, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control[num_columns - 1], (HANDLE) ppt2, NULL);
  return dlg;
}

static void 
EditOrgModColumn 
(CharPtr            mod_name,
 SourceAssistantPtr sap, 
 SeqEntryPtr        seq_list,
 Int2               seqPackage)
{
  WindoW                 w;
  GrouP                  h, k, c;
  DialoG                 dlg;
  TagListPtr             tlp;
  ModalAcceptCancelData  acd;
  ButtoN                 b;
  Boolean                is_nontext;
  Int4                   mod_type;
  GrouP                  instr_grp = NULL;
  SetColumnValuesData    scvd;
  ButtoN                 clear_btn;
  CharPtr                mod_label;
  ButtoN                 apply_all_btn;
  ValNodePtr             row_list = NULL;
  Int4                   num_columns;
  IDAndTitleEditPtr      iatep;
  Boolean                allow_multi;
  
  if (StringHasNoText (mod_name) || (sap == NULL && seq_list == NULL))
  {
    return;
  }
  
  if (StringICmp (mod_name, "moltype") == 0)
  {
    mod_label = StringSave ("molecule type");
  }
  else
  {
    mod_label = StringSave (mod_name);
  }

  is_nontext = IsNonTextModifier (mod_name);
  mod_type = GetModifierType (mod_name);
  
  if (mod_type == eModifierType_Organism
      || mod_type == eModifierType_Location
      || mod_type == eModifierType_GeneticCode)
  {
    EditOrganismColumn (sap, seq_list);
    return;
  }
  
  allow_multi = AllowMultipleValues (mod_name);
  if (sap == NULL)
  {
    iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  }
  else
  {
    iatep = SourceAssistantToIDAndTitleEdit (sap);
  }
  
  row_list = IDAndTitleEditToModifierColumnTagList (iatep, mod_name, 
                                                    is_nontext, mod_type, 
                                                    allow_multi);
                                                    
  num_columns = 2;
  
  w = MovableModalWindow (-20, -13, -10, -10, mod_label, NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  mod_label = MemFree (mod_label);

  instr_grp = MakeInstructionGroup (h, is_nontext, mod_type);  
  
  scvd.all_val_dlg = SingleModValDialog (h, is_nontext, mod_type, seqPackage);
  PointerToDialog (scvd.all_val_dlg, NULL);
  scvd.is_nontext = is_nontext;
  scvd.mod_type = mod_type;
  apply_all_btn = PushButton (h, "Apply above value to all sequences", EditOrgModApplyAll);
  SetObjectExtra (apply_all_btn, &scvd, NULL);

  k = HiddenGroup (h, -1, 0, NULL);

  dlg = CreateValueListDialog (k, mod_name, seqPackage, iatep, -1);
  
  scvd.tlp = (TagListPtr) GetObjectExtra (dlg);
  
  if (mod_type == eModifierType_MolType)
  {
    clear_btn = PushButton (h, "Reset All to Genomic DNA", ClearColumnValues);
  }
  else if (mod_type == eModifierType_Topology)
  {
    clear_btn = PushButton (h, "Reset All to Linear", ClearColumnValues);
  }
  else if (mod_type == eModifierType_Molecule)
  {
    clear_btn = PushButton (h, "Reset All to DNA", ClearColumnValues);
  }
  else
  {
    clear_btn = PushButton (h, "Clear All Values", ClearColumnValues);
  }
  SetObjectExtra (clear_btn, &scvd, NULL);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) scvd.all_val_dlg,
                              (HANDLE) apply_all_btn,
                              (HANDLE) k, 
                              (HANDLE) clear_btn, 
                              (HANDLE) c, 
                              (HANDLE) instr_grp,
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
  if (! acd.cancelled)
  {
    tlp = GetObjectExtra (dlg);
    
    ApplyModifierColumnTagListToIDAndTitleEdit (mod_name, is_nontext,
                                                mod_type, seqPackage,
                                                taglist_alists [num_columns - 1], 
                                                tlp->vnp,
                                                num_columns, iatep, -1);
    if (sap == NULL)
    {
      ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
    }
    else
    {
      ApplyIDAndTitleEditToSourceAssistant (sap, iatep);
    }
    
    UpdateOrgModDlg (sap);  
  }
  Remove (w);

  if (mod_type == eModifierType_NucGeneticCode
      || mod_type == eModifierType_MitoGeneticCode)
  {
    taglist_alists [num_columns - 1] = FreeGeneticCodeEnum (taglist_alists [num_columns - 1]);
  }
  iatep = IDAndTitleEditFree (iatep);
}

static void 
EditModsForOneSequence 
(CharPtr            mod_name,
 SourceAssistantPtr sap, 
 SeqEntryPtr        seq_list,
 Int2               seqPackage,
 Int4               seq_num)
{
  WindoW                 w;
  GrouP                  h, k, c;
  DialoG                 dlg = NULL;
  TagListPtr             tlp;
  ModalAcceptCancelData  acd;
  ButtoN                 b;
  Boolean                is_nontext;
  CharPtr                new_value = NULL;
  Int4                   mod_type;
  GrouP                  instr_grp = NULL;
  CharPtr                mod_label;
  Int4                   num_columns = 1;
  IDAndTitleEditPtr      iatep;
  Boolean                allow_multi;
  DialoG                 all_val_dlg = NULL;
  
  if (StringHasNoText (mod_name) || (sap == NULL && seq_list == NULL))
  {
    return;
  }
  
  if (StringICmp (mod_name, "moltype") == 0)
  {
    mod_label = StringSave ("molecule type");
  }
  else
  {
    mod_label = StringSave (mod_name);
  }

  is_nontext = IsNonTextModifier (mod_name);
  mod_type = GetModifierType (mod_name);
  
  allow_multi = AllowMultipleValues (mod_name);
  if (sap == NULL)
  {
    iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  }
  else
  {
    iatep = SourceAssistantToIDAndTitleEdit (sap);
  }
  
  /* make sure sequence number is in range */
  if (seq_num > iatep->num_sequences || seq_num < 0)
  {
    iatep = IDAndTitleEditFree (iatep);
    return;
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, mod_label, NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  instr_grp = MakeInstructionGroup (h, is_nontext, mod_type);  

  k = HiddenGroup (h, -1, 0, NULL);
  
  if (allow_multi)
  {
    dlg = CreateValueListDialog (k, mod_name, seqPackage, iatep, seq_num);
  }
  else
  {
    all_val_dlg = SingleModValDialog (h, is_nontext, mod_type, seqPackage);
    /* set initial value */
    new_value = FindValueFromPairInDefline (mod_name, iatep->title_list [seq_num]);
    PointerToDialog (all_val_dlg, new_value);
    new_value = MemFree (new_value);
  }
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) k, 
                              (HANDLE) c, 
                              (HANDLE) instr_grp,
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
  if (! acd.cancelled)
  {
    if (allow_multi)
    {
      tlp = GetObjectExtra (dlg);
    
      ApplyModifierColumnTagListToIDAndTitleEdit (mod_name, is_nontext,
                                                  mod_type, seqPackage,
                                                  taglist_alists [num_columns - 1], 
                                                  tlp->vnp,
                                                  num_columns, iatep, seq_num);
    }
    else
    {
      new_value = DialogToPointer (all_val_dlg);
      iatep->title_list [seq_num] = ReplaceValueInOneDefLine (iatep->title_list [seq_num],
                                                              mod_name, new_value);
      new_value = MemFree (new_value);
    }
    if (sap == NULL)
    {
      ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
    }
    else
    {
      ApplyIDAndTitleEditToSourceAssistant (sap, iatep);
    }
    
    UpdateOrgModDlg (sap);  
  }
  Remove (w);

  if (mod_type == eModifierType_NucGeneticCode
      || mod_type == eModifierType_MitoGeneticCode)
  {
    taglist_alists [num_columns - 1] = FreeGeneticCodeEnum (taglist_alists [num_columns - 1]);
  }
  iatep = IDAndTitleEditFree (iatep);
}

static void UpdateGeneticCodesForIDAndTitleEdit (IDAndTitleEditPtr iatep)
{
  Int4 seq_num;
  CharPtr     taxname, location, gcode_name;
  Int4        gcode;
  ValNodePtr  gencodelist;
  
  if (iatep == NULL)
  {
    return;
  }
  gencodelist = GetGeneticCodeValNodeList ();

  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    taxname = FindValueFromPairInDefline ("organism", iatep->title_list [seq_num]);
    location = FindValueFromPairInDefline ("location", iatep->title_list [seq_num]);
    if (StringHasNoText (location))
    {
      location = MemFree (location);
      location = StringSave ("genomic");
    }
    gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
    taxname = MemFree (taxname);
    location = MemFree (location);
    if (gcode > 0)
    {
      gcode_name = GeneticCodeStringFromIntAndList (gcode, gencodelist);
      iatep->title_list [seq_num] = ReplaceValueInOneDefLine (iatep->title_list [seq_num], "genetic_code", gcode_name);
    }
  }
  gencodelist = ValNodeFreeData (gencodelist);  
}

static Boolean 
ContinueWithAutopopulatedGeneticCodes 
(SeqEntryPtr        seq_list,
 SourceAssistantPtr sap,
 ValNodePtr         row_list,
 Int4               affected_row)
{
  ValNodePtr  autopop_list = NULL, already_have = NULL;
  SeqEntryPtr sep, nuc_sep;
  BioseqPtr   bsp;
  Int4        j;
  CharPtr     list_msg = NULL;
  Int4        num_sequences = 0;
  Boolean     rval = TRUE;
  CharPtr     taxname, location, gcode_name;
  Int4        gcode;
  ValNodePtr  row_vnp, col_vnp;
  
  if (seq_list == NULL && sap == NULL && row_list == NULL)
  {
    return FALSE;
  }
  
  if (seq_list != NULL)
  {
    for (sep = seq_list, j = 0; sep != NULL; sep = sep->next, j++)
    {
      if (affected_row != -1 && affected_row != j)
      {
        continue;
      }
      bsp = NULL;
      if (IS_Bioseq (sep))
      {
        bsp = (BioseqPtr) sep->data.ptrvalue;
      }
      else if (IS_Bioseq_set (sep))
      {
        nuc_sep = FindNucSeqEntry (sep);
        if (nuc_sep != NULL && IS_Bioseq (nuc_sep))
        {
          bsp = (BioseqPtr) nuc_sep->data.ptrvalue;
        }
      }
      if (bsp == NULL)
      {
        continue;
      }
      
      taxname = GetModValueFromSeqEntry (sep, "organism");
      location = GetModValueFromSeqEntry (sep, "location");
      if (StringHasNoText (location))
      {
        location = MemFree (location);
        location = StringSave ("genomic");
      }
      gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
      taxname = MemFree (taxname);
      location = MemFree (location);
      if (gcode > 0)
      {
        if (bsp != NULL)
        {
          ValNodeAddPointer (&autopop_list, 0, SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT));
        }
      }
      else
      {
        gcode_name = GetModValueFromSeqEntry (sep, "genetic_code");
        if (!StringHasNoText (gcode_name))
        {
          ValNodeAddPointer (&already_have, 0, SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT));
        }
        gcode_name = MemFree (gcode_name);
      }
      num_sequences++;
    }
  }
  else if (sap != NULL)
  {
    for (j = 0; j < sap->num_deflines; j++)
    {
      if (affected_row != -1 && affected_row != j)
      {
        continue;
      }
      taxname = GetValueFromTitle ("organism", sap->defline_list [j]);
      location = GetValueFromTitle ("location", sap->defline_list [j]);
      gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
      taxname = MemFree (taxname);
      location = MemFree (location);
      if (gcode > 0)
      {
        ValNodeAddPointer (&autopop_list, 0, StringSave (sap->id_list[j]));
      }
      else
      {
        gcode_name = GetValueFromTitle ("genetic_code", sap->defline_list [j]);
        if (!StringHasNoText (gcode_name))
        {
          ValNodeAddPointer (&already_have, 0, StringSave (sap->id_list[j]));
        }
        gcode_name = MemFree (gcode_name);
      }
    }
    num_sequences = sap->num_deflines;
  }
  else if (row_list != NULL)
  {
    num_sequences = 0;
    for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
    {
      col_vnp = row_vnp->data.ptrvalue;
      if (col_vnp != NULL && col_vnp->next != NULL 
          && col_vnp->next->next != NULL
          && col_vnp->next->next->next != NULL)
      {
        taxname = col_vnp->next->data.ptrvalue;
        location = col_vnp->next->next->data.ptrvalue;
        gcode = GetGeneticCodeForTaxNameAndLocation (taxname, location);
        if (gcode > 0)
        {
          ValNodeAddPointer (&autopop_list, 0, StringSave (col_vnp->data.ptrvalue));
        }
      }
      num_sequences++;
    }
  }
  
  if (autopop_list != NULL)
  {
    list_msg = CreateListMessage ("Sequence", 
               (autopop_list->next == NULL ? 
               " has a genetic code determined by the location and scientific name.  The genetic code for this sequence cannot be edited."
               : " have genetic codes determined by the location and scientific name.  The genetic code for these sequences cannot be edited."),
               autopop_list);
    if (ValNodeLen (autopop_list) == num_sequences || affected_row != -1)
    {
      Message (MSG_ERROR, list_msg);
      rval = FALSE;
    }
    else
    {
      if (ANS_NO == Message (MSG_YN, 
                     "%s  Do you want to edit the genetic code for the remaining sequences?",
                     list_msg))
      {
        rval = FALSE;
      }
    }
    list_msg = MemFree (list_msg);
  }
  autopop_list = ValNodeFreeData (autopop_list);
  
  if (rval && already_have != NULL && affected_row == -1 && num_sequences > 1)
  {
    list_msg = CreateListMessage ("Sequence",
                     (already_have->next == NULL ?
                     " already has a genetic code.  Do you wish to overwrite this value?"
                     : "already have genetic codes.  Do you wish to overwrite these values?"),
                     already_have);                     
    if (ANS_NO == Message (MSG_YN, list_msg))
    {
      rval = FALSE;
    }
    list_msg = MemFree (list_msg);
  }
  return rval;
}

static Boolean 
IDAndTitleEditHasAllDefaultValues 
(Int2              mod_type, 
 CharPtr           mod_name,
 IDAndTitleEditPtr iatep)
{
  Boolean               orig_all_default = TRUE;
  Int4                  seq_num;
  CharPtr               mod_value;
  
  if (mod_type != eModifierType_MolType
      && mod_type != eModifierType_Topology
      && mod_type != eModifierType_Location)
  {
    return FALSE;
  }
  if (iatep == NULL)
  {
    return FALSE;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences && orig_all_default; seq_num++)
  {
    if (iatep->is_seg && iatep->is_seg [seq_num])
    {
      continue;
    }
    mod_value = FindValueFromPairInDefline (mod_name, iatep->title_list [seq_num]);
    if ((mod_type == eModifierType_MolType && StringICmp (mod_value, "Genomic DNA") != 0)
        || (mod_type == eModifierType_Topology && StringICmp (mod_value, "Linear") != 0)
        || (mod_type == eModifierType_Location && StringICmp (mod_value, "Genomic") != 0))
    {
      orig_all_default = FALSE;
    }
    mod_value = MemFree (mod_value);
  }
  return orig_all_default;
}

static Boolean 
RowListHasAllDefaultValues 
(Int2              mod_type, 
 CharPtr           mod_name,
 ValNodePtr        row_list,
 Int4              row_list_column)
{
  Boolean               orig_all_default = TRUE;
  Int4                  seq_num;
  CharPtr               mod_value;
  Int4                  num_rows;
  
  
  if (mod_type != eModifierType_MolType
      && mod_type != eModifierType_Topology
      && mod_type != eModifierType_Location)
  {
    return FALSE;
  }
  if (row_list == NULL)
  {
    return FALSE;
  }
  
  num_rows = ValNodeLen (row_list);

  for (seq_num = 0; seq_num < num_rows && orig_all_default; seq_num++)
  {
    mod_value = GetRowListCellText (row_list, seq_num, row_list_column);
    if ((mod_type == eModifierType_MolType && StringICmp (mod_value, "Genomic DNA") != 0)
        || (mod_type == eModifierType_Topology && StringICmp (mod_value, "Linear") != 0)
        || (mod_type == eModifierType_Location && StringICmp (mod_value, "Genomic") != 0))
    {
      orig_all_default = FALSE;
    }
    mod_value = MemFree (mod_value);
  }
  return orig_all_default;
}

static void 
ApplyOrgModColumnOrCell 
(CharPtr            mod_name,
 CharPtr            suggested_value,
 Int4               row, 
 SourceAssistantPtr sap,
 SeqEntryPtr        seq_list,
 ValNodePtr         row_list,
 Int4               row_list_column,
 Int2               seqPackage)
{
  WindoW                w;
  GrouP                 h, c;
  GrouP                 instr_grp = NULL;
  Int4                  j;
  ModalAcceptCancelData acd;
  ButtoN                b;
  Boolean               is_nontext;
  CharPtr               new_value = NULL;
  CharPtr               title;
  CharPtr               all_seq_fmt = "%s (all sequences)";
  CharPtr               one_seq_fmt = "%s (Seq_ID %s)";
  Int4                  num_sequences = 0;
  Char                  id_txt[128];
  Int2                  mod_type;
  ValNodePtr            row_vnp = NULL, col_vnp;
  Int4                  row_num, col_num;
  CharPtr               mod_label;
  Boolean               done;
  DialoG                val_dlg;
  PrompT                apply_to_all_prompt = NULL;
  Boolean               orig_all_default = FALSE;
  IDAndTitleEditPtr     iatep = NULL;
  Int4                  seq_num;
  
  if (StringHasNoText (mod_name) || row < -1)
  {
    return;
  }
  if (sap == NULL && seq_list == NULL && row_list == NULL)
  {
    return;
  }

  if (sap != NULL)
  {
    iatep = SourceAssistantToIDAndTitleEdit (sap);
  }
  else if (seq_list != NULL)
  {
    iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  }
  
  if (iatep != NULL)
  {
    num_sequences = 0;
    for (j = 0; j < iatep->num_sequences; j++)
    {
      if (iatep->is_seg != NULL && iatep->is_seg [j])
      {
        continue;
      }
      num_sequences++;
    }
  }
  else
  {
    num_sequences = ValNodeLen (row_list);
  }
    
  if (row >= num_sequences)
  {
    return;
  }
  
  is_nontext = IsNonTextModifier (mod_name);
  mod_type = GetModifierType (mod_name);
  
  if (row < 0)
  {
    if (iatep != NULL)
    {
      orig_all_default = IDAndTitleEditHasAllDefaultValues (mod_type, mod_name, iatep);
    }
    else
    {
      orig_all_default = RowListHasAllDefaultValues (mod_type, mod_name, row_list, row_list_column);
    }
  }

  /* get label to use in window */
  if (mod_type == eModifierType_MolType)
  {
    mod_label = StringSave ("molecule type");
  }
  else
  {
    mod_label = StringSave (mod_name);
  }

  if (row == -1)
  {
    title = (CharPtr) MemNew ((StringLen (mod_label) + StringLen (all_seq_fmt)) * sizeof (Char));
    sprintf (title, all_seq_fmt, mod_label);
  }
  else
  {
    if (iatep != NULL)
    {
      StringNCpy (id_txt, iatep->id_list [row], sizeof (id_txt) - 1);
      id_txt[sizeof(id_txt) - 1] = 0;
    }
    else if (row_list != NULL)
    {
      for (row_vnp = row_list, row_num = 0;
           row_vnp != NULL && row_num < row;
           row_vnp = row_vnp->next, row_num++)
      {
        
      }
      if (row_vnp != NULL)
      {
        col_vnp = row_vnp->data.ptrvalue;
        if (col_vnp != NULL)
        {
          StringNCpy (id_txt, col_vnp->data.ptrvalue, sizeof (id_txt) - 1);
          id_txt[sizeof(id_txt) - 1] = 0;
        }
        else
        {
          return;
        }
      }
      else
      {
        return;
      }
    }
    title = (CharPtr) MemNew ((StringLen (mod_label) 
                               + StringLen (one_seq_fmt)
                               + StringLen (id_txt)) * sizeof (Char));
    sprintf (title, one_seq_fmt, mod_label, id_txt);
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, title, NULL);
  title = MemFree (title);
  
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  if (row == -1)
  {
    apply_to_all_prompt = StaticPrompt (h, "This will apply to all sequences in the record.",
                                        0, 0, programFont, 'l');
  }
  
  instr_grp = MakeInstructionGroup (h, is_nontext, mod_type);
  
  val_dlg = SingleModValDialog (h, is_nontext, mod_type, seqPackage);
  PointerToDialog (val_dlg, suggested_value);
    
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  if (instr_grp == NULL)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) val_dlg,
                                (HANDLE) c, 
                                (HANDLE) apply_to_all_prompt, 
                                NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) val_dlg,
                                (HANDLE) instr_grp, 
                                (HANDLE) c, 
                                (HANDLE) apply_to_all_prompt,
                                NULL);
  }

  mod_label = MemFree (mod_label);

  Show (w);
  Select (w);
  
  done = FALSE;
  while (!done)
  {
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
      done = TRUE;
    }
    else if (row < 0 
             && ! orig_all_default
             && ANS_NO == Message (MSG_YN, "Are you sure you want to apply this value to all of your sequences?"))
    {
      /* do nothing - they'll be able to cancel from the dialog if they want to */
    }
    else
    {
      /* prepare value */
      new_value = DialogToPointer (val_dlg);
      
      /* apply values */
      if (iatep != NULL)
      {
        for (j = 0, seq_num = 0; j < iatep->num_sequences; j++)
        {
          /* don't apply modifiers to segments */
          if (iatep->is_seg != NULL && iatep->is_seg [j])
          {
            continue;
          }
          
          if (seq_num == row || row == -1)
          {
            iatep->title_list [j] = ReplaceValueInOneDefLine (iatep->title_list [j],
                                                              mod_name,
                                                              new_value);
          }
          seq_num++;
        }
        if (mod_type == eModifierType_Organism 
            || mod_type == eModifierType_Location
            || mod_type == eModifierType_GeneticCode)
        {
          UpdateGeneticCodesForIDAndTitleEdit (iatep);
        }
        if (seq_list != NULL)
        {
          ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
        }
        else if (sap != NULL)
        {
          ApplyIDAndTitleEditToSourceAssistant (sap, iatep);
        }
      }
      else if (row_list != NULL)
      {
        if (row < 0)
        {
          SetRowListColumn (row_list, row_list_column, new_value);
        }
        else if (row_vnp != NULL)
        {
          for (col_vnp = row_vnp->data.ptrvalue, col_num = 0;
               col_vnp != NULL && col_num < row_list_column;
               col_vnp = col_vnp->next, col_num++)
          {
          }
          if (col_vnp != NULL)
          {
            col_vnp->data.ptrvalue = MemFree (col_vnp->data.ptrvalue);
            col_vnp->data.ptrvalue = StringSave (new_value);
          }
        }
      }
      new_value = MemFree (new_value);
      UpdateOrgModDlg (sap);
      done = TRUE;
    } 
  }
  Remove (w);
  iatep = IDAndTitleEditFree (iatep);
}

static void 
ApplyOrgModColumn 
(CharPtr mod_name,
 CharPtr suggested_value,
 SourceAssistantPtr sap)
{
  Int4    mod_type;
  Boolean is_nontext;
  IDAndTitleEditPtr iatep;
  Boolean           all_default;
  
  mod_type = GetModifierType (mod_name);
  is_nontext = IsNonTextModifier (mod_name);
  
  if (mod_type == eModifierType_GeneticCode)
  {
    if (! ContinueWithAutopopulatedGeneticCodes (NULL, sap, NULL, -1))
    {
      return;
    }
  }
  else if (!StringHasNoText (suggested_value)  && ! IsNonTextModifier (mod_name)
      && sap->num_deflines > 1)
  {
    iatep = SourceAssistantToIDAndTitleEdit (sap);
    all_default = IDAndTitleEditHasAllDefaultValues (mod_type, mod_name, iatep);
    iatep = IDAndTitleEditFree (iatep);
  
    if (!all_default 
        && ANS_YES != Message (MSG_YN, "Warning!  Some sequences already contain "
                            "a value for %s.  Are you sure you want to "
                            "overwrite these values?", mod_name))
    {
      return;
    }
  }
  ApplyOrgModColumnOrCell (mod_name, suggested_value, -1, sap, NULL, NULL, 0, sap->seqPackage);
}

static void ApplyEditOrgModColumnBtn (ButtoN b, Boolean apply)
{
  SourceAssistantPtr sap;
  CharPtr            mod_name = NULL;
  CharPtr            suggested_value;
  ValNodePtr         vnp;
  SourceQualDescPtr  sqdp;

  sap = (SourceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL) return;
  vnp = DialogToPointer (sap->mod_type_dlg);
  if (vnp != NULL)
  {
    if (vnp->choice == eModifierType_Organism)
    {
      mod_name = "organism";
    }
    else if (vnp->choice == eModifierType_Location)
    {
      mod_name = "location";
    }
    else if (vnp->choice == 0 && vnp->data.ptrvalue != NULL)
    {
      sqdp = (SourceQualDescPtr) vnp->data.ptrvalue;
      mod_name = sqdp->name;
    }
  }
  if (!StringHasNoText (mod_name)) {
    if (apply)
    {
      if (IsNonTextModifier (mod_name))
      {
        suggested_value = StringSave ("TRUE");
      }
      else
      {
        suggested_value = GetFirstDeflineValue (sap, mod_name);
      }
      ApplyOrgModColumn (mod_name, suggested_value, sap);
      suggested_value = MemFree (suggested_value);
    }
    else
    {
      EditOrgModColumn (mod_name, sap, NULL, sap->seqPackage);
    }
  }
  vnp = ValNodeFreeData (vnp);
}

static void ApplyOrgModColBtn (ButtoN b)
{
  ApplyEditOrgModColumnBtn (b, TRUE);
}

static void EditOrgModColBtn (ButtoN b)
{
  ApplyEditOrgModColumnBtn (b, FALSE);
}

static void SourceAssistantImportModsBtn (ButtoN b)
{
  SourceAssistantPtr sap;
  IDAndTitleEditPtr  iatep;
  Boolean            rval;

  sap = (SourceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL) return;

  iatep = SourceAssistantToIDAndTitleEdit (sap);
  rval = ImportModifiersToIDAndTitleEdit (iatep);
  if (rval)
  {
    ApplyIDAndTitleEditToSourceAssistant (sap, iatep);
    UpdateOrgModDlg (sap);
  }
  iatep = IDAndTitleEditFree (iatep);

}

static void SourceAssistantClearAllModifiers (ButtoN b)
{
  SourceAssistantPtr sap;
  Int4               j;
  ValNodePtr         found_modifiers = NULL, vnp;

  sap = (SourceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL) return;
  
  if (Message (MSG_YN, 
      "Are you sure you want to remove all of the source qualifiers from all of your sequences?")
      == ANS_NO)
  {
    return;
  }

  for (j = 0; j < sap->num_deflines; j++)
  {
    found_modifiers = BuildModifierTypeList (found_modifiers, 
                                             sap->defline_list[j],
                                             FALSE);
  }
  
  for (j = 0; j < sap->num_deflines; j++)
  {
    for (vnp = found_modifiers; vnp != NULL; vnp = vnp->next)
    { 
      if (StringICmp (vnp->data.ptrvalue, "genetic_code") == 0
          || StringICmp (vnp->data.ptrvalue, "organism") == 0
          || StringICmp (vnp->data.ptrvalue, "location") == 0
          || StringICmp (vnp->data.ptrvalue, "gencode_comment") == 0
          || StringICmp (vnp->data.ptrvalue, "moltype") == 0
          || StringICmp (vnp->data.ptrvalue, "topology") == 0)
      {
        continue;
      }
      RemoveValueFromDefline (vnp->data.ptrvalue, sap->defline_list [j]);
    }
  }
  
  found_modifiers = ValNodeFreeData (found_modifiers);
  UpdateOrgModDlg (sap);
}

static void OrgModDblClick (PoinT cell_coord, CharPtr header_text, CharPtr cell_text, Pointer userdata)
{
  SourceAssistantPtr sap;
  Int4               mod_type;
  
  sap = (SourceAssistantPtr) userdata;
  if (sap == NULL)
  {
    return;
  }
  
  mod_type = GetModifierType (header_text);
  
  if (cell_coord.x == 0)
  {
    return;
  }
  else if (cell_coord.y == 0)
  {
    EditOrgModColumn (cell_text, sap, NULL, sap->seqPackage);
  }
  else if (mod_type != eModifierType_GeneticCode 
           || ContinueWithAutopopulatedGeneticCodes (NULL, sap, NULL, cell_coord.y - 1))
  {
    EditModsForOneSequence (header_text, sap, NULL, sap->seqPackage, cell_coord.y - 1);
  }
}

static void SeqEntryListToSourceAssistant (SeqEntryPtr seq_list, SourceAssistantPtr sap)
{
  IDAndTitleEditPtr iatep;
  Int4              seq_num, sap_num;
  
  if (sap == NULL)
  {
    return;
  }
  sap->num_deflines = 0;
  
  iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  if (iatep == NULL || iatep->num_sequences < 1)
  {
    iatep = IDAndTitleEditFree (iatep);
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    else
    {
      sap->num_deflines ++;
    }
  }
  
  if (sap->num_deflines < 1)
  {
    iatep = IDAndTitleEditFree (iatep);
    return;
  }
  
  sap->defline_list = (CharPtr PNTR) MemNew (sizeof (CharPtr) * sap->num_deflines);
  sap->id_list = (CharPtr PNTR) MemNew (sizeof (CharPtr) * sap->num_deflines);
  if (sap->defline_list == NULL || sap->id_list == NULL)
  {
    sap->defline_list = MemFree (sap->defline_list);
    sap->id_list = MemFree (sap->id_list);
    return;
  }
  
  sap_num = 0;
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    sap->id_list [sap_num] = StringSave (iatep->id_list [seq_num]);
    sap->defline_list [sap_num] = StringSave (iatep->title_list [seq_num]);
    sap_num++;
  }
  iatep = IDAndTitleEditFree (iatep);
}

static void ApplySourceAssistantToSeqEntryList (SourceAssistantPtr sap, SeqEntryPtr seq_list)
{
  IDAndTitleEditPtr iatep;
  Int4              seq_num, sap_num;
  
  if (sap == NULL || seq_list == NULL)
  {
    return;
  }
  
  iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  if (iatep == NULL || iatep->num_sequences < 1)
  {
    return;
  }
  
  sap_num = 0;
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (iatep->is_seg != NULL && iatep->is_seg [seq_num])
    {
      continue;
    }
    iatep->title_list [seq_num] = MemFree (iatep->title_list [seq_num]);
    iatep->title_list [seq_num] = StringSave (sap->defline_list [sap_num]);
    sap_num++;
  }
  ApplyIDAndTitleEditToSeqEntryList (seq_list, iatep);
  iatep = IDAndTitleEditFree (iatep);
}

static Boolean ShowRedSeqID (Int4 row, ValNodePtr row_list, Pointer userdata)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       row_num;
  
  if (OrganismMatchesAnotherRow (row, row_list, userdata))
  {
    return TRUE;
  }
  
  /* find the row we're interested in */
  for (row_vnp = row_list->next, row_num = 1;
       row_vnp != NULL && row_num != row; 
       row_vnp = row_vnp->next, row_num++)
  {
  }
  if (row_vnp == NULL || row_vnp->data.ptrvalue == NULL)
  {
    return TRUE;
  }
  
  /* the second column contains the organism names */
  col_vnp = row_vnp->data.ptrvalue;
  if (col_vnp == NULL 
      || col_vnp->next == NULL 
      || StringHasNoText (col_vnp->next->data.ptrvalue))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static Int4 GetStandardTableDisplayDialogWidth (SequencesFormPtr sqfp)
{
  Int4 doc_width = 0;
  RecT r;
  
  if (sqfp == NULL || sqfp->annot_tbs == NULL)
  {
    SelectFont (GetTableDisplayDefaultFont ());
    doc_width = CharWidth ('0') * 40;
  }
  else
  {
    GetPosition (sqfp->annot_tbs, &r);
    doc_width = r.right - r.left;
  }
  return doc_width;
}


NLM_EXTERN Boolean SourceAssistantForDeflines (SeqEntryPtr seq_list, Int4 doc_width, Int2 seqPackage)
{
  SourceAssistantData sad;
  WindoW              w;
  GrouP               h, k, g, g2, mod_btn_grp, c;
  Int4                i;
  ButtoN              export_btn;
  PrompT              ppt1, ppt2;
  ValNodePtr          modifier_choice_list = NULL, qual_list = NULL;
  Boolean             rval = FALSE;
  ButtoN              b;
    
  SeqEntryListToSourceAssistant (seq_list, &sad);
  if (sad.num_deflines < 1)
  {
    return rval;
  }
  
  SendHelpScrollMessage (helpForm, "Organism Page", "Add Source Modifiers");  
  
  sad.seqPackage = seqPackage;
    
  sad.done = FALSE;
  sad.cancelled = FALSE;
  modedit_widths [0] = 7;
  modedit_widths [1] = 18;
    
  w = MovableModalWindow (-20, -13, -10, -10, "Specify Source Modifiers", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  k = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (k, "Import source modifiers table", 0, popupMenuHeight, programFont, 'l');
  b = PushButton (k, "Select File", SourceAssistantImportModsBtn);
  SetObjectExtra (b, &sad, NULL);
  
  g = NormalGroup (h, -1, 0, "", programFont, NULL);
  g2 = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (g2, "Select Modifier", 0, popupMenuHeight, programFont, 'l');
  
  ValNodeAddPointer (&modifier_choice_list, eModifierType_Organism, StringSave ("Organism"));
  ValNodeAddPointer (&modifier_choice_list, eModifierType_Location, StringSave ("Location"));
  qual_list = GetSourceQualDescList (TRUE, TRUE, FALSE, FALSE);
  AddOneSourceQualDesc (&qual_list, "Fwd-PCR-primer-name", FALSE, SUBSRC_fwd_primer_name, 0);
  AddOneSourceQualDesc (&qual_list, "Fwd-PCR-primer-seq", FALSE, SUBSRC_fwd_primer_seq, 0);
  AddOneSourceQualDesc (&qual_list, "Rev-PCR-primer-name", FALSE, SUBSRC_rev_primer_name, 0);
  AddOneSourceQualDesc (&qual_list, "Rev-PCR-primer-seq", FALSE, SUBSRC_rev_primer_seq, 0);
  qual_list = ValNodeSort (qual_list, SortVnpBySourceQualDesc);
  ValNodeLink (&modifier_choice_list, qual_list);

  /* note - ValNodeSelectionDialog cleans up modifier_choice_list */
  sad.mod_type_dlg = ValNodeSelectionDialog (g2, modifier_choice_list, 6,
                                             SourceQualValNodeName,
                                             ValNodeSimpleDataFree,
                                             SourceQualValNodeDataCopy,
                                             SourceQualValNodeMatch,
                                             "modifier",
                                             NULL, NULL, FALSE);
  modifier_choice_list = NULL;
  qual_list = NULL;
  
  /* set default value for mod_type_dlg */
  ValNodeAddPointer (&qual_list, eModifierType_Organism, StringSave ("Organism"));
  PointerToDialog (sad.mod_type_dlg, qual_list);
  qual_list = ValNodeFreeData (qual_list);
    
  mod_btn_grp = HiddenGroup (g, 2, 0, NULL);
  b = PushButton (mod_btn_grp, "Apply One Value to All", ApplyOrgModColBtn);
  SetObjectExtra (b, &sad, NULL);
  b = PushButton (mod_btn_grp, "Edit Individual Values", EditOrgModColBtn);
  SetObjectExtra (b, &sad, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g2, (HANDLE) mod_btn_grp, NULL);
    
  sad.mod_doc = DocumentPanel (h, doc_width, stdLineHeight * 4);
  SetDocAutoAdjust (sad.mod_doc, TRUE);

  ppt1 = StaticPrompt (h, "Sequence IDs in red have missing organism names", 0, 0, programFont, 'l');
  ppt2 = StaticPrompt (h, "or have source information that matches at least one other sequence.", 0, 0, programFont, 'l');

  sad.orgmod_dlg = TableDisplayDialog (h, doc_width, stdLineHeight * 8, 1, 1,
                                       OrgModDblClick, &sad,
                                       ShowRedSeqID, NULL);
  UpdateOrgModDlg (&sad);  

  export_btn = PushButton (h, "Export Source Modifier Table", SourceAssistantExport);
  SetObjectExtra (export_btn, &sad, NULL);
  
  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton(c, "OK", SourceAssistantOk);
  SetObjectExtra (b, &sad, NULL);
  b = PushButton (c, "Clear All Source Modifiers", SourceAssistantClearAllModifiers);
  SetObjectExtra (b, &sad, NULL);
  b = PushButton(c, "Cancel", SourceAssistantCancel);
  SetObjectExtra (b, &sad, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) k,
                              (HANDLE) g, 
                              (HANDLE) sad.mod_doc,
                              (HANDLE) sad.orgmod_dlg, 
                              (HANDLE) ppt1,
                              (HANDLE) ppt2,
                              (HANDLE) export_btn, 
                              (HANDLE) c, 
                              NULL);
  
  Show(w); 
  Select (w);
  while (!sad.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (!sad.cancelled)
  {
    ApplySourceAssistantToSeqEntryList (&sad, seq_list);
    rval = TRUE;
  }
  
  for (i = 0; i < sad.num_deflines; i++)
  {
    sad.defline_list[i] = MemFree (sad.defline_list[i]);
    sad.id_list[i] = MemFree (sad.id_list[i]);
  }
  sad.defline_list = MemFree (sad.defline_list);
  sad.id_list = MemFree (sad.id_list);
  
  Remove (w);
  return rval;
}


static void SourceAssistant (ButtoN b)
{
  SequencesFormPtr    sqfp;
  SeqEntryPtr         seq_list;
  FastaPagePtr        fpp;
  Int4                doc_width = stdCharWidth * 40;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL) return;
  
  doc_width = GetStandardTableDisplayDialogWidth (sqfp);
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);

  if (SourceAssistantForDeflines (seq_list, doc_width, sqfp->seqPackage)) {
    SeqEntryPtrToSourceTab (sqfp);
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) 
    {
      Reset (fpp->doc);
      FormatFastaDoc (fpp);
    }
  }
  
}

static void ApplyOneValueToAllSequencesDialog (ButtoN b, CharPtr mod_name)
{
  WindoW               w;
  SequencesFormPtr     sqfp;
  CharPtr              value;
  SeqEntryPtr          seq_list;
  FastaPagePtr         fpp;
  Int4                 mod_type;
  Boolean              is_nontext;
  IDAndTitleEditPtr    iatep;
  Boolean              all_default;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL)
  {
    return;
  }
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    return;
  }
  
  is_nontext = IsNonTextModifier (mod_name);
  mod_type = GetModifierType (mod_name);
  value = GetFirstModValueFromSeqEntryTitles (seq_list, mod_name);
 
  if (mod_type == eModifierType_GeneticCode)
  {
    if (! ContinueWithAutopopulatedGeneticCodes (seq_list, NULL, NULL, -1))
    {
      value = MemFree (value);
      return;
    }
  }
  else if (!StringHasNoText (value)  && ! is_nontext
      && seq_list->next != NULL)
  {
    iatep = SeqEntryListToIDAndTitleEdit (seq_list);
    all_default = IDAndTitleEditHasAllDefaultValues (mod_type, mod_name, iatep);
    iatep = IDAndTitleEditFree (iatep);
  
    if (!all_default
        && ANS_YES != Message (MSG_YN, "Warning!  Some sequences already contain "
                            "a value for %s.  Are you sure you want to "
                            "overwrite these values?", mod_name))
    {
      value = MemFree (value);
      return;
    }
  }


  w = ParentWindow ((Nlm_GraphiC) b);
  if (w != (WindoW)sqfp->form)
  {
    Remove (w);
  }

  ApplyOrgModColumnOrCell (mod_name, value, -1, NULL, seq_list, NULL, 0, sqfp->seqPackage);
  value = MemFree (value);
  
  SeqEntryPtrToSourceTab (sqfp);
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp != NULL) 
  {
    Reset (fpp->doc);
    FormatFastaDoc (fpp);
  } 
}

static void ApplyModValuesIndividually (ButtoN b, CharPtr mod_name)
{
  SequencesFormPtr     sqfp;
  SeqEntryPtr          seq_list;
  FastaPagePtr         fpp;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL)
  {
    return;
  }
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    return;
  }
  
  EditOrgModColumn (mod_name, NULL, seq_list, sqfp->seqPackage);
  SeqEntryPtrToSourceTab (sqfp);  
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp != NULL) 
  {
    Reset (fpp->doc);
    FormatFastaDoc (fpp);
  }
}

static void ApplyModValuesAllBtn (ButtoN b)
{
  WindoW w;
  CharPtr mod_name;
  
  w = ParentWindow ((Nlm_GraphiC) b);
  mod_name = GetObjectExtra (w);
  if (StringHasNoText (mod_name))
  {
    return;
  }
  ApplyOneValueToAllSequencesDialog (b, mod_name);
  Remove (w);  
}

static void ApplyModValuesIndividuallyBtn (ButtoN b)
{
  WindoW w;
  CharPtr mod_name;
  
  w = ParentWindow ((Nlm_GraphiC) b);
  mod_name = GetObjectExtra (w);
  if (StringHasNoText (mod_name))
  {
    return;
  }
  ApplyModValuesIndividually (b, mod_name);
  Remove (w);  
}

static void SpecifyModValueButton (ButtoN b, CharPtr mod_name)
{
  SequencesFormPtr sqfp;
  WindoW           w;
  GrouP            h;
  SeqEntryPtr      seq_list;
  Char             title[255];
  ButtoN           apply_one_btn;
  ButtoN           apply_all_btn;
  ButtoN           cancel_btn;
  CharPtr          mod_label = NULL;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL || StringHasNoText (mod_name))
  {
    return;
  }
  
  if (StringICmp (mod_name, "moltype") == 0)
  {
    mod_label = StringSave ("molecule type");
  }
  else
  {
    mod_label = StringSave (mod_name);
  }
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    Message (MSG_ERROR, "You must add sequences before you can add %s information!", mod_label);
    mod_label = MemFree (mod_label);
    return;
  }

  sprintf (title, "Edit %s Information", mod_label);
  title [5] = TO_UPPER (title [5]);
  
  w = MovableModalWindow (-20, -13, -10, -10, title, NULL);
  SetObjectExtra (w, StringSave (mod_name), StdCleanupExtraProc);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  sprintf (title, "Apply one %s to all sequences", mod_label);
  apply_one_btn = PushButton (h, title, ApplyModValuesAllBtn);
  SetObjectExtra (apply_one_btn, sqfp, NULL);
  sprintf (title, "Apply %s to sequences individually", mod_label);
  apply_all_btn = PushButton (h, title, ApplyModValuesIndividuallyBtn);
  SetObjectExtra (apply_all_btn, sqfp, NULL);
  cancel_btn = PushButton (h, "Cancel", StdCancelButtonProc);
  mod_label = MemFree (mod_label);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) apply_one_btn,
                             (HANDLE) apply_all_btn,
                             (HANDLE) cancel_btn,
                             NULL);
  Show(w); 
  Select (w);  
  
}

static void SpecifyOrganismLocationGeneticCodeButton (ButtoN b)
{
  SequencesFormPtr sqfp;
  SeqEntryPtr      seq_list;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL)
  {
    return;
  }
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    Message (MSG_ERROR, "You must add sequences before you can add organisms, locations, or genetic codes!");
  }
  else
  {
    ApplyModValuesIndividually (b, "organism");
  }
}

static void ClearAllSequenceModifiers (ButtoN b)
{
  SequencesFormPtr sqfp;
  SeqEntryPtr      seq_list, sep, nsep;
  CharPtr          ttl;
  ValNodePtr       found_modifiers = NULL, mod_vnp;
  FastaPagePtr     fpp;
  Int4             mod_type;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL)
  {
    return;
  }
  
  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    return;
  }  
  
  if (ANS_YES != Message (MSG_YN, "Are you sure you want to clear all of the modifiers on all of your sequences?"))
  {
    return;
  }
  
  for (sep = seq_list; sep != NULL; sep = sep->next)
  {
    ttl = NULL;
    nsep = FindNucSeqEntry (sep);
    SeqEntryExplore (nsep, (Pointer) (&ttl), FindFirstTitle);
    found_modifiers = BuildModifierTypeList (found_modifiers, ttl, FALSE);
    for (mod_vnp = found_modifiers; mod_vnp != NULL; mod_vnp = mod_vnp->next)
    {
      mod_type = GetModifierType (mod_vnp->data.ptrvalue);
      if (mod_type != eModifierType_Protein
          && mod_type != eModifierType_Location
          && mod_type != eModifierType_Origin
          && mod_type != eModifierType_MolType
          && mod_type != eModifierType_Molecule
          && mod_type != eModifierType_Topology
          && mod_type != eModifierType_GeneticCode
          && mod_type != eModifierType_GeneticCodeComment
          && mod_type != eModifierType_Organism)
      {
        ApplyOneModToSeqEntry (nsep, mod_vnp->data.ptrvalue, NULL);
      }
    }
    found_modifiers = ValNodeFreeData (found_modifiers);
  }
  SeqEntryPtrToSourceTab (sqfp);
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp != NULL) 
  {
    Reset (fpp->doc);
    FormatFastaDoc (fpp);
  }
}

static void SummaryDblClick (PoinT cell_coord, CharPtr header_text, CharPtr cell_text, Pointer userdata)
{
  ValNodePtr        found_modifiers = NULL, vnp;
  SequencesFormPtr  sqfp;
  SeqEntryPtr       seq_list;
  Int4              mod_num;
  FastaPagePtr      fpp;
  IDAndTitleEditPtr iatep;
  
  if (cell_coord.y < 1 || userdata == NULL)
  {
    return;
  }
  
  sqfp = (SequencesFormPtr) userdata;

  seq_list = GetSeqEntryFromSequencesForm (sqfp);
  if (seq_list == NULL)
  {
    return;
  }  
  
  iatep = SeqEntryListToIDAndTitleEdit (seq_list);
  
  /* get list of modifiers */ 
  found_modifiers = GetListOfCurrentSourceModifiers (iatep);
   
  for (vnp = found_modifiers, mod_num = 1;
       vnp != NULL && mod_num < cell_coord.y;
       vnp = vnp->next, mod_num++)
  {
  }
  if (vnp != NULL)
  {
    EditOrgModColumn (vnp->data.ptrvalue, NULL, seq_list, sqfp->seqPackage);
    SeqEntryPtrToSourceTab (sqfp);  
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) 
    {
      Reset (fpp->doc);
      FormatFastaDoc (fpp);
    }
  }
  ValNodeFreeData (found_modifiers);
  iatep = IDAndTitleEditFree (iatep);
}

static GrouP CreateSourceTab (GrouP h, SequencesFormPtr sqfp)
{
  GrouP              mod_grp;
  GrouP              src_btns_grp;
  GrouP              k;
  Int4               doc_width;

  if (h == NULL || sqfp == NULL)
  {
    return NULL;
  }
 
  modedit_widths [0] = 7;
  modedit_widths [1] = 18;

  mod_grp = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (mod_grp, 10, 20);

  doc_width = GetStandardTableDisplayDialogWidth (sqfp);

  sqfp->org_doc = DocumentPanel (mod_grp, doc_width, stdLineHeight * 3);
  SetDocAutoAdjust (sqfp->org_doc, TRUE);
  sqfp->ident_org_grp = HiddenGroup (mod_grp, -1, 0, NULL);
  StaticPrompt (sqfp->ident_org_grp, "Some sequences have identical source information.", 0, 0, programFont, 'c');
  StaticPrompt (sqfp->ident_org_grp, "Edit the source information using the", 0, 0, programFont, 'c');
  StaticPrompt (sqfp->ident_org_grp, "'Add Source Modifiers' button below.", 0, 0, programFont, 'c');
  Disable (sqfp->source_assist_btn);
  
  sqfp->summary_dlg = TableDisplayDialog (mod_grp, doc_width, stdLineHeight * 8, 1, 1,
                                       SummaryDblClick, sqfp,
                                       NULL, NULL);
 

  sqfp->specify_orgs_btn = PushButton (mod_grp, 
                                       "Add Organisms, Locations, and Genetic Codes", 
                                       SpecifyOrganismLocationGeneticCodeButton);
  SetObjectExtra (sqfp->specify_orgs_btn, sqfp, NULL);
  Disable (sqfp->specify_orgs_btn);
  k = HiddenGroup (mod_grp, -1, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  src_btns_grp = HiddenGroup (k, 2, 0, NULL);
  SetGroupSpacing (src_btns_grp, 10, 10);
  sqfp->import_mod_btn = PushButton (src_btns_grp, "Import Source Modifiers", ImportModifiersButtonProc);
  SetObjectExtra (sqfp->import_mod_btn, sqfp, NULL);
  Disable (sqfp->import_mod_btn);
  sqfp->source_assist_btn = PushButton (src_btns_grp, "Add Source Modifiers", SourceAssistant);
  SetObjectExtra (sqfp->source_assist_btn, sqfp, NULL);
  Disable (sqfp->source_assist_btn);
  sqfp->specify_locs_btn = NULL;
  sqfp->specify_gcode_btn = NULL;
  sqfp->specify_mgcode_btn = NULL;
    
  sqfp->clear_mods_btn = PushButton (k, "Clear All Source Modifiers", ClearAllSequenceModifiers);
  SetObjectExtra (sqfp->clear_mods_btn, sqfp, NULL);
  Disable (sqfp->clear_mods_btn);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) src_btns_grp, (HANDLE) sqfp->clear_mods_btn, NULL);

  SeqEntryPtrToSourceTab (sqfp);
  AlignObjects (ALIGN_CENTER, (HANDLE) sqfp->org_doc,
                              (HANDLE) sqfp->ident_org_grp,
                              (HANDLE) sqfp->summary_dlg,
                              (HANDLE) sqfp->specify_orgs_btn,
                              (HANDLE) k,
                              NULL);
  
  return mod_grp;
}

typedef struct fastasummary 
{
  DIALOG_MESSAGE_BLOCK
  PrompT summary_ppt;
} FastaSummaryData, PNTR FastaSummaryPtr;

static void SequencesToFastaSummary (DialoG d, Pointer userdata)
{
  FastaSummaryPtr dlg;
  SeqEntryPtr     seq_list, sep, nsep;
  Int4            num_sequences = 0, tot_len = 0;
  BioseqPtr       bsp;
  CharPtr         str_format = "%d sequences, total length %d";
  Char            str[255];
  
  dlg = (FastaSummaryPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    return;
  }
  
  seq_list = (SeqEntryPtr) userdata;
  for (sep = seq_list; sep != NULL; sep = sep->next)
  {
    num_sequences++;
    if (IS_Bioseq (sep)) 
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) 
      {
        tot_len += bsp->length;
      }
    } 
    else if (IS_Bioseq_set (sep)) 
    {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL && IS_Bioseq (nsep)) 
      {
        bsp = (BioseqPtr) nsep->data.ptrvalue;
        if (bsp != NULL) 
        {
          tot_len += bsp->length;
        }
      }
    }
  }
  sprintf (str, str_format, num_sequences, tot_len);
  SetTitle (dlg->summary_ppt, str);
}

static DialoG FastaSummaryDialog (GrouP parent)
{
  FastaSummaryPtr dlg;
  GrouP           p;
  
  dlg = (FastaSummaryPtr) MemNew (sizeof (FastaSummaryData));
  p = HiddenGroup (parent, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SequencesToFastaSummary;
  dlg->fromdialog = NULL;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;

  dlg->summary_ppt = StaticPrompt (p, NULL, stdCharWidth * 20,
                                   popupMenuHeight, programFont, 'l');  
  return (DialoG) p;
}

/* The Sequence Assistant will allow users to paste in FASTA or create sequences one at a time.
 * The list control for selecting a sequence to view and edit should be a DocPanel, so that we
 * can add and remove sequences without needing to destroy and recreate the dialog.
 * The sequences created should be stored as a SeqEntry list.
 */
typedef struct sequenceassistant 
{
  DoC            sequence_selector;
  DialoG         summary_dlg;
  DialoG         sequence_table;
  ButtoN         edit_btn;
  ButtoN         delete_btn;
  ButtoN         delete_all_btn;
  ButtoN         import_btn;
  
  Int2           sequence_row;  
  
  SeqEntryPtr    seq_list;
  
  Int2           seqPackage;

  Boolean        done;
  Boolean        cancelled;  
} SequenceAssistantData, PNTR SequenceAssistantPtr;

static ValNodePtr PrepareSequenceAssistantTableData (SequenceAssistantPtr sap)
{
  ValNodePtr         column_list = NULL, row_list = NULL;
  ValNodePtr         header_list, header_vnp;
  Int4               column_width, num_columns = 0;
  Int4               max_column_width = 20;
  SeqEntryPtr        sep, nsep;
  BioseqPtr          bsp;
  Char               tmp[128];
  CharPtr            ttl = NULL, id;
  CharPtr            valstr;
  
  if (sap == NULL)
  {
    return NULL;
  }
  AddDefaultModifierValues (sap->seq_list);
      
  /* create header line for table */
  /* store max column width in choice */
  ValNodeAddPointer (&column_list, 6, StringSave ("Seq ID"));
  ValNodeAddPointer (&column_list, 6, StringSave ("Length"));
  ValNodeAddPointer (&column_list, 8, StringSave ("Molecule"));
  ValNodeAddPointer (&column_list, 8, StringSave ("Topology"));
  ValNodeAddPointer (&column_list, 11, StringSave ("Title"));
    
  ValNodeAddPointer (&row_list, 0, column_list);
  header_list = column_list;
  
  num_columns = 3;

  /* create data lines for table */
  for (sep = sap->seq_list; sep != NULL; sep = sep->next)
  {
    bsp = NULL;
    
    if (IS_Bioseq (sep))
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
    }
    else if (IS_Bioseq_set (sep))
    {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL && IS_Bioseq (nsep))
      {
        bsp = (BioseqPtr) nsep->data.ptrvalue;
      }
    }
    
    column_list = NULL;
    
    /* add Sequence ID */
    header_vnp = header_list;
    id = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);
    column_width = MAX (StringLen (id), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, id);
      
    /* add length */
    header_vnp = header_vnp->next;
    sprintf (tmp, "%d", bsp->length);
    column_width = MAX (StringLen (tmp), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, StringSave (tmp));
    
    ttl = NULL;
    SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);

    /* add molecule */
    header_vnp = header_vnp->next;
    valstr = FindValueFromPairInDefline ("moltype", ttl);
    if (StringHasNoText (valstr))
    {
      valstr = MemFree (valstr);
      valstr = StringSave ("Genomic DNA");
    }
    column_width = MAX (StringLen (tmp), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, valstr);

    /* add topology */
    header_vnp = header_vnp->next;
    valstr = FindValueFromPairInDefline ("topology", ttl);
    if (StringHasNoText (valstr))
    {
      valstr = MemFree (valstr);
      valstr = StringSave ("Linear");
    }
    column_width = MAX (StringLen (tmp), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, valstr);

    /* add title */  
    header_vnp = header_vnp->next;
    column_width = MAX (StringLen (ttl), header_vnp->choice);
    column_width = MIN (column_width, max_column_width);
    header_vnp->choice = column_width;
    ValNodeAddPointer (&column_list, 0, StringSave (ttl));
  
    ValNodeAddPointer (&row_list, 0, column_list);
  }
  return row_list;
}

static void PopulateSequenceSelector (SequenceAssistantPtr sap, DoC selector)
{
  SeqEntryPtr sep, nsep;
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  Char        tmp[128];
  Int2        seq_num;
  
  if (sap == NULL || selector == NULL)
  {
    return;
  }
  Reset (selector);
  for (sep = sap->seq_list, seq_num = 0; sep != NULL; sep = sep->next, seq_num++)
  {
    bsp = NULL;
    if (IS_Bioseq (sep))
    {
      bsp = sep->data.ptrvalue;
    }
    else if (IS_Bioseq_set (sep))
    {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL && IS_Bioseq (nsep))
      {
        bsp = nsep->data.ptrvalue;
      }
    }
    
    if (bsp != NULL)
    {
      /* add to sequence_selector doc */
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp) - 1);
  	  AppendText (selector, tmp, &faParFmt, &faColFmt, programFont);  	  
    }
  }
  InvalDocRows (selector, 0, 0, 0);
  sap->sequence_row = -1;
  Disable (sap->edit_btn);
  Disable (sap->delete_btn);
  
}

static void UpdateSequenceAssistant (SequenceAssistantPtr sap)
{
  ValNodePtr row_list;
  
  if (sap == NULL)
  {
    return;
  }
  
  row_list = PrepareSequenceAssistantTableData (sap);
  PointerToDialog (sap->sequence_table, row_list);
  FreeTableDisplayRowList (row_list);

  PopulateSequenceSelector (sap, sap->sequence_selector);
  PointerToDialog (sap->summary_dlg, sap->seq_list);
  
  /* set title for import button */
  if (sap->seq_list != NULL && sap->seq_list->next != NULL)
  {
    SetTitle (sap->import_btn, "Import Additional Nucleotide FASTA");
  }
  else
  {
    SetTitle (sap->import_btn, "Import Nucleotide FASTA");
  }
  
  if (sap->seq_list == NULL)
  {
    Disable (sap->delete_all_btn);
  }
  else
  {
    Enable (sap->delete_all_btn);
  }
}


static void ImportSequenceAssistantEditData (SequenceAssistantPtr sap, CharPtr seq_str)
{
  Char         path [PATH_MAX];
  SeqEntryPtr  new_sep_list;

  if (sap == NULL || StringHasNoText (seq_str))
  {
    return;
  }
  
  TmpNam (path);
  new_sep_list = GetSequencesFromFile (path, sap->seq_list);
  if (new_sep_list != NULL)
  {
    ValNodeLink (&sap->seq_list, new_sep_list);
  }
  UpdateSequenceAssistant (sap);

  FileRemove (path);
}

static SeqEntryPtr CopySeqEntryList (SeqEntryPtr seq_list)
{
  SeqEntryPtr new_seq_list, new_seq, last_seq;
  ErrSev      oldsev;
  
  if (seq_list == NULL)
  {
    return NULL;
  }
  
  oldsev = ErrSetMessageLevel (SEV_MAX);
  
  new_seq = AsnIoMemCopy ((Pointer) seq_list,
                          (AsnReadFunc) SeqEntryAsnRead,
                          (AsnWriteFunc) SeqEntryAsnWrite);
  new_seq_list = new_seq;
  last_seq = new_seq;
  
  seq_list = seq_list->next;
  
  while (last_seq != NULL && seq_list != NULL)
  {   
    new_seq = AsnIoMemCopy ((Pointer) seq_list,
                            (AsnReadFunc) SeqEntryAsnRead,
                            (AsnWriteFunc) SeqEntryAsnWrite);
    last_seq->next = new_seq;
    last_seq = last_seq->next;
    seq_list = seq_list->next;
  }
  
  ErrSetMessageLevel (oldsev);
                          
  return new_seq_list;
}


/* This section of code is used for correcting errors in the sequence IDs and
 * titles.
 */

/* This structure contains the data used for editing and updating the new
 * and existing sequence lists.
 */
typedef struct seqidedit 
{
  WindoW            w;
  IDAndTitleEditPtr iatep_new;
  IDAndTitleEditPtr iatep_current;
  DialoG            new_dlg;
  DialoG            current_dlg;
  GrouP             show_all_grp;
  DoC               auto_correct_doc;
  DialoG            bracket_dlg;
  PaneL             badvalue_pnl;
  PaneL             unrec_mod_pnl;
  ButtoN            auto_correct_btn;
  Boolean           auto_correct_ids;
  Boolean           auto_correct_bracketing;
  Boolean           auto_correct_modnames;
  ButtoN            accept_btn;
  ButtoN            refresh_err_btn;
  ButtoN            refresh_error_list_btn;
  
  Boolean           seqid_edit_phase;
  Boolean           is_nuc;
} SeqIdEditData, PNTR SeqIdEditPtr;
 
/* This section of code is used for detecting errors in lists of sequence IDs
 * and titles.
 */
static Uint2 idedit_types [] = {
  TAGLIST_PROMPT, TAGLIST_PROMPT, TAGLIST_TEXT, TAGLIST_TEXT
};

static Uint2 idedit_widths [] = {
  6, 5, 10, 40,
};

static Boolean HasMissingIDs (IDAndTitleEditPtr iatep)
{
  Int4 i;
  
  if (iatep == NULL)
  {
    return FALSE;
  }
  for (i = 0; i < iatep->num_sequences; i++)
  {
    if (StringHasNoText (iatep->id_list [i]))
    {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean IsDuplicateEditID (IDAndTitleEditPtr iatep_new, Int4 new_pos, IDAndTitleEditPtr iatep_current)
{
  Int4 j;
  
  if (iatep_new == NULL || iatep_new->num_sequences == 0
      || new_pos < 0 || new_pos >= iatep_new->num_sequences)
  {
    return FALSE;
  }
  
  for (j = 0; j < iatep_new->num_sequences; j++)
  {
    if (j == new_pos)
    {
      continue;
    }
    if (StringICmp (iatep_new->id_list [new_pos], iatep_new->id_list [j]) == 0)
    {
      return TRUE;
    }
  }
  if (iatep_current != NULL)
  {
    for (j = 0; j < iatep_current->num_sequences; j++)
    {
      if (StringICmp (iatep_new->id_list [new_pos], iatep_current->id_list [j]) == 0)
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean EditHasDuplicateIDs (IDAndTitleEditPtr iatep_new, IDAndTitleEditPtr iatep_current)
{
  Int4 i;
  ValNodeBlock block;
  ValNodePtr vnp;
  CharPtr prev = NULL;
  Boolean found_duplicate = FALSE;
  
  if (iatep_new == NULL || iatep_new->num_sequences == 0)
  {
    return FALSE;
  }

  InitValNodeBlock (&block, NULL);
  for (i = 0; i < iatep_new->num_sequences; i++)
  {
    ValNodeAddPointerToEnd (&block, 0, iatep_new->id_list[i]);
  }
  if (iatep_current != NULL)
  {
    for (i = 0; i < iatep_current->num_sequences; i++)
    {
      ValNodeAddPointerToEnd (&block, 0, iatep_current->id_list[i]);
    }
  }
  if (block.head == NULL || block.head->next == NULL) 
  {
    block.head = ValNodeFree (block.head);
    return FALSE;
  }

  block.head = ValNodeSort (block.head, SortVnpByString);
  prev = block.head->data.ptrvalue;
  for (vnp = block.head->next; vnp != NULL && !found_duplicate; vnp = vnp->next) 
  {
    if (StringICmp (prev, vnp->data.ptrvalue) == 0) 
    {
      found_duplicate = TRUE;
    }
    prev = vnp->data.ptrvalue;
  }
  block.head = ValNodeFree (block.head);
  return found_duplicate;
}


static Boolean StringHasSpace (CharPtr str)
{
  CharPtr cp;
  Boolean has_space = FALSE;

  /* look for space in ID */
  for (cp = str; cp != NULL && *cp != 0 && !has_space; cp++) 
  {
    if (isspace (*cp)) {
      has_space = TRUE;
    }
  }
  return has_space;
}


static Boolean EditHasSpaceInIDs (IDAndTitleEditPtr iatep)
{
  Int4 i;
  Boolean has_space = FALSE;
  
  if (iatep != NULL)
  {
    for (i = 0; i < iatep->num_sequences && !has_space; i++)
    {
      has_space = StringHasSpace (iatep->id_list[i]);
    }
  }
  return has_space;
}


/* This function creates a list of suggested IDs and titles 
 * based on errors in the IDs of the original lists.
 * The errors that can be auto-corrected are:
 *     * spaces in the sequence IDs
 *     * brackets in the sequence IDs 
 */
static IDAndTitleEditPtr SuggestCorrectionForLocalIDs (IDAndTitleEditPtr iatep_new, IDAndTitleEditPtr iatep_current)
{
  CharPtr           add_str, cp;
  Int4              len, add_str_len;
  Int4              seq_num;
  IDAndTitleEditPtr iatep_corrected = NULL;
  CharPtr           new_title_start;
  
  if (iatep_new == NULL || iatep_new->num_sequences < 1)
  {
    return NULL;
  }
  
  iatep_corrected = IDAndTitleEditCopy (iatep_new);

  for (seq_num = 0; seq_num < iatep_corrected->num_sequences; seq_num++)
  {
    if (!StringHasNoText (iatep_new->id_list [seq_num])
        && IsDuplicateEditID (iatep_new, seq_num, iatep_current)
        && !StringHasNoText (iatep_new->title_list [seq_num])
        && iatep_new->title_list [seq_num][0] != '[')
    {
      len = StringCSpn (iatep_corrected->title_list [seq_num], " \t[]");
      if (len > 0)
      {
        add_str_len = len + StringLen (iatep_corrected->id_list [seq_num]);
        add_str = (CharPtr) MemNew ((add_str_len + 1) * sizeof (Char));
        if (add_str != NULL)
        {
          StringCpy (add_str, iatep_corrected->id_list [seq_num]);
          StringNCat (add_str, iatep_corrected->title_list [seq_num], len);
          add_str [add_str_len] = 0;
          iatep_corrected->id_list [seq_num] = MemFree (iatep_corrected->id_list [seq_num]);
          iatep_corrected->id_list [seq_num] = add_str;
          
          new_title_start = iatep_corrected->title_list [seq_num] + len;
          new_title_start += StringSpn (new_title_start, " \t");
          new_title_start = StringSave (new_title_start);
          iatep_corrected->title_list [seq_num] = MemFree (iatep_corrected->title_list [seq_num]);
          iatep_corrected->title_list [seq_num] = new_title_start;
        }
      }
    }
    
    /* suggest correction for brackets in sequence IDs */
    if ((cp = StringChr (iatep_corrected->id_list [seq_num], '[')) != NULL)
    {
      len = StringLen (cp);
      if (StringNCmp (cp, iatep_corrected->title_list [seq_num], len) != 0)
      {
        add_str_len = len + StringLen (iatep_corrected->title_list [seq_num]);
        add_str = (CharPtr) MemNew ((add_str_len + 2) * sizeof (Char));
        if (add_str != NULL)
        {
          StringCpy (add_str, cp);
          StringCat (add_str, " ");
          StringCat (add_str, iatep_corrected->title_list [seq_num]);
          iatep_corrected->title_list [seq_num] = MemFree (iatep_corrected->title_list [seq_num]);
          iatep_corrected->title_list [seq_num] = add_str;
          *cp = 0;
        }
      }
      *cp = 0;
    }

    /* suggest correction for spaces in sequence IDs */
    cp = iatep_corrected->id_list [seq_num];
    while (*cp != 0) {
      if (isspace (*cp)) {
        *cp = '_';
      }
      cp++;
    }

  }
  return iatep_corrected;  
}

/* This function creates a list of suggested titles based on bracketing errors
 * in the original list.
 */
static IDAndTitleEditPtr SuggestCorrectionForTitleBracketing (IDAndTitleEditPtr iatep_orig)
{
  IDAndTitleEditPtr iatep_corrected;
  Int4              seq_num, msg_num;
  
  if (iatep_orig == NULL || iatep_orig->num_sequences < 1)
  {
    return NULL;
  }
  
  iatep_corrected = IDAndTitleEditCopy (iatep_orig);

  for (seq_num = 0; seq_num < iatep_corrected->num_sequences; seq_num++)
  {
    msg_num = DetectBadBracketing (iatep_corrected->title_list[seq_num]);
    if (msg_num != 0)
    {
      iatep_corrected->title_list[seq_num] = MemFree (iatep_corrected->title_list[seq_num]);
      iatep_corrected->title_list[seq_num] = SuggestCorrectBracketing(iatep_orig->title_list[seq_num]);
    }
  }
  return iatep_corrected;
}

/* This function indicates whether there are bracketing errors in the supplied list
 * of sequence IDs and titles. 
 */
static Boolean EditNeedsBracketingFixes (IDAndTitleEditPtr iatep_orig)
{
  Int4    seq_num;
  Boolean needs_fix = FALSE;

  if (iatep_orig == NULL || iatep_orig->num_sequences < 1)
  {
    return FALSE;
  }
  
  for (seq_num = 0; seq_num < iatep_orig->num_sequences && ! needs_fix; seq_num++)
  {
    if (DetectBadBracketing (iatep_orig->title_list[seq_num]) != 0)
    {
      needs_fix = TRUE;
    }
  }
  return needs_fix;
}

/* These functions are used to find and list unrecognized modifier names in
 * definition lines.
 */

static Boolean IsUnrecognizedModifierName (ModifierInfoPtr mip, Boolean is_nuc)
{
  if (mip == NULL
      || (mip->modtype == eModifierType_SourceQual
  	      && mip->subtype == 255
  	      && StringICmp (mip->name, "note-subsrc") != 0 
  	      && StringICmp (mip->name, "note-orgmod") != 0)
  	  || (!is_nuc && mip->modtype != eModifierType_Protein
  	      && StringICmp (mip->name, "note-orgmod") != 0))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}
 
/* This function searches a single definition line for unrecognized modifiers
 * and adds them to the list if they are not already on the list.
 */
static void 
AddUnrecognizedModifiersForOneDefinitionLine 
(CharPtr         defline,
 ValNodePtr PNTR unrecognized_list,
 Boolean         is_nuc)
{
  ValNodePtr      modifier_info_list;
  ValNodePtr      info_vnp, type_vnp;
  ModifierInfoPtr mip;
  
  if (StringHasNoText (defline) || unrecognized_list == NULL)
  {
    return;
  }
  
  modifier_info_list = ParseAllBracketedModifiers (defline);
  for (info_vnp = modifier_info_list; info_vnp != NULL; info_vnp = info_vnp->next)
  {
    mip = (ModifierInfoPtr)info_vnp->data.ptrvalue;
    if (mip == NULL || !IsUnrecognizedModifierName (mip, is_nuc))
    {
      mip = ModifierInfoFree (mip);
      continue;
    }
    for (type_vnp = *unrecognized_list;
         type_vnp != NULL && StringICmp (mip->name, type_vnp->data.ptrvalue) != 0;
         type_vnp = type_vnp->next)
    {
    }
    if (type_vnp == NULL)
    {
      ValNodeAddPointer (unrecognized_list, 0, StringSave (mip->name));
    }
  }
}

/* This function searches all of the titles in the supplied list of sequence IDs and titles
 * for unrecognized modifier names and adds them to the unrecognized_list if they are
 * not already in the list.
 */
static void 
AddUnrecognizedModifiers 
(IDAndTitleEditPtr iatep,
 ValNodePtr PNTR   unrecognized_list,
 Boolean           is_nuc)
{
  Int4            seq_num;
  
  if (iatep == NULL || unrecognized_list == NULL)
  {
    return;
  }

  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    AddUnrecognizedModifiersForOneDefinitionLine (iatep->title_list [seq_num],
                                                  unrecognized_list,
                                                  is_nuc);
  }  
}

/* This function searches all of the titles in the new and existing sets of sequences
 * for unrecognized modifier names and generates a list of unique unrecognized modifier names.
 */
static ValNodePtr 
ListUnrecognizedModifiers 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 Boolean           is_nuc)
{
  ValNodePtr      unrecognized_list = NULL;

  AddUnrecognizedModifiers (iatep_new, &unrecognized_list, is_nuc);
  AddUnrecognizedModifiers (iatep_current, &unrecognized_list, is_nuc);
  
  return unrecognized_list;
}

/* This section of code will look for inappropriate values in definition line pairs */
typedef struct badvalue
{
  CharPtr seq_id;
  CharPtr mod_name;
  CharPtr value;
} BadValueData, PNTR BadValuePtr;

static BadValuePtr BadValueFree (BadValuePtr bvp)
{
  if (bvp != NULL)
  {
    bvp->seq_id = MemFree (bvp->seq_id);
    bvp->mod_name = MemFree (bvp->mod_name);
    bvp->value = MemFree (bvp->value);
    bvp = MemFree (bvp);
  }
  return bvp;
}

static ValNodePtr BadValueListFree (ValNodePtr list)
{
  if (list != NULL)
  {
    list->next = BadValueListFree (list->next);
    list->data.ptrvalue = BadValueFree (list->data.ptrvalue);
    list = ValNodeFree (list);
  }
  return list;
}

static BadValuePtr BadValueNew (CharPtr seq_id, CharPtr mod_name, CharPtr value)
{
  BadValuePtr bvp;
  
  bvp = (BadValuePtr) MemNew (sizeof (BadValueData));
  if (bvp != NULL)
  {
    bvp->seq_id = StringSave (seq_id);
    bvp->mod_name = StringSave (mod_name);
    if (StringHasNoText (value))
    {
      bvp->value = NULL;
    }
    else
    {
      bvp->value = StringSave (value);
    }
  }
  return bvp;
}

/* The FixModName structure and the functions SetFixModNameAccept
 * and FixOneModifierName are used to present a dialog that allows
 * a user to replace a modifier name, either for one sequence or
 * for all sequences.
 * SetFixModNameAccept is used to prevent the user from clicking on
 * Accept before choosing a new modifier name.
 */
typedef struct fixmodname
{
  DialoG  name_list;
  ButtoN  accept_btn;
} FixModNameData, PNTR FixModNamePtr;

static void SetFixModNameAccept (Pointer userdata)
{
  FixModNamePtr fmp;
  ValNodePtr    vnp;
  Boolean      ok_to_accept = TRUE;

  fmp = (FixModNamePtr) userdata;
  if (fmp == NULL)
  {
    return;
  }
  
  vnp = DialogToPointer (fmp->name_list);
  if (vnp == NULL)
  {
    ok_to_accept = FALSE;
  }
  vnp = ValNodeFreeData (vnp);
  
  if (ok_to_accept)
  {
    Enable (fmp->accept_btn);
  }
  else
  {
    Disable (fmp->accept_btn);
  }
  
}

static CharPtr 
ReplaceOneModifierValue 
(CharPtr title,
 CharPtr orig_name, 
 CharPtr orig_value,
 CharPtr repl_value,
 Boolean is_nontext,
 Boolean copy_to_note);
static void 
UpdateIdAndTitleEditDialog 
(DialoG            d,
 IDAndTitleEditPtr iatep_new, 
 IDAndTitleEditPtr iatep_current, 
 Boolean           seqid_edit_phase,
 Boolean           show_all,
 Boolean           is_nuc);
static void ShowErrorInstructions (Pointer userdata);
static void ScrollTagListToSeqId (DialoG d, CharPtr seq_id);


static Boolean 
FixOneModifierName 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 CharPtr           seq_id,
 CharPtr           orig_mod_name,
 Boolean           is_nuc)
{
  ValNodePtr vnp;
  WindoW     w;
  ValNodePtr mod_choices;
  GrouP      h, action_type, c;
  PrompT       p;
  FixModNamePtr fmp;
  ModalAcceptCancelData acd;
  Boolean               rval = FALSE;
  ButtoN                b;
  CharPtr               prompt_txt;
  CharPtr               prompt_fmt = "Please choose a valid modifier name to replace %s:";
  CharPtr               radio_txt;
  CharPtr               radio_fmt = "For sequence '%s' only";
  CharPtr               repl_name;
  Int4                  action_type_val, seq_num;
  
  
  if ((iatep_new == NULL && iatep_current == NULL)
      || StringHasNoText (orig_mod_name))
  {
    return FALSE;
  }
  
  fmp = (FixModNamePtr) MemNew (sizeof (FixModNameData));
  if (fmp == NULL)
  {
    return FALSE;
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, "Replace Modifier Name", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  SetObjectExtra (w, fmp, StdCleanupExtraProc);
  
  prompt_txt = (CharPtr) MemNew ((StringLen (prompt_fmt) + StringLen (orig_mod_name)) * sizeof (Char));
  if (prompt_txt != NULL)
  {
    sprintf (prompt_txt, prompt_fmt, orig_mod_name);
  }
  p = StaticPrompt (h, prompt_txt,
                    0, 0, programFont, 'l');
  prompt_txt = MemFree (prompt_txt);

  mod_choices = GetFastaModifierList (is_nuc, !is_nuc);
  fmp->name_list = ValNodeSelectionDialog (h, mod_choices, 6,
                                          SourceQualValNodeName,
                                          ValNodeSimpleDataFree,
                                          SourceQualValNodeDataCopy,
                                          SourceQualValNodeMatch,
                                          "modifier",
                                          SetFixModNameAccept, fmp, FALSE);

  if (StringHasNoText (seq_id) 
      || (iatep_new == NULL && iatep_current->num_sequences == 1)
      || (iatep_current == NULL && iatep_new->num_sequences == 1))
  {
    action_type = NULL;
  }
  else
  {
    action_type = HiddenGroup (h, 0, 2, NULL);
    SetGroupSpacing (action_type, 10, 10);
    RadioButton (action_type, "For all sequences");
    radio_txt = (CharPtr) MemNew ((StringLen (radio_fmt) + StringLen (seq_id)) * sizeof (Char));
    if (radio_txt != NULL)
    {
      sprintf (radio_txt, radio_fmt, seq_id);
    }
    RadioButton (action_type, radio_txt);
    radio_txt = MemFree (radio_txt);
    SetValue (action_type, 1);
  }
  
  c = HiddenGroup (h, 2, 0, NULL);
  fmp->accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (fmp->accept_btn, &acd, NULL);
  Disable (fmp->accept_btn);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p,
                              (HANDLE) fmp->name_list,
                              (HANDLE) c, 
                              (HANDLE) action_type, 
                              NULL);

  Show (w);
  Select (w);
  
  acd.cancelled = FALSE;
  acd.accepted = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.cancelled)
  {
    rval = FALSE;
  }
  else
  {
    vnp = DialogToPointer (fmp->name_list);
    repl_name = SourceQualValNodeName (vnp);
    if (action_type == NULL)
    {
      action_type_val = 1;
    }
    else
    {
      action_type_val = GetValue (action_type);
    }
    
    for (seq_num = 0; iatep_new != NULL && seq_num < iatep_new->num_sequences; seq_num++)
    {
      if (action_type_val == 1 /* replace value for all sequences */
          || StringCmp (iatep_new->id_list [seq_num], seq_id) == 0)
      {
        iatep_new->title_list [seq_num] = ReplaceOneModifierName (iatep_new->title_list [seq_num],
                                                                  orig_mod_name, 
                                                                  repl_name);
      }
    }
    for (seq_num = 0; iatep_current != NULL && seq_num < iatep_current->num_sequences; seq_num++)
    {
      if (action_type_val == 1 /* replace value for all sequences */
          || StringCmp (iatep_current->id_list [seq_num], seq_id) == 0)
      {
        iatep_current->title_list [seq_num] = ReplaceOneModifierName (iatep_current->title_list [seq_num],
                                                                      orig_mod_name, 
                                                                      repl_name);
      }
    }
    repl_name = MemFree (repl_name);
    vnp = ValNodeFreeData (vnp);
    
    rval = TRUE;
  }
    
  Remove (w);
  
  return rval;    
}

/* This function presents a dialog that allows a user to replace a value for a
 * modifier, either in a single sequence or for every sequence.
 * The user also has the option to copy the original value into a note.
 */
static Boolean 
FixOneModifierValue 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 CharPtr           seq_id,
 CharPtr           orig_mod_name,
 CharPtr           orig_mod_value,
 Int4              mod_type)
{
  WindoW     w;
  GrouP      h, action_type, c, instr_grp;
  PrompT       p;
  ModalAcceptCancelData acd;
  Boolean               rval = FALSE;
  ButtoN                b, accept_btn;
  CharPtr               prompt_txt;
  CharPtr               prompt_fmt = "Please choose a valid value for %s to replace %s where %s=%s:";
  CharPtr               radio_txt;
  CharPtr               radio_fmt = "For sequence '%s' only";
  Int4                  action_type_val, seq_num;
  ButtoN                copy_to_note_btn;
  DialoG                new_value_dlg;
  Boolean               is_nontext = FALSE;
  CharPtr               new_value;
  Boolean               copy_to_note;
  
  if ((iatep_new == NULL && iatep_current == NULL)
      || StringHasNoText (orig_mod_name)
      || StringHasNoText (seq_id))
  {
    return FALSE;
  }
    
  w = MovableModalWindow (-20, -13, -10, -10, "Replace Modifier Value", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  prompt_txt = (CharPtr) MemNew ((StringLen (prompt_fmt) 
                                  + 2 * StringLen (orig_mod_name) 
                                  + 2 * StringLen (orig_mod_value)) * sizeof (Char));
  if (prompt_txt != NULL)
  {
    sprintf (prompt_txt, prompt_fmt, orig_mod_name, 
                                     orig_mod_value == NULL ? "" : orig_mod_value,
                                     orig_mod_name, 
                                     orig_mod_value == NULL ? "" : orig_mod_value);
  }
  p = StaticPrompt (h, prompt_txt,
                    0, 0, programFont, 'l');
  prompt_txt = MemFree (prompt_txt);
  
  if (mod_type == eModifierType_SourceQual)
  {
    is_nontext = IsNonTextModifier (orig_mod_name);
  }
  
  instr_grp = MakeInstructionGroup (h, is_nontext, mod_type);  
  
  new_value_dlg = SingleModValDialog (h, is_nontext, mod_type, 0);

  action_type = HiddenGroup (h, 0, 2, NULL);
  SetGroupSpacing (action_type, 10, 10);
  RadioButton (action_type, "For all sequences");
  radio_txt = (CharPtr) MemNew ((StringLen (radio_fmt) + StringLen (seq_id)) * sizeof (Char));
  if (radio_txt != NULL)
  {
    sprintf (radio_txt, radio_fmt, seq_id);
  }
  RadioButton (action_type, radio_txt);
  radio_txt = MemFree (radio_txt);
  SetValue (action_type, 1);
  
  copy_to_note_btn = CheckBox (h, "Copy original value to note", NULL);
  SetStatus (copy_to_note_btn, TRUE);
  
  c = HiddenGroup (h, 2, 0, NULL);
  accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p,
                              (HANDLE) new_value_dlg,
                              (HANDLE) action_type, 
                              (HANDLE) copy_to_note_btn,
                              (HANDLE) c, 
                              (HANDLE) instr_grp,
                              NULL);

  Show (w);
  Select (w);
  
  acd.cancelled = FALSE;
  acd.accepted = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.cancelled)
  {
    rval = FALSE;
  }
  else
  {
    new_value = DialogToPointer (new_value_dlg);
    action_type_val = GetValue (action_type);
    copy_to_note = GetStatus (copy_to_note_btn);
    
    for (seq_num = 0; iatep_new != NULL && seq_num < iatep_new->num_sequences; seq_num++)
    {
      if (action_type_val == 1 /* replace value for all sequences */
          || StringCmp (iatep_new->id_list [seq_num], seq_id) == 0)
      {
        iatep_new->title_list [seq_num] = ReplaceOneModifierValue (iatep_new->title_list [seq_num],
                                                                   orig_mod_name,
                                                                   orig_mod_value,
                                                                   new_value,
                                                                   is_nontext,
                                                                   copy_to_note); 
      }
    }
    for (seq_num = 0; iatep_current != NULL && seq_num < iatep_current->num_sequences; seq_num++)
    {
      if (action_type_val == 1 /* replace value for all sequences */
          || StringCmp (iatep_current->id_list [seq_num], seq_id) == 0)
      {
        iatep_current->title_list [seq_num] = ReplaceOneModifierValue (iatep_current->title_list [seq_num],
                                                                       orig_mod_name, 
                                                                       orig_mod_value,
                                                                       new_value,
                                                                       is_nontext,
                                                                       copy_to_note); 
      }
    }
    new_value = MemFree (new_value);
    
    rval = TRUE;
  }
    
  Remove (w);
  
  return rval;    
}

static void FindBadLocationInTitle (CharPtr seq_id, CharPtr title, ValNodePtr PNTR badlist)
{
  CharPtr value;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) || badlist == NULL
      || FindValuePairInDefLine ("location", title, NULL) == NULL)
  {
    return;
  }
  
  value = FindValueFromPairInDefline ("location", title);
  if (StringHasNoText (value))
  {
    ValNodeAddPointer (badlist, eModifierType_Location, BadValueNew (seq_id, "location", NULL));
  }
  else if (!IsValueInEnumAssoc (value, biosource_genome_simple_alist))
  {
    ValNodeAddPointer (badlist, eModifierType_Location, BadValueNew (seq_id, "location", value));
  }
  value = MemFree (value);
}

static void FindBadOriginInTitle (CharPtr seq_id, CharPtr title, ValNodePtr PNTR badlist)
{
  CharPtr value;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) || badlist == NULL
      || FindValuePairInDefLine ("origin", title, NULL) == NULL)
  {
    return;
  }
  
  value = FindValueFromPairInDefline ("origin", title);
  if (StringHasNoText (value))
  {
    ValNodeAddPointer (badlist, eModifierType_Origin, BadValueNew (seq_id, "origin", NULL));
  }
  else if (!IsValueInEnumAssoc (value, biosource_origin_alist))
  {
    ValNodeAddPointer (badlist, eModifierType_Origin, BadValueNew (seq_id, "origin", value));
  }
  value = MemFree (value);
}

static void FindBadTopologyInTitle (CharPtr seq_id, CharPtr title, ValNodePtr PNTR badlist)
{
  CharPtr value;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) || badlist == NULL
      || FindValuePairInDefLine ("topology", title, NULL) == NULL)
  {
    return;
  }
  
  value = FindValueFromPairInDefline ("topology", title);
  if (StringHasNoText (value))
  {
    ValNodeAddPointer (badlist, eModifierType_Topology, BadValueNew (seq_id, "topology", NULL));
  }
  else if (!IsValueInEnumAssoc (value, topology_nuc_alist))
  {
    ValNodeAddPointer (badlist, eModifierType_Topology, BadValueNew (seq_id, "topology", value));
  }
  value = MemFree (value);
}

static void FindBadMolTypeInTitle (CharPtr seq_id, CharPtr title, ValNodePtr PNTR badlist)
{
  CharPtr value;
  Int4    moltype;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) || badlist == NULL
      || FindValuePairInDefLine ("moltype", title, NULL) == NULL)
  {
    return;
  }

  value = FindValueFromPairInDefline ("moltype", title);
  moltype = MolTypeFromString (value);
  if (moltype == 0)
  {
    ValNodeAddPointer (badlist, eModifierType_MolType, BadValueNew (seq_id, "moltype", value));
  }
  value = MemFree (value);  
}

static void FindBadMoleculeInTitle (CharPtr seq_id, CharPtr title, ValNodePtr PNTR badlist)
{
  CharPtr value;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) || badlist == NULL
      || FindValuePairInDefLine ("molecule", title, NULL) == NULL)
  {
    return;
  }

  value = FindValueFromPairInDefline ("molecule", title);
  if (StringICmp (value, "dna") != 0 && StringICmp (value, "rna") != 0)
  {
    ValNodeAddPointer (badlist, eModifierType_Molecule, BadValueNew (seq_id, "molecule", value));
  }
  value = MemFree (value);    
}

static void 
FindBadGeneticCodeInTitle 
(CharPtr seq_id,
 CharPtr title, 
 CharPtr mod_name, 
 ValNodePtr PNTR badlist)
{
  CharPtr value;
  Int4    gcode;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) 
      || StringHasNoText (mod_name)
      || badlist == NULL
      || FindValuePairInDefLine (mod_name, title, NULL) == NULL)
  {
    return;
  }

  value = FindValueFromPairInDefline (mod_name, title);
  gcode = GeneticCodeFromString (value);
  if (gcode == 0)
  {
    ValNodeAddPointer (badlist, GetModifierType (mod_name), BadValueNew (seq_id, mod_name, value));
  }
  value = MemFree (value);
}

static void FindBadNonTextValueInTitle (CharPtr seq_id, CharPtr title, CharPtr mod_name, ValNodePtr PNTR badlist)
{
  CharPtr value;
  
  if (StringHasNoText (seq_id) || StringHasNoText (title) 
      || StringHasNoText (mod_name)
      || badlist == NULL
      || FindValuePairInDefLine (mod_name, title, NULL) == NULL)
  {
    return;
  }

  value = FindValueFromPairInDefline (mod_name, title);
  if (!StringHasNoText (value) && StringICmp (value, "TRUE") != 0 && StringICmp (value, "FALSE") != 0)
  {
    ValNodeAddPointer (badlist, eModifierType_SourceQual, BadValueNew (seq_id, mod_name, value));
  }
  value = MemFree (value);
}

static void FindBadValuesInTitle (CharPtr seq_id, CharPtr title, ValNodePtr PNTR badlist)
{  
  if (StringHasNoText (seq_id) || StringHasNoText (title) || badlist == NULL)
  {
    return;
  }
  
  FindBadLocationInTitle (seq_id, title, badlist);
  FindBadOriginInTitle (seq_id, title, badlist);
  FindBadGeneticCodeInTitle (seq_id, title, "gcode", badlist);
  FindBadGeneticCodeInTitle (seq_id, title, "mgcode", badlist);
  FindBadGeneticCodeInTitle (seq_id, title, "genetic_code", badlist);
  FindBadMolTypeInTitle (seq_id, title, badlist);
  FindBadMoleculeInTitle (seq_id, title, badlist);
  FindBadTopologyInTitle (seq_id, title, badlist);
  
  /* check nontext modifiers */
  FindBadNonTextValueInTitle (seq_id, title, "transgenic", badlist);
  FindBadNonTextValueInTitle (seq_id, title, "germline", badlist);
  FindBadNonTextValueInTitle (seq_id, title, "metagenomic", badlist);
  FindBadNonTextValueInTitle (seq_id, title, "environmental-sample", badlist);
  FindBadNonTextValueInTitle (seq_id, title, "rearranged", badlist);
}

static void FindBadValuesInIDsAndTitles (IDAndTitleEditPtr iatep, ValNodePtr PNTR bad_list)
{
  Int4 seq_num;
  
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    FindBadValuesInTitle (iatep->id_list [seq_num], iatep->title_list [seq_num], bad_list);
  }
}

static CharPtr 
GetIDAndTitleErrorMessage 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current, 
 Int4              seq_num,
 Boolean           has_dups,
 Boolean           has_missing,
 Boolean           seqid_edit_phase,
 Boolean           has_bracket,
 Boolean           has_unrec_mods,
 Boolean           is_nuc)
{
  ValNodePtr unrec_mod_list = NULL, bad_value_list;  
  Boolean    is_dup, has_id_bracket = FALSE, has_space = FALSE;
  Int4       msg_num = 0;
  CharPtr    err_msg = "";
  BadValuePtr bvp;

  /* get appropriate error message */
  unrec_mod_list = NULL;
  bad_value_list = NULL;
    
  /* determine whether this is a duplicate */
  if (is_nuc)
  {
    is_dup = IsDuplicateEditID (iatep_new, seq_num, iatep_current);
  }
  else
  {
    is_dup = FALSE;
  }
  
  /* look for bracket in ID */
  if (StringChr (iatep_new->id_list [seq_num], '['))
  {
    has_id_bracket = TRUE;   
  }

  /* look for space in ID */
  has_space = StringHasSpace (iatep_new->id_list [seq_num]);
  
  if (has_dups || has_missing || has_id_bracket || has_space)
  {
    if (StringHasNoText (iatep_new->id_list [seq_num]))
    {
      err_msg = "Missing ID";
    }
    else if (is_dup)
    {
      err_msg = "Duplicate ID";
    }
    else if (has_id_bracket)
    {
      err_msg = "Bracket in ID";
    }
    else if (has_space)
    {
      err_msg = "Space in ID";
    }
  }
  else if (seqid_edit_phase)
  {
    err_msg = "";
  }
  else if (has_bracket)
  {
    msg_num = DetectBadBracketing (iatep_new->title_list [seq_num]);
    /* we have bracketing problems */
    switch (msg_num)
    {
      case BRACKET_ERR_MISMATCHED_BRACKETS:
        err_msg = "Mismatched []";
        break;
      case BRACKET_ERR_MISSING_EQUALS:
        err_msg = "Missing '='";
        break;
      case BRACKET_ERR_MULT_EQUALS:
        err_msg = "Too many '='";
        break;
      case BRACKET_ERR_NO_MOD_NAME:
        err_msg = "Missing name";
        break;
      case BRACKET_ERR_MISMATCHED_QUOTES:
        err_msg = "Mismatched \" or '";
        break;
    }
  }
  else if (has_unrec_mods)
  {
    AddUnrecognizedModifiersForOneDefinitionLine (iatep_new->title_list [seq_num],
                                                  &unrec_mod_list,
                                                  is_nuc);
    if (unrec_mod_list != NULL)
    {
      err_msg = unrec_mod_list->data.ptrvalue;
    }
  }
  else
  {
    FindBadValuesInTitle (iatep_new->id_list [seq_num], iatep_new->title_list [seq_num], &bad_value_list);
    if (bad_value_list != NULL && (bvp = bad_value_list->data.ptrvalue) != NULL)
    {
      err_msg = bvp->mod_name;
    }
  }
  
  err_msg = StringSave (err_msg);
  ValNodeFreeData (unrec_mod_list);
  BadValueListFree (bad_value_list);
  return err_msg;
}

static CharPtr GetTagListErrValueForSeqNum (TagListPtr tlp, Int4 seq_num)
{
  Char       seq_str [15];
  ValNodePtr vnp;
  Int4       row_num;
  CharPtr    pos_str;
  
  if (tlp == NULL)
  {
    return NULL;
  }
  
  sprintf (seq_str, "%d", seq_num + 1);
  for (vnp = tlp->vnp, row_num = 0; vnp != NULL; vnp = vnp->next, row_num++)
  {
    pos_str = GetTagListValueEx (tlp, row_num, 1);
    if (StringCmp (pos_str, seq_str) == 0)
    {
      pos_str = MemFree (pos_str);
      return GetTagListValueEx (tlp, row_num, 0);
    }
    pos_str = MemFree (pos_str);
  }
  return NULL;
}

/* This function displays information from a list of sequence IDs and titles
 * in a TagList dialog.
 * The TagList dialog has four columns: Error, Position, Sequence ID, and Title.
 */
static void 
UpdateIdAndTitleEditDialog 
(DialoG            d,
 IDAndTitleEditPtr iatep_new, 
 IDAndTitleEditPtr iatep_current,
 Boolean           seqid_edit_phase, 
 Boolean           show_all,
 Boolean           is_nuc)
{
  TagListPtr tlp;
  Int4       seq_num, len;
  CharPtr    str;
  CharPtr    str_format = "%s\t%d\t%s\t%s\n";
  Int4       num_shown;
  CharPtr    err_msg, old_err_msg;
  ValNodePtr taglist_data = NULL;
  Int4       row_to_show, row_to_hide;
  Boolean    has_dups, has_missing, has_bracket, has_space, has_unrec_mods = FALSE;
  ValNodePtr unrec_mods = NULL;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL || iatep_new == NULL || iatep_new->num_sequences == 0)
  {
    return;
  }
  
  has_dups = EditHasDuplicateIDs (iatep_new, iatep_current);
  has_missing = HasMissingIDs (iatep_new) || HasMissingIDs (iatep_current);
  has_bracket = EditNeedsBracketingFixes (iatep_new) 
                || EditNeedsBracketingFixes (iatep_current);
  has_space = EditHasSpaceInIDs (iatep_new) || EditHasSpaceInIDs (iatep_current);

  if (!has_bracket)
  {
    unrec_mods = ListUnrecognizedModifiers (iatep_new, iatep_current, is_nuc);
  }
  if (unrec_mods != NULL)
  {
    has_unrec_mods = TRUE;
    unrec_mods = ValNodeFreeData (unrec_mods);
  }
  
  num_shown = 0;
  for (seq_num = 0; seq_num < iatep_new->num_sequences; seq_num++)
  {
    err_msg = GetIDAndTitleErrorMessage (iatep_new, iatep_current, 
                                         seq_num, has_dups,
                                         has_missing, seqid_edit_phase,
                                         has_bracket,
                                         has_unrec_mods, is_nuc);
    if (seqid_edit_phase && StringHasNoText (err_msg))
    {
      old_err_msg = GetTagListErrValueForSeqNum (tlp, seq_num);
      if (StringCmp (old_err_msg, "Duplicate ID") == 0
          || StringCmp (old_err_msg, "Missing ID") == 0
          || StringCmp (old_err_msg, "Space in ID") == 0
          || StringCmp (old_err_msg, "Fixed") == 0)
      {
        err_msg = MemFree (err_msg);
        err_msg = StringSave ("Fixed");
      }
      old_err_msg = MemFree (old_err_msg);
    }
                                         
    if (StringHasNoText (err_msg) && !show_all)
    {
      err_msg = MemFree (err_msg);
      continue;
    }

    len = StringLen (str_format) + StringLen (err_msg) + 20
                     + StringLen (iatep_new->id_list [seq_num]) 
                     + StringLen (iatep_new->title_list [seq_num]);
    str = MemNew (len * sizeof (Char));
    if (str != NULL) 
    {
      sprintf (str, str_format, 
               err_msg,
               seq_num + 1,
               StringHasNoText (iatep_new->id_list [seq_num]) ? "" : iatep_new->id_list [seq_num],
               StringHasNoText (iatep_new->title_list [seq_num]) ? "" : iatep_new->title_list [seq_num]);               
      ValNodeAddPointer (&taglist_data, 0, StringSave (str));
    }
    err_msg = MemFree (err_msg);
    num_shown ++;    
  }
  
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET); 
  tlp->vnp = taglist_data;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  tlp->max = MAX ((Int2) 0, (Int2) (num_shown - tlp->rows));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  CorrectBarMax (tlp->left_bar, tlp->max);
  CorrectBarPage (tlp->left_bar, tlp->rows - 1, tlp->rows - 1);
  for (row_to_show = 0; row_to_show < MIN (num_shown, tlp->rows); row_to_show ++)
  {
    SafeShow (tlp->control [row_to_show * MAX_TAGLIST_COLS + 2]);
    SafeShow (tlp->control [row_to_show * MAX_TAGLIST_COLS + 3]);
  }
  if (tlp->max > 0) {
    SafeShow (tlp->bar);
    SafeShow (tlp->left_bar);
  } else {
    SafeHide (tlp->bar);
    SafeHide (tlp->left_bar);
    for (row_to_hide = num_shown; row_to_hide < tlp->rows; row_to_hide ++)
    {
      SafeHide (tlp->control [row_to_hide * MAX_TAGLIST_COLS + 2]);
      SafeHide (tlp->control [row_to_hide * MAX_TAGLIST_COLS + 3]);
    }    
  }
}

/* This function copies the contents of a TagList dialog into a list
 * of Sequence IDs and titles.  The first two columns of the TagList
 * dialog, Error and Position, are ignored.
 */
static void UpdateIdAndTitleData (DialoG d, IDAndTitleEditPtr iatep)
{
  CharPtr      str;
  Int4         num_rows, row_num, seq_pos;
  TagListPtr   tlp;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL || iatep == NULL)
  {
    return;
  }
  
  num_rows = ValNodeLen (tlp->vnp);
  for (row_num = 0; row_num < num_rows; row_num++)
  {
    /* get position for this sequence */
    str = GetTagListValueEx (tlp, row_num, 1);
    seq_pos = atoi (str);
    str = MemFree (str);
    if (seq_pos < 1 || seq_pos > iatep->num_sequences)
    {
      continue;
    }
    seq_pos --;
    
    /* collect ID */
    iatep->id_list [seq_pos] = MemFree (iatep->id_list [seq_pos]);
    iatep->id_list [seq_pos] = GetTagListValueEx (tlp, row_num, 2);
    TrimSpacesAroundString (iatep->id_list [seq_pos]);
    
    /* collect title */
    iatep->title_list [seq_pos] = MemFree (iatep->title_list [seq_pos]);
    iatep->title_list [seq_pos] = GetTagListValueEx (tlp, row_num, 3);
  }    
}

static void SetIDAndTitleEditDialogErrorColumn
(DialoG            d,
 IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 Boolean           seqid_edit_phase,
 Boolean           is_nuc)
{
  TagListPtr tlp;
  ValNodePtr vnp;
  Int4       seq_num, seq_pos;
  Int4       row_num;
  CharPtr    err_msg, old_err_msg;
  Boolean    has_dups, has_missing, has_bracket, has_unrec_mods = FALSE, has_space;
  CharPtr    str;
  ValNodePtr unrec_mods = NULL;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL || iatep_new == NULL || iatep_new->num_sequences == 0)
  {
    return;
  }
  
  has_dups = EditHasDuplicateIDs (iatep_new, iatep_current);
  has_missing = HasMissingIDs (iatep_new) || HasMissingIDs (iatep_current);
  has_bracket = EditNeedsBracketingFixes (iatep_new) 
                || EditNeedsBracketingFixes (iatep_current);
  has_space = EditHasSpaceInIDs (iatep_new) || EditHasSpaceInIDs (iatep_current);
  
  if (!has_bracket)
  {
    unrec_mods = ListUnrecognizedModifiers (iatep_new, iatep_current, is_nuc);
  }
  if (unrec_mods != NULL)
  {
    has_unrec_mods = TRUE;
    unrec_mods = ValNodeFreeData (unrec_mods);
  }
  
  for (row_num = 0, vnp = tlp->vnp; vnp != NULL; row_num++, vnp = vnp->next)
  {
    /* get position for this sequence */
    str = GetTagListValueEx (tlp, row_num, 1);
    if (str == NULL)
    {
      continue;
    }
    seq_pos = atoi (str);
    str = MemFree (str);
    if (seq_pos < 1 || seq_pos > iatep_new->num_sequences)
    {
      continue;
    }
    seq_num = seq_pos - 1;

    /* get appropriate error message */
    err_msg = GetIDAndTitleErrorMessage (iatep_new, iatep_current, 
                                         seq_num, has_dups,
                                         has_missing, seqid_edit_phase,
                                         has_bracket,
                                         has_unrec_mods,
                                         is_nuc);
                                         
    if (seqid_edit_phase && StringHasNoText (err_msg))
    {
      old_err_msg = GetTagListValueEx (tlp, row_num, 0);
      if (StringCmp (old_err_msg, "Duplicate ID") == 0
          || StringCmp (old_err_msg, "Missing ID") == 0
          || StringCmp (old_err_msg, "Space in ID") == 0
          || StringCmp (old_err_msg, "Fixed") == 0)
      {
        err_msg = MemFree (err_msg);
        err_msg = StringSave ("Fixed");
      }
      old_err_msg = MemFree (old_err_msg);
    }                                                                                  
    
    SetTagListValue (tlp, row_num, 0, err_msg);

    err_msg = MemFree (err_msg);    
  }
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
}

static void ClearIDAndTitleEditDialogErrorColumn (DialoG d)
{
  TagListPtr tlp;
  ValNodePtr vnp;
  Int4       row_num;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL)
  {
    return;
  }
  
  for (row_num = 0, vnp = tlp->vnp; vnp != NULL; row_num++, vnp = vnp->next)
  {    
    SetTagListValue (tlp, row_num, 0, "");
  }
}


static void 
UpdateIDAndTitleEditDialogErrorColumns
(DialoG            d_new,
 DialoG            d_current,
 IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 Boolean           seqid_edit_phase,
 Boolean           is_nuc)
{
  if (d_new == NULL)
  {
    iatep_new = NULL;
  }
  else
  {
    iatep_new = IDAndTitleEditCopy (iatep_new);
    UpdateIdAndTitleData (d_new, iatep_new);
  }
  
  if (d_current == NULL)
  {
    iatep_current = NULL;
  }
  else
  {
    iatep_current = IDAndTitleEditCopy (iatep_current);
    UpdateIdAndTitleData (d_current, iatep_current);
  }
  
  SetIDAndTitleEditDialogErrorColumn (d_new, iatep_new, iatep_current, 
                                      seqid_edit_phase, is_nuc);
  SetIDAndTitleEditDialogErrorColumn (d_current, iatep_current, iatep_new,
                                      seqid_edit_phase, is_nuc);
  iatep_new = IDAndTitleEditFree (iatep_new);
  iatep_current = IDAndTitleEditFree (iatep_current);
}


typedef struct unrecmods
{
  DialoG PNTR unrec_dlg;
  ValNodePtr unrecognized_list;
  Int4        num_unrecognized;  
  ButtoN      accept_btn;
} UnrecModsData, PNTR UnrecModsPtr;

static void CleanupUnrecMods (GraphiC g, VoidPtr data)

{
  UnrecModsPtr ump;

  ump = (UnrecModsPtr) data;
  if (ump != NULL)
  {
    ump->unrecognized_list = ValNodeFreeData (ump->unrecognized_list);
  }
  MemFree (ump);
}

static void SetUnrecAccept (Pointer userdata)
{
  UnrecModsPtr ump;
  ValNodePtr   vnp;
  Int4       repl_num;
  Boolean      ok_to_accept = TRUE;

  ump = (UnrecModsPtr) userdata;
  if (ump == NULL)
  {
    return;
  }
  
  for (repl_num = 0; repl_num < ump->num_unrecognized && ok_to_accept && repl_num < 3; repl_num++)
  {
    vnp = DialogToPointer (ump->unrec_dlg [repl_num]);
    if (vnp == NULL)
    {
      ok_to_accept = FALSE;
    }
    vnp = ValNodeFreeData (vnp);
  }
  if (ok_to_accept)
  {
    Enable (ump->accept_btn);
  }
  else
  {
    Disable (ump->accept_btn);
  }
  
}

static void ReplaceThreeUnrecognizedModifiers (IDAndTitleEditPtr iatep, UnrecModsPtr ump)
{
  ValNodePtr vnp, repl_vnp;
  Int4       repl_num, seq_num;
  CharPtr    repl_name;

  if (iatep == NULL || ump == NULL)
  {
    return;
  }
  
  for (repl_num = 0, vnp = ump->unrecognized_list;
       repl_num < ump->num_unrecognized && vnp != NULL && repl_num < 3;
       repl_num++, vnp = vnp->next)
  {
    if (StringHasNoText (vnp->data.ptrvalue))
    {
      continue;
    }
    repl_vnp = DialogToPointer (ump->unrec_dlg [repl_num]);
    if (repl_vnp == NULL)
    {
      continue;
    }
    repl_name = SourceQualValNodeName (repl_vnp);
    for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
    {
      iatep->title_list [seq_num] = ReplaceOneModifierName (iatep->title_list [seq_num],
                                                            vnp->data.ptrvalue, 
                                                            repl_name);
    }
    repl_name = MemFree (repl_name);
    repl_vnp = ValNodeFreeData (repl_vnp);
  }  
}

static Boolean 
FixThreeUnrecognizedModifiers 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 ValNodePtr        unrecognized_list)
{
  ValNodePtr vnp;
  Int4       repl_num;
  WindoW     w;
  ValNodePtr mod_choices;
  GrouP      h, g, k, c;
  PrompT       p;
  UnrecModsPtr ump;
  ModalAcceptCancelData acd;
  Boolean               rval = FALSE;
  ButtoN                b;
  
  if (unrecognized_list == NULL || (iatep_new == NULL && iatep_current == NULL))
  {
    return FALSE;
  }
  
  ump = (UnrecModsPtr) MemNew (sizeof (UnrecModsData));
  if (ump == NULL)
  {
    return FALSE;
  }
  
  ump->unrecognized_list = unrecognized_list;
  ump->num_unrecognized = ValNodeLen(ump->unrecognized_list);

  ump->unrec_dlg = (DialoG PNTR) MemNew (sizeof (DialoG) * ump->num_unrecognized);

  w = MovableModalWindow (-20, -13, -10, -10, "Choose Valid Modifiers", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  SetObjectExtra (w, ump, CleanupUnrecMods);
  
  p = StaticPrompt (h, "Please choose a valid modifier name to replace these invalid names.",
                    0, 0, programFont, 'l');

  g = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  for (repl_num = 0, vnp = ump->unrecognized_list;
       repl_num < ump->num_unrecognized && repl_num < 3 && vnp != NULL;
       repl_num++, vnp = vnp->next)
  {
    k = HiddenGroup (g, 2, 0, NULL);
    SetGroupSpacing (k, 2, 2);
    mod_choices = GetFastaModifierList (TRUE, TRUE);
    StaticPrompt (k, vnp->data.ptrvalue, 0, 0, programFont, 'l');
    ump->unrec_dlg [repl_num] = ValNodeSelectionDialog (k, mod_choices, 6,
                                          SourceQualValNodeName,
                                          ValNodeSimpleDataFree,
                                          SourceQualValNodeDataCopy,
                                          SourceQualValNodeMatch,
                                          "modifier",
                                          SetUnrecAccept, ump, FALSE);

  }
  
  c = HiddenGroup (h, 2, 0, NULL);
  ump->accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (ump->accept_btn, &acd, NULL);
  Disable (ump->accept_btn);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p,
                              (HANDLE) g, 
                              (HANDLE) c, 
                              NULL);

  Show (w);
  Select (w);
  
  acd.cancelled = FALSE;
  acd.accepted = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.cancelled)
  {
    rval = FALSE;
  }
  else
  {
    ReplaceThreeUnrecognizedModifiers (iatep_new, ump);
    ReplaceThreeUnrecognizedModifiers (iatep_current, ump);
    
    rval = TRUE;
  }
    
  Remove (w);
  
  return rval;    
}

static Boolean 
FixAllUnrecognizedModifiers 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 SeqIdEditPtr siep)
{
  ValNodePtr unrecognized_list;
  Boolean    rval = TRUE;
  Boolean    show_all;
  
  if (siep == NULL)
  {
    return FALSE;
  }
  if (GetValue (siep->show_all_grp) == 2)
  {
    show_all = TRUE;
  }
  else
  {
    show_all = FALSE;
  }
  
  unrecognized_list = ListUnrecognizedModifiers (iatep_new, iatep_current, siep->is_nuc);
  while (unrecognized_list != NULL && rval)
  {
    rval = FixThreeUnrecognizedModifiers (iatep_new, iatep_current, unrecognized_list);
    unrecognized_list = ListUnrecognizedModifiers (iatep_new, iatep_current, siep->is_nuc);
    UpdateIdAndTitleEditDialog (siep->new_dlg, siep->iatep_new, siep->iatep_current, 
                                siep->seqid_edit_phase, show_all, siep->is_nuc);
    UpdateIdAndTitleEditDialog (siep->current_dlg, siep->iatep_current, siep->iatep_new, 
                                siep->seqid_edit_phase, show_all, siep->is_nuc);
    ShowErrorInstructions (siep);
  }  
  return rval;
}


static ParData     extendedIDParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData     extendedIDColFmt[] = 
  {
    {0, 0, 40, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 40, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 40, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 40, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };

static Boolean AnyBracketsInIDs (IDAndTitleEditPtr iatep)
{
  Int4 seq_num;
  Boolean rval = FALSE;
  
  if (iatep == NULL)
  {
    return FALSE;
  }
  for (seq_num = 0; seq_num < iatep->num_sequences && !rval; seq_num++)
  {
    if (StringChr (iatep->id_list [seq_num], '['))
    {
      rval = TRUE;
    }
  }
  return rval;
}

static Boolean AnyIDCorrectionsToList 
(IDAndTitleEditPtr iatep,
 IDAndTitleEditPtr iatep_current,
 BoolPtr           space_corr,
 BoolPtr           bracket_corr)
{
  IDAndTitleEditPtr suggested;
  Boolean           any_to_show = FALSE;
  Int4              seq_num;

  if (iatep == NULL || space_corr == NULL || bracket_corr == NULL)
  {
    return FALSE;
  }

  suggested = SuggestCorrectionForLocalIDs (iatep, iatep_current);
  if (suggested == NULL || iatep->num_sequences != suggested->num_sequences)
  {
    suggested = IDAndTitleEditFree (suggested);
    return FALSE;
  }
  
  for (seq_num = 0;
       seq_num < iatep->num_sequences && (!any_to_show || !*space_corr || !*bracket_corr);
       seq_num++)
  {
    if (! StringHasNoText (suggested->id_list [seq_num])
        && ! StringHasNoText (iatep->title_list [seq_num])
        && StringCmp (iatep->id_list [seq_num], suggested->id_list [seq_num]) != 0)
    {
      any_to_show = TRUE;
      if (StringChr (iatep->id_list [seq_num], '[') != NULL)
      {
        *bracket_corr = TRUE;
      }
      else
      {
        *space_corr = TRUE;
      }
    }
  }
  suggested = IDAndTitleEditFree (suggested);
  
  return any_to_show;  
}

static Boolean ListIDCorrections
(IDAndTitleEditPtr iatep,
 IDAndTitleEditPtr iatep_current,
 CharPtr           str,
 DoC               doc)
{
  IDAndTitleEditPtr suggested;
  CharPtr           doc_txt;
  CharPtr           doc_txt_fmt = "%s%d\t%s\t%s\t%s\n";
  Int4              len;
  Boolean           any_to_show = FALSE;
  Int4              seq_num;
  
  
  if (iatep == NULL || doc == NULL)
  {
    return FALSE;
  }

  suggested = SuggestCorrectionForLocalIDs (iatep, iatep_current);
  if (suggested == NULL || iatep->num_sequences != suggested->num_sequences)
  {
    suggested = IDAndTitleEditFree (suggested);
    return FALSE;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    if (StringHasNoText (suggested->id_list [seq_num])
        || StringHasNoText (iatep->title_list [seq_num])
        || StringCmp (iatep->id_list [seq_num], suggested->id_list [seq_num]) == 0)
    {
      continue;
    }
   
    len = StringLen (doc_txt_fmt) 
                     + StringLen (suggested->id_list [seq_num]) 
                     + StringLen (iatep->id_list[seq_num])
                     + StringLen (str)
                     + 15
                     + StringLen (iatep->title_list [seq_num]);
    doc_txt = (CharPtr) MemNew (len * sizeof (Char));
    if (doc_txt != NULL)
    {
      sprintf (doc_txt, doc_txt_fmt, 
               str == NULL ? "" : str, 
               seq_num + 1,
               suggested->id_list [seq_num], 
               iatep->id_list [seq_num] == NULL ? "" : iatep->id_list [seq_num],
               iatep->title_list [seq_num]);
      AppendText (doc, doc_txt, &extendedIDParFmt, extendedIDColFmt, programFont);
      doc_txt = MemFree (doc_txt);
      any_to_show = TRUE;
    }
  }
  
  suggested = IDAndTitleEditFree (suggested);
  
  return any_to_show;
}

static Boolean 
ShowExtendedIDCorrections 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current, 
 DoC               doc)
{
  Boolean     any_to_show = FALSE, space_corr = FALSE, bracket_corr = FALSE;
  RecT        r;
  
  if (doc == NULL || iatep_new == NULL)
  {
    return FALSE;
  }
  
  if (! AnyIDCorrectionsToList (iatep_new, iatep_current, &space_corr, &bracket_corr)
      && ! AnyIDCorrectionsToList (iatep_current, iatep_new, &space_corr, &bracket_corr))
  {
    return FALSE;
  }

  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  extendedIDColFmt[0].pixWidth = (r.right - r.left) / 10;
  extendedIDColFmt[1].pixWidth = (r.right - r.left) / 4;
  extendedIDColFmt[2].pixWidth = (r.right - r.left) / 4;
  extendedIDColFmt[3].pixWidth = (r.right - r.left) 
                                  - extendedIDColFmt[0].pixWidth
                                  - extendedIDColFmt[1].pixWidth
                                  - extendedIDColFmt[2].pixWidth;
  
  if (space_corr)
  {
    AppendText (doc, "Your sequence IDs are not unique.  Did you try to put spaces in your sequence IDs?  This is not allowed.\n", NULL, NULL, programFont);
  }
  if (AnyBracketsInIDs (iatep_new) || AnyBracketsInIDs (iatep_current))
  {
    AppendText (doc, "Did you forget to put spaces between your sequence IDs and your titles?  This is not allowed.\n", NULL, NULL, programFont);
  }
    
  AppendText (doc, "\nPosition\tSuggested ID\tOriginal ID\tOriginal Title\n", &extendedIDParFmt, extendedIDColFmt, programFont);

  if (iatep_current == NULL)
  {
    any_to_show = ListIDCorrections (iatep_new, iatep_current, "", doc);
  }
  else
  {
    any_to_show = ListIDCorrections (iatep_new, iatep_current, "new:", doc);
    any_to_show |= ListIDCorrections (iatep_current, iatep_new, "existing:", doc);    
  }
    
  return any_to_show;
}

static Boolean SequenceIDsHaveNonFixableBrackets 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current, 
 DoC               doc)
{
  if (doc == NULL)
  {
    return FALSE;
  }
  
  if (! AnyBracketsInIDs (iatep_new) && ! AnyBracketsInIDs (iatep_current))
  {
    return FALSE;
  }

  Reset (doc);
  AppendText (doc, "Your sequence IDs contain brackets ('[' and/or ']').  This is not allowed.\n", NULL, NULL, programFont);
    
  return TRUE;  
}

/* This section of code is used to display bracketing errors in an AutonomousPanel.
 * The panel has a frozen title row and a frozen column with sequence IDs in it.
 * The panel scrolling affects only the sequence titles.
 * There are two rows for each sequence.  The top row displays the original title;
 * the bottom row displays the suggested bracketing corrections.  The differences
 * between the original and suggested titles will be colored in red.
 * Both rows in a pair will have the same background color which alternates between gray
 * and white for each pair.
 * Clicking on a sequence will scroll to the next difference for that pair.
 */
typedef struct diffdlg
{
  DIALOG_MESSAGE_BLOCK
  PaneL             pnl;
  IDAndTitleEditPtr new_original;
  IDAndTitleEditPtr new_suggested;
  IDAndTitleEditPtr existing_original;
  IDAndTitleEditPtr existing_suggested;
  FonT              display_font;
  Int4              char_width;
  Int4              descent;
  Int4              num_header_rows;
  Int4              max_title_length;
  Int4              max_id_length;
  Int4              table_inset;
  
} DiffDlgData, PNTR DiffDlgPtr;

static void DrawDiffDlgExplanation (Int4 x, Int4 y, Int4 descent, Int4 win_width)
{
  RecT rct;
  
  /* draw explanation rows */
  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  PaintStringEx ("Some of your titles have bracketing errors.", x, y);
  y += stdLineHeight;
  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  PaintStringEx ("Double-click on 'original' to scroll to the next error for that pair.", x, y);  
}

static void DrawDiffDlgTitle (Int4 x, Int4 y, Int4 char_width, Int4 descent, Int4 max_id_length, Int4 win_width)
{
  RecT rct;
  
  /* draw title row */
  DkGray ();
  InvertColors ();
  White ();
  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  PaintStringEx ("Sequence ID", x, y);
  x += (max_id_length + 1) * char_width;
  x += 10 * char_width;
  PaintStringEx ("Title", x, y);
  InvertColors ();
  Black ();
}

static void 
PaintColorizedString 
(Int4       x,
 Int4       y,
 Int4       char_width,
 CharPtr    paintstring, 
 Int4       string_offset,
 Boolean    shade,
 ValNodePtr diff_list,
 Int4       diff_choice)
{
  CharPtr    cp;
  Char       buf [2];
  ValNodePtr diff_vnp;
  
  if (shade)
  {
    White ();
  }
  else
  {
    Black ();
  }
  if (paintstring == NULL || (Int4) StringLen (paintstring) <= string_offset)
  {
    PaintStringEx (" ", x, y);
  }
  else
  {
    cp = paintstring + string_offset;
    diff_vnp = diff_list;
    buf [1] = 0;
    
    while (*cp != 0)
    {
      while (diff_vnp != NULL && (diff_vnp->choice != diff_choice || diff_vnp->data.intvalue < string_offset))
      {
        diff_vnp = diff_vnp->next;
      }
      if (diff_vnp != NULL && diff_vnp->choice == diff_choice && string_offset == diff_vnp->data.intvalue)
      {
        Red ();
      }
    
      buf [0] = *cp;
      PaintStringEx (buf, x, y);
      x += char_width;
      string_offset++;
      cp++;
    
      if (shade)
      {
        White ();
      }
      else
      {
        Black ();
      }
    }
  }  
}

static void 
DrawDiffDlgRow 
(Int4 x, 
 Int4 y, 
 Int4 char_width, 
 Int4 descent, 
 Int4 max_id_length, 
 Int4 win_width,
 CharPtr id_str,
 CharPtr title_str,
 Int4       offset,
 ValNodePtr diff_list,
 Int4       choice_num,
 Boolean    shade)
{
  RecT rct;
  PoinT      pt1, pt2;

  if (shade)
  {
    Gray ();
    InvertColors ();
    White ();
  }
  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  if (id_str != NULL)
  {
    PaintStringEx (id_str, x, y);
  }
  x += (max_id_length + 1) * char_width;
  if (choice_num == 1)
  {
    PaintStringEx (" original", x, y);
  }
  else
  {
    PaintStringEx ("suggested", x, y);
  }
  x += 9 * char_width;
  pt1.x = x + 2;
  pt2.x = x + 2;
  pt1.y = y + descent;
  pt2.y = y - stdLineHeight + descent;
  Black ();
  DrawLine (pt1, pt2);

  x += char_width;
  
  PaintColorizedString (x, y, char_width, title_str, offset, shade, diff_list,
                        choice_num);
  
  if (shade)
  {
    InvertColors ();
    Black ();  
  }
}

/* This function produces a ValNode list of integers.
 * The choice value (1 or 2) indicates whether the difference is in string 1
 * or string 2; the integer value indicates the offset in that string of the
 * difference.
 * Space characters are ignored when computing differences.
 * The output from this function is used for displaying the differences between
 * the original definition line and a suggested bracketing correction.
 */
static ValNodePtr GetTextDifferences (CharPtr str1, CharPtr str2)
{
  ValNodePtr diff_list = NULL;
  CharPtr    cp1, cp2, diff_end1, diff_end2;
  Int4       offset1, offset2, j;
  
  if (str1 == NULL && str2 == NULL)
  {
    return NULL;
  }
  cp1 = str1;
  cp2 = str2;
  
  offset1 = 0;
  offset2 = 0;
  
  while (cp1 != NULL || cp2 != NULL)
  {
    /* skip over spaces in cp1 */
    while (cp1 != NULL && *cp1 != 0 && isspace (*cp1))
    {
      cp1 ++;
      offset1 ++;
    }
    if (cp1 != NULL && *cp1 == 0)
    {
      cp1 = NULL;
    }
    /* skip over spaces in cp2 */
    while (cp2 != NULL && *cp2 != 0 && isspace (*cp2))
    {
      cp2 ++;
      offset2 ++;
    }
    if (cp2 != NULL && *cp2 == 0)
    {
      cp2 = NULL;
    }
    
    if (cp1 == NULL && cp2 == NULL)
    {
      /* both NULL, do nothing */
    }
    else if (cp1 == NULL)
    {
      ValNodeAddInt (&diff_list, 2, offset2);
    }
    else if (cp2 == NULL)
    {
      ValNodeAddInt (&diff_list, 1, offset1);
    }
    else
    {
      if (*cp1 != *cp2)
      {
        if ((diff_end1 = StringSearch (cp1, cp2)) != NULL)
        {
          while (diff_end1 != cp1)
          {
            if (!isspace (*cp1))
            {
              ValNodeAddInt (&diff_list, 1, offset1);
            }
            offset1++;
            cp1++;
          }
        }
        else if ((diff_end2 = StringSearch (cp2, cp1))!= NULL)
        {
          while (diff_end2 != cp2)
          {
            if (!isspace (*cp2))
            {
              ValNodeAddInt (&diff_list, 2, offset2);
            }
            offset2++;
            cp2++;
          }
        }
        else if (*(cp1 + 1) == *cp2)
        {
          if (!isspace (*cp1))
          {
            ValNodeAddInt (&diff_list, 1, offset1);
          }
          offset1++;
          cp1++;
        }
        else if (*(cp2 + 1) == *cp1)
        {
          if (!isspace (*cp2))
          {
            ValNodeAddInt (&diff_list, 2, offset2);
          }
          offset2++;
          cp2++;
        }
        else if ((Int4)StringLen (cp1) > len_fake_modifier_name
                 && StringNCmp (cp1, fake_modifier_name, len_fake_modifier_name) == 0
                 && *(cp1 + len_fake_modifier_name) == *cp2)
        {
          /* show all of fake modifier name in red */
          for (j = 0; j < len_fake_modifier_name; j++)
          {
            ValNodeAddInt (&diff_list, 1, offset1);
            cp1++;
            offset1++;
          }
        }
        else if ((Int4) StringLen (cp2) > len_fake_modifier_name
                 && StringNCmp (cp2, fake_modifier_name, len_fake_modifier_name) == 0
                 && *(cp2 + len_fake_modifier_name) == *cp1)
        {
          /* show all of fake modifier name in red */
          for (j = 0; j < len_fake_modifier_name; j++)
          {
            ValNodeAddInt (&diff_list, 2, offset2);
            cp2++;
            offset2++;
          }
        }        
        else
        {
          diff_end1 = StringChr (cp1, *cp2);
          diff_end2 = StringChr (cp2, *cp1);
          if (diff_end1 == NULL || diff_end2 == NULL)
          {
            ValNodeAddInt (&diff_list, 1, offset1);
            ValNodeAddInt (&diff_list, 2, offset2);
          }
          else if (diff_end1 - cp1 < diff_end2 - cp2)
          {
            while (diff_end1 != cp1)
            {
              if (!isspace (*cp1))
              {
                ValNodeAddInt (&diff_list, 1, offset1);
              }
              offset1++;
              cp1++;
            }
          }
          else
          {
            while (diff_end2 != cp2)
            {
              if (!isspace (*cp2))
              {
                ValNodeAddInt (&diff_list, 2, offset2);
              }
              offset2++;
              cp2++;
            }
          }
        }
      }
    }
    
    if (cp1 != NULL)
    {
      cp1++;
      offset1++;
    }
    
    if (cp2 != NULL)
    {
      cp2++;
      offset2++;
    }
  }
  return diff_list;
}

static Int4 
CountTitleCorrectionRows 
(IDAndTitleEditPtr original, IDAndTitleEditPtr suggested)
{
  Int4 num_rows = 0, seq_num;
  
  if (original == NULL || suggested == NULL 
      || original->num_sequences != suggested->num_sequences)
  {
    return 0;
  }
  
  for (seq_num = 0; seq_num < original->num_sequences; seq_num++)
  {
    if (StringCmp (original->title_list [seq_num], suggested->title_list [seq_num]) != 0)
    {
      num_rows += 2;
    }
  }
  return num_rows;
}

static Int4 
DrawDiffPair 
(Int4              x,
 Int4              y,
 Int4              last_y,
 DiffDlgPtr        dlg,
 IDAndTitleEditPtr original,
 IDAndTitleEditPtr suggested,
 Int4              row_length,
 Int4Ptr           start_row,
 Int4              start_col,
 BoolPtr           shade)
{
  ValNodePtr diff_list;
  Int4       row, seq_num, visible_row;
  
  if (dlg == NULL || start_row == NULL || shade == NULL
      || y > last_y
      || original == NULL || suggested == NULL 
      || original->num_sequences != suggested->num_sequences)
  {
    return y;
  }
  
  SelectFont (dlg->display_font);  
  
  /* draw difference rows */
  diff_list = NULL;
  visible_row = 0;
  for (row = 0; 
       row < original->num_sequences * 2 && y <= last_y;
       row++)
  {
    seq_num = row / 2;
    
    if (row % 2 == 0)
    {
      /* draw original */
      if (StringCmp (original->title_list [seq_num], suggested->title_list [seq_num]) != 0)
      {
        diff_list = ValNodeFree (diff_list);
        diff_list = GetTextDifferences (original->title_list [seq_num],
                                    suggested->title_list [seq_num]);
        if (visible_row == *start_row)
        {
          DrawDiffDlgRow (x, y, dlg->char_width, dlg->descent, dlg->max_id_length, row_length,
                          original->id_list [seq_num],
                          original->title_list [seq_num],
                          start_col, diff_list, 1, *shade);
          y += stdLineHeight;
          (*start_row) ++;
        }
        visible_row++;
      }
    }
    else
    {
      /* draw suggested */
      if (StringCmp (original->title_list [seq_num], suggested->title_list [seq_num]) != 0)
      {
        if (diff_list == NULL)
        {
          /* only calculate if it was NULL, otherwise use same as previous diff_list */
          diff_list = GetTextDifferences (original->title_list [seq_num],
                                          suggested->title_list [seq_num]);
        }
        if (visible_row == *start_row)
        {
          DrawDiffDlgRow (x, y, dlg->char_width, dlg->descent, dlg->max_id_length, row_length,
                          suggested->id_list [seq_num],
                          suggested->title_list [seq_num],
                          start_col, diff_list, 2, *shade);
          y += stdLineHeight;
          (*start_row) ++;
        }
        visible_row++;
        /* toggle the shading */
        *shade = !(*shade);
      }
    }
  }
  diff_list = ValNodeFree (diff_list);
  return y;
}

static void OnDrawDiffDlg (PaneL p)
{
  DiffDlgPtr dlg;
  BaR          sb_vert, sb_horiz;
  Int4         start_row, start_col;
  RecT         r;
  Int4         x, y, row_length, last_y;
  Int4         num_new_rows, num_existing_rows, num_rows, visible_rows;
  Int4         new_vmax, new_hmax, old_vmax, old_hmax;
  Boolean      shade = TRUE;

  dlg = (DiffDlgPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  
  num_new_rows = CountTitleCorrectionRows (dlg->new_original, dlg->new_suggested);
  num_existing_rows = CountTitleCorrectionRows (dlg->existing_original, dlg->existing_suggested);
  num_rows = num_new_rows + num_existing_rows;
  
  if (!EditNeedsBracketingFixes (dlg->new_original) && ! EditNeedsBracketingFixes (dlg->existing_original))
  {
    return;
  }
  
  SelectFont (dlg->display_font);
  
  sb_vert  = GetSlateVScrollBar ((SlatE) p);
  Enable (sb_vert);
  sb_horiz = GetSlateHScrollBar ((SlatE) p);
  Enable (sb_horiz);
  
  start_row = GetBarValue (sb_vert);
  start_col = GetBarValue (sb_horiz);
  
  ObjectRect (p, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = r.left + 1;
  y = r.top + stdLineHeight;
  SelectFont (programFont); 
  
  row_length = r.right - r.left - 2;
  
  visible_rows = (r.bottom - r.top - 2 * dlg->table_inset) / stdLineHeight - dlg->num_header_rows;
  new_vmax = num_rows - visible_rows;
  new_hmax = dlg->max_title_length - 1;
  if (new_vmax < 0)
  {
    new_vmax = 0;
  }
  if (new_hmax < 0)
  {
    new_hmax = 0;
  }
  old_vmax = GetBarMax (sb_vert);
  old_hmax = GetBarMax (sb_horiz);
  
  if (old_vmax != new_vmax)
  {
    CorrectBarMax (sb_vert, new_vmax);
    if (start_row > new_vmax)
    {
      start_row = new_vmax;
    }
    CorrectBarValue (sb_vert, start_row);
    CorrectBarPage (sb_vert, 1, 1);
  }
  
  if (old_hmax != new_hmax)
  {
    CorrectBarMax (sb_horiz, new_hmax);
    if (start_col > new_hmax)
    {
      start_col = new_hmax;
    }
    CorrectBarValue (sb_horiz, start_col);
    CorrectBarPage (sb_horiz, 1, 1);
  }
  
  last_y = r.bottom - 2 * dlg->table_inset;

  /* draw explanatory text */
  DrawDiffDlgExplanation (x, y, dlg->descent, row_length);
  y+= 2 * stdLineHeight;
  
  /* draw title row */
  DrawDiffDlgTitle (x, y, dlg->char_width, dlg->descent, dlg->max_id_length, row_length);
  y+= stdLineHeight;
  
  y = DrawDiffPair (x, y, last_y, dlg, dlg->new_original, dlg->new_suggested,
                    row_length, &start_row, start_col, &shade);
  
  start_row -= num_new_rows;

  y = DrawDiffPair (x, y, last_y, dlg, dlg->existing_original, dlg->existing_suggested,
                    row_length, &start_row, start_col, &shade);
  
}

static void OnVScrollDiffDlg (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  RecT r;
  
  ObjectRect (s, &r);
  InvalRect (&r);
}

static void OnHScrollDiffDlg (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  RecT r;
  
  ObjectRect (s, &r);
  InvalRect (&r);
}

static PoinT GetDiffDlgCoord (DiffDlgPtr dlg, PoinT pt)
{
  BaR sb_horiz;
  BaR sb_vert;
  Int4 start_row, start_col;
  RecT r;
  PoinT cell_coord;
  Int4  x, y;
  
  cell_coord.x = -1;
  cell_coord.y = -1;
  
  if (dlg == NULL)
  {
    return cell_coord;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) dlg->pnl);
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->pnl);
  
  start_row = GetBarValue (sb_vert);
  start_col = GetBarValue (sb_horiz);
  
  ObjectRect (dlg->pnl, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = pt.x - r.left;
  y = pt.y - r.top;
  
  cell_coord.y = y / stdLineHeight;
  
  if (cell_coord.y >= dlg->num_header_rows)
  {
    cell_coord.y += GetBarValue (sb_vert);
  }
  
  cell_coord.x = x / dlg->char_width;
  if (cell_coord.x >= dlg->max_id_length + 10)
  {
    cell_coord.x += GetBarValue (sb_horiz);
  }

  return cell_coord;
}

static void 
ScrollForDiffInRow 
(Int4              row, 
 IDAndTitleEditPtr new_original, 
 IDAndTitleEditPtr new_suggested,
 IDAndTitleEditPtr existing_original,
 IDAndTitleEditPtr existing_suggested,
 BaR               sb_horiz)
{
  Int4 seq_num, displayed_row;
  ValNodePtr diff_list = NULL, vnp;
  Boolean    found_row = FALSE;
  Int4       scroll_val = 0;
  
  if (sb_horiz == NULL
      || (new_original == NULL && existing_original == NULL)
      || (new_original == NULL && new_suggested != NULL) 
      || (new_original != NULL && new_suggested == NULL) 
      || (new_original != NULL && new_original->num_sequences != new_suggested->num_sequences)
      || (existing_original == NULL && existing_suggested != NULL) 
      || (existing_original != NULL && existing_suggested == NULL) 
      || (existing_original != NULL && existing_original->num_sequences != existing_suggested->num_sequences))
  {
    return;
  }
  
  displayed_row = 0;
  
  if (new_original != NULL)
  {
    for (seq_num = 0;
         seq_num < new_original->num_sequences && ! found_row; 
         seq_num++)
    {
      if (StringCmp (new_original->title_list [seq_num],
                     new_suggested->title_list [seq_num]) != 0)
      {
        if (displayed_row == row)
        {
          found_row = TRUE;
          diff_list = GetTextDifferences (new_original->title_list [seq_num],
                                          new_suggested->title_list [seq_num]);
        }
        else
        {
          displayed_row += 2;
        }
      }
    }
  }

  if (existing_original != NULL)
  {
    for (seq_num = 0;
         seq_num < existing_original->num_sequences && ! found_row; 
         seq_num++)
    {
      if (StringCmp (existing_original->title_list [seq_num],
                     existing_suggested->title_list [seq_num]) != 0)
      {
        if (row == displayed_row || row == displayed_row + 1 )
        {
          found_row = TRUE;
          diff_list = GetTextDifferences (existing_original->title_list [seq_num],
                                          existing_suggested->title_list [seq_num]);
        }
        else
        {
          displayed_row += 2;
        }
      }
    }
  }
  
  if (diff_list == NULL)
  {
    scroll_val = 0;
  }
  else
  {
    scroll_val = GetBarValue (sb_horiz);
    vnp = diff_list;
    while (vnp != NULL && vnp->data.intvalue <= scroll_val)
    {
      vnp = vnp->next;
    }
    if (vnp != NULL)
    {
      scroll_val = vnp->data.intvalue;
    }
    else
    {
      scroll_val = diff_list->data.intvalue;
    }
  }
  SetBarValue (sb_horiz, scroll_val);
  diff_list = ValNodeFree (diff_list);
}

static void OnClickDiffDlg (PaneL p, PoinT pt)
{
  DiffDlgPtr dlg;
  Boolean    dbl_click;
  PoinT      cell_coord;

  dlg = (DiffDlgPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  dbl_click = dblClick;
  if (dbl_click)
  {
    cell_coord = GetDiffDlgCoord (dlg, pt);
    if (cell_coord.y >= dlg->num_header_rows)
    {
      cell_coord.y -= dlg->num_header_rows;
      ScrollForDiffInRow (cell_coord.y, dlg->new_original, dlg->new_suggested,
                          dlg->existing_original, dlg->existing_suggested,
                          GetSlateHScrollBar ((SlatE) dlg->pnl)); 
    }
  }
}

typedef struct diffset
{
  IDAndTitleEditPtr new_original;
  IDAndTitleEditPtr new_suggested;
  IDAndTitleEditPtr existing_original;
  IDAndTitleEditPtr existing_suggested;
} DiffSetData, PNTR DiffSetPtr;

static void SetToDiffDlg (DialoG d, Pointer userdata)
{
  DiffDlgPtr  dlg;
  DiffSetPtr  dsp;
  Int4        seq_num;
  RecT        r;
  
  dlg = (DiffDlgPtr) GetObjectExtra (d);
  dsp = (DiffSetPtr) userdata;
  if (dlg == NULL || dsp == NULL 
      || (dsp->new_original == NULL && dsp->existing_original == NULL)
      || (dsp->new_original == NULL && dsp->new_suggested != NULL)
      || (dsp->new_original != NULL && dsp->new_suggested == NULL)
      || (dsp->new_original != NULL && dsp->new_original->num_sequences != dsp->new_suggested->num_sequences)
      || (dsp->existing_original == NULL && dsp->existing_suggested != NULL)
      || (dsp->existing_original != NULL && dsp->existing_suggested == NULL)
      || (dsp->existing_original != NULL && dsp->existing_original->num_sequences != dsp->existing_suggested->num_sequences))
  {
    return;
  }
  
  dlg->new_original = IDAndTitleEditCopy (dsp->new_original);
  dlg->new_suggested = IDAndTitleEditCopy (dsp->new_suggested);
  dlg->existing_original = IDAndTitleEditCopy (dsp->existing_original);
  dlg->existing_suggested = IDAndTitleEditCopy (dsp->existing_suggested);
  
  dlg->max_id_length = 0;
  dlg->max_title_length = 0;
  
  /* get max lengths from new set */
  if (dsp->new_original != NULL)
  {
    for (seq_num = 0; seq_num < dsp->new_original->num_sequences; seq_num++)
    {
      /* we want the maximum length only for those rows we'll actually display */
      if (StringCmp (dsp->new_original->title_list [seq_num],
                     dsp->new_suggested->title_list [seq_num]) == 0)
      {
        continue;
      }
      /* max ID length */
      dlg->max_id_length = MAX (dlg->max_id_length, (Int4) StringLen (dsp->new_original->id_list [seq_num]));
      dlg->max_id_length = MAX (dlg->max_id_length, (Int4) StringLen (dsp->new_suggested->id_list [seq_num]));
    
      /* max title length */
      dlg->max_title_length = MAX (dlg->max_title_length, (Int4) StringLen (dsp->new_original->title_list [seq_num]));
      dlg->max_title_length = MAX (dlg->max_title_length, (Int4) StringLen (dsp->new_suggested->title_list [seq_num]));
    }
  }
  /* get max lengths from existing set */
  if (dsp->existing_original != NULL)
  {
    for (seq_num = 0; seq_num < dsp->existing_original->num_sequences; seq_num++)
    {
      /* we want the maximum length only for those rows we'll actually display */
      if (StringCmp (dsp->existing_original->title_list [seq_num],
                     dsp->existing_suggested->title_list [seq_num]) == 0)
      {
        continue;
      }
      /* max ID length */
      dlg->max_id_length = MAX (dlg->max_id_length, (Int4) StringLen (dsp->existing_original->id_list [seq_num]));
      dlg->max_id_length = MAX (dlg->max_id_length, (Int4) StringLen (dsp->existing_suggested->id_list [seq_num]));
    
      /* max title length */
      dlg->max_title_length = MAX (dlg->max_title_length, (Int4) StringLen (dsp->existing_original->title_list [seq_num]));
      dlg->max_title_length = MAX (dlg->max_title_length, (Int4) StringLen (dsp->existing_suggested->title_list [seq_num]));
    }
  }
  ObjectRect (dlg->pnl, &r);
  InvalRect (&r);   
}

static void CleanupDifferenceDialog (GraphiC g, Pointer data)
{
  DiffDlgPtr dlg;
  
  dlg = (DiffDlgPtr) data;
  if (dlg != NULL)
  {
    dlg->new_original = IDAndTitleEditFree (dlg->new_original);
    dlg->new_suggested = IDAndTitleEditFree (dlg->new_suggested);
    dlg->existing_original = IDAndTitleEditFree (dlg->existing_original);
    dlg->existing_suggested = IDAndTitleEditFree (dlg->existing_suggested);
    dlg = MemFree (dlg);
  }
}

static DialoG 
ShowDifferenceDialog 
(GrouP parent, 
 Int4  width,
 Int4  height)
{
  DiffDlgPtr dlg;
  GrouP        p;
  
  dlg = (DiffDlgPtr) MemNew (sizeof (DiffDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (parent, 1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupDifferenceDialog);
  
  dlg->dialog = (DialoG) p;
  dlg->todialog = SetToDiffDlg;
  dlg->fromdialog = NULL;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;

  dlg->new_original = NULL;
  dlg->new_suggested = NULL;  
  dlg->existing_original = NULL;
  dlg->existing_suggested = NULL;  

#ifdef WIN_MAC
  dlg->display_font = ParseFont ("Monaco, 9");
#endif
#ifdef WIN_MSWIN
  dlg->display_font = ParseFont ("Courier, 9");
#endif
#ifdef WIN_MOTIF
  dlg->display_font = ParseFont ("fixed, 12");
#endif
  SelectFont (dlg->display_font);
  dlg->char_width  = CharWidth ('0');
  dlg->descent = Descent ();
  dlg->table_inset = 4;
  
  dlg->max_title_length = 0;
  
  dlg->num_header_rows = 3;
  
  dlg->pnl = AutonomousPanel4 (p, width, height, OnDrawDiffDlg,
                               OnVScrollDiffDlg, OnHScrollDiffDlg,
                               sizeof (DiffDlgData), NULL, NULL); 
  SetObjectExtra (dlg->pnl, dlg, NULL);
  SetPanelClick(dlg->pnl, OnClickDiffDlg, NULL, NULL, NULL);
  
  return (DialoG) p;    
}

static Boolean 
ShowBracketingCorrections 
(IDAndTitleEditPtr iatep_new, 
 IDAndTitleEditPtr iatep_existing,
 DialoG dlg)
{
  IDAndTitleEditPtr  suggested_new, suggested_existing;
  Boolean            any_to_show = FALSE;
  DiffSetData        dsd;
  Int4               num_to_show;
  
  if (dlg == NULL)
  {
    return FALSE;
  }

  suggested_new = SuggestCorrectionForTitleBracketing (iatep_new);
  suggested_existing = SuggestCorrectionForTitleBracketing (iatep_existing);
  
  num_to_show = CountTitleCorrectionRows (iatep_new, suggested_new)
                + CountTitleCorrectionRows (iatep_existing, suggested_existing);
                
  if (num_to_show > 0)
  {
    dsd.new_original = iatep_new;
    dsd.new_suggested = suggested_new;
    dsd.existing_original = iatep_existing;
    dsd.existing_suggested = suggested_existing;
    PointerToDialog (dlg, &dsd);
    any_to_show = TRUE;
  }
  suggested_new = IDAndTitleEditFree (suggested_new);
  suggested_new = IDAndTitleEditFree (suggested_new);
  
  return any_to_show;
}

typedef Int4 (*DrawExplanationFunc) PROTO ((Int4, Int4, Int4, Int4, Int4));

typedef ValNodePtr (*ColorizeStringFunc) PROTO ((CharPtr, Pointer, Boolean));

typedef void (*UpdateColorizedPanelParentProc) PROTO ((Pointer));

typedef void (*ScrollParentProc) PROTO ((Int4, Pointer));

typedef struct colorizeddeflinedlg
{
  PaneL                          pnl;
  IDAndTitleEditPtr              iatep_new;
  IDAndTitleEditPtr              iatep_current;
  Boolean                        is_nuc;
  FonT                           display_font;
  Int4                           char_width;
  Int4                           descent;
  Int4                           max_title_length;
  Int4                           max_id_length;
  Int4                           table_inset;
  Int4                           num_header_rows;
  DrawExplanationFunc            draw_explanation; 
  Boolean                        edit_values; 
  ColorizeStringFunc             colorize_title;
  Pointer                        colorize_data;
  UpdateColorizedPanelParentProc update_parent;
  Pointer                        update_parent_data;
  ScrollParentProc               scroll_parent;
  Pointer                        scroll_parent_data;
} ColorizedDeflineDlgData, PNTR ColorizedDeflineDlgPtr;

static Int4 CountRowsWithColor 
(IDAndTitleEditPtr  iatep, 
 Int4Ptr            max_id_length,
 Int4Ptr            max_title_length,
 ColorizeStringFunc colorize_title,
 Pointer            colorize_data,
 Boolean            is_nuc)
{
  ValNodePtr diff_list;
  Int4       seq_num, num_rows_with_color = 0;
  
  if (iatep == NULL || colorize_title == NULL)
  {
    return 0;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    diff_list = colorize_title (iatep->title_list [seq_num], colorize_data, is_nuc);
    if (diff_list != NULL)
    {
      num_rows_with_color ++;
      diff_list = ValNodeFree (diff_list);
      if (max_id_length != NULL)
      {
        *max_id_length = MAX (*max_id_length, (Int4) StringLen (iatep->id_list [seq_num]));
      }
      if (max_title_length != NULL)
      {
        *max_title_length = MAX (*max_title_length, (Int4) StringLen (iatep->title_list [seq_num]));
      }
    }
  }
  return num_rows_with_color;
}

static void DrawDeflineDlgTitle (Int4 x, Int4 y, Int4 char_width, Int4 descent, Int4 max_id_length, Int4 win_width)
{
  RecT rct;
  
  /* draw title row */
  DkGray ();
  InvertColors ();
  White ();
  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  PaintStringEx ("Sequence ID", x, y);
  x += (max_id_length + 2) * char_width;
  PaintStringEx ("Title", x, y);
  InvertColors ();
  Black ();
}

static void 
DrawDeflineDlgRow 
(Int4 x, 
 Int4 y, 
 Int4 char_width, 
 Int4 descent, 
 Int4 max_id_length, 
 Int4 win_width,
 CharPtr id_str,
 CharPtr title_str,
 Int4       offset,
 ValNodePtr diff_list,
 Int4       choice_num)
{
  RecT rct;
  PoinT      pt1, pt2;

  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  if (id_str != NULL)
  {
    PaintStringEx (id_str, x, y);
  }
  x += (max_id_length + 1) * char_width;
  pt1.x = x + 2;
  pt2.x = x + 2;
  pt1.y = y + descent;
  pt2.y = y - stdLineHeight + descent;
  Black ();
  DrawLine (pt1, pt2);

  x += char_width;
  
  PaintColorizedString (x, y, char_width, title_str, offset, FALSE, diff_list, 
                        choice_num);
}

static void 
DrawColorizedDeflinesInSet 
(Int4                x,
 Int4Ptr             y,
 Int4                last_y,
 Int4                char_width,
 Int4                descent,
 Int4                max_id_length,
 Int4                row_length,
 Int4Ptr             start_row,
 Int4                start_col,
 IDAndTitleEditPtr   iatep,
 ColorizeStringFunc  colorize_title,
 Pointer             colorize_data,
 Boolean             is_nuc)
{
  Int4       row, visible_row;
  ValNodePtr diff_list;
  
  if (iatep == NULL || y == NULL || start_row == NULL || colorize_title == NULL)
  {
    return;
  }
  
  visible_row = 0;
  for (row = 0; 
       row < iatep->num_sequences && *y <= last_y;
       row++)
  {
    diff_list = NULL;
    diff_list = colorize_title (iatep->title_list [row], colorize_data, is_nuc);
    if (diff_list != NULL)
    {
      if (visible_row == *start_row)
      {
        DrawDeflineDlgRow (x, *y, char_width, descent, max_id_length, row_length,
                          iatep->id_list [row],
                          iatep->title_list [row],
                          start_col, diff_list, 1);
        (*y) += stdLineHeight;
        (*start_row) ++;
      }
      visible_row++;
      diff_list = ValNodeFree (diff_list);
    }
  }
  
} 

static void OnDrawColorizedDeflineDlg (PaneL p)
{
  ColorizedDeflineDlgPtr dlg;
  BaR                    sb_vert, sb_horiz;
  Int4                   start_row, start_col;
  RecT                   r;
  Int4                   x, y, row_length, last_y;
  Int4                   num_new_rows, num_existing_rows, num_rows;
  Int4                   visible_rows;
  Int4                   new_vmax, new_hmax, old_vmax, old_hmax;

  dlg = (ColorizedDeflineDlgPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  
  num_rows = 0;
  dlg->max_id_length = 10;
  dlg->max_title_length = 5;
  num_new_rows = CountRowsWithColor (dlg->iatep_new, 
                                     &(dlg->max_id_length), &(dlg->max_title_length),
                                     dlg->colorize_title, dlg->colorize_data, dlg->is_nuc);
  num_existing_rows = CountRowsWithColor (dlg->iatep_current, 
                                          &(dlg->max_id_length), &(dlg->max_title_length),
                                          dlg->colorize_title, dlg->colorize_data, dlg->is_nuc);
  num_rows = num_new_rows + num_existing_rows;
  
  SelectFont (dlg->display_font);
  
  sb_vert  = GetSlateVScrollBar ((SlatE) p);
  Enable (sb_vert);
  sb_horiz = GetSlateHScrollBar ((SlatE) p);
  Enable (sb_horiz);
  
  start_row = GetBarValue (sb_vert);
  start_col = GetBarValue (sb_horiz);
  
  ObjectRect (p, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = r.left + 1;
  y = r.top + stdLineHeight;
  SelectFont (programFont); 
  
  row_length = r.right - r.left - 2;
  
  dlg->num_header_rows = 0;
  /* draw explanatory text */
  if (dlg->draw_explanation != NULL)
  {
    dlg->num_header_rows = (dlg->draw_explanation) (x, y, 
                                                    dlg->char_width, 
                                                    dlg->descent, 
                                                    row_length);
    y += stdLineHeight * dlg->num_header_rows;
  }
  
  /* draw title row */
  DrawDeflineDlgTitle (x, y, dlg->char_width, dlg->descent, dlg->max_id_length, row_length);
  y+= stdLineHeight;
  dlg->num_header_rows ++;  
  
  visible_rows = (r.bottom - r.top - 2 * dlg->table_inset) / stdLineHeight - dlg->num_header_rows;
  new_vmax = num_rows - visible_rows;
  new_hmax = dlg->max_title_length - 1;
  if (new_vmax < 0)
  {
    new_vmax = 0;
  }
  if (new_hmax < 0)
  {
    new_hmax = 0;
  }
  old_vmax = GetBarMax (sb_vert);
  old_hmax = GetBarMax (sb_horiz);
  
  if (old_vmax != new_vmax)
  {
    CorrectBarMax (sb_vert, new_vmax);
    if (start_row > new_vmax)
    {
      start_row = new_vmax;
    }
    CorrectBarValue (sb_vert, start_row);
    CorrectBarPage (sb_vert, 1, 1);
  }
  
  if (old_hmax != new_hmax)
  {
    CorrectBarMax (sb_horiz, new_hmax);
    if (start_col > new_hmax)
    {
      start_col = new_hmax;
    }
    CorrectBarValue (sb_horiz, start_col);
    CorrectBarPage (sb_horiz, 1, 1);
  }
  
  last_y = r.bottom - 2 * dlg->table_inset;

  DrawColorizedDeflinesInSet (x, &y, last_y, dlg->char_width, dlg->descent,
                              dlg->max_id_length, row_length, &start_row, start_col,
                              dlg->iatep_new, 
                              dlg->colorize_title, dlg->colorize_data,
                              dlg->is_nuc);

  start_row -= num_new_rows;

  DrawColorizedDeflinesInSet (x, &y, last_y, dlg->char_width, dlg->descent,
                              dlg->max_id_length, row_length, &start_row, start_col,
                              dlg->iatep_current, 
                              dlg->colorize_title, dlg->colorize_data,
                              dlg->is_nuc);
}

static PoinT GetColorizedDeflineCoord (ColorizedDeflineDlgPtr dlg, PoinT pt)
{
  BaR sb_horiz;
  BaR sb_vert;
  Int4 start_row, start_col;
  RecT r;
  PoinT cell_coord;
  Int4  x, y, apparent_x, apparent_y, vis_row, new_rows = 0;
  ValNodePtr diff_list;
  
  cell_coord.x = -1;
  cell_coord.y = -1;
  
  if (dlg == NULL)
  {
    return cell_coord;
  }
  
  sb_vert  = GetSlateVScrollBar ((SlatE) dlg->pnl);
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->pnl);
  
  start_row = GetBarValue (sb_vert);
  start_col = GetBarValue (sb_horiz);
  
  ObjectRect (dlg->pnl, &r);
  InsetRect (&r, dlg->table_inset, dlg->table_inset);
  x = pt.x - r.left;
  y = pt.y - r.top;
  
  apparent_y = y / stdLineHeight;
  
  if (apparent_y < dlg->num_header_rows)
  {
    cell_coord.y = -1;
  }
  else
  {
    apparent_y = apparent_y - dlg->num_header_rows + start_row;
    cell_coord.y = 0;
    vis_row = -1;
    while (vis_row < apparent_y && dlg->iatep_new != NULL && cell_coord.y < dlg->iatep_new->num_sequences)
    {
      diff_list = NULL;
      diff_list = dlg->colorize_title (dlg->iatep_new->title_list [cell_coord.y], 
                                       dlg->colorize_data,
                                       dlg->is_nuc);
      if (diff_list != NULL)
      {
        vis_row++;
        new_rows ++;
        if (vis_row < apparent_y)
        {
          cell_coord.y ++;
        }
      }
      else
      {
        cell_coord.y ++;
      }
      diff_list = ValNodeFree (diff_list);
    }
    while (vis_row < apparent_y && dlg->iatep_current != NULL 
           && cell_coord.y - new_rows < dlg->iatep_current->num_sequences)
    {
      diff_list = NULL;
      diff_list = dlg->colorize_title (dlg->iatep_current->title_list [cell_coord.y - new_rows],
                                       dlg->colorize_data,
                                       dlg->is_nuc);
      if (diff_list != NULL)
      {
        vis_row++;
        if (vis_row < apparent_y)
        {
          cell_coord.y ++;
        }
      }
      else
      {
        cell_coord.y ++;
      }
      diff_list = ValNodeFree (diff_list);
    }
  }
  
  apparent_x = x / dlg->char_width;
  if (apparent_x <= dlg->max_id_length + 1)
  {
    cell_coord.x = -1;
  }
  else 
  {
    cell_coord.x = apparent_x - dlg->max_id_length - 2 + start_col;
  }

  return cell_coord;
}

static BadValuePtr GetModValuePairForCoord (ColorizedDeflineDlgPtr dlg, PoinT coord, BoolPtr is_value)
{
  Int4                   current_num;
  CharPtr                seq_id = NULL, title = NULL, mod_name = NULL, mod_value = NULL;
  CharPtr                cp, eq_loc, start_bracket, end_bracket;
  BadValuePtr            bvp = NULL;
  ModifierInfoPtr        mip;

  if (dlg == NULL || coord.x < 0 || coord.y < 0)
  {
    return NULL;
  }
  
  if (dlg->iatep_new != NULL && coord.y < dlg->iatep_new->num_sequences)
  {
    seq_id = dlg->iatep_new->id_list [coord.y];
    title = dlg->iatep_new->title_list [coord.y];
  }
  else if (dlg->iatep_current != NULL)
  {
    current_num = coord.y;
    if (dlg->iatep_new != NULL)
    {
      current_num -= dlg->iatep_new->num_sequences;
    }
    if (current_num < dlg->iatep_current->num_sequences)
    {
      seq_id = dlg->iatep_current->id_list [current_num];
      title = dlg->iatep_current->id_list [current_num];
    }
  }
  if (seq_id == NULL || title == NULL || (Int4) StringLen (title) < coord.x)
  {
    return NULL;
  }
  
  cp = title + coord.x;
  
  mip = ParseOneBracketedModifier (title, &start_bracket, &end_bracket);
  while (mip != NULL && start_bracket != NULL && end_bracket != NULL
         && cp > end_bracket)
  {
    mip = ModifierInfoFree (mip);
    mip = ParseOneBracketedModifier (end_bracket + 1, &start_bracket, &end_bracket);
  }
  mip = ModifierInfoFree (mip);
  
  mod_name = NULL;
  mod_value = NULL;
  
  if (start_bracket <= cp && end_bracket >= cp)
  {
    eq_loc = NextBracketToken (start_bracket + 1);
    if (eq_loc != NULL && *eq_loc == '=')
    {
      mod_name = (CharPtr) MemNew ((eq_loc - start_bracket) * sizeof (Char));
      if (mod_name != NULL)
      {
        StringNCpy (mod_name, start_bracket + 1, eq_loc - start_bracket - 1);
        mod_name [eq_loc - start_bracket - 1] = 0;
        TrimSpacesAroundString (mod_name);
      }
      
      mod_value = (CharPtr) MemNew (end_bracket - eq_loc);
      if (mod_value != NULL)
      {
        StringNCpy (mod_value, eq_loc + 1, end_bracket - eq_loc - 1);
        mod_value [end_bracket - eq_loc - 1] = 0;
        TrimSpacesAroundString (mod_value);
      }
    }
  }
  
  if (mod_name == NULL || mod_value == NULL)
  {
    mod_name = MemFree (mod_name);
    mod_value = MemFree (mod_value);
    return NULL;
  }
  
  bvp = BadValueNew (seq_id, mod_name, mod_value);
  mod_name = MemFree (mod_name);
  mod_value = MemFree (mod_value);
  
  if (is_value != NULL)
  {
    if (title + coord.x < eq_loc)
    {
      *is_value = FALSE;
    }
    else if (title + coord.x >= eq_loc)
    {
      *is_value = TRUE;
    }
  }
  return bvp;
}

static void ScrollToColor (ColorizedDeflineDlgPtr dlg, PoinT coord)
{
  BaR  sb_horiz;
  Int4 current_scroll_pos, current_row = 0, new_scroll_pos = 0;
  ValNodePtr diff_list, vnp;
  CharPtr    title = NULL;
  Boolean    found_in_new = FALSE;
  
  if (dlg == NULL || dlg->colorize_title == NULL || coord.y < 0)
  {
    return;
  }
  
  sb_horiz = GetSlateHScrollBar ((SlatE) dlg->pnl);
  current_scroll_pos = GetBarValue (sb_horiz);
  
  current_row = coord.y;
  
  if (dlg->iatep_new != NULL)
  {
    if (coord.y < dlg->iatep_new->num_sequences)
    {
      title = dlg->iatep_new->title_list [coord.y];
      found_in_new = TRUE;
    }
    else
    {
      current_row = coord.y - dlg->iatep_new->num_sequences;
    }
  }
  
  if (!found_in_new && dlg->iatep_current != NULL && current_row < dlg->iatep_current->num_sequences)
  {
    title = dlg->iatep_current->title_list [current_row];
  }
  
  if (title == NULL)
  {
    return;
  }
  
  diff_list = (dlg->colorize_title) (title, dlg->colorize_data, dlg->is_nuc);
  vnp = diff_list;
  while (vnp != NULL && vnp->data.intvalue <= current_scroll_pos)
  {
    vnp = vnp->next;
  }
  while (vnp != NULL && vnp->data.intvalue == current_scroll_pos + 1)
  {
    vnp = vnp->next;
    current_scroll_pos ++;
  }
  if (vnp != NULL)
  {
    new_scroll_pos = vnp->data.intvalue;
  }
  if (new_scroll_pos > GetBarMax (sb_horiz))
  {
    new_scroll_pos = GetBarMax (sb_horiz);
  }
  SetBarValue (sb_horiz, new_scroll_pos); 
   
}

static void OnClickColorizedDeflinePanel (PaneL p, PoinT pt)
{
  ColorizedDeflineDlgPtr dlg;
  PoinT                  cell_coord;
  BadValuePtr            bvp;
  Boolean                is_value = FALSE, rval = FALSE;

  dlg = (ColorizedDeflineDlgPtr) GetObjectExtra (p);
  if (dlg == NULL || ! dblClick)
  {
    return;
  }

  cell_coord = GetColorizedDeflineCoord (dlg, pt);
  if (cell_coord.y < 0)
  {
    return;
  }
  else if (cell_coord.x < 0)
  {
    if (dlg->scroll_parent != NULL)
    {
      (dlg->scroll_parent) (cell_coord.y, dlg->scroll_parent_data);
    }
    ScrollToColor (dlg, cell_coord);
  }
  else
  {
    bvp = GetModValuePairForCoord (dlg, cell_coord, &is_value);
    if (bvp != NULL)
    {
      if (is_value)
      {
        if (dlg->edit_values)
        {
          rval = FixOneModifierValue (dlg->iatep_new,
                                      dlg->iatep_current,
                                      bvp->seq_id,
                                      bvp->mod_name,
                                      bvp->value,
                                      GetModifierType (bvp->mod_name));
        }
      }
      else
      {
        rval = FixOneModifierName (dlg->iatep_new,
                                   dlg->iatep_current,
                                   bvp->seq_id,
                                   bvp->mod_name,
                                   dlg->is_nuc);
      }
    }
    bvp = BadValueFree (bvp);
    if (rval && dlg->update_parent != NULL)
    {
      (dlg->update_parent) (dlg->update_parent_data);
    }
  }
}

static void 
UpdateColorizedDeflinePanelData 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 PaneL             pnl)
{
  ColorizedDeflineDlgPtr dlg;
  RecT                   r;
  
  dlg = (ColorizedDeflineDlgPtr) GetObjectExtra (pnl);
  if (dlg == NULL)
  {
    return;
  }
  dlg->iatep_new = iatep_new;
  dlg->iatep_current = iatep_current;

  
  ObjectRect ((SlatE) pnl, &r);
  InvalRect (&r);
}

static PaneL 
ColorizedDeflinePanel 
(GrouP parent, 
 Int4  width,
 Int4  height,
 IDAndTitleEditPtr   iatep_new,
 IDAndTitleEditPtr   iatep_current,
 Boolean             is_nuc,
 DrawExplanationFunc draw_explanation,
 Boolean             edit_values,
 ColorizeStringFunc  colorize_title,
 Pointer             colorize_data,
 UpdateColorizedPanelParentProc update_parent,
 Pointer                        update_parent_data,
 ScrollParentProc               scroll_parent,
 Pointer                        scroll_parent_data)
{
  ColorizedDeflineDlgPtr dlg;
  
  dlg = (ColorizedDeflineDlgPtr) MemNew (sizeof (ColorizedDeflineDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
#ifdef WIN_MAC
  dlg->display_font = ParseFont ("Monaco, 9");
#endif
#ifdef WIN_MSWIN
  dlg->display_font = ParseFont ("Courier, 9");
#endif
#ifdef WIN_MOTIF
  dlg->display_font = ParseFont ("fixed, 12");
#endif
  SelectFont (dlg->display_font);
  dlg->char_width  = CharWidth ('0');
  dlg->descent = Descent ();
  dlg->table_inset = 4;
  
  dlg->max_title_length = 0;
  dlg->max_id_length = 0;
  
  dlg->draw_explanation = draw_explanation;
  dlg->edit_values = edit_values;
  dlg->colorize_title = colorize_title;
  dlg->colorize_data = colorize_data;
  dlg->iatep_new = iatep_new;
  dlg->iatep_current = iatep_current;
  dlg->is_nuc = is_nuc;
  dlg->update_parent = update_parent;
  dlg->update_parent_data = update_parent_data;
  dlg->scroll_parent = scroll_parent;
  dlg->scroll_parent_data = scroll_parent_data;
  
  dlg->pnl = AutonomousPanel4 (parent, width, height, OnDrawColorizedDeflineDlg,
                               OnVScrollDiffDlg, OnHScrollDiffDlg,
                               sizeof (ColorizedDeflineDlgData), NULL, NULL); 
  SetObjectExtra (dlg->pnl, dlg, NULL);

  SetPanelClick(dlg->pnl, OnClickColorizedDeflinePanel, NULL, NULL, NULL); 
  
  return dlg->pnl;    
}

static void UpdateSeqIdEditForColorizedPanel (Pointer userdata)
{
  Boolean      show_all;
  SeqIdEditPtr siep;
  
  siep = (SeqIdEditPtr) userdata;
  if (siep == NULL)
  {
    return;
  }
  
  if (GetValue (siep->show_all_grp) == 2)
  {
    show_all = TRUE;
  }
  else
  {
    show_all = FALSE;
  }
  
  UpdateIdAndTitleEditDialog (siep->new_dlg, 
                              siep->iatep_new, 
                              siep->iatep_current, 
                              siep->seqid_edit_phase,
                              show_all,
                              siep->is_nuc);
  UpdateIdAndTitleEditDialog (siep->current_dlg, 
                              siep->iatep_current, 
                              siep->iatep_new, 
                              siep->seqid_edit_phase,
                              show_all,
                              siep->is_nuc);
  ShowErrorInstructions (siep);
}

static void ScrollSeqIdEditForColorizedPanel (Int4 seq_num, Pointer userdata)
{
  SeqIdEditPtr siep;
  Int4         current_num;
  
  if (seq_num < 0)
  {
    return;
  }
  
  siep = (SeqIdEditPtr) userdata;
  if (siep == NULL)
  {
    return;
  }
  
  if (siep->iatep_new != NULL && seq_num < siep->iatep_new->num_sequences)
  {
    ScrollTagListToSeqId (siep->new_dlg, siep->iatep_new->id_list [seq_num]);
  }
  else if (siep->iatep_current != NULL)
  {
    current_num = seq_num;
    if (siep->iatep_new != NULL)
    {
      current_num -= siep->iatep_new->num_sequences;
    }
    if (current_num < siep->iatep_current->num_sequences)
    {
      ScrollTagListToSeqId (siep->current_dlg,
                            siep->iatep_current->id_list [current_num]);
    }
  }
}

static Int4 
DrawExplanation 
(Int4 x,
 Int4 y, 
 Int4 char_width, 
 Int4 descent, 
 Int4 win_width,
 CharPtr line1, 
 CharPtr exp_part1, 
 CharPtr exp_red, 
 CharPtr exp_part2)
{
  RecT rct;
  Int4 num_lines = 0, tmp_x;
  CharPtr  dbl_click = "Double-click on a sequence ID to scroll to the next invalid ";
  
  
  /* draw first explanation row */
  LoadRect (&rct, x, y + descent,
            x + win_width, 
            y - stdLineHeight + descent);
  EraseRect (&rct);

  if (!StringHasNoText (line1))
  {
    LoadRect (&rct, x, y + descent,
              x + win_width, 
              y - stdLineHeight + descent);
    EraseRect (&rct);
    PaintStringEx (line1, x, y);
    y += stdLineHeight;
    num_lines ++;
  }
  
  if (!StringHasNoText (exp_part1) || !StringHasNoText (exp_red) || !StringHasNoText (exp_part2))
  {
    LoadRect (&rct, x, y + descent,
              x + win_width, 
              y - stdLineHeight + descent);
    EraseRect (&rct);

    tmp_x = x;
    PaintStringEx (exp_part1, tmp_x, y);
    tmp_x += StringLen (exp_part1) * char_width;
    
    Red ();
    PaintStringEx (exp_red, tmp_x, y);
    Black ();
    tmp_x += StringLen (exp_red) * char_width;
  
    PaintStringEx (exp_part2, tmp_x, y);
    
    y += stdLineHeight;
    num_lines ++;
  }
  
  if (!StringHasNoText (exp_red))
  {
    /* draw second row */
    LoadRect (&rct, x, y + descent,
              x + win_width, 
              y - stdLineHeight + descent);
    EraseRect (&rct);
    PaintStringEx (dbl_click, x, y);
    tmp_x = x + StringLen (dbl_click) * char_width;
    Red ();
    PaintStringEx (exp_red, tmp_x, y);
    Black ();
    tmp_x += StringLen (exp_red) * char_width;
    PaintStringEx (".", tmp_x, y);
    y += stdLineHeight;
    num_lines ++;
  }
  
  
  return num_lines;
}

static Int4 
DrawInvalidNameExplanation 
(Int4 x,
 Int4 y, 
 Int4 char_width, 
 Int4 descent, 
 Int4 win_width)
{
  /* draw explanation rows */
  return DrawExplanation (x, y, char_width, descent, win_width, 
                          "Some of your modifiers have invalid names.",
                          "Double-click on a ",
                          "modifier name",
                          " to change the name.");
}

static ValNodePtr ColorizeUnrecognizedNames (CharPtr title, Pointer userdata, Boolean is_nuc)
{
  Int4            offset;
  CharPtr         stop, cp;
  ModifierInfoPtr mip;
  ValNodePtr      diff_list = NULL;
    
  if (StringHasNoText (title))
  {
    return NULL;
  }

  cp = StringChr (title, '[');
  mip = ParseOneBracketedModifier (cp, NULL, &stop);
  while (mip != NULL && stop != NULL)
  {
  	if (IsUnrecognizedModifierName (mip, is_nuc))
  	{
  	  cp++;
      offset = cp - title;
      while (*cp != 0 && *cp != '=')
      {
        if (!isspace (*cp))
        {
          ValNodeAddInt (&diff_list, 1, offset);
        }
        cp++;
        offset++;
      }
  	}
  	mip = ModifierInfoFree (mip);
  	cp = StringChr (stop + 1, '[');
  	mip = ParseOneBracketedModifier (cp, NULL, &stop);
  }
  mip = ModifierInfoFree (mip);
  return diff_list;
}

static PaneL 
UnrecognizedModifiersPanel 
(GrouP                          parent,
 Int4                           width,
 Int4                           height,
 IDAndTitleEditPtr              iatep_new,
 IDAndTitleEditPtr              iatep_current,
 Boolean                        is_nuc,
 UpdateColorizedPanelParentProc update_parent,
 Pointer                        update_parent_data,
 ScrollParentProc               scroll_parent,
 Pointer                        scroll_parent_data)
{
  return ColorizedDeflinePanel (parent, width, height, 
                                iatep_new,
                                iatep_current,
                                is_nuc,
                                DrawInvalidNameExplanation,
                                FALSE,
                                ColorizeUnrecognizedNames,
                                NULL,
                                update_parent,
                                update_parent_data,
                                scroll_parent,
                                scroll_parent_data);
}

static Int4 
DrawBadValuesExplanation 
(Int4 x, 
 Int4 y, 
 Int4 char_width, 
 Int4 descent, 
 Int4 win_width)
{
  /* draw explanation rows */
  return DrawExplanation (x, y, char_width, descent, win_width,
                          "Some of your modifiers have inappropriate values.",
                          "Double-click on a modifier name to replace the modifier name, double-click on a ",
                          "value",
                          " to replace the value.");
}

static ValNodePtr ColorizeInvalidValues (CharPtr title, Pointer userdata, Boolean is_nuc)
{
  Int4            offset;
  CharPtr         start, stop, cp;
  ModifierInfoPtr mip;
  ValNodePtr      diff_list = NULL;
    
  if (StringHasNoText (title))
  {
    return NULL;
  }

  cp = StringChr (title, '[');
  mip = ParseOneBracketedModifier (cp, &start, &stop);
  while (mip != NULL && stop != NULL)
  {
  	if (ModifierHasInvalidValue (mip))
  	{
  	  cp = NextBracketToken (start + 1);
      offset = cp - title;
      while (*cp != 0 && cp != stop)
      {
        if (!isspace (*cp))
        {
          ValNodeAddInt (&diff_list, 1, offset);
        }
        cp++;
        offset++;
      }
  	}
  	mip = ModifierInfoFree (mip);
  	cp = StringChr (stop + 1, '[');
  	mip = ParseOneBracketedModifier (cp, &start, &stop);
  }
  mip = ModifierInfoFree (mip);
  return diff_list;
}

static PaneL
InvalidValuesPanel
(GrouP                          parent,
 Int4                           width,
 Int4                           height,
 IDAndTitleEditPtr              iatep_new,
 IDAndTitleEditPtr              iatep_current,
 Boolean                        is_nuc,
 UpdateColorizedPanelParentProc update_parent,
 Pointer                        update_parent_data,
 ScrollParentProc               scroll_parent,
 Pointer                        scroll_parent_data)
{
  return ColorizedDeflinePanel (parent, width, height,
                                iatep_new,
                                iatep_current, 
                                is_nuc, 
                                DrawBadValuesExplanation,
                                TRUE,
                                ColorizeInvalidValues,
                                NULL,
                                update_parent,
                                update_parent_data,
                                scroll_parent,
                                scroll_parent_data);                                
}
  
static Boolean 
ShowUnrecognizedModifiers 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 Boolean           is_nuc)
{
  ValNodePtr   unrecognized_list;
  
  unrecognized_list = ListUnrecognizedModifiers (iatep_new, iatep_current, is_nuc);
  
  if (unrecognized_list == NULL)
  {
    return FALSE;
  }
  else
  {
    unrecognized_list = ValNodeFreeData (unrecognized_list);
    return TRUE;
  }
}

static Boolean ShowBadValues (IDAndTitleEditPtr iatep_new, IDAndTitleEditPtr iatep_current)
{
  ValNodePtr badlist = NULL;
  Boolean    rval = FALSE;
  
  FindBadValuesInIDsAndTitles (iatep_new, &badlist);
  FindBadValuesInIDsAndTitles (iatep_current, &badlist);

  if (badlist != NULL)
  {
    rval = TRUE;
  }

  badlist = BadValueListFree (badlist);
  return rval;
}

static Boolean 
IDAndTitleEditIDsNeedFix 
(IDAndTitleEditPtr iatep_new,
 IDAndTitleEditPtr iatep_current,
 Boolean           is_nuc)
{
  if (HasMissingIDs (iatep_new) || HasMissingIDs (iatep_current)
      || (EditHasDuplicateIDs (iatep_new, iatep_current) && is_nuc)
      || AnyBracketsInIDs (iatep_new) || AnyBracketsInIDs (iatep_current)
      || EditHasSpaceInIDs (iatep_new) || EditHasSpaceInIDs (iatep_current))
  {
    return TRUE;
  }
  else 
  {
    return FALSE;
  }
}

static Boolean IDAndTitleEditTitlesNeedFix (IDAndTitleEditPtr iatep, Boolean is_nuc)
{
  Boolean    need_fix = FALSE;
  ValNodePtr unrec_list = NULL;
 
  if (EditNeedsBracketingFixes (iatep)
      || (unrec_list = ListUnrecognizedModifiers (iatep, NULL, is_nuc)) != NULL
      || ShowBadValues (iatep, NULL))
  {
    need_fix = TRUE;
  }
  unrec_list = ValNodeFreeData (unrec_list);
  return need_fix;
}

static Boolean EditIDsNeedFix (SeqIdEditPtr siep)
{
  Boolean      need_fix = FALSE;
  
  if (siep == NULL || siep->iatep_new == NULL || siep->iatep_new->num_sequences < 1)
  {
    return FALSE;
  }
    
  if (IDAndTitleEditIDsNeedFix (siep->iatep_new, siep->iatep_current, siep->is_nuc)
      || IDAndTitleEditTitlesNeedFix (siep->iatep_new, siep->is_nuc)
      || IDAndTitleEditTitlesNeedFix (siep->iatep_current, siep->is_nuc))
  {
    need_fix = TRUE;
  }
  return need_fix;
}

static Boolean IsSeqNumHidden (DialoG d, Int4 seq_num)
{
  TagListPtr tlp;
  Char       seq_str [15];
  ValNodePtr vnp;
  Int4       row_num;
  CharPtr    pos_str;
  Boolean    rval = TRUE;
  
  tlp = (TagListPtr) GetObjectExtra (d);
  
  if (tlp == NULL)
  {
    return FALSE;
  }
  
  sprintf (seq_str, "%d", seq_num + 1);
  for (vnp = tlp->vnp, row_num = 0; 
       vnp != NULL && rval; 
       vnp = vnp->next, row_num++)
  {
    pos_str = GetTagListValueEx (tlp, row_num, 1);
    if (StringCmp (pos_str, seq_str) == 0)
    {
      rval = FALSE;
    }
    pos_str = MemFree (pos_str);
  }
  return rval;
}

static Boolean SomeErrorMessagesAreHidden (SeqIdEditPtr siep)
{
  Boolean has_dups, has_missing, has_bracket, has_unrec_mods = FALSE;
  ValNodePtr unrec_mods = NULL;
  Int4    seq_num;
  CharPtr err_msg;
  Boolean some_hidden = FALSE;
  
  if (siep == NULL)
  {
    return FALSE;
  }

  has_dups = EditHasDuplicateIDs (siep->iatep_new, siep->iatep_current);
  has_missing = HasMissingIDs (siep->iatep_new) || HasMissingIDs (siep->iatep_current);
  has_bracket = EditNeedsBracketingFixes (siep->iatep_new) 
                || EditNeedsBracketingFixes (siep->iatep_current);

  if (!has_bracket)
  {
    unrec_mods = ListUnrecognizedModifiers (siep->iatep_new, siep->iatep_current, siep->is_nuc);
  }
  if (unrec_mods != NULL)
  {
    has_unrec_mods = TRUE;
    unrec_mods = ValNodeFreeData (unrec_mods);
  }
  
  
  if (siep->iatep_new != NULL)
  {
    for (seq_num = 0; seq_num < siep->iatep_new->num_sequences && !some_hidden; seq_num++)
    {
      if (IsSeqNumHidden (siep->new_dlg, seq_num))
      {
        err_msg = GetIDAndTitleErrorMessage (siep->iatep_new, siep->iatep_current, 
                                         seq_num, has_dups,
                                         has_missing, siep->seqid_edit_phase,
                                         has_bracket,
                                         has_unrec_mods,
                                         siep->is_nuc);
        if (!StringHasNoText (err_msg))
        {
          some_hidden = TRUE;
        }
        err_msg = MemFree (err_msg);
      }
    }
  }
  if (siep->iatep_current != NULL)
  {
    for (seq_num = 0; seq_num < siep->iatep_current->num_sequences && !some_hidden; seq_num++)
    {
      if (IsSeqNumHidden (siep->current_dlg, seq_num))
      {
        err_msg = GetIDAndTitleErrorMessage (siep->iatep_current, siep->iatep_new, 
                                         seq_num, has_dups,
                                         has_missing, siep->seqid_edit_phase,
                                         has_bracket,
                                         has_unrec_mods,
                                         siep->is_nuc);
        if (!StringHasNoText (err_msg))
        {
          some_hidden = TRUE;
        }
        err_msg = MemFree (err_msg);
      }
    }
  }
  return some_hidden;
}

static void ShowErrorInstructions (Pointer userdata)
{
  SeqIdEditPtr      siep;
  Boolean           has_missing, has_dups, has_bracket_ids, has_space;
  Boolean           old_seqid_edit_phase;
  Boolean           show_all;
    
  siep = (SeqIdEditPtr) userdata;
  if (siep == NULL)
  {
    return;
  }
  
  UpdateIdAndTitleData (siep->new_dlg, siep->iatep_new);
  UpdateIdAndTitleData (siep->current_dlg, siep->iatep_current);
  
  has_missing = HasMissingIDs (siep->iatep_new) || HasMissingIDs (siep->iatep_current);
  has_dups = EditHasDuplicateIDs (siep->iatep_new, siep->iatep_current);
  has_bracket_ids = AnyBracketsInIDs (siep->iatep_new) || AnyBracketsInIDs (siep->iatep_current);
  has_space = EditHasSpaceInIDs (siep->iatep_new) || EditHasSpaceInIDs (siep->iatep_current);
  
  old_seqid_edit_phase = siep->seqid_edit_phase;
  siep->seqid_edit_phase |= has_missing | (has_dups && siep->is_nuc) | has_bracket_ids;
    
  if (siep->seqid_edit_phase)
  {
    SetTitle (siep->w, "Provide Sequence IDs For Your Sequences");
    Show (siep->refresh_err_btn);
  }
  else
  {
    SetTitle (siep->w, "Provide Correctly Formatted Titles For Your Sequences");
    Hide (siep->refresh_err_btn);
  }
  
  if (siep->seqid_edit_phase && !old_seqid_edit_phase)
  {
    if (GetValue (siep->show_all_grp) == 2)
    {
      show_all = TRUE;
    }
    else
    {
      show_all = FALSE;
    }  

    UpdateIdAndTitleEditDialog (siep->new_dlg, 
                                siep->iatep_new, siep->iatep_current,
                                siep->seqid_edit_phase,
                                show_all,
                                siep->is_nuc);
    UpdateIdAndTitleEditDialog (siep->current_dlg, 
                                siep->iatep_current, siep->iatep_new,
                                siep->seqid_edit_phase,
                                show_all,
                                siep->is_nuc);
                                  
  }

  
  Reset (siep->auto_correct_doc);
  Show (siep->auto_correct_doc);
  Hide (siep->bracket_dlg);
  Hide (siep->badvalue_pnl);
  Hide (siep->unrec_mod_pnl);
  if (has_missing || (has_dups && siep->is_nuc) || has_space)
  {
    if (has_missing)
    {
      AppendText (siep->auto_correct_doc, "Some of your sequences lack sequence IDs.", NULL, NULL, programFont);
    }
    if (has_dups && siep->is_nuc)
    {
      AppendText (siep->auto_correct_doc, "Some of your sequence IDs are duplicated.", NULL, NULL, programFont);
    }
    AppendText (siep->auto_correct_doc, "Please provide unique sequence IDs for every sequence.", NULL, NULL, programFont);

    if (has_space)
    {
      AppendText (siep->auto_correct_doc, "Some of your sequence IDs contain spaces.", NULL, NULL, programFont);
    }
    
    if (SomeErrorMessagesAreHidden (siep))
    {
      AppendText (siep->auto_correct_doc, "Press 'Refresh Error List' to see the complete list of sequences with errors.", NULL, NULL, programFont);
    }
    
    if (ShowExtendedIDCorrections (siep->iatep_new, siep->iatep_current, siep->auto_correct_doc))
    {
      siep->auto_correct_ids = TRUE;
      SetTitle (siep->auto_correct_btn, "Autocorrect Sequence IDs");
      Show (siep->auto_correct_btn);
    }
    else
    {
      siep->auto_correct_ids = FALSE;
      Hide (siep->auto_correct_btn);
    }
    siep->auto_correct_bracketing = FALSE;
    siep->auto_correct_modnames = FALSE;
  }
  else if (ShowExtendedIDCorrections (siep->iatep_new, siep->iatep_current, siep->auto_correct_doc))
  {
    siep->auto_correct_ids = TRUE;
    SetTitle (siep->auto_correct_btn, "Autocorrect Sequence IDs");
    Show (siep->auto_correct_btn);
    siep->auto_correct_bracketing = FALSE;
    siep->auto_correct_modnames = FALSE;
  }
  else if (SequenceIDsHaveNonFixableBrackets (siep->iatep_new, siep->iatep_current, siep->auto_correct_doc))
  {
    Hide (siep->auto_correct_btn);
  }
  else if (siep->seqid_edit_phase && EditIDsNeedFix (siep))
  {
    siep->auto_correct_ids = FALSE;
    siep->auto_correct_bracketing = FALSE;
    siep->auto_correct_modnames = FALSE;
    Reset (siep->auto_correct_doc);
    AppendText (siep->auto_correct_doc, "Sequence ID errors have been corrected.\nErrors are present within your sequence titles.\nClick 'Proceed to Title Correction' to correct sequence title errors.\n", NULL, NULL, programFont);
    SetTitle (siep->auto_correct_btn, "Proceed to Title Correction");
    Show (siep->auto_correct_btn);
  }
  else if (ShowBracketingCorrections (siep->iatep_new, siep->iatep_current, siep->bracket_dlg))
  {
    Show (siep->bracket_dlg);
    Hide (siep->auto_correct_doc);
    SetTitle (siep->auto_correct_btn, "Autocorrect Bracketing");
    Show (siep->auto_correct_btn);
    siep->auto_correct_ids = FALSE;
    siep->auto_correct_bracketing = TRUE;
    siep->auto_correct_modnames = FALSE;
    UpdateIDAndTitleEditDialogErrorColumns (siep->new_dlg, siep->current_dlg, 
                                            siep->iatep_new, siep->iatep_current,
                                            siep->seqid_edit_phase, siep->is_nuc);

  }
  else if (ShowUnrecognizedModifiers (siep->iatep_new, siep->iatep_current, siep->is_nuc))
  {
    Hide (siep->auto_correct_doc);
    Show (siep->unrec_mod_pnl);

    Show (siep->auto_correct_btn);
    SetTitle (siep->auto_correct_btn, "Correct Modifier Names");
    siep->auto_correct_ids = FALSE;
    siep->auto_correct_bracketing = FALSE;
    siep->auto_correct_modnames = TRUE;
    UpdateIDAndTitleEditDialogErrorColumns (siep->new_dlg, siep->current_dlg, 
                                            siep->iatep_new, siep->iatep_current,
                                            siep->seqid_edit_phase, siep->is_nuc);
  }
  else if (ShowBadValues (siep->iatep_new, siep->iatep_current))
  {
    Show (siep->badvalue_pnl);
    Hide (siep->auto_correct_doc);
    Hide (siep->auto_correct_btn);
  }
  else
  {
    AppendText (siep->auto_correct_doc, "Sequence ID and title errors have been corrected.", NULL, NULL, programFont);
    
    Hide (siep->auto_correct_btn);
    siep->auto_correct_ids = FALSE;
    siep->auto_correct_bracketing = FALSE;
    siep->auto_correct_modnames = TRUE;
  }
  UpdateDocument (siep->auto_correct_doc, 0, 0);
  UpdateIDAndTitleEditDialogErrorColumns (siep->new_dlg, siep->current_dlg, 
                                          siep->iatep_new, siep->iatep_current,
                                          siep->seqid_edit_phase, siep->is_nuc); 

  if (IDAndTitleEditIDsNeedFix (siep->iatep_new, siep->iatep_current, siep->is_nuc))
  {
    Disable (siep->accept_btn);
  }
  else
  {
    Enable (siep->accept_btn);
  }
}
 
static TaglistCallback callback_list[4] = 
 { ShowErrorInstructions, ShowErrorInstructions, ShowErrorInstructions, ShowErrorInstructions };

static GrouP CreateIDsAndTitlesDialog (GrouP parent, SeqIdEditPtr siep, Boolean is_new)
{
  GrouP      g, k, ppt;
  PrompT     p_desc = NULL, p1, p2, p3, p4;
  Int4       num_sequences;
  TagListPtr tlp;
  DialoG     dlg;
  
  if (siep == NULL)
  {
    return NULL;
  }
  
  g = HiddenGroup (parent, -1, 0, NULL);

  if (is_new)
  {
    if (siep->iatep_current != NULL)
    {
      p_desc = StaticPrompt (g, "New Sequences", 0, 0, programFont, 'l');
    }
  }
  else
  {
    p_desc = StaticPrompt (g, "Existing Sequences", 0, 0, programFont, 'l');
  }
  
  ppt = HiddenGroup (g, 4, 0, NULL);
  p1 = StaticPrompt (ppt, "Error", 6 * stdCharWidth, 0, programFont, 'l');
  p2 = StaticPrompt (ppt, "Position", 5 * stdCharWidth, 0, programFont, 'l');
  p3 = StaticPrompt (ppt, "Sequence ID", 10 * stdCharWidth, 0, programFont, 'l');
  p4 = StaticPrompt (ppt, "Title", 40 * stdCharWidth, 0, programFont, 'l');

  if (is_new)
  {
    num_sequences = siep->iatep_new->num_sequences;
  }
  else
  {
    num_sequences = siep->iatep_current->num_sequences;
  }


  k = NormalGroup (g, 1, 0, "", programFont, NULL);
  dlg = CreateTagListDialogExEx (k, MIN (num_sequences, 4), 4, 2,
                                 idedit_types, idedit_widths,
                                 NULL, TRUE, TRUE, NULL, NULL,
                                 callback_list, siep, TRUE);

  tlp = (TagListPtr) GetObjectExtra (dlg);  
  if (tlp == NULL) return NULL;
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg, (HANDLE) p_desc, NULL);
  
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [0], (HANDLE) p1, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [1], (HANDLE) p2, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [2], (HANDLE) p3, NULL);
  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [3], (HANDLE) p4, NULL);
  
  if (is_new)
  {
    siep->new_dlg = dlg;
  }
  else
  {
    siep->current_dlg = dlg;
  }

  return g;
}

static void ScrollTagListToSeqId (DialoG d, CharPtr seq_id)
{
  TagListPtr tlp;
  ValNodePtr vnp;
  Int4       row_num, sb_max;
  Int4       scroll_value = -1;
  CharPtr    id_from_tag;
  
  tlp = (TagListPtr) GetObjectExtra (d); 
  
  if (tlp == NULL || StringHasNoText (seq_id))
  {
    return;
  }
  
  for (row_num = 0, vnp = tlp->vnp;
       vnp != NULL && scroll_value < 0;
       vnp = vnp->next, row_num++)
  {
    id_from_tag = GetTagListValueEx (tlp, row_num, 2);
    if (StringCmp (id_from_tag, seq_id) == 0)
    {
      scroll_value = row_num;
    }
    id_from_tag = MemFree (id_from_tag);
  }
  
  if (scroll_value < 0)
  {
    return;
  }
  
  sb_max = GetBarMax (tlp->bar);
  
  if (scroll_value >= sb_max)
  {
    SetBarValue (tlp->bar, sb_max);
  }
  else
  {
    SetBarValue (tlp->bar, scroll_value);
  }
  SendMessageToDialog (d, VIB_MSG_REDRAW);
}

static void ShowAllSequences (GrouP g)
{
  SeqIdEditPtr siep;
  Boolean      show_all = FALSE;
  
  siep = (SeqIdEditPtr) GetObjectExtra (g);
  if (siep == NULL)
  {
    return;
  }
  
  show_all = HasMissingIDs (siep->iatep_new) || HasMissingIDs (siep->iatep_current);
  
  if (show_all)
  {
    SetValue (siep->show_all_grp, 2);
    Disable (siep->show_all_grp);
  }
  else
  {
    Enable (siep->show_all_grp);
  }
  
  if (GetValue (siep->show_all_grp) == 2)
  {
    show_all = TRUE;
  }

  if (siep->new_dlg != NULL)
  {
    UpdateIdAndTitleData (siep->new_dlg, siep->iatep_new);
    UpdateIdAndTitleEditDialog (siep->new_dlg, siep->iatep_new, siep->iatep_current, 
                                siep->seqid_edit_phase, show_all,
                                siep->is_nuc);
  }
  
  if (siep->current_dlg != NULL)
  {
    UpdateIdAndTitleData (siep->current_dlg, siep->iatep_current);
    UpdateIdAndTitleEditDialog (siep->current_dlg, siep->iatep_current, siep->iatep_new,
                                siep->seqid_edit_phase, show_all,
                                siep->is_nuc);
  }

  ShowErrorInstructions (siep);
}

static void RefreshErrorList (ButtoN b)
{
  SeqIdEditPtr siep;
  
  siep = (SeqIdEditPtr) GetObjectExtra (b);
  if (siep == NULL)
  {
    return;
  }
  ShowAllSequences (siep->show_all_grp);
}

static void AutoCorrectIDsAndTitles (ButtoN b)
{
  SeqIdEditPtr      siep;
  IDAndTitleEditPtr suggested, suggested_new, suggested_current;
  Boolean           show_all;
  
  siep = (SeqIdEditPtr) GetObjectExtra (b);
  if (siep == NULL)
  {
    return;
  }
  
  if (siep->auto_correct_ids)
  {
    suggested_new = SuggestCorrectionForLocalIDs (siep->iatep_new, siep->iatep_current);
    suggested_current = SuggestCorrectionForLocalIDs (siep->iatep_current, siep->iatep_new);
    UpdateColorizedDeflinePanelData (suggested_new, suggested_current, siep->badvalue_pnl);
    UpdateColorizedDeflinePanelData (suggested_new, suggested_current, siep->unrec_mod_pnl);

    siep->iatep_new = IDAndTitleEditFree (siep->iatep_new);
    siep->iatep_new = suggested_new;
    
    if (siep->iatep_current == NULL)
    {
      suggested_current = IDAndTitleEditFree (siep->iatep_current);
    }
    else
    {
      siep->iatep_current = IDAndTitleEditFree (siep->iatep_current);
      siep->iatep_current = suggested_current;
    }
  }
  else if (siep->seqid_edit_phase)
  {
    siep->seqid_edit_phase = FALSE;
  }
  else if (siep->auto_correct_bracketing)
  {
    suggested = SuggestCorrectionForTitleBracketing (siep->iatep_new);
    UpdateColorizedDeflinePanelData (suggested, siep->iatep_current, siep->badvalue_pnl);
    UpdateColorizedDeflinePanelData (suggested, siep->iatep_current, siep->unrec_mod_pnl);

    siep->iatep_new = IDAndTitleEditFree (siep->iatep_new);
    siep->iatep_new = suggested;
    if (siep->iatep_current != NULL)
    {
      suggested = SuggestCorrectionForTitleBracketing (siep->iatep_current);
      UpdateColorizedDeflinePanelData (siep->iatep_new, suggested, siep->badvalue_pnl);
      UpdateColorizedDeflinePanelData (siep->iatep_new, suggested, siep->unrec_mod_pnl);

      siep->iatep_current = IDAndTitleEditFree (siep->iatep_current);
      siep->iatep_current = suggested;      
    }
  }
  else if (siep->auto_correct_modnames)
  {
    FixAllUnrecognizedModifiers (siep->iatep_new, siep->iatep_current, siep);
  }

  show_all = HasMissingIDs (siep->iatep_new) || HasMissingIDs (siep->iatep_current);
  if (show_all)
  {
    Disable (siep->show_all_grp);
  }
  else
  {
    Enable (siep->show_all_grp);
  }

  if (GetValue (siep->show_all_grp) == 2)
  {
    show_all = TRUE;
  }

  UpdateIdAndTitleEditDialog (siep->new_dlg, siep->iatep_new, siep->iatep_current,
                              siep->seqid_edit_phase, show_all, siep->is_nuc);
  UpdateIdAndTitleEditDialog (siep->current_dlg, siep->iatep_current, siep->iatep_new, 
                              siep->seqid_edit_phase, show_all, siep->is_nuc);

  ShowErrorInstructions (siep);
}

static void RefreshErrorButton (ButtoN b)
{
  SeqIdEditPtr siep;
  
  siep = (SeqIdEditPtr) GetObjectExtra (b);
  if (siep == NULL)
  {
    return;
  }
  
  UpdateIdAndTitleData (siep->new_dlg, siep->iatep_new);
  ClearIDAndTitleEditDialogErrorColumn (siep->new_dlg);
  UpdateIdAndTitleData (siep->current_dlg, siep->iatep_current);
  ClearIDAndTitleEditDialogErrorColumn (siep->current_dlg);
  UpdateIDAndTitleEditDialogErrorColumns (siep->new_dlg, siep->current_dlg,
                                          siep->iatep_new, siep->iatep_current,
                                          siep->seqid_edit_phase, siep->is_nuc);
}

static void CleanupDuplicateAndEmptyPairs (IDAndTitleEditPtr iatep)
{
  Int4 seq_num;
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    iatep->title_list [seq_num] = RemoveAllDuplicatePairsFromOneTitle (iatep->title_list [seq_num]);
    RemoveMeaninglessEmptyPairsFromOneTitle (iatep->title_list [seq_num]);
  }
}

/* This function will insert the escape character before each double-quotation mark
 * in a string, starting at the position offset and ending at the portion of the
 * string that corresponds to the position stop_offset in the original string - 
 * this offset is adjusted for the additional escape characters inserted.
 */
static CharPtr EscapeQuotesBetweenPositions (CharPtr orig_title, Int4 start_offset, Int4 stop_offset)
{
  CharPtr next_quote;
  
  if (StringHasNoText (orig_title))
  {
    return orig_title;
  }
  
  next_quote = NextUnescapedQuote (orig_title + start_offset);
  while (next_quote != NULL && next_quote - orig_title < stop_offset)
  {
    orig_title = InsertStringAtOffset (orig_title, "\\", next_quote - orig_title);
    stop_offset ++;
    next_quote = NextUnescapedQuote (orig_title + start_offset);
  }
  return orig_title;
  
}

/* This function will insert the escape character before each double-quotation mark in a
 * string, starting at position offset.  This is useful when the first portion of a 
 * title appears to be parseable, but the remainder is not.
 */
static CharPtr EscapeQuotesAfterOffset (CharPtr orig_title, Int4 offset)
{
  CharPtr next_quote;
  
  if (StringHasNoText (orig_title))
  {
    return orig_title;
  }
  
  next_quote = NextUnescapedQuote (orig_title + offset);
  while (next_quote != NULL)
  {
    orig_title = InsertStringAtOffset (orig_title, "\\", next_quote - orig_title);
    next_quote = NextUnescapedQuote (orig_title + offset);
  }
  return orig_title;
}

/* This function inserts escape characters before each double-quotation mark in
 * a string starting at position offset, and then puts the portion of the string
 * starting at position offset inside a pair of double-quotes.
 */
static CharPtr QuoteToEndFromOffset (CharPtr orig_title, Int4 offset)
{
  if (StringHasNoText (orig_title))
  {
    return orig_title;
  }
  
  orig_title = EscapeQuotesAfterOffset (orig_title, offset);
  orig_title = InsertStringAtOffset (orig_title, "\"", offset);
  orig_title = InsertStringAtOffset (orig_title, "\"", StringLen (orig_title));
  return orig_title;
}

/* This function finds the first position where bracketing is incorrect
 * and puts the remainder of the string inside quotation marks.
 */
static CharPtr FixUnparseableBracketing (CharPtr orig_title)
{
  CharPtr next_token, next_next_token, last_start = NULL;
  Char    ch_expected;
  
  next_token = NextBracketToken (orig_title);
  if (next_token == NULL)
  {
    return orig_title;
  }
  
  if (*next_token == '[')
  {
    last_start = next_token;
  }
  else
  {
    orig_title = QuoteToEndFromOffset (orig_title, next_token - orig_title);
    return orig_title;
  }
  
  while (next_token != NULL)
  {
    ch_expected = ExpectToken (next_token);
    next_next_token = NextBracketToken (next_token + 1);
    if ((next_next_token == NULL && ch_expected != '[')
         || (next_next_token != NULL && *next_next_token != ch_expected))
    {
      orig_title = QuoteToEndFromOffset (orig_title, last_start - orig_title);
      return orig_title;
    }
    else
    {
      next_token = next_next_token;
      if (next_token != NULL && *next_token == '[')
      {
        last_start = next_token;
      }
    }
  }
  return orig_title;  
}

/* This function finds individual bracketed name-value pairs that either have
 * unrecognizable names or invalid values and puts them in quotation marks, so
 * that they will not be parsed.
 */
static CharPtr 
QuoteUnrecognizedModifierNamesAndValues 
(CharPtr orig_title,
 Boolean is_nuc)
{
  CharPtr         start, stop, cp;
  ModifierInfoPtr mip;
  Int4            stop_offset, start_offset;
  
  cp = orig_title;
  mip = ParseOneBracketedModifier (cp, &start, &stop);
  while (mip != NULL && stop != NULL && start != NULL)
  {
    if (IsUnrecognizedModifierName (mip, is_nuc))
    {
      /* note - put stop quote in first, because position of stop will change
       * after putting in quote for start */
      stop_offset = stop - orig_title + 1;
      start_offset = start - orig_title;
      if (is_nuc)
      {
        orig_title = InsertStringAtOffset (orig_title, "\"", stop_offset);
      }
      else
      {
        orig_title = InsertStringAtOffset (orig_title, "\"]", stop_offset);
      }
      orig_title = EscapeQuotesBetweenPositions (orig_title, start_offset, stop_offset);
      if (is_nuc)
      {
        orig_title = InsertStringAtOffset (orig_title, "\"", start_offset);
      }
      else
      {
        orig_title = InsertStringAtOffset (orig_title, "[comment=\"", start_offset);
      }
      cp = orig_title + start_offset;
    }
    else
    {
      cp = stop + 1;
    }
    mip = ModifierInfoFree (mip);    
    mip = ParseOneBracketedModifier (cp, &start, &stop);  
  }
  return orig_title;
}

/* This section will insert double-quotation marks around the sections of a title
 * that cannot be parsed.
 */
static void QuoteUnparseableSections (IDAndTitleEditPtr iatep, Boolean is_nuc)
{
  Int4 seq_num;
  
  if (iatep == NULL)
  {
    return;
  }
  
  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    iatep->title_list [seq_num] = FixUnparseableBracketing (iatep->title_list [seq_num]);
    iatep->title_list [seq_num] = QuoteUnrecognizedModifierNamesAndValues (iatep->title_list [seq_num], is_nuc);
  }  
}

static void PositionRefreshErrorListBtn (SeqIdEditPtr siep)
{
  RecT a, b;
  
  if (siep == NULL || siep->refresh_err_btn == NULL 
      || siep->refresh_error_list_btn == NULL)
  {
    return;
  }
  
  ObjectRect (siep->refresh_err_btn, &a);
  ObjectRect (siep->refresh_error_list_btn, &b);
  b.left += a.right - a.left + 10;
  b.right += a.right - a.left + 10;
  SetPosition (siep->refresh_error_list_btn, &b);
}

/* This function checks to see if all of the sequences in new_list and current_list
 * have non-empty, unique sequences and if any of their titles have bracketing errors,
 * invalid modifier names, or invalid modifier values.
 * Problems with sequence IDs must be fixed; problems with titles that are left
 * unfixed will be cordoned off with double-quotation marks.
 */
NLM_EXTERN Boolean FixIDsAndTitles (SeqEntryPtr new_list, SeqEntryPtr current_list, Boolean is_nuc)
{
  WindoW                w;
  GrouP                 h, instr_grp, k, j, c, new_grp, current_grp = NULL;
  ButtoN                b;
  Boolean               need_fix = FALSE;
  ModalAcceptCancelData acd;
  Boolean               show_all;
  SeqIdEditData         sied;
  Boolean               rval;
  
  sied.iatep_new = SeqEntryListToIDAndTitleEdit (new_list);
  sied.iatep_current = SeqEntryListToIDAndTitleEdit (current_list);
  sied.is_nuc = is_nuc;
    
  /* check for unique IDs - don't need to present dialog if they
   * are all present and unique */
  need_fix = EditIDsNeedFix (&sied);
  
  /* if no fixes are needed, do not present dialog */
  if (!need_fix)
  {
    sied.iatep_new = IDAndTitleEditFree (sied.iatep_new);
    sied.iatep_current = IDAndTitleEditFree (sied.iatep_current);
    return TRUE;
  }
  
  show_all = HasMissingIDs (sied.iatep_new) || HasMissingIDs (sied.iatep_current);
  
  if (!show_all 
      && (! EditHasDuplicateIDs (sied.iatep_new, sied.iatep_current) || !sied.is_nuc)
      && !AnyBracketsInIDs (sied.iatep_new) 
      && !AnyBracketsInIDs (sied.iatep_current))
  {
    sied.seqid_edit_phase = FALSE;
  }
  else
  {
    sied.seqid_edit_phase = TRUE;
  }

  
  w = MovableModalWindow (-20, -13, -10, -10, "Provide Sequence IDs for your Sequences", NULL);
  sied.w = w;
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  instr_grp = HiddenGroup (h, -1, 0, NULL);
  k = HiddenGroup (instr_grp, 0, 0, NULL);
  sied.auto_correct_doc = DocumentPanel (k, stdCharWidth * 63, stdLineHeight * 12);
  SetDocAutoAdjust (sied.auto_correct_doc, TRUE);
  SetObjectExtra (sied.auto_correct_doc, &sied, NULL);
  sied.bracket_dlg = ShowDifferenceDialog (k, stdCharWidth * 63, stdLineHeight * 12);
  sied.badvalue_pnl = InvalidValuesPanel (k, stdCharWidth * 63, stdLineHeight * 12, 
                                          sied.iatep_new, sied.iatep_current, is_nuc,
                                          UpdateSeqIdEditForColorizedPanel,
                                          &sied,
                                          ScrollSeqIdEditForColorizedPanel,
                                          &sied);
  sied.unrec_mod_pnl = UnrecognizedModifiersPanel (k, stdCharWidth * 63, stdLineHeight * 12, 
                                                   sied.iatep_new, sied.iatep_current, is_nuc,
                                                   UpdateSeqIdEditForColorizedPanel,
                                                   &sied,
                                                   ScrollSeqIdEditForColorizedPanel,
                                                   &sied);

  AlignObjects (ALIGN_CENTER, (HANDLE) sied.auto_correct_doc, 
                              (HANDLE) sied.bracket_dlg, 
                              (HANDLE) sied.badvalue_pnl,
                              (HANDLE) sied.unrec_mod_pnl,
                              NULL);

  j = HiddenGroup (h, 0, 0, NULL);
  sied.refresh_err_btn = PushButton (j, "Clear Fixed Errors", RefreshErrorButton);
  SetObjectExtra (sied.refresh_err_btn, &sied, NULL);
  
  sied.refresh_error_list_btn = PushButton (j, "Refresh Error List", RefreshErrorList);
  SetObjectExtra (sied.refresh_error_list_btn, &sied, NULL);

  sied.auto_correct_btn = PushButton (j, "Make automatic corrections", AutoCorrectIDsAndTitles);
  SetObjectExtra (sied.auto_correct_btn, &sied, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) sied.auto_correct_btn, NULL);

  sied.new_dlg = NULL;
  sied.current_dlg = NULL;
  
  new_grp = CreateIDsAndTitlesDialog (h, &sied, TRUE);
  if (sied.iatep_current != NULL)
  {
    current_grp = CreateIDsAndTitlesDialog (h, &sied, FALSE);
  }
  
  UpdateIdAndTitleEditDialog (sied.new_dlg, sied.iatep_new, sied.iatep_current,
                              sied.seqid_edit_phase, show_all, sied.is_nuc);
  UpdateIdAndTitleEditDialog (sied.current_dlg, sied.iatep_current, sied.iatep_new,
                              sied.seqid_edit_phase, show_all, sied.is_nuc);

  k = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  sied.show_all_grp = HiddenGroup (k, 2, 0, ShowAllSequences);
  RadioButton (sied.show_all_grp, "Show only sequences with errors");
  RadioButton (sied.show_all_grp, "Show all sequences in set");
  SetValue (sied.show_all_grp, 1);
  SetObjectExtra (sied.show_all_grp, &sied, NULL);
  if (show_all)
  {
    Disable (sied.show_all_grp);
  }
  else
  {
    Enable (sied.show_all_grp);
  }
  
  
  c = HiddenGroup (h, 2, 0, NULL);
  sied.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (sied.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  ShowErrorInstructions (&sied);
    
  AlignObjects (ALIGN_CENTER, (HANDLE) instr_grp,
                              (HANDLE) k, 
                              (HANDLE) c, 
                              (HANDLE) new_grp,
                              (HANDLE) current_grp,
                              NULL);
                              
  AlignObjects (ALIGN_LEFT, (HANDLE) sied.refresh_err_btn,
                            (HANDLE) sied.new_dlg,
                            NULL);                              

  PositionRefreshErrorListBtn (&sied);

  Show (w);
  Select (w);
  
  acd.cancelled = FALSE;
  while (need_fix && ! acd.cancelled)
  {
    acd.accepted = FALSE;
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
    if (! acd.cancelled)
    {
      UpdateIdAndTitleData (sied.new_dlg, sied.iatep_new);
      UpdateIdAndTitleData (sied.current_dlg, sied.iatep_current);
     
      if (IDAndTitleEditTitlesNeedFix (sied.iatep_new, sied.is_nuc)
          || IDAndTitleEditTitlesNeedFix (sied.iatep_current, sied.is_nuc))
      {
        if (ANS_YES == Message (MSG_YN, "Your titles contain unparseable elements.  If you choose to continue, the unparseable sections will be placed in quotes and no data will be parsed from these sections.  You will have to add information about the sequence manually.  Are you sure you want to continue?"))
        {
          QuoteUnparseableSections (sied.iatep_new, is_nuc);
          QuoteUnparseableSections (sied.iatep_current, is_nuc);
          need_fix = FALSE;
        }
      }
      else
      {
        need_fix = FALSE;
      }
             
      show_all = HasMissingIDs (sied.iatep_new) || HasMissingIDs (sied.iatep_current);
      if (show_all)
      {
        Disable (sied.show_all_grp);
      }
      else
      {
        Enable (sied.show_all_grp);
        if (GetValue (sied.show_all_grp) == 2)
        {
          show_all = TRUE;
        }
        else
        {
          show_all = FALSE;
        }
      }
      UpdateIdAndTitleEditDialog (sied.new_dlg, sied.iatep_new, sied.iatep_current,
                                  sied.seqid_edit_phase, show_all, sied.is_nuc);
      UpdateIdAndTitleEditDialog (sied.current_dlg, sied.iatep_current, sied.iatep_new,
                                  sied.seqid_edit_phase, show_all, sied.is_nuc);
      ShowErrorInstructions (&sied);
    }
  }
    
  if (acd.cancelled)
  {
    rval = FALSE;
  }
  else
  {
    CleanupDuplicateAndEmptyPairs (sied.iatep_new);
    CleanupDuplicateAndEmptyPairs (sied.iatep_current);
    ApplyIDAndTitleEditToSeqEntryList (new_list, sied.iatep_new);
    ApplyIDAndTitleEditToSeqEntryList (current_list, sied.iatep_current);
    rval = TRUE;
  }
  Remove (w);
  
  sied.iatep_new = IDAndTitleEditFree (sied.iatep_new);
  sied.iatep_current = IDAndTitleEditFree (sied.iatep_current);
  
  return rval;  
}

static Boolean CollectIDsAndTitles (SeqEntryPtr new_list, SeqEntryPtr current_list, Boolean is_nuc)
{

  ArrowCursor ();
  
  return FixIDsAndTitles (new_list, current_list, is_nuc);
}

/* The following section of code is used for editing titles with the Sequence Assistant.
 * A user may edit a single title or the list of titles from the Sequence Assistant.
 * This code provides the same bracketing, modifier name, and modifier value checks
 * that are used when sequences are imported or created.
 */
 
typedef struct titleedit
{
  DialoG            bracket_dlg;
  PaneL             unrec_mod_pnl;
  PaneL             badvalue_pnl;
  TexT              title_txt;
  DialoG            multi_title;
  IDAndTitleEditPtr iatep;
  ButtoN            accept_btn;
} TitleEditData, PNTR TitleEditPtr;

static void ShowTitleEditErrors (TitleEditPtr tep)
{
  if (tep == NULL)
  {
    return;
  }
  
  Hide (tep->bracket_dlg);
  Hide (tep->unrec_mod_pnl);
  Hide (tep->badvalue_pnl);
  
  if (ShowBracketingCorrections (tep->iatep, NULL, tep->bracket_dlg))
  {
    Show (tep->bracket_dlg);
    Disable (tep->accept_btn);
  }
  else if (ShowUnrecognizedModifiers (tep->iatep, NULL, TRUE))
  {
    Show (tep->unrec_mod_pnl);
    Disable (tep->accept_btn);
  }
  else if (ShowBadValues (tep->iatep, NULL))
  {
    Show (tep->badvalue_pnl);
    Disable (tep->accept_btn);
  }
  else
  {
    Enable (tep->accept_btn);
  }  
}

static void OnTitleEditChange (TexT t)
{
  TitleEditPtr tep;
  
  tep = (TitleEditPtr) GetObjectExtra (t);
  if (tep == NULL)
  {
    return;
  }
  
  tep->iatep->title_list [0] = MemFree (tep->iatep->title_list [0]);
  tep->iatep->title_list [0] = SaveStringFromText (tep->title_txt);

  ShowTitleEditErrors (tep);
}

static void UpdateTitleEditForColorizedPanel (Pointer userdata)
{
  TitleEditPtr tep;
  
  tep = (TitleEditPtr) userdata;
  if (tep == NULL)
  {
    return;
  }

  SetTitle (tep->title_txt, tep->iatep->title_list [0]);  
  ShowTitleEditErrors (tep);
}

static void EditOneSequenceTitle (SequenceAssistantPtr sap, Int4 seq_num)
{
  Int4                  seq_pos;
  SeqEntryPtr           sep, nsep;
  WindoW                w;
  CharPtr               title = NULL;
  BioseqPtr             bsp = NULL;
  CharPtr               title_fmt = "Title for %s";
  GrouP                 h, g, c, err_grp;
  ButtoN                b;
  ModalAcceptCancelData acd;
  SeqDescrPtr           sdp;
  TitleEditData         ted;
  
  if (sap == NULL)
  {
    return;
  }
  
  for (seq_pos = 0, sep = sap->seq_list;
       seq_pos != seq_num && sep != NULL;
       seq_pos ++, sep = sep->next)
  {
  }
  if (sep == NULL)
  {
    return;
  }
  
  if (IS_Bioseq (sep)) 
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  } 
  else if (IS_Bioseq_set (sep)) 
  {
    nsep = FindNucSeqEntry (sep);
    if (nsep != NULL && IS_Bioseq (nsep)) 
    {
      bsp = (BioseqPtr) nsep->data.ptrvalue;
    }
  }
  if (bsp == NULL)
  {
    return;
  }
  
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  
  ted.iatep = IDAndTitleEditNew ();
  if (ted.iatep == NULL)
  {
    return;
  }
  
  /* set up IDAndTitleEdit to use for single sequence */
  ted.iatep->num_sequences = 1;
  ted.iatep->id_list = (CharPtr PNTR) MemNew (ted.iatep->num_sequences * sizeof (CharPtr));
  ted.iatep->title_list = (CharPtr PNTR) MemNew (ted.iatep->num_sequences * sizeof (CharPtr));
  ted.iatep->id_list [0] = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);
  if (sdp != NULL && !StringHasNoText (sdp->data.ptrvalue))
  {
    ted.iatep->title_list [0] = StringSave (sdp->data.ptrvalue);
  }
  else
  {
    ted.iatep->title_list [0] = StringSave ("");
  }  
  
  title = (CharPtr) MemNew ((StringLen (title_fmt) 
                               + StringLen (ted.iatep->id_list [0])) * sizeof (Char));
  sprintf (title, title_fmt, ted.iatep->id_list [0]);
  w = MovableModalWindow (-20, -13, -10, -10, title, NULL);
  title = MemFree (title);
  
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  err_grp = HiddenGroup (h, 0, 0, NULL);
  ted.bracket_dlg = ShowDifferenceDialog (err_grp, stdCharWidth * 63, stdLineHeight * 5);
  ted.badvalue_pnl = InvalidValuesPanel (err_grp, stdCharWidth * 63, stdLineHeight * 5, 
                                         ted.iatep, NULL, TRUE,
                                         UpdateTitleEditForColorizedPanel,
                                         &ted,
                                         NULL, NULL);
  ted.unrec_mod_pnl = UnrecognizedModifiersPanel (err_grp, stdCharWidth * 63, stdLineHeight * 5, 
                                                  ted.iatep, NULL, TRUE,
                                                  UpdateTitleEditForColorizedPanel,
                                                  &ted,
                                                  NULL, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ted.bracket_dlg,
                              (HANDLE) ted.badvalue_pnl,
                              (HANDLE) ted.unrec_mod_pnl,
                              NULL);
  
  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Title", 0, popupMenuHeight, programFont, 'l');
  ted.title_txt = DialogText (g, "", 40, OnTitleEditChange);
  SetObjectExtra (ted.title_txt, &ted, NULL);
  if (sdp != NULL && !StringHasNoText (sdp->data.ptrvalue))
  {
    SetTitle (ted.title_txt, sdp->data.ptrvalue);
  }
  
  c = HiddenGroup (h, 2, 0, NULL);
  ted.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (ted.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) err_grp, (HANDLE) g, (HANDLE) c, NULL);

  Show(w); 
  Select (w);
  
  ShowTitleEditErrors (&ted);
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  while (!acd.cancelled && ! acd.accepted)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.accepted)
  {
    CleanupDuplicateAndEmptyPairs (ted.iatep);
    SetTitle (ted.title_txt, ted.iatep->title_list [0]);
    if (sdp == NULL)
    {
      sdp = SeqDescrNew (bsp->descr);
      if (bsp->descr == NULL)
      {
        bsp->descr = sdp;
      }
      sdp->choice = Seq_descr_title;
    }
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = SaveStringFromText (ted.title_txt);
  }
  
  Remove (w);
}

#define SEQUENCE_ASSISTANT_MOLECULE_COLUMN 2
#define SEQUENCE_ASSISTANT_TOPOLOGY_COLUMN 3
#define SEQUENCE_ASSISTANT_TITLE_COLUMN    4

static Uint2 titleedit_types [] = {
  TAGLIST_PROMPT, TAGLIST_TEXT
};

static Uint2 titleedit_widths [] = {
  10, 40,
};

static void IDAndTitleEditToTagData (DialoG d, Pointer userdata) 
{
  Int4              len;
  CharPtr           str;
  IDAndTitleEditPtr iatep;
  TagListPtr        tlp;
  Int4              seq_num;
  ValNodePtr        list = NULL;
  
  tlp =(TagListPtr) GetObjectExtra (d);
  if (tlp == NULL)
  {
    return;
  }
  
  iatep = (IDAndTitleEditPtr) userdata;
  if (iatep == NULL)
  {
    return;
  }

  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    len = StringLen (iatep->id_list [seq_num]) + StringLen (iatep->title_list [seq_num]) + 4;
    str = (CharPtr) MemNew (len * sizeof (Char));
    if (str != NULL) {
      sprintf (str, "%s\t%s\n", iatep->id_list [seq_num], 
                                iatep->title_list [seq_num] == NULL ? "" : iatep->title_list [seq_num]);
    }
    ValNodeAddPointer (&list, 0, str);
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = list;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  tlp->max = MAX ((Int2) 0, (Int2) (iatep->num_sequences - tlp->rows));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  if (tlp->max > 0) {
    SafeShow (tlp->bar);
  } else {
    SafeHide (tlp->bar);
  }  
}

static Pointer TagDataToIDAndTitleEdit (DialoG d)
{
  IDAndTitleEditPtr iatep;
  TagListPtr        tlp;
  Int4              seq_num;
  
  tlp =(TagListPtr) GetObjectExtra (d);
  if (tlp == NULL)
  {
    return NULL;
  }
  
  iatep = IDAndTitleEditNew ();
  if (iatep == NULL)
  {
    return NULL;
  }
  
  iatep->num_sequences = ValNodeLen (tlp->vnp);
  iatep->id_list = (CharPtr PNTR) MemNew (iatep->num_sequences * sizeof (CharPtr));
  iatep->title_list = (CharPtr PNTR) MemNew (iatep->num_sequences * sizeof (CharPtr));

  for (seq_num = 0; seq_num < iatep->num_sequences; seq_num++)
  {
    iatep->id_list [seq_num] = GetTagListValueEx (tlp, seq_num, 0);
    iatep->title_list [seq_num] = GetTagListValueEx (tlp, seq_num, 1);
  }
  return iatep;  
}

static void UpdateMultiTitleEditForColorizedPanel (Pointer userdata)
{
  TitleEditPtr tep;
  
  tep = (TitleEditPtr) userdata;
  if (tep == NULL)
  {
    return;
  }
  
  PointerToDialog (tep->multi_title, tep->iatep);

  ShowTitleEditErrors (tep);
}

static void ShowMultiTitleErrors (Pointer userdata)
{
  TitleEditPtr tep;
  
  tep = (TitleEditPtr) userdata;
  
  if (tep == NULL)
  {
    return;
  }
  
  tep->iatep = IDAndTitleEditFree (tep->iatep);
  tep->iatep = DialogToPointer (tep->multi_title);
  UpdateColorizedDeflinePanelData (tep->iatep, NULL, tep->badvalue_pnl);
  UpdateColorizedDeflinePanelData (tep->iatep, NULL, tep->unrec_mod_pnl);
  
  ShowTitleEditErrors (tep);
}  

static TaglistCallback title_callback_list[2] = 
 { ShowMultiTitleErrors, ShowMultiTitleErrors };
 
static void EditSequenceTitleColumns (SequenceAssistantPtr sap)
{
  WindoW                w;
  GrouP                 h, err_grp, c;
  PrompT                ppt;
  ButtoN                b;
  Int4                  rows_shown = 0;
  TagListPtr            tlp;
  ModalAcceptCancelData acd;
  TitleEditData         ted;
  
  if (sap == NULL || sap->seq_list == NULL)
  {
    return;
  }

  ted.iatep = SeqEntryListToIDAndTitleEdit (sap->seq_list);
  if (ted.iatep == NULL)
  {
    return;
  }

  rows_shown = MIN (ted.iatep->num_sequences, 5);
  
  w = MovableModalWindow (-20, -13, -10, -10, "Sequence Titles", NULL);
  
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  err_grp = HiddenGroup (h, 0, 0, NULL);
  ted.bracket_dlg = ShowDifferenceDialog (err_grp, stdCharWidth * 63, stdLineHeight * 5);
  ted.badvalue_pnl = InvalidValuesPanel (err_grp, stdCharWidth * 63, stdLineHeight * 5, 
                                         ted.iatep, NULL, TRUE,
                                         UpdateMultiTitleEditForColorizedPanel,
                                         &ted,
                                         NULL, NULL);
  ted.unrec_mod_pnl = UnrecognizedModifiersPanel (err_grp, stdCharWidth * 63, stdLineHeight * 5, 
                                                  ted.iatep, NULL, TRUE,
                                                  UpdateMultiTitleEditForColorizedPanel,
                                                  &ted,
                                                  NULL, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ted.bracket_dlg,
                              (HANDLE) ted.badvalue_pnl,
                              (HANDLE) ted.unrec_mod_pnl,
                              NULL);  
  
  ppt = StaticPrompt (h, "Title", 18 * stdCharWidth, 0, programFont, 'l');    
  
  ted.multi_title = CreateTagListDialogExEx (h, rows_shown, 2, 2,
                                           titleedit_types, titleedit_widths,
                                           NULL, TRUE, TRUE, 
                                           IDAndTitleEditToTagData, 
                                           TagDataToIDAndTitleEdit,
                                           title_callback_list, &ted, FALSE);

  PointerToDialog (ted.multi_title, ted.iatep);
  
  tlp = (TagListPtr) GetObjectExtra (ted.multi_title);  
  if (tlp == NULL) return;
    

  c = HiddenGroup (h, 2, 0, NULL);
  ted.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (ted.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) ted.multi_title, (HANDLE) c, (HANDLE) NULL);

  AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [1], (HANDLE) ppt, NULL);

  ShowMultiTitleErrors (&ted);

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
  Hide (w);
  if (acd.accepted)
  {
    ted.iatep = IDAndTitleEditFree (ted.iatep);
    ted.iatep = DialogToPointer (ted.multi_title);
    ApplyIDAndTitleEditToSeqEntryList (sap->seq_list, ted.iatep);
    UpdateSequenceAssistant (sap);
  }
  ted.iatep = IDAndTitleEditFree (ted.iatep);
  Remove (w);
}


/* This section of code is for importing sequences from a file or creating new sequences.
 */
 
NLM_EXTERN SeqEntryPtr GetSequencesFromFile (CharPtr path, SeqEntryPtr current_list) 
{
  FILE         *fp;
  SeqEntryPtr  new_sep_list, new_sep, test_sep;
  Boolean      cancelled = FALSE;
  Boolean      chars_stripped = FALSE;

  fp = FileOpen (path, "r");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return NULL;
  }
      
  new_sep_list = ImportSequencesFromFile (fp, NULL, TRUE, TRUE, NULL, NULL, &chars_stripped);
  if (chars_stripped && new_sep_list != NULL)
  {
    if (ANS_CANCEL == Message (MSG_OKC, "Illegal characters will be stripped from your sequence data.  Do you want to continue?"))
    {
      new_sep_list = SeqEntryFree (new_sep_list);
      return NULL;
    }
  }

  if (new_sep_list == NULL)
  {
    Message (MSG_ERROR, "Unable to read sequences");
    return NULL;
  }
  else if (! RejectZeroLengthSequences (&new_sep_list))
  {
    return NULL;
  }
  else if (!CollectIDsAndTitles (new_sep_list, current_list, TRUE))
  { 
    new_sep = new_sep_list;   
    while (new_sep != NULL)
    {
      test_sep = new_sep->next;
      SeqEntryFree (new_sep);
      new_sep = test_sep;
    }
    FileClose (fp);
    return NULL;
  }
  
  if (cancelled)
  {
    new_sep = new_sep_list;
    while (new_sep != NULL)
    {
      test_sep = new_sep->next;
      SeqEntryFree (new_sep);
      new_sep = test_sep;
    }
    FileClose (fp);
    return NULL;
  }
     
  FileClose (fp);

  return new_sep_list;
}

static SeqEntryPtr GetSequencesFromText (TexT t, SeqEntryPtr current_list)
{
  CharPtr      seq_str;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep_list;
  FILE *fp;

  seq_str = SaveStringFromText (t);

  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return NULL;
  fprintf (fp, "%s", seq_str);
  FileClose (fp);
  
  seq_str = MemFree (seq_str);
  
  sep_list = GetSequencesFromFile (path, current_list);
  FileRemove (path);
  return sep_list;
}

static ValNodePtr ReadLinesOfFile (CharPtr path)
{
  ReadBufferData    rbd;
  CharPtr           line;
  FILE              *f;
  ValNodePtr        line_list = NULL;
  
  f = FileOpen (path, "r");
  if (f == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return NULL;
  }
  
  rbd.fp = f;
  rbd.current_data = NULL;
  line = AbstractReadFunction (&rbd);
  while (line != NULL)
  {
    ValNodeAddPointer (&line_list, 0, line);
    line = AbstractReadFunction (&rbd);
  }
  FileClose (f);
  return line_list;
}

static void AddLineListToText (ValNodePtr line_list, TexT t)
{
  CharPtr              old_seqstr, new_seqstr;
  Int4                 len;
  ValNodePtr           vnp;

  if (line_list == NULL || t == NULL)
  {
    return;
  }

  old_seqstr = SaveStringFromText (t);
  len = StringLen (old_seqstr) + 1;
  for (vnp = line_list; vnp != NULL; vnp = vnp->next)
  {
    len += StringLen (vnp->data.ptrvalue) + 3;
  }
  
  new_seqstr = (CharPtr) MemNew (len * sizeof (Char));
  if (new_seqstr != NULL)
  {
    StringCpy (new_seqstr, old_seqstr);
    StringCat (new_seqstr, "\n");
    for (vnp = line_list; vnp != NULL; vnp = vnp->next)
    {
      StringCat (new_seqstr, vnp->data.ptrvalue);
      StringCat (new_seqstr, "\n");
    }
    SetTitle (t, new_seqstr);
    new_seqstr = MemFree (new_seqstr);
  }
  old_seqstr = MemFree (old_seqstr);
}

static void AddSequenceImportFasta (ButtoN b)
{
  CharPtr              extension;
  Char                 path [PATH_MAX];
  TexT                 t;
  ValNodePtr           line_list;
  
  t = (TexT) GetObjectExtra (b);
  if (t == NULL)
  {
    return;
  }
  
  /* get filename from user */
  extension = GetAppProperty ("FastaNucExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  
  line_list = ReadLinesOfFile (path);
  AddLineListToText (line_list, t);
  ValNodeFreeData (line_list);
}

static void CheckSequenceAssistantCharInput (TexT t)
{
  SequenceAssistantPtr sap;
  CharPtr              seq_str;
  CharPtr              found_bracket, found_next_bracket;
  CharPtr              found_ret;
  Int4                 num_seq = 0;
  MsgAnswer            ans;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (t);
  if (sap == NULL) return;
  
  seq_str = SaveStringFromText (t);
  found_bracket = StringChr (seq_str, '>');
  if (found_bracket == NULL || *(found_bracket + 1) == 0
      || found_bracket != seq_str)
  {
    MemFree (seq_str);
    return;
  }
  found_ret = StringChr (found_bracket, '\n');
  if (found_ret == NULL)
  {
    MemFree (seq_str);
    return;
  }
  if (*(found_ret + 1) == 0)
  {
    MemFree (seq_str);
    return;
  }
  
  while (found_bracket != NULL)
  {
    found_ret = StringChr (found_bracket, '\n');
    if (found_ret == NULL)
    {
      /* last line is a defline */
      found_next_bracket = NULL;
    }
    else
    {
      found_next_bracket = StringChr (found_ret, '>');
    }
    num_seq++;
    found_bracket = found_next_bracket;
  }
  ans = Message (MSG_YN, "You are pasting in %d sequences, correct?", num_seq);
  if (ans == ANS_YES)
  {
    SetTitle (t, StringSave (""));

    ImportSequenceAssistantEditData (sap, seq_str);
    UpdateSequenceAssistant (sap);
  }
  MemFree (seq_str); 
}

static CharPtr GetSequenceString (SeqEntryPtr sep)
{
  CharPtr     seqbuf;
  BioseqPtr   bsp;
  SeqPortPtr  spp;
  Int2        ctr;
  Int4        read_len;
  
  if (sep == NULL || ! IS_Bioseq (sep))
  {
    return NULL;
  }
  
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || bsp->length < 1) 
  {
    return NULL;
  }

  spp = SeqPortNew (bsp, 0, bsp->length - 1, Seq_strand_plus, Seq_code_iupacna);
  seqbuf = (CharPtr) MemNew ((bsp->length + 1) * sizeof (Char));
  if (seqbuf != NULL)
  {
    SeqPortSeek (spp, 0, SEEK_SET);
    read_len = 0;
    while (read_len < bsp->length)
    {
      ctr = SeqPortRead (spp, (UcharPtr)(seqbuf + read_len), INT2_MAX);
      seqbuf[ctr + read_len] = 0;
      read_len += INT2_MAX;
    }
  }
  spp = SeqPortFree (spp);
  
  return seqbuf;
}

static CharPtr ReformatSequenceText (CharPtr seq_text)
{
  CharPtr src, dst;
  CharPtr new_text;
  Int4    num_lines;
  Int4    len;
  Int4    line_len = 80;
  Int4    counter;

  if (StringHasNoText (seq_text))
  {
  	MemFree (seq_text);
  	return NULL;
  }
  len = StringLen (seq_text);
  num_lines = len / line_len;
  len += num_lines + 2;
  new_text = (CharPtr) MemNew (len * sizeof (Char));
  if (new_text == NULL)
  {
  	return seq_text;
  }
  dst = new_text;
  counter = 0;
  for (src = seq_text; *src != 0; src++)
  {
  	if (!isspace ((Int4)(*src)))
  	{
  	  *dst = *src;
  	  dst++;
  	  counter++;
  	  if (counter == line_len)
  	  {
  	  	*dst = '\n';
  	  	dst++;
  	  	counter = 0;
  	  }
  	}
  }
  *dst = 0;
  MemFree (seq_text);
  return new_text;
}

static void FixStringForByteStore (CharPtr seq_str)
{
  CharPtr cp_src, cp_dst;
  
  if (seq_str == NULL)
  {
    return;
  }
  
  cp_src = seq_str;
  cp_dst = seq_str;
  while (*cp_src != 0)
  {
    if (isalpha (*cp_src))
    {
      *cp_dst = TO_UPPER (*cp_src);
      cp_dst++;
    }
    cp_src++;
  }
  *cp_dst = 0;
}

static Boolean 
IsDuplicateID (SeqEntryPtr seq_list, BioseqPtr edit_bsp, SeqIdPtr sip)
{
  SeqEntryPtr  sep;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqIdPtr     tmp_sip;
  Boolean      is_dup = FALSE;
  
  if (seq_list == NULL || sip == NULL)
  {
    return FALSE;
  }
  
  for (sep = seq_list; sep != NULL && !is_dup; sep = sep->next)
  {
    if (sep->data.ptrvalue == NULL)
    {
      continue;
    }
    if (IS_Bioseq (sep)) 
    {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != edit_bsp)
      {
        for (tmp_sip = sip; tmp_sip != NULL; tmp_sip = tmp_sip->next)
        {
          if (SeqIdIn (tmp_sip, bsp->id))
          {
            is_dup = TRUE;
          }  
        }
      }
    } 
    else if (IS_Bioseq_set (sep)) 
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      is_dup |= IsDuplicateID (bssp->seq_set, edit_bsp, sip);
    }
  }
  return is_dup;
}

static void PasteSequenceAssistant (IteM i)
{
  TexT    txt;
  CharPtr str;
  
  txt = (TexT) GetObjectExtra (i);
  if (txt == NULL)
  {
    return;
  }
  
  str = ClipboardToString ();
  SetTitle (txt, str);
  str = MemFree (str);
}

static void SequenceAssistantAddSequence (SequenceAssistantPtr sap)
{
  ModalAcceptCancelData acd;
  WindoW                w;
  GrouP                 h, sequence_grp, c;
  ButtoN                import_fasta_btn;
  TexT                  sequence_txt;
  Char                  str [200];
  ButtoN                b;
  SeqEntryPtr           new_sep;
  Boolean               done = FALSE;
  MenU                  edit_menu;
  IteM                  local_item;
  
  if (sap == NULL)
  {
    return;
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, "Add New Sequence", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  import_fasta_btn = PushButton (h, "Import Nucleotide FASTA", AddSequenceImportFasta);
  
  sequence_grp = NormalGroup (h, 1, 0, "Sequence Characters", programFont, NULL);
  SetGroupSpacing (sequence_grp, 10, 10);
  StaticPrompt (sequence_grp, "Paste or type the nucleotide sequence.", 0, 
                popupMenuHeight, programFont, 'l');
  sequence_txt = ScrollText (sequence_grp, 60, 10, programFont, FALSE, NULL);
  SetObjectExtra (sequence_txt, sap, NULL);
  
  SetObjectExtra (import_fasta_btn, sequence_txt, NULL);
  
  sprintf (str, "You may only use the valid IUPAC characters (%s).", 
                valid_iupac_characters);
  MultiLinePrompt (sequence_grp, str, 60 * stdCharWidth, programFont);
    
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton(c, "Save", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton(c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) import_fasta_btn,
                              (HANDLE) sequence_grp,
                              (HANDLE) c, 
                              NULL);

  /* Edit Menu */
  edit_menu = PulldownMenu (w, "Edit");
  local_item = CommandItem (edit_menu, "Paste", PasteSequenceAssistant);
  SetObjectExtra (local_item, sequence_txt, NULL);
                              
  Show(w); 
  Select (w);

  while (!done)
  {
    acd.accepted = FALSE;
    acd.cancelled = FALSE;  
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
  
    if (acd.accepted)
    {  
      new_sep = GetSequencesFromText (sequence_txt, sap->seq_list);
      if (new_sep != NULL)
      {
        ValNodeLink (&(sap->seq_list), new_sep);
        UpdateSequenceAssistant (sap);
        done = TRUE;
      }
    }
    else if (acd.cancelled)
    {
      done = TRUE;
    }
    acd.accepted = FALSE;
    acd.cancelled = TRUE;
  }
  
  Remove (w);
}

static void SequenceAssistantEditSequence (SequenceAssistantPtr sap, Int4 seq_num)
{
  ModalAcceptCancelData acd;
  WindoW                w;
  GrouP                 h, g1 = NULL, g2, sequence_grp, err_grp, c;
  SeqEntryPtr           sep;
  Int4                  i;
  Char                  window_title [150];
  BioseqPtr             bsp = NULL;
  SeqIdPtr              sip;
  Char                  id_txt [128];
  TexT                  sequence_id_txt, sequence_txt;
  CharPtr               ttl = NULL;
  Char                  str [200];
  ButtoN                b;
  CharPtr               seqbuf;
  SeqEntryPtr           last_sep = NULL, new_sep;
  CharPtr               id_str, title_str, seq_str;
  ValNodePtr            err_list = NULL;
  Int2                  seq_num_for_error = seq_num;
  Boolean               done = FALSE;
  ByteStorePtr          new_bs = NULL;
  SeqIdPtr              new_sip = NULL;
  SeqDescrPtr           sdp;
  MenU                  edit_menu;
  IteM                  local_item;
  TitleEditData         ted;
  
  if (sap == NULL)
  {
    return;
  }
  
  for (i = 0, sep = sap->seq_list; i != seq_num && sep != NULL; i++, sep = sep->next)
  {
  }
  
  if (sep != NULL)
  {
    if (sep->data.ptrvalue == NULL)
    {
      return;
    }
    else if (IS_Bioseq_set (sep))
    {
      Message (MSG_ERROR, "Can't edit segmented set!");
      return;
    }
    else if (! IS_Bioseq (sep))
    {
      return;
    }
    bsp = sep->data.ptrvalue;
    if (bsp->repr == Seq_repr_delta)
    {
      Message (MSG_ERROR, "Can't edit gapped sequence!");
      return;
    }
    sip = SeqIdFindWorst (bsp->id);
    SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);    
    sprintf (window_title, "Edit %s", id_txt);
  }
  else
  {
    sprintf (window_title, "Add new sequence");
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, window_title, NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  /* set up IDAndTitleEdit to use for single sequence */
  ted.iatep = IDAndTitleEditNew ();
  if (ted.iatep == NULL)
  {
    return;
  }
  
  ted.iatep->num_sequences = 1;
  ted.iatep->id_list = (CharPtr PNTR) MemNew (ted.iatep->num_sequences * sizeof (CharPtr));
  ted.iatep->title_list = (CharPtr PNTR) MemNew (ted.iatep->num_sequences * sizeof (CharPtr));
  ted.iatep->id_list[0] = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);
  
  ttl = NULL;
  SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
  if (StringHasNoText (ttl))
  {
    ted.iatep->title_list [0] = StringSave ("");
  }
  else
  {
    ted.iatep->title_list [0] = StringSave (ttl);
  }
  
  err_grp = HiddenGroup (h, 0, 0, NULL);
  ted.bracket_dlg = ShowDifferenceDialog (err_grp, stdCharWidth * 63, stdLineHeight * 5);
  ted.badvalue_pnl = InvalidValuesPanel (err_grp, stdCharWidth * 63, stdLineHeight * 5, 
                                         ted.iatep, NULL, TRUE,
                                         UpdateTitleEditForColorizedPanel,
                                         &ted,
                                         NULL, NULL);
  ted.unrec_mod_pnl = UnrecognizedModifiersPanel (err_grp, stdCharWidth * 63, stdLineHeight * 5, 
                                                  ted.iatep, NULL, TRUE,
                                                  UpdateTitleEditForColorizedPanel,
                                                  &ted,
                                                  NULL, NULL);
  
  /* users can enter titles and IDs for individual sequences */
  g2 = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g2, "Sequence ID", 0, popupMenuHeight, programFont, 'l');
  sequence_id_txt = DialogText (g2, "", 20, NULL);
  if (bsp != NULL)
  {
    SetTitle (sequence_id_txt, id_txt);
  }
  StaticPrompt (g2, "Sequence Title", 0, popupMenuHeight, programFont, 'l');
  ted.title_txt = DialogText (g2, "", 20, OnTitleEditChange);
  SetObjectExtra (ted.title_txt, &ted, NULL);
  if (!StringHasNoText (ttl))
  {
    SetTitle (ted.title_txt, ttl);
  }
  
  sequence_grp = NormalGroup (h, 1, 0, "Sequence Characters", programFont, NULL);
  SetGroupSpacing (sequence_grp, 10, 10);
  StaticPrompt (sequence_grp, "Paste or type the nucleotide sequence.", 0, 
                popupMenuHeight, programFont, 'l');
  sequence_txt = ScrollText (sequence_grp, 60, 10, programFont, FALSE, CheckSequenceAssistantCharInput);
  SetObjectExtra (sequence_txt, sap, NULL);
  seqbuf = GetSequenceString (sep);
  seqbuf = ReformatSequenceText (seqbuf);

  SetTitle (sequence_txt, seqbuf);
  MemFree (seqbuf);
  
  sprintf (str, "You may only use the valid IUPAC characters (%s).", 
                valid_iupac_characters);
  MultiLinePrompt (sequence_grp, str, 60 * stdCharWidth, programFont);
    
  c = HiddenGroup (h, 2, 0, NULL);
  ted.accept_btn = PushButton(c, "Save", ModalAcceptButton);
  SetObjectExtra (ted.accept_btn, &acd, NULL);
  b = PushButton(c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) err_grp, (HANDLE) g2, (HANDLE) sequence_grp,
                                (HANDLE) c, (HANDLE) g1, NULL);
  
  /* Edit Menu */
  edit_menu = PulldownMenu (w, "Edit");
  local_item = CommandItem (edit_menu, "Paste", PasteSequenceAssistant);
  SetObjectExtra (local_item, sequence_txt, NULL);
                
  ShowTitleEditErrors (&ted);                
                              
  Show(w); 
  Select (w);

  while (!done)
  {
    acd.accepted = FALSE;
    acd.cancelled = FALSE;  
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
  
    if (acd.accepted)
    {  
      id_str = SaveStringFromText (sequence_id_txt);
      id_str = ReformatLocalId (id_str);
      title_str = SaveStringFromText (ted.title_txt);
      seq_str = SaveStringFromText (sequence_txt);
      if (seq_num_for_error < 0)
      {
        seq_num_for_error = 0;
      }
  
      if (StringHasNoText (seq_str))
      {
        Message (MSG_ERROR, "You must supply sequence characters!");
      }
      else if (StringHasNoText (id_str))
      {
        Message (MSG_ERROR, "You must supply a sequence ID!");
      }
      else
      {
        done = TRUE;
        if (!SeqCharsOk (seq_str, seq_num_for_error, id_str, &err_list))
        {
          if (!ContinueWithErrorList (err_list, TRUE))
          {
            done = FALSE;
          }
        }
        
        
        if (done)
        { 
          new_sip = MakeSeqID (id_str);
          if (IsDuplicateID (sap->seq_list, bsp, new_sip))
          {
            Message (MSG_ERROR,
                   "Sequence IDs must be unique within the record!  %s is a duplicate ID",
                   id_str);
            done = FALSE;
            new_sip = SeqIdFree (new_sip);
          }
        }
        
        if (done)
        {
          FixStringForByteStore (seq_str);
          new_bs = BSNew (1000);
          if (new_bs != NULL)
          {
            BSWrite (new_bs, (VoidPtr) seq_str, (Int4) StringLen (seq_str));
          }
          
          if (bsp == NULL)
          {
            /* create new Bioseq and add to list */
            bsp = BioseqNew ();
            
            new_sep = SeqEntryNew ();
            new_sep->choice = 1;
            new_sep->data.ptrvalue = bsp;
            for (sep = sap->seq_list;
                 sep != NULL;
                 sep = sep->next)
            {
              last_sep = sep;
            }
            if (last_sep == NULL)
            {
              sap->seq_list = new_sep;
            }
            else
            {
              last_sep->next = new_sep;
            }
          }

          /* replace ID */
          bsp->id = SeqIdFree (bsp->id);
          bsp->id = new_sip;
            
          /* replace title */
          for (sdp = bsp->descr;
               sdp != NULL && sdp->choice != Seq_descr_title; 
               sdp = sdp->next)
          {
          }
          if (sdp == NULL)
          {
            sdp = SeqDescrNew (bsp->descr);
            if (bsp->descr == NULL)
            {
              bsp->descr = sdp;
            }
            sdp->choice = Seq_descr_title;
          }
          sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
          sdp->data.ptrvalue = StringSave (title_str);
          
          /* replace sequence data */
          bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
          bsp->repr = Seq_repr_raw;
          bsp->mol = Seq_mol_na;
          bsp->seq_data_type = Seq_code_iupacna;            
          bsp->seq_data = (SeqDataPtr) new_bs;
          bsp->length = BSLen (new_bs);
          
          BioseqPack (bsp);

          UpdateSequenceAssistant (sap);          
        }
      }
      id_str = MemFree (id_str);
      title_str = MemFree (title_str);
      seq_str = MemFree (seq_str);
      err_list = ValNodeFreeData (err_list);
    }
    else if (acd.cancelled)
    {
      done = TRUE;
    }
    acd.accepted = FALSE;
    acd.cancelled = TRUE;
  }
  
  ted.iatep = IDAndTitleEditFree (ted.iatep);
  Remove (w);
}


static void SequenceAssistantOk (SequenceAssistantPtr sap)
{
  if (sap == NULL) return;
  sap->cancelled = FALSE;
  sap->done = TRUE;
}

static void SequenceAssistantOkButton (ButtoN b)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (b);
  SequenceAssistantOk (sap);
}

static void SequenceAssistantOkItem (IteM i)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (i);
  SequenceAssistantOk (sap);
}


static void SequenceAssistantCancel (SequenceAssistantPtr sap)
{  
  if (sap == NULL) return;
  
  if (Message (MSG_YN, 
      "Are you sure you want to cancel (and lose all your editing)?") != ANS_YES)
  {
    return;
  }
  sap->cancelled = TRUE;
  sap->done = TRUE;
}

static void SequenceAssistantCancelButton (ButtoN b)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (b);
  SequenceAssistantCancel (sap);
}

static void SequenceAssistantCancelItem (IteM i)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (i);
  SequenceAssistantCancel (sap);
}

static void DeleteAllSequences (SequenceAssistantPtr sap)
{
  SeqEntryPtr sep, next_sep;
  if (sap == NULL || sap->seq_list == NULL)
  {
    return;
  }
  if (Message (MSG_YN, "Are you sure you want to delete all of your sequences?") != ANS_YES)
  {
    return;
  }
  
  sep = sap->seq_list;
  while (sep != NULL) 
  {
    next_sep = sep->next;
    sep->next = NULL;
    SeqEntryFree (sep);
    sep = next_sep;
  }
  sap->seq_list = NULL;
  UpdateSequenceAssistant (sap);
}

static void DeleteAllSequencesButton (ButtoN b)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (b);
  DeleteAllSequences (sap);
}

static void DeleteAllSequencesItem (IteM i)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (i);
  DeleteAllSequences (sap);
}

static void SelectSequenceDoc (DoC d, PoinT pt)
{
  Int2      item, row, prevrow;
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (d);
  if (sap == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, NULL, NULL);
  if (item > 0 && row > 0) {
    prevrow = sap->sequence_row;
    sap->sequence_row = item;
    if (item != prevrow)
    {
      if (prevrow != -1)
      {
        InvalDocRows (d, prevrow, 1, 1);
      }
      InvalDocRows (d, item, 1, 1);
    }  	
    Enable (sap->edit_btn);
    Enable (sap->delete_btn);
  }
}

static Boolean SequenceHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (doc);
  if (sap == NULL) return FALSE;
  
  if (item == sap->sequence_row) return TRUE;
  return FALSE;
}

static void AddSequence (SequenceAssistantPtr sap)
{
  if (sap == NULL) return;
  SequenceAssistantAddSequence (sap);
}

static void AddSequenceButton (ButtoN b)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (b);

  AddSequence (sap);
}

static void AddSequenceItem (IteM i)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (i);
  AddSequence (sap);  
}

static Boolean ConfirmSequenceDelete (SeqEntryPtr sep)
{
  MsgAnswer            ans;
  BioseqPtr            bsp;
  SeqIdPtr             sip;
  Char                 tmp[128];
  SeqDescrPtr          sdp;
  
  if (sep == NULL)
  {
    return FALSE;
  }
  
  bsp = (BioseqPtr) sep->data.ptrvalue;
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  sip = SeqIdFindWorst (bsp->id);
  SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
  if (sdp != NULL && ! StringHasNoText (sdp->data.ptrvalue))
  {
    ans = Message (MSG_YN, "Are you sure you want to delete %s (%s)?", tmp, sdp->data.ptrvalue);
  }
  else
  {
    ans = Message (MSG_YN, "Are you sure you want to delete %s?", tmp);
  }
  if (ans == ANS_YES)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

static void DeleteSequence (SequenceAssistantPtr sap)
{
  Int2                 seq_num;
  SeqEntryPtr          sep, last_sep = NULL;
  
  if (sap == NULL || sap->seq_list == NULL) return;
  
  if (sap->seq_list != NULL)
  {
    for (sep = sap->seq_list, seq_num = 1;
         sep != NULL && seq_num < sap->sequence_row;
         sep = sep->next, seq_num++)
    {
      last_sep = sep;
    }
    if (!ConfirmSequenceDelete (sep))
    {
      return;
    }
    
    if (sep == NULL)
    {
      /* do nothing, deleted non-existent sequence */
    }
    else if (last_sep == NULL)
    {
      /* remove first in list */
      sap->seq_list = sap->seq_list->next;
      sep->next = NULL;
      sep = SeqEntryFree (sep);
    }
    else
    {
      last_sep->next = sep->next;
      sep->next = NULL;
      sep = SeqEntryFree (sep);
    }
    if (sap->sequence_row > 1)
    {
      sap->sequence_row --;
    }
  }
  UpdateSequenceAssistant (sap);
}

static void DeleteSequenceButton (ButtoN b)
{
  SequenceAssistantPtr sap;

  sap = (SequenceAssistantPtr) GetObjectExtra (b);
  DeleteSequence (sap);  
}

static void DeleteSequenceItem (IteM i)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (i);
  DeleteSequence (sap);  
}

static void ImportFastaFileItem (IteM i)
{
  SequenceAssistantPtr sap;
  CharPtr              extension;
  Char                 path [PATH_MAX];
  SeqEntryPtr          new_sep_list;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (i);
  if (sap == NULL)
  {
    return;
  }

  extension = GetAppProperty ("FastaNucExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  
  new_sep_list = GetSequencesFromFile (path, sap->seq_list);
  if (new_sep_list != NULL
      && ImportedSequenceTypeOk (new_sep_list, sap->seqPackage))
  {
    ValNodeLink (&sap->seq_list, new_sep_list);
  }
  else
  {
    new_sep_list = SeqEntryFree (new_sep_list);
  }

  UpdateSequenceAssistant (sap);
}

static void ImportFastaFileButton (ButtoN b)
{
  SequenceAssistantPtr sap;
  CharPtr              extension;
  Char                 path [PATH_MAX];
  SeqEntryPtr          new_sep_list;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL)
  {
    return;
  }

  extension = GetAppProperty ("FastaNucExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;

  new_sep_list = GetSequencesFromFile (path, sap->seq_list);
  if (new_sep_list != NULL
      && ImportedSequenceTypeOk (new_sep_list, sap->seqPackage))
  {
    ValNodeLink (&sap->seq_list, new_sep_list);
  }
  else
  {
    new_sep_list = SeqEntryFree (new_sep_list);
  }

  UpdateSequenceAssistant (sap);
}

static Boolean 
SequenceAssistantValidateSegments 
(SeqEntryPtr seq_list,
 ValNodePtr PNTR err_list)
{
  SeqEntryPtr sep;
  Boolean     all_one_char = TRUE;
  CharPtr     seqbuf;
  Char        first_seg_char = 0;
  Char        first_seg_char_this = 0;
  Boolean     rval = TRUE;
  Int4        total_len = 0;
  Boolean     non_N = FALSE;
  
  if (seq_list == NULL)
  {
    return TRUE;
  }
  
  for (sep = seq_list; sep != NULL; sep = sep->next)
  {
    seqbuf = GetSequenceString (sep);
    if (!StringHasNoText (seqbuf))
    {
      if (all_one_char)
      {
        if (IsSequenceAllOneCharacter(seqbuf))
        {
          if (first_seg_char == 0)
          {
            first_seg_char = seqbuf[StringSpn(seqbuf, " \t\n")];
          }
          else
          {
            first_seg_char_this = seqbuf[StringSpn(seqbuf, " \t\n")];
            if (first_seg_char_this != first_seg_char)
            {
              all_one_char = FALSE;
            }
          }
        }
        else
        {
          all_one_char = FALSE;
        }
      }
      
      total_len += CountSeqChars (seqbuf);
      non_N |= ! IsSequenceAllNs (seqbuf);
    }
    seqbuf = MemFree (seqbuf);
  }
  if (all_one_char && err_list != NULL)
  {
    ValNodeAddPointer (err_list, CREATE_FASTA_WARNING,
             StringSave ("Your segmented set sequences all consist entirely of the same character."));
  }
  if (!non_N)
  {
    if (err_list != NULL)
    {
      ValNodeAddPointer (err_list, CREATE_FASTA_REQUIRED, 
                         StringSave ("Your segmented set consists entirely of Ns. "
                         "This is not a valid sequence.  Please edit."));

    }
    rval = FALSE;
  }
  if (total_len < 50)
  {
    if (err_list != NULL)
    {
      /* Note - this is a required error because only small molecules can have
       * less than 50 base pairs, and the small molecules should never have
       * been sequenced in segments.
       */
      ValNodeAddPointer (err_list, CREATE_FASTA_REQUIRED, StringSave (
                     "You have fewer than 50 total base pairs in this "
                     "segmented set. GenBank will not accept segmented sets with "
                     "fewer than 50 base pairs. Please edit your sequence."));
    }
    rval = FALSE;
  }
  return rval;
}

static Boolean 
SequenceAssistantValidateOneBioseqContentAndLength 
(BioseqPtr       bsp,
 ValNodePtr PNTR all_N_list,
 ValNodePtr PNTR all_one_char_list,
 ValNodePtr PNTR too_short_list)
{
  CharPtr     seqbuf;
  Boolean     rval = TRUE;
  Int2        seq_num = 0;
  CharPtr     id_str;
  SeqEntryPtr sep;
  
  if (bsp == NULL || too_short_list == NULL || all_N_list == NULL || all_one_char_list == NULL)
  {
    return FALSE;
  }
  
  sep = SeqMgrGetSeqEntryForData (bsp);
  
  seqbuf = GetSequenceString (sep);
  
  id_str = SeqIdWholeLabel (SeqIdFindWorst (bsp->id), PRINTID_REPORT);

  if (CountSeqChars (seqbuf) < 50)
  {
    ValNodeAddPointer (too_short_list, seq_num, StringSave (id_str));
    rval = FALSE;
  }
  if (IsSequenceAllNs (seqbuf))
  {
    ValNodeAddPointer (all_N_list, seq_num, id_str);
    rval = FALSE;
  }
  else if (IsSequenceAllOneCharacter(seqbuf))
  {
    ValNodeAddPointer (all_one_char_list, seq_num, id_str);
    rval = FALSE;
  }
  else 
  {
    id_str = MemFree (id_str);
  }
  seqbuf = MemFree (seqbuf);
  return rval;
}

static Boolean 
AddBioseqErrors 
(ValNodePtr all_N_list,
 ValNodePtr all_one_char_list,
 ValNodePtr too_short_list,
 ValNodePtr PNTR err_list)
{
  CharPtr err_msg = NULL;
  Boolean rval = TRUE;
  
  if (all_N_list != NULL)
  {
    if (err_list != NULL)
    {
      if (all_N_list->next == NULL)
      {
        err_msg = CreateListMessage ("In sequence", 
                     " there are only Ns. This not a valid sequence, please edit it.",
                     all_N_list);
      }
      else
      {
        err_msg = CreateListMessage ("In sequence", 
                   " there are only Ns. These are not valid sequences, please edit them.",
                     all_N_list);
      }
      ValNodeAddPointer (err_list, CREATE_FASTA_REQUIRED, err_msg);
    }
    rval = FALSE;
  }
  if (all_one_char_list != NULL)
  {
    if (err_list != NULL)
    {
      if (all_one_char_list->next == NULL)
      {
        err_msg = CreateListMessage ("In sequence", 
                   " one character is repeated for the entire sequence. This is not a valid sequence, please edit it.",
                   all_one_char_list);
      }
      else
      {
        err_msg = CreateListMessage ("In sequence", 
                     " one character is repeated for the entire sequence. These are not valid sequences, please edit them.",
                     all_one_char_list);
      }
      ValNodeAddPointer (err_list, CREATE_FASTA_WARNING, err_msg);
    }
    rval = FALSE;
  }
  if (too_short_list != NULL && err_list != NULL)
  {
    if (too_short_list->next == NULL)
    {
      err_msg = CreateListMessage ("Sequence", " is shorter than 50 base pairs. "
                     "GenBank will not accept sequences with "
                     "fewer than 50 base pairs. Please edit your sequence or "
                     "make sure that your comment explains why your sequence "
                     "is so short.",
                     too_short_list);
    }
    else
    {
      err_msg = CreateListMessage ("Sequence", " are shorter than 50 base pairs. "
                     "GenBank will not accept sequences with "
                     "fewer than 50 base pairs. Please edit your sequences or "
                     "make sure that your comments explain why your sequences "
                     "are so short.",
                     too_short_list);
    }
    ValNodeAddPointer (err_list, CREATE_FASTA_WARNING, err_msg);
  }
  return rval;
}

/* This function will add to a list of errors.  Any errors in the list with a choice of 0
 * cause the sequence to be unacceptable.  Any errors in the list with a choice of 1 are
 * a yes-no question - yes means that the user wants to go back and correct the problems.
 * no means the user would like to continue anyway.
 */
static Boolean 
SequenceAssistantValidateContentAndLength 
(SeqEntryPtr     seq_list,
 ValNodePtr PNTR all_N_list,
 ValNodePtr PNTR all_one_char_list,
 ValNodePtr PNTR too_short_list,
 ValNodePtr PNTR err_list)
{
  Boolean           rval = TRUE;
  BioseqSetPtr      bssp;
  
  if (seq_list == NULL)
  {
    return TRUE;
  }
  
  if (IS_Bioseq_set (seq_list))
  {
    bssp = (BioseqSetPtr) seq_list->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_parts)
    {
/*      rval = SequenceAssistantValidateSegments (bssp->seq_set, err_list); */
    }
    else
    {
      rval = SequenceAssistantValidateContentAndLength (bssp->seq_set, 
                                                        all_N_list,
                                                        all_one_char_list,
                                                        too_short_list,
                                                        err_list);
    }
  }
  else if (IS_Bioseq (seq_list))
  {
    rval = SequenceAssistantValidateOneBioseqContentAndLength (seq_list->data.ptrvalue,
                                                               all_N_list,
                                                               all_one_char_list,
                                                               too_short_list);
  }
  rval &= SequenceAssistantValidateContentAndLength (seq_list->next, 
                                                            all_N_list,
                                                            all_one_char_list,
                                                            too_short_list,
                                                           err_list);
  return rval;
}

static Boolean SequenceAssistantValidate (SeqEntryPtr seq_list)
{
  ValNodePtr err_list = NULL;
  Boolean    rval;
  ValNodePtr all_N_list = NULL;
  ValNodePtr all_one_char_list = NULL;
  ValNodePtr too_short_list = NULL;
  
  rval = SequenceAssistantValidateContentAndLength (seq_list, 
                                                    &all_N_list,
                                                    &all_one_char_list,
                                                    &too_short_list,
                                                    &err_list);
  
  rval &= AddBioseqErrors (all_N_list, all_one_char_list, too_short_list, &err_list);

  all_N_list = ValNodeFreeData (all_N_list);
  too_short_list = ValNodeFreeData (too_short_list);
  all_one_char_list = ValNodeFreeData (all_one_char_list);

  if (err_list != NULL)
  {
    rval = ContinueWithErrorList (err_list, TRUE);
  }
  err_list = ValNodeFreeData (err_list);
  return rval;
}

static void SequenceDblClick (PoinT cell_coord, CharPtr header_text, CharPtr cell_text, Pointer userdata)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) userdata;
  if (sap == NULL)
  {
    return;
  }
  if (cell_coord.x < 2)
  {
    return;
  }
  if (cell_coord.y == 0)
  {
    if (cell_coord.x == SEQUENCE_ASSISTANT_MOLECULE_COLUMN)
    {
      /* edit all molecule types */
      EditOrgModColumn ("moltype", NULL, sap->seq_list, sap->seqPackage);
    }
    else if (cell_coord.x == SEQUENCE_ASSISTANT_TOPOLOGY_COLUMN)
    {
      /* edit all topologies */
      EditOrgModColumn ("topology", NULL, sap->seq_list, sap->seqPackage);
    }
    else
    {
      /* edit all titles */
      EditSequenceTitleColumns (sap);
    }
    UpdateSequenceAssistant (sap);
  }
  else
  {
    if (cell_coord.x == SEQUENCE_ASSISTANT_MOLECULE_COLUMN)
    {
      /* edit one molecule type */
      ApplyOrgModColumnOrCell ("moltype", cell_text, cell_coord.y - 1, NULL, sap->seq_list,
                               NULL, 0, sap->seqPackage);
    }
    else if (cell_coord.x == SEQUENCE_ASSISTANT_TOPOLOGY_COLUMN)
    {
      /* edit one topology */
      ApplyOrgModColumnOrCell ("topology", cell_text, cell_coord.y - 1, NULL, sap->seq_list,
                               NULL, 0, sap->seqPackage);
    }
    else
    {
      /* edit one title */
      EditOneSequenceTitle (sap, cell_coord.y - 1);
    }
    UpdateSequenceAssistant (sap);
  }
}

static void SequenceAssistantEditButton (ButtoN b)
{
  SequenceAssistantPtr sap;
  
  sap = (SequenceAssistantPtr) GetObjectExtra (b);
  if (sap == NULL)
  {
    return;
  }
  
  SequenceAssistantEditSequence (sap, sap->sequence_row - 1);
}

static void SequenceAssistant (ButtoN b)
{
  SequencesFormPtr      sqfp;
  SequenceAssistantData sad;
  WindoW                w;
  MenU                  m;
  IteM                  i;
  GrouP                 h, k, edit_grp, selector_grp, c;
  FastaPagePtr          fpp;
  RecT                  r;
  SeqEntryPtr           sep, next_sep;
  ButtoN                add_btn;
  Int4                  doc_width;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL) 
  {
    return;
  }
  
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp == NULL) 
  {
    return;
  }
  
  sad.done = FALSE;
  sad.cancelled = FALSE;
  
  sad.seq_list = CopySeqEntryList (fpp->list);
  sad.sequence_row = -1;
  
  sad.seqPackage = sqfp->seqPackage;
    
  w = MovableModalWindow (-20, -13, -10, -10, "Specify Sequences", NULL);
  /* add menus */  
  m = PulldownMenu (w, "File");
  i = CommandItem (m, "Import FASTA file", ImportFastaFileItem);
  SetObjectExtra (i, &sad, NULL);
  i = CommandItem (m, "Done", SequenceAssistantOkItem);
  SetObjectExtra (i, &sad, NULL);
  i = CommandItem (m, "Cancel", SequenceAssistantCancelItem);
  SetObjectExtra (i, &sad, NULL);
  
  /* edit menu */
  m = PulldownMenu (w, "Edit");
  i = CommandItem (m, "Add Sequence", AddSequenceItem);
  SetObjectExtra (i, &sad, NULL);
  i = CommandItem (m, "Delete Sequence", DeleteSequenceItem);
  SetObjectExtra (i, &sad, NULL);
  i = CommandItem (m, "Delete All Sequences", DeleteAllSequencesItem);
  SetObjectExtra (i, &sad, NULL);
    
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  k = HiddenGroup (h, 2, 0, NULL);  
  sad.import_btn = PushButton (k, "Import Additional Nucleotide FASTA", ImportFastaFileButton);
  SetObjectExtra (sad.import_btn, &sad, NULL);
  
  add_btn = PushButton (k, "Add New Sequence", AddSequenceButton);
  SetObjectExtra (add_btn, &sad, NULL);
  
  edit_grp = HiddenGroup (h, 2, 0, NULL);
  sad.sequence_selector = DocumentPanel (edit_grp, stdCharWidth * 10, stdLineHeight * 6);
  SetObjectExtra (sad.sequence_selector, &sad, NULL);
  SetDocProcs (sad.sequence_selector, SelectSequenceDoc, NULL, NULL, NULL);
  SetDocShade (sad.sequence_selector, NULL, NULL, SequenceHighlight, NULL);
  selector_grp = HiddenGroup (edit_grp, 0, 4, NULL);
  sad.edit_btn = PushButton (selector_grp, "Edit Sequence", SequenceAssistantEditButton);
  SetObjectExtra (sad.edit_btn, &sad, NULL);
  StaticPrompt (selector_grp, "", 0, popupMenuHeight, programFont, 'l');
  sad.delete_btn = PushButton (selector_grp, "Delete Sequence", DeleteSequenceButton);
  SetObjectExtra (sad.delete_btn, &sad, NULL);
  sad.delete_all_btn = PushButton (selector_grp, "Delete All Sequences", DeleteAllSequencesButton);
  SetObjectExtra (sad.delete_all_btn, &sad, NULL);

  sad.summary_dlg = FastaSummaryDialog (h);

  doc_width = GetStandardTableDisplayDialogWidth (sqfp);
  
  sad.sequence_table = TableDisplayDialog (h, doc_width, stdLineHeight * 8, 1, 2,
                                       SequenceDblClick, &sad,
                                       NULL, NULL);
  
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton(c, "Done", SequenceAssistantOkButton);
  SetObjectExtra (b, &sad, NULL);
  b = PushButton(c, "Cancel", SequenceAssistantCancelButton);
  SetObjectExtra (b, &sad, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) k,
                              (HANDLE) edit_grp,
                              (HANDLE) sad.summary_dlg,
                              (HANDLE) sad.sequence_table,
                              (HANDLE) c, NULL);

  UpdateSequenceAssistant (&sad);
                                
  Show(w); 
  Select (w);
  
  while (!sad.done)
  {
    while (!sad.done)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
    if (!sad.cancelled)
    {

      if (SequenceAssistantValidate (sad.seq_list))
      {
        /* check for number of sequences */
        if (sad.seq_list != NULL && sad.seq_list->next != NULL
            && PackageTypeIsSingle (sqfp->seqPackage))
        {
          if (Message (MSG_YN, "You are importing multiple sequences - did you intend to create a batch submission?") == ANS_YES)
          {
            sqfp->seqPackage = SEQ_PKG_GENBANK;
            fpp->single = FALSE;
            SafeHide (fpp->singleIdGrp);
          }
          else
          {
            sad.done = FALSE;
            sad.cancelled = FALSE;
          }
        }
      }
      else
      {
        sad.done = FALSE;
        sad.cancelled = FALSE;
      }
      
      if (sad.done)
      {
        ResetFastaPage (fpp);
        Reset (fpp->doc);
        if (sad.seq_list == NULL)
        {          
          SafeHide (fpp->have_seq_instr_grp);
          SafeShow (fpp->instructions);
          Update ();
          Enable (fpp->import_btn);
          SetTitle (fpp->import_btn, "Import Nucleotide FASTA");
          Disable (fpp->clear_btn);
          ClearOrganismModifiers (sqfp);
          Disable (sqfp->molecule_btn);
          Disable (sqfp->topology_btn);
        }
        else
        {    
          /* these statements make sure the column width is large enough */
          ObjectRect (fpp->doc, &r);
          InsetRect (&r, 4, 4);
          faColFmt.pixWidth = r.right - r.left;
 
          fpp->list = sad.seq_list;
          SafeHide (fpp->instructions);
          Update ();
          
          if (PackageTypeIsSingle (sqfp->seqPackage) || sqfp->seqPackage == SEQ_PKG_GENOMICCDNA)
          {
            Disable (fpp->import_btn);
          }
          else
          {
            Enable (fpp->import_btn);
            SetTitle (fpp->import_btn, "Import Additional Nucleotide FASTA");
          }
          Enable (fpp->clear_btn);
          FormatFastaDoc (fpp);
          SafeShow (fpp->have_seq_instr_grp);
          NucleotideImportFinish (sqfp);
          if (sqfp->seqFormat == SEQ_FMT_FASTA
              && SuggestJumpingToWizard (fpp->list)) {
            fpp->list = NULL;
            Remove (sqfp->form);
          }
        }      
      }
      else
      {
        sad.cancelled = FALSE;
        sad.done = FALSE;
      }  
    }
    else
    {
      /* clean up list of sequences from form, since they will not be used */
      sep = sad.seq_list;
      while (sep != NULL) 
      {
        next_sep = sep->next;
        sep->next = NULL;
        SeqEntryFree (sep);
        sep = next_sep;
      }
      sad.seq_list = NULL;
    }
  }
  Remove (w);
}

static void SpecifyMolecule (ButtoN b)
{
  SpecifyModValueButton (b, "moltype");
}

static void SpecifyTopology (ButtoN b)
{
  SpecifyModValueButton (b, "topology");
}

static void CustomAlnSettings (ButtoN b)
{
  SequencesFormPtr sqfp;
  TSequenceInfoPtr new_settings = NULL;
  PhylipPagePtr    ppp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL) return;

  ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (ppp == NULL) return;

  new_settings = GetAlignmentOptions (NULL, ppp->aln_settings);
  if (new_settings != NULL) 
  {
    SequenceInfoFree (ppp->aln_settings);
    ppp->aln_settings = new_settings;
  }
}

static GrouP CreateNucleotideTab (GrouP h, SequencesFormPtr sqfp)
{
  GrouP              q, g, x, y, k;
  ButtoN             b = NULL;
  Handle             h1 = NULL, h2 = NULL;
  Boolean            single;
  GrouP              import_btn_grp = NULL;
  FastaPagePtr       fpp;
  
  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 10);
  g = HiddenGroup (q, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);

  y = HiddenGroup (g, -2, 0, NULL);
  SetGroupSpacing (y, 10, 2);

  if (PackageTypeIsSet (sqfp->seqPackage) 
      && sqfp->seqPackage != SEQ_PKG_GENBANK /* exclude batch submissions */
      && sqfp->seqPackage != SEQ_PKG_TSA /* exclude TSA */
      && sqfp->seqFormat == SEQ_FMT_FASTA) 
	{
        sqfp->makeAlign = CheckBox (g, "Create Alignment", NULL);
      /*
      if (sqfp->seqPackage < SEQ_PKG_GENBANK) {
        SetStatus (sqfp->makeAlign, TRUE);
      }
      */
  }

  k = HiddenGroup (g, 0, 2, NULL);
  if (sqfp->seqFormat == SEQ_FMT_FASTA) {
    single = PackageTypeIsSingle (sqfp->seqPackage);
    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      sqfp->dnaseq = CreateFastaDialog (k, "", TRUE, FALSE, fastaGenMsg, single, NULL);
      fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
      import_btn_grp = HiddenGroup (g, 4, 0, NULL);
      fpp->import_btn = PushButton (import_btn_grp, "Import Genomic FASTA", ImportBtnProc);
    } else {
      sqfp->dnaseq = CreateFastaDialog (k, "", TRUE, FALSE, fastaNucMsg, single, &(sqfp->seqPackage));
      fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
      import_btn_grp = HiddenGroup (g, 4, 0, NULL);
      fpp->import_btn = PushButton (import_btn_grp, "Import Additional Nucleotide FASTA", ImportBtnProc);
      SetTitle (fpp->import_btn, "Import Nucleotide FASTA");
      if (sqfp->seqPackage == SEQ_PKG_GAPPED) 
      {
        fpp->is_delta = TRUE;
      }
      else
      {
        fpp->is_delta = FALSE;
      }
    }
    SetObjectExtra (fpp->import_btn, sqfp, NULL);
    
    if (sqfp->seqPackage != SEQ_PKG_GAPPED 
        && sqfp->seqPackage != SEQ_PKG_GENOMICCDNA)
    {
      b = PushButton (import_btn_grp, "Add/Modify Sequences", SequenceAssistant);
      SetObjectExtra (b, sqfp, NULL);
    }
    fpp->clear_btn = PushButton (import_btn_grp, "Clear Sequences", ClearSequencesButton);
    SetObjectExtra (fpp->clear_btn, sqfp, NULL);
    Disable (fpp->clear_btn);
  } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
    sqfp->dnaseq = CreatePhylipDialog (k, "", phylipNucMsg, sqfp->seqFormat, "",
                                       sqfp->seqPackage);
    import_btn_grp = HiddenGroup (g, 4, 0, NULL);
    b = PushButton (import_btn_grp, "Import Nucleotide Alignment", ImportBtnProc);
    SetObjectExtra (b, sqfp, NULL);
    b = PushButton (import_btn_grp, "Optional Alignment Settings", CustomAlnSettings);
    SetObjectExtra (b, sqfp, NULL);
  }
  
  x = HiddenGroup (g, -4, 0, NULL);  
  
  sqfp->molecule_btn = PushButton (x, "Specify Molecule", SpecifyMolecule);
  SetObjectExtra (sqfp->molecule_btn, sqfp, NULL);
  Disable (sqfp->molecule_btn);
  sqfp->topology_btn = PushButton (x, "Specify Topology", SpecifyTopology);
  SetObjectExtra (sqfp->topology_btn, sqfp, NULL);
  Disable (sqfp->topology_btn);

  if (sqfp->makeAlign != NULL) {
    h1 = (Handle) sqfp->makeAlign;
    h2 = (Handle) import_btn_grp;
  } else {
    h1 = import_btn_grp;
    h2 = NULL;
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, (HANDLE) k,
                  (HANDLE) h1, (HANDLE) h2, NULL);
  return q;
}

static GrouP CreateSequencingMethodTab (GrouP h, SequencesFormPtr sqfp)
{
  sqfp->sequencing_method_dlg = SequencingMethodDialog (h);
  return (GrouP) sqfp->sequencing_method_dlg;
}

static GrouP CreateTranscriptsTab (GrouP h, SequencesFormPtr sqfp)
{
  GrouP   q, g, y, k;
  ButtoN  b;
  
  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 20);
  g = HiddenGroup (q, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  y = HiddenGroup (g, -2, 0, NULL);
  SetGroupSpacing (y, 10, 2);
  sqfp->partialmRNA5 = CheckBox (y, "Incomplete at 5' end", NULL);
  sqfp->partialmRNA3 = CheckBox (y, "Incomplete at 3' end", NULL);

  k = HiddenGroup (g, 0, 2, NULL);
  sqfp->mrnaseq = CreateFastaDialog (k, "", TRUE, TRUE, fastaMrnaMsg, FALSE, NULL);
  b = PushButton (g, "Import Transcript FASTA", ImportBtnProc);
  SetObjectExtra (b, sqfp, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) y, (HANDLE) k, (HANDLE) b, NULL);
  return q;
}


static GrouP CreateProteinTab (GrouP h, SequencesFormPtr sqfp)
{
  GrouP        q, g, x, y, k, prot_btns;
  ButtoN       mrna = NULL, b;
  Char         str [32];
  FastaPagePtr fpp;
  
  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 20);
  g = HiddenGroup (q, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  x = HiddenGroup (g, -1, 0, NULL);
  sqfp->protTechBoth = CheckBox (x,
             "Conceptual translation confirmed by peptide sequencing", NULL);
  y = HiddenGroup (g, -2, 0, NULL);
  SetGroupSpacing (y, 10, 2);
  sqfp->partialN = CheckBox (y, "Incomplete at NH2 end", NULL);
  sqfp->partialC = CheckBox (y, "Incomplete at CO2H end", NULL);

  sqfp->makeMRNA = FALSE;
  if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA) {
    mrna = CheckBox (g, "Create initial mRNA with CDS intervals", ChangeMrnaFlag);
    SetObjectExtra (mrna, sqfp, NULL);
    if (GetAppParam ("SEQUIN", "PREFERENCES", "CREATEMRNA", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        sqfp->makeMRNA = TRUE;
      }
    }
  }
  SafeSetStatus (mrna, sqfp->makeMRNA);
  k = HiddenGroup (g, 0, 2, NULL);
  sqfp->protseq = CreateFastaDialog (k, "", FALSE, FALSE, fastaProtMsg, FALSE, NULL);
  prot_btns = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (prot_btns, 10, 10);
  b = PushButton (prot_btns, "Import Protein FASTA", ImportBtnProc);
  SetObjectExtra (b, sqfp, NULL);
  
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
  if (fpp != NULL)
  {
    fpp->clear_btn = PushButton (prot_btns, "Clear Protein Sequences", ClearSequencesButton);
    SetObjectExtra (fpp->clear_btn, sqfp, NULL);
    Disable (fpp->clear_btn);
  }
  
  AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, (HANDLE) k,
                (HANDLE) prot_btns, (HANDLE) mrna, NULL);
    

  return q;
}


static GrouP CreateAnnotTab (GrouP h, SequencesFormPtr sqfp)
{
  GrouP  q, z;
  PrompT ppt1, ppt2;
  RecT   comment_rect, defline_rect;
  
  q = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (q, 10, 10);
  ppt1 = StaticPrompt (q, "Add feature across full length of all sequences",
                       0, 0, programFont, 'l');
  sqfp->feature_info = SubmissionFeatureInfoDialog (q);
  ppt2 = StaticPrompt (q, "Add title to all sequences if not in definition line",
                       0, 0, programFont, 'c');
  z = HiddenGroup (q, 2, 0, NULL);
  StaticPrompt (z, "Title       ", 0, 3 * Nlm_stdLineHeight, programFont, 'c');
  sqfp->defline = ScrollText (z, 20, 3, programFont, TRUE, NULL);
  sqfp->orgPrefix = CheckBox (q, "Prefix title with organism name", NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) sqfp->feature_info,
                (HANDLE) ppt2,
                (HANDLE) sqfp->orgPrefix, NULL);
  
  /* fix title scroll to be the same size as and aligned with the comment scroll */
  GetSubmissionFeatureInfoFeatCommentRect(sqfp->feature_info, &comment_rect);
  ObjectRect (sqfp->defline, &defline_rect);
  /* because the doc has a vertical scroll bar and the set position subtracts the
   * width of the scroll bar before positioning the list, must add the width of
   * the scroll bar to the rightt.
   */
  defline_rect.left = comment_rect.left;
  defline_rect.right = comment_rect.right + Nlm_vScrollBarWidth;

  SetPosition (sqfp->defline, &defline_rect);
  
  return q;
}


static void RemoveSequencesFromSubmission (IteM i)
{
  SequencesFormPtr sqfp;
  SeqEntryPtr  sep;

  sqfp = (SequencesFormPtr) GetObjectExtra (i);
  if (sqfp == NULL) {
    return;
  }

  sep = GetSeqEntryFromSequencesForm (sqfp);
  if (RemoveSequencesFromWizardList (&sep)) {
    /* redraw */
    SetSequencesForSubmissionForm ((WindoW)sqfp->form, sep, 0);
    Update();
  }
}


static void SetSequencesFormPage (SequencesFormPtr sqfp, Int4 page);
static void RemoveAllSequencesFromSubmission (IteM i)
{
  SequencesFormPtr sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (i);
  if (sqfp == NULL) {
    return;
  }
  /* delete the sequences and organism modifiers */
  DeleteAllSequencesFromForm (sqfp);
  /* put the user back on the sequences page */
  SetSequencesFormPage (sqfp, 0);
  Update();
}


static void CleanupSequencesForm (GraphiC g, Pointer data)
{
  SequencesFormPtr sqfp;
  
  if (data != NULL)
  {
    sqfp = (SequencesFormPtr) data;
    sqfp->nuc_prot_assoc_list = FreeAssociationList (sqfp->nuc_prot_assoc_list);
  }
  StdCleanupFormProc (g, data);
}

extern ForM CreateInitOrgNucProtForm (Int2 left, Int2 top, CharPtr title,
                                      FormatBlockPtr format,
                                      BtnActnProc goToNext,
                                      BtnActnProc goBack,
                                      WndActnProc activateForm)

{
  GrouP              c;
  GrouP              h;
  GrouP              j;
  GrouP              tbs_grp;
  Int2               page;
  StdEditorProcsPtr  sepp;
  SequencesFormPtr   sqfp;
  WindoW             w;
  MenU               edit_menu = NULL;
  IteM               i;

  w = NULL;
  sqfp = MemNew (sizeof (SequencesForm));
  if (sqfp != NULL) {

    if (format != NULL) {
      sqfp->seqPackage = format->seqPackage;
      sqfp->seqFormat = format->seqFormat;
      sqfp->numSeqs = format->numSeqs;
      sqfp->submType = format->submType;
    } else {
      sqfp->seqPackage = SEQ_PKG_SINGLE;
      sqfp->seqFormat = SEQ_FMT_FASTA;
      sqfp->numSeqs = 0;
      sqfp->submType = SEQ_ORIG_SUBMISSION;
    }
    sqfp->nuc_prot_assoc_list = NULL;

    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, sqfp, CleanupSequencesForm);
    sqfp->form = (ForM) w;
    sqfp->toform = NULL;
    if (sqfp->seqFormat == SEQ_FMT_FASTA) {
      sqfp->fromform = FastaSequencesFormToSeqEntryPtr;
    } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
      sqfp->fromform = PhylipSequencesFormToSeqEntryPtr;
    }
    sqfp->testform = NULL;
    sqfp->importform = ImportSequencesForm;
    sqfp->exportform = ExportSequencesForm;
    sqfp->formmessage = SequencesFormMessage;

#ifndef WIN_MAC
    edit_menu = CreateSqnInitialFormMenus (w);
    if (edit_menu && sqfp->seqFormat == SEQ_FMT_FASTA) {
      i = CommandItem (edit_menu, "Sequence Deletion Tool", RemoveSequencesFromSubmission);
      SetObjectExtra (i, sqfp, NULL);
      i = CommandItem (edit_menu, "Delete All Sequences", RemoveAllSequencesFromSubmission);
      SetObjectExtra (i, sqfp, NULL);
    }
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      sqfp->appmessage = sepp->handleMessages;
    }

    SetGroupSpacing (w, 10, 10);

    j = HiddenGroup (w, 10, 0, NULL);

    tbs_grp = HiddenGroup (j, 0, 0, NULL);
    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      sqfp->seq_tbs = CreateFolderTabs (tbs_grp, cdnaGenSeqFormTabs, NUCLEOTIDE_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
      sqfp->annot_tbs = CreateFolderTabs (tbs_grp, cdnaGenAnnotFormTabs, 1,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
    } else if (sqfp->seqPackage == SEQ_PKG_GAPPED) {
      sqfp->seq_tbs = CreateFolderTabs (tbs_grp, seqSegSeqFormTabs, NUCLEOTIDE_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
      sqfp->annot_tbs = CreateFolderTabs (tbs_grp, seqSegAnnotFormTabs, 1,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
    } else {
      sqfp->seq_tbs = CreateFolderTabs (tbs_grp, popPhyMutSeqFormTabs, NUCLEOTIDE_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
      sqfp->annot_tbs = CreateFolderTabs (tbs_grp, popPhyMutAnnotFormTabs, 1,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) sqfp->seq_tbs, (HANDLE) sqfp->annot_tbs, NULL);
    Hide (sqfp->annot_tbs);
    sqfp->currentPage = 0;
    page = 0;

    h = HiddenGroup (w, 0, 0, NULL);

    sqfp->pages [page] = CreateNucleotideTab (h, sqfp);
    Hide (sqfp->pages [page]);
    sqfp->tagFromPage [page] = NUCLEOTIDE_PAGE;
    page++;

    sqfp->pages[page] = CreateSequencingMethodTab (h, sqfp);
    Hide (sqfp->pages [page]);
    sqfp->tagFromPage [page] = SEQUENCING_METHOD_PAGE;
    page++;

    sqfp->pages [page] = CreateSourceTab (h, sqfp);
    Hide (sqfp->pages [page]);
    sqfp->tagFromPage [page] = ORGANISM_PAGE;
    page++;
    

    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      sqfp->pages [page] = CreateTranscriptsTab (h, sqfp);
      Hide (sqfp->pages [page]);
      sqfp->tagFromPage [page] = MRNA_PAGE;
      page++;
    }
    
    sqfp->pages [page] = CreateProteinTab (h, sqfp);
    Hide (sqfp->pages [page]);
    sqfp->tagFromPage [page] = PROTEIN_PAGE;
    page++;

    if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA 
        && sqfp->seqPackage != SEQ_PKG_GAPPED)
    {
      sqfp->pages [page] = CreateAnnotTab (h, sqfp);
      Hide (sqfp->pages [page]);
      sqfp->tagFromPage [page] = ANNOTATE_PAGE;
      page++;
    }

    sqfp->numPages = page;

    if (sqfp->seqPackage == SEQ_PKG_TSA) {
      ChooseFeatureTypeForSubmissionFeatureInfoDialog (sqfp->feature_info, eSubmitAnnotType_None);
    }
    
    c = HiddenGroup (w, 3, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    sqfp->goToPrev = goBack;
    sqfp->prevBtn = PushButton (c, " << Prev Form ", PrevSequencesFormBtn);
    SetObjectExtra (sqfp->prevBtn, sqfp, NULL);
    sqfp->goToNext = goToNext;
    sqfp->nextBtn = PushButton (c, " Next Page >> ", NextSequencesFormBtn);
    SetObjectExtra (sqfp->nextBtn, sqfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) j, (HANDLE) c,
                  (HANDLE) sqfp->pages [0], (HANDLE) sqfp->pages [1],
                  (HANDLE) sqfp->pages [2], (HANDLE) sqfp->pages [3], NULL);

    RealizeWindow (w);

    SafeSetTitle (sqfp->prevBtn, "<< Prev Form");
    SafeSetTitle (sqfp->nextBtn, "Next Page >>");

    sqfp->activate = activateForm;
    SetActivate (w, InitOrgNucProtFormActivate);

    SendMessageToDialog (sqfp->seq_tbs, VIB_MSG_INIT);
    SendMessageToDialog (sqfp->annot_tbs, VIB_MSG_INIT);
    SendMessageToDialog (sqfp->dnaseq, VIB_MSG_INIT);
    SendMessageToDialog (sqfp->protseq, VIB_MSG_INIT);

    Show (sqfp->pages [sqfp->currentPage]);
  }
  return (ForM) w;
}


static void SetSequencesFormPage (SequencesFormPtr sqfp, Int4 page)
{
  if (sqfp == NULL) {
    return;
  }
  if (page > -1 && page < 4 && page != sqfp->currentPage) {
    Hide (sqfp->pages[sqfp->currentPage]);
    if (page > 1) {
      sqfp->show_annot = TRUE;
      Hide (sqfp->seq_tbs);
      Show (sqfp->annot_tbs);
      SetValue (sqfp->annot_tbs, page - 2);
    } else {
      sqfp->show_annot = FALSE;
      Show (sqfp->seq_tbs);
      Hide (sqfp->annot_tbs);
      SetValue (sqfp->seq_tbs, page);
    }
  }
}

NLM_EXTERN void SetSequencesForSubmissionForm (WindoW w, SeqEntryPtr sep_list, Int4 page)
{
  SequencesFormPtr sqfp;
  FastaPagePtr     fpp;

  sqfp = GetObjectExtra (w);
  if (sqfp == NULL || sqfp->seqFormat != SEQ_FMT_FASTA || (fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq)) == NULL) {
    return;
  }

  fpp->list = sep_list;
  fpp->errmsgs = ValNodeFreeData (fpp->errmsgs);
  if (fpp->list == NULL) {
    fpp->path [0] = '\0';
    SafeHide (fpp->have_seq_instr_grp);
    Reset (fpp->doc);
    SafeShow (fpp->instructions);
    if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA)
    {
      SetTitle (fpp->import_btn, "Import Nucleotide FASTA");
    }
    Enable (fpp->import_btn);
    Disable (fpp->clear_btn);
  } else {
    SafeHide (fpp->instructions);
    Update ();
    SetTitle (fpp->import_btn, "Import Additional Nucleotide FASTA");
    Enable (fpp->clear_btn);
    Reset (fpp->doc);
    FormatFastaDoc (fpp);
    SafeShow (fpp->have_seq_instr_grp);
  }
  NucleotideImportFinish (sqfp);
  SetSequencesFormPage (sqfp, page);
}

static void MakePubAndDefLine (SequinBlockPtr sbp, SeqEntryPtr sep)

{
  AffilPtr     affil;
  AuthListPtr  alp;
  CitGenPtr    cgp;
  PubdescPtr   pdp;
  ValNodePtr   pep;
  ValNodePtr   vnp;
  /*
  BioseqSetPtr  bssp;
  Char          str [256];
  CharPtr       ttl;
  */

  if (sep == NULL) return;
  /*
  if (SeqEntryGetTitle (sep) != NULL) return;
  ttl = NULL;
  SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
  if (ttl != NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_title);
    if (vnp != NULL) {
      StringNCpy_0 (str, ttl, sizeof (str) - 32);
      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && bssp->_class == 1) {
          StringCat (str, ", and translated products");
        }
      }
      vnp->data.ptrvalue = StringSave (str);
    }
  }
  */
  if (sbp == NULL || sbp->citsubauthors == NULL) return;
  pdp = PubdescNew ();
  if (pdp != NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_pub);
    if (vnp != NULL) {
      vnp->data.ptrvalue = (Pointer) pdp;
      pdp->reftype = 0;
      pep = ValNodeNew (NULL);
      pdp->pub = pep;
      if (pep != NULL) {
        cgp = CitGenNew ();
        if (cgp != NULL) {
          pep->choice = PUB_Gen;
          pep->data.ptrvalue = cgp;
          cgp->cit = StringSave ("unpublished");
          alp = AsnIoMemCopy ((Pointer) sbp->citsubauthors,
                              (AsnReadFunc) AuthListAsnRead,
                              (AsnWriteFunc) AuthListAsnWrite);
          cgp->authors = alp;
          if (alp != NULL) {
            affil = AsnIoMemCopy ((Pointer) sbp->citsubaffil,
                                  (AsnReadFunc) AffilAsnRead,
                                  (AsnWriteFunc) AffilAsnWrite);
            alp->affil = affil;
            if (affil != NULL) {
              affil->phone = MemFree (affil->phone);
              affil->fax = MemFree (affil->fax);
              affil->email = MemFree (affil->email);
            }
          }
          cgp->title = sbp->citsubtitle;
          sbp->citsubtitle = NULL;
        }
      }
    }
  }
}

extern SubmitBlockPtr ConvertSequinBlockToSubmitBlock (SequinBlockPtr sqp);

extern SubmitBlockPtr ConvertSequinBlockToSubmitBlock (SequinBlockPtr sqp)

{
  AffilPtr        affil;
  AuthorPtr       ap;
  AuthListPtr     authors;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  DatePtr         dp;
  CharPtr         os;
  SubmitBlockPtr  sbp;
  Char            str [64];

  sbp = NULL;
  if (sqp != NULL) {
    sbp = SubmitBlockNew ();
    if (sbp != NULL) {
      sbp->subtype = 1;
      os = GetOpSysString ();
      if (os != NULL) {
        sprintf (str, "Sequin %s - %s", SEQUIN_APPLICATION, os);
      } else {
        sprintf (str, "Sequin %s", SEQUIN_APPLICATION);
      }
      sbp->tool = StringSave (str);
      MemFree (os);
      sbp->reldate = sqp->releasedate;
      dp = sbp->reldate;
      if (dp != NULL && dp->data [0] == 1 && dp->data [1] > 0) {
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
      cip = ContactInfoNew ();
      if (cip != NULL) {
        ap = sqp->contactperson;
        cip->contact = ap;
        if (ap != NULL) {
          affil = sqp->citsubaffil;
          if (affil != NULL) {
            if (ap->affil != NULL) {
              affil->phone = MemFree (affil->phone);
              affil->fax = MemFree (affil->fax);
              affil->email = MemFree (affil->email);
              affil->phone = StringSave (ap->affil->phone);
              affil->fax = StringSave (ap->affil->fax);
              affil->email = StringSave (ap->affil->email);
              ap->affil = AffilFree (ap->affil);
            }
            ap->affil = affil;
          }
        }
      }
      sbp->contact = cip;
      csp = CitSubFromContactInfo (cip);
      sbp->cit = csp;
      if (csp != NULL) {
        authors = csp->authors;
        if (authors != NULL) {
          affil = authors->affil;
          authors->affil = NULL;
          csp->authors = AuthListFree (csp->authors);
          csp->authors = sqp->citsubauthors;
          authors = csp->authors;
          if (authors != NULL) {
            authors->affil = affil;
            if (affil != NULL) {
              affil->phone = MemFree (affil->phone);
              affil->fax = MemFree (affil->fax);
              affil->email = MemFree (affil->email);
            }
          }
        }
      }
      sbp->hup = sqp->holduntilpublished;
    }
    MemFree (sqp);
  }
  return sbp;
}

extern Uint2 PackageFormResults (SequinBlockPtr sbp, SeqEntryPtr sep, Boolean makePubAndDefLine)

{
  Uint2         entityID;
  SeqSubmitPtr  ssp;

  entityID = 0;
  if (sep != NULL) {
    if (sbp != NULL) {
      ssp = SeqSubmitNew ();
      if (ssp != NULL) {
        ssp->datatype = 1;
        ssp->data = (Pointer) sep;
        if (makePubAndDefLine) {
          MakePubAndDefLine (sbp, sep);
        }
        sbp->citsubtitle = MemFree (sbp->citsubtitle);
        ssp->sub = ConvertSequinBlockToSubmitBlock (sbp);
        ObjMgrConnect (OBJ_SEQENTRY, sep->data.ptrvalue, OBJ_SEQSUB, (Pointer) ssp);
        if (! ObjMgrRegister (OBJ_SEQSUB, (Pointer) ssp)) {
          ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
        }
      } else {
        if (! ObjMgrRegister (OBJ_SEQENTRY, (Pointer) sep)) {
          ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
        }
      }
    } else {
      if (! ObjMgrRegister (OBJ_SEQENTRY, (Pointer) sep)) {
        ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
      }
    }
    if (EntrezASN1Detected (sep)) {
      ErrPostEx (SEV_WARNING, 0, 0, "This record was retrieved from Entrez.");
    }
    entityID = ObjMgrGetEntityIDForChoice (sep);
  }
  return entityID;
}

static void GetRawBsps (BioseqPtr bsp, Pointer userdata)

{
  ValNodePtr PNTR  head;

  if (bsp->repr != Seq_repr_raw) return;
  head = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (head, 0, (Pointer) bsp);
}

extern void ParseInMoreProteins (IteM i)
{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  NucProtAssocPtr   nuc_prot_assoc_list;
  CharPtr      extension;
  Char         path [PATH_MAX];
  FILE         *fp;
  Boolean      parseSeqId;
  SeqEntryPtr  sep, prot_list = NULL, nuc_list = NULL;
  Boolean      chars_stripped = FALSE;
  Boolean      makeMRNA = FALSE;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  extension = GetAppProperty ("FastaProtExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  ans = Message (MSG_YN, "Do FASTA definition lines start with seqID?");
  parseSeqId = (Boolean) (ans == ANS_YES);

  WatchCursor();
  Update();
  /* read in proteins */
  prot_list = ImportSequencesFromFile (fp, NULL, FALSE, parseSeqId, NULL, NULL, &chars_stripped);
  FileClose (fp);  
  if (chars_stripped && prot_list != NULL)
  {
    if (ANS_CANCEL == Message (MSG_OKC, "Illegal characters will be stripped from your sequence data.  Do you want to continue?"))
    {
      prot_list = SeqEntryFree (prot_list);
      ArrowCursor ();
      Update ();
      return;
    }
  }
  
  /* get list of nucleotide sequences */
  nuc_list = sep;

  ArrowCursor();
  Update();
  nuc_prot_assoc_list = AssignProteinsForSequenceSet (nuc_list, prot_list, TRUE);
  if (nuc_prot_assoc_list == NULL) {
    prot_list = SeqEntryFree (prot_list);
    return;
  } 

  ans = Message (MSG_YN, "Do you wish to make default mRNAs?");
  makeMRNA = (Boolean) (ans == ANS_YES);
  
  ReplaceDuplicateProteinIDs (nuc_list, prot_list);
  AssignProteinsToSelectedNucleotides (nuc_prot_assoc_list,
                                       nuc_list,
                                       prot_list,
                                       NULL,
                                       1, makeMRNA);

  nuc_prot_assoc_list = FreeAssociationList (nuc_prot_assoc_list);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();  
}

extern void ParseInMoreMRNAs (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  Int4         count;
  CharPtr      errormsg;
  CharPtr      extension;
  FILE         *fp;
  ValNodePtr   head;
  Boolean      isLocalUnknownID;
  SeqEntryPtr  last;
  SeqEntryPtr  list;
  MonitorPtr   mon;
  SeqEntryPtr  nextsep;
  SeqEntryPtr  nucsep;
  ObjectIdPtr  oid;
  Boolean      parseSeqId;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  Char         str [32];
  Char         tmp [128];
  ValNodePtr   vnp;

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

  ans = Message (MSG_YN, "Do FASTA definition lines start with seqID?");
  parseSeqId = (Boolean) (ans == ANS_YES);

  WatchCursor ();
  Update ();

  count = 0;
  list = NULL;
  last = NULL;
  head = NULL;
  errormsg = NULL;

  nucsep = FindNucSeqEntry (sep);
  slp = CreateWholeInterval (sep);

  nextsep = SequinFastaToSeqEntryEx (fp, TRUE, &errormsg, parseSeqId, NULL);
  while (nextsep != NULL) {
    count++;
    if (IS_Bioseq (nextsep) && nextsep->data.ptrvalue != NULL) {
      bsp = (BioseqPtr) nextsep->data.ptrvalue;
      isLocalUnknownID = FALSE;
      sip = bsp->id;
      if (sip != NULL && sip->choice == SEQID_LOCAL) {
        oid = (ObjectIdPtr) sip->data.ptrvalue;
        if (oid != NULL && oid->str != NULL) {
          isLocalUnknownID = (Boolean) (StringICmp (oid->str, "unknown") == 0);
        }
      }
      if ((! parseSeqId) || isLocalUnknownID) {
        sip = MakeNewProteinSeqId (slp, NULL);
        if (sip != NULL) {
          bsp->id = SeqIdFree (bsp->id);
          bsp->id = sip;
          SeqMgrReplaceInBioseqIndex (bsp);
        }
      }
    }
    SeqEntryPack (nextsep);
    if (last != NULL) {
      last->next = nextsep;
      last = nextsep;
    } else {
      list = nextsep;
      last = list;
    }
    if (! StringHasNoText (errormsg)) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->data.ptrvalue = errormsg;
        errormsg = NULL;
      }
    }
    nextsep = SequinFastaToSeqEntryEx (fp, TRUE, &errormsg, parseSeqId, NULL);
  }

  SeqLocFree (slp);
  FileClose (fp);

  ArrowCursor ();
  Update ();

  if (head != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      Message (MSG_POSTERR, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }
    ValNodeFreeData (head);
    ans = Message (MSG_YN, "Errors detected - do you wish to proceed?");
    if (ans == ANS_NO) {
      sep = list;
      while (sep != NULL) {
        nextsep = sep->next;
        sep->next = NULL;
        SeqEntryFree (sep);
        sep = nextsep;
      }
      return;
    }
  }

  if (list == NULL) return;
  
  WatchCursor ();
  Update ();

  nucsep = FindNucSeqEntry (sep);
  if (nucsep == NULL) return;

  mon = MonitorStrNewEx ("Reading mRNA sequences", 20, FALSE);
  count = 0;
  while (list != NULL) {
    nextsep = list->next;
    list->next = NULL;
    count++;
    if (mon != NULL) {
      str [0] = '\0';
      tmp [0] = '\0';
      bsp = (BioseqPtr) list->data.ptrvalue;
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
        SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      }
      sprintf (str, "Processing sequence %d [%s]", (int) count, tmp);
      MonitorStrValue (mon, str);
      Update ();
    }
    AutomaticMrnaProcess (nucsep, list, FALSE, FALSE);
    SeqEntryFree (list);
    list = nextsep;
  }
  mon = MonitorFree (mon);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

/*#ifdef ALLOW_DOWNLOAD*/
typedef struct fetchform {
  FORM_MESSAGE_BLOCK
  GrouP           accntype;
  TexT            accession;
  ButtoN          accept;
} FetchForm, PNTR FetchFormPtr;

static void FetchFormMessage (ForM f, Int2 mssg)

{
  FetchFormPtr  ffp;

  ffp = (FetchFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    switch (mssg) {
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
        if (ffp->appmessage != NULL) {
          ffp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void ExamineIdProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  Int2       i;
  BoolPtr    idTypes;
  SeqIdPtr   sip;

  if (sep == NULL || sep->data.ptrvalue == NULL || mydata == NULL) return;
  idTypes = (BoolPtr) mydata;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sip = bsp->id;
    while (sip != NULL) {
      i = (Int2) sip->choice;
      if (i >= 0 && i < NUM_SEQID) {
        (idTypes [i])++;
      }
      sip = sip->next;
    }
  }
}

static Boolean OwnedByOtherDatabase (SeqEntryPtr sep, BoolPtr idTypes)

{
  Int2  i;

  if (sep == NULL || idTypes == NULL) return FALSE;
  for (i = 0; i < NUM_SEQID; i++) {
    idTypes [i] = FALSE;
  }
  BioseqExplore (sep, (Pointer) idTypes, ExamineIdProc);
  if (! (idTypes [SEQID_GENBANK])) return TRUE;
  if (idTypes [SEQID_EMBL] || idTypes [SEQID_DDBJ]) return TRUE;
  if (! FindNucSeqEntry (sep)) return TRUE;
  return FALSE;
}

static Int4 AccessionToGi (CharPtr string)
{
   /*
   CharPtr str;
   LinkSetPtr lsp;
   Int4 gi;

   str = MemNew (StringLen (string) + 10);
   sprintf (str, "\"%s\" [ACCN]", string);
   lsp = EntrezTLEvalString (str, TYP_NT, -1, NULL, NULL);
   MemFree (str);
   if (lsp == NULL) return 0;
   if (lsp->num <= 0) {
       LinkSetFree (lsp);
       return 0;
   }
   gi = lsp->uids [0];
   LinkSetFree (lsp);
   return gi;
   */
   Int4      gi;
   SeqIdPtr  sip;

   sip = SeqIdFromAccessionDotVersion (string);
   if (sip == NULL) return 0;
   gi = GetGIForSeqId (sip);
   SeqIdFree (sip);
   return gi;
}

static void LookForReplacedByCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr   bsp;
  SeqHistPtr  hist;
  BoolPtr     rsult;

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  hist = bsp->hist;
  if (hist == NULL) return;
  if (hist->replaced_by_ids != NULL) {
    rsult = (BoolPtr) mydata;
    if (rsult == NULL) return;
    *rsult = TRUE;
  }
}

#ifdef USE_SMARTNET
extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
#endif


extern void LaunchDisplay (Uint2 entityID)
{
  Int2          handled;
  Char          str [32];

  seqviewprocs.filepath = str;
  seqviewprocs.forceSeparateViewer = TRUE;
  handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
  seqviewprocs.filepath = NULL;
  ArrowCursor ();
  if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
    Message (MSG_FATAL, "Unable to launch viewer.");
    return;
  } else {
    SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
  }
  ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
  ObjMgrSetDirtyFlag (entityID, TRUE);
}

extern void DownloadAndDisplay (Int4 uid)
{
  CharPtr       accn = NULL;
  MsgAnswer     ans;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype;
  CharPtr       dbname;
  Uint2         entityID;
  Int2          handled;
  Boolean       idTypes [NUM_SEQID];
  Boolean       isReplaced = FALSE;
  SeqEntryPtr   sep;
  Char          str [32];

  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, /* -1 */ 0);
    /* EntrezFini (); */
    if (sep == NULL) {
      ArrowCursor ();
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
    if (IS_Bioseq (sep)) {
      datatype = OBJ_BIOSEQ;
    } else if (IS_Bioseq_set (sep)) {
      datatype = OBJ_BIOSEQSET;
    } else {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
    SeqEntryExplore (sep, (Pointer) (&isReplaced), LookForReplacedByCallback);
    if (isReplaced) {
      ans = Message (MSG_YN, "This record has been replaced.  Are you sure you want to edit it?");
      if (ans == ANS_NO) {
        SeqEntryFree (sep);
        ArrowCursor ();
        return;
      }
    }
    dataptr = (Pointer) sep->data.ptrvalue;
  } else if (! StringHasNoText (accn)) {
#ifdef USE_SMARTNET
    if (accn != NULL) {
      dataptr = ReadFromTPASmart (accn, &datatype, NULL);
      if (dataptr == NULL) {
        dataptr = ReadFromSmart (accn, &datatype, NULL);
        if (dataptr == NULL) {
          dataptr = ReadFromDirSub (accn, &datatype, NULL);
        }
      }
    }
#endif
  }
  if (dataptr != NULL) {
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
        if (sep != NULL && OwnedByOtherDatabase (sep, idTypes)) {
          dbname = NULL;
          if (idTypes [SEQID_EMBL]) {
            dbname = "EMBL";
          } else if (idTypes [SEQID_DDBJ]) {
            dbname = "DDBJ";
          }
        }
        seqviewprocs.filepath = str;
        seqviewprocs.forceSeparateViewer = TRUE;
        handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                    OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
        seqviewprocs.filepath = NULL;
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
        ArrowCursor ();
      }
    } else {
      ArrowCursor ();
    }
  } else {
    /* EntrezFini (); */
    ArrowCursor ();
    Message (MSG_OK, "Unable to find this record in the database.");
  }
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


static Boolean StrToULong (CharPtr str, Uint4Ptr longval)

{
  Char           ch;
  Int2           i;
  Int2           len;
  Char           local [64];
  Boolean        nodigits;
  Boolean        rsult;
  unsigned long val;

  rsult = FALSE;
  if (longval != NULL) {
    *longval = (Uint4) 0;
  }
  len = (Int2) StringLen (str);
  if (len != 0) {
    rsult = TRUE;
    nodigits = TRUE;
    for (i = 0; i < len; i++) {
      ch = str [i];
      if (ch == ' ' || ch == '+' || ch == '-') {
      } else if (ch < '0' || ch > '9') {
        rsult = FALSE;
      } else {
        nodigits = FALSE;
      }
    }
    if (nodigits) {
      rsult = FALSE;
    }
    if (rsult && longval != NULL) {
      StringNCpy_0 (local, str, sizeof (local));
      if (sscanf (local, "%lu", &val) == 1) {
        *longval = val;
      }
    }
  }
  return rsult;
}


static void DownloadProc (ButtoN b)

{
  CharPtr       accn = NULL;
  MsgAnswer     ans;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype;
  CharPtr       dbname;
  Uint2         entityID;
  FetchFormPtr  ffp;
  Int4          flags = 0;
  Int2          handled;
  Boolean       idTypes [NUM_SEQID];
  Boolean       isReplaced = FALSE;
  Int2          retcode = 0;
  SeqEntryPtr   sep;
  Char          str [64];
  Uint4         tid;
  CharPtr       tmp1 = NULL;
  CharPtr       tmp2 = NULL;
  Int4          uid;
  long int      val;
  ForM          w;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  w = ffp->form;
  Hide (w);
  WatchCursor ();
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  TrimSpacesAroundString (str);
  if (StringHasNoText (str)) {
    Message (MSG_OK, "Please enter an accession number or gi");
    Show (w);
    Select (w);
    Select (ffp->accession);
    return;
  }

  tmp1 = StringChr (str, ',');
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

  sep = NULL;
  uid = 0;
  tid = 0;
  /*
  if (! EntrezIsInited ()) {
    if (! SequinEntrezInit ("Sequin", FALSE, NULL)) {
      Remove (w);
      Show (startupForm);
      Select (startupForm);
      ArrowCursor ();
      return;
    }
  }
  */
  if (GetValue (ffp->accntype) == 1) {
    /*
    sip = ValNodeNew (NULL);
    if (sip != NULL) {
      tsip = TextSeqIdNew ();
      if (tsip != NULL) {
        tsip->accession = StringSave (str);
        sip->choice = SEQID_GENBANK;
        sip->data.ptrvalue = (Pointer) tsip;
        uid = EntrezFindSeqId (sip);
        if (uid == 0) {
          sip->choice = SEQID_EMBL;
          uid = EntrezFindSeqId (sip);
        }
        if (uid == 0) {
          sip->choice = SEQID_DDBJ;
          uid = EntrezFindSeqId (sip);
        }
      }
    }
    SeqIdFree (sip);
    */
    if (StringNICmp (str, "ti|", 3) == 0 && IsAllDigits (str + 3)) {
      if (! StrToULong (str + 3, &tid)) {
        tid = 0;
      }
      if (tid > 0) {
        sep = PubSeqSynchronousQueryTI (tid, 0, /* -1 */ 0);
      }
    } else {      
      uid = AccessionToGi (str);
      accn = str;
    }
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, retcode, flags);
  }
  if (uid > 0 || tid > 0) {
    /* EntrezFini (); */
    if (sep == NULL) {
      ArrowCursor ();
      Message (MSG_OK, "Unable to find this record in the database.");
      Show (w);
      Select (w);
      Select (ffp->accession);
      return;
    }
    if (IS_Bioseq (sep)) {
      datatype = OBJ_BIOSEQ;
    } else if (IS_Bioseq_set (sep)) {
      datatype = OBJ_BIOSEQSET;
    } else {
      ArrowCursor ();
      Message (MSG_OK, "Unable to find this record in the database.");
      Show (w);
      Select (w);
      Select (ffp->accession);
      return;
    }
    Remove (w);
    SeqEntryExplore (sep, (Pointer) (&isReplaced), LookForReplacedByCallback);
    if (isReplaced) {
      ans = Message (MSG_YN, "This record has been replaced.  Are you sure you want to edit it?");
      if (ans == ANS_NO) {
        SeqEntryFree (sep);
        Show (startupForm);
        Select (startupForm);
        ArrowCursor ();
        return;
      }
    }
    dataptr = (Pointer) sep->data.ptrvalue;
  } else if (! StringHasNoText (accn)) {
#ifdef USE_SMARTNET
    if (accn != NULL) {
      dataptr = ReadFromTPASmart (accn, &datatype, NULL);
      if (dataptr == NULL) {
        dataptr = ReadFromSmart (accn, &datatype, NULL);
        if (dataptr == NULL) {
          dataptr = ReadFromDirSub (accn, &datatype, NULL);
        }
      }
    }
#endif
  }
  if (dataptr != NULL) {
    entityID = ObjMgrRegister (datatype, dataptr);
    SeqEntrySetScope (NULL);
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
        if (sep != NULL && OwnedByOtherDatabase (sep, idTypes)) {
          dbname = NULL;
          if (idTypes [SEQID_EMBL]) {
            dbname = "EMBL";
          } else if (idTypes [SEQID_DDBJ]) {
            dbname = "DDBJ";
          }
        }
        if (datatype != OBJ_SEQSUB && uid > 0) {
          ArrowCursor ();
          Update ();
          if (!indexerVersion && Message (MSG_YN, repackageMsg) == ANS_YES) {
            globalEntityID = entityID;
            globalsep = sep;
            StringNCpy_0 (globalPath, str, sizeof (globalPath));
            WatchCursor ();
            Update ();
            w = CreateSubmitBlockForm (-50, -33, "Submitting Authors",
                                       FALSE, TRUE, NULL, JustRegisterSeqEntryBtn,
                                       AddSubmitBlockToSeqEntry);
            ArrowCursor ();
            if (w != NULL) {
              Show (w);
              Select (w);
              SendHelpScrollMessage (helpForm, "Submitting Authors Form", NULL);
              return;
            } else {
              Message (MSG_FATAL, "Unable to create window.");
              SeqEntryFree (sep);
              Show (startupForm);
              Select (startupForm);
              return;
            }
          }
        }
        seqviewprocs.filepath = str;
        seqviewprocs.forceSeparateViewer = TRUE;
        handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                    OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
        seqviewprocs.filepath = NULL;
        ArrowCursor ();
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
          Message (MSG_FATAL, "Unable to launch viewer.");
          SeqEntryFree (sep);
          Show (startupForm);
          Select (startupForm);
          return;
        } else {
          SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
        }
        ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
        ObjMgrSetDirtyFlag (entityID, TRUE);
      } else {
        Message (MSG_ERROR, "Unable to process object type %d.", (int) datatype);
        ObjMgrDelete (datatype, dataptr);
        Show (startupForm);
        Select (startupForm);
        ArrowCursor ();
      }
    } else {
      Show (startupForm);
      Select (startupForm);
      ArrowCursor ();
    }
  } else {
    /* EntrezFini (); */
    ArrowCursor ();
    Message (MSG_OK, "Unable to find this record in the database.");
    Show (w);
    Select (w);
    Select (ffp->accession);
  }
}

extern void DownloadAndUpdateProc (ButtoN b)

{
  FetchFormPtr  ffp;
  SeqEntryPtr   sep;
  Char          str [32];
  Int4          uid;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ParentWindow (b));
  WatchCursor ();
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  if (StringHasNoText (str)) {
    Remove (ParentWindow (b));
    ArrowCursor ();
    return;
  }
  sep = NULL;
  uid = 0;
  /*
  if (! EntrezIsInited ()) {
    if (! SequinEntrezInit ("Sequin", FALSE, NULL)) {
      Remove (ParentWindow (b));
      ArrowCursor ();
      return;
    }
  }
  */
  if (GetValue (ffp->accntype) == 1) {
    uid = AccessionToGi (str);
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  Remove (ParentWindow (b));
  ArrowCursor ();
  Update ();
  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, /* -1 */ 0);
    /* EntrezFini (); */
    if (sep == NULL) {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
    if (IS_Bioseq (sep)) {
    } else if (IS_Bioseq_set (sep)) {
    } else {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
  }

  SqnReadAlignView ((BaseFormPtr) ffp, updateTargetBspKludge, sep, TRUE);
}

extern void DownloadAndExtendProc (ButtoN b)

{
  FetchFormPtr  ffp;
  SeqEntryPtr   sep;
  Char          str [32];
  Int4          uid;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ParentWindow (b));
  WatchCursor ();
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  if (StringHasNoText (str)) {
    Remove (ParentWindow (b));
    ArrowCursor ();
    return;
  }
  sep = NULL;
  uid = 0;
  /*
  if (! EntrezIsInited ()) {
    if (! SequinEntrezInit ("Sequin", FALSE, NULL)) {
      Remove (ParentWindow (b));
      ArrowCursor ();
      return;
    }
  }
  */
  if (GetValue (ffp->accntype) == 1) {
    uid = AccessionToGi (str);
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  Remove (ParentWindow (b));
  ArrowCursor ();
  Update ();
  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, /* -1 */ 0);
    /* EntrezFini (); */
    if (sep == NULL) {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
    if (IS_Bioseq (sep)) {
    } else if (IS_Bioseq_set (sep)) {
    } else {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
  }

  SqnReadAlignView ((BaseFormPtr) ffp, updateTargetBspKludge, sep, FALSE);
}

static void CancelFetchProc (ButtoN b)

{
  StdCancelButtonProc (b);
  Show (startupForm);
  Select (startupForm);
}


static void FetchTextProc (TexT t)

{
  Boolean       alldigits;
  FetchFormPtr  ffp;
  CharPtr       ptr;
  Char          str [64];

  ffp = (FetchFormPtr) GetObjectExtra (t);
  if (ffp == NULL) return;
  GetTitle (t, str, sizeof (str));
  if (StringHasNoText (str)) {
    SafeDisable (ffp->accept);
  } else {
    SafeEnable (ffp->accept);
    ptr = StringChr (str, ',');
    if (ptr != NULL) {
      *ptr = '\0';
    }
    TrimSpacesAroundString (str);
    alldigits = IsAllDigits (str);
    if (alldigits) {
      SafeSetValue (ffp->accntype, 2);
    } else {
      SafeSetValue (ffp->accntype, 1);
    }
  }
}

extern void CommonFetchFromNet (BtnActnProc actn, BtnActnProc cancel)

{
  GrouP              c;
  FetchFormPtr       ffp;
  GrouP              g;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  Hide (startupForm);
  Update ();
  w = NULL;
  ffp = MemNew (sizeof (FetchForm));
  if (ffp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Download From Entrez", NULL);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->formmessage = FetchFormMessage;

#ifndef WIN_MAC
    CreateSqnInitialFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      ffp->appmessage = sepp->handleMessages;
    }
    SetGroupSpacing (w, 10, 10);

    g = HiddenGroup (w, -3, 0, NULL);
    StaticPrompt (g, "Type", 0, stdLineHeight, programFont, 'l');
    ffp->accntype = HiddenGroup (g, 4, 0, NULL);
    RadioButton (ffp->accntype, "Accession");
    RadioButton (ffp->accntype, "GI");
    SetValue (ffp->accntype, 1);
    ffp->accession = DialogText (g, "", 6, FetchTextProc);
    SetObjectExtra (ffp->accession, ffp, NULL);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    ffp->accept = DefaultButton (c, "Retrieve", actn);
    SetObjectExtra (ffp->accept, ffp, NULL);
    Disable (ffp->accept);
    PushButton (c, "Cancel", cancel);

    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
    }
    Select (ffp->accession);
    Show (w);
    Select (w);
    Update ();
  } else {
    Show (startupForm);
    Select (startupForm);
  }
}

extern void FetchFromNet (ButtoN b)

{
  CommonFetchFromNet (DownloadProc, CancelFetchProc);
}

/*#else
#define FetchFromNet NULL
#endif*/

/*
static Boolean FindPerfectSubMatch (CharPtr prot, CharPtr trans, Int4 start,
                                    Int4 len, Uint1 frame, Int2 strand,
                                    Int4Ptr fromPtr, Int4Ptr toPtr)

{
  int      ch;
  Int2     d [256];
  Int4     from;
  int      i;
  int      j;
  int      k;
  size_t   protLen;
  Boolean  rsult;
  Int4     to;
  size_t   transLen;

  rsult = FALSE;
  from = 0;
  to = 0;
  if (prot != NULL && trans != NULL) {
    protLen = StringLen (prot);
    transLen = StringLen (trans);
    if (protLen <= transLen) {
      for (ch = 0; ch < 256; ch++) {
        d [ch] = protLen;
      }
      for (j = 0; j < protLen - 1; j++) {
        d [(int) prot [j]] = protLen - j - 1;
      }
      i = protLen;
      do {
        j = protLen;
        k = i;
        do {
          k--;
          j--;
        } while (j >= 0 && prot [j] == trans [k]);
        if (j >= 0) {
          i += d [(int) trans [i - 1]];
        }
      } while (j >= 0 && i <= transLen);
      if (j < 0) {
        i -= protLen;
        from = (long) (i * 3 + (frame - 1));
        to = from + 3 * protLen;
        if (trans [i + protLen] == '*') {
          to += 3;
        }
        if (strand == Seq_strand_plus) {
          from += 1;
        } else if (strand == Seq_strand_minus) {
          from = len - from;
          to = len - to + 1;
        }
        rsult = TRUE;
      }
    }
  }
  if (fromPtr != NULL) {
    *fromPtr = from + start;
  }
  if (toPtr != NULL) {
    *toPtr = to + start;
  }
  return rsult;
}

static Boolean CheckOneFrame (BioseqPtr bsp, Int4 start, Int4 len,
                              CharPtr prot, Int2 gencode,
                              Uint1 frame, Int2 strand,
                              Int4Ptr fromPtr, Int4Ptr toPtr)

{
  ByteStorePtr  bs;
  Char          ch;
  ValNodePtr    code;
  CdRegionPtr   crp;
  CharPtr       ptr;
  Boolean       rsult;
  SeqFeatPtr    sfp;
  CharPtr       trans;
  ValNodePtr    vnp;

  rsult = FALSE;
  if (bsp != NULL && gencode > 0) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_CDREGION;
      crp = CdRegionNew ();
      sfp->data.value.ptrvalue = (Pointer) crp;
      if (crp != NULL) {
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
        AddIntToSeqFeat (sfp, start, start + len - 1, bsp, -1, -1, strand);
        trans = NULL;
        bs = ProteinFromCdRegion (sfp, TRUE);
        if (bs != NULL) {
          trans = BSMerge (bs, NULL);
          BSFree (bs);
        }
        if (trans != NULL) {
          ptr = trans;
          ch = *ptr;
          while (ch != '\0') {
            *ptr = TO_UPPER (ch);
            ptr++;
            ch = *ptr;
          }
          if (trans [0] == '-') {
            trans [0] = prot [0];
          }
          rsult = FindPerfectSubMatch (prot, trans, start, len,
                                       frame, strand, fromPtr, toPtr);
          MemFree (trans);
        }
      }
      SeqFeatFree (sfp);
    }
  }
  return rsult;
}

#define PREDICT_BLOCK_SIZE 30000L

static SeqLocPtr FindSingleCodingInterval (BioseqPtr nuc, BioseqPtr prot, Int2 genCode)

{
  Int4        cdsFrom;
  Int4        cdsTo;
  Char        ch;
  Int4        cntr;
  Uint1       frame;
  Int4        from;
  Int4        incr;
  Int4        len;
  Boolean     matched;
  size_t      protLen;
  CharPtr     protstr;
  CharPtr     ptr;
  SeqFeatPtr  sfp;
  SeqLocPtr   slp;
  Int4        start;
  Int2        strand;
  Int4        tmp;
  Int4        to;

  slp = NULL;
  if (nuc != NULL && prot != NULL) {
    cdsFrom = 0;
    cdsTo = 0;
    strand = Seq_strand_unknown;
    protstr = NULL;
    if (prot->length > 0) {
      protstr = BSMerge (prot->seq_data, NULL);
      if (protstr != NULL) {
        ptr = protstr;
        ch = *ptr;
        while (ch != '\0') {
          *ptr = TO_UPPER (ch);
          ptr++;
          ch = *ptr;
        }
        protLen = StringLen (protstr);
        matched = FALSE;
        for (frame = 1; frame <= 3 && (! matched); frame++) {
          strand = Seq_strand_plus;
          start = 0;
          cntr = nuc->length;
          len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          while (len > 0 && (! matched)) {
            incr = MIN (cntr, PREDICT_BLOCK_SIZE);
            matched = CheckOneFrame (nuc, start, len, protstr, genCode, frame,
                                     strand, &cdsFrom, &cdsTo);
            start += incr;
            cntr -= incr;
            len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          }
        }
        for (frame = 1; frame <= 3 && (! matched); frame++) {
          strand = Seq_strand_minus;
          start = 0;
          cntr = nuc->length;
          len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          while (len > 0 && (! matched)) {
            incr = MIN (cntr, PREDICT_BLOCK_SIZE);
            matched = CheckOneFrame (nuc, start, len, protstr, genCode, frame,
                                     strand, &cdsFrom, &cdsTo);
            start += incr;
            cntr -= incr;
            len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          }
        }
        if (matched) {
          sfp = SeqFeatNew ();
          if (sfp != NULL) {
            from = cdsFrom - 1;
            to = cdsTo - 1;
            if (from > to) {
              tmp = from;
              from = to;
              to = tmp;
            }
            AddIntToSeqFeat (sfp, from, to, nuc, -1, -1, strand);
            slp = sfp->location;
            sfp->location = NULL;
          }
          SeqFeatFree (sfp);
        }
      }
      MemFree (protstr);
    }
  }
  return slp;
}
*/

static BioseqPtr FindBioseqByIDString (
  CharPtr    seqid,
  SeqEntryPtr sep
)
{
  BioseqSetPtr      bssp;
  BioseqPtr         bsp;
  CharPtr           str;

  if (sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    /* this also delves into nuc-prot sets */
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)) ||
                         bssp->_class == 1)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        bsp = FindBioseqByIDString (seqid, sep);
        if (bsp != NULL) return bsp;
      }
      return NULL;
    }
  }
  if (!IS_Bioseq (sep)) return FALSE;

  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return NULL;
  str = SeqIdWholeLabel (bsp->id, PRINTID_REPORT);
  if (StringCmp (str, seqid) != 0) {
    bsp = NULL;
  }
  str = MemFree (str);
  return bsp;
}


typedef struct tpafromfile {
  BioseqPtr bsp;
  ValNodePtr acc_list;
} TPAFromFileData, PNTR TPAFromFilePtr;


static TPAFromFilePtr TPAFromFileNew (BioseqPtr bsp, ValNodePtr acc_list)
{
  TPAFromFilePtr tpfp;

  tpfp = (TPAFromFilePtr) MemNew (sizeof (TPAFromFileData));
  tpfp->bsp = bsp;
  tpfp->acc_list = acc_list;
  return tpfp;
}

static TPAFromFilePtr TPAFromFileFree (TPAFromFilePtr tpfp)
{
  if (tpfp != NULL) 
  {
    tpfp->acc_list = ValNodeFreeData (tpfp->acc_list);
    tpfp = MemFree (tpfp);
  }
  return tpfp;
}


static ValNodePtr TPAFromFileListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) 
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = TPAFromFileFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;  
}


NLM_EXTERN UserObjectPtr GetTPAAssembly (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  UserObjectPtr     uop = NULL;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL && uop == NULL)
  {
    uop = (UserObjectPtr)sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringICmp (uop->type->str, "TpaAssembly") == 0)
    {
      return uop;
    }
    else
    { 
      uop = NULL;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }
  return uop;
}


static Boolean HasTPAAccessions (BioseqPtr bsp)
{
  UserObjectPtr uop;
  Boolean       found = FALSE;
  UserFieldPtr  ufp;

  uop = GetTPAAssembly (bsp);
  if (uop != NULL)
  {
    ufp = uop->data;
    while (ufp != NULL && ufp->choice != 11) 
    {
      ufp = ufp->next;
    }
    if (ufp != NULL)
    {
      found = TRUE;
    }
  }
  return found;
} 


static ValNodePtr ReadTPAAccessionNumbersFile (FILE *fp, SeqEntryPtr sep, BoolPtr replace)
{
  LogInfoPtr notfound_lip, alreadyhas_lip;
  Boolean    need_seqid;
  ValNodePtr tpalist = NULL;
  BioseqPtr      bsp;
  Char          str [8192];
  size_t        len = 8192;
  Char          seqid[100];
  Int4          seqid_len;
  CharPtr       cp;
  CharPtr       acc_end;
  Boolean       found_end;
  ValNodePtr    acc_list;
  Boolean       found_other_acc = FALSE;
  MsgAnswer     ans;
  SeqIdPtr      sip;

  notfound_lip = OpenLog ("Accessions in Table Not Found in Record");
  alreadyhas_lip = OpenLog ("Accessions that Already Have TPA Acccesion Numbers");

  need_seqid = TRUE;
  acc_list = NULL;
  ReadLine (fp, str, len);
  while (Nlm_fileDone) 
  {
    cp = str;
    if (strlen (str) == 0)
    {
      ReadLine (fp, str, len);
      continue;
    }
    if (need_seqid)
    {
      seqid_len = StringCSpn (str, " \t");
      if (seqid_len > 0)
      {
        StringNCpy (seqid, str, seqid_len);
        seqid [seqid_len] = 0;
        need_seqid = FALSE;
      }
      cp = str + seqid_len + 1;
    }
    if (need_seqid)
    {
      ReadLine (fp, str, len);
      continue;
    }
    if (str [strlen (str) - 1] != ',')
    {
      need_seqid = TRUE;
    }
    
    found_end = FALSE;
    while (*cp != 0)
    {
      if (*cp == ' ' || *cp == ',' || *cp == '\t')
      {
        cp++;
      }
      else
      {
        acc_end = cp + 1;
        while (*acc_end != 0 && *acc_end != ',')
        {
          acc_end++;
        }
        if (*acc_end == 0)
        {
          found_end = TRUE;
        }
        else
        {
          *acc_end = 0;
        }
        ValNodeAddStr (&acc_list, 0, StringSave (cp));
        if (found_end)
        {
          cp = acc_end;
        }
        else
        {
          cp = acc_end + 1;
        }
      }
    }

    if (need_seqid == TRUE)
    {
      sip = CreateSeqIdFromText (seqid, sep);
      bsp = BioseqFind (sip);
      sip = SeqIdFree (sip);
      if (bsp == NULL) {
        bsp = FindBioseqByIDString (seqid, sep);
      }
      if (bsp == NULL) 
      {
        fprintf (notfound_lip->fp, 
                 "Unable to update accession numbers for %s (not found)\n",
                 seqid);
        notfound_lip->data_in_log = TRUE;
        acc_list = ValNodeFreeData (acc_list);
      }
      else
      {
        ValNodeAddPointer (&tpalist, 0, TPAFromFileNew (bsp, acc_list));
        if (HasTPAAccessions (bsp)) 
        {
          fprintf (alreadyhas_lip->fp,
                   "Already have accessions for %s\n",
                   seqid);
          alreadyhas_lip->data_in_log = TRUE;
        }
      }
      acc_list = NULL;
    }
      
    ReadLine (fp, str, len);
  }
  if (acc_list != NULL)
  {
    bsp = FindBioseqByIDString (seqid, sep);
    if (bsp == NULL)
    {
      fprintf (notfound_lip->fp,
               "Unable to update accession numbers for %s (not found)\n",
               seqid);
      notfound_lip->data_in_log = TRUE;
      acc_list = ValNodeFreeData (acc_list);
    }
    else
    {
      ValNodeAddPointer (&tpalist, 0, TPAFromFileNew (bsp, acc_list));
      if (GetTPAAssembly (bsp) != NULL) 
      {
        fprintf (alreadyhas_lip->fp,
                 "Already have accessions for %s\n",
                 seqid);
        alreadyhas_lip->data_in_log = TRUE;
      }
    }
    acc_list = NULL;
  }

  CloseLog (notfound_lip);
  notfound_lip = FreeLog (notfound_lip);
  found_other_acc = alreadyhas_lip->data_in_log;
  CloseLog (alreadyhas_lip);
  alreadyhas_lip = FreeLog (alreadyhas_lip);

  if (tpalist == NULL) 
  {
    Message (MSG_ERROR, "Unable to read accession lists from file!");
  }
  else if (replace != NULL)
  {
    if (found_other_acc) 
    {
      ans = Message (MSG_YNC, "Some accessions already have TPA accession numbers.  Remove these first?");
      if (ans == ANS_CANCEL)
      {
        tpalist = TPAFromFileListFree (tpalist);
      }
      else if (ans == ANS_YES) 
      {
        *replace = TRUE;
      }
      else
      {
        *replace = FALSE;
      }
    } 
    else
    {
      *replace = FALSE;
    }
  }

  return tpalist;
}


NLM_EXTERN void RemoveOldTPAAccessions (UserObjectPtr uop)
{
  UserFieldPtr  ufp, ufp_next, ufp_prev = NULL;
  ObjectIdPtr   oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "TpaAssembly") != 0) return;

  ufp = uop->data;
  while (ufp != NULL) 
  {
    ufp_next = ufp->next;
    if (ufp->choice == 11)
    {
      if (ufp_prev == NULL)
      {
        uop->data = ufp_next;
      }
      else
      {
        ufp_prev->next = ufp_next;
      }
      ufp->next = NULL;
      ufp = UserFieldFree (ufp);
    }
    else 
    {
      ufp_prev = ufp;
    }
    ufp = ufp_next;
  }
}

/* Need to get list of accessions from file, determine which ones are in the table,
 * and determine if any already have TPA accessions.
 * If any already have TPA accessions, ask if this should add or replace.
 */
extern void LoadTPAAccessionNumbersFromFile (
  IteM i
)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;
  Char          path [PATH_MAX];
  FILE          *fp;

  ValNodePtr    tpalist, vnp, vnp_acc;
  Boolean       replace = FALSE;
  UserObjectPtr uop;
  SeqDescrPtr   sdp;
  TPAFromFilePtr tpfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;
  
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  tpalist = ReadTPAAccessionNumbersFile (fp, sep, &replace);
  FileClose (fp);

  if (tpalist != NULL) 
  {
    for (vnp = tpalist; vnp != NULL; vnp = vnp->next)
    {
      tpfp = (TPAFromFilePtr) vnp->data.ptrvalue;
      uop = GetTPAAssembly(tpfp->bsp);
      if (uop == NULL) 
      {
        sdp = CreateNewDescriptor (SeqMgrGetSeqEntryForData (tpfp->bsp), Seq_descr_user);
        uop = CreateTpaAssemblyUserObject ();
        sdp->data.ptrvalue = uop;
      }

      if (replace)
      {
        RemoveOldTPAAccessions (uop);
      }
      for (vnp_acc = tpfp->acc_list; vnp_acc != NULL; vnp_acc = vnp_acc->next)
      {
        AddAccessionToTpaAssemblyUserObject (uop, vnp_acc->data.ptrvalue, 0, 0);
      }
    }
    tpalist = TPAFromFileListFree (tpalist);

    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    ArrowCursor ();
    Update ();
  }
}

static void AddHistory (
  BioseqPtr  bsp,
  ValNodePtr acc_list
)
{
  SeqHistPtr      hist;
  ValNodePtr      vnp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  Uint4           whichdb;
  Char            prefix [20];

  if (bsp == NULL || acc_list == NULL) return;
  hist = bsp->hist;
  if (hist == NULL)
  {
    hist = SeqHistNew ();
    if (hist == NULL) return;
    bsp->hist = hist;
  }
  for (vnp = acc_list; vnp != NULL; vnp = vnp->next) {
    tsip = TextSeqIdNew ();
    if (tsip == NULL) return;
    tsip->accession = StringSave (vnp->data.ptrvalue);

    sip = ValNodeNew (hist->replace_ids);
    if (hist->replace_ids == NULL) {
      hist->replace_ids = sip;
    }
    if (sip == NULL) return;

    sip->data.ptrvalue = (Pointer) tsip;

    StringNCpy_0 (prefix, (CharPtr) vnp->data.ptrvalue, sizeof (prefix));
    whichdb = WHICH_db_accession (prefix);
    if (ACCN_IS_EMBL (whichdb)) {
      sip->choice = SEQID_EMBL;
    } else if (ACCN_IS_DDBJ (whichdb)) {
      sip->choice = SEQID_DDBJ;
    } else {
      sip->choice = SEQID_GENBANK;
    }
  }
  if (hist != NULL
    && hist->assembly == NULL
    && hist->replace_date == NULL
    && hist->replace_ids == NULL
    && hist->replaced_by_date == NULL
    && hist->replaced_by_ids == NULL
    && hist->deleted_date == NULL
    && ! hist->deleted)
  {
      bsp->hist = SeqHistFree (bsp->hist);
  }
}

static Boolean DoIDsMatch (CharPtr seqid, BioseqPtr bsp, Boolean AllowLocal)
{
  CharPtr      str;
  Int4         seqid_len;
  SeqIdPtr     sip;
  Boolean      rval = FALSE;

  if (bsp == NULL) return FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next)
  {
    if (sip->choice != SEQID_LOCAL || AllowLocal)
    {
      str = SeqIdWholeLabel (sip, PRINTID_REPORT);
      seqid_len = StringCSpn (str, ".");
      if (seqid_len > 0)
      {
        str [ seqid_len ] = 0;
      }
      if (StringCmp (str, seqid) == 0) rval = TRUE;
      str = MemFree (str);
    }
  }
  return rval;
}

static Boolean AddAccessionToGenbankBlock (
  CharPtr     seqid,
  ValNodePtr  acc_list,
  SeqEntryPtr sep,
  Boolean     add_hist
)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  GBBlockPtr   gbp;
  ValNodePtr   last_one;
  SeqDescrPtr       sdp;

  if (seqid == NULL || acc_list == NULL
    || sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    /* this also delves into nuc-prot sets */
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)) ||
                         bssp->_class == 1)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        if (AddAccessionToGenbankBlock (seqid, acc_list, sep, add_hist))
        {
          return TRUE;
        }
      }
      return FALSE;
    }
  }
  if (!IS_Bioseq (sep)) return FALSE;

  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return FALSE;
  if (! DoIDsMatch (seqid, bsp, TRUE)) return FALSE;

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);

  if (sdp == NULL)
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_genbank);
    if (sdp == NULL) return FALSE;
  }
 
  if (sdp->data.ptrvalue == NULL)
  {
    sdp->data.ptrvalue = GBBlockNew ();
    if (sdp->data.ptrvalue == NULL) return FALSE;
  }
 
  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  
  for (last_one = gbp->extra_accessions;
       last_one != NULL && last_one->next != NULL;
       last_one = last_one->next)
  {}
  if (last_one == NULL)
  {
    gbp->extra_accessions = acc_list;
  }
  else
  {
    last_one->next = acc_list;
  }
  if (add_hist)
  {
    AddHistory (bsp, acc_list);
  }
  return TRUE;
}
 
static void LoadSecondaryAccessionNumbersPlusHistFromFile (
  IteM    i,
  Boolean add_hist
)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;
  Char          path [PATH_MAX];
  FILE          *fp;
  Char          str [8192];
  size_t        len = 8192;
  Boolean       need_seqid;
  Char          seqid[100];
  Int4          seqid_len;
  CharPtr       cp;
  CharPtr       acc_end;
  Boolean       found_end;
  ValNodePtr    acc_list;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;
  
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  need_seqid = TRUE;
  acc_list = NULL;
  ReadLine (fp, str, len);
  while (Nlm_fileDone || str[0] != 0) 
  {
    cp = str;
    if (strlen (str) == 0)
    {
      ReadLine (fp, str, len);
      continue;
    }
    seqid_len = StringCSpn (str, " \t");
    if (seqid_len > 0)
    {
      StringNCpy (seqid, str, seqid_len);
      seqid [seqid_len] = 0;
      cp = str + seqid_len + 1;
    }
    else
    {
      ReadLine (fp, str, len);
      continue;
    }
    
    found_end = FALSE;
    while (*cp != 0)
    {
      if (*cp == ' ' || *cp == ' ')
      {
        cp++;
      }
      else
      {
        acc_end = cp + 1;
        while (*acc_end != 0 && *acc_end != ' ')
        {
          acc_end++;
        }
        if (*acc_end == 0)
        {
          found_end = TRUE;
        }
        else
        {
          *acc_end = 0;
        }
        ValNodeAddStr (&acc_list, 0, StringSave (cp));
        if (found_end)
        {
          cp = acc_end;
        }
        else
        {
          cp = acc_end + 1;
        }
      }
    }

    /* do something with accession list */
    if ( ! AddAccessionToGenbankBlock (seqid, acc_list, sep, add_hist))
    {
      Message (MSG_ERROR,
               "Unable to update accession numbers for %s (not found)",
               seqid);
    }
    acc_list = NULL;
      
    ReadLine (fp, str, len);
  }
  if (acc_list != NULL
    && ! AddAccessionToGenbankBlock (seqid, acc_list, sep, add_hist))
  {
    Message (MSG_ERROR,
             "Unable to update accession numbers for %s (not found)",
             seqid);
  }

  FileClose (fp);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  return;
}

extern void LoadSecondaryAccessionNumbersFromFile (
  IteM i
)
{
  LoadSecondaryAccessionNumbersPlusHistFromFile (i, FALSE);
}

extern void LoadHistoryAccessionNumbersFromFile (
  IteM i
)
{
  LoadSecondaryAccessionNumbersPlusHistFromFile (i, TRUE);
}


static void RemoveEmptyGenomeProjectIDCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovn;
  UserObjectPtr uop;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp;
  
  if (sdp == NULL || sdp->extended == 0 || sdp->choice != Seq_descr_user) return;
    
  if (sdp->data.ptrvalue != NULL) {       
    uop = (UserObjectPtr) sdp->data.ptrvalue;    
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "GenomeProjectsDB") == 0) {
      ufp = uop->data;
      while (ufp != NULL) {
        oip = ufp->label;
        if (oip != NULL
            && StringCmp (oip->str, "ProjectID") == 0 || StringCmp (oip->str, "ParentID") == 0) {
          if(ufp->choice == 2 && ufp->data.intvalue != 0) {
            /* found nonempty ID */
            return;
          }          
        } else if (oip != NULL) {
          /* found nonempty field of other type */
          return;
        }
        ufp = ufp->next;
      }
      /* found no nonempty fields */
      ovn = (ObjValNodePtr) sdp;
      ovn->idx.deleteme = TRUE;
    }
  }      
}


extern void RemoveEmptyGenomeProjectIDs (IteM i)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  VisitDescriptorsInSep (sep, NULL, RemoveEmptyGenomeProjectIDCallback);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  return;
}


static void RemoveAllGenomeProjectIDCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjValNodePtr ovn;
  UserObjectPtr uop;
  ObjectIdPtr   oip;
  
  if (sdp == NULL || sdp->extended == 0 || sdp->choice != Seq_descr_user) return;
    
  if (sdp->data.ptrvalue != NULL) {       
    uop = (UserObjectPtr) sdp->data.ptrvalue;    
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "GenomeProjectsDB") == 0) {
      ovn = (ObjValNodePtr) sdp;
      ovn->idx.deleteme = TRUE;
    }
  }      
}


extern void RemoveGenomeProjectIdsBaseForm (BaseFormPtr bfp)
{
  SeqEntryPtr   sep;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  
  VisitDescriptorsInSep (sep, NULL, RemoveAllGenomeProjectIDCallback);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  return;
}

extern void RemoveGenomeProjectIDs  (IteM i)
{
  BaseFormPtr   bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  RemoveGenomeProjectIdsBaseForm (bfp);
}


CharPtr MostUsedFeatureList[] = { 
  "CDS",
  "exon",
  "Gene",
  "intron",
  "mRNA",
  "rRNA",
  "RNA"
};

extern ValNodePtr InsertMostUsedFeatureValNodes (ValNodePtr old_list)
{
  ValNodePtr new_list, new_item, old_item;
  Int4       index;

  new_list = NULL;
  for (index = 0;
       index < sizeof (MostUsedFeatureList) / sizeof (CharPtr);
       index ++)
  {
    old_item = FindExactStringListMatch ( old_list, MostUsedFeatureList [index])
;
    if (old_item == NULL) continue;
    new_item = ValNodeNew ( new_list);
    if (new_item == NULL) return old_list;
    new_item->choice = old_item->choice;
    new_item->data.ptrvalue = StringSave (MostUsedFeatureList [index]);
    if (new_list == NULL) new_list = new_item;
  }
  if (new_item != NULL)
  {
    if (old_list != NULL &&
      ( StringCmp (old_list->data.ptrvalue, "All") == 0
       || StringCmp (old_list->data.ptrvalue, "[ALL FEATURES]") == 0))
    {
      new_item->next = old_list->next;
      old_list->next = new_list;
      new_list = old_list;
    }
    else
    {
      new_item->next = old_list;
    }
  }
  else
  {
    new_list = old_list;
  }
  return new_list;
}

static EnumFieldAssocPtr FindEnumFieldAssoc (
  EnumFieldAssocPtr alist,
  CharPtr findStr
)
{
  EnumFieldAssocPtr ap;
  
  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    if (StringCmp (ap->name, findStr) == 0) return ap;
  }
  return NULL;
}

static void CopyEnumFieldAssoc (EnumFieldAssocPtr ap1, EnumFieldAssocPtr ap2)
{
  if (ap1 == NULL || ap2 == NULL) return;

  ap1->name = StringSave (ap2->name);
  ap1->value = ap2->value;
}

extern EnumFieldAssocPtr InsertMostUsedFeatureEnumFieldAssoc (
  EnumFieldAssocPtr alist
)
{
  Int4              num_total_fields, index, new_index;
  EnumFieldAssocPtr ap, new_alist, old_ap;

  num_total_fields = sizeof (MostUsedFeatureList) / sizeof (CharPtr);

  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    num_total_fields ++;
  }
  /* need the last null field */
  num_total_fields ++;

  new_alist = MemNew (num_total_fields * sizeof (EnumFieldAssoc));
  if (new_alist == NULL) return alist;

  /* copy the first item if wildcard */
  if (StringCmp (alist->name, "[ALL FEATURES]") == 0)
  {
    CopyEnumFieldAssoc (new_alist, alist);
    new_index = 1;
  }
  else
  {
    new_index = 0;
  }

  for (index = 0;
       index < sizeof (MostUsedFeatureList) / sizeof (CharPtr);
       index ++)
  {
    old_ap = FindEnumFieldAssoc (alist, MostUsedFeatureList [index]);
    if (old_ap == NULL) continue;
    CopyEnumFieldAssoc (new_alist + new_index++, old_ap);
  }

  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    CopyEnumFieldAssoc (new_alist + new_index ++, ap);
  }
  /* copy over the last null field */
  if (ap != NULL)
  {
    CopyEnumFieldAssoc (new_alist + new_index ++, ap);
  }
  return new_alist;
  
}

static Uint2 UnusualFeatureTypes [] = {
  FEATDEF_ORG,
  FEATDEF_mutation,
  FEATDEF_site_ref,
  FEATDEF_gap,
  FEATDEF_NON_STD_RESIDUE,
  FEATDEF_NUM
};
 
extern ValNodePtr BuildFeatureValNodeList (
  Boolean prefer_most_used,
  CharPtr wild_card_name,
  Int4    wild_card_value,
  Boolean skip_unusual,
  Boolean skip_import
)
{
  FeatDefPtr  curr;
  ValNodePtr  head, vnp;
  Uint1       key;
  CharPtr     label = NULL;
  Uint1       subtype;
  Int4        index;
  Boolean     skip;
  Char        str [256];

  head = NULL;
  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    skip = FALSE;
    if (skip_unusual)
    {
      for (index = 0;
           ! skip && index < sizeof ( UnusualFeatureTypes ) / sizeof (Uint2);
           index ++)
      {
        if (curr->featdef_key == UnusualFeatureTypes [ index ]) skip = TRUE;
      }
    }
    if (key != FEATDEF_BAD && ! skip && !IsUnwantedFeatureType(key)) {
      
      subtype = curr->featdef_key;
	  if (subtype == FEATDEF_PUB)
	  {
        StringNCpy_0 (str, curr->typelabel, sizeof (str) - 15);
        StringCat (str, " (Publication)");
	  }
	  else if (subtype != FEATDEF_misc_RNA &&
          subtype != FEATDEF_precursor_RNA &&
          subtype != FEATDEF_mat_peptide &&
          subtype != FEATDEF_sig_peptide &&
          subtype != FEATDEF_transit_peptide &&
          subtype != FEATDEF_Imp_CDS)
      {
        StringNCpy_0 (str, curr->typelabel, sizeof (str) - 1);
      }
      else if (! skip_import)
      {
        StringNCpy_0 (str, curr->typelabel, sizeof (str) - 10);
        StringCat (str, "_imp");
      }
      else
      {
        skip = TRUE;
      }
      if (! skip)
      {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          vnp->data.ptrvalue = StringSave (str);
        }
      }
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  if (head != NULL) {
    head = SortValNode (head, CompareFeatureValNodeStrings);
    head = InsertMostUsedFeatureValNodes (head);
    if (wild_card_name != NULL)
    {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->choice = wild_card_value;
        vnp->data.ptrvalue = StringSave (wild_card_name);
        vnp->next = head;
        head = vnp;
      }
    }
  }
  return head;
}

extern void RemoveOldName (OrgRefPtr orp)
{
  OrgModPtr prev = NULL, curr, next_mod;
  
  if (orp == NULL || orp->orgname == NULL) return;
  
  curr = orp->orgname->mod;
  while (curr != NULL)
  {
    next_mod = curr->next;
    if (curr->subtype == ORGMOD_old_name)
    {
      if (prev == NULL)
      {
        orp->orgname->mod = curr->next;
      }
      else
      {
        prev->next = curr->next;
      }
      curr->next = NULL;
      OrgModFree (curr);
    }
    else
    {
      prev = curr;
    }
    
    curr = next_mod;
  }
  
}

extern void SetTaxNameAndRemoveTaxRef (OrgRefPtr orp, CharPtr taxname)
{
  Boolean         remove_taxrefs = FALSE;

  if (orp == NULL) return;

  if ( taxname == NULL || orp->taxname == NULL
    || StringCmp (taxname, orp->taxname) != 0)
  {
    remove_taxrefs = TRUE;
  }
  MemFree (orp->taxname);
  orp->taxname = taxname;

  if (! remove_taxrefs) return;

  orp->common = MemFree (orp->common);

  RemoveTaxRef (orp);
}

static Boolean
FindMatchingProprotein 
(SeqFeatPtr sfp,
 SeqMgrFeatContextPtr fcontext,
 BioseqPtr prot_bsp)
{
  SeqFeatPtr        prot_sfp;
  SeqMgrFeatContext pcontext;
  CharPtr           start;

  if (prot_bsp == NULL || fcontext == NULL) return FALSE;
  if (StringNICmp (fcontext->label, "encodes ", 8) == 0) {
    start = fcontext->label + 8;
  } else {
    start = fcontext->label;
  }
  prot_sfp = NULL;
  while ((prot_sfp = SeqMgrGetNextFeature (prot_bsp, prot_sfp, 
                                           0, 0, &pcontext)) != NULL) {
    if (StringCmp (pcontext.label, start) == 0) {
      return TRUE;
    } 
  }
  return FALSE;
}


static void 
RemoveRedundantProproteinMiscFeatsOnBioseq
(BioseqPtr bsp,
 Pointer userdata)
{
  SeqFeatPtr        sfp, cds;
  SeqMgrFeatContext fcontext, cds_context;
  BioseqPtr         bsp_prot;

  sfp = NULL;

  /* list misc feats */
  while ((sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext)) != NULL) {
    if (fcontext.featdeftype == FEATDEF_misc_feature
        &&  StringStr(fcontext.label, "proprotein") != NULL) {
      cds = NULL;
      while ((cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &cds_context)) != NULL) {
        if (cds_context.left <= fcontext.left
            &&  cds_context.right >= fcontext.right) {
          /* Get Protein sequence, look for matching proprotein feat */
          bsp_prot = BioseqFind (SeqLocId(cds->product));
          if (FindMatchingProprotein (sfp, &fcontext, bsp_prot)) {
            sfp->idx.deleteme = TRUE;
          }
        }
      }
    }
  }

}


extern void RemoveRedundantProproteinMiscFeats (IteM i)
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

  /* Visit each bioseq to remove redundant proprotein misc feats */
  VisitBioseqsInSep (sep, NULL, RemoveRedundantProproteinMiscFeatsOnBioseq);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

typedef struct typestraindata
{
  FORM_MESSAGE_BLOCK

  GrouP  strain_or_comment_grp;
  DialoG string_constraint_dlg;
  ButtoN remove_found_btn;

  Boolean             search_strain;
  StringConstraintPtr string_constraint;
  Boolean             remove_found_text;
} TypeStrainData, PNTR TypeStrainPtr;


static void AddTypeStrainCommentsProc (BioSourcePtr biop, Pointer userdata)
{
  TypeStrainPtr      tsp;
  CharPtr            orig_note;
  CharPtr            format = "type strain of %s";
  CharPtr            match = NULL, orig_match, tmp, taxname;
  ValNode            vn, vn_comment, vn_taxname;

  if (biop == NULL || biop->org == NULL || biop->org->taxname == NULL) return;

  tsp = (TypeStrainPtr) userdata;

  vn_comment.choice = SourceQualChoice_textqual;
  vn_comment.data.intvalue = Source_qual_orgmod_note;
  vn_comment.next = NULL;
  vn_taxname.choice = SourceQualChoice_textqual;
  vn_taxname.data.intvalue = Source_qual_taxname;
  vn_taxname.next = NULL;

  orig_note = GetSourceQualFromBioSource (biop, &vn, NULL);
  if (orig_note != NULL && StringStr (orig_note, "type strain of") != NULL) {
    orig_note = MemFree (orig_note);
    return;
  }
  orig_note = MemFree (orig_note);

  taxname = GetSourceQualFromBioSource (biop, &vn_taxname, NULL);
  
  if (tsp != NULL && !IsStringConstraintEmpty (tsp->string_constraint)) {
    if (tsp->search_strain) {
      vn.choice = SourceQualChoice_textqual;
      vn.data.intvalue = Source_qual_strain;
      vn.next = NULL;  
      match = GetSourceQualFromBioSource (biop, &vn, NULL);
    } else {
      vn.choice = SourceQualChoice_textqual;
      vn.data.intvalue = Source_qual_orgmod_note;
      vn.next = NULL;  
      match = GetSourceQualFromBioSource (biop, &vn, NULL);
      if (match == NULL) {
        vn.choice = Source_qual_subsource_note;
        match = GetSourceQualFromBioSource (biop, &vn, NULL);
      }
    }
    if (!DoesStringMatchConstraint (match, tsp->string_constraint)) {
      match = MemFree (match);
      taxname = MemFree (taxname);
      return;
    }

    orig_match = StringSave (match);
    if (tsp->remove_found_text && RemoveStringConstraintPortionFromString (&match, tsp->string_constraint)) {
      if (StringHasNoText (match)) {
        RemoveSourceQualFromBioSource (biop, &vn, tsp->string_constraint);
      } else {
        TrimSpacesAroundString (match);
        SetSourceQualInBioSource (biop, &vn, tsp->string_constraint, match, ExistingTextOption_replace_old);
      }
      FindReplaceString (&taxname, orig_match, match, TRUE, TRUE);
      SetSourceQualInBioSource (biop, &vn_taxname, NULL, taxname, ExistingTextOption_replace_old);
    }
    orig_match = MemFree (orig_match);
    match = MemFree (match);
  }

  vn.choice = SourceQualChoice_textqual;
  vn.data.intvalue = Source_qual_strain;
  vn.next = NULL;  

  if (!StringHasNoText (taxname)) {
    tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (format) + StringLen (taxname)));
    sprintf (tmp, format, taxname);
    SetSourceQualInBioSource (biop, &vn_comment, NULL, tmp, ExistingTextOption_append_semi);
    tmp = MemFree (tmp);
  }
  taxname = MemFree (taxname);
}


static void AddTypeStrainCommentsWithConstraintProc (ButtoN b)
{
  TypeStrainPtr tsp;
  SeqEntryPtr   sep;
  
  tsp = (TypeStrainPtr) GetObjectExtra (b);
  if (tsp == NULL) return;
  sep = GetTopSeqEntryForEntityID (tsp->input_entityID);
  if (sep == NULL) return;

  tsp->string_constraint = DialogToPointer (tsp->string_constraint_dlg);
  tsp->remove_found_text = GetStatus (tsp->remove_found_btn);
  if (GetValue (tsp->strain_or_comment_grp) == 1) {
    tsp->search_strain = TRUE;
  } else {
    tsp->search_strain = FALSE;
  }
  
  /* Visit each bioseq to remove redundant proprotein misc feats */
  VisitBioSourcesInSep (sep, tsp, AddTypeStrainCommentsProc);

  tsp->string_constraint = StringConstraintFree (tsp->string_constraint);

  ObjMgrSetDirtyFlag (tsp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, tsp->input_entityID, 0, 0);
  Remove (tsp->form);
  ArrowCursor ();
  Update ();
}

extern void AddTypeStrainCommentsWithConstraint (IteM i)
{
  BaseFormPtr    bfp;
  TypeStrainPtr  tsp;
  WindoW         w;
  GrouP          h, c;
  PrompT         p;
  ButtoN         b;
  StringConstraintPtr scp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
	
  tsp = (TypeStrainPtr) MemNew (sizeof (TypeStrainData));
  if (tsp == NULL) return;
  tsp->input_entityID = bfp->input_entityID;

  w = FixedWindow (-50, -33, -10, -10, "Add Type Strain Comments", StdCloseWindowProc);
  if (w == NULL) {
	MemFree (tsp);
	return;
  }
  tsp->form = (ForM) w;
  SetObjectExtra (w, tsp, StdCleanupFormProc);
  
  h = HiddenGroup (w, -1, 0, NULL);

  p = StaticPrompt (h, "When", 0, dialogTextHeight, systemFont, 'c');
  tsp->strain_or_comment_grp = HiddenGroup (h, 2, 0, NULL);
  RadioButton (tsp->strain_or_comment_grp, "Strain");
  RadioButton (tsp->strain_or_comment_grp, "Comment");
  SetValue (tsp->strain_or_comment_grp, 1);

  tsp->string_constraint_dlg = StringConstraintDialog (h, NULL, FALSE, NULL, NULL);
  scp = StringConstraintNew ();
  scp->match_location = String_location_ends;
  scp->case_sensitive = TRUE;
  PointerToDialog (tsp->string_constraint_dlg, scp);
  scp = StringConstraintFree (scp);

  tsp->remove_found_btn = CheckBox (h, "Remove found text", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", AddTypeStrainCommentsWithConstraintProc);
  SetObjectExtra (b, tsp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, tsp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) tsp->strain_or_comment_grp, (HANDLE) tsp->string_constraint_dlg, 
                (HANDLE) tsp->remove_found_btn, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

extern void AddTypeStrainCommentsToAll (IteM i)
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

  /* Visit each bioseq to remove redundant proprotein misc feats */
  VisitBioSourcesInSep (sep, NULL, AddTypeStrainCommentsProc);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

extern void SqnNewAlign (BioseqPtr bsp1, BioseqPtr bsp2, SeqAlignPtr PNTR salp)
{
  BLAST_SummaryOptions *options = NULL;
  Uint1 mol_was;

  if (bsp1 == NULL || bsp2 == NULL || salp == NULL) return;

  *salp = NULL;
  if (ISA_na (bsp1->mol) != ISA_na (bsp2->mol)) return;

  mol_was = bsp2->mol;
  bsp2->mol = bsp1->mol;
  BLAST_SummaryOptionsInit(&options);

  options->cutoff_evalue = 0.001;
  if (bsp1->length > 10000 || bsp2->length > 10000)
  {
    options->filter_string = StringSave ("m L");
    options->word_size = 20;
    options->cutoff_evalue = act_get_eval (60);
    if (ISA_na (bsp1->mol))
    {
      options->program = eBlastn;
    }
    else
    {
      options->program = eBlastp;
    }
    options->hint = eNone;
  }

  BLAST_TwoSequencesSearch(options, bsp1, bsp2, salp);
  bsp2->mol = mol_was;
  BLAST_SummaryOptionsFree(options);
  
}

/* This section of code is for the Remove Sequences From Alignments function. */

typedef struct alignmentsequencelist {
  SeqIdPtr sip;
  Char     descr[255];
} AlignmentSequenceListData, PNTR AlignmentSequenceListPtr;

static void 
ListSequencesInSeqEntry 
(SeqEntryPtr sep,
 ValNodePtr PNTR list, 
 Boolean show_nucs, 
 Boolean show_prots)
{
  BioseqPtr                bsp;
  BioseqSetPtr             bssp;
  ValNodePtr               vnp;
  AlignmentSequenceListPtr aslp;
  Int4                     offset;
  SeqIdPtr                 bsp_sip;
  
  if (sep == NULL) return;
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    if (!show_nucs && ISA_na (bsp->mol))
    {
      return;
    }
    if (!show_prots && ISA_aa (bsp->mol))
    {
      return;
    }
    aslp = (AlignmentSequenceListPtr) MemNew (sizeof (AlignmentSequenceListData));
    if (aslp == NULL) return;
    aslp->sip = bsp->id;
    aslp->descr[0] = 0;
	  aslp->descr[253] = 0;
    offset = 0;
    for (bsp_sip = bsp->id; bsp_sip != NULL && offset < 250; bsp_sip = bsp_sip->next) {
	  if (aslp->descr[0] != 0) {
	    aslp->descr[offset] = ':';
	    offset ++;
	  }
      SeqIdWrite (bsp_sip, aslp->descr + offset, PRINTID_TEXTID_ACCESSION, 254 - offset);
      offset = StringLen (aslp->descr);
	}
    vnp = ValNodeNew (*list);
    if (vnp != NULL)
    {
      vnp->data.ptrvalue = aslp;
    }
    if (*list == NULL)
    {
      *list = vnp;
    }
  }
  else
  {
  	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
    {
      ListSequencesInSeqEntry (sep, list, show_nucs, show_prots);
    }
  }
}

typedef struct sequencelistctrl
{
  ValNodePtr      sequence_list;
  Nlm_LstActnProc actn;
  Pointer         userdata;
  
} SequenceListCtrlData, PNTR SequenceListCtrlPtr;

static void CleanupSequenceListCtrl (
  GraphiC g,
  VoidPtr data
)

{
  SequenceListCtrlPtr slcp;

  slcp = (SequenceListCtrlPtr) data;
  if (slcp != NULL) {
	  slcp->sequence_list = ValNodeFreeData (slcp->sequence_list);
  }
  MemFree (slcp);
}


static void SequenceListCtrlAction (LisT l)
{
  SequenceListCtrlPtr slcp;
  
  slcp = (SequenceListCtrlPtr) GetObjectExtra (l);
  if (slcp == NULL) return;
  
  if (slcp->actn != NULL)
  {
    SetObjectExtra (l, slcp->userdata, NULL);
    (slcp->actn) (l);
    SetObjectExtra (l, slcp, CleanupSequenceListCtrl);
  }
}

extern LisT 
MakeSequenceListControl 
(GrouP g,
 SeqEntryPtr sep,
 Nlm_LstActnProc actn,
 Pointer userdata,
 Boolean show_nucs,
 Boolean show_prots)
{
  LisT                     list_ctrl;
  SequenceListCtrlPtr      slcp;
  ValNodePtr               vnp;
  AlignmentSequenceListPtr aslp;
  
  slcp = (SequenceListCtrlPtr) MemNew (sizeof (SequenceListCtrlData));
  slcp->actn = actn;
  slcp->userdata = userdata;
  ListSequencesInSeqEntry (sep, &slcp->sequence_list, show_nucs, show_prots);
  
  list_ctrl = MultiList (g, 20, 8, SequenceListCtrlAction);
  SetObjectExtra (list_ctrl, slcp, CleanupSequenceListCtrl);
  
  for (vnp = slcp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	  if (aslp != NULL) 
	  {
      ListItem (list_ctrl, aslp->descr);
	  }
  }

  return list_ctrl;
  
}


extern void SelectAllSequencesInListCtrl (LisT l)
{
  SequenceListCtrlPtr   slcp;
  ValNodePtr            vnp;
  Int2                  val;
  
  
  slcp = (SequenceListCtrlPtr) GetObjectExtra (l);
  if (slcp == NULL) return;
  
  for (val = 1, vnp = slcp->sequence_list; vnp != NULL; vnp = vnp->next, val++)
  {
    SetItemStatus (l, val, TRUE);
  }  
}


extern void UnSelectAllSequencesInListCtrl (LisT l)
{
  SequenceListCtrlPtr   slcp;
  ValNodePtr            vnp;
  Int2                  val;
  
  
  slcp = (SequenceListCtrlPtr) GetObjectExtra (l);
  if (slcp == NULL) return;
  
  for (val = 1, vnp = slcp->sequence_list; vnp != NULL; vnp = vnp->next, val++)
  {
    SetItemStatus (l, val, FALSE);
  }  
}


extern ValNodePtr GetSelectedSequenceList (LisT l)
{
  SequenceListCtrlPtr      slcp;
  ValNodePtr               sip_list = NULL, vnp;
  Int2                     val;
  AlignmentSequenceListPtr aslp;

  slcp = (SequenceListCtrlPtr) GetObjectExtra (l);
  if (slcp == NULL) return NULL;
  
  val = 1;
  for (vnp = slcp->sequence_list; vnp != NULL; vnp = vnp->next) 
  {
    aslp = vnp->data.ptrvalue;
	  if (aslp == NULL) continue;
	  if (GetItemStatus (l, val)) 
	  {
	    ValNodeAddPointer (&sip_list, 0, aslp->sip);
	  }
	  val++;
  }
  
  return sip_list;
}

/* This function is used so that a sequence ID will only appear once in the list,
 * even if it appears in more than one alignment or subalignment.
 */
static Boolean IsIDAlreadyInList (SeqIdPtr sip, ValNodePtr list)
{
  ValNodePtr vnp;
  AlignmentSequenceListPtr aslp;
  
  if (sip == NULL) return FALSE;
  
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    aslp = (AlignmentSequenceListPtr) vnp->data.ptrvalue;
    if (aslp != NULL && SeqIdComp (aslp->sip, sip) == SIC_YES)
    {
      return TRUE;
    }
  }
  return FALSE;
}

/* This function creates the list of sequence IDs and descriptions to use in 
 * the Remove Sequences From Alignments dialog.
 */
static void ListSequencesInAlignmentsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip_list, sip, bsp_sip;
  ValNodePtr PNTR list;
  ValNodePtr  vnp; 
  AlignmentSequenceListPtr aslp;
  BioseqPtr                bsp;
  Int4                     offset;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  while (salp != NULL) 
  {
    list = (ValNodePtr PNTR)userdata;
    sip_list = SeqAlignIDList (salp);
    if (sip_list == NULL) return;
    for (sip = sip_list; sip != NULL; sip = sip->next) {
      if (IsIDAlreadyInList (sip, *list)) continue;
      aslp = (AlignmentSequenceListPtr) MemNew (sizeof (AlignmentSequenceListData));
	  if (aslp == NULL) return;
	  aslp->sip = sip;
	  bsp = BioseqFindCore (sip);
	  if (bsp != NULL) {
		  aslp->descr[0] = 0;
		  aslp->descr[253] = 0;
		  offset = 0;
		  for (bsp_sip = bsp->id; bsp_sip != NULL && offset < 250; bsp_sip = bsp_sip->next) {
			if (aslp->descr[0] != 0) {
			  aslp->descr[offset] = '\t';
			  offset ++;
			}
		    SeqIdWrite (bsp_sip, aslp->descr + offset, PRINTID_TEXTID_ACCESSION, 254 - offset);
			offset = StringLen (aslp->descr);
		  }
	  } else {
        SeqIdWrite (sip, aslp->descr, PRINTID_TEXTID_ACCESSION, 254);	    
	  }
	  vnp = ValNodeNew (*list);
	  vnp->data.ptrvalue = aslp;
	  if (*list == NULL) {
		  *list = vnp;
	  }	  
    }
    salp = salp->next;
  }
}

static ValNodePtr ListSequencesInAlignments (SeqEntryPtr sep)
{
	ValNodePtr list = NULL;
    VisitAnnotsInSep (sep, (Pointer) &list, ListSequencesInAlignmentsCallback);
    return list;
}

static LisT MakeAlignmentSequenceListControl (GrouP g, SeqEntryPtr sep, Nlm_LstActnProc actn, Pointer userdata)
{
  LisT                     list_ctrl;
  SequenceListCtrlPtr      slcp;
  ValNodePtr               vnp;
  AlignmentSequenceListPtr aslp;
  
  slcp = (SequenceListCtrlPtr) MemNew (sizeof (SequenceListCtrlData));
  slcp->actn = actn;
  slcp->userdata = userdata;
  slcp->sequence_list = ListSequencesInAlignments (sep);
  
  list_ctrl = MultiList (g, 16, 16, SequenceListCtrlAction);
  SetObjectExtra (list_ctrl, slcp, CleanupSequenceListCtrl);
  
  for (vnp = slcp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	  if (aslp != NULL) 
	  {
      ListItem (list_ctrl, aslp->descr);
	  }
  }

  return list_ctrl;
  
}

typedef struct removeseqfromaligndata {
  FORM_MESSAGE_BLOCK
  DialoG      clickable_list_dlg;
  DialoG      constraint_dlg;
  ValNodePtr  seq_list;
  LisT        sequence_list_ctrl;
  TexT        mark_seqid_list;
  TexT        mark_seqlen;

  SeqEntryPtr sep;
  Boolean     remove_all_from_alignments;
  Boolean     no_remove_all_from_alignments;
  Boolean     remove_all_products;
  Boolean     no_remove_all_products;
} RemoveSeqFromAlignData, PNTR RemoveSeqFromAlignPtr;

static void DoRemoveSequencesFromAlignment (ButtoN b)
{
  RemoveSeqFromAlignPtr    rp;
  WindoW                   w;
  ValNodePtr               vnp, sip_list;
  
  if (b == NULL) return;
  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  
  w = (WindoW) rp->form;
  Hide (w);
  /* first, check for pairwise alignments */
  sip_list = GetSelectedSequenceList (rp->sequence_list_ctrl);
  for (vnp = sip_list; vnp != NULL; vnp = vnp->next)
  {
    if (IsSequenceFirstInPairwise (rp->sep, (SeqIdPtr) vnp->data.ptrvalue))
	  {
	  	Message (MSG_ERROR, "One of the selected sequences is the first in a pairwise alignment."
	  	"  You must convert the alignment to a multiple alignment before trying to remove this sequence.");
      Remove (rp->form);  
      ValNodeFree (sip_list);
      return;
	  }
  }

  for (vnp = sip_list; vnp != NULL; vnp = vnp->next)
  {
    RemoveSequenceFromAlignments (rp->sep, (SeqIdPtr) vnp->data.ptrvalue);
  }
 
  ValNodeFree (sip_list);
  DeleteMarkedObjects (rp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (rp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rp->input_entityID, 0, 0);
  Remove (rp->form);  
}


extern void RemoveSequencesFromAlignment (IteM i)
{
  BaseFormPtr              bfp;
  WindoW                   w;
  RemoveSeqFromAlignPtr    rp;
  GrouP                    h, k, c;
  ButtoN                   b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;

  rp = (RemoveSeqFromAlignPtr) MemNew (sizeof (RemoveSeqFromAlignData));
  if (rp == NULL) return;
  rp->input_entityID = bfp->input_entityID;
  rp->sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (rp->sep == NULL) 
  {
	  MemFree (rp);
	  return;
  }

  w = FixedWindow (-50, -33, -10, -10, "Remove Sequences From Alignment", StdCloseWindowProc);
  if (w == NULL) {
	MemFree (rp);
	return;
  }
  rp->form = (ForM) w;
  SetObjectExtra (w, rp, StdCleanupFormProc);
  
  h = HiddenGroup (w, -1, 0, NULL);
  k = HiddenGroup (h, 2, 0, NULL);

  rp->sequence_list_ctrl = MakeAlignmentSequenceListControl (k, rp->sep, NULL, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveSequencesFromAlignment);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, rp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

/* End of Remove Sequences From Alignments function code. */

/* This section of code is used for removing sequences from the record. */

static void DoesBioseqHaveFeaturesWithProductsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR list;
  ValNodePtr vnp;
  
  if (sfp == NULL || userdata == NULL) return;
  list = (ValNodePtr PNTR) userdata;
  
  if (sfp->product != NULL)
  {
  	vnp = ValNodeNew (*list);
  	if (vnp != NULL)
  	{
  	  vnp->data.ptrvalue = sfp;
  	}
  	if (*list == NULL)
  	{
  	  *list = vnp;
  	}
  }
}

static void RemoveBioseq (BioseqPtr bsp, RemoveSeqFromAlignPtr rp);

static void RemoveBioseqProducts (ValNodePtr product_feature_list, RemoveSeqFromAlignPtr rp)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  BioseqPtr  bsp;
  
  for (vnp = product_feature_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
  	  bsp = BioseqFindFromSeqLoc (sfp->product);
  	  sfp->product = SeqLocFree (sfp->product);
  	  RemoveBioseq (bsp, rp);
    }
  }
}

static void RemoveNucProtSet (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;

  if (sep == NULL || !IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;

  bssp->idx.deleteme = TRUE;
}

typedef struct removealnorproductans 
{
  WindoW  w;
  Boolean ans;
  Boolean do_all;
  Boolean done;
} RemoveAlnOrProductAnsData, PNTR RemoveAlnOrProductAnsPtr;

static void RemoveAlnOrProductYes (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = TRUE;
  rp->do_all = FALSE;
  Remove (rp->w);
  rp->done = TRUE;
}

static void RemoveAlnOrProductYesAll (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = TRUE;
  rp->do_all = TRUE;
  Remove (rp->w);
  rp->done = TRUE;
}

static void RemoveAlnOrProductNo (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = FALSE;
  rp->do_all = FALSE;
  Remove (rp->w);
  rp->done = TRUE;
}

static void RemoveAlnOrProductNoAll (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = FALSE;
  rp->do_all = TRUE;
  Remove (rp->w);
  rp->done = TRUE;
}

static Boolean GetRemoveProducts (RemoveSeqFromAlignPtr rp, CharPtr idstr)
{
  RemoveAlnOrProductAnsData rd;

  GrouP                    g, h, c;
  ButtoN                   b;
  CharPtr                  prompt_fmt = "%s contains features that have products (proteins, etc.).  Would you like to remove the product sequences?";
  CharPtr                  prompt_str = NULL;
  
  if (rp == NULL || idstr == NULL) return FALSE;
  if (rp->remove_all_products) return TRUE;
  if (rp->no_remove_all_products) return FALSE;

  prompt_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (prompt_fmt) + StringLen (idstr)));
  if (prompt_str == NULL) return FALSE;
  sprintf (prompt_str, prompt_fmt, idstr);
  rd.w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup(rd.w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  rd.done = FALSE;
  g = HiddenGroup (h, 1, 0, NULL);
  StaticPrompt (g, prompt_str, 0, popupMenuHeight, programFont, 'l');
  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton(c, "Yes", RemoveAlnOrProductYes);
  SetObjectExtra (b, &rd, NULL);
  b = PushButton(c, "Remove All", RemoveAlnOrProductYesAll);
  SetObjectExtra (b, &rd, NULL);
  b = DefaultButton(c, "No", RemoveAlnOrProductNo);
  SetObjectExtra (b, &rd, NULL);
  b = DefaultButton(c, "Remove None", RemoveAlnOrProductNoAll);
  SetObjectExtra (b, &rd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  prompt_str = MemFree (prompt_str);
  
  Show(rd.w); 
  Select (rd.w);
  rd.done = FALSE;
  while (!rd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (rd.do_all)
  {
    if (rd.ans)
    {
  	  rp->remove_all_products = TRUE;
  	  rp->no_remove_all_products = FALSE;
    }
    else
    {
  	  rp->remove_all_products = FALSE;
  	  rp->no_remove_all_products = TRUE;
    }
  }
  return rd.ans;
}


static void RemoveBioseq (BioseqPtr bsp, RemoveSeqFromAlignPtr rp)
{
  ValNodePtr   product_feature_list = NULL;
  CharPtr      str = NULL;
  SeqEntryPtr  sep;
  
  if (bsp == NULL || rp == NULL) return;	
  
  if (IsBioseqInAnyAlignment (bsp, rp->input_entityID))
  {
    if (!rp->remove_all_from_alignments && !rp->no_remove_all_from_alignments) 
    {
	    if (ANS_YES == Message (MSG_YN, "This sequence is part of an alignment.  Would you like to remove it from the alignment?")) 
      {
        rp->remove_all_from_alignments = TRUE;
      }
      else
      {
        rp->no_remove_all_from_alignments = TRUE;
      }
    }

    if (rp->remove_all_from_alignments) {
      RemoveSequenceFromAlignments (rp->sep, bsp->id);
    }

  }
  VisitFeaturesOnBsp (bsp, &product_feature_list, DoesBioseqHaveFeaturesWithProductsCallback);
  if (product_feature_list != NULL)
  {
    str = SeqIdWholeLabel (bsp->id, PRINTID_REPORT);
    if (GetRemoveProducts (rp, str))
    {
      RemoveBioseqProducts (product_feature_list, rp);
    }
    str = MemFree (str);
  }
        
  bsp->idx.deleteme = TRUE;
  if (ISA_na (bsp->mol)) {
    /* remove nuc-prot set if we are deleting the nucleotide */
    sep = GetBestTopParentForData (rp->input_entityID, bsp);
    RemoveNucProtSet (sep);
  }
  ValNodeFree (product_feature_list);
  
}


static void DoRemoveSequencesFromRecord (ButtoN b)
{
  RemoveSeqFromAlignPtr    rp;
  WindoW                   w;
  ValNodePtr               vnp, item_vnp;
  BioseqPtr                bsp;
  ClickableItemPtr         cip;
  Boolean                  any_removed = FALSE;
  
  if (b == NULL) return;
  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  
  w = (WindoW) rp->form;
  Hide (w);

  for (vnp = rp->seq_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (!cip->chosen) continue;
    if (cip != NULL && cip->item_list != NULL) {
      for (item_vnp = cip->item_list; item_vnp != NULL; item_vnp = item_vnp->next) {
        if (item_vnp->choice == OBJ_BIOSEQ) {
          bsp = item_vnp->data.ptrvalue;
	        if (bsp != NULL) {
	          RemoveBioseq (bsp, rp);
            any_removed = TRUE;
	        }
        }
      }
    }
  }
  if (!any_removed) {
    Show (w);
    Message (MSG_ERROR, "You have not selected any sequences for removal!\n(Mark the sequences you want to remove by checking the box next to the sequence or by using the Mark button next to the constraint.)");
    return;
  }
  
  DeleteMarkedObjects (rp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (rp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rp->input_entityID, 0, 0);
  Remove (rp->form);  
}


static void ChooseCategoriesByMacroSequenceConstraint (ValNodePtr value_list, SequenceConstraintPtr scp, Boolean do_choose)
{
  ClickableItemPtr cip;

  while (value_list != NULL) {
    cip = (ClickableItemPtr) value_list->data.ptrvalue;
    if (cip != NULL) {
      if (cip->item_list != NULL 
          && cip->item_list->choice == OBJ_BIOSEQ 
          && cip->item_list->data.ptrvalue != NULL
          && DoesSequenceMatchSequenceConstraint ((BioseqPtr) cip->item_list->data.ptrvalue, scp)) {
        cip->chosen = do_choose;
      } else {
          ChooseCategoriesByMacroSequenceConstraint (cip->subcategories, scp, do_choose);
      }
    }
    value_list = value_list->next;
  }
}


static void SelectSequenceIDsForRemoval (ButtoN b)
{
  RemoveSeqFromAlignPtr rp;
  SequenceConstraintPtr  scp;

  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp == NULL) return;

  scp = DialogToPointer (rp->constraint_dlg);
  ChooseCategoriesByMacroSequenceConstraint (rp->seq_list, scp, TRUE);
  PointerToDialog (rp->clickable_list_dlg, rp->seq_list);
  scp = SequenceConstraintFree (scp);
}


static void UnselectAllSequences (ButtoN b)
{
  RemoveSeqFromAlignPtr rp;

  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp != NULL) {
    ChooseCategories (rp->seq_list, FALSE);
    PointerToDialog (rp->clickable_list_dlg, rp->seq_list);
  }
}


static void ResortMarkedSegCategories (ButtoN b)
{
  RemoveSeqFromAlignPtr rp;

  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp != NULL) {
    rp->seq_list = ValNodeSort (rp->seq_list, SortVnpByClickableItemChosen);
    PointerToDialog (rp->clickable_list_dlg, rp->seq_list);
  }
}


static void RemoveSequencesFromRecordBaseForm (BaseFormPtr bfp)
{
  WindoW                   w;
  RemoveSeqFromAlignPtr    rp;
  GrouP                    h, c;
  ButtoN                   mark_btn, b;
  SeqEntryPtr              sep;
  SequenceConstraint       seqd;
  StringConstraint         scd;
  PrompT                   ppt, ppt2;

  if (bfp == NULL) return;

  rp = (RemoveSeqFromAlignPtr) MemNew (sizeof (RemoveSeqFromAlignData));
  if (rp == NULL) return;
  rp->input_entityID = bfp->input_entityID;
  rp->sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (rp->sep == NULL) 
  {
	  MemFree (rp);
	  return;
  }
  
  rp->remove_all_from_alignments = FALSE;
  rp->remove_all_products = FALSE;
  rp->no_remove_all_from_alignments = FALSE;
  rp->no_remove_all_products = FALSE;
  
  w = FixedWindow (-50, -33, -10, -10, "Remove Sequences From Record", StdCloseWindowProc);
  if (w == NULL) 
  {
	  MemFree (rp);
	  return;
  }
  rp->form = (ForM) w;
  SetObjectExtra (w, rp, StdCleanupFormProc);
  
  h = HiddenGroup (w, -1, 0, NULL);

  rp->clickable_list_dlg = CreateClickableListDialogExEx (h, "Sequences to Remove", "",
                                                      "Use checkbox to mark sequences to remove",
                                                      "Single click to navigate to sequence in record",
                                                      ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                      GetDiscrepancyItemText,
                                                      stdCharWidth * 30,
                                                      stdCharWidth * 30 + 5,
                                                      TRUE, FALSE, TRUE);

  ppt2 = StaticPrompt (h, "Mark sequences where:", 0, dialogTextHeight, programFont, 'c');
  rp->constraint_dlg = SequenceConstraintDialog (h, NULL, NULL);

  /* set up default constraint */
  MemSet (&seqd, 0, sizeof (SequenceConstraint));
  MemSet (&scd, 0, sizeof (StringConstraint));
  scd.match_location = String_location_inlist;
  seqd.id = &scd;
  PointerToDialog (rp->constraint_dlg, &seqd);

  mark_btn = PushButton (h, "Mark", SelectSequenceIDsForRemoval);
  SetObjectExtra (mark_btn, rp, NULL);

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, &(rp->seq_list), ListAllSequences);

  PointerToDialog (rp->clickable_list_dlg, rp->seq_list);

  ppt = StaticPrompt (h, "Warning - this will remove sequence from Sequin and SMART", 0, dialogTextHeight, programFont, 'c');

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveSequencesFromRecord);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Unmark All Sequences", UnselectAllSequences);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Resort Marked", ResortMarkedSegCategories);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, rp, NULL);
  AlignObjects (ALIGN_CENTER, 
                (HANDLE) rp->clickable_list_dlg,
                (HANDLE) ppt2, 
                (HANDLE) rp->constraint_dlg, 
                (HANDLE) mark_btn, 
                (HANDLE) ppt,
                (HANDLE) c,
                NULL);

  RealizeWindow (w);
  Show (w);
  Update ();  
}


static void MarkBFPRemoveListBtn (ButtoN b)
{
  RemoveSeqFromAlignPtr    rp;
  SequenceConstraintPtr      scp;

  rp = (RemoveSeqFromAlignPtr)GetObjectExtra (b);
  if (rp == NULL) {
    return;
  }

  scp = SequenceConstraintNew ();
  scp->id = StringConstraintNew ();
  scp->id->match_location = String_location_inlist;
  scp->id->case_sensitive = FALSE;
  scp->id->match_text = SaveStringFromText (rp->mark_seqid_list);
  ChooseCategoriesByMacroSequenceConstraint (rp->seq_list, scp, TRUE);
  PointerToDialog (rp->clickable_list_dlg, rp->seq_list);
  scp = SequenceConstraintFree (scp);
}


static void MarkBFPRemoveLenBtn (ButtoN b)
{
  RemoveSeqFromAlignPtr    rp;
  SequenceConstraintPtr      scp;
  CharPtr txt;
  Int4    num;

  rp = (RemoveSeqFromAlignPtr)GetObjectExtra (b);
  if (rp == NULL) {
    return;
  }

  txt = SaveStringFromText (rp->mark_seqlen);
  num = atoi (txt);
  txt = MemFree (txt);
  if (num < 1) {
    return;
  }
  scp = SequenceConstraintNew ();
  scp->length = ValNodeNew (NULL);
  scp->length->choice = QuantityConstraint_less_than;
  scp->length->data.intvalue = num;
  ChooseCategoriesByMacroSequenceConstraint (rp->seq_list, scp, TRUE);
  PointerToDialog (rp->clickable_list_dlg, rp->seq_list);
  scp = SequenceConstraintFree (scp);
}


NLM_EXTERN void SubmitterRemoveSequencesFromRecordBaseForm (BaseFormPtr bfp)
{
  WindoW                   w;
  RemoveSeqFromAlignPtr    rp;
  GrouP                    h, mark_grp, sg, c;
  ButtoN                   b, unselect;
  SeqEntryPtr              sep;

  if (bfp == NULL) return;

  rp = (RemoveSeqFromAlignPtr) MemNew (sizeof (RemoveSeqFromAlignData));
  if (rp == NULL) return;
  rp->input_entityID = bfp->input_entityID;
  rp->sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (rp->sep == NULL) 
  {
	  MemFree (rp);
	  return;
  }
  
  rp->remove_all_from_alignments = FALSE;
  rp->remove_all_products = FALSE;
  rp->no_remove_all_from_alignments = FALSE;
  rp->no_remove_all_products = FALSE;
  
  w = FixedWindow (-50, -33, -10, -10, "Delete Sequences", StdCloseWindowProc);
  if (w == NULL) 
  {
	  MemFree (rp);
	  return;
  }
  rp->form = (ForM) w;
  SetObjectExtra (w, rp, StdCleanupFormProc);
  
  h = HiddenGroup (w, -1, 0, NULL);

  rp->clickable_list_dlg = CreateClickableListDialogExEx (h, "Sequences in your file", "Sequences selected for deletion",
                                                      "Use checkbox to select sequences for deletion",
                                                      "Single click to navigate to sequence in record",
                                                      ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                      GetDiscrepancyItemText,
                                                      stdCharWidth * 30,
                                                      stdCharWidth * 30 + 5,
                                                      TRUE, FALSE, TRUE);

  mark_grp = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (mark_grp, 10, 10);
  b = PushButton (mark_grp, "Select sequences in this list:", MarkBFPRemoveListBtn);
  SetObjectExtra (b, rp, NULL);
  rp->mark_seqid_list = DialogText (mark_grp, "", 30, NULL);
  b = PushButton (mark_grp, "Select sequences less than:", MarkBFPRemoveLenBtn);
  SetObjectExtra (b, rp, NULL);
  sg = HiddenGroup (mark_grp, 2, 0, NULL);
  SetGroupSpacing (sg, 10, 10);
  rp->mark_seqlen = DialogText (sg, "200", 10, NULL);
  StaticPrompt (sg, "nucleotides in length", 0, 0, programFont, 'l');

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, &(rp->seq_list), ListAllSequences);

  PointerToDialog (rp->clickable_list_dlg, rp->seq_list);

  unselect = PushButton (h, "Unselect All Sequences", UnselectAllSequences);
  SetObjectExtra (unselect, rp, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveSequencesFromRecord);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, rp, NULL);
  AlignObjects (ALIGN_CENTER, 
                (HANDLE) rp->clickable_list_dlg,
                (HANDLE) mark_grp, 
                (HANDLE) unselect,
                (HANDLE) c,
                NULL);

  RealizeWindow (w);
  Show (w);
  Update ();  
}



NLM_EXTERN void ButtonOrMenuItemTemplate (IteM i, MenuOrButtonFunc func)
{
  BaseFormPtr              bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL || func == NULL) return;
  func(bfp);
}


NLM_EXTERN void ButtonOrMenuButtonTemplate (ButtoN b, MenuOrButtonFunc func)
{
  BaseFormPtr              bfp;

  bfp = GetObjectExtra (b);

  if (bfp == NULL || func == NULL) return;
  func(bfp);
}


extern void RemoveSequencesFromRecord (IteM i)
{
  ButtonOrMenuItemTemplate(i, RemoveSequencesFromRecordBaseForm);
}


extern void RemoveSequencesFromRecordBtn (ButtoN b)
{
  ButtonOrMenuButtonTemplate(b, RemoveSequencesFromRecordBaseForm);
}


typedef struct removeseqfromwizardform {
  FORM_MESSAGE_BLOCK
  DialoG clickable_list_dlg;
  TexT   mark_seqid_list;
  TexT   mark_seqlen;
  ValNodePtr seq_list;
} RemoveSeqFromWizardFormData, PNTR RemoveSeqFromWizardFormPtr;


static void UnselectAllSequencesForRemove (ButtoN b)
{
  RemoveSeqFromWizardFormPtr frm;

  frm = (RemoveSeqFromWizardFormPtr) GetObjectExtra (b);
  if (frm != NULL) {
    ChooseCategories (frm->seq_list, FALSE);
    PointerToDialog (frm->clickable_list_dlg, frm->seq_list);
  }
}


static void MarkRemoveListBtn (ButtoN b)
{
  RemoveSeqFromWizardFormPtr frm;
  SequenceConstraintPtr      scp;

  frm = (RemoveSeqFromWizardFormPtr) GetObjectExtra (b);
  if (frm == NULL || TextHasNoText (frm->mark_seqid_list)) {
    return;
  }

  scp = SequenceConstraintNew ();
  scp->id = StringConstraintNew ();
  scp->id->match_location = String_location_inlist;
  scp->id->case_sensitive = FALSE;
  scp->id->match_text = SaveStringFromText (frm->mark_seqid_list);
  ChooseCategoriesByMacroSequenceConstraint (frm->seq_list, scp, TRUE);
  PointerToDialog (frm->clickable_list_dlg, frm->seq_list);
  scp = SequenceConstraintFree (scp);
}


static void MarkRemoveLenBtn (ButtoN b)
{
  RemoveSeqFromWizardFormPtr frm;
  SequenceConstraintPtr      scp;
  CharPtr txt;
  Int4    num;

  frm = (RemoveSeqFromWizardFormPtr) GetObjectExtra (b);
  if (frm == NULL || TextHasNoText (frm->mark_seqlen)) {
    return;
  }
  txt = SaveStringFromText (frm->mark_seqlen);
  num = atoi (txt);
  txt = MemFree (txt);
  if (num < 1) {
    return;
  }
  scp = SequenceConstraintNew ();
  scp->length = ValNodeNew (NULL);
  scp->length->choice = QuantityConstraint_less_than;
  scp->length->data.intvalue = num;
  ChooseCategoriesByMacroSequenceConstraint (frm->seq_list, scp, TRUE);
  PointerToDialog (frm->clickable_list_dlg, frm->seq_list);
  scp = SequenceConstraintFree (scp);
}


NLM_EXTERN Boolean RemoveSequencesFromWizardList (ValNodePtr PNTR sequences)
{
  WindoW w;
  RemoveSeqFromWizardFormPtr frm;
  GrouP h, mark_grp, sg, c;
  ButtoN b, unselect;
  ModalAcceptCancelData acd;
  ValNodeBlock removal_list;
  ValNodePtr   vnp;
  SeqEntryPtr  sep, sep_prev = NULL, sep_next;
  ClickableItemPtr cip;
  Boolean          rval = FALSE;

  if (sequences == NULL || *sequences == NULL) {
    return FALSE;
  }

  InitValNodeBlock (&removal_list, NULL);
  for (sep = *sequences; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      cip = ClickableItemForBioseq (sep->data.ptrvalue);
      ValNodeAddPointerToEnd (&removal_list, 0, cip);
    }
  }

  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  acd.third_option = FALSE;

  frm = (RemoveSeqFromWizardFormPtr) MemNew (sizeof (RemoveSeqFromWizardFormData));
  w = ModalWindow(-20, -13, -10, -10, NULL);
  frm->form = (ForM) w;
  SetObjectExtra (w, frm, StdCleanupFormProc);
  frm->seq_list = removal_list.head;
  
  h = HiddenGroup (w, -1, 0, NULL);


  frm->clickable_list_dlg = CreateClickableListDialogExEx (h, "Sequences in your file", "Sequences selected for deletion",
                                                      "Use checkbox to select sequences to delete",
                                                      "",
                                                      NULL, NULL, NULL,
                                                      GetDiscrepancyItemText,
                                                      stdCharWidth * 30,
                                                      stdCharWidth * 30 + 5,
                                                      TRUE, FALSE, TRUE);
  PointerToDialog (frm->clickable_list_dlg, frm->seq_list);
  
  mark_grp = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (mark_grp, 10, 10);
  b = PushButton (mark_grp, "Select sequences in this list:", MarkRemoveListBtn);
  SetObjectExtra (b, frm, NULL);
  frm->mark_seqid_list = DialogText (mark_grp, "", 30, NULL);
  b = PushButton (mark_grp, "Select sequences less than:", MarkRemoveLenBtn);
  SetObjectExtra (b, frm, NULL);
  sg = HiddenGroup (mark_grp, 2, 0, NULL);
  SetGroupSpacing (sg, 10, 10);
  frm->mark_seqlen = DialogText (sg, "200", 10, NULL);
  StaticPrompt (sg, "nucleotides in length", 0, 0, programFont, 'l');
  
  unselect = PushButton (h, "Unselect All Sequences", UnselectAllSequencesForRemove);
  SetObjectExtra (unselect, frm, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Delete and Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  AlignObjects (ALIGN_CENTER, 
                (HANDLE) frm->clickable_list_dlg,
                (HANDLE) mark_grp,
                (HANDLE) unselect,
                (HANDLE) c,
                NULL);

  RealizeWindow (w);
  Show (w);
  Update ();  
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (acd.accepted) 
  {
    /* actually remove the sequences */
    for (sep = *sequences, vnp = frm->seq_list;
         sep != NULL && vnp != NULL;
         sep = sep_next, vnp = vnp->next) 
    {
      sep_next = sep->next;
      if (IS_Bioseq (sep) 
          && (cip = (ClickableItemPtr) (vnp->data.ptrvalue)) != NULL
          && cip->chosen)
      {
        if (sep_prev == NULL) 
        {
          *sequences = sep_next;
        }
        else
        {
          sep_prev->next = sep_next;
        }
        sep->next = NULL;
        sep = SeqEntryFree (sep);
        rval = TRUE;
      }
      else 
      {
        sep_prev = sep;
      }
    }
  }
  frm->seq_list = FreeClickableList (frm->seq_list);
  Remove (w);
  return rval;
}


