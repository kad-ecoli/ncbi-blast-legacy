/*   sqnutil2.c
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
* File Name:  sqnutil2.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   9/2/97
*
* $Revision: 6.515 $
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
#include <simple.h>
#include <validerr.h>
#include <findrepl.h>
#include <alignmgr2.h>
#include <alignval.h>
#include <objvalid.h>
#include <valapi.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

#include <utilpub.h>

static CharPtr SqnTrimSpacesAroundString (CharPtr str)

{
  Uchar    ch;  /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
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
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
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

static CharPtr SqnStringSave (CharPtr from)

{
  size_t  len;
  CharPtr to;

  to = NULL;
  len = StringLen (from);
  if (len > 0) {
    to = (CharPtr) MemGet (len + 1, FALSE);
    if (to != NULL) {
      MemCpy (to, from, len + 1);
      SqnTrimSpacesAroundString (to);
    }
  }
  return to;
}

NLM_EXTERN void UpdateLocalId (BioseqPtr bsp, CharPtr localId)

{
  Char         ch;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  SeqIdPtr     sip;
  long         val;

  if (bsp != NULL) {
    if (localId != NULL) {
      sip = bsp->id;
      while (sip != NULL && sip->choice != SEQID_LOCAL) {
        sip = sip->next;
      }
      oip = NULL;
      if (sip != NULL) {
        oip = (ObjectIdPtr) sip->data.ptrvalue;
      } else {
        sip = ValNodeNew (bsp->id);
        if (bsp->id == NULL) {
          bsp->id = sip;
        }
        if (sip != NULL) {
          oip = ObjectIdNew ();
          sip->choice = SEQID_LOCAL;
          sip->data.ptrvalue = (Pointer) oip;
        }
      }
      if (oip != NULL) {
        oip->str = MemFree (oip->str);
        if (sscanf (localId, "%ld", &val) == 1) {
          oip->id = (Int4) val;
        } else {
          oip->str = SqnStringSave (localId);
          ptr = oip->str;
          ch = *ptr;
          while (ch != '\0') {
            if (ch == '|') {
              *ptr = '~';
            }
            ptr++;
            ch = *ptr;
          }
        }
      }
      SeqMgrReplaceInBioseqIndex (bsp);
    }
  }
}

NLM_EXTERN void UpdateTitle (BioseqPtr bsp, CharPtr title)

{
  ValNodePtr  vnp;

  if (bsp != NULL) {
    if (title != NULL) {
      vnp = NULL;
      if (bsp->descr != NULL) {
        vnp = ValNodeFindNext (bsp->descr, NULL, Seq_descr_title);
      }
      if (vnp == NULL) {
        vnp = ValNodeNew (bsp->descr);
        if (vnp != NULL) {
          vnp->choice = Seq_descr_title;
        }
        if (bsp->descr == NULL) {
          bsp->descr = vnp;
        }
      }
      if (vnp != NULL) {
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        vnp->data.ptrvalue = SqnStringSave (title);
      }
    }
  }
}

NLM_EXTERN GeneRefPtr CreateNewGeneRef (CharPtr locus, CharPtr allele,
                             CharPtr desc, Boolean pseudo)

{
  GeneRefPtr  geneRef;

  geneRef = GeneRefNew ();
  if (geneRef != NULL) {
    geneRef->locus = SqnStringSave (locus);
    geneRef->allele = SqnStringSave (allele);
    geneRef->desc = SqnStringSave (desc);
    geneRef->pseudo = pseudo;
    if (geneRef->locus == NULL && geneRef->allele == NULL && geneRef->desc == NULL) {
      geneRef = GeneRefFree (geneRef);
    }
  }
  return geneRef;
}

NLM_EXTERN ProtRefPtr CreateNewProtRef (CharPtr name, CharPtr desc,
                             CharPtr ec, CharPtr activity)

{
  ProtRefPtr  protRef;
  ValNodePtr  vnp;

  protRef = ProtRefNew ();
  if (protRef != NULL) {
    if (name != NULL && *name != '\0') {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->data.ptrvalue = SqnStringSave (name);
        protRef->name = vnp;
      }
    }
    protRef->desc = SqnStringSave (desc);
    if (ec != NULL && *ec != '\0') {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->data.ptrvalue = SqnStringSave (ec);
        protRef->ec = vnp;
      }
    }
    if (activity != NULL && *activity != '\0') {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->data.ptrvalue = SqnStringSave (activity);
        protRef->activity = vnp;
      }
    }
    if (protRef->name == NULL && protRef->desc == NULL &&
        protRef->ec == NULL && protRef->activity == NULL) {
      protRef = ProtRefFree (protRef);
    }
  }
  return protRef;
}

NLM_EXTERN CdRegionPtr CreateNewCdRgn (Uint1 frame, Boolean orf, Int2 genCode)

{
  CdRegionPtr  cdRgn;
  ValNodePtr   code;
  ValNodePtr   vnp;

  cdRgn = CdRegionNew ();
  if (cdRgn != NULL) {
    cdRgn->orf = orf;
    cdRgn->conflict = FALSE;
    cdRgn->frame = frame;
    cdRgn->gaps = 0;
    cdRgn->mismatch = 0;
    cdRgn->stops = 0;
    code = ValNodeNew (NULL);
    if (code != NULL) {
      code->choice = 254;
      vnp = ValNodeNew (NULL);
      code->data.ptrvalue = vnp;
      if (vnp != NULL) {
        vnp->choice = 2;
        vnp->data.intvalue = (Int4) genCode;
      }
    }
    cdRgn->genetic_code = code;
    cdRgn->code_break = NULL;
  }
  return cdRgn;
}

NLM_EXTERN void SetSeqFeatData (SeqFeatPtr sfp, Pointer data)

{
  if (sfp != NULL) {
    sfp->data.value.ptrvalue = (Pointer) data;
  }
}

NLM_EXTERN void SetSeqFeatProduct (SeqFeatPtr sfp, BioseqPtr bsp)

{
  ValNodePtr  slp;

  if (sfp != NULL) {
    sfp->product = SeqLocFree (sfp->product);
    if (bsp != NULL && bsp->id != NULL) {
      slp = ValNodeNew (NULL);
      if (slp != NULL) {
        slp->choice = 3;
        slp->data.ptrvalue = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
      }
      sfp->product = slp;
    }
  }
}

NLM_EXTERN void ResetSeqFeatInterval (SeqFeatPtr sfp)

{
  if (sfp != NULL) {
    sfp->location = SeqLocFree (sfp->location);
  }
}

NLM_EXTERN void AddSeqFeatInterval (SeqFeatPtr sfp, BioseqPtr bsp, Int4 from,
                         Int4 to, Boolean partial5, Boolean partial3)

{
  Int2  fuzz_from;
  Int2  fuzz_to;
  Int2  strand;
  Int4  tmp;

  if (sfp != NULL && bsp != NULL) {
    strand = Seq_strand_plus;
    if (from > to) {
      tmp = from;
      from = to;
      to = tmp;
      strand = Seq_strand_minus;
    }
    fuzz_from = -1;
    fuzz_to = -1;
    if (partial5) {
      fuzz_from = 2;
    }
    if (partial3) {
      fuzz_to = 1;
    }
    AddIntToSeqFeat (sfp, from - 1, to - 1, bsp, fuzz_from, fuzz_to, strand);
  }
}

NLM_EXTERN void AddSeqLocPoint (SeqLocPtr PNTR old_slp, SeqIdPtr sip, Int4 location,
                      Boolean fuzz_before, Boolean fuzz_after, Int2 strand)

{
    SeqLocPtr slp, tmp, tmp2;
    SeqPntPtr spp;
    IntFuzzPtr ifp;
    Int2 fuzz;

  if (old_slp == NULL)
  {
    return;
  }
    spp = SeqPntNew();
    spp->point = location - 1;
    spp->id = SeqIdDup(sip);
    spp->strand = (Uint1)strand;

    fuzz = -1;
    if (fuzz_before) {
      fuzz = 4;        /* tl */
    } else if (fuzz_after) {
      fuzz = 3;        /* tr */
    }
    if (fuzz >= 0)
    {
        ifp = IntFuzzNew();
        ifp->choice = 4;   /* lim */
        ifp->a = (Int4)fuzz;
        spp->fuzz = ifp;
    }

    slp = ValNodeNew(NULL);
    slp->choice = SEQLOC_PNT;
    slp->data.ptrvalue = (Pointer)spp;

    if (*old_slp == NULL)
    {
        *old_slp = slp;
        return;
    }

    tmp = *old_slp;
    if (tmp->choice == SEQLOC_MIX)   /* second one already */
    {
        tmp2 = (ValNodePtr)(tmp->data.ptrvalue);
        while (tmp2->next != NULL)
            tmp2 = tmp2->next;
        tmp2->next = slp;
    }
    else                             /* create a chain */
    {
        tmp2 = ValNodeNew(NULL);
        tmp2->choice = SEQLOC_MIX;
        tmp2->data.ptrvalue = (Pointer)tmp;
        tmp->next = slp;
        *old_slp = tmp2;
    }
}

NLM_EXTERN void AddSeqFeatPoint (SeqFeatPtr sfp, BioseqPtr bsp, Int4 location,
                      Boolean fuzz_before, Boolean fuzz_after, Int2 strand)

{
  AddSeqLocPoint (&(sfp->location), SeqIdFindBest(bsp->id, 0), location,
                  fuzz_before, fuzz_after, strand);
}

typedef struct seqlocrange {
  Int4        left;
  Int4        right;
  Uint1        strand;
  Uint1     choice;
  struct seqlocrange PNTR next;
 } SeqLocRange, PNTR SeqLocRangePtr;
 
static SeqLocRangePtr SeqLocRangeFree (SeqLocRangePtr slrp)

{
  SeqLocRangePtr  next;

  while (slrp != NULL) {
    next = slrp->next;
    MemFree (slrp);
    slrp = next;
  }
  return NULL;
}


static Boolean IsLocationOnCircularBioseq (SeqLocPtr slp)
{
  BioseqPtr bsp;
  Boolean is_circular = FALSE;

  bsp = BioseqFind (SeqLocId(slp));
  if (bsp != NULL && bsp->topology == TOPOLOGY_CIRCULAR) {
      is_circular = TRUE;
  }
  return is_circular;
}


static SeqLocRangePtr CollectRanges (BioseqPtr target, SeqLocPtr slp)

{
  SeqLocRangePtr  change;
  SeqLocPtr       curr;
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;
  Int4            left;
  Int4            right;
  SeqLocRangePtr  slrp;
  Uint1           strand;
  Boolean         is_circular;
  Boolean         left_flip = FALSE, right_flip = FALSE;

  if (target == NULL) return NULL;
  head = NULL;
  last = NULL;
  curr = SeqLocFindNext (slp, NULL);
  while (curr != NULL) {
    if (curr->choice != SEQLOC_NULL) {
      is_circular = IsLocationOnCircularBioseq (curr);
#if 1
      GetLeftAndRightOffsetsInBioseq (curr, target, &left, &right, is_circular, &left_flip, &right_flip);
#else
      left = GetOffsetInBioseqEx (curr, target, SEQLOC_LEFT_END, is_circular);
      right = GetOffsetInBioseqEx (curr, target, SEQLOC_RIGHT_END, is_circular);
#endif
      strand = SeqLocStrand (curr);
      /* left > right if within a minus strand delta seq component, flip strand here */
      if (left > right || left_flip || right_flip) {
        if (strand == Seq_strand_minus) {
          strand = Seq_strand_plus;
        } else {
          strand = Seq_strand_minus;
        }
      }
      if (left != -1 && right != -1) {
        slrp = MemNew (sizeof (SeqLocRange));
        if (slrp != NULL) {
          slrp->left = left;
          slrp->right = right;
          slrp->strand = strand;
          slrp->choice = curr->choice;
          if (head == NULL) {
            head = slrp;
          } else if (last != NULL) {
            last->next = slrp;
          } else {
            ErrPostEx (SEV_ERROR, 0, 0, "SeqLocMerge list problem");
            SeqLocRangeFree (head);
            return NULL;
          }
          last = slrp;
        }
      }
    }
    curr = SeqLocFindNext (slp, curr);
  }
  if (head == NULL || target->topology != TOPOLOGY_CIRCULAR) return head;
#if 1
  GetLeftAndRightOffsetsInBioseq (slp, target, &left, &right, target->topology == TOPOLOGY_CIRCULAR, &left_flip, &right_flip);
#else
  left = GetOffsetInBioseqEx (slp, target, SEQLOC_LEFT_END, target->topology == TOPOLOGY_CIRCULAR);
  right = GetOffsetInBioseqEx (slp, target, SEQLOC_RIGHT_END, target->topology == TOPOLOGY_CIRCULAR);
#endif
  if (left == -1 || right == -1 || left <= right) return head;
  /* feature spans origin */
  change = NULL;
  left = head->left;
  strand = SeqLocStrand (slp);
  if (strand == Seq_strand_minus) {
    for (slrp = head->next; slrp != NULL; slrp = slrp->next) {
      if (slrp->left > left && change == NULL) {
        change = slrp;
      }
    }
  } else {
    for (slrp = head->next; slrp != NULL; slrp = slrp->next) {
      if (slrp->left < left && change == NULL) {
        change = slrp;
      }
    }
  }
  if (change == NULL) return head;
  if (strand == Seq_strand_minus) {
    for (slrp = change; slrp != NULL; slrp = slrp->next) {
      slrp->left -= target->length;
      slrp->right -= target->length;
    }
  } else {
    for (slrp = head; slrp != NULL && slrp != change; slrp = slrp->next) {
      slrp->left -= target->length;
      slrp->right -= target->length;
    }
  }
  return head;
}

static int LIBCALLBACK CompareRanges (VoidPtr ptr1, VoidPtr ptr2)

{
  SeqLocRangePtr   slrp1;
  SeqLocRangePtr   slrp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    slrp1 = *((SeqLocRangePtr PNTR) ptr1);
    slrp2 = *((SeqLocRangePtr PNTR) ptr2);
    if (slrp1 != NULL && slrp2 != NULL) {
      if (slrp1->left > slrp2->left) {
        return 1;
      } else if (slrp1->left < slrp2->left) {
        return -1;
      } else if (slrp1->right > slrp2->right) {
        return 1;
      } else if (slrp1->right < slrp2->right) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static int LIBCALLBACK CompareReverseRanges (VoidPtr ptr1, VoidPtr ptr2)

{
  return (0 - CompareRanges (ptr1, ptr2));
}

static SeqLocRangePtr SortRanges (SeqLocRangePtr list, Boolean reverse)

{
  size_t          count;
  SeqLocRangePtr  PNTR head;
  size_t          i;
  SeqLocRangePtr  tmp;

  if (list != NULL) {
    count = 0;
    tmp = list;
    while (tmp != NULL) {
      count++;
      tmp = tmp->next;
    }
    if (count > 0) {
      head = MemNew ((count + 1) * sizeof (SeqLocRangePtr));
      if (head != NULL) {
        tmp = list;
        i = 0;
        while (tmp != NULL && i < count) {
          head [i] = tmp;
          tmp = tmp->next;
          i++;
        }
        if (reverse) {
          StableMergeSort (head, count, sizeof (SeqLocRangePtr), CompareReverseRanges);
        } else {
          StableMergeSort (head, count, sizeof (SeqLocRangePtr), CompareRanges);
        }
        for (i = 0; i < count; i++) {
          tmp = head [i];
          tmp->next = head [i + 1];
        }
        list = head [0];
        MemFree (head);
      }
    }
  }
  return list;
}

static SeqLocRangePtr MergeOverlaps (SeqLocRangePtr list, Boolean fuse_joints, Boolean merge_overlaps)

{
  SeqLocRangePtr  last;
  SeqLocRangePtr  next;
  SeqLocRangePtr  this;

  if (list != NULL) {
    this = list->next;
    last = list;
    while (this != NULL) {
      next = this->next;
      if (merge_overlaps && this->left <= last->right) {
        last->right = MAX (this->right, last->right);
        MemFree (this);
        last->next = next;
      } else if (fuse_joints &&
                 (this->left == last->right + 1 && last->right != -1)) {
        last->right = MAX (this->right, last->right);
        MemFree (this);
        last->next = next;
      } else {
        last = this;
      }
      this = next;
    }
  }
  return list;
}

static SeqLocPtr SeqLocFromRange (SeqLocRangePtr head, BioseqPtr target,
                                  Boolean partial5, Boolean partial3,
                                  Boolean add_null)

{
  SeqLocPtr   firstSlp;
  Int4        from;
  Int2        fuzz_from;
  Int2        fuzz_to;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Boolean     notFirst;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  Int2        strand;
  Int4        tmp;
  SeqLocPtr   tmploc1;
  SeqLocPtr   tmploc2;
  Int4        to;
  SeqLocPtr   master_loc = NULL;

  if (head == NULL) return NULL;
  slp = NULL;
  notFirst = FALSE;
  while (head != NULL) {
    fuzz_from = -1;
    fuzz_to = -1;
    from = head->left;
    to = head->right;
    strand = head->strand;
    if (from > to) {
      tmp = from;
      from = to;
      to = tmp;
    }
    if (add_null && notFirst) {
      slp = ValNodeNew (NULL);
      if (slp != NULL) {
        slp->choice = SEQLOC_NULL;
        tmploc1 = master_loc;
        if (tmploc1 != NULL) {
          if (tmploc1->choice == SEQLOC_MIX) {
            tmploc2 = (ValNodePtr) (tmploc1->data.ptrvalue);
            if (tmploc2 != NULL) {
              while (tmploc2->next != NULL) {
                tmploc2 = tmploc2->next;
              }
              tmploc2->next = slp;
            }
          } else {
            tmploc2 = ValNodeNew (NULL);
            if (tmploc2 != NULL) {
              tmploc2->choice = SEQLOC_MIX;
              tmploc2->data.ptrvalue = (Pointer) tmploc1;
              tmploc1->next = slp;
              master_loc = tmploc2;
            }
          }
        }
      }
    }
    if (head->choice == SEQLOC_PNT) {
      AddPntToSeqLoc (&master_loc, from, target, fuzz_from, strand);
    } else {
      AddIntToSeqLoc (&master_loc, from, to, SeqIdFindBest(target->id, 0),
                      fuzz_from, fuzz_to, strand);
    }
    notFirst = TRUE;
    head = head->next;
  }
  firstSlp = NULL;
  lastSlp = NULL;
  slp = SeqLocFindNext (master_loc, NULL);
  while (slp != NULL) {
    if (firstSlp == NULL) {
      firstSlp = slp;
    }
    lastSlp = slp;
    slp = SeqLocFindNext (master_loc, slp);
  }
  if (firstSlp != NULL && firstSlp->choice == SEQLOC_INT &&
      firstSlp->data.ptrvalue != NULL && partial5) {
    sip = (SeqIntPtr) firstSlp->data.ptrvalue;
    ifp = IntFuzzNew ();
    if (ifp != NULL) {
      ifp->choice = 4;
      if (sip->strand == Seq_strand_minus ||
          sip->strand == Seq_strand_both_rev) {
        sip->if_to = ifp;
        ifp->a = 1;
      } else {
        sip->if_from = ifp;
        ifp->a = 2;
      }
    }
  }
  if (lastSlp != NULL && lastSlp->choice == SEQLOC_INT &&
      lastSlp->data.ptrvalue != NULL && partial3) {
    sip = (SeqIntPtr) lastSlp->data.ptrvalue;
    ifp = IntFuzzNew ();
    if (ifp != NULL) {
      ifp->choice = 4;
      if (sip->strand == Seq_strand_minus ||
          sip->strand == Seq_strand_both_rev) {
        sip->if_from = ifp;
        ifp->a = 2;
      } else {
        sip->if_to = ifp;
        ifp->a = 1;
      }
    }
  }
  return master_loc;
}


static SeqLocPtr SimpleMerge (BioseqPtr target, SeqLocPtr to, SeqLocPtr from)
{
  SeqIntPtr sint;

  if (target != NULL && to != NULL && from == NULL
      && to->choice == SEQLOC_INT && (sint = (SeqIntPtr) to->data.ptrvalue) != NULL
      && SeqIdIn (sint->id, target->id)) {
    return SeqLocCopy (to);
  } else {
    return NULL;
  }
}

   
NLM_EXTERN SeqLocPtr SeqLocMergeExEx (BioseqPtr target, SeqLocPtr to, SeqLocPtr from,
                       Boolean single_interval, Boolean fuse_joints,
                       Boolean merge_overlaps, Boolean add_null,
                       Boolean ignore_mixed, Boolean ignore_out_of_order)

{
  SeqLocRangePtr  curr;
  SeqLocRangePtr  slrp;
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;
  SeqLocRangePtr  tmp;
  Boolean         mixed;
  Boolean         partial5;
  Boolean         partial3;
  SeqLocPtr       slp;
  Uint1           strand;
  Boolean         unordered;

  if (target == NULL) return NULL;
  if (to == NULL && from == NULL) return NULL;

  if ((slp = SimpleMerge (target, to, from)) != NULL) {
    return slp;
  }

  slp = NULL;
  partial5 = FALSE;
  partial3 = FALSE;
  head = CollectRanges (target, to);
  if (head == NULL) {
    head = CollectRanges (target, from);
  } else {
    last = head;
    while (last->next != NULL) {
      last = last->next;
    }
    last->next = CollectRanges (target, from);
  }
  if (head != NULL) {

    /* test for mixed strands */
    mixed = FALSE;
    unordered = FALSE;
    last = head;
    strand = head->strand;
    curr = head->next;
    while (curr != NULL) {
      if (curr->strand == Seq_strand_minus) {
        if (strand == Seq_strand_plus || strand == Seq_strand_unknown) {
          mixed = TRUE;
        }
        if (last->right < curr->right) {
          unordered = TRUE;
        }
      } else {
        if (strand == Seq_strand_minus) {
          mixed = TRUE;
        }
        if (last->left > curr->left) {
          unordered = TRUE;
        }
      }
      last = curr;
      curr = curr->next;
    }

    /* but can override mixed strands behavior */
    if (ignore_mixed) {
      mixed = FALSE;
    }
    if (ignore_out_of_order) {
      unordered = FALSE;
    }

    if ((! mixed) && (! unordered)) {
      strand = head->strand;
      head = SortRanges (head, FALSE);
      head = MergeOverlaps (head, fuse_joints, merge_overlaps);
      if (single_interval) {
        last = head;
        while (last->next != NULL) {
          last = last->next;
        }
        if (head->left < 0 && last->left >= 0)
        {
          head->right = -1;
          last->left = 0;
          if (head->next != last)
          {
            /* remove intervening intervals */
            tmp = head->next;
            head->next = last;
            for (last = tmp; last->next != head->next; last = last->next)
            {
            }
            last->next = NULL;
            SeqLocRangeFree (tmp);
          }
        }
        else
        {
          head->left = MIN (head->left, last->left);
          head->right = MAX (head->right, last->right);
          head->next = SeqLocRangeFree (head->next);
        }
      }
      last = head;
      while (last != NULL) {
        last->strand = strand;
        last = last->next;
      }
      if (strand == Seq_strand_minus) {
        head = SortRanges (head, TRUE);
      }
    }

    for (slrp = head; slrp != NULL; slrp = slrp->next) {
      if (slrp->left < 0) {
        slrp->left += target->length;
      }
      if (slrp->right < 0) {
        slrp->right += target->length;
      }
    }
    slp = SeqLocFromRange (head, target, partial5, partial3, add_null);
    head = SeqLocRangeFree (head);
  }
  return slp;
}

NLM_EXTERN SeqLocPtr SeqLocMergeEx (BioseqPtr target, SeqLocPtr to, SeqLocPtr from,
                       Boolean single_interval, Boolean fuse_joints,
                       Boolean merge_overlaps, Boolean add_null)

{
  return SeqLocMergeExEx (target, to, from, single_interval, fuse_joints, merge_overlaps, add_null, FALSE, TRUE);
}

NLM_EXTERN SeqLocPtr SeqLocMerge (BioseqPtr target, SeqLocPtr to, SeqLocPtr from,
                       Boolean single_interval, Boolean fuse_joints,
                       Boolean add_null)

{
  return SeqLocMergeExEx (target, to, from, single_interval, fuse_joints, TRUE, add_null, FALSE, TRUE);
}

NLM_EXTERN Boolean SeqLocBadSortOrder (BioseqPtr bsp, SeqLocPtr slp)

{
  SeqLocRangePtr  curr;
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;

  if (bsp == NULL || slp == NULL) return FALSE;
  if (SeqLocCheck (slp) == SEQLOCCHECK_WARNING) return FALSE;
  /*
  if (SeqLocId (slp) == NULL) return FALSE;
  */
  head = CollectRanges (bsp, slp);
  if (head == NULL) return FALSE;
  if (head->next == NULL) {
    SeqLocRangeFree (head);
    return FALSE;
  }
  last = head;
  curr = head->next;
  while (curr != NULL) {
    if (curr->strand == Seq_strand_minus) {
      if (last->right < curr->right) {
        SeqLocRangeFree (head);
        return TRUE;
      }
    } else {
      if (last->left > curr->left) {
        SeqLocRangeFree (head);
        return TRUE;
      }
    }
    last = curr;
    curr = curr->next;
  }
  SeqLocRangeFree (head);
  return FALSE;
}

NLM_EXTERN Boolean SeqLocMixedStrands (BioseqPtr bsp, SeqLocPtr slp)

{
  SeqLocRangePtr  curr;
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;
  Uint1           strand;

  if (bsp == NULL || slp == NULL) return FALSE;
  if (SeqLocCheck (slp) == SEQLOCCHECK_WARNING) return FALSE;
  /*
  if (SeqLocId (slp) == NULL) return FALSE;
  */
  head = CollectRanges (bsp, slp);
  if (head == NULL) return FALSE;
  if (head->next == NULL) {
    SeqLocRangeFree (head);
    return FALSE;
  }
  last = head;
  strand = last->strand;
  curr = head->next;
  while (curr != NULL) {
    if (curr->strand == Seq_strand_minus) {
      if (strand == Seq_strand_plus || strand == Seq_strand_unknown) {
        SeqLocRangeFree (head);
        return TRUE;
      }
    } else {
      if (strand == Seq_strand_minus) {
        SeqLocRangeFree (head);
        return TRUE;
      }
    }
    last = curr;
    curr = curr->next;
  }
  SeqLocRangeFree (head);
  return FALSE;
}

static void ConvertToFeatsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent, Boolean toProts)

{
  BioseqPtr   bsp;
  SeqFeatPtr  sfp;
  ValNodePtr  vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol) && (! toProts)) return;
  vnp = (ValNodePtr) mydata;
  if (vnp == NULL) return;
  while (vnp != NULL) {
    switch (vnp->choice) {
      case Seq_descr_pub :
        sfp = CreateNewFeature (sep, NULL, SEQFEAT_PUB, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = AsnIoMemCopy ((Pointer) vnp->data.ptrvalue,
                                                   (AsnReadFunc) PubdescAsnRead,
                                                   (AsnWriteFunc) PubdescAsnWrite);
        }
        break;
      case Seq_descr_source :
        sfp = CreateNewFeature (sep, NULL, SEQFEAT_BIOSRC, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = AsnIoMemCopy ((Pointer) vnp->data.ptrvalue,
                                                   (AsnReadFunc) BioSourceAsnRead,
                                                   (AsnWriteFunc) BioSourceAsnWrite);
        }
        break;
      case Seq_descr_comment :
        sfp = CreateNewFeature (sep, NULL, SEQFEAT_COMMENT, NULL);
        if (sfp != NULL) {
          sfp->comment = StringSave ((CharPtr) vnp->data.ptrvalue);
        }
        break;
      default :
        break;
    }
    vnp = vnp->next;
  }
}

static void ConvertToFeatsOnNucsAndProts (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  ConvertToFeatsCallback (sep, mydata, index, indent, TRUE);
}

static void ConvertToFeatsOnNucsOnly (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  ConvertToFeatsCallback (sep, mydata, index, indent, FALSE);
}

typedef struct descstringcheck {
    CharPtr  findstring;
    Boolean  stringfound;
    AsnIoPtr aip;
} DescStringCheckData, PNTR DescStringCheckPtr;

static void LIBCALLBACK AsnWriteConvertForDCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr            pchFind;
  CharPtr            pchSource;
  DescStringCheckPtr dscp;

  dscp = (DescStringCheckPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) {
    pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
    pchFind = dscp->findstring;
    if (StringSearch (pchSource, pchFind) != NULL) {
      dscp->stringfound = TRUE;
    }
  }
}

static Boolean PubSrcComHasSubstring (ObjMgrTypePtr omtp, Pointer ptr, DescStringCheckPtr dscp)
{
  if (omtp == NULL || dscp == NULL || dscp->findstring == NULL || StringHasNoText (dscp->findstring)) {
    return TRUE;
  }
  if (ptr == NULL || dscp->aip == NULL) {
    return FALSE;
  }
  dscp->stringfound = FALSE;
  (omtp->asnwrite) (ptr, dscp->aip, NULL);
  return dscp->stringfound;
}

static void 
ExtractPubSrcComDescs 
(SeqEntryPtr sep,
 ValNodePtr PNTR head,
 Boolean pub,
 Boolean src,
 Boolean com,
 ObjMgrTypePtr omtp,
 DescStringCheckPtr dscp
 )

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;
  Boolean       ok_to_extract;

  if (sep == NULL || head == NULL) return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    ok_to_extract = FALSE;
    if ((sdp->choice == Seq_descr_pub && pub) ||
        (sdp->choice == Seq_descr_source && src) ||
        (sdp->choice == Seq_descr_comment && com)) {
      if (PubSrcComHasSubstring (omtp, sdp, dscp)) {
        ok_to_extract = TRUE;
      }
    }
    if (ok_to_extract) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      ValNodeLink (head, sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}




static Boolean 
HasPubSrcComDescriptors 
(SeqEntryPtr sep,
 Boolean pub,
 Boolean src,
 Boolean com,
 ObjMgrTypePtr omtp,
 DescStringCheckPtr dscp
)
{
  BioseqPtr           bsp;
  BioseqSetPtr        bssp;
  ValNodePtr          list = NULL;
  Boolean             rval = FALSE;

  if (sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (IS_Bioseq (sep)) {
    bsp = sep->data.ptrvalue;
    list = bsp->descr;
  } else if (IS_Bioseq_set(sep)){
    bssp = sep->data.ptrvalue;
    list = bssp->descr;
  }
  if (list == NULL) return FALSE;
  while (list != NULL && !rval) {
    if ((list->choice == Seq_descr_pub && pub) ||
        (list->choice == Seq_descr_source && src) ||
        (list->choice == Seq_descr_comment && com)) {
       if (PubSrcComHasSubstring (omtp, list, dscp)) {
        rval = TRUE;
      }
    }
    list = list->next;
  }
  return rval;
}

static void 
PropagatePubSrcComDescriptors 
(SeqEntryPtr sep,
 Boolean pub,
 Boolean src,
 Boolean com,
 ObjMgrTypePtr omtp,
 DescStringCheckPtr dscp
)
{
  ValNodePtr   sdp_list = NULL;
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
  SeqEntryPtr  seqentry;

  if (sep == NULL || ! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr)sep->data.ptrvalue;
  if (bssp == NULL || bssp->descr == NULL) return;
 
  ExtractPubSrcComDescs (sep, &sdp_list, pub, src, com, omtp, dscp);
  
  if (sdp_list == NULL) return;
  for (seqentry = bssp->seq_set; seqentry != NULL; seqentry =seqentry->next) {
    if (IS_Bioseq_set (seqentry)) {
      bssp = seqentry->data.ptrvalue;
      ValNodeLink (&(bssp->descr),
                   AsnIoMemCopy ((Pointer) sdp_list,
                                 (AsnReadFunc) SeqDescrAsnRead,
                                 (AsnWriteFunc) SeqDescrAsnWrite));
    } else if (IS_Bioseq (seqentry)){
      bsp = (BioseqPtr) seqentry->data.ptrvalue;
      ValNodeLink (&(bsp->descr),
                   AsnIoMemCopy ((Pointer) sdp_list,
                                 (AsnReadFunc) SeqDescrAsnRead,
                                 (AsnWriteFunc) SeqDescrAsnWrite));
    }
  }
  SeqDescrFree (sdp_list);
}

NLM_EXTERN Boolean 
ConvertPubSrcComDescsToFeats 
(SeqEntryPtr sep,
 Boolean pub,
 Boolean src,
 Boolean com,
 Boolean toProts,
 Boolean PNTR asked_about_prop,
 Boolean PNTR propagate_descriptors,
 CharPtr findstring
 )

{
  BioseqSetPtr        bssp;
  ValNodePtr          head;
  Boolean             rsult;
  ValNodePtr          sdp;
  SeqEntryPtr         set_sep;
  MsgAnswer           ans;
  DescStringCheckData dscd;
  AsnExpOptPtr        aeop;
  ObjMgrPtr           omp;
  ObjMgrTypePtr       omtp = NULL;


  rsult = FALSE;
  if (! (pub || src || com)) return FALSE;
  if (sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (findstring != NULL && ! StringHasNoText (findstring)) {
    omp = ObjMgrGet();
    if (omp != NULL) {
      omtp = ObjMgrTypeFind(omp, OBJ_SEQDESC, NULL, NULL);
    }
  }

  if (findstring != NULL && ! StringHasNoText (findstring) && omtp != NULL) {
    dscd.aip = AsnIoNullOpen ();
    aeop = AsnExpOptNew (dscd.aip, NULL, NULL, AsnWriteConvertForDCallBack);
    if (aeop != NULL) {
      aeop->user_data = (Pointer) &dscd;
    }
    dscd.findstring = findstring;
  } else {
    dscd.aip = NULL;
    dscd.findstring = NULL;
  }
  dscd.stringfound = FALSE;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (HasPubSrcComDescriptors (sep, pub, src, com, omtp, &dscd)){
      if (asked_about_prop != NULL 
          && *asked_about_prop == FALSE
          && propagate_descriptors != NULL
          && *propagate_descriptors == FALSE) {
        ans = Message (MSG_YN, "Do you want to propagate descriptors on sets so that they can be converted?");
        *asked_about_prop = TRUE;
        if (ans == ANS_YES) {
          *propagate_descriptors = TRUE;
        }
      }
      if (propagate_descriptors != NULL && *propagate_descriptors) {
        PropagatePubSrcComDescriptors (sep, pub, src, com, omtp, &dscd);
      }
    }
    for (set_sep = bssp->seq_set; set_sep != NULL; set_sep = set_sep->next) {
      if (ConvertPubSrcComDescsToFeats (set_sep, pub, src, com, toProts && ! pub, asked_about_prop, propagate_descriptors, findstring)) {
        rsult = TRUE;
      }
    }
    if (dscd.aip != NULL) {
      AsnIoClose (dscd.aip);  
      dscd.aip = NULL;
    }
    return rsult;
  }
  head = NULL;
  ExtractPubSrcComDescs (sep, &head, pub, src, com, omtp, &dscd);
  rsult = (head != NULL);
  if (toProts) {
    BioseqExplore (sep, head, ConvertToFeatsOnNucsAndProts);
  } else {
    BioseqExplore (sep, head, ConvertToFeatsOnNucsOnly);
  }
  for (sdp = head; sdp != NULL; sdp = sdp->next) {
    switch (sdp->choice) {
      case Seq_descr_pub :
        PubdescFree ((PubdescPtr) sdp->data.ptrvalue);
        break;
      case Seq_descr_source :
        BioSourceFree ((BioSourcePtr) sdp->data.ptrvalue);
        break;
      case Seq_descr_comment :
        MemFree (sdp->data.ptrvalue);
        break;
      default :
        break;
    }
  }
  ValNodeFree (head);
  if (dscd.aip != NULL) {
    AsnIoClose (dscd.aip);  
    dscd.aip = NULL;
  }
  return rsult;
}

/* from Colombe */

static CharPtr sqn_string_complement (CharPtr str)
{
  CharPtr strp;

  for (strp = str; *strp != '\0'; strp++) {
         if (*strp == 'A') *strp = 'T';
         else if (*strp == 'T') *strp = 'A';
         else if (*strp == 'C') *strp = 'G';
         else if (*strp == 'G') *strp = 'C';
  }
  *strp = '\0';
  return str;
}

static CharPtr sqn_string_reverse (CharPtr str)
{
  Char    car;
  Int4    j;
  Int4    k;

  j = 0;
  k = StringLen (str) - 1;
  while (j < k) {
    car = str[j]; str[j] = str[k]; str[k] = car;
    j++;
    k--;
  }
  return str;
}

static Int4 sqn_ReadBufferFromSep (SeqPortPtr spp, CharPtr buffer, Int4 from, Int4 to, Int4 buffsegstart)
{
  Uint1    residue;
  Int4     k;
  Int4     pos;

  SeqPortSeek (spp, from, SEEK_SET);
  k = buffsegstart;
  pos = from;
  residue = SeqPortGetResidue(spp);
  while (pos < to && residue != SEQPORT_EOF)
  {
    if ( ! IS_residue(residue)) {
      /*
      switch (residue)
      {  
           case SEQPORT_VIRT:
              Message(MSG_OK,"SEQPORT_VIRT [%d=%c] at %ld\n", (int)residue, (char)residue, (long)pos);
              break;
           case SEQPORT_EOS:
              Message(MSG_OK,"[EOS]\n");
              break;
           default:
              Message(MSG_OK,"unknown char\n");
              break;
      }  
      pos++;
      */
    } else {
      buffer[k] = (Char) residue;
      k++;  
      pos++;
    }
    residue = SeqPortGetResidue(spp);
  }
  buffer[k] = '\0';
  return k;
}
 
static CharPtr sqn_load_seq_data (SeqIdPtr sip, Int4 from, Int4 to, Boolean is_prot, Int4 *lenp)
{
  BioseqPtr        bsp;
  SeqLocPtr        slp;
  SeqPortPtr       spp;
  CharPtr          str = NULL;
  Int4             lens;

  if (from > -1 && to > -1 && from >= to)
     return NULL;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     if (from < 0 || from > bsp->length -1)
        from = 0;
     if (to < 0 || to > bsp->length -1)
        to = bsp->length -1;
     BioseqUnlock (bsp);
     slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
     if (is_prot)
        spp = SeqPortNewByLoc (slp, Seq_code_ncbistdaa);
     else
        spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
     if (spp != NULL) {
        str = MemNew ((to-from+4) * sizeof(Char));
        lens = sqn_ReadBufferFromSep (spp, str, 0, to -from +1, 0);
        SeqPortFree (spp);
        if (lenp != NULL)
           *lenp = lens;
     }   
     SeqLocFree (slp);
  }
  return str;
}

static Int4 getlengthforid (SeqIdPtr sip)
{
  BioseqPtr        bsp;
  Int4             lens=0;

  if (sip==NULL)
     return 0;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     lens = bsp->length;
     BioseqUnlock (bsp);
  }
  return lens;
}

NLM_EXTERN SeqLocPtr StringSearchInBioseq (SeqIdPtr sip, CharPtr sub)
{
  SeqLocPtr slp=NULL;
  CharPtr   strdb,
            strtmp;
  Int4      lenbsp,
            fromp, top;
  Int4      lensub,
            lens,
            offset,
            shiftRange,
            maxRange;
  Boolean   firstpass=TRUE;

  lensub = StringLen (sub);
  maxRange = (Int4) MAX ((Int4)1000, lensub);
  lenbsp = getlengthforid (sip);
  while (slp == NULL) 
  {
   fromp = 0;
   top = MIN ((Int4)(fromp+maxRange), lenbsp) -1;
   while (fromp <= lenbsp && slp == NULL)
   {
     strdb = sqn_load_seq_data (sip, fromp, top, FALSE, &lens); 
     if (strdb != NULL)
     {
        offset = 0;
        strtmp = StringISearch (strdb, sub);
        if (strtmp != NULL) {
           offset =(Int4)abs(abs((long)strdb)-abs((long)strtmp));
           offset += fromp;
           if (offset > 0) {
              if (firstpass)
                 slp = SeqLocIntNew (offset, offset+lensub-1, Seq_strand_plus, sip);
              else 
                 slp = SeqLocIntNew (offset, offset+lensub-1, Seq_strand_minus, sip);
           }
        }
        MemFree (strdb);
     }
     shiftRange = maxRange - lensub;
     fromp = fromp + shiftRange;
     top = MIN ((Int4)(fromp+maxRange), lenbsp);
   }
   if (!firstpass) {
      sub = sqn_string_complement (sub);
      sub = sqn_string_reverse (sub);
      break;
   }
   firstpass=FALSE;
   sub = sqn_string_complement (sub);
   sub = sqn_string_reverse (sub);
  }
  return slp;
}

/*****************************************************************************
*
*   SequinEntryList (sep, mydata, mycallback, index, indent)
*       traverses all Seq-entry nodes beginning with sep
*       calls mycallback () at each node
*       Does enter BioseqSets of _class "parts", but ignores the
*       parts set itself
*
*****************************************************************************/

NLM_EXTERN Int4 SequinEntryList (SeqEntryPtr sep, Pointer mydata,
                                 SeqEntryFunc mycallback,
                                 Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return index;
  if (IS_Bioseq (sep)) {
    if (mycallback != NULL)
      (*mycallback) (sep, mydata, index, indent);
    return index + 1;
  }
  /*
  if (Bioseq_set_class (sep) == 4) return index;
  index++;
  */
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  sep = bssp->seq_set;
  indent++;
  while (sep != NULL) {
    index = SequinEntryList (sep, mydata, mycallback, index, indent);
    sep = sep->next;
  }
  return index;
}

/* functions to parse [org=Drosophila melanogaster] and [gene=lacZ] from titles */

NLM_EXTERN SqnTagPtr SqnTagParse (CharPtr ttl)

{
  Int2       num_tags;
  CharPtr    ptr;
  CharPtr    query;
  SqnTagPtr  stp;
  CharPtr    str;
  CharPtr    tmp;

  if (StringHasNoText (ttl)) return NULL;
  stp = (SqnTagPtr) MemNew (sizeof (SqnTag));
  if (stp == NULL) return NULL;
  query = StringSave (ttl);

  str = StringChr (query, '[');
  for (num_tags = 0; num_tags < MAX_SQN_TAGS && str != NULL; num_tags++) {
    str++;
    stp->tag [num_tags] = str;
    ptr = StringChr (str, '=');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      stp->val [num_tags] = ptr;
      if (*ptr == '"') {
        ptr++;
        stp->val [num_tags] = ptr;
        tmp = StringChr (ptr, '"');
        if (tmp != NULL) {
          *tmp = '\0';
          tmp++;
          ptr = tmp;
        }
      }
      tmp = StringChr (ptr, ']');
      if (tmp != NULL) {
        *tmp = '\0';
        tmp++;
        str = StringChr (tmp, '[');
      }
      TrimSpacesAroundString (stp->tag [num_tags]);
      TrimSpacesAroundString (stp->val [num_tags]);
    }
  }

  stp->query = query;
  stp->num_tags = num_tags;

  return stp;
}

NLM_EXTERN SqnTagPtr SqnTagFree (SqnTagPtr stp)

{
  if (stp == NULL) return NULL;
  MemFree (stp->query);
  return MemFree (stp);
}

extern Boolean StringsAreEquivalent (CharPtr str1, CharPtr str2)

{
  Char  ch1, ch2;

  if (StringHasNoText (str1) && StringHasNoText (str2)) return TRUE;
  if (StringHasNoText (str1) || StringHasNoText (str2)) return FALSE;

  ch1 = *str1;
  ch2 = *str2;
  while (ch1 != '\0' && ch2 != '\0') {
    if (TO_LOWER (ch1) != TO_LOWER (ch2)) {
      if ((ch1 != '-' && ch1 != '_' && ch1 != ' ') || (ch2 != '_' && ch2 != '-' && ch2 != ' ')) return FALSE;
    }
    str1++;
    str2++;
    ch1 = *str1;
    ch2 = *str2;
  }

  if (TO_LOWER (ch1) != TO_LOWER (ch2)) {
    if ((ch1 != '-' && ch1 != '_' && ch1 != ' ') || (ch2 != '_' && ch2 != '-' && ch2 != ' ')) return FALSE;
  }

  return TRUE;
}

NLM_EXTERN Uint1 EquivalentSubSourceEx (CharPtr str, Boolean allow_discouraged_and_discontinued)
{
  Int4    i;
  Uint1   subtype = 0;

  for (i = 0; current_subsource_subtype_alist [i].name != NULL && subtype == 0; i++) {
    if (StringsAreEquivalent (str, current_subsource_subtype_alist [i].name)) {
      subtype = current_subsource_subtype_alist [i].value;
    }
  }
  for (i = 0; subsource_aliases[i].name != NULL && subtype == 0; i++) {
    if (StringsAreEquivalent (str, subsource_aliases [i].alias)) {
      subtype = subsource_aliases [i].value;
    }
  }
  if (allow_discouraged_and_discontinued && subtype == 0) {
    for (i = 0; discouraged_subsource_subtype_alist [i].name != NULL && subtype == 0; i++) {
      if (StringsAreEquivalent (str, discouraged_subsource_subtype_alist [i].name)) {
        subtype = discouraged_subsource_subtype_alist [i].value;
      }
    }
    for (i = 0; discontinued_subsource_subtype_alist [i].name != NULL && subtype == 0; i++) {
      if (StringsAreEquivalent (str, discontinued_subsource_subtype_alist [i].name)) {
        subtype = discontinued_subsource_subtype_alist [i].value;
      }
    }
  }

  return subtype;
}


NLM_EXTERN Uint1 EquivalentSubSource (CharPtr str)
{
  return EquivalentSubSourceEx (str, FALSE);
}


NLM_EXTERN Uint1 EquivalentOrgModEx (CharPtr str, Boolean allow_discouraged_and_discontinued)
{
  Int4    i;
  Uint1   subtype = 0;

  for (i = 0; current_orgmod_subtype_alist [i].name != NULL && subtype == 0; i++) {
    if (StringsAreEquivalent (str, current_orgmod_subtype_alist [i].name)) {
      subtype = current_orgmod_subtype_alist [i].value;
    }
  }
  for (i = 0; orgmod_aliases[i].name != NULL && subtype == 0; i++) {
    if (StringsAreEquivalent (str, orgmod_aliases [i].alias)) {
      subtype = orgmod_aliases [i].value;
    }
  }
  if (allow_discouraged_and_discontinued && subtype == 0) {
    for (i = 0; discouraged_orgmod_subtype_alist [i].name != NULL && subtype == 0; i++) {
      if (StringsAreEquivalent (str, discouraged_orgmod_subtype_alist [i].name)) {
        subtype = discouraged_orgmod_subtype_alist [i].value;
      }
    }
    for (i = 0; discontinued_orgmod_subtype_alist [i].name != NULL && subtype == 0; i++) {
      if (StringsAreEquivalent (str, discontinued_orgmod_subtype_alist [i].name)) {
        subtype = discontinued_orgmod_subtype_alist [i].value;
      }
    }
  }

  return subtype;
}


NLM_EXTERN Uint1 EquivalentOrgMod (CharPtr str)
{
  return EquivalentOrgModEx (str, FALSE);
}


NLM_EXTERN CharPtr SqnTagFind (SqnTagPtr stp, CharPtr tag)

{
  Int2  i;

  if (stp == NULL || StringHasNoText (tag)) return NULL;
  for (i = 0; i < stp->num_tags; i++) {
    if (stp->tag [i] != NULL && StringsAreEquivalent (stp->tag [i], tag)) {
      stp->used [i] = TRUE;
      return stp->val [i];
    }
  }
  return NULL;
}


NLM_EXTERN ValNodePtr SqnTagFindMultiple (SqnTagPtr stp, CharPtr tag)

{
  Int2  i;
  ValNodePtr list = NULL;

  if (stp == NULL || StringHasNoText (tag)) return NULL;
  for (i = 0; i < stp->num_tags; i++) {
    if (stp->tag [i] != NULL && StringsAreEquivalent (stp->tag [i], tag)) {
      stp->used [i] = TRUE;
      ValNodeAddPointer (&list, 0, stp->val[i]);
    }
  }
  return list;
}


NLM_EXTERN CharPtr SqnTagFindUnused (SqnTagPtr stp, CharPtr tag)

{
  Int2  i;

  if (stp == NULL || StringHasNoText (tag)) return NULL;
  for (i = 0; i < stp->num_tags; i++) {
    if (stp->tag [i] != NULL && StringsAreEquivalent (stp->tag [i], tag) && !stp->used[i]) {
      stp->used [i] = TRUE;
      return stp->val [i];
    }
  }
  return NULL;
}


static void AddSqnTagToSubSource (BioSourcePtr biop, Uint1 subtype, CharPtr subname)
{
  SubSourcePtr ssp;
  if (biop == NULL || subtype == 0 || subname == NULL) return;

  for (ssp = biop->subtype;
       ssp != NULL && ssp->subtype != subtype;
       ssp = ssp->next) continue;
  if (ssp != NULL) {
    ssp->name = MemFree (ssp->name);
    ssp->name = StringSave (subname);
  } else {
    ssp = SubSourceNew ();
    if (ssp != NULL) {
      ssp->subtype = subtype;
      ssp->name = StringSave (subname);
      ssp->next = biop->subtype;
      biop->subtype = ssp;
    }
  }
}

static void SqnTagFindSubSourceQuals (SqnTagPtr stp, BioSourcePtr biop)
{
  Int4    i, j;
  CharPtr str = NULL;

  for (i = 0; current_subsource_subtype_alist [i].name != NULL; i++) {
    str = SqnTagFind (stp, current_subsource_subtype_alist [i].name);
    if (str == NULL) {
      for (j = 0; subsource_aliases[j].name != NULL && str == NULL; j++) {
        if (subsource_aliases[j].value == current_subsource_subtype_alist[i].value) {
          str = SqnTagFind (stp, subsource_aliases [j].alias);
        }
      }
    }
    if (str != NULL) {
      AddSqnTagToSubSource (biop, current_subsource_subtype_alist[i].value, str);
    }
  }
}

static void AddSqnTagToOrgMod (OrgNamePtr onp, Uint1 subtype, CharPtr subname)
{
  OrgModPtr omp;
  if (onp == NULL || subname == NULL || subtype == 0) return;

  for (omp = onp->mod;
       omp != NULL && omp->subtype != subtype;
       omp = omp->next) continue;
  if (omp != NULL) {
    omp->subname = MemFree (omp->subname);
    omp->subname = StringSave (subname);
  } else {
    omp = OrgModNew ();
    if (omp != NULL) {
      omp->subtype = subtype;
      omp->subname = StringSave (subname);
      omp->next = onp->mod;
      onp->mod = omp;
    }
  }
}

static void SqnTagFindOrgModQuals (SqnTagPtr stp, OrgNamePtr onp)
{
  Int4    i, j;
  CharPtr str;

  for (i = 0; current_orgmod_subtype_alist [i].name != NULL; i++) {
    str = SqnTagFind (stp, current_orgmod_subtype_alist [i].name);
    if (str == NULL) {
      for (j = 0; orgmod_aliases[j].name != NULL && str == NULL; j++) {
        if (orgmod_aliases[j].value == current_orgmod_subtype_alist[i].value) {
          str = SqnTagFind (stp, orgmod_aliases [j].alias);
        }
      }
    }
    if (str != NULL) {
      AddSqnTagToOrgMod (onp, current_orgmod_subtype_alist[i].value, str);
    }
  }
}

/* functions to extract BioSource, MolInfo, and Bioseq information from parsed titles */

static CharPtr sqntag_biosrc_genome_list [] = {
  "?", "genomic", "chloroplast", "chromoplast", "kinetoplast",
  "mitochondrion", "plastid", "macronuclear", "extrachromosomal",
  "plasmid", "transposon", "insertion sequence", "cyanelle",
  "proviral", "virion", "nucleomorph", "apicoplast", "leucoplast",
  "proplastid", "endogenous-virus", "hydrogenosome", "chromosome",
  "chromatophore", NULL
};

static CharPtr sqntag_biosrc_origin_list [] = {
  "?", "natural", "natural mutant", "mutant", "artificial",
  "synthetic", "other", NULL
};


static void SqnTagParsePrimers (SqnTagPtr stp, BioSourcePtr biop)
{
  ValNode quals[4];
  Int4    qual_types[] = { SUBSRC_fwd_primer_name, SUBSRC_fwd_primer_seq, SUBSRC_rev_primer_name, SUBSRC_rev_primer_seq};
  Int4    qual_defs[] = { Source_qual_fwd_primer_name, Source_qual_fwd_primer_seq, Source_qual_rev_primer_name, Source_qual_rev_primer_seq};
  Int4 num_quals = 4, qual;
  Int4 i, j;

  if (stp == NULL || stp->num_tags == 0 || biop == NULL) return;

  for (i = 0; i < num_quals; i++) {
    MemSet (quals + i, 0, sizeof (ValNode));
    quals[i].choice = SourceQualChoice_textqual;
    quals[i].data.intvalue = qual_defs[i];
  }

  for (i = 0; i < stp->num_tags; i++) {
    if (stp->tag [i] != NULL) {
      qual = EquivalentSubSourceEx (stp->tag[i], TRUE);
      for (j = 0; j < num_quals; j++) {
        if (qual == qual_types[j]) {
          stp->used [i] = TRUE;
          SetSourceQualInBioSource (biop, quals + j, NULL, stp->val[i], ExistingTextOption_add_qual);
          break;
        }
      }
    }
  }

}


NLM_EXTERN BioSourcePtr ParseTitleIntoBioSource (
  SqnTagPtr stp,
  CharPtr organism,
  BioSourcePtr biop
)

{
  DbtagPtr      db;
  Int2          i;
  ObjectIdPtr   oip;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  CharPtr       ptr;
  SubSourcePtr  ssp;
  CharPtr       str;
  int           val;
  ValNodePtr    vnp, list, list_vnp;


  if ((stp == NULL || stp->num_tags == 0) && StringHasNoText (organism)) return biop;

  if (biop == NULL) {
    biop = BioSourceNew ();
    if (biop == NULL) return biop;
  }
  if (biop->org == NULL) {
    biop->org = OrgRefNew ();
  }
  orp = biop->org;
  if (orp->orgname == NULL) {
    orp->orgname = OrgNameNew ();
  }
  onp = orp->orgname;

  str = SqnTagFind (stp, "organism");
  if (str == NULL) {
    str = SqnTagFind (stp, "org");
  }
  if (organism == NULL) {
    organism = str;
  }
  if (! StringHasNoText (organism)) {
    if (StringICmp (orp->taxname, organism) != 0) {

      /* if command line or fasta defline organism doesn't match, clear template */

      biop->org = OrgRefFree (biop->org);
      biop->subtype = SubSourceFree (biop->subtype);

      /* then recreate orgref and orgname structures, save organism name */

      biop->org = OrgRefNew ();
      orp = biop->org;
      orp->orgname = OrgNameNew ();
      onp = orp->orgname;

      orp->taxname = StringSave (organism);
    }
  }

  if (stp == NULL) return biop;

  str = SqnTagFind (stp, "location");
  if (str != NULL) {
    if (StringICmp (str, "mitochondrial") == 0) {
      str = "mitochondrion";
    } else if (StringICmp (str, "provirus") == 0) {
      str = "proviral";
    }
    for (i = 0; sqntag_biosrc_genome_list [i] != NULL; i++) {
      if (StringsAreEquivalent (str, sqntag_biosrc_genome_list [i])) {
        biop->genome = (Uint1) i;
      }
    }
  }

  str = SqnTagFind (stp, "origin");
  if (str != NULL) {
    for (i = 0; sqntag_biosrc_origin_list [i] != NULL; i++) {
      if (StringsAreEquivalent (str, sqntag_biosrc_origin_list [i])) {
        biop->origin = (Uint1) i;
      }
    }
    if (biop->origin == 6) {
      biop->origin = 255;
    }
  }

  SqnTagFindOrgModQuals (stp, onp);

  SqnTagFindSubSourceQuals (stp, biop);

  SqnTagParsePrimers (stp, biop);

  list = SqnTagFindMultiple (stp, "db_xref");
  for (list_vnp = list; list_vnp != NULL; list_vnp = list_vnp->next) {
    str = list_vnp->data.ptrvalue;
    vnp = ValNodeNew (NULL);
    db = DbtagNew ();
    vnp->data.ptrvalue = db;
    ptr = StringChr (str, ':');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      db->db = StringSave (str);
      oip = ObjectIdNew ();
      oip->str = StringSave (ptr);
      db->tag = oip;
    } else {
      db->db = StringSave ("?");
      oip = ObjectIdNew ();
      oip->str = StringSave (str);
      db->tag = oip;
    }
    vnp->next = orp->db;
    orp->db = vnp;
  }
  list = ValNodeFree (list);

  str = SqnTagFind (stp, "division");
  if (str == NULL) {
    str = SqnTagFind (stp, "div");
  }
  if (str != NULL) {
    onp->div = MemFree (onp->div);
    onp->div = StringSave (str);
  }

  str = SqnTagFind (stp, "lineage");
  if (str != NULL) {
    onp->lineage = MemFree (onp->lineage);
    onp->lineage = StringSave (str);
  }

  str = SqnTagFind (stp, "gcode");
  if (str != NULL && sscanf (str, "%d", &val) == 1) {
    onp->gcode = (Uint1) val; /* cytoplasmic */
  }

  str = SqnTagFind (stp, "mgcode");
  if (str != NULL && sscanf (str, "%d", &val) == 1) {
    onp->mgcode = (Uint1) val; /* mitochondrial */
  }

  str = SqnTagFind (stp, "pgcode");
  if (str != NULL && sscanf (str, "%d", &val) == 1) {
    onp->pgcode = (Uint1) val; /* plastid */
  }

  str = SqnTagFind (stp, "note");
  if (str == NULL) {
    str = SqnTagFind (stp, "notes");
  }
  if (str != NULL) {
    ssp = SubSourceNew ();
    if (ssp != NULL) {
      ssp->subtype = (Uint1) SUBSRC_other;
      ssp->name = StringSave (str);
      ssp->next = biop->subtype;
      biop->subtype = ssp;
    }
  }

  str = SqnTagFind (stp, "focus");
  if (str != NULL) {
    if (StringICmp (str, "TRUE") == 0) {
      biop->is_focus = TRUE;
    }
  }

  return biop;
}

static CharPtr molinfo_biomol_list [] = {
  "?", "genomic", "precursor RNA", "mRNA", "rRNA", "tRNA", "snRNA",
  "scRNA", "peptide", "other-genetic", "genomic-mRNA", "cRNA", "snoRNA",
  "transcribed RNA", "non-coding RNA", "transfer-messenger RNA", NULL
};

static CharPtr molinfo_completeness_list [] = {
  "unknown", "complete", "partial", "no-left", "no-right", "no-ends", "has-left", "has-right", NULL
};

NLM_EXTERN void ReadTechFromString (CharPtr str, MolInfoPtr mip)
{
  Int4 i;
  
  if (mip == NULL || str == NULL)
  {
    return;
  }

  i = TechFromTechName (str);
  if (i > -1) {
    mip->tech = (Uint1) i;
  }
}

NLM_EXTERN void ReadCompletenessFromString (CharPtr str, MolInfoPtr mip)
{
  Int4 i;
  
  if (mip == NULL || str == NULL)
  {
    return;
  }
  
  for (i = 0; molinfo_completeness_list [i] != NULL; i++) {
    if (StringsAreEquivalent (str, molinfo_completeness_list [i])) {
      mip->completeness = (Uint1) i;
    }
  }
}

NLM_EXTERN MolInfoPtr ParseTitleIntoMolInfo (
  SqnTagPtr stp,
  MolInfoPtr mip
)

{
  Int2     i;
  CharPtr  str;

  if (stp == NULL) return mip;

  if (mip == NULL) {
    mip = MolInfoNew ();
    if (mip == NULL) return mip;
  }

  str = SqnTagFind (stp, "moltype");
  if (str == NULL) {
    str = SqnTagFind (stp, "mol-type");
  }
  if (str == NULL) {
    str = SqnTagFind (stp, "mol_type");
  }
  if (str != NULL) {
    for (i = 0; molinfo_biomol_list [i] != NULL; i++) {
      if (StringsAreEquivalent (str, molinfo_biomol_list [i])) {
        mip->biomol = (Uint1) i;
      }
    }
  }

  str = SqnTagFind (stp, "tech");
  ReadTechFromString (str, mip);

  str = SqnTagFind (stp, "completeness");
  if (str == NULL) {
    str = SqnTagFind (stp, "completedness");
  }
  ReadCompletenessFromString (str, mip);

  return mip;
}

NLM_EXTERN BioseqPtr ParseTitleIntoBioseq (
  SqnTagPtr stp,
  BioseqPtr bsp
)

{
  CharPtr  str;

  if (stp == NULL || bsp == NULL) return bsp;

  str = SqnTagFind (stp, "topology");
  if (str == NULL) {
    str = SqnTagFind (stp, "top");
  }
  if (str != NULL) {
    if (StringICmp (str, "linear") == 0) {
      bsp->topology = TOPOLOGY_LINEAR;
    } else if (StringICmp (str, "circular") == 0) {
      bsp->topology = TOPOLOGY_CIRCULAR;
    }
  }

  str = SqnTagFind (stp, "molecule");
  if (str == NULL) {
    str = SqnTagFind (stp, "mol");
  }
  if (str != NULL) {
    if (StringICmp (str, "dna") == 0) {
      bsp->mol = Seq_mol_dna;
    } else if (StringICmp (str, "rna") == 0) {
      bsp->mol = Seq_mol_rna;
    }
  }

  str = SqnTagFind (stp, "strand");
  if (str != NULL) {
    if (StringICmp (str, "single") == 0) {
      bsp->strand = 1;
    } else if (StringICmp (str, "double") == 0) {
      bsp->strand = 2;
    } else if (StringICmp (str, "mixed") == 0) {
      bsp->strand = 3;
    }
  }

  return bsp;
}

NLM_EXTERN GeneRefPtr ParseTitleIntoGeneRef (
  SqnTagPtr stp,
  GeneRefPtr grp
)

{
  CharPtr  str;

  if (stp == NULL || grp == NULL) return grp;

  str = SqnTagFind (stp, "gene");
  if (str != NULL) {
    grp->locus = StringSave (str);
  }

  str = SqnTagFind (stp, "allele");
  if (str != NULL) {
    grp->allele = StringSave (str);
  }

  str = SqnTagFind (stp, "gene_syn");
  if (str == NULL) {
    str = SqnTagFind (stp, "gene_synonym");
  }
  if (str != NULL) {
    ValNodeCopyStr (&(grp->syn), 0, str);
  }

  str = SqnTagFind (stp, "locus_tag");
  if (str != NULL) {
    grp->locus_tag = StringSave (str);
  }

  return grp;
}

NLM_EXTERN ProtRefPtr ParseTitleIntoProtRef (
  SqnTagPtr stp,
  ProtRefPtr prp
)

{
  CharPtr  str;

  if (stp == NULL || prp == NULL) return prp;

  str = SqnTagFind (stp, "protein");
  if (str == NULL) {
    str = SqnTagFind (stp, "prot");
  }
  if (str != NULL) {
    ValNodeCopyStr (&(prp->name), 0, str);
  }

  str = SqnTagFind (stp, "prot_desc");
  if (str != NULL) {
    prp->desc = StringSave (str);
  }

  str = SqnTagFind (stp, "EC_number");
  if (str != NULL) {
    ValNodeCopyStr (&(prp->ec), 0, str);
  }

  str = SqnTagFind (stp, "activity");
  if (str == NULL) {
    str = SqnTagFind (stp, "function");
  }
  if (str != NULL) {
    ValNodeCopyStr (&(prp->activity), 0, str);
  }

  return prp;
}

static Boolean ParseAccessionRange (
  CharPtr accn,
  CharPtr prefix,
  Int4Ptr startp,
  Int4Ptr stopp,
  Int2Ptr digitsp
)

{
  Char      ch;
  Int2      digits;
  CharPtr   ptr, tmp;
  Int4      start, stop;
  long int  val;

  if (StringHasNoText (accn)) return FALSE;
  if (prefix == NULL || startp == NULL || stopp == NULL || digitsp == NULL) return FALSE;

  ptr = accn;
  ch = *ptr;
  while (IS_ALPHA (ch)) {
    *prefix = ch;
    prefix++;
    ptr++;
    ch = *ptr;
  }
  *prefix = '\0';

  tmp = StringChr (ptr, '-');
  if (tmp == NULL) return FALSE;
  *tmp = '\0';
  tmp++;

  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  start = (Int4) val;

  digits = 0;
  while (IS_DIGIT (ch)) {
    digits++;
    ptr++;
    ch = *ptr;
  }

  ptr = tmp;
  ch = *ptr;
  while (IS_ALPHA (ch)) {
    ptr++;
    ch = *ptr;
  }

  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  stop = (Int4) val;

  *startp = start;
  *stopp = stop;
  *digitsp = digits;

  return TRUE;
}

static void DoAddToSecAccn (
  GBBlockPtr gbp,
  CharPtr accn
)

{
  Int2  digits, j;
  Int4  idx;
  Char  numbers [32];
  Char  prefix [16];
  Int4  start, stop;
  Char  tmp [64];

  if (StringChr (accn, '-') != NULL) {
    if (ParseAccessionRange (accn, prefix, &start, &stop, &digits)) {
      for (idx = start; idx <= stop; idx++) {
        sprintf (numbers, "%*ld", digits, (long) idx);
        for (j = 0; j < digits && numbers [j] != '\0'; j++) {
          if (numbers [j] == ' ') {
            numbers [j] = '0';
          }
        }
        StringCpy (tmp, prefix);
        StringCat (tmp, numbers);
        ValNodeCopyStr (&(gbp->extra_accessions), 0, tmp);
      }
    }
  } else {
    ValNodeCopyStr (&(gbp->extra_accessions), 0, accn);
  }
}

NLM_EXTERN GBBlockPtr ParseTitleIntoGenBank (
  SqnTagPtr stp,
  GBBlockPtr gbp
)

{
  Char     ch;
  CharPtr  last;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (stp == NULL) return gbp;

  if (gbp == NULL) {
    gbp = GBBlockNew ();
    if (gbp == NULL) return gbp;
  }

  str = SqnTagFind (stp, "secondary-accession");
  if (str == NULL) {
    str = SqnTagFind (stp, "secondary-accessions");
  }
  if (str != NULL) {
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',') {
        *ptr = '\0';
        if (! StringHasNoText (last)) {
          TrimSpacesAroundString (last);
          DoAddToSecAccn (gbp, last);
        }
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      DoAddToSecAccn (gbp, last);
    }
    MemFree (tmp);
  }

  str = SqnTagFind (stp, "keyword");
  if (str == NULL) {
    str = SqnTagFind (stp, "keywords");
  }
  if (str != NULL) {
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',' || ch == ';') {
        *ptr = '\0';
        if (! StringHasNoText (last)) {
          TrimSpacesAroundString (last);
          ValNodeCopyStr (&(gbp->keywords), 0, last);
        }
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      ValNodeCopyStr (&(gbp->keywords), 0, last);
    }
    MemFree (tmp);
  }

  return gbp;
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


static void AddStringToSeqHist (
  SeqHistPtr shp,
  CharPtr str
)

{
  Char          prefix [20];
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;
  Uint4         whichdb;

  if (shp == NULL || StringHasNoText (str)) return;
  sip = ValNodeAdd (&(shp->replace_ids));
  if (sip == NULL) return;

  if (StringNICmp (str, "other|", 6) == 0) {
    tsip = TextSeqIdNew ();
    tsip->accession = StringSave (str + 6);
    sip->data.ptrvalue = tsip;
    sip->choice = SEQID_OTHER;
  } else if (StringNICmp (str, "ref_seq|", 8) == 0) {
    tsip = TextSeqIdNew ();
    tsip->accession = StringSave (str + 8);
    sip->data.ptrvalue = tsip;
    sip->choice = SEQID_OTHER;
  } else if (StringNICmp (str, "gi|", 3) == 0 && IsAllDigits (str + 3)) {
    sip->data.intvalue = atoi (str + 3);
    sip->choice = SEQID_GI;
  } else if (IsAllDigits (str)) {
    sip->data.intvalue = atoi (str);
    sip->choice = SEQID_GI;
  } else {
    tsip = TextSeqIdNew ();
    StringNCpy_0 (prefix, str, sizeof (prefix));
    whichdb = WHICH_db_accession (prefix);
    if (ACCN_IS_EMBL (whichdb)) {
      sip->choice = SEQID_EMBL;
    } else if (ACCN_IS_DDBJ (whichdb)) {
      sip->choice = SEQID_DDBJ;
    } else {
      sip->choice = SEQID_GENBANK;
    }
    sip->data.ptrvalue = (Pointer) tsip;
    tsip->accession = StringSave (str);
  }
}

static void DoAddToSeqHist (
  SeqHistPtr shp,
  CharPtr accn
)

{
  Int2  digits, j;
  Int4  idx;
  Char  numbers [32];
  Char  prefix [16];
  Int4  start, stop;
  Char  tmp [64];

  if (StringChr (accn, '-') != NULL) {
    if (ParseAccessionRange (accn, prefix, &start, &stop, &digits)) {
      for (idx = start; idx <= stop; idx++) {
        sprintf (numbers, "%*ld", digits, (long) idx);
        for (j = 0; j < digits && numbers [j] != '\0'; j++) {
          if (numbers [j] == ' ') {
            numbers [j] = '0';
          }
        }
        StringCpy (tmp, prefix);
        StringCat (tmp, numbers);
        AddStringToSeqHist (shp, tmp);
      }
    }
  } else {
    AddStringToSeqHist (shp, accn);
  }
}

NLM_EXTERN SeqHistPtr ParseStringIntoSeqHist (
  SeqHistPtr shp,
  CharPtr str
)

{
  Char     ch;
  CharPtr  last;
  CharPtr  ptr;
  CharPtr  tmp;

  if (shp == NULL) {
    shp = SeqHistNew ();
    if (shp == NULL) return shp;
  }

  if (str != NULL) {
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',') {
        *ptr = '\0';
        if (! StringHasNoText (last)) {
          TrimSpacesAroundString (last);
          DoAddToSeqHist (shp, last);
        }
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      DoAddToSeqHist (shp, last);
    }
    MemFree (tmp);
  }

  return shp;
}

NLM_EXTERN SeqHistPtr ParseTitleIntoSeqHist (
  SqnTagPtr stp,
  SeqHistPtr shp
)

{
  Char     ch;
  CharPtr  last;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (stp == NULL) return shp;

  if (shp == NULL) {
    shp = SeqHistNew ();
    if (shp == NULL) return shp;
  }

  str = SqnTagFind (stp, "secondary-accession");
  if (str == NULL) {
    str = SqnTagFind (stp, "secondary-accessions");
  }
  if (str != NULL) {
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',') {
        *ptr = '\0';
        if (! StringHasNoText (last)) {
          TrimSpacesAroundString (last);
          DoAddToSeqHist (shp, last);
        }
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      DoAddToSeqHist (shp, last);
    }
    MemFree (tmp);
  }

  return shp;
}

NLM_EXTERN void ParseTitleIntoSubmitBlock (
  SqnTagPtr stp,
  SubmitBlockPtr sbp
)

{
  DatePtr  dp;
  CharPtr  str;

  if (stp == NULL || sbp == NULL) return;

  str = SqnTagFind (stp, "hup");
  if (str != NULL) {
    sbp->hup = FALSE;
    sbp->reldate = DateFree (sbp->reldate);
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "y") == 0) {
        sbp->hup = TRUE;
        dp = DateCurr ();
        sbp->reldate = dp;
        if (dp != NULL) {
          if (dp->data [0] == 1) {
            (dp->data [1])++;
          }
        }
      } else {
        dp = DateParse (str);
        if (dp != NULL) {
          sbp->hup = TRUE;
          sbp->reldate = dp;
        }
      }
    }
  }
}

NLM_EXTERN UserObjectPtr ParseTitleIntoTpaAssembly (
  SqnTagPtr stp,
  UserObjectPtr uop
)

{
  Char     ch;
  CharPtr  last;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (stp == NULL) return uop;

  if (uop == NULL) {
    uop = CreateTpaAssemblyUserObject ();
    if (uop == NULL) return uop;
  }

  str = SqnTagFind (stp, "primary");
  if (str == NULL) {
    str = SqnTagFind (stp, "primary_accessions");
  }
  if (str != NULL) {
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',') {
        *ptr = '\0';
        if (! StringHasNoText (last)) {
          TrimSpacesAroundString (last);
          AddAccessionToTpaAssemblyUserObject (uop, last, 0, 0);
        }
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (! StringHasNoText (last)) {
      TrimSpacesAroundString (last);
      AddAccessionToTpaAssemblyUserObject (uop, last, 0, 0);
    }
    MemFree (tmp);
  }

  return uop;
}

NLM_EXTERN UserObjectPtr ParseTitleIntoGenomeProjectsDB (
  SqnTagPtr stp,
  UserObjectPtr uop
)

{
  Char      ch;
  CharPtr   last;
  CharPtr   ptr;
  Int4      projectID;
  CharPtr   str;
  CharPtr   tmp;
  long int  val;

  if (stp == NULL) return uop;

  if (uop == NULL) {
    uop = CreateGenomeProjectsDBUserObject ();
    if (uop == NULL) return uop;
  }

  str = SqnTagFind (stp, "project");
  if (str == NULL) {
    str = SqnTagFind (stp, "projects");
  }
  if (str != NULL) {
    tmp = StringSave (str);
    last = tmp;
    ptr = last;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ',' || ch == ';') {
        *ptr = '\0';
        if (StringDoesHaveText (last)) {
          TrimSpacesAroundString (last);
          if (sscanf (last, "%ld", &val) == 1 && val > 0) {
            projectID = (Int4) val;
            AddIDsToGenomeProjectsDBUserObject (uop, projectID, 0);
          }
        }
        ptr++;
        last = ptr;
        ch = *ptr;
      } else {
        ptr++;
        ch = *ptr;
      }
    }
    if (StringDoesHaveText (last)) {
      TrimSpacesAroundString (last);
      if (sscanf (last, "%ld", &val) == 1 && val > 0) {
        projectID = (Int4) val;
        AddIDsToGenomeProjectsDBUserObject (uop, projectID, 0);
      }
    }
    MemFree (tmp);
  }

  return uop;
}

static ValNodePtr ParseCommaStringList (
  CharPtr str
)

{
  Char        ch;
  ValNodePtr  head = NULL, tail = NULL;
  CharPtr     last;
  CharPtr     ptr;
  CharPtr     tmp;

  if (StringHasNoText (str)) return NULL;

  tmp = StringSave (str);
  if (tmp == NULL) return NULL;

  last = tmp;
  ptr = last;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == ',' || ch == ';') {
      *ptr = '\0';
      if (StringDoesHaveText (last)) {
        TrimSpacesAroundString (last);
        ValNodeCopyStrEx (&head, &tail, 0, last);
      }
      ptr++;
      last = ptr;
      ch = *ptr;
    } else {
      ptr++;
      ch = *ptr;
    }
  }
  if (StringDoesHaveText (last)) {
    TrimSpacesAroundString (last);
    ValNodeCopyStrEx (&head, &tail, 0, last);
  }

  MemFree (tmp);

  return head;
}


NLM_EXTERN void AddFieldStringToDbLinkUserObject (
  CharPtr str,
  CharPtr field_name,
  UserObjectPtr uop
)
{
  CharPtr PNTR     cpp;
  ValNodePtr       head, vnp;
  Int4             i, num;

  head = ParseCommaStringList (str);
  if (head != NULL) {
    num = 0;
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      num++;
    }
    if (num > 0) {
      cpp = (CharPtr PNTR) MemNew (sizeof (CharPtr) * num);
      if (cpp != NULL) {
        i = 0;
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          cpp [i] = str;
          i++;
        }
        if (i > 0) {
          AddFieldsToDBLinkUserObject (uop, field_name, i, cpp);
        }
      }
      MemFree (cpp);
    }
  }
  ValNodeFreeData (head);
}


NLM_EXTERN UserObjectPtr ParseTitleIntoDBLinkBioProject (
  SqnTagPtr stp,
  UserObjectPtr uop
)

{
  CharPtr          str;

  if (stp == NULL) return uop;

  if (uop == NULL) {
    uop = CreateDBLinkUserObject ();
    if (uop == NULL) return uop;
  }

  str = SqnTagFind (stp, "bioproject");
  if (str == NULL) {
    str = SqnTagFind (stp, "bioprojects");
  }
  if (str == NULL) return uop;

  AddFieldStringToDbLinkUserObject (str, "BioProject", uop);

  return uop;
}

NLM_EXTERN UserObjectPtr ParseTitleIntoDBLinkBioSample (
  SqnTagPtr stp,
  UserObjectPtr uop
)

{
  CharPtr PNTR     cpp;
  ValNodePtr       head, vnp;
  Int4             i, num;
  CharPtr          str;

  if (stp == NULL) return uop;

  if (uop == NULL) {
    uop = CreateDBLinkUserObject ();
    if (uop == NULL) return uop;
  }

  str = SqnTagFind (stp, "biosample");
  if (str == NULL) {
    str = SqnTagFind (stp, "biosamples");
  }
  if (str == NULL) return uop;

  head = ParseCommaStringList (str);
  if (head != NULL) {
    num = 0;
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      num++;
    }
    if (num > 0) {
      cpp = (CharPtr PNTR) MemNew (sizeof (CharPtr) * num);
      if (cpp != NULL) {
        i = 0;
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          cpp [i] = str;
          i++;
        }
        if (i > 0) {
          AddBioSampleIDsToDBLinkUserObject (uop, i, cpp);
        }
      }
      MemFree (cpp);
    }
  }
  ValNodeFreeData (head);

  return uop;
}

NLM_EXTERN void AddPubsFromTitle (
  SqnTagPtr stp,
  SeqDescrPtr PNTR desc_list
)

{
  CharPtr         str;
  PubdescPtr      pdp;

  if (stp == NULL || desc_list == NULL) return;
  

  str = SqnTagFindUnused (stp, "PubMed");
  if (str == NULL) 
  {
    str = SqnTagFindUnused (stp, "PMID");
  }
  while (str != NULL)
  {
    pdp = PubdescNew();
    ValNodeAddInt (&(pdp->pub), PUB_PMid, atoi(str));
    SeqDescrAddPointer (desc_list, Seq_descr_pub, (Pointer) pdp);
    str = SqnTagFindUnused (stp, "PubMed");
    if (str == NULL) 
    {
      str = SqnTagFindUnused (stp, "PMID");
    }
  }
}

NLM_EXTERN UserObjectPtr ParseStringIntoStructuredComment (
  UserObjectPtr uop,
  CharPtr str,
  CharPtr prefix,
  CharPtr suffix
)

{
  Char     ch;
  CharPtr  field;
  CharPtr  item;
  CharPtr  last;
  CharPtr  ptr;
  CharPtr  tmp;

  if (uop == NULL) {
    uop = CreateStructuredCommentUserObject (prefix, suffix);
    if (uop == NULL) return uop;
  }
  if (str == NULL) return uop;

  tmp = StringSave (str);
  if (tmp == NULL) return uop;

  last = tmp;
  if (StringDoesHaveText (prefix)) {
    ptr = StringStr (last, prefix);
    if (ptr != NULL) {
      last = ptr + StringLen (prefix);
    }
  }
  if (StringDoesHaveText (suffix)) {
    ptr = StringStr (last, suffix);
    if (ptr != NULL) {
      *ptr = '\0';
    }
  }

  ptr = last;
  ch = *ptr;
  while (ch != '\0') {
    field = last;
    ptr = StringChr (last, '=');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      item = ptr;
      last = StringChr (ptr, ';');
      if (last != NULL) {
        *last = '\0';
        last++;
        ch = *last;
      } else {
        ch = '\0';
      }
      TrimSpacesAroundString (field);
      TrimSpacesAroundString (item);
      AddItemStructuredCommentUserObject (uop, field, item);
    } else {
      ch = '\0';
    }
  }

  MemFree (tmp);

  return uop;
}

/* PHRAP file reading functions */

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

static SeqEntryPtr ReadPhrapDNA (FileCachePtr fcp, CharPtr id)

{
  ByteStorePtr  bs = NULL;
  BioseqPtr     bsp = NULL;
  Char          buf [256];
  Char          ch;
  Boolean       goOn = TRUE;
  Boolean       nonewline;
  CharPtr       p;
  CharPtr       q;
  SeqEntryPtr   sep = NULL;
  CharPtr       str;

  if (fcp == NULL || HasNoText (id)) return NULL;
  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;
  bs = BSNew (1000);
  if (bs == NULL) return NULL;

  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

  bsp->mol = Seq_mol_na;
  bsp->seq_data_type = Seq_code_iupacna;
  bsp->repr = Seq_repr_raw;
  bsp->length = 0;
  bsp->id = MakeSeqID (id);
  SeqMgrAddToBioseqIndex (bsp);

  goOn = TRUE;
  while (goOn) {
    str = FileCacheReadLine (fcp, buf, sizeof (buf), &nonewline);
    if (HasNoText (str)) {
      goOn = FALSE;
    } else {
      p = str;
      q = str;
      ch = *p;
      while (ch != '\0') {
        if (! (IS_ALPHA (ch))) {
          p++;
        } else {
          ch = TO_UPPER (ch);
          if (ch == 'X') {
            ch = 'N';
          }
          *q = ch;
          p++;
          q++;
        }
        ch = *p;
      }
      *q = '\0';
      BSWrite (bs, (VoidPtr) str, (Int4) StringLen (str));
    }
  }

  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);

  BioseqPack (bsp);
  return sep;
}

NLM_EXTERN SeqGraphPtr ReadPhrapQualityFC (FileCachePtr fcp, BioseqPtr bsp)

{
  ByteStorePtr  bs = NULL;
  Char          buf [256];
  Uint1         bytes [128];
  Char          ch;
  Boolean       goOn = TRUE;
  Int2          i;
  Int2          max = INT2_MIN;
  Int2          min = INT2_MAX;
  Boolean       nonewline;
  CharPtr       p;
  Int4          pos;
  CharPtr       q;
  Char          prefix [256];
  size_t        prefixlen;
  SeqGraphPtr   sgp = NULL;
  SeqIntPtr     sintp;
  CharPtr       str;
  int           val;

  if (fcp == NULL || bsp == NULL) return NULL;
  sgp = SeqGraphNew ();
  if (sgp == NULL) return NULL;
  bs = BSNew (1000);
  if (bs == NULL) return NULL;

  goOn = TRUE;
  buf [0] = '\0';
  prefix [0] = '\0';
  while (goOn) {
    StringCpy (buf, prefix);
    prefix [0] = '\0';
    prefixlen = StringLen (buf);
    pos = FileCacheTell (fcp);
    if (NULL == FileCacheReadLine (fcp, buf + prefixlen, sizeof (buf) - prefixlen, &nonewline)) {
      goOn = FALSE;
    } else {
      /* above function returned prefix characters past buf start */
      str = buf;
      if (str [0] == '>') {
        goOn = FALSE;
        if (str [0] == '>') {
          FileCacheSeek (fcp, pos);
        }
      } else {
        i = 0;
        p = str;
        ch = *p;
        while (ch != '\0') {
          while (IS_WHITESP (ch)) {
            p++;
            ch = *p;
          }
          q = p;
          ch = *q;
          while (IS_DIGIT (ch)) {
            q++;
            ch = *q;
          }
          *q = '\0';
          q++;
  
          if (ch == '\0' && nonewline) {
            StringCpy (prefix, p);
          } else {
            if (*p != '\0') {
              if (sscanf (p, "%d", &val) == 1) {
                if (val < 0 || val > 255) {
                  /* error */
                  val = 0;
                }
                bytes [i] = (Uint1) val;
                i++;
                max = MAX (max, (Int2) val);
                min = MIN (min, (Int2) val);
              }
            }
            p = q;
            ch = *p;
          }
        }
        if (i > 0) {
          BSWrite (bs, (Pointer) bytes, (Int4) i);
        }
      }
    }
  }

  sgp->numval = BSLen (bs);
  if (sgp->numval == 0) {
    sgp = SeqGraphFree (sgp);
    return sgp;
  }
  
  BSPutByte (bs, EOF);
  sgp->title = StringSave ("Phrap Quality");
  if (bsp->length != sgp->numval) {
    sgp->flags [0] = 1;
    sgp->compr = (bsp->length) / sgp->numval;
  } else {
    sgp->flags [0] = 0;
    sgp->compr = 1;
  }
  sgp->flags [1] = 0;
  sgp->flags [2] = 3;
  sgp->axis.intvalue = 0;
  sgp->min.intvalue = min;
  sgp->max.intvalue = max;
  sgp->a = 1.0;
  sgp->b = 0;
  sgp->values = (Pointer) bs;

  sintp = SeqIntNew ();
  sintp->from = 0;
  sintp->to = bsp->length - 1;
  sintp->id = SeqIdDup (bsp->id);
  ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);

  return sgp;
}

NLM_EXTERN SeqGraphPtr ReadPhrapQuality (FILE *fp, BioseqPtr bsp)

{
  FileCache    fc;
  Int4         pos;
  SeqGraphPtr  sgp;

  if (fp == NULL || bsp == NULL) return NULL;
  FileCacheSetup (&fc, fp);
  sgp = ReadPhrapQualityFC (&fc, bsp);
  pos = FileCacheTell (&fc);
  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, pos);
  fseek (fp, pos, SEEK_SET);
  return sgp;
}

static Boolean PhrapSequenceHasClipping (FileCachePtr fcp)

{
  Char     buf [256];
  Boolean  goOn = TRUE;
  Boolean  nonewline;
  Boolean  rsult = FALSE;
  CharPtr  str;

  if (fcp == NULL) return FALSE;
  goOn = TRUE;
  while (goOn) {
    str = FileCacheReadLine (fcp, buf, sizeof (buf), &nonewline);
    if (HasNoText (str)) {
      goOn = FALSE;
    } else {
      if (StringNCmp (str, "Clipping", 8) == 0) {
        rsult = TRUE;
      }
    }
  }
  return rsult;
}

static CharPtr BioseqGetLocalIdStr (BioseqPtr bsp)

{
  ObjectIdPtr  oip;
  SeqIdPtr     sip;

  if (bsp == NULL) return NULL;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip != NULL && oip->str != NULL) {
        return oip->str;
      }
    }
  }
  return NULL;
}

static SeqAnnotPtr NewGraphSeqAnnot (CharPtr name, SeqGraphPtr sgp)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (! HasNoText (name)) {
    SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

static CharPtr taglist [] = {
  "", "DNA", "CO", "BaseQuality", "BQ", "Sequence", NULL
};

/* Phrap reading function based on sample code supplied by C. Magness */
NLM_EXTERN SeqEntryPtr ReadPhrapFile (FILE *fp)

{
  BioseqPtr    bsp;
  Char         buf [256];
  FileCache    fc;
  Boolean      goOn = TRUE;
  SeqEntryPtr  head = NULL;
  Int2         i;
  SeqEntryPtr  lastsep;
  SeqGraphPtr  lastsgp;
  Boolean      nonewline;
  CharPtr      p;
  Int4         pos;
  CharPtr      q;
  SeqAnnotPtr  sap;
  SeqEntryPtr  sep = NULL;
  SeqGraphPtr  sgp;
  CharPtr      str;
  Int2         tag;

  if (fp == NULL) return NULL;

  FileCacheSetup (&fc, fp);

  goOn = TRUE;
  while (goOn) {
    str = FileCacheReadLine (&fc, buf, sizeof (buf), &nonewline);
    if (str == NULL) {
      goOn = FALSE;
    } else if (! HasNoText (str)) {
      p = StringChr (str, ' ');
      if (p != NULL) {
        *p = '\0';
        p++;
      }
      tag = 0;
      for (i = 0; taglist [i] != NULL; i++) {
        if (StringCmp (str, taglist [i]) == 0) {
          tag = i;
        }
      }
      if (tag != 0) {
        if (p != NULL) {
          q = StringChr (p, ' ');
          if (q != NULL) {
            *q = '\0';
          }
        }
        switch (tag) {
          case 1 :
          case 2 :
            if (p != NULL) {
              sep = ReadPhrapDNA (&fc, p);
              ValNodeLink (&head, sep);
            }
            /* for new format, sep points to current sequence */
            break;
          case 3 :
            if (p != NULL) {
              sep = head;
              while (sep != NULL && StringCmp (p, BioseqGetLocalIdStr ((BioseqPtr) sep->data.ptrvalue)) != 0) {
                sep = sep->next;
              }
            }
            /* and flow through to case 4 */
          case 4 :
            if (sep != NULL) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              sgp = ReadPhrapQualityFC (&fc, bsp);
              if (sgp != NULL) {
                for (sap = bsp->annot; sap != NULL; sap = sap->next) {
                  if (sap->type == 3) {
                    for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
                      continue;
                    }
                    lastsgp->next = sgp;
                    break;
                  }
                }
                if (sap == NULL) {
                  if (bsp->annot != NULL) {
                    for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
                      continue;
                    }
                    sap->next = NewGraphSeqAnnot ("Graphs", sgp);
                  } else {
                    bsp->annot = NewGraphSeqAnnot ("Graphs", sgp);
                  }
                }
              }
            }
            break;
          case 5 :
            /* unlinkes and removes sep if Clipping line present */
            if (p != NULL) {
              if (PhrapSequenceHasClipping (&fc)) {
                sep = head;
                lastsep = NULL;
                while (sep != NULL && StringCmp (p, BioseqGetLocalIdStr ((BioseqPtr) sep->data.ptrvalue)) != 0) {
                  lastsep = sep;
                  sep = sep->next;
                }
                if (sep != NULL) {
                  if (lastsep != NULL) {
                    lastsep->next = sep->next;
                    sep->next = NULL;
                    SeqEntryFree (sep);
                  } else {
                    head = sep->next;
                    sep->next = NULL;
                    SeqEntryFree (sep);
                  }
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

  pos = FileCacheTell (&fc);
  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, pos);
  fseek (fp, pos, SEEK_SET);

  return head;
}

static ValNodePtr ParseContigOrFeatureTableString (CharPtr contigs, Boolean tabDelimited)

{
  Char        ch;
  Int4        i, j, k;
  CharPtr     str;
  Char        tmp [2048];
  ValNodePtr  vnp;

  vnp = NULL;
  i = 0;
  while (StringLen (contigs + i) > 0) {
    str = contigs + i;
    k = 0;
    ch = str [k];
    while (ch == ' ') {
      k++;
      ch = str [k];
    }
    j = 0;
    if (tabDelimited) {
      while (ch != '\0' && ch != '\t') {
        j++;
        ch = str [j + k];
      }
    } else {
      while (ch != '\0' && ch != ',' && (! (IS_WHITESP (ch)))) {
        j++;
        ch = str [j + k];
      }
    }
    if (ch == '\0') {
      i += j + k;
    } else {
      str [j + k] = '\0';
      i += j + k + 1;
    }
    if (StringLen (str + k) < sizeof (tmp)) {
      StringNCpy_0 (tmp, str + k, sizeof (tmp));
      SqnTrimSpacesAroundString (tmp);
      if (HasNoText (tmp)) {
        ValNodeAdd (&vnp);
      } else {
        ValNodeCopyStr (&vnp, 0, tmp);
      }
    } else {
      ValNodeAddPointer (&vnp, 0, StringSave (str));
    }
  }
  if (vnp != NULL) {
    vnp->choice = (Uint1) ValNodeLen (vnp);
  }
  return vnp;
}

/* ReversePhrap coerces BioseqReverse to work on the SeqGraph byte store */

static void ReversePhrap (SeqGraphPtr sgp, Pointer userdata)

{
  ByteStorePtr  bs;
  Bioseq        bsq;

  if (sgp == NULL || sgp->values == NULL) return;
  if (StringICmp (sgp->title, "Phrap Quality") != 0) return;
  if (sgp->flags [1] != 0 || sgp->flags [2] != 3) return;

  bs = (ByteStorePtr) sgp->values;

  MemSet ((Pointer) &bsq, 0, sizeof (Bioseq));
  bsq.repr = Seq_repr_raw;
  bsq.mol = Seq_mol_na;
  bsq.length = BSLen (bs);
  bsq.seq_data_type = Seq_code_iupacna;
  bsq.seq_data = (SeqDataPtr) bs;

  BioseqReverse (&bsq);
}

NLM_EXTERN SeqEntryPtr SetPhrapContigOrder (SeqEntryPtr head, CharPtr contigs)

{
  BioseqPtr    bsp;
  Char         ch;
  CharPtr      id;
  size_t       len;
  Boolean      minus;
  SeqEntryPtr  sep, lastsep, nextsep, newhead;
  ValNodePtr   vnphead, vnp;

  if (head == NULL || contigs == NULL) return head;
  vnphead = ParseContigOrFeatureTableString (contigs, FALSE);
  if (vnphead == NULL) return head;
  newhead = NULL;
  for (vnp = vnphead; vnp != NULL; vnp = vnp->next) {
    sep = head;
    lastsep = NULL;
    id = (CharPtr) vnp->data.ptrvalue;
    len = StringLen (id);
    minus = FALSE;

    /* look for + or - after accession, indicating orientation */

    if (len > 1) {
      ch = id [len - 1];
      if (ch == '+') {
        id [len - 1] = '\0';
      } else if (ch == '-') {
        id [len - 1] = '\0';
        minus = TRUE;
      }
    }
    while (sep != NULL &&
           StringCmp (id, BioseqGetLocalIdStr ((BioseqPtr) sep->data.ptrvalue)) != 0) {
      lastsep = sep;
      sep = sep->next;
    }
    if (sep != NULL) {
      if (lastsep != NULL) {
        lastsep->next = sep->next;
        sep->next = NULL;
        ValNodeLink (&newhead, sep);
      } else {
        head = sep->next;
        sep->next = NULL;
        ValNodeLink (&newhead, sep);
      }

      /* if - orientation, reverse complement sequence */

      if (minus) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          BioseqRevComp (bsp);

          /* and then reverse phrap scores */

          VisitGraphsOnBsp (bsp, NULL, ReversePhrap);
        }
      }
    }
  }
  for (sep = head; sep != NULL; sep = nextsep) {
    nextsep = sep->next;
    sep->next = NULL;
    SeqEntryFree (sep);
    sep = nextsep;
  }
  ValNodeFreeData (vnphead);
  return newhead;
}

/* More automatic version of ReadAsnFastaOrFlatFile */

NLM_EXTERN CharPtr AsnIoGets (AsnIoPtr aip);
NLM_EXTERN Int2 AsnIoReadBlock (AsnIoPtr aip);

NLM_EXTERN Int4 ReadSequenceAsnFile (
  CharPtr inputFile,
  Boolean binary,
  Boolean compressed,
  Pointer userdata,
  ScanBioseqSetFunc callback
)

{
  AsnIoPtr       aip;
  AsnModulePtr   amp;
  AsnTypePtr     atp, atp_bss, atp_se;
  BytePtr        bp;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Char           buf [128];
  Char           ch;
  Pointer        dataptr = NULL;
  Uint2          datatype;
  Uint2          entityID;
  FILE           *fp = NULL;
  Int4           index = 0;
  Boolean        is_bioseq_set = FALSE;
  size_t         len;
  ObjMgrPtr      omp;
  ObjMgrDataPtr  omdp = NULL;
  ObjMgrTypePtr  omtp = NULL;
  SeqEntryPtr    sep;
  CharPtr        tag;
  CharPtr        tmp;
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
          Message (MSG_FATAL, "Unable to find zcat or gzcat in ScanBioseqSetRelease");
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

  if (binary) {
    AsnIoReadBlock (aip);
  } else {
    AsnIoGets (aip);
  }

  buf [0] = '\0';

  if (aip->bytes > 3 && aip->buf != NULL) {
    bp = aip->buf;
    if (bp [0] > 127 || bp [1] > 127) {
      if (bp [0] == 48 && bp [1] == 128) {
        is_bioseq_set = TRUE;
      }

    } else {

      len = MIN ((size_t) aip->bytes, sizeof (buf));
      StringNCpy_0 (buf, (CharPtr) bp, len);
      if (StringStr (buf, "::=") != NULL) {
        if (StringStr (buf, "Bioseq-set ::=") != NULL) {
          is_bioseq_set = TRUE;
        } else {
          /* first skip past empty space at start of line */

          tag = buf;
          ch = *tag;
          while (ch != '\0' && IS_WHITESP (ch)) {
            tag++;
            ch = *tag;
          }

          /* now find ASN tag */

          tmp = tag;
          ch = *tmp;
          while (ch != '\0' && (! IS_WHITESP (ch))) {
            tmp++;
            ch = *tmp;
          }
          *tmp = '\0';

          omp = ObjMgrReadLock ();
          omtp = ObjMgrTypeFind (omp, 0, tag, NULL);
          ObjMgrUnlock ();
        }
      }
    }
  }


  if (is_bioseq_set) {
    atp = atp_bss;

    while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
      if (atp == atp_se) {
        SeqMgrHoldIndexing (TRUE);
        sep = SeqEntryAsnRead (aip, atp);
        SeqMgrHoldIndexing (FALSE);
        if (sep != NULL) {
          callback (sep, userdata);
          index++;
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
    }
  }

  if (omtp != NULL) {
    aip->scan_for_start = TRUE;
    SeqMgrHoldIndexing (TRUE);
    dataptr = (*(omtp->asnread)) (aip, NULL);
    SeqMgrHoldIndexing (FALSE);

    if (dataptr == NULL) {
      ErrPostEx (SEV_ERROR, 0, 0, "Couldn't read type [%s]", omtp->asnname);
    } else {
      datatype = omtp->datatype;
      if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
        omp = ObjMgrReadLock ();
        omdp = ObjMgrFindByData (omp, dataptr);
        ObjMgrUnlock ();
        if (omdp != NULL && omdp->choice == NULL) {
          /* always want sep above bsp or bssp */
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
        }
      }

      entityID = ObjMgrRegister (datatype, dataptr);
      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL) {
        callback (sep, userdata);
        index++;
      }

      ObjMgrFree (datatype, dataptr);

      omp = ObjMgrGet ();
      ObjMgrReapOne (omp);
      SeqMgrClearBioseqIndex ();
      ObjMgrFreeCache (0);
      FreeSeqIdGiCache ();

      SeqEntrySetScope (NULL);
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

/* ReadAsnFastaOrFlatFile section */

/* GetSeqId skips past LOCUS or ID, or starts past >, skips any white space, then
takes the next token as the seqID.  The return value points to the remaining
copied text, which for feature tables may contain a desired Seq-annot name. */

static CharPtr GetSeqId (CharPtr seqid, CharPtr str, size_t max, Boolean skiptag, Boolean trimwhite)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (seqid != NULL) {
    *seqid = '\0';
  }
  if (str == NULL || seqid == NULL) return FALSE;
  if (skiptag) {
    ch = *str;
    while (ch != '\0' && (! IS_WHITESP (ch))) {
      str++;
      ch = *str;
    }
  }
  ch = *str;
  while (IS_WHITESP (ch)) {
    str++;
    ch = *str;
  }

  StringNCpy_0 (seqid, str, max);
  str = seqid;

  /* find first token, or anything within quotation marks */

  while (ch != '\0' && (! IS_WHITESP (ch))) {
    if (ch == '"') {
      str++;
      ch = *str;
      while (ch != '\0' && ch != '"') {
        str++;
        ch = *str;
      }
      if (ch == '"') {
        str++;
        ch = *str;
      }
    } else {
      str++;
      ch = *str;
    }
  }
  *str = '\0';

  if (ch != 0) {    
    /* trim optional annot name */

    *str = '\0';
    str++;
    ch = *str;
    while (ch != '\0' && (IS_WHITESP (ch))) {
      str++;
      ch = *str;
    }
    if (trimwhite) {
      ptr = str;
      while (ch != '\0' && (! IS_WHITESP (ch))) {
        ptr++;
        ch = *ptr;
      }
      *ptr = '\0';
    }
  }


  /* remove quotation marks in seqid */

  dst = seqid;
  ptr = seqid;
  ch = *ptr;
  while (ch != '\0') {
    if (ch != '"') {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

/* Build contig section */
static void  AddNucToContig (CharPtr accnString, Int4 from, Int4 to,
                             Int4 size, Int2 strand, BioseqPtr segseq,
                             BoolPtr hasgaps, Boolean isgap)

{
  Boolean       allDigits;
  Char          ch;
  DbtagPtr      dp;
  CharPtr       ptr;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  TextSeqIdPtr  tsip;
  long int      val;

  slp = ValNodeNew ((ValNodePtr) segseq->seq_ext);
  if (slp == NULL) return;
  if (segseq->seq_ext == NULL) {
    segseq->seq_ext = (Pointer) slp;
  }

  sintp = SeqIntNew ();
  sintp->from = from;
  sintp->to = to;
  sintp->strand = (Uint1) strand;

  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = (Pointer) sintp;

  if (isgap) {
    sip = ValNodeNew (NULL);
    /* sip = MakeUniqueSeqID ("gap_"); */
    dp = DbtagNew ();
    dp->db = StringSave ("SeqLit");
    dp->tag = ObjectIdNew ();
    dp->tag->id = 0;
    dp->tag->str = NULL;
    sip->choice = SEQID_GENERAL;
    sip->data.ptrvalue = dp;
    if (hasgaps != NULL) {
      *hasgaps = TRUE;
    }
  } else {
    allDigits = TRUE;
    ptr = accnString;
    ch = *ptr;
    while (ch != '\0' && allDigits) {
      if (! IS_DIGIT (ch)) {
        allDigits = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (allDigits && sscanf (accnString, "%ld", &val) == 1) {
      sip = ValNodeNew (NULL);
      sip->choice = (Uint1) SEQID_GI;
      sip->data.intvalue = val;
    } else {
      sip = SeqIdFromAccessionDotVersion (accnString);
      if (sip == NULL) {
        sip = ValNodeNew (NULL);
        tsip = TextSeqIdNew ();
        tsip->accession = StringSave (accnString);
        sip->choice = (Uint1) SEQID_GENBANK;
        sip->data.ptrvalue = tsip;
      }
    }
  }

  sintp->id = sip;

  segseq->length += size;
}

#define accnString field [0]
#define startString field [1]
#define stopString field [2]
#define sizeString field [3]
#define strandString field [4]

static void AdjustContigValues (ValNodePtr line)

{
  Int2        i;
  ValNodePtr  nextAccn;
  ValNodePtr  nextStart;
  long int    num;
  ValNodePtr  thisStop;
  Char        tmp [32];
  Int4        val;

  if (line == NULL) return;
  for (i = 0, thisStop = line->data.ptrvalue; i < 2 && thisStop != NULL; i++, thisStop = thisStop->next) {
    continue;
  }
  line = line->next;
  while (line != NULL && line->data.ptrvalue == NULL) {
    line = line->next;
  }
  if (line == NULL) {
    if (thisStop != NULL) {
      if (sscanf ((CharPtr) thisStop->data.ptrvalue, "%ld", &num) == 1) {
        val = (Int4) num;
        val++;
        sprintf (tmp, "%ld", (long) val);
        thisStop->data.ptrvalue = MemFree (thisStop->data.ptrvalue);
        thisStop->data.ptrvalue = StringSave (tmp);
      }
    }
    return;
  }
  nextAccn = line->data.ptrvalue;
  if (nextAccn != NULL && StringICmp (nextAccn->data.ptrvalue, "gap") == 0) return;
  for (i = 0, nextStart = line->data.ptrvalue; i < 1 && nextStart != NULL; i++, nextStart = nextStart->next) {
    continue;
  }
  if (thisStop != NULL && nextStart != NULL) {
    thisStop->data.ptrvalue = MemFree (thisStop->data.ptrvalue);
    thisStop->data.ptrvalue = StringSave ((CharPtr) nextStart->data.ptrvalue);
  }
}

static void ProcessOneContigLine (ValNodePtr line, BioseqPtr segseq, Int4 lineNum,
                                  BoolPtr hasgaps, Boolean coordsOnMaster)

{
  Boolean     badNumber;
  CharPtr     field [5];
  Int2        i;
  Boolean     isgap;
  long int    num;
  Int4        size;
  Int4        start;
  Int4        stop;
  Int2        strand = Seq_strand_unknown;
  Int4        tmp;
  ValNodePtr  vnp;

  if (line == NULL || segseq == NULL) return;
  vnp = line->data.ptrvalue;
  if (vnp != NULL) {
    for (i = 0; i < 5; i++) {
      field [i] = NULL;
    }
    start = -1;
    stop = -1;
    size = -1;
    for (i = 0, vnp = line->data.ptrvalue; i < 5 && vnp != NULL; i++, vnp = vnp->next) {
      if (field [i] == NULL && (! HasNoText ((CharPtr) vnp->data.ptrvalue))) {
        field [i] = (CharPtr) vnp->data.ptrvalue;
      }
    }
  }

  if (HasNoText (accnString)) return;

  badNumber = FALSE;
  if (sizeString != NULL && sscanf (sizeString, "%ld", &num) == 1) {
    size = num;
  } else {
    size = -1;
  }
  if (startString != NULL && sscanf (startString, "%ld", &num) == 1) {
    start = num;
  } else {
    start = -1;
    badNumber = TRUE;
  }
  if (stopString != NULL && sscanf (stopString, "%ld", &num) == 1) {
    stop = num;
  } else {
    stop = -1;
    badNumber = TRUE;
  }
  if (start < 1 || stop < 1) {
    badNumber = TRUE;
  }
  isgap = FALSE;
  if (StringICmp (accnString, "gap") == 0) {
    if (size >= 0) {
      isgap = TRUE;
      badNumber = FALSE;
      start = 1;
      stop = size;
    }
  }

  if (badNumber) {
    if (startString == NULL) startString = "";
    if (stopString == NULL) stopString = "";
    if (start < 1 && stop < 1) {
      Message (MSG_POST, "Bad number in line %ld - start '%s', stop '%s'",
               (long) lineNum, startString, stopString);
    } else if (start < 1) {
      Message (MSG_POST, "Bad number in line %ld - start '%s'", (long) lineNum, startString);
    } else if (stop < 1) {
      Message (MSG_POST, "Bad number in line %ld - stop '%s'", (long) lineNum, stopString);
    } else {
      Message (MSG_POST, "Bad number in line %ld", (long) lineNum);
    }
    return;
  }

  if (isgap) {
    start = 0;
    stop = size - 1;
  } else {
    if (coordsOnMaster && start == stop) {
      Message (MSG_POST, "Ignoring accession %s", accnString);
      return;
    }

    start--;
    stop--;

    strand = Seq_strand_plus;
    if (strandString != NULL) {
      if (StringStr (strandString, "minus") ||
          StringChr (strandString, '-') ||
          StringStr (strandString, "complement")) {
        strand = Seq_strand_minus;
      }
    }
    if (start > stop) {
      tmp = start;
      start = stop;
      stop = tmp;
      strand = Seq_strand_minus;
    }
    if (strandString != NULL) {
      if (StringStr (strandString, "plus") || StringChr (strandString, '+')) {
        strand = Seq_strand_plus;
      }
    }

    if (coordsOnMaster) {
      stop -= (start + 1);
      start = 0;
    }

    size = ABS (stop - start) + 1;
  }

  AddNucToContig (accnString, start, stop, size, strand, segseq, hasgaps, isgap);
}

static void FreeFeatureTable (ValNodePtr head)

{
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    vnp->data.ptrvalue = ValNodeFreeData (vnp->data.ptrvalue);
  }
  ValNodeFreeData (head);
}

static SeqEntryPtr ReadContigListExEx (FileCachePtr fcp, Boolean coordinatesOnMaster, CharPtr seqid, CharPtr title)

{
  BioseqPtr    bsp;
  DeltaSeqPtr  dsp;
  Boolean      hasgaps;
  ValNodePtr   head = NULL;
  Char         line [1023];
  Int4         lineNum;
  Int4         pos;
  SeqEntryPtr  sep;
  SeqIdPtr     sip = NULL;
  CharPtr      str;
  Char         tmp [128];
  ValNodePtr   vnp;

  if (fcp == NULL) return NULL;

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  if (str != NULL && StringNICmp (line, ">Assembly", 9) == 0) {
    if (seqid == NULL && title == NULL) {
      title = GetSeqId (tmp, line, sizeof (tmp), TRUE, FALSE);
      seqid = tmp;
    }
    str = FileCacheGetString (fcp, line, sizeof (line));
  }
  while (str != NULL) {
    if (! HasNoText (line)) {
      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringNCmp (line, "//", 2) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        break;
      }
      vnp = ParseContigOrFeatureTableString (line, TRUE);
      if (vnp != NULL) {
        ValNodeAddPointer (&head, 0, (Pointer) vnp);
      }
    }
    pos = FileCacheTell (fcp);
    str = FileCacheGetString (fcp, line, sizeof (line));
  }
  if (head == NULL) return NULL;

  bsp = BioseqNew ();
  if (bsp == NULL) {
    FreeFeatureTable (head);
    return NULL;
  }
  bsp->mol = Seq_mol_dna;
  bsp->repr = Seq_repr_seg;
  bsp->seq_ext_type = 1;
  bsp->length = 0;
  if (! StringHasNoText (seqid)) {
    sip = SeqIdFindBest (MakeSeqID (seqid), 0);
  }
  if (sip == NULL) {
    sip = MakeUniqueSeqID ("contig_");
  }
  bsp->id = sip;

  if (! StringHasNoText (title)) {
    str = StringSaveNoNull (title);
    if (str != NULL) {
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
    }
  }

  if (coordinatesOnMaster) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      if (vnp->data.ptrvalue != NULL) {
        AdjustContigValues (vnp);
      }
    }
  }

  lineNum = 0;
  hasgaps = FALSE;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    lineNum++;
    if (vnp->data.ptrvalue != NULL) {
      ProcessOneContigLine (vnp, bsp, lineNum, &hasgaps, coordinatesOnMaster);
    }
  }

  FreeFeatureTable (head);

  if (bsp->seq_ext == NULL) return NULL;

  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

  if (hasgaps) {
    dsp = GappedSeqLocsToDeltaSeqs (bsp->seq_ext);
    if (dsp != NULL) {
      bsp->seq_ext = SeqLocSetFree ((ValNodePtr) bsp->seq_ext);
      bsp->repr = Seq_repr_delta;
      bsp->seq_ext_type = 4;
      bsp->seq_ext = (Pointer) dsp;
    }
  }

  return sep;
}

NLM_EXTERN SeqEntryPtr ReadContigListEx (FILE *fp, Boolean coordinatesOnMaster, CharPtr seqid, CharPtr title)

{
  FileCache    fc;
  Int4         pos;
  SeqEntryPtr  sep;

  if (fp == NULL) return NULL;

  FileCacheSetup (&fc, fp);

  sep = ReadContigListExEx (&fc, coordinatesOnMaster, seqid, title);

  pos = FileCacheTell (&fc);
  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, pos);
  fseek (fp, pos, SEEK_SET);

  return sep;
}

NLM_EXTERN SeqEntryPtr ReadContigList (FILE *fp, Boolean coordinatesOnMaster)

{
  FileCache    fc;
  Int4         pos;
  SeqEntryPtr  sep;

  if (fp == NULL) return NULL;

  FileCacheSetup (&fc, fp);

  sep = ReadContigListExEx (&fc, coordinatesOnMaster, NULL, NULL);

  pos = FileCacheTell (&fc);
  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, pos);
  fseek (fp, pos, SEEK_SET);

  return sep;
}

/* PreCheckSeqForProteinType saves the current file position, then reads lines of
sequence, checking each character for letters that appear only in proteins.  It then
restores the file position, and returns true if it thinks it found a protein. */

static Boolean PreCheckSeqForProteinType (FileCachePtr fcp, Boolean PNTR non_prot_char)

{
  Char     ch;
  Boolean  isProt = FALSE;
  Char     line [1023];
  CharPtr  p;
  Int4     pos;
  CharPtr  str;

  if (fcp == NULL || non_prot_char == NULL) return FALSE;

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "[", 1) == 0 ||
          StringNCmp (line, "]", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringNCmp (line, "//", 2) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        return isProt;
      }

      p = line;
      ch = *p;
      while (ch != '\0') {
        if (! (IS_ALPHA (ch))) {
          p++;
        } else {
          /*
          ch = TO_UPPER (ch);
          if (StringChr ("EFILPQZ", ch) != NULL) {
            isProt = TRUE;
          }
          */
          if (non_prot_char [(int) ch]) {
            isProt = TRUE;
          }
          p++;
        }
        ch = *p;
      }

    }

    str = FileCacheGetString (fcp, line, sizeof (line));
  }

  FileCacheSetup (fcp, fcp->fp);
  FileCacheSeek (fcp, pos);
  fseek (fcp->fp, pos, SEEK_SET);
  return isProt;
}

static Int4 CountBadChars (Int4Ptr bad_chars, Boolean include_warn)
{
  Int4 num_bad_chars = 0, i;
  
  for (i = 0; i < 255; i++)
  {  
      if (bad_chars [i] > 0)
      {
        if (include_warn || (i != '-' && i != '?' && i != ':' && !isdigit(i)))
        {
          num_bad_chars ++;
        }
      }
  }
  return num_bad_chars;
}

static Int4 ReportFlatFileDNABadChars (Int4Ptr bad_chars, CharPtr origline, CharPtr id_str)
{
  Int4    num_bad_chars = 0;
  Int4    i;
  CharPtr msg = NULL;
  CharPtr msg_format = "%d different illegal characters were found:";
  CharPtr msg_format_single = "One illegal character (%c) was found %d times.";
  CharPtr msg_format_one = "'%c' (%d),";
  CharPtr msg_ptr;
  Boolean   indexerVersion;
  
  /* don't warn indexers about numbers */
  indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
  if (indexerVersion)
  {
    for (i = 0; i < 255; i++)
    {
      if (isdigit (i))
      {
        /* don't report numbers as bad characters */
        bad_chars [i] = 0;
      }
    }
  }
  
  num_bad_chars = CountBadChars(bad_chars, TRUE);
  
  if (num_bad_chars == 0) return 0;
  
  if (num_bad_chars == 1)
  {
      msg = (CharPtr) MemNew ((StringLen (msg_format_single) + 15) * sizeof (Char));
      if (msg != NULL)
      {
        for (i = 0; i < 255; i++)
        {
            if (bad_chars [i] > 0)
            {
            sprintf (msg, msg_format_single, i, bad_chars [i]);
            }
        }
      }
  }
  else
  {
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (msg_format) + 15
                   + num_bad_chars * (StringLen (msg_format_one) + 18)));
      if (msg != NULL)
      {
        sprintf (msg, msg_format, num_bad_chars);
        msg_ptr = msg + StringLen (msg);
        for (i = 0; i < 255; i ++)
        {
            if (bad_chars [i] > 0)
            {
              sprintf (msg_ptr, msg_format_one, i, bad_chars [i]);
              msg_ptr += StringLen (msg_ptr);
            }
        }
        msg_ptr --;
        *msg_ptr = 0;
      }
  }
  StringCpy (origline + 60, "...");
  origline [63] = '\0';
  if (indexerVersion) {
    Message (MSG_POSTERR, "%s", msg);
    Message (MSG_POSTERR, "Offending line started at: %s", origline);
    if (!StringHasNoText (id_str))
    {
      Message (MSG_POSTERR, "Sequence %s:", id_str);
    }
  } else {
    Message (MSG_POSTERR, "Sequence %s: %s. Offending line started at: %s.",
                        StringHasNoText (id_str) ? "no id provided" : id_str,
                        msg,
                        origline);
  }

  
  num_bad_chars = CountBadChars (bad_chars, FALSE);
  
  return num_bad_chars;  
}

static Boolean IsNonSeqChar (Char ch, Boolean is_prot)
{
  if (! isalpha (ch))
  {
    return TRUE;
  }
  
  ch = TO_LOWER (ch);
  if (StringChr ("atgcbdhkmnrsuvwy", ch) != NULL)
  {
    return FALSE;
  }
  else if (is_prot && StringChr ("efilpqxz", ch) != NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

/* ReadFlatFileDNA reads lines of sequence into a byte store.  Unless it is forced to be
treated as a nucleotide or a protein, it first calls PreCheckSeqForProteinType to look at
the sequence in advance, checking for protein-specific letters. If it encounters a non-
printing character, it completes the read but returns NULL. */

static ByteStorePtr ReadFlatFileDNA (FileCachePtr fcp, BoolPtr protPtr, Boolean forceNuc,
                                     Boolean forceProt, Boolean fastaAsSimpleSeq,
                                     Boolean strictCheck, BoolPtr perr,
                                     CharPtr id_str)

{
  Char           ch;
  ByteStorePtr   bs = NULL;
  Boolean        isProt = FALSE;
  Char           line [1023];
  Boolean        noErrors = TRUE;
  Char           origline [1023];
  CharPtr        p;
  CharPtr        q;
  Int4           pos;
  CharPtr        str;
  Int4           bad_char [256];
  Boolean        non_prot_char [256];
  Int4           num_bad = 0;
  Boolean        is_nuc_char [256];
  Boolean        is_prot_char [256];
  CharPtr        nuc_list = "atgcbdhkmnrsuvwy";
  CharPtr        prot_list = "abcdefghijklmnopqrstuvwxyz";
  CharPtr        ptr;

  if (fcp == NULL) return NULL;
  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  
  if (perr != NULL)
  {
    *perr = FALSE;
  }

  MemSet (is_nuc_char, 0, sizeof (is_nuc_char));

  ptr = nuc_list;
  ch = *ptr;
  while (ch != '\0') {
    is_nuc_char [(int) ch] = TRUE;
    ch = TO_UPPER (ch);
    is_nuc_char [(int) ch] = TRUE;
    ptr++;
    ch = *ptr;
  }

  MemSet (is_prot_char, 0, sizeof (is_prot_char));

  ptr = prot_list;
  ch = *ptr;
  while (ch != '\0') {
    is_prot_char [(int) ch] = TRUE;
    ch = TO_UPPER (ch);
    is_prot_char [(int) ch] = TRUE;
    ptr++;
    ch = *ptr;
  }

  if (forceNuc) {
    isProt = FALSE;
  } else if (forceProt) {
    isProt = TRUE;
  } else if (protPtr != NULL) {
    MemSet (non_prot_char, 0, sizeof (non_prot_char));
    non_prot_char [(int) 'E'] = TRUE;
    non_prot_char [(int) 'F'] = TRUE;
    non_prot_char [(int) 'I'] = TRUE;
    non_prot_char [(int) 'L'] = TRUE;
    non_prot_char [(int) 'P'] = TRUE;
    non_prot_char [(int) 'Q'] = TRUE;
    non_prot_char [(int) 'Z'] = TRUE;
    non_prot_char [(int) 'e'] = TRUE;
    non_prot_char [(int) 'f'] = TRUE;
    non_prot_char [(int) 'i'] = TRUE;
    non_prot_char [(int) 'l'] = TRUE;
    non_prot_char [(int) 'p'] = TRUE;
    non_prot_char [(int) 'q'] = TRUE;
    non_prot_char [(int) 'z'] = TRUE;
    isProt = PreCheckSeqForProteinType (fcp, non_prot_char);
  }
  if (protPtr != NULL) {
    *protPtr = isProt;
  }

  MemSet (bad_char, 0, sizeof (bad_char));

  origline [0] = '\0';

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "[", 1) == 0 ||
          StringNCmp (line, "]", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        num_bad = ReportFlatFileDNABadChars (bad_char, origline, id_str);
        if (perr != NULL && num_bad > 0)
        {
          *perr = TRUE;
        }
        return bs;
      } else if (StringNCmp (line, "//", 2) == 0) {
        num_bad = ReportFlatFileDNABadChars (bad_char, origline, id_str);
        if (perr != NULL && num_bad > 0)
        {
          *perr = TRUE;
        }
        return bs;
      }

      if (noErrors) {
        StringNCpy_0 (origline, line, sizeof (origline));
      }
      p = line;
      q = line;
      ch = *p;
      while (ch != '\0') {
        ch = TO_UPPER (ch);
        if (IS_WHITESP (ch)) {
        } else if (! (IS_ALPHA (ch))) {
          if (isProt && (ch == '*' || ch == '-')) {
            *q = ch;
            q++;
          } else if (! IS_PRINT (ch)) {
            bs = BSFree (bs);
          }
          else if (IS_DIGIT (ch))
          {
              /* only report for strictCheck */
              if (strictCheck)
              {
                bad_char [(int) ch] ++;
                noErrors = FALSE;
              }
          }
/*        Do not convert question marks to Ns  
          else if (ch == '?')
          {
              bad_char [(int) ch] ++;
              *q = 'N';
              q++;
              noErrors = FALSE;
          } */
          else
          {
              bad_char [(int) ch] ++;
              noErrors = FALSE;
          }
        } else {
          /* if (IsNonSeqChar (ch, isProt)) */
          if ((isProt && (! is_prot_char [(int) ch])) || ((! isProt) && (! is_nuc_char [(int) ch]))) {
            bad_char [(int) ch] ++;
            noErrors = FALSE;
          }
          else
          {
            if (! fastaAsSimpleSeq) {
              ch = TO_UPPER (ch);
            }
            if (! isProt) {
              if (ch == 'U') {
                ch = 'T';
              } else if (ch == 'u') {
                ch = 't';
              } else if (ch == 'X') {
                ch = 'N';
              } else if (ch == 'x') {
                ch = 'n';
              }
            }
            *q = ch;
            q++;
          }
        }
        p++;
        ch = *p;
      }
      *q = '\0';
      if (bs != NULL) {
        BSWrite (bs, (VoidPtr) line, (Int4) StringLen (line));
      }

    }

    pos = FileCacheTell (fcp);
    str = FileCacheGetString (fcp, line, sizeof (line));
  }

  if (bs != NULL && BSLen (bs) < 1) {
    bs = BSFree (bs);
  }
  num_bad = ReportFlatFileDNABadChars (bad_char, origline, id_str);
  if (perr != NULL && num_bad > 0)
  {
    *perr = TRUE;
  }

  return bs;
}

static SimpleSeqPtr ByteStoreToSimpleSeq (ByteStorePtr bs, CharPtr seqid, CharPtr title)

{
  SimpleSeqPtr  ssp;

  if (bs == NULL) return NULL;
  ssp = SimpleSeqNew ();
  if (ssp == NULL) return NULL;

  ssp->seq = bs;
  ssp->seqlen = BSLen (bs);

  if (! HasNoText (seqid)) {
    ssp->id [0] = StringSave (seqid);
    ssp->numid = 1;
  }
  if (! HasNoText (title)) {
    ssp->title = StringSave (title);
  }

  return ssp;
}

/* ReadFeatureTable reads lines of feature intervals and qualifiers into a Seq-annot. */

#define NUM_FTABLE_COLUMNS  6

#define START_TAG           0
#define STOP_TAG            1
#define FEAT_TYPE_TAG       2
#define QUAL_TYPE_TAG       3
#define QUAL_VAL_TAG        4
#define STRAND_TAG          5

#define startStr   field [START_TAG]
#define stopStr    field [STOP_TAG]
#define featType   field [FEAT_TYPE_TAG]
#define qualType   field [QUAL_TYPE_TAG]
#define qualVal    field [QUAL_VAL_TAG]
#define strandStr  field [STRAND_TAG]


static Char UnexpectedCharInPositionString (CharPtr str)
{
  CharPtr cp;

  if (str == NULL) {
    return 0;
  }

  cp = str;
  while (*cp == '<' || *cp == '>' || *cp == '^' || isdigit (*cp) || *cp == '-') {
    cp++;
  }
  return *cp;
}


static Boolean ParseFeatTableLine (CharPtr line, Int4Ptr startP, Int4Ptr stopP,
                                   BoolPtr partial5P, BoolPtr partial3P, BoolPtr ispointP,
                                   BoolPtr isminusP, CharPtr PNTR featP, CharPtr PNTR qualP,
                                   CharPtr PNTR valP, Int4 offset, Int4 lin_num)

{
  Boolean     badNumber;
  CharPtr     field [NUM_FTABLE_COLUMNS];
  Int2        i;
  Boolean     isminus = FALSE;
  Boolean     ispoint = FALSE;
  size_t      len;
  ValNodePtr  parsed;
  Boolean     partial5 = FALSE;
  Boolean     partial3 = FALSE;
  Int4        start;
  Int4        stop;
  CharPtr     str;
  Int4        tmp;
  long int    val;
  ValNodePtr  vnp;
  Char        badch;

  if (line == NULL || HasNoText (line)) return FALSE;
  if (*line == '[') return FALSE; /* offset and other instructions encoded in brackets */
  parsed = ParseContigOrFeatureTableString (line, TRUE);
  if (parsed == NULL) return FALSE;

  for (i = 0; i < NUM_FTABLE_COLUMNS; i++) {
    field [i] = NULL;
  }
  start = -1;
  stop = -1;
  vnp = parsed;
  for (i = 0; i < NUM_FTABLE_COLUMNS && vnp != NULL; i++) {
    if (field [i] == NULL) {
      if (! HasNoText ((CharPtr) vnp->data.ptrvalue)) {
        field [i] = (CharPtr) vnp->data.ptrvalue;
      }
    }
    vnp = vnp->next;
  }

  badNumber = FALSE;
  str = startStr;
  badch = UnexpectedCharInPositionString (str);
  if (badch != 0) {
    Message (MSG_POSTERR, "Unexpected characters in from column of line %d - first bad character is '%c'", lin_num, badch);
  }
  if (str != NULL && *str == '<') {
    partial5 = TRUE;
    str++;
  }
  len = StringLen (str);
  if (len > 1 && str [len - 1] == '^') {
    ispoint = TRUE;
    str [len - 1] = '\0';
  }
  if (str != NULL && sscanf (str, "%ld", &val) == 1) {
    start = val;
  } else {
    start = -1;
    badNumber = TRUE;
  }
  str = stopStr;
  badch = UnexpectedCharInPositionString (str);
  if (badch != 0) {
    Message (MSG_POSTERR, "Unexpected characters in to column of line %d - first bad character is '%c'", lin_num, badch);
  }
  if (str != NULL && *str == '>') {
    partial3 = TRUE;
    str++;
  }
  if (str != NULL && sscanf (str, "%ld", &val) == 1) {
    stop = val;
  } else {
    stop = -1;
    badNumber = TRUE;
  }

  if (badNumber) {
    start = -1;
    stop = -1;
  } else {
    start--;
    stop--;
    if (strandStr != NULL) {
      if (StringStr (strandStr, "minus") ||
          StringChr (strandStr, '-') ||
          StringStr (strandStr, "complement")) {
        if (start < stop) {
          tmp = start;
          start = stop;
          stop = tmp;
        }
        isminus = TRUE;
      }
    }
  }

  *startP = start + offset;
  *stopP = stop + offset;
  *partial5P = partial5;
  *partial3P = partial3;
  *ispointP = ispoint;
  *isminusP = isminus;
  *featP = StringSaveNoNull (featType);
  *qualP = StringSaveNoNull (qualType);
  *valP = StringSaveNoNull (qualVal);

  ValNodeFreeData (parsed);
  return TRUE;
}

static CharPtr aaList [] = {
  "-", "Gap", "Gap",        /* cannot be recognized because we split tRNA-xxx */
  "A", "Ala", "Alanine",
  "B", "Asx", "Asp or Asn",
  "C", "Cys", "Cysteine",
  "D", "Asp", "Aspartic Acid",
  "E", "Glu", "Glutamic Acid",
  "F", "Phe", "Phenylalanine",
  "G", "Gly", "Glycine",
  "H", "His", "Histidine",
  "I", "Ile", "Isoleucine",
  "K", "Lys", "Lysine",
  "L", "Leu", "Leucine",
  "M", "Met", "Methionine",
  "N", "Asn", "Asparagine",
  "P", "Pro", "Proline",
  "Q", "Gln", "Glutamine",
  "R", "Arg", "Arginine",
  "S", "Ser", "Serine",
  "T", "Thr", "Threonine",
  "V", "Val", "Valine",
  "W", "Trp", "Tryptophan",
  "X", "Xxx", "Undetermined or atypical",
  "Y", "Tyr", "Tyrosine",
  "Z", "Glx", "Glu or Gln",
  "U", "Sec", "Selenocysteine",
  "*", "Ter", "Termination",
  "O", "Pyl", "Pyrrolysine",
  "J", "Xle", "Leu or Ile",
  NULL, NULL, NULL
};

NLM_EXTERN Uint1 FindTrnaAA3 (CharPtr str)

{
  Uint1    aa;
  Int2     i;
  CharPtr  ptr;
  Char     tmp [128];

  if (HasNoText (str)) return 0;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  SqnTrimSpacesAroundString (tmp);
  for (i = 0; aaList [i] != NULL; i += 3) {
    if (StringICmp (aaList [i + 1], tmp) == 0) {
      ptr = aaList [i];
      aa = (Uint1) ptr [0];
      return aa;
    }
  }
  if (StringICmp ("fMet", tmp) == 0) return (Uint1) 'M';
  if (StringICmp ("OTHER", tmp) == 0) return (Uint1) 'X';
  if (StringICmp ("Aspartate", tmp) == 0) return (Uint1) 'D';
  if (StringICmp ("Glutamate", tmp) == 0) return (Uint1) 'E';
  return 0;
}

NLM_EXTERN Uint1 FindTrnaAA (CharPtr str)

{
  Uint1    aa;
  Int2     i;
  Int2     j;
  CharPtr  ptr;
  Char     tmp [128];

  if (HasNoText (str)) return 0;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  SqnTrimSpacesAroundString (tmp);
  for (i = 0; aaList [i] != NULL; i += 3) {
    for (j = 0; j < 3; j++) {
      if (StringICmp (aaList [i + j], tmp) == 0) {
        ptr = aaList [i];
        aa = (Uint1) ptr [0];
        return aa;
      }
    }
  }
  if (StringICmp ("fMet", tmp) == 0) return (Uint1) 'M';
  if (StringICmp ("OTHER", tmp) == 0) return (Uint1) 'X';
  if (StringICmp ("Aspartate", tmp) == 0) return (Uint1) 'D';
  if (StringICmp ("Glutamate", tmp) == 0) return (Uint1) 'E';
  return 0;
}

NLM_EXTERN CharPtr FindTrnaAAIndex (CharPtr str)

{
  Int2  i;
  Int2  j;
  Char  tmp [128];

  if (StringHasNoText (str)) return 0;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  TrimSpacesAroundString (tmp);
  for (i = 0; aaList [i] != NULL; i += 3) {
    for (j = 0; j < 3; j++) {
      if (StringICmp (aaList [i + j], tmp) == 0) {
        return aaList [i + 1];
      }
    }
  }
  if (StringICmp ("fMet", tmp) == 0) return "Methionine";
  if (StringICmp ("OTHER", tmp) == 0) return "Selenocysteine";
  if (StringICmp ("Aspartate", tmp) == 0) return "Aspartic Acid";
  if (StringICmp ("Glutamate", tmp) == 0) return "Glutamic Acid";
  return NULL;
}

NLM_EXTERN Char FindResidueByName (CharPtr res_name, SeqCodeTablePtr sctp)
{
  Int2    i;
  Uint1   last;
  Int2    res = INVALID_RESIDUE;
  Int4    len;
  
  if (res_name == NULL || sctp == NULL) return INVALID_RESIDUE;
  last = LastResidueInCode (sctp);
  len = StringLen (res_name);
  if (len == 1)
  {
    res = GetResidueForSymbol (sctp, res_name[0]);
  }
  else
  {
    for (i = 0; aaList [3 * i] != NULL && res == INVALID_RESIDUE; i++)
    {
      if (StringICmp (res_name, aaList [3 * i + 1]) == 0
          || StringICmp (res_name, aaList [3 * i + 2]) == 0)
      {
        res = GetResidueForSymbol (sctp, aaList [3 * i][0]);
      }
    }
  }
  return (Char) res;
}


NLM_EXTERN Uint1 ParseTRnaString (CharPtr strx, BoolPtr justTrnaText, Uint1Ptr cdP, Boolean noSingleLetter)

{
  Uint1       aa;
  Char        ch;
  Char        codon [16];
  Uint1       curraa;
  ValNodePtr  head;
  Int2        i;
  Boolean     justt = TRUE;
  CharPtr     str;
  tRNA        tr;
  ValNodePtr  vnp;

  if (justTrnaText != NULL) {
    *justTrnaText = FALSE;
  }
  if (cdP != NULL) {
    for (i = 0; i < 6; i++) { 
      cdP [i] = 255;
    }
  }
  if (StringHasNoText (strx)) return 0;

  MemSet ((Pointer) &tr, 0, sizeof (tRNA));
  for (i = 0; i < 6; i++) {
    tr.codon [i] = 255;
  }

  aa = 0;
  head = TokenizeTRnaString (strx);
  for (vnp = head; (aa == 0 || aa == 'A') && vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    curraa = FindTrnaAA (str);
    if (noSingleLetter && StringLen (str) == 1) {
      curraa = 0;
    }
    if (curraa != 0) {
      if (aa == 0 || aa == 'A') {
        aa = curraa;
      }
    } else if (StringICmp ("tRNA", str) != 0 &&
               StringICmp ("transfer", str) != 0 &&
               StringICmp ("RNA", str) != 0 &&
               StringICmp ("product", str) != 0) {
      if (cdP != NULL && StringLen (str) == 3) {
        StringCpy (codon, str);
        for (i = 0; i < 3; i++) {
          if (codon [i] == 'U') {
            codon [i] = 'T';
          }
        }
        if (ParseDegenerateCodon (&tr, (Uint1Ptr) codon)) {
          for (i = 0; i < 6; i++) { 
            cdP [i] = tr.codon [i];
          }
        } else {
          justt = FALSE;
        }
      } else {
        justt = FALSE;
      }
    }
  }
  for (; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    curraa = FindTrnaAA (str);
    if (curraa != 0) {
    } else if (StringICmp ("tRNA", str) != 0 &&
               StringICmp ("transfer", str) != 0 &&
               StringICmp ("RNA", str) != 0 &&
               StringICmp ("product", str) != 0) {
      if (cdP != NULL && StringLen (str) == 3) {
        StringCpy (codon, str);
        for (i = 0; i < 3; i++) {
          if (codon [i] == 'U') {
            codon [i] = 'T';
          }
        }
        if (ParseDegenerateCodon (&tr, (Uint1Ptr) codon)) {
          for (i = 0; i < 6; i++) { 
            cdP [i] = tr.codon [i];
          }
        } else {
          justt = FALSE;
        }
      } else {
        justt = FALSE;
      }
    }
  }
  ValNodeFreeData (head);
  if (justt) {
    str = strx;
    ch = *str;
    while (ch != '\0') {
      if (IS_DIGIT (ch)) {
        justt = FALSE;
      }
      str++;
      ch = *str;
    }
  }
  if (justTrnaText != NULL) {
    *justTrnaText = justt;
  }
  return aa;
}

static Boolean ThreeLettersPlusDigits (CharPtr str)

{
  Char    ch;
  Uint2    i;
  size_t  len;

  if (StringHasNoText (str)) return FALSE;
  len = StringLen (str);
  if (len < 4) return FALSE;
  for (i = 0; i < 3; i++) {
    ch = str [i];
    if (! IS_ALPHA (ch)) return FALSE;
  }
  for (i = 3; i < len; i++) {
    ch = str [i];
    if (! IS_DIGIT (ch)) return FALSE;
  }
  return TRUE;
}

NLM_EXTERN ValNodePtr TokenizeTRnaString (CharPtr strx)

{
  Char        ch;
  ValNodePtr  head;
  Int2        i, j, k;
  size_t      len;
  CharPtr     ptr;
  Char        str [256];
  CharPtr     strs;
  Char        tmp [128];

  if (HasNoText (strx)) return NULL;
  strs = StringSave (strx);
  head = NULL;
  /* SGD Tx(NNN)c or Tx(NNN)c#, where x is the amino acid, c is the chromosome (A-P, Q for mito),
     and optional # is presumably for individual tRNAs with different anticodons and the same
     amino acid */
  len = StringLen (strs);
  if (len >= 8 && len <= 10) {
    if (strs [0] == 'T' || strs [0] == 't') {
      if (IS_ALPHA (strs [1]) && strs [2] == '('
          && strs [6] == ')' && IS_ALPHA (strs [7])) {
        if (len == 8 ||
            (len == 9 && IS_DIGIT (strs [8])) ||
            (len == 10 && IS_DIGIT (strs [8]) && IS_DIGIT (strs [9]))) {
          tmp [0] = '('; /* parse SGD tRNA anticodon */
          tmp [1] = strs [5]; /* reverse */
          tmp [2] = strs [4];
          tmp [3] = strs [3];
          tmp [4] = ')';
          tmp [5] = '\0';
          for (i = 1; i < 4; i++) {
            ch = tmp [i];
            ch = TO_UPPER (ch);
            switch (ch) {
              case 'A' :
                ch = 'T';
                break;
              case 'C' :
                ch = 'G';
                break;
              case 'G' :
                ch = 'C';
                break;
              case 'T' :
                ch = 'A';
                break;
              case 'U' :
                ch = 'A';
                break;
              default :
                ch = '?';
                break;
            }
            tmp [i] = ch; /* and complement */
          }
          ValNodeCopyStr (&head, 0, tmp);
          tmp [0] = strs [1]; /* parse SGD tRNA amino acid */
          tmp [1] = '\0';
          ValNodeCopyStr (&head, 0, tmp);
          MemFree (strs);
          return head;
        }
      }
    }
  }
  ptr = strs;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '*') {  /* keep possible terminator tRNA symbol */
    } else if (IS_WHITESP (ch) ||
               ch == '-' || ch == ',' || ch == ';' ||
               ch == ':' || ch == '(' || ch == ')' ||
               ch == '=' || ch == '\'' || ch == '_' ||
               ch == '~') {
     *ptr = ' ';
    }
    ptr++;
    ch = *ptr;
  }
  i = 0;
  while (StringLen (strs + i) > 0) {
    StringNCpy_0 (str, strs + i, sizeof (str));
    k = 0;
    ch = str [k];
    while (ch == ' ') {
      k++;
      ch = str [k];
    }
    j = 0;
    while (ch != '\0' && ch != ' ') {
      j++;
      ch = str [j + k];
    }
    if (ch == ' ') {
      str [j + k] = '\0';
      i += j + k + 1;
    } else {
      i += j + k;
    }
    StringNCpy_0 (tmp, str + k, sizeof (tmp));
    if (StringNICmp (tmp, "tRNA", 4) == 0) {
      tmp [0] = ' ';
      tmp [1] = ' ';
      tmp [2] = ' ';
      tmp [3] = ' ';
    }
    SqnTrimSpacesAroundString (tmp);
    if (! HasNoText (tmp)) {
      if (ThreeLettersPlusDigits (tmp)) {
        tmp [3] = '\0';
      }
      ValNodeCopyStr (&head, 0, tmp);
    }
  }
  MemFree (strs);
  return head;
}

static CharPtr bondList [] = {
  "", "disulfide", "thiolester", "xlink", "thioether", NULL
};

static CharPtr siteList [] = {
  "", "active", "binding", "cleavage", "inhibit", "modified",
  "glycosylation", "myristoylation", "mutagenized", "metal binding",
  "phosphorylation", "acetylation", "amidation", "methylation",
  "hydroxylation", "sulfatation", "oxidative deamination",
  "pyrrolidone carboxylic acid", "gamma carboxyglutamic acid",
  "blocked", "lipid binding", "np binding", "DNA binding",
  "signal peptide", "transit peptide", "transmembrane region",
  "nitrosylation", NULL
};

static void StripHyphens (CharPtr str)

{
  Char     ch;
  CharPtr  ptr;

  if (str == NULL) return;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '-') {
      *ptr = ' ';
    }
    ptr++;
    ch = *ptr;
  }
}

/*
static Uint1 ParseCodon (CharPtr str)

{
  Char   ch;
  Uint1  codon [4];
  Int2   i, j;

  if (str == NULL) return 255;
  for (i = 0, j = 1; i < 3; i++, j++) {
    ch = TO_UPPER (str [j]);
    codon [i] = (Uint1) ch;
  }
  codon [3] = 0;
  return IndexForCodon (codon, Seq_code_iupacna);
}
*/

extern Boolean ParseAnticodon (SeqFeatPtr sfp, CharPtr val, Int4 offset);

static CharPtr orgRefList [] = {
  "", "organism", "mitochondrion", "div", "lineage", "gcode", "mgcode", "pgcode", NULL
};

static OrgNamePtr GetOrMakeOnp (OrgRefPtr orp)

{
  OrgNamePtr  onp;

  if (orp == NULL) return NULL;
  if (orp->orgname != NULL) return orp->orgname;
  onp = OrgNameNew ();
  orp->orgname = onp;
  return onp;
}

static Boolean ParseQualIntoBioSource (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  BioSourcePtr  biop;
  Int2          found, j;
  int           num;
  OrgModPtr     omp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SubSourcePtr  ssp;
  CharPtr       str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return FALSE;
  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;

  found = 0;
  for (j = 0, str = orgRefList [j]; str != NULL; j++, str = orgRefList [j]) {
    if (StringsAreEquivalent (qual, str)) {
      found = j;
    }
  }
  if (found > 0) {
    switch (found) {
      case 1 :
        orp->taxname = MemFree (orp->taxname);
        orp->taxname = StringSave (val);
        break;
      case 2 :
        biop->genome = GENOME_mitochondrion;
        break;
      case 3 :
        onp = GetOrMakeOnp (orp);
        if (onp == NULL) return FALSE;
        onp->div = MemFree (onp->div);
        onp->div = StringSave (val);
        break;
      case 4 :
        onp = GetOrMakeOnp (orp);
        if (onp == NULL) return FALSE;
        onp->lineage = MemFree (onp->lineage);
        onp->lineage = StringSave (val);
        break;
      case 5 :
        onp = GetOrMakeOnp (orp);
        if (onp == NULL) return FALSE;
        if (sscanf (val, "%d", &num) == 1) {
          onp->gcode = (Uint1) num;
        }
        break;
      case 6 :
        onp = GetOrMakeOnp (orp);
        if (onp == NULL) return FALSE;
        if (sscanf (val, "%d", &num) == 1) {
          onp->mgcode = (Uint1) num;
        }
        break;
      case 7 :
        onp = GetOrMakeOnp (orp);
        if (onp == NULL) return FALSE;
        if (sscanf (val, "%d", &num) == 1) {
          onp->pgcode = (Uint1) num;
        }
        break;
      default :
        break;
    }
    return TRUE;
  }

  found = EquivalentOrgMod (qual);
  if (found > 0) {
    if (found == 32) {
      found = 253;
    } else if (found == 33) {
      found = 254;
    } else if (found == 34) {
      found = 255;
    }
    onp = GetOrMakeOnp (orp);
    if (onp == NULL) return FALSE;
    omp = OrgModNew ();
    if (omp == NULL) return FALSE;
    omp->subtype = (Uint1) found;
    omp->subname = StringSave (val);
    omp->next = onp->mod;
    onp->mod = omp;
    return TRUE;
  }

  found = EquivalentSubSource (qual);

  if (found > 0) {
    ssp = SubSourceNew ();
    if (ssp == NULL) return FALSE;
    ssp->subtype = (Uint1) found;
    ssp->name = StringSave (val);
    ssp->next = biop->subtype;
    biop->subtype = ssp;
    return TRUE;
  }

  return FALSE;
}

static Boolean IS_real (CharPtr str)

{
  Char     ch;
  Boolean  nodigits = TRUE;
  Boolean  isinteger = TRUE;

  if (StringHasNoText (str)) return FALSE;
  ch = *str;
  while (ch != '\0') {
    if (ch == '+' || ch == '-' || ch == '.' || ch == 'E' || ch == 'e') {
      isinteger = FALSE;
    } else if (ch < '0' || ch > '9') {
      return FALSE;
    } else {
      nodigits = FALSE;
    }
    str++;
    ch = *str;
  }
  if (nodigits) return FALSE;
  if (isinteger) return FALSE;
  return TRUE;
}

static void AddFieldToSnpStsCloneUserObject (UserObjectPtr uop, CharPtr qual, CharPtr val)

{
  UserFieldPtr  curr;
  UserFieldPtr  prev = NULL;
  long int      num;
  ObjectIdPtr   oip;
  double        dbl;

  if (uop == NULL || StringHasNoText (qual) || StringHasNoText (val)) return;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, qual) == 0) {
      break;
    }
    prev = curr;
  }

  if (curr == NULL) {
    curr = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->str = StringSave (qual);
    curr->label = oip;
    if (IS_real (val) && sscanf (val, "%lf", &dbl) == 1) {
      curr->choice = 3; /* real */
      curr->data.realvalue = (FloatHi) dbl;
    } else if (sscanf (val, "%ld", &num) == 1) {
      curr->choice = 2; /* integer */
      curr->data.intvalue = (Int4) num;
    } else {
      curr->choice = 1; /* visible string */
      curr->data.ptrvalue = StringSave (val);
    }

    /* link at end of list */

    if (prev != NULL) {
      prev->next = curr;
    } else {
      uop->data = curr;
    }
  }
}

static CharPtr snpQualList [] = {
  "", "snp_class", "weight", "chrcnt", "ctgcnt", "loccnt", "snp_het", "snp_het_se",
  "snp_maxrate", "snp_gtype", "snp_linkout", "snp_valid", NULL
};

static UserObjectPtr CreateSnpUserObject (void)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  uop = UserObjectNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("dbSnpSynonymyData");
  uop->type = oip;

  return uop;
}

static UserObjectPtr GetSnpUserObject (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (sfp == NULL) return NULL;
  if (sfp->ext == NULL) {
    sfp->ext = CreateSnpUserObject ();
  }
  uop = sfp->ext;
  if (uop == NULL) return NULL;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "dbSnpSynonymyData") != 0) return NULL;
  return uop;
}

static Boolean ParseQualIntoSnpUserObject (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  Int2           found, j;
  CharPtr        str;
  UserObjectPtr  uop;

  found = 0;
  for (j = 0, str = snpQualList [j]; str != NULL; j++, str = snpQualList [j]) {
    if (StringICmp (qual, str) == 0) {
      found = j;
    }
  }

  if (found > 0) {
    uop = GetSnpUserObject (sfp);
    if (uop == NULL) return FALSE;
    AddFieldToSnpStsCloneUserObject (uop, qual, val);
    return TRUE;
  }

  return FALSE;
}

static CharPtr stsQualList [] = {
  "", "sts_dsegs", "sts_aliases", "weight", NULL
};

static UserObjectPtr CreateStsUserObject (void)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  uop = UserObjectNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("stsUserObject");
  uop->type = oip;

  return uop;
}

static UserObjectPtr GetStsUserObject (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (sfp == NULL) return NULL;
  if (sfp->ext == NULL) {
    sfp->ext = CreateStsUserObject ();
  }
  uop = sfp->ext;
  if (uop == NULL) return NULL;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "stsUserObject") != 0) return NULL;
  return uop;
}

static Boolean ParseQualIntoStsUserObject (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  Int2           found, j;
  CharPtr        str;
  UserObjectPtr  uop;

  found = 0;
  for (j = 0, str = stsQualList [j]; str != NULL; j++, str = stsQualList [j]) {
    if (StringICmp (qual, str) == 0) {
      found = j;
    }
  }

  if (found > 0) {
    uop = GetStsUserObject (sfp);
    if (uop == NULL) return FALSE;
    AddFieldToSnpStsCloneUserObject (uop, qual, val);
    return TRUE;
  }
  return FALSE;
}

static CharPtr cloneQualList [] = {
  "", "clone_id", "method", "sequence", "bac_ends", "STS", "weight", NULL
};

static UserObjectPtr CreateCloneUserObject (void)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  uop = UserObjectNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("cloneUserObject");
  uop->type = oip;

  return uop;
}

static UserObjectPtr GetCloneUserObject (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (sfp == NULL) return NULL;
  if (sfp->ext == NULL) {
    sfp->ext = CreateCloneUserObject ();
  }
  uop = sfp->ext;
  if (uop == NULL) return NULL;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "cloneUserObject") != 0) return NULL;
  return uop;
}

static Boolean ParseQualIntoCloneUserObject (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  Int2           found, j;
  CharPtr        str;
  UserObjectPtr  uop;

  found = 0;
  for (j = 0, str = cloneQualList [j]; str != NULL; j++, str = cloneQualList [j]) {
    if (StringICmp (qual, str) == 0) {
      found = j;
    }
  }

  if (found > 0) {
    uop = GetCloneUserObject (sfp);
    if (uop == NULL) return FALSE;
    AddFieldToSnpStsCloneUserObject (uop, qual, val);
    return TRUE;
  }
  return FALSE;
}

/* gene ontology user object parsing */

static CharPtr goQualList [] = {
  "", "go_process", "go_component", "go_function", NULL
};

static CharPtr goQualType [] = {
  "", "Process", "Component", "Function", NULL
};

/* later will need to be able to deal with CombinedFeatureUserObjects */

static UserObjectPtr GetGeneOntologyUserObject (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (sfp == NULL) return NULL;
  if (sfp->ext == NULL) {
    sfp->ext = CreateGeneOntologyUserObject ();
  }
  uop = sfp->ext;
  if (uop == NULL) return NULL;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "GeneOntology") != 0) return NULL;
  return uop;
}

static Boolean ParseQualIntoGeneOntologyUserObject (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  CharPtr        fields [4];
  Int2           found, j;
  long int       num;
  Int4           pmid = 0;
  CharPtr        str, ptr, tmp, goref = NULL;
  UserObjectPtr  uop;

  found = 0;
  for (j = 0, str = goQualList [j]; str != NULL; j++, str = goQualList [j]) {
    if (StringICmp (qual, str) == 0) {
      found = j;
    }
  }

  if (found > 0) {
    uop = GetGeneOntologyUserObject (sfp);
    if (uop == NULL) return FALSE;
    str = StringSave (val);
    for (j = 0; j < 4; j++) {
      fields [j] = NULL;
    }
    ptr = str;
    for (j = 0; j < 4 && ptr != NULL; j++) {
      fields [j] = ptr;
      TrimSpacesAroundString (ptr);
      ptr = StringChr (ptr, '|');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
    }
    tmp = fields [1];
    if (tmp != NULL && *tmp != '\0') {
      if (StringNICmp (tmp, "GO:", 3) == 0) {
        fields [1] = tmp + 3;
      }
    }
    tmp = fields [2];
    if (tmp != NULL && *tmp != '\0') {
      if (StringNICmp (tmp, "GO_REF:", 7) == 0) {
        fields [2] = tmp + 7;
      }
    }
    if (fields [2] != NULL && sscanf (fields [2], "%ld", &num) == 1) {
      pmid = (Int4) num;
      tmp = fields [2];
      if (*tmp == '0') {
        pmid = 0;
        goref = tmp;
      }
    }
    AddToGeneOntologyUserObject (uop, goQualType [found], fields [0],
                                 fields [1], pmid, goref, fields [3]);
    MemFree (str);
    return TRUE;
  }
  return FALSE;
}

static CharPtr okayCategoryPrefixes [] = {
  "",
  "COORDINATES:",
  "DESCRIPTION:",
  "EXISTENCE:",
  NULL
};

static CharPtr okayInferencePrefixes [] = {
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

static Boolean InvalidInference (CharPtr str)

{
  Int2    best, j;
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return TRUE;

  for (j = 0; okayCategoryPrefixes [j] != NULL; j++) {
    len = StringLen (okayCategoryPrefixes [j]);
    if (StringNICmp (str, okayCategoryPrefixes [j], len) != 0) continue;
    str += len;
    ch = *str;
    while (ch == ' ') {
      str++;
      ch = *str;
    }
    break;
  }

  if (StringHasNoText (str)) return TRUE;

  best = -1;
  for (j = 0; okayInferencePrefixes [j] != NULL; j++) {
    len = StringLen (okayInferencePrefixes [j]);
    if (StringNICmp (str, okayInferencePrefixes [j], len) != 0) continue;
    best = j;
  }
  if (best >= 0 && okayInferencePrefixes [best] != NULL) return FALSE;

  return TRUE;
}

static void ParseCodonRecognized (CharPtr val, tRNAPtr trp)

{
  Char        buf [256];
  Char        codon [16];
  ValNodePtr  head = NULL;
  Int2        i;
  Int2        j;
  CharPtr     ptr;
  CharPtr     str;
  tRNA        tr;
  ValNodePtr  vnp;

  if (trp == NULL) return;
  for (j = 0; j < 6; j++) {
    trp->codon [j] = 255;
  }
  if (StringHasNoText (val)) return;

  MemSet ((Pointer) &tr, 0, sizeof (tRNA));

  StringNCpy_0 (buf, val, sizeof (buf));
  str = buf;
  while (StringDoesHaveText (str)) {
    ptr = StringChr (str, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    if (StringDoesHaveText (str)) {
      for (j = 0; j < 6; j++) {
        tr.codon [j] = 255;
      }
      StringCpy (codon, str);
      for (i = 0; i < 3; i++) {
        if (codon [i] == 'U') {
          codon [i] = 'T';
        }
      }
      ParseDegenerateCodon (&tr, (Uint1Ptr) codon);
      for (i = 0; i < 6; i++) {
        if (tr.codon [i] == 255) continue;
        ValNodeAddInt (&head, 0, (long) tr.codon [i]);
      }
    }
    str = ptr;
  }
  if (head == NULL) return;

  head = ValNodeSort (head, SortByIntvalue);
  head = UniqueIntValNode (head);
  for (vnp = head, j = 0; vnp != NULL && j < 6; vnp = vnp->next, j++) {
    trp->codon [j] = (Uint1) vnp->data.intvalue;
  }
}

static UserObjectPtr CreateNomenclatureUserObject (
  void
)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  uop = UserObjectNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("OfficialNomenclature");
  uop->type = oip;

  return uop;
}

static UserObjectPtr GetNomenclatureUserObject (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (sfp == NULL) return NULL;
  if (sfp->ext == NULL) {
    sfp->ext = CreateNomenclatureUserObject ();
  }
  uop = sfp->ext;
  if (uop == NULL) return NULL;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "OfficialNomenclature") != 0) return NULL;
  return uop;
}

static void AddToNomenclatureUserObject (
  UserObjectPtr uop,
  CharPtr status,
  CharPtr symbol,
  CharPtr name,
  CharPtr source
)

{
  UserFieldPtr  last = NULL;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp;

  if (uop == NULL) return;
  if (StringHasNoText (symbol)) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "OfficialNomenclature") != 0) return;

  ufp = UserFieldNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("Symbol");
  ufp->label = oip;
  ufp->choice = 1; /* visible string */
  ufp->data.ptrvalue = (Pointer) StringSave (symbol);

  uop->data = ufp;
  last = ufp;

  if (StringDoesHaveText (name)) {
    ufp = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->str = StringSave ("Name");
    ufp->label = oip;
    ufp->choice = 1; /* visible string */
    ufp->data.ptrvalue = (Pointer) StringSave (name);
    last->next = ufp;
    last = ufp;
  }

  if (StringDoesHaveText (source)) {
    ufp = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->str = StringSave ("DataSource");
    ufp->label = oip;
    ufp->choice = 1; /* visible string */
    ufp->data.ptrvalue = (Pointer) StringSave (source);
    last->next = ufp;
    last = ufp;
  }

  if (StringDoesHaveText (status)) {
    ufp = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->str = StringSave ("Status");
    ufp->label = oip;
    ufp->choice = 1; /* visible string */
    if (StringICmp (status, "Official") == 0) {
      ufp->data.ptrvalue = (Pointer) StringSave ("Official");
    } else if (StringICmp (status, "Interim") == 0) {
      ufp->data.ptrvalue = (Pointer) StringSave ("Interim");
    } else {
      ufp->data.ptrvalue = (Pointer) StringSave ("?");
    }
    last->next = ufp;
    last = ufp;
  }
}

static void ParseQualIntoNomenclatureUserObject (SeqFeatPtr sfp, CharPtr val)

{
  CharPtr        fields [4];
  Int2           j;
  CharPtr        str, ptr;
  UserObjectPtr  uop;

  if (sfp == NULL) return;
  if (StringHasNoText (val)) return;

  str = StringSave (val);
  for (j = 0; j < 4; j++) {
    fields [j] = NULL;
  }
  ptr = str;
  for (j = 0; j < 4 && ptr != NULL; j++) {
    fields [j] = ptr;
    TrimSpacesAroundString (ptr);
    ptr = StringChr (ptr, '|');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
  }

  uop = GetNomenclatureUserObject (sfp);
  AddToNomenclatureUserObject (uop, fields [0], fields [1], fields [2], fields [3]);

  MemFree (str);
}

static void TrailingCommaFix (CharPtr str)

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
    str [len - 1] = '_';
    str [len] = '\0';
  }
}

static CharPtr singletonList [] = {
  "artificial location",
  "artificial-location",
  "artificial_location",
  "exception",
  "mitochondrion",
  "order",
  "pseudo",
  "ribosomal slippage",
  "ribosomal-slippage",
  "ribosomal_slippage",
  "trans splicing",
  "trans-splicing",
  "trans_splicing"
};

static Boolean IsSingletonQual (CharPtr str)

{
  Int2  L, R, mid;

  if (str == NULL || *str == '\0') return FALSE;

  L = 0;
  R = (sizeof (singletonList) / sizeof (CharPtr)) - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringCmp (singletonList [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringCmp (singletonList [R], str) == 0) {
    return TRUE;
  }

  return FALSE;
}

static void AddQualifierToFeatureEx (SeqFeatPtr sfp, CharPtr qual, CharPtr val, Int4 offset, Int4 lin_num)

{
  Uint1           aa;
  AffilPtr        affil;
  AuthListPtr     alp;
  Boolean         bail;
  Uint1           codon [6];
  CdRegionPtr     crp;
  CitSubPtr       csp;
  DbtagPtr        db;
  GBQualPtr       gbq;
  GeneRefPtr      grp;
  ImpFeatPtr      ifp = NULL;
  Boolean         isGeneDesc = FALSE;
  Boolean         isGeneSyn = FALSE;
  Boolean         isLocusTag = FALSE;
  Boolean         isNomenclature = FALSE;
  Boolean         isCytMap = FALSE;
  Boolean         isGenMap = FALSE;
  Boolean         isRadMap = FALSE;
  Boolean         isAuthor = FALSE;
  Boolean         isAffil = FALSE;
  Boolean         isMuid = FALSE;
  Boolean         isPmid = FALSE;
  Int2            j;
  Boolean         justTrnaText;
  GBQualPtr       last;
  size_t          len;
  int             num;
  ObjectIdPtr     oip;
  PubdescPtr      pdp;
  ProtRefPtr      prp = NULL;
  CharPtr         ptr;
  Int2            qnum;
  RnaRefPtr       rrp;
  CharPtr         str;
  CharPtr         tag;
  tRNAPtr         trna;
  long            uid;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || HasNoText (qual)) return;
  if (HasNoText (val)) {
    if (! IsSingletonQual (qual)) return;
    
    if (StringICmp (qual, "artificial location") == 0 || StringICmp (qual, "artificial-location") == 0) {
      qual = "artificial_location";
    } else if (StringICmp (qual, "ribosomal slippage") == 0 || StringICmp (qual, "ribosomal-slippage") == 0) {
      qual = "ribosomal_slippage";
    } else if (StringICmp (qual, "trans splicing") == 0 || StringICmp (qual, "trans-splicing") == 0) {
      qual = "trans_splicing";
    }
  }
  qnum = GBQualNameValid (qual);
  if (qnum <= -1) {
    if (StringNCmp (qual, "gene_syn", 8) == 0) {
      qnum = GBQUAL_gene;
      isGeneSyn = TRUE;
    } else if (StringNCmp (qual, "gene_desc", 9) == 0) {
      qnum = GBQUAL_gene;
      isGeneDesc = TRUE;
    } else if (StringNCmp (qual, "locus_tag", 9) == 0) {
      qnum = GBQUAL_gene;
      isLocusTag = TRUE;
    } else if (StringNCmp (qual, "nomenclature", 12) == 0) {
      qnum = GBQUAL_gene;
      isNomenclature = TRUE;
    } else if (StringNCmp (qual, "gen_map", 7) == 0) {
      qnum = GBQUAL_gene;
      isGenMap = TRUE;
    } else if (StringNCmp (qual, "cyt_map", 7) == 0) {
      qnum = GBQUAL_gene;
      isCytMap = TRUE;
    } else if (StringNCmp (qual, "rad_map", 7) == 0) {
      qnum = GBQUAL_gene;
      isRadMap = TRUE;
    } else if (sfp->data.choice == SEQFEAT_PUB) {
      if (StringICmp (qual, "pmid") == 0 || StringICmp (qual, "PubMed") == 0) {
        isPmid = TRUE;
      } else if (StringICmp (qual, "muid") == 0 || StringICmp (qual, "MEDLINE") == 0) {
        isMuid = TRUE;
      } else if (StringICmp (qual, "Author") == 0) {
        isAuthor = TRUE;
      } else if (StringICmp (qual, "Affil") == 0 || StringICmp (qual, "Affiliation") == 0) {
        isAffil = TRUE;
      }
    }
  }
  if (qnum == GBQUAL_evidence) {
    qnum = -1; /* no longer legal */
  }
  if (qnum <= -1) {
    bail = TRUE;
    if (sfp->data.choice == SEQFEAT_IMP) {
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue; /* for variation user object */
    }
    if (sfp->data.choice == SEQFEAT_REGION && (StringCmp (qual, "region") == 0 || StringCmp (qual, "region_name") == 0)) {
      sfp->data.value.ptrvalue = MemFree (sfp->data.value.ptrvalue);
      sfp->data.value.ptrvalue = StringSave (val);
    } else if (sfp->data.choice == SEQFEAT_BOND && StringCmp (qual, "bond_type") == 0) {
      StripHyphens (val);
      sfp->data.value.intvalue = 255;
      for (j = 0; bondList [j] != NULL; j++) {
        if (StringNICmp (val, bondList [j], StringLen (bondList [j])) == 0) {
          sfp->data.value.intvalue = j;
        }
      }
    } else if (sfp->data.choice == SEQFEAT_SITE && StringCmp (qual, "site_type") == 0) {
      StripHyphens (val);
      sfp->data.value.intvalue = 255;
      for (j = 0; siteList [j] != NULL; j++) {
        if (StringNICmp (val, siteList [j], StringLen (siteList [j])) == 0) {
          sfp->data.value.intvalue = j;
        }
      }
    } else if (sfp->data.choice == SEQFEAT_PUB) {
      if (isPmid) {
        if (sscanf (val, "%ld", &uid) == 1) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          if (pdp != NULL) {
            ValNodeAddInt (&(pdp->pub), PUB_PMid, (Int4) uid);
          }
        }
      } else if (isMuid) {
        if (sscanf (val, "%ld", &uid) == 1) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          if (pdp != NULL) {
            ValNodeAddInt (&(pdp->pub), PUB_Muid, (Int4) uid);
          }
        }
      } else if (isAuthor || isAffil) {
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        csp = NULL;
        if (pdp != NULL) {
          for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
            if (vnp->choice == PUB_Sub) {
              csp = (CitSubPtr) vnp->data.ptrvalue;
              break;
            }
          }
          if (csp == NULL) {
            csp = CitSubNew ();
            if (csp != NULL) {
              csp->date = DateCurr ();
              ValNodeAddPointer (&(pdp->pub), PUB_Sub, (Pointer) csp);
            }
          }
          if (csp != NULL) {
            alp = csp->authors;
            if (alp == NULL) {
              alp = AuthListNew ();
              if (alp != NULL) {
                alp->choice = 3;
                csp->authors = alp;
              }
            }
            if (alp != NULL) {
              if (isAuthor) {
                ValNodeCopyStr (&(alp->names), 3, val);
              } else if (isAffil) {
                affil = alp->affil;
                if (affil == NULL) {
                  affil = AffilNew ();
                  alp->affil = affil;
                }
                if (affil != NULL) {
                  affil->choice = 1;
                  affil->affil = StringSave (val);
                }
              }
            }
          }
        }
      }
    } else if (sfp->data.choice == SEQFEAT_PUB && (StringICmp (qual, "muid") == 0 || StringICmp (qual, "MEDLINE") == 0)) {
    } else if (sfp->data.choice == SEQFEAT_BIOSRC && ParseQualIntoBioSource (sfp, qual, val)) {
    } else if (sfp->data.choice == SEQFEAT_CDREGION && StringCmp (qual, "prot_desc") == 0) {
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
        xref = xref->next;
      }
      if (xref == NULL) {
        prp = ProtRefNew ();
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          xref->data.choice = SEQFEAT_PROT;
          xref->data.value.ptrvalue = (Pointer) prp;
          xref->next = sfp->xref;
          sfp->xref = xref;
        }
      }
      if (xref != NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
      }
      if (prp == NULL) return;
      prp->desc = MemFree (prp->desc);
      prp->desc = StringSaveNoNull (val);
    } else if (sfp->data.choice == SEQFEAT_CDREGION && StringCmp (qual, "prot_note") == 0) {
      bail = FALSE;
    } else if (sfp->data.choice == SEQFEAT_PROT && StringCmp (qual, "prot_desc") == 0) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp != NULL) {
        prp->desc = MemFree (prp->desc);
        prp->desc = StringSaveNoNull (val);
      }
    } else if (sfp->data.choice == SEQFEAT_CDREGION && StringCmp (qual, "secondary_accession") == 0) {
      bail = FALSE;
    } else if (sfp->data.choice == SEQFEAT_RNA &&
               (StringCmp (qual, "codon_recognized") == 0 || StringCmp (qual, "codons_recognized") == 0)) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 3) {
        if (rrp->ext.choice == 0 && rrp->ext.value.ptrvalue == NULL) {
          rrp->ext.choice = 2;
          trna = (tRNAPtr) MemNew (sizeof (tRNA));
          rrp->ext.value.ptrvalue = (Pointer) trna;
          if (trna != NULL) {
            trna->aatype = 2;
            for (j = 0; j < 6; j++) {
              trna->codon [j] = 255;
            }
          }
        }
        trna = (tRNAPtr) rrp->ext.value.ptrvalue;
        ParseCodonRecognized (val, trna);
        /*
        StringNCpy_0 ((CharPtr) codon, val, sizeof (codon));
        if (StringLen ((CharPtr) codon) == 3) {
          for (j = 0; j < 3; j++) {
            if (codon [j] == 'U') {
              codon [j] = 'T';
            }
          }
          if (trna != NULL) {
            ParseDegenerateCodon (trna, (Uint1Ptr) codon);
          }
        }
        */
      }
    } else if (ifp != NULL && StringICmp (ifp->key, "variation") == 0 && ParseQualIntoSnpUserObject (sfp, qual, val)) {
    } else if (ifp != NULL && StringICmp (ifp->key, "STS") == 0 && ParseQualIntoStsUserObject (sfp, qual, val)) {
    } else if (ifp != NULL && StringICmp (ifp->key, "misc_feature") == 0 && ParseQualIntoCloneUserObject (sfp, qual, val)) {
    } else if ((sfp->data.choice == SEQFEAT_GENE ||
                sfp->data.choice == SEQFEAT_CDREGION ||
                sfp->data.choice == SEQFEAT_RNA) &&
               ParseQualIntoGeneOntologyUserObject (sfp, qual, val)) {
    } else if (sfp->data.choice == SEQFEAT_RNA && StringCmp (qual, "comment") == 0) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 2) {
        bail = FALSE;
      }
    } else {
      if (lin_num > 0) {
        ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_UnknownImpFeatQual, "Unknown qualifier '%s', relative line %ld", qual, (long) lin_num);
      } else {
        ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_UnknownImpFeatQual, "Unknown qualifier '%s'", qual);
      }
    }
    if (bail) return;
  }
  if (qnum == GBQUAL_note) {
    if (sfp->comment == NULL) {
      sfp->comment = StringSave (val);
    } else {
      len = StringLen (sfp->comment) + StringLen (val) + 5;
      str = MemNew (sizeof (Char) * len);
      StringCpy (str, sfp->comment);
      /*
      StringCat (str, "; ");
      */
      StringCat (str, "~");
      StringCat (str, val);
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = str;
    }
    return;
  } else if (qnum == GBQUAL_pseudo) {
    sfp->pseudo = TRUE;
    return;
  } else if ((qnum == GBQUAL_gene || qnum == GBQUAL_locus_tag) && sfp->data.choice != SEQFEAT_GENE) {
    if (StringCmp (val, "-") == 0) {
      val = NULL;
    }
    xref = sfp->xref;
    while (xref != NULL && xref->data.choice != SEQFEAT_GENE) {
      xref = xref->next;
    }
    if (xref == NULL) {
      grp = GeneRefNew ();
      xref = SeqFeatXrefNew ();
      if (xref != NULL) {
        xref->data.choice = SEQFEAT_GENE;
        xref->data.value.ptrvalue = (Pointer) grp;
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
    }
    if (xref != NULL) {
      grp = (GeneRefPtr) xref->data.value.ptrvalue;
      if (grp == NULL) return;
      if (isGeneSyn) {
        ValNodeCopyStr (&(grp->syn), 0, val);
      } else if (isGeneDesc) {
        grp->desc = StringSave (val);
      } else if (isLocusTag || qnum == GBQUAL_locus_tag) {
        grp->locus_tag = StringSave (val);
      } else if (grp->locus == NULL) {
        grp->locus = StringSave (val);
      } else {
        ValNodeCopyStr (&(grp->syn), 0, val);
      }
    }
    return;
  } else if (qnum == GBQUAL_db_xref) {
    if (StringICmp (val, "GI") == 0) {
      ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Reserved db_xref value %s", val);
      return;
    }
    if (StringICmp (val, "NID") == 0 ||
        StringICmp (val, "PID") == 0 ||
        StringICmp (val, "PIDg") == 0 ||
        StringICmp (val, "PIDe") == 0 ||
        StringICmp (val, "PIDd") == 0) {
      ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Obsolete db_xref value %s", val);
      return;
    }

    vnp = ValNodeNew (NULL);
    db = DbtagNew ();
    vnp->data.ptrvalue = db;
    tag = val;
    ptr = StringChr (tag, ':');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      db->db = StringSave (tag);
      oip = ObjectIdNew ();
      oip->str = StringSave (ptr);
      db->tag = oip;
    } else {
      db->db = StringSave ("?");
      oip = ObjectIdNew ();
      oip->str = StringSave (tag);
      db->tag = oip;
    }
    if (sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      vnp->next = grp->db;
      grp->db = vnp;
    } else {
      vnp->next = sfp->dbxref;
      sfp->dbxref = vnp;
    }
    return;
  } else if (qnum == GBQUAL_replace && StringCmp (val, "-") == 0) {
    val = "";
  } else if (qnum == GBQUAL_evidence) {
    /*
    if (StringICmp (val, "experimental") == 0) {
      sfp->exp_ev = 1;
    } else if (StringICmp (val, "not_experimental") == 0 ||
               StringICmp (val, "non_experimental") == 0 ||
               StringICmp (val, "not-experimental") == 0 ||
               StringICmp (val, "non-experimental") == 0) {
      sfp->exp_ev = 2;
    }
    */
    return;
  } else if (qnum == GBQUAL_exception) {
    sfp->excpt = TRUE;
    if (! HasNoText (val)) {
      sfp->except_text = StringSave (val);
    }
    return;
  }

  if (qnum == GBQUAL_old_locus_tag || qnum == GBQUAL_experiment) {

    /* fall through to add as gbqual */

  } else if (qnum == GBQUAL_inference) {

    if (InvalidInference (val)) {
      ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_InvalidQualifierValue, "Invalid inference value %s", val);
      return;
    }

  } else if (sfp->data.choice == SEQFEAT_GENE) {
    if (qnum == GBQUAL_gene || qnum == GBQUAL_allele || qnum == GBQUAL_map || qnum == GBQUAL_locus_tag) {
      if (qnum == GBQUAL_gene) {
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        if (grp != NULL) {
          if (isGeneSyn) {
            ValNodeCopyStr (&(grp->syn), 0, val);
          } else if (isGeneDesc) {
            grp->desc = StringSave (val);
          } else if (isLocusTag) {
            grp->locus_tag = StringSave (val);
          } else if (isGenMap || isCytMap || isRadMap) {
            /* fall through to add as gbqual */
          } else if (isNomenclature) {
            ParseQualIntoNomenclatureUserObject (sfp, val);
          } else if (grp->locus == NULL) {
            grp->locus = StringSave (val);
          } else {
            ValNodeCopyStr (&(grp->syn), 0, val);
          }
        }
      } else if (qnum == GBQUAL_allele) {
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        if (grp != NULL) {
          grp->allele = StringSave (val);
        }
      } else if (qnum == GBQUAL_map) {
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        if (grp != NULL) {
          grp->maploc = StringSave (val);
        }
      } else if (qnum == GBQUAL_locus_tag) {
        grp = (GeneRefPtr) sfp->data.value.ptrvalue;
        if (grp != NULL) {
          grp->locus_tag = StringSave (val);
        }
      }
      if (isGenMap || isCytMap || isRadMap) {
        /* fall through to add as gbqual */
      } else {
        return;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    if (qnum == GBQUAL_function || qnum == GBQUAL_EC_number || qnum == GBQUAL_product) {
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
        xref = xref->next;
      }
      if (xref == NULL) {
        prp = ProtRefNew ();
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          xref->data.choice = SEQFEAT_PROT;
          xref->data.value.ptrvalue = (Pointer) prp;
          xref->next = sfp->xref;
          sfp->xref = xref;
        }
      }
      if (xref != NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
      }
      if (prp == NULL) return;
      if (qnum == GBQUAL_function) {
        ValNodeCopyStr (&(prp->activity), 0, val);
      } else if (qnum == GBQUAL_EC_number) {
        ValNodeCopyStr (&(prp->ec), 0, val);
      } else if (qnum == GBQUAL_product) {
        TrailingCommaFix (val);
        ValNodeCopyStr (&(prp->name), 0, val);
      }
      return;
    } else if (qnum == GBQUAL_transl_except) {
      if (ParseCodeBreak (sfp, val, offset)) return;
    } else if (qnum == GBQUAL_codon_start) {
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (sscanf (val, "%d", &num) == 1 && crp != NULL) {
        if (num > 0 && num < 4) {
          crp->frame = (Uint1) num;
        }
      }
      return;
    }
  } else if (sfp->data.choice == SEQFEAT_PROT) {
    if (qnum == GBQUAL_function || qnum == GBQUAL_EC_number || qnum == GBQUAL_product) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp != NULL) {
        if (qnum == GBQUAL_function) {
          ValNodeCopyStr (&(prp->activity), 0, val);
        } else if (qnum == GBQUAL_EC_number) {
          ValNodeCopyStr (&(prp->ec), 0, val);
        } else if (qnum == GBQUAL_product) {
          TrailingCommaFix (val);
          ValNodeCopyStr (&(prp->name), 0, val);
        }
        return;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    if (qnum == GBQUAL_product) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) return;
      if (rrp->type == 3) {
        aa = ParseTRnaString (val, &justTrnaText, codon, FALSE);
        if (aa != 0) {
          if (rrp->ext.choice == 0 && rrp->ext.value.ptrvalue == NULL) {
            rrp->ext.choice = 2;
            trna = (tRNAPtr) MemNew (sizeof (tRNA));
            rrp->ext.value.ptrvalue = (Pointer) trna;
            if (trna != NULL) {
              trna->aatype = 2;
              for (j = 0; j < 6; j++) {
                trna->codon [j] = 255;
              }
            }
          }
          trna = (tRNAPtr) rrp->ext.value.ptrvalue;
          if (trna != NULL) {
            if (justTrnaText) {
              for (j = 0; j < 6; j++) {
                trna->codon [j] = codon [j];
              }
            }
            trna->aa = aa;
          }
          if (aa == 'M') {
            if (StringStr (val, "fMet") != NULL) {
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
        } else {
          if (sfp->comment == NULL) {
            sfp->comment = StringSave (val);
          } else {
            len = StringLen (sfp->comment) + StringLen (val) + 5;
            str = MemNew (sizeof (Char) * len);
            StringCpy (str, sfp->comment);
            StringCat (str, "; ");
            StringCat (str, val);
            sfp->comment = MemFree (sfp->comment);
            sfp->comment = str;
          }
        }
        return;
      } else if (rrp->type != 255) {
        if (rrp->ext.choice == 1) {
          /*
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          */
          if (sfp->comment == NULL) {
            sfp->comment = StringSave (val);
          } else {
            len = StringLen (sfp->comment) + StringLen (val) + 5;
            str = MemNew (sizeof (Char) * len);
            StringCpy (str, sfp->comment);
            StringCat (str, "; ");
            StringCat (str, val);
            sfp->comment = MemFree (sfp->comment);
            sfp->comment = str;
          }
        } else {
          rrp->ext.choice = 1;
          TrailingCommaFix (val);
          rrp->ext.value.ptrvalue = StringSave (val);
        }
        return;
      }
    } else if (qnum == GBQUAL_anticodon) {
      if (ParseAnticodon (sfp, val, offset)) return;
    }
  } else if (sfp->data.choice == SEQFEAT_BIOSRC) {
    if (ParseQualIntoBioSource (sfp, qual, val)) return;
  }

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

NLM_EXTERN void AddQualifierToFeature (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  AddQualifierToFeatureEx (sfp, qual, val, 0, 0);
}

NLM_EXTERN SeqLocPtr AddIntervalToLocation (SeqLocPtr loc, SeqIdPtr sip,
                                            Int4 start, Int4 stop,
                                            Boolean partial5, Boolean partial3)

{
  Int4        flip;
  IntFuzzPtr  ifp;
  Boolean     is_first;
  SeqLocPtr   rsult = NULL;
  SeqIntPtr   sintp;
  SeqLocPtr   slp;
  Uint1       strand;
  SeqLocPtr   tmp;

  if (sip == NULL) return NULL;

  sintp = SeqIntNew ();
  strand = Seq_strand_plus;
  if (start > stop) {
    flip = start;
    start = stop;
    stop = flip;
    strand = Seq_strand_minus;
  }
  sintp->from = start;
  sintp->to = stop;
  sintp->strand = strand;
  sintp->id = SeqIdDup (sip);

  if (partial5) {
    ifp = IntFuzzNew ();
    if (ifp != NULL) {
      ifp->choice = 4;
      if (strand == Seq_strand_minus || strand == Seq_strand_both_rev) {
        sintp->if_to = ifp;
        ifp->a = 1;
      } else {
        sintp->if_from = ifp;
        ifp->a = 2;
      }
    }
  }

  if (partial3) {
    ifp = IntFuzzNew ();
    if (ifp != NULL) {
      ifp->choice = 4;
      if (strand == Seq_strand_minus || strand == Seq_strand_both_rev) {
        sintp->if_from = ifp;
        ifp->a = 2;
      } else {
        sintp->if_to = ifp;
        ifp->a = 1;
      }
    }
  }

  slp = ValNodeAddPointer (NULL, SEQLOC_INT, (Pointer) sintp);

  if (loc == NULL) return slp;

  if (loc->choice == SEQLOC_MIX) {
    tmp = (ValNodePtr) (loc->data.ptrvalue);
    while (tmp->next != NULL) {
      tmp = tmp->next;
    }
    tmp->next = slp;
    rsult = loc;
  } else {
    tmp = ValNodeNew (NULL);
    tmp->choice = SEQLOC_MIX;
    tmp->data.ptrvalue = (Pointer) loc;
    loc->next = slp;
    rsult = tmp;
  }

  if (SeqLocStrand (rsult) == Seq_strand_other) {
    is_first = TRUE;
    slp = SeqLocFindNext (rsult, NULL);
    while (slp != NULL) {
      if (slp->choice == SEQLOC_INT) {
        sintp = (SeqIntPtr) slp->data.ptrvalue;
        if (sintp != NULL) {
          /* exon of one base in feature */
          if (sintp->from == sintp->to) {
            sintp->strand = Seq_strand_minus;
            if (is_first) {
              ifp = sintp->if_from;
              if (ifp != NULL && ifp->choice == 4 && ifp->a == 2 && sintp->if_to == NULL) {
                sintp->if_from = NULL;
                sintp->if_to = ifp;
                ifp->a = 1;
              }
            } else if (partial3) {
              ifp = sintp->if_to;
              if (ifp != NULL && ifp->choice == 4 && ifp->a == 1 && sintp->if_from == NULL) {
                sintp->if_to = NULL;
                sintp->if_from = ifp;
                ifp->a = 2;
              }
            }
          }
        }
      }
      slp = SeqLocFindNext (rsult, slp);
      is_first = FALSE;
    }
  }

  return rsult;
}

static void PutNullsBetween (SeqLocPtr loc)

{
  SeqLocPtr  next;
  SeqLocPtr  tmp;
  SeqLocPtr  vnp;

  if (loc == NULL) return;
  if (loc->choice != SEQLOC_MIX) return;

  vnp = (ValNodePtr) (loc->data.ptrvalue);
  while (vnp != NULL && vnp->next != NULL) {
    next = vnp->next;
    tmp = ValNodeNew (NULL);
    if (tmp != NULL) {
      tmp->choice = SEQLOC_NULL;
      tmp->next = vnp->next;
      vnp->next = tmp;
    }
    vnp = next;
  }
}

static CharPtr TokenizeAtWhiteSpace (CharPtr str)

{
  Char     ch;
  CharPtr  ptr;

  if (str == NULL) return NULL;
  ptr = str;
  ch = *ptr;

  while (ch != '\0' && (IS_WHITESP (ch))) {
    ptr++;
    ch = *ptr;
  }
  while (ch != '\0' && (! IS_WHITESP (ch))) {
    ptr++;
    ch = *ptr;
  }
  if (ch != '\0') {
    *ptr = '\0';
    ptr++;
  }

  return ptr;
}

static void ParseWhitespaceIntoTabs (CharPtr line)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (StringHasNoText (line)) return;
  len = StringLen (line) + 10;

  str = MemNew (len);
  if (str == NULL) return;

  ptr = line;
  ch = *ptr;
  if (IS_WHITESP (ch)) {
    /* qualifier value line */
    StringCat (str, "\t\t\t");
    TrimSpacesAroundString (ptr);
    tmp = TokenizeAtWhiteSpace (ptr);
    if (tmp != NULL) {
      while (isspace (*tmp)) {
        tmp++;
      }
    }
    StringCat (str, ptr);
    StringCat (str, "\t");
    StringCat (str, tmp);
  } else {
    /* location and possible feature key line */
    TrimSpacesAroundString (ptr);
    tmp = TokenizeAtWhiteSpace (ptr);
    StringCat (str, ptr);
    StringCat (str, "\t");
    ptr = tmp;
    tmp = TokenizeAtWhiteSpace (ptr);
    StringCat (str, ptr);
    ptr = tmp;
    if (! StringHasNoText (ptr)) {
      tmp = TokenizeAtWhiteSpace (ptr);
      StringCat (str, "\t");
      StringCat (str, ptr);
    }
  }

  /* replace original with tab-delimited table */
  StringCpy (line, str);

  MemFree (str);
}


static CharPtr ReadTheRestOfTheLine (FileCachePtr fcp, CharPtr original_buffer)
{
  Char         line [2047];
  CharPtr      str;
  Boolean      nonewline = TRUE;
  ValNodeBlock extra;
  Int4         len = 1;
  ValNodePtr   vnp;

  InitValNodeBlock(&extra, NULL);
  ValNodeAddPointerToEnd (&extra, 0, StringSave(original_buffer));
  len += StringLen (original_buffer);
  while (nonewline) {
    nonewline = FALSE;
    str = FileCacheReadLine (fcp, line, sizeof (line), &nonewline);
    if (str == NULL) {
      nonewline = FALSE;
    } else {
      ValNodeAddPointerToEnd (&extra, 0, StringSave (line));
      len += StringLen (line);
    }
  }
  str = (CharPtr) MemNew (sizeof (Char) * len);
  *str = 0;
  for (vnp = extra.head; vnp != NULL; vnp = vnp->next) {
    StringCat(str, vnp->data.ptrvalue);
  }
  str[len - 1] = 0;
  return str;
}


static SeqAnnotPtr ReadFeatureTableEx (FileCachePtr fcp, CharPtr seqid, CharPtr annotname, Int4Ptr p_line)

{
  Boolean        allowWhitesp  = TRUE;
  BioSourcePtr   biop;
  Char           buf [128];
  CdRegionPtr    crp;
  AnnotDescrPtr  desc;
  Boolean        endsinspace;
  CharPtr        feat;
  IntFuzzPtr     fuzz;
  GeneRefPtr     grp;
  Int2           idx;
  ImpFeatPtr     ifp;
  Boolean        inquals = FALSE;
  Boolean        isminus;
  Boolean        isnote;
  Boolean        ispoint;
  Int2           j;
  CharPtr        label;
  size_t         len;
  Char           line [2047];
  Int4           lin_num = 1;
  CharPtr        loc;
  Boolean        nonewline;
  long int       num;
  Int4           offset = 0;
  OrgRefPtr      orp;
  Boolean        partial5;
  Boolean        partial3;
  PubdescPtr     pdp;
  Int4           pos;
  SeqFeatPtr     prev = NULL;
  ProtRefPtr     prp;
  CharPtr        qual;
  Uint1          rnatype;
  RnaRefPtr      rrp;
  SeqAnnotPtr    sap = NULL;
  SeqFeatPtr     sfp = NULL;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  SeqPntPtr      spp;
  Int4           start;
  Int4           stop;
  SqnTagPtr      stp;
  CharPtr        str;
  CharPtr        tmp;
  CharPtr        val;
  ValNodePtr     vslp;
  Boolean        free_str = FALSE;

  if (fcp == NULL || fcp->fp == NULL || seqid == NULL) return NULL;
  sip = SeqIdFindBest (MakeSeqID (seqid), 0);
  if (sip == NULL) return NULL;

  pos = FileCacheTell (fcp);
  str = FileCacheReadLine (fcp, line, sizeof (line), &nonewline);
  if (nonewline) {
    str = ReadTheRestOfTheLine (fcp, line);
    if (StringDoesHaveText (str)) {
      free_str = TRUE;
    } else {
      str = MemFree (str);
    }
  }

  if (p_line != NULL) {
    lin_num = *p_line;
  }
  lin_num++;

  while (str != NULL) {

    isnote = FALSE;
    endsinspace = FALSE;
    len = StringLen (str);
    if (len > 2 && str [len - 1] == ' ') {
      endsinspace = TRUE;
    }

    if (! HasNoText (str)) {

      if (StringNCmp (str, ">", 1) == 0 ||
          StringNCmp (str, "LOCUS ", 6) == 0 ||
          StringNCmp (str, "ID ", 3) == 0 ||
          StringStr (str, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        SeqIdFree (sip);
        if (p_line != NULL) {
          *p_line = lin_num;
        }
        if (free_str) {
          str = MemFree (str);
        }
        return sap;
      } else if (StringNCmp (str, "//", 2) == 0) {
        SeqIdFree (sip);
        if (p_line != NULL) {
          *p_line = lin_num;
        }
        if (free_str) {
          str = MemFree (str);
        }
        return sap;
      }

      if (allowWhitesp) {
        ParseWhitespaceIntoTabs (str);
      }

      feat = NULL;
      qual = NULL;
      val = NULL;

      if (*str == '[') {
        stp = SqnTagParse (str);
        if (stp != NULL) {
          tmp = SqnTagFind (stp, "offset");
          if (tmp != NULL) {
            if (sscanf (tmp, "%ld", &num) == 1) {
              offset = (Int4) num;
            }
          }
        }
        SqnTagFree (stp);

      } else if (StringNICmp (str, "ORDER", 5) == 0) {

        if (sfp != NULL) {
          PutNullsBetween (sfp->location);
        }

      } else if (ParseFeatTableLine (str, &start, &stop, &partial5, &partial3, &ispoint,
                                     &isminus, &feat, &qual, &val, offset, lin_num)) {
        if (feat != NULL && start >= 0 && stop >= 0) {

          if (sap == NULL) {
            sap = SeqAnnotNew ();
            if (sap != NULL) {
              sap->type = 1;
              if (! HasNoText (annotname)) {
                desc = AnnotDescrNew (NULL);
                if (desc != NULL) {
                  desc->choice = Annot_descr_name;
                  desc->data.ptrvalue = StringSave (annotname);
                  sap->desc = desc;
                }
              }
            }
          }

          sfp = SeqFeatNew ();
          if (sfp != NULL && sap != NULL) {
            if (sap->data != NULL) {
              if (prev == NULL) {
                prev = sap->data;
                while (prev->next != NULL) {
                  prev = prev->next;
                }
              }
              prev->next = sfp;
              prev = sfp;
            } else {
              sap->data = (Pointer) sfp;
              prev = sfp;
            }

            if (StringCmp (feat, "gene") == 0) {

              sfp->data.choice = SEQFEAT_GENE;
              grp = GeneRefNew ();
              if (grp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) grp;
              }

            } else if (StringCmp (feat, "CDS") == 0) {

              sfp->data.choice = SEQFEAT_CDREGION;
              crp = CreateNewCdRgn (1, FALSE, 0);
              if (crp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) crp;
              }

            } else if (StringStr (feat, "RNA") != NULL) {

              sfp->data.choice = SEQFEAT_RNA;
              rrp = RnaRefNew ();
              if (rrp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) rrp;
                rnatype = 255;
                if (StringCmp (feat, "precursor_RNA") == 0) {
                  rnatype = 1;
                } else if (StringCmp (feat, "mRNA") == 0) {
                  rnatype = 2;
                } else if (StringCmp (feat, "tRNA") == 0) {
                  rnatype = 3;
                } else if (StringCmp (feat, "rRNA") == 0) {
                  rnatype = 4;
                } else if (StringCmp (feat, "snRNA") == 0) {
                  rnatype = 255;
                  rrp->ext.choice = 1;
                  rrp->ext.value.ptrvalue = StringSave ("ncRNA");
                  AddQualifierToFeatureEx (sfp, "ncRNA_class", feat, offset, lin_num);
                } else if (StringCmp (feat, "scRNA") == 0) {
                  rnatype = 255;
                  rrp->ext.choice = 1;
                  rrp->ext.value.ptrvalue = StringSave ("ncRNA");
                  AddQualifierToFeatureEx (sfp, "ncRNA_class", feat, offset, lin_num);
                } else if (StringCmp (feat, "snoRNA") == 0) {
                  rnatype = 255;
                  rrp->ext.choice = 1;
                  rrp->ext.value.ptrvalue = StringSave ("ncRNA");
                  AddQualifierToFeatureEx (sfp, "ncRNA_class", feat, offset, lin_num);
                } else if (StringCmp (feat, "misc_RNA") == 0) {
                  rnatype = 255;
                  rrp->ext.choice = 1;
                  rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
                } else if (StringCmp (feat, "ncRNA") == 0) {
                  rnatype = 255;
                  rrp->ext.choice = 1;
                  rrp->ext.value.ptrvalue = StringSave ("ncRNA");
                } else if (StringCmp (feat, "tmRNA") == 0) {
                  rnatype = 255;
                  rrp->ext.choice = 1;
                  rrp->ext.value.ptrvalue = StringSave ("tmRNA");
                }
                rrp->type = rnatype;
              }

            } else if (StringCmp (feat, "Protein") == 0) {

              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              if (prp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
              }

            } else if (StringCmp (feat, "proprotein") == 0 || StringCmp (feat, "preprotein") == 0) {

              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              if (prp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
                prp->processed = 1;
              }

            } else if (StringCmp (feat, "mat_peptide") == 0) {

              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              if (prp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
                prp->processed = 2;
              }

            } else if (StringCmp (feat, "sig_peptide") == 0) {

              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              if (prp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
                prp->processed = 3;
              }

            } else if (StringCmp (feat, "transit_peptide") == 0) {

              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              if (prp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
                prp->processed = 4;
              }

            } else if (StringCmp (feat, "source") == 0) {

              sfp->data.choice = SEQFEAT_BIOSRC;
              biop = BioSourceNew ();
              if (biop != NULL) {
                orp = OrgRefNew ();
                biop->org = orp;
                sfp->data.value.ptrvalue = (Pointer) biop;
              }

            } else if (StringCmp (feat, "Region") == 0) {

              sfp->data.choice = SEQFEAT_REGION;
              sfp->data.value.ptrvalue = StringSave ("?");

            } else if (StringCmp (feat, "Bond") == 0) {

              sfp->data.choice = SEQFEAT_BOND;
              sfp->data.value.intvalue = 255;

            } else if (StringCmp (feat, "Site") == 0) {

              sfp->data.choice = SEQFEAT_SITE;
              sfp->data.value.intvalue = 255;

            } else if (StringICmp (feat, "REFERENCE") == 0 || StringICmp (feat, "CITSUB") == 0) {

              sfp->data.choice = SEQFEAT_PUB;
              pdp = PubdescNew ();
              if (pdp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) pdp;
              }

            } else {
              sfp->data.choice = SEQFEAT_IMP;
              ifp = ImpFeatNew ();
              if (ifp != NULL) {
                ifp->key = StringSave (feat);
                sfp->data.value.ptrvalue = (Pointer) ifp;
              }

              idx = -1;
              for (j = 0; j < ParFlat_TOTAL_GBFEAT; j++) {
                if (StringCmp (ParFlat_GBFeat [j].key, feat) == 0) {
                  idx = j;
                }
              }
              if (idx == -1) {
                ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_UnknownImpFeatKey, "Unknown feature %s", feat);
                if (ifp != NULL) {
                  ifp->key = MemFree (ifp->key);
                  ifp->key = StringSave ("misc_feature");
                  StringCpy (buf, "UNKNOWN FEATURE KEY ");
                  if (StringLen (feat) < 100) {
                    StringCat (buf, feat);
                  }
                  AddQualifierToFeatureEx (sfp, "note", buf, offset, lin_num);
                }
              }
            }

            if (ispoint) {
              spp = SeqPntNew ();
              if (spp != NULL) {
                spp->point = start;
                if (isminus) {
                  spp->strand = Seq_strand_minus;
                }
                spp->id = SeqIdDup (sip);
                fuzz = IntFuzzNew ();
                if (fuzz != NULL) {
                  fuzz->choice = 4;
                  fuzz->a = 3;
                  spp->fuzz = fuzz;
                }
                slp = ValNodeNew (NULL);
                if (slp != NULL) {
                  slp->choice = SEQLOC_PNT;
                  slp->data.ptrvalue = (Pointer) spp;
                  sfp->location = slp;
                }
              }
            } else {
              sfp->location = AddIntervalToLocation (NULL, sip, start, stop, partial5, partial3);
            }

            if (partial5 || partial3) {
              sfp->partial = TRUE;
            }
          }

          inquals = FALSE;

        } else if (start >= 0 && stop >= 0 && feat == NULL && qual == NULL && val == NULL && sfp != NULL) {

          if (inquals) {

            ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_ImpFeatBadLoc, "Unexpected intervals after qualifiers (start %ld, stop %ld)", (long) start, (long) stop);

          } else {

            sfp->location = AddIntervalToLocation (sfp->location, sip, start, stop, partial5, partial3);

            if (partial5 || partial3) {
              sfp->partial = TRUE;
            }
          }

        } else if (sfp != NULL && qual != NULL && (val != NULL || IsSingletonQual (qual))) {

          if (StringICmp (qual, "order") == 0) {
            if (! LocationHasNullsBetween (sfp->location)) {
              slp = SeqLocFindNext (sfp->location, NULL);
              if (slp != NULL) {
                vslp = ValNodeNew (NULL);
                if (vslp != NULL) {
                  vslp->choice = SEQLOC_NULL;
                  vslp->next = slp->next;
                  slp->next = vslp;
                }
              }
            }
            NormalizeNullsBetween (sfp->location);
          } else {
            if (StringICmp (qual, "note") == 0) {
              isnote = TRUE;
            }
            AddQualifierToFeatureEx (sfp, qual, val, offset, lin_num);
          }

          inquals = TRUE;

        } else if (sfp != NULL && qual != NULL && val == NULL) {

          label = (CharPtr) FeatDefTypeLabel (sfp);
          if (label == NULL) {
            label = "?";
          }
          loc = SeqLocPrint (sfp->location);
          if (loc == NULL) {
            loc = StringSave ("?");
          }
          if (lin_num > 0) {
            ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_WrongQualOnImpFeat, "Qualifier '%s' has no value on %s feature at %s, relative line %ld", qual, label, loc, (long) lin_num);
          } else {
            ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_WrongQualOnImpFeat, "Qualifier '%s' has no value on %s feature at %s", qual, label, loc);
          }
          MemFree (loc);

        } else if (feat != NULL) {

          ErrPostEx (SEV_ERROR, ERR_SEQ_FEAT_ImpFeatBadLoc, "Bad location on feature %s (start %ld, stop %ld)", feat, (long) start, (long) stop);
        }
      } else {
        Message (MSG_POSTERR, "Unrecognized line in feature table: %s", str);
      }

      /* ParseFeatTableLine copies these three strings, so free here */

      feat = MemFree (feat);
      qual = MemFree (qual);
      val = MemFree (val);

    }

#if 0
    /* commented out - always read in entire line now */
    /* if humongously long line /note, now extends by concatenation */

    while (nonewline && str != NULL) {
      str = FileCacheReadLine (fcp, line, sizeof (line), &nonewline);
      lin_num++;
      if (isnote && sfp != NULL && StringDoesHaveText (str)) {
        if (sfp->comment == NULL) {
          sfp->comment = StringSave (val);
        } else {
          len = StringLen (sfp->comment) + StringLen (str) + 5;
          tmp = MemNew (sizeof (Char) * len);
          StringCpy (tmp, sfp->comment);
          if (endsinspace) {
            StringCat (tmp, " ");
            endsinspace = FALSE;
          }
          StringCat (tmp, str);
          sfp->comment = MemFree (sfp->comment);
          sfp->comment = tmp;
        }
      }
    }
#endif

    pos = FileCacheTell (fcp);
    if (free_str) {
      str = MemFree (str);
      free_str = FALSE;
    }

    str = FileCacheReadLine (fcp, line, sizeof (line), &nonewline);
    if (nonewline) {
      str = ReadTheRestOfTheLine (fcp, line);
      if (StringDoesHaveText (str)) {
        free_str = TRUE;
      } else {
        str = MemFree (str);
      }
    } else {
      free_str = FALSE;
    }

    lin_num++;
  }

  if (free_str) {
    str = MemFree (str);
  }

  SeqIdFree (sip);
  if (p_line != NULL) {
    *p_line = lin_num;
  }
  return sap;
}

static SeqAnnotPtr ReadFeatureTable (FileCachePtr fcp, CharPtr seqid, CharPtr annotname)
{
  return ReadFeatureTableEx (fcp, seqid, annotname, NULL);
}

/* ReadVecScreenTable reads lines of vector screen output into a Seq-annot. */

static SeqAnnotPtr ReadVecScreenTable (FileCachePtr fcp, CharPtr seqid, CharPtr annotname)

{
  Char            ch;
  CharPtr         database = NULL;
  Char            date [32];
  DatePtr         dp;
  AnnotDescrPtr   desc;
  GeneRefPtr      grp;
  ImpFeatPtr      ifp;
  Char            line [1023];
  Char            matchtype [64];
  Char            note [128];
  Int4            pos;
  SeqFeatPtr      prev;
  CharPtr         ptr;
  SeqAnnotPtr     sap = NULL;
  CharPtr         screen = NULL;
  SeqFeatPtr      sfp = NULL;
  SeqIdPtr        sip;
  long int        start;
  long int        stop;
  CharPtr         str;
  SeqFeatXrefPtr  xref;

  if (fcp == NULL || seqid == NULL) return NULL;
  sip = SeqIdFindBest (MakeSeqID (seqid), 0);
  if (sip == NULL) return NULL;
  matchtype [0] = '\0';

  date [0] = '\0';
  dp = DateCurr ();
  DatePrint (dp, date);
  DateFree (dp);

  ptr = StringStr (annotname, "Database:");
  if (ptr != NULL) {
    ptr += 9;
    ch = *ptr;
    while (ch == ' ') {
      ptr++;
      ch = *ptr;
    }
    database = ptr;
  }

  ptr = StringStr (annotname, "Screen:");
  if (ptr != NULL) {
    ptr += 7;
    ch = *ptr;
    while (ch == ' ') {
      ptr++;
      ch = *ptr;
    }
    screen = ptr;
    while (ch != '\0' && ch != ' ') {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';
  }

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        SeqIdFree (sip);
        return sap;
      } else if (StringNCmp (line, "//", 2) == 0) {
        SeqIdFree (sip);
        return sap;
      }

      if (sscanf (line, "%ld\t%ld", &start, &stop) == 2) {
        start--;
        stop--;

        if (start >= 0 && stop >= 0) {
          if (! HasNoText (matchtype)) {

            if (sap == NULL) {
              sap = SeqAnnotNew ();
              if (sap != NULL) {
                sap->type = 1;
                if (! HasNoText (annotname)) {
                  desc = AnnotDescrNew (NULL);
                  if (desc != NULL) {
                    desc->choice = Annot_descr_name;
                    desc->data.ptrvalue = StringSave ("VecScreen");
                    sap->desc = desc;
                  }
                }
              }
            }

            if (sfp == NULL) {
              sfp = SeqFeatNew ();
              if (sfp != NULL) {

                /* make misc_feature for now */

                sfp->data.choice = SEQFEAT_IMP;
                ifp = ImpFeatNew ();
                if (ifp != NULL) {
                  ifp->key = StringSave ("misc_feature");
                }
                AddQualifierToFeature (sfp, "standard_name", "Vector Contamination");
                AddQualifierToFeature (sfp, "phenotype", matchtype);

                if ((! StringHasNoText (database)) && (! StringHasNoText (screen))) {
                  sprintf (note, "Screened against %s using %s on %s", database, screen, date);
                  sfp->comment = StringSave (note);
                }

                /* suppress /gene */

                grp = GeneRefNew ();
                if (grp != NULL) {
                  xref = SeqFeatXrefNew ();
                  sfp->xref = xref;
                  if (xref != NULL) {
                    xref->data.choice = SEQFEAT_GENE;
                    xref->data.value.ptrvalue = (Pointer) grp;
                  }
                }

                sfp->data.value.ptrvalue = (Pointer) ifp;

                if (sap != NULL) {
                  if (sap->data != NULL) {
                    prev = sap->data;
                    while (prev->next != NULL) {
                      prev = prev->next;
                    }
                    prev->next = sfp;
                  } else {
                    sap->data = (Pointer) sfp;
                  }
                }

                sfp->location = AddIntervalToLocation (NULL, sip, (Int4) start, (Int4) stop, FALSE, FALSE);
              }

            } else {

              sfp->location = AddIntervalToLocation (sfp->location, sip, (Int4) start, (Int4) stop, FALSE, FALSE);

            }
          }
        }

      } else {
        StringNCpy_0 (matchtype, line, sizeof (matchtype));
        sfp = NULL;
        if (StringCmp (matchtype, "No hits found") == 0) {
          sprintf (note, "No vector hits found for %s", seqid);
          Message (MSG_POST, "%s\n", note);
        }
      }

    }

    pos = FileCacheTell (fcp);
    str = FileCacheGetString (fcp, line, sizeof (line));
  }

  SeqIdFree (sip);
  return sap;
}

/* ReadRestrictionSiteTable reads lines of restriction enzyme names or cut sites into a Seq-annot. */

static SeqLocPtr AddPointToLocation (SeqLocPtr loc, SeqIdPtr sip, Int4 pt)

{
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;

  if (sip == NULL) return NULL;

  if (loc == NULL) {
    pspp = PackSeqPntNew ();
    pspp->id = SeqIdDup (sip);
    slp = ValNodeNew (NULL);
    slp->choice = SEQLOC_PACKED_PNT;
    slp->data.ptrvalue = (Pointer) pspp;
    loc = slp;
  }

  if (loc != NULL && loc->choice == SEQLOC_PACKED_PNT) {
    pspp = (PackSeqPntPtr) loc->data.ptrvalue;
    if (pspp != NULL) {
      PackSeqPntPut (pspp, pt);
    }
  }

  return loc;
}

static SeqAnnotPtr ReadRestrictionSiteTable (FileCachePtr fcp, CharPtr seqid, CharPtr annotname)

{
  DbtagPtr       dbt;
  AnnotDescrPtr  desc;
  Char           line [1023];
  Char           name [64];
  ObjectIdPtr    oip;
  Int4           pos;
  SeqFeatPtr     prev;
  Int4           pt;
  RsiteRefPtr    rrp;
  SeqAnnotPtr    sap = NULL;
  SeqFeatPtr     sfp = NULL;
  SeqIdPtr       sip;
  CharPtr        str;
  long int       val;

  if (fcp == NULL || seqid == NULL) return NULL;
  sip = SeqIdFindBest (MakeSeqID (seqid), 0);
  if (sip == NULL) return NULL;
  name [0] = '\0';

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        SeqIdFree (sip);
        return sap;
      } else if (StringNCmp (line, "//", 2) == 0) {
        SeqIdFree (sip);
        return sap;
      }

      if (sscanf (line, "%ld", &val) == 1) {
        pt = (Int4) val;

        if (! HasNoText (name)) {

          if (sap == NULL) {
            sap = SeqAnnotNew ();
            if (sap != NULL) {
              sap->type = 1;
              if (! HasNoText (annotname)) {
                desc = AnnotDescrNew (NULL);
                if (desc != NULL) {
                  desc->choice = Annot_descr_name;
                  desc->data.ptrvalue = StringSave (annotname);
                  sap->desc = desc;
                }
              }
            }
          }

          if (sfp == NULL) {
            sfp = SeqFeatNew ();
            if (sfp != NULL) {
              sfp->data.choice = SEQFEAT_RSITE;
              dbt = DbtagNew ();
              if (dbt != NULL) {
                dbt->db = StringSave ("REBASE");
                oip = ObjectIdNew ();
                if (oip != NULL) {
                  oip->str = StringSave (name);
                }
                dbt->tag = oip;
              }
              rrp = ValNodeNew (NULL);
              if (rrp != NULL) {
                rrp->choice = 2;
                rrp->data.ptrvalue = dbt;
              }
              sfp->data.value.ptrvalue = (Pointer) rrp;

              if (sap != NULL) {
                if (sap->data != NULL) {
                  prev = sap->data;
                  while (prev->next != NULL) {
                    prev = prev->next;
                  }
                  prev->next = sfp;
                } else {
                  sap->data = (Pointer) sfp;
                }
              }
            }
          }

          if (sfp != NULL) {
            sfp->location = AddPointToLocation (sfp->location, sip, pt);
          }

        }

      } else {
        StringNCpy_0 (name, line, sizeof (name));
        sfp = NULL;
      }

    }

    pos = FileCacheTell (fcp);
    str = FileCacheGetString (fcp, line, sizeof (line));
  }

  SeqIdFree (sip);
  return sap;
}

/* ReadMessageStrings allows retired services to announce replacement URLs. */

static void ReadMessageStrings (FileCachePtr fcp)

{
  Boolean     done = FALSE;
  ValNodePtr  head = NULL;
  size_t      len;
  Char        line [1023];
  Int4        pos;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  if (fcp == NULL) return;

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  while (str != NULL && (! done)) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        done = TRUE;
      } else if (StringNCmp (line, "//", 2) == 0) {
        done = TRUE;
      }

      if (! done) {
        ValNodeCopyStr (&head, 0, line);
      }
      /* Message (MSG_POST, "%s\n", line); */
    }

    if (! done) {
      pos = FileCacheTell (fcp);
      str = FileCacheGetString (fcp, line, sizeof (line));
    }
  }

  for (vnp = head, len = 0; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      len += StringLen (str) + 1;
    }
  }
  if (len > 0) {
    ptr = MemNew (sizeof (Char) * (len + 2));
    if (ptr != NULL) {
      for (vnp = head, tmp = NULL; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (str != NULL) {
          if (tmp == NULL) {
            tmp = ptr;
          } else {
            tmp = StringMove (tmp, "\n");
          }
          tmp = StringMove (tmp, str);
        }
      }
      Message (MSG_POST, "%s\n", ptr);
      MemFree (ptr);
    }
  }

  ValNodeFreeData (head);
}

/* ReadUidList reads lines of uids (or accessions) into a byte store. */

static ByteStorePtr ReadUidList (FileCachePtr fcp, Boolean nucdb, Boolean lastResortSeqIDs)

{
  Boolean       allDigits;
  Boolean       abort = FALSE;
  ByteStorePtr  bs;
  Char          ch;
  Char          line [1023];
  Int4          pos;
  CharPtr       ptr;
  CharPtr       str;
  TextSeqId     tsid;
  Int4          uid;
  long int      val;
  ValNode       vn;

  if (fcp == NULL) return NULL;
  bs = BSNew (128);
  if (bs == NULL) return NULL;

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        FileCacheSeek (fcp, pos);
        if (abort) {
          bs = BSFree (bs);
        }
        return bs;
      } else if (StringNCmp (line, "//", 2) == 0) {
        if (abort) {
          bs = BSFree (bs);
        }
        return bs;
      }

      allDigits = TRUE;
      ptr = line;
      ch = *ptr;
      while (ch != '\0' && allDigits) {
        if (! IS_DIGIT (ch)) {
          allDigits = FALSE;
        }
        ptr++;
        ch = *ptr;
      }
      if (allDigits && sscanf (line, "%ld", &val) == 1) {
        uid = (Int4) val;
        BSWrite (bs, &uid, sizeof (Int4));
      } else if (nucdb) {
        tsid.name = NULL;
        tsid.accession = line;
        tsid.release = NULL;
        tsid.version = INT2_MIN;
        vn.choice = (Uint1) SEQID_GENBANK;
        vn.data.ptrvalue = (Pointer) (&tsid);
        uid = GetGIForSeqId (&vn);
        if (uid > 0) {
          BSWrite (bs, &uid, sizeof (Int4));
        } else if (lastResortSeqIDs) {
          abort = TRUE;
        }
      }

    }

    pos = FileCacheTell (fcp);
    str = FileCacheGetString (fcp, line, sizeof (line));
  }

  if (abort) {
    bs = BSFree (bs);
  }
  return bs;
}


static Boolean DoesBioseqAccessionMatchList (BioseqPtr bsp, ValNodePtr accn_list)
{
  ValNodePtr vnp, vnp_m;
  ValNodePtr match_list = NULL;
  SeqIdPtr   sip, sip_next;
  CharPtr    id, cp;
  Boolean    found_match = FALSE;
  DbtagPtr   dbtag;

  if (bsp == NULL) {
    return FALSE;
  }

  /* note - in match_list, 1 indicates that memory needs to be freed, 0 not */
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    sip_next = sip->next;
    sip->next = NULL;
    id = SeqIdWholeLabel (sip, PRINTID_FASTA_LONG);
    sip->next = sip_next;
    if (id != NULL) {
      /* remove terminating pipe character */
      if (id[StringLen(id) - 1] == '|') 
      {
        id[StringLen(id) - 1] = 0;
      }
      ValNodeAddPointer(&match_list, 1, id);

      /* remove leading pipe identifier */
      cp = StringChr (id, '|');
      if (cp != NULL)
      {
        cp = cp + 1;
        ValNodeAddPointer (&match_list, 0, cp);
      } else {
        cp = id;
      }

      /* try ID without version */
      id = StringSave (cp);
      cp = StringChr (id, '.');
      if (cp != NULL) 
      {
        *cp = 0;
        ValNodeAddPointer (&match_list, 1, id);
      } else {
        id = MemFree (id);
      }

      /* just bankit number */
      if (sip->choice == SEQID_GENERAL 
          && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL) {
        if (StringCmp (dbtag->db, "BankIt") == 0) {
          if (dbtag->tag->str != NULL) {
            ValNodeAddPointer (&match_list, 0, dbtag->tag->str);
          }
        } else if (StringCmp (dbtag->db, "NCBIFILE") == 0 && dbtag->tag != NULL) {
          ValNodeAddPointer (&match_list, 0, dbtag->tag->str);
          if ((cp = StringRChr (dbtag->tag->str, '/')) != NULL) {
            ValNodeAddPointer (&match_list, 0, cp + 1);
          }
        }
      }
    }
  }

  for (vnp = accn_list; vnp != NULL && !found_match; vnp = vnp->next) {
    for (vnp_m = match_list; vnp_m != NULL && !found_match; vnp_m = vnp_m->next) {
      if (StringICmp (vnp->data.ptrvalue, vnp_m->data.ptrvalue) == 0) {
        found_match = TRUE;
      }
    }
  }

  /* special free for match_list */
  vnp = ValNodeExtractList (&match_list, 1);
  vnp = ValNodeFreeData (vnp);
  match_list = ValNodeFree (match_list);

  return found_match;
}


static Boolean DoesBioseqSetAccessionMatchList (BioseqSetPtr bssp, ValNodePtr accn_list)
{
  BioseqPtr bsp;
  Boolean   rval = FALSE;

  if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot
    && bssp->seq_set != NULL && IS_Bioseq (bssp->seq_set)) {
    bsp = bssp->seq_set->data.ptrvalue;
    rval = DoesBioseqAccessionMatchList(bsp, accn_list);
  }
  return rval;
}


static Boolean DoesSeqEntryAccessionMatchList (SeqEntryPtr sep, ValNodePtr accn_list)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  Boolean      rval = FALSE;

  if (sep == NULL) {
    return FALSE;
  }

  if (IS_Bioseq (sep)) {
    bsp = sep->data.ptrvalue;
    rval = DoesBioseqAccessionMatchList (bsp, accn_list);
  } else {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    rval = DoesBioseqSetAccessionMatchList (bssp, accn_list);
  }

  return rval;
}


static Boolean s_IsDelimiter (Char ch)
{
  if (isspace (ch) || ch == ',' || ch == ';') {
    return TRUE;
  } else {
    return FALSE;
  }
}


static ValNodePtr ListFromString (CharPtr accn_list)
{
  CharPtr start, stop, id;
  ValNodePtr list = NULL;
  Int4 len;

  start = accn_list;
  while (*start != 0) {
    while (*start != 0 && s_IsDelimiter(*start)) {
      start++;
    }
    if (*start != 0) {
      stop = start + 1;
      while (*stop != 0 && !s_IsDelimiter (*stop)) {
        stop++;
      }
      len = stop - start;
      id = (CharPtr) MemNew (sizeof (Char) * (len + 1));
      StringNCpy (id, start, len);
      id[len] = 0;
      ValNodeAddPointer (&list, 0, id);
      start = stop;
    }
  }

  return list;
}

typedef struct setatp {
  AsnModulePtr amp;
  AsnTypePtr atp_class;
  AsnTypePtr atp_seqset;
  AsnTypePtr atp_se;
  AsnTypePtr atp_descr;
  AsnTypePtr atp_descr_e;
  AsnTypePtr atp_set_desc;
  AsnTypePtr atp_bioseq_desc;
  AsnTypePtr atp_desc;
  AsnTypePtr atp_annot;
  AsnTypePtr atp_bioseq_annot;
  AsnTypePtr atp_annot_e;
  AsnTypePtr atp_bioseq_annot_e;
  AsnTypePtr atp_id;
  AsnTypePtr atp_coll;
  AsnTypePtr atp_date;
  AsnTypePtr atp_level;
  AsnTypePtr atp_release;
  AsnTypePtr atp_bss;
  AsnTypePtr atp_bioseq;
  AsnTypePtr atp_seqentry;
  AsnTypePtr atp_seq;
  AsnTypePtr atp_set;
  AsnTypePtr atp_seqsubmit;
  AsnTypePtr atp_sub;
  AsnTypePtr atp_seqsubmit_data;
  AsnTypePtr atp_seqsubmit_data_entries_E;
  AsnTypePtr atp_seqsubmit_data_entries;
  AsnTypePtr atp_seqsubmit_data_entries_set;
  AsnTypePtr atp_bioseq_id_E;
  AsnTypePtr atp_seqdesc_pub;
} SetAtpData, PNTR SetAtpPtr;


static SetAtpPtr GetSetAtp (void)
{
  AsnModulePtr amp;
  AsnTypePtr atp_class;
  AsnTypePtr atp_seqset;
  AsnTypePtr atp_se;
  AsnTypePtr atp_descr;
  AsnTypePtr atp_descr_e;
  AsnTypePtr atp_set_desc;
  AsnTypePtr atp_bioseq_desc;
  AsnTypePtr atp_desc;
  AsnTypePtr atp_annot;
  AsnTypePtr atp_bioseq_annot;
  AsnTypePtr atp_annot_e;
  AsnTypePtr atp_bioseq_annot_e;
  AsnTypePtr atp_id;
  AsnTypePtr atp_coll;
  AsnTypePtr atp_date;
  AsnTypePtr atp_level;
  AsnTypePtr atp_release;
  AsnTypePtr atp_bss;
  AsnTypePtr atp_seqentry;
  AsnTypePtr atp_seq;
  AsnTypePtr atp_set;
  AsnTypePtr atp_seqsubmit;
  AsnTypePtr atp_sub;
  AsnTypePtr atp_seqsubmit_data;
  AsnTypePtr atp_seqsubmit_data_entries_E;
  AsnTypePtr atp_seqsubmit_data_entries;
  AsnTypePtr atp_seqsubmit_data_entries_set;
  AsnTypePtr atp_bioseq;
  AsnTypePtr atp_bioseq_id_E;
  AsnTypePtr atp_seqdesc_pub;
  SetAtpPtr  sp;

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_POSTERR, "Unable to load AsnAllModPtr");
    return NULL;
  }

  atp_seqset = AsnFind ("Bioseq-set.seq-set");
  if (atp_seqset == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set");
    return NULL;
  }

  atp_se = AsnFind ("Bioseq-set.seq-set.E");
  if (atp_se == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
    return NULL;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set");
    return NULL;
  }

  atp_bioseq = AsnFind ("Bioseq");
  if (atp_bioseq == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq");
    return NULL;
  }

  atp_class = AsnFind ("Bioseq-set.class");
  if (atp_class == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.class");
    return NULL;
  }

  atp_descr = AsnFind ("Seq-descr");
  if (atp_descr == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-descr");
    return NULL;
  }

  atp_descr_e = AsnFind ("Seq-descr.E");
  if (atp_descr_e == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-descr.E");
    return NULL;
  }

  atp_set_desc = AsnFind ("Bioseq-set.descr");
  if (atp_set_desc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.descr");
    return NULL;
  }

  atp_bioseq_desc = AsnFind ("Bioseq.descr");
  if (atp_bioseq_desc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq.descr");
    return NULL;
  }

  atp_desc = AsnFind ("Seqdesc");
  if (atp_desc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seqdesc");
    return NULL;
  }

  atp_annot = AsnFind ("Bioseq-set.annot");
  if (atp_annot == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.annot");
    return NULL;
  }

  atp_bioseq_annot = AsnFind ("Bioseq.annot");
  if (atp_bioseq_annot == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq.annot");
    return NULL;
  }

  atp_annot_e = AsnFind ("Bioseq-set.annot.E");
  if (atp_annot_e == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.annot.E");
    return NULL;
  }

  atp_bioseq_annot_e = AsnFind ("Bioseq.annot.E");
  if (atp_bioseq_annot_e == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq.annot.E");
    return NULL;
  }

  atp_id = AsnFind ("Bioseq-set.id");
  if (atp_id == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.id");
    return NULL;
  }

  atp_coll = AsnFind ("Bioseq-set.coll");
  if (atp_coll == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.coll");
    return NULL;
  }
  atp_date = AsnFind ("Bioseq-set.date");
  if (atp_date == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.date");
    return NULL;
  }
  atp_level = AsnFind ("Bioseq-set.level");
  if (atp_level == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.level");
    return NULL;
  }
  atp_release = AsnFind ("Bioseq-set.release");
  if (atp_release == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.release");
    return NULL;
  }

  atp_seqentry = AsnFind ("Seq-entry");
  if (atp_seqentry == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry");
    return NULL;
  }

  atp_seq = AsnFind ("Seq-entry.seq");
  if (atp_seq == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry.seq");
    return NULL;
  }

  atp_set = AsnFind ("Seq-entry.set");
  if (atp_set == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-entry.set");
    return NULL;
  }

  atp_seqsubmit = AsnFind ("Seq-submit");
  if (atp_seqsubmit == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    return NULL;
  }

  atp_sub = AsnFind ("Seq-submit.sub");
  if (atp_sub == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    return NULL;
  }

  atp_seqsubmit_data = AsnFind ("Seq-submit.data");
  if (atp_seqsubmit_data == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data");
    return NULL;
  }

  atp_seqsubmit_data_entries_E = AsnFind ("Seq-submit.data.entrys.E");
  if (atp_seqsubmit_data_entries_E == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys.E");
    return NULL;
  }

  atp_seqsubmit_data_entries = AsnFind ("Seq-submit.data.entrys");
  if (atp_seqsubmit_data_entries == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys");
    return NULL;
  }

  atp_seqsubmit_data_entries_set = AsnFind ("Seq-submit.data.entrys.E.set");
  if (atp_seqsubmit_data_entries == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys.E.set");
    return NULL;
  }

  atp_bioseq_id_E = AsnFind ("Bioseq.id.E");
  if (atp_bioseq_id_E == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq.id.E");
    return NULL;
  }

  atp_seqdesc_pub = AsnFind ("Seqdesc.pub");
  if (atp_seqdesc_pub == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seqdesc.pub");
    return NULL;
  }

  sp = (SetAtpPtr) MemNew (sizeof(SetAtpData));
  sp->amp = amp;
  sp->atp_class = atp_class;
  sp->atp_seqset = atp_seqset;
  sp->atp_se = atp_se;
  sp->atp_descr = atp_descr;
  sp->atp_descr_e = atp_descr_e;
  sp->atp_set_desc = atp_set_desc;
  sp->atp_bioseq_desc = atp_bioseq_desc;
  sp->atp_desc = atp_desc;
  sp->atp_annot = atp_annot;
  sp->atp_bioseq_annot = atp_bioseq_annot;
  sp->atp_annot_e = atp_annot_e;
  sp->atp_bioseq_annot_e = atp_bioseq_annot_e;
  sp->atp_id = atp_id;
  sp->atp_coll = atp_coll;
  sp->atp_date = atp_date;
  sp->atp_level = atp_level;
  sp->atp_release = atp_release;
  sp->atp_bss = atp_bss;
  sp->atp_bioseq = atp_bioseq;
  sp->atp_seqentry = atp_seqentry;
  sp->atp_seq = atp_seq;
  sp->atp_set = atp_set;
  sp->atp_seqsubmit = atp_seqsubmit;
  sp->atp_sub = atp_sub;
  sp->atp_seqsubmit_data = atp_seqsubmit_data;
  sp->atp_seqsubmit_data_entries_E = atp_seqsubmit_data_entries_E;
  sp->atp_seqsubmit_data_entries = atp_seqsubmit_data_entries;
  sp->atp_seqsubmit_data_entries_set = atp_seqsubmit_data_entries_set;
  sp->atp_bioseq_id_E = atp_bioseq_id_E;
  sp->atp_seqdesc_pub = atp_seqdesc_pub;

  return sp;
}


static BioseqSetPtr BioseqSetPartialRead (AsnIoPtr aip, AsnTypePtr PNTR orig, SetAtpPtr sp)
{
  DataVal av;
  AsnTypePtr atp, oldatp;
  BioseqSetPtr bsp=NULL;
  SeqEntryPtr curr, next;


	if (aip == NULL)
		return bsp;

	if (orig == NULL || *orig == NULL)           /* BioseqSet ::= (self contained) */
		atp = AsnReadId(aip, sp->amp, sp->atp_bss);
	else
		atp = AsnLinkType(*orig, sp->atp_bss);    /* link in local tree */

  oldatp = atp;
  if (atp == NULL) {
    if (orig != NULL) {
      *orig = atp;
    }
    return bsp;
  }

	bsp = BioseqSetNew();
	if (bsp == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    curr = NULL;

    while ((atp = AsnReadId(aip, sp->amp, atp)) != oldatp)
    {
		  if (atp == NULL) goto erret;
      if (atp == sp->atp_id)
		  {
        bsp->id = ObjectIdAsnRead(aip, atp);
			  if (bsp->id == NULL) goto erret;
		  }
      else if (atp == sp->atp_coll)
		  {
        bsp->coll = DbtagAsnRead(aip, atp);
			  if (bsp->coll == NULL) goto erret;
		  }
      else if (atp == sp->atp_date)
		  {
        bsp->date = DateAsnRead(aip, atp);
			  if (bsp->date == NULL) goto erret;
		  }
      else if (atp == sp->atp_set_desc)
		  {
        bsp->descr = SeqDescrAsnRead(aip, atp);
			  if (bsp->descr == NULL) goto erret;
		  }
      else if (atp == sp->atp_se)
      {
		  	if ((next = SeqEntryAsnRead(aip, atp)) != NULL)
			  {
				  if (IS_Bioseq(next))
					  SeqMgrConnect(SM_BIOSEQ, next->data.ptrvalue,
						              SM_BIOSEQSET, (Pointer) bsp);
				  else
					  SeqMgrConnect(SM_BIOSEQSET, next->data.ptrvalue,
						              SM_BIOSEQSET, (Pointer) bsp);


          if (curr == NULL)
  			    bsp->seq_set = next;
      		else
           curr->next = next;
          curr = next;
			  }
      }
      else if (atp == sp->atp_annot)
      {
        bsp->annot = SeqAnnotSetAsnRead(aip, atp, sp->atp_annot_e);
				if (bsp->annot == NULL) goto erret;
      }
      else
      {
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* takes care of everything else */
        if (atp == sp->atp_level)
          bsp->level = (Int2)av.intvalue;
        else if (atp == sp->atp_class)
			  {
          bsp->_class = (Uint1)av.intvalue;
          if (bsp->_class != BioseqseqSet_class_nuc_prot) {
            if (orig != NULL) {
              *orig = atp;
            }
            bsp->descr = NULL;
            bsp = BioseqSetFree (bsp);
            return NULL;
          }

			  }
        else if (atp == sp->atp_release)
          bsp->release = (CharPtr)av.ptrvalue;
      }
    }
  if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end BioseqSet */



ret:
  AsnUnlinkType(*orig);     /*  unlink local tree */
  return bsp;
erret:
  aip->io_failure = TRUE;
  bsp = BioseqSetFree(bsp);
  goto ret;
}


NLM_EXTERN SeqEntryPtr ReadFilteredAsn (FILE *fp, Boolean is_binary, CharPtr accn_list, Uint2Ptr entityIDptr)
{
  AsnIoPtr       aip;
  SetAtpPtr      sp;
  AsnTypePtr     atp, atp_ssp;
  SeqEntryPtr    sep = NULL, inner_sep, last_sep = NULL;
  BioseqSetPtr   bssp = NULL;
  BioseqSetPtr   nuc_set;
  BioseqPtr      bsp;
  ValNodePtr     id_match_list;
  SeqDescrPtr    sdp = NULL;
  Uint1          holding_set_class = BioseqseqSet_class_genbank;

  if (fp == NULL) return NULL;

  sp = GetSetAtp ();
  if (sp == NULL) {
    return NULL;
  }

  atp_ssp = AsnFind ("Seq-submit");
  if (atp_ssp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    sp = MemFree (sp);
    return NULL;
  }

  aip = AsnIoNew (is_binary ? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file");
    sp = MemFree (sp);
    return NULL;
  }

  if ((atp = AsnReadId (aip, sp->amp, atp_ssp)) != NULL) {
    AsnReadVal (aip, atp, NULL);
    atp = AsnReadId (aip, sp->amp, atp);
  } else {
    AsnIoFree (aip, FALSE);
    rewind (fp);
    aip = AsnIoNew (is_binary ? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
    atp = AsnReadId (aip, sp->amp, sp->atp_seqentry);
    if (atp == NULL) {
      AsnIoFree (aip, FALSE);
      rewind (fp);
      aip = AsnIoNew (is_binary ? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
      atp = AsnReadId (aip, sp->amp, sp->atp_bss);
    } else {
      AsnReadVal (aip, atp, NULL);
      atp = AsnReadId(aip, sp->amp, atp);
    }
  }
  if (atp == NULL) {
    AsnIoFree (aip, FALSE);
    sp = MemFree (sp);
    return NULL;
  }

  id_match_list = ListFromString (accn_list);

  bssp = BioseqSetNew ();
  bssp->_class = holding_set_class;
  sep = SeqEntryNew ();
  sep->choice = 2;
  sep->data.ptrvalue = bssp;

  SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);


  while (! aip->io_failure && atp != NULL) {
    if (atp == sp->atp_set) {
      nuc_set = BioseqSetPartialRead (aip, &atp, sp);
      if (nuc_set != NULL) {
        if (DoesBioseqSetAccessionMatchList (nuc_set, id_match_list)) {
          inner_sep = SeqEntryNew();
          inner_sep->choice = 2;
          inner_sep->data.ptrvalue = nuc_set;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) nuc_set, inner_sep);

          if (last_sep == NULL) {
            bssp->seq_set = inner_sep;
          } else {
            last_sep->next = inner_sep;
          }
          last_sep = inner_sep;
          SeqMgrLinkSeqEntry (inner_sep, OBJ_BIOSEQ, bssp);
        } else {
          nuc_set = BioseqSetFree (nuc_set);
        }
      }
    } else if (atp == sp->atp_seq) {
      bsp = BioseqAsnRead (aip, atp);
      if (DoesBioseqAccessionMatchList (bsp, id_match_list)) {
        inner_sep = SeqEntryNew();
        inner_sep->choice = 1;
        inner_sep->data.ptrvalue = bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, inner_sep);

        if (last_sep == NULL) {
          bssp->seq_set = inner_sep;
        } else {
          last_sep->next = inner_sep;
        }
        last_sep = inner_sep;
        SeqMgrLinkSeqEntry (inner_sep, OBJ_BIOSEQ, bssp);
      } else {
        bsp = BioseqFree (bsp);
      }
    } else if (atp == sp->atp_set_desc) {
      sdp = SeqDescrAsnRead (aip, atp);      
      ValNodeLink (&(bssp->descr), ValNodeExtractList (&sdp, Seq_descr_pub));
      sdp = SeqDescrFree (sdp);
    } else {
      AsnReadVal (aip, atp, NULL);
    }
    atp = AsnReadId (aip, sp->amp, atp);
  }

  id_match_list = ValNodeFreeData (id_match_list);

  AsnIoFree (aip, FALSE);
  if (bssp != NULL) {
    *entityIDptr = ObjMgrRegister (OBJ_SEQENTRY, sep);
  }
  sp = MemFree (sp);
  return sep;
}


static Boolean IdListsCoincide (SeqIdPtr list1, SeqIdPtr list2)
{
  SeqIdPtr sip1, sip2;
  Boolean  found_match = FALSE;
  Boolean  found_mismatch = FALSE;
  Uint1    comp;

  for (sip1 = list1; sip1 != NULL && !found_mismatch; sip1 = sip1->next) {
    for (sip2 = list2; sip2 != NULL && !found_mismatch; sip2 = sip2->next) {
      comp = SeqIdComp (sip1, sip2);
      if (comp == SIC_YES) {
        found_match = TRUE;
      } else if (comp == SIC_NO) {
        found_mismatch = TRUE;
      }
    }
  }
  if (found_match && !found_mismatch) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean DoesSeqReplaceSeq (BioseqPtr seq1, BioseqPtr seq2)
{
  if (seq1 == NULL || seq2 == NULL) {
    return FALSE;
  }
  return IdListsCoincide (seq1->id, seq2->id);
}


static Boolean DoesSetReplaceSeq (BioseqSetPtr set, BioseqPtr seq)
{
  BioseqPtr nuc;
  if (set == NULL || seq == NULL || set->_class != BioseqseqSet_class_nuc_prot
      || set->seq_set == NULL || !IS_Bioseq (set->seq_set)
      || (nuc = (BioseqPtr) set->seq_set->data.ptrvalue) == NULL) {
    return FALSE;
  }
  return DoesSeqReplaceSeq (nuc, seq);
}


static Boolean DoesSetReplaceSet (BioseqSetPtr set1, BioseqSetPtr set2)
{
  if (set1 == NULL || set2 == NULL
      || set1->_class != BioseqseqSet_class_nuc_prot
      || set2->_class != BioseqseqSet_class_nuc_prot
      || set1->seq_set == NULL
      || set2->seq_set == NULL
      || !IS_Bioseq (set1->seq_set)
      || !IS_Bioseq (set2->seq_set)) {
    return FALSE;
  } else {
    return DoesSeqReplaceSeq (set1->seq_set->data.ptrvalue, set2->seq_set->data.ptrvalue);
  }
}


static SeqEntryPtr FindReplacementSeqEntry (SeqEntryPtr edited, SeqEntryPtr orig)
{
  BioseqPtr bsp_e, bsp_o;
  BioseqSetPtr bssp;
  SeqEntryPtr  sep_replace = NULL, tmp;

  if (edited == NULL || orig == NULL) {
    return NULL;
  }

  if (IS_Bioseq (edited)) {
    bsp_e = (BioseqPtr) edited->data.ptrvalue;
    if (IS_Bioseq (orig)) {
      bsp_o = (BioseqPtr) orig->data.ptrvalue;
      if (IdListsCoincide(bsp_e->id, bsp_o->id)) {
        sep_replace = edited;
      }
    } else {
      if (DoesSetReplaceSeq (orig->data.ptrvalue, bsp_e)) {
        sep_replace = edited;
      }
    }
  } else if (IS_Bioseq_set (edited)) {
    bssp = (BioseqSetPtr) edited->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_nuc_prot) {
      if (IS_Bioseq (orig)) {
        if (DoesSetReplaceSeq (bssp, orig->data.ptrvalue)) {
          sep_replace = edited;
        }
      } else {
        if (DoesSetReplaceSet (bssp, orig->data.ptrvalue)) {
          sep_replace = edited;
        }
      }
    } else {
      for (tmp = bssp->seq_set; tmp != NULL && sep_replace == NULL; tmp = tmp->next) {
        sep_replace = FindReplacementSeqEntry (tmp, orig);
      }
    }
  }

  return sep_replace;
}


static Boolean BioseqSetWriteBefore (BioseqSetPtr bsp, AsnIoPtr aip, AsnTypePtr orig, SetAtpPtr sp)
{
	DataVal av;
	AsnTypePtr atp;
	Boolean retval = FALSE;

	if (aip == NULL)
		return FALSE;

  /* first write Seq-entry lead-in */
  if (!AsnWriteChoice(aip, sp->atp_se, (Int2)2, &av)) goto erret;

	atp = AsnLinkType(orig, sp->atp_bss);   /* link local tree */
	if (atp == NULL) return FALSE;

	if (bsp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

	if (! AsnOpenStruct(aip, atp, (Pointer)bsp)) goto erret;
    
  if (bsp->id != NULL)
	{
        if (! ObjectIdAsnWrite(bsp->id, aip, sp->atp_id)) goto erret;
	}
  if (bsp->coll != NULL)
	{
        if (! DbtagAsnWrite(bsp->coll, aip, sp->atp_coll)) goto erret;
	}
  if (bsp->level != INT2_MIN)
  {
    av.intvalue = bsp->level;
    if (! AsnWrite(aip, sp->atp_level, &av)) goto erret;
  }
  if (bsp->_class != 0)
  {
    av.intvalue = bsp->_class;
    if (! AsnWrite(aip, sp->atp_class, &av)) goto erret;
  }
  if (bsp->release != NULL)
  {
    av.ptrvalue = bsp->release;
    if (! AsnWrite(aip, sp->atp_release, &av)) goto erret;
  }
  if (bsp->date != NULL)
	{
      if (! DateAsnWrite(bsp->date, aip, sp->atp_date)) goto erret;
	}
  if (bsp->descr != NULL)              /* Seq-descr optional */
	{
    if (! SeqDescrAsnWrite(bsp->descr, aip, sp->atp_set_desc)) goto erret;
	}

  if (! AsnOpenStruct(aip, sp->atp_seqset, (Pointer)bsp->seq_set)) goto erret;

	retval = TRUE;
erret:
	return retval;

}


static Boolean BioseqSetWriteAfter (BioseqSetPtr bsp, AsnIoPtr aip, AsnTypePtr orig, SetAtpPtr sp)
{
	Boolean retval = FALSE;

	if (aip == NULL)
		return FALSE;

  if (! AsnCloseStruct(aip, sp->atp_seqset, (Pointer)bsp->seq_set)) goto erret;
    if (bsp->annot != NULL)              /* annotation optional */
	{
        if (! SeqAnnotSetAsnWrite(bsp->annot, aip, sp->atp_annot, sp->atp_annot_e)) goto erret;
	}

    if (! AsnCloseStruct(aip, orig, (Pointer)bsp)) goto erret;
	retval = TRUE;
erret:

  return retval;
}


static void SeqEntryCopyReplace (AsnIoPtr aip_in, AsnIoPtr aip_out, SeqEntryPtr edited, AsnTypePtr orig, SetAtpPtr sp);

static SeqEntryPtr BioseqSetCopyReplace (AsnIoPtr aip_in, AsnIoPtr aip_out, SeqEntryPtr edited,
                                          AsnTypePtr PNTR orig, SetAtpPtr sp)
{
  DataVal      av;
  AsnTypePtr   atp, oldatp;
  BioseqSetPtr bsp=NULL, edited_set;
  SeqEntryPtr  curr, next;
  Boolean      wrote_front = FALSE;
  SeqDescrPtr  tmp;
  SeqEntryPtr  tmp_sep, replace;
  SeqEntryPtr  sep_return = NULL;

	if (aip_in == NULL)
		return NULL;

	if (orig == NULL || *orig == NULL)           /* BioseqSet ::= (self contained) */
		atp = AsnReadId(aip_in, sp->amp, sp->atp_bss);
	else
		atp = AsnLinkType(*orig, sp->atp_bss);    /* link in local tree */

  oldatp = atp;
  if (atp == NULL) {
    if (orig != NULL) {
      *orig = atp;
    }
    return NULL;
  }

	bsp = BioseqSetNew();
	if (bsp == NULL) goto erret;

  edited_set = (BioseqSetPtr) edited->data.ptrvalue;

	if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;    /* read the start struct */
    curr = NULL;

    while ((atp = AsnReadId(aip_in, sp->amp, atp)) != oldatp)
    {
		  if (atp == NULL) goto erret;
      if (atp == sp->atp_id)
		  {
        bsp->id = ObjectIdAsnRead(aip_in, atp);
			  if (bsp->id == NULL) goto erret;
		  }
      else if (atp == sp->atp_coll)
		  {
        bsp->coll = DbtagAsnRead(aip_in, atp);
			  if (bsp->coll == NULL) goto erret;
		  }
      else if (atp == sp->atp_date)
		  {
        bsp->date = DateAsnRead(aip_in, atp);
			  if (bsp->date == NULL) goto erret;
		  }
      else if (atp == sp->atp_set_desc)
		  {
        bsp->descr = SeqDescrAsnRead(aip_in, atp);
			  if (bsp->descr == NULL) goto erret;
		  }
      else if (atp == sp->atp_se)
      {
        /* if this is a nuc prot set, read in the entries so we can see if this is
         * a candidate for replacement.
         * otherwise, write out the first part of the bioseq set here,
         * then write/replace the individual seq-entries.
         */

        if (bsp->_class == BioseqseqSet_class_nuc_prot) {
		  	  if ((next = SeqEntryAsnRead(aip_in, atp)) != NULL)
			    {
				    if (IS_Bioseq(next))
					    SeqMgrConnect(SM_BIOSEQ, next->data.ptrvalue,
						                SM_BIOSEQSET, (Pointer) bsp);
				    else
					    SeqMgrConnect(SM_BIOSEQSET, next->data.ptrvalue,
						                SM_BIOSEQSET, (Pointer) bsp);


            if (curr == NULL)
  			      bsp->seq_set = next;
      		  else
             curr->next = next;
            curr = next;
			    }
        } else {
          /* write front, loop through lower items */
          if (!wrote_front) {
            /* remove old pubs */
            tmp = ValNodeExtractList (&(bsp->descr), Seq_descr_pub);
            tmp = SeqDescrFree (tmp);
            ValNodeLink (&(bsp->descr), edited_set->descr);
            edited_set->descr = NULL;
            BioseqSetWriteBefore (bsp, aip_out, *orig, sp);
            wrote_front = TRUE;
          }
          SeqEntryCopyReplace (aip_in, aip_out, edited, atp, sp);
        }
      }
      else if (atp == sp->atp_annot)
      {
        bsp->annot = SeqAnnotSetAsnRead(aip_in, atp, sp->atp_annot_e);
				if (bsp->annot == NULL) goto erret;
      }
      else
      {
        if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;    /* takes care of everything else */
        if (atp == sp->atp_level)
          bsp->level = (Int2)av.intvalue;
        else if (atp == sp->atp_class)
			  {
          bsp->_class = (Uint1)av.intvalue;
			  }
        else if (atp == sp->atp_release)
          bsp->release = (CharPtr)av.ptrvalue;
      }
    }
  if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;   /* end BioseqSet */

  if (bsp->_class == BioseqseqSet_class_nuc_prot) {
    tmp_sep = SeqEntryNew ();
    tmp_sep->choice = 2;
    tmp_sep->data.ptrvalue = bsp;
    replace = FindReplacementSeqEntry (edited, tmp_sep);
    if (replace != NULL) {
      sep_return = AsnIoMemCopy (replace, (AsnReadFunc) SeqEntryAsnRead, (AsnWriteFunc) SeqEntryAsnWrite);
    } else {
      sep_return = SeqEntryNew();
      sep_return->choice = 2;
      sep_return->data.ptrvalue = bsp;
    }
  } else {
    BioseqSetWriteAfter (bsp, aip_out, atp, sp);
    bsp = BioseqSetFree (bsp);
  }
  
ret:
  AsnUnlinkType(*orig);     /*  unlink local tree */
  return sep_return;
erret:
  aip_in->io_failure = TRUE;
  aip_out->io_failure = TRUE;
  bsp = BioseqSetFree(bsp);
  goto ret;
}

static void SeqEntryCopyReplace (AsnIoPtr aip_in, AsnIoPtr aip_out, SeqEntryPtr edited, AsnTypePtr orig, SetAtpPtr sp)
{
	DataVal av;
	AsnTypePtr atp;
  SeqEntryPtr sep=NULL, tmp_sep;
  SeqEntryPtr replacement_sep;
	Uint1 type = 0;
  BioseqSetPtr bssp;

	if (aip_in == NULL || aip_out == NULL || sp == NULL)
		return;

	if (orig == NULL)           /* SeqEntry ::= (self contained) */
		atp = AsnReadId(aip_in, sp->amp, sp->atp_seqentry);
	else
		atp = AsnLinkType(orig, sp->atp_seqentry);    /* link in local tree */
	if (atp == NULL) return;

	sep = SeqEntryNew();
	if (sep == NULL) goto erret;

	if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;    /* read the CHOICE */

  atp = AsnReadId(aip_in, sp->amp, atp); 
  if (atp == NULL) goto erret;   /* which choice? */
  if (atp == sp->atp_seq)
  {
    sep->choice = 1;
    sep->data.ptrvalue = (Pointer) BioseqAsnRead(aip_in, atp);
    type = (Uint1)SM_BIOSEQ;
    replacement_sep = FindReplacementSeqEntry (edited, sep);
    if (replacement_sep != NULL) {
      sep = SeqEntryFree (sep);
      sep = replacement_sep;
    }
    av.ptrvalue = (Pointer)sep;
    if (!AsnWriteChoice(aip_out, orig, (Int2)sep->choice, &av)) goto erret;
    if (sep->choice == 1)
    {
      if (! BioseqAsnWrite((BioseqPtr)sep->data.ptrvalue, aip_out, sp->atp_seq))
		    goto erret;
    }
    else if (sep->choice == 2)
    {
      if (! BioseqSetAsnWrite((BioseqSetPtr)sep->data.ptrvalue, aip_out, sp->atp_set))
		      goto erret;
    }
    /* need to do this so that we don't free part of the edited set before we're done with it. */
    if (replacement_sep != NULL) {
      sep = NULL;
    }

  }
  else if (atp == sp->atp_set)
  {
    /* HERE, we need to read in the first part of the set, determine if it's a nuc-prot set.
     * if nuc-prot, read the whole set and look for the replacement seq-entry, then write
     * otherwise, write the front part, then loop through the seq-set recursively
     */

    tmp_sep = BioseqSetCopyReplace (aip_in, aip_out, edited, &atp, sp);

    if (tmp_sep != NULL) {
      av.ptrvalue = (Pointer)tmp_sep;
      AsnIoFlush (aip_out);
      if (!AsnWriteChoice(aip_out, sp->atp_se, (Int2)sep->choice, &av)) goto erret;
      if (tmp_sep->choice == 1)
      {
        if (! BioseqAsnWrite((BioseqPtr)tmp_sep->data.ptrvalue, aip_out, sp->atp_seq))
		      goto erret;
      }
      else if (tmp_sep->choice == 2)
      {
        if (! BioseqSetAsnWrite((BioseqSetPtr)tmp_sep->data.ptrvalue, aip_out, sp->atp_set))
		        goto erret;
      }
      tmp_sep = SeqEntryFree (tmp_sep);
    }
  }
  else if (atp == sp->atp_set_desc) 
  {
    /* write out descriptors from holding set instead */
    bssp = edited->data.ptrvalue;
    SeqDescrAsnWrite (bssp->descr, aip_out, sp->atp_set_desc);
  }

  sep = SeqEntryFree (sep);

ret:
  AsnUnlinkType(orig);      /*  unlink local tree */
	return;
erret:
  aip_in->io_failure = TRUE;
  aip_out->io_failure = TRUE;
	sep = SeqEntryFree(sep);
	goto ret;
}


NLM_EXTERN void ReintegrateFilteredAsn (SeqEntryPtr sep, FILE *orig_file, FILE *output, Boolean is_binary)
{
  AsnIoPtr       aip_in, aip_out;
  SetAtpPtr      sp;
  AsnTypePtr     atp, atp_ssp;
  DataVal        dv;

  if (orig_file == NULL || output == NULL) {
    return;
  }

  sp = GetSetAtp ();
  if (sp == NULL) {
    return;
  }

  atp_ssp = AsnFind ("Seq-submit");
  if (atp_ssp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    sp = MemFree (sp);
    return;
  }

  aip_in = AsnIoNew (is_binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, orig_file, NULL, NULL, NULL);
  if (aip_in == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file");
    sp = MemFree (sp);
    return;
  }

  aip_out = AsnIoNew (is_binary ? ASNIO_BIN_OUT : ASNIO_TEXT_OUT, output, NULL, NULL, NULL);

  if ((atp = AsnReadId (aip_in, sp->amp, atp_ssp)) != NULL) {
    AsnReadVal (aip_in, atp, &dv);
    AsnWrite (aip_out, atp, &dv);
    atp = AsnReadId (aip_in, sp->amp, atp);
  } else {
    AsnIoFree (aip_in, FALSE);
    rewind (orig_file);
    aip_in = AsnIoNew (is_binary ? ASNIO_BIN_IN : ASNIO_TEXT_IN, orig_file, NULL, NULL, NULL);
    atp = AsnReadId (aip_in, sp->amp, sp->atp_seqentry);
    if (atp == NULL) {
      AsnIoFree (aip_in, FALSE);
      rewind (orig_file);
      aip_in = AsnIoNew (is_binary ? ASNIO_BIN_IN : ASNIO_TEXT_IN, orig_file, NULL, NULL, NULL);
      atp = AsnReadId (aip_in, sp->amp, sp->atp_bss);
    } else {
      AsnReadVal (aip_in, atp, NULL);
      AsnWrite (aip_out, atp, &dv);
      atp = AsnReadId(aip_in, sp->amp, atp);
    }
  }
  if (atp == NULL) {
    AsnIoFree (aip_in, FALSE);
    AsnIoFree (aip_out, FALSE);
    sp = MemFree (sp);
    return;
  }

  while (! aip_in->io_failure && atp != NULL) {
    if (atp == sp->atp_se) {
      SeqEntryCopyReplace (aip_in, aip_out, sep, atp, sp);
    } else {
      AsnReadVal (aip_in, atp, &dv);
      AsnWrite (aip_out, atp, &dv);
    }
    atp = AsnReadId (aip_in, sp->amp, atp);
  }

  AsnIoFree (aip_in, FALSE);
  AsnIoFree (aip_out, FALSE);
  sp = MemFree (sp);
}


static CharPtr MakePubLabelString (PubdescPtr pdp)

{
  Char        buf [521];
  CitGenPtr   cgp;
  ValNodePtr  vnp;

  if (pdp == NULL) return NULL;

  vnp = pdp->pub;

  /* skip over just serial number */

  if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp != NULL) {
      if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
        if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
          vnp = vnp->next;
        }
      }
    }
  }

  if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
    return StringSaveNoNull (buf);
  }

  return NULL;
}


static CharPtr GetDescriptorLabel (SeqDescrPtr sdp)
{
  if (sdp == NULL) {
    return NULL;
  } else if (sdp->choice == Seq_descr_pub) {
    return MakePubLabelString (sdp->data.ptrvalue);
  } else {
    return NULL;
  }
}


NLM_EXTERN DescStreamPtr DescStreamNew (SeqDescPtr sdp, BioseqPtr parent)
{
  DescStreamPtr ds;


  ds = (DescStreamPtr) MemNew (sizeof (DescStreamData));
  if (sdp != NULL) {
    ds->orig = AsnIoMemCopy (sdp, (AsnReadFunc) SeqDescAsnRead, (AsnWriteFunc) SeqDescAsnWrite);
    ds->replace = AsnIoMemCopy (sdp, (AsnReadFunc) SeqDescAsnRead, (AsnWriteFunc) SeqDescAsnWrite);
    ds->text = GetDescriptorLabel(ds->orig);
  }
  if (parent != NULL) {
    ds->owners = SeqIdDup (SeqIdFindBest (parent->id, SEQID_GENBANK));
    ds->last_owner = ds->owners;
  }


  return ds;
}


NLM_EXTERN DescStreamPtr DescStreamFree (DescStreamPtr ds)
{
  if (ds != NULL) {
    ds->orig = SeqDescFree (ds->orig);
    ds->replace = SeqDescFree (ds->replace);
    ds->owners = SeqIdSetFree (ds->owners);
    ds = MemFree (ds);
  }
  return ds;
}


NLM_EXTERN ValNodePtr DescStreamListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->data.ptrvalue = DescStreamFree (vnp->data.ptrvalue);
    vnp->next = NULL;
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static Boolean DoDescriptorsMatch (SeqDescPtr sdp1, SeqDescPtr sdp2)
{
  if (sdp1 == NULL && sdp2 == NULL) {
    return TRUE;
  } else if (sdp1 == NULL || sdp2 == NULL) {
    return FALSE;
  } else if (sdp1->choice != sdp2->choice) {
    return FALSE;
  } else if (sdp1->choice == Seq_descr_pub) {
    return PubdescContentMatch (sdp1->data.ptrvalue, sdp2->data.ptrvalue);
  } else {
    return AsnIoMemComp (sdp1, sdp2, (AsnWriteFunc) SeqDescAsnWrite);
  }
}


static void AddToDescStream (ValNodeBlockPtr vb, SeqDescPtr sdp, BioseqPtr parent)
{
  DescStreamPtr dsp_new, dsp;
  CharPtr txt;
  ValNodePtr vnp, prev = NULL, vnp_new;
  Boolean add_to_prev = FALSE;

  if (vb == NULL) {
    return;
  }
  if (vb->head == NULL) {
    ValNodeAddPointerToEnd (vb, 0, DescStreamNew (sdp, parent));
  } else {
    txt = GetDescriptorLabel(sdp);
    vnp = vb->head;
    dsp = vnp->data.ptrvalue;
    while (vnp != NULL && StringCmp (txt, dsp->text) < 0) {
      prev = vnp;
      vnp = vnp->next;
      if (vnp != NULL) {
        dsp = vnp->data.ptrvalue;
      }
    }
    if (vnp == NULL) {
      ValNodeAddPointerToEnd (vb, 0, DescStreamNew (sdp, parent));
    } else {
      while (vnp != NULL && StringCmp (txt, dsp->text) == 0
        && !(add_to_prev = DoDescriptorsMatch (sdp, dsp->orig)) ) {
        prev = vnp;
        vnp = vnp->next;
        if (vnp != NULL) {
          dsp = vnp->data.ptrvalue;
        }
      }
      if (add_to_prev) {
        dsp->last_owner->next = SeqIdDup (SeqIdFindBest (parent->id, SEQID_GENBANK));
        dsp->last_owner = dsp->last_owner->next;
      } else {
        dsp_new = DescStreamNew (sdp, parent);
        vnp_new = ValNodeNew (NULL);
        vnp_new->data.ptrvalue = dsp_new;
        if (prev == NULL) {
          vb->head = vnp_new;
          vb->tail = vnp_new;
        } else {
          vnp_new->next = prev->next;
          prev->next = vnp_new;
          if (vnp_new->next == NULL) {
            vb->tail = vnp_new;
          }
        }
      }
      txt = MemFree (txt);
    }
  }
}


static int DescStreamCompare (DescStreamPtr ds1, DescStreamPtr ds2)
{
  if (ds1 == NULL && ds2 == NULL) {
    return 0;
  } else if (ds1 == NULL) {
    return -1;
  } else if (ds2 == NULL) {
    return 1;
  } else if (ds1->text == NULL) {
    return -1;
  } else if (ds2->text == NULL) {
    return 1;
  } else {
    return StringCmp (ds1->text, ds2->text);
  }
}


static int LIBCALLBACK SortVnpByDescStream (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  return DescStreamCompare(vnp1->data.ptrvalue, vnp2->data.ptrvalue);
}


static void RecombineDescStreamList (ValNodePtr PNTR p_list)
{
  ValNodePtr vnp, cmp, tmp;
  DescStreamPtr d1, d2;

  if (p_list == NULL || *p_list == NULL) {
    return;
  }

  *p_list = ValNodeSort (*p_list, SortVnpByDescStream);

  for (vnp = *p_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      d1 = (DescStreamPtr) vnp->data.ptrvalue;
      for (cmp = vnp->next; 
           cmp != NULL && (d2 = (DescStreamPtr) cmp->data.ptrvalue) != NULL && StringCmp (d1->text, d2->text) == 0;
           cmp = cmp->next) {
        if (cmp->choice == 0 && DoDescriptorsMatch (d1->orig, d2->orig)) {
          /* combine owner lists */
          if (d1->last_owner == NULL) {
            d1->owners = d2->owners;
            d1->last_owner = d1->owners;
          } else {
            d1->last_owner->next = d2->owners;
          }
          d2->owners = NULL;
          if (d1->last_owner != NULL) {
            while (d1->last_owner->next != NULL) {
              d1->last_owner = d1->last_owner->next;
            }
          }
            
          /* add dependencies */
          d1->num_dependent += d2->num_dependent;
          /* mark choice for later extraction and deletion */
          cmp->choice = 1;
        }
      }
    }
  }

  tmp = ValNodeExtractList (p_list, 1);
  tmp = DescStreamListFree (tmp);
}
  

static void AddPubCitationsFromFeat (SeqFeatPtr sfp, ValNodePtr desc_stream_list)
{
  ValNodePtr repl_v;
  DescStreamPtr d;
  PubdescPtr    pdp;
  ValNodePtr    vnp;
  ValNode       vn_p, vn_c;
  Boolean       found = FALSE;


  if (sfp == NULL || sfp->cit == NULL || sfp->cit->choice != 1 || sfp->cit->data.ptrvalue == NULL) 
  {
    return;
  }

  MemSet (&vn_p, 0, sizeof (ValNode));
  MemSet (&vn_c, 0, sizeof (ValNode));

  /* note - there could be multiple identical copies of a pub in the list,
   * we only need to count the match once - we will combine the totals
   * in RecombineDescStreamList.
   */
  for (repl_v = desc_stream_list; repl_v != NULL && !found; repl_v = repl_v->next) 
  {
    d = (DescStreamPtr) repl_v->data.ptrvalue;
    if (d->orig != NULL 
        && d->orig->choice == Seq_descr_pub
        && (pdp = (PubdescPtr) d->orig->data.ptrvalue) != NULL)
    {
      for (vnp = sfp->cit->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
        /* each vnp is a pub */
        vn_p.choice = PUB_Equiv;
        vn_p.data.ptrvalue = pdp->pub;
        vn_c.choice = PUB_Equiv;
        vn_c.data.ptrvalue = vnp;

        if (PubLabelMatch (&vn_p, &vn_c) == 0) 
        {
          d->num_dependent ++;
          found = TRUE;
        }
      }
    }
  }
}


static void AddPubCitationsFromAnnot (SeqAnnotPtr annot, ValNodePtr desc_stream_list)
{
  SeqFeatPtr sfp;

  if (annot == NULL || annot->type != 1) 
  {
    return;
  }
  for (sfp = annot->data; sfp != NULL; sfp = sfp->next) 
  {
    AddPubCitationsFromFeat (sfp, desc_stream_list);
  }
}


static void AddPubCitationsFromAnnotSet (SeqAnnotPtr annot, ValNodePtr desc_stream_list)
{
  while (annot != NULL)
  {
    AddPubCitationsFromAnnot (annot, desc_stream_list);
    annot = annot->next;
  }
}


static void AddPubCitationsFromSet (BioseqSetPtr bssp, ValNodePtr desc_stream_list)
{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;

  if (bssp == NULL || desc_stream_list == NULL) 
  {
    return;
  }

  AddPubCitationsFromAnnotSet (bssp->annot, desc_stream_list);
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
  {
    if (IS_Bioseq (sep) && (bsp = (BioseqPtr) sep->data.ptrvalue) != NULL) 
    {
      AddPubCitationsFromAnnotSet (bsp->annot, desc_stream_list);
    }
    else if (IS_Bioseq_set (sep))
    {
      AddPubCitationsFromSet (sep->data.ptrvalue, desc_stream_list);
    }
  }
}


static void FixCitationsInAnnot (SeqAnnotPtr annot, ValNodePtr desc_stream_list)
{
  SeqFeatPtr sfp;
  ValNodePtr repl_v;
  DescStreamPtr d;
  PubdescPtr    pdp, pdp_r;
  ValNodePtr    ppr, vnp, vnp_prev, vnp_next;
  ValNode       vn, vn_p, vn_c, vn_tmp;

  if (annot == NULL || annot->type != 1) 
  {
    return;
  }
  MemSet (&vn_p, 0, sizeof (ValNode));
  MemSet (&vn_c, 0, sizeof (ValNode));
  for (sfp = annot->data; sfp != NULL; sfp = sfp->next) 
  {
    if (sfp->cit == NULL || sfp->cit->choice != 1 || sfp->cit->data.ptrvalue == NULL) 
    {
      continue;
    }
    for (repl_v = desc_stream_list; repl_v != NULL; repl_v = repl_v->next) 
    {
      d = (DescStreamPtr) repl_v->data.ptrvalue;
      if (d->orig != NULL 
          && d->orig->choice == Seq_descr_pub
          && (pdp = (PubdescPtr) d->orig->data.ptrvalue) != NULL)
      {
        vnp_prev = NULL;
        for (vnp = sfp->cit->data.ptrvalue; vnp != NULL; vnp = vnp_next) {
          vnp_next = vnp->next;
          /* each vnp is a pub */
          vn_p.choice = PUB_Equiv;
          vn_p.data.ptrvalue = pdp->pub;
          vn_c.choice = PUB_Equiv;
          vn_c.data.ptrvalue = vnp;

          if (PubLabelMatch (&vn_p, &vn_c) == 0) 
          {
            if (d->replace != NULL
                && d->replace->choice == Seq_descr_pub
                && (pdp_r = (PubdescPtr) d->replace->data.ptrvalue) != NULL) 
            {
              /* update Seq-feat cit */
              MemSet ((Pointer) &vn, 0, sizeof (ValNode));
              MemCopy (&vn, sfp->cit, sizeof (ValNode));
              vn_p.choice = PUB_Equiv;
              vn_p.data.ptrvalue = pdp_r->pub;
              ppr = MinimizePub (&vn_p);
              
              vn_tmp.choice = vnp->choice;
              vn_tmp.data.ptrvalue = vnp->data.ptrvalue;
              vnp->choice = ppr->choice;
              vnp->data.ptrvalue = ppr->data.ptrvalue;
              ppr->choice = vn_tmp.choice;
              ppr->data.ptrvalue = vn_tmp.data.ptrvalue;
              ppr = PubFree (ppr);
              vnp_prev = vnp;
            } 
            else
            {
              /* remove Seq-feat Cit */
              if (vnp_prev == NULL) 
              {
                sfp->cit->data.ptrvalue = vnp->next;
              }
              else
              {
                vnp_prev->next = vnp->next;
              }
              vnp->next = NULL;
              vnp = PubFree (vnp);
            }
          }
          else
          {
            vnp_prev = vnp;
          }
        }
      }
    }
  }
}


static void FixCitationsInAnnotSet (SeqAnnotPtr sap, ValNodePtr desc_stream_list)
{
  while (sap != NULL) 
  {
    FixCitationsInAnnot (sap, desc_stream_list);
    sap = sap->next;
  }
}


static void FixCitationsInSet (BioseqSetPtr bssp, ValNodePtr desc_stream_list)
{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;

  if (bssp == NULL || desc_stream_list == NULL) 
  {
    return;
  }

  FixCitationsInAnnotSet (bssp->annot, desc_stream_list);
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
  {
    if (IS_Bioseq (sep) && (bsp = (BioseqPtr) sep->data.ptrvalue) != NULL) 
    {
      FixCitationsInAnnotSet (bsp->annot, desc_stream_list);
    }
    else if (IS_Bioseq_set (sep))
    {
      FixCitationsInSet (sep->data.ptrvalue, desc_stream_list);
    }
  }
}


typedef struct streamreader {
  ValNodeBlock desc_list;
  SeqDescrPtr parent_list;
  ValNodeBlock seqid_list;
} StreamReaderData, PNTR StreamReaderPtr;


static AsnTypePtr StreamingSkipElement (AsnIoPtr aip, AsnTypePtr orig, SetAtpPtr sp)
{
  AsnTypePtr   atp;
  DataVal av;

  if (AsnReadVal(aip, orig, &av) <= 0) return NULL;

  atp = AsnReadId(aip, sp->amp, orig); if (atp == NULL) return NULL;
  while (atp != orig && atp != NULL) {
    AsnReadVal(aip, atp, &av);
    AsnKillValue (atp, &av);
    atp = AsnReadId(aip, sp->amp, atp);
  }

  /* close structure */
  if (atp == orig) {
    AsnReadVal (aip, atp, &av);
    AsnKillValue (atp, &av);
  }
  return atp;
}


static void StreamingReadAny (AsnIoPtr aip, AsnTypePtr atp, SetAtpPtr sp, StreamReaderPtr sr);

static AsnTypePtr StreamingReadBioseqSet (AsnIoPtr aip, AsnTypePtr orig, SetAtpPtr sp, StreamReaderPtr sr)
{
  DataVal      av;
  AsnTypePtr   atp, oldatp;
  BioseqSetPtr bsp=NULL;
  SeqEntryPtr  curr, next;
  BioseqPtr    nuc_bsp;
  SeqDescPtr   sdp = NULL;
  SeqAnnotPtr  annot;
  SeqDescPtr   last_parent, first_parent;

	if (aip == NULL || sp == NULL || sr == NULL)
		return orig;

	if (orig == NULL)           /* BioseqSet ::= (self contained) */
		atp = AsnReadId(aip, sp->amp, sp->atp_bss);
	else
		atp = AsnLinkType(orig, sp->atp_bss);    /* link in local tree */

  oldatp = atp;
  if (atp == NULL) {
    return atp;
  }

	bsp = BioseqSetNew();
	if (bsp == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

  curr = NULL;

  while ((atp = AsnReadId(aip, sp->amp, atp)) != oldatp)
  {
	  if (atp == NULL) goto erret;
    if (atp == sp->atp_id)
	  {
      bsp->id = ObjectIdAsnRead(aip, atp);
		  if (bsp->id == NULL) goto erret;
	  }
    else if (atp == sp->atp_coll)
	  {
      bsp->coll = DbtagAsnRead(aip, atp);
		  if (bsp->coll == NULL) goto erret;
	  }
    else if (atp == sp->atp_date)
	  {
      bsp->date = DateAsnRead(aip, atp);
		  if (bsp->date == NULL) goto erret;
	  }
    else if (atp == sp->atp_set_desc)
	  {
      bsp->descr = SeqDescrAsnRead (aip, atp);
		  if (bsp->descr == NULL) goto erret;
	  }
    else if (atp == sp->atp_seqset && bsp->_class != BioseqseqSet_class_nuc_prot) 
    {
      first_parent = bsp->descr;
      bsp->descr = NULL;
      last_parent = sr->parent_list;
      if (last_parent == NULL) {
        sr->parent_list = first_parent;
      } else {
        while (last_parent->next != NULL) {
          last_parent = last_parent->next;
        }
        last_parent->next = first_parent;
      }
      /* reading members of set that is not nuc-prot */
      StreamingReadAny (aip, atp, sp, sr);
      if (last_parent == NULL) {
        sr->parent_list = NULL;
      } else {
        last_parent->next = NULL;
      }
      first_parent = SeqDescrFree (first_parent);
    }
    else if (atp == sp->atp_se)
    {
      /* reading members of set that is nuc-prot */
  	  if ((next = SeqEntryAsnRead(aip, atp)) != NULL)
	    {
		    if (IS_Bioseq(next))
			    SeqMgrConnect(SM_BIOSEQ, next->data.ptrvalue,
				                SM_BIOSEQSET, (Pointer) bsp);
		    else
			    SeqMgrConnect(SM_BIOSEQSET, next->data.ptrvalue,
				                SM_BIOSEQSET, (Pointer) bsp);


        if (curr == NULL)
		      bsp->seq_set = next;
  		  else
         curr->next = next;
        curr = next;
	    }
    }
    else if (atp == sp->atp_annot)
    {
      annot = SeqAnnotSetAsnRead(aip, atp, sp->atp_annot_e);
			if (annot == NULL) goto erret;
      if (bsp != NULL) {
        bsp->annot = annot;
      }
    }
    else
    {
      if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* takes care of everything else */
      if (atp == sp->atp_level)
        bsp->level = (Int2)av.intvalue;
      else if (atp == sp->atp_class)
		  {
        bsp->_class = (Uint1)av.intvalue;
		  }
      else if (atp == sp->atp_release)
      {
        bsp->release = (CharPtr)av.ptrvalue;
      }
    }
  }
  if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end BioseqSet */



ret:
  AsnUnlinkType(orig);     /*  unlink local tree */

  /* if this was a nuc-prot set, add it to the list now */
  if (bsp->_class == BioseqseqSet_class_nuc_prot) {
    if (bsp->seq_set != NULL && IS_Bioseq (bsp->seq_set)) {
      nuc_bsp = bsp->seq_set->data.ptrvalue;
      if (nuc_bsp != NULL) {
        ValNodeLinkToEnd (&(sr->seqid_list), SeqIdDup (SeqIdFindBest (nuc_bsp->id, SEQID_GENBANK)));
      }
      for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
        if (sdp->choice == Seq_descr_pub) {
          ValNodeAddPointerToEnd (&(sr->desc_list), 0, DescStreamNew (sdp, nuc_bsp));
        }
      }
      for (sdp = sr->parent_list; sdp != NULL; sdp = sdp->next) {
        if (sdp->choice == Seq_descr_pub) {
          ValNodeAddPointerToEnd (&(sr->desc_list), 0, DescStreamNew (sdp, nuc_bsp));
        }
      }
    }
    /* count feature citations */
    AddPubCitationsFromSet (bsp, sr->desc_list.head);
  }

  bsp = BioseqSetFree (bsp);

  return atp;
erret:
  aip->io_failure = TRUE;
  bsp = BioseqSetFree(bsp);
  goto ret;
}


static BioseqPtr LIBCALL StreamingReadBioseq (AsnIoPtr aip, AsnTypePtr orig, SetAtpPtr sp)
{
    DataVal av;
    AsnTypePtr atp;
    BioseqPtr bsp=NULL;
    Int2 level;

    if (aip == NULL)
        return bsp;

    if (! ProgMon("Read Bioseq"))
        return bsp;

    if (orig == NULL)           /* Bioseq ::= (self contained) */
        atp = AsnReadId(aip, sp->amp, sp->atp_bioseq);
    else
        atp = AsnLinkType(orig, sp->atp_bioseq);    /* link in local tree */
    if (atp == NULL) return bsp;

    bsp = BioseqNew();
    if (bsp == NULL) goto erret;

    level = AsnGetLevel(aip);     /* for skipping */

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    atp = AsnReadId(aip, sp->amp, atp); if (atp == NULL) goto erret;  /* id required, start struct */
    bsp->id = SeqIdSetAsnRead(aip, atp, sp->atp_bioseq_id_E);
    if (bsp->id == NULL) goto erret;

    atp = AsnReadId(aip, sp->amp, atp); if (atp == NULL) goto erret;
    if (atp == sp->atp_bioseq_desc)           /* descr optional */
    {
        bsp->descr = SeqDescrAsnRead (aip, atp);
		    if (bsp->descr == NULL) goto erret;
        atp = AsnReadId(aip, sp->amp, atp); if (atp == NULL) goto erret;
    }

    atp = StreamingSkipElement(aip, atp, sp);
    if (atp == NULL) goto erret;

    atp = AsnReadId(aip, sp->amp, atp); if (atp == NULL) goto erret;

    if (atp == sp->atp_bioseq_annot)
    {
        bsp->annot = SeqAnnotSetAsnRead(aip, atp, sp->atp_bioseq_annot_e);
        if (bsp->annot == NULL) goto erret;
        atp = AsnReadId(aip, sp->amp, atp); if (atp == NULL) goto erret;
    }

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end Bioseq */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return bsp;
erret:
    aip->io_failure = TRUE;
    bsp = BioseqFree(bsp);
    goto ret;
}


static void StreamingReadAny (AsnIoPtr aip, AsnTypePtr atp, SetAtpPtr sp, StreamReaderPtr sr)
{
  BioseqPtr      nuc_bsp;
  SeqDescrPtr    sdp = NULL;
  AsnTypePtr     atp_orig;
  Boolean        first = TRUE;
  DataVal        av;

  if (aip == NULL || sp == NULL || sr == NULL) {
    return;
  }
  atp_orig = atp;

  while (! aip->io_failure && atp != NULL && (first || atp != atp_orig)) {
    first = FALSE;
    if (atp == sp->atp_set) {
      atp = StreamingReadBioseqSet (aip, atp, sp, sr);
    } else {
      if (atp == sp->atp_seq) {
        nuc_bsp = StreamingReadBioseq (aip, atp, sp);
        if (nuc_bsp != NULL) {
          ValNodeLinkToEnd (&(sr->seqid_list), SeqIdDup (SeqIdFindBest (nuc_bsp->id, SEQID_GENBANK)));
        }
        for (sdp = nuc_bsp->descr; sdp != NULL; sdp = sdp->next) {
          if (sdp->choice == Seq_descr_pub) {
            AddToDescStream (&(sr->desc_list), sdp, nuc_bsp);
          }
        }
        for (sdp = sr->parent_list; sdp != NULL; sdp = sdp->next) {
          if (sdp->choice == Seq_descr_pub) {
            AddToDescStream (&(sr->desc_list), sdp, nuc_bsp);
          }
        }
        AddPubCitationsFromAnnotSet(nuc_bsp->annot, sr->desc_list.head);
        nuc_bsp = BioseqFree (nuc_bsp);
      } else if (atp == sp->atp_set_desc) {
        ValNodeLink (&(sr->parent_list), SeqDescrAsnRead (aip, atp));
      } else {
        AsnReadVal (aip, atp, &av);
        AsnKillValue (atp, &av);
      }
    }
    atp = AsnReadId (aip, sp->amp, atp);
  }
  if (atp == atp_orig) {
    AsnReadVal (aip, atp, NULL);
  }
}

static Boolean StreamingReadSeqEntry (AsnIoPtr aip, SetAtpPtr sp, StreamReaderPtr sr)
{
  AsnTypePtr atp;
  BioseqPtr  nuc_bsp;
  SeqDescPtr sdp;
  DataVal    av;

  atp = AsnReadId (aip, sp->amp, sp->atp_seqentry);
  if (atp == NULL) {
    return FALSE;
  }

  AsnReadVal (aip, atp, NULL);
  atp = AsnReadId(aip, sp->amp, atp);

  if (atp == sp->atp_set) {
    atp = StreamingReadBioseqSet (aip, atp, sp, sr);
  } else {
    if (atp == sp->atp_seq) {
      nuc_bsp = BioseqAsnRead (aip, atp);
      if (nuc_bsp != NULL) {
        ValNodeLinkToEnd (&(sr->seqid_list), SeqIdDup (SeqIdFindBest (nuc_bsp->id, SEQID_GENBANK)));
      }
      for (sdp = nuc_bsp->descr; sdp != NULL; sdp = sdp->next) {
        if (sdp->choice == Seq_descr_pub) {
          AddToDescStream (&(sr->desc_list), sdp, nuc_bsp);
        }
      }
      for (sdp = sr->parent_list; sdp != NULL; sdp = sdp->next) {
        if (sdp->choice == Seq_descr_pub) {
          AddToDescStream (&(sr->desc_list), sdp, nuc_bsp);
        }
      }
      nuc_bsp = BioseqFree (nuc_bsp);
    } else {
      AsnReadVal (aip, atp, &av);
      AsnKillValue (atp, &av);
    }
  }
  return TRUE;
}


static Boolean StreamingReadSeqSubmit (AsnIoPtr aip, SetAtpPtr sp, StreamReaderPtr sr)
{
  AsnTypePtr atp;
  BioseqPtr  nuc_bsp;
  SeqDescPtr sdp;
  SubmitBlockPtr sbp;

  atp = AsnReadId (aip, sp->amp, sp->atp_seqsubmit);
  if (atp == NULL) {
    return FALSE;
  }

  AsnReadVal (aip, atp, NULL);
  atp = AsnReadId(aip, sp->amp, atp);
  
  if (atp == sp->atp_sub) {
    sbp = SubmitBlockAsnRead (aip, atp);
    sbp = SubmitBlockFree (sbp);
    atp = AsnReadId (aip, sp->amp, atp);
  }

  while (atp != NULL) {
    if (atp == sp->atp_set) {
      atp = StreamingReadBioseqSet (aip, atp, sp, sr);
    } else {
      if (atp == sp->atp_seq) {
        nuc_bsp = BioseqAsnRead (aip, atp);
        if (nuc_bsp != NULL) {
          ValNodeLinkToEnd (&(sr->seqid_list), SeqIdDup (SeqIdFindBest (nuc_bsp->id, SEQID_GENBANK)));
        }
        for (sdp = nuc_bsp->descr; sdp != NULL; sdp = sdp->next) {
          if (sdp->choice == Seq_descr_pub) {
            AddToDescStream (&(sr->desc_list), sdp, nuc_bsp);
          }
        }
        for (sdp = sr->parent_list; sdp != NULL; sdp = sdp->next) {
          if (sdp->choice == Seq_descr_pub) {
            AddToDescStream (&(sr->desc_list), sdp, nuc_bsp);
          }
        }
        nuc_bsp = BioseqFree (nuc_bsp);
      } else {
        AsnReadVal (aip, atp, NULL);
      }
    }
    atp = AsnReadId (aip, sp->amp, atp);
  }
  return TRUE;
}


/* note - for now, we're just doing pubs.  later I'll find a way to create a text label
 * for other descriptors so that we can sort them also.
 * note - we also need to update seqfeatcits when we write out.
 */
NLM_EXTERN ValNodePtr StreamAsnForDescriptors (FILE *fp, Boolean is_binary, Boolean is_batch, Boolean is_submit, SeqIdPtr PNTR sip_list)
{
  AsnIoPtr       aip;
  SetAtpPtr      sp;
  StreamReaderData sr;
  Boolean          rval;
  ValNodePtr       tmp;

  if (fp == NULL) return NULL;

  sp = GetSetAtp ();
  if (sp == NULL) {
    return NULL;
  }

  aip = AsnIoNew (is_binary ? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file");
    sp = MemFree (sp);
    return NULL;
  }

  MemSet (&sr, 0, sizeof (StreamReaderData));
  if (sip_list != NULL) {
    InitValNodeBlock (&(sr.seqid_list), *sip_list);
  }

  if (is_submit) {
    StreamingReadSeqSubmit (aip, sp, &sr);
  } else if (is_batch) {
    while (rval = StreamingReadSeqEntry (aip, sp, &sr)) {
    }
  } else {
    StreamingReadSeqEntry (aip, sp, &sr);
  }

  AsnIoFree (aip, FALSE);
  sp = MemFree (sp);
  sr.parent_list = SeqDescrFree (sr.parent_list);

  /* combine list items */
  RecombineDescStreamList(&(sr.desc_list.head));

  if (sip_list == NULL) {
    sr.seqid_list.head = SeqIdSetFree(sr.seqid_list.head);
  } else {
    *sip_list = sr.seqid_list.head;
  }

  /* set up on-all */
  tmp = SeqIdListToValNodeSeqIdList (*sip_list);
  SetOnAllValsForDescStreamList(sr.desc_list.head, tmp);
  tmp = ValNodeSeqIdListFree (tmp);

  return sr.desc_list.head;
}


static SeqDescrPtr GetDescriptorsForBioseq (BioseqPtr bsp, ValNodePtr desc_stream_list)
{
  ValNodePtr vnp;
  DescStreamPtr d;
  SeqIdPtr      sip, sip_tmp;
  Boolean       found;
  SeqDescrPtr   sdp = NULL;

  if (bsp == NULL || desc_stream_list == NULL) {
    return NULL;
  }

  for (vnp = desc_stream_list; vnp != NULL; vnp = vnp->next) {
    d = (DescStreamPtr) vnp->data.ptrvalue;
    if (d->replace != NULL) {
      found = FALSE;
      if (d->on_all) {
        found = TRUE;
      } else {
        /* note - we can use just the best one, because that's the one that was copied */
        sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
        found = FALSE;
        for (sip_tmp = d->owners; sip_tmp != NULL && !found; sip_tmp = sip_tmp->next) {
          if (SeqIdComp(sip, sip_tmp) == SIC_YES) {
            found = TRUE;
          }
        }
      }
      if (found) {
        ValNodeLink (&sdp, AsnIoMemCopy (d->replace, (AsnReadFunc) SeqDescAsnRead, (AsnWriteFunc) SeqDescAsnWrite));
      }
    }
  }
  return sdp;
}


static SeqDescrPtr GetDescriptorsForBioseqSet (BioseqSetPtr bssp, ValNodePtr desc_stream_list)
{
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_nuc_prot 
      || bssp->seq_set == NULL
      || desc_stream_list == NULL) 
  {
    return NULL;
  }

  if (IS_Bioseq(bssp->seq_set)) 
  {
    return GetDescriptorsForBioseq (bssp->seq_set->data.ptrvalue, desc_stream_list);
  } 
  else 
  {
    /* this had better be a segset */
    bssp = bssp->seq_set->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset 
        && bssp->seq_set != NULL
        && IS_Bioseq (bssp->seq_set))
    {
      return GetDescriptorsForBioseq (bssp->seq_set->data.ptrvalue, desc_stream_list);
    }
  }
  return NULL;
}


static void StreamingReadWriteAny (AsnIoPtr aip_in, AsnIoPtr aip_out, AsnTypePtr atp, SetAtpPtr sp, ValNodePtr desc_stream_list);

static AsnTypePtr 
StreamingReadWriteBioseqSet 
(AsnIoPtr   aip_in, 
 AsnIoPtr   aip_out, 
 ValNodePtr desc_stream_list, 
 AsnTypePtr orig, 
 AsnTypePtr set_top,
 SetAtpPtr  sp)
{
  DataVal      av;
  AsnTypePtr   atp, oldatp;
  BioseqSetPtr bsp=NULL;
  SeqEntryPtr  curr, next;
  SeqDescrPtr  tmp;
  SeqAnnotPtr  annot;

	if (aip_in == NULL || aip_out == NULL || sp == NULL)
		return orig;

	if (orig == NULL)           /* BioseqSet ::= (self contained) */
		atp = AsnReadId(aip_in, sp->amp, sp->atp_bss);
	else
		atp = AsnLinkType(orig, sp->atp_bss);    /* link in local tree */

  oldatp = atp;
  if (atp == NULL) {
    return atp;
  }

	bsp = BioseqSetNew();
	if (bsp == NULL) goto erret;

	if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;    /* read the start struct */

  curr = NULL;

  while ((atp = AsnReadId(aip_in, sp->amp, atp)) != oldatp)
  {
	  if (atp == NULL) goto erret;
    if (atp == sp->atp_id)
	  {
      bsp->id = ObjectIdAsnRead(aip_in, atp);
		  if (bsp->id == NULL) goto erret;
	  }
    else if (atp == sp->atp_coll)
	  {
      bsp->coll = DbtagAsnRead(aip_in, atp);
		  if (bsp->coll == NULL) goto erret;
	  }
    else if (atp == sp->atp_date)
	  {
      bsp->date = DateAsnRead(aip_in, atp);
		  if (bsp->date == NULL) goto erret;
	  }
    else if (atp == sp->atp_set_desc)
	  {
      bsp->descr = SeqDescrAsnRead(aip_in, atp);
		  if (bsp->descr == NULL) goto erret;
      /* remove descriptors that are being streamed and edited.
       * for now, this is just pubs.
       */
      tmp = ValNodeExtractList (&(bsp->descr), Seq_descr_pub);
      tmp = SeqDescrFree (tmp);
	  }
    else if (atp == sp->atp_seqset && bsp->_class != BioseqseqSet_class_nuc_prot) 
    {
      /* reading members of set that is not nuc-prot */
      StreamingReadWriteAny (aip_in, aip_out, atp, sp, desc_stream_list);
    }
    else if (atp == sp->atp_se)
    {
      if (bsp == NULL) 
      {
      } 
      else 
      {
        /* reading members of set that is nuc-prot */
	  	  if ((next = SeqEntryAsnRead(aip_in, atp)) != NULL)
		    {
			    if (IS_Bioseq(next))
				    SeqMgrConnect(SM_BIOSEQ, next->data.ptrvalue,
					                SM_BIOSEQSET, (Pointer) bsp);
			    else
				    SeqMgrConnect(SM_BIOSEQSET, next->data.ptrvalue,
					                SM_BIOSEQSET, (Pointer) bsp);


          if (curr == NULL)
			      bsp->seq_set = next;
    		  else
           curr->next = next;
          curr = next;
		    }
      }
    }
    else if (atp == sp->atp_annot)
    {
      annot = SeqAnnotSetAsnRead(aip_in, atp, sp->atp_annot_e);
			if (annot == NULL) goto erret;
      bsp->annot = annot;
    }
    else
    {
      if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;    /* takes care of everything else */
      if (atp == sp->atp_level)
        bsp->level = (Int2)av.intvalue;
      else if (atp == sp->atp_class)
		  {
        bsp->_class = (Uint1)av.intvalue;
        if (bsp->_class != BioseqseqSet_class_nuc_prot) {
          /* remove descriptors that are being streamed and replaced */
          tmp = ValNodeExtractList (&(bsp->descr), Seq_descr_pub);
          tmp = SeqDescrFree (tmp);
          BioseqSetWriteBefore (bsp, aip_out, orig, sp);
        }
		  }
      else if (atp == sp->atp_release)
      {
        if (bsp != NULL) {
          bsp->release = (CharPtr)av.ptrvalue;
        }
      }
    }
  }
  if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;   /* end BioseqSet */



ret:

  AsnUnlinkType(orig);     /*  unlink local tree */

  FixCitationsInSet (bsp, desc_stream_list);
  /* if this was a nuc-prot set, change descriptors and write out now */
  if (bsp->_class == BioseqseqSet_class_nuc_prot) {
    ValNodeLink (&(bsp->descr), GetDescriptorsForBioseqSet(bsp, desc_stream_list));
    AsnWriteChoice(aip_out, set_top, 2, &av);
    BioseqSetAsnWrite(bsp, aip_out, atp);  
  } else {
    BioseqSetWriteAfter (bsp, aip_out, atp, sp);
  }

  bsp = BioseqSetFree (bsp);

  return atp;
erret:
  aip_in->io_failure = TRUE;
  bsp = BioseqSetFree(bsp);
  goto ret;
}


static void StreamingReadWriteBioseq (AsnIoPtr aip_in, AsnIoPtr aip_out, AsnTypePtr atp, ValNodePtr desc_stream_list)
{
  BioseqPtr   nuc_bsp;
  SeqDescrPtr tmp;

  nuc_bsp = BioseqAsnRead (aip_in, atp);
  tmp = ValNodeExtractList (&(nuc_bsp->descr), Seq_descr_pub);
  tmp = SeqDescrFree (tmp);
  ValNodeLink (&(nuc_bsp->descr), GetDescriptorsForBioseq (nuc_bsp, desc_stream_list));
  FixCitationsInAnnotSet (nuc_bsp->annot, desc_stream_list);
  BioseqAsnWrite (nuc_bsp, aip_out, atp);
  nuc_bsp = BioseqFree (nuc_bsp);
}


static void StreamingReadWriteAny (AsnIoPtr aip_in, AsnIoPtr aip_out, AsnTypePtr atp, SetAtpPtr sp, ValNodePtr desc_stream_list)
{
  AsnTypePtr     atp_orig;
  Boolean        first = TRUE;
  DataVal       av;

  if (aip_in == NULL || aip_out == NULL || sp == NULL) {
    return;
  }
  atp_orig = atp;

  while (! aip_in->io_failure && atp != NULL && (first || atp != atp_orig)) {
    first = FALSE;
    if (atp == sp->atp_set) {
      atp = StreamingReadWriteBioseqSet (aip_in, aip_out, desc_stream_list, atp, sp->atp_se, sp);
    } else {
      if (atp == sp->atp_seq) {
        AsnWriteChoice(aip_out, sp->atp_se, (Int2)1, &av);
        StreamingReadWriteBioseq (aip_in, aip_out, atp, desc_stream_list);
      } else {
        AsnReadVal (aip_in, atp, NULL);
      }
    }
    atp = AsnReadId (aip_in, sp->amp, atp);
  }
  if (atp == atp_orig) {
    AsnReadVal (aip_in, atp, NULL);
  }
}


static Boolean StreamingReadWriteSeqEntry (ValNodePtr desc_stream_list, AsnIoPtr aip_in, AsnIoPtr aip_out, SetAtpPtr sp)
{
  AsnTypePtr atp;
  DataVal av;

  atp = AsnReadId (aip_in, sp->amp, sp->atp_seqentry);
  if (atp == NULL) {
    return FALSE;
  }

  AsnReadVal (aip_in, atp, NULL);
  atp = AsnReadId(aip_in, sp->amp, atp);

  if (atp == sp->atp_set) {
    atp = StreamingReadWriteBioseqSet (aip_in, aip_out, desc_stream_list, atp, sp->atp_se, sp);
  } else if (atp == sp->atp_seq) {
    /* first write Seq-entry lead-in */
    AsnWriteChoice(aip_out, sp->atp_se, (Int2)1, &av);
    StreamingReadWriteBioseq (aip_in, aip_out, atp, desc_stream_list);
  } else {
    AsnReadVal (aip_in, atp, NULL);
  }
  return TRUE;
}


static Boolean StreamingReadWriteSeqSubmit (ValNodePtr desc_stream_list, AsnIoPtr aip_in, AsnIoPtr aip_out, SetAtpPtr sp)
{
  AsnTypePtr atp, oldatp;
  DataVal av;
  SeqSubmitPtr ssp;
  SubmitBlockPtr sbp;

  atp = AsnReadId (aip_in, sp->amp, sp->atp_seqsubmit);
  if (atp == NULL) {
    return FALSE;
  }

  AsnReadVal (aip_in, atp, NULL);

  ssp = SeqSubmitNew ();
  if (! AsnOpenStruct(aip_out, atp, (Pointer)ssp)) {
    ssp = SeqSubmitFree (ssp);
    return FALSE;
  }

  atp = AsnReadId(aip_in, sp->amp, atp);
  
  if (atp == sp->atp_sub) 
  {
    sbp = SubmitBlockAsnRead (aip_in, atp);
    SubmitBlockAsnWrite (sbp, aip_out, atp);
    sbp = SubmitBlockFree (sbp);
    atp = AsnReadId (aip_in, sp->amp, atp);
  }

  if (atp == NULL) goto erret;
  if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;

  atp = AsnReadId(aip_in, sp->amp, atp);  /* read the data */
  if (atp == NULL) goto erret;
	oldatp = atp;     /* the SET OF */
  if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;

  while ((atp = AsnReadId(aip_in, sp->amp, atp)) != oldatp && atp != NULL)
  {
	  if (atp == sp->atp_seqsubmit_data_entries_E)
	  {
      if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;
      if (! AsnWriteChoice(aip_out, sp->atp_seqsubmit_data, (Int2)1, &av)) goto erret;
      if (! AsnOpenStruct(aip_out, sp->atp_seqsubmit_data_entries, ssp->data)) goto erret; 
    } else if (atp == sp->atp_set) {
      atp = StreamingReadWriteBioseqSet (aip_in, aip_out, desc_stream_list, atp, sp->atp_seqsubmit_data_entries_E, sp);
    } else if (atp == sp->atp_seq) {
      /* first write Seq-entry lead-in */
      AsnWriteChoice(aip_out, sp->atp_se, (Int2)1, &av);
      StreamingReadWriteBioseq (aip_in, aip_out, atp, desc_stream_list);
    }
  }
  if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;   /* end set of */
  if (! AsnCloseStruct(aip_out, atp, (Pointer) ssp)) goto erret;

  atp = AsnReadId(aip_in, sp->amp, atp);
  if (atp == NULL) goto erret;
  if (AsnReadVal(aip_in, atp, &av) <= 0) goto erret;  /* end struct */
  if (! AsnCloseStruct(aip_out, atp, (Pointer)ssp)) goto erret;

ret:
  ssp = SeqSubmitFree (ssp);
	return TRUE;
erret:
  aip_in->io_failure = TRUE;
  ssp = SeqSubmitFree(ssp);
  return FALSE;

}


NLM_EXTERN void WriteAsnWithReplacedDescriptors (ValNodePtr desc_stream_list, FILE *orig_file, FILE *output, Boolean is_binary, Boolean is_batch, Boolean is_submit)
{
  AsnIoPtr       aip_in, aip_out;
  SetAtpPtr      sp;
  Boolean        rval;

  if (orig_file == NULL || output == NULL) {
    return;
  }

  sp = GetSetAtp ();
  if (sp == NULL) {
    return;
  }

  aip_in = AsnIoNew (is_binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, orig_file, NULL, NULL, NULL);
  if (aip_in == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file");
    sp = MemFree (sp);
    return;
  }

  aip_out = AsnIoNew (is_binary ? ASNIO_BIN_OUT : ASNIO_TEXT_OUT, output, NULL, NULL, NULL);

  if (is_submit) {
    StreamingReadWriteSeqSubmit (desc_stream_list, aip_in, aip_out, sp);
    AsnIoFlush (aip_out);
  } else if (is_batch) {
    while (rval = StreamingReadWriteSeqEntry(desc_stream_list, aip_in, aip_out, sp)) {
      AsnIoReset (aip_out);
    }
  } else {
    rval = StreamingReadWriteSeqEntry(desc_stream_list, aip_in, aip_out, sp);
    AsnIoFlush (aip_out);
  }
  AsnIoFree (aip_in, FALSE);
  AsnIoFree (aip_out, FALSE);

  sp = MemFree (sp);
}


NLM_EXTERN Boolean IdListsMatch (SeqIdPtr sip_list, ValNodePtr all_sip)
{
  Boolean found = FALSE, any_missing = FALSE;
  ValNodePtr vnp;

  if (sip_list == NULL || all_sip == NULL) {
    return FALSE;
  }

  if (ValNodeLen (sip_list) != ValNodeLen (all_sip)) {
    return FALSE;
  }

  while (sip_list != NULL) {
    found = FALSE;
    for (vnp = all_sip; vnp != NULL && !found; vnp = vnp->next) {
      if (vnp->choice == 0 && SeqIdComp (vnp->data.ptrvalue, sip_list) == SIC_YES) {
        vnp->choice = 1;
        found = TRUE;
      }
    }
    sip_list = sip_list->next;
  }
  for (vnp = all_sip; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      any_missing = TRUE;
    }
    vnp->choice = 0;
  }
  return !any_missing;
}


NLM_EXTERN void SetOnAllValsForDescStreamList (ValNodePtr desc_list, ValNodePtr all_sip)
{
  ValNodePtr vnp;
  DescStreamPtr d;

  for (vnp = desc_list; vnp != NULL; vnp = vnp->next) {
    d = (DescStreamPtr) vnp->data.ptrvalue;
    d->on_all = IdListsMatch(d->owners, all_sip);
  }
}


/* ReadAsnFastaOrFlatFileEx reads lines, looking for starts of ASN.1, FASTA, GenBank, EMBL,
or GenPept files.  It then calls the appropriate read function, which is responsible for
reading the sequence (or object) and restoring the file pointer to the beginning of the
next record. */

NLM_EXTERN Pointer ReadAsnFastaOrFlatFileEx (
  FILE *fp,
  Uint2Ptr datatypeptr,
  Uint2Ptr entityIDptr,
  Boolean forceNuc,
  Boolean forceProt,
  Boolean parseFastaSeqId,
  Boolean fastaAsSimpleSeq,
  BoolPtr chars_stripped
)

{
  AsnIoPtr       aip;
  CharPtr        annotname;
  Int4           begin;
  ByteStorePtr   bs = NULL;
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp;
  Char           ch;
  Uint1          choice = 0;
  Boolean        coordinatesOnMaster;
  Uint2          datatype;
  Int2           db = -1;
  FileCache      fc;
  Boolean        inLetters;
  Boolean        isProt = FALSE;
  Int4           j;
  long           len;
  Char           line [4096];
  Boolean        mayBeAccessionList = TRUE;
  Boolean        mayBePlainFasta = TRUE;
  SeqFeatPtr     nextsfp;
  Int2           numDigits;
  Int2           numLetters;
  Int4           numLinks;
  ObjMgrDataPtr  omdp;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp = NULL;
  PubdescPtr     pdp;
  Int4           pos;
  ValNodePtr     pip;
  Pointer PNTR   prevsfp;
  ProjectPtr     proj = NULL;
  BoolPtr        protPtr;
  Pointer        ptr = NULL;
  SeqAnnotPtr    sap = NULL;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp;
  Char           seqid [2048];
  SimpleSeqPtr   ssp = NULL;
  CharPtr        str;
  CharPtr        tag;
  CharPtr        title = NULL;
  CharPtr        tmp;
  Int4           uid;
  long int       val;
  ValNodePtr     vnp;
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (fp == NULL) return NULL;

  if (datatypeptr != NULL) *datatypeptr = 0;
  if (entityIDptr != NULL) *entityIDptr = 0;

  if (forceNuc) {
    isProt = FALSE;
    protPtr = NULL;
  } else if (forceProt) {
    isProt = TRUE;
    protPtr = NULL;
  } else {
    protPtr = &isProt;
  }

  seqid [0] = '\0';

  FileCacheSetup (&fc, fp);

  pos = FileCacheTell (&fc);
  begin = pos;
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);

  if (str == NULL) return NULL; /* already at end of file */

  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringStr (line, "::=") != NULL) {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;

        /* first skip past empty space at start of line */

        tag = line;
        ch = *tag;
        while (ch != '\0' && IS_WHITESP (ch)) {
          tag++;
          ch = *tag;
        }

        /* now find ASN tag */

        tmp = tag;
        ch = *tmp;
        while (ch != '\0' && (! IS_WHITESP (ch))) {
          tmp++;
          ch = *tmp;
        }
        *tmp = '\0';

        omp = ObjMgrReadLock ();
        omtp = ObjMgrTypeFind (omp, 0, tag, NULL);
        ObjMgrUnlock ();

        if (omtp != NULL) {
          FileCacheFree (&fc, FALSE);
          fseek (fp, pos, SEEK_SET);
          aip = AsnIoNew (ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
          aip->scan_for_start = TRUE;
          SeqMgrHoldIndexing (TRUE);
          ptr = (*(omtp->asnread)) (aip, NULL);
          SeqMgrHoldIndexing (FALSE);
          pos = AsnIoTell (aip);
          AsnIoFree (aip, FALSE);
          FileCacheSetup (&fc, fp);
          FileCacheSeek (&fc, pos);
          fseek (fp, pos, SEEK_SET);

          if (ptr == NULL) {
            ErrPostEx (SEV_ERROR, 0, 0, "Couldn't read type [%s]", omtp->asnname);
          } else {
            datatype = omtp->datatype;
            if (datatypeptr != NULL) {
              *datatypeptr = datatype;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (datatype, ptr);
            }
            if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
              omp = ObjMgrReadLock ();
              omdp = ObjMgrFindByData (omp, ptr);
              ObjMgrUnlock ();
              if (omdp != NULL && omdp->choice == NULL) {
                /* always want sep above bsp or bssp */
                sep = SeqEntryNew ();
                if (sep != NULL) {
                  if (datatype == OBJ_BIOSEQ) {
                    bsp = (BioseqPtr) ptr;
                    sep->choice = 1;
                    sep->data.ptrvalue = bsp;
                    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
                  } else if (datatype == OBJ_BIOSEQSET) {
                    bssp = (BioseqSetPtr) ptr;
                    sep->choice = 2;
                    sep->data.ptrvalue = bssp;
                    SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
                  } else {
                    sep = SeqEntryFree (sep);
                  }
                }
              }
            }
          }
        } else {
          ErrPostEx (SEV_ERROR, 0, 0, "Couldn't read type [%s]", line);
        }
        return ptr;

      } else if (line [0] == '>') {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        db = -1;
        if (StringNCmp (line, ">PubMed", 7) == 0) {
          db = 0;
        } else if (StringNCmp (line, ">Protein", 8) == 0) {
          db = 1;
        } else if (StringNCmp (line, ">Nucleotide", 11) == 0) {
          db = 2;
        } else if (StringNCmp (line, ">Structure", 10) == 0) {
          db = 3;
        } else if (StringNCmp (line, ">Genome", 7) == 0) {
          db = 4;
        }
        if (db != -1) {

          bs = ReadUidList (&fc, (Boolean) (db == 2), FALSE);
          if (bs != NULL) {
            proj = ProjectNew ();
            if (proj != NULL) {
              pip = ValNodeNew (NULL);
              if (pip != NULL) {
                switch (db) {
                  case 0 :
                    choice = ProjectItem_pmuid;
                    break;
                  case 1 :
                    choice = ProjectItem_protuid;
                    break;
                  case 2 :
                    choice = ProjectItem_nucuid;
                    break;
                  case 3 :
                    choice = ProjectItem_genomeuid;
                    break;
                  case 4 :
                    choice = ProjectItem_structuid;
                    break;
                  default :
                    choice = 0;
                    break;
                }
                pip->choice = choice;
                proj->data = pip;
                numLinks = BSLen (bs) / sizeof (Int4);
                BSSeek (bs, 0L, 0);
                for (j = 0; j < numLinks; j++) {
                  BSRead (bs, &uid, sizeof (Int4));
                  ValNodeAddInt ((ValNodePtr PNTR) &(pip->data.ptrvalue), choice, uid);
                  /*
                  switch (db) {
                    case 0 :
                      ValNodeAddInt (&(pip->data.ptrvalue), ProjectItem_pmid, uid);
                      break;
                    case 1 :
                    case 2 :
                    case 3 :
                      sip = ValNodeNew (NULL);
                      if (sip != NULL) {
                        sip->choice = SEQID_GI;
                        sip->data.intvalue = uid;
                      }
                      break;
                    case 4 :
                      break;
                    default :
                      break;
                  }
                  */
                }
              }
            }
            bs = BSFree (bs);

            if (datatypeptr != NULL) {
              *datatypeptr = OBJ_PROJECT;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_PROJECT, (Pointer) proj);
            }

            pos = FileCacheTell (&fc);
            FileCacheSetup (&fc, fp);
            FileCacheSeek (&fc, pos);
            fseek (fp, pos, SEEK_SET);

            return (Pointer) proj;
          }

        } else if (StringNICmp (line, ">Feature", 8) == 0) {

          annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE, FALSE);
          if (! HasNoText (seqid)) {
            sap = ReadFeatureTable (&fc, seqid, annotname);
            if (sap != NULL && sap->type == 1) {
              sfp = (SeqFeatPtr) sap->data;
              prevsfp = (Pointer PNTR) &(sap->data);
              while (sfp != NULL) {
                nextsfp = sfp->next;
                if (sfp->data.choice == SEQFEAT_PUB) {
                  pdp = (PubdescPtr) sfp->data.value.ptrvalue;
                  if (pdp != NULL && pdp->pub == NULL) {
                    *(prevsfp) = sfp->next;
                    sfp->next = NULL;
                    SeqFeatFree (sfp);
                  } else {
                    prevsfp = (Pointer PNTR) &(sfp->next);
                  }
                } else {
                  prevsfp = (Pointer PNTR) &(sfp->next);
                }
                sfp = nextsfp;
              }
              if (sap->data == NULL) {
                sap = SeqAnnotFree (sap);
              }
            }
            if (sap != NULL) {
              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_SEQANNOT;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
              }

              pos = FileCacheTell (&fc);
              FileCacheSetup (&fc, fp);
              FileCacheSeek (&fc, pos);
              fseek (fp, pos, SEEK_SET);

              return (Pointer) sap;
            }
          }

        } else if (StringNICmp (line, ">Vector", 7) == 0) {

          annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE, FALSE);
          if (! HasNoText (seqid)) {
            sap = ReadVecScreenTable (&fc, seqid, annotname);
            if (sap != NULL) {
              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_SEQANNOT;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
              }

              pos = FileCacheTell (&fc);
              FileCacheSetup (&fc, fp);
              FileCacheSeek (&fc, pos);
              fseek (fp, pos, SEEK_SET);

              return (Pointer) sap;
            }
          }

        } else if (StringNICmp (line, ">Restriction", 12) == 0) {

          annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE, TRUE);
          if (! HasNoText (seqid)) {
            sap = ReadRestrictionSiteTable (&fc, seqid, annotname);
            if (sap != NULL) {
              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_SEQANNOT;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
              }

              pos = FileCacheTell (&fc);
              FileCacheSetup (&fc, fp);
              FileCacheSeek (&fc, pos);
              fseek (fp, pos, SEEK_SET);

              return (Pointer) sap;
            }
          }

        } else if (StringNICmp (line, ">Assembly", 9) == 0) {

          coordinatesOnMaster = FALSE;
          if (StringISearch (line, "Master") != NULL) {
            coordinatesOnMaster = TRUE;
          }
          annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE, FALSE);
          sep = ReadContigListExEx (&fc, coordinatesOnMaster, seqid, annotname);
          if (sep != NULL && IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {

              oip = ObjectIdNew ();
              oip->str = StringSave ("info");
              uop = UserObjectNew ();
              uop->type = oip;
              uop->_class = StringSave ("Genomes");

              oip = ObjectIdNew ();
              oip->id = 0;
              ufp = UserFieldNew ();
              ufp->choice = 2;
              ufp->data.intvalue = 0;
              ufp->label = oip;

              uop->data = ufp;

              vnp = SeqDescrNew (NULL);
              vnp->choice = Seq_descr_user;
              vnp->data.ptrvalue = (Pointer) uop;
              vnp->next = bsp->descr;
              bsp->descr = vnp;

              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_BIOSEQ;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
              }
            }

            pos = FileCacheTell (&fc);
            FileCacheSetup (&fc, fp);
            FileCacheSeek (&fc, pos);
            fseek (fp, pos, SEEK_SET);

            return (Pointer) bsp;
          }

        } else if (StringNICmp (line, ">Virtual", 8) == 0) {

          tmp = GetSeqId (seqid, line, sizeof (seqid), TRUE, TRUE);
          if (! HasNoText (seqid)) {
            TrimSpacesAroundString (tmp);
            if (tmp != NULL && sscanf (tmp, "%ld", &val) == 1) {
              sep = SeqEntryNew ();
              if (sep != NULL) {
                bsp = BioseqNew ();
                if (bsp != NULL) {
                  sep->choice = 1;
                  sep->data.ptrvalue = (Pointer) bsp;
                  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

                  bsp->mol = Seq_mol_na;

                  bsp->repr = Seq_repr_virtual;
                  bsp->length = (Int4) val;
                  bsp->id = MakeSeqID (seqid);
                  SeqMgrAddToBioseqIndex (bsp);

                  if (datatypeptr != NULL) {
                    *datatypeptr = OBJ_BIOSEQ;
                  }
                  if (entityIDptr != NULL) {
                    *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
                  }
                }
              }

              pos = FileCacheTell (&fc);
              FileCacheSetup (&fc, fp);
              FileCacheSeek (&fc, pos);
              fseek (fp, pos, SEEK_SET);

              return (Pointer) bsp;
            }
          }

        } else if (StringNICmp (line, ">Message", 8) == 0) {

          ReadMessageStrings (&fc);

        } else if (line [1] == '?') {

          sep = SeqEntryNew ();
          if (sep != NULL) {
            bsp = BioseqNew ();
            if (bsp != NULL) {
              sep->choice = 1;
              sep->data.ptrvalue = (Pointer) bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

              tmp = line + 2;
              ch = *tmp;
              while (IS_WHITESP (ch)) {
                tmp++;
                ch = *tmp;
              }
              if (StringNCmp (tmp, "unk100", 6) == 0) {
                bsp->id = MakeSeqID ("lcl|unk100");
                tmp += 3;
              } else {
                bsp->id = MakeSeqID ("lcl|gap");
              }
              SeqMgrAddToBioseqIndex (bsp);

              bsp->repr = Seq_repr_virtual;
              if (*tmp != '\0' && sscanf (tmp, "%ld", &len) == 1 && len > 0) {
                bsp->length = (Int4) len;
              } else {
                bsp->length = -1;
              }
              if (isProt) {
                bsp->mol = Seq_mol_aa;
                bsp->seq_data_type = Seq_code_ncbieaa;
              } else {
                bsp->mol = Seq_mol_na;
                bsp->seq_data_type = Seq_code_iupacna;
              }

              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_BIOSEQ;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
              }
            }
          }

          pos = FileCacheTell (&fc);
          FileCacheSetup (&fc, fp);
          FileCacheSeek (&fc, pos);
          fseek (fp, pos, SEEK_SET);

          return (Pointer) bsp;

        } else {

          title = NULL;
          tmp = StringChr (line + 1, '[');
          if (tmp != NULL) {
            if (parseFastaSeqId)
            {
              if (StringStr (tmp, "[") != NULL && StringStr (tmp, "=") != NULL) {
                TrimSpacesAroundString (tmp);
                title = StringSave (tmp);
              }
            }
            else
            {
              title = StringSaveNoNull (line + 1);
              TrimSpacesAroundString (title);
            }
          } else if (fastaAsSimpleSeq) {
            if (parseFastaSeqId)
            {
              tmp = StringChr (line + 1, ' ');
            }
            else
            {
              tmp = line + 1;
            }
            if (tmp != NULL) {
              tmp++;
              TrimSpacesAroundString (tmp);
              title = StringSaveNoNull (tmp);
            }
          }
          else if (!parseFastaSeqId)
          {
            title = StringSaveNoNull (line + 1);
          }
          if (parseFastaSeqId) {
            tmp = line + 1;
            ch = *tmp;
            while (IS_WHITESP (ch)) {
              tmp++;
              ch = *tmp;
            }
            if (ch == '[') {
              parseFastaSeqId = FALSE;
            }
          }
          if (parseFastaSeqId) {
            GetSeqId (seqid, line + 1, sizeof (seqid), FALSE, TRUE);
            if (! HasNoText (seqid)) {
              tmp = StringStr (line + 1, seqid);
              if (tmp != NULL) {
                tmp += StringLen (seqid);
                if (! StringHasNoText (tmp)) {
                  TrimSpacesAroundString (tmp);
                  title = MemFree (title);
                  title = StringSaveNoNull (tmp);
                }
              }
              bs = ReadFlatFileDNA (&fc, protPtr, forceNuc, forceProt, fastaAsSimpleSeq,
                                    FALSE, chars_stripped, seqid);
            }
          } else {
            bs = ReadFlatFileDNA (&fc, protPtr, forceNuc, forceProt, fastaAsSimpleSeq,
                                  FALSE, chars_stripped, NULL);
          }
          if (bs == NULL && title != NULL) {
            title = MemFree (title);
          }
        }

      } else if (StringNCmp (line, "LOCUS ", 6) == 0 || StringNCmp (line, "ID ", 3) == 0) {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        GetSeqId (seqid, line, sizeof (seqid), TRUE, TRUE);

      } else if (StringNCmp (line, "ACCESSION ", 10) == 0) {

        if (StringStr (line + 10, "unknown") == NULL) {
          mayBePlainFasta = FALSE;
          mayBeAccessionList = FALSE;
          /* locus may not be unique, but accession should be, so it overrides locus */
          GetSeqId (seqid, line, sizeof (seqid), TRUE, TRUE);
        }

      } else if (StringNCmp (line, "ORIGIN", 6) == 0 || StringNCmp (line, "SQ ", 3) == 0) {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        if (! HasNoText (seqid)) {
          bs = ReadFlatFileDNA (&fc, protPtr, forceNuc, forceProt, fastaAsSimpleSeq,
                                FALSE, chars_stripped, seqid);
        }

      } else if (line [0] == '[' || line [0] == ']') {

        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);

        return NULL;

      } else {

        if (mayBePlainFasta) {
          tmp = line;
          ch = *tmp;
          while (ch != '\0') {
            if (IS_WHITESP (ch)) {
            } else if (! (IS_ALPHA (ch) || IS_DIGIT (ch) || ch == '*' || ch == '-')) {
              mayBePlainFasta = FALSE;
            } else if (protPtr != NULL) {
              ch = TO_UPPER (ch);
              if (StringChr ("EFILPQZ", ch) != NULL) {
                isProt = TRUE;
              }
            }
            tmp++;
            ch = *tmp;
          }
        }
        if (mayBeAccessionList) {
          inLetters = TRUE;
          numLetters = 0;
          numDigits = 0;
          tmp = line;
          ch = *tmp;
          while (ch != '\0') {
            if (IS_WHITESP (ch)) {
            } else if (IS_ALPHA (ch)) {
              if (! inLetters) {
                mayBeAccessionList = FALSE;
                numLetters++;
              }
            } else if (IS_DIGIT (ch)) {
              inLetters = FALSE;
              numDigits++;
            } else {
              mayBeAccessionList = FALSE;
            }
            tmp++;
            ch = *tmp;
          }
          if (numLetters == 1 && numDigits == 5) {
          } else if (numLetters == 2 && numDigits == 6) {
          } else {
            mayBeAccessionList = FALSE;
          }
        }
      }

      if (bs != NULL) {
        if (fastaAsSimpleSeq) {
          ssp = ByteStoreToSimpleSeq (bs, seqid, title);
          if (ssp != NULL) {
            if (datatypeptr != NULL) {
              *datatypeptr = OBJ_FASTA;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_FASTA, (Pointer) ssp);
            }
          }
          return (Pointer) ssp;
        }

        sep = SeqEntryNew ();
        if (sep != NULL) {
          bsp = BioseqNew ();
          if (bsp != NULL) {
            sep->choice = 1;
            sep->data.ptrvalue = (Pointer) bsp;
            SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

            if (isProt) {
              bsp->mol = Seq_mol_aa;
              bsp->seq_data_type = Seq_code_ncbieaa;
            } else {
              bsp->mol = Seq_mol_na;
              bsp->seq_data_type = Seq_code_iupacna;
            }

            bsp->repr = Seq_repr_raw;
            bsp->length = 0;
            if (parseFastaSeqId) {
              bsp->id = MakeSeqID (seqid);
            } else {
              bsp->id = MakeNewProteinSeqId (NULL, NULL);
            }
            SeqMgrAddToBioseqIndex (bsp);

            bsp->seq_data = (SeqDataPtr) bs;
            bsp->length = BSLen (bs);

            BioseqPack (bsp);

            if (title != NULL) {
              vnp = CreateNewDescriptor (sep, Seq_descr_title);
              if (vnp != NULL) {
                vnp->data.ptrvalue = (Pointer) title;
                title = NULL;
              }
              bsp->descr = vnp;
            }

            if (datatypeptr != NULL) {
              *datatypeptr = OBJ_BIOSEQ;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
            }
          }
        }

        pos = FileCacheTell (&fc);
        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);

        return (Pointer) bsp;
      }

    }

    pos = FileCacheTell (&fc);
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  if (mayBePlainFasta) {

    FileCacheSetup (&fc, fp);
    FileCacheSeek (&fc, begin);
    fseek (fp, begin, SEEK_SET);
    if (fastaAsSimpleSeq) {
      bs = ReadFlatFileDNA (&fc, NULL, (Boolean) (! isProt), (Boolean) (isProt), 
                            fastaAsSimpleSeq, FALSE, chars_stripped, NULL);
      if (bs != NULL) {
        ssp = ByteStoreToSimpleSeq (bs, NULL, NULL);
        if (ssp != NULL) {
          if (datatypeptr != NULL) {
            *datatypeptr = OBJ_FASTA;
          }
          if (entityIDptr != NULL) {
            *entityIDptr = ObjMgrRegister (OBJ_FASTA, (Pointer) ssp);
          }
        }

        pos = FileCacheTell (&fc);
        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);

        return (Pointer) ssp;
      }
    }

    /*
    sep = FastaToSeqEntryEx (fp, (Boolean) (! isProt), NULL, FALSE);
    if (sep != NULL && IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (datatypeptr != NULL) {
          *datatypeptr = OBJ_BIOSEQ;
        }
        if (entityIDptr != NULL) {
          *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
        }
      }
      return (Pointer) bsp;
    }
    */

    bs = ReadFlatFileDNA (&fc, NULL, (Boolean) (! isProt), (Boolean) (isProt),
                          FALSE, FALSE, chars_stripped, NULL);
    if (bs != NULL) {

      sep = SeqEntryNew ();
      if (sep != NULL) {
        bsp = BioseqNew ();
        if (bsp != NULL) {
          sep->choice = 1;
          sep->data.ptrvalue = (Pointer) bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

          if (isProt) {
            bsp->mol = Seq_mol_aa;
            bsp->seq_data_type = Seq_code_ncbieaa;
          } else {
            bsp->mol = Seq_mol_na;
            bsp->seq_data_type = Seq_code_iupacna;
          }

          bsp->repr = Seq_repr_raw;
          bsp->length = 0;
          bsp->id = MakeUniqueSeqID (NULL);
          SeqMgrAddToBioseqIndex (bsp);

          bsp->seq_data = (SeqDataPtr) bs;
          bsp->length = BSLen (bs);

          BioseqPack (bsp);

          if (datatypeptr != NULL) {
            *datatypeptr = OBJ_BIOSEQ;
          }
          if (entityIDptr != NULL) {
            *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
          }
        }
      }

      pos = FileCacheTell (&fc);
      FileCacheSetup (&fc, fp);
      FileCacheSeek (&fc, pos);
      fseek (fp, pos, SEEK_SET);

      return (Pointer) bsp;
    }
  }

  if (mayBeAccessionList) {

    FileCacheSetup (&fc, fp);
    FileCacheSeek (&fc, begin);
    fseek (fp, begin, SEEK_SET);
    bs = ReadUidList (&fc, TRUE, TRUE);
    if (bs != NULL) {
      numLinks = BSLen (bs) / sizeof (Int4);
      if (numLinks < 1) {
        bs = BSFree (bs);
        return NULL;
      }
      proj = ProjectNew ();
      if (proj != NULL) {
        pip = ValNodeNew (NULL);
        if (pip != NULL) {
          choice = ProjectItem_nucuid;
          pip->choice = choice;
          proj->data = pip;
          BSSeek (bs, 0L, 0);
          for (j = 0; j < numLinks; j++) {
            BSRead (bs, &uid, sizeof (Int4));
            ValNodeAddInt ((ValNodePtr PNTR) &(pip->data.ptrvalue), choice, uid);
          }
        }
      }
      bs = BSFree (bs);

      if (datatypeptr != NULL) {
        *datatypeptr = OBJ_PROJECT;
      }
      if (entityIDptr != NULL) {
        *entityIDptr = ObjMgrRegister (OBJ_PROJECT, (Pointer) proj);
      }

      pos = FileCacheTell (&fc);
      FileCacheSetup (&fc, fp);
      FileCacheSeek (&fc, pos);
      fseek (fp, pos, SEEK_SET);

      return (Pointer) proj;
    }

  }

  return NULL;
}

NLM_EXTERN Pointer ReadAsnFastaOrFlatFile (FILE *fp, Uint2Ptr datatypeptr, Uint2Ptr entityIDptr,
                                           Boolean forceNuc, Boolean forceProt,
                                           Boolean parseFastaSeqId, Boolean fastaAsSimpleSeq)
{
  return ReadAsnFastaOrFlatFileEx (fp, datatypeptr, entityIDptr,
                                  forceNuc, forceProt, 
                                  parseFastaSeqId, fastaAsSimpleSeq,
                                  NULL);
}

NLM_EXTERN Pointer ReadFeatureTableFile (
  FILE *fp,
  Uint2Ptr datatypeptr,
  Uint2Ptr entityIDptr,
  Int4Ptr lineP,
  BoolPtr failP
)

{
  CharPtr        annotname;
  Int4           begin;
  FileCache      fc;
  Char           line [4096];
  SeqFeatPtr     nextsfp;
  PubdescPtr     pdp;
  Int4           pos;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap = NULL;
  SeqFeatPtr     sfp;
  Char           seqid [2048];
  CharPtr        str;

  if (failP != NULL) *failP = FALSE;

  if (fp == NULL) return NULL;

  if (datatypeptr != NULL) *datatypeptr = 0;
  if (entityIDptr != NULL) *entityIDptr = 0;

  seqid [0] = '\0';

  FileCacheSetup (&fc, fp);

  pos = FileCacheTell (&fc);
  begin = pos;
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);

  if (str == NULL) return NULL; /* already at end of file */

  while (str != NULL) {

    if (lineP != NULL) {
      (*lineP)++;
    }

    if (! HasNoText (line)) {

      if (StringNICmp (line, ">Feature", 8) == 0) {

        annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE, FALSE);
        if (! HasNoText (seqid)) {
          sap = ReadFeatureTableEx (&fc, seqid, annotname, lineP);
          if (sap != NULL && sap->type == 1) {
            sfp = (SeqFeatPtr) sap->data;
            prevsfp = (Pointer PNTR) &(sap->data);
            while (sfp != NULL) {
              nextsfp = sfp->next;
              if (sfp->data.choice == SEQFEAT_PUB) {
                pdp = (PubdescPtr) sfp->data.value.ptrvalue;
                if (pdp != NULL && pdp->pub == NULL) {
                  *(prevsfp) = sfp->next;
                  sfp->next = NULL;
                  SeqFeatFree (sfp);
                } else {
                  prevsfp = (Pointer PNTR) &(sfp->next);
                }
              } else {
                prevsfp = (Pointer PNTR) &(sfp->next);
              }
              sfp = nextsfp;
            }
            if (sap->data == NULL) {
              sap = SeqAnnotFree (sap);
            }
          }
          if (sap != NULL) {
            if (datatypeptr != NULL) {
              *datatypeptr = OBJ_SEQANNOT;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
            }

            pos = FileCacheTell (&fc);
            FileCacheSetup (&fc, fp);
            FileCacheSeek (&fc, pos);
            fseek (fp, pos, SEEK_SET);

            return (Pointer) sap;
          }
        }

      } else {
        if (failP != NULL) {
          *failP = TRUE;
        }
        return NULL;
      }
    }

    pos = FileCacheTell (&fc);
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  return NULL;
}

extern BioseqPtr ReadFastaOnly (FILE *fp, 
                              Boolean forceNuc, Boolean forceProt,
                              BoolPtr chars_stripped,
                              CharPtr lastchar)

{
  Int4           begin;
  ByteStorePtr   bs = NULL;
  BioseqPtr      bsp = NULL;
  Char           ch;
  FileCache      fc;
  Boolean        isProt = FALSE;
  Char           line [4096];
  Boolean        mayBePlainFasta = TRUE;
  Int4           pos;
  BoolPtr        protPtr;
  SeqEntryPtr    sep;
  Char           seqid [2048];
  CharPtr        str;
  CharPtr        title = NULL;
  CharPtr        tmp;
  ValNodePtr     vnp;

  if (fp == NULL) return NULL;
 
  if (lastchar != NULL) {
    *lastchar = 0;
  }

  if (forceNuc) {
    isProt = FALSE;
    protPtr = NULL;
  } else if (forceProt) {
    isProt = TRUE;
    protPtr = NULL;
  } else {
    protPtr = &isProt;
  }

  seqid [0] = '\0';

  FileCacheSetup (&fc, fp);

  pos = FileCacheTell (&fc);
  begin = pos;
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);

  if (str == NULL) return NULL; /* already at end of file */

  while (str != NULL) {

    if (! StringHasNoText (line)) {

      if (StringStr (line, "::=") != NULL) {
        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);
        return NULL;
      } else if (line [0] == '>') {
        title = NULL;
        tmp = line + 1;
        ch = *tmp;
        while (IS_WHITESP (ch)) {
          tmp++;
          ch = *tmp;
        }
        title = StringSaveNoNull (tmp);

        bs = ReadFlatFileDNA (&fc, protPtr, forceNuc, forceProt, FALSE,
                              FALSE, chars_stripped, NULL);
        if (bs == NULL && title != NULL) {
          title = MemFree (title);
        }

      } else if (StringNCmp (line, "LOCUS ", 6) == 0 
                 || StringNCmp (line, "ID ", 3) == 0
                 || StringNCmp (line, "ACCESSION ", 10) == 0
                 || StringNCmp (line, "ORIGIN", 6) == 0 
                 || StringNCmp (line, "SQ ", 3) == 0
                 || line [0] == '[' || line [0] == ']'
                 ) {
        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);
        return NULL;
      } else {
        tmp = line;
        ch = *tmp;
        while (ch != '\0') {
          if (IS_WHITESP (ch)) {
          } else if (! (IS_ALPHA (ch) || IS_DIGIT (ch) || ch == '*' || ch == '-')) {
            FileCacheSetup (&fc, fp);
            FileCacheSeek (&fc, pos);
            fseek (fp, pos, SEEK_SET);
            if (lastchar != NULL) {
              *lastchar = ch;
            }
            return NULL;
          } else if (protPtr != NULL) {
            ch = TO_UPPER (ch);
            if (StringChr ("EFILPQZ", ch) != NULL) {
              isProt = TRUE;
            }
          }
          tmp++;
          ch = *tmp;
        }
      }

      if (bs != NULL) {
        sep = SeqEntryNew ();
        if (sep != NULL) {
          bsp = BioseqNew ();
          if (bsp != NULL) {
            sep->choice = 1;
            sep->data.ptrvalue = (Pointer) bsp;
            SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

            if (isProt) {
              bsp->mol = Seq_mol_aa;
              bsp->seq_data_type = Seq_code_ncbieaa;
            } else {
              bsp->mol = Seq_mol_na;
              bsp->seq_data_type = Seq_code_iupacna;
            }

            bsp->repr = Seq_repr_raw;
            bsp->length = 0;
            bsp->id = MakeNewProteinSeqId (NULL, NULL);
            SeqMgrAddToBioseqIndex (bsp);

            bsp->seq_data = (SeqDataPtr) bs;
            bsp->length = BSLen (bs);

            BioseqPack (bsp);

            if (title != NULL) {
              vnp = CreateNewDescriptor (sep, Seq_descr_title);
              if (vnp != NULL) {
                vnp->data.ptrvalue = (Pointer) title;
                title = NULL;
              }
              bsp->descr = vnp;
            }
          }
        }

        pos = FileCacheTell (&fc);
        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);

        return bsp;
      }

    }

    pos = FileCacheTell (&fc);
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  if (mayBePlainFasta) {

    FileCacheSetup (&fc, fp);
    FileCacheSeek (&fc, begin);
    fseek (fp, begin, SEEK_SET);

    bs = ReadFlatFileDNA (&fc, NULL, (Boolean) (! isProt), (Boolean) (isProt),
                          FALSE, FALSE, chars_stripped, NULL);
    if (bs != NULL) {

      sep = SeqEntryNew ();
      if (sep != NULL) {
        bsp = BioseqNew ();
        if (bsp != NULL) {
          sep->choice = 1;
          sep->data.ptrvalue = (Pointer) bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

          if (isProt) {
            bsp->mol = Seq_mol_aa;
            bsp->seq_data_type = Seq_code_ncbieaa;
          } else {
            bsp->mol = Seq_mol_na;
            bsp->seq_data_type = Seq_code_iupacna;
          }

          bsp->repr = Seq_repr_raw;
          bsp->length = 0;
          bsp->id = MakeUniqueSeqID (NULL);
          SeqMgrAddToBioseqIndex (bsp);

          bsp->seq_data = (SeqDataPtr) bs;
          bsp->length = BSLen (bs);

          BioseqPack (bsp);
        }
      }

      pos = FileCacheTell (&fc);
      FileCacheSetup (&fc, fp);
      FileCacheSeek (&fc, pos);
      fseek (fp, pos, SEEK_SET);

      return bsp;
    }
  }

  return NULL;
}


/* ReadDeltaFasta reads a FASTA file, combining raw sequence and >?unk100 lines into
 * a delta Bioseq.  The file pointer stops at the next FASTA with a real SeqID. 
 * The contents of perr are set to TRUE if characters were stripped from the
 * FASTA other than numbers.
 */

static ValNodePtr ReadDeltaLits (FileCachePtr fcp, BoolPtr perr, CharPtr idstr)

{
  ByteStorePtr  bs = NULL;
  Char          ch;
  Uint1         choice;
  ValNodePtr    head = NULL;
  long          len;
  Char          line [1023];
  CharPtr       str, tmp;
  Int4          pos;
  Boolean       error_flag = FALSE;

  if (fcp == NULL) return NULL;

  pos = FileCacheTell (fcp);
  str = FileCacheGetString (fcp, line, sizeof (line));

  while (str != NULL) {

    if (StringDoesHaveText (line)) {
      TrimSpacesAroundString (line);

      if ((line [0] == '>' && line [1] != '?') || line [0] == '[') {
        if (bs != NULL) {
          ValNodeAddPointer (&head, 1, (Pointer) bs);
        }

        FileCacheSeek (fcp, pos);
        return head;
      }

      if (line [0] == ']') {
        ErrPostEx (SEV_ERROR, 0, 0, "Unbalanced square bracket in ReadDeltaLits");
      
        if (bs != NULL) {
          ValNodeAddPointer (&head, 1, (Pointer) bs);
        }

        return head;
      }

      if (line [0] == '>' && line [1] == '?') {
        if (bs != NULL) {
          ValNodeAddPointer (&head, 1, (Pointer) bs);
          bs = NULL;
        }

        tmp = line + 2;
        ch = *tmp;
        while (IS_WHITESP (ch)) {
          tmp++;
          ch = *tmp;
        }
        choice = 2;
        if (StringNCmp (tmp, "unk100", 6) == 0) {
          choice = 3;
          tmp += 3;
        }
        if (*tmp != '\0' && sscanf (tmp, "%ld", &len) == 1 && len > 0) {
          ValNodeAddInt (&head, choice, (Int4) len);
        } else {
          ValNodeAddInt (&head, choice, 0);
        }

      } else {
        FileCacheSeek (fcp, pos);
        error_flag = FALSE;
        bs = ReadFlatFileDNA (fcp, NULL, TRUE, FALSE, TRUE, TRUE, &error_flag, idstr);
        if (perr != NULL)
        {
          *perr |= error_flag;
        }
      }
    }

    pos = FileCacheTell (fcp);
    str = FileCacheGetString (fcp, line, sizeof (line));
  }

  if (bs != NULL) {
    ValNodeAddPointer (&head, 1, (Pointer) bs);
  }

  return head;
}

/* perrors is set to TRUE if characters other than digits had to be stripped
 * from the FASTA sequence characters.
 */
static BioseqPtr ReadDeltaSet (FileCachePtr fcp, BoolPtr perrors, CharPtr idstr)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp = NULL;
  ValNodePtr    head, vnp;
  IntFuzzPtr    ifp;
  Boolean       is_unk100;
  SeqLitPtr     slitp;

  if (fcp == NULL) return NULL;

  head = ReadDeltaLits (fcp, perrors, idstr);
  if (head == NULL) return NULL;

  if (head->next == NULL && head->choice == 1) {
    bs = (ByteStorePtr) head->data.ptrvalue;
    if (bs == NULL) return NULL;
  
    bsp = BioseqNew ();
    if (bsp == NULL) return NULL;

    bsp->repr = Seq_repr_raw;
    bsp->seq_data_type = Seq_code_iupacna;
    bsp->mol = Seq_mol_dna;

    bsp->seq_data = (SeqDataPtr) bs;
    bsp->length = BSLen (bs);

    ValNodeFree (head);

    return bsp;
  }

  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;

  bsp->repr = Seq_repr_delta;
  bsp->seq_ext_type = 4;
  bsp->mol = Seq_mol_dna;
  bsp->length = 0;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    slitp = (SeqLitPtr) MemNew (sizeof (SeqLit));
    if (slitp == NULL) continue;

    if (vnp->choice == 1) {
      bs = (ByteStorePtr) vnp->data.ptrvalue;
      if (bs == NULL) continue;

      slitp->length = BSLen (bs);
      slitp->seq_data_type = Seq_code_iupacna;
      slitp->seq_data = (SeqDataPtr) bs;

    } else if (vnp->choice == 2 || vnp->choice == 3) {
      is_unk100 = (Boolean) vnp->choice == 3;
    
      slitp->length = vnp->data.intvalue;
      if (slitp->length < 1 || is_unk100) {
        if (slitp->length < 1) {
          slitp->length = 0;
        }
        ifp = IntFuzzNew ();
        if (ifp != NULL) {
          ifp->choice = 4;
          slitp->fuzz = ifp;
        }
      }
    }

    bsp->length += slitp->length;
    ValNodeAddPointer ((ValNodePtr PNTR) &(bsp->seq_ext), (Int2) 2, (Pointer) slitp);
  }

  ValNodeFree (head);

  return bsp;
}

NLM_EXTERN BioseqPtr ReadDeltaFastaEx (FILE *fp, Uint2Ptr entityIDptr, BoolPtr chars_stripped)

{
  Int4         begin, pos;
  BioseqPtr    bsp = NULL;
  FileCache    fc;
  Char         line [4096], seqid [2048];
  CharPtr      msg = NULL, str, title = NULL, tmp;
  SeqEntryPtr  sep;

  if (fp == NULL) return NULL;
  
  if (chars_stripped != NULL)
  {
    *chars_stripped = FALSE;
  }

  if (entityIDptr != NULL) *entityIDptr = 0;

  seqid [0] = '\0';

  FileCacheSetup (&fc, fp);

  pos = FileCacheTell (&fc);
  begin = pos;
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);

  while (str != NULL) {

    if (StringDoesHaveText (line)) {
      TrimSpacesAroundString (line);

      if (StringStr (line, "::=") != NULL) {
        msg = "ReadDeltaFasta does not read ASN.1";
      } else if (StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringNCmp (line, "ACCESSION ", 10) == 0 ||
          StringNCmp (line, "ORIGIN", 6) == 0 ||
          StringNCmp (line, "SQ ", 3) == 0) {
        msg = "ReadDeltaFasta does not read flatfiles";
      } else if (StringNCmp (line, ">PubMed", 7) == 0 ||
          StringNCmp (line, ">Protein", 8) == 0 ||
          StringNCmp (line, ">Nucleotide", 11) == 0 ||
          StringNCmp (line, ">Structure", 10) == 0 ||
          StringNCmp (line, ">Genome", 7) == 0) {
        msg = "ReadDeltaFasta does not read uid lists";
      } else if (StringNICmp (line, ">Feature", 8) == 0 ||
          StringNICmp (line, ">Vector", 7) == 0 ||
          StringNICmp (line, ">Restriction", 12) == 0 ||
          StringNICmp (line, ">Assembly", 9) == 0 ||
          StringNICmp (line, ">Virtual", 8) == 0 ||
          StringNICmp (line, ">Message", 8) == 0) {
        msg = "ReadDeltaFasta does not read special lists";
      } else if (line [0] == '[') {
        msg = "ReadDeltaFasta does not read bracketed sets";
      } else if (line [0] == '>' && StringHasNoText (line + 1)) {
        msg = "ReadDeltaFasta does not read empty deflines";
      } else if (line [0] != '>') {
        msg = "ReadDeltaFasta needs a defline";
      }

      if (msg != NULL) {
        ErrPostEx (SEV_ERROR, 0, 0, "%s", msg);

        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);

        return NULL;
      }

      if (line [0] == '>') {

        title = NULL;
        tmp = StringChr (line + 1, '[');
        if (tmp != NULL) {
          if (StringStr (tmp, "[") != NULL && StringStr (tmp, "=") != NULL) {
            TrimSpacesAroundString (tmp);
            title = StringSave (tmp);
          } else {
            title = StringSaveNoNull (line + 1);
            TrimSpacesAroundString (title);
          }
        }

        tmp = GetSeqId (seqid, line + 1, sizeof (seqid), FALSE, FALSE);

        if (StringDoesHaveText (seqid)) {

          if (StringDoesHaveText (tmp)) {
            TrimSpacesAroundString (tmp);
            title = MemFree (title);
            title = StringSaveNoNull (tmp);
          }

          bsp = ReadDeltaSet (&fc, chars_stripped, seqid);

          if (bsp != NULL) {
        
            sep = SeqEntryNew ();
            if (sep != NULL) {
              sep->choice = 1;
              sep->data.ptrvalue = (Pointer) bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
            }

            if (title != NULL) {
              SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) title);
            }

            if (StringNICmp (seqid, "lcl|", 4) == 0 
                ||StringNICmp (seqid, "gnl|", 4) == 0) {
                bsp->id = SeqIdParse (seqid);
            }
            if (bsp->id == NULL) {
              bsp->id = MakeSeqID (seqid);
            }
            SeqMgrAddToBioseqIndex (bsp);

            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
            }

            pos = FileCacheTell (&fc);
            FileCacheSetup (&fc, fp);
            FileCacheSeek (&fc, pos);
            fseek (fp, pos, SEEK_SET);

            return bsp;
          }
        }

        MemFree (title);
      }
    }

    pos = FileCacheTell (&fc);
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, begin);
  fseek (fp, begin, SEEK_SET);

  return NULL;
}

NLM_EXTERN BioseqPtr ReadDeltaFasta (FILE *fp, Uint2Ptr entityIDptr)

{
  Boolean chars_stripped = FALSE;
  
  return ReadDeltaFastaEx (fp, entityIDptr, &chars_stripped);
}

NLM_EXTERN BioseqPtr ReadDeltaFastaWithEmptyDefline (FILE *fp, Uint2Ptr entityIDptr, BoolPtr chars_stripped)

{
  Int4         begin, pos;
  BioseqPtr    bsp = NULL;
  FileCache    fc;
  Char         line [4096], seqid [2048];
  SeqEntryPtr  sep;
  CharPtr      str;

  if (fp == NULL) return NULL;
  
  if (chars_stripped != NULL)
  {
    *chars_stripped = FALSE;
  }

  if (entityIDptr != NULL) *entityIDptr = 0;

  seqid [0] = '\0';

  FileCacheSetup (&fc, fp);

  pos = FileCacheTell (&fc);
  begin = pos;
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  if (str != NULL && StringDoesHaveText (line))
  {
    TrimSpacesAroundString (line);
    if (line [0] == '>' && line [1] == 0)
    {
      bsp = ReadDeltaSet (&fc, chars_stripped, NULL);

      if (bsp != NULL) {
      
        sep = SeqEntryNew ();
        if (sep != NULL) {
          sep->choice = 1;
          sep->data.ptrvalue = (Pointer) bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        }

        bsp->id = MakeUniqueSeqID ("delta_");
        SeqMgrAddToBioseqIndex (bsp);

        if (entityIDptr != NULL) {
          *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
        }

        pos = FileCacheTell (&fc);
        FileCacheSetup (&fc, fp);
        FileCacheSeek (&fc, pos);
        fseek (fp, pos, SEEK_SET);

        return bsp;
      }
    }
  }

  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, begin);
  fseek (fp, begin, SEEK_SET);

  return NULL;
}

/* general purpose text finite state machine */
/* based on Practical Algorithms for Programmers by Binstock and Rex */

typedef struct fsagoto {
  Char             ch;
  Int4             newstate;
  struct fsagoto * next;
} GotoItem, PNTR GotoPtr;

typedef struct fsastate {
  GotoPtr       transition;
  ValNodePtr    matchfound;
  Int4          onfailure;
} StateItem, PNTR StatePtr;

#define FAIL_STATE -1

static StatePtr GetState (
  StatePtr PNTR stateTable,
  Int4 state
)

{
  StatePtr  sp;

  if (state < 0) return NULL;

  sp = stateTable [state];
  if (sp == NULL) {
    sp = (StatePtr) MemNew (sizeof (StateItem));
    stateTable [state] = sp;
  }

  return sp;
}

static Int4 GotoState (StatePtr PNTR stateTable, Int4 state,
                       Char ch, Boolean zeroFailureReturnsZero)

{
  GotoPtr   gp;
  StatePtr  sp;

  sp = GetState (stateTable, state);
  if (sp == NULL) return 0;

  for (gp = sp->transition; gp != NULL; gp = gp->next) {
    if (gp->ch == ch) return gp->newstate;
  }

  if (state == 0 && zeroFailureReturnsZero) return 0;

  return FAIL_STATE;
}

/*
#define FailState(stateTable,state) stateTable [state].onfailure
*/

static Int4 FailState (
  StatePtr PNTR stateTable,
  Int4 state
)

{
  StatePtr  sp;

  sp = GetState (stateTable, state);
  if (sp == NULL) return 0;

  return sp->onfailure;
}

static void AddTransition (StatePtr PNTR stateTable, Int4 oldState,
                           Char ch, Int4 newState)

{
  GotoPtr   gp;
  GotoPtr   prev;
  StatePtr  sp;

  gp = (GotoPtr) MemNew (sizeof (GotoItem));
  if (gp == NULL) return;

  gp->ch = ch;
  gp->newstate = newState;

  sp = GetState (stateTable, oldState);
  if (sp == NULL) return;

  prev = sp->transition;
  if (prev == NULL) {
    sp->transition = gp;
  } else {
    while (prev->next != NULL) {
      prev = prev->next;
    }
    prev->next = gp;
  }
}

static void AddOutput (StatePtr PNTR stateTable, Int4 state, CharPtr word)

{
  StatePtr    sp;
  ValNodePtr  vnp;

  sp = GetState (stateTable, state);
  if (sp == NULL) return;

  for (vnp = sp->matchfound; vnp != NULL; vnp = vnp->next) {
    if (StringCmp (word, (CharPtr) vnp->data.ptrvalue) == 0) return;
  }

  ValNodeCopyStr (&(sp->matchfound), 0, word);
}

static Int4 EnterWord (StatePtr PNTR stateTable, CharPtr word,
                       Int4 highState, Int4 maxState)

{
  Char     ch;
  Int4     next;
  CharPtr  ptr;
  Int4     state;

  state = 0;
  next = 0;

  /* try to overlay beginning of word onto existing table */

  for (ptr = word, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
    next = GotoState (stateTable, state, ch, FALSE);
    if (next == FAIL_STATE) break;
    state = next;
  }

  /* now create new states for remaining characters in word */

  for ( ; ch != '\0'; ptr++, ch = *ptr) {
    highState++;
    if (highState >= maxState) return highState;

    AddTransition (stateTable, state, ch, highState);
    state = highState;
  }

  /* at end of word record match information */

  AddOutput (stateTable, state, word);

  return highState;
}

static void QueueAdd (Int4Ptr queue, Int4 qbeg, Int4 val)

{
  Int4  q;

  q = queue [qbeg];
  if (q == 0) {
    queue [qbeg] = val;
  } else {
    for ( ; queue [q] != 0; q = queue [q]) continue;
    queue [q] = val;
  }
  queue [val] = 0;
}

static void FindFail (StatePtr PNTR stateTable, Int4 state,
                      Int4 newState, Char ch)

{
  Int4        next;
  StatePtr    sp;
  ValNodePtr  vnp;

  /* traverse existing failure path */

  next = GotoState (stateTable, state, ch, TRUE);

  while ((next = GotoState (stateTable, state, ch, TRUE)) == FAIL_STATE) {
    state = FailState (stateTable, state);
  }

  /* add new failure state */

  sp = GetState (stateTable, newState);
  if (sp == NULL) return;

  sp->onfailure = next;

  /* add matches of substring at new state */

  sp = GetState (stateTable, next);
  if (sp == NULL) return;

  for (vnp = sp->matchfound; vnp != NULL; vnp = vnp->next) {
    AddOutput (stateTable, newState, (CharPtr) vnp->data.ptrvalue);
  }
}

static void ComputeFail (StatePtr PNTR stateTable, Int4Ptr queue, Int4 highState)

{
  GotoPtr   gp;
  Int4      qbeg, r, s, state;
  StatePtr  sp;

  qbeg = 0;
  queue [0] = 0;

  /* queue up states reached directly from state 0 (depth 1) */

  sp = GetState (stateTable, 0);
  if (sp == NULL) return;

  for (gp = sp->transition; gp != NULL; gp = gp->next) {
    s = gp->newstate;

    sp = GetState (stateTable, s);
    if (sp == NULL) return;

    sp->onfailure = 0;
    QueueAdd (queue, qbeg, s);
  }

  while (queue [qbeg] != 0) {
    r = queue [qbeg];
    qbeg = r;

    /* depth 1 states beget depth 2 states, etc. */

    sp = GetState (stateTable, r);
    if (sp == NULL) return;

    for (gp = sp->transition; gp != NULL; gp = gp->next) {
      s = gp->newstate;
      QueueAdd (queue, qbeg, s);

      /*
         State   Substring   Transitions   Failure
           2       st          a ->   3       6
           3       sta         l ->   4
           6       t           a ->   7       0
           7       ta          p ->   8

         For example, r = 2 (st), if 'a' would go to s = 3 (sta).
         From previous computation, 2 (st) fails to 6 (t).
         Thus, check state 6 (t) for any transitions using 'a'.
         Since 6 (t) 'a' -> 7 (ta), therefore set fail [3] -> 7.
      */

      state = FailState (stateTable, r);
      FindFail (stateTable, state, s, gp->ch);
    }
  }
}

typedef struct TextFsa {
  StatePtr PNTR  stateTable;
  ValNodePtr     siteList;
  Int4           highState;
  Int4           numWords;
  Int4           longestWord;
  Boolean        primed;
} TextFsaData;

static void PrimeStateTable (TextFsaPtr tbl)

{
  Int4           highState;
  Int4           maxState;
  Int4Ptr        queue;
  StatePtr PNTR  stateTable;
  ValNodePtr     vnp;
  CharPtr        word;

  if (tbl == NULL || tbl->siteList == NULL || tbl->primed) return;

  for (maxState = 1, vnp = tbl->siteList; vnp != NULL; vnp = vnp->next) {
    word = (CharPtr) vnp->data.ptrvalue;
    maxState += StringLen (word);
  }

  maxState++;

  stateTable = (StatePtr PNTR) MemNew (sizeof (StatePtr) * (size_t) maxState);
  queue = (Int4Ptr) MemNew (sizeof (Int4) * maxState);

  if (stateTable == NULL || queue == NULL) {
    MemFree (stateTable);
    MemFree (queue);
    Message (MSG_POST, "FiniteStateSearch unable to allocate buffers");
    return;
  }

  for (highState = 0, vnp = tbl->siteList; vnp != NULL; vnp = vnp->next) {
    word = (CharPtr) vnp->data.ptrvalue;
    highState = EnterWord (stateTable, word, highState, maxState);
  }

  if (highState >= maxState) {
    ErrPostEx (SEV_ERROR, 0, 0, "FiniteStateSearch cannot handle more than %d states", (int) highState);
  }

  ComputeFail (stateTable, queue, highState);

  MemFree (queue);

  tbl->stateTable = stateTable;
  tbl->highState = highState;
  tbl->primed = TRUE;
}

NLM_EXTERN TextFsaPtr TextFsaNew (void)

{
  TextFsaPtr  tbl;

  tbl = (TextFsaPtr) MemNew (sizeof (TextFsaData));
  if (tbl == NULL) return NULL;
  tbl->stateTable = NULL;
  tbl->siteList = NULL;
  tbl->primed = FALSE;
  return tbl;
}

NLM_EXTERN void TextFsaAdd (TextFsaPtr tbl, CharPtr word)

{
  Int4  len;

  if (tbl == NULL) return;
  len = (Int4) StringLen (word);
  if (len < 1) return;
  ValNodeCopyStr (&(tbl->siteList), 0, word);
  (tbl->numWords)++;
  if (len > tbl->longestWord) {
    tbl->longestWord = len;
  }
}

NLM_EXTERN Int4 TextFsaNext (TextFsaPtr tbl, Int4 currState,
                             Char ch, ValNodePtr PNTR matches)

{
  Int4           next;
  StatePtr       sp;
  StatePtr PNTR  stateTable;

  if (matches != NULL) {
    *matches = NULL;
  }
  if (tbl == NULL) return 0;
  if (! tbl->primed) {
    PrimeStateTable (tbl);
  }
  stateTable = tbl->stateTable;
  if (stateTable == NULL) return 0;

  while ((next = GotoState (stateTable, currState, ch, TRUE)) == FAIL_STATE) {
    currState = FailState (stateTable, currState);
  }

  if (matches != NULL) {

    sp = GetState (stateTable, next);
    if (sp == NULL) return next;

    *matches = sp->matchfound;
  }

  return next;
}

NLM_EXTERN Boolean TextFsaGetStats (
  TextFsaPtr tbl,
  Int4Ptr highStateP,
  Int4Ptr numWordsP,
  Int4Ptr longestWordP
)

{
  if (tbl == NULL) return FALSE;
  if (highStateP != NULL) {
    *highStateP = tbl->highState;
  }
  if (numWordsP != NULL) {
    *numWordsP = tbl->numWords;
  }
  if (longestWordP != NULL) {
    *longestWordP = tbl->longestWord;
  }
  return TRUE;
}

NLM_EXTERN TextFsaPtr TextFsaFree (TextFsaPtr tbl)

{
  GotoPtr        gp;
  Int4           highState;
  GotoPtr        nxtgp;
  StatePtr       sp;
  Int4           state;
  StatePtr PNTR  stateTable;

  if (tbl == NULL) return NULL;

  stateTable = tbl->stateTable;
  if (stateTable != NULL) {
    highState = tbl->highState;

    for (state = 0; state <= highState; state++) {
      sp = stateTable [state];
      if (sp == NULL) continue;

      gp = sp->transition;
      while (gp != NULL) {
        nxtgp = gp->next;
        MemFree (gp);
        gp = nxtgp;
      }

      sp->matchfound = ValNodeFreeData (sp->matchfound);
      sp = MemFree (sp);
      stateTable[state] = NULL;
    }

    stateTable = MemFree (stateTable);
  }

  tbl->siteList = ValNodeFreeData (tbl->siteList);

  return MemFree (tbl);
}

/* sequence quality exchange */

typedef struct gphgetdata {
  ValNodePtr  vnp;
  BioseqPtr   bsp;
} GphGetData, PNTR GphGetPtr;

typedef struct gphitem {
  SeqGraphPtr  sgp;
  Int4         left;
  Int4         right;
  Int2         index;
} GphItem, PNTR GphItemPtr;

static void GetGraphsProc (SeqGraphPtr sgp, Pointer userdata)

{
  GphGetPtr   ggp;
  GphItemPtr  gip;

  ggp = (GphGetPtr) userdata;
  if (ggp == NULL || sgp == NULL) return;
  /* only phrap or gap4 currently allowed */
  if (StringICmp (sgp->title, "Phrap Quality") == 0 ||
      StringICmp (sgp->title, "Phred Quality") == 0 ||
      StringICmp (sgp->title, "Gap4") == 0) {
    /* data type must be bytes */
    if (sgp->flags[2] == 3) {
      if (SeqIdIn (SeqLocId (sgp->loc), ggp->bsp->id)) {
        gip = (GphItemPtr) MemNew (sizeof (GphItem));
        if (gip == NULL) return;
        gip->sgp = sgp;
        gip->left = GetOffsetInBioseq (sgp->loc, ggp->bsp, SEQLOC_LEFT_END);
        gip->right = GetOffsetInBioseq (sgp->loc, ggp->bsp, SEQLOC_RIGHT_END);
        ValNodeAddPointer (&(ggp->vnp), 0, (Pointer) gip);
      }
    }
  }
}

static int LIBCALLBACK SortSeqGraphProc (VoidPtr ptr1, VoidPtr ptr2)

{
  GphItemPtr  gip1, gip2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gip1 = (GphItemPtr) vnp1->data.ptrvalue;
  gip2 = (GphItemPtr) vnp2->data.ptrvalue;
  if (gip1 == NULL || gip2 == NULL) return 0;
  if (gip1->left > gip2->left) {
    return 1;
  } else if (gip1->left < gip2->left) {
    return -1;
  } else if (gip1->right > gip2->right) {
    return -1;
  } else if (gip2->right < gip2->right) {
    return 1;
  }
  return 0;
}

/* gets valnode list of sorted graphs in GphItem structures */

static ValNodePtr GetSeqGraphsOnBioseq (BioseqPtr bsp)

{
  GphGetData  ggd;
  GphItemPtr  gip;
  Int2        index;
  ValNodePtr  vnp;

  ggd.vnp = NULL;
  ggd.bsp = bsp;
  VisitGraphsOnBsp (bsp, (Pointer) &ggd, GetGraphsProc);
  for (vnp = ggd.vnp, index = 1; vnp != NULL; vnp = vnp->next, index++) {
    gip = (GphItemPtr) vnp->data.ptrvalue;
    if (gip != NULL) {
      gip->index = index;
    }
  }
  ggd.vnp = ValNodeSort (ggd.vnp, SortSeqGraphProc);
  return ggd.vnp;
}

static void PrintQualProc (CharPtr buf, Uint4 buflen, Pointer userdata)

{
  FILE  *fp;

  fp = (FILE*) userdata;
  fprintf (fp, "%s", buf);
}

NLM_EXTERN void PrintQualityScores (BioseqPtr bsp, FILE *fp)

{
  PrintQualityScoresToBuffer (bsp, TRUE, (Pointer) fp, PrintQualProc);
}

NLM_EXTERN void PrintQualityScoresToBuffer (BioseqPtr bsp, Boolean gapIsZero, Pointer userdata, QualityWriteFunc callback)

{
  ByteStorePtr  bs;
  Char          id [41], buf [84], tmp [16];
  Int4          curpos = 0, c80, i, len = 0, min = INT4_MAX, max = INT4_MIN;
  Uint2         entityID;
  Int2          gap;
  GphItemPtr    gip;
  ValNodePtr    head, vnp;
  SeqGraphPtr   sgp;
  SeqIdPtr      sip, sip2;
  CharPtr       title = NULL, ptr;
  Int2          val;

  if (bsp == NULL || callback == NULL) return;
  entityID = ObjMgrGetEntityIDForPointer (bsp);
  head = GetSeqGraphsOnBioseq (bsp);

  /* skip bioseqs with no quality graphs */

  if (head == NULL) return;

  /* find accession */

  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  if (sip->choice == SEQID_GI) {
    sip2 = GetSeqIdForGI (sip->data.intvalue);
    if (sip2 != NULL) {
      sip = sip2;
    }
  }
  SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);

  if (gapIsZero) {
    gap = 0;
  } else {
    gap = -1;
  }

  /* get min, max, title, but currently won't use cumulative length */

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gip = (GphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL) continue;
    sgp = gip->sgp;
    min = MIN ((Int4) min, (Int4) sgp->min.intvalue);
    max = MAX ((Int4) max, (Int4) sgp->max.intvalue);
    len += sgp->numval;
    if (title == NULL) {
      title = sgp->title;
    }
  }
  if (title == NULL) {
    title = "?";
  }
  len = bsp->length; /* report full length of bioseq */
  if (min == INT4_MAX) {
    min = 0;
  }
  if (max == INT4_MIN) {
    max = 0;
  }
  sprintf (buf, ">%s %s (Length:%ld, Min: %ld, Max: %ld)\n", id, title,
           (long) len, (long) min, (long) max);
  callback (buf, sizeof (buf), userdata);

  c80 = 0;
  ptr = buf;
  buf [0] = '\0';

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gip = (GphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL) continue;
    sgp = gip->sgp;

    /* expand gaps by padding with 0s (optionally -1) */

    while (curpos < gip->left) {
      if (c80 == 20) {
        c80 = 0;
        ptr = StringMove (ptr, "\n");
        callback (buf, sizeof (buf), userdata);
        ptr = buf;
        buf [0] = '\0';
      }
      sprintf (tmp, "%3d", (int) gap);
      ptr = StringMove (ptr, tmp);
      curpos++;
      c80++;
    }

    /* now at proper position, write actual scores */

    bs = (ByteStorePtr) sgp->values;
    BSSeek (bs, 0, SEEK_SET);
    for (i = 0; i < sgp->numval; i++) {
      val = (Int2) BSGetByte (bs);
      if (c80 == 20) {
        c80 = 0;
        ptr = StringMove (ptr, "\n");
        callback (buf, sizeof (buf), userdata);
        ptr = buf;
        buf [0] = '\0';
      }
      if (val < 100) {
        sprintf (tmp, "%3d", (int) val);
      } else {
        sprintf (tmp, "%4d", (int) val);
      }
      ptr = StringMove (ptr, tmp);
      curpos++;
      c80++;
    }
  }

  /* expand any remaining space at end by padding with 0s (optionally -1) */

  while (curpos < bsp->length) {
    if (c80 == 20) {
      c80 = 0;
      ptr = StringMove (ptr, "\n");
      callback (buf, sizeof (buf), userdata);
      ptr = buf;
      buf [0] = '\0';
    }
    sprintf (tmp, "%3d", (int) gap);
    ptr = StringMove (ptr, tmp);
    curpos++;
    c80++;
  }

  ptr = StringMove (ptr, "\n");
  callback (buf, sizeof (buf), userdata);

  ValNodeFreeData (head);
}


NLM_EXTERN void TrimSeqGraph (SeqGraphPtr sgp, Int4 num_to_trim, Boolean from_left)
{
  FloatHiPtr   new_flvalues = NULL, old_flvalues;
  Int4Ptr      new_intvalues = NULL, old_intvalues;
  ByteStorePtr new_bytevalues = NULL, old_bytevalues;
  Int4         new_len;
  Int4         start_pos;
  FloatHi      fhmax = 0.0, fhmin = 0.0;
  Int4         intmax = 0, intmin = 0;
  Int2         bs_max = 0, bs_min = 0;
  Int4         new_pos, old_pos;
  Int2         val;
  Int4         loc_stop;
  Boolean      changed = FALSE;
  
  if (sgp == NULL || num_to_trim < 1)
  {
    return;
  }
  
  new_len = sgp->numval - num_to_trim;
  if (from_left)
  {
    start_pos = num_to_trim;
  }
  else
  {
    start_pos = 0;
  }
  
  if (sgp->flags[2] == 1)
  {
    new_flvalues = (FloatHiPtr) MemNew (new_len * sizeof (FloatHi));
    old_flvalues = (FloatHiPtr) sgp->values;
    new_pos = 0;
    old_pos = start_pos;
    while (old_pos < sgp->numval && new_pos < new_len)
    {
      new_flvalues [new_pos] = old_flvalues[start_pos];
      if (old_pos == start_pos)
      {
        fhmax = new_flvalues[new_pos];
        fhmin = new_flvalues[new_pos];
      }
      else
      {
        if (fhmax < new_flvalues[new_pos])
        {
          fhmax = new_flvalues[new_pos];
        }
        
        if (fhmin > new_flvalues[new_pos])
        {
          fhmin = new_flvalues[new_pos];
        }
      }
      new_pos++;
      old_pos++;
    }
    old_flvalues = MemFree (old_flvalues);
    sgp->values = new_flvalues;
    sgp->numval = new_len;
    sgp->max.realvalue = fhmax;
    sgp->min.realvalue = fhmin;
    changed = TRUE;
  }
  else if (sgp->flags[2] == 2)
  {
    new_intvalues = (Int4Ptr) MemNew (new_len * sizeof (FloatHi));
    old_intvalues = (Int4Ptr) sgp->values;
    new_pos = 0;
    old_pos = start_pos;
    while (old_pos < sgp->numval && new_pos < new_len)
    {
      new_intvalues [new_pos] = old_intvalues[start_pos];
      if (old_pos == start_pos)
      {
        intmax = new_intvalues[new_pos];
        intmin = new_intvalues[new_pos];
      }
      else
      {
        if (intmax < new_intvalues[new_pos])
        {
          intmax = new_intvalues[new_pos];
        }
        
        if (intmin > new_intvalues[new_pos])
        {
          intmin = new_intvalues[new_pos];
        }
      }
      new_pos++;
      old_pos++;
    }
    old_intvalues = MemFree (old_intvalues);
    sgp->values = new_intvalues;
    sgp->numval = new_len;
    sgp->max.intvalue = intmax;
    sgp->min.intvalue = intmin;
    changed = TRUE;
  }
  else if (sgp->flags[2] == 3)
  {
    new_bytevalues = BSNew(new_len + 1);
    old_bytevalues = (ByteStorePtr) sgp->values;
    new_pos = 0;
    old_pos = start_pos;
    while (old_pos < sgp->numval && new_pos < new_len)
    {
      BSSeek (old_bytevalues, old_pos, SEEK_SET);
      BSSeek (new_bytevalues, new_pos, SEEK_SET);
      val = (Int2) BSGetByte (old_bytevalues);
      BSPutByte (new_bytevalues, val);

      if (old_pos == start_pos)
      {
        bs_max = val;
        bs_min = val;
      }
      else
      {
        if (bs_max < val)
        {
          bs_max = val;
        }
        
        if (bs_min > val)
        {
          bs_min = val;
        }
      }
      new_pos++;
      old_pos++;
    }
    old_bytevalues = BSFree (old_bytevalues);
    sgp->values = new_bytevalues;
    sgp->numval = new_len;
    sgp->max.intvalue = bs_max;
    sgp->min.intvalue = bs_min;
    changed = TRUE;
  }
  if (changed) 
  {
    loc_stop = SeqLocStop (sgp->loc);  
    sgp->loc = SeqLocDelete (sgp->loc, SeqLocId (sgp->loc), 
                             loc_stop - num_to_trim + 1, 
                             loc_stop, FALSE, &changed);
  }
}


NLM_EXTERN void TrimQualityScores (BioseqPtr bsp, Int4 num_to_trim, Boolean from_left)
{
  ValNodePtr    qual_scores, vnp;
  GphItemPtr    gip;

  if (bsp == NULL) return;
  qual_scores = GetSeqGraphsOnBioseq (bsp);
  for (vnp = qual_scores; vnp != NULL; vnp = vnp->next)
  {
    gip = (GphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL) continue;
    TrimSeqGraph (gip->sgp, num_to_trim, from_left);
  }
  
}


NLM_EXTERN void ReverseSeqGraph (SeqGraphPtr sgp)
{
  FloatHiPtr   flvalues;
  Int4Ptr      intvalues;
  ByteStorePtr new_bytevalues = NULL, old_bytevalues;
  Int4         pos, mid, antipos;
  FloatHi      fswap;
  Int4         iswap;
  Int2         val;
  Int4         loc_start, loc_stop, diff_left, diff_right;
  Boolean      changed = FALSE;
  BioseqPtr    bsp;
  
  if (sgp == NULL)
  {
    return;
  }
    
  if (sgp->flags[2] == 1)
  {
    flvalues = (FloatHiPtr) sgp->values;
    mid = sgp->numval / 2 - 1;
    for (pos = 0; mid; pos++) {
      fswap = flvalues[pos];
      flvalues[pos] = flvalues[sgp->numval - pos - 1];
      flvalues[sgp->numval - pos] = fswap;
    }
    changed = TRUE;
  }
  else if (sgp->flags[2] == 2)
  {
    intvalues = (Int4Ptr) sgp->values;
    mid = sgp->numval / 2 - 1;
    for (pos = 0; mid; pos++) {
      iswap = intvalues[pos];
      intvalues[pos] = intvalues[sgp->numval - pos - 1];
      intvalues[sgp->numval - pos] = iswap;
    }
    changed = TRUE;
  }
  else if (sgp->flags[2] == 3)
  {
    new_bytevalues = BSNew(sgp->numval + 1);
    old_bytevalues = (ByteStorePtr) sgp->values;
    pos = 0;
    antipos = sgp->numval - 1;
    while (pos < sgp->numval)
    {
      BSSeek (old_bytevalues, antipos, SEEK_SET);
      BSSeek (new_bytevalues, pos, SEEK_SET);
      val = (Int2) BSGetByte (old_bytevalues);
      BSPutByte (new_bytevalues, val);
      BSSeek (new_bytevalues, pos, SEEK_SET);
      val = (Int2) BSGetByte (new_bytevalues);
      pos++;
      antipos--;
    }
    old_bytevalues = BSFree (old_bytevalues);
    sgp->values = new_bytevalues;
    changed = TRUE;
  }
  if (changed) 
  {
    bsp = BioseqLockById (SeqLocId (sgp->loc));
    if (bsp != NULL) {
      loc_start = SeqLocStart (sgp->loc);
      loc_stop = SeqLocStop (sgp->loc);
      if (loc_start < loc_stop) {
        diff_left = loc_start;
        diff_right = bsp->length - loc_stop - 1;
      } else {
        diff_left = loc_stop;
        diff_right = bsp->length - loc_start - 1;
      }
      if (diff_right != diff_left) {
        SeqEdAdjustFeatureInterval (sgp->loc, diff_right - diff_left, eSlide, 0, bsp);
      }
      BioseqUnlock (bsp);
    }
  }
}


NLM_EXTERN void ReverseQualityScores (BioseqPtr bsp)
{
  ValNodePtr    qual_scores, vnp;
  GphItemPtr    gip;

  if (bsp == NULL) return;
  qual_scores = GetSeqGraphsOnBioseq (bsp);
  for (vnp = qual_scores; vnp != NULL; vnp = vnp->next)
  {
    gip = (GphItemPtr) vnp->data.ptrvalue;
    if (gip == NULL) continue;
    ReverseSeqGraph (gip->sgp);
  }
  
}


NLM_EXTERN BytePtr GetScoresbySeqId (SeqIdPtr sip, Int4Ptr bsplength)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Int4          curpos = 0, i;
  Uint2         entityID;
  GphItemPtr    gip;
  ValNodePtr    head, vnp;
  Int4          len;
  SeqGraphPtr   sgp;
  BytePtr       str = NULL;

  if (bsplength != NULL) {
    *bsplength = 0;
  }
  if (sip == NULL) return NULL;

  bsp = BioseqLockById (sip);
  if (bsp == NULL) return NULL;

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  head = GetSeqGraphsOnBioseq (bsp);

  if (head != NULL && ISA_na (bsp->mol) && bsp->length < MAXALLOC) {
    str = MemNew (sizeof (Byte) * (bsp->length + 2));
    if (str != NULL) {

      len = bsp->length;
      if (bsplength != NULL) {
        *bsplength = len;
      }

      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        gip = (GphItemPtr) vnp->data.ptrvalue;
        if (gip == NULL) continue;
        sgp = gip->sgp;

        /* expand gaps by padding with 0s (now 255) */

        while (curpos < gip->left && curpos < len) {
          str [curpos] = 255;
          curpos++;
        }

        /* now at proper position, write actual scores */

        bs = (ByteStorePtr) sgp->values;
        BSSeek (bs, 0, SEEK_SET);
        for (i = 0; i < sgp->numval && curpos < len; i++) {
          str [curpos] = (Byte) BSGetByte (bs);
          curpos++;
        }
      }

      /* expand any remaining space at end by padding with 0s (now 255) */

      while (curpos < len) {
        str [curpos] = 255;
        curpos++;
      }

    }
  }

  ValNodeFreeData (head);
  BioseqUnlock (bsp);
  return str;
}

NLM_EXTERN BytePtr GetScoresbyAccessionDotVersion (CharPtr accession, Int4Ptr bsplength)

{
  BytePtr   bs;
  SeqIdPtr  sip;

  if (bsplength != NULL) {
    *bsplength = 0;
  }
  if (StringHasNoText (accession)) return NULL;
  sip = SeqIdFromAccessionDotVersion (accession);
  if (sip == NULL) return NULL;

  bs = GetScoresbySeqId (sip, bsplength);
  sip = SeqIdFree (sip);
  return bs;
}

static void PrintAScore (
  FILE* fp,
  Int2 val,
  Int2Ptr linepos
)

{
  if (*linepos >= 20) {
    fprintf (fp, "\n");
    *linepos = 0;
  }
  if (val == 255) {
    val = -1;
  }
  fprintf (fp, "%3d", (int) val);
  (*linepos)++;
}

NLM_EXTERN void PrintQualityScoresForContig (
  BioseqPtr bsp,
  Boolean gapIsZero,
  FILE* fp
)

{
  Char         accn [41];
  BytePtr      bp;
  DeltaSeqPtr  dsp;
  Int2         gap;
  Int4         i;
  Int4         len;
  Int2         linepos = 0;
  SeqIdPtr     sip, sip2;
  SeqLitPtr    slitp;
  SeqLocPtr    slp;
  Int4         tstart, tstop;

  if (bsp == NULL || fp == NULL) return;
  if (bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL) return;

  /* find accession */

  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return;
  if (sip->choice == SEQID_GI) {
    sip2 = GetSeqIdForGI (sip->data.intvalue);
    if (sip2 != NULL) {
      sip = sip2;
    }
  }
  SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn) - 1);
  fprintf (fp, ">%s\n", accn);

  if (gapIsZero) {
    gap = 0;
  } else {
    gap = -1;
  }

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {

      slp = (SeqLocPtr) dsp->data.ptrvalue;
      if (slp == NULL || slp->choice == SEQLOC_NULL) continue;

      sip = SeqLocId (slp);
      if (sip == NULL) continue;

      /*
      if (sip->choice == SEQID_GI) {
        gi = sip->data.intvalue;
        accn [0] = '\0';
      } else {
        SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn) - 1);
        gi = 0;
      }
      */
      bp = GetScoresbySeqId (sip, &len);
      if (bp == NULL) {
        len = SeqLocLen (slp);
        for (i = 0; i < len; i++) {
          PrintAScore (fp, gap, &linepos);
        }
        continue;
      }

      tstart = SeqLocStart (slp);
      tstop = SeqLocStop (slp);

      len = tstop - tstart + 1;
      if (len == SeqLocLen (slp)) {
        if (SeqLocStrand (slp) == Seq_strand_minus) {
          for (i = tstop; i >= tstart; i--) {
            PrintAScore (fp, bp [i], &linepos);
          }
        } else {
          for (i = tstart; i <= tstop; i++) {
            PrintAScore (fp, bp [i], &linepos);
          }
        }
      }

      MemFree (bp);

    } else if (dsp->choice == 2) {

      slitp = (SeqLitPtr) dsp->data.ptrvalue;
      if (slitp == NULL /* || slitp->seq_data != NULL */) continue;
      for (i = 0; i < slitp->length; i++) {
        PrintAScore (fp, gap, &linepos);
      }
    }
  }

  fprintf (fp, "\n");
}

typedef struct phrapdata {
  BioseqPtr  bsp;
  Int4       length;
  BytePtr    scores;
} PhrapData, PNTR PhrapDataPtr;

NLM_EXTERN SeqAnnotPtr PhrapGraphForContig (
  BioseqPtr bsp
)

{
  ByteStorePtr  bs;
  BytePtr       bp, str = NULL, ptr, tmp;
  Byte          by;
  DeltaSeqPtr   dsp;
  Int4          i, len, tstart, tstop;
  Int2          max = INT2_MIN;
  Int2          min = INT2_MAX;
  BioseqPtr     pbsp;
  PhrapDataPtr  pdp;
  ValNodePtr    phplist = NULL, vnp;
  SeqAnnotPtr   sap = NULL;
  SeqGraphPtr   sgp, lastsgp = NULL;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLitPtr     slitp;
  SeqLocPtr     slp;

  if (bsp == NULL) return NULL;
  if (bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL) return NULL;
  if ((! ISA_na (bsp->mol)) || bsp->length >= MAXALLOC) return NULL;

  str = MemNew (sizeof (Byte) * (bsp->length + 2));
  if (str == NULL) return NULL;

  /* initialize every byte to 255 (gap) so only real regions will be kept */

  for (ptr = str, i = 0; i < bsp->length; ptr++, i++) {
    *ptr = 255;
  }

  ptr = str;

  /* lock all components once, get uniqued list of component Bioseqs and scores */

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {

      slp = (SeqLocPtr) dsp->data.ptrvalue;
      if (slp == NULL || slp->choice == SEQLOC_NULL) continue;

      sip = SeqLocId (slp);
      if (sip == NULL) continue;

      pbsp = BioseqLockById (sip);
      if (pbsp == NULL) continue;

      for (vnp = phplist; vnp != NULL; vnp = vnp->next) {
        pdp = (PhrapDataPtr) vnp->data.ptrvalue;
        if (SeqIdIn (sip, pdp->bsp->id)) break;
      }
      if (vnp == NULL) {
        pdp = (PhrapDataPtr) MemNew (sizeof (PhrapData));
        if (pdp == NULL) continue;
        pdp->bsp = pbsp;
        pdp->scores = GetScoresbySeqId (pbsp->id, &(pdp->length));
        ValNodeAddPointer (&phplist, 0, (Pointer) pdp);
      } else {
        BioseqUnlock (pbsp);
      }
    }
  }

  /* build master byte array of scores */

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {

      slp = (SeqLocPtr) dsp->data.ptrvalue;
      if (slp == NULL || slp->choice == SEQLOC_NULL) continue;

      sip = SeqLocId (slp);
      if (sip == NULL) continue;

      bp = NULL;
      for (vnp = phplist; vnp != NULL; vnp = vnp->next) {
        pdp = (PhrapDataPtr) vnp->data.ptrvalue;
        if (SeqIdIn (sip, pdp->bsp->id)) {
          bp = pdp->scores;
          break;
        }
      }
      if (bp == NULL) {
        len = SeqLocLen (slp);
        for (i = 0; i < len; i++) {
          *ptr = 255;
          ptr++;
        }
        continue;
      }

      tstart = SeqLocStart (slp);
      tstop = SeqLocStop (slp);

      len = tstop - tstart + 1;
      if (len == SeqLocLen (slp)) {
        if (SeqLocStrand (slp) == Seq_strand_minus) {
          for (i = tstop; i >= tstart; i--) {
            *ptr = bp [i];
            ptr++;
          }
        } else {
          for (i = tstart; i <= tstop; i++) {
            *ptr = bp [i];
            ptr++;
          }
        }
      }

    } else if (dsp->choice == 2) {

      slitp = (SeqLitPtr) dsp->data.ptrvalue;
      if (slitp == NULL || slitp->seq_data != NULL) continue;
      for (i = 0; i < slitp->length; i++) {
        *ptr = 255;
        ptr++;
      }
    }
  }

  /* now make graphs */

  i = 0;
  ptr = str;
  while (i < bsp->length) {
    by = *ptr;
    while (by == 255 && i < bsp->length) {
      i++;
      ptr++;
      by = *ptr;
    }
    if (i < bsp->length) {
      tstart = i;
      tmp = ptr;
      len = 0;
      max = INT2_MIN;
      min = INT2_MAX;
      while (by != 255 && i < bsp->length) {
        max = MAX (max, (Int2) by);
        min = MIN (min, (Int2) by);
        len++;
        i++;
        ptr++;
        by = *ptr;
      }
      tstop = i;
      sgp = SeqGraphNew ();
      if (sgp != NULL) {
        bs = BSNew (len + 1);
        if (bs != NULL) {
          BSWrite (bs, (Pointer) tmp, (Int4) len);
          sgp->numval = BSLen (bs);
          BSPutByte (bs, EOF);
          sgp->title = StringSave ("Phrap Quality");
          if (bsp->length != sgp->numval) {
            sgp->flags [0] = 1;
            sgp->compr = (len) / sgp->numval;
          } else {
            sgp->flags [0] = 0;
            sgp->compr = 1;
          }
          sgp->flags [1] = 0;
          sgp->flags [2] = 3;
          sgp->axis.intvalue = 0;
          sgp->min.intvalue = min;
          sgp->max.intvalue = max;
          sgp->a = 1.0;
          sgp->b = 0;
          sgp->values = (Pointer) bs;

          sintp = SeqIntNew ();
          sintp->from = tstart;
          sintp->to = tstop - 1;
          sintp->id = SeqIdDup (bsp->id);
          ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);

          if (lastsgp != NULL) {
            lastsgp->next = sgp;
          }
          lastsgp = sgp;

          if (sap == NULL) {
            sap = SeqAnnotNew ();
            if (sap != NULL) {
              SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave ("Graphs"));
              sap->type = 3;
              sap->data = (Pointer) sgp;
            }
          }
        }
      }
    }
  }

  /*
  sgp = SeqGraphNew ();
  if (sgp != NULL) {
    bs = BSNew (bsp->length);
    if (bs != NULL) {
      BSWrite (bs, (Pointer) str, (Int4) bsp->length);
      sgp->numval = BSLen (bs);
      BSPutByte (bs, EOF);
      sgp->title = StringSave ("Phrap Quality");
      if (bsp->length != sgp->numval) {
        sgp->flags [0] = 1;
        sgp->compr = (bsp->length) / sgp->numval;
      } else {
        sgp->flags [0] = 0;
        sgp->compr = 1;
      }
      sgp->flags [1] = 0;
      sgp->flags [2] = 3;
      sgp->axis.intvalue = 0;
      sgp->min.intvalue = min;
      sgp->max.intvalue = max;
      sgp->a = 1.0;
      sgp->b = 0;
      sgp->values = (Pointer) bs;

      sintp = SeqIntNew ();
      sintp->from = 0;
      sintp->to = bsp->length - 1;
      sintp->id = SeqIdDup (bsp->id);
      ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);

      if (lastsgp != NULL) {
        lastsgp->next = sgp;
      }
      lastsgp = sgp;

      if (sap == NULL) {
        sap = SeqAnnotNew ();
        if (sap != NULL) {
          SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave ("Graphs"));
          sap->type = 3;
          sap->data = (Pointer) sgp;
        }
      }
    }
  }
  */

  /* remove remaining single lock from components */

  for (vnp = phplist; vnp != NULL; vnp = vnp->next) {
    pdp = (PhrapDataPtr) vnp->data.ptrvalue;
    BioseqUnlock (pdp->bsp);
    MemFree (pdp->scores);
  }
  ValNodeFreeData (phplist);

  /* free master byte array */

  MemFree (str);

  return sap;
}

/* Functions for correcting capitalization */

typedef struct replaceitempair {
  CharPtr FindString;
  CharPtr ReplaceString;
} ReplaceItemPair, PNTR ReplaceItemPairPtr;


ReplaceItemPair AbbreviationList[] = {
 { "arabidopsis thaliana", "Arabidopsis thaliana" },
 { "adp", "ADP" },
 { "adp-", "ADP-" },
 { "atp", "ATP" },
 { "atp-", "ATP-" },
 { "bac", "BAC" },
 { "caenorhabditis elegans", "Caenorhabditis elegans" },
 { "cdna", "cDNA" },
 { "cdnas", "cDNAs" },
 { "coa", "CoA" },
 { "coi", "COI" },
 { "coii", "COII" },
 { "danio rerio", "Danio rerio" },
 { "dna", "DNA" },
 { "dna-", "DNA-" },
 { "drosophila melanogaster", "Drosophila melanogaster" },
 { "dsrna", "dsRNA" },
 { "escherichia coli", "Escherichia coli" },
 { "hiv", "HIV" },
 { "hiv-1", "HIV-1" },
 { "hiv-2", "HIV-2" },
 { "hnrna", "hnRNA" },
 { "homo sapiens", "Homo sapiens" },
 { "mhc", "MHC" },
 { "mrna", "mRNA" },
 { "mtdna", "mtDNA" },
 { "mus musculus", "Mus musculus" },
 { "nadh", "NADH" },
 { "nov.", "nov." },
 { "nov..", "nov.." },
 { "pcr", "PCR" },
 { "pcr-", "PCR-" },
 { "rattus norvegicus", "Rattus norvegicus" },
 { "rapd", "RAPD" },
 { "rdna", "rDNA" },
 { "rna", "RNA" },
 { "rna-", "RNA-" },
 { "rrna", "rRNA" },
 { "rt-pcr", "RT-PCR" },
 { "saccharomyces cerevisiae", "Saccharomyces cerevisiae" },
 { "scrna", "scRNA" },
 { "siv-1", "SIV-1" },
 { "snp", "SNP"     },
 { "snps", "SNPs"   },
 { "snrna", "snRNA" },
 { "sp.", "sp." },
 { "sp..", "sp.." },
 { "ssp.", "ssp." },
 { "ssp..", "ssp.." },
 { "ssrna", "ssRNA" },
 { "subsp.", "subsp." },
 { "subsp..", "subsp.." },
 { "trna", "tRNA" },
 { "trna-", "tRNA-" },
 { "var.", "var." },
 { "var..", "var.." },
 { "uk", "UK" },
 { "usa", "USA" },
 { "U.S.A.", "USA" },
 { "U.S.A", "USA" },
 { "United States of America", "USA" },
 {"(hiv)", "(HIV)" },
 {"(hiv1)", "(HIV1)" },
 {"(hiv-1)", "(HIV-1)" }
};

ReplaceItemPair SpecialAbbreviationList[] = {
 { "sp.", "sp." },
 { "nov.", "nov." },
 { "ssp.", "ssp." },
 { "var.", "var." },
 { "subsp.", "subsp." }
};

NLM_EXTERN void FixAbbreviationsInElement (CharPtr PNTR pEl)
{
  int i;
  CharPtr NewPtr;
  Boolean whole_word;

  if (pEl == NULL) return;
  if (*pEl == NULL) return;

  for (i = 0; i < sizeof (AbbreviationList) / sizeof (ReplaceItemPair); i++)
  {
    if (AbbreviationList[i].FindString[StringLen (AbbreviationList[i].FindString) - 1] == '-')
    {
      whole_word = FALSE;
    }
    else
    {
      whole_word = TRUE;
    }
    FindReplaceString (pEl, AbbreviationList[i].FindString,
            AbbreviationList[i].ReplaceString, FALSE, whole_word);
  }
  for (i = 0; i < sizeof (SpecialAbbreviationList) / sizeof (ReplaceItemPair); i++)
  {
    FindReplaceString (pEl, SpecialAbbreviationList[i].FindString,
            SpecialAbbreviationList[i].ReplaceString, FALSE, TRUE);
    if (StringLen (*pEl) >= StringLen (SpecialAbbreviationList[i].ReplaceString)
      && StringCmp ((*pEl) + StringLen (*pEl) - StringLen (SpecialAbbreviationList[i].ReplaceString), SpecialAbbreviationList[i].ReplaceString) == 0)
    {
      NewPtr = MemNew (StringLen (*pEl) + 2);
      if (NewPtr == NULL) return;
      StringCpy (NewPtr, *pEl);
      StringCat (NewPtr, ".");
      MemFree (*pEl);
      *pEl = NewPtr;
    }
  }
}

ReplaceItemPair ShortWordList[] = {
 { "A", "a" },
 { "About", "about" },
 { "And", "and" },
 { "At", "at" },
 { "But", "but" },
 { "By", "by" },
 { "For", "for" },
 { "In", "in" },
 { "Is", "is" },
 { "Of", "of" },
 { "On", "on" },
 { "Or", "or" },
 { "The", "the" },
 { "To", "to" },
 { "With", "with" }
};

static void FixShortWordsInElement (CharPtr PNTR pEl)
{
  Int2 i;

  if (pEl == NULL) return;
  if (*pEl == NULL) return;

  for (i = 0; i < sizeof (ShortWordList) / sizeof (ReplaceItemPair); i++)
  {
    FindReplaceString (pEl, ShortWordList[i].FindString,
            ShortWordList[i].ReplaceString, FALSE, TRUE);
  }
  if (isalpha ((Int4)((*pEl)[0])))
  {
    (*pEl)[0] = toupper ((*pEl)[0]);
  }
}

NLM_EXTERN void 
FixCapitalizationInElement 
(CharPtr PNTR pEl,
 Boolean      bAbbrev, 
 Boolean      bShortWords,
 Boolean      bApostrophes)
{
  CharPtr pCh;
  Boolean bSendToLower;
 
  if(pEl == NULL) return; 
  if(*pEl == NULL) return; 

  bSendToLower = FALSE;
  for(pCh = *pEl; *pCh != 0; pCh++)
  {
    if(isalpha((Int4)(*pCh)))
    {
      if(bSendToLower)
      {
        *pCh = tolower(*pCh);
      }
      else
      {
        *pCh = toupper(*pCh);
        bSendToLower = TRUE;
      }
    }
    else if (bApostrophes || *pCh != '\'')
    {
      bSendToLower = FALSE;
    }
  }
  if (bShortWords)
    FixShortWordsInElement (pEl);
  if (bAbbrev)
    FixAbbreviationsInElement (pEl);
}


static ReplaceItemPair s_CountryFixes[] = {
  { "chnia", "China" },
  { "pr china", "P.R. China" },
  { "prchina", "P.R. China" },
  { "p.r.china", "P.R. China" },
  { "p.r china", "P.R. China" },
  { "p, r, china", "P.R. China" },
  { "rok", "ROK" },
  { "rsa", "RSA" },
  { "roc", "ROC" },
  { "uae", "UAE" }
};

#define NUM_CountryFixes sizeof (s_CountryFixes) / sizeof (ReplaceItemPair)


static void InsertMissingSpacesAfterCommas (CharPtr PNTR pString)
{
  Int4 num_new_spaces = 0;
  CharPtr str, cp, new_str, src, dst;

  if (pString == NULL || *pString == NULL) {
    return;
  }

  str = *pString;
  cp = StringChr (str, ',');
  while (cp != NULL) {
    if (*(cp + 1) != 0 && !isspace (*(cp + 1))) {
      num_new_spaces++;
    }
    cp = StringChr (cp + 1, ',');
  }

  if (num_new_spaces == 0) {
    return;
  }

  new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + num_new_spaces + 1));
  src = str;
  dst = new_str;
  while (*src != 0) {
    *dst = *src;
    ++dst;
    if (*src == ',' && *(src + 1) != 0 && !isspace (*(src + 1))) {
      *dst = ' ';
      ++dst;
    }
    ++src;
  }
  *dst = 0;
  str = MemFree (str);
  *pString = new_str;
}


static void InsertMissingSpacesAfterNo (CharPtr PNTR pString)
{
  Int4 num_new_spaces = 0;
  CharPtr str, cp, new_str, src;

  if (pString == NULL || *pString == NULL) {
    return;
  }

  str = *pString;
  cp = StringISearch (str, "No.");
  while (cp != NULL) {
    if (isalpha(*(cp + 3)) || isdigit(*(cp + 3))) {
      num_new_spaces++;
    }
    cp = StringISearch (cp + 3, "No.");
  }

  new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + num_new_spaces + 1));
  new_str[0] = 0;

  src = str;
  cp = StringISearch (src, "No.");
  while (cp != NULL) {
    StringNCat (new_str, src, cp - src);
    StringCat (new_str, "No.");
    if (isalpha(*(cp + 3)) || isdigit(*(cp + 3))) {
      StringCat (new_str, " ");
    }
    src = cp + 3;
    cp = StringISearch (src, "No.");
  }
  StringCat (new_str, src);

  str = MemFree (str);
  *pString = new_str;
}


NLM_EXTERN void FixCapitalizationInCountryStringEx (CharPtr PNTR pCountry, Boolean punct_only)
{
  Int4 i;

  if (pCountry == NULL || StringICmp (*pCountry, "US") == 0) {
    return;
  }
  InsertMissingSpacesAfterCommas (pCountry);
  InsertMissingSpacesAfterNo (pCountry);
  if (!punct_only) {
    FixCapitalizationInElement (pCountry, TRUE, TRUE, FALSE);

  }
  for (i = 0; i < NUM_CountryFixes; i++) {
    FindReplaceString (pCountry, s_CountryFixes[i].FindString,
            s_CountryFixes[i].ReplaceString, FALSE, TRUE);
  }
}

NLM_EXTERN void FixCapitalizationInCountryString (CharPtr PNTR pCountry)
{
  FixCapitalizationInCountryStringEx (pCountry, FALSE);
}


NLM_EXTERN void FixCapitalizationInAuthor (AuthorPtr pAuthor)
{
  NameStdPtr pNameStandard;
  CharPtr    cp;
  
  if (pAuthor == NULL)
    return;
  else if(pAuthor->name->choice != 2)
    return;
  pNameStandard = pAuthor->name->data;
  if (pNameStandard != NULL)
  {
    FixCapitalizationInElement (&(pNameStandard->names[0]), FALSE, FALSE, TRUE);
    FixCapitalizationInElement (&(pNameStandard->names[1]), FALSE, FALSE, FALSE);
    /* Set initials to all caps */
    for (cp = pNameStandard->names[4]; cp != NULL && *cp != 0; cp++)
    {
      *cp = toupper (*cp);
    }
  }
}


NLM_EXTERN void FixStateAbbreviationsInAffil (AffilPtr affil, LogInfoPtr lip)
{
  CharPtr abbrev;

  if (affil == NULL) {
    return;
  }
  if (StringCmp (affil->country, "USA") == 0) {
    abbrev = GetStateAbbreviation (affil->sub);
    if (abbrev != NULL) {
      if (lip != NULL) {
        if (lip->fp != NULL) {
          fprintf (lip->fp, "Changed %s to %s\n", affil->sub, abbrev);
        }
        lip->data_in_log = TRUE;
      }
      affil->sub = MemFree (affil->sub);
      affil->sub = StringSave (abbrev);
    }
  }
}


NLM_EXTERN void FixCapsInPubAffilEx (AffilPtr affil, Boolean punct_only)
{
  if (affil == NULL) return;
  if (!punct_only) {
    FixCapitalizationInElement (&(affil->affil), TRUE, TRUE, FALSE);
    FixAffiliationShortWordsInElement (&(affil->affil));
    FixCapitalizationInElement (&(affil->div), TRUE, TRUE, FALSE);
    FixAffiliationShortWordsInElement (&(affil->div));
    FixCapitalizationInElement (&(affil->city), FALSE, TRUE, FALSE);
    FixAffiliationShortWordsInElement (&(affil->city));
  }
  FixKnownAbbreviationsInElement (&(affil->affil));
  FixKnownAbbreviationsInElement (&(affil->street));
  FixKnownAbbreviationsInElement (&(affil->div));
  FixKnownAbbreviationsInElement (&(affil->city));
  
  InsertMissingSpacesAfterCommas (&(affil->affil));
  InsertMissingSpacesAfterNo (&(affil->affil));
  InsertMissingSpacesAfterCommas (&(affil->div));
  InsertMissingSpacesAfterNo (&(affil->div));
  InsertMissingSpacesAfterCommas (&(affil->city));
  InsertMissingSpacesAfterNo (&(affil->city));

  /* special handling for states */
  if (punct_only) {
    InsertMissingSpacesAfterCommas (&(affil->sub));
  } else {
    if (affil->sub != NULL && StringLen (affil->sub) == 2
      && isalpha((Int4)(affil->sub[0]))	&& isalpha((Int4)(affil->sub[1])))
    {
      affil->sub[0] = toupper(affil->sub[0]);
      affil->sub[1] = toupper(affil->sub[1]);
    } else {
      FixCapitalizationInElement (&(affil->sub), FALSE, TRUE, FALSE);
      FixAffiliationShortWordsInElement (&(affil->sub));
      InsertMissingSpacesAfterCommas (&(affil->sub));
    }
  }

  if (!punct_only) {
    FixCapitalizationInCountryString (&(affil->country));
    FixCapitalizationInElement (&(affil->street), FALSE, TRUE, FALSE);
    FixAffiliationShortWordsInElement (&(affil->street));
    FixStateAbbreviationsInAffil (affil, NULL);
  }
  if (StringCmp (affil->country, "USA") == 0) {
    FixStateAbbreviationsInAffil (affil, NULL);
  }
  InsertMissingSpacesAfterCommas (&(affil->street));
  InsertMissingSpacesAfterNo (&(affil->street));
}


NLM_EXTERN void FixCapsInPubAffil (AffilPtr affil)
{
  FixCapsInPubAffilEx (affil, FALSE);
}

ReplaceItemPair AffiliationShortWordList[] = {
 { "Au", "au" } ,
 { "Aux", "aux" } ,
 { "A La", "a la" } ,
 { "De La", "de la" } ,
 { "De", "de" } ,
 { "Del", "del" } ,
 { "Des", "des" } ,
 { "Du", "du" } ,
 { "Et", "et" } ,
 { "La", "la" },
 { "Le", "le" },
 { "Les", "les" },
 { "Rue", "rue" },
 { "Po Box", "PO Box" },
 { "Pobox", "PO Box" },
 { "P.O box", "P.O. Box" },
 { "P.Obox", "P.O. Box" },
 { "Y", "y" }
};

NLM_EXTERN void FixAffiliationShortWordsInElement (CharPtr PNTR pEl)
{
  Int2 i;
  CharPtr cp;

  if (pEl == NULL) return;
  if (*pEl == NULL) return;

  for (i = 0; i < sizeof (AffiliationShortWordList) / sizeof (ReplaceItemPair); i++)
  {
    FindReplaceString (pEl, AffiliationShortWordList[i].FindString,
            AffiliationShortWordList[i].ReplaceString, FALSE, TRUE);
  }
  if (isalpha ((Int4)((*pEl)[0])))
  {
    (*pEl)[0] = toupper ((*pEl)[0]);
  }

  /* fix d' */
  cp = StringStr (*pEl, "D'");
  while (cp != NULL) {
    if (cp == *pEl || !isalpha(*(cp - 1))) {
      *cp = 'd';
      if (isalpha (*(cp + 2))) {
        *(cp + 2) = toupper(*(cp + 2));
      }
    }    
    cp = StringStr (cp + 1, "D'");
  }
}


ReplaceItemPair KnownAbbreviationList[] = {
 { "po box", "PO Box" },
 { "Pobox", "PO Box" },
 { "P.O box", "P.O. Box" },
 { "P.Obox", "P.O. Box" },
 { "PO.Box", "P.O. Box" },
 { "PO. Box", "P.O. Box" },
 { "pr china", "P.R. China" },
 { "prchina", "P.R. China" },
 { "p.r.china", "P.R. China" },
 { "p.r china", "P.R. China" },
 { "p, r, china", "P.R. China" },
 { "p,r, china", "P.R. China" },
 { "p,r,china", "P.R. China" }
};

NLM_EXTERN void FixKnownAbbreviationsInElement (CharPtr PNTR pEl)
{
  Int2 i;

  if (pEl == NULL) return;
  if (*pEl == NULL) return;

  for (i = 0; i < sizeof (KnownAbbreviationList) / sizeof (ReplaceItemPair); i++)
  {
    FindReplaceString (pEl, KnownAbbreviationList[i].FindString,
            KnownAbbreviationList[i].ReplaceString, FALSE, TRUE);
  }
}


NLM_EXTERN void FixOrgNamesInString (CharPtr str, ValNodePtr org_names)
{
  ValNodePtr vnp;
  CharPtr    cp, taxname;
  Int4       taxname_len;
  
  if (StringHasNoText (str) || org_names == NULL) return;
  for (vnp = org_names; vnp != NULL; vnp = vnp->next)
  {
    taxname = (CharPtr) org_names->data.ptrvalue;
    taxname_len = StringLen (taxname);
    cp = StringISearch (str, taxname);
    while (cp != NULL)
    {
      StringNCpy (cp, taxname, taxname_len);
      cp = StringISearch (cp + taxname_len, taxname);
    }
  }
}


NLM_EXTERN void ResetCapitalization (Boolean first_is_upper, CharPtr pString)
{
  CharPtr pCh;
  Boolean was_digit = FALSE;

  pCh = pString;
  if (pCh == NULL) return;
  if (*pCh == '\0') return;
  
  if (first_is_upper)
  {
    /* Set first character to upper */
    *pCh = toupper (*pCh);    
  }
  else
  {
    /* set first character to lower */
    *pCh = tolower (*pCh);    
  }
  
  if (isdigit ((Int4)(*pCh)))
  {
      was_digit = TRUE;
  }
  pCh++;
  /* Set rest of characters to lower */
  while (*pCh != '\0')
  {
    if (was_digit 
        && (*pCh == 'S' || *pCh == 's') 
        && (isspace ((Int4)(*(pCh + 1))) || *(pCh + 1) == 0))
    {
      *pCh = toupper (*pCh);
      was_digit = FALSE;
    }
    else if (isdigit ((Int4)(*pCh)))
    {
      was_digit = TRUE;
    }
    else
    {
      was_digit = FALSE;
      *pCh = tolower (*pCh);
    }
    pCh++;
  }  
}


NLM_EXTERN SeqIdPtr CreateSeqIdFromText (CharPtr id_str, SeqEntryPtr sep)
{
  BioseqPtr   bsp = NULL;
  CharPtr     tmpstr;
  SeqIdPtr    sip = NULL;
  SeqEntryPtr scope;
  
  if (StringStr (id_str, "|") == NULL) {
    tmpstr = (CharPtr) MemNew (sizeof (Char) * (StringLen (id_str) + 20));
    sprintf (tmpstr, "lcl|%s", id_str);
    sip = SeqIdParse (tmpstr);
    if (sip != NULL) {
      scope = SeqEntrySetScope (sep);
      bsp = BioseqFind (sip);
      SeqEntrySetScope (scope);
      if (bsp == NULL) {
        sip = SeqIdFree (sip);
      }
    } 
    if (bsp == NULL) {
      sprintf (tmpstr, "gb|%s", id_str);
      sip = SeqIdParse (tmpstr);
      if (sip != NULL) {
        scope = SeqEntrySetScope (sep);
        bsp = BioseqFind (sip);
        SeqEntrySetScope (scope);
        if (bsp == NULL) {
          sip = SeqIdFree (sip);
        }
      }
    }
    if (bsp == NULL) {
      sprintf (tmpstr, "gnl|%s", id_str);
      sip = SeqIdParse (tmpstr);
      if (sip != NULL) {
        scope = SeqEntrySetScope (sep);
        bsp = BioseqFind (sip);
        SeqEntrySetScope (scope);
        if (bsp == NULL) {
          sip = SeqIdFree (sip);
        }
      }
    }
    if (bsp == NULL) {
      if (StringNICmp (id_str, "bankit", 6) == 0) {
        sprintf (tmpstr, "gnl|BankIt|%s", id_str + 6);
      } else {
        sprintf (tmpstr, "gnl|BankIt|%s", id_str);
      }
      sip = SeqIdParse (tmpstr);
      if (sip != NULL) {
        scope = SeqEntrySetScope (sep);
        bsp = BioseqFind (sip);
        SeqEntrySetScope (scope);
        if (bsp == NULL) {
          sip = SeqIdFree (sip);
        }
      }
    }

    if (bsp == NULL) {
      sprintf (tmpstr, "gnl|NCBIFILE|%s", id_str);
      sip = SeqIdParse (tmpstr);
      if (sip != NULL) {
        scope = SeqEntrySetScope (sep);
        bsp = BioseqFind (sip);
        SeqEntrySetScope (scope);
        if (bsp == NULL) {
          sip = SeqIdFree (sip);
        }
      }
    }

    if (bsp == NULL) {
      sprintf (tmpstr, "ref|%s", id_str);
      sip = SeqIdParse (tmpstr);
      if (sip != NULL) {
        scope = SeqEntrySetScope (sep);
        bsp = BioseqFind (sip);
        SeqEntrySetScope (scope);
        if (bsp == NULL) {
          sip = SeqIdFree (sip);
        }
      }
    }

    if (bsp == NULL && IsAllDigits (id_str)) {
      sprintf (tmpstr, "gi|%s", id_str);
      sip = SeqIdParse (tmpstr);
      if (sip != NULL) {
        scope = SeqEntrySetScope (sep);
        bsp = BioseqFind (sip);
        SeqEntrySetScope (scope);
        if (bsp == NULL) {
          sip = SeqIdFree (sip);
        }
      }
    }
    MemFree (tmpstr);
  } else {
    sip = SeqIdParse (id_str);
    if (sip != NULL) {
      scope = SeqEntrySetScope (sep);
      bsp = BioseqFind (sip);
      SeqEntrySetScope (scope);
      if (bsp == NULL) {
        sip = SeqIdFree (sip);
      }
    }
  }        
  return sip;
}


NLM_EXTERN SeqLocPtr SeqLocWholeNew (BioseqPtr bsp)
{
  ValNodePtr vnp;

  if (bsp == NULL) return NULL;

  vnp = ValNodeNew (NULL);

  if (vnp == NULL) return NULL;

  vnp->choice = SEQLOC_WHOLE;
  vnp->data.ptrvalue = (Pointer) SeqIdDup (SeqIdFindBest (bsp->id, 0));
  return (SeqLocPtr)vnp;
}


NLM_EXTERN Int4 GetDeltaSeqLen (DeltaSeqPtr dsp)
{
  Int4 len = 0;
  SeqLitPtr slp;

  if (dsp == NULL || dsp->data.ptrvalue == NULL) {
    /* do nothing, empty */
  } else if (dsp->choice == 1) {
    len = SeqLocLen ((SeqLocPtr)(dsp->data.ptrvalue));
  } else if (dsp->choice == 2) {
    slp = (SeqLitPtr) dsp->data.ptrvalue;
    len = slp->length;
  }
  return len;
}


/* The following section of code is used for retranslating a CDS and updating
 * the protein features based on an alignment between the old and new protein
 * sequences.
 */
static Int4 SeqEdRemapCoord (SeqAlignPtr salp, Int4 coord, Boolean move_up, Int4 len)

{
  Int4 aln_pos;
  
  if (salp == NULL) return -1;
  aln_pos = AlnMgr2MapBioseqToSeqAlign (salp, coord, 1);
  while (aln_pos == -1)
  {
      if (move_up)
      {
        if (coord >= len - 1)
        {
            return len - 1;
        }
        else
        {
            coord ++;
        }
      }
      else
      {
        if (coord <= 0)
        {
            return 0;
        }
        else
        {
            coord --;
        }
      }
      aln_pos = AlnMgr2MapBioseqToSeqAlign (salp, coord, 1);
  }
  return AlnMgr2MapSeqAlignToBioseq (salp, aln_pos, 2);
}


static void SeqEdRemapSeqIntLoc (SeqAlignPtr salp, SeqIntPtr sintp, Int4 seq_len)

{  
  if (salp == NULL || sintp == NULL) return;
  sintp->from = SeqEdRemapCoord (salp, sintp->from, TRUE, seq_len);
  sintp->to = SeqEdRemapCoord (salp, sintp->to, FALSE, seq_len); 
}

static void SeqEdRemapSeqPntLoc (SeqAlignPtr salp, SeqPntPtr spp, Int4 seq_len)

{  
  if (salp == NULL || spp == NULL) return;
  
  spp->point = SeqEdRemapCoord (salp, spp->point, FALSE, seq_len);
}


static void SeqEdRemapPackSeqPnt (SeqAlignPtr salp, PackSeqPntPtr pspp, Int4 seq_len)

{
  Uint1          used;

  if (salp == NULL || pspp == NULL) return;
  for (used = 0; used < pspp->used; used++) 
  {
    pspp->pnts [used] = SeqEdRemapCoord (salp, pspp->pnts [used], FALSE, seq_len);
  }
}


NLM_EXTERN void SeqEdRemapLocation (SeqAlignPtr salp, SeqLocPtr slp, Int4 seq_len)

{

  if (slp == NULL) return;
  switch (slp->choice) {
    case SEQLOC_INT :
      SeqEdRemapSeqIntLoc (salp, slp->data.ptrvalue, seq_len);
      break;
    case SEQLOC_PNT :
      SeqEdRemapSeqPntLoc (salp, slp->data.ptrvalue, seq_len);
      break;
    case SEQLOC_PACKED_PNT :
      SeqEdRemapPackSeqPnt (salp, slp->data.ptrvalue, seq_len);
      break;
    default :
      break;
  }
}


static void MakeLocationMatchEntireSequence (SeqLocPtr slp, BioseqPtr bsp)
{
  SeqIntPtr sip;
  
  if (slp == NULL || bsp == NULL) return;
  
  if (slp->choice == SEQLOC_WHOLE)
  {
      SeqIdFree (slp->data.ptrvalue);
      slp->data.ptrvalue = SeqIdDup (bsp->id);
  }
  else if (slp->choice == SEQLOC_INT)
  {
    sip = (SeqIntPtr) slp->data.ptrvalue;
    if (sip == NULL)
    {
      sip = SeqIntNew ();
      slp->data.ptrvalue = sip;
    }
    if (sip != NULL)
    {
      sip->from = 0;
      sip->to = bsp->length - 1;
    }
  }
}


NLM_EXTERN Boolean SeqEdFixProteinFeatures (BioseqPtr oldbsp, BioseqPtr newbsp, Boolean force_fix, GlobalAlignFunc align_func)
{
  SeqAlignPtr       salp = NULL;
  Boolean           revcomp = FALSE;
  SeqFeatPtr        sfp;
  Boolean           tried_to_get_alignment = FALSE;
  Boolean           unmappable_feats = FALSE;
  SeqLocPtr         slp_tmp = NULL;
  SeqAnnotPtr       sap;
  ProtRefPtr        prp;
  
  if (oldbsp == NULL || newbsp == NULL) return FALSE;

  /* get alignment between old and new proteins */

  if (ISA_na (oldbsp->mol) != ISA_na (newbsp->mol)) return FALSE;

  /* iterate through the features on the old protein.  Full length features
   * should be set to the new length.  Other features should be mapped through
   * the alignment (if possible), otherwise warn the user that they could not
   * be remapped. */      

  if (!force_fix)
  {
    for (sap = oldbsp->annot; sap != NULL && !unmappable_feats; sap = sap->next) {
      if (sap->type == 1) {
        for (sfp = sap->data; sfp != NULL && !unmappable_feats; sfp = sfp->next) {
          if (sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL) {
            continue;
          }
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp->processed != 0)
          {
            if (salp == NULL)
            {
              if (align_func != NULL)
              {
                salp = align_func (oldbsp, newbsp, &revcomp);
              }
            }
            if (salp == NULL)
            {
              unmappable_feats = TRUE;
            } 
            else
            {
              slp_tmp = (SeqLocPtr) AsnIoMemCopy (sfp->location,
                                                  (AsnReadFunc) SeqLocAsnRead,
                                                  (AsnWriteFunc) SeqLocAsnWrite);
              SeqEdRemapLocation (salp, slp_tmp, newbsp->length);                
              if (slp_tmp == NULL)
              {
                unmappable_feats = TRUE;
              }
              else
              {
                slp_tmp = SeqLocFree (slp_tmp);
              }
            }
          }
        }
      }
    }
    if (unmappable_feats)
    {
      return FALSE;
    }
  }

  for (sap = oldbsp->annot; sap != NULL && !unmappable_feats; sap = sap->next) {
    if (sap->type == 1) {
      for (sfp = sap->data; sfp != NULL && !unmappable_feats; sfp = sfp->next) {
        if (sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL) {
          continue;
        }
        prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        if (prp->processed == 0)
        {
          /* make new location match new sequence length */
          MakeLocationMatchEntireSequence (sfp->location, newbsp);
        }
        else
        {
          if (salp == NULL && !tried_to_get_alignment && align_func != NULL)
          {
            salp = align_func (oldbsp, newbsp, &revcomp);
            tried_to_get_alignment = TRUE;
          }
          if (salp != NULL)
          {
            SeqEdRemapLocation (salp, sfp->location, newbsp->length);                
          }
          else
          {
            unmappable_feats = TRUE;
          }
        }
      }
    }
  } 
        
  if (salp != NULL)
  {
    SeqAlignFree (salp);
  }
  if (unmappable_feats)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


NLM_EXTERN void SeqEdTranslateOneCDS (SeqFeatPtr sfp, BioseqPtr featbsp, Uint2 entityID, GlobalAlignFunc align_func)
{
  ByteStorePtr  bs;
  Char          ch;
  CharPtr       prot;
  CharPtr       ptr;
  Int4          star_at_end = 0;
  BioseqPtr     old_prot;
  SeqIdPtr      new_prot_id;
  SeqEntryPtr   parent, new_prot_sep;
  SeqLocPtr     slp;
  Uint1         seq_data_type;
  Int4          old_length;
  BioseqPtr     newbsp;
  ProtRefPtr    prp;
  SeqFeatPtr    prot_sfp;
  SeqDataPtr    sdp;
  
  if (featbsp == NULL || sfp == NULL || sfp->location == NULL 
      || sfp->data.choice != SEQFEAT_CDREGION) 
  {
      return;
  }
  
  old_prot = BioseqFindFromSeqLoc (sfp->product);
  new_prot_id = SeqIdDup (SeqLocId (sfp->product));
  if (new_prot_id == NULL)
  {
      new_prot_id = MakeNewProteinSeqId (sfp->location, featbsp->id);
  }
      
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs != NULL) {
    prot = BSMerge (bs, NULL);
    bs = BSFree (bs);
    if (prot != NULL) {
      ptr = prot;
      ch = *ptr;
      while (ch != '\0') {
        *ptr = TO_UPPER (ch);
        if (ch == '*') {
          star_at_end = 1;
        } else {
          star_at_end = 0;
        }
        ptr++;
        ch = *ptr;
      }
      if (star_at_end)
      {
          *(ptr - 1) = 0;
      }
      bs = BSNew (1000);
      if (bs != NULL) {
        ptr = prot;
        BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
      }
      MemFree (prot);
    }
    newbsp = BioseqNew ();
    if (newbsp != NULL) {
      newbsp->id = SeqIdParse ("lcl|CdRgnTransl");
      newbsp->repr = Seq_repr_raw;
      newbsp->mol = Seq_mol_aa;
      newbsp->seq_data_type = Seq_code_ncbieaa;
      newbsp->seq_data = (SeqDataPtr) bs;
      newbsp->length = BSLen (bs);
      
      if (old_prot == NULL)
      {
          /* need to create a new protein sequence */
          SeqIdFree (newbsp->id);
          newbsp->id = new_prot_id;
          new_prot_sep = SeqEntryNew ();
          new_prot_sep->choice = 1;
          new_prot_sep->data.ptrvalue = newbsp;
          parent = GetBestTopParentForData (entityID, featbsp);
          if (parent != NULL)
          {
            AddSeqEntryToSeqEntry (parent, new_prot_sep, TRUE);
          }
          slp = ValNodeNew (NULL);
          if (slp != NULL)
          {
            slp->choice = SEQLOC_WHOLE;
            slp->data.ptrvalue = SeqIdDup (new_prot_id);
          }
          sfp->product = slp;
          
          /* create full length protein feature */
        prp = ProtRefNew ();
        prot_sfp = CreateNewFeature (new_prot_sep, NULL, SEQFEAT_PROT, NULL);
        if (prot_sfp != NULL) {
          prot_sfp->data.value.ptrvalue = (Pointer) prp;
        }
      }
      else
      {
        /* propagate features to new protein */
        if (!SeqEdFixProteinFeatures (old_prot, newbsp, TRUE, align_func))
        {
          Message (MSG_ERROR, "Unable to construct alignment between old and new "
                  "proteins - you will need to adjust the protein features "
                  "manually.");
        }
                 
          /* then replace old protein with new */
          seq_data_type = old_prot->seq_data_type;
          sdp = old_prot->seq_data;
          old_length = old_prot->length;
          old_prot->seq_data_type = newbsp->seq_data_type;
          old_prot->seq_data = newbsp->seq_data;
          old_prot->length = newbsp->length;
          newbsp->seq_data_type = seq_data_type;
          newbsp->seq_data = sdp;
          newbsp->length = old_length;
          BioseqFree (newbsp);
      }
    }
  }
}


static SeqLocPtr RemoveGapsFromSegmentedLocation (SeqLocPtr slp, BioseqPtr bsp)
{
  SeqLocPtr loc_slp, slp_new, slp_tmp, loc_list = NULL, loc_last = NULL;
  SeqIdPtr  sip;
  Uint1     strand;
  Int4      seq_offset, start_pos, end_pos, piece_len;
  
  if (slp == NULL || bsp == NULL || bsp->repr != Seq_repr_seg)
  {
    return slp;
  }
  
  loc_slp = SeqLocFindNext (slp, NULL);

  while (loc_slp != NULL)
  {
    strand = SeqLocStrand (loc_slp);
    start_pos = SeqLocStart (loc_slp);
    end_pos = SeqLocStop (loc_slp);
  
    /* create list of locations */
    slp_tmp = bsp->seq_ext;
    seq_offset = 0;
    while (slp_tmp != NULL)
    {
      piece_len = SeqLocLen (slp_tmp);
      
      if (seq_offset < end_pos
          && seq_offset + piece_len >= start_pos)
      {
        sip = SeqLocId (slp_tmp);
      
        slp_new = SeqLocIntNew (MAX (0, start_pos - seq_offset),
                                MIN (piece_len, end_pos - seq_offset),
                                strand, sip);

        if (slp_new != NULL)
        {
          if (loc_last == NULL)
          {
            loc_list = slp_new;
          }
          else
          {
            loc_last->next = slp_new;
          }
          loc_last = slp_new;
        }
      }
      seq_offset += piece_len;
      slp_tmp = slp_tmp->next;
    }
    loc_slp = SeqLocFindNext (slp, loc_slp);
  }
  
  if (loc_list == NULL)
  {
    /* failed to convert - do not change */
  }
  else if (loc_list->next == NULL)
  {
    /* only found one piece */
    slp = SeqLocFree (slp);
    slp = loc_list;
  }
  else
  {
    /* make mixed location */
    slp_new = ValNodeNew (NULL);
    if (slp_new != NULL)
    {
      slp_new->choice = SEQLOC_MIX;
      slp_new->data.ptrvalue = loc_list;
      slp = SeqLocFree (slp);
      slp = slp_new;
    }
  }
  return slp;
}

static Boolean PointInInterval (Int4 interval_start, Int4 interval_length, Int4 point)
{
  if (point >= interval_start && point < interval_start + interval_length) 
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

NLM_EXTERN Boolean GapInLocation (Int4 seq_offset, Int4 length, SeqLocPtr loc)
{
  SeqLocPtr slp;
  Int4      start, stop;
  Boolean   gap_in_location = FALSE;
  
  slp = SeqLocFindNext (loc, NULL);
  
  while (slp != NULL && ! gap_in_location)
  {
    start = SeqLocStart (slp);
    stop = SeqLocStop (slp);
    
    if (PointInInterval (seq_offset, length, start)
        || PointInInterval (seq_offset, length, stop)
        || PointInInterval (start, stop - start + 1, seq_offset)
        || PointInInterval (start, stop - start + 1, seq_offset + length - 1))
    {
      gap_in_location = TRUE;
    }
    
    slp = SeqLocFindNext (loc, slp);
  }
  return gap_in_location;
}


static Boolean AdjustForThisGap (DeltaSeqPtr dsp, Uint4 options, Int4 seq_offset, SeqLocPtr before)
{
  Int4 dsp_len;
  Int4 range_start, range_stop, tmp;

  dsp_len = GetDeltaSeqLen(dsp);
  if (!IsDeltaSeqGap(dsp)) {
    return FALSE;
  }

  if (!(options & eAdjustFeatForGap_unknown_gaps) && IsDeltaSeqUnknownGap(dsp)) {
    return FALSE;
  }

  if (!(options & eAdjustFeatForGap_known_gaps) && IsDeltaSeqKnownGap (dsp)) {
    return FALSE;
  }

  if (options & eAdjustFeatForGap_split_in_intron) {
    range_start = SeqLocStart (before);
    range_stop = SeqLocStop (before);
    if (range_stop < range_start) {
      tmp = range_stop;
      range_stop = range_start;
      range_start = tmp;
    }
    if (seq_offset >= range_start && seq_offset <= range_stop) {
      return TRUE;
    } else if (seq_offset + dsp_len >= range_start && seq_offset + dsp_len < range_stop) {
      return TRUE;
    } else if (seq_offset <= range_start && seq_offset + dsp_len > range_stop) {
      return TRUE;
    }
  }

  if (GapInLocation (seq_offset, dsp_len, before)) {
    return TRUE;
  }

  return FALSE;
}


NLM_EXTERN void
LocationContainsGaps
(SeqLocPtr slp,
 BioseqPtr bsp,
 Uint4     options,
 BoolPtr   terminal_gaps,
 BoolPtr   internal_gaps,
 BoolPtr   entirely_in_gap)
{
  DeltaSeqPtr dsp;
  Int4        seq_offset = 0;
  Int4        dsp_len;
  Boolean     has_terminal_gaps = FALSE;
  Boolean     has_internal_gaps = FALSE;
  Boolean     all_sublocs_in_gap = TRUE, this_subloc_in_gap = FALSE, this_subloc_start_gap = FALSE, this_subloc_stop_gap = FALSE;
  Boolean     right_sublocs_in_gap = FALSE;
  SeqLocPtr   tmp_slp;
  Int4        start, stop;

  if (slp == NULL || bsp == NULL 
      || bsp->repr != Seq_repr_delta
      || bsp->seq_ext_type != 4
      || bsp->seq_ext == NULL)
  {
    return;
  }
    
  for (tmp_slp = SeqLocFindNext (slp, NULL); 
       tmp_slp != NULL;
       tmp_slp = SeqLocFindNext (slp, tmp_slp))
  {
    seq_offset = 0;
    start = SeqLocStart (tmp_slp);
    stop = SeqLocStop (tmp_slp);
    this_subloc_in_gap = FALSE;
    this_subloc_start_gap = FALSE;
    this_subloc_stop_gap = FALSE;
    for (dsp = (DeltaSeqPtr) bsp->seq_ext;
         dsp != NULL && seq_offset <= stop && !this_subloc_in_gap;
         dsp = dsp->next)
    {
      dsp_len = GetDeltaSeqLen(dsp);
      if (AdjustForThisGap (dsp, options, seq_offset, tmp_slp))
      {
        if (PointInInterval (seq_offset, dsp_len, start)
            && PointInInterval (seq_offset, dsp_len, stop))
        {
          this_subloc_in_gap = TRUE;
        }
        else if (PointInInterval (seq_offset, dsp_len, start))
        {
          this_subloc_start_gap = TRUE;
        }
        else if (PointInInterval (seq_offset, dsp_len, stop)) 
        {
          this_subloc_stop_gap = TRUE;
        }
        else 
        {
          has_internal_gaps = TRUE;
        }
      }
      seq_offset += dsp_len;
    }

    if (this_subloc_in_gap)
    {
      /* all sublocs up to this point have been in the gap, so still part of left gap */
      if (all_sublocs_in_gap)
      {
        has_terminal_gaps = TRUE;
      }
      /* could be part of chain of sublocs on the right in gap */
      right_sublocs_in_gap = TRUE;
    }
    else
    {
      if (right_sublocs_in_gap && !all_sublocs_in_gap)
      {
        /* a chain of prior gapped sublocs has ended that did not start at the left end */
        has_internal_gaps = TRUE;
      }

      if (this_subloc_start_gap)
      { 
        if (all_sublocs_in_gap)
        {
          has_terminal_gaps = TRUE;
        }
        else
        {
          /* gap on left, not part of chain of sublocs on left in gap */
          has_internal_gaps = TRUE;
        }
      } 
      
      if (this_subloc_stop_gap)
      {
        /* could be start of chain of sublocs on the right in gap */
        right_sublocs_in_gap = TRUE;
      }
      else
      {
        right_sublocs_in_gap = FALSE;
      }
      /* at least this subloc is not completely contained in a gap */
      all_sublocs_in_gap = FALSE;
    }
  }

  if (right_sublocs_in_gap)
  {
    has_terminal_gaps = TRUE;
  }

  if (all_sublocs_in_gap)
  {
    has_terminal_gaps = FALSE;
    has_internal_gaps = FALSE;
  }

  if (terminal_gaps != NULL)
  {
    *terminal_gaps = has_terminal_gaps;
  }

  if (internal_gaps != NULL)
  {
    *internal_gaps = has_internal_gaps;
  }

  if (entirely_in_gap != NULL)
  {
    *entirely_in_gap = all_sublocs_in_gap;
  }
}


static SeqLocPtr 
RemoveGapsFromDeltaLocation 
(SeqLocPtr slp,
 BioseqPtr bsp,
 Uint4     options,
 Boolean   set_partial_ends,
 BoolPtr   split)
{
  DeltaSeqPtr dsp;
  Int4        seq_offset = 0;
  Int4        dsp_len;
  SeqLocPtr   loc_list = NULL, prev_loc = NULL;
  SeqIdPtr    sip, before_sip, after_sip;
  Boolean     changed, partial5, partial3;
  SeqLocPtr   before = NULL, after = NULL;
  Int4        start, stop;
  Boolean     delete_for_this_gap;
  Uint1       strand;

  if (slp == NULL || bsp == NULL 
      || bsp->repr != Seq_repr_delta
      || bsp->seq_ext_type != 4
      || bsp->seq_ext == NULL)
  {
    return slp;
  }
  
  sip = SeqLocId (slp);
  if (sip == NULL)
  {
    return slp;
  }
  
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  
  seq_offset = 0;
  before = SeqLocCopy (slp);
  loc_list = before;
  
  for (dsp = (DeltaSeqPtr) bsp->seq_ext;
       dsp != NULL;
       dsp = dsp->next)
  {
    dsp_len = GetDeltaSeqLen(dsp);
    if (AdjustForThisGap (dsp, options, seq_offset, before))
    {
      delete_for_this_gap = TRUE;
      start = SeqLocStart (before);
      stop = SeqLocStop (before);
      if (PointInInterval (seq_offset, dsp_len, start)
          && PointInInterval (seq_offset, dsp_len, stop))
      {
        loc_list = SeqLocFree (loc_list);
        before = NULL;
        after = NULL;
        break;
      }
      else if (!PointInInterval (seq_offset, dsp_len, start) 
          && !PointInInterval (seq_offset, dsp_len, stop)) 
      {
        if (!(options & eAdjustFeatForGap_split_internal)) 
        {
          delete_for_this_gap = FALSE;
        }
        else
        {
          if (split != NULL)
          {
            *split = TRUE;
          }
        }
      }
      else if (!(options & eAdjustFeatForGap_trim_ends))
      {
        delete_for_this_gap = FALSE;
      }

      if (delete_for_this_gap)
      {
        /* get strand - determines direction of partials */
        strand = SeqLocStrand (before);

        /* we make a copy of the original location */
        after = SeqLocCopy (before);
        
        /* note - we need to use duplicates of the SeqID returned by
         * SeqLocId, just in case the first location in a mixed location 
         * is deleted, which would free the result from SeqLocId 
         */
        after_sip = SeqIdDup (SeqLocId (after));
        before_sip = SeqIdDup (SeqLocId (before));
        /* in the "after" location, we free everything before the 
         * end of the gap.
         */
        after = SeqLocDeleteEx (after, after_sip,
                          0, seq_offset + dsp_len - 1,
                          FALSE, &changed, &partial5, &partial3);
                          
        /* in the "before" location, we free everything after the
         * beginning of the gap.
         */
        before = SeqLocDeleteEx (before, before_sip, 
                          seq_offset, bsp->length,
                          FALSE, &changed, &partial5, &partial3);
        
        /* handle partialness for ends */
        CheckSeqLocForPartial (after, &partial5, &partial3);
        if (strand == Seq_strand_minus) 
        {
          if (before == NULL)
          {
            /* truncated at 3' end */
            SetSeqLocPartial (after, partial5, set_partial_ends);
          }
          else
          {
            SetSeqLocPartial (after, partial5, TRUE);
          }

        }
        else
        {        
          if (before == NULL)
          {
            /* truncated at 5' end*/
            SetSeqLocPartial (after, set_partial_ends, partial3);
          }
          else
          {
            SetSeqLocPartial (after, TRUE, partial3);
          }
        }
       
        CheckSeqLocForPartial (before, &partial5, &partial3);
        if (strand == Seq_strand_minus)
        {
          if (after == NULL) 
          {
            /* truncated at 5' end*/
            SetSeqLocPartial (before, set_partial_ends, partial3);
          } else {
            SetSeqLocPartial (before, TRUE, partial3);
          }
        }
        else
        {
          if (after == NULL) 
          {
            /* truncated */
            SetSeqLocPartial (before, partial5, set_partial_ends);
          } else {
            SetSeqLocPartial (before, partial5, TRUE);
          }
        }
        
        /* we're done with these IDs now */                  
        after_sip = SeqIdFree (after_sip);
        before_sip = SeqIdFree (before_sip);                          
                          
        if (before == NULL)
        {
          if (prev_loc == NULL)
          {
            loc_list = after;
          }
          else
          {
            prev_loc->next = after;
          }
        }
        else
        {
          before->next = after;
          prev_loc = before;
        }        
        before = after;
      }
    }
    seq_offset += dsp_len;
  }

  slp = SeqLocFree (slp);
  return loc_list;  
}


static SeqLocPtr 
RemoveGapsFromLocation 
(SeqLocPtr slp, 
 Uint4     options,
 Boolean   set_partial_ends,
 BoolPtr   split)
{
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  
  if (slp == NULL)
  {
    return slp;
  }
  
  sip = SeqLocId (slp);
  if (sip == NULL)
  {
    return slp;
  }
  
  bsp = BioseqFind (sip);
  if (bsp == NULL)
  {
    return slp;
  }
  else if (bsp->repr == Seq_repr_seg)
  {
    if (!(options & eAdjustFeatForGap_split_internal)) {
      return slp;
    } else {
      return RemoveGapsFromSegmentedLocation (slp, bsp);
    }
  }
  else if (bsp->repr == Seq_repr_delta
          && bsp->seq_ext_type == 4
          && bsp->seq_ext != NULL)
  {
    return RemoveGapsFromDeltaLocation (slp, bsp, options, set_partial_ends, split);
  }
  else
  {
    return slp;
  }
}

NLM_EXTERN void AdjustFrame (SeqFeatPtr sfp, BioseqPtr oldprot)
{
  ByteStorePtr bs;
  CdRegionPtr  crp;
  CharPtr      oldprot_str, newprot_str;
  Uint1        orig_frame, best_frame = 0;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || sfp->data.value.ptrvalue == NULL
      || oldprot == NULL)
  {
    return;
  }
  
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  
  oldprot_str = GetSequenceByBsp (oldprot);
  
  orig_frame = crp->frame;
  for (crp->frame = 1; crp->frame <= 3 && best_frame == 0; crp->frame++)
  {
    newprot_str = NULL;
    bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    if (bs != NULL) 
    {
      newprot_str = BSMerge (bs, NULL);
      bs = BSFree (bs);
    }
    
    if (StringHasNoText (newprot_str))
    {
      newprot_str = MemFree (newprot_str);
      continue;
    }
    if (StringSearch (oldprot_str, newprot_str) != NULL
        || StringSearch (oldprot_str, newprot_str + 1) != NULL)
    {
      best_frame = crp->frame;
    }
    else
    {
      newprot_str [StringLen (newprot_str) - 1] = 0;
      if (StringSearch (oldprot_str, newprot_str) != NULL
          || StringSearch (oldprot_str, newprot_str + 1) != NULL)
      {
        best_frame = crp->frame;
      }
    }
    newprot_str = MemFree (newprot_str);
  }
  oldprot_str = MemFree (oldprot_str);
  if (best_frame > 0)
  {
    crp->frame = best_frame;
  }
  else
  {
    crp->frame = orig_frame;
  }

}


NLM_EXTERN AdjustFeatForGapPtr AdjustFeatForGapFree (AdjustFeatForGapPtr agp)
{
  if (agp != NULL) {
    agp->feature_list = ValNodeFree (agp->feature_list);
    agp->features_in_gap = ValNodeFree (agp->features_in_gap);
    agp = MemFree (agp);
  }
  return agp;
}
 

NLM_EXTERN Boolean FeatureOkForFeatureList (SeqFeatPtr sfp, ValNodePtr feature_list)
{
  ValNodePtr vnp;
  Boolean    rval = FALSE;

  if (sfp == NULL) return FALSE;
  if (feature_list == NULL) return TRUE;
  for (vnp = feature_list; vnp != NULL && !rval; vnp = vnp->next) {
    if (vnp->choice == 255 || vnp->choice == sfp->idx.subtype) {
      rval = TRUE;
    }
  }
  return rval;
}


NLM_EXTERN SeqFeatPtr GetGeneForFeature (SeqFeatPtr sfp)
{
  BioseqPtr bsp;
  GeneRefPtr grp;
  SeqFeatPtr overlap_gene = NULL;
  Boolean    is_suppressed;
  SeqMgrFeatContext fcontext;

  grp = SeqMgrGetGeneXref (sfp);
  is_suppressed = SeqMgrGeneIsSuppressed (grp);
  if (is_suppressed) return NULL;

  if (grp != NULL) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp == NULL) return NULL;
    if (StringDoesHaveText (grp->locus_tag)) {
      overlap_gene = SeqMgrGetGeneByLocusTag (bsp, grp->locus_tag, &fcontext);
    } else if (StringDoesHaveText (grp->locus)) {
      overlap_gene = SeqMgrGetFeatureByLabel (bsp, grp->locus, SEQFEAT_GENE, 0, &fcontext);
    }
  } else {
    overlap_gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  }
  return overlap_gene;
}


const char *cds_gap_comment = "coding region disrupted by sequencing gap";

NLM_EXTERN void AddCDSGapComment (SeqFeatPtr sfp)
{
  CharPtr new_comment = NULL;
  
  if (sfp == NULL || StringSearch (sfp->comment, cds_gap_comment) != NULL)
  {
    return;
  }
  
  if (StringHasNoText (sfp->comment))
  {
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = StringSave (cds_gap_comment);
  }
  else
  {
    new_comment = (CharPtr) MemNew ((StringLen (sfp->comment) 
                                     + StringLen (cds_gap_comment)
                                     + 4) * sizeof (Char));
    StringCpy (new_comment, sfp->comment);
    StringCat (new_comment, "; ");
    StringCat (new_comment, cds_gap_comment);
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;
  }
}


NLM_EXTERN BioseqPtr 
AddProteinSequenceCopy 
(BioseqPtr  protbsp, 
 BioseqPtr  featbsp,
 SeqFeatPtr new_sfp,
 Uint2      entityID)
{
  Char        str [128];
  SeqIdPtr    new_id;
  BioseqPtr   new_protbsp;
  SeqEntryPtr prot_sep, parent_sep;
  SeqFeatPtr  prot_sfp, prot_cpy;
  SeqAnnotPtr sap;
  SeqLocPtr   slp;

  if (protbsp == NULL || featbsp == NULL || new_sfp == NULL)
  {
    return NULL;
  }
  
  parent_sep = GetBestTopParentForData (entityID, featbsp);
  if (parent_sep == NULL 
      || !IS_Bioseq_set (parent_sep) 
      || parent_sep->data.ptrvalue == NULL)
  {
    return NULL;
  }
  
  SeqIdWrite (protbsp->id, str, PRINTID_REPORT, sizeof (str) - 1);    
  new_id = MakeUniqueSeqID (str);                                            
  new_protbsp = BioseqCopyEx (new_id, protbsp, 0, protbsp->length - 1,
                                      Seq_strand_plus, TRUE);
  ValNodeLink (&(new_protbsp->descr),
               AsnIoMemCopy ((Pointer) protbsp->descr,
                             (AsnReadFunc) SeqDescrAsnRead,
                             (AsnWriteFunc) SeqDescrAsnWrite));
                                      
  prot_sep = SeqEntryNew ();
  if (prot_sep != NULL)
  {
    prot_sep->choice = 1;
    prot_sep->data.ptrvalue = new_protbsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) new_protbsp, prot_sep);
    AddSeqEntryToSeqEntry (parent_sep, prot_sep, TRUE);
  }
  
  SeqMgrAddToBioseqIndex (new_protbsp);           
  new_sfp->product = SeqLocWholeNew (new_protbsp); 

  if (new_protbsp->annot == NULL) {
    for (sap = protbsp->annot; sap != NULL; sap = sap->next) {
      if (sap->type == 1) {
        for (prot_sfp = sap->data; prot_sfp != NULL; prot_sfp = prot_sfp->next) {
          slp = SeqLocCopy (prot_sfp->location);
          slp = SeqLocReplaceID (slp, new_protbsp->id);
          prot_cpy = CreateNewFeatureOnBioseq (new_protbsp, SEQFEAT_PROT, slp);
          prot_cpy->data.value.ptrvalue = AsnIoMemCopy (prot_sfp->data.value.ptrvalue, (AsnReadFunc) ProtRefAsnRead, (AsnWriteFunc) ProtRefAsnWrite);
        }
      }
    }
  }
  SeqMgrIndexFeatures (entityID, NULL);
  return new_protbsp;
}


NLM_EXTERN void SetProductSequencePartials (BioseqPtr protbsp, Boolean partial5, Boolean partial3)
{
  SeqFeatPtr  prot_sfp;
  SeqMgrFeatContext context;
  SeqMgrDescContext dcontext;
  SeqDescrPtr       sdp;
  MolInfoPtr        mip;

  if (protbsp == NULL)
  {
    return;
  }
        
  /* set partials on product */
  prot_sfp = SeqMgrGetNextFeature (protbsp, NULL, 
                                          SEQFEAT_PROT, FEATDEF_PROT, &context);
  if (prot_sfp != NULL)
  {
    SetSeqLocPartial (prot_sfp->location, partial5, partial3);
    prot_sfp->partial = partial5 || partial3;
  }
        
  sdp = SeqMgrGetNextDescriptor (protbsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL)
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (partial5 && partial3) {
      mip->completeness = 5;
    } else if (partial5) {
      mip->completeness = 3;
    } else if (partial3) {
      mip->completeness = 4;
    } else if (prot_sfp != NULL && prot_sfp->partial) {
      mip->completeness = 2;
    } else {
      mip->completeness = 0;
    }
  }  
}


static SeqLocPtr MakeMixedLocFromLocList (SeqLocPtr loc_list)
{
  SeqLocPtr   slp_mix, tmp_slp, prev_slp = NULL, next_slp;

  /* make location mixed */    
  slp_mix = ValNodeNew(NULL);
  slp_mix->choice = SEQLOC_MIX;
  slp_mix->data.ptrvalue = loc_list;
  tmp_slp = loc_list;
  while (tmp_slp != NULL) {
    next_slp = tmp_slp->next;
    if (tmp_slp->choice == SEQLOC_MIX) {
      if (tmp_slp->data.ptrvalue == NULL) {
        /* empty mixed loc, just remove it */
        if (prev_slp == NULL) {
          slp_mix->data.ptrvalue = tmp_slp->next;
        } else {
          prev_slp->next = tmp_slp->next;
        }
        tmp_slp->next = NULL;
        tmp_slp = SeqLocFree (tmp_slp);
      } else {
        /* take sublocations and promote them */
        if (prev_slp == NULL) {
          slp_mix->data.ptrvalue = tmp_slp->data.ptrvalue;
        } else {
          prev_slp->next = tmp_slp->data.ptrvalue;
        }
        prev_slp = tmp_slp->data.ptrvalue;
        while (prev_slp->next != NULL) {
          prev_slp = prev_slp->next;
        }
        prev_slp->next = next_slp;
        tmp_slp->next = NULL;
        tmp_slp->data.ptrvalue = NULL;
        tmp_slp = SeqLocFree (tmp_slp);
      }
    } else {
      prev_slp = tmp_slp;
    }
    tmp_slp = next_slp;
  }
  return slp_mix;
}


NLM_EXTERN void AdjustFeatureForGapsCallback (SeqFeatPtr sfp, Pointer data)
{
  AdjustFeatForGapPtr afgp;
  BioseqPtr   protbsp = NULL, new_protbsp;
  SeqFeatPtr  new_sfp, tmp, gene, mrna;
  Boolean     partial5, partial3;
  Uint2       entityID;
  SeqLocPtr   slp_new;
  Boolean     split;
  Boolean     split_internal;
  ValNodePtr  tmp_list;
  BioseqPtr   gapped_bioseq;
  SeqMgrFeatContext fcontext;
  Boolean           set_partial_ends;

  if (sfp == NULL || data == NULL || sfp->location == NULL || sfp->idx.deleteme) return;

  afgp = (AdjustFeatForGapPtr) data;

  if (!FeatureOkForFeatureList(sfp, afgp->feature_list)) return;

  gapped_bioseq = BioseqFind (SeqLocId (sfp->location));
  
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  slp_new = (SeqLocPtr) AsnIoMemCopy (sfp->location, 
                                      (AsnReadFunc) SeqLocAsnRead, 
                                      (AsnWriteFunc) SeqLocAsnWrite);                                         
  split = FALSE;
  set_partial_ends = afgp->options & eAdjustFeatForGap_make_partial;
  if (set_partial_ends && ! (afgp->options & eAdjustFeatForGap_partial_for_pseudo)) {
    if (sfp->pseudo) {
      set_partial_ends = FALSE;
    } else {
      gene = GetGeneForFeature (sfp);
      if (gene != NULL && gene->pseudo) {
        set_partial_ends = FALSE;
      }
    }
  }

  /* handle overlapping features if coding region */
  if (sfp->data.choice == SEQFEAT_CDREGION) 
  {    
    split_internal = afgp->options & eAdjustFeatForGap_split_internal;
    afgp->options &= ~eAdjustFeatForGap_split_internal;
    tmp_list = afgp->feature_list;
    afgp->feature_list = NULL;
    gene = GetGeneForFeature (sfp);
    AdjustFeatureForGapsCallback (gene, afgp);
    mrna = SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
    afgp->options |= split_internal;
    AdjustFeatureForGapsCallback (mrna, afgp);
    afgp->feature_list = tmp_list;
  }

  slp_new = RemoveGapsFromLocation (slp_new, afgp->options, set_partial_ends, &split);
  if (slp_new == NULL) {
    ValNodeAddPointer (&(afgp->features_in_gap), OBJ_SEQFEAT, sfp);
    return;
  } else if (SeqLocCompare (slp_new, sfp->location) == SLC_A_EQ_B) {
    slp_new = SeqLocFree (slp_new);
    return;
  }

  if (split) {
    AddCDSGapComment(sfp);
  }

  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp_new;

  if (gapped_bioseq != NULL && gapped_bioseq->repr == Seq_repr_delta)
  {
    entityID = gapped_bioseq->idx.entityID;
   
    while (sfp->location->next != NULL)
    {
      /* create a copy of the feature for each interval between gaps */
      tmp = sfp->next;
      sfp->next = NULL;
      new_sfp = (SeqFeatPtr) AsnIoMemCopy (sfp, 
                                          (AsnReadFunc) SeqFeatAsnRead, 
                                          (AsnWriteFunc) SeqFeatAsnWrite);
      sfp->next = new_sfp;
      new_sfp->next = tmp;

      new_sfp->location = sfp->location->next;
      new_sfp->comment = StringSave (sfp->comment);

      sfp->location->next = NULL;
        
      /* fix partials */
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      sfp->partial = partial5 || partial3;

      if (sfp->data.choice == SEQFEAT_CDREGION) {
        protbsp = BioseqFindFromSeqLoc (sfp->product);        
        if (protbsp != NULL)
        {
          new_protbsp = AddProteinSequenceCopy (protbsp, gapped_bioseq, new_sfp, entityID);                                                  
        }
                
        /* adjust frame */
        AdjustFrame (sfp, protbsp);

        /* retranslate coding region */
        SeqEdTranslateOneCDS (sfp, gapped_bioseq, entityID, afgp->align_func);
        
        /* set partials on product */
        if (protbsp == NULL)
        {
          protbsp = BioseqFindFromSeqLoc (sfp->product);
        }
        SetProductSequencePartials (protbsp, partial5, partial3);
        protbsp = new_protbsp;
      }    
      partial5 = TRUE;
      sfp = new_sfp;
    }
    /* fix partials for last feature */
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    sfp->partial = partial5 || partial3;

    if (sfp->data.choice == SEQFEAT_CDREGION) 
    {    
      /* adjust frame */
      protbsp = BioseqFindFromSeqLoc (sfp->product);
      AdjustFrame (sfp, protbsp);

      /* retranslate coding region */
      SeqEdTranslateOneCDS (sfp, gapped_bioseq, entityID, afgp->align_func);
      /* set partials on product */
      SetProductSequencePartials (protbsp, partial5, partial3);
    }
  }
  else 
  {
    /* make location mixed */
    sfp->location = MakeMixedLocFromLocList (sfp->location);
  }
}


NLM_EXTERN void MarkFeaturesInGapsForDeletion (AdjustFeatForGapPtr afgp)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  BioseqPtr  product_bsp;

  if (afgp == NULL || afgp->features_in_gap == NULL)
  {
    return;
  }
  
  for (vnp = afgp->features_in_gap; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    sfp->idx.deleteme = TRUE;
    if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL)
    {
      product_bsp = BioseqFindFromSeqLoc (sfp->product);
      if (product_bsp != NULL)
      {
        product_bsp->idx.deleteme = TRUE;
      }
    }
  }
  afgp->features_in_gap = ValNodeFree (afgp->features_in_gap);
}

NLM_EXTERN void AdjustCDSLocationsForUnknownGapsCallback (SeqFeatPtr sfp, Pointer data)
{
  AdjustFeatForGapData agd;

  agd.feature_list = NULL;
 
  agd.options = eAdjustFeatForGap_unknown_gaps | eAdjustFeatForGap_make_partial
                | eAdjustFeatForGap_split_internal | eAdjustFeatForGap_trim_ends;

  agd.align_func = (GlobalAlignFunc) data;

  agd.features_in_gap = NULL;

  AdjustFeatureForGapsCallback (sfp, &agd);
  MarkFeaturesInGapsForDeletion (&agd);
}


NLM_EXTERN void RevCompOneFeatForBioseq (SeqFeatPtr sfp, BioseqPtr bsp)
{
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  CodeBreakPtr cbp;
  CdRegionPtr  crp;
  RnaRefPtr    rrp;
  tRNAPtr      trp;
  Boolean      split;

  if (sfp == NULL || bsp == NULL) return;

  sip = SeqLocId (sfp->location);
  if (sip != NULL) {
    if (SeqIdIn (sip, bsp->id)) {
      slp = SeqLocCopyRegion (sip, sfp->location, bsp, 0,
                              bsp->length - 1, Seq_strand_minus, &split);
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp;
      switch (sfp->data.choice) {
        case SEQFEAT_CDREGION :
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              sip = SeqLocId (cbp->loc);
              slp = SeqLocCopyRegion (sip, cbp->loc, bsp, 0,
                                      bsp->length - 1, Seq_strand_minus, &split);
              cbp->loc = SeqLocFree (cbp->loc);
              cbp->loc = slp;
            }
          }
          break;
        case SEQFEAT_RNA :
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->ext.choice == 2) {
            trp = (tRNAPtr) rrp->ext.value.ptrvalue;
            if (trp != NULL && trp->anticodon != NULL) {
              sip = SeqLocId (trp->anticodon);
              slp = SeqLocCopyRegion (sip, trp->anticodon, bsp, 0,
                                      bsp->length - 1, Seq_strand_minus, &split);
              trp->anticodon = SeqLocFree (trp->anticodon);
              trp->anticodon = slp;
            }
          }
          break;
        default :
          break;
      }
    }
  }
}


static void RevCompFeatsOnBioseq (BioseqPtr bsp)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    RevCompOneFeatForBioseq (sfp, bsp);
  }
}


static Boolean s_IsSkippable (Char ch)
{
  if (isspace (ch) || ch == ',' || ch == '"') {
    return TRUE;
  } else {
    return FALSE;
  }
}


static ValNodePtr MakeTokensFromLineAnySeparator (CharPtr line)
{
  CharPtr token_start, token_end, token;
  ValNodePtr tokens = NULL;
  Int4       len;

  if (StringHasNoText (line)) {
    return NULL;
  }

  token_start = line;
  while (*token_start != 0) {
    while (s_IsSkippable (*token_start)) {
      token_start ++;
    }
    if (*token_start != 0) {
      token_end = token_start + 1;
      while (*token_end != 0 && !s_IsSkippable (*token_end)) {
        token_end ++;
      }
      len = token_end - token_start + 1;
      token = (CharPtr) MemNew (sizeof (Char) * len);
      StringNCpy (token, token_start, len - 1);
      token[len - 1] = 0;
      ValNodeAddPointer (&tokens, 0, token);
      token_start = token_end;
    }
  }
  return tokens;
}


static ValNodePtr MakeTokensFromLineTab (CharPtr line)
{
  CharPtr token_start, token;
  ValNodePtr tokens = NULL;
  Int4       len;

  if (StringHasNoText (line)) {
    return NULL;
  }

  token_start = line;
  while (*token_start != 0) {
    len = StringCSpn (token_start, "\t");
    token = (CharPtr) MemNew (sizeof (Char) * (len + 1));
    StringNCpy (token, token_start, len);
    token[len] = 0;
    ValNodeAddPointer (&tokens, 0, token);
    token_start += len;
    if (*token_start == '\t') {
      token_start++;
    }
  }
  return tokens;
}


NLM_EXTERN ValNodePtr MakeTokensFromLine (CharPtr line)
{
  ValNodePtr tokens = NULL;

  if (StringChr (line, '\t') == NULL) {
    tokens = MakeTokensFromLineAnySeparator (line);
  } else {
    tokens = MakeTokensFromLineTab (line);
  }
  return tokens;
}


static ValNodePtr MakeTranscriptomeIDTokensFromLine (CharPtr line)
{
  ValNodePtr tokens = NULL, vnp, vnp_next, tmp, prev;

  if (StringChr (line, '\t') == NULL) {
    tokens = MakeTokensFromLineAnySeparator (line);
  } else {
    tokens = MakeTokensFromLineTab (line);
    if (tokens != NULL && tokens->next != NULL) {
      prev = tokens;
      for (vnp = tokens->next; vnp != NULL; vnp = vnp_next) {
        vnp_next = vnp->next;
        tmp = MakeTokensFromLineAnySeparator (vnp->data.ptrvalue);
        if (tmp != NULL) {
          ValNodeLink (&tmp, vnp->next);
          prev->next = tmp;
          vnp->next = NULL;
          vnp = ValNodeFreeData (vnp);
        } else {
          prev = vnp;
        }
      }
    }
  }
  return tokens;

}


NLM_EXTERN Boolean HasExistingSeqHistAssembly (ValNodePtr list)
{
  TranscriptomeIdsPtr t;
  Boolean has_tables = FALSE;

  while (list != NULL && !has_tables) {
    t = list->data.ptrvalue;
    if (t != NULL && t->consensus_bsp != NULL 
        && t->consensus_bsp->hist != NULL 
        && t->consensus_bsp->hist->assembly != NULL) {
      has_tables = TRUE;
    }
    list = list->next;
  }
  return has_tables;
}


NLM_EXTERN void DeleteSeqHistAssembliesForList (ValNodePtr list)
{
  TranscriptomeIdsPtr t;

  while (list != NULL) {
    t = (TranscriptomeIdsPtr) list->data.ptrvalue;
    if (t != NULL && t->consensus_bsp != NULL && t->consensus_bsp->hist != NULL
        && t->consensus_bsp->hist->assembly != NULL) {
      t->consensus_bsp->hist->assembly = SeqAlignFree (t->consensus_bsp->hist->assembly);
    }
    list = list->next;
  }
}


static void AddTSARangeError (ValNodePtr PNTR range_list, CharPtr id, Int4 start, Int4 stop)
{
  CharPtr      big_range_fmt = "%s: Large gap in coverage (>50) from %d to %d";
  CharPtr      med_range_fmt = "%s: Medium gap in coverage (10-50) from %d to %d";
  CharPtr      small_range_fmt = "%s: Small gap in coverage (<10) from %d to %d";
  CharPtr      fmt;
  CharPtr      range; 
  Int4         diff;

  diff = stop - start + 1;
  if (diff > 50) {
    fmt = big_range_fmt;
  } else if (diff < 10) {
    fmt = small_range_fmt;
  } else {
    fmt = med_range_fmt;
  }

  range = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (id) + 30));
  sprintf (range, fmt, id, start, stop);
  ValNodeAddPointer (range_list, 0, range);
}


NLM_EXTERN ValNodePtr ReportCoverageForBioseqSeqHist (BioseqPtr bsp)
{
  ValNodePtr   range_list = NULL;
  SeqAlignPtr  salp;
  Int4         assembly_from, assembly_to, tmp, zero_start, primary_from, primary_to;
  Int4         aln_pos, i;
  Int4Ptr      coverage;
  Char         id_buf[255];
  Char         id_buf2[255];
  CharPtr      err_msg;
  CharPtr      no_assembly_fmt = "Consensus sequence %s has no assembly";
  CharPtr      gaps_fmt = "Too many gaps in alignment between %s and %s";
  Int4         assem_len, prim_len;
  CharPtr      range;

  if (bsp == NULL) return NULL;

  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
  if (bsp->hist == NULL || bsp->hist->assembly == NULL) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_assembly_fmt) + StringLen (id_buf)));
    sprintf (err_msg, no_assembly_fmt, id_buf);
    ValNodeAddPointer (&range_list, 0, err_msg);
  } else {
    coverage = (Int4Ptr) MemNew (sizeof (Int4) * bsp->length);
    MemSet (coverage, 0, sizeof (Int4) * bsp->length);

    for (salp = bsp->hist->assembly; salp != NULL; salp = salp->next) {
      AlnMgr2GetNthSeqRangeInSA(salp, 1, &assembly_from, &assembly_to);
      for (i = assembly_from; i <= assembly_to; i++) {
        if (coverage[i] == 0) {
          aln_pos = AlnMgr2MapBioseqToSeqAlign(salp, i, 1);
          if (AlnMgr2MapSeqAlignToBioseq (salp, aln_pos, 1) > -1) {
            coverage[i] = 1;
          }
        }
      }
      AlnMgr2GetNthSeqRangeInSA (salp, 2, &primary_from, &primary_to);
      if (assembly_to > assembly_from) {
        assem_len = assembly_to - assembly_from + 1;
      } else {
        assem_len = assembly_from - assembly_to + 1;
      }
      if (primary_to > primary_from) {
        prim_len = primary_to - primary_from + 1;
      } else {
        prim_len = primary_from - primary_to + 1;
      }
      if (ABS(assem_len - prim_len) >= .1 * assem_len || ABS(assem_len - prim_len) >= .1 * prim_len) {
        SeqIdWrite (AlnMgr2GetNthSeqIdPtr (salp, 2), id_buf2, PRINTID_REPORT, sizeof (id_buf2) - 1);
        range = (CharPtr) MemNew (sizeof (Char) * (StringLen (gaps_fmt) + StringLen (id_buf) + StringLen (id_buf2)));
        sprintf (range, gaps_fmt, id_buf, id_buf2);
        ValNodeAddPointer (&range_list, 0, range);
      } else if (ABS(assem_len - prim_len) > 50) {
        SeqIdWrite (AlnMgr2GetNthSeqIdPtr (salp, 2), id_buf2, PRINTID_REPORT, sizeof (id_buf2) - 1);
        range = (CharPtr) MemNew (sizeof (Char) * (StringLen (gaps_fmt) + StringLen (id_buf) + StringLen (id_buf2)));
        sprintf (range, gaps_fmt, id_buf, id_buf2);
        ValNodeAddPointer (&range_list, 0, range);
      }
    }

    zero_start = -1;
    for (tmp = 0; tmp < bsp->length; tmp++) {
      if (coverage[tmp] == 0) {
        if (zero_start == -1) {
          zero_start = tmp;
        }
      } else if (zero_start > -1) {
        /* note - print values as 1-based rather than 0-based coordinates.  Second value is actually tmp - 1 + 1 */
        AddTSARangeError (&range_list, id_buf, zero_start + 1, tmp);
        zero_start = -1;
      }
    }
    if (coverage[bsp->length - 1] == 0) {
      /* note - print values as 1-based rather than 0-based coordinates.  Second value is actually bsp->length - 1 + 1 */
      AddTSARangeError (&range_list, id_buf, zero_start + 1, bsp->length);
    }
    coverage = MemFree (coverage);
  }
  return range_list;
}


static Boolean IsSingleEstCoverage (BioseqPtr bsp)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip;

  if (bsp == NULL || bsp->hist == NULL || bsp->hist->assembly == NULL) {
    return FALSE;
  }
  if (bsp->hist->assembly->next == NULL) {
    return TRUE;
  } else {
    sip = AlnMgr2GetNthSeqIdPtr (bsp->hist->assembly, 2);
    for (salp = bsp->hist->assembly->next; salp != NULL; salp = salp->next) {
      if (SeqIdComp (AlnMgr2GetNthSeqIdPtr (salp, 2), sip) == SIC_NO) {
        return FALSE;
      }
    }
    return TRUE;
  }
}


NLM_EXTERN ValNodePtr ReportCoverageForTranscriptomeIdsListSeqHist (ValNodePtr list)
{
  ValNodePtr range_list = NULL, new_list;
  TranscriptomeIdsPtr t;
  Char                id_str[255];
  CharPtr             good_fmt = "Coverage is complete for %s";
  CharPtr             msg;
  Int4                num_single = 0, num_mult = 0;
  CharPtr             single_fmt = "%d records are covered by a single EST";
  CharPtr             mult_fmt = "%d records are covered by multiple ESTs";

  while (list != NULL) {
    t = (TranscriptomeIdsPtr) list->data.ptrvalue;
    if (t != NULL && t->consensus_bsp != NULL) {
      new_list = ReportCoverageForBioseqSeqHist (t->consensus_bsp);
      if (new_list == NULL) {
        SeqIdWrite (SeqIdFindBest (t->consensus_bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
        msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (good_fmt) + StringLen (id_str)));
        sprintf (msg, good_fmt, id_str);
        ValNodeAddPointer (&range_list, 0, msg);
      } else {
        ValNodeLink (&range_list, new_list);
      }
      if (IsSingleEstCoverage(t->consensus_bsp)) {
        num_single++;
      } else {
        num_mult++;
      }
    }
    list = list->next;
  }
  msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (single_fmt) + 15));
  sprintf (msg, single_fmt, num_single);
  ValNodeAddPointer (&range_list, 0, msg);
  msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (mult_fmt) + 15));
  sprintf (msg, mult_fmt, num_mult);
  ValNodeAddPointer (&range_list, 0, msg);
  
  return range_list;
}


NLM_EXTERN ValNodePtr ReportGapsInSeqHistAlignmentForBioseq (BioseqPtr bsp)
{
  ValNodePtr   err_list = NULL;
  SeqAlignPtr  salp;
  DenseSegPtr  dsp;
  Char         id_buf[255];
  Char         id_buf2[255];
  CharPtr      err_msg;
  CharPtr      gaps_fmt = "Alignment between %s and %s contains gaps";
  Int4         j, max_gap;
  Boolean      has_gaps;

  if (bsp == NULL || bsp->hist == NULL || bsp->hist->assembly == NULL) return NULL;

  max_gap = bsp->length / 10;
  if (max_gap > 50) {
    max_gap = 50;
  }
  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
  for (salp = bsp->hist->assembly; salp != NULL; salp = salp->next) {
    if (salp->segtype == SAS_DENSEG && (dsp = (DenseSegPtr) salp->segs) != NULL && dsp->numseg > 1) {
      has_gaps = FALSE;
      for (j = 0; j < dsp->numseg && !has_gaps; j++) {
        if (dsp->lens[j] >= max_gap) {
          if (dsp->starts[j * dsp->dim] == -1 || dsp->starts[(j * dsp->dim) + 1] == -1) {
            has_gaps = TRUE;
          }
        }
      }
      if (has_gaps) {
        SeqIdWrite (AlnMgr2GetNthSeqIdPtr (salp, 2), id_buf2, PRINTID_REPORT, sizeof (id_buf2) - 1);          
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (gaps_fmt) + StringLen (id_buf) + StringLen (id_buf2)));
        sprintf (err_msg, gaps_fmt, id_buf, id_buf2);
        ValNodeAddPointer (&err_list, 0, err_msg);
      }
    }        
  }
  return err_list;
}


NLM_EXTERN ValNodePtr ReportGapsInSeqHistAlignmentsForIdsList (ValNodePtr list)
{
  TranscriptomeIdsPtr t;
  ValNodePtr          errs = NULL;

  while (list != NULL) {
    t = (TranscriptomeIdsPtr) list->data.ptrvalue;
    if (t != NULL && t->consensus_bsp != NULL) {
      ValNodeLink (&errs, ReportGapsInSeqHistAlignmentForBioseq (t->consensus_bsp));
    }
    list = list->next;
  }
  return errs;
}


static ValNodePtr GetSeqHistAlignmentSummaryForRange (Int4 cons_start, Int4 cons_stop, SeqAlignPtr salp_list)
{
  Int4 aln_start, aln_stop;
  Uint1Ptr buf1 = NULL, buf2 = NULL;
  Int4 buf_size = -1;
  Int4 aln_len = 0;
  SeqIdPtr sip;
  Char id1[255], id2[255];
  CharPtr aln_msg;
  ValNodePtr summary = NULL;
  Boolean    show_consensus = TRUE;

  if (salp_list == NULL || cons_start < 0 || cons_stop < cons_start) {
    return NULL;
  }

  while (salp_list != NULL) {
    sip = AlnMgr2GetNthSeqIdPtr (salp_list, 1);
    SeqIdWrite (sip, id1, PRINTID_REPORT, sizeof (id1) - 1);
    sip = SeqIdFree (sip);
    sip = AlnMgr2GetNthSeqIdPtr (salp_list, 2);
    SeqIdWrite (sip, id2, PRINTID_REPORT, sizeof (id2) - 1);
    sip = SeqIdFree (sip);

    aln_start = AlnMgr2MapBioseqToSeqAlign(salp_list, cons_start, 1);
    aln_stop = AlnMgr2MapBioseqToSeqAlign(salp_list, cons_stop, 1);

    if (aln_start >= 0 && aln_stop >= 0) {    
      if (buf_size < aln_stop - aln_start + 2) {
        buf1 = MemFree (buf1);
        buf2 = MemFree (buf2);
        buf_size = aln_stop - aln_start + 2;
        buf1 = (Uint1Ptr) MemNew (sizeof (Uint1) * buf_size);
        buf2 = (Uint1Ptr) MemNew (sizeof (Uint1) * buf_size);
      }
      if (show_consensus) {
        AlignmentIntervalToString (salp_list, 1, aln_start, aln_stop, 1, FALSE, 
                                    buf1, buf2, &aln_len, FALSE);
        aln_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (id1) + aln_len + 3));
        StringCpy (aln_msg, id1);
        StringCat (aln_msg, " ");
        StringNCat (aln_msg, (CharPtr) buf2, aln_len);
        ValNodeAddPointer (&summary, 0, aln_msg);
        if (aln_start == aln_stop) {
          show_consensus = FALSE;
        }
      }
      AlignmentIntervalToString (salp_list, 2, aln_start, aln_stop, 2, FALSE, 
                                  buf1, buf2, &aln_len, FALSE);
      aln_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (id2) + aln_len + 3));
      StringCpy (aln_msg, id2);
      StringCat (aln_msg, " ");
      StringNCat (aln_msg, (CharPtr) buf2, aln_len);
      if (show_consensus) {
        StringCat (aln_msg, "\n");
      }
      ValNodeAddPointer (&summary, 0, aln_msg);
    }
    salp_list = salp_list->next;
  }
  buf1 = MemFree (buf1);
  buf2 = MemFree (buf2);
  return summary;  
}


NLM_EXTERN ValNodePtr ReportConsensusMatchForBioseqSeqHist (BioseqPtr bsp)
{
  ValNodePtr   err_list = NULL;
  SeqAlignPtr  salp;
  Int4         assembly_from = 0, assembly_to = 0;
  Int4         aln_pos, i, read_pos, pct;
  Int4Ptr      coverage;
  Int4Ptr      match;
  Char         id_buf[255];
  CharPtr      err_msg;
  CharPtr      err_fmt = "Consensus sequence %s matches less than half of reads at position %d";
  CharPtr      err_range_fmt = "Consensus sequence %s matches less than half of reads at positions %d-%d";
  SeqIdPtr     sip;
  BioseqPtr    read_bsp;
  Char         buf1[2], buf2[2];
  Uint1        read_strand;
  Int4         start_range = -1;

  if (bsp == NULL || bsp->hist == NULL || bsp->hist->assembly == NULL) return NULL;

  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
  coverage = (Int4Ptr) MemNew (sizeof (Int4) * bsp->length);
  MemSet (coverage, 0, sizeof (Int4) * bsp->length);
  match = (Int4Ptr) MemNew (sizeof (Int4) * bsp->length);
  MemSet (match, 0, sizeof (Int4) * bsp->length);

  for (salp = bsp->hist->assembly; salp != NULL; salp = salp->next) {
    sip = AlnMgr2GetNthSeqIdPtr (salp, 2);
    read_bsp = BioseqLockById (sip);
    sip = SeqIdFree (sip);
    read_strand = SeqAlignStrand (salp, 1);
    AlnMgr2GetNthSeqRangeInSA(salp, 1, &assembly_from, &assembly_to);
    for (i = assembly_from; i <= assembly_to; i++) {
      aln_pos = AlnMgr2MapBioseqToSeqAlign(salp, i, 1);
      if ((read_pos = AlnMgr2MapSeqAlignToBioseq (salp, aln_pos, 2)) > -1) {
        coverage[i] ++;
        SeqPortStreamInt (bsp, i, i, Seq_strand_plus, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf1, NULL);
        SeqPortStreamInt (read_bsp, read_pos, read_pos, read_strand, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) buf2, NULL);
        if (buf1[0] == buf2[0]) {
          match[i] ++;
        }
      }
    }
    BioseqUnlock (read_bsp);
  }

  for (i = assembly_from; i <= assembly_to; i++) {
    if (coverage[i] > 0) {
      pct = (100 * match[i]) / coverage[i];
      if (pct < 50) {
        if (start_range < 0) {
          start_range = i;
        }
      } else {
        if (start_range > -1) {
          if (i > start_range + 1) {
            err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (err_range_fmt) + StringLen (id_buf) + 30));
            sprintf (err_msg, err_range_fmt, id_buf, start_range + 1, i);
          } else {
            err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (err_fmt) + StringLen (id_buf) + 15));
            sprintf (err_msg, err_fmt, id_buf, i);
          }
          ValNodeAddPointer (&err_list, 0, err_msg);
          ValNodeLink (&err_list, GetSeqHistAlignmentSummaryForRange (start_range, i - 1, bsp->hist->assembly));
          start_range = -1;
        }
      }
    }
  }

  if (start_range > -1) {
    if (i > start_range + 1) {
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (err_range_fmt) + StringLen (id_buf) + 30));
      sprintf (err_msg, err_range_fmt, id_buf, start_range + 1, i);
    } else {
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (err_fmt) + StringLen (id_buf) + 15));
      sprintf (err_msg, err_fmt, id_buf, i);
    }
    ValNodeAddPointer (&err_list, 0, err_msg);
    ValNodeLink (&err_list, GetSeqHistAlignmentSummaryForRange (start_range, i - 1, bsp->hist->assembly));
    start_range = -1;
  }

  coverage = MemFree (coverage);
  match = MemFree (match);
  return err_list;
}


static int LIBCALLBACK SortAlignmentByRange (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  SeqAlignPtr salp1, salp2;
  Int4        from1 = -1, from2 = -1, to1 = -1, to2 = -1;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  salp1 = vnp1->data.ptrvalue;
  salp2 = vnp2->data.ptrvalue;
  
  AlnMgr2GetNthSeqRangeInSA(salp1, 1, &from1, &to1);
  AlnMgr2GetNthSeqRangeInSA(salp2, 1, &from2, &to2);

  if (from1 < from2) {
    return -1;
  } else if (from1 > from2) {
    return 1;
  } else if (to1 < to2) {
    return -1;
  } else if (to1 > to2) {
    return 1;
  } else {
    return 0;
  }
}


NLM_EXTERN SeqAlignPtr SortPairwiseAlignmentsByFirstSeqRange (SeqAlignPtr salp)
{
  ValNodePtr list = NULL, vnp;
  SeqAlignPtr salp_tmp, salp_next, salp_prev = NULL;

  if (salp == NULL || salp->next == NULL) {
    return salp;
  }

  for (salp_tmp = salp; salp_tmp != NULL; salp_tmp = salp_next) {
    salp_next = salp_tmp->next;
    salp_tmp->next = NULL;
    ValNodeAddPointer (&list, 0, salp_tmp);
  }
  list = ValNodeSort (list, SortAlignmentByRange);
  salp = list->data.ptrvalue;
  salp_prev = salp;
  for (vnp = list->next; vnp != NULL; vnp = vnp->next) {
    salp_prev->next = vnp->data.ptrvalue;
    salp_prev = salp_prev->next;
  }
  list = ValNodeFree (list);
  return salp;
}
 
 
/* nth is the sequence in the alignment to reverse (1 is first, 2 is second) */
extern void ReverseAlignmentStrand (SeqAlignPtr salp, Int4 nth)
{
  DenseSegPtr dsp;
  SeqIdPtr    sip;
  BioseqPtr   bsp;
  Int4        i, j;
  
  if (salp == NULL || salp->segtype != SAS_DENSEG || salp->segs == NULL)
  {
    return;
  }
  
  dsp = (DenseSegPtr) salp->segs;

  if (dsp->strands == NULL) {
    dsp->strands = (Uint1Ptr) MemNew (dsp->numseg * dsp->dim * sizeof (Uint1));
    MemSet (dsp->strands, Seq_strand_plus, dsp->numseg * dsp->dim * sizeof (Uint1));
  }
  
  sip = dsp->ids;
  i = 1;
  while (sip != NULL && i < nth) {
    sip = sip->next;
    i++;
  }
  bsp = BioseqFind (sip);
  if (bsp == NULL)
  {
    return;
  }
  for (i = 0; i < dsp->numseg; i++)
  {
    j = (i * dsp->dim) + nth - 1;
    
    if (dsp->starts[j] > -1)
    {
      dsp->starts[j] = bsp->length - dsp->starts[j] - dsp->lens[i];
    }
    if (dsp->strands [j] == Seq_strand_minus)
    {
      dsp->strands [j] = Seq_strand_plus;
    }
    else
    {
      dsp->strands [j] = Seq_strand_minus;
    }
  }
  
}


NLM_EXTERN ValNodePtr 
MakeTranscriptomeAssemblySeqHist 
(TranscriptomeIdsPtr t,
 LocalAlignFunc aln_func,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  BioseqPtr   read_bsp;
  ValNodePtr  salp_list = NULL;
  SeqAlignPtr salp, salp_prev;
  ValNodePtr  vnp;
  SeqIdPtr    sip;
  Boolean     dirty;
  ValNodePtr  err_list = NULL;
  CharPtr     err_msg;
  CharPtr     no_aln_fmt = "No alignment between %s and consensus sequence %s";
  CharPtr     invalid_aln_fmt = "Alignment between %s and consensus sequence %s is invalid";
  CharPtr     not_replaced_fmt = "Existing assembly for %s was not replaced";
  CharPtr     no_download_fmt = "Unable to download %s";
  Char        id_buf[255];
  Char        consensus_id_buf[255];
  ErrSev      old_sev;
  
  if (t == NULL || t->consensus_bsp == NULL || t->token_list == NULL) return NULL;

  SeqIdWrite (SeqIdFindBest (t->consensus_bsp->id, SEQID_GENBANK),
              consensus_id_buf, PRINTID_REPORT, sizeof (consensus_id_buf) - 1);

  if (t->consensus_bsp->hist != NULL && t->consensus_bsp->hist->assembly != NULL) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (not_replaced_fmt) + StringLen (consensus_id_buf)));
    sprintf (err_msg, not_replaced_fmt, id_buf);
    ValNodeAddPointer (&err_list, 0, err_msg);
    return err_list;
  }    

 
  for (vnp = t->token_list; vnp != NULL; vnp = vnp->next) {
    if (StringChr (vnp->data.ptrvalue, '|') == NULL) {
      sprintf (id_buf, "gb|%s", (char *) vnp->data.ptrvalue);
    } else {
      sprintf (id_buf, "%s", (char *) vnp->data.ptrvalue);
    }
    sip = MakeSeqID (id_buf);
    read_bsp = BioseqLockById (sip);
    if (read_bsp == NULL) {
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_download_fmt) + StringLen (id_buf)));
      sprintf (err_msg, no_download_fmt, id_buf);
      ValNodeAddPointer (&err_list, 0, err_msg);
    }
    sip = SeqIdFree (sip);

    if (read_bsp == NULL) continue;
    salp = aln_func (t->consensus_bsp, read_bsp);
    
    old_sev = ErrSetMessageLevel (SEV_INFO);
    SeqIdWrite (SeqIdFindBest (read_bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_LONG, sizeof (id_buf) - 1);
    if (salp == NULL) 
    {
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_aln_fmt) + StringLen (id_buf) + StringLen (consensus_id_buf)));
      sprintf (err_msg, no_aln_fmt, id_buf, consensus_id_buf);
      ValNodeAddPointer (&err_list, 0, err_msg);
    }
    else if (! ValidateSeqAlign (salp, t->consensus_bsp->idx.entityID, FALSE, FALSE, TRUE, FALSE, FALSE, &dirty))
    {
      /* if the new alignment wasn't valid, don't add to multiple alignment */
      salp = SeqAlignFree (salp);
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (invalid_aln_fmt) + StringLen (id_buf) + StringLen (consensus_id_buf)));
      sprintf (err_msg, invalid_aln_fmt, id_buf, consensus_id_buf);
      ValNodeAddPointer (&err_list, 0, err_msg);
    }
    else  
    {
      ValNodeAddPointer (&salp_list, OBJ_SEQALIGN, salp);
    }
    ErrSetMessageLevel (old_sev);

    if (change_notify) {
      change_notify (change_userdata);
    }
  }

  salp_list = ValNodeSort (salp_list, SortAlignmentByRange);

  if (salp_list != NULL) {
    if (t->consensus_bsp->hist == NULL) {
      t->consensus_bsp->hist = SeqHistNew ();
    }
    if (t->consensus_bsp->hist->assembly != NULL) {
      t->consensus_bsp->hist->assembly = SeqAlignFree (t->consensus_bsp->hist->assembly);
    }
    t->consensus_bsp->hist->assembly = salp_list->data.ptrvalue;
    salp_prev = t->consensus_bsp->hist->assembly;
    for (vnp = salp_list->next; vnp != NULL; vnp = vnp->next) {
      salp_prev->next = vnp->data.ptrvalue;
      salp_prev = salp_prev->next;
    }
  }

  return err_list;
}


NLM_EXTERN ValNodePtr 
ApplyTranscriptomeIdsListToSeqEntrySeqHist 
(ValNodePtr           list,
 LocalAlignFunc       aln_func,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata)
{
  ValNodePtr err_list = NULL;

  while (list != NULL) {
    ValNodeLink (&err_list, MakeTranscriptomeAssemblySeqHist (list->data.ptrvalue, aln_func, change_notify, change_userdata));
    list = list->next;
  }
  return err_list;
}


NLM_EXTERN TranscriptomeIdsPtr TranscriptomeIdsNew (BioseqPtr bsp, ValNodePtr token_list)
{
  TranscriptomeIdsPtr t;

  t = (TranscriptomeIdsPtr) MemNew (sizeof (TranscriptomeIdsData));
  t->consensus_bsp = bsp;
  t->token_list = token_list;
  return t;
}


NLM_EXTERN TranscriptomeIdsPtr TranscriptomeIdsFree (TranscriptomeIdsPtr t)
{
  if (t != NULL) {
    t->token_list = ValNodeFreeData (t->token_list);
    t = MemFree (t);
  }
  return t;
}


NLM_EXTERN ValNodePtr TranscriptomeIdsListFree (ValNodePtr list)
{
  ValNodePtr list_next;
  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = TranscriptomeIdsFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static BioseqPtr GetTranscriptomeBioseqFromStringId (CharPtr str)
{
  BioseqPtr bsp = NULL, tbsp;
  CharPtr   rev_id;
  Int4      len;
  SeqMgrPtr smp;
  Int4      i, j, num = -1, imin = 0, imax;
  SeqIdIndexElementPtr PNTR sipp;
  SeqEntryPtr               scope;
  Boolean   found;

  if (StringHasNoText (str)) {
    return NULL;
  }

  scope = SeqEntryGetScope();

/*  SeqMgrProcessNonIndexedBioseq(FALSE); */

  /* reverse the string we're looking for */
  len = StringLen (str);
  rev_id = (CharPtr) MemNew (sizeof (Char) * (len + 1));
  for (i = 0; i < len; i++) {
    rev_id[i] = toupper (str[len - i - 1]);
  }
  rev_id[len] = 0;

  smp = SeqMgrReadLock();
  imax = smp->BioseqIndexCnt - 1;
  sipp = smp->BioseqIndex;

  while (imax >= imin) {
    i = (imax + imin)/2;
    if (StringLen (sipp[i]->str) > len && StringNCmp (sipp[i]->str, rev_id, len) == 0) {
      if (sipp[i]->str[len] == '|') {
        num = i;
        break;
      } else {
        /* search down the list */
        j = i - 1;
        found = FALSE;
        while (!found && j > -1 && StringLen (sipp[j]->str) > len 
               && StringNCmp (sipp[j]->str, rev_id, len) == 0) {
          if (sipp[j]->str[len] == '|') {
            found = TRUE;
          } else {
            j--;
          }
        }
        if (found) {
          num = j;
          break;
        }
        /* search up the list */
        j = i + 1;
        while (!found && j < smp->BioseqIndexCnt && StringLen (sipp[j]->str) > len
               && StringNCmp (sipp[j]->str, rev_id, len) == 0) {
          if (sipp[j]->str[len] == '|') {
            found = TRUE;
          } else {
            j++;
          }
        }
        if (found) {
          num = j;
          break;
        } else {
          break;
        }
      }
    } else if ((j = StringCmp (sipp[i]->str, rev_id)) > 0) {
      imax = i - 1;
    } else if (j < 0) {
      imin = i + 1;
    } else {
      num = i;
      break;
    }
  }

  if (num > -1) {
    if (scope == NULL) {
      /* no scope set, take the first one found */
      bsp = sipp[num]->omdp->dataptr;
    } else {
      /* check in scope */
        tbsp = (BioseqPtr)(sipp[num]->omdp->dataptr);
        if (ObjMgrIsChild(scope->data.ptrvalue, tbsp))
        {
            bsp = tbsp;
        }
        else
        {                  /* not in scope, could be duplicate SeqId */
            i = num-1;
            while ((i >= 0) && (bsp == NULL) 
                   && (StringLen (sipp[i]->str) > len
                       && sipp[i]->str[len] == '|'
                       && StringNCmp(sipp[i]->str, rev_id, len) == 0))  /* back up */
            {
               tbsp = (BioseqPtr)(sipp[i]->omdp->dataptr);
               if (ObjMgrIsChild(scope->data.ptrvalue, tbsp))
               {
                   bsp = tbsp;
               }
               i--;
            }
            i = num + 1;
            imax = smp->BioseqIndexCnt - 1;
            while ((bsp == NULL) && (i <= imax)
                   && (StringLen (sipp[i]->str) > len
                       && sipp[i]->str[len] == '|'
                       && StringNCmp(sipp[i]->str, rev_id, len) == 0)) 
            {
               tbsp = (BioseqPtr)(sipp[i]->omdp->dataptr);
               if (ObjMgrIsChild(scope->data.ptrvalue, tbsp))
               {
                   bsp = tbsp;
               }
               i++;
            }
        }
    }
  }

  SeqMgrUnlock();
  return bsp;
}


NLM_EXTERN ValNodePtr GetTranscriptomeIdsList (FILE *fp, SeqEntryPtr sep, ValNodePtr PNTR err_list)
{
  ReadBufferData    rbd;
  ValNodePtr        list = NULL, token_list;
  SeqIdPtr          sip;
  BioseqPtr         bsp;
  CharPtr           line, err_str;
  CharPtr           bad_id_fmt = "Unable to make SeqId from %s";

  rbd.current_data = NULL;
  rbd.fp = fp;

  line = AbstractReadFunction (&rbd); 

  while (line != NULL && line[0] != EOF) {
    if (!StringHasNoText (line)) {
      token_list = MakeTranscriptomeIDTokensFromLine (line);
      if (token_list != NULL && token_list->next != NULL) {     
        sip = CreateSeqIdFromText (token_list->data.ptrvalue, sep);
        bsp = BioseqFind (sip);
        sip = SeqIdFree (sip);
        if (bsp == NULL) {
          bsp = GetTranscriptomeBioseqFromStringId (token_list->data.ptrvalue);
        }
        if (bsp == NULL) {
          err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (token_list->data.ptrvalue) + StringLen (bad_id_fmt)));
          sprintf (err_str, bad_id_fmt, token_list->data.ptrvalue);
          ValNodeAddPointer (err_list, 0, err_str);
        } else {
          ValNodeAddPointer (&list, 0, TranscriptomeIdsNew(bsp, token_list->next));
          token_list->next = NULL;
        }
      }
      token_list = ValNodeFreeData (token_list);
    }
    line = MemFree (line);
    line = AbstractReadFunction (&rbd);
  }
  return list;      
}


static void GetExistingTSATableIdsCallback (BioseqPtr bsp, Pointer userdata)
{
  TranscriptomeIdsPtr t;
  ValNodePtr        token_list = NULL;
  SeqAlignPtr       salp;
  Char              buf[255];
  DenseSegPtr       dsp;

  if (bsp == NULL || userdata == NULL || ISA_aa (bsp->mol) || bsp->hist == NULL || bsp->hist->assembly == NULL) {
    return;
  }

  salp = bsp->hist->assembly;
  while (salp != NULL) {
    if (salp->segtype == SAS_DENSEG && salp->segs != NULL) {
      dsp = (DenseSegPtr) salp->segs;
      if (dsp->dim == 2 && dsp->ids != NULL && dsp->ids->next != NULL) {
        SeqIdWrite (dsp->ids->next, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
        ValNodeAddPointer (&token_list, 0, StringSave (buf));
      }
    }
    salp = salp->next;
  }

  if (token_list != NULL) {
    t = TranscriptomeIdsNew (bsp, token_list);
    if (t != NULL) {
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, t);
    }
  }
}


NLM_EXTERN ValNodePtr GetExistingTSATableIds (SeqEntryPtr sep)
{
  ValNodePtr ids_list = NULL;

  VisitBioseqsInSep (sep, &ids_list, GetExistingTSATableIdsCallback);
  return ids_list;
}


static Int4 ReadNumberFromPortionOfString (CharPtr str, Int4 len)
{
  CharPtr num_buf;
  Int4    val;

  if (str == NULL) return -1;
  num_buf = (CharPtr) MemNew (sizeof (Char) * (len + 1));
  StringNCpy (num_buf, str, len);
  num_buf[len] = 0;
  val = atoi (num_buf);
  num_buf = MemFree (num_buf);
  return val;
}


static Boolean IsStringInSpan (CharPtr str, CharPtr first, CharPtr second)
{
  Int4 prefix_len = 0;
  Boolean rval = FALSE;
  Int4 first_num, second_num, str_num;
  CharPtr cp, cp1, cp2, suf1, suf2, suf_str;

  if (StringHasNoText (str)) {
    return FALSE;
  } else if (StringCmp (str, first) == 0 || StringCmp (str, second) == 0) {
    return TRUE;
  } else if (StringHasNoText (first) || StringHasNoText (second)) {
    return FALSE;
  }

  if (IsAllDigits (first)) {
    if (IsAllDigits (str) && IsAllDigits (second)) {
      str_num = atoi (str);
      first_num = atoi (first);
      second_num = atoi (second);
      if ((str_num > first_num && str_num < second_num)
          || (str_num > second_num && str_num < first_num)) {
        rval = TRUE;
      }
    }
  } else if (IsAllDigits(second)) {
    cp = first;
    while (!isdigit (*cp)) {
      prefix_len ++;
      cp++;
    }
    if (StringNCmp (str, first, prefix_len) == 0
        && IsAllDigits (str + prefix_len)
        && IsAllDigits (first + prefix_len)) {
      first_num = atoi (cp);
      second_num = atoi (second);
      str_num = atoi (str);
      if ((str_num > first_num && str_num < second_num)
          || (str_num > second_num && str_num < first_num)) {
        rval = TRUE;
      }
    }
  } else {
    /* determine length of prefix */
    cp1 = first;
    cp2 = second;
    while (*cp1 != 0 && *cp2 != 0 && *cp1 == *cp2) {
      prefix_len++;
      cp1++;
      cp2++;
    }
    if (*cp1 != 0 && *cp2 != 0 
        && isdigit (*cp1) && isdigit (*cp2)
        && StringNCmp (str, first, prefix_len) == 0) {
      if (IsAllDigits (cp1) && IsAllDigits (cp2) && IsAllDigits (str + prefix_len)) {
        first_num = atoi (cp1);
        second_num = atoi (cp2);
        str_num = atoi (str + prefix_len);
        if ((str_num > first_num && str_num < second_num)
            || (str_num > second_num && str_num < first_num)) {
          rval = TRUE;
        }
      } else {
        /* determine whether there is a suffix */
        suf1 = cp1 + StringSpn (cp1, "0123456789");
        suf2 = cp2 + StringSpn (cp2, "0123456789");
        suf_str = str + prefix_len + StringSpn (str + prefix_len, "0123456789");
        if (StringCmp (suf1, suf2) == 0 && StringCmp (suf1, suf_str) == 0) {
          /* suffixes match */
          first_num = ReadNumberFromPortionOfString (cp1, suf1 - cp1);
          second_num = ReadNumberFromPortionOfString (cp2, suf2 - cp2);
          str_num = ReadNumberFromPortionOfString (str + prefix_len, suf_str - str - prefix_len);
          if ((str_num > first_num && str_num < second_num)
              || (str_num > second_num && str_num < first_num)) {
            rval = TRUE;
          }
        }
      }
    }
  }
  return rval;
}


static Boolean GetSpanFromHyphenInString (CharPtr str, CharPtr hyphen, CharPtr PNTR first, CharPtr PNTR second)
{
  CharPtr cp;
  Int4    len;

  *first = NULL;
  *second = NULL;

  if (hyphen == str) {
    return FALSE;
  }

  /* find range start */
  cp = hyphen - 1;
  while (isspace (*cp) && cp != str) {
    cp--;
  }
  
  while (!isspace (*cp) && *cp != ',' && *cp != ';' && cp != str) {
    cp--;
  }
  
  len = hyphen - cp;
  *first = (CharPtr) MemNew (sizeof (Char) * (len + 1));
  StringNCpy (*first, cp, len);
  (*first)[len] = 0;
  TrimSpacesAroundString (*first);
  
  /* find range end */
  cp = hyphen + 1;
  while (isspace (*cp)) {
    cp++;
  }
  while (*cp != 0 && !isspace (*cp) && *cp != ',' && *cp != ';') {
    cp++;
  }

  len = cp - hyphen;
  if (*cp != 0 && !isspace (*cp)) {
    len--;
  }
  *second = (CharPtr) MemNew (sizeof (Char) * (len + 1));
  StringNCpy (*second, hyphen + 1, len);
  (*second)[len] = 0;
  TrimSpacesAroundString (*second);  
  return TRUE;
}


NLM_EXTERN Boolean IsStringInSpanInList (CharPtr str, CharPtr list)
{
  CharPtr cp, hyphen;
  Int4    prefix_len = 0, suffix_len;
  CharPtr num_start, range_start = NULL, range_end = NULL;
  Int4    str_val;
  Boolean rval = FALSE;

  if (StringHasNoText (list) || StringHasNoText (str)) {
    return FALSE;
  }

  cp = str;
  while (isalpha (*cp)) {
    prefix_len++;
    cp++;
  }
  if (*cp == 0) {
    return FALSE;
  }

  num_start = cp;
  while (isdigit (*cp)) {
    cp++;
  }
  suffix_len = StringLen (cp);

  str_val = ReadNumberFromPortionOfString (num_start, cp - num_start);

  /* find ranges */

  hyphen = StringChr (list, '-');
  while (hyphen != NULL && !rval) {
    if (hyphen == list) {
      hyphen = StringChr (hyphen + 1, '-');
    } else {
      if (GetSpanFromHyphenInString (list, hyphen, &range_start, &range_end)) {
        if (IsStringInSpan (str, range_start, range_end)) {
          rval = TRUE;
        }
        range_start = MemFree (range_start);
        range_end = MemFree (range_end);
      }
      hyphen = StringChr (hyphen + 1, '-');
    }
  }
  
  return rval;
}


static void ParseGoTermsFromFieldsSfp (UserObjectPtr uop, Pointer userdata)
{
  ObjectIdPtr  oip;
  UserFieldPtr ufp_process, ufp_new, ufp2, ufp, ufp_last = NULL, new_list = NULL, last_new = NULL;
  CharPtr      cp;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    for (ufp_process = uop->data; ufp_process != NULL; ufp_process = ufp_process->next) {
      if (ufp_process->label == NULL 
          || (StringCmp (ufp_process->label->str, "Process") != 0
              && StringCmp (ufp_process->label->str, "Function") != 0
              && StringCmp (ufp_process->label->str, "Component") != 0)) {
        continue;
      }
      for (ufp2 = ufp_process->data.ptrvalue; ufp2 != NULL; ufp2 = ufp2->next) {
        if (ufp2->choice != 11) continue;
        new_list = NULL;
        ufp_last = NULL;
        last_new = NULL;
        for (ufp = ufp2->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
          ufp_last = ufp;
          if (ufp->label != NULL && StringCmp (ufp->label->str, "text string") == 0
              && (cp = StringISearch (ufp->data.ptrvalue, "GO:")) != NULL) {
            ufp_new = UserFieldNew ();
            ufp_new->label = ObjectIdNew ();
            ufp_new->label->str = StringSave ("go id");
            ufp_new->choice = 1; /* visible string - need to keep leading zeroes */
            ufp_new->data.ptrvalue = StringSave (cp + 3);
            *cp = 0;
            cp--;
            while (cp > (CharPtr) ufp->data.ptrvalue && isspace (*cp)) {
              *cp = 0;
              cp--;
            }
            if (last_new == NULL) {
              new_list = ufp_new;
            } else {
              last_new->next = ufp_new;
            }
            last_new = ufp_new;
          }
        }
        if (new_list != NULL) {
          if (ufp_last == NULL) {
            uop->data = new_list;
          } else {
            ufp_last->next = new_list;
          }
        }
      }
    }
  }
}


static void ParseGoTermsFromFieldsCallback (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->ext == NULL) return;

  VisitUserObjectsInUop (sfp->ext, NULL, ParseGoTermsFromFieldsSfp);

}


NLM_EXTERN void ParseGoTermsFromFields (SeqEntryPtr sep)
{
  VisitFeaturesInSep (sep, NULL, ParseGoTermsFromFieldsCallback);
}


typedef  Int4  (*Nlm_LenToToken) PROTO ((CharPtr));
typedef  Int4  (*Nlm_TokenLen) PROTO ((CharPtr));
typedef  Boolean (*Nlm_IsToken) PROTO ((CharPtr));

static Int4 LenToTabToken (CharPtr cp) 
{
  if (cp == NULL) {
    return 0;
  } else {
    return StringCSpn (cp, "\t\n");
  }
}


static Int4 LenTabToken (CharPtr cp)
{
  if (cp == NULL) {
    return 0;
  } else {
    return 1;
  }
}


static Boolean IsTab (CharPtr cp)
{
  if (cp == NULL) {
    return FALSE;
  } else if (*cp == '\t') {
    return TRUE;
  } else {
    return FALSE;
  }
}


static ValNodePtr 
ReadOneColumnListEx 
(CharPtr line, Nlm_LenToToken len_to_token_func, Nlm_TokenLen token_len_func, Nlm_IsToken is_token_func)
{
  CharPtr p_start, p_end, p_quote, zap;
  Char    tmp_end;
  Int4    plen, quote_len;
  Boolean found_end;
  ValNodePtr col_list = NULL;
  Char       term[2];
  ValNodeBlock col_block;

  if (StringHasNoText (line) || len_to_token_func == NULL || token_len_func == NULL || is_token_func == NULL) return NULL;
  term[1] = 0;
  p_start = line;
  found_end = FALSE;
  col_block.head = NULL;
  col_block.tail = NULL;
  while (*p_start != 0 && !found_end)
  {
    quote_len = StringCSpn (p_start, "\"");
    plen = len_to_token_func (p_start);
    if (quote_len < plen) {
      p_quote = p_start + quote_len;
      term[0] = *p_quote;
      plen = quote_len + StringCSpn (p_start + quote_len + 1, term);
      plen += len_to_token_func (p_start + plen);
    }
    if (plen == 0)
    {
      if (is_token_func(p_start))
      {
        ValNodeAddPointerToEnd (&col_block, 0, StringSave (""));
      }
      p_start+=token_len_func(p_start);
      if (*p_start == 0) {
        if (col_list != NULL)
        {
          ValNodeAddPointerToEnd (&col_block, 0, StringSave (""));
        }
      }
      continue;
    }
    if (plen == StringLen (p_start))
    {
      found_end = TRUE;
      p_end = p_start + plen;
    }
    else
    {
      p_end = p_start + plen;
      tmp_end = *p_end;
      *p_end = 0;
    }
    while (*p_start == ' ') {
      ++p_start;
    }
    zap = p_end - 1;
    while (zap > p_start && *zap == ' ') {
      *zap = 0;
      zap--;
    }
    
    ValNodeAddPointerToEnd (&col_block, 0, StringSave (p_start));
    if (!found_end)
    {
      *p_end = tmp_end;
      p_start = p_end + token_len_func(p_end);
    }
  }
  return col_block.head; 
}


NLM_EXTERN ValNodePtr ReadOneColumnList (CharPtr line)
{
  return ReadOneColumnListEx (line, LenToTabToken, LenTabToken, IsTab);
}


NLM_EXTERN ValNodePtr ExtractNthValNode (ValNodePtr PNTR list, Int4 nth)
{
  ValNodePtr prev = NULL, this_vnp;
  if (nth < 0)
  {
    return NULL;
  }
  
  this_vnp = *list;
  while (nth > 0 && this_vnp != NULL)
  {
    prev = this_vnp;
    this_vnp = this_vnp->next;
    nth --;
  }
  
  if (this_vnp != NULL)
  {
    if (prev == NULL)
    {
      *list = (*list)->next;
    }
    else
    {
      prev->next = this_vnp->next;
    }
    this_vnp->next = NULL;
  }
  return this_vnp;
  
}


static void RemoveEmptyRowsFromTabTable (ValNodePtr PNTR line_list)
{
  ValNodePtr vnp_prev = NULL, vnp_next, vnp;

  if (line_list == NULL || *line_list == NULL) {
    return;
  }
  for (vnp = *line_list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (vnp->data.ptrvalue == NULL) {
      if (vnp_prev == NULL) {
        *line_list = vnp_next;
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


static void RemoveEmptyColumnsFromTabTable (ValNodePtr PNTR line_list)
{
  ValNodePtr row_vnp, col_vnp, del_vnp;
  Int4       num_col, max_col = 0, i;
  BoolPtr    col_empty;

  if (line_list == NULL || *line_list == NULL) {
    return;
  }

  for (row_vnp = *line_list; row_vnp != NULL; row_vnp = row_vnp->next) {
    num_col = ValNodeLen (row_vnp->data.ptrvalue);
    if (num_col > max_col) {
      max_col = num_col;
    }
  }

  col_empty = (BoolPtr) MemNew (sizeof (Boolean) * max_col);
  for (i = 0; i < max_col; i++) {
    col_empty[i] = TRUE;
  }

  for (row_vnp = *line_list; row_vnp != NULL; row_vnp = row_vnp->next) {
    for (col_vnp = row_vnp->data.ptrvalue, num_col = 0;
         col_vnp != NULL;
         col_vnp = col_vnp->next, num_col++) {
      if (!StringHasNoText (col_vnp->data.ptrvalue)) {
        col_empty[num_col] = FALSE;
      }
    }
  }

  for (i = max_col - 1; i >= 0; i--) {
    if (col_empty[i]) {
      for (row_vnp = *line_list; row_vnp != NULL; row_vnp = row_vnp->next) {  
        col_vnp = row_vnp->data.ptrvalue;
        del_vnp = ExtractNthValNode (&col_vnp, i);
        row_vnp->data.ptrvalue = col_vnp;
        del_vnp = ValNodeFreeData (del_vnp);
      }
    }
  }
  col_empty = MemFree (col_empty);

  RemoveEmptyRowsFromTabTable (line_list);
}


NLM_EXTERN ValNodePtr ReadTabTableFromFile (FILE *fp)
{
  Int4          max_columns, num_cols, num_discarded = 0;
  ValNodePtr    header_line;
  ValNodePtr    line_list, column_list;
  ValNodePtr    vnp, last_line = NULL;
  ReadBufferData rbd;
  CharPtr        line;

  if (fp == NULL) return NULL;
  rbd.fp = fp;
  rbd.current_data = NULL;

  line_list = NULL;
  max_columns = 0;
  header_line = NULL;
  line = AbstractReadFunction (&rbd);
  while (line != NULL) 
  {
    column_list = ReadOneColumnList (line);
    if (column_list != NULL)
    {
      vnp = ValNodeAddPointer (&last_line, 0, column_list);
      if (line_list == NULL) {
        line_list = last_line;
      }
      last_line = vnp;
      num_cols = ValNodeLen (column_list);
      if (num_cols > max_columns)
      {
        max_columns = num_cols;
        header_line = vnp;
      }
    }
    line = MemFree (line);
    line = AbstractReadFunction (&rbd);
  }
  /* throw out all lines before header line */
  if (header_line != line_list)
  {
    num_discarded = 1;
    vnp = line_list;
    while (vnp != NULL && vnp->next != header_line)
    {
      num_discarded++;
      vnp = vnp->next;
    }
    if (vnp != NULL) {
      vnp->next = NULL;
    }
    ValNodeFreeData (line_list);
    line_list = NULL;
    Message (MSG_OKC, "Warning - the first row of the table did not have enough columns for headers for the following rows - %d rows were discarded before a row with enough columns to provide headers was found.", num_discarded);
  }

  RemoveEmptyColumnsFromTabTable (&header_line);
  return header_line;
}


NLM_EXTERN ValNodePtr FlipTabTableAxes (ValNodePtr row_list)
{
  ValNodePtr vnp, vnp_c;
  ValNodePtr new_table = NULL, vnp_new_row = NULL, vnp_new;

  if (row_list == NULL) {
    return NULL;
  }

  new_table = ValNodeNew (NULL);
  for (vnp = row_list; vnp != NULL; vnp = vnp->next) {
    vnp_c = vnp->data.ptrvalue;
    vnp_new_row = new_table;
    while (vnp_c != NULL) {
      if (vnp_new_row == NULL) {
        vnp_new_row = ValNodeNew (new_table);
      }
      vnp_new = vnp_new_row->data.ptrvalue;
      ValNodeAddPointer (&vnp_new, 0, StringSave (vnp_c->data.ptrvalue));
      vnp_new_row->data.ptrvalue = vnp_new;
      vnp_new_row = vnp_new_row->next;
      vnp_c = vnp_c->next;
    }
  }

  return new_table;
}

NLM_EXTERN ValNodePtr FreeTabTable (ValNodePtr row_list)
{
  ValNodePtr row_vnp, column_list;
  
  if (row_list != NULL)
  {
    /* free table text */
    for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next)
    {
      column_list = (ValNodePtr) row_vnp->data.ptrvalue;
      row_vnp->data.ptrvalue = ValNodeFreeData (column_list);
    }
    row_list = ValNodeFree (row_list);
  }
  return row_list;
}


NLM_EXTERN ValNodePtr CopyTabTable (ValNodePtr row_list)
{
  ValNodeBlock row_block;
  ValNodeBlock col_block;
  ValNodePtr row, col;

  InitValNodeBlock(&row_block, NULL);
  for (row = row_list; row != NULL; row = row->next) {
    InitValNodeBlock(&col_block, NULL);
    for (col = row->data.ptrvalue; col != NULL; col = col->next) {
      ValNodeAddPointerToEnd (&col_block, col->choice, StringSave(col->data.ptrvalue));
    }
    ValNodeAddPointerToEnd (&row_block, 0, col_block.head);
  }
  return row_block.head;
}


NLM_EXTERN void WriteTabTableToFile (ValNodePtr table, FILE *fp)
{
  ValNodePtr line, vnp;

  for (line = table; line != NULL; line = line->next) {
    for (vnp = line->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      fprintf (fp, "%s%s", vnp->data.ptrvalue == NULL ? "" : (CharPtr) vnp->data.ptrvalue, vnp->next == NULL ? "\n" : "\t");
    }
  }
}


NLM_EXTERN ValNodePtr CountTabTableBlanks (ValNodePtr row_list)
{
  ValNodePtr   line_vnp, col_vnp, blank_vnp;
  Int4         num_rows = 0;
  ValNodePtr   blank_list = NULL;
  
  if (row_list == NULL) return NULL;  

  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next)
  {
    col_vnp = line_vnp->data.ptrvalue;
    blank_vnp = blank_list;
    while (col_vnp != NULL || blank_vnp != NULL) {
      if (blank_vnp == NULL) {
        /* for all rows prior to this one, this column was blank */
        blank_vnp = ValNodeAddInt (&blank_list, 0, num_rows);
      }
      if (col_vnp == NULL || StringHasNoText (col_vnp->data.ptrvalue)) {
        blank_vnp->data.intvalue ++;
      }
      if (col_vnp != NULL) {
        col_vnp = col_vnp->next;
      }
      blank_vnp = blank_vnp->next;
    }
    num_rows ++;
  }
  return blank_list;
}


NLM_EXTERN void RemoveQuotesFromTabTable (ValNodePtr row_list)
{
  ValNodePtr line_vnp, col_vnp;
  CharPtr    val;
  Int4       len, i;

  if (row_list == NULL) return;

  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next)
  {
    col_vnp = line_vnp->data.ptrvalue;
    while (col_vnp != NULL) {
      val = col_vnp->data.ptrvalue;
      len = StringLen (val);
      /* remove double quotes */
      if (val != NULL && val[0] == '"' && val[len - 1] == '"') {
        for (i = 1; i < len - 1; i++) {
          val[i - 1] = val[i];
        }
        val[i - 1] = 0;
      }
      col_vnp = col_vnp->next;
    }
  }
}


NLM_EXTERN void ReparseTabTableConvertFirstSpaceToTab (ValNodePtr row_list)
{
  ValNodePtr line_vnp, col_vnp, new_vnp;
  CharPtr    first_text, second_text, first_space;

  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next) 
  {
    col_vnp = line_vnp->data.ptrvalue;
    first_text = col_vnp->data.ptrvalue;
    if ((first_space = StringChr (first_text, ' ')) != NULL) {
      second_text = first_space + StringSpn (first_space, " ");
      if (*second_text != 0) {
        /* terminate first text at first space */
        *first_space = 0;
        /* create new column with text after first space */
        second_text = StringSave (second_text);
        new_vnp = ValNodeNew (NULL);
        new_vnp->data.ptrvalue = second_text;
        /* insert new column */
        new_vnp->next = col_vnp->next;
        col_vnp->next = new_vnp;
      }
    }
  }
}


NLM_EXTERN void ReparseTabTableSeparateColumnAtDelimiter (ValNodePtr row_list, Char delimiter, Int4 col, Boolean stop_after_first)
{
  ValNodePtr line_vnp, col_vnp, new_vnp, next_col;
  CharPtr    first_text, second_text, first_space;
  Int4       col_num;

  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next) 
  {
    col_vnp = line_vnp->data.ptrvalue;
    col_num = 0;
    while (col_num < col && col_vnp != NULL) {
      col_num++;
      col_vnp = col_vnp->next;
    }
    if (col_vnp != NULL) {
      next_col = col_vnp->next;
      while (col_vnp != next_col) {
        first_text = col_vnp->data.ptrvalue;
        if ((first_space = StringChr (first_text, delimiter)) != NULL) {
          second_text = first_space + 1;
          if (*second_text != 0) {
            /* terminate first text at first delimiter */
            *first_space = 0;
            /* create new column with text after first delimiter */
            second_text = StringSave (second_text);
            new_vnp = ValNodeNew (NULL);
            new_vnp->data.ptrvalue = second_text;
            /* insert new column */
            new_vnp->next = col_vnp->next;
            col_vnp->next = new_vnp;
          }
        }
        if (stop_after_first) {
          col_vnp = next_col;
        } else {
          col_vnp = col_vnp->next;
        }
      }
    }
  }
}


static Int4 LenToNextTabOrMultispace (CharPtr cp)
{
  Int4 len = 0;
  Boolean found = FALSE;

  if (StringHasNoText (cp)) {
    return 0;
  }

  while (!found && *cp != 0) {
    if (*cp == '\t' || *cp == '\n' || (*cp == ' ' && *(cp + 1) == ' ')) {
      found = TRUE;
    } else {
      ++len;
      ++cp;
    }
  }
  return len;
}


static Int4 LenTabOrMultispace (CharPtr cp)
{
  Int4 len = 0;

  if (StringHasNoText (cp)) {
    len = 0;
  } else if (*cp == '\t' || *cp == '\n') {
    len = 1;
  } else {
    len = StringSpn (cp, " ");
  }
  return len;
}


static Boolean IsTabOrMultiSpace (CharPtr cp)
{
  if (cp == NULL || *cp == 0) {
    return FALSE;
  } else if (*cp == '\t' || *cp == '\n') {
    return TRUE;
  } else if (*cp == ' ' && *(cp + 1) == ' ') {
    return TRUE;
  } else {
    return FALSE;
  }
}


NLM_EXTERN void ReparseTabTableConvertMultiSpaceToTab (ValNodePtr row_list)
{
  ValNodePtr line_vnp, col_vnp, new_cols, col_prev, col_next, last_vnp;

  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next) 
  {
    col_prev = NULL;
    for (col_vnp = line_vnp->data.ptrvalue; col_vnp != NULL; col_vnp = col_next) 
    {
      col_next = col_vnp->next;
      new_cols = ReadOneColumnListEx (col_vnp->data.ptrvalue, LenToNextTabOrMultispace, LenTabOrMultispace, IsTabOrMultiSpace);
      if (new_cols != NULL) {
        /* insert new columns */
        last_vnp = new_cols;
        while (last_vnp->next != NULL) {
          last_vnp = last_vnp->next;
        }
        last_vnp->next = col_vnp->next;
        col_vnp->next = NULL;
        col_vnp = ValNodeFreeData (col_vnp);
        if (col_prev == NULL) {
          line_vnp->data.ptrvalue = new_cols;
        } else {
          col_prev->next = new_cols;
        }
        col_prev = last_vnp;
      } else {
        col_prev = col_vnp;
      }
    }
  }
}


/* first intended use is to create file ID columns from first two columns of table.
 * second intended use is to combine columns to make ID list.
 * Note that column_pos needs to be a sorted list of integers giving a zero-based column offset.
 */
NLM_EXTERN void CombineTabTableColumns (ValNodePtr row_list, ValNodePtr column_pos, CharPtr delimiter)
{
  ValNodePtr line_vnp, col_vnp, col_prev, col_next, offset_vnp, add_vnp;
  Int4       col_num, len;
  CharPtr    tmp;

  if (row_list == NULL || column_pos == NULL || column_pos->next == NULL) {
    return;
  }


  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next) 
  {
    col_prev = NULL;
    add_vnp = NULL;
    offset_vnp = column_pos;
    for (col_vnp = line_vnp->data.ptrvalue, col_num = 0; col_vnp != NULL && offset_vnp != NULL; col_vnp = col_next, col_num++) 
    {
      col_next = col_vnp->next;
      if (col_num == offset_vnp->data.intvalue) {
        if (add_vnp == NULL) {
          add_vnp = col_vnp;
          col_prev = col_vnp;
        } else {
          if (StringHasNoText (col_vnp->data.ptrvalue)) {
            /* do nothing - no need to add blank to blank */
          } else if (StringHasNoText (add_vnp->data.ptrvalue)) {
            /* move from col_vnp */
            add_vnp->data.ptrvalue = MemFree (add_vnp->data.ptrvalue);
            add_vnp->data.ptrvalue = col_vnp->data.ptrvalue;
            col_vnp->data.ptrvalue = NULL;
          } else {
            /* combine with delimiter */
            len = StringLen (add_vnp->data.ptrvalue) + StringLen (delimiter) + StringLen (col_vnp->data.ptrvalue) + 1;
            tmp = (CharPtr) MemNew (sizeof (Char) * len);
            sprintf (tmp, "%s%s%s", (char *) add_vnp->data.ptrvalue, delimiter == NULL ? "" : delimiter, (char *) col_vnp->data.ptrvalue);
            add_vnp->data.ptrvalue = MemFree (add_vnp->data.ptrvalue);
            add_vnp->data.ptrvalue = tmp;
          }
          col_prev->next = col_vnp->next;
          col_vnp->next = NULL;
          col_vnp = ValNodeFreeData (col_vnp);
        }
        offset_vnp = offset_vnp->next;
      } else {
        col_prev = col_vnp;
      }
    }
  }
}


NLM_EXTERN void AddTextToTabTableColumn (ValNodePtr row_list, Int4 col, CharPtr text, Uint2 existing_text)
{
  ValNodePtr row_vnp, col_vnp;
  Int4       i;
  CharPtr    str;

  if (text == NULL) {
    return;
  }

  for (row_vnp = row_list; row_vnp != NULL; row_vnp = row_vnp->next) {
    for (col_vnp = row_vnp->data.ptrvalue, i = 0;
         col_vnp != NULL;
         col_vnp = col_vnp->next, i++) {
      if (i == col) {
        str = col_vnp->data.ptrvalue;
        SetStringValue (&str, text, existing_text);
        col_vnp->data.ptrvalue = str;
        break;
      }
    }
  }
}


static int LIBCALLBACK SortTableRowByColumn (VoidPtr ptr1, VoidPtr ptr2, Int4 column)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ValNodePtr  col1, col2;
  Int4        colpos = 1;
  int         rval = 0;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  col1 = vnp1->data.ptrvalue;
  col2 = vnp2->data.ptrvalue;
  while (col1 != NULL && col2 != NULL && colpos < column) {
    col1 = col1->next;
    col2 = col2->next;
    colpos++;
  }
  if (col1 == NULL && col2 == NULL) {
    rval = 0;
  } else if (col1 == NULL) {
    rval = -1;
  } else if (col2 == NULL) {
    rval = 1;
  } else {
    rval = StringCmp (col1->data.ptrvalue, col2->data.ptrvalue);
  }
  return rval;
}


static Int4 s_TableRowSortColumn = 0;

static int LIBCALLBACK SortTableRowByColumnStatic (VoidPtr ptr1, VoidPtr ptr2)
{
  return SortTableRowByColumn (ptr1, ptr2, s_TableRowSortColumn);
}


NLM_EXTERN ValNodePtr SortTableRowByAnyColumn (ValNodePtr table, Int4 column)
{
  s_TableRowSortColumn = column;
  table = ValNodeSort (table, SortTableRowByColumnStatic);
  return table;
}


NLM_EXTERN TwoStringHashPtr TwoStringHashFree (TwoStringHashPtr tsh)
{
  Int4 i;

  if (tsh != NULL) {
    for (i = 0; i < tsh->num_lines; i++) {
      tsh->table[2 * i] = MemFree (tsh->table[2 * i]);
      tsh->table[2 * i + 1] = MemFree (tsh->table[2 * i + 1]);
    }
    tsh->table = MemFree (tsh->table);
    tsh = MemFree (tsh);
  }
  return tsh;
}


NLM_EXTERN ValNodePtr GetNthValNode (ValNodePtr list, Int4 n)
{
  Int4 pos = 1;
  ValNodePtr vnp;

  if (n < 1) {
    return NULL;
  }
  for (vnp = list; vnp != NULL && pos < n; vnp = vnp->next) 
  {
    pos++;
  }
  return vnp;
}


NLM_EXTERN TwoStringHashPtr MakeTwoStringHashFromTabTable (ValNodePtr line_list, Int4 column1, Int4 column2)
{
  ValNodePtr tmp, vnp, col1, col2;
  Int4       len, i;
  TwoStringHashPtr tsh;

  tmp = CopyTabTable(line_list);
  tmp = SortTableRowByAnyColumn (tmp, column1);
  len = ValNodeLen (tmp);

  tsh = (TwoStringHashPtr) MemNew (sizeof (TwoStringHashData));
  tsh->table = (CharPtr PNTR) MemNew (sizeof (CharPtr) * len * 2);
  for (i = 0, vnp = tmp; vnp != NULL; vnp = vnp->next) {
    col1 = GetNthValNode (vnp->data.ptrvalue, column1);
    col2 = GetNthValNode (vnp->data.ptrvalue, column2);
    if (col1 != NULL && col2 != NULL && !StringHasNoText (col1->data.ptrvalue) && !StringHasNoText (col2->data.ptrvalue)) {
      tsh->table[2 * i] = StringSave (col1->data.ptrvalue);
      tsh->table[2 * i + 1] = StringSave (col2->data.ptrvalue);
      i++;
    }
  }
  tsh->num_lines = i;
  tmp = FreeTabTable(tmp);
  return tsh;
}


NLM_EXTERN CharPtr GetValueFromTwoStringHash (CharPtr key, TwoStringHashPtr tsh)
{
  Int4    min = 0, num = -1, i, j;
  Int4    max;
  CharPtr tmp;

  if (StringHasNoText (key) || tsh == NULL) {
    return NULL;
  }
  max = tsh->num_lines - 1;

  while (max >= min)
  {
    i = (max + min)/2;
    tmp = tsh->table[2 * i];
    if ((j = StringCmp(tmp, key)) > 0)
    {
      max = i - 1;
    }
    else if (j < 0)
    {
      min = i + 1;
    }
    else
    {
      num = i;
      break;
    }
  }
  if (num == -1) {
    return NULL;
  } else {
    return tsh->table[2 * num + 1];
  }
}


static void AddToContextList (Char ch, CharPtr PNTR strp, ValNodePtr PNTR search_list)
{
  ValNodePtr vnp, vnp_last = NULL, vnp2, clist;

  if (strp == NULL || search_list == NULL) return;

  /* group contexts for the same character together */
  vnp = *search_list;
  while (vnp != NULL && vnp->choice != (Uint1)ch)
  {
    vnp_last = vnp;
    vnp = vnp->next;
  }
  if (vnp == NULL)
  {
    vnp = ValNodeNew(NULL);
    if (vnp_last == NULL) 
    {
      *search_list = vnp;
    }
    else
    {
      vnp_last->next = vnp;
    }
  }
  vnp->choice = (Uint1) ch;
  clist = vnp->data.ptrvalue;
  /* don't add the same string twice for the same character */
  vnp2 = clist;
  vnp_last = NULL;
  while (vnp2 != NULL && vnp2->data.ptrvalue != strp)
  {
    vnp_last = vnp2;
    vnp2 = vnp2->next;
  }
  if (vnp2 == NULL)
  {
    vnp2 = ValNodeNew (NULL);
    vnp2->data.ptrvalue = strp;
    if (vnp_last == NULL)
    {
      clist = vnp2;
    }
    else
    {
      vnp_last->next = vnp2;
    }
  }
  vnp->data.ptrvalue = clist;
}


NLM_EXTERN ValNodePtr FreeContextList (ValNodePtr context_list)
{
  ValNodePtr vnp;

  for (vnp = context_list; vnp != NULL; vnp = vnp->next)
  {
    vnp->data.ptrvalue = ValNodeFree (vnp->data.ptrvalue);
  }
  context_list = ValNodeFree (context_list);
  return context_list;
}


NLM_EXTERN void SpecialCharFindWithContext (CharPtr PNTR strp, Pointer userdata, BoolPtr did_find, BoolPtr did_change)
{
  CharPtr cp;
  Boolean found_any = FALSE;

  if (strp == NULL || *strp == NULL || userdata == NULL) return;

  cp = *strp;
  while (*cp != 0)
  {
    if (*cp < ' ' || *cp > '~')
    {
      found_any = TRUE;
      AddToContextList (*cp, strp, (ValNodePtr PNTR) userdata);
    }
    cp++;
  }
  if (found_any && did_find != NULL)
  {
    *did_find = TRUE;
  }
}


NLM_EXTERN ValNodePtr ScanTabTableForSpecialCharacters (ValNodePtr row_list)
{
  ValNodePtr special_list = NULL, line_vnp, col_vnp;

  if (row_list == NULL) return NULL;  

  for (line_vnp = row_list; line_vnp != NULL; line_vnp = line_vnp->next)
  {
    col_vnp = line_vnp->data.ptrvalue;
    while (col_vnp != NULL) {
      SpecialCharFindWithContext ((CharPtr PNTR)&(col_vnp->data.ptrvalue), &special_list, NULL, NULL);
      col_vnp = col_vnp->next;
    }
  }
  return special_list;
}


NLM_EXTERN ValNodePtr AutoReplaceSpecialCharactersInText (CharPtr PNTR text)
{
  CharPtr cp, str, new_str, cp_dst;
  Int4    len;
  Int4    extra_len = 0;
  Boolean any = FALSE;
  CharPtr replace_fmt = "Replaced '%c' with '%s'";
  ValNodePtr repl_list = NULL;
  CharPtr repl_str;

  if (text == NULL || (cp = *text) == NULL) {
    return NULL;
  }

  while (*cp != 0) {
    if (*cp < ' ' || *cp > '~') {
#ifdef OS_WINNT
      str = GetSpecialWinCharacterReplacement ((unsigned char) *cp);
#else
      str = GetSpecialMacCharacterReplacement ((unsigned char) *cp);
#endif
      len = StringLen (str);
      if (len > 1) {
        extra_len += len - 1;
      }
      any = TRUE;
    }
    ++cp;
  }
  if (any) {
    new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (*text) + extra_len + 1));
    cp = *text;
    cp_dst = new_str;
    while (*cp != 0) {
      if (*cp < ' ' || *cp > '~') {
#ifdef OS_WINNT
        str = GetSpecialWinCharacterReplacement ((unsigned char) *cp);
#else
        str = GetSpecialMacCharacterReplacement ((unsigned char) *cp);
#endif
        repl_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (replace_fmt) + StringLen (str)));
        sprintf (repl_str, replace_fmt, *cp, str == NULL ? "" : str);
        ValNodeAddPointer (&repl_list, 0, repl_str);
        if (str != NULL) {
          while (*str != 0) {
            *cp_dst = *str;
            cp_dst++;
            str++;
          }
        }
      } else {
        *cp_dst = *cp;
        cp_dst++;
      }
      cp++;
    }
    *text = MemFree (*text);
    *text = new_str;
  }
  return repl_list;
}


NLM_EXTERN void AutoReplaceSpecialCharactersWithMessage (CharPtr PNTR text)
{
  ValNodePtr list, vnp;

  list = AutoReplaceSpecialCharactersInText(text);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    Message (MSG_POSTERR, "%s", vnp->data.ptrvalue);
  }
  list = ValNodeFreeData (list);
}


NLM_EXTERN ValNodePtr AutoReplaceSpecialCharactersInTabTable (ValNodePtr row_list)
{
  ValNodePtr repl_list = NULL, col;
  CharPtr cp;

  while (row_list != NULL) {
    for (col = row_list->data.ptrvalue; col != NULL; col = col->next) {
      cp = col->data.ptrvalue;
      ValNodeLink (&repl_list, AutoReplaceSpecialCharactersInText(&cp));
      col->data.ptrvalue = cp;
    }
    row_list = row_list->next;
  }
  return repl_list;
}


NLM_EXTERN void AutoFixSpecialCharactersInEntity (Uint2 entityID)
{
  ValNodePtr bad_list = NULL, vnp, vnp_c;
  Char       label[2];
  CharPtr    repl;

  label[1] = 0;
  StringActionInEntity (entityID, FALSE, UPDATE_NEVER, NULL, NULL, NULL, TRUE,
                        SpecialCharFindWithContext, NULL, &bad_list);
  for (vnp = bad_list; vnp != NULL; vnp = vnp->next) 
  {
#ifdef OS_WINNT
    repl = GetSpecialWinCharacterReplacement ((unsigned char) vnp->choice);
#else
    repl = GetSpecialMacCharacterReplacement ((unsigned char) vnp->choice);
#endif
    label[0] = vnp->choice;
    Message (MSG_POSTERR, "Replaced '%s' with '%s'", label, repl == NULL ? "" : repl);
    for (vnp_c = vnp->data.ptrvalue; vnp_c != NULL; vnp_c = vnp_c->next)
    {
      FindReplaceString (vnp_c->data.ptrvalue, label, repl, TRUE, FALSE);
    }
  }
  bad_list = FreeContextList (bad_list);
}


/* Functions for reassigning affiliations of authors for Flu sequences */
typedef struct authaffil {
  CharPtr affil;
  ValNodePtr authors;
} AuthAffilData, PNTR AuthAffilPtr;


static AuthAffilPtr AuthAffilNew (CharPtr affil)
{
  AuthAffilPtr a;

  a = (AuthAffilPtr) MemNew (sizeof (AuthAffilData));
  a->affil = StringSave (affil);
  a->authors = NULL;
  return a;
}

static AuthAffilPtr AuthAffilFree (AuthAffilPtr a)
{
  if (a != NULL) {
    a->affil = MemFree (a->affil);
    a->authors = ValNodeFreeData (a->authors);
    a = MemFree (a);
  }
  return a;
}


static ValNodePtr AuthAffilListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = AuthAffilFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


typedef struct splitpub {
  BioseqPtr bsp;
  ValNodePtr auth_affil_list;
} SplitPubData, PNTR SplitPubPtr;


static SplitPubPtr SplitPubNew (BioseqPtr bsp)
{
  SplitPubPtr s;

  s = (SplitPubPtr) MemNew (sizeof (SplitPubData));
  s->bsp = bsp;
  s->auth_affil_list = NULL;
  return s;
}


static SplitPubPtr SplitPubFree (SplitPubPtr s)
{
  if (s != NULL) {
    s->auth_affil_list = AuthAffilListFree (s->auth_affil_list);
    s = MemFree (s);
  }
  return s;
}


NLM_EXTERN ValNodePtr SplitPubListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = SplitPubFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static int LIBCALLBACK SortBySplitPubTabTableRow (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ValNodePtr  col1, col2;
  int         rval = 0;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  col1 = vnp1->data.ptrvalue;
  col2 = vnp2->data.ptrvalue;
  while (col1 != NULL && col2 != NULL && rval == 0) {
    rval = StringCmp (col1->data.ptrvalue, col2->data.ptrvalue);
    col1 = col1->next;
    col2 = col2->next;
  }
  if (rval == 0) {
    if (col1 == NULL && col2 != NULL) {
      rval = -1;
    } else if (col1 != NULL && col2 == NULL) {
      rval = 1;
    }
  }
  return rval;
}


NLM_EXTERN ValNodePtr MakeSplitPubListFromTabList (ValNodePtr PNTR tab_table, SeqEntryPtr sep, ValNodePtr PNTR err_list)
{
  ValNodePtr  row, col;
  SeqIdPtr    sip;
  BioseqPtr   bsp = NULL;
  CharPtr     last_id = NULL;
  CharPtr     not_found_fmt = "Unable to locate sequence for %s";
  CharPtr     err_msg;
  ValNodePtr   split_pub_list = NULL;
  SplitPubPtr  s = NULL;
  AuthAffilPtr a = NULL;
  CharPtr      cp;

  if (tab_table == NULL || sep == NULL) {
    return NULL;
  }

  *tab_table = ValNodeSort (*tab_table, SortBySplitPubTabTableRow);

  for (row = *tab_table; row != NULL; row = row->next) {
    col = row->data.ptrvalue;    
    if (col == NULL) {
      continue;
    }
    if (last_id == NULL || StringCmp (last_id, col->data.ptrvalue) != 0) {
      sip = CreateSeqIdFromText (col->data.ptrvalue, sep);
      bsp = BioseqFind (sip);
      sip = SeqIdFree (sip);
      if (bsp == NULL) {
        if (err_list != NULL) {
          err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (not_found_fmt) + StringLen (col->data.ptrvalue)));
          sprintf (err_msg, not_found_fmt, col->data.ptrvalue);
          ValNodeAddPointer (err_list, 0, err_msg);
        }
        s = NULL;
      } else {
        s = SplitPubNew (bsp);
        ValNodeAddPointer (&split_pub_list, 0, s);
      }
      a = NULL;
      last_id = col->data.ptrvalue;
    }
    if (s == NULL) {
      continue;
    }
    col = col->next;
    if (col == NULL) {
      continue;
    }
    if (a == NULL || StringCmp (a->affil, col->data.ptrvalue) != 0) {
      a = AuthAffilNew (col->data.ptrvalue);
      ValNodeAddPointer (&(s->auth_affil_list), 0, a);
    }
    col = col->next;
    if (col != NULL && col->data.ptrvalue != NULL) {
      /* skip dr, mr, ms, mrs title */
      cp = col->data.ptrvalue;
      if (StringNICmp (cp, "dr ", 3) == 0
          || StringNICmp (cp, "dr.", 3) == 0
          || StringNICmp (cp, "mr ", 3) == 0
          || StringNICmp (cp, "mr.", 3) == 0
          || StringNICmp (cp, "ms ", 3) == 0
          || StringNICmp (cp, "ms.", 3) == 0) {
        cp += 3;
        cp += StringSpn (cp, " ");
      } else if (StringNICmp (cp, "mrs ", 4) == 0
                 || StringNICmp (cp, "mrs.", 4) == 0) {
        cp += 4;
        cp += StringSpn (cp, " ");
      } else if (StringNICmp (cp, "prof ", 5) == 0
                 || StringNICmp (cp, "prof.", 5) == 0) {
        cp += 4;
        cp += StringSpn (cp, " ");
      }
      ValNodeAddPointer (&a->authors, 0, StringSave (cp));
    }
  }
  return split_pub_list;
}


NLM_EXTERN AuthListPtr GetAuthorListForPub (PubPtr the_pub)
{
  CitGenPtr  cgp;
  CitSubPtr  csp;
  CitArtPtr  cap;
  CitBookPtr cbp;
  CitPatPtr  cpp;
  AuthListPtr alp = NULL;

  if (the_pub == NULL) return NULL;
  
  switch (the_pub->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      alp = cgp->authors;
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      alp = csp->authors;
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      alp = cap->authors;
      break;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      alp = cbp->authors;
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) the_pub->data.ptrvalue;
      alp = cpp->authors;
      break;
    default :
      break;
  }
  return alp;
}


static void ReplaceAuthorListForPub (PubPtr the_pub, AuthListPtr new_alp)
{
  CitGenPtr  cgp;
  CitSubPtr  csp;
  CitArtPtr  cap;
  CitBookPtr cbp;
  CitPatPtr  cpp;
  AuthListPtr alp = NULL;

  if (the_pub == NULL) return;
  
  switch (the_pub->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) the_pub->data.ptrvalue;
      alp = cgp->authors;
      cgp->authors = new_alp;
      break;
    case PUB_Sub :
      csp = (CitSubPtr) the_pub->data.ptrvalue;
      alp = csp->authors;
      csp->authors = new_alp;
      break;
    case PUB_Article :
      cap = (CitArtPtr) the_pub->data.ptrvalue;
      alp = cap->authors;
      cap->authors = new_alp;
      break;
    case PUB_Book :
    case PUB_Man :
      cbp = (CitBookPtr) the_pub->data.ptrvalue;
      alp = cbp->authors;
      cbp->authors = new_alp;
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) the_pub->data.ptrvalue;
      alp = cpp->authors;
      cpp->authors = new_alp;
      break;
    default :
      break;
  }
  alp = AuthListFree (alp);
}


static Boolean DoesStringMatchAuthor (CharPtr str, AuthorPtr author)
{
  Boolean    rval = FALSE;
  NameStdPtr nsp;
  CharPtr    cp;
  Int4       len;

  if (StringHasNoText (str) || author == NULL || author->name == NULL) {
    return FALSE;
  }
  switch (author->name->choice) {
    case 2:
      nsp = (NameStdPtr) author->name->data;
      if (nsp != NULL) {
        cp = StringRChr (str, ' ');
        if (cp == NULL) {
          cp = StringRChr (str, '.');
        }
        if (cp != NULL && StringCmp (cp + 1, nsp->names[0]) == 0) {
          len = StringCSpn (str, " .");
          if (StringNCmp (str, nsp->names[1], len) == 0) {
            rval = TRUE;
          }
        }
      }
      break;
    case 4:
    case 5:
      if (StringCmp (author->name->data, str) == 0) {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static AuthListPtr RemoveAuthorsFromAuthList (AuthListPtr alp, ValNodePtr auth_list)
{
  AuthListPtr removed = NULL;
  ValNodePtr  vnp_a, vnp_next, vnp_remove, vnp_prev;
  ValNodePtr  extract_list = NULL;

  if (alp == NULL || auth_list == NULL || alp->names == NULL || alp->choice != 1) {
    return NULL;
  }

  for (vnp_remove = auth_list; vnp_remove != NULL; vnp_remove = vnp_remove->next) {
    vnp_prev = NULL;
    for (vnp_a = alp->names; vnp_a != NULL; vnp_a = vnp_next) {
      vnp_next = vnp_a->next;
      if ((vnp_a->choice == 1 && DoesStringMatchAuthor (vnp_remove->data.ptrvalue, vnp_a->data.ptrvalue))
          || (vnp_a->choice == 2 && StringCmp (vnp_remove->data.ptrvalue, vnp_a->data.ptrvalue))) {
        if (vnp_prev == NULL) {
          alp->names = vnp_a->next;
        } else {
          vnp_prev->next = vnp_a->next;
        }
        vnp_a->next = NULL;
        ValNodeLink (&extract_list, vnp_a);
      } else {
        vnp_prev = vnp_a;
      }
    }
  }

  if (extract_list != NULL) {
    removed = AuthListNew ();
    removed->choice = 1;
    removed->names = extract_list;
  }
  return removed;
}


static AuthListPtr RemoveAuthorsFromPub (PubdescPtr pub, ValNodePtr auth_list)
{
  AuthListPtr alp, removed_alp;

  if (pub == NULL || auth_list == NULL) return NULL;
  alp = GetAuthorListForPub (pub->pub);
  removed_alp = RemoveAuthorsFromAuthList (alp, auth_list);
  return removed_alp;
}


typedef struct newpub {
  SeqDescrPtr PNTR descr;
  PubdescPtr pub;
} NewPubData, PNTR NewPubPtr;


static NewPubPtr NewPubNew (SeqDescrPtr PNTR descr, PubdescPtr pub)
{
  NewPubPtr n;

  n = (NewPubPtr) MemNew (sizeof (NewPubData));
  n->descr = descr;
  n->pub = pub;
  return n;
}


static void SplitPubsForBioseq (SplitPubPtr s)
{
  SeqDescrPtr       sdp, new_sdp;
  SeqDescrPtr PNTR  descr;
  SeqMgrDescContext context;
  AuthListPtr       alp;
  PubdescPtr        pubdesc, pubdesc_new;
  ValNodePtr        vnp;
  AuthAffilPtr      aa;
  ValNodePtr        new_pubs = NULL;
  ObjValNodePtr     ovp;
  BioseqPtr         p_bsp;
  BioseqSetPtr      p_bssp;
  NewPubPtr         n;

  if (s == NULL || s->bsp == NULL || s->auth_affil_list == NULL) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (s->bsp, NULL, Seq_descr_pub, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (s->bsp, sdp, Seq_descr_pub, &context)) {
    pubdesc = sdp->data.ptrvalue;
    for (vnp = s->auth_affil_list; vnp != NULL; vnp = vnp->next) {
      aa = (AuthAffilPtr) vnp->data.ptrvalue;
      if (aa != NULL) {
        alp = RemoveAuthorsFromPub (pubdesc, aa->authors);
        if (alp != NULL) {          
          pubdesc_new = AsnIoMemCopy (pubdesc, (AsnReadFunc) PubdescAsnRead, (AsnWriteFunc) PubdescAsnWrite);
          alp->affil = AffilNew ();
          alp->affil->affil = StringSave (aa->affil);
          ReplaceAuthorListForPub (pubdesc_new->pub, alp);
          descr = &(s->bsp->descr);
          if (sdp->extended == 1) {
            ovp = (ObjValNodePtr) sdp;
            if (ovp->idx.parenttype == OBJ_BIOSEQ) {
              p_bsp = (BioseqPtr) ovp->idx.parentptr;
              descr = &(p_bsp->descr);
            } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
              p_bssp = (BioseqSetPtr) ovp->idx.parentptr;
              descr = &(p_bssp->descr);
            }
          }
          ValNodeAddPointer (&new_pubs, 0, NewPubNew (descr, pubdesc_new));
        }
      }
    }
  }
  for (vnp = new_pubs; vnp != NULL; vnp = vnp->next) {
    n = (NewPubPtr) vnp->data.ptrvalue;

    new_sdp = SeqDescrNew (*(n->descr));
    new_sdp->choice = Seq_descr_pub;
    new_sdp->data.ptrvalue = n->pub;
    if (n->descr == NULL) {
      *(n->descr) = new_sdp;
    }
  }
  new_pubs = ValNodeFreeData (new_pubs);

}


NLM_EXTERN void SplitPubsByList (ValNodePtr split_list)
{
  ValNodePtr vnp;

  for (vnp = split_list; vnp != NULL; vnp = vnp->next) {
    SplitPubsForBioseq (vnp->data.ptrvalue);
  }
}


static void AddStructuredCommentCallback (BioseqPtr bsp, Pointer data)
{
  UserObjectPtr uop;
  SeqDescrPtr   sdp;

  if (bsp == NULL || ISA_aa (bsp->mol) || (uop = (UserObjectPtr) data) == NULL) {
    return;
  }

  sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_user);
  sdp->data.ptrvalue = AsnIoMemCopy (uop, (AsnReadFunc) UserObjectAsnRead, (AsnWriteFunc) UserObjectAsnWrite);
}


static CharPtr official_prefix_list[] = {
  "HIVDataBaseData",
  "MIGS-Data",
  "MIMS-Data",
  "MIENS-Data",
  "MIGS:3.0-Data",
  "GISAID_EpiFlu(TM)Data",
  "FluData",
  "EpifluData",
  "International Barcode of Life (iBOL)Data", 
  "Assembly-Data",
  "Genome-Assembly-Data",
  NULL
};


NLM_EXTERN ValNodePtr GetStructuredCommentPrefixList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; official_prefix_list[i] != NULL; i++) {
    ValNodeAddPointer (&list, 0, StringSave (official_prefix_list[i]));
  }
  return list;
}


static Int4 GetDbnameCoreLen (CharPtr dbname)
{
  Int4 len = StringLen (dbname);
  if (len > 4 && StringICmp (dbname + len - 4, "Data") == 0) {
    len -= 4;
  }
  if (len > 1 && StringNICmp (dbname + len - 1, "-", 1) == 0) {
    len -= 1;
  }
  return len;
}


static CharPtr MatchesOfficialStructuredCommentDbname (CharPtr dbname)
{
  Int4 i;
  Int4 len_orig;
  Int4 len_can;

  len_orig = GetDbnameCoreLen (dbname);
  for (i = 0; official_prefix_list[i] != NULL; i++) {
    len_can = GetDbnameCoreLen (official_prefix_list[i]);
    if (len_orig == len_can && StringNICmp (dbname, official_prefix_list[i], len_orig) == 0) {
      return official_prefix_list[i];
    }
  }
  if (StringNICmp (dbname, "HIV-Database", len_orig) == 0) {
    return "HIVDatabase";
  }
  return NULL;
}


NLM_EXTERN CharPtr StructuredCommentDbnameFromString (CharPtr string)
{
  CharPtr dbname, tmp;
  Int4    len;

  if (StringHasNoText (string)) {
    return NULL;
  }

  dbname = StringSave (string + StringSpn (string, "##"));
  len = StringLen (dbname);
  if (len > 2 && StringCmp (dbname + len - 2, "##") == 0) {
    dbname[len - 2] = 0;
    len -= 2;
  }
  if (len > 6 && StringCmp (dbname + len - 6, "-START") == 0) {
    dbname[len - 6] = 0;
    len -= 6;
  }
  if (len > 6 && StringCmp (dbname + len - 4, "-END") == 0) {
    dbname[len - 4] = 0;
    len -= 4;
  }

  /* correct for weirdnesses with -data for recognizable prefixes */
  tmp = MatchesOfficialStructuredCommentDbname (dbname);
  if (tmp != NULL) {
    dbname = MemFree (dbname);
    dbname = StringSave (tmp);
  }
  return dbname;
}


static CharPtr MakeStructuredCommentPrefixFromString (CharPtr orig)
{
  CharPtr    core, new_prefix;
  Int4       core_len;

  if (StringHasNoText (orig)) {
    return StringSave ("##Metadata-START##");
  }

  core = StructuredCommentDbnameFromString(orig);
  core_len = StringLen (core);
  
  new_prefix = (CharPtr) MemNew (sizeof (Char) * (11 + core_len));
  StringCpy (new_prefix, "##");
  StringNCat (new_prefix, core, core_len);
  StringCat (new_prefix, "-START##");
  core = MemFree (core);
  return new_prefix;
}


static CharPtr MakeStructuredCommentSuffixFromString (CharPtr orig)
{
  CharPtr    core, new_suffix;
  Int4       core_len;

  if (StringHasNoText (orig)) {
    return StringSave ("##Metadata-END##");
  }

  core = StructuredCommentDbnameFromString(orig);
  core_len = StringLen (core);
  
  new_suffix = (CharPtr) MemNew (sizeof (Char) * (9 + core_len));
  StringCpy (new_suffix, "##");
  StringNCat (new_suffix, core, core_len);
  StringCat (new_suffix, "-END##");
  core = MemFree (core);
  return new_suffix;
}


/* This function reads in a tab-delimited table.  The first line is a header.
 * If not apply_to_all, the first column must contain sequence IDs.  The remaining cells in the first
 * line are the names of fields to create in structured comments.
 */
NLM_EXTERN ValNodePtr CreateStructuredCommentsFromRow (ValNodePtr header, ValNodePtr values, CharPtr id_str, ValNodePtr PNTR err_list)
{
  ValNodePtr comment_list = NULL;
  ValNodePtr vnp_h, vnp_l;
  UserObjectPtr uop = NULL;
  CharPtr       suffix = NULL, fix, msg;
  CharPtr       extra_data_fmt = "Too many fields for sequence %s";

  if (header == NULL || values == NULL) {
    return NULL;
  }

  vnp_h = header;
  vnp_l = values;

  while (vnp_h != NULL && vnp_l != NULL) {
    if (!StringHasNoText (vnp_l->data.ptrvalue)) {
      if (StringICmp (vnp_h->data.ptrvalue, "StructuredCommentPrefix") == 0) {
        if (suffix != NULL) {
          fix = MakeStructuredCommentSuffixFromString (suffix);
          AddItemStructuredCommentUserObject (uop, "StructuredCommentSuffix", fix);
          fix = MemFree (fix);
          suffix = MemFree (suffix);
        }
        uop = CreateStructuredCommentUserObject (NULL, NULL);
        ValNodeAddPointer (&comment_list, 0, uop);
        suffix = StringSave (vnp_l->data.ptrvalue);
        fix = MakeStructuredCommentPrefixFromString (suffix);
        AddItemStructuredCommentUserObject (uop, vnp_h->data.ptrvalue, fix);
        fix = MemFree (fix);
      } else if (StringICmp (vnp_h->data.ptrvalue, "StructuredCommentSuffix") == 0) {
        fix = MakeStructuredCommentSuffixFromString (vnp_l->data.ptrvalue);
        AddItemStructuredCommentUserObject (uop, vnp_h->data.ptrvalue, fix);
        fix = MemFree (fix);
        suffix = MemFree (suffix);
        uop = NULL;
      } else {
        if (uop == NULL) {
          uop = CreateStructuredCommentUserObject (NULL, NULL);
          ValNodeAddPointer (&comment_list, 0, uop);
        }
        AddItemStructuredCommentUserObject (uop, vnp_h->data.ptrvalue, vnp_l->data.ptrvalue);
      }
    }
    vnp_h = vnp_h->next;
    vnp_l = vnp_l->next;
  }
  if (uop != NULL && suffix != NULL) {
    fix = MakeStructuredCommentSuffixFromString (suffix);
    AddItemStructuredCommentUserObject (uop, "StructuredCommentSuffix", fix);
    fix = MemFree (fix);
  }
  suffix = MemFree (suffix);

  if (err_list != NULL && id_str != NULL) {
    while (vnp_l != NULL && StringHasNoText (vnp_l->data.ptrvalue)) {
      vnp_l = vnp_l->next;
    }
    if (vnp_l != NULL) {
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (extra_data_fmt) + StringLen (id_str)));
      sprintf (msg, extra_data_fmt, id_str);
      ValNodeAddPointer (err_list, 0, msg);
    }
  }

  return comment_list;
}


NLM_EXTERN void CreateStructuredCommentsForAllFromTable (SeqEntryPtr sep, ValNodePtr header, ValNodePtr line, ValNodePtr PNTR err_list)
{
  ValNodePtr tmp, vnp_l;
  UserObjectPtr uop;

  while (line != NULL) {
    tmp = CreateStructuredCommentsFromRow (header, line->data.ptrvalue, NULL, err_list);
    for (vnp_l = tmp; vnp_l != NULL; vnp_l = vnp_l->next) {
      uop = (UserObjectPtr) vnp_l->data.ptrvalue;
      VisitBioseqsInSep (sep, uop, AddStructuredCommentCallback);
      uop = UserObjectFree (uop);
    }
    tmp = ValNodeFree (tmp);
    line = line->next;
  }
}


NLM_EXTERN ValNodePtr CreateStructuredCommentsFromFile (FILE *fp, SeqEntryPtr sep, Boolean apply_to_all)
{
  ValNodePtr err_list = NULL;
  ValNodePtr table, header, line, vnp_h, vnp_l, tmp;
  SeqIdPtr   sip;
  CharPtr    id_str;
  CharPtr    bad_id_fmt = "Unable to find sequence for %s";
  CharPtr    msg;
  BioseqPtr  bsp;
  SeqDescrPtr sdp;

  if (fp == NULL || sep == NULL) {
    return NULL;
  }
  table = ReadTabTableFromFile (fp);
  if (table == NULL || table->next == NULL || table->data.ptrvalue == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("Unable to read table from file"));
    table = FreeTabTable (table);
    return err_list;
  }
  if (apply_to_all) {
    tmp = FlipTabTableAxes (table);
    table = FreeTabTable (table);
    table = tmp;
  }

  header = table->data.ptrvalue;
  if (header == NULL || header->data.ptrvalue == NULL || header->next == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("Bad header line"));
    table = FreeTabTable (table);
    return err_list;
  }
  line = table->next;

  if (apply_to_all) {
    CreateStructuredCommentsForAllFromTable (sep, header, line, &err_list);
  } else {
    while (line != NULL) {
      vnp_h = header;
      vnp_l = line->data.ptrvalue;
      if (vnp_l != NULL)  {
        id_str = vnp_l->data.ptrvalue;
        sip = CreateSeqIdFromText (id_str, sep);
        if (sip == NULL || (bsp = BioseqFind (sip)) == NULL) {
          msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_id_fmt) + StringLen (id_str)));
          sprintf (msg, bad_id_fmt, id_str);
          ValNodeAddPointer (&err_list, 0, msg);
        } else {
          tmp = CreateStructuredCommentsFromRow (header->next, vnp_l->next, id_str, &err_list);
          for (vnp_l = tmp; vnp_l != NULL; vnp_l = vnp_l->next) {
            sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_user);
            sdp->data.ptrvalue = vnp_l->data.ptrvalue;
          }
          tmp = ValNodeFree (tmp);
        }
      }
      line = line->next;
    }
  }
  return err_list;
}


NLM_EXTERN void AddDatabaseNameToStructuredComment (UserObjectPtr uop, CharPtr dbname)
{
  UserFieldPtr   curr;
  Boolean        hasPrefix = FALSE;
  Boolean        hasSuffix = FALSE;
  ObjectIdPtr    oip;
  CharPtr        prefix_fmt = "##%sData-START##";
  CharPtr        suffix_fmt = "##%sData-END##";
  CharPtr        prefix, suffix;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "StructuredComment") != 0) return;

  if (StringHasNoText (dbname)) {
    dbname = "Meta";
  } 

  prefix = (CharPtr) MemNew (sizeof (Char) * (StringLen (prefix_fmt) + StringLen(dbname)));
  sprintf (prefix, prefix_fmt, dbname);
  suffix = (CharPtr) MemNew (sizeof (Char) * (StringLen (suffix_fmt) + StringLen(dbname)));
  sprintf (suffix, suffix_fmt, dbname);

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "StructuredCommentPrefix") == 0) {
      hasPrefix = TRUE;
      if (curr->choice == 1) {
        MemFree (curr->data.ptrvalue);
        curr->data.ptrvalue = (Pointer) StringSave (prefix);
      }
    }
  }
  if (! hasPrefix) {
    AddItemStructuredCommentUserObject (uop, "StructuredCommentPrefix", prefix);
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "StructuredCommentSuffix") == 0) {
      hasSuffix = TRUE;
      if (curr->choice == 1) {
        MemFree (curr->data.ptrvalue);
        curr->data.ptrvalue = (Pointer) StringSave (suffix);
      }
    }
  }
  if (! hasSuffix) {
    AddItemStructuredCommentUserObject (uop, "StructuredCommentSuffix", suffix);
  }
  prefix = MemFree (prefix);
  suffix = MemFree (suffix);
}


static ValNodePtr RowFromStructuredComment (UserObjectPtr uop, ValNodePtr PNTR header)
{
  UserFieldPtr   curr;
  ObjectIdPtr    oip;
  ValNodePtr     vnp_h, vnp_v;
  ValNodePtr     values = NULL;
  CharPtr        label;

  if (uop == NULL || uop->type == NULL 
      || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return NULL;
  }

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL) {
      label = GetObjectIdString(oip);
      for (vnp_h = *header, vnp_v = values; 
           vnp_h != NULL && (StringCmp (oip->str, vnp_h->data.ptrvalue) != 0 || (vnp_v != NULL && vnp_v->data.ptrvalue != NULL)); 
           vnp_h = vnp_h->next, vnp_v = vnp_v->next) {
        if (vnp_v == NULL) {
          vnp_v = ValNodeNew (values);
          if (values == NULL) {
            values = vnp_v;
          }
        }        
      }
      if (vnp_h == NULL) {
        vnp_h = ValNodeNew (*header);
        if (*header == NULL) {
          *header = vnp_h;
        }
        vnp_h->data.ptrvalue = label;
        label = NULL;
      }
      label = MemFree (label);
      if (vnp_v == NULL) {
        vnp_v = ValNodeNew (values);
        if (values == NULL) {
          values = vnp_v;
        }
      }
      vnp_v->data.ptrvalue = StringSave (curr->data.ptrvalue);
    }
  }
  return values;
}


static void GetStructuredCommentsForBioseq(BioseqPtr bsp, Pointer data)
{
  SeqMgrDescContext context;
  SeqDescPtr sdp;
  ValNodePtr header = NULL;
  ValNodePtr list = NULL;
  ValNodePtr PNTR table = NULL;
  ValNodePtr vnp;
  Char       id_txt[200];

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  table = (ValNodePtr PNTR) data;
  if (*table != NULL) {
    header = (*table)->data.ptrvalue;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) {
    ValNodeLink (&list, RowFromStructuredComment (sdp->data.ptrvalue, &header));
  }

  if (list != NULL) {
    if (*table == NULL) {
      vnp = ValNodeNew (NULL);
      vnp->data.ptrvalue = StringSave ("SeqId");
      vnp->next = header;
      header = vnp;
      ValNodeAddPointer (table, 0, header);
      vnp = ValNodeNew (NULL);
      SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
      vnp->data.ptrvalue = StringSave (id_txt);
      vnp->next = list;
      list = vnp;
    } else {      
      (*table)->data.ptrvalue = header;
      /* placeholder for SeqId already exists */
      SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
      list->data.ptrvalue = StringSave (id_txt);
    }
    ValNodeAddPointer (table, 0, list);
  }
}


NLM_EXTERN ValNodePtr CreateStructuredCommentTableFromSeqEntry (SeqEntryPtr sep)
{
  ValNodePtr table = NULL;

  VisitBioseqsInSep (sep, &table, GetStructuredCommentsForBioseq);
  return table;
}


static SeqPortPtr SeqPortFromAlignmentInterval (Int4 seqstart, Int4 seqstop, Uint1 strand, BioseqPtr bsp)
{
  SeqIntPtr  sinp;
  SeqLocPtr  slp;
  SeqPortPtr spp;

  if (bsp == NULL || seqstart >= bsp->length - 1) return NULL;
  seqstop = MIN (bsp->length -1, seqstop);
  sinp = SeqIntNew();
  if (sinp == NULL) return NULL;
  sinp->from = seqstart;
  sinp->to = seqstop;
  sinp->strand = strand;
  sinp->id = SeqIdDup (SeqIdFindBest (bsp->id, 0));
  slp = ValNodeNew (NULL);
  if (slp == NULL) {
    SeqIntFree (sinp);
    return NULL;
  }
  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = (Pointer) sinp;
  spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
  SeqLocFree (slp);
  return spp;
}


static void SetSequenceIntervalBuf
(SeqAlignPtr salp,
 BioseqPtr   bsp,
 Int4        row,
 Int4        start,                             
 Int4        stop,
 Int4Ptr     seqstart,
 Int4Ptr     seqstop,
 Int4        aln_len,                             
 Uint1Ptr    target_buf)
{
  Int4       buf_len = stop - start + 1;
  Uint1      strand;
  Int4       i;
  SeqPortPtr spp;

  if (seqstart == NULL || seqstop == NULL)
  {
    return;
  }
  
  *seqstart = ALNMGR_GAP;
  *seqstop = ALNMGR_GAP;
  if (bsp == NULL)
  {
    return;
  }
  strand = SeqAlignStrand (salp, row - 1);
  MemSet (target_buf, 0, buf_len);
  /* if this is a minus strand sequence, start is stop and stop is start */
  if (strand == Seq_strand_minus) {
    *seqstop = AlnMgr2MapSeqAlignToBioseq(salp, start, row);
    *seqstart  = AlnMgr2MapSeqAlignToBioseq(salp, stop, row);
  } else {
    *seqstart = AlnMgr2MapSeqAlignToBioseq(salp, start, row);
    *seqstop  = AlnMgr2MapSeqAlignToBioseq(salp, stop, row);
  }
  
  if (strand == Seq_strand_minus) {
    i = stop;
    while ((*seqstart == ALNMGR_GAP || *seqstart == ALNMGR_ROW_UNDEFINED) && i > 0) { /* count backward if we are in the gap */
      i--;
      *seqstart = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
  } else {
    i = start;
    while ((*seqstart == ALNMGR_GAP || *seqstart == ALNMGR_ROW_UNDEFINED) && i < aln_len) { /* count forward if we in the gap */
      i++;
      *seqstart = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
  }
  
  if (*seqstop < 0 || *seqstop>=bsp->length) *seqstop = bsp->length - 1;  /* -1 means exeed sequence length */
  
  if (*seqstop > -1 && *seqstart > -1 && *seqstop - *seqstart > stop - start) {
    *seqstop = *seqstart + stop - start;
  }


  if (strand == Seq_strand_minus) {
    i = start;
    while (*seqstop == ALNMGR_GAP && i > 0) { /* count backward if we are in the gap */
      i--;
      *seqstop = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
  } else {
    i = stop;
    while (*seqstop == ALNMGR_GAP && i < aln_len) { /* count forward if we are in the gap */
      i++;
      *seqstop = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
  }
  
  if (*seqstart == ALNMGR_GAP  &&  *seqstop == ALNMGR_GAP) {
    return;
  }
  if (*seqstop  < 0) *seqstop  = bsp->length - 1;
  if (*seqstart < 0) *seqstart = *seqstop;
  if (*seqstop < *seqstart) {
    *seqstop = *seqstart = 0;
  }
  if (strand == Seq_strand_minus) {
    if (*seqstop - *seqstart > buf_len) 
      *seqstart = *seqstop - buf_len;
  } else {
    if (*seqstop - *seqstart > buf_len) *seqstop = *seqstart + buf_len;  /* not to exeed the current line */
  }

  spp = SeqPortFromAlignmentInterval (*seqstart, *seqstop, strand, bsp);
  SeqPortRead  (spp, target_buf, *seqstop - *seqstart + 1);
  SeqPortFree  (spp);
}


NLM_EXTERN void 
AlignmentIntervalToString 
(SeqAlignPtr salp,
 Int4        row,
 Int4        start,
 Int4        stop,
 Int4        target_row,
 Boolean     view_whole_entity,
 Uint1Ptr    seqbuf,
 Uint1Ptr    alnbuf,
 Int4 PNTR   alnbuffer_len,
 Boolean     show_substitutions)
{
  Int4       aln_len = AlnMgr2GetAlnLength(salp, FALSE);
  SeqIdPtr   sip     = AlnMgr2GetNthSeqIdPtr(salp, row);
  BioseqPtr  bsp     = BioseqLockById(sip);
  Int4       alnbuf_len = stop - start + 1;
  Uint1      strand;
  Int4       seqstart, seqstop;
  Int4       i, k;
  SeqPortPtr spp;
  Int4       seq_len;
  Uint1      target_strand;
  SeqIdPtr   sip_target;
  BioseqPtr  bsp_target;
  Int4       target_start;
  Int4       target_stop;
  Uint1Ptr   target_buf;
  Int4       aln_pos;

  MemSet(alnbuf, '-', alnbuf_len); /* assume all gaps and fill the sequence later */
  MemSet(seqbuf, 0, alnbuf_len);
  if (target_row < 0 || bsp == NULL)
  {
    BioseqUnlock (bsp);
    SeqIdFree    (sip);
    return;
  }
  
  if (stop > aln_len && start > aln_len)
  {
    BioseqUnlock (bsp);
    SeqIdFree    (sip);
    return;
  }

  if (stop > aln_len) {
    MemSet (alnbuf + aln_len - start, 0, stop - aln_len);
    stop = aln_len - 1;
    alnbuf_len = stop - start + 1;
  }

  if (alnbuffer_len != NULL) {
    *alnbuffer_len = alnbuf_len;
  }

  strand = SeqAlignStrand (salp, row - 1);
  target_strand = SeqAlignStrand (salp, target_row - 1);
  /* if this is a minus strand sequence, start is stop and stop is start */
  if (strand == Seq_strand_minus) {
    seqstop = AlnMgr2MapSeqAlignToBioseq(salp, start, row);
    seqstart  = AlnMgr2MapSeqAlignToBioseq(salp, stop,  row);
  } else {
    seqstart = AlnMgr2MapSeqAlignToBioseq(salp, start, row);
    seqstop  = AlnMgr2MapSeqAlignToBioseq(salp, stop,  row);
  }
  
  if (strand == Seq_strand_minus) {
    i = stop;
    while ((seqstart == ALNMGR_GAP || seqstart == ALNMGR_ROW_UNDEFINED) && i > 0) { /* count backward if we are in the gap */
      i--;
      seqstart = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
  } else {
    i = start;
    while ((seqstart == ALNMGR_GAP || seqstart == ALNMGR_ROW_UNDEFINED) && i < aln_len) { /* count forward if we in the gap */
      i++;
      seqstart = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
  }
  
  if (seqstop == -1 || seqstop>=bsp->length)
  {
    seqstop = bsp->length - 1;  /* -1 means exeed sequence length */
  }
  
  if (strand == Seq_strand_minus) {
    i = start;
    while (seqstop == ALNMGR_GAP && i > 0) { /* count backward if we are in the gap */
      i--;
      seqstop = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
    if (i == 0) {
      /* gap goes to beginning of sequence, count forward until we are no longer in the gap */
      i = start;
      while (seqstop < 0 && i < stop) {
        i++;
        seqstop = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
      }
    }
  } else {
    i = stop;
    while (seqstop < 0 && i < aln_len) { /* count forward if we are in the gap */
      i++;
      seqstop = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
    }
    if (i == aln_len) {
      /* gap goes to end of sequence, count backwards until we are no longer in the gap */
      i = stop;
      while (seqstop < 0 && i > start) {
        i--;
        seqstop = AlnMgr2MapSeqAlignToBioseq(salp, i, row);
      }
    }
  }
  
  if (seqstart == ALNMGR_GAP  &&  seqstop == ALNMGR_GAP) seqstart = seqstop = 0;  /* whole line are gaps */
  if (seqstop < seqstart) {
    seqstart = seqstop = 0; /* treat whole line as gap */
  }
  if (seqstop  < 0) seqstop  = bsp->length - 1;
  if (seqstart < 0) seqstart = seqstop;
  if (strand == Seq_strand_minus) {
    if (seqstop - seqstart > alnbuf_len)
    {
      seqstart = seqstop - alnbuf_len;
    }
  } else {
    if (seqstop - seqstart > alnbuf_len) 
    {
      seqstop = seqstart + alnbuf_len;  /* not to exeed the current line */
    }
  }

  spp = SeqPortFromAlignmentInterval (seqstart, seqstop, strand, bsp);
  SeqPortRead  (spp, seqbuf, seqstop - seqstart + 1);
  if (seqbuf [stop - start] == 0) {
    seq_len = StringLen ((CharPtr) seqbuf);
  } else {
    seq_len = stop - start + 1;
  }
  SeqPortFree  (spp);
  BioseqUnlock (bsp);
  SeqIdFree    (sip);

  if (row != target_row  &&  ! view_whole_entity  &&  target_row != ALNMGR_ROW_UNDEFINED)  {
    sip_target = AlnMgr2GetNthSeqIdPtr(salp, target_row);
    bsp_target = BioseqLockById(sip_target);

    target_buf = (Uint1Ptr) MemNew (stop - start + 1);
    MemSet (target_buf, 0, stop - start + 1);
    SetSequenceIntervalBuf (salp, bsp_target, target_row, start, stop, 
                              &target_start, &target_stop, aln_len, target_buf);
  } else {
    sip_target = NULL;
    bsp_target = NULL;
    target_buf = NULL;
  }

  k = 0;
  i = 0;

  for (aln_pos = start; aln_pos <= stop; aln_pos ++) {
    Int4 seq_pos = AlnMgr2MapSeqAlignToBioseq(salp, aln_pos, row);
    Int4 target_pos = AlnMgr2MapSeqAlignToBioseq(salp, aln_pos, target_row);

    if (seq_pos >= 0 && (seq_pos < seqstart || seq_pos > seqstop)) {
      seq_pos = -1;
    }
    if (seq_pos >= 0) {
      alnbuf [aln_pos - start] = TO_LOWER (seqbuf[k]);
      if (show_substitutions)
      {
        /* Handle mismatches (insert dots when matched) */
        if (row != target_row  &&  ! view_whole_entity  &&  target_row != ALNMGR_ROW_UNDEFINED)  {
          if(target_pos >= 0  && target_pos < bsp_target->length) { /* no gap in the target sequence */
            if (seqbuf[k] == target_buf[i]) {
              alnbuf[aln_pos - start] = '.';
            }
          }
        }   /* mismatches */
      }
      k++;
    }
    if (target_pos >= 0) {
      i++;
    }
  }    

  if (alnbuf[alnbuf_len] == 0 && alnbuffer_len != NULL) {
    *alnbuffer_len = StringLen ((CharPtr) alnbuf);
  }

  if (bsp_target != NULL) {
    BioseqUnlock (bsp_target);
  }
  if (sip_target != NULL) {
    SeqIdFree (sip_target);
  }
  if (target_buf != NULL) {
    MemFree (target_buf);
  }
}


NLM_EXTERN void SetDescriptorPropagate (BioseqSetPtr bssp)
{
  BioseqPtr         bsp;
  SeqEntryPtr       seqentry;
  ValNodePtr        sourcedescr;

  if (bssp != NULL) {
    sourcedescr = bssp->descr;
    if (sourcedescr != NULL) {
      bssp->descr = NULL;
      seqentry = bssp->seq_set;
      while (seqentry != NULL) {
        if (seqentry->data.ptrvalue != NULL) {
          if (seqentry->choice == 1) {
            bsp = (BioseqPtr) seqentry->data.ptrvalue;
            ValNodeLink (&(bsp->descr),
                         AsnIoMemCopy ((Pointer) sourcedescr,
                                       (AsnReadFunc) SeqDescrAsnRead,
                                       (AsnWriteFunc) SeqDescrAsnWrite));
          } else if (seqentry->choice == 2) {
            bssp = (BioseqSetPtr) seqentry->data.ptrvalue;
            ValNodeLink (&(bssp->descr),
                         AsnIoMemCopy ((Pointer) sourcedescr,
                                       (AsnReadFunc) SeqDescrAsnRead,
                                       (AsnWriteFunc) SeqDescrAsnWrite));
          }
        }
        seqentry = seqentry->next;
      }
      SeqDescrFree (sourcedescr);
    }
  }
}


/* This function will look for nested sets of the same type
 * and remove the inner set.
 */
static Boolean RemoveDuplicateNestedSetsInSeqEntry (SeqEntryPtr top_sep)
{
  BioseqSetPtr bssp, lower_bssp;
  SeqEntryPtr  sep, sep_next, sep_tmp, sep_prev = NULL;
  SeqDescrPtr  last_sdp;
  SeqAnnotPtr  last_sap;
  Boolean      rval = FALSE;

  if (top_sep == NULL || !IS_Bioseq_set (top_sep) 
      || (bssp = (BioseqSetPtr) top_sep->data.ptrvalue) == NULL
      || bssp->seq_set == NULL) {
    return FALSE;
  }
  
  sep = bssp->seq_set;
  while (sep != NULL) {
    sep_next = sep->next;
    rval |= RemoveDuplicateNestedSetsInSeqEntry (sep);
    if (IS_Bioseq_set (sep) 
        && (lower_bssp = (BioseqSetPtr) sep->data.ptrvalue) != NULL
        && bssp->_class == lower_bssp->_class) {
      /* if this is the only set, move the descriptors up, otherwise
       * propagate the descriptors down.
       */
      if (sep->next == NULL && sep == bssp->seq_set) {
        if (bssp->descr == NULL) {
          bssp->descr = lower_bssp->descr;
        } else {
          last_sdp = bssp->descr;
          while (last_sdp->next != NULL) {
            last_sdp = last_sdp->next;
          }
          last_sdp->next = lower_bssp->descr;
        }
        lower_bssp->descr = NULL;
      } else {
        SetDescriptorPropagate (lower_bssp);
      }
      /* copy annotations to parent */
      if (bssp->annot == NULL) {
        bssp->annot = lower_bssp->annot;
      } else {
        last_sap = bssp->annot;
        while (last_sap->next != NULL) {
          last_sap = last_sap->next;
        }
        last_sap->next = lower_bssp->annot;
      }
      lower_bssp->annot = NULL;

      /* insert members of lower set in this position in upper set */
      if (lower_bssp->seq_set == NULL) {
        if (sep_prev == NULL) {
          bssp->seq_set = sep_next;
        } else {
          sep_prev->next = sep_next;
        }
      } else {
        if (sep_prev == NULL) {
          bssp->seq_set = lower_bssp->seq_set;
        } else {
          sep_prev->next = lower_bssp->seq_set;
        }
        sep_tmp = lower_bssp->seq_set;
        sep_prev = sep_tmp;
        while (sep_tmp->next != NULL) {
          sep_prev = sep_tmp;
          sep_tmp = sep_tmp->next;
        }
        sep_tmp->next = sep_next;
        lower_bssp->seq_set = NULL;
      }
      sep->next = NULL;
      sep = SeqEntryFree (sep);
      rval = TRUE;
    } else {
      sep_prev = sep;
    }
    sep = sep_next;
  }
  return rval;
}


NLM_EXTERN Boolean RemoveDuplicateNestedSetsForEntityID (Uint2 entityID)
{
  SeqEntryPtr       top_sep;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             top_parenttype;
  Pointer           top_parentptr;
  Boolean           rval = FALSE;

  top_sep = GetTopSeqEntryForEntityID (entityID);
  if (top_sep == NULL) return FALSE;

  SaveSeqEntryObjMgrData (top_sep, &omdptop, &omdata);
  GetSeqEntryParent (top_sep, &top_parentptr, &top_parenttype);

  rval = RemoveDuplicateNestedSetsInSeqEntry(top_sep);

  SeqMgrLinkSeqEntry (top_sep, top_parenttype, top_parentptr);
    
  SeqMgrClearFeatureIndexes (entityID, NULL);
  SeqMgrIndexFeatures (entityID, NULL);

  RestoreSeqEntryObjMgrData (top_sep, omdptop, &omdata);
  NormalizeDescriptorOrder (top_sep);
  
  SeqMgrClearFeatureIndexes (entityID, NULL);
  SeqMgrIndexFeatures (entityID, NULL);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);    
  return rval;
}


NLM_EXTERN CharPtr KeywordForStructuredCommentName (UserObjectPtr uop)
{
  UserFieldPtr ufp;
  CharPtr prefix = NULL;
  CharPtr keyword = NULL;

  if (uop == NULL) {
    return NULL;
  }

  for (ufp = uop->data; ufp != NULL && prefix == NULL; ufp = ufp->next) {
    if (ufp->label != NULL 
        && StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
        && ufp->choice == 1) {
      prefix = ufp->data.ptrvalue;
    }
  }

  if (prefix == NULL) {
    keyword = NULL;
  } else if (StringCmp (prefix, "##MIGS-Data-START##") == 0) {
    keyword = StringSave ("GSC:MIGS:2.1");
  } else if (StringCmp (prefix, "##MIMS-Data-START##") == 0) {
    keyword = StringSave ("GSC:MIMS:2.1");
  } else if (StringCmp (prefix, "##MIENS-Data-START##") == 0) {
    keyword = StringSave ("GSC:MIENS:2.1");
  } else if (StringCmp (prefix, "##MIGS:3.0-Data-START##") == 0) {
    keyword = StringSave ("GSC:MIxS;MIGS:3.0");
  } else if (StringCmp (prefix, "##MIMS:3.0-Data-START##") == 0) {
    keyword = StringSave ("GSC:MIxS;MIMS:3.0");
  } else if (StringCmp (prefix, "##MIMARKS:3.0-Data-START##") == 0) {
    keyword = StringSave ("GSC:MIxS;MIMARKS:3.0");
  }

  return keyword;
}


static Boolean HasKeyword (BioseqPtr bsp, CharPtr keyword)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  GBBlockPtr gb;
  Boolean    has_keyword = FALSE;
  ValNodePtr vnp;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
       sdp != NULL && !has_keyword;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &context)) {
    if ((gb = (GBBlockPtr) sdp->data.ptrvalue) != NULL) {
      for (vnp = gb->keywords; vnp != NULL && !has_keyword; vnp = vnp->next) {
        if (StringCmp (vnp->data.ptrvalue, keyword) == 0) {
          has_keyword = TRUE;
        }
      }
    }
  }
  return has_keyword;
}


static void AddKeywordToBioseq (BioseqPtr bsp, CharPtr keyword)
{
  SeqDescPtr sdp;
  SeqEntryPtr sep;
  GBBlockPtr gb;
  ValNodePtr vnp;

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    return;
  }
  sdp = GetDescrOnSeqEntry (sep, Seq_descr_genbank);
  if (sdp == NULL) {
    sdp = NewDescrOnSeqEntry (sep, Seq_descr_genbank);
    if (sdp != NULL) {
      sdp->data.ptrvalue = (Pointer) GBBlockNew ();
    }
  }
  if (sdp == NULL) return;
  gb = (GBBlockPtr) sdp->data.ptrvalue;
  if (gb == NULL) {
    gb = GBBlockNew ();
    sdp->data.ptrvalue = gb;
  }
  if (gb == NULL) return;

  for (vnp = gb->keywords; vnp; vnp = vnp->next) {
    if (StringCmp((CharPtr)vnp->data.ptrvalue, keyword) == 0) {
      return;
    }
  }
  ValNodeAddPointer (&(gb->keywords), 0, StringSave (keyword));
}


static void RemoveKeywordFromBioseq (BioseqPtr bsp, CharPtr keyword)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  GBBlockPtr gb;
  ValNodePtr vnp, vnp_next, prev;
  ObjValNodePtr ovn;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &context)) {
    gb = (GBBlockPtr) sdp->data.ptrvalue;
    if (gb != NULL) {
      prev = NULL;
      for (vnp = gb->keywords; vnp; vnp = vnp_next) {       
        vnp_next = vnp->next;
        if (StringCmp((CharPtr)vnp->data.ptrvalue, keyword) == 0) {
          if (prev == NULL) {
            gb->keywords = vnp->next;
          } else {
            prev->next = vnp->next;
          }
          vnp->next = NULL;
          vnp = ValNodeFreeData (vnp);
        } else {
          prev = vnp;
        }
      }
      if (GBBlockIsCompletelyEmpty(gb)) {
        if (sdp->extended) {
          ovn = (ObjValNodePtr) sdp;
          ovn->idx.deleteme = TRUE;
        }
      }
    }
  }
}


NLM_EXTERN Boolean HasKeywordForStructuredCommentName (BioseqPtr bsp, UserObjectPtr uop)
{
  CharPtr keyword = NULL;
  Boolean    has_keyword = FALSE;

  if (bsp == NULL || uop == NULL || (keyword = KeywordForStructuredCommentName(uop)) == NULL) {
    return FALSE;
  }

  has_keyword = HasKeyword(bsp, keyword);
  
  keyword = MemFree (keyword);
  return has_keyword;
}


static void AddStructuredCommentKeywordsCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  UserObjectPtr     uop;
  CharPtr           keyword;

  if (bsp == NULL) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    keyword = KeywordForStructuredCommentName(uop);
    if (keyword != NULL) {
      if (IsStructuredCommentValid (uop, NULL, NULL) == eFieldValid_Valid) {
        AddKeywordToBioseq(bsp, keyword);
      }
    }
    keyword = MemFree (keyword);
  }
}


NLM_EXTERN void AddStructuredCommentKeywords (Uint2 entityID)
{
  SeqEntryPtr sep;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) {
    return;
  }

  VisitBioseqsInSep (sep, NULL, AddStructuredCommentKeywordsCallback);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);    
}


static void RemoveStructuredCommentKeywordsCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  UserObjectPtr     uop;
  CharPtr           keyword;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    keyword = KeywordForStructuredCommentName(uop);
    if (keyword != NULL) {
      if (IsStructuredCommentValid (uop, NULL, NULL) != eFieldValid_Valid) {
        RemoveKeywordFromBioseq(bsp, keyword);
      }
    }
    keyword = MemFree (keyword);
  }
}


NLM_EXTERN void RemoveStructuredCommentKeywords (Uint2 entityID)
{
  SeqEntryPtr sep;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) {
    return;
  }

  VisitBioseqsInSep (sep, NULL, RemoveStructuredCommentKeywordsCallback);
  DeleteMarkedObjects (entityID, 0, NULL);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);    
}


static Boolean StartsWith(CharPtr str, CharPtr start)
{
  Int4 str_len, start_len;

  str_len = StringLen (str);
  start_len = StringLen (start);

  if (str_len < start_len || StringNICmp(str, start, start_len) != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static Boolean EndsWith(CharPtr str, CharPtr end)
{
  Int4 str_len, end_len;

  str_len = StringLen (str);
  end_len = StringLen (end);

  if (str_len < end_len || StringICmp(str + str_len - end_len, end) != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static void TrimPrimerSeqJunkFromString (CharPtr str)
{
  Int4 len, start_len = 0, end_len = 0;
  CharPtr src, dst;
  
  if (StringHasNoText (str)) {
    return;
  }
  len = StringLen (str);

  if (StartsWith (str, "5'-") || StartsWith (str, "5`-")) {
    start_len = 3;
  } else if (StartsWith (str, "5-") || StartsWith (str, "5'") || StartsWith (str, "5`")) {
    start_len = 2;
  } else if (StartsWith (str, "-")) {
    start_len = 1;
  }

  if (EndsWith (str, "-3'") || EndsWith (str, "-3`")) {
    end_len = 3;
  } else if (EndsWith (str, "-3") || EndsWith(str, "3'") || EndsWith(str, "3`")) {
    end_len = 2;
  } else if (EndsWith (str, "-")) {
    end_len = 1;
  }
  
  if (end_len > 0 || start_len > 0) {
    src = str + start_len;
    dst = str;
    len -= (end_len + start_len);

    while (len > 0) {
      *dst = *src;
      src++;
      dst++;
      len--;
    }
    *dst = 0;
  }

}


static Boolean TrimJunkFromPrimer (PCRPrimerPtr pp, FILE *log_fp)
{
  CharPtr orig = NULL;
  Boolean rval = FALSE;

  if (pp == NULL || StringHasNoText (pp->seq)) {
    return FALSE;
  }
  if (log_fp != NULL) {
    orig = StringSave (pp->seq);
  }
  TrimPrimerSeqJunkFromString (pp->seq);
  if (log_fp != NULL && StringCmp (orig, pp->seq) != 0) {
    fprintf (log_fp, "Changed primer seq from %s to %s\n", orig, pp->seq);
    rval = TRUE;
  }
  orig = MemFree (orig);
  return rval;
}


static Boolean TrimPrimerSeqJunkOnBioSource (BioSourcePtr biop, FILE *log_fp)
{
  PCRReactionSetPtr ps;
  PCRPrimerPtr      pp;
  Boolean           rval = FALSE;

  if (biop == NULL) {
    return FALSE;
  }

  for (ps = biop->pcr_primers; ps != NULL; ps = ps->next) {
    for (pp = ps->forward; pp != NULL; pp = pp->next) {
      rval |= TrimJunkFromPrimer(pp, log_fp);
    }
    for (pp = ps->reverse; pp != NULL; pp = pp->next) {
      rval |= TrimJunkFromPrimer(pp, log_fp);
    }
  }

  return rval;
}


static void TrimPrimerSeqJunkDescrCallback (SeqDescrPtr sdp, Pointer data)
{
  LogInfoPtr lip;


  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    lip = (LogInfoPtr) data;
    if (TrimPrimerSeqJunkOnBioSource (sdp->data.ptrvalue, lip == NULL ? NULL : lip->fp) && lip != NULL) {
      lip->data_in_log = TRUE;
    }
  }
}


static void TrimPrimerSeqJunkFeatCallback (SeqFeatPtr sfp, Pointer data)
{
  LogInfoPtr lip;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
    lip = (LogInfoPtr) data;
    if (TrimPrimerSeqJunkOnBioSource (sfp->data.value.ptrvalue, lip == NULL ? NULL : lip->fp) && lip != NULL) {
      lip->data_in_log = TRUE;
    }
  }
}


NLM_EXTERN Boolean TrimPrimerSeqJunkInSeqEntry (SeqEntryPtr sep, FILE *log_fp)
{
  LogInfoData lid;

  MemSet (&lid, 0, sizeof (LogInfoData));
  lid.fp = log_fp;
  VisitDescriptorsInSep (sep, &lid, TrimPrimerSeqJunkDescrCallback);
  VisitFeaturesInSep (sep, &lid, TrimPrimerSeqJunkFeatCallback);
  return lid.data_in_log;
}


static Boolean IsUSA (CharPtr country)
{
  if (StringICmp (country, "USA") == 0
      || StringICmp (country, "United States of America") == 0
      || StringICmp (country, "United States") == 0
      || StringICmp (country, "U.S.A.") == 0
      || StringICmp (country, "U S A") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void FixStateAbbreviationsInCitSub (CitSubPtr csp, LogInfoPtr lip)
{
  if (csp != NULL && csp->authors != NULL 
      && csp->authors->affil != NULL
      && IsUSA(csp->authors->affil->country)) {
    if (StringCmp (csp->authors->affil->country, "USA") != 0) {
      if (lip != NULL) {
        if (lip->fp != NULL) {
          fprintf (lip->fp, "Changed %s to USA\n", csp->authors->affil->country);
        }
        lip->data_in_log = TRUE;
      }
      csp->authors->affil->country = MemFree (csp->authors->affil->country);
      csp->authors->affil->country = StringSave ("USA");
    }
    FixStateAbbreviationsInAffil (csp->authors->affil, NULL);
  }
}


static void AbbreviateCitSubAffilStatesCallback (PubdescPtr pdp, Pointer data)
{
  ValNodePtr vnp;
  CitSubPtr  csp;
  LogInfoPtr lip;

  if (pdp == NULL) return;
  lip = (LogInfoPtr)data;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Sub) {
      csp = (CitSubPtr) vnp->data.ptrvalue;  
      FixStateAbbreviationsInCitSub (csp, lip);
    }
  }  
}


NLM_EXTERN Boolean FixUsaAndStateAbbreviations (Uint2 entityID, FILE *log_fp)
{
  SeqEntryPtr sep;
  LogInfoData lid;
  SeqSubmitPtr ssp;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL)
    return FALSE;

  MemSet (&lid, 0, sizeof (LogInfoData));
  lid.fp = log_fp;
  VisitPubdescsInSep (sep, &lid, AbbreviateCitSubAffilStatesCallback);

  ssp = FindSeqSubmitForSeqEntry (sep);
  if (ssp != NULL && ssp->sub != NULL && ssp->sub->cit != NULL) {
    FixStateAbbreviationsInCitSub (ssp->sub->cit, &lid);
  }
  return lid.data_in_log;
}


static ValNodePtr FindExonForInterval (BioseqPtr bsp, SeqLocPtr slp, Boolean match_from_exactly, Boolean match_to_exactly)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        sfp;
  ValNodePtr        list = NULL;
  Int4              from, to, feat_from, feat_to;
  Uint1             strand;
  SeqPntPtr         spp;
  SeqIntPtr         sint;

  if (slp == NULL) {
    return NULL;
  } else if (slp->choice == SEQLOC_PNT) {
    spp = (SeqPntPtr) slp->data.ptrvalue;
    from = spp->point;
    to = spp->point;
    strand = spp->strand;
  } else if (slp->choice == SEQLOC_INT) {
    sint = (SeqIntPtr) slp->data.ptrvalue;
    from = sint->from;
    to = sint->to;
    strand = sint->strand;
  } else {
    return NULL;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_exon, &context);
       sfp != NULL && context.left <= to;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_exon, &context)) 
  {
    /* note - have to use location values, rather than context.left and context.right,
     * because exon may already have been altered for another mRNA/CDS
     */
    if (sfp->location == NULL) {
      /* no location */
      continue;
    } else if (sfp->location->choice == SEQLOC_PNT) {
      spp = (SeqPntPtr) sfp->location->data.ptrvalue;
      feat_from = spp->point;
      feat_to = spp->point;
    } else if (sfp->location->choice == SEQLOC_INT) {
      sint = (SeqIntPtr) sfp->location->data.ptrvalue;
      feat_from = sint->from;
      feat_to = sint->to;
    } else {
      /* not handling other types of locations */
      continue;
    }
    if (context.numivals != 1) {
      /* not going to match multi-interval exons */
    } else if (match_from_exactly && feat_from != from) {
      /* no match on from */
    } else if (!match_from_exactly && (feat_from < from || feat_from > to)) {
      /* less restrictive match fails for from */
    } else if (match_to_exactly && feat_to != to) {
      /* no match on to */
    } else if (!match_to_exactly && (feat_to > to || feat_to < from)) {
      /* less restrictive match fails for to */
    } else if ((strand == Seq_strand_minus && context.strand != Seq_strand_minus)
      || (strand != Seq_strand_minus && context.strand == Seq_strand_minus)) {
      /* strand match fails */
    } else {
      ValNodeAddPointer (&list, OBJ_SEQFEAT, sfp);
    }
  }
  return list;
}


static ValNodePtr SaveOrigExonPositions (ValNodePtr exon_list)
{
  ValNodePtr vnp;
  SeqFeatPtr exon;
  CharPtr    orig_loc;
  ValNodePtr loc_list = NULL;

  for (vnp = exon_list; vnp != NULL; vnp = vnp->next) 
  {
    exon = (SeqFeatPtr) vnp->data.ptrvalue;
    orig_loc = SeqLocPrintUseBestID (exon->location);
    ValNodeAddPointer (&loc_list, 0, orig_loc);
  }
  return loc_list;
}


static void FixExonsForInterval (ValNodePtr list, Int4 from_diff, Int4 to_diff)
{
  ValNodePtr vnp;
  SeqFeatPtr exon;
  SeqPntPtr  spp;
  SeqIntPtr  sint;

  if (list == NULL) {
    return;
  }
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    exon = vnp->data.ptrvalue;
    if (exon != NULL && exon->location != NULL) {
      if (exon->location->choice == SEQLOC_PNT) {
        spp = (SeqPntPtr) exon->location->data.ptrvalue;
        sint = SeqIntNew ();
        sint->id = spp->id;
        spp->id = NULL;
        sint->strand = spp->strand;
        sint->to = spp->point;
        sint->from = spp->point;
        spp = SeqPntFree (spp);
        exon->location->data.ptrvalue = sint;
      }
      sint = (SeqIntPtr) exon->location->data.ptrvalue;
      sint->from += from_diff;
      sint->to += to_diff;
    }
  }
}

typedef struct exonloclist {
  ValNodePtr feature_list;
  ValNodePtr orig_loc_list;
} ExonLocListData, PNTR ExonLocListPtr;


static ExonLocListPtr ExonLocListNew (BioseqPtr bsp, SeqLocPtr slp, Boolean match_from_exactly, Boolean match_to_exactly)
{
  ExonLocListPtr el = (ExonLocListPtr) MemNew (sizeof (ExonLocListData));
  el->feature_list = FindExonForInterval(bsp, slp, match_from_exactly, match_to_exactly);
  if (el->feature_list == NULL) {
    el = MemFree (el);
  } else {
    el->orig_loc_list = SaveOrigExonPositions(el->feature_list);
  }
  return el;
}


static ExonLocListPtr ExonLocListFree (ExonLocListPtr el)
{
  if (el != NULL) {
    el->feature_list = ValNodeFree (el->feature_list);
    el->orig_loc_list = ValNodeFreeData (el->orig_loc_list);
    el = MemFree (el);
  }
  return el;
}


static void ReportExonLocationChanges (ExonLocListPtr el, LogInfoPtr lip)
{
  ValNodePtr exon_v, orig;
  SeqFeatPtr exon;
  CharPtr    new_loc;

  if (lip == NULL || el == NULL) {
    return;
  }
  for (exon_v = el->feature_list, orig = el->orig_loc_list; exon_v != NULL && orig != NULL; exon_v = exon_v->next, orig = orig->next) {
    exon = (SeqFeatPtr) exon_v->data.ptrvalue;
    new_loc = SeqLocPrintUseBestID (exon->location);
    if (StringCmp (orig->data.ptrvalue, new_loc) != 0) {
      if (lip->fp != NULL) {
        fprintf (lip->fp, "Adjusted location for splice consensus: %s became %s\n", (char*) orig->data.ptrvalue, new_loc);
      }
      lip->data_in_log = TRUE;
    }
    new_loc = MemFree (new_loc);
  }
}


static void AdjustForConsensusSpliceCallback (SeqFeatPtr sfp, Pointer data)
{
  SeqLocPtr slp, slp_last = NULL;
  SeqIdPtr  sip;
  BioseqPtr bsp;
  Boolean   partial5, partial3, partial5_last, partial3_last, first = TRUE;
  Uint1     strand = Seq_strand_unknown, this_strand;
  Int4      prev_pos = -1, this_pos, exon_len, exon_len_last = -1;
  CharPtr   buf;
  Int4      len, start, stop, diff;
  Boolean   match;
  SeqIntPtr sint;
  SeqPntPtr spp;
  ExonLocListPtr last_exon_list = NULL, this_exon_list = NULL;
  /* variables used for logging change */
  CharPtr orig_loc = NULL, new_loc;
  Boolean changed = FALSE;
  LogInfoPtr lip;

  if (sfp == NULL 
      || (sfp->data.choice != SEQFEAT_CDREGION && sfp->idx.subtype != FEATDEF_mRNA)
      || sfp->location == NULL 
      || (sfp->location->choice != SEQLOC_MIX && sfp->location->choice != SEQLOC_PACKED_INT)
      || (sip = SeqLocId (sfp->location)) == NULL
      || (bsp = BioseqLockById (sip)) == NULL) {
    return;
  }

  /* we're not going to handle mixed-strand exons */
  for (slp = sfp->location->data.ptrvalue; slp != NULL && strand != Seq_strand_other; slp = slp->next) {
    this_strand = SeqLocStrand (slp);
    if (this_strand == Seq_strand_minus) {
      if (first) {
        strand = Seq_strand_minus;
      } else if (strand != Seq_strand_minus) {
        strand = Seq_strand_other;
      }
    } else {
      if (strand == Seq_strand_minus) {
        strand = Seq_strand_other;
      }
    }
    first = FALSE;
  }

  if (strand == Seq_strand_other) {
    BioseqUnlock (bsp);
    return;
  }

  if ((lip = (LogInfoPtr)data) != NULL && lip->fp != NULL) {
    orig_loc = SeqLocPrintUseBestID (sfp->location);
  }

  first = TRUE;
  for (slp = sfp->location->data.ptrvalue; slp != NULL; slp = slp->next) {
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    exon_len = SeqLocLen (slp);
    /* record underlying exon features */
    this_exon_list = ExonLocListNew (bsp, slp, TRUE, TRUE);

    if (!first && !partial5 && !partial3_last
        && (slp_last->choice == SEQLOC_INT || slp_last->choice == SEQLOC_PNT)
        && (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_PNT)) { 

      /* check for donor and acceptor pair */
      /* maximum search space is beginning of previous exon to end of current exon */
      exon_len_last = SeqLocLen (slp_last);
      if (strand == Seq_strand_minus) {
        this_pos = SeqLocStart (slp);
      } else {
        this_pos = SeqLocStop (slp);
      }
      start = MIN (this_pos, prev_pos);
      stop = MAX (this_pos, prev_pos);
      len = stop - start + 1;
      buf = (CharPtr) MemNew (sizeof (Char) * (len + 1));
      SeqPortStreamInt (bsp, start, stop, strand, EXPAND_GAPS_TO_DASHES, (Pointer) buf, NULL);
      if (buf[exon_len_last] == 'G' && (buf[exon_len_last + 1] == 'T' || buf[exon_len_last + 1] == 'C')
          && buf[len - exon_len - 1] == 'G' && buf[len - exon_len - 2] == 'A') {
        /* already have donor acceptor pair */
      } else {
        match = FALSE;
        /* search forward */
        if ((slp_last->choice == SEQLOC_INT || slp_last->choice == SEQLOC_PNT)
            && slp->choice == SEQLOC_INT) {
          diff = 1;
          while (diff < exon_len && !match) {
            if (buf[exon_len_last + diff - 1] != buf[len - exon_len + diff - 1]) {
              diff = exon_len;
            } else if (buf[exon_len_last + diff] == 'G' && buf[exon_len_last + 1 + diff] == 'T'
                && buf[len - exon_len - 1 + diff] == 'G' && buf[len - exon_len - 2 + diff] == 'A') {
              match = TRUE;
            } else {
              diff++;
            }
          }
          if (match) {
            if (slp_last->choice == SEQLOC_PNT) {
              spp = (SeqPntPtr) slp_last->data.ptrvalue;
              sint = SeqIntNew ();
              sint->id = spp->id;
              spp->id = NULL;
              sint->strand = spp->strand;
              sint->to = spp->point;
              sint->from = spp->point;
              spp = SeqPntFree (spp);
              slp_last->data.ptrvalue = sint;
            }
            sint = (SeqIntPtr) slp_last->data.ptrvalue;
            if (sint->strand == Seq_strand_minus) {
              sint->from -= diff;
            } else {
              sint->to += diff;
            }
            sint = (SeqIntPtr) slp->data.ptrvalue;
            if (sint->strand == Seq_strand_minus) {
              sint->to -= diff;
            } else {
              sint->from += diff;
            }
            changed = TRUE;
          }
        }
        /* search backward */
        if (!match && slp_last->choice == SEQLOC_INT
            && (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_PNT)) {
          diff = 1;
          while (diff < exon_len_last && !match) {
            if (buf[exon_len_last - diff] != buf[len - exon_len - diff]) {
              diff = exon_len_last;
            } else if (buf[exon_len_last - diff] == 'G' && buf[exon_len_last + 1 - diff] == 'T'
                && buf[len - exon_len - 1 - diff] == 'G' && buf[len - exon_len - 2 - diff] == 'A') {
              match = TRUE;
            } else {
              diff++;
            }
          }
          if (match) {
            if (slp->choice == SEQLOC_PNT) {
              spp = (SeqPntPtr) slp->data.ptrvalue;
              sint = SeqIntNew ();
              sint->id = spp->id;
              spp->id = NULL;
              sint->strand = spp->strand;
              sint->to = spp->point;
              sint->from = spp->point;
              spp = SeqPntFree (spp);
              slp->data.ptrvalue = sint;
            }
            sint = (SeqIntPtr) slp->data.ptrvalue;
            if (sint->strand == Seq_strand_minus) {
              sint->to += diff;
              if (this_exon_list != NULL) {
                FixExonsForInterval (this_exon_list->feature_list, 0, diff);
              }
            } else {
              sint->from -= diff;
              if (this_exon_list != NULL) {
                FixExonsForInterval (this_exon_list->feature_list, -diff, 0);
              }
            }
            sint = (SeqIntPtr) slp_last->data.ptrvalue;
            if (sint->strand == Seq_strand_minus) {
              sint->from += diff;
              if (last_exon_list != NULL) {
                FixExonsForInterval (last_exon_list->feature_list, diff, 0);
              }
            } else {
              sint->to -= diff;
              if (last_exon_list != NULL) {
                FixExonsForInterval (last_exon_list->feature_list, 0, -diff);
              }
            }
            changed = TRUE;
          }
        }
      }

      buf = MemFree (buf);
    }

    if (strand == Seq_strand_minus) {
      prev_pos = SeqLocStop (slp);
    } else {
      prev_pos = SeqLocStart (slp);
    }

    partial5_last = partial5;
    partial3_last = partial3;
    slp_last = slp;
    ReportExonLocationChanges (last_exon_list, lip);
    last_exon_list = ExonLocListFree (last_exon_list);
    last_exon_list = this_exon_list;
    first = FALSE;
  }  

  ReportExonLocationChanges (last_exon_list, lip);
  last_exon_list = ExonLocListFree (last_exon_list);

  BioseqUnlock (bsp);

  if (changed) {
    if (lip->fp != NULL) {
      new_loc = SeqLocPrintUseBestID (sfp->location);
      fprintf (lip->fp, "Adjusted location for splice consensus: %s became %s\n", orig_loc, new_loc);
      new_loc = MemFree (new_loc);
    }
    lip->data_in_log = TRUE;
  }
  orig_loc = MemFree (orig_loc);
}


static void AdjustSeqEntryForConsensusSpliceBioseqCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext dcontext;
  BioSourcePtr biop;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || (biop = (BioSourcePtr)sdp->data.ptrvalue) == NULL
      || (biop->genome != GENOME_genomic && biop->genome != GENOME_unknown)
      || (biop->org != NULL && biop->org->orgname != NULL && StringISearch (biop->org->orgname->lineage, "viruses") != NULL)
      || !HasTaxonomyID(biop)) 
  {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext)) 
  {
    AdjustForConsensusSpliceCallback (sfp, data);
  }
}


NLM_EXTERN Boolean AdjustSeqEntryForConsensusSpliceEx (SeqEntryPtr sep, FILE *log_fp)
{
  LogInfoData lid;

  if (sep == NULL) {
    return FALSE;
  }
  MemSet (&lid, 0, sizeof (LogInfoData));
  lid.fp = log_fp;

  VisitBioseqsInSep (sep, &lid, AdjustSeqEntryForConsensusSpliceBioseqCallback);
  return lid.data_in_log;
}

NLM_EXTERN void AdjustSeqEntryForConsensusSplice (SeqEntryPtr sep)
{
  AdjustSeqEntryForConsensusSpliceEx (sep, NULL);
}


NLM_EXTERN CharPtr ValNodeSeqIdName (ValNodePtr vnp)
{
  Char buf[100];

  if (vnp == NULL || vnp->data.ptrvalue == NULL)
  {
    return NULL;
  }
  else
  {
    SeqIdWrite (vnp->data.ptrvalue, buf, PRINTID_FASTA_SHORT, sizeof (buf) - 1);
    return StringSave (buf);
  }
}


NLM_EXTERN void ValNodeSeqIdFree (ValNodePtr vnp)
{
  if (vnp != NULL && vnp->data.ptrvalue != NULL)
  {
    vnp->data.ptrvalue = SeqIdFree (vnp->data.ptrvalue);
  }
}


NLM_EXTERN ValNodePtr ValNodeSeqIdCopy (ValNodePtr vnp)
{
  ValNodePtr vnp_copy = NULL;
  if (vnp != NULL)
  {
    ValNodeAddPointer (&vnp_copy, vnp->choice, SeqIdDup (vnp->data.ptrvalue));
  }
  return vnp_copy;
}

NLM_EXTERN Boolean ValNodeSeqIdMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  if (SeqIdComp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) == SIC_YES) 
  {
    return TRUE;
  } 
  else 
  {
    return FALSE;
  }
}


NLM_EXTERN ValNodePtr ValNodeSeqIdListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = SeqIdFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


NLM_EXTERN ValNodePtr ValNodeSeqIdListCopy (ValNodePtr list)
{
  ValNodePtr vnp, list_copy = NULL, list_prev = NULL;

  while (list != NULL) {
    vnp = ValNodeNew (list_prev);
    vnp->data.ptrvalue = SeqIdDup (list->data.ptrvalue);
    if (list_copy == NULL) {
      list_copy = vnp;
    }
    list_prev = vnp;
    list = list->next;
  }
  return list_copy;
}


NLM_EXTERN ValNodePtr SeqIdListToValNodeSeqIdList (SeqIdPtr sip_list)
{
  SeqIdPtr sip;
  ValNodePtr list = NULL, vnp_p = NULL, vnp;

  for (sip = sip_list; sip != NULL; sip = sip->next) {
    vnp = ValNodeNew (vnp_p);
    if (vnp_p == NULL) {
      list = vnp;
    }
    vnp->data.ptrvalue = SeqIdDup (sip);
    vnp_p = vnp;
  }
  return list;
}


NLM_EXTERN SeqIdPtr ValNodeSeqIdListToSeqIdList (ValNodePtr vnp_list)
{
  ValNodePtr vnp;
  SeqIdPtr sip_list = NULL, sip_prev = NULL, sip;

  for (vnp = vnp_list; vnp != NULL; vnp = vnp->next) {
    sip = SeqIdDup (vnp->data.ptrvalue);
    if (sip_prev == NULL) {
      sip_list = sip;
    } else {
      sip_prev->next = sip;
    }
    sip_prev = sip;
  }
  return sip_list;
}


