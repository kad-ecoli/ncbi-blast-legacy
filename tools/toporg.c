static char const rcsid[] = "$Id: toporg.c,v 6.189 2011/10/28 15:15:23 kans Exp $";

#include <stdio.h>
#include <ncbi.h>
#include <sequtil.h>
#include <toasn3.h>
#include <toporg.h>
#include <tfuns.h>
#include <utilpub.h>

/* includes for new cleanup functions from Sequin */
#include <sqnutils.h>
#include <gather.h>
#include <explore.h>
#include <subutil.h>
#include <tofasta.h>
#include <objfdef.h>
#include <valid.h>

static ValNodePtr GetDescrNoTitles (ValNodePtr PNTR descr);

SeqDescrPtr remove_descr PROTO((SeqDescrPtr head, SeqDescrPtr x));
/****************************************************************************
 *   move org-ref, modif, mol_type, date, title and pubs
 *   to seg-set level in segmented set, if not there
 *   move org-ref, modif, date, title to nuc-prot level, if not there
 *   notice that mol_type and pub are not moved to nucprot level
 *****************************************************************************/
void toporg(SeqEntryPtr sep)
{

  SeqEntryExplore(sep, NULL, ChkSegset);
  SeqEntryExplore(sep, NULL, ChkNucProt);

  return;
}

/****************************************************************************
 *     check for backbone entry and expand (to whole) OrgRef feature
 *    if it's not whole and the only one OrgRef in the entry
 *     whole fetures would be converted to descr later in FindOrg function
 *
 *     check if org-ref, modif, mol_type, date and title
 *     are the same for all segments
 *     move them to seg-set level in segmented set, if not already there
 *****************************************************************************/
void ChkSegset (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

  BioseqSetPtr  bssp, tmp;
  BioseqPtr     bsp;
  SeqEntryPtr   segsep, parts = NULL;
  ValNodePtr    vnp = NULL, set_vnp = NULL, upd_date_vnp = NULL;
  ValNodePtr    org, modif, mol, date, v /*, title */;
  SeqAnnotPtr   sap = NULL;
  SeqFeatPtr    tmp_sfp, sfp0, sfp;
  SeqIdPtr      sidp;
  Pointer       pnt;
  SeqLocPtr     slp;
  SeqIntPtr     sip;
  OrgRefPtr     orp;
  Boolean       is_org = FALSE, is_modif = FALSE, is_title = FALSE;
  Boolean       is_date = FALSE, is_mol = FALSE;
  Boolean       is_na = FALSE, is_bb = FALSE, whole = FALSE;
  Int2          count = 0;
  Int4          len;
  static Char   msg[51];

  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    vnp = bsp->descr;
    sap = bsp->annot;
    len = bsp->length;
    if (bsp->mol != Seq_mol_aa) {
      is_na = TRUE;
    }
    if (bsp == NULL || !is_na) {
      return;
    }
    for (sidp = bsp->id; sidp != NULL; sidp = sidp->next) {
      pnt = sidp->data.ptrvalue;
      switch (sidp->choice)
      {
        case SEQID_LOCAL:      /* local */
        case SEQID_GIIM:       /* giimid */
        case SEQID_PATENT:     /* patent seq id */
        case SEQID_GENERAL:    /* general */
        case SEQID_GI:         /* gi */
        case SEQID_PDB:
        case SEQID_GENBANK:    /* genbank */
        case SEQID_EMBL:       /* embl */
        case SEQID_PIR:        /* pir   */
        case SEQID_SWISSPROT:  /* swissprot */
        case SEQID_OTHER:      /* other */
        case SEQID_DDBJ:
        case SEQID_PRF:
        case SEQID_TPG:
        case SEQID_TPE:
        case SEQID_TPD:
          continue;
        case SEQID_GIBBSQ:     /* gibbseq */
        case SEQID_GIBBMT:     /* gibbmt */
          is_bb = TRUE;
          break;
        default:
          continue;
      }
    }
    if (!is_bb) {
      return;
    }
    if (sap != NULL && sap->type == 1)   /* feature table */
    {
      tmp_sfp = (SeqFeatPtr) (sap->data);
      sfp0 = SeqFeatExtractList(&(tmp_sfp), SEQFEAT_ORG);
      for (sfp = sfp0; sfp != NULL; sfp = sfp->next) {
        orp = (OrgRefPtr)(sfp->data.value.ptrvalue);
        count++;
        if ((whole = check_whole(sfp, bsp->length)) == TRUE) {
          break;
        }
      }
      if (!whole && count == 1) {
        StringNCpy(msg, SeqLocPrint(sfp0->location), 50);
        ErrPostEx(SEV_WARNING, 0, 2,
                  "Backbone entry source with bad OrgRef feature: %s", msg);
        slp = sfp0->location;
        sip = slp->data.ptrvalue;
        sip->from = 0;
        sip->to = len-1;
      }
      tmp_sfp = tie_feat(tmp_sfp, sfp0);
      sap->data = tmp_sfp;
      bsp->annot = sap;
    }
    return;
  }

  bssp = (BioseqSetPtr)(sep->data.ptrvalue);
  if (bssp->_class != 2) {    /*  do the rest for segset only */
    return;
  }
  segsep = bssp->seq_set;
  set_vnp = bssp->descr;
  if (segsep->next == NULL) {
    return;
  }
  if (!IS_Bioseq(segsep->next)) {
    tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
    parts = tmp->seq_set;
    if (parts == NULL) {
      return;
    }
  }
  /*    title = SrchSegChoice(parts, Seq_descr_title);*/
  org = SrchSegChoice(parts, Seq_descr_org);
  mol = SrchSegChoice(parts, Seq_descr_mol_type);
  modif = SrchSegChoice(parts, Seq_descr_modif);
  date = SrchSegChoice(parts, Seq_descr_update_date);
  for(v=set_vnp; v != NULL; v = v->next) {
    if (v->choice == Seq_descr_org) {
      is_org = TRUE;
    }
    if (v->choice == Seq_descr_org) {
      is_title = TRUE;
    }
    if (v->choice == Seq_descr_mol_type) {
      if (mol && mol->data.intvalue == v->data.intvalue) {
        is_mol = TRUE;
        ValNodeFree(mol);
      }
    }
    if (v->choice == Seq_descr_modif) {
      is_modif = TRUE;
    }
    if (v->choice == Seq_descr_update_date) {
      is_date = TRUE;
      upd_date_vnp = v;
    }
  }
  /*
   if (!is_title) {
   set_vnp = tie_next(set_vnp, title);
   }
   */
  if (!is_modif) {
    if (set_vnp != NULL) {
      set_vnp = tie_next(set_vnp, modif);
    } else {
      ValNodeLink (&(bssp->descr), modif);
    }
  }
  if (!is_org) {
    if (set_vnp != NULL) {
      set_vnp = tie_next(set_vnp, org);
    } else {
      ValNodeLink (&(bssp->descr), org);
    }
  }
  if (!is_mol) {
    if (set_vnp != NULL) {
      set_vnp = tie_next(set_vnp, mol);
    } else {
      ValNodeLink (&(bssp->descr), mol);
    }
  }
  if (!is_date) {
    if (set_vnp != NULL) {
      set_vnp = tie_next(set_vnp, date);
    } else {
      ValNodeLink (&(bssp->descr), date);
    }
  } else if (upd_date_vnp != NULL && date != NULL) {
    upd_date_vnp->data.ptrvalue = DateFree ((DatePtr) upd_date_vnp->data.ptrvalue);
    upd_date_vnp->data.ptrvalue = DateDup ((DatePtr) date->data.ptrvalue);
  }
  SrchSegSeqMol(parts);

  return;
}

static Int2 PubLabelMatchEx (ValNodePtr vnp1, ValNodePtr vnp2, Boolean checkdates)

{
  AffilPtr     afp1, afp2;
  AuthListPtr  alp1, alp2;
  CitGenPtr    cgp1 = NULL, cgp2 = NULL;
  CitSubPtr    csp1 = NULL, csp2 = NULL;
  Int2         ret;

  if (vnp1 == NULL || vnp2 == NULL) return -1;
  ret = PubLabelMatch (vnp1, vnp2);
  if (ret != 0) return ret;

  while (vnp1 != NULL) {
    if (vnp1->choice == PUB_Sub) {
      csp1 = (CitSubPtr) vnp1->data.ptrvalue;
    } else if (vnp1->choice == PUB_Gen) {
      cgp1 = (CitGenPtr) vnp1->data.ptrvalue;
    }
    vnp1 = vnp1->next;
  }

  while (vnp2 != NULL) {
    if (vnp2->choice == PUB_Sub) {
      csp2 = (CitSubPtr) vnp2->data.ptrvalue;
    } else if (vnp2->choice == PUB_Gen) {
      cgp2 = (CitGenPtr) vnp2->data.ptrvalue;
    }
    vnp2 = vnp2->next;
  }

  if (csp1 != NULL && csp2 != NULL) {
    if (checkdates && csp1->date != NULL && csp2->date != NULL && DateMatch (csp1->date, csp2->date, FALSE) != 0) return -1;
    if (StringICmp (csp1->descr, csp2->descr) != 0) return -1;
    alp1 = csp1->authors;
    alp2 = csp2->authors;
    if (alp1 == NULL || alp2 == NULL) return 0;
    if (AuthListMatch (alp1, alp2, TRUE) != 0) return -1;
    afp1 = alp1->affil;
    afp2 = alp2->affil;
    if (afp1 != NULL && afp2 != NULL) {
      if (! AsnIoMemComp (afp1, afp2, (AsnWriteFunc) AffilAsnWrite)) return -1;
    } else if (afp1 != NULL || afp2 != NULL) {
      return -1;
    }
    return 0;
  }

  if (cgp1 != NULL && cgp2 != NULL) {
    if (checkdates && cgp1->date != NULL && cgp2->date != NULL && DateMatch (cgp1->date, cgp2->date, FALSE) != 0) return -1;
    if (StringICmp (cgp1->cit, cgp2->cit) != 0) return -1;
    alp1 = cgp1->authors;
    alp2 = cgp2->authors;
    if (alp1 == NULL || alp2 == NULL) return 0;
    if (AuthListMatch (alp1, alp2, TRUE) != 0) return -1;
    afp1 = alp1->affil;
    afp2 = alp2->affil;
    if (afp1 != NULL && afp2 != NULL) {
      if (! AsnIoMemComp (afp1, afp2, (AsnWriteFunc) AffilAsnWrite)) return -1;
    } else if (afp1 != NULL || afp2 != NULL) {
      return -1;
    }
    return 0;
  }

  return 0;
}

static void RemovePubFromParts(SeqEntryPtr sep, ValNodePtr pub)
{
  BioseqPtr    b;
  ValNodePtr   v, vnp, next;
  SeqEntryPtr  s;
  PubdescPtr   pdp, p;

  for (vnp = pub; vnp; vnp = vnp->next) {
    pdp = vnp->data.ptrvalue;
    for (s= sep; s; s=s->next) {
      b = (BioseqPtr)(s->data.ptrvalue);
      for (v=b->descr; v; v=next) {
        next = v->next;
        if (v->choice != Seq_descr_pub)
          continue;
        p = v->data.ptrvalue;
        if (PubLabelMatchEx (pdp->pub, p->pub, TRUE) == 0) {
          if (pdp->name != NULL || pdp->fig != NULL
              || pdp->num != NULL || pdp->maploc != NULL
              || pdp->comment != NULL) {
            continue;
          } else {
            b->descr = remove_descr(b->descr, v);
          }
        }
      }
    }
  }
  return;
}
/***************************************************************************
 *    0  match
 *    1 no match
 *    2 type unknown or not implemented
 ***************************************************************************/
static Int2 NumberingMatch(ValNodePtr num1, ValNodePtr num2)
{
  NumContPtr nc1, nc2;
  NumEnumPtr ne1, ne2;
  NumRealPtr nr1, nr2;

  if (num1 == NULL || num2 == NULL)
    return 0;
  if (num1->choice != num2->choice)
    return 1;
  switch(num1->choice)
  {
    case Numbering_cont:
      nc1 = (NumContPtr) num1->data.ptrvalue;
      nc2 = (NumContPtr) num2->data.ptrvalue;
      if (nc1->refnum == nc2->refnum) {
        return 0;
      }
      break;
    case Numbering_enum:
      ne1 = (NumEnumPtr) num1->data.ptrvalue;
      ne2 = (NumEnumPtr) num2->data.ptrvalue;
      if (ne1->num != ne2->num)
        return 1;
      return (StringCmp(ne1->buf, ne2->buf));
    case Numbering_ref_source:
      break;
    case Numbering_ref_align:
      break;
    case Numbering_real:
      nr1 = (NumRealPtr) num1->data.ptrvalue;
      nr2 = (NumRealPtr) num2->data.ptrvalue;
      return (StringCmp(nr1->units, nr2->units));
    default:
      break;
  }
  return 2;
}

static AuthListPtr AlpFromPdp (PubdescPtr pdp, CitSubPtr csp)

{
  AuthListPtr  alp = NULL;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  ValNodePtr   vnp;

  if (csp != NULL) {
    alp = csp->authors;
    if (alp != NULL) return alp;
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          alp = cgp->authors;
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          alp = csp->authors;
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          alp = cap->authors;
        }
        break;
      case PUB_Book :
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          alp = cbp->authors;
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          alp = cpp->authors;
        }
        break;
      default :
        break;
    }

    if (alp != NULL) return alp;
  }

  return NULL;
}

static Boolean PubdescMatch (PubdescPtr p1, PubdescPtr p2)
{
  AuthListPtr alp1, alp2;

  if (p1 == NULL || p2 == NULL)
    return TRUE;
  if (p1->name && p2->name) {
    if (StringCmp(p1->name, p2->name) != 0)
      return FALSE;
  }
  if (p1->fig && p2->fig) {
    if (StringCmp(p1->fig, p2->fig) != 0)
      return FALSE;
  }
  if (p1->maploc && p2->maploc) {
    if (StringCmp(p1->maploc, p2->maploc) != 0)
      return FALSE;
  }
  if (p1->comment && p2->comment) {
    if (StringCmp(p1->comment, p2->comment) != 0)
      return FALSE;
  }
  if (p1->num && p2->num) {
    if (NumberingMatch(p1->num, p2->num) != 0)
      return FALSE;
  }
  /* do full author match */
  alp1 = AlpFromPdp (p1, NULL);
  alp2 = AlpFromPdp (p2, NULL);
  if (alp1 != NULL && alp2 != NULL) {
    if (AuthListMatch (alp1, alp2, TRUE) != 0)
      return FALSE;
  }
  return TRUE;
}

/* return list of pubs that are the same in all segmets */
static ValNodePtr CheckSegsForPub(SeqEntryPtr sep)
{
  BioseqPtr    bsp, b;
  ValNodePtr   vnp, v, list = NULL, vnpnext, new, next;
  PubdescPtr   pdp, p, new_p;
  Boolean      same;
  SeqEntryPtr  s;

  if (sep == NULL) {
    return NULL;
  }
  if (!IS_Bioseq(sep)) {
    return NULL;
  }
  bsp = (BioseqPtr)(sep->data.ptrvalue);
  /* first bioseq from parts */
  for (vnp=bsp->descr; vnp; vnp=vnpnext) {
    vnpnext = vnp->next;
    if (vnp->choice != Seq_descr_pub)
      continue;
    pdp = vnp->data.ptrvalue;
    for (s= sep->next, same = FALSE; s; s=s->next) {
      b = (BioseqPtr)(s->data.ptrvalue);
      /* added to skip virtual spacers */
      if (b->repr == Seq_repr_virtual && b->descr == NULL) continue;
      for (v=b->descr; v; v=next) {
        next = v->next;
        if (v->choice != Seq_descr_pub)
          continue;
        p = v->data.ptrvalue;
        if (PubLabelMatchEx (pdp->pub, p->pub, FALSE) == 0) {
          if (PubdescMatch(pdp, p) == TRUE) {
            same = TRUE;
            break;
          }
        }
      }
      if (v == NULL) {
        same = FALSE;
        break;
      }
    }
    if (same == TRUE) {
      new = SeqDescrNew(NULL);
      new->choice = Seq_descr_pub;
      new_p = AsnIoMemCopy(pdp, (AsnReadFunc) PubdescAsnRead,
                           (AsnWriteFunc) PubdescAsnWrite);
      new->data.ptrvalue = new_p;
      list = tie_next(list, new);
    }
  }
  return list;
}
/* move identical pubs in segmented set to the set level */
void MoveSegmPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

  BioseqSetPtr  bssp, tmp;
  SeqEntryPtr   segsep, parts = NULL;
  ValNodePtr    v, pub, vv, next;
  PubdescPtr    pdp, pdpv;

  if (IS_Bioseq(sep)) {
    return;
  }
  bssp = (BioseqSetPtr)(sep->data.ptrvalue);
  if (bssp->_class != 2) {    /*  do the rest for segset only */
    return;
  }
  segsep = bssp->seq_set;
  if (segsep->next == NULL) {
    return;
  }
  if (!IS_Bioseq(segsep->next)) {
    tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
    parts = tmp->seq_set;
    if (parts == NULL) {
      return;
    }
  }
  /* find identical pubs in all segments */
  pub = CheckSegsForPub(parts);
  if (pub) {
    RemovePubFromParts(parts, pub);
  }
  /* check if pub is already on the set descr */
  for(v=bssp->descr; v != NULL; v = v->next) {
    if (v->choice != Seq_descr_pub)
      continue;
    for (vv = pub; vv; vv = next) {
      next = vv->next;
      pdp = vv->data.ptrvalue;
      pdpv = v->data.ptrvalue;
      if (PubLabelMatchEx (pdp->pub, pdpv->pub, FALSE) == 0) {
        PubdescFree(pdp);
        pub = remove_node(pub, vv);
      }
    }
  }

  bssp->descr = tie_next(bssp->descr, pub);

  return;
}

static Boolean AllPartsHaveTitles (BioseqSetPtr bssp)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  CharPtr      str;

  if (bssp == NULL) return FALSE;
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (sep->choice != 1) return FALSE;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return FALSE;
    str = BioseqGetTitle (bsp);
    if (StringHasNoText (str)) return FALSE;
  }

  return TRUE;
}

void ChkNucProt (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqSetPtr  bssp, tmp, bssp2, bssp3;
  BioseqPtr     bsp = NULL, bsp2, bsp3, bsp4;
  SeqEntryPtr   seqsep, sep2, sep3, sep4;
  ValNodePtr    descr = NULL, head = NULL, vnp, next;
  Int2Ptr       ip;
  Boolean       is_org = FALSE, is_modif = FALSE, is_title = FALSE;
  Boolean       is_date = FALSE, is_nuc /* , is_pub = FALSE */;
  CharPtr       npstitle = NULL, seqtitle = NULL;
  ValNodePtr    PNTR prev;
  Char          ch;

  if (IS_Bioseq(sep)) {
    return;
  }
  bssp = (BioseqSetPtr)(sep->data.ptrvalue);
  if (bssp->_class != BioseqseqSet_class_nuc_prot) {    /*  do the rest for nuc-prot only */
    return;
  }
  seqsep = bssp->seq_set;
  if (seqsep == NULL || seqsep->data.ptrvalue == NULL) return;

  /* if protein first, extract nucleotide, put first */
  if (seqsep->choice == 1) {
    bsp3 = seqsep->data.ptrvalue;
    if (bsp3 != NULL && ISA_aa (bsp3->mol)) {
      prev = (ValNodePtr PNTR) &(bssp->seq_set);
      vnp = bssp->seq_set;
      while (vnp != NULL) {
        next = vnp->next;
        is_nuc = FALSE;
        if (IS_Bioseq (vnp)) {
          bsp4 = (BioseqPtr) vnp->data.ptrvalue;
          if (bsp4 != NULL && ISA_na (bsp4->mol)) {
            is_nuc = TRUE;
          }
        }
        if (is_nuc) {
          *prev = vnp->next;
          ValNodeLink (&head, vnp);
          vnp->next = NULL;
        } else {
          prev = (ValNodePtr PNTR) &(vnp->next);
        }
        vnp = next;
      }
      if (head != NULL) {
        vnp = bssp->seq_set;
        bssp->seq_set = head;
        ValNodeLink (&(bssp->seq_set), vnp);
      }
    }

    seqsep = bssp->seq_set;
    if (seqsep == NULL || seqsep->data.ptrvalue == NULL) return;
  }

  if (seqsep->choice == 1) {
    bsp = seqsep->data.ptrvalue;
    descr = bsp->descr;
  } else if (seqsep->choice == 2) {
    tmp = seqsep->data.ptrvalue;
    descr = tmp->descr;
  }
  ip = (Int2Ptr) data;
  if (bssp->descr == NULL) {
    bssp->descr = GetDescrNoTitles(&descr);
  } else {
    for (vnp = bssp->descr; vnp!= NULL; vnp= vnp->next) {
      if (vnp->choice == Seq_descr_title) {
        is_title = TRUE;
        npstitle = (CharPtr) vnp->data.ptrvalue;
      }
      if (vnp->choice == Seq_descr_org) {
        is_org = TRUE;
      }
      if (vnp->choice == Seq_descr_modif) {
        is_modif = TRUE;
      }
      if (vnp->choice == Seq_descr_update_date) {
        is_date = TRUE;
      }
      /*    if (vnp->choice == Seq_descr_pub) {
       is_pub = TRUE;
       }
       */
    }
    /* look for old style nps title, remove if based on nucleotide title, also remove exact duplicate */
    if (npstitle != NULL) {
      seqtitle = SeqEntryGetTitle (seqsep);
      if (seqtitle != NULL) {
        ch = *npstitle;
        while (ch != '\0' && ch == *seqtitle) {
          npstitle++;
          ch = *npstitle;
          seqtitle++;
        }
        if (ch == '\0' || StringCmp (npstitle, ", and translated products") == 0) {
          vnp = ValNodeExtractList (&(bssp->descr), Seq_descr_title);
          ValNodeFreeData (vnp);
          if (ip != NULL) {
            *ip = 1;
          }
        } else {
          /* now removing any unrelated title */
          vnp = ValNodeExtractList (&(bssp->descr), Seq_descr_title);
          ValNodeFreeData (vnp);
          if (ip != NULL) {
            *ip = 2;
          }
        }
      } else if (bsp != NULL) {
        /* if no nucleotide title, move nps title to it */
        vnp = ValNodeExtractList (&(bssp->descr), Seq_descr_title);
        bsp->descr = ValNodeLink (&(bsp->descr), vnp);
        descr = bsp->descr;
        if (ip != NULL) {
          *ip = 3;
        }
      } else if (seqsep->choice == 2) {
        if (ip != NULL) {
          *ip = 4;
        }
        /* get segmented sequence in segset */
        sep2 = bssp->seq_set;
        if (sep2 != NULL && sep2->choice == 2 && SeqEntryGetTitle (sep2) == NULL) {
          bssp2 = sep2->data.ptrvalue;
          if (bssp2 != NULL && bssp2->_class == BioseqseqSet_class_segset) {
            sep3 = bssp2->seq_set;
            if (sep3 != NULL && sep3->choice == 1) {
              bsp2 = sep3->data.ptrvalue;
              if (bsp2 != NULL && BioseqGetTitle (bsp2) == NULL) {
                sep4 = sep3->next;
                if (sep4 != NULL && sep4->choice == 2) {
                  bssp3 = sep4->data.ptrvalue;
                  if (bssp3 != NULL && bssp3->_class == BioseqseqSet_class_parts) {
                    if (AllPartsHaveTitles (bssp3)) {
                      /* if no segmented nucleotide bioseq title, move nps title to it */
                      vnp = ValNodeExtractList (&(bssp->descr), Seq_descr_title);
                      bsp2->descr = ValNodeLink (&(bsp2->descr), vnp);
                      if (ip != NULL) {
                        *ip = 5;
                      }
                    } else {
                      if (ip != NULL) {
                        *ip = 6;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        if (ip != NULL) {
          *ip = 7;
        }
      }
    }
    /*
     if (!is_title) {
     vnp = ValNodeExtractList(&descr, Seq_descr_title);
     bssp->descr = ValNodeLink(&(bssp->descr), vnp);
     }
     */
    if (!is_modif && check_GIBB(descr)) {
      vnp = ValNodeExtractList(&descr, Seq_descr_modif);
      if (vnp != NULL) {
        bssp->descr = ValNodeLink(&(bssp->descr), vnp);
      }
    }
    if (!is_org) {
      vnp = ValNodeExtractList(&descr, Seq_descr_org);
      if (vnp != NULL) {
        bssp->descr = ValNodeLink(&(bssp->descr), vnp);
      }
    }
    if (!is_date) {
      vnp = ValNodeExtractList(&descr, Seq_descr_update_date);
      if (vnp != NULL) {
        bssp->descr = ValNodeLink(&(bssp->descr), vnp);
      }
    }
    /*    vnp = ValNodeExtractList(&descr, Seq_descr_pub);
     if (!is_pub)
     bssp->descr = ValNodeLink(&(bssp->descr), vnp);
     */
  }
  if (seqsep->choice == 1) {
    bsp = seqsep->data.ptrvalue;
    bsp->descr = descr;
  }
  if (seqsep->choice == 2) {
    tmp = seqsep->data.ptrvalue;
    tmp->descr = descr;
  }
  return;
}

void MoveNPPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

  BioseqSetPtr  bssp, tmp;
  BioseqPtr     bsp;
  SeqEntryPtr   seqsep;
  ValNodePtr    descr = NULL, vnp = NULL, vnext, v, v_copy;
  PubdescPtr    pdp, pdp_copy;

  if (IS_Bioseq(sep)) {
    return;
  }
  bssp = (BioseqSetPtr)(sep->data.ptrvalue);
  if (bssp->_class != 1) {    /*  do the rest for nuc-prot only */
    return;
  }
  seqsep = bssp->seq_set;
  if (seqsep == NULL) {
    return;
  }
  if (seqsep->choice == 1) {
    bsp = seqsep->data.ptrvalue;
    descr = bsp->descr;
  }
  if (seqsep->choice == 2) {
    tmp = seqsep->data.ptrvalue;
    descr = tmp->descr;
  }
  if (bssp->descr == NULL) {
    bssp->descr = GetDescrNoTitles(&descr);
  } else {
    /* move pubs to nuc-prot level */
    vnp = ValNodeExtractList(&descr, Seq_descr_pub);
    for (v=vnp; v; v=vnext) {
      vnext = v->next;
      pdp = (PubdescPtr) v->data.ptrvalue;
      if (pdp->num != NULL || pdp->name != NULL || pdp->fig != NULL
          || pdp->comment !=NULL) {
        pdp_copy = AsnIoMemCopy(pdp, (AsnReadFunc) PubdescAsnRead,
                                (AsnWriteFunc) PubdescAsnWrite);
        v_copy = SeqDescrNew(NULL);
        v_copy->choice = Seq_descr_pub;
        v_copy->data.ptrvalue = pdp_copy;
        descr = ValNodeLink(&(descr), v_copy);
        PubdescFree (pdp);
        vnp = remove_node(vnp, v);
      }
    }
    if (vnp != NULL) {
      bssp->descr = ValNodeLink(&(bssp->descr), vnp);
    }
  }
  if (seqsep->choice == 1) {
    bsp = seqsep->data.ptrvalue;
    bsp->descr = descr;
  }
  if (seqsep->choice == 2) {
    tmp = seqsep->data.ptrvalue;
    tmp->descr = descr;
  }
  return;
}

static SeqDescrPtr GetSeqDescFromSeqEntry (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sep == NULL) return NULL;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return NULL;
    return bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return NULL;
    return bssp->descr;
  }

  return NULL;
}

/* return list of pubs that are the same on all pop/phy/mut components */
static SeqDescrPtr CheckSegsForPopPhyMut (SeqEntryPtr sep)

{
  SeqDescrPtr    descr;
  ValNodePtr     head;
  ValNodePtr     last;
  SeqDescrPtr    list = NULL;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp1, pdp2;
  Boolean        same;
  SeqDescrPtr    sdp1, sdp2;
  SeqEntryPtr    tmp;
  ValNodePtr     vnp, vnp1, vnp2;

  for (sdp1 = GetSeqDescFromSeqEntry (sep); sdp1 != NULL; sdp1 = sdp1->next) {
    if (sdp1->choice != Seq_descr_pub) continue;
    pdp1 = (PubdescPtr) sdp1->data.ptrvalue;
    if (pdp1 == NULL) continue;
    head = NULL;
    last = NULL;
    for (tmp = sep->next, same = FALSE; tmp != NULL; tmp = tmp->next) {
      for (sdp2 = GetSeqDescFromSeqEntry (tmp); sdp2 != NULL; sdp2 = sdp2->next) {
        if (sdp2->choice != Seq_descr_pub) continue;
        pdp2 = (PubdescPtr) sdp2->data.ptrvalue;
        if (pdp2 == NULL) continue;
        if (PubLabelMatchEx (pdp1->pub, pdp2->pub, FALSE) == 0) {
          if (PubdescMatch (pdp1, pdp2)) {
            same = TRUE;
            vnp = ValNodeAddPointer (&last, 0, (Pointer) sdp2);
            if (head == NULL) {
              head = vnp;
            }
            last = vnp;
            break;
          }
        }
      }
      if (sdp2 == NULL) {
        same = FALSE;
        break;
      }
    }
    if (same) {
      descr = SeqDescrNew (NULL);
      descr->choice = Seq_descr_pub;
      descr->data.ptrvalue = AsnIoMemCopy (pdp1,
                                           (AsnReadFunc) PubdescAsnRead,
                                           (AsnWriteFunc) PubdescAsnWrite);
      list = tie_next (list, descr);
      /* mark original pubs for deletion */
      if (sdp1->extended) {
        ovp = (ObjValNodePtr) sdp1;
        ovp->idx.deleteme = 1;
      }
      for (vnp1 = head; vnp1 != NULL; vnp1 = vnp1->next) {
        vnp2 = (ValNodePtr) vnp1->data.ptrvalue;
        if (vnp2->extended) {
          ovp = (ObjValNodePtr) vnp2;
          ovp->idx.deleteme = 1;
        }
      }
    }
    ValNodeFree (head);
  }

  return list;
}

/* move identical pubs in pop/phy/mut components to the set level */
static void MovePopPhyMutPubsProc (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;
  ValNodePtr    v, pub, vv, next;
  PubdescPtr    pdp, pdpv;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if ((bssp->_class < BioseqseqSet_class_mut_set ||
       bssp->_class > BioseqseqSet_class_eco_set) &&
      bssp->_class != BioseqseqSet_class_wgs_set &&
      bssp->_class != BioseqseqSet_class_small_genome_set) return;
  pub = CheckSegsForPopPhyMut (bssp->seq_set);
  if (pub == NULL) return;
  /* check if pub is already on the set descr */
  for(v=bssp->descr; v != NULL; v = v->next) {
    if (v->choice != Seq_descr_pub)
      continue;
    for (vv = pub; vv; vv = next) {
      next = vv->next;
      pdp = vv->data.ptrvalue;
      pdpv = v->data.ptrvalue;
      if (PubLabelMatchEx (pdp->pub, pdpv->pub, FALSE) == 0) {
        PubdescFree(pdp);
        pub = remove_node(pub, vv);
      }
    }
  }

  bssp->descr = tie_next(bssp->descr, pub);
}

void MovePopPhyMutPubs (SeqEntryPtr sep)

{
  if (sep == NULL) return;
  SeqEntryExplore(sep, (Pointer) NULL, MovePopPhyMutPubsProc);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, sep);
}

static void AddFeatToBioseq (SeqFeatPtr sfp, BioseqPtr bsp)

{
  SeqFeatPtr   prev;
  SeqAnnotPtr  sap;

  if (sfp == NULL || bsp == NULL) return;
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
  }
  sap = bsp->annot;
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

static void MoveFeatsOnPartsProc (BioseqSetPtr bssp, Pointer userdata)

{
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  BioseqPtr      target;

  if (bssp == NULL || bssp->_class != BioseqseqSet_class_parts) return;

  sap = bssp->annot;
  prevsap = (Pointer PNTR) &(bssp->annot);

  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        /* target = GetBioseqGivenSeqLoc (sfp->location, sfp->idx.entityID); */
        target = BioseqFindFromSeqLoc (sfp->location);
        if (target != NULL) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          AddFeatToBioseq (sfp, target);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void MoveFeatsFromPartsSet (SeqEntryPtr sep)

{
  VisitSetsInSep (sep, NULL, MoveFeatsOnPartsProc);
}

Boolean CmpOrgById(BioSourcePtr b1, BioSourcePtr b2)
{
  DbtagPtr    d1 = NULL, d2 = NULL;
  ValNodePtr  vnp;

  if (b1 == NULL || b2 == NULL) {
    return FALSE;
  }
  if (b1->org ==  NULL || b2->org == NULL) {
    return FALSE;
  }
  for (vnp = b1->org->db; vnp; vnp = vnp->next) {
    d1 = (DbtagPtr) vnp->data.ptrvalue;
    if (StringCmp(d1->db, "taxon") == 0) {
      break;
    }
  }
  for (vnp = b2->org->db; vnp; vnp = vnp->next) {
    d2 = (DbtagPtr) vnp->data.ptrvalue;
    if (StringCmp(d2->db, "taxon") == 0) {
      break;
    }
  }
  if (d1 && d2) {
    if (d1->tag->id == d2->tag->id) {
      return TRUE;
    } else {
    }
  } else if (StringICmp(b1->org->taxname, b2->org->taxname) == 0) {
    return TRUE;
  }
  return FALSE;
}

BioSourcePtr BioSourceMerge(BioSourcePtr host, BioSourcePtr guest)
{
  SubSourcePtr  ssp, sp;
  OrgModPtr     omp, homp;
  OrgNamePtr    onp;

  if (host == NULL && guest == NULL) {
    return NULL;
  }
  if (host == NULL && guest != NULL) {
    host = AsnIoMemCopy(guest, (AsnReadFunc) BioSourceAsnRead,
                        (AsnWriteFunc) BioSourceAsnWrite);
    return host;
  }
  if (host != NULL && guest == NULL) {
    return host;
  }
  if (host->genome == 0 && guest->genome != 0) {
    host->genome = guest->genome;
  }
  if (host->origin == 0 && guest->origin != 0) {
    host->origin = guest->origin;
  }
  for (ssp = guest->subtype; ssp; ssp = ssp->next) {
    sp = AsnIoMemCopy(ssp, (AsnReadFunc) SubSourceAsnRead,
                      (AsnWriteFunc) SubSourceAsnWrite);
    host->subtype = tie_next_subtype(host->subtype, sp);
  }
  if (guest->org->orgname) {
    for (omp = guest->org->orgname->mod; omp; omp = omp->next) {
      homp = AsnIoMemCopy(omp, (AsnReadFunc) OrgModAsnRead,
                          (AsnWriteFunc) OrgModAsnWrite);
      if ((onp = host->org->orgname)    == NULL) {
        onp = OrgNameNew();
        host->org->orgname = onp;
      }
      onp->mod = tie_next_OrgMod(onp->mod, homp);
    }
  }
  return host;
}

BioSourcePtr BioSourceCommon(BioSourcePtr host, BioSourcePtr guest)
{
  SubSourcePtr  ssp, sp, spnext;
  OrgModPtr     omp, om, ompnext;

  if (host->genome != guest->genome) {
    host->genome = 0;
  }
  if (host->origin != guest->origin) {
    host->origin = 0;
  }
  for (sp = host->subtype; sp; sp = spnext) {
    spnext = sp->next;
    for (ssp = guest->subtype; ssp; ssp = ssp->next) {
      if (sp->subtype == ssp->subtype &&
          StringCmp(sp->name, ssp->name) == 0) {
        break;
      }
    }
    if (ssp == NULL) {
      host->subtype = remove_subtype(host->subtype, sp);
    }
  }
  if (CmpOrgById(host, guest) == FALSE) {
    OrgRefFree(host->org);
    host->org = NULL;
    return host;
  }
  if (StringExtCmp(host->org->common, guest->org->common) != 0) {
    MemFree(host->org->common);
    host->org->common = NULL;
  }
  if (guest->org->orgname == NULL) {
    MemFree(host->org->orgname);
    host->org->orgname = NULL;
  } else {
    if (host->org->orgname) {
      for (omp = host->org->orgname->mod; omp; omp = ompnext) {
        ompnext = omp->next;
        for (om = guest->org->orgname->mod; om; om = om->next) {
          if (om->subtype == omp->subtype &&
              StringCmp(om->subname, omp->subname) == 0) {
            break;
          }
        }
      }
      if (omp == NULL) {
        host->org->orgname->mod =
        remove_OrgMod(host->org->orgname->mod, omp);
      }
    }
  }
  return host;
}

static Boolean EmptyBioSource(BioSourcePtr bio)
{
  if (bio == NULL || bio->org == NULL) return TRUE;
  if (bio->org->taxname == NULL && bio->org->common == NULL && bio->org->db == NULL) return TRUE;
  return FALSE;
}

void StripBSfromTop (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqSetPtr  bssp, tmp;
  ValNodePtr    vnp, bio_vnp = NULL;
  SeqEntryPtr   segsep, parts, cur;
  BioseqPtr     bsp;
  BioSourcePtr  biotop = NULL, bio = NULL;
  Boolean       first = TRUE;

  if (IS_Bioseq(sep)) {
    return;
  }
  bssp = (BioseqSetPtr)(sep->data.ptrvalue);
  if (bssp->_class != 2) {    /*  do the rest for segset only */
    return;
  }
  bio_vnp = ValNodeExtractList(&(bssp->descr), Seq_descr_source);
  if (bio_vnp) {
    biotop = ((BioSourcePtr)bio_vnp->data.ptrvalue);
  } else {
    bio_vnp = SeqDescrNew(NULL);
    bio_vnp->choice = Seq_descr_source;
  }
  segsep = bssp->seq_set;
  if (segsep->next == NULL) {
    return;
  }
  if (!IS_Bioseq(segsep->next)) {
    tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
    parts = tmp->seq_set;
    if (parts == NULL) {
      return;
    }
    for (cur = parts; cur; cur = cur->next) {
      bsp = cur->data.ptrvalue;
      for (vnp = bsp->descr; vnp; vnp=vnp->next) {
        if (vnp->choice == Seq_descr_source) {
          bio = ((BioSourcePtr)vnp->data.ptrvalue);
          break;
        }
      }
      if (bio) {
        if (biotop == NULL && first == TRUE) {
          biotop = AsnIoMemCopy(bio, (AsnReadFunc) BioSourceAsnRead,
                                (AsnWriteFunc) BioSourceAsnWrite);
          bio_vnp->data.ptrvalue = biotop;
          first = FALSE;
        } else {
          biotop = BioSourceCommon(biotop, bio);
        }
      }
    }
    if (biotop != NULL) {
      if (EmptyBioSource(biotop)) {
        BioSourceFree(biotop);
        ValNodeFree(bio_vnp);
      } else {
        bssp->descr = tie_next(bssp->descr, bio_vnp);
      }
    }
  }
  return;
}

void StripBSfromParts (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

  BioseqSetPtr  bssp, tmp;
  BioseqPtr     bsp;
  SeqEntryPtr   segsep, parts = NULL, cur;
  ValNodePtr    descr, vnp, set_vnp;
  BSMapPtr      bsmp;
  Boolean       empty;
  BioSourcePtr  biosp = NULL;

  bsmp = (BSMapPtr) data;
  if (IS_Bioseq(sep)) {
    return;
  }
  bssp = (BioseqSetPtr)(sep->data.ptrvalue);
  if (bssp->_class != 2) {    /*  do the rest for segset only */
    return;
  }
  vnp = ValNodeExtractList(&(bssp->descr), Seq_descr_source);
  if (vnp) {
    biosp = (BioSourcePtr) vnp->data.ptrvalue;
    ValNodeFree(vnp);
    vnp=NULL;
  }
  segsep = bssp->seq_set;
  if (segsep->next == NULL) {
    return;
  }
  set_vnp = bssp->descr;
  if (!IS_Bioseq(segsep->next)) {
    tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
    parts = tmp->seq_set;
  }
  for (cur = parts, empty = TRUE; cur; cur = cur->next) {
    bsp = cur->data.ptrvalue;
    vnp = ValNodeExtractList(&(bsp->descr), Seq_descr_source);
    if (vnp) {
      biosp = BioSourceMerge(biosp, vnp->data.ptrvalue);
      BioSourceFree((BioSourcePtr) (vnp->data.ptrvalue));
      ValNodeFree(vnp);
    }
  }
  if (biosp) {
    descr = SeqDescrNew(NULL);
    descr->choice = Seq_descr_source;
    descr->data.ptrvalue = biosp;
    bssp->descr = ValNodeLink(&(bssp->descr), descr);
  }

  return;
}

/*------------------------ GetDescr() --------------------------*/
/*****************************************************************************
 *  GetDescr:
 *                                                                    8-12-93
 ******************************************************************************/
ValNodePtr GetDescr(ValNodePtr PNTR descr)
{
  ValNodePtr  vnp, hvnp = NULL;

  vnp = ValNodeExtractList(descr, Seq_descr_title);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  vnp = ValNodeExtractList(descr, Seq_descr_org);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  if ( check_GIBB(*descr)) {
    vnp = ValNodeExtractList(descr, Seq_descr_modif);
    if (vnp != NULL) {
      hvnp = ValNodeLink(&hvnp, vnp);
    }
  }
  vnp = ValNodeExtractList(descr, Seq_descr_comment);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  vnp = ValNodeExtractList(descr, Seq_descr_pub);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  vnp = ValNodeExtractList(descr, Seq_descr_update_date);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  return (hvnp);

} /* GetDescr */

static ValNodePtr GetDescrNoTitles (ValNodePtr PNTR descr)
{
  ValNodePtr  vnp, hvnp = NULL;

  vnp = ValNodeExtractList(descr, Seq_descr_org);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  if ( check_GIBB(*descr)) {
    vnp = ValNodeExtractList(descr, Seq_descr_modif);
    if (vnp != NULL) {
      hvnp = ValNodeLink(&hvnp, vnp);
    }
  }
  vnp = ValNodeExtractList(descr, Seq_descr_comment);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  vnp = ValNodeExtractList(descr, Seq_descr_pub);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  vnp = ValNodeExtractList(descr, Seq_descr_update_date);
  if (vnp != NULL) {
    hvnp = ValNodeLink(&hvnp, vnp);
  }

  return (hvnp);

} /* GetDescrNoTitles */

/*------------------------ check_GIBB() --------------------------*/
/*****************************************************************************
 *  check_GIBB:
 *                                                                    8-12-93
 ******************************************************************************/
Boolean check_GIBB(ValNodePtr descr)
{
  ValNodePtr  vnp, modif;
  Int4        gmod;

  if (descr == NULL) {
    return FALSE;
  }
  for (vnp = descr; vnp && vnp->choice != Seq_descr_modif; vnp = vnp->next)
    continue;
  if (vnp == NULL) {
    return FALSE;
  }
  modif = (ValNodePtr) vnp->data.ptrvalue;
  if (modif == NULL) {
    return FALSE;
  }
  gmod = modif->data.intvalue;
  if (gmod == Seq_descr_GIBB_mod_dna || gmod == Seq_descr_GIBB_mod_rna ||
      gmod == Seq_descr_GIBB_mod_est || gmod == Seq_descr_GIBB_mod_complete
      || gmod == Seq_descr_GIBB_mod_partial) {
    return FALSE;
  }
  return TRUE;
}

/*----------------------------- SrchSegChoice() --------------------------*/
/*****************************************************************************
 *  SrchSegChoice:
 ******************************************************************************/
ValNodePtr SrchSegChoice(SeqEntryPtr sep, Uint1 choice)
{
  BioseqPtr   bsp;
  ValNodePtr  hvnp = NULL;

  if (sep == NULL) {
    return NULL;
  }
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    /* first bioseq from parts */
    if (CheckSegDescrChoice(sep, choice)) {     /*identical */
      hvnp = ValNodeExtractList(&(bsp->descr), choice);
      CleanUpSeqDescrChoice(sep->next, choice);
    }
  }
  return (hvnp);

} /* SrchSegChoice */

/*---------------------------- SrchSegSeqMol() --------------------------*/
/*************************************************************************
 *  SrchSegSeqMol:
 *                                                             5-14-93
 **************************************************************************/
void SrchSegSeqMol(SeqEntryPtr sep)
{
  BioseqPtr    bsp = NULL;
  SeqEntryPtr  cursep;
  Uint1        mol;
  /*
  CharPtr     str1, str2;
  */

  if (sep == NULL || sep->next) {
    return;
  }
  if (IS_Bioseq(sep)) {
    bsp = sep->data.ptrvalue;
  }
  if (bsp == NULL) {
    return;
  }
  mol = bsp->mol;

  for (cursep = sep->next; cursep != NULL; cursep = cursep->next) {
    if (IS_Bioseq(sep)) {
      bsp = cursep->data.ptrvalue;
    } else {
      continue;
    }
    if (mol != bsp->mol) {
      break;
    }
  }

  return;

} /* SrchSegSeqMol */

/*------------------------ CheckSegDescrChoice() -------------------------*/
/*****************************************************************************
 *  CheckSegDescrChoice:
 *                                                                  5-18-93
 ******************************************************************************/
Boolean CheckSegDescrChoice(SeqEntryPtr sep, Uint1 choice)
{
  BioseqPtr     bsp;
  SeqEntryPtr   cursep;
  ValNodePtr    vnp, mvnp;
  Boolean       same;
  Boolean       no_choice = TRUE;
  BioSourcePtr  biosp = NULL;
  OrgRefPtr     orp;
  CharPtr       title = NULL;
  DatePtr       dp = NULL;
  CharPtr       org = NULL;
  Int4          modif = -1, mol = -1;
  PubdescPtr    pdp = NULL;

  for (cursep = sep, same = TRUE;
       cursep != NULL && same; cursep = cursep->next) {
    bsp = cursep->data.ptrvalue;
    for (vnp = bsp->descr; vnp != NULL && vnp->choice != choice;
         vnp = vnp->next)
      continue;

    if (vnp == NULL) {
      same = FALSE;
      no_choice = TRUE;
    } else if (choice == Seq_descr_org) {
      no_choice = FALSE;
      orp = vnp->data.ptrvalue;

      if (org == NULL)
        org = orp->taxname;
      else if (StringCmp(org, orp->taxname) != 0)
        same = FALSE;
    } else if (choice == Seq_descr_source) {
      no_choice = FALSE;
      if (biosp == NULL) {
        biosp = vnp->data.ptrvalue;
      } else if (BSComparison(biosp,
                              (BioSourcePtr) vnp->data.ptrvalue) != 0) {
        same = FALSE;
      }
    } else if (choice == Seq_descr_mol_type) {
      no_choice = FALSE;
      if (mol == -1)
        mol = vnp->data.intvalue;
      else if (mol != vnp->data.intvalue)
        same = FALSE;
    } else if (choice == Seq_descr_modif) {
      no_choice = FALSE;
      mvnp = vnp->data.ptrvalue;

      if (modif == -1)
        modif = mvnp->data.intvalue;
      else if (modif != mvnp->data.intvalue)
        same = FALSE;
    }else if (choice == Seq_descr_update_date) {
      no_choice = FALSE;
      if (dp == NULL)
        dp = vnp->data.ptrvalue;
      else if (DateMatch(dp, vnp->data.ptrvalue, TRUE) != 0)
        same = FALSE;
    } else if (choice == Seq_descr_pub) {
      no_choice = FALSE;
      if (pdp == NULL)
        pdp = vnp->data.ptrvalue;
      else if (PubMatch(pdp->pub,
                        ((PubdescPtr)(vnp->data.ptrvalue))->pub) != 0)
        same = FALSE;
    } else if (choice == Seq_descr_title) {
      no_choice = FALSE;
      if (title == NULL)
        title = vnp->data.ptrvalue;
      else if (StringCmp(title, (CharPtr) vnp->data.ptrvalue) != 0)
        same = FALSE;
    } else {
      no_choice = FALSE;
      ErrPostEx(SEV_WARNING, 0, 2,
                "Unrecognized choice: %d", choice);
      same = FALSE;
    }
  }
  if (same == FALSE && no_choice == TRUE && choice != Seq_descr_update_date) {
    same = TRUE;
  }
  return (same);

} /* CheckSegDescrChoice */


void StripProtXref (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr       bsp, prot = NULL;
  BioseqSetPtr    bssp;
  ValNodePtr      vnp = NULL;
  SeqAnnotPtr     sap = NULL, ap, pap = NULL;
  SeqFeatPtr      sfp, psfp = NULL, head;
  ProtRefPtr      prp, pprp = NULL;
  SeqFeatXrefPtr  xrp, xrpnext;
  SeqIdPtr        sid;
  SeqLocPtr       slp;
  TextSeqIdPtr    tsip;

  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    sap = bsp->annot;
    if (bsp != NULL) {
      for (sid = bsp->id; sid != NULL; sid = sid->next) {
        if (sid->choice != SEQID_OTHER) continue;
        tsip = (TextSeqIdPtr) sid->data.ptrvalue;
        if (tsip == NULL) continue;
        if (StringNCmp (tsip->accession, "NC_", 3) == 0) return;
        if (StringNCmp (tsip->accession, "NG_", 3) == 0) return;
        if (StringNCmp (tsip->accession, "NT_", 3) == 0) return;
        if (StringNCmp (tsip->accession, "NW_", 3) == 0) return;
        if (StringNCmp (tsip->accession, "AC_", 3) == 0) return;
      }
    }
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    sap = bssp->annot;
  }
  for (ap = sap; ap != NULL; ap = ap->next) {
    if (ap->type != 1) continue;
    head = (SeqFeatPtr)(ap->data);
    for (sfp = head; sfp; sfp = sfp->next) {
      if (sfp->data.choice != SEQFEAT_CDREGION) continue;
      prot = NULL;
      pap = NULL;
      psfp = NULL;
      pprp = NULL;
      vnp = sfp->product;
      if (vnp != NULL) {
        if (vnp->choice == SEQLOC_WHOLE) {
          sid = vnp->data.ptrvalue;
          prot = BioseqFind(sid);
        }
      }
      if (prot != NULL) {
        pap = prot->annot;
        if (pap != NULL) {
          for (psfp = pap->data; psfp; psfp=psfp->next) {
            if (psfp->data.choice == SEQFEAT_PROT) {
              pprp = psfp->data.value.ptrvalue;
              break;
            }
          }
        }
      }
      if (vnp != NULL) {         /* sfp->product != NULL */
        for (xrp = sfp->xref; xrp != NULL; xrp = xrpnext) {
          xrpnext = xrp->next;
          if (xrp->data.choice == SEQFEAT_PROT) {
            prp = xrp->data.value.ptrvalue;
            if (pap != NULL && pprp == NULL) {
              if (psfp == NULL) {
                psfp = SeqFeatNew();
                psfp->data.choice = SEQFEAT_PROT;
                slp = ValNodeNew(NULL);
                slp->choice = SEQLOC_WHOLE;
                slp->data.ptrvalue = SeqIdDup(sid);
                psfp->location = slp;
                pap->data = tie_feat(pap->data, psfp);
              }
              psfp->data.value.ptrvalue = AsnIoMemCopy(prp,
                                                       (AsnReadFunc) ProtRefAsnRead,
                                                       (AsnWriteFunc) ProtRefAsnWrite);
            }
            sfp->xref = remove_xref(sfp->xref, xrp);
          }
        }
      }
    }
    ap->data = head;
  }
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    bsp->annot = sap;
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    bssp->annot = sap;
  }
  return;
}

static SeqLocPtr GetAnticodonFromObject(SeqFeatPtr sfp)
{
  UserObjectPtr  usop;
  UserFieldPtr   ufp;
  Int4Ptr        ints;
  SeqLocPtr      slp;
  SeqIntPtr      sip;
  Int4           from, to;

  if (sfp == NULL) {
    return NULL;
  }
  if ((usop = sfp->ext) == NULL) {
    return NULL;
  }
  if (StringICmp (usop->_class, "NCBI") != 0) {
    return NULL;
  }
  ufp = usop->data;
  if (ufp && ufp->choice == 8) {  /* ints */
    ints = (Int4Ptr) ufp->data.ptrvalue;
    from = ints[0];
    to = ints[1];
  }
  sip = SeqIntNew();
  sip->from = from;
  sip->to = to;
  sip->id = SeqIdDup(SeqLocId(sfp->location));
  slp = ValNodeNew(NULL);
  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = sip;
  sfp->ext = usop->next /* NULL */;
  UserObjectFree (usop);
  return slp;

}

/*--------------------------- CheckMaps() --------------------------*/
/***************************************************************************
 *   CheckMaps:
 *   -- find all /map and Gene-ref
 *        if all maps are the same put it to Biosource and remove quals
 *   -- change User-Object anticodon in tRNA to SeqLoc
 ****************************************************************************/
void CheckMaps (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioSourcePtr   biop;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     descr;
  SeqAnnotPtr    sap = NULL, ap;
  CharPtr        qval;
  SeqFeatPtr     sfp;
  GeneRefPtr     grp;
  QualMapPtr     qmp;
  RnaRefPtr      rrp;
  SubSourcePtr   ssp;
  tRNAPtr        trna;
  GBQualPtr      q, qnext;

  qmp = data;
  if (qmp->same == FALSE) {
    return;
  }
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
    sap = bsp->annot;
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
    sap = bssp->annot;
  }

  while (descr != NULL) {
    if (descr->choice == Seq_descr_source) {
      biop = (BioSourcePtr) descr->data.ptrvalue;
      if (biop != NULL) {
        for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
          if (ssp->subtype == SUBSRC_map && ssp->name != NULL) {
            if (qmp->name == NULL) {
              qmp->name = StringSave (ssp->name);
            } else if (StringCmp (qmp->name, ssp->name) != 0) {
              qmp->same = FALSE;
            }
          }
        }
      }
    }
    descr = descr->next;
  }

  /* look for all the same maploc and place it to SubSource*/
  for (ap = sap; ap != NULL; ap = ap->next) {
    if (ap->type != 1) {
      continue;
    }
    for (sfp = (SeqFeatPtr)(ap->data); sfp; sfp = sfp->next) {
      if ((qval = get_qvalue(sfp->qual, "map")) != NULL) {
        if (qmp->name == NULL) {
          qmp->name = StringSave(qval);
        } else if (StringCmp(qmp->name, qval) != 0) {
          qmp->same = FALSE;
          break;
        }
      }
      if (sfp->data.choice == SEQFEAT_GENE) {
        grp = sfp->data.value.ptrvalue;
        if ((qval = grp->maploc) != NULL) {
          if (qmp->name == NULL) {
            qmp->name = StringSave(qval);
          } else if (StringCmp(qmp->name, qval) != 0) {
            qmp->same = FALSE;
            break;
          }
        }
      }
      if (sfp->data.choice == SEQFEAT_RNA) {
        rrp = sfp->data.value.ptrvalue;
        if (rrp->type == 3 && rrp->ext.choice == 2) {
          trna = rrp->ext.value.ptrvalue;
          if (sfp->ext != NULL && trna->anticodon == NULL) {
            trna->anticodon = GetAnticodonFromObject(sfp);
            for (q = sfp->qual; q; q = qnext) {
              qnext = q->next;
              if (StringCmp(q->qual, "anticodon") == 0) {
                sfp->qual = remove_qual(sfp->qual, q);
              }
            }
          }
        }
      }
    }
  }
  return;
}

void StripMaps (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  BioSourcePtr    biosp = NULL;
  SubSourcePtr    ssp;
  ValNodePtr      descr = NULL, vnp;
  SeqAnnotPtr     sap = NULL, ap;
  CharPtr         qval = NULL;
  SeqFeatPtr      sfp;
  GeneRefPtr      grp;
  QualMapPtr      qmp;
  SeqFeatXrefPtr  xrp;

  qmp = data;
  if (qmp->same == FALSE || qmp->name == NULL) {
    return;
  }
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
    sap = bsp->annot;
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
    sap = bssp->annot;
  }
  for ( vnp = descr; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Seq_descr_source) {
      biosp = vnp->data.ptrvalue;
      break;
    }
  }
  for (ap = sap; ap != NULL; ap = ap->next) {
    if (ap->type == 1) {
      for (sfp = (SeqFeatPtr)(ap->data); sfp; sfp = sfp->next) {
        qval = qvalue_extract(&(sfp->qual), "map");
        if(qval != NULL)
        {
          MemFree(qval);
          qval = NULL;
        }
        if (sfp->data.choice == SEQFEAT_GENE) {
          grp = sfp->data.value.ptrvalue;
          if(grp->maploc != NULL)
          {
            MemFree(grp->maploc);
            grp->maploc = NULL;
          }
        }
        for (xrp = sfp->xref; xrp != NULL; xrp = xrp->next) {
          if (xrp->data.choice == SEQFEAT_GENE) {
            grp = xrp->data.value.ptrvalue;
            if(grp->maploc != NULL)
            {
              MemFree(grp->maploc);
              grp->maploc = NULL;
            }
          }
        }
      }
    }
  }
  if (biosp != NULL) {  /* has biosource */
    ssp = SubSourceNew();
    ssp->subtype = 2;    /*map */
    ssp->name = StringSave(qmp->name);
    biosp->subtype = tie_next_subtype(biosp->subtype, ssp);
  }
  /*    if (qval) {
   MemFree(qval);
   }*/
  return;
}

static Boolean GBQualPresent(CharPtr ptr, GBQualPtr gbqual)

{
  Boolean    present=FALSE;
  GBQualPtr  qual;

  for (qual=gbqual; qual != NULL; qual=qual->next)
    if (StringCmp(ptr, qual->qual) == 0)
    {
      present = TRUE;
      break;
    }

  return present;
}

static CharPtr ExamineGBQual (CharPtr ptr, GBQualPtr gbqual)

{
  GBQualPtr  qual;

  for (qual = gbqual; qual != NULL; qual = qual->next) {
    if (StringCmp (ptr, qual->qual) == 0) return qual->val;
  }

  return NULL;
}

void MapsToGenref (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    descr = NULL, head, last, vnp;
  SeqAnnotPtr   sap = NULL, ap;
  CharPtr       qval= NULL, name;
  Boolean       same;
  SeqFeatPtr    sfp, cur;
  SeqLocPtr     loc;
  GeneRefPtr    grp;

  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
    sap = bsp->annot;
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
    sap = bssp->annot;
  }

  /* maps are different */
  for (ap = sap; ap != NULL; ap = ap->next) {
    if (ap->type == 1) {
      for (sfp = (SeqFeatPtr)(ap->data); sfp; sfp = sfp->next) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          grp = sfp->data.value.ptrvalue;
          name = NULL;
          head = NULL;
          last = NULL;
          same = TRUE;
          for (cur = (SeqFeatPtr)(ap->data); cur; cur = cur->next) {
            if ((GBQualPresent("map", cur->qual)) == FALSE) {
              continue;
            }
            if (SeqLocAinB(cur->location, sfp->location) < 0 ) {
              continue;
            }
            vnp = ValNodeAddPointer (&last, 0, cur);
            if (head == NULL) {
              head = vnp;
            }
            last = vnp;
          }
          for (vnp = head; vnp != NULL; vnp = vnp->next) {
            cur = (SeqFeatPtr) vnp->data.ptrvalue;
            qval = ExamineGBQual("map", (cur->qual));
            if (name == NULL) {
              name = qval;
            } else if (StringICmp (name, qval) != 0) {
              same = FALSE;
            }
          }
          if (same && name != NULL) {
            if (grp->maploc == NULL && name != NULL) {
              grp->maploc = StringSave(name);
            }
            name = NULL;
            loc = NULL;
            for (vnp = head; vnp != NULL; vnp = vnp->next) {
              cur = (SeqFeatPtr) vnp->data.ptrvalue;
              if (cur == NULL) continue;
              qval = qvalue_extract(&(cur->qual), "map");
              MemFree (qval);
            }
          }
          ValNodeFree (head);
        } /* if SEQFEAT_GENE */
      }
    } /* if ftable */
  }
  return;
}

static Boolean CheckMinPub(ValNodePtr pub, Boolean is_ref_seq_prot)
{
  CitGenPtr  gen;

  if (pub == NULL) {
    return TRUE;
  }
  if (pub->choice == PUB_Muid || pub->choice == PUB_PMid) {
    if (pub->next == NULL) {
      if (is_ref_seq_prot) return FALSE;
      return TRUE;
    } else {
      return (CheckMinPub(pub->next, is_ref_seq_prot));
    }
  }
  if (pub->choice == PUB_Gen) {
    gen = pub->data.ptrvalue;
    if (gen->cit != NULL && gen->journal == NULL && gen->authors == NULL
        && gen->volume == NULL && gen->pages == NULL) {
      if (pub->next == NULL) {
        return TRUE;
      } else {
        return (CheckMinPub(pub->next, FALSE));
      }
    }
  }
  return FALSE;
}

static Boolean OkayToFuseRemarks (CharPtr com1, CharPtr com2)

{
  if (com1 != NULL && com2 != NULL) {
    if (StringICmp (com1, com2) != 0) return FALSE;
  }

  return TRUE;
}

static ValNodePtr AddToListEx (ValNodePtr list, ValNodePtr check, PubdescPtr pdp, Boolean is_ref_seq_prot)
{
  ValNodePtr    v, vnext;
  PubdescPtr    vpdp;
  PubStructPtr  psp;
  ValNodePtr    pubequ1 = NULL, pubequ2 = NULL;
  Boolean       is_1;

  if (pdp == NULL) {
    return NULL;
  }
  for (v = check; v != NULL; v = v->next) {
    psp = v->data.ptrvalue;
    if (psp->start != 2) {
      continue;
    }
    if (PubLabelMatchEx (psp->pub, pdp->pub, TRUE) == 0) {
      return list;
    }
  }
  if (pdp->name == NULL && pdp->fig == NULL && pdp->fig == NULL) {
    if (CheckMinPub(pdp->pub, is_ref_seq_prot) == TRUE) {   /* do not add minimum pub */
      return list;
    }
  }
  for (v = list; v != NULL; v = vnext) {
    vnext = v->next;
    vpdp = v->data.ptrvalue;
    if (vpdp->pub->next != NULL) {
      pubequ1 = SeqDescrNew(NULL);
      is_1 = TRUE;
      pubequ1->choice = PUB_Equiv;
      pubequ1->data.ptrvalue = vpdp->pub;
    } else {
      is_1 = FALSE;
      pubequ1 = vpdp->pub;
    }
    if (pdp->pub->next != NULL) {
      pubequ2 = SeqDescrNew(NULL);
      pubequ2->choice = PUB_Equiv;
      pubequ2->data.ptrvalue = pdp->pub;
    } else {
      pubequ2 = pdp->pub;
    }
    if (PubLabelMatchEx (pubequ1, pubequ2, TRUE) == 0 && OkayToFuseRemarks (pdp->comment, vpdp->comment)) {
      if (pdp->reftype == 2 && vpdp->reftype == 1) {
        vpdp->reftype = 2;
      }
      if (pdp->reftype == 1 && vpdp->reftype == 2) {
        pdp->reftype = 2;
      }
      if (vpdp->comment != NULL && pdp->comment == NULL) {
        pdp->comment = StringSave (vpdp->comment);
      } else if (vpdp->comment == NULL && pdp->comment != NULL) {
        vpdp->comment = StringSave (pdp->comment);
      }
      if (SelectBestPub(pubequ1, pubequ2) >= 0) {
        if (is_1) {
          ValNodeFree(pubequ1);
        }
        if (pdp->pub->next != NULL) {
          ValNodeFree(pubequ2);
        }
        return list;
      }
      PubdescFree((PubdescPtr) (v->data.ptrvalue));
      list = remove_node(list, v);
    }
    if (is_1) {
      ValNodeFree(pubequ1);
    }
    if (pdp->pub->next != NULL) {
      ValNodeFree(pubequ2);
    }
  }
  if (pdp && pdp->pub) {
    v = SeqDescrNew(NULL);
    v->choice = Seq_descr_pub;
    v->data.ptrvalue = AsnIoMemCopy(pdp, (AsnReadFunc) PubdescAsnRead,
                                    (AsnWriteFunc) PubdescAsnWrite);
  }
  list = ValNodeLink(&list, v);
  /*  may be sort ???? */
  return list;
}

ValNodePtr AddToList(ValNodePtr list, ValNodePtr check, PubdescPtr pdp)
{
  return AddToListEx (list, check, pdp, TRUE);
}

void CheckCitSubNew(ValNodePtr vnp)
{
  CitSubPtr    csp;
  AuthListPtr  alp;
  ImprintPtr   imp;

  if (vnp == NULL)
    return;
  if (vnp->choice != PUB_Sub)
    return;
  csp = (CitSubPtr) vnp->data.ptrvalue;
  if (csp != NULL) {
    alp = csp->authors;
    imp = csp->imp;
    if (alp != NULL && alp->affil == NULL &&
        imp != NULL && imp->pub != NULL) {
      alp->affil = imp->pub;
      imp->pub = NULL;
    }
    if (csp->date == NULL && imp != NULL && imp->date != NULL) {
      csp->date = imp->date;
      imp->date = NULL;
    }
    if (imp != NULL && imp->pub == NULL) {
      csp->imp = ImprintFree (csp->imp);
    }
  }
  return;
}

void ChangeCitSub (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  ValNodePtr    descr = NULL, vnp, v;
  PubdescPtr    pdp;

  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
  }
  for (vnp = descr; vnp; vnp=vnp->next) {
    if (vnp->choice != Seq_descr_pub) {
      continue;
    }
    if ((pdp = vnp->data.ptrvalue) == NULL) {
      continue;
    }
    for (v = pdp->pub; v; v=v->next) {
      if (v->choice == PUB_Sub) {
        CheckCitSubNew(v);
      }
    }
  }
}

/***************************************************************************
 *   NewPubs:
 *   -- find all ImpFeat "sites"
 *        change to pubdesc with reftype 'sites'
 *    -- find all other sfp->cit
 *        change to pubdesc with reftype 'feats'
 *    -- pubs are moved from SeqAnnot to Seqdescr on the same level
 ****************************************************************************/
void NewPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  ValNodePtr    descr = NULL, pubset, tmp, pubequ;
  ValNodePtr    next_pubequ, pub, min_pub;
  SeqAnnotPtr   sap = NULL, ap, apnext;
  SeqFeatPtr    sfp, cur, curnext;
  ImpFeatPtr    ifp;
  PubdescPtr    pubdesc;
  ValNodePtr    check = NULL, np_list = NULL;
  SeqIdPtr      sip;
  Boolean       is_ref_seq_prot = FALSE;

  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
    sap = bsp->annot;
    if (ISA_aa (bsp->mol)) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_OTHER) {
          is_ref_seq_prot = TRUE;
        }
      }
    }
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
    sap = bssp->annot;
  }
  tmp = ValNodeExtractList(&descr, Seq_descr_pub);
  if (tmp != NULL) {
    np_list = ValNodeLink(&np_list, tmp);
  }
  for (ap = sap; ap != NULL; ap = apnext) {
    apnext = ap->next;
    if (ap->type != 1) {
      continue;
    }
    sfp = (SeqFeatPtr)(ap->data);
    for (cur = sfp; cur; cur = curnext) {
      curnext = cur->next;
      if (cur->cit == NULL) {
        continue;
      }
      pubset = cur->cit;
      pub = NULL;
      min_pub = NULL;
      pubequ = pubset->data.ptrvalue;
      while (pubequ) {
        next_pubequ = pubequ->next;
        pubdesc = PubdescNew();
        if (pubequ->choice == PUB_Equiv) {
          pubdesc->pub = pubequ->data.ptrvalue;
        } else {
          pubdesc->pub = pubequ;
        }
        if (cur->data.choice == SEQFEAT_IMP) {
          ifp = cur->data.value.ptrvalue;
          if (StringCmp(ifp->key, "Site-ref") == 0) {
            pubdesc->reftype = 1; /* sites */
            np_list = AddToListEx (np_list, check, pubdesc, is_ref_seq_prot);
            min_pub = MinimizePub(pubequ);
            pub = tie_next(pub, min_pub);
            MemFree(pubdesc);
          } else {
            pubdesc->reftype = 2;
            np_list = AddToListEx (np_list, check, pubdesc, is_ref_seq_prot);
            min_pub = MinimizePub(pubequ);
            pub = tie_next(pub, min_pub);
            MemFree(pubdesc);
          }
        } else {
          pubdesc->reftype = 2;
          np_list = AddToListEx (np_list, check, pubdesc, is_ref_seq_prot);
          min_pub = MinimizePub(pubequ);
          pub = tie_next(pub, min_pub);
          MemFree(pubdesc);
        }
        PubFree(pubequ);
        pubequ = next_pubequ;
      }
      if (pub && pubset) {
        pubset->data.ptrvalue = pub;
      }
    }
    ap->data = sfp;
  }

  for (ap = sap; ap != NULL; ap = apnext) {
    apnext = ap->next;
    if (ap->data == NULL) {
      sap = remove_annot(sap, ap);
    }
  }

  if (bssp != NULL) {
    descr = tie_next(descr, np_list);
    bssp->descr = descr;
    bssp->annot = sap;
  } else {
    descr = tie_next(descr, np_list);
    bsp->descr = descr;
    bsp->annot = sap;
  }
  return;
}

void CmpPub (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr        bsp;
  PubdescPtr       pdp, vpdp;
  ValNodePtr       pubequ1, pubequ2, v;
  PubdescPtr PNTR  ppdp;
  Boolean          is_1;

  ppdp = data;
  pdp = *ppdp;
  if (pdp == NULL) {
    return;
  }
  if (!IS_Bioseq(sep)) {
    return;
  }
  bsp = (BioseqPtr)(sep->data.ptrvalue);
  if (bsp->mol == Seq_mol_aa) {
    return;
  }
  for(v = bsp->descr; v; v=v->next) {
    if (v->choice != Seq_descr_pub) {
      continue;
    }
    vpdp = v->data.ptrvalue;
    if (vpdp->pub->next != NULL) {
      pubequ1 = SeqDescrNew(NULL);
      is_1 = TRUE;
      pubequ1->choice = PUB_Equiv;
      pubequ1->data.ptrvalue = vpdp->pub;
    } else {
      is_1 = FALSE;
      pubequ1 = vpdp->pub;
    }
    if (pdp->pub->next != NULL) {
      pubequ2 = (NULL);
      pubequ2->choice = PUB_Equiv;
      pubequ2->data.ptrvalue = pdp->pub;
    } else {
      pubequ2 = pdp->pub;
    }
    if (PubMatch(pubequ1, pubequ2) == 0) {
      if (is_1) {
        ValNodeFree(pubequ1);
      }
      if (pdp->pub->next != NULL) {
        ValNodeFree(pubequ2);
      }
      break;
    }
    if (is_1) {
      ValNodeFree(pubequ1);
    }
    if (pdp->pub->next != NULL) {
      ValNodeFree(pubequ2);
    }
  }
  if (v == NULL) {
    PubdescFree(*ppdp);
    *ppdp = NULL;
  }

  return;
}

/***********************************************************************
 *    delete pubs from Bioseqs if they are already on the top level
 *    don't delete Pubdesc if additional info (name fig, num etc) is present
 ************************************************************************/
void DeletePubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr   bsp;
  PubdescPtr  pdp, vpdp;
  ValNodePtr  pubequ1, pubequ2, v, vnext, descr = NULL;
  Boolean     is_1, is_2;


  pdp = data;
  if (pdp == NULL) {
    return;
  }
  if (!IS_Bioseq(sep)) {
    return;
  }
  bsp = (BioseqPtr)(sep->data.ptrvalue);
  if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
    return;
  /*
   if (bsp->mol == Seq_mol_aa) {
   return;
   }
   */
  descr = bsp->descr;
  for(v = descr; v; v = vnext) {
    vnext = v->next;
    if (v->choice != Seq_descr_pub) {
      continue;
    }
    vpdp = v->data.ptrvalue;   /* from the Bioseq */
    if (vpdp->name != NULL || vpdp->fig != NULL
        || vpdp->num != NULL || vpdp->maploc != NULL
        || vpdp->comment != NULL) {
      continue;
    }
    if (vpdp->pub->next != NULL) {
      is_1 = TRUE;
      pubequ1 = ValNodeNew(NULL);
      pubequ1->choice = PUB_Equiv;
      pubequ1->data.ptrvalue = vpdp->pub;
    } else {
      is_1 = FALSE;
      pubequ1 = vpdp->pub;
    }
    if (pdp->pub->next != NULL) {  /* from the set */
      is_2 = TRUE;
      pubequ2 = ValNodeNew(NULL);
      pubequ2->choice = PUB_Equiv;
      pubequ2->data.ptrvalue = pdp->pub;
    } else {
      is_2 = FALSE;
      pubequ2 = pdp->pub;
    }
    if (PubMatch(pubequ1, pubequ2) == 0) {
      PubdescFree((PubdescPtr) (v->data.ptrvalue));
      descr = remove_node(descr, v);
    }
    if (is_1) {
      ValNodeFree(pubequ1);
    }
    if (is_2) {
      ValNodeFree(pubequ2);
    }
  }
  bsp->descr = descr;
  return;
}

void MoveSetPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp;
  Boolean        first;
  PubSetListPtr  psp;
  ValNodePtr     descr = NULL,v, vnext, tmp, set_list;
  PubdescPtr     tmp_pdp, pdp;

  psp = data;
  set_list = psp->list;
  first = psp->first;
  if (IS_Bioseq(sep) && (first == TRUE)) {
    bsp = sep->data.ptrvalue;
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
  } else {
    bssp = sep->data.ptrvalue;
    if (bssp->_class == 4 && bssp->seq_set != NULL) {
      bsp = (BioseqPtr) bssp->seq_set->data.ptrvalue;
      descr = bsp->descr;
    }
  }
  if (bsp == NULL) {
    return;
  }
  if (first == FALSE) {
    return;
  }
  for (v = descr; v; v=vnext) {
    vnext = v->next;
    if (v->choice == Seq_descr_pub) {
      pdp = v->data.ptrvalue;
      tmp_pdp = AsnIoMemCopy(pdp,
                             (AsnReadFunc) PubdescAsnRead, (AsnWriteFunc) PubdescAsnWrite);
      SeqEntryExplore(sep, &tmp_pdp, CmpPub);
      if (tmp_pdp != NULL) {
        tmp = SeqDescrNew(NULL);
        tmp->choice = Seq_descr_pub;
        tmp->data.ptrvalue = tmp_pdp;
        set_list = tie_next(set_list, tmp);
        SeqEntryExplore(sep, tmp_pdp, DeletePubs);
        first = FALSE;
      }
    }
  }
  psp->list = set_list;
  psp->first = first;
  data = psp;
  return;
}

void FindOldLineage (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    descr = NULL, vnp;
  GBBlockPtr    gb;
  CharPtr PNTR  linp;
  CharPtr       lineage;

  linp = (CharPtr PNTR) data;
  lineage = *linp;
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
  } else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
  }
  for (vnp = descr; vnp; vnp= vnp->next) {
    if (vnp->choice == Seq_descr_genbank) {
      gb = vnp->data.ptrvalue;
      if (gb->taxonomy) {
        if (*linp) {
          MemFree(*linp);
        }
        *linp = gb->taxonomy;
        gb->taxonomy = NULL;
      }
      break;
    }
  }
}

void NewLineage (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    descr = NULL, vnp;
  BioSourcePtr  biosp;
  OrgRefPtr     orp = NULL;
  OrgNamePtr    omp;
  CharPtr PNTR  linp;
  CharPtr       lineage;

  linp = (CharPtr PNTR) data;
  lineage = *linp;
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr)(sep->data.ptrvalue);
    if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
      return;
    descr = bsp->descr;
  }
  else {
    bssp = (BioseqSetPtr)(sep->data.ptrvalue);
    descr = bssp->descr;
  }
  for (vnp = descr; vnp; vnp= vnp->next) {
    if (vnp->choice == Seq_descr_source) {
      biosp = vnp->data.ptrvalue;
      orp = (OrgRefPtr) biosp->org;
      break;
    }
  }
  if (orp && lineage) {
    if (orp->orgname == NULL) {
      omp = OrgNameNew();
      orp->orgname = omp;
    }
    if(orp->orgname->lineage != NULL)
    {
      MemFree(orp->orgname->lineage);
    }
    orp->orgname->lineage = StringSave(lineage);
  }

}

/****************************************************************************
 *  delete_valnode:
 *****************************************************************************/
static ValNodePtr delete_valnode(ValNodePtr host, Uint1 choice)
{
  Boolean     first;
  ValNodePtr  curvnp, prevnp;

  for (curvnp = host, first = TRUE; curvnp != NULL
       && curvnp->choice != choice; curvnp = curvnp->next) {

    if (first) {
      prevnp = curvnp;
      first = FALSE;
    }
    else
      prevnp = prevnp->next;
  }

  if (curvnp == NULL) {
    return host;
  }
  if (first) {
    host = curvnp->next;
  } else {
    prevnp->next = curvnp->next;
  }
  curvnp->next = NULL;

  switch (choice) {
    case Seq_descr_org:
      OrgRefFree(curvnp->data.ptrvalue);
      break;
    case Seq_descr_modif:
      ValNodeFree(curvnp->data.ptrvalue);
      break;
    case Seq_descr_update_date:
      DateFree(curvnp->data.ptrvalue);
      break;
    case Seq_descr_mol_type:
      break;
    default:
      break;
  }

  ValNodeFree(curvnp);

  return host;
}

/*------------------- CleanUpSeqDescrChoice() -------------------------*/
/****************************************************************************
 *  CleanUpSeqDescrChoice:
 *                                                                    5-21-93
 *****************************************************************************/
void CleanUpSeqDescrChoice(SeqEntryPtr sep, Uint1 choice)
{
  BioseqPtr    bsp;
  SeqEntryPtr  cursep;

  for (cursep = sep; cursep != NULL; cursep = cursep->next) {
    bsp = cursep->data.ptrvalue;

    bsp->descr = delete_valnode(bsp->descr, choice);
  }

} /* CleanUpSeqDescrChoice */

/**********************************************************/
SeqDescrPtr remove_descr(SeqDescrPtr head, SeqDescrPtr x)
{
  SeqDescrPtr  v;
  SeqDescrPtr  p;

  if(head == NULL)
    return(NULL);

  if(x == head)
  {
    head = x->next;
    x->next = NULL;
    SeqDescFree(x);
    return(head);
  }
  for(v = head; v != NULL && v != x; v = v->next)
    p = v;

  if(v != NULL)
  {
    p->next = x->next;
    x->next = NULL;
    SeqDescFree(x);
  }
  return(head);
}


/* Cleanup functions originally from Sequin */

static void FindConsistentMolInfo (SeqEntryPtr sep, MolInfoPtr PNTR mipp, BoolPtr consist)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  MolInfoPtr    mip;
  ValNodePtr    sdp = NULL;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (*mipp == NULL) {
          *mipp = mip;
        } else {
          if ((*mipp)->biomol != mip->biomol ||
              (*mipp)->tech != mip->tech ||
              (*mipp)->completeness != mip->completeness ||
              StringICmp ((*mipp)->techexp, mip->techexp) != 0) {
            *consist = FALSE;
          }
        }
      }
    }
    sdp = sdp->next;
  }
  if (bssp == NULL) return;
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    FindConsistentMolInfo (sep, mipp, consist);
  }
}

static void RemoveMolInfoCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
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
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_molinfo) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

extern void NormalizeSegSeqMolInfo (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  Boolean       consistent;
  MolInfoPtr    master;
  MolInfoPtr    mip;
  ValNodePtr    sdp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp != NULL && bssp->_class != BioseqseqSet_class_segset) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        NormalizeSegSeqMolInfo (sep);
      }
      return;
    }
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
      mip = NULL;
      consistent = TRUE;
      FindConsistentMolInfo (sep, &mip, &consistent);
      if (mip != NULL && consistent) {
        master = MolInfoNew ();
        if (master == NULL) return;
        master->biomol = mip->biomol;
        master->tech = mip->tech;
        master->completeness = mip->completeness;
        master->techexp = StringSaveNoNull (mip->techexp);
        SeqEntryExplore (sep, NULL, RemoveMolInfoCallback);
        sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
        if (sdp != NULL) {
          sdp->data.ptrvalue = (Pointer) master;
        }
      }
    }
  }
}

static void CollectPseudoCdsProducts (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqContextPtr  bcp;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  CharPtr           label;
  size_t            len;
  BioseqPtr         product;
  SeqFeatPtr        prot;
  ProtRefPtr        prp;
  Boolean           pseudo;
  GBQualPtr         gbqual;
  SeqAnnotPtr       sap;
  SeqFeatPtr        sfp;
  CharPtr           str;
  ValNodePtr PNTR   vnpp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  vnpp = (ValNodePtr PNTR) mydata;
  if (vnpp == NULL) return;
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
        if (sfp->data.choice == SEQFEAT_CDREGION) {
          pseudo = sfp->pseudo;
          if (! pseudo) {
            gbqual = sfp->qual;
            while (gbqual != NULL) {
              if (StringICmp (gbqual->qual, "pseudo") == 0) {
                pseudo = TRUE;
              }
              gbqual = gbqual->next;
            }
          }
          if (pseudo) {
            product = BioseqFind (SeqLocId (sfp->product));
            if (product != NULL) {
              ValNodeAddPointer (vnpp, 0, (Pointer) product);
              sfp->product = SeqLocFree (sfp->product);
              prot = SeqMgrGetBestProteinFeature (product, NULL);
              if (prot == NULL) {
                bcp = BioseqContextNew (product);
                prot = BioseqContextGetSeqFeat (bcp, SEQFEAT_PROT, NULL, NULL, 0);
                BioseqContextFree (bcp);
              }
              if (prot != NULL) {
                prp = (ProtRefPtr) prot->data.value.ptrvalue;
                if (prp != NULL) {
                  label = NULL;
                  if (prp->name != NULL) {
                    label = prp->name->data.ptrvalue;
                  } else if (prp->desc != NULL) {
                    label = prp->desc;
                  }
                  if (label != NULL) {
                    if (sfp->comment == NULL) {
                      sfp->comment = StringSaveNoNull (label);
                    } else {
                      len = StringLen (sfp->comment) + StringLen (label) + 5;
                      str = MemNew (sizeof (Char) * len);
                      StringCpy (str, sfp->comment);
                      StringCat (str, "; ");
                      StringCat (str, label);
                      sfp->comment = MemFree (sfp->comment);
                      sfp->comment = str;
                    }
                  }
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
        case SEQID_TPE :
        case SEQID_TPD :
          *isEmblOrDdbj = TRUE;
          break;
          break;
        default :
          break;
      }
    }
  }
}

static void CleanUpPseudoProductsEx (Uint2 entityID, SeqEntryPtr sep, Boolean doPseudo)

{
  BioseqPtr      bsp;
  Char           id [41];
  Boolean        isEmblOrDdbj = FALSE;
  Uint4          itemID;
  ValNodePtr     list;
  OMProcControl  ompc;
  ValNodePtr     vnp;

  if (entityID == 0 || sep == NULL) return;
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);
  if (isEmblOrDdbj) return;
  list = NULL;
  SeqEntryExplore (sep, &list, CollectPseudoCdsProducts);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
    if (itemID > 0) {
      if (doPseudo) {
        MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
        ompc.do_not_reload_from_cache = TRUE;
        ompc.input_entityID = entityID;
        ompc.input_itemID = itemID;
        ompc.input_itemtype = OBJ_BIOSEQ;
        if (! DetachDataForProc (&ompc, FALSE)) {
          Message (MSG_POSTERR, "DetachDataForProc failed");
        }
      } else {
        SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
        ErrPostEx (SEV_WARNING, 0, 2, "Accession %s is product of pseudo CDS", id);
      }
    }
  }
  ValNodeFree (list);
}

extern void CleanUpPseudoProducts (Uint2 entityID, SeqEntryPtr sep)

{
  CleanUpPseudoProductsEx (entityID, sep, TRUE);
}

extern void CleanupGenbankCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Boolean        empty;
  GBBlockPtr     gbp;
  ValNodePtr     nextsdp;
  Pointer PNTR   prevsdp;
  ValNodePtr     sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
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
    empty = FALSE;
    if (sdp->choice == Seq_descr_genbank && sdp->data.ptrvalue != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      /* gbp->source = MemFree (gbp->source); */
      /* gbp->origin = MemFree (gbp->origin); */
      gbp->taxonomy = MemFree (gbp->taxonomy);
      if (gbp->extra_accessions == NULL && gbp->source == NULL &&
          gbp->keywords == NULL && gbp->origin == NULL &&
          gbp->date == NULL && gbp->entry_date == NULL &&
          gbp->div == NULL && gbp->taxonomy == NULL) {
        empty = TRUE;
      }
    }
    if (empty) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

static void BarCodeTechToKeyword (BioseqPtr bsp, Pointer userdata)

{
  GBBlockPtr   gbp;
  MolInfoPtr   mip;
  SeqDescrPtr  sdp;
  ValNodePtr   vnp;

  if (bsp == NULL) return;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;

  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL || mip->tech != MI_TECH_barcode) return;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_genbank, NULL);
  if (sdp == NULL) {
    gbp = GBBlockNew ();
    if (gbp != NULL) {
      sdp = SeqDescrAddPointer (&(bsp->descr), Seq_descr_genbank, (Pointer) gbp);
    }
  }
  if (sdp == NULL || sdp->choice != Seq_descr_genbank) return;

  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, "BARCODE") == 0) return;
  }

  vnp = ValNodeCopyStr (NULL, 0, "BARCODE");
  if (vnp == NULL) return;

  vnp->next = gbp->keywords;
  gbp->keywords = vnp;
}

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

extern void MergeAdjacentAnnotsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   nextsap;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
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

static Boolean HasEvidenceOrInferenceQual (SeqFeatPtr sfp)

{
  GBQualPtr  gbq;

  if (sfp == NULL) return FALSE;
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "experiment") == 0) return TRUE;
    if (StringICmp (gbq->qual, "inference") == 0) return TRUE;
  }
  return FALSE;
}

extern void CleanupEmptyFeatCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Boolean        empty;
  GeneRefPtr     grp;
  ImpFeatPtr     ifp;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  ProtRefPtr     prp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  ValNodePtr     vnp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        empty = FALSE;
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
              grp->db == NULL && grp->syn == NULL &&
              (! sfp->pseudo) && (! grp->pseudo) &&
              (sfp->exp_ev == 0) &&
              (! HasEvidenceOrInferenceQual (sfp))) {
            empty = TRUE;
          }
          if (empty) {
            /* if it has a comment, convert to misc_feature */
            if (! EmptyOrNullString (sfp->comment)) {
              ifp = (ImpFeatPtr) MemNew (sizeof (ImpFeat));
              if (ifp != NULL) {
                ifp->key = StringSave ("misc_feature");
                sfp->data.choice = SEQFEAT_IMP;
                sfp->data.value.ptrvalue = (Pointer) ifp;
                sfp->idx.subtype = 0;
                GeneRefFree (grp);
                empty = FALSE;
              }
            }
          }
        } else if (sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp->processed != 3 && prp->processed != 4 &&
              prp->name == NULL && sfp->comment != NULL) {
            if (StringICmp (sfp->comment, "putative") != 0) {
              ValNodeAddStr (&(prp->name), 0, sfp->comment);
              sfp->comment = NULL;
            }
          }
          if (prp->processed == 2 && prp->name == NULL) {
            ValNodeCopyStr (&(prp->name), 0, "unnamed");
          }
          if (prp->processed != 3 && prp->processed != 4) {
            vnp = prp->name;
            if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
                EmptyOrNullString (prp->desc) &&
                prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
              empty = TRUE;
            }
          }
        } else if (sfp->data.choice == SEQFEAT_COMMENT && EmptyOrNullString (sfp->comment)) {
          empty = TRUE;
        }
        if (empty) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          SeqFeatFree (sfp);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveBioSourceOnPopSet (SeqEntryPtr sep, OrgRefPtr master)

{
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  OrgRefPtr     orp;
  ValNodePtr    sdp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 ||
        (bssp->_class >= 13 && bssp->_class <= 16) ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_small_genome_set) { /* now on phy and mut sets */
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
      if (sdp == NULL) return;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop == NULL) return;
      orp = biop->org;
      if (orp == NULL) return;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RemoveBioSourceOnPopSet (sep, orp);
      }
      sdp = ValNodeExtract (&(bssp->descr), Seq_descr_source);
      SeqDescrFree (sdp);
      return;
    }
    /* if (bssp->_class == 7 || bssp->_class == 13 || bssp->_class == 15) return; */
    if (bssp->_class == 7) { /* also handle genbank supersets */
      orp = NULL;
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
        if (biop != NULL) {
          orp = biop->org;
        }
      }
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RemoveBioSourceOnPopSet (sep, orp);
      }
      sdp = ValNodeExtract (&(bssp->descr), Seq_descr_source);
      SeqDescrFree (sdp);
      return;
    }
  }
  if (master == NULL) return;
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

extern void RemoveMolInfoOnPopSet (SeqEntryPtr sep, MolInfoPtr master)

{
  BioseqSetPtr  bssp;
  MolInfoPtr    mip;
  ValNodePtr    sdp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 ||
        (bssp->_class >= 13 && bssp->_class <= 16) ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_small_genome_set) {
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
      if (sdp == NULL) return;
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip == NULL) return;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RemoveMolInfoOnPopSet (sep, mip);
      }
      sdp = ValNodeExtract (&(bssp->descr), Seq_descr_molinfo);
      SeqDescrFree (sdp);
      return;
    }
  }
  if (master == NULL) return;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  if (sdp != NULL) return;
  mip = MolInfoNew ();
  if (mip == NULL) return;
  mip->biomol = master->biomol;
  mip->tech = master->tech;
  if (StringDoesHaveText (master->techexp)) {
    mip->techexp = StringSave (master->techexp);
  }
  mip->completeness = master->completeness;
  if (StringDoesHaveText (master->gbmoltype)) {
    mip->gbmoltype = StringSave (master->gbmoltype);
  }
  sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = (Pointer) mip;
}

/* NoBiosourceOrTaxonId also looks for lineage and division */

extern Boolean NoBiosourceOrTaxonId (SeqEntryPtr sep)

{
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  DbtagPtr      dbt;
  Boolean       notaxid;
  ObjectIdPtr   oid;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  ValNodePtr    sdp;
  ValNodePtr    vnp;

  if (sep == NULL) return TRUE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (bssp->_class >= 13 && bssp->_class <= 16) ||
                         bssp->_class == BioseqseqSet_class_wgs_set ||
                         bssp->_class == BioseqseqSet_class_small_genome_set)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        if (NoBiosourceOrTaxonId (sep)) return TRUE;
      }
      return FALSE;
    }
  }
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
  if (sdp == NULL) return TRUE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return TRUE;
  orp = biop->org;
  if (orp == NULL) return TRUE;
  vnp = orp->db;
  if (vnp == NULL) return TRUE;
  notaxid = TRUE;
  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringCmp (dbt->db, "taxon") == 0) {
        oid = dbt->tag;
        if (oid != NULL) {
          if (oid->str == NULL && oid->id > 0) {
            notaxid = FALSE;
          }
        }
      }
    }
    vnp = vnp->next;
  }
  if (notaxid) return TRUE;
  onp = orp->orgname;
  if (onp == NULL) return TRUE;
  if (StringHasNoText (onp->lineage) || StringHasNoText (onp->div)) return TRUE;
  return FALSE;
}

static void CollectGeneFeatures (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  SeqAnnotPtr      sap;
  SeqFeatPtr       sfp;
  ValNodePtr PNTR  vnpp;

  if (sep == NULL || sep->data.ptrvalue == NULL || mydata == NULL) return;
  vnpp = (ValNodePtr PNTR) mydata;
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
        if (sfp->data.choice == SEQFEAT_GENE) {
          ValNodeAddPointer (vnpp, 0, (Pointer) sfp);
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static void ExtendGeneWithinNucProt (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 ||
        (bssp->_class >= 13 && bssp->_class <= 16) ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_small_genome_set) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        ExtendGeneWithinNucProt (sep);
      }
      return;
    }
  }
  vnp = NULL;
  SeqEntryExplore (sep, (Pointer) (&vnp), CollectGeneFeatures);
  if (vnp != NULL && vnp->next == NULL) {
    SeqEntryExplore (sep, NULL, CorrectGeneFeatLocation);
  }
  ValNodeFree (vnp);
}

extern void ExtendGeneFeatIfOnMRNA (Uint2 entityID, SeqEntryPtr sep)

{
  if (entityID < 1 && sep == NULL) return;
  if (entityID > 0 && sep == NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
  }
  if (sep == NULL) return;
  ExtendGeneWithinNucProt (sep);
}

static Boolean ConvertPubFeatDescProc (GatherObjectPtr gop)

{
  BioseqPtr    bsp;
  size_t       len;
  PubdescPtr   pdp;
  SeqDescPtr   sdp;
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp;
  SeqIdPtr     sip;
  CharPtr      str;
  ValNode      vn;

  if (gop->itemtype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gop->dataptr;
  /* look for publication features */
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) return TRUE;
  /* get bioseq by feature location */
  sip = SeqLocId (sfp->location);
  bsp = BioseqFind (sip);
  if (bsp == NULL) return TRUE;
  sip = SeqIdFindBest(bsp->id, 0);
  if (sip == NULL) return TRUE;
  vn.choice = SEQLOC_WHOLE;
  vn.extended = 0;
  vn.data.ptrvalue = (Pointer) sip;
  vn.next = NULL;
  /* is feature full length? */
  if (SeqLocCompare (sfp->location, &vn) != SLC_A_EQ_B) return TRUE;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return TRUE;
  sdp = CreateNewDescriptor (sep, Seq_descr_pub);
  if (sdp == NULL) return TRUE;
  /* move publication from feature to descriptor */
  sdp->data.ptrvalue = sfp->data.value.ptrvalue;
  sfp->data.value.ptrvalue = NULL;
  /* flag old feature for removal */
  sfp->idx.deleteme = TRUE;
  /* move or append comment to pubdesc comment */
  if (sfp->comment == NULL) return TRUE;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return TRUE;
  if (pdp->comment == NULL) {
    pdp->comment = sfp->comment;
  } else {
    len = StringLen (sfp->comment) + StringLen (pdp->comment) + 5;
    str = MemNew (sizeof (Char) * len);
    StringCpy (str, pdp->comment);
    StringCat (str, "; ");
    StringCat (str, sfp->comment);
    pdp->comment = MemFree (pdp->comment);
    pdp->comment = str;
  }
  sfp->comment = NULL;
  return TRUE;
}

extern void ConvertFullLenPubFeatToDesc (SeqEntryPtr sep)

{
  Boolean      objMgrFilter [OBJ_MAX];
  SeqEntryPtr  oldscope;

  if (sep == NULL) return;
  oldscope = SeqEntrySetScope (sep);

  MemSet ((Pointer) objMgrFilter, FALSE, sizeof (objMgrFilter));
  objMgrFilter [OBJ_SEQFEAT] = TRUE;

  GatherObjectsInEntity (0, OBJ_SEQENTRY, (Pointer) sep,
                         ConvertPubFeatDescProc, NULL, objMgrFilter);

  SeqEntrySetScope (oldscope);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
}

static Boolean GatherConvertSourceFeatDescProc (GatherObjectPtr gop)

{
  SeqFeatPtr  sfp;

  if (gop->itemtype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gop->dataptr;
  ConvertSourceFeatDescProc (sfp, NULL);
  return TRUE;
}

static void LookForTransgenic (SeqDescrPtr sdp, Pointer userdata)

{
  BioSourcePtr  biop;
  BoolPtr       is_trans_or_focus;
  SubSourcePtr  ssp;

  if (sdp == NULL || sdp->choice != Seq_descr_source) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  if (biop->is_focus) {
    is_trans_or_focus = (BoolPtr) userdata;
    *is_trans_or_focus = TRUE;
    return;
  }
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_transgenic) {
      is_trans_or_focus = (BoolPtr) userdata;
      *is_trans_or_focus = TRUE;
      return;
    }
  }
}

extern void ConvertFullLenSourceFeatToDesc (SeqEntryPtr sep)

{
  Boolean      is_trans_or_focus = FALSE;
  Boolean      objMgrFilter [OBJ_MAX];
  SeqEntryPtr  oldscope;

  if (sep == NULL) return;
  VisitDescriptorsInSep (sep, (Pointer) &is_trans_or_focus, LookForTransgenic);
  if (is_trans_or_focus) return;

  oldscope = SeqEntrySetScope (sep);

  MemSet ((Pointer) objMgrFilter, FALSE, sizeof (objMgrFilter));
  objMgrFilter [OBJ_SEQFEAT] = TRUE;

  GatherObjectsInEntity (0, OBJ_SEQENTRY, (Pointer) sep,
                         GatherConvertSourceFeatDescProc, NULL, objMgrFilter);

  SeqEntrySetScope (oldscope);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
}

static Int4 LoopSeqEntryToAsn3 (
  SeqEntryPtr sep,
  Boolean strip,
  Boolean correct,
  SeqEntryFunc taxfun,
  SeqEntryFunc taxmerge,
  Boolean gpipeMode,
  Boolean isEmblOrDdbj
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Boolean       taxserver;

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (bssp->_class >= 13 && bssp->_class <= 16) ||
                         bssp->_class == BioseqseqSet_class_wgs_set ||
                         bssp->_class == BioseqseqSet_class_small_genome_set)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += LoopSeqEntryToAsn3 (sep, strip, correct, taxfun, taxmerge, gpipeMode, isEmblOrDdbj);
      }
      return rsult;
    }
  }
  oldscope = SeqEntrySetScope (sep);
  taxserver = (Boolean) (taxfun != NULL || taxmerge != NULL);
  rsult = SeqEntryToAsn3Ex (sep, strip, correct, taxserver, taxfun, taxmerge, gpipeMode, isEmblOrDdbj);
  SeqEntrySetScope (oldscope);
  return rsult;
}

static Boolean DeleteBadMarkedGeneXrefs (GatherObjectPtr gop)

{
  GeneRefPtr           grp;
  SeqFeatXrefPtr       nextxref;
  SeqFeatXrefPtr PNTR  prevxref;
  SeqFeatPtr           sfp;
  SeqFeatPtr           sfpx;
  Boolean              unlink;
  SeqFeatXrefPtr       xref;

  if (gop->itemtype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gop->dataptr;
  xref = sfp->xref;
  prevxref = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  while (xref != NULL) {
    nextxref = xref->next;
    unlink = FALSE;
    if (xref->specialCleanupFlag && xref->data.choice == SEQFEAT_GENE) {
      grp = (GeneRefPtr) xref->data.value.ptrvalue;
      if (grp != NULL) {
        sfpx = SeqMgrGetOverlappingGene (sfp->location, NULL);
        if (sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE) {
          unlink = TRUE;
        }
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
  return TRUE;
}

static Boolean MarkMovedGeneGbquals (GatherObjectPtr gop)

{
  GBQualPtr       gbq;
  GeneRefPtr      grp;
  BoolPtr         hasMarkedGenesP;
  GBQualPtr       nextqual;
  GBQualPtr PNTR  prevqual;
  SeqFeatPtr      sfp;
  SeqFeatXrefPtr  xref;

  if (gop->itemtype != OBJ_SEQFEAT) return TRUE;
  hasMarkedGenesP = (BoolPtr) gop->userdata;
  if (hasMarkedGenesP == NULL) return TRUE;
  sfp = (SeqFeatPtr) gop->dataptr;
  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    TrimSpacesAroundString (gbq->qual);
    TrimSpacesAroundString (gbq->val);
    nextqual = gbq->next;
    if (StringICmp (gbq->qual, "gene") == 0 && (! StringHasNoText (gbq->val))) {
      grp = GeneRefNew ();
      grp->locus = StringSave (gbq->val);
      xref = SeqFeatXrefNew ();
      xref->data.choice = SEQFEAT_GENE;
      xref->data.value.ptrvalue = (Pointer) grp;
      xref->specialCleanupFlag = TRUE; /* flag to test for overlapping gene later */
      xref->next = sfp->xref;
      sfp->xref = xref;
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      gbq->qual = MemFree (gbq->qual);
      gbq->val = MemFree (gbq->val);
      GBQualFree (gbq);
      *hasMarkedGenesP = TRUE;
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
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
    if (sdp->choice != Seq_descr_title) continue;
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

static CitSubPtr GetCspFromPdp (PubdescPtr pdp)

{
  ValNodePtr  vnp;

  if (pdp == NULL) return NULL;
  vnp = pdp->pub;
  if (vnp == NULL) return NULL;
  if (vnp->choice != PUB_Sub) return NULL;
  return (CitSubPtr) vnp->data.ptrvalue;
}

static Boolean CitSubsMatch (CitSubPtr csp1, CitSubPtr csp2)

{
  AffilPtr     afp1, afp2;
  AuthListPtr  alp1, alp2;

  if (csp1 == NULL || csp2 == NULL) return FALSE;
  if (DateMatch (csp1->date, csp2->date, FALSE) != 0) return FALSE;
  if (StringICmp (csp1->descr, csp2->descr) != 0) return FALSE;
  alp1 = csp1->authors;
  alp2 = csp2->authors;
  if (alp1 == NULL || alp2 == NULL) return FALSE;
  if (AuthListMatch (alp1, alp2, TRUE) != 0) return FALSE;
  afp1 = alp1->affil;
  afp2 = alp2->affil;
  if (afp1 != NULL && afp2 != NULL) {
    if (! AsnIoMemComp (afp1, afp2, (AsnWriteFunc) AffilAsnWrite)) return FALSE;
  }
  return TRUE;
}

static void MergeEquivCitSubs (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AuthListPtr    alp1, alp2;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  CitSubPtr      csp1, csp2;
  SeqDescrPtr    descr = NULL;
  SeqDescrPtr    lastcit;
  PubdescPtr     lastpdp;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;
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

  lastcit = NULL;
  lastpdp = NULL;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_pub) continue;
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp == NULL) continue;
    if (lastcit != NULL && lastpdp != NULL) {
      csp2 = GetCspFromPdp (pdp);
      if (csp2 != NULL) {
        if (CitSubsMatch (csp1, csp2) &&
            OkayToFuseRemarks (pdp->comment, lastpdp->comment)) {
          alp1 = csp1->authors;
          alp2 = csp2->authors;
          if (alp1 != NULL && alp2 != NULL) {
            if (alp1->affil == NULL && alp2->affil != NULL) {
              alp1->affil = alp2->affil;
              alp2->affil = NULL;
            }
          }
          if (lastpdp->comment == NULL && pdp->comment != NULL) {
            lastpdp->comment = pdp->comment;
            pdp->comment = NULL;
          }
          if (sdp->extended != 0) {
            ovp = (ObjValNodePtr) sdp;
            ovp->idx.deleteme = TRUE;
          }
        } else {
          lastcit = sdp;
          csp1 = csp2;
        }
      }
    } else {
      csp1 = GetCspFromPdp (pdp);
      if (csp1 != NULL) {
        lastcit = sdp;
      }
    }
    lastpdp = pdp;
  }
}

static void MergeMultipleDates (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  DatePtr        dp1, dp2;
  SeqDescrPtr    descr = NULL;
  SeqDescrPtr    lastdate;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;
  Int2           status;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    descr = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    descr = bssp->descr;
  } else return;

  lastdate = NULL;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_create_date) continue;
    if (lastdate != NULL) {
      dp1 = (DatePtr) lastdate->data.ptrvalue;
      dp2 = (DatePtr) sdp->data.ptrvalue;
      status = DateMatch (dp1, dp2, FALSE);
      if (status == 1) {
        if (sdp->extended != 0) {
          ovp = (ObjValNodePtr) sdp;
          ovp->idx.deleteme = TRUE;
        }
      } else {
        if (lastdate->extended != 0) {
          ovp = (ObjValNodePtr) lastdate;
          ovp->idx.deleteme = TRUE;
        }
        lastdate = sdp;
      }
    } else {
      lastdate = sdp;
    }
  }

  lastdate = NULL;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_update_date) continue;
    if (lastdate != NULL) {
      dp1 = (DatePtr) lastdate->data.ptrvalue;
      dp2 = (DatePtr) sdp->data.ptrvalue;
      status = DateMatch (dp1, dp2, FALSE);
      if (status == 1) {
        if (sdp->extended != 0) {
          ovp = (ObjValNodePtr) sdp;
          ovp->idx.deleteme = TRUE;
        }
      } else {
        if (lastdate->extended != 0) {
          ovp = (ObjValNodePtr) lastdate;
          ovp->idx.deleteme = TRUE;
        }
        lastdate = sdp;
      }
    } else {
      lastdate = sdp;
    }
  }
}

static void MolInfoUpdate (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqDescrPtr    descr = NULL;
  Uint1          meth;
  MolInfoPtr     mip = NULL;
  Uint1          mol;
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
    if (sdp->choice != Seq_descr_molinfo) continue;
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }
  if (mip == NULL) return;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    switch (sdp->choice) {
      case Seq_descr_mol_type :
        mol = sdp->data.intvalue;
        if (mol != 0 && mip->biomol == 0) {
          switch (mol) {
            case MOLECULE_TYPE_GENOMIC :
              mip->biomol = MOLECULE_TYPE_GENOMIC;
              break;
            case MOLECULE_TYPE_PRE_MRNA :
              mip->biomol = MOLECULE_TYPE_PRE_MRNA;
              break;
            case MOLECULE_TYPE_MRNA :
              mip->biomol = MOLECULE_TYPE_MRNA;
              break;
            case MOLECULE_TYPE_RRNA :
              mip->biomol = MOLECULE_TYPE_RRNA;
              break;
            case MOLECULE_TYPE_TRNA :
              mip->biomol = MOLECULE_TYPE_TRNA;
              break;
            case MOLECULE_TYPE_SNRNA :
              mip->biomol = MOLECULE_TYPE_SNRNA;
              break;
            case MOLECULE_TYPE_SCRNA :
              mip->biomol = MOLECULE_TYPE_SCRNA;
              break;
            case MOLECULE_TYPE_PEPTIDE :
              mip->biomol = MOLECULE_TYPE_PEPTIDE;
              break;
            case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL :
              mip->biomol = MOLECULE_TYPE_OTHER_GENETIC_MATERIAL;
              break;
            case MOLECULE_TYPE_GENOMIC_MRNA_MIX :
              mip->biomol = MOLECULE_TYPE_GENOMIC_MRNA_MIX;
              break;
            case 255 :
              mip->biomol = 255;
              break;
            default :
              break;
          }
        }
        if (sdp->extended != 0) {
          ovp = (ObjValNodePtr) sdp;
          ovp->idx.deleteme = TRUE;
        }
        break;
      case Seq_descr_modif :
        break;
      case Seq_descr_method :
        meth = sdp->data.intvalue;
        if (meth != 0 && mip->tech == 0) {
          switch (meth) {
            case METHOD_concept_transl :
              mip->tech = MI_TECH_concept_trans;
              break;
            case METHOD_seq_pept :
              mip->tech = MI_TECH_seq_pept;
              break;
            case METHOD_both :
              mip->tech = MI_TECH_both;
              break;
            case METHOD_seq_pept_overlap :
              mip->tech = MI_TECH_seq_pept_overlap;
              break;
            case METHOD_seq_pept_homol :
              mip->tech = MI_TECH_seq_pept_homol;
              break;
            case METHOD_concept_transl_a :
              mip->tech = MI_TECH_concept_trans_a;
              break;
            case METHOD_other :
              mip->tech = MI_TECH_other;
              break;
            default :
              break;
          }
        }
        if (sdp->extended != 0) {
          ovp = (ObjValNodePtr) sdp;
          ovp->idx.deleteme = TRUE;
        }
        break;
      default :
        break;
    }
  }
}

static AuthListPtr SSECGetAuthListPtr (PubdescPtr pdp)

{
  AuthListPtr  alp = NULL;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  ValNodePtr   vnp;

  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          alp = cgp->authors;
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          alp = csp->authors;
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          alp = cap->authors;
        }
        break;
      case PUB_Book :
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          alp = cbp->authors;
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          alp = cpp->authors;
        }
        break;
      default :
        break;
    }

    if (alp != NULL) return alp;
  }

  return NULL;
}

static Boolean JustMuid (ValNodePtr ppr)

{
  ValNodePtr  vnp;

  if (ppr == NULL) return FALSE;
  if (ppr->choice == PUB_Muid) return TRUE;
  if (ppr->choice == PUB_Equiv) {
    for (vnp = (ValNodePtr) ppr->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      if (JustMuid (vnp)) return TRUE;
    }
  }
  return FALSE;
}

static void FixZeroMuid (ValNodePtr ppr, Int4 muid)

{
  ValNodePtr  vnp;

  if (ppr == NULL) return;
  if (ppr->choice == PUB_Muid) {
    ppr->data.intvalue = muid;
  }
  if (ppr->choice == PUB_Equiv) {
    for (vnp = (ValNodePtr) ppr->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      FixZeroMuid (vnp, muid);
    }
  }
}

static void RepairBadBackbonePub (PubdescPtr pdp, Int4 muid)

{
  ValNodePtr  vnp;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (! JustMuid (vnp)) return;
  }
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    FixZeroMuid (vnp, muid);
  }
}

static void RemoveZeroMuids (ValNodePtr  ppr, ValNodePtr PNTR prev)

{
  ValNodePtr  next;

  /* if only muid 0, cannot leave empty pdp->pub */

  if (ppr == NULL || ppr->next == NULL) return;

  while (ppr != NULL) {
    next = ppr->next;
    if (ppr->choice == PUB_Muid && ppr->data.intvalue == 0) {
      *prev = ppr->next;
      ppr->next = NULL;
      ValNodeFree (ppr);
    } else {
      prev = (ValNodePtr PNTR) &(ppr->next);
    }
    ppr = next;
  }
}

static Boolean IsPatent (PubdescPtr pdp)

{
  ValNodePtr  vnp;

  if (pdp == NULL) return FALSE;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Patent) return TRUE;
  }
  return FALSE;
}

extern Boolean IsPubContentBad (PubdescPtr pdp)
{
  AuthListPtr  alp = NULL;
  AuthorPtr    ap;
  CitArtPtr    cap;
  CitJourPtr   cjp;
  ImprintPtr   imp;
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      title = NULL;
  ValNodePtr   ttl, vnp;

  if (pdp == NULL) return FALSE;
  /* keep anything with a figure - backbone entry */

  if (! StringHasNoText (pdp->fig)) return FALSE;

  /* look for at least one author name */

  alp = SSECGetAuthListPtr (pdp);
  if (alp == NULL) return TRUE;
  if (IsPatent (pdp)) {
    /* patents can get away with no authors */
  } else if (alp->choice == 1) {
    names = alp->names;
    if (names == NULL) return TRUE;
    ap = (AuthorPtr) names->data.ptrvalue;
    if (ap == NULL) return TRUE;
    pid = ap->name;
    if (pid == NULL) return TRUE;
    if (pid->choice == 2) {
      nsp = (NameStdPtr) pid->data;
      if (nsp == NULL) return TRUE;
      if (StringHasNoText (nsp->names [0])) return TRUE;
    } else if (pid->choice == 3 || pid->choice == 4) {
      if (StringHasNoText ((CharPtr) pid->data)) return TRUE;
    }
  } else if (alp->choice == 2 || alp->choice == 3) {
    names = alp->names;
    if (names == NULL) return TRUE;
    if (StringHasNoText ((CharPtr) names->data.ptrvalue)) return TRUE;
  }

  /* look for CitArt journal */

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          switch (cap->from) {
            case 1 :
              cjp = (CitJourPtr) cap->fromptr;
              if (cjp == NULL) return TRUE;
              if (cjp->title == NULL) return TRUE;
              for (ttl = cjp->title; ttl != NULL; ttl = ttl->next) {
                if (! StringHasNoText ((CharPtr) ttl->data.ptrvalue)) {
                  title = (CharPtr) ttl->data.ptrvalue;
                }
              }
              if (title == NULL) return TRUE;
              imp = cjp->imp;
              if (imp == NULL) return TRUE;
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

  return FALSE;
}

static Boolean IsPubBad (PubdescPtr pdp, Pointer userdata)

{
  Int4        muid;
  Int4Ptr     muidp;
  ValNodePtr  vnp;

  if (pdp == NULL) return FALSE;

  /* single pmid not cleared here, left for CheckMinPub with RefSeq protein exception */

  vnp = pdp->pub;
  if (vnp != NULL && vnp->next == NULL && vnp->choice == PUB_PMid) {
    /* but first check for 0 pmid, mark for removal */
    if (vnp->data.intvalue == 0) return TRUE;
    return FALSE;
  }

  /* if single real muid, repair 0 muid backbone references */

  muidp = (Int4Ptr) userdata;
  if (muidp != NULL) {
    muid = *muidp;
    if (muid != 0 && muid != -1) {
      RepairBadBackbonePub (pdp, muid);
    }
  }

  /* remove remaining 0 muids */

  RemoveZeroMuids (pdp->pub, &(pdp->pub));

  return IsPubContentBad (pdp);
}

static void RemoveBadPubFeat (SeqFeatPtr sfp, Pointer userdata)

{
  PubdescPtr  pdp;

  if (sfp->data.choice != SEQFEAT_PUB) return;
  pdp = (PubdescPtr) sfp->data.value.ptrvalue;
  if (IsPubBad (pdp, userdata)) {
    sfp->idx.deleteme = TRUE;
  }
}

static void RemoveBadPubDescr (SeqDescrPtr sdp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;

  if (sdp->choice != Seq_descr_pub) return;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (IsPubBad (pdp, userdata)) {
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}

static void LookForUniqMuidProc (ValNodePtr ppr, Int4Ptr muidp)

{
  Int4        muid;
  ValNodePtr  vnp;

  if (ppr == NULL || muidp == NULL) return;
  if (*muidp == -1) return;
  switch (ppr->choice) {
    case PUB_Muid :
      muid = ppr->data.intvalue;
      if (muid == 0) {
      } else if (*muidp == 0) {
        *muidp = muid;
      } else if (*muidp != muid) {
        *muidp = -1;
      }
      break;
    case PUB_Equiv :
      for (vnp = (ValNodePtr) ppr->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
        LookForUniqMuidProc (vnp, muidp);
      }
      break;
    default :
      break;
  }
}

static void LookForUniqueMuid (PubdescPtr pdp, Pointer userdata)

{
  Int4Ptr     muidp;
  ValNodePtr  vnp;

  muidp = (Int4Ptr) userdata;
  if (*muidp == -1) return;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    LookForUniqMuidProc (vnp, muidp);
  }
}

static Boolean EmptyExceptSerialNumber (CitGenPtr cgp)
{
  if (cgp == NULL) return FALSE;

  if (!StringHasNoText (cgp->cit)
      || cgp->authors != NULL
      || cgp->muid != 0
      || cgp->journal != NULL
      || !StringHasNoText (cgp->volume)
      || !StringHasNoText (cgp->issue)
      || !StringHasNoText (cgp->pages)
      || cgp->date != NULL
      || !StringHasNoText (cgp->title)
      || cgp->pmid != 0) {
    return FALSE;
  } else {
    return TRUE;
  }
}

static Boolean ArePubsMergeableForFig(PubdescPtr fig, PubdescPtr nofig)
{
  CitGenPtr  cgp_fig, cgp_nofig;

  if (fig == NULL || nofig == NULL) return FALSE;

  /* name */
  if (!StringHasNoText (nofig->name) && !StringHasNoText (fig->name)
      && !StringCmp (nofig->name, fig->name)) {
    return FALSE;
  }
  /* fig */
  if (!StringHasNoText (nofig->fig) && !StringHasNoText (fig->fig)
      && !StringCmp (nofig->fig, fig->fig)) {
    return FALSE;
  }
  /* num */
  if (NumberingMatch(nofig->num, fig->num) != 0) {
    return FALSE;
  }

  if ((nofig->numexc && !fig->numexc)
      || (!nofig->numexc && fig->numexc)) {
    return FALSE;
  }
  if ((nofig->poly_a && !fig->poly_a)
      || (!nofig->poly_a && fig->poly_a)) {
    return FALSE;
  }

  /* maploc */
  if (!StringHasNoText (nofig->maploc) && !StringHasNoText (fig->maploc)
      && !StringCmp (nofig->maploc, fig->maploc)) {
    return FALSE;
  }

  /* seq-raw */
  if (!StringHasNoText (nofig->seq_raw) && !StringHasNoText (fig->seq_raw)
      && !StringCmp (nofig->seq_raw, fig->seq_raw)) {
    return FALSE;
  }

  /* align-group */
  if (nofig->align_group > 0 && fig->align_group > 0
      && nofig->align_group != fig->align_group) {
    return FALSE;
  }

  /* comment */
  if (!StringHasNoText (nofig->comment) && !StringHasNoText (fig->comment)
      && !StringCmp (nofig->comment, fig->comment)) {
    return FALSE;
  }

  /* reftype */
  if (nofig->reftype > 0 && fig->reftype > 0
      && nofig->reftype != fig->reftype) {
    return FALSE;
  }

  if (nofig->pub != NULL && fig->pub != NULL) {
    if (nofig->pub->next != NULL
        || fig->pub->next != NULL
        || nofig->pub->choice != PUB_Gen
        || fig->pub->choice != PUB_Gen
        || nofig->pub->data.ptrvalue == NULL
        || fig->pub->data.ptrvalue == NULL) {
      return FALSE;
    }
    cgp_fig = fig->pub->data.ptrvalue;
    cgp_nofig = nofig->pub->data.ptrvalue;
    if (!EmptyExceptSerialNumber (cgp_fig)
        || !EmptyExceptSerialNumber (cgp_nofig)
        || !cgp_fig->serial_number != cgp_nofig->serial_number) {
      return FALSE;
    }
  }

  return TRUE;
}

/* rescue pub with just fig that is in same chain as real pub by merging data */

static void MergePubFigInChain (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqDescrPtr    descr = NULL;
  PubdescPtr     hasfig = NULL;
  Int4Ptr        muidp;
  PubdescPtr     nofig = NULL;
  ObjValNodePtr  ovp = NULL;
  PubdescPtr     pdp;
  SeqDescrPtr    sdp_fig, sdp_nofig;

  muidp = (Int4Ptr) mydata;
  if (*muidp == -1) return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    descr = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    descr = bssp->descr;
  } else return;

  /* attempt to combine any descriptor with a fig with
   * a different publication without a fig, but only if
   * they do not have conflicting data.
   */
  for (sdp_fig = descr; sdp_fig != NULL; sdp_fig = sdp_fig->next) {
    if (sdp_fig->choice != Seq_descr_pub
        || sdp_fig->extended == 0) continue;
    pdp = (PubdescPtr) sdp_fig->data.ptrvalue;
    if (pdp == NULL) continue;
    if (! StringHasNoText (pdp->fig)) {
      hasfig = pdp;
      ovp = (ObjValNodePtr) sdp_fig;

      for (sdp_nofig = descr; sdp_nofig != NULL; sdp_nofig = sdp_nofig->next) {
        if (sdp_nofig->choice != Seq_descr_pub) continue;
        pdp = (PubdescPtr) sdp_nofig->data.ptrvalue;
        if (pdp != NULL && StringHasNoText (pdp->fig)
            && ArePubsMergeableForFig(hasfig, pdp)) {
          nofig = pdp;
          if (StringHasNoText (nofig->name)) {
            nofig->name = MemFree (nofig->name);
            nofig->name = hasfig->name;
            hasfig->name = NULL;
          }
          if (StringHasNoText (nofig->fig)) {
            nofig->fig = MemFree (nofig->fig);
            nofig->fig = hasfig->fig;
            hasfig->fig = NULL;
          }
          if (nofig->num == NULL) {
            nofig->num = hasfig->num;
            hasfig->num = NULL;
          }
          if (hasfig->numexc) {
            nofig->numexc = hasfig->numexc;
          }
          if (hasfig->poly_a) {
            nofig->poly_a = hasfig->poly_a;
          }
          if (hasfig->align_group > 0) {
            nofig->align_group = hasfig->align_group;
          }
          if (StringHasNoText (nofig->maploc)) {
            nofig->maploc = MemFree (nofig->maploc);
            nofig->maploc = hasfig->maploc;
            hasfig->maploc = NULL;
          }
          if (StringHasNoText (nofig->seq_raw)) {
            nofig->seq_raw = MemFree (nofig->seq_raw);
            nofig->seq_raw = hasfig->seq_raw;
            hasfig->seq_raw = NULL;
          }
          if (StringHasNoText (nofig->comment)) {
            nofig->comment = MemFree (nofig->comment);
            nofig->comment = hasfig->comment;
            hasfig->comment = NULL;
          }
          if (hasfig->reftype > 0) {
            nofig->reftype = hasfig->reftype;
          }

          if (nofig->pub == NULL) {
            nofig->pub = hasfig->pub;
          }

          ovp->idx.deleteme = TRUE;
        }
      }
    }
  }
}

static void CorrectSfpExceptText (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL || StringHasNoText (sfp->except_text)) return;
  if (StringICmp (sfp->except_text, "reasons cited in publication") == 0) {
    sfp->except_text = MemFree (sfp->except_text);
    sfp->except_text = StringSave ("reasons given in citation");
  }
}

static Boolean IsCodonCorrect (tRNAPtr trp, Uint1 taa)

{
  Uint1           aa;
  Uint1           from;
  SeqMapTablePtr  smtp;

  if (trp == NULL) return TRUE;
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
  if (aa > 0 && aa != 255) {
    if (taa != aa && aa != 'U') {
      return FALSE;
    }
  }
  return TRUE;
}

static    Uint1 codon_xref [4] = {   /* mapping from NCBI2na to codon codes */
2,  /* A */
1,  /* C */
3,  /* G */
0 }; /* T */

static void CorrectTrnaCodons (SeqFeatPtr sfp, Pointer userdata)

{
  Uint1           alt [4];
  Char            ch;
  CharPtr         codes;
  Uint1           codon [4];
  Int2            i, j, k;
  Uint1           index;
  RnaRefPtr       rrp;
  Uint1           residue;
  SeqMapTablePtr  smtp;
  Uint1           taa;
  tRNAPtr         trp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->type != 3 || rrp->ext.choice != 2) return;
  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) return;
  if (trp->codon [0] == 255 || trp->codon [1] != 255 || trp->codon [0] >= 64) return;
  codes = (CharPtr) userdata;
  if (codes == NULL) return;

  /* see if existing codon is correct in current genetic code */

  index = trp->codon [0];
  taa = codes [index];
  if (IsCodonCorrect (trp, taa)) return;

  /* convert codon to string */

  smtp = SeqMapTableFind (Seq_code_iupacna, Seq_code_ncbi2na);
  if (smtp == NULL) return;
  for (i = 0, j = 16; i < 3; i++, j /= 4) {
    residue = (Uint1) ((Int2) index / j);
    index -= (Uint1) (residue * j);
    for (k = 0; k < 4; k++) {
      if (codon_xref[k] == residue) {
        residue = (Uint1) k;
        break;
      }
    }
    residue = SeqMapTableConvert (smtp, residue);
    codon [i] = residue;
  }
  codon [3] = 0;

  /* reverse complement */

  for (i = 0; i < 3; i++) {
    ch = (Char) codon [2 - i];
    if (ch == 'A') {
      ch = 'T';
    } else if (ch == 'C') {
      ch = 'G';
    } else if (ch == 'G') {
      ch = 'C';
    } else if (ch == 'T') {
      ch = 'A';
    }
    alt [i] = ch;
  }
  alt [3] = 0;

  /* see if revcomp codon is correct in current genetic code */

  index = IndexForCodon (alt, Seq_code_iupacna);
  taa = codes [index];
  if (IsCodonCorrect (trp, taa)) {
    trp->codon [0] = index;
    return;
  }

  /* just complement */

  for (i = 0; i < 3; i++) {
    ch = (Char) codon [i];
    if (ch == 'A') {
      ch = 'T';
    } else if (ch == 'C') {
      ch = 'G';
    } else if (ch == 'G') {
      ch = 'C';
    } else if (ch == 'T') {
      ch = 'A';
    }
    alt [i] = ch;
  }
  alt [3] = 0;

  /* see if complement codon is correct in current genetic code */

  index = IndexForCodon (alt, Seq_code_iupacna);
  taa = codes [index];
  if (IsCodonCorrect (trp, taa)) {
    trp->codon [0] = index;
    return;
  }

  /* just reverse */

  for (i = 0; i < 3; i++) {
    ch = (Char) codon [2 - i];
    alt [i] = ch;
  }
  alt [3] = 0;

  /* see if reverse codon is correct in current genetic code */

  index = IndexForCodon (alt, Seq_code_iupacna);
  taa = codes [index];
  if (IsCodonCorrect (trp, taa)) {
    trp->codon [0] = index;
    return;
  }
}

static void FindSingleBioSource (BioSourcePtr biop, Pointer userdata)

{
  BioSourcePtr PNTR  biopp;

  biopp = (BioSourcePtr PNTR) userdata;
  if (biop == NULL || biopp == NULL) return;
  *biopp = biop;
}

static void CleanupOldName (BioSourcePtr biop, Pointer userdata)

{
  OrgModPtr       next, omp;
  OrgNamePtr      onp;
  OrgRefPtr       orp;
  OrgModPtr PNTR  prev;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL || StringHasNoText (orp->taxname)) return;
  onp = orp->orgname;
  if (onp == NULL || onp->mod == NULL) return;
  prev = &(onp->mod);
  omp = *prev;
  while (omp != NULL) {
    next = omp->next;
    if (omp->subtype == ORGMOD_old_name && StringCmp (orp->taxname, omp->subname) == 0 && StringHasNoText (omp->attrib)) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      prev = &(omp->next);
    }
    omp = next;
  }
}

static void CleanupOrgModNote (BioSourcePtr biop, Pointer userdata)

{
  CharPtr         gbacr = NULL, gbana = NULL, gbsyn = NULL, taxname = NULL;
  OrgModPtr       next, omp;
  OrgNamePtr      onp;
  OrgRefPtr       orp;
  OrgModPtr PNTR  prev;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  taxname = orp->taxname;
  onp = orp->orgname;
  if (onp == NULL || onp->mod == NULL) return;
  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == ORGMOD_gb_acronym) {
      gbacr = omp->subname;
    } else if (omp->subtype == ORGMOD_gb_anamorph) {
      gbana = omp->subname;
    } else if (omp->subtype == ORGMOD_gb_synonym) {
      gbsyn = omp->subname;
    }
  }
  prev = &(onp->mod);
  omp = *prev;
  while (omp != NULL) {
    next = omp->next;
    if (omp->subtype == ORGMOD_other &&
        StringDoesHaveText (omp->subname) &&
        (StringCmp (taxname, omp->subname) == 0 ||
         StringCmp (gbacr, omp->subname) == 0 ||
         StringCmp (gbana, omp->subname) == 0 ||
         StringCmp (gbsyn, omp->subname) == 0)) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      prev = &(omp->next);
    }
    omp = next;
  }
}

static SeqFeatPtr GetUnambigOverlappingGene (BioseqPtr bsp, SeqLocPtr slp)

{
  SeqMgrFeatContext  context;
  SeqFeatPtr         gene;
  Int2               i;
  Int4Ptr            ivals;
  Int2               j;
  SeqFeatPtr         next;
  Int2               numivals;

  gene = SeqMgrGetOverlappingGene (slp, &context);
  if (gene == NULL) return NULL;
  numivals = context.numivals;
  ivals = context.ivals;
  next = SeqMgrGetNextFeature (bsp, gene, SEQFEAT_GENE, 0, &context);
  if (next == NULL) return gene;
  if (numivals != context.numivals) return gene;
  for (i = 0, j = 0; i < numivals; i++) {
    if (ivals [j] != context.ivals [j]) return gene;
    j++;
    if (ivals [j] != context.ivals [j]) return gene;
    j++;
  }
  return NULL;
}

static void LookForPeptides (SeqFeatPtr sfp, Pointer userdata)

{
  BoolPtr     hasPepP;
  ProtRefPtr  prp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  hasPepP = (BoolPtr) userdata;
  if (hasPepP == NULL) return;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;
  if (prp->processed > 0) {
    *hasPepP = TRUE;
  }
}

static void RemoveUnneededGeneXrefs (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr            bp;
  SeqFeatXrefPtr     curr;
  SeqMgrFeatContext  fcontext;
  GeneRefPtr         grp;
  GeneRefPtr         grpx;
  Boolean            hasPeptide;
  Boolean            isEmblDdbj = FALSE;
  SeqFeatXrefPtr     PNTR last;
  SeqFeatXrefPtr     next;
  BioseqPtr          pbsp;
  Boolean            redundantgenexref;
  SeqFeatPtr         sfp;
  SeqFeatPtr         sfpx;
  CharPtr            syn1;
  CharPtr            syn2;

  bp = (BoolPtr) userdata;
  if (bp != NULL) {
    isEmblDdbj = *bp;
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->data.choice != SEQFEAT_GENE) {
      grp = SeqMgrGetGeneXref (sfp);
      if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
        sfpx = GetUnambigOverlappingGene (bsp, sfp->location);
        if (sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE) {
          grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
          if (grpx != NULL) {
            redundantgenexref = FALSE;
            if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
              if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
                redundantgenexref = TRUE;
              }
            } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
              if (StringICmp (grp->locus, grpx->locus) == 0) {
                redundantgenexref = TRUE;
              }
            } else if (grp->syn != NULL && grpx->syn != NULL) {
              syn1 = (CharPtr) grp->syn->data.ptrvalue;
              syn2 = (CharPtr) grpx->syn->data.ptrvalue;
              if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
                if (StringICmp (syn1, syn2) == 0) {
                  redundantgenexref = TRUE;
                }
              }
            }
            if (redundantgenexref && isEmblDdbj && sfp->data.choice == SEQFEAT_CDREGION) {
              hasPeptide = FALSE;
              pbsp = BioseqFindFromSeqLoc (sfp->product);
              if (pbsp != NULL) {
                VisitFeaturesOnBsp (pbsp, (Pointer) &hasPeptide, LookForPeptides);
                if (hasPeptide) {
                  redundantgenexref = FALSE;
                }
              }
            }
            if (redundantgenexref) {
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
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}

static void StripBadTitleFromProteinProducts (BioseqPtr bsp, Pointer userdata)

{
  BioseqSetPtr       bssp;
  CharPtr            buf;
  size_t             buflen = 1001;
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  ObjValNodePtr      ovp;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  Uint1              tech;
  CharPtr            title;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) return;
  }
  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;

  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    while (bssp != NULL && bssp->_class != BioseqseqSet_class_nuc_prot) {
      if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) bssp->idx.parentptr;
      } else {
        bssp = NULL;
      }
    }
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_nuc_prot) {
      title = (CharPtr) vnp->data.ptrvalue;
      tech = 0;
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
      if (sdp != NULL) {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip != NULL) {
          tech = mip->tech;
        }
      }
      buf = MemNew (sizeof (Char) * (buflen + 1));
      if (buf != NULL && NewCreateDefLineBuf (NULL, bsp, buf, buflen, TRUE, FALSE)) {
        if (StringICmp (buf, title) != 0) {
          if (vnp->extended != 0) {
            ovp = (ObjValNodePtr) vnp;
            ovp->idx.deleteme = TRUE;
          }
        }
      }
      MemFree (buf);
    }
  }
}

static void MarkBadProtTitlesInNucProts (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class == BioseqseqSet_class_genbank ||
      (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) ||
      bssp->_class == BioseqseqSet_class_wgs_set ||
      bssp->_class == BioseqseqSet_class_small_genome_set) {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      /* StripTitleFromProtsInNucProts (sep); */
      MarkBadProtTitlesInNucProts (sep);
    }
    return;
  }
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;
  VisitBioseqsInSep (sep, NULL, StripBadTitleFromProteinProducts);
}

static void MakeNcbiCleanupObject (SeqEntryPtr sep, Boolean gpipeMode)

{
  DatePtr        dp;
  ValNodePtr     sdp;
  UserObjectPtr  uop;

  dp = DateCurr ();
  if (dp == NULL) return;

  uop = CreateNcbiCleanupUserObject ();
  if (uop == NULL) return;

  if (gpipeMode) {
    AddStringToNcbiCleanupUserObject (uop, "method", "GpipeSeqEntryCleanup");
  } else {
    AddStringToNcbiCleanupUserObject (uop, "method", "SeriousSeqEntryCleanup");
  }
  AddIntegerToNcbiCleanupUserObject (uop, "version", NCBI_CLEANUP_VERSION);

  AddIntegerToNcbiCleanupUserObject (uop, "month", dp->data [2]);
  AddIntegerToNcbiCleanupUserObject (uop, "day", dp->data [3]);
  AddIntegerToNcbiCleanupUserObject (uop, "year", dp->data [1] + 1900);

  sdp = NewDescrOnSeqEntry (sep, Seq_descr_user);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = uop;
}

static void SeriousSeqEntryCleanupEx (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge, Boolean doPseudo, Boolean gpipeMode)

{
  BioSourcePtr    biop;
  BioseqSetPtr    bssp;
  Int2            code;
  CharPtr         codes;
  Uint2           entityID;
  GeneticCodePtr  gncp;
  Boolean         hasMarkedGenes = FALSE;
  Boolean         isEmblOrDdbj = FALSE;
  ErrSev          lsev;
  ErrSev          msev;
  Int4            muid = 0;
  Boolean         objMgrFilter [OBJ_MAX];
  SeqEntryPtr     oldscope;
  SeqEntryPtr     tmp;
  ValNodePtr      vnp;

  if (sep == NULL) return;
  oldscope = SeqEntrySetScope (sep);
  msev = ErrSetMessageLevel (SEV_MAX);
  lsev = ErrSetLogLevel (SEV_MAX);
  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  /* clear indexes, since CleanupEmptyFeatCallback removes genes, etc. */
  SeqMgrClearFeatureIndexes (entityID, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);
  RemoveDuplicateNestedSetsForEntityID (entityID);
  SeqMgrClearFeatureIndexes (entityID, NULL);
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
      tmp = bssp->seq_set;
      if (tmp != NULL && tmp->next == NULL && (IS_Bioseq (tmp))) {
        /* coerce genbank set on top of single sequence to nuc-prot set for unnecessary set removal */
        bssp->_class = BioseqseqSet_class_nuc_prot;
        RenormalizeNucProtSets (sep, TRUE);
      }
    }
  }
  MemSet ((Pointer) objMgrFilter, FALSE, sizeof (objMgrFilter));
  objMgrFilter [OBJ_SEQFEAT] = TRUE;
  GatherObjectsInEntity (entityID, 0, NULL, MarkMovedGeneGbquals, (Pointer) &hasMarkedGenes, objMgrFilter);
  BasicSeqEntryCleanup (sep);
  SeqEntryExplore (sep, NULL, CleanupGenbankCallback);
  ConvertFullLenSourceFeatToDesc (sep);
  ConvertFullLenPubFeatToDesc (sep);
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);
  /* reindex, since PseudoGeneOverlap gets best overlapping gene */
  SeqMgrIndexFeatures (entityID, NULL);
  EntryChangeImpFeat(sep);     /* change any CDS ImpFeat to real CdRegion */
  /* MoveRnaGBQualProductToName (sep); */ /* move rna gbqual product to rna-ref.ext.name */
  /* MoveProtGBQualProductToName (sep); */ /* move prot gbqual product to prot-ref.name */
  /* MoveCdsGBQualProductToName (sep); */ /* move cds gbqual product to prot-ref.name */
  /* MoveFeatGBQualsToFields (sep); */ /* move feature partial, exception to fields */
  /* ExtendGeneFeatIfOnMRNA (0, sep); */ /* gene on mRNA is full length */

  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);
  VisitBioseqsInSep (sep, NULL, ExtendSingleGeneOnMRNA);

  RemoveBioSourceOnPopSet (sep, NULL);
  RemoveMolInfoOnPopSet (sep, NULL);
  /*
  SeqEntryExplore (sep, NULL, DeleteMultipleTitles);
  */
  SeqEntryExplore (sep, NULL, RemoveMultipleTitles);
  SeqEntryExplore (sep, NULL, MergeMultipleDates);
  VisitPubdescsInSep (sep, (Pointer) &muid, LookForUniqueMuid);
  VisitDescriptorsInSep (sep, (Pointer) &muid, RemoveBadPubDescr);
  VisitFeaturesInSep (sep, (Pointer) &muid, RemoveBadPubFeat);
  SeqEntryExplore (sep, (Pointer) &muid, MergePubFigInChain);
  VisitFeaturesInSep (sep, NULL, CorrectSfpExceptText);
  if (! isEmblOrDdbj) {
    SeqEntryExplore (sep, NULL, MergeEquivCitSubs);
  }
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
  EntryMergeDupBioSources (sep); /* do before and after SE2A3 */
  LoopSeqEntryToAsn3 (sep, TRUE, FALSE, taxfun, taxmerge, gpipeMode, isEmblOrDdbj);
  /* EntryStripSerialNumber(sep); */ /* strip citation serial numbers */
  MovePopPhyMutPubs (sep);
  EntryChangeGBSource (sep);   /* at least remove redundant information in GBBlocks */
  EntryCheckGBBlock (sep);
  SeqEntryMoveDbxrefs (sep); /* db_xref gbqual to sfp->dbxref */
  EntryMergeDupBioSources (sep);
  SeqEntryExplore (sep, NULL, GetRidOfEmptyFeatsDescCallback);
  CleanUpPseudoProductsEx (entityID, sep, doPseudo);
  RenormalizeNucProtSets (sep, TRUE);
  /*
  StripTitleFromProtsInNucProts (sep);
  */
  MarkBadProtTitlesInNucProts (sep);
  MoveFeatsFromPartsSet (sep);
  move_cds_ex (sep, doPseudo);
  SeqEntryExplore (sep, NULL, MolInfoUpdate);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
  /* do these again, since SE2A3 can create full length source feature */
  SeqEntryExplore (sep, NULL, CleanupGenbankCallback);
  ConvertFullLenSourceFeatToDesc (sep);
  ConvertFullLenPubFeatToDesc (sep);
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);
  /* VisitBioseqsInSep (sep, NULL, BarCodeTechToKeyword); */
  UpdateReplacedECNumbers (sep);
  if (GetAppProperty ("NcbiTbl2Asn") == NULL) {
    DeleteBadECNumbers (sep);
  }
  /* reindex, since CdEndCheck (from CdCheck) gets best overlapping gene */
  SeqMgrIndexFeatures (entityID, NULL);
  biop = NULL;
  if (VisitBioSourcesInSep (sep, (Pointer) &biop, FindSingleBioSource) == 1) {
    code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
    gncp = GeneticCodeFind (code, NULL);
    if (gncp == NULL) {
      gncp = GeneticCodeFind (1, NULL);
    }
    if (gncp != NULL) {
      codes = NULL;
      for (vnp = (ValNodePtr) gncp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == 3) {
          codes = (CharPtr) vnp->data.ptrvalue;
        }
      }
      if (codes != NULL) {
        VisitFeaturesInSep (sep, (Pointer) codes, CorrectTrnaCodons);
      }
    }
  }
  VisitBioSourcesInSep (sep, NULL, CleanupOldName);
  VisitBioSourcesInSep (sep, NULL, CleanupOrgModNote);
  CdCheck (sep, NULL);
  /* do again to catch occasional duplicate pub on some set components */
  SeqEntryPubsAsn4 (sep, isEmblOrDdbj);
  SeqMgrIndexFeatures (entityID, NULL);
  if (hasMarkedGenes) {
    MemSet ((Pointer) objMgrFilter, FALSE, sizeof (objMgrFilter));
    objMgrFilter [OBJ_SEQFEAT] = TRUE;
    GatherObjectsInEntity (entityID, 0, NULL, DeleteBadMarkedGeneXrefs, NULL, objMgrFilter);
  }
  if (! gpipeMode) {
    if (! isEmblOrDdbj) {   /* for now leave gene xrefs on EMBL and DDBJ */
      VisitBioseqsInSep (sep, (Pointer) &isEmblOrDdbj, RemoveUnneededGeneXrefs);
    }
  }
  ResynchCodingRegionPartials (sep);
  ResynchMessengerRNAPartials (sep);
  ResynchProteinPartials (sep);
  InstantiateProteinTitles (entityID, NULL);
  SeqMgrClearFeatureIndexes (entityID, NULL);
  MakeNcbiCleanupObject (sep, gpipeMode);
  BasicSeqEntryCleanup (sep);
  if (! gpipeMode) {
    NormalizeDescriptorOrder (sep);
  }
  TransTableFreeAll ();
  ErrSetMessageLevel (msev);
  ErrSetLogLevel (lsev);
  SeqEntrySetScope (oldscope);
}

extern void SeriousSeqEntryCleanup (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  SeriousSeqEntryCleanupEx (sep, taxfun, taxmerge, TRUE, FALSE);
}

extern void SeriousSeqEntryCleanupBulk (SeqEntryPtr sep)

{
  SeriousSeqEntryCleanupEx (sep, NULL, NULL, FALSE, FALSE);
}

extern void GpipeSeqEntryCleanup (SeqEntryPtr sep)

{
  SeriousSeqEntryCleanupEx (sep, NULL, NULL, TRUE, TRUE);
}

extern void SeriousSeqAnnotCleanup (SeqAnnotPtr sap)

{
  AnnotDescrPtr  adp;
  DatePtr        dp;
  AnnotDescrPtr  last;
  UserObjectPtr  uop;

  if (sap == NULL) return;

  RemoveAllSeqAnnotCleanupUserObjs (sap);

  BasicSeqAnnotCleanup (sap);

  dp = DateCurr ();
  if (dp == NULL) return;

  uop = CreateNcbiCleanupUserObject ();
  if (uop == NULL) return;

  AddStringToNcbiCleanupUserObject (uop, "method", "SeriousSeqAnnotCleanup");
  AddIntegerToNcbiCleanupUserObject (uop, "version", NCBI_CLEANUP_VERSION);

  AddIntegerToNcbiCleanupUserObject (uop, "month", dp->data [2]);
  AddIntegerToNcbiCleanupUserObject (uop, "day", dp->data [3]);
  AddIntegerToNcbiCleanupUserObject (uop, "year", dp->data [1] + 1900);

  adp = AnnotDescrNew (NULL);
  if (adp == NULL) return;

  adp->choice = Annot_descr_user;
  adp->data.ptrvalue = uop;

  if (sap->desc == NULL) {
    sap->desc = adp;
  } else {
    last = sap->desc;
    while (last->next != NULL) {
      last = last->next;
    }
    last->next = adp;
  }
}

static void CorrectGenCodeIndexedCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  Int2Ptr         pGenCode;
  ValNodePtr      vnp;
  Boolean         need_replacement = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION 
      || sfp->data.value.ptrvalue == NULL
      || userdata == NULL) return;
 
  pGenCode = (Int2Ptr) userdata;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp->genetic_code != NULL
      && crp->genetic_code->choice == 254) {
    if (crp->genetic_code->data.ptrvalue == NULL) {
      vnp = ValNodeNew (NULL);
      vnp->choice = 2;
      vnp->data.intvalue = (Int4) *pGenCode;
    } else {
      vnp = crp->genetic_code->data.ptrvalue;
      if (vnp->next == NULL && vnp->choice == 2) {
        vnp->data.intvalue = (Int4) *pGenCode;
      } else {
        need_replacement = TRUE;
      }
    }
  } else {
    need_replacement = TRUE;
  }
  if (need_replacement) {
    gc = GeneticCodeNew ();
    if (gc == NULL) return;
    crp->genetic_code = GeneticCodeFree (crp->genetic_code);
    vnp = ValNodeNew (NULL);
    gc->data.ptrvalue = vnp;
    if (vnp != NULL) {
      vnp->choice = 2;
      vnp->data.intvalue = (Int4) *pGenCode;
    }
    crp->genetic_code = gc;
  }
}

static void CorrectGenCodesBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;

  if (bsp == NULL || userdata == NULL) return;
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext)) {
    CorrectGenCodeIndexedCallback (sfp, userdata);
  }

}

typedef struct gencodescan {
  Boolean mito;
  Boolean plastid;
  Int2    nuclCode;
  Int2    mitoCode;
  Int2    pstdCode;
  Boolean already_found;
} GenCodeScanData, PNTR GenCodeScanPtr;

static void JustGetGenCodeFromOrgRef (OrgRefPtr orp, GenCodeScanPtr gp)
{
  OrgNamePtr onp;

  if (orp == NULL || orp->orgname == NULL || gp == NULL || gp->already_found) return;
  onp = orp->orgname;

  gp->nuclCode = onp->gcode;
  gp->mitoCode = onp->mgcode;
  gp->pstdCode = onp->pgcode;
}

static void JustGetGenCodeFromBiop (BioSourcePtr biop, GenCodeScanPtr gp)
{
  if (biop == NULL || gp == NULL) return;
  if (gp->already_found && !biop->is_focus) return;

  gp->mito = (Boolean) (biop->genome == GENOME_kinetoplast ||
                        biop->genome == GENOME_mitochondrion ||
                        biop->genome == GENOME_hydrogenosome);

  gp->plastid = (Boolean) (biop->genome == GENOME_chloroplast ||
                                biop->genome == GENOME_chromoplast ||
                                biop->genome == GENOME_plastid ||
                                biop->genome == GENOME_cyanelle ||
                                biop->genome == GENOME_apicoplast ||
                                biop->genome == GENOME_leucoplast ||
                                biop->genome == GENOME_proplastid ||
                                biop->genome == GENOME_chromatophore);

  JustGetGenCodeFromOrgRef (biop->org, gp);
  gp->already_found = TRUE;
}


static void JustGetGenCodeFromFeat (SeqFeatPtr sfp, Pointer userdata) 
{
  GenCodeScanPtr gp;

  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;

  gp = (GenCodeScanPtr) userdata;

  JustGetGenCodeFromBiop (sfp->data.value.ptrvalue, gp);
}

static void JustGetGenCodeFromDesc (SeqDescrPtr sdp, Pointer userdata)
{
  GenCodeScanPtr gp;

  if (sdp == NULL || userdata == NULL || sdp->choice != Seq_descr_source) return;

  gp = (GenCodeScanPtr) userdata;

  JustGetGenCodeFromBiop (sdp->data.ptrvalue, gp);
}

static Int2 JustGetGenCodeForSeqEntry (SeqEntryPtr sep) 
{
  GenCodeScanData gd;

  gd.already_found = FALSE;
  gd.mito = FALSE;
  gd.mitoCode = 0;
  gd.nuclCode = 0;
  gd.pstdCode = 0;
  gd.plastid = FALSE;

  VisitDescriptorsInSep (sep, &gd, JustGetGenCodeFromDesc);
  VisitFeaturesInSep (sep, &gd, JustGetGenCodeFromFeat);

  if (gd.plastid) {
    if (gd.pstdCode > 0) {
      return gd.pstdCode;
    } else {
      return 11;
    }
  } else if (gd.mito) {
    return gd.mitoCode;
  } else {
    return gd.nuclCode;
  }
}

extern void CorrectGenCodes (SeqEntryPtr sep, Uint2 entityID)

{
  BioseqSetPtr  bssp;
  Int2          genCode;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        CorrectGenCodes (sep, entityID);
      }
      return;
    }
  }

  genCode = JustGetGenCodeForSeqEntry(sep);
  VisitFeaturesInSep (sep, &genCode, CorrectGenCodeIndexedCallback);
  VisitBioseqsInSep (sep, &genCode, CorrectGenCodesBioseqCallback);
}



/* ConvertSegSetToDeltaSeq section */

static void FindSegSet (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqSetPtr      bssp;
  SeqEntryPtr PNTR  sepp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_segset) return;
  sepp = (SeqEntryPtr PNTR) mydata;
  if (sepp == NULL) return;
  *sepp = sep;
}

static void LookForMixedMols (BioseqPtr bsp, Pointer userdata)

{
  Uint1Ptr  molp;

  if (bsp == NULL) return;
  molp = (Uint1Ptr) userdata;
  if (molp == NULL) return;

  /* Boolean OR bsp->mols within segset to look for mixtures */
  *molp |= bsp->mol;
}

static void CheckForMissingMolInfo (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr            bp;
  SeqMgrDescContext  context;
  SeqDescrPtr        sdp;

  if (bsp == NULL) return;
  bp = (BoolPtr) userdata;
  if (bp == NULL) return;

  if (bsp->repr == Seq_repr_virtual) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL) {
    *bp = TRUE;
  }
}

typedef struct gbdata {
  GBBlockPtr  gbp;
  CharPtr     source;
  CharPtr     origin;
  Boolean     mixedsources;
  Boolean     mixedorigins;
} GBData, PNTR GBDataPtr;

static void FindGenBankDiffs (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  GBDataPtr   gbdp;
  GBBlockPtr  gbp;

  if (sdp == NULL || sdp->choice != Seq_descr_genbank) return;
  gbdp = (GBDataPtr) userdata;
  if (gbdp == NULL) return;

  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  if (gbp == NULL) return;

  if (gbdp->gbp == NULL) {
    gbdp->gbp = gbp;
    gbdp->source = gbp->source;
    gbdp->origin = gbp->origin;
    return;
  }

  if (StringICmp (gbp->source, gbdp->source) != 0) {
    gbdp->mixedsources = TRUE;
  }
  if (StringICmp (gbp->origin, gbdp->origin) != 0) {
    gbdp->mixedorigins = TRUE;
  }
}

static void CopyFirstGBBlock (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  GBBlockPtr       gbp;
  GBBlockPtr PNTR  gbpp;


  if (sdp == NULL || sdp->choice != Seq_descr_genbank) return;
  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  if (gbp == NULL) return;

  gbpp = (GBBlockPtr PNTR) userdata;
  if (gbpp == NULL) return;

  if (*gbpp != NULL) return;
  *gbpp = (GBBlockPtr) AsnIoMemCopy (gbp, (AsnReadFunc) GBBlockAsnRead, (AsnWriteFunc) GBBlockAsnWrite);
}

static void CopyFirstMolInfo (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  MolInfoPtr       mip;
  MolInfoPtr PNTR  mipp;


  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;

  mipp = (MolInfoPtr PNTR) userdata;
  if (mipp == NULL) return;

  if (*mipp != NULL) return;
  *mipp = (MolInfoPtr) AsnIoMemCopy (mip, (AsnReadFunc) MolInfoAsnRead, (AsnWriteFunc) MolInfoAsnWrite);
}

static void CopyFirstTitle (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  CharPtr       title;
  CharPtr PNTR  titlep;


  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  title = (CharPtr) sdp->data.ptrvalue;
  if (title == NULL) return;

  titlep = (CharPtr PNTR) userdata;
  if (titlep == NULL) return;

  if (*titlep != NULL) return;
  *titlep = (CharPtr) StringSave (title);
}

static void AddPartAccns (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Char        buf [64];
  GBBlockPtr  gbp;
  SeqIdPtr    sip;
 
  if (bsp == NULL) return;
  gbp = (GBBlockPtr) userdata;
  if (gbp == NULL) return;

  if (bsp->repr == Seq_repr_virtual) return;

  sip = SeqIdFindBestAccession (bsp->id);
  if (sip == NULL) return;

  SeqIdWrite (sip, buf, PRINTID_TEXTID_ACCESSION, sizeof (buf));
  if (StringHasNoText (buf)) return;

  ValNodeCopyStr (&(gbp->extra_accessions), 0, buf);
}

static void AddPartHist (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Char        buf [64];
  BioseqPtr   deltabsp;
  SeqHistPtr  shp;
  SeqIdPtr    sip;

  if (bsp == NULL) return;
  deltabsp = (BioseqPtr) userdata;
  if (deltabsp == NULL) return;

  if (bsp->repr == Seq_repr_virtual) return;

  sip = SeqIdFindBestAccession (bsp->id);
  if (sip == NULL) return;

  SeqIdWrite (sip, buf, PRINTID_TEXTID_ACCESSION, sizeof (buf));
  if (StringHasNoText (buf)) return;

  shp = ParseStringIntoSeqHist (deltabsp->hist, buf);
  if (deltabsp->hist == NULL) {
    deltabsp->hist = shp;
  }
}

static Boolean LIBCALLBACK AddSegToDeltaSeq (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  CharPtr     bases;
  BioseqPtr   bsp;
  BioseqPtr   deltabsp;
  IntFuzzPtr  ifp;
  SeqLocPtr   loc;
  SeqIdPtr    sip;
  SeqLitPtr   slitp;

  if (slp == NULL || context == NULL) return FALSE;
  deltabsp = (BioseqPtr) context->userdata;
  if (deltabsp == NULL) return FALSE;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc == NULL) return TRUE;
    sip = SeqLocId (loc);
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqFind (sip);
  if (bsp == NULL) return TRUE;

  if (bsp->repr == Seq_repr_virtual) {
    if (deltabsp->seq_ext != NULL) {
      /* insert gap of unknown length (by convention, 100 bases) between the previous segment and this one. */
      slitp = (SeqLitPtr) MemNew (sizeof (SeqLit));
      if (slitp != NULL) {
        slitp->length = bsp->length;
        if (slitp->length == 100) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            slitp->fuzz = ifp;
          }
        }
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slitp);
        deltabsp->length += slitp->length;
      }
    }
    return TRUE;
  }

  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) return TRUE;

  slitp = (SeqLitPtr) MemNew (sizeof (SeqLit));
  if (slitp != NULL) {
    slitp->length = StringLen (bases);
    ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slitp);
    slitp->seq_data = (SeqDataPtr) BSNew (slitp->length);
    slitp->seq_data_type = Seq_code_iupacna;
    AddBasesToByteStore ((ByteStorePtr) slitp->seq_data, bases);
    deltabsp->length += slitp->length;
  }

  return TRUE;
}

static void MoveAnnotsToDeltaSeq (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  BioseqPtr     deltabsp;
  SeqAnnotPtr   nextsap;
  SeqFeatPtr    nextsfp;
  Pointer PNTR  prevsap;
  Pointer PNTR  prevsfp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL) return;
  deltabsp = (BioseqPtr) mydata;
  if (deltabsp == NULL) return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;

  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        *(prevsfp) = sfp->next;
        sfp->next = NULL;
        AddFeatToBioseq (sfp, deltabsp);
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static Boolean RptUnitInBaseRange (CharPtr str, Int4Ptr fromP, Int4Ptr toP)

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
  if (fromP != NULL) {
    *fromP = val - 1;
  }
  ptr += 2;
  if (StringHasNoText (ptr)) return FALSE;
  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  if (toP != NULL) {
    *toP = val - 1;
  }
  return TRUE;
}

static BioseqPtr FindFirstLocalBioseq (SeqLocPtr loc)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;
  SeqLocPtr  slp = NULL;

  if (loc == NULL) return NULL;

  while ((slp = SeqLocFindNext (loc, slp)) != NULL) {
    sip = SeqLocId (slp);
    if (sip != NULL) {
      bsp = BioseqFindCore (sip);
      if (bsp != NULL) return bsp;
    }
  }

  return NULL;
}

static void MapSegFeatToMaster (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr     bsp, ptbsp;
  Char          buf [64];
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Int4          from, to;
  GBQualPtr     gbq;
  Boolean       hasNulls;
  Int4          lim;
  Boolean       noLeft;
  Boolean       noRight;
  ValNodePtr    partiallist = NULL, emptypartials = NULL;
  RnaRefPtr     rrp;
  SeqInt        sint;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLocPtr     slp = NULL;
  tRNAPtr       trna;
  ValNode       vn;

  if (sfp == NULL || sfp->location == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol)) return;
  if (bsp->repr != Seq_repr_seg) return;

  partiallist = GetSeqLocPartialSet (sfp->location);
  CheckSeqLocForPartialEx (sfp->location, &noLeft, &noRight, &lim);
  hasNulls = LocationHasNullsBetween (sfp->location);

  if (sfp->data.choice == SEQFEAT_GENE) {
    slp = SeqLocMergeExEx (bsp, sfp->location, NULL, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE);
    hasNulls = FALSE;
    sfp->partial = FALSE;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    slp = SeqLocMergeExEx (bsp, sfp->location, NULL, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE);
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    slp = SeqLocMergeExEx (bsp, sfp->location, NULL, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE);
  } else {
    slp = SeqLocMergeExEx (bsp, sfp->location, NULL, FALSE, TRUE, FALSE, hasNulls, FALSE, FALSE);
  }
  if (slp == NULL) {
    ValNodeFree (partiallist);
    return;
  }

  ptbsp = FindFirstLocalBioseq (sfp->location);

  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  emptypartials = GetSeqLocPartialSet (sfp->location);
  FreeAllFuzz (sfp->location);
  SetSeqLocPartialEx (sfp->location, noLeft, noRight, lim);
  if (ValNodeLen (partiallist) == ValNodeLen (emptypartials)) {
    SetSeqLocPartialSet (sfp->location, partiallist);
  }
  ValNodeFree (partiallist);
  ValNodeFree (emptypartials);
  sfp->partial = (sfp->partial || noLeft || noRight || hasNulls);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL && crp->code_break != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          slp = SeqLocMerge (bsp, cbp->loc, NULL, FALSE, TRUE, FALSE);
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
          slp = SeqLocMerge (bsp, trna->anticodon, NULL, FALSE, TRUE, FALSE);
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

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "rpt_unit_range") != 0) continue;
    if (! RptUnitInBaseRange (gbq->val, &from, &to)) continue;
     if (ptbsp == NULL || ptbsp->repr != Seq_repr_raw) continue;
    sip = SeqIdFindBest (ptbsp->id, 0);
    if (sip == NULL) continue;
    MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
    MemSet ((Pointer) &vn, 0, sizeof (ValNode));
    if (from < to) {
      sint.from = from;
      sint.to = to;
      sint.strand = Seq_strand_plus;
    } else {
      sint.to = from;
      sint.from = to;
      sint.strand = Seq_strand_minus;
    }
    sint.id = sip;
    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
    slp = SeqLocMerge (bsp, &vn, NULL, FALSE, TRUE, FALSE);
    if (slp != NULL) {
      if (slp->choice == SEQLOC_INT) {
        sintp = (SeqIntPtr) slp->data.ptrvalue;
        if (sintp != NULL) {
          buf [0] = '\0';
          if (sintp->strand == Seq_strand_minus) {
            sprintf (buf, "%ld..%ld", (long) sintp->to + 1, (long) sintp->from + 1);
          } else {
            sprintf (buf, "%ld..%ld", (long) sintp->from + 1, (long) sintp->to + 1);
          }
          if (StringDoesHaveText (buf)) {
            gbq->val = MemFree (gbq->val);
            gbq->val = StringSave (buf);
          }
        }
      }
      SeqLocFree (slp);
    }
  }
}

static void PartDescToFeatProc (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  PubdescPtr    pdp;
  SeqFeatPtr    sfp;

  if (sdp == NULL) return;
  bsp = (BioseqPtr) userdata;
  if (bsp == NULL) return;

  if (sdp->choice == Seq_descr_pub) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp == NULL) return;
    sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_PUB, NULL);
    if (sfp != NULL) {
      sfp->data.value.ptrvalue = AsnIoMemCopy (pdp, (AsnReadFunc) PubdescAsnRead, (AsnWriteFunc) PubdescAsnWrite);
    }
  } else if (sdp->choice == Seq_descr_source) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop == NULL) return;
    sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_BIOSRC, NULL);
    if (sfp != NULL) {
      sfp->data.value.ptrvalue = AsnIoMemCopy (biop, (AsnReadFunc) BioSourceAsnRead, (AsnWriteFunc) BioSourceAsnWrite);
    }
  }
}

static void MovePartDescToFeat (
  BioseqPtr bsp,
  Pointer userdata
)

{
  if (bsp == NULL) return;
  if (bsp->repr == Seq_repr_virtual) return;
  VisitDescriptorsOnBsp (bsp, (Pointer) bsp, PartDescToFeatProc);
}

static void PartCitSubDatesCompProc (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  Int2          compare;
  CitSubPtr     csp;
  DatePtr PNTR  dpp;
  PubdescPtr    pdp;
  ValNodePtr    vnp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub) return;
  dpp = (DatePtr PNTR) userdata;
  if (dpp == NULL) return;

  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Sub) continue;
    csp = (CitSubPtr) vnp->data.ptrvalue;
    if (csp == NULL || csp->date == NULL) continue;
    if (*dpp == NULL) {
      *dpp = csp->date;
    } else {
      compare = DateMatch (*dpp, csp->date, TRUE);
      if (compare == 1) {
        *dpp = csp->date;
      }
    }
  }
}

static void CompPartCitSubDates (
  BioseqPtr bsp,
  Pointer userdata
)

{
  if (bsp == NULL) return;
  VisitDescriptorsOnBsp (bsp, userdata, PartCitSubDatesCompProc);
}

static void PartCitSubDatesSyncProc (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  CitSubPtr   csp;
  DatePtr     dp;
  PubdescPtr  pdp;
  ValNodePtr  vnp;

  if (sdp == NULL || sdp->choice != Seq_descr_pub) return;
  dp = (DatePtr) userdata;
  if (dp == NULL) return;

  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Sub) continue;
    csp = (CitSubPtr) vnp->data.ptrvalue;
    if (csp == NULL || csp->date == NULL) continue;
    if (DateMatch (csp->date, dp, TRUE) != 0) {
      csp->date = DateFree (csp->date);
      csp->date = DateDup (dp);
    }
  }
}

static void SyncPartCitSubDates (
  BioseqPtr bsp,
  Pointer userdata
)

{
  VisitDescriptorsOnBsp (bsp, userdata, PartCitSubDatesSyncProc);
}

extern void SegSeqNullToVirtual (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  Uint2         entityID;
  GBData        gbd;
  Boolean       hasnulls = FALSE;
  Boolean       hasvirt = FALSE;
  Boolean       missingMolInfo = FALSE;
  Uint1         mol = 0;
  SeqEntryPtr   partssep = NULL;
  BioseqSetPtr  partsset = NULL;
  BioseqPtr     segbsp = NULL;
  SeqEntryPtr   segsep = NULL;
  SeqEntryPtr   segseq = NULL;
  BioseqSetPtr  segset = NULL;
  SeqLocPtr     slp;
  BioseqPtr     virtbsp;
  ValNode       vn;
  SeqLocPtr     vslp;

  if (sep == NULL) return;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;

  /* skip pop/phy/mut/eco/wgs sets for now */
  if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) return;
  if (bssp->_class == BioseqseqSet_class_wgs_set) return;
  if (bssp->_class == BioseqseqSet_class_small_genome_set) return;

  /* find SeqEntryPtr parent of single seg set */
  SeqEntryExplore (sep, (Pointer) &segsep, FindSegSet);
  if (segsep == NULL) return;

  /* do not handle just segset without nucprot set wrapper for now */
  if (segsep == sep) return;

  /* skip the few cases of mixed molecule types */
  VisitBioseqsInSep (segsep, (Pointer) &mol, LookForMixedMols);
  if (mol != Seq_mol_dna && mol != Seq_mol_rna && mol != Seq_mol_aa && mol != Seq_mol_na && mol != Seq_mol_other) return;

  /* skip mixed gbblock source or origin fields for now */
  MemSet ((Pointer) &gbd, 0, sizeof (gbd));
  VisitDescriptorsInSep (sep, (Pointer) &gbd, FindGenBankDiffs);
  if (gbd.mixedsources || gbd.mixedorigins) return;

  /* avoid copying protein molinfo if nucleotide molinfo is missing */
  VisitBioseqsInSep (sep, (Pointer) &missingMolInfo, CheckForMissingMolInfo);
  if (missingMolInfo) return;

  if (! IS_Bioseq_set (segsep)) return;
  segset = (BioseqSetPtr) segsep->data.ptrvalue;
  if (segset == NULL) return;
  segseq = segset->seq_set;
  if (segseq == NULL) return;
  if (! IS_Bioseq (segseq)) return;
  segbsp = (BioseqPtr) segseq->data.ptrvalue;
  if (segbsp == NULL) return;
  if (segbsp->repr != Seq_repr_seg) return;
  partssep = segseq->next;
  if (partssep == NULL) return;
  if (! IS_Bioseq_set (partssep)) return;
  partsset = (BioseqSetPtr) partssep->data.ptrvalue;
  if (partsset == NULL) return;

  entityID = segbsp->idx.entityID;

  /* check to see if it doesn't need conversion */
  MemSet ((Pointer) &vn, 0, sizeof (vn));
  vn.choice = SEQLOC_MIX;
  vn.extended = 0;
  vn.data.ptrvalue = segbsp->seq_ext;
  vn.next = NULL;
  slp = SeqLocFindNext (&vn, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_NULL) {
      hasnulls = TRUE;
    }
    slp = SeqLocFindNext (&vn, slp);
  }

  /* if no nulls or virtuals, add a null and then fill in between all components */
  if (! hasnulls) {
    if (! hasvirt) {
      slp = SeqLocFindNext (&vn, NULL);
      if (slp != NULL) {
        vslp = ValNodeNew (NULL);
        if (vslp != NULL) {
          vslp->choice = SEQLOC_NULL;
          vslp->next = slp->next;
          slp->next = vslp;
          NormalizeNullsBetween (&vn);
          hasnulls = TRUE;
        }
      }
    }
  }

  if (! hasnulls) return;

  /* virtual part of conventional 100 base gap length */
  virtbsp = BioseqNew ();
  if (virtbsp == NULL) return;
  virtbsp->id = MakeUniqueSeqID ("virtual_");
  virtbsp->repr = Seq_repr_virtual;
  virtbsp->mol = segbsp->mol;
  virtbsp->seq_data_type = 0;
  virtbsp->seq_ext_type = 0;
  virtbsp->seq_data = NULL;
  virtbsp->seq_ext = NULL;
  virtbsp->length = 100;

  /* put virtual segments between real segments */
  MemSet ((Pointer) &vn, 0, sizeof (vn));
  vn.choice = SEQLOC_MIX;
  vn.extended = 0;
  vn.data.ptrvalue = segbsp->seq_ext;
  vn.next = NULL;
  slp = SeqLocFindNext (&vn, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_NULL) {
      slp->choice = SEQLOC_WHOLE;
      slp->data.ptrvalue = SeqIdDup (virtbsp->id);
      segbsp->length += virtbsp->length;
    }
    slp = SeqLocFindNext (&vn, slp);
  }

  /* package virtual bioseq in parts */
  ValNodeAddPointer (&(partsset->seq_set), 1, (Pointer) virtbsp);

  /* reindex for new segmented bioseq length */
  SeqMgrIndexFeatures (entityID, NULL);
}

static void ForceSegSeqNullToVirtual (SeqEntryPtr sep)

{
  Uint2         entityID;
  Boolean       hasnulls = FALSE;
  Boolean       hasvirt = FALSE;
  SeqEntryPtr   partssep = NULL;
  BioseqSetPtr  partsset = NULL;
  BioseqPtr     segbsp = NULL;
  SeqEntryPtr   segsep = NULL;
  SeqEntryPtr   segseq = NULL;
  BioseqSetPtr  segset = NULL;
  SeqLocPtr     slp;
  BioseqPtr     vbsp;
  BioseqPtr     virtbsp;
  ValNode       vn;
  SeqLocPtr     vslp;

  if (sep == NULL) return;

  SeqEntryExplore (sep, (Pointer) &segsep, FindSegSet);
  if (segsep == NULL) return;

  if (! IS_Bioseq_set (segsep)) return;
  segset = (BioseqSetPtr) segsep->data.ptrvalue;
  if (segset == NULL) return;
  segseq = segset->seq_set;
  if (segseq == NULL) return;
  if (! IS_Bioseq (segseq)) return;
  segbsp = (BioseqPtr) segseq->data.ptrvalue;
  if (segbsp == NULL) return;
  if (segbsp->repr != Seq_repr_seg) return;
  partssep = segseq->next;
  if (partssep == NULL) return;
  if (! IS_Bioseq_set (partssep)) return;
  partsset = (BioseqSetPtr) partssep->data.ptrvalue;
  if (partsset == NULL) return;

  entityID = segbsp->idx.entityID;

  /* check to see if it doesn't need conversion */
  MemSet ((Pointer) &vn, 0, sizeof (vn));
  vn.choice = SEQLOC_MIX;
  vn.extended = 0;
  vn.data.ptrvalue = segbsp->seq_ext;
  vn.next = NULL;
  slp = SeqLocFindNext (&vn, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_NULL) {
      hasnulls = TRUE;
    } else {
      vbsp = BioseqFindFromSeqLoc (slp);
      if (vbsp != NULL && vbsp->repr == Seq_repr_virtual) {
        hasvirt = TRUE;
      }
    }
    slp = SeqLocFindNext (&vn, slp);
  }

  /* if no nulls or virtuals, add a null and then fill in between all components */
  if (! hasnulls) {
    if (! hasvirt) {
      slp = SeqLocFindNext (&vn, NULL);
      if (slp != NULL) {
        vslp = ValNodeNew (NULL);
        if (vslp != NULL) {
          vslp->choice = SEQLOC_NULL;
          vslp->next = slp->next;
          slp->next = vslp;
          NormalizeNullsBetween (&vn);
          hasnulls = TRUE;
        }
      }
    }
  }

  if (! hasnulls) return;

  /* virtual part of conventional 100 base gap length */
  virtbsp = BioseqNew ();
  if (virtbsp == NULL) return;
  virtbsp->id = MakeUniqueSeqID ("virtual_");
  virtbsp->repr = Seq_repr_virtual;
  virtbsp->mol = segbsp->mol;
  virtbsp->seq_data_type = 0;
  virtbsp->seq_ext_type = 0;
  virtbsp->seq_data = NULL;
  virtbsp->seq_ext = NULL;
  virtbsp->length = 100;

  /* put virtual segments between real segments */
  MemSet ((Pointer) &vn, 0, sizeof (vn));
  vn.choice = SEQLOC_MIX;
  vn.extended = 0;
  vn.data.ptrvalue = segbsp->seq_ext;
  vn.next = NULL;
  slp = SeqLocFindNext (&vn, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_NULL) {
      slp->choice = SEQLOC_WHOLE;
      slp->data.ptrvalue = SeqIdDup (virtbsp->id);
      segbsp->length += virtbsp->length;
    }
    slp = SeqLocFindNext (&vn, slp);
  }

  /* package virtual bioseq in parts */
  ValNodeAddPointer (&(partsset->seq_set), 1, (Pointer) virtbsp);

  /* reindex for new segmented bioseq length */
  SeqMgrIndexFeatures (entityID, NULL);
}

extern Boolean ConvertSegSetToDeltaSeq (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  GBData        gbd;
  DatePtr       cp = NULL;
  DatePtr       dp = NULL;
  BioseqPtr     deltabsp;
  Uint2         entityID;
  GBBlockPtr    gbp = NULL;
  MolInfoPtr    mip = NULL;
  Boolean       missingMolInfo = FALSE;
  Uint1         mol = 0;
  SeqEntryPtr   partssep = NULL;
  BioseqPtr     segbsp = NULL;
  SeqEntryPtr   segsep = NULL;
  SeqEntryPtr   segseq = NULL;
  BioseqSetPtr  segset = NULL;
  SeqIdPtr      tmpid;
  CharPtr       ttl = NULL;
  ValNodePtr    vnp;
  /*
  SeqLocPtr     slp;
  BioseqPtr     virtbsp;
  ValNode       vn;
  */

  if (sep == NULL) return FALSE;
  if (! IS_Bioseq_set (sep)) return FALSE;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return FALSE;

  /* skip pop/phy/mut/eco/wgs sets for now */
  if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) return FALSE;
  if (bssp->_class == BioseqseqSet_class_wgs_set) return FALSE;
  if (bssp->_class == BioseqseqSet_class_small_genome_set) return FALSE;

  /* find SeqEntryPtr parent of single seg set */
  SeqEntryExplore (sep, (Pointer) &segsep, FindSegSet);
  if (segsep == NULL) return FALSE;

  /* do not handle just segset without nucprot set wrapper for now */
  if (segsep == sep) return FALSE;

  /* skip the few cases of mixed molecule types */
  VisitBioseqsInSep (segsep, (Pointer) &mol, LookForMixedMols);
  if (mol != Seq_mol_dna && mol != Seq_mol_rna && mol != Seq_mol_aa && mol != Seq_mol_na && mol != Seq_mol_other) return FALSE;

  /* skip mixed gbblock source or origin fields for now */
  MemSet ((Pointer) &gbd, 0, sizeof (gbd));
  VisitDescriptorsInSep (sep, (Pointer) &gbd, FindGenBankDiffs);
  if (gbd.mixedsources || gbd.mixedorigins) return FALSE;

  /* avoid copying protein molinfo if nucleotide molinfo is missing */
  VisitBioseqsInSep (sep, (Pointer) &missingMolInfo, CheckForMissingMolInfo);
  if (missingMolInfo) return FALSE;

  /*
  ConvertSegSetsToDeltaSequences (sep);
  */

  if (! IS_Bioseq_set (segsep)) return FALSE;
  segset = (BioseqSetPtr) segsep->data.ptrvalue;
  if (segset == NULL) return FALSE;
  segseq = segset->seq_set;
  if (segseq == NULL) return FALSE;
  if (! IS_Bioseq (segseq)) return FALSE;
  segbsp = (BioseqPtr) segseq->data.ptrvalue;
  if (segbsp == NULL) return FALSE;
  if (segbsp->repr != Seq_repr_seg) return FALSE;
  partssep = segseq->next;
  if (partssep == NULL) return FALSE;

  /* synchronize dates of cit-sub descriptors on parts */
  VisitBioseqsInSep (partssep, (Pointer) &dp, CompPartCitSubDates);
  cp = DateDup (dp);
  VisitBioseqsInSep (partssep, (Pointer) cp, SyncPartCitSubDates);
  DateFree (cp);
  /* move pub descriptors up from parts before adding virtual sequence */
  SeqEntryPubsAsn4 (sep, FALSE);

  /* then convert any remaining pub and source descriptors on parts to features */
  VisitBioseqsInSep (partssep, NULL, MovePartDescToFeat);

  /* put virtual segments of conventional 100 bases between real segments */
  ForceSegSeqNullToVirtual (sep);

  /* map feature locations to segmented bioseq */
  VisitFeaturesInSep (sep, NULL, MapSegFeatToMaster);

  VisitDescriptorsInSep (segsep, (Pointer) &gbp, CopyFirstGBBlock);
  VisitDescriptorsInSep (segsep, (Pointer) &mip, CopyFirstMolInfo);
  VisitDescriptorsInSep (segsep, (Pointer) &ttl, CopyFirstTitle);

  if (gbp == NULL) {
    gbp = GBBlockNew ();
  }
  /* populate secondary accessions */
  VisitBioseqsInSep (partssep, (Pointer) gbp, AddPartAccns);

  /*
  deltabsp = GetDeltaSeqFromMasterSeg (segbsp);
  */

  deltabsp = BioseqNew ();
  if (deltabsp == NULL) return FALSE;
  deltabsp->id = MakeUniqueSeqID ("delta_");
  deltabsp->repr = Seq_repr_delta;
  deltabsp->mol = segbsp->mol;
  deltabsp->seq_data_type = 0;
  deltabsp->seq_ext_type = 4;
  deltabsp->seq_data = NULL;
  deltabsp->seq_ext = NULL;
  deltabsp->length = 0;

  /* populate Seq-hist.replaces */
  VisitBioseqsInSep (partssep, (Pointer) deltabsp, AddPartHist);

  /* construct delta seq from segmented parts */
  SeqMgrExploreSegments (segbsp, (Pointer) deltabsp, AddSegToDeltaSeq);

  /* move features in segset to delta seq */
  SeqEntryExplore (segsep, (Pointer) deltabsp, MoveAnnotsToDeltaSeq);

  /* insert delta sequence into chain */
  vnp = SeqDescrNew (NULL);
  if (vnp == NULL) return FALSE;
  vnp->choice = 1;
  vnp->data.ptrvalue = (Pointer) deltabsp;
  vnp->next = segsep->next;
  segsep->next = vnp;

  /* keep segmented bioseq IDs */
  tmpid = segbsp->id;
  segbsp->id = deltabsp->id;
  deltabsp->id = tmpid;

  SeqMgrDeleteFromBioseqIndex (segbsp);
  SeqMgrReplaceInBioseqIndex (deltabsp);

  /* remove old segset */
  segset->idx.deleteme = TRUE;
  DeleteMarkedObjects (0, OBJ_SEQENTRY, sep);

  /*
  BioseqFree (virtbsp);
  */

  /* add descriptors */
  if (ttl != NULL) {
    SeqDescrAddPointer (&(deltabsp->descr), Seq_descr_title, (Pointer) ttl);
  }
  if (mip != NULL) {
    SeqDescrAddPointer (&(deltabsp->descr), Seq_descr_molinfo, (Pointer) mip);
  }
  if (gbp != NULL) {
    SeqDescrAddPointer (&(deltabsp->descr), Seq_descr_genbank, (Pointer) gbp);
  }

  BioseqPack (deltabsp);

  /* indexes are out of date */
  SeqMgrClearFeatureIndexes (entityID, NULL);

  SeriousSeqEntryCleanup (sep, NULL, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);

  return TRUE;
}

static void ConvPartDescToFeatCallback (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_parts) return;

  VisitBioseqsInSet (bssp, NULL, MovePartDescToFeat);
}

extern void ConvertPartDescToFeat (SeqEntryPtr sep)

{
  VisitSetsInSep (sep, NULL, ConvPartDescToFeatCallback);
}

extern void SimpleAutoDef (SeqEntryPtr sep)

{
  BioseqPtr                     bsp = NULL;
  BioseqSetPtr                  bssp = NULL;
  ValNodePtr                    defline_clauses = NULL;
  Uint2                         entityID = 0;
  DeflineFeatureRequestList     feature_requests;
  Int4                          index;
  ValNodePtr                    modifier_indices = NULL;
  ModifierItemLocalPtr          modList;
  OrganismDescriptionModifiers  odmp;
  SeqEntryPtr                   oldscope;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    entityID = bsp->idx.entityID;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    entityID = bssp->idx.entityID;
  } else return;

  if (entityID < 1) return;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;

  InitFeatureRequests (&feature_requests);

  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  InitOrganismDescriptionModifiers (&odmp, sep);

  RemoveNucProtSetTitles (sep);
  oldscope = SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID, &feature_requests,
                              DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                              &defline_clauses);
  if (AreFeatureClausesUnique (defline_clauses)) {
    modifier_indices = GetModifierIndicesFromModList (modList);
  } else {
    modifier_indices = FindBestModifiers (sep, modList);
  }

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, &odmp);
  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL) {
    for (index = 0; index < NumDefLineModifiers (); index++) {
      ValNodeFree (modList [index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
  /*
  RemovePopsetTitles (sep);
  */
  AddPopsetTitles (sep, &feature_requests, DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);

  SeqEntrySetScope (oldscope);
}

