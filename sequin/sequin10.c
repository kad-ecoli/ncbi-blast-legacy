/*   sequin10.c
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
* File Name:  sequin10.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   9/3/2003
*
* $Revision: 1.562 $
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
#include <sqnutils.h>
#include <subutil.h>
#include <explore.h>
#include <edutil.h>
#include <tofasta.h>
#include <gbftdef.h>
#include <gbfeat.h>
#include <biosrc.h>
#include <findrepl.h>
#include <asnenbin.h>
#include <cdrgn.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macrodlg.h>
#include <macroapi.h>
#include <seqpanel.h>
#include <salpanel.h>
#include <alignmgr2.h>
#include <salpedit.h>
#include <objvalid.h>
#include <valapi.h>
#include <seqsub.h>
#include <tax3api.h>
#include <valid.h>
#include <vsmutil.h>

typedef struct deflineformdata {
  FEATURE_FORM_BLOCK

  ModifierItemLocalPtr modList;
  ButtoN PNTR modifier_btns;
  DeflineFeatureRequestList feature_requests;
  ButtoN    feature_btns[NumRemovableItems];
  GrouP     sourceListGrp;
  PopuP     customGrp;
  ButtoN    use_labels;
  ButtoN    keep_paren;
  ButtoN    exclude_sp;
  ButtoN    exclude_cf;
  ButtoN    exclude_aff;
  ButtoN    exclude_nr;
  ButtoN    allow_mod_at_end_of_taxname;
  ButtoN    include_country_extra;
  ButtoN    allow_semicolon_in_modifier;
  GrouP     clone_isolate_HIV_rule_num;
  PopuP     modLimit;
  PopuP     organelle_popup;
  BioseqPtr target_bsp;
  ButtoN    modify_only_target;
  ButtoN    suppress_alt_splice_phrase;
  ButtoN    remove_subfeatures;
  PopuP     featurePopup;
  GrouP     featureOptsGrp;
  GrouP     optional_features_grp;
  PopuP     misc_feat_parse_rule;
  GrouP     promoter_type;
  ButtoN    alternate_splice_flag;
  ButtoN    use_ncrna_note;
  ButtoN    suppress_locus_tags;
  ButtoN    gene_cluster_opp_strand;
  DialoG    suppressed_feature_list;
  GrouP     suppressed_feature_grp;
  Boolean   popset_only;
} DefLineFormData, PNTR DefLineFormPtr;

static void DefLineFormMessageProc (ForM f, Int2 mssg)

{
  DefLineFormPtr  dlfp;

  dlfp = (DefLineFormPtr) GetObjectExtra (f);
  if (dlfp != NULL) {
    if (dlfp->appmessage != NULL) {
      dlfp->appmessage (f, mssg);
    }
  }
}

static void CleanupDefLineForm (
  GraphiC g,
  VoidPtr data
)

{
  DefLineFormPtr  dlfp;
  Int4            i;

  dlfp = (DefLineFormPtr) data;
  if (dlfp != NULL) {
    if (dlfp->modList != NULL)
    {
      for (i=0; i < NumDefLineModifiers (); i++)
      {
        ValNodeFree (dlfp->modList[i].values_seen);
      }
      MemFree (dlfp->modList);
    }
    dlfp->modifier_btns = MemFree (dlfp->modifier_btns);
  }
  StdCleanupFormProc (g, data);
}

static void ChangeCustomPopup (PopuP p)

{
  DefLineFormPtr  dlfp;

  dlfp = (DefLineFormPtr) GetObjectExtra (p);
  if (dlfp == NULL) return;
  if (GetValue (p) == 1) {
    SafeDisable (dlfp->sourceListGrp);
  } else {
    SafeEnable (dlfp->sourceListGrp);
  }
}

static void ChangeFeaturePopup (PopuP p)
{
  DefLineFormPtr  dlfp;

  dlfp = (DefLineFormPtr) GetObjectExtra (p);
  if (dlfp == NULL) return;
  if (GetValue (p) == 2 || GetValue (p) == 3)
  {
    SafeDisable (dlfp->featureOptsGrp);
    SafeDisable (dlfp->suppressed_feature_grp);
    SafeDisable (dlfp->optional_features_grp);
    SafeDisable (dlfp->organelle_popup);
    SafeDisable (dlfp->alternate_splice_flag);
  }
  else
  {
    SafeEnable (dlfp->featureOptsGrp);
    SafeEnable (dlfp->suppressed_feature_grp);
    SafeEnable (dlfp->optional_features_grp);
    SafeEnable (dlfp->organelle_popup);
    SafeEnable (dlfp->alternate_splice_flag);
  }
}


static void DoAutoDefLine (ButtoN b)
{
  DefLineFormPtr dlfp;
  SeqEntryPtr sep;
  ValNodePtr modifier_indices = NULL;
  Int2 feature_index;
  OrganismDescriptionModifiers odmp;
  Int2 product_flag, feature_list_type;
  Int4 i;
  Boolean alternate_splice_flag;
  Boolean gene_cluster_opp_strand;

  dlfp = GetObjectExtra (b);
  if (b == NULL) return;
  Hide (dlfp->form);
  WatchCursor ();
  Update ();

  if (!dlfp->popset_only) {
    InitOrganismDescriptionModifiers (&odmp, NULL);
    odmp.use_labels = GetStatus (dlfp->use_labels);
    odmp.keep_paren = GetStatus (dlfp->keep_paren);
    odmp.exclude_sp = GetStatus (dlfp->exclude_sp);
    odmp.exclude_cf = GetStatus (dlfp->exclude_cf);
    odmp.exclude_aff = GetStatus (dlfp->exclude_aff);
    odmp.exclude_nr = GetStatus (dlfp->exclude_nr);
    odmp.include_country_extra = GetStatus (dlfp->include_country_extra);
    odmp.use_modifiers = TRUE;
    odmp.allow_semicolon_in_modifier = GetStatus (dlfp->allow_semicolon_in_modifier);
    odmp.allow_mod_at_end_of_taxname = !GetStatus (dlfp->allow_mod_at_end_of_taxname);

    if (dlfp->clone_isolate_HIV_rule_num == NULL)
    {
      odmp.clone_isolate_HIV_rule_num = clone_isolate_HIV_rule_want_both;
    }
    else
    {
      odmp.clone_isolate_HIV_rule_num = GetValue (dlfp->clone_isolate_HIV_rule_num);
    }

    odmp.max_mods = GetValue (dlfp->modLimit);
    if (odmp.max_mods > 1)
    {
      odmp.max_mods = odmp.max_mods - 1;
    }
    else
    {
      odmp.max_mods = -99;
    }

  }

  if (GetValue (dlfp->customGrp) == 1)
  {
    /* take all features */
    for (feature_index = 0; feature_index < NumDefLineModifiers (); feature_index++)
    {
      if (dlfp->modList[feature_index].any_present)
      {
        ValNodeAddInt (&modifier_indices, 0, feature_index);
      }
    }
  }
  else
  {
    /* take selected features */
    for (feature_index = 0; feature_index < NumDefLineModifiers (); feature_index++)
    {
      if (GetStatus (dlfp->modifier_btns[feature_index]))
      {
        ValNodeAddInt (&modifier_indices, 0, feature_index);
      }
    }
  }

  feature_list_type = GetValue (dlfp->featurePopup);
  switch (feature_list_type)
  {
    case DEFLINE_USE_FEATURES :
      dlfp->feature_requests.feature_list_type = DEFLINE_USE_FEATURES;
      break;
    case DEFLINE_PARTIAL_SEQUENCE :
      dlfp->feature_requests.feature_list_type = DEFLINE_PARTIAL_SEQUENCE;
      break;
    case DEFLINE_COMPLETE_SEQUENCE :
      dlfp->feature_requests.feature_list_type = DEFLINE_COMPLETE_SEQUENCE;
      break;
    case DEFLINE_COMPLETE_GENOME :
      dlfp->feature_requests.feature_list_type = DEFLINE_COMPLETE_GENOME;
      break;
    case DEFLINE_PARTIAL_GENOME :
      dlfp->feature_requests.feature_list_type = DEFLINE_PARTIAL_GENOME;
      break;
    case DEFLINE_SEQUENCE :
      dlfp->feature_requests.feature_list_type = DEFLINE_SEQUENCE;
      break;
    default:
      dlfp->feature_requests.feature_list_type = DEFLINE_USE_FEATURES;
      break;
  }

  for (i=0; i< NumRemovableItems; i++)
  {
    dlfp->feature_requests.keep_items[i] 
              = GetStatus (dlfp->feature_btns[i]);
  }
  if (GetValue (dlfp->promoter_type) == 1) 
  {
    dlfp->feature_requests.add_fake_promoters = TRUE;
  }
  else
  {
    dlfp->feature_requests.add_fake_promoters = FALSE;
  }

  dlfp->feature_requests.suppress_alt_splice_phrase = 
                 GetStatus (dlfp->suppress_alt_splice_phrase);

  dlfp->feature_requests.use_ncrna_note = 
                 GetStatus (dlfp->use_ncrna_note);
  
  dlfp->feature_requests.remove_subfeatures = 
                 GetStatus (dlfp->remove_subfeatures);

  dlfp->feature_requests.suppress_locus_tags = 
                 GetStatus (dlfp->suppress_locus_tags);

  dlfp->feature_requests.misc_feat_parse_rule = 
                 GetValue (dlfp->misc_feat_parse_rule);
                 
  dlfp->feature_requests.suppressed_feature_list = DialogToPointer (dlfp->suppressed_feature_list);                 

  product_flag = GetValue (dlfp->organelle_popup) - 1;
  alternate_splice_flag = GetStatus (dlfp->alternate_splice_flag);
  gene_cluster_opp_strand = GetStatus (dlfp->gene_cluster_opp_strand);
 
  if (dlfp->target_bsp != NULL && GetStatus (dlfp->modify_only_target))
  {
    sep = GetBestTopParentForData (dlfp->input_entityID, dlfp->target_bsp);
  }
  else
  { 
    sep = GetTopSeqEntryForEntityID (dlfp->input_entityID);
  }
  if (sep == NULL) return;

  RemovePopsetTitles (sep);
  if (dlfp->popset_only) {
    AddPopsetTitles (sep, &dlfp->feature_requests, product_flag, alternate_splice_flag, gene_cluster_opp_strand);
  } else {
    AutoDefForSeqEntry (sep, dlfp->input_entityID, &odmp, dlfp->modList, modifier_indices,
                        &dlfp->feature_requests, product_flag, alternate_splice_flag, gene_cluster_opp_strand);
  }
                              
  dlfp->feature_requests.suppressed_feature_list = ValNodeFree (dlfp->feature_requests.suppressed_feature_list);
  modifier_indices = ValNodeFree (modifier_indices);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (dlfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlfp->input_entityID, 0, 0);
  Remove (dlfp->form);
  SeqEntrySetScope (NULL);
}

static void IsTricky (
  BioSourcePtr biop,
  Pointer userdata
)
{
  Boolean PNTR pTricky;
  static CharPtr long_name = "Human immunodeficiency virus";
  static CharPtr short_name = "HIV-";
  OrgModPtr     mod;
  OrgNamePtr    onp;
  SubSourcePtr  ssp;

  pTricky = (Boolean PNTR) userdata;
  if (pTricky == NULL
    || *pTricky
    || biop == NULL
    || biop->org == NULL
    || biop->org->taxname == NULL
    || biop->org->orgname == NULL
    || biop->subtype == NULL)
  {
    return;
  }
  if (StringNICmp (biop->org->taxname, long_name, StringLen (long_name)) != 0
    && StringNICmp (biop->org->taxname, short_name, StringLen (short_name)) != 0)
  {
    return;
  }

  onp = biop->org->orgname;
  if (onp == NULL) return;
  mod = onp->mod;
  while (mod != NULL
    && mod->subtype != ORGMOD_isolate)
  {
    mod = mod->next;
  }
  if (mod == NULL || mod->subname == NULL)
  {
    return;
  }

  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != SUBSRC_clone)
  {
    ssp = ssp->next;
  }
  if (ssp == NULL || ssp->name == NULL)
  {
    return;
  }
  *pTricky = TRUE;
  return;
}

static Boolean HasTrickyHIVRecords (
  SeqEntryPtr sep
)
{
  Boolean     has_tricky;

  if (sep == NULL) return FALSE;
  
  has_tricky = FALSE;
  VisitBioSourcesInSep (sep, &has_tricky, IsTricky);
  return has_tricky;
}


static GrouP CreateDefLineFormHIVRule (
  GrouP h
)
{
  GrouP hiv_rule;

  hiv_rule = NormalGroup (h, 3, 0, "HIV rule", programFont, NULL);
  RadioButton (hiv_rule, "Prefer Clone");
  RadioButton (hiv_rule, "Prefer Isolate");
  RadioButton (hiv_rule, "Want Both Isolate and Clone");
  SetValue (hiv_rule, clone_isolate_HIV_rule_want_both);
  return hiv_rule;
}

static PopuP CreateDefLineFormModLimitPopup (
  GrouP q,
  Int2 count
)
{
  PopuP limit;
  Int2  i;
  Char  str[10];

  StaticPrompt (q, "Max modifiers per line", 0, popupMenuHeight, programFont, 'l');
  limit = PopupList (q, TRUE, NULL);
  PopupItem (limit, "no limit");
  for (i = 1; i <= count; i++) {
    sprintf (str, "%d", (int) i);
    PopupItem (limit, str);
  }

  SetValue (limit, 1);
  
  return limit;
}

static PopuP CreateDefLineFormModifierListPopuP (
  GrouP          q,
  DefLineFormPtr dlfp
)
{
  PopuP popup;
  StaticPrompt (q, "Modifier List", 0, popupMenuHeight, programFont, 'l');
  popup = PopupList (q, TRUE, ChangeCustomPopup);
  SetObjectExtra (popup, dlfp, NULL);
  PopupItem (popup, "All");
  PopupItem (popup, "Custom");
  SetValue (popup, 2);
  
  return popup;
}
  
static void SetHIVRuleEnable (Handle a)
{
  DefLineFormPtr dlfp;
 
  dlfp = GetObjectExtra (a);
  if (dlfp == NULL) return;
  if (dlfp->clone_isolate_HIV_rule_num == NULL) return;
  if (dlfp->modList == NULL) return;
  if ((dlfp->modifier_btns [DEFLINE_POS_Clone] != NULL
       && GetStatus (dlfp->modifier_btns [DEFLINE_POS_Clone]))
      || (dlfp->modifier_btns [DEFLINE_POS_Isolate] != NULL
          && GetStatus (dlfp->modifier_btns [DEFLINE_POS_Isolate])))
  {
    Disable (dlfp->clone_isolate_HIV_rule_num);
  }
  else
  {
    Enable (dlfp->clone_isolate_HIV_rule_num);
  }
}

static GrouP CreateDefLineFormSourceGroup (
  GrouP          h,
  DefLineFormPtr dlfp,
  SeqEntryPtr    sep
)
{
  GrouP p, g, g1, r;
  Int2    item_index, num_label_columns, num_item_rows;
  Int2    num_item_columns, count = 0;
  Char    ch;
  Boolean has_tricky;

  p = NormalGroup (h, -1, 0, "SOURCE", programFont, NULL);
  SetGroupSpacing (p, 10, 10);

  count = 0;
  for (item_index=0; item_index < NumDefLineModifiers (); item_index++)
  {
    if (dlfp->modList[item_index].any_present)
    {
      count++;
    }
  }
  
  g = HiddenGroup (p, 4, 0, NULL);
  SetGroupSpacing (g, 20, 10);
  dlfp->customGrp = CreateDefLineFormModifierListPopuP (g, dlfp);
  dlfp->modLimit = CreateDefLineFormModLimitPopup (g, count);
  dlfp->use_labels = CheckBox (p, "Use labels", NULL);
  SetStatus (dlfp->use_labels, TRUE);
  
  num_item_rows = 6;
  if (count > 18)
    num_item_rows = 8;
  if (count > 24)
    num_item_rows = 10;
  if (count > 30)
    num_item_rows = 12;

  num_item_columns = count / num_item_rows;
  if (count % num_item_rows != 0) num_item_columns ++;

  num_label_columns = 3;
  if (count > 6)
    num_label_columns --;
  if (count > 12)
    num_label_columns --;

  dlfp->sourceListGrp = NormalGroup (p,
                                     0,
                                     4,
                                     "Available Modifiers",
                                     programFont, NULL);
  SetGroupSpacing (dlfp->sourceListGrp, 10, 10);                                     

  for (item_index=0; item_index < NumDefLineModifiers (); item_index++)
  {
    if (dlfp->modList[item_index].any_present)
    {
      g1 = HiddenGroup (dlfp->sourceListGrp, num_label_columns, 0, NULL);
      SetGroupSpacing (g1, 10, 10);
      dlfp->modifier_btns[item_index] = CheckBox (g1,
                                             DefLineModifiers[item_index].name,
                                             (BtnActnProc) SetHIVRuleEnable);
      SetObjectExtra (dlfp->modifier_btns [item_index], dlfp, NULL);
      SetStatus (dlfp->modifier_btns [item_index],
                 dlfp->modList[item_index].required);
      
      if (num_label_columns > 1)
      {
        StaticPrompt (g1,
                      dlfp->modList[item_index].status,
                      0, popupMenuHeight, programFont, 'l');
        if (num_label_columns > 2)
        {
		  if (StringLen (dlfp->modList[item_index].first_value_seen) > 50)
		  {
		    ch = dlfp->modList[item_index].first_value_seen [50];
            dlfp->modList[item_index].first_value_seen [50] = 0;
		  }
		  else
		  {
		    ch = 0;
		  }
          StaticPrompt (g1,
                        dlfp->modList[item_index].first_value_seen,
                        0, popupMenuHeight, programFont, 'l');
	      if (ch != 0)
		  {
		    dlfp->modList[item_index].first_value_seen [50] = ch;
		  }
        }
      }
    }
    else
    {
      dlfp->modifier_btns[item_index] = NULL;
      dlfp->modList[item_index].required = FALSE;
    }
  }

  Enable (dlfp->sourceListGrp);


  has_tricky = HasTrickyHIVRecords (sep);
  if (has_tricky)
  {
    dlfp->clone_isolate_HIV_rule_num = CreateDefLineFormHIVRule (p);
    SetObjectExtra (dlfp->clone_isolate_HIV_rule_num, dlfp, NULL);
    SetHIVRuleEnable (dlfp->clone_isolate_HIV_rule_num);
  }
  else
  {
    dlfp->clone_isolate_HIV_rule_num = NULL;
  }
  r = NormalGroup (p, 3, 0, "Other Options", programFont, NULL);
  SetGroupSpacing (r, 10, 10);
  dlfp->keep_paren = CheckBox (r, "Leave in parenthetical organism info", NULL);
  SetStatus (dlfp->keep_paren, TRUE);
  dlfp->include_country_extra = CheckBox (r, "Include text after colon in country", NULL);
  SetStatus (dlfp->include_country_extra, FALSE);
  dlfp->allow_semicolon_in_modifier = CheckBox (r, "Include text after semicolon in modifiers", NULL);
  SetStatus (dlfp->allow_semicolon_in_modifier, FALSE);
  dlfp->exclude_sp = CheckBox (r, "Do not apply modifier to 'sp.' organisms",
         NULL);
  SetStatus (dlfp->exclude_sp, ShouldExcludeSp (sep));
  dlfp->exclude_cf = CheckBox (r, "Do not apply modifier to 'cf.' organisms",
         NULL);
  SetStatus (dlfp->exclude_cf, FALSE);
  
  dlfp->exclude_aff = CheckBox (r, "Do not apply modifier to 'aff.' organisms",
         NULL);
  SetStatus (dlfp->exclude_aff, FALSE);
  dlfp->exclude_nr = CheckBox (r, "Do not apply modifier to 'nr.' organisms", NULL);
  SetStatus (dlfp->exclude_nr, FALSE);
  
  dlfp->allow_mod_at_end_of_taxname = CheckBox (r, "Do not apply modifier to organisms with matching tax name value", NULL);
  SetStatus (dlfp->allow_mod_at_end_of_taxname, TRUE);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlfp->use_labels,
                              (HANDLE) dlfp->sourceListGrp,
                              (HANDLE) r,
                              (HANDLE) dlfp->clone_isolate_HIV_rule_num,
                              NULL);
  return p;
}

static PopuP CreateDefLineFormOrganellePopup (
  GrouP h
)
{
  PopuP o;

  o = PopupList (h, FALSE, NULL);
  PopupItem (o, "No mitochondrial or chloroplast suffix");
  PopupItem (o, "Nuclear gene(s) for mitochondrial product(s)");
  PopupItem (o, "Nuclear gene(s) for chloroplast product(s)");
  PopupItem (o, "Nuclear gene(s) for kinetoplast product(s)");
  PopupItem (o, "Nuclear gene(s) for plastid product(s)");
  PopupItem (o, "Nuclear gene(s) for chromoplast product(s)");
  PopupItem (o, "Nuclear gene(s) for cyanelle product(s)");
  PopupItem (o, "Nuclear gene(s) for apicoplast product(s)");
  PopupItem (o, "Nuclear gene(s) for leucoplast product(s)");
  PopupItem (o, "Nuclear gene(s) for proplastid product(s)");
  PopupItem (o, "Nuclear genes based on CDS products");
  SetValue (o, DEFAULT_ORGANELLE_CLAUSE + 1);
  return o;
}

static GrouP CreateDefLineFormFeatureListPopuP (
  GrouP          g,
  DefLineFormPtr dlfp
)
{
  GrouP q;

  q = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (q, "Features or Complete", 0, popupMenuHeight, programFont, 'l');
  dlfp->featurePopup = PopupList (q, TRUE, ChangeFeaturePopup);
  SetObjectExtra (dlfp->featurePopup, dlfp, NULL);
  PopupItem (dlfp->featurePopup, "List Features");
  PopupItem (dlfp->featurePopup, "Complete Sequence");
  PopupItem (dlfp->featurePopup, "Partial Sequence");
  PopupItem (dlfp->featurePopup, "Complete Genome");
  PopupItem (dlfp->featurePopup, "Partial Genome");
  PopupItem (dlfp->featurePopup, "Sequence");
  SetValue (dlfp->featurePopup, 1);
  
  return q;
}

static void SetMiscFeatRuleEnable (Handle a)
{
  DefLineFormPtr dlfp;

  if (a == NULL || ( dlfp = (DefLineFormPtr) GetObjectExtra (a)) == NULL)
  {
    return;
  }
  if (GetStatus (dlfp->feature_btns[RemovableNoncodingProductFeat]))
  {
    Enable (dlfp->misc_feat_parse_rule);
  }
  else
  {
    Disable (dlfp->misc_feat_parse_rule);
  }
  if (GetStatus (dlfp->feature_btns[RemovablePromoter]))
  {
    Enable (dlfp->promoter_type);
  }
  else
  {
    Disable (dlfp->promoter_type);
  }
}
  
static GrouP CreateDefLineFormFeatureOptionsGroup (
  GrouP          h,
  DefLineFormPtr dlfp
)
{
  GrouP p, g1, s1, s2, k, k2;
  Int4  i;
  SeqEntryPtr sep;

  p = NormalGroup (h, -1, 0, "FEATURES", programFont, NULL);
  SetGroupSpacing (p, 10, 10);

  g1 = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g1, 10, 10);
  s1 = HiddenGroup (g1, -1, 0, NULL);
  SetGroupSpacing (s1, 10, 10);
  CreateDefLineFormFeatureListPopuP (s1, dlfp);
  dlfp->organelle_popup = CreateDefLineFormOrganellePopup (s1);
  dlfp->alternate_splice_flag = CheckBox (s1,
	                "Append 'alternatively spliced'", NULL);
  dlfp->use_ncrna_note = CheckBox (s1, "Use ncRNA note if no class or product", NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlfp->featurePopup, 
                              (HANDLE) dlfp->organelle_popup, 
                              (HANDLE) dlfp->use_ncrna_note, 
                              (HANDLE) dlfp->alternate_splice_flag, NULL);
  
  dlfp->optional_features_grp = NormalGroup (g1, 3, 0,
                   "Optional Features", programFont, NULL);
                   SetGroupSpacing (dlfp->optional_features_grp, 10, 10);
  for (i=0; i < NumRemovableItems; i++)
  {
    if (i == 0 || i == 4 || i == 7)
    {
      k = HiddenGroup (dlfp->optional_features_grp, 0, 4, NULL);
      SetGroupSpacing (k, 10, 10);
    }
    if (i != RemovableCDS)
    {
      dlfp->feature_btns[i] = 
                     CheckBox (k, GetRemovableItemName (i), 
                               (BtnActnProc) SetMiscFeatRuleEnable);
      SetObjectExtra (dlfp->feature_btns[i], dlfp, NULL);
      SetStatus (dlfp->feature_btns[i], FALSE);
      if (i == RemovableNoncodingProductFeat) {
        s2 = HiddenGroup (k, 2, 0, NULL);
        StaticPrompt (s2, "        ", 8, popupMenuHeight, programFont, 'l');
        dlfp->misc_feat_parse_rule = PopupList (s2, TRUE, NULL);
        PopupItem (dlfp->misc_feat_parse_rule, "Use comment before first semicolon");
        PopupItem (dlfp->misc_feat_parse_rule, "Look for Noncoding Products");
        SetValue (dlfp->misc_feat_parse_rule, 2);
        Disable (dlfp->misc_feat_parse_rule);
      } else if (i == RemovablePromoter) {
        k2 = HiddenGroup (k, 2, 0, NULL);
        SetGroupSpacing (k2, 10, 10);
        StaticPrompt (k2, "        ", 8, popupMenuHeight, programFont, 'l');
        dlfp->promoter_type = HiddenGroup (k2, 0, 2, NULL);
        RadioButton (dlfp->promoter_type, "All");
        RadioButton (dlfp->promoter_type, "If present");
        SetValue (dlfp->promoter_type, 1);
        Disable (dlfp->promoter_type);
        AlignObjects (ALIGN_CENTER, (HANDLE) dlfp->feature_btns[i],
                                    (HANDLE) dlfp->promoter_type,
                                    NULL);
      }
    }
  }
                                      
  dlfp->suppressed_feature_grp = NormalGroup (g1, -1, 0,
                    "Suppress Features", programFont, NULL);
  sep = GetTopSeqEntryForEntityID(dlfp->input_entityID);
  dlfp->suppressed_feature_list = FeatureSelectionDialogEx (dlfp->suppressed_feature_grp, TRUE, sep, NULL, NULL);
  
  dlfp->featureOptsGrp = NormalGroup (g1, 0, 4, "Suppress", programFont, NULL);
  SetGroupSpacing (dlfp->featureOptsGrp, 10, 10);
  dlfp->remove_subfeatures = CheckBox (dlfp->featureOptsGrp, "Mobile element subfeatures", NULL);
  SetStatus (dlfp->remove_subfeatures, FALSE);
  
  dlfp->gene_cluster_opp_strand = CheckBox (dlfp->featureOptsGrp,
             "Gene cluster/locus subfeatures (both strands)", NULL);
  SetStatus (dlfp->gene_cluster_opp_strand, FALSE);

  dlfp->suppress_locus_tags = CheckBox (dlfp->featureOptsGrp, "Locus tags", NULL);
  SetStatus (dlfp->suppress_locus_tags, FALSE);

  dlfp->suppress_alt_splice_phrase = CheckBox (dlfp->featureOptsGrp, "Alternative splice phrase", NULL);
  SetStatus (dlfp->suppress_alt_splice_phrase, FALSE);

  return p;
}

static void CreateDefLineForm (
  BaseFormPtr bfp,
  ModifierItemLocalPtr modList,
  DeflineFeatureRequestList *feature_requests,
  Boolean popset_only
)
{
  DefLineFormPtr dlfp;
  WindoW         w;
  GrouP          h, k, g1 = NULL, c, feat_opts;
  ButtoN         b;
  SeqEntryPtr    sep;

  dlfp = MemNew (sizeof (DefLineFormData));
  if (dlfp == NULL) return;
  dlfp->input_entityID = bfp->input_entityID;
  dlfp->input_itemID = bfp->input_itemID;
  dlfp->input_itemtype = bfp->input_itemtype;
  dlfp->popset_only = popset_only;

  dlfp->modifier_btns = (ButtoN PNTR) MemNew (sizeof (ButtoN) * NumDefLineModifiers());

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  w = FixedWindow (-50, -33, -10, -10, "Automatic Definition Line",
  StdCloseWindowProc);
  SetObjectExtra (w, dlfp, CleanupDefLineForm);
  dlfp->form = (ForM) w;
  dlfp->formmessage = DefLineFormMessageProc;

  dlfp->modList = modList;
  dlfp->feature_requests = *feature_requests;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  k = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  if (!popset_only) {
    g1 = CreateDefLineFormSourceGroup (k, dlfp, sep);
  }

  feat_opts = CreateDefLineFormFeatureOptionsGroup (k, dlfp);
  AlignObjects (ALIGN_CENTER, (HANDLE) feat_opts, (HANDLE) g1, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoAutoDefLine);
  SetObjectExtra (b, dlfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  dlfp->target_bsp = GetBioseqGivenIDs (dlfp->input_entityID,
          dlfp->input_itemID,
          dlfp->input_itemtype);
  dlfp->modify_only_target = CheckBox (c, "Only modify targeted record", NULL);
  SetStatus (dlfp->modify_only_target, FALSE);
  if (dlfp->target_bsp == NULL)
  {
    Disable (dlfp->modify_only_target);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) k,
                (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


/* The AddBestComboModifiersToModList function sets the modifiers
 * selected by the FindBestCombo process as selected.
 */
static void AddBestModifiersToModList (
  ValNodePtr modifier_indices,
  ModifierItemLocalPtr modList
)
{
  ValNodePtr vnp;

  for (vnp = modifier_indices; vnp != NULL; vnp = vnp->next)
  {
    modList[vnp->data.intvalue].required = TRUE;
  }
}


extern void AutoDefEntityIDNoOptions (Uint2 entityID, Boolean use_modifiers)
{
  SeqEntryPtr sep;
  DeflineFeatureRequestList feature_requests;
  ModifierItemLocalPtr modList;
  ValNodePtr modifier_indices = NULL;
  OrganismDescriptionModifiers odmp;
  Int4 index;
  ValNodePtr defline_clauses = NULL;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;

  InitFeatureRequests (&feature_requests);

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;
  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  InitOrganismDescriptionModifiers (&odmp, sep);
  odmp.use_modifiers = use_modifiers;

  RemoveNucProtSetTitles (sep);  
  SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID, &feature_requests,
                              DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                              &defline_clauses);
  if ( AreFeatureClausesUnique (defline_clauses))
  {
    modifier_indices = GetModifierIndicesFromModList (modList);
  }
  else
  {
    modifier_indices = FindBestModifiersForDeflineClauseList (defline_clauses, modList);
  }

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, &odmp);

  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL)
  {
    for (index=0; index < NumDefLineModifiers (); index++)
    {
      ValNodeFree (modList[index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
  RemovePopsetTitles (sep);
  AddPopsetTitles (sep, &feature_requests, DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);
}


extern void AutoDefBaseFormCommon (
  BaseFormPtr bfp,
  Boolean use_form,
  Boolean use_modifiers,
  Boolean popset_only
)
{
  SeqEntryPtr sep;
  DeflineFeatureRequestList feature_requests;
  ModifierItemLocalPtr modList;
  ValNodePtr modifier_indices = NULL;
  OrganismDescriptionModifiers odmp;
  Int4 index;
  ValNodePtr defline_clauses = NULL;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  InitFeatureRequests (&feature_requests);

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;
  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);
  if (use_form)
  {
    modifier_indices = FindBestModifiersEx (sep, modList, TRUE);
    AddBestModifiersToModList ( modifier_indices, modList);
    modifier_indices = ValNodeFree (modifier_indices);
    CreateDefLineForm (bfp, modList, &feature_requests, popset_only);
  }
  else
  {
    WatchCursor ();
    Update ();
    InitOrganismDescriptionModifiers (&odmp, sep);
    odmp.use_modifiers = use_modifiers;

    RemoveNucProtSetTitles (sep);  
    SeqEntrySetScope (sep);

    BuildDefLineFeatClauseList (sep, bfp->input_entityID, &feature_requests,
                                DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                                &defline_clauses);
    if ( AreFeatureClausesUnique (defline_clauses))
    {
      modifier_indices = GetModifierIndicesFromModList (modList);
    }
    else
    {
      modifier_indices = FindBestModifiersForDeflineClauseList (defline_clauses, modList);
    }

    if (!popset_only) {
      BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                                  modifier_indices, &odmp);
    }

    DefLineFeatClauseListFree (defline_clauses);
    if (modList != NULL)
    {
      for (index=0; index < NumDefLineModifiers (); index++)
      {
        ValNodeFree (modList[index].values_seen);
      }
      MemFree (modList);
    }
    modifier_indices = ValNodeFree (modifier_indices);

    if (!popset_only) {
      ClearProteinTitlesInNucProts (bfp->input_entityID, NULL);
      InstantiateProteinTitles (bfp->input_entityID, NULL);
    }
    RemovePopsetTitles (sep);
    AddPopsetTitles (sep, &feature_requests, DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);

    ArrowCursor ();
    Update ();
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);  
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    SeqEntrySetScope (NULL);
  }
}


static const Int4 s_auto_def_id_preferred_quals[] = {
  DEFLINE_POS_Strain,
  DEFLINE_POS_Clone,
  DEFLINE_POS_Isolate,
  DEFLINE_POS_Cultivar,
  DEFLINE_POS_Specimen_voucher,
};

static const Int4 k_num_auto_def_id_preferred_quals = sizeof (s_auto_def_id_preferred_quals) / sizeof (Int4);

extern void AutoDefIdEx (Uint2 entityID, DefLineType feature_list_type)
{
  SeqEntryPtr sep;
  DeflineFeatureRequestList feature_requests;
  ModifierItemLocalPtr modList;
  ValNodePtr modifier_indices = NULL;
  OrganismDescriptionModifiers odmp;
  Int4 index;
  Boolean added_required = FALSE;
  ValNodePtr defline_clauses = NULL;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;

  InitFeatureRequests (&feature_requests);
  feature_requests.feature_list_type = feature_list_type;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;
  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);
  /* first look for first modifier in list that is present on all sources */
  for (index = 0; index < k_num_auto_def_id_preferred_quals && !added_required; index++) {
    if (modList[s_auto_def_id_preferred_quals[index]].all_present) {
      modList[s_auto_def_id_preferred_quals[index]].required = TRUE;
      added_required = TRUE;
    }
  }
  /* if not found, then look for first modifier in list that is present on any sources */
  for (index = 0; index < k_num_auto_def_id_preferred_quals && !added_required; index++) {
    if (modList[s_auto_def_id_preferred_quals[index]].any_present) {
      modList[s_auto_def_id_preferred_quals[index]].required = TRUE;
      added_required = TRUE;
    }
  }

  WatchCursor ();
  Update ();
  InitOrganismDescriptionModifiers (&odmp, sep);
  odmp.use_modifiers = TRUE;

  RemoveNucProtSetTitles (sep);  
  SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID, &feature_requests,
                              DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                              &defline_clauses);
  if ( AreFeatureClausesUnique (defline_clauses))
  {
    modifier_indices = GetModifierIndicesFromModList (modList);
  }
  else
  {
    modifier_indices = FindBestModifiersForDeflineClauseList (defline_clauses, modList);
  }

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, &odmp);

  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL)
  {
    for (index=0; index < NumDefLineModifiers (); index++)
    {
      ValNodeFree (modList[index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
  RemovePopsetTitles (sep);
  AddPopsetTitles (sep, &feature_requests, DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);

}


extern void AutoDefId (Uint2 entityID)
{
  AutoDefIdEx (entityID, DEFLINE_USE_FEATURES);
}


extern void AutoDefStrain (
  BaseFormPtr bfp
)
{

  if (bfp == NULL) return;

  WatchCursor ();
  Update ();
  AutoDefId (bfp->input_entityID);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);  
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  SeqEntrySetScope (NULL);
}


extern void AutoDefMiscFeat (
  BaseFormPtr bfp
)
{
  SeqEntryPtr sep;
  DeflineFeatureRequestList feature_requests;
  ModifierItemLocalPtr modList;
  ValNodePtr modifier_indices = NULL;
  OrganismDescriptionModifiers odmp;
  Int4 index;
  ValNodePtr defline_clauses = NULL;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  InitFeatureRequests (&feature_requests);
  feature_requests.keep_items[RemovableNoncodingProductFeat] = TRUE;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;
  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  WatchCursor ();
  Update ();
  InitOrganismDescriptionModifiers (&odmp, sep);
  odmp.use_modifiers = TRUE;

  RemoveNucProtSetTitles (sep);  
  SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, bfp->input_entityID, &feature_requests,
                              DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                              &defline_clauses);
  if ( AreFeatureClausesUnique (defline_clauses))
  {
    modifier_indices = GetModifierIndicesFromModList (modList);
  }
  else
  {
    modifier_indices = FindBestModifiersForDeflineClauseList (defline_clauses, modList);
  }

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, &odmp);

  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL)
  {
    for (index=0; index < NumDefLineModifiers (); index++)
    {
      ValNodeFree (modList[index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  ClearProteinTitlesInNucProts (bfp->input_entityID, NULL);
  InstantiateProteinTitles (bfp->input_entityID, NULL);
  RemovePopsetTitles (sep);
  AddPopsetTitles (sep, &feature_requests, DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);  
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  SeqEntrySetScope (NULL);
}


static void AutoDefCommon (IteM i, Boolean use_form, Boolean use_modifiers)
{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;
  AutoDefBaseFormCommon (bfp, use_form, use_modifiers, FALSE);
}

extern void AutoDef (IteM i)
{
  AutoDefCommon (i, FALSE, TRUE);
}

extern void AutoDefWithOptions (IteM i)
{
  AutoDefCommon (i, TRUE, TRUE);
}


extern void AutoDefWithoutModifiers (IteM i)
{
  AutoDefCommon (i, FALSE, FALSE);
}

extern void AutoDefToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  AutoDefBaseFormCommon (bfp, FALSE, TRUE, FALSE);
}

extern void AutoDefOptionsToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  AutoDefBaseFormCommon (bfp, TRUE, TRUE, FALSE);
}


extern void AutoDefStrainToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  AutoDefStrain (bfp);
}


extern void AutoDefMiscFeatToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  AutoDefMiscFeat (bfp);
}


extern void AddPopsetTitlesItem (IteM i)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (i);
  if (bfp == NULL) return;
  AutoDefBaseFormCommon (bfp, TRUE, TRUE, TRUE);
}


static void TrimQuotesAroundString (CharPtr str)
{
  Int4 len;
  CharPtr src, dst;

  if (StringHasNoText (str) || *str != '"') return;
  len = StringLen (str);
  if (str[len - 1] != '"') return;
  src = str + 1;
  dst = str;
  while (*src != '"')
  {
    *dst = *src;
    dst++;
    src++;
  }
  *dst = 0;
}


typedef struct tablelinedata {
  Int4 num_parts;
  ValNodePtr parts;
} TableLineData, PNTR TableLinePtr;

static TableLinePtr MakeTableData (CharPtr line, ValNodePtr PNTR special_list)
{
  TableLinePtr tlp;
  CharPtr p_start, p_end;
  Int4    plen;
  Boolean found_end;
  CharPtr val;
  ValNodePtr vnp;

  if (line == NULL) return NULL;
  tlp = MemNew (sizeof (TableLineData));
  if (tlp == NULL) return NULL;
  p_start = line;
  found_end = FALSE;
  while (*p_start != 0 && !found_end)
  {
    plen = StringCSpn (p_start, "\t\n");
    if (plen == 0)
    {
      ValNodeAddStr (&tlp->parts, 0, StringSave (""));
      tlp->num_parts ++;
      p_start++;
      if (*p_start == 0) {
        if (tlp->num_parts > 0)
        {
          ValNodeAddStr (&tlp->parts, 0, StringSave (""));
          tlp->num_parts ++;
        }
      }
      continue;
    }
    if (plen == StringLen (p_start))
    {
      found_end = TRUE;
    }
    else
    {
      p_end = p_start + plen;
      *p_end = 0;
    }
    TrimSpacesAroundString (p_start);
    TrimQuotesAroundString (p_start);
    val = StringSave (p_start);

    vnp = ValNodeAddStr (&tlp->parts, 0, val);
    SpecialCharFindWithContext ((CharPtr PNTR)&(vnp->data.ptrvalue), special_list, NULL, NULL);
    tlp->num_parts ++;
    if (!found_end)
    {
      p_start = p_end + 1;
    }
  }
  if (tlp->num_parts == 0)
  {
    MemFree (tlp);
    return NULL;
  } 
  return tlp;  
}

static void CleanUpTableData (ValNodePtr vnp)
{
  TableLinePtr tlp;
  ValNodePtr   vnp_next;

  while (vnp != NULL)
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    tlp = vnp->data.ptrvalue;
    if (tlp != NULL)
    {
      ValNodeFreeData (tlp->parts);
    }
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
}

static void FormatPopupWithTableDataColumns (PopuP p, ValNodePtr header_line)
{
  TableLinePtr tlp;
  ValNodePtr   vnp;
  CharPtr      val_fmt = "Column %d (%s)";
  CharPtr      val;
  Int4         col = 1;

  if (p == NULL || header_line == NULL) return;
  tlp = (TableLinePtr) header_line->data.ptrvalue;
  while ((tlp == NULL || tlp->num_parts < 1) && header_line->next != NULL) {
    header_line = header_line->next;
    tlp = header_line->data.ptrvalue;
  }
  if (tlp == NULL || tlp->num_parts < 1) return;
  vnp = tlp->parts;
  while (vnp != NULL) {
    val = MemNew ((StringLen(val_fmt) + StringLen (vnp->data.ptrvalue) + 15) * sizeof (Char));
    sprintf (val, val_fmt, col, vnp->data.ptrvalue == NULL ? "" : vnp->data.ptrvalue);
    PopupItem (p, val);
    MemFree (val);
    vnp = vnp->next;
    col++;
  }
}


enum orgmodmatchtype 
{
  eMatchAccession = 1,
  eMatchLocalID,
  eMatchTaxName,
  eMatchTMSMART,
  eMatchBankIt,
  eMatchGeneral,
  eMatchNone
};


typedef struct orgmodtablecolumn {
  Int4       match_choice;
  ValNodePtr apply_choice;
} OrgModTableColumnData, PNTR OrgModTableColumnPtr;


static OrgModTableColumnPtr OrgModTableColumnFree (OrgModTableColumnPtr t)
{
  if (t != NULL) {
    t->apply_choice = ValNodeFreeData (t->apply_choice);
    t = MemFree (t);
  }
  return t;
}


static ValNodePtr OrgModTableColumnListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) 
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = OrgModTableColumnFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


typedef struct orgmodtablecolumndlg {
  DIALOG_MESSAGE_BLOCK

  GrouP  action_choice;
  PopuP  match_choice;
  DialoG apply_choice;
 
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} OrgModTableColumnDlgData, PNTR OrgModTableColumnDlgPtr;

static Pointer OrgModTableColumnDialogToData (DialoG d)
{
  OrgModTableColumnPtr o;
  OrgModTableColumnDlgPtr dlg;
  
  dlg = (OrgModTableColumnDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  o = (OrgModTableColumnPtr) MemNew (sizeof (OrgModTableColumnData));

  if (GetValue (dlg->action_choice) == 1) 
  {
    o->match_choice = GetValue (dlg->match_choice);
    o->apply_choice = NULL;
  }
  else
  {
    o->match_choice = 0;
    o->apply_choice = DialogToPointer (dlg->apply_choice);
  }
  return o;
}


static void ChangeOrgModTableColumnMatchChoice (PopuP p)
{
  OrgModTableColumnDlgPtr dlg;

  dlg = (OrgModTableColumnDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;
  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify) (dlg->change_userdata);
  }  
}
  

static void ChangeOrgModTableColumnActionChoice (GrouP g)
{
  OrgModTableColumnDlgPtr dlg;

  dlg = (OrgModTableColumnDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;

  if (GetValue (dlg->action_choice) == 1)
  {
    Enable (dlg->match_choice);
    Disable (dlg->apply_choice);
  }
  else
  {
    Enable (dlg->apply_choice);
    Disable (dlg->match_choice);
  }  
  if (dlg->change_notify != NULL) 
  {
    (dlg->change_notify) (dlg->change_userdata);
  }  
}


static void OrgModTableColumnDataToDialog (DialoG d, Pointer data)
{
  OrgModTableColumnPtr o;
  OrgModTableColumnDlgPtr dlg;
  
  dlg = (OrgModTableColumnDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  o = (OrgModTableColumnPtr) data;

  if (o == NULL)
  {
    SetValue (dlg->action_choice, 1);
    SetValue (dlg->match_choice, eMatchNone);
    PointerToDialog (dlg->apply_choice, NULL);
  }
  else
  {
    if (o->match_choice == 0)
    {
      SetValue (dlg->action_choice, 2);
      PointerToDialog (dlg->apply_choice, o->apply_choice);
    } 
    else
    {
      SetValue (dlg->action_choice, 1);
      SetValue (dlg->match_choice, o->match_choice);
    }
  }
  ChangeOrgModTableColumnActionChoice (dlg->action_choice);
}


static DialoG 
OrgModTableColumnDialog 
(GrouP                h,
 CharPtr              value_string,
 Int4                 num_blank,
 Nlm_ChangeNotifyProc change_notify,
 Pointer              change_userdata)
{
  OrgModTableColumnDlgPtr dlg;
  GrouP            p, g;
  ValNodePtr       qual_selection_list;
  CharPtr real_title;
  CharPtr title_fmt = "%s (%d are blank)";

  dlg = (OrgModTableColumnDlgPtr) MemNew (sizeof (OrgModTableColumnDlgData));

  if (num_blank > 0) {
    real_title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (value_string) + 15));
    sprintf (real_title, title_fmt, value_string, num_blank);
  } else {
    real_title = value_string;
  }

  p = NormalGroup (h, 3, 0, real_title, programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  if (real_title != value_string) {
    real_title = MemFree (real_title);
  }

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = OrgModTableColumnDialogToData;
  dlg->todialog = OrgModTableColumnDataToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->action_choice = HiddenGroup (p, 0, 2, ChangeOrgModTableColumnActionChoice);
  SetObjectExtra (dlg->action_choice, dlg, NULL);
  SetGroupSpacing (dlg->action_choice, 10, 10);
  RadioButton (dlg->action_choice, "Match to");
  RadioButton (dlg->action_choice, "Apply to");
  SetValue (dlg->action_choice, 1);

  g = HiddenGroup (p, 0, 3, NULL);
  /* list of matchable items */
  dlg->match_choice = PopupList (g, TRUE, ChangeOrgModTableColumnMatchChoice);
  SetObjectExtra (dlg->match_choice, dlg, NULL);
  PopupItem (dlg->match_choice, "Accession");
  PopupItem (dlg->match_choice, "Local ID");
  PopupItem (dlg->match_choice, "Organism Taxonomy Name");
  PopupItem (dlg->match_choice, "TMSMART ID");
  PopupItem (dlg->match_choice, "BankIt ID");
  PopupItem (dlg->match_choice, "General ID");
  PopupItem (dlg->match_choice, "None");
  SetValue (dlg->match_choice, eMatchNone);


  /* list of organism modifiers */
  qual_selection_list = GetSourceQualDescList (TRUE, TRUE, FALSE, FALSE);
  
  ValNodeAddPointer (&qual_selection_list, 1, StringSave ("Tax Name"));
  
    /* note - the ValNodeSelectionDialog will free the qual_selection_list when done */
  dlg->apply_choice = ValNodeSelectionDialog (g, qual_selection_list, 8, SourceQualValNodeName,
                                ValNodeSimpleDataFree, SourceQualValNodeDataCopy,
                                SourceQualValNodeMatch, "qual choice", 
                                change_notify, change_userdata, FALSE);

  ChangeOrgModTableColumnActionChoice(dlg->action_choice);
  return (DialoG) p;  
}


typedef struct orgmodloadformdata {
  FEATURE_FORM_BLOCK

  ValNodePtr        line_list;
  BioSourcePtr PNTR biop_list;
  DialoG PNTR       columns;
  Int4              num_columns;
  Uint2             entityID;
  ButtoN            replace_with_blank_btn;
  ButtoN            accept_button;
  ButtoN            leave_dlg_up_btn;
  DialoG            taxname_options;

  /* These members are for importing modifiers in the
   * submission dialogs.
   */
  SeqEntryPtr       seq_list;
  Boolean           done;
  Boolean           mods_added;
  
} OrgModLoadFormData, PNTR OrgModLoadFormPtr;

static void SetFormModsEnable (Pointer userdata)
{
  OrgModLoadFormPtr form_data;
  OrgModTableColumnPtr o;
  Boolean           have_apply = FALSE;
  Boolean           have_action = FALSE;
  Boolean           have_apply_taxname = FALSE;
  Int4              column_index;

  form_data = (OrgModLoadFormPtr) userdata;
  if (form_data == NULL) return;
 
  for (column_index = 0;
       column_index < form_data->num_columns;
       column_index ++)
  {
    o = (OrgModTableColumnPtr) DialogToPointer (form_data->columns[column_index]);
    if (o == NULL) continue;
    if (o->match_choice > 0 && o->match_choice != eMatchNone)
    {
      have_action = TRUE;
    }
    else if (o->apply_choice != NULL) 
    {
      have_apply = TRUE;
      if (o->apply_choice->choice > 0 && StringCmp (o->apply_choice->data.ptrvalue, "Tax Name") == 0)
      {
        have_apply_taxname = TRUE;
      }
    }
    o = OrgModTableColumnFree (o);
  }

  if (have_apply_taxname)
  {
    Enable (form_data->taxname_options);
  }
  else
  {
    Disable (form_data->taxname_options);
  }

  /* must have an action selected if there is more than one line */
  if ( ! have_action
    && form_data->line_list != NULL
    && form_data->line_list->next != NULL)
  {
    Disable (form_data->accept_button);
    return;
  }

  if (have_apply)
  {
    Enable (form_data->accept_button);
  }
  else
  {
    Disable (form_data->accept_button);
  }
}

static void CleanupOrgModLoadForm (
  GraphiC g,
  VoidPtr data
)
{
  OrgModLoadFormPtr form_data;

  form_data = (OrgModLoadFormPtr)data;
  if (form_data == NULL) return;
  form_data->columns = MemFree (form_data->columns);
  CleanUpTableData (form_data->line_list);
  StdCleanupFormProc (g, data);
} 

static CharPtr genome_names[] = 
{
  "",
  "genomic",
  "chloroplast",
  "chromoplast",
  "kinteoplast",
  "mitochondrion",
  "plastid",
  "macronuclear",
  "extrachrom",
  "plasmid",
  "transposon",
  "insertion_seq",
  "cyanelle",
  "proviral",
  "virion",
  "nucleomorph",
  "apicoplast",
  "leucoplast",
  "proplastid",
  "endogenous_virus",
  "hydrogenosome",
  "chromosome",
  "chromatophore"
};

static Uint1 GetGenomeValFromString (CharPtr value_string)
{
  Uint1 genome_val;
  
  if (StringHasNoText (value_string))
  {
    return 0;
  }
  for (genome_val = 1; genome_val < sizeof (genome_names); genome_val++)
  {
    if (StringICmp (value_string, genome_names[genome_val]) == 0)
    {
      return genome_val;
    }
  }
  return 0;
}


static void AddOneQualToOrg 
( BioSourcePtr biop,
  CharPtr      value_string,
  SourceQualDescPtr sqdp,
  ExistingTextPtr etp)
{
  OrgRefPtr     orp;
  OrgModPtr     mod, last_mod;
  OrgNamePtr    onp;
  SubSourcePtr  ssp, last_ssp;

  if (biop == NULL || value_string == NULL || sqdp == NULL) return;

  orp = biop->org;
  if (orp == NULL) return;

  if (sqdp->isOrgMod)
  {
    onp = orp->orgname;
    if (onp == NULL) {
      onp = OrgNameNew ();
      if (onp == NULL) return;
      orp->orgname = onp;
    }
    mod = onp->mod;
    last_mod = NULL;
    while (mod != NULL
      && mod->subtype != sqdp->subtype)
    {
      last_mod = mod;
      mod = mod->next;
    }
    if (mod != NULL)
    {
      if (StringHasNoText (value_string)) {
        if (last_mod == NULL) {
          onp->mod = mod->next;
        } else {
          last_mod->next = mod->next;
        }
        mod->next = NULL;
        mod = OrgModFree (mod);
      } else {
        mod->subname = HandleExistingText (mod->subname,
                                           StringSave (value_string),
                                           etp); 
      }
    }
    else if (!StringHasNoText (value_string))
    {
      mod = OrgModNew ();
      mod->subtype = sqdp->subtype;
      mod->subname = StringSave (value_string);
      if (last_mod == NULL)
      {
        onp->mod = mod;
      }
      else
      {
        last_mod->next = mod;
      }
    }
  } else {
    ssp = biop->subtype;
    last_ssp = NULL;
    while (ssp != NULL && ssp->subtype != sqdp->subtype)
    {
      last_ssp = ssp;
      ssp = ssp->next;
    }
    if (ssp != NULL)
    {
      if (StringHasNoText (value_string)) {
        if (last_ssp == NULL) {
          biop->subtype = ssp->next;
        } else {
          last_ssp->next = ssp->next;
        }
        ssp->next = NULL;
        ssp = SubSourceFree (ssp);
      } else {
        ssp->name = HandleExistingText (ssp->name,
                                        StringSave (value_string),
                                        etp);
      }
    }
    else if (!StringHasNoText (value_string))
    {
      ssp = SubSourceNew ();
      ssp->subtype = sqdp->subtype;
      ssp->name = StringSave (value_string);
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

static void ApplyTaxNameToOrg (BioSourcePtr biop, CharPtr value_string, ExistingTextPtr etp)
{
  if (biop == NULL) return;
  if (biop->org == NULL)
  {
  	biop->org = OrgRefNew ();
  }
  if (biop->org == NULL) return;
  biop->org->taxname = HandleExistingText (biop->org->taxname,
                                           StringSave (value_string),
                                           etp);
}

static void ApplyLocationToOrg (BioSourcePtr biop, CharPtr value_string)
{
  if (biop != NULL)
  {
    biop->genome = GetGenomeValFromString (value_string);
  }
}


static Boolean HasExtraAccession (
  CharPtr acc_str,
  GBBlockPtr gbp)
{
  ValNodePtr vnp;
  if (acc_str == NULL || gbp == NULL) return FALSE;
  for (vnp = gbp->extra_accessions;
       vnp != NULL && StringCmp (acc_str, vnp->data.ptrvalue) != 0;
       vnp = vnp->next)
  {}
  if (vnp != NULL)
    return TRUE;
  else
    return FALSE;
}

static Boolean
IDIsInTextList 
(CharPtr id,
 CharPtr acc_str,
 Boolean look_for_tmsmart)
{
  CharPtr wanted, found;
  Int4    len;

  if (id == NULL || acc_str == NULL) return FALSE;

  if (StringNCmp (acc_str, "TMSMART:", 8) == 0) {
    if (look_for_tmsmart) {
      wanted = acc_str + 8;
    } else {
      return FALSE;
    }
  } else {
    if (look_for_tmsmart) {
      return FALSE;
    } else {
      wanted = acc_str;
    }
  }
  len = StringLen (wanted);
  found = StringStr (id, wanted);
  if (found == NULL) {
    return FALSE;
  } else if (*(found + len) != 0 
             && *(found + len) != ','
             && ! isspace ((Int4)*(found + len))) {
    return FALSE;
  } else if (found != id
             && * (found - 1) != ','
             && ! isspace ((Int4)*(found - 1))) {
    return FALSE;
  } else {
    return TRUE;
  }
}

static Boolean IDListHasValue (
  CharPtr id,
  SeqIdPtr list,
  Boolean only_local,
  Boolean look_for_tmsmart,
  Boolean look_for_general)
{
  SeqIdPtr sip;
  Char     acc_str [256];
  Int4     match_len, match_len2, db_len;
  DbtagPtr gnl_tag; 

  for (sip = list; sip != NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_GENERAL 
        && look_for_general 
        && StringNICmp (id, "gnl|", 3) == 0
        && (db_len = StringCSpn (id + 4, "|")) > 0)
    {
      gnl_tag = (DbtagPtr) sip->data.ptrvalue;
      if (gnl_tag != NULL
          && StringNCmp (id + 4, gnl_tag->db, db_len) == 0
          && ((gnl_tag->tag->id > 0 && atoi (id + 5 + db_len) == gnl_tag->tag->id)
             || StringCmp (gnl_tag->tag->str, id + 5 + db_len) == 0))
      {
        return TRUE;
      }
    }
    if ((! only_local && sip->choice != SEQID_LOCAL)
      || (only_local && sip->choice == SEQID_LOCAL))
    {
      SeqIdWrite (sip, acc_str, PRINTID_REPORT, sizeof (acc_str));
      if (look_for_tmsmart)
      {
        if (StringNCmp (acc_str, "TMSMART:", 8) == 0
          && StringCmp (id, acc_str + 8) == 0)
        {
          return TRUE;
        } else if (IDIsInTextList (id, acc_str, TRUE)) {
          return TRUE;
        }
      } else {
        if (only_local)
        {
          if (StringCmp (id, acc_str) == 0)
          {
          	return TRUE;
          }
        }
        else 
        {
          match_len = StringCSpn (acc_str, ".");
          match_len2 = StringCSpn (id, ".");
          if (match_len == match_len2
            && match_len > 0 && StringNCmp (id, acc_str, match_len) == 0)
          {
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}


static CharPtr FindBankitNumberInTableString (CharPtr table_string)
{
  CharPtr bankit_start;
  
  if (StringHasNoText (table_string))
  {
    return table_string;
  }
  bankit_start = table_string + StringSpn (table_string, " ");
  if (StringNICmp (bankit_start, "bankit", 6) == 0)
  {
    bankit_start += 6;
    bankit_start += StringSpn (bankit_start, " ");
  }
  return bankit_start;
}



typedef struct applycolumn {
  ValNodePtr apply_choice;
  CharPtr    apply_value;
} ApplyColumnData, PNTR ApplyColumnPtr;


static ApplyColumnPtr ApplyColumnFree (ApplyColumnPtr p)
{
  if (p != NULL) {
    p->apply_value = MemFree (p->apply_value);
    p->apply_choice = ValNodeFreeData (p->apply_choice);
    p = MemFree (p);
  }
  return p;
}


static ValNodePtr ApplyColumnListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = ApplyColumnFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


typedef struct tabletobiosource {
  ValNodePtr identifiers;
  ValNodePtr columns_to_apply;
  ValNodePtr biop_list;
} TableToBioSourceData, PNTR TableToBioSourcePtr;


static TableToBioSourcePtr TableToBioSourceFree (TableToBioSourcePtr t)
{
  if (t != NULL) {
    t->identifiers = ValNodeFreeData (t->identifiers);
    t->columns_to_apply = ApplyColumnListFree (t->columns_to_apply);
    t->biop_list = ValNodeFree (t->biop_list);
    t = MemFree (t);
  }
  return t;
}


static ValNodePtr TableToBioSourceListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) 
  {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = TableToBioSourceFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;;
  }
  return vnp;
}


static TableToBioSourcePtr TableToBioSourceFromTableLine (OrgModTableColumnPtr PNTR o_list, Int4 num_columns, TableLinePtr tlp)
{
  Int4       column_index;
  ValNodePtr part;
  TableToBioSourcePtr t;
  ApplyColumnPtr      a;
  
  if (o_list == NULL || tlp == NULL || tlp->parts == NULL) return NULL;

  t = (TableToBioSourcePtr) MemNew (sizeof (TableToBioSourceData));
  for (column_index = 0, part = tlp->parts;
       column_index < num_columns;
       column_index ++)
  {
    if (o_list[column_index] == NULL) continue;
    if (o_list[column_index]->match_choice > 0 
        && o_list[column_index]->match_choice != eMatchNone
        && part != NULL)
    {
      ValNodeAddPointer (&(t->identifiers), o_list[column_index]->match_choice, StringSave (part->data.ptrvalue));
    }
    else if (o_list[column_index]->apply_choice != NULL)
    {
      a = (ApplyColumnPtr) MemNew (sizeof (ApplyColumnData));
      a->apply_choice = SourceQualValNodeDataCopy (o_list[column_index]->apply_choice);
      if (part == NULL)
      {
        a->apply_value = StringSave ("");
      }
      else
      {   
        a->apply_value = StringSave (part->data.ptrvalue);
      }

      ValNodeAddPointer (&(t->columns_to_apply), 0, a);
    }

    if (column_index < tlp->num_parts && part != NULL) 
    {
      part = part->next;
    }
  }
  return t;
}


static ValNodePtr OrgModLoadFormToApplyList (OrgModLoadFormPtr f)
{
  ValNodePtr line;
  ValNodePtr apply_list = NULL;
  Int4       column_index;
  OrgModTableColumnPtr PNTR o_list;

  if (f == NULL) return NULL;

  o_list = (OrgModTableColumnPtr PNTR) MemNew (f->num_columns * sizeof (OrgModTableColumnPtr));
  for (column_index = 0; column_index < f->num_columns; column_index++)
  {
    o_list[column_index] = (OrgModTableColumnPtr) DialogToPointer (f->columns[column_index]);
  }

  for (line = f->line_list; line != NULL; line = line->next)
  {
    ValNodeAddPointer (&apply_list, 0, TableToBioSourceFromTableLine (o_list, f->num_columns, line->data.ptrvalue));
  }
  
  for (column_index = 0; column_index < f->num_columns; column_index++)
  {
    o_list[column_index] = OrgModTableColumnFree (o_list[column_index]);
  }
  o_list = MemFree (o_list);
  return apply_list;
}


static void AddBioSourceToValNodeListIfNotAlreadyPresent (ValNodePtr PNTR list, BioSourcePtr biop)
{
  ValNodePtr vnp;

  if (list == NULL || biop == NULL) return;

  vnp = *list;
  if (vnp == NULL) {
    ValNodeAddPointer (list, 0, biop);
  } else {
    while (vnp->next != NULL) {
      if (biop == vnp->data.ptrvalue) {
        return;
      }
      vnp = vnp->next;
    }
    vnp->next = ValNodeNew(NULL);
    vnp->next->choice = 0;
    vnp->next->data.ptrvalue = biop;
  }
}


static void FindOrganismsForTableToBioSourceByTaxName (BioSourcePtr biop, Pointer userdata)
{
  ValNodePtr apply_list, vnp, vnp_id;
  TableToBioSourcePtr    t;  

  if (biop == NULL || biop->org == NULL || userdata == NULL) return;

  apply_list = (ValNodePtr) userdata;
  for (vnp = apply_list; vnp != NULL; vnp = vnp->next)
  {
    t = (TableToBioSourcePtr) vnp->data.ptrvalue;
    for (vnp_id = t->identifiers; vnp_id != NULL; vnp_id = vnp_id->next)
    {
      if (vnp_id->choice == eMatchTaxName && StringCmp (biop->org->taxname, vnp_id->data.ptrvalue) == 0) {
        AddBioSourceToValNodeListIfNotAlreadyPresent (&(t->biop_list), biop);
      }
    }
  }
}


/* look for matches based on accession numbers, local IDs, or bankit IDs */
static void FindOrganismsForTableToBioSourceByBioseq (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr             apply_list, vnp, vnp_id;
  TableToBioSourcePtr    t;  
  SeqDescrPtr            sdp;
  GBBlockPtr             gbp;
  BioSourcePtr           biop;
  Boolean                use_local_id, look_for_tmsmart, look_for_general;
  Char                   match_str [30];
  SeqIdPtr               sip;
  DbtagPtr               dp;
  ObjectIdPtr            oip;

  biop = GetBiopForBsp (bsp);
  if (biop == NULL || biop->org == NULL || userdata == NULL) return;

  apply_list = (ValNodePtr) userdata;

  /* Find GenBankBlockPtr for extra accessions */
  gbp = NULL;
  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);
  if (sdp != NULL)
  {
    gbp = sdp->data.ptrvalue;
  }

  /* Find Bankit ID */
  match_str [0] = 0;
  for(sip = bsp->id; sip != NULL; sip = sip->next) {
    if(sip->choice == SEQID_GENERAL) {
      dp = (DbtagPtr) sip->data.ptrvalue;
      if(StringICmp(dp->db, "BankIt") == 0) {
        oip = dp->tag;
        sprintf (match_str, "%d", oip->id);
        break;
      }
    }
  }

  for (vnp = apply_list; vnp != NULL; vnp = vnp->next)
  {
    t = (TableToBioSourcePtr) vnp->data.ptrvalue;
    for (vnp_id = t->identifiers; vnp_id != NULL; vnp_id = vnp_id->next)
    {
      if (vnp_id->choice == eMatchAccession 
          || vnp_id->choice == eMatchLocalID 
          || vnp_id->choice == eMatchTMSMART
          || vnp_id->choice == eMatchGeneral)
      {
        if (vnp_id->choice == eMatchAccession 
            || vnp_id->choice == eMatchTMSMART
            || vnp_id->choice == eMatchGeneral) {
          use_local_id = FALSE;
        } else {
          use_local_id = TRUE;
        }
        if (vnp_id->choice == eMatchTMSMART) {
          look_for_tmsmart = TRUE;
        } else {
          look_for_tmsmart = FALSE;
        }

        if (vnp_id->choice == eMatchGeneral)
        {
          look_for_general = TRUE;
        }
        else
        {
          look_for_general = FALSE;
        }

        if (IDListHasValue ( vnp_id->data.ptrvalue,
                            bsp->id, use_local_id, look_for_tmsmart, look_for_general)
              || (! use_local_id && 
                  HasExtraAccession ( vnp_id->data.ptrvalue, gbp)))
        {
          AddBioSourceToValNodeListIfNotAlreadyPresent (&(t->biop_list), biop);
        }
      }        
      else if (vnp_id->choice == eMatchBankIt) 
      {
        if (StringCmp (FindBankitNumberInTableString(vnp_id->data.ptrvalue), 
                          match_str) == 0) {
          AddBioSourceToValNodeListIfNotAlreadyPresent (&(t->biop_list), biop);
        }
      }
    }
  }
}

  
static void FindOrganismsForTableToBioSourceList (ValNodePtr apply_list, SeqEntryPtr sep)
{
  VisitBioSourcesInSep (sep, apply_list, FindOrganismsForTableToBioSourceByTaxName);
  VisitBioseqsInSep (sep, apply_list, FindOrganismsForTableToBioSourceByBioseq);

}


static Boolean ListOrganismsWithMultipleRows (ValNodePtr apply_list)
{
  LogInfoPtr          lip;
  ValNodePtr          vnp, vnp_compare, vnp_biop, vnp_biop2;
  BioSourcePtr        biop;
  TableToBioSourcePtr t, t2;
  Int4                row1_num, row2_num;
  Boolean             rval;

  lip = OpenLog ("Organisms Affected by More Than One Row");
  for (vnp = apply_list, row1_num = 0; vnp != NULL; vnp = vnp->next, row1_num++) 
  {
    t = (TableToBioSourcePtr) vnp->data.ptrvalue;
    if (t == NULL || t->biop_list == NULL) continue;
    for (vnp_biop = t->biop_list; vnp_biop != NULL; vnp_biop = vnp_biop->next)
    {
      biop = vnp_biop->data.ptrvalue;
      if (biop == NULL) continue;
      for (vnp_compare = vnp->next, row2_num = row1_num + 1; vnp_compare != NULL; vnp_compare = vnp_compare->next, row2_num++)
      {
        t2 = (TableToBioSourcePtr) vnp_compare->data.ptrvalue;
        for (vnp_biop2 = t2->biop_list; vnp_biop2 != NULL; vnp_biop2 = vnp_biop2->next)
        {
          if (vnp_biop2->data.ptrvalue == biop)
          {
            fprintf (lip->fp, "Two rows (%d and %d) match the same organism", row1_num, row2_num);
            if (t->identifiers != NULL && t->identifiers->data.ptrvalue != NULL
                && t2->identifiers != NULL && t2->identifiers->data.ptrvalue != NULL)
            {
              if (StringCmp (t->identifiers->data.ptrvalue, t2->identifiers->data.ptrvalue) == 0)
              {
                fprintf (lip->fp, " (both are %s).\n", (CharPtr) t->identifiers->data.ptrvalue);
              }
              else
              {
                fprintf (lip->fp, " (%s and %s).\n", (CharPtr) t->identifiers->data.ptrvalue, (CharPtr) t2->identifiers->data.ptrvalue);
              }
            }
            else if (t->identifiers != NULL && t->identifiers->data.ptrvalue != NULL)
            {
              fprintf (lip->fp, " (%s).\n", (CharPtr) t->identifiers->data.ptrvalue);
            }
            else if (t2->identifiers != NULL && t2->identifiers->data.ptrvalue != NULL)
            {
              fprintf (lip->fp, " (%s).\n", (CharPtr) t2->identifiers->data.ptrvalue);
            }
            else
            {
              fprintf (lip->fp, ".\n");
            }
                
            lip->data_in_log = TRUE;
          }
        }
      }
    }
  }
  rval = lip->data_in_log;
  CloseLog (lip);
  FreeLog (lip);
  return rval;
}


static Boolean ListRowsWithMultipleOrganisms (ValNodePtr apply_list)
{
  LogInfoPtr          lip;
  ValNodePtr          vnp, vnp_id;
  TableToBioSourcePtr t;
  Boolean             rval;

  lip = OpenLog ("Rows with Multiple Organisms");
  for (vnp = apply_list; vnp != NULL; vnp = vnp->next) 
  {
    t = (TableToBioSourcePtr) vnp->data.ptrvalue;
    if (t == NULL || t->biop_list == NULL || t->biop_list->next == NULL || t->identifiers == NULL) continue;
    fprintf (lip->fp, "Row with the following match columns applies to %d organisms\n", ValNodeLen (t->biop_list));
    for (vnp_id = t->identifiers; vnp_id != NULL; vnp_id = vnp_id->next) 
    {
      fprintf (lip->fp, "\t%s\n", (CharPtr) vnp_id->data.ptrvalue);
    }
    lip->data_in_log = TRUE;
  }

  rval = lip->data_in_log;
  CloseLog (lip);
  FreeLog (lip);
  return rval;
}


static void ListRowsWithoutOrganisms (ValNodePtr apply_list)
{
  LogInfoPtr          lip;
  ValNodePtr          vnp, vnp_id;
  TableToBioSourcePtr t;

  lip = OpenLog ("Rows without Organisms");
  for (vnp = apply_list; vnp != NULL; vnp = vnp->next) 
  {
    t = (TableToBioSourcePtr) vnp->data.ptrvalue;
    if (t == NULL || t->biop_list != NULL || t->identifiers == NULL) continue;
    for (vnp_id = t->identifiers; vnp_id != NULL; vnp_id = vnp_id->next) 
    {
      fprintf (lip->fp, "No match found for %s\n", (CharPtr) vnp_id->data.ptrvalue);
      lip->data_in_log = TRUE;
    }
  }

  CloseLog (lip);
  FreeLog (lip);
}


/* This code is used to detect existing modifiers before a change is made */
static void FindOneQualOnOrg 
( BioSourcePtr biop,
  SourceQualDescPtr sqdp,
  GetSamplePtr gsp)
{
  OrgModPtr     mod;
  SubSourcePtr  ssp;

  if (biop == NULL || gsp == NULL || sqdp == NULL) return;

  if (sqdp->isOrgMod)
  {
    if (biop->org == NULL
        || biop->org->orgname == NULL)
    {
      return;
    }
    
    mod = biop->org->orgname->mod;
    while (mod != NULL
           && mod->subtype != sqdp->subtype)
    {
      mod = mod->next;
    }
    if (mod != NULL)
    {
      gsp->num_found ++;
      if (gsp->sample_text == NULL)
      {
        gsp->sample_text = StringSave (mod->subname);
      }
      else if (StringCmp (gsp->sample_text, mod->subname) != 0)
      {
        gsp->all_same = FALSE;
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
    if (ssp != NULL)
    {
      gsp->num_found ++;
      if (gsp->sample_text == NULL)
      {
        gsp->sample_text = StringSave (ssp->name);
      }
      else if (StringCmp (gsp->sample_text, ssp->name) != 0)
      {
        gsp->all_same = FALSE;
      }
    }
  }
}


static GetSamplePtr DetectExistingModifiersForTableToBioSourceList (ValNodePtr apply_list, Boolean replace_with_blank)
{
  ValNodePtr vnp, vnp_col, vnp_biop;
  TableToBioSourcePtr    t; 
  BioSourcePtr           biop; 
  ApplyColumnPtr         a;
  GetSamplePtr           gsp;

  gsp = GetSampleNew ();
  for (vnp = apply_list; vnp != NULL; vnp = vnp->next) 
  {
    t = (TableToBioSourcePtr) vnp->data.ptrvalue;
    for (vnp_biop = t->biop_list; vnp_biop != NULL; vnp_biop = vnp_biop->next) 
    {
      biop = vnp_biop->data.ptrvalue;
      for (vnp_col = t->columns_to_apply; vnp_col != NULL; vnp_col = vnp_col->next) 
      {
        a = (ApplyColumnPtr) vnp_col->data.ptrvalue;

        if (replace_with_blank || ! StringHasNoText (a->apply_value))
        {
          if (a->apply_choice->choice == 0 && a->apply_choice->data.ptrvalue != NULL)
          {
            FindOneQualOnOrg (biop, a->apply_choice->data.ptrvalue, gsp);
          }
          else if (a->apply_choice->choice > 0
                 && (StringICmp (a->apply_choice->data.ptrvalue, "Tax Name") == 0
                     || StringICmp (a->apply_choice->data.ptrvalue, "Organism") == 0))

          {  
            if (biop->org != NULL && !StringHasNoText (biop->org->taxname))
            {
              gsp->num_found ++;
              if (gsp->sample_text == NULL)
              {
                gsp->sample_text = StringSave (biop->org->taxname);
              }
              else if (StringCmp (gsp->sample_text, biop->org->taxname) != 0)
              {
                gsp->all_same = FALSE;
              }
            }
          }
          else if (a->apply_choice->choice > 0 && StringICmp (a->apply_choice->data.ptrvalue, "Location") == 0)
          {
            if (biop->genome > 0)
            {
              gsp->num_found ++;
              if (gsp->sample_text == NULL)
              {
                gsp->sample_text = StringSave (genome_names [biop->genome]);
              }
              else if (StringCmp (gsp->sample_text, genome_names [biop->genome]) != 0)
              {
                gsp->all_same = FALSE;
              }
            }
          }
        }        
      }  
    }
  }
  return gsp;
}


static void 
ApplyOneTableToBioSource 
(TableToBioSourcePtr t,
 Boolean             replace_with_blank,
 ExistingTextPtr     etp,
 TaxnameOptionsPtr   taxname_options)
{
  ValNodePtr vnp_col, vnp_biop;
  BioSourcePtr biop;
  ApplyColumnPtr a;

  if (t == NULL || t->columns_to_apply == NULL || t->biop_list == NULL) return;

  for (vnp_biop = t->biop_list; vnp_biop != NULL; vnp_biop = vnp_biop->next) 
  {
    biop = (BioSourcePtr) vnp_biop->data.ptrvalue;
    if (biop == NULL) continue;
    for (vnp_col = t->columns_to_apply; vnp_col != NULL; vnp_col = vnp_col->next)
    {
      a = (ApplyColumnPtr) vnp_col->data.ptrvalue;
      if (a == NULL 
          || a->apply_choice == NULL
          || (!replace_with_blank && StringHasNoText (a->apply_value))) continue;

      if (a->apply_choice->choice == 0 && a->apply_choice->data.ptrvalue != NULL)
      {
        AddOneQualToOrg (biop, a->apply_value, a->apply_choice->data.ptrvalue,
                          etp);
      }
      else if (a->apply_choice->choice > 0 
                && (StringICmp (a->apply_choice->data.ptrvalue, "Tax Name") == 0
                    || StringICmp (a->apply_choice->data.ptrvalue, "Organism") == 0))
      {  
        ApplyTaxNameToOrg (biop, a->apply_value, etp);
        ApplyTaxnameOptionsToBioSource (biop, taxname_options);
      }
      else if (a->apply_choice->choice > 0 && StringICmp (a->apply_choice->data.ptrvalue, "Location") == 0)
      {
        ApplyLocationToOrg (biop, a->apply_value);
      }
    }
  }
}


static void 
ApplyTableToBioSourceList 
(ValNodePtr        apply_list,
 Boolean           replace_with_blank,
 ExistingTextPtr   etp,
 TaxnameOptionsPtr taxname_options)
{
  ValNodePtr vnp;

  for (vnp = apply_list; vnp != NULL; vnp = vnp->next)
  {
    ApplyOneTableToBioSource (vnp->data.ptrvalue, replace_with_blank, etp, taxname_options);
  }
}


static void DoAcceptFormMods (ButtoN b)
{
  OrgModLoadFormPtr form_data;
  SeqEntryPtr   sep;
  GetSamplePtr gsp;
  ValNodePtr apply_list;
  Boolean    replace_with_blank;
  ExistingTextPtr etp;
  TaxnameOptionsPtr taxname_options;

  form_data = GetObjectExtra (b);
  if (form_data == NULL) return;

  sep = GetTopSeqEntryForEntityID (form_data->entityID);
  if (sep == NULL) return;
  
  replace_with_blank = GetStatus (form_data->replace_with_blank_btn);
  apply_list = OrgModLoadFormToApplyList (form_data);
  FindOrganismsForTableToBioSourceList (apply_list, sep);
  if (ListOrganismsWithMultipleRows (apply_list)) 
  {
    Message (MSG_ERROR, "Multiple rows in the table apply to the same organism!  Table cannot be imported.");
    apply_list = TableToBioSourceListFree (apply_list);
    return;
  }
  ListRowsWithoutOrganisms (apply_list);
  if (ListRowsWithMultipleOrganisms (apply_list)) {
    if (ANS_CANCEL == Message (MSG_OKC, "Some rows in the table apply to more than one organism - continue?"))
    {
      apply_list = TableToBioSourceListFree (apply_list);
      return;
    }
  }
  gsp = DetectExistingModifiersForTableToBioSourceList (apply_list, replace_with_blank);
  etp = GetExistingTextHandlerInfo (gsp == NULL ? 0 : gsp->num_found, FALSE);
  gsp = GetSampleFree (gsp);
  if (etp != NULL && etp->existing_text_choice == eExistingTextChoiceCancel)
  {
    etp = MemFree (etp);
    apply_list = TableToBioSourceListFree (apply_list);
    return;
  }
  taxname_options = DialogToPointer (form_data->taxname_options);
  ApplyTableToBioSourceList (apply_list, replace_with_blank, etp, taxname_options);
  taxname_options = MemFree (taxname_options);
  etp = MemFree (etp);
  apply_list = TableToBioSourceListFree (apply_list);

  Update ();
  ObjMgrSetDirtyFlag (form_data->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, form_data->entityID, 0, 0);
  if (!GetStatus (form_data->leave_dlg_up_btn))
  {
    Remove (form_data->form);
  }
}


static ValNodePtr ReadTableData (void)
{
  Char          path [PATH_MAX];
  Int4          max_columns;
  ValNodePtr    header_line;
  ValNodePtr    line_list;
  TableLinePtr  tlp;
  ValNodePtr    vnp;
  ReadBufferData rbd;
  CharPtr        line;
  ValNodePtr     special_list = NULL;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return NULL;
  
  rbd.fp = FileOpen (path, "r");
  if (rbd.fp == NULL) return NULL;
  rbd.current_data = NULL;

  line_list = NULL;
  max_columns = 0;
  header_line = NULL;
  line = AbstractReadFunction (&rbd);
  while (line != NULL) 
  {
    tlp = MakeTableData (line, &special_list);
    if (tlp != NULL)
    {
      vnp = ValNodeNew (line_list);
      if (vnp == NULL) return NULL;
      if (line_list == NULL) line_list = vnp;
      vnp->data.ptrvalue = tlp;
      if (tlp->num_parts > max_columns)
      {
        max_columns = tlp->num_parts;
        header_line = vnp;
      }
    }
    line = MemFree (line);
    line = AbstractReadFunction (&rbd);
  }
  FileClose (rbd.fp);
  if (special_list != NULL && !FixSpecialCharactersForStringsInList (special_list, "The table contains special characters\nand cannot be used until they are replaced.", FALSE))
  {
    line_list = ValNodeFreeData (line_list);
    header_line = NULL;
  }
  /* throw out all lines before header line */
  else if (header_line != line_list)
  {
    vnp = line_list;
    while (vnp != NULL && vnp->next != header_line)
    {
      vnp = vnp->next;
    }
    vnp->next = NULL;
    ValNodeFreeData (line_list);
    line_list = NULL;
  }
  special_list = FreeContextList (special_list);
  return header_line;
}


static Int4Ptr GetColumnBlankCounts (ValNodePtr header_line, Int4 num_columns)
{
  ValNodePtr   line_vnp, part_vnp;
  TableLinePtr tlp;
  Int4         i;
  Int4Ptr      blank_list;
  
  if (header_line == NULL || num_columns < 1) return NULL;

  blank_list = (Int4Ptr) MemNew (sizeof (Int4) * num_columns);

  for (line_vnp = header_line; line_vnp != NULL; line_vnp = line_vnp->next)
  {
    if (line_vnp->data.ptrvalue == NULL) continue;
    tlp = (TableLinePtr) line_vnp->data.ptrvalue;
    part_vnp = tlp->parts;
    for (i = 0; i < num_columns; i++) 
    {
      if (part_vnp == NULL || StringHasNoText (part_vnp->data.ptrvalue)) 
      {
        blank_list[i]++;
      }
      if (part_vnp != NULL) 
      {
        part_vnp = part_vnp->next;
      }
    }
  }
  return blank_list;
}


NLM_EXTERN void LoadOrganismModifierTableBaseForm (BaseFormPtr bfp, Boolean IsTaxConsult)
{
  ValNodePtr    header_line;
  ValNodePtr    vnp;
  TableLinePtr  tlp;
  WindoW        w;
  GrouP         h, g, k, c;
  OrgModLoadFormPtr form_data;
  Int4          index;
  Int4          max_columns;
  Int4Ptr       blank_list = NULL;
  OrgModTableColumnData o;

  if (bfp == NULL) return;

  header_line = ReadTableData ();
  if (header_line == NULL || header_line->data.ptrvalue == NULL) return;
  tlp = header_line->data.ptrvalue;
  max_columns = 0;
  for (vnp = tlp->parts; vnp != NULL; vnp = vnp->next)
  {
    max_columns ++;
  }
  
  form_data = MemNew (sizeof (OrgModLoadFormData));
  if (form_data == NULL) return;
  form_data->entityID = bfp->input_entityID;
  form_data->line_list = header_line;

  form_data->num_columns = max_columns;
  form_data->columns = (DialoG PNTR) MemNew (max_columns * sizeof (DialoG));

  /* now create a dialog to display values */
  w = FixedWindow (-50, -33, -10, -10, "Table Conversion", StdCloseWindowProc);
  SetObjectExtra (w, form_data, CleanupOrgModLoadForm);
  form_data->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 3, 0, NULL);

  /* pre-analyze table for blanks, so we can display the information in the column dialogs */
  blank_list = GetColumnBlankCounts (header_line, form_data->num_columns);

  tlp = header_line->data.ptrvalue;
  for (index = 0, vnp = tlp->parts; index < max_columns; index++)
  {
    form_data->columns[index] = OrgModTableColumnDialog (g, vnp == NULL ? NULL : vnp->data.ptrvalue,
                                                         blank_list == NULL ? 0 : blank_list[index],
                                                         SetFormModsEnable, form_data);
    if (vnp != NULL) vnp = vnp->next;
  }
  blank_list = MemFree (blank_list);
  if (max_columns > 1 && IsTaxConsult)
  {
    o.match_choice = eMatchAccession;
    o.apply_choice = NULL;
    PointerToDialog (form_data->columns[0], &o);
    o.match_choice = 0;
    ValNodeAddPointer (&o.apply_choice, 1, "Tax Name");
    PointerToDialog (form_data->columns[1], &o);
    o.apply_choice = ValNodeFree (o.apply_choice);
  }

  k = HiddenGroup (h, 1, 0, NULL);
  form_data->replace_with_blank_btn = CheckBox (k, "Erase current value when blank found in table", NULL);

  form_data->taxname_options = TaxnameOptionsDialog (h, NULL, NULL);
  Disable (form_data->taxname_options);

  c = HiddenGroup (h, 4, 0, NULL);
  form_data->accept_button = DefaultButton (c, "Accept", DoAcceptFormMods);
  SetObjectExtra (form_data->accept_button, form_data, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  form_data->leave_dlg_up_btn = CheckBox (c, "Leave Dialog Up", NULL);

  SetFormModsEnable(form_data);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) form_data->taxname_options, (HANDLE) c, NULL);

  RealizeWindow (w);
  Show (w);
  Update ();
}

static void LoadOrganismModifierTableEx (IteM i, Boolean IsTaxConsult)
{
  BaseFormPtr   bfp;
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  LoadOrganismModifierTableBaseForm (bfp, IsTaxConsult);
}

extern void LoadOrganismModifierTable (IteM i)
{
  LoadOrganismModifierTableEx (i, FALSE);
}

extern void LoadTaxConsult (IteM i)
{
  LoadOrganismModifierTableEx (i, TRUE);
}


static void LoadOrganismModifierTableButtonEx (ButtoN b, Boolean IsTaxConsult)
{
  BaseFormPtr   bfp;

  bfp = (BaseFormPtr) GetObjectExtra(b);
  if (b == NULL) {
    return;
  }
  LoadOrganismModifierTableBaseForm (bfp, IsTaxConsult);
}


extern void LoadTaxConsultBtn (ButtoN b)
{
  LoadOrganismModifierTableButtonEx (b, TRUE);
}



extern Boolean ReplaceImportModifierName (CharPtr PNTR orig_name, Int4 col_num, ValNodePtr preferred_list)
{
  ModalAcceptCancelData acd;
  WindoW                w;
  Char                  title[128];
  GrouP                 g, k, c;
  ButtoN                b;
  DialoG                modifier_name_choice;
  ValNodePtr            mod_name_choice_list = NULL, default_vnp, choice_vnp;
  
  if (orig_name == NULL)
  {
    return FALSE;
  }
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  sprintf (title, "Choose import modifier for column %d", col_num);
  w = MovableModalWindow(-20, -13, -10, -10, title, NULL);
  g = HiddenGroup (w, -1, 0, NULL);
  
  k = NULL;
  if (!StringHasNoText (*orig_name))
  {
    k = HiddenGroup (g, 2, 0, NULL);
    StaticPrompt (k, "Unrecognized column name:", 0, popupMenuHeight, programFont, 'r');
    StaticPrompt (k, *orig_name, 0, popupMenuHeight, programFont, 'l');
  }
  
  ValNodeAddPointer (&mod_name_choice_list, 1, StringSave ("Ignore Column"));
  ValNodeLink (&mod_name_choice_list, ValNodeDupList(preferred_list, SourceQualValNodeCopy));
  ValNodeLink(&mod_name_choice_list, GetSourceQualDescList (TRUE, TRUE, FALSE, FALSE));
  ValNodeAddPointer (&mod_name_choice_list, 2, StringSave ("Organism"));
  ValNodeAddPointer (&mod_name_choice_list, 3, StringSave ("Location"));
  
  modifier_name_choice = ValNodeSelectionDialog (g, mod_name_choice_list, 6, SourceQualValNodeName,
                                ValNodeSimpleDataFree, SourceQualValNodeDataCopy,
                                SourceQualValNodeMatch, "feature list", 
                                NULL, NULL, FALSE);
  default_vnp = ValNodeNew (NULL);
  default_vnp->choice = 1;
  default_vnp->data.ptrvalue = StringSave ("Ignore Column");
  default_vnp->next = NULL;
  PointerToDialog (modifier_name_choice, default_vnp);
  default_vnp = ValNodeFreeData (default_vnp);

  c = HiddenGroup (g, 2, 0, NULL);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  if (k == NULL)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) modifier_name_choice, (HANDLE) c, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) modifier_name_choice, (HANDLE) c, NULL);
  }
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
  choice_vnp = (ValNodePtr) DialogToPointer (modifier_name_choice);
  
  Remove (w);
  if (acd.cancelled)
  {
    return FALSE;
  }
  else
  {
    *orig_name = MemFree (*orig_name);
    if (choice_vnp != NULL && choice_vnp->choice != 1)
    {
      *orig_name = SourceQualValNodeName (choice_vnp);
    }
    return TRUE;
  }
}


static ValNodePtr FreeTableData (ValNodePtr lines)
{
  TableLinePtr tlp;
  
  if (lines == NULL)
  {
    return NULL;
  }
  lines->next = FreeTableData (lines->next);
  tlp = (lines->data.ptrvalue);
  if (tlp != NULL)
  {
    tlp->parts = ValNodeFreeData (tlp->parts);
    tlp = MemFree (tlp);
    lines->data.ptrvalue = NULL;
  }
  return ValNodeFree (lines);
}


typedef struct exportorgtable
{
  ModifierItemLocalPtr modList;
  BioseqPtr            bsp;
  Boolean              list_tax_name;
  Boolean              list_accession;
  Boolean              list_local;
  Boolean              list_general;
  Boolean              launch_excel;

  ValNodePtr           organism_id_profile;
  FILE *fp;  
} ExportOrgTableData, PNTR ExportOrgTablePtr;

#define ORGANISM_ID_PROFILE_HAS_ACCESSION 1
#define ORGANISM_ID_PROFILE_HAS_LOCAL     2
#define ORGANISM_ID_PROFILE_HAS_GENERAL   4
#define ORGANISM_ID_PROFILE_HAS_TAX_NAME  8

static void ExportOneOrganism (BioSourcePtr biop, Pointer userdata)
{
  ExportOrgTablePtr eotp;
  Char              acc_str [256];
  SeqIdPtr          sip;
  SeqIdPtr          acc_sip = NULL, local_sip = NULL, gen_sip = NULL;
  Int4              i;
  OrgModPtr         mod;
  SubSourcePtr      ssp;
  Boolean           need_tab = FALSE;
  
  if (biop == NULL || userdata == NULL) return;
  eotp = (ExportOrgTablePtr) userdata;
  
  if (eotp->bsp == NULL || eotp->modList == NULL || eotp->fp == NULL) return;
  
  for (sip = eotp->bsp->id;
       sip != NULL && (acc_sip == NULL || local_sip == NULL);
       sip = sip->next)
  {
    if (acc_sip == NULL && sip->choice == SEQID_GENBANK)
    {
      acc_sip = sip;
    }
    else if (local_sip == NULL && sip->choice == SEQID_LOCAL)
    {
      local_sip = sip;
    }
    else if (gen_sip == NULL && sip->choice == SEQID_GENERAL)
    {
      gen_sip = sip;
    }
  }

  if (eotp->list_accession)
  {
    /* get accession number and print to column */
    if (acc_sip == NULL)
    {
      fprintf (eotp->fp, " ");
    }
    else
    {
      sip = acc_sip->next;
      acc_sip->next = NULL;
      SeqIdWrite (acc_sip, acc_str, PRINTID_TEXTID_ACC_VER, sizeof (acc_str));
      acc_sip->next = sip;
      fprintf (eotp->fp, "%s", acc_str);
    }
    need_tab = TRUE;
  }
  
  if (eotp->list_local)
  {
    if (need_tab) {
      fprintf (eotp->fp, "\t");
    }
    /* get local ID and print to column */
    if (local_sip == NULL)
    {
      fprintf (eotp->fp, " ");
    }
    else
    {
      sip = local_sip->next;
      local_sip->next = NULL;
      SeqIdWrite (local_sip, acc_str, PRINTID_TEXTID_ACCESSION, sizeof (acc_str));
      local_sip->next = sip;
      fprintf (eotp->fp, "%s", acc_str);
    }
    need_tab = TRUE;
  }
  
  if (eotp->list_general)
  {
    if (need_tab) {
      fprintf (eotp->fp, "\t");
    }
    /* get general ID and print to column */
    if (gen_sip == NULL)
    {
      fprintf (eotp->fp, " ");
    }
    else
    {
      sip = gen_sip->next;
      gen_sip->next = NULL;
      SeqIdWrite (gen_sip, acc_str, PRINTID_FASTA_GENERAL, sizeof (acc_str));
      gen_sip->next = sip;
      fprintf (eotp->fp, "%s", acc_str);
    }
    need_tab = TRUE;
  }
  
  if (eotp->list_tax_name)
  {
    if (need_tab) {
      fprintf (eotp->fp, "\t");
    }
    /* get tax name and print to column */
    if (biop->org != NULL && ! StringHasNoText (biop->org->taxname))
    {
      fprintf (eotp->fp, "%s", biop->org->taxname);
    }
    else
    {
      fprintf (eotp->fp, " ");
    }
    need_tab = TRUE;
  }
  
  /* print modifiers for each available column */
  for (i= 0; i < NumDefLineModifiers (); i++)
  {
    if (eotp->modList[i].any_present)
    {
      mod = NULL;
      ssp = NULL;
      if (DefLineModifiers[i].isOrgMod)
      {
        if ( biop->org != NULL && biop->org->orgname != NULL) 
        {
          mod = biop->org->orgname->mod;
          while (mod != NULL
                 && mod->subtype != DefLineModifiers[i].subtype)
          {
            mod = mod->next;
          }
        }
      }
      else
      {
        ssp = biop->subtype;
        while (ssp != NULL && ssp->subtype != DefLineModifiers[i].subtype)
        {
          ssp = ssp->next;
        }
      }
      if (need_tab) {
        fprintf (eotp->fp, "\t");
      }
	  if (IsNonTextModifier (DefLineModifiers[i].name))
	  {
	    if (mod == NULL && ssp == NULL)
	    {
	      fprintf (eotp->fp, "FALSE");
		}
		else
		{
		  fprintf (eotp->fp, "TRUE");
		}
	  }
      else if (mod != NULL && !StringHasNoText (mod->subname))
      {
        fprintf (eotp->fp, "%s", mod->subname);
      }
      else if (ssp != NULL && !StringHasNoText (ssp->name))
      {
        fprintf (eotp->fp, "%s", ssp->name);
      }
      else
      {
        fprintf (eotp->fp, " ");
      }
      need_tab = TRUE;
    }
  }
  fprintf (eotp->fp, "\n");
}

static void ExportOrganisms (SeqEntryPtr sep, ExportOrgTablePtr eotp)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  nsep;
  BioseqPtr    bsp;
  
  if (sep == NULL || eotp == NULL || sep->data.ptrvalue == NULL) return;
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && !ISA_aa (bsp->mol)) {
      eotp->bsp = bsp;
      VisitBioSourcesOnBsp (eotp->bsp, eotp, ExportOneOrganism);
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_nuc_prot)
    {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL)
      {
        eotp->bsp = (BioseqPtr) nsep->data.ptrvalue;
        VisitBioSourcesOnSep (sep, eotp, ExportOneOrganism);
      }
    }
    else
    {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        ExportOrganisms (sep, eotp);
      }
    }
  }
}

static void GetOneOrganismIdProfile (BioSourcePtr biop, Pointer userdata)
{
  ExportOrgTablePtr eotp;
  SeqIdPtr          sip;
  Int4              profile_value = 0;

  if (biop == NULL || userdata == NULL) return;
  eotp = (ExportOrgTablePtr) userdata;
  
  if (eotp->bsp == NULL) return;

  for (sip = eotp->bsp->id;
       sip != NULL;
       sip = sip->next)
  {
    if (sip->choice == SEQID_GENBANK)
    {
	  profile_value |= ORGANISM_ID_PROFILE_HAS_ACCESSION;
    }
    else if (sip->choice == SEQID_LOCAL)
    {
      profile_value |= ORGANISM_ID_PROFILE_HAS_LOCAL;
    }
    else if (sip->choice == SEQID_GENERAL)
    {
      profile_value |= ORGANISM_ID_PROFILE_HAS_GENERAL;
    }
  }
  
  if (biop->org != NULL && ! StringHasNoText (biop->org->taxname))
  {
    profile_value |= ORGANISM_ID_PROFILE_HAS_TAX_NAME;
  }

  ValNodeAddInt (&(eotp->organism_id_profile), 0, profile_value);
}

static void GetOrganismIDProfile (SeqEntryPtr sep, ExportOrgTablePtr eotp)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  nsep;
  
  if (sep == NULL || eotp == NULL || sep->data.ptrvalue == NULL) return;
  
  if (IS_Bioseq (sep))
  {
    eotp->bsp = (BioseqPtr) sep->data.ptrvalue;
    VisitBioSourcesOnBsp (eotp->bsp, eotp, GetOneOrganismIdProfile);  
  }
  else if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->_class == BioseqseqSet_class_nuc_prot)
    {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL)
      {
        eotp->bsp = (BioseqPtr) nsep->data.ptrvalue;
        VisitBioSourcesOnSep (sep, eotp, GetOneOrganismIdProfile);
      }
    }
    else
    {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        GetOrganismIDProfile (sep, eotp);
      }
    }
  }
}

static void GetOrganismIDProfileSummary (ValNodePtr organism_id_profile, Int4Ptr available, Int4Ptr recommended)
{
  ValNodePtr vnp;
  Int4       has_vals = 0, recommended_vals = 0, missing_vals = 0;

  if (available != NULL)
  {
    *available = 0;
  }
  if (recommended != NULL)
  {
    *recommended = 0;
  }

  for (vnp = organism_id_profile;
	   vnp != NULL;
	   vnp = vnp->next)
  {
    has_vals |= vnp->data.intvalue;
	missing_vals |= !(vnp->data.intvalue);
	if (vnp->data.intvalue & ORGANISM_ID_PROFILE_HAS_ACCESSION)
	{
	  recommended_vals |= ORGANISM_ID_PROFILE_HAS_ACCESSION;
	}
    else if (vnp->data.intvalue & ORGANISM_ID_PROFILE_HAS_LOCAL)
	{
	  recommended_vals |= ORGANISM_ID_PROFILE_HAS_LOCAL;
	}
	else if (vnp->data.intvalue & ORGANISM_ID_PROFILE_HAS_GENERAL)
	{
	  recommended_vals |= ORGANISM_ID_PROFILE_HAS_GENERAL;
	}
	else if (vnp->data.intvalue & ORGANISM_ID_PROFILE_HAS_TAX_NAME)
	{
	  recommended_vals |= ORGANISM_ID_PROFILE_HAS_TAX_NAME;
	}
  }

  if (recommended_vals & ORGANISM_ID_PROFILE_HAS_TAX_NAME)
  {
    if (!(missing_vals & ORGANISM_ID_PROFILE_HAS_TAX_NAME))
	{
	  recommended_vals = ORGANISM_ID_PROFILE_HAS_TAX_NAME;
	}
  }
  else if (recommended_vals & ORGANISM_ID_PROFILE_HAS_GENERAL)
  {
    if (!(missing_vals & ORGANISM_ID_PROFILE_HAS_GENERAL))
	{
	  recommended_vals = ORGANISM_ID_PROFILE_HAS_GENERAL;
	}
  }
  else if (recommended_vals & ORGANISM_ID_PROFILE_HAS_LOCAL)
  {
    if (!(missing_vals & ORGANISM_ID_PROFILE_HAS_LOCAL))
	{
	  recommended_vals = ORGANISM_ID_PROFILE_HAS_LOCAL;
	}
  }
  
  if (available != NULL)
  {
    *available = has_vals;
  }

  if (recommended != NULL)
  {
    *recommended = recommended_vals;
  }
}

typedef struct exportmodform
{
  FEATURE_FORM_BLOCK

  DialoG selected_mods;
  ButtoN list_tax_name;
  ButtoN list_accession;
  ButtoN list_local;
  ButtoN list_general;
  ButtoN list_additional_mods;
  ButtoN launch_excel;
  ButtoN accept_btn;
} ExportModFormData, PNTR ExportModFormPtr;

static void SetExportFormModsAccept (Pointer userdata)
{
  ExportModFormPtr form_data;
  ValNodePtr       selected_mods;

  form_data = (ExportModFormPtr) userdata;
  if (form_data == NULL) return;
  
  if (GetStatus (form_data->list_additional_mods))
  {
    Enable (form_data->selected_mods); 
    selected_mods = (ValNodePtr) DialogToPointer (form_data->selected_mods);
    if (selected_mods == NULL)
    {
      Disable (form_data->accept_btn);
      return;
    }
    else
    {
	    selected_mods = ValNodeFreeData (selected_mods);
    }
  }
  else
  {
    Disable (form_data->selected_mods);
  }

	if (!GetStatus (form_data->list_tax_name) && ! GetStatus (form_data->list_accession)
		&& !GetStatus (form_data->list_local) && ! GetStatus (form_data->list_general))
	{
    Disable (form_data->accept_btn);
	}
	else
	{
	  Enable (form_data->accept_btn);
	}
}

static void SetExportFormModsAcceptBtn (ButtoN b)
{
  SetExportFormModsAccept (GetObjectExtra (b));
}

static Boolean SelectModifiersForExport (ExportOrgTablePtr eotp)
{
  Int4                  idx;
  ValNodePtr            available_mods_list = NULL, vnp;
  ModalAcceptCancelData acd;
  ExportModFormPtr      form_data;
  WindoW                w;
  GrouP                 h, g, c;
  ButtoN                b;
  ValNodePtr            selected_mods = NULL;
  Boolean               found_mod;
  Int4                  has_vals = 0, recommended_vals = 0;

  if (eotp == NULL)
  {
    return FALSE;
  }

  for (idx = 0; idx < NumDefLineModifiers (); idx++)
  {
    if (eotp->modList[idx].any_present)
    {
      ValNodeAddPointer (&available_mods_list, idx, StringSave (DefLineModifiers [idx].name));
    }
  }

  GetOrganismIDProfileSummary (eotp->organism_id_profile, &has_vals, &recommended_vals);

  form_data = (ExportModFormPtr) MemNew (sizeof (ExportModFormData));

  w = FixedWindow (-50, -33, -10, -10, "Choose Modifiers for Export", StdCloseWindowProc);
  SetObjectExtra (w, form_data, StdCleanupFormProc);
  form_data->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 5, NULL);
  form_data->list_tax_name = CheckBox (g, "Tax Name", SetExportFormModsAcceptBtn);
  SetObjectExtra (form_data->list_tax_name, form_data, NULL);
  if (!(has_vals & ORGANISM_ID_PROFILE_HAS_TAX_NAME))
  {
    Disable (form_data->list_tax_name);
  }
  else if (recommended_vals & ORGANISM_ID_PROFILE_HAS_TAX_NAME)
  {
    SetStatus (form_data->list_tax_name, TRUE);
  }
  form_data->list_accession = CheckBox (g, "Accession", SetExportFormModsAcceptBtn);
  SetObjectExtra (form_data->list_accession, form_data, NULL);
  if (!(has_vals & ORGANISM_ID_PROFILE_HAS_ACCESSION))
  {
    Disable (form_data->list_accession);
  }
  else if (recommended_vals & ORGANISM_ID_PROFILE_HAS_ACCESSION)
  {
    SetStatus (form_data->list_accession, TRUE);
  }
  form_data->list_local = CheckBox (g, "Local ID", SetExportFormModsAcceptBtn);
  SetObjectExtra (form_data->list_local, form_data, NULL);
  if (!(has_vals & ORGANISM_ID_PROFILE_HAS_LOCAL))
  {
    Disable (form_data->list_local);
  }
  else if (recommended_vals & ORGANISM_ID_PROFILE_HAS_LOCAL)
  {
    SetStatus (form_data->list_local, TRUE);
  }
  form_data->list_general = CheckBox (g, "General ID", SetExportFormModsAcceptBtn);
  SetObjectExtra (form_data->list_general, form_data, NULL);
  if (!(has_vals & ORGANISM_ID_PROFILE_HAS_GENERAL))
  {
    Disable (form_data->list_general);
  }
  else if (recommended_vals & ORGANISM_ID_PROFILE_HAS_GENERAL)
  {
    SetStatus (form_data->list_general, TRUE);
  }
  
  form_data->list_additional_mods = CheckBox (g, "Additional Modifiers", SetExportFormModsAcceptBtn);
  SetObjectExtra (form_data->list_additional_mods, form_data, NULL);
  SetStatus (form_data->list_additional_mods, TRUE);

  form_data->selected_mods = ValNodeSelectionDialog (h, available_mods_list, TALL_SELECTION_LIST,
                                ValNodeStringName,
                                ValNodeSimpleDataFree, 
                                ValNodeStringCopy,
                                ValNodeStringMatch,
                                "qualifiers", 
                                SetExportFormModsAccept, form_data,
                                TRUE);

  /* initialize dialog to "All" */
  SendMessageToDialog (form_data->selected_mods, NUM_VIB_MSG + 1);

#ifdef WIN_MSWIN
  form_data->launch_excel = CheckBox (h, "Launch Excel", NULL);
#endif

  c = HiddenGroup (h, 2, 0, NULL);
  form_data->accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (form_data->accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) form_data->selected_mods, (HANDLE) g, (HANDLE) c, NULL);

  SetExportFormModsAccept (form_data);

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
  
  Remove (w);
  if (acd.cancelled)
  {
    return FALSE;
  }
  else
  {
    if (GetStatus (form_data->list_additional_mods))
    {
      selected_mods = (ValNodePtr) DialogToPointer (form_data->selected_mods);
    }
    else
    {
      selected_mods = NULL;
    }

    for (idx = 0; idx < NumDefLineModifiers (); idx++)
    {
      if (eotp->modList[idx].any_present)
      {
        found_mod = FALSE;
   	    for (vnp = selected_mods; vnp != NULL && ! found_mod; vnp = vnp->next)
        {
          if (StringCmp (vnp->data.ptrvalue, DefLineModifiers [idx].name) == 0)
          {
            found_mod = TRUE;
          }
        }
        if (!found_mod)
        {
          eotp->modList[idx].any_present = FALSE;
        }
      }
    }
    selected_mods = ValNodeFreeData (selected_mods);

    eotp->list_tax_name = GetStatus (form_data->list_tax_name);
    eotp->list_accession = GetStatus (form_data->list_accession);
    eotp->list_local = GetStatus (form_data->list_local);
    eotp->list_general = GetStatus (form_data->list_general);
#ifdef WIN_MSWIN
    eotp->launch_excel = GetStatus (form_data->launch_excel);
#endif

    return TRUE;
  }  
}

extern void ExportOrganismTable (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  ExportOrgTableData eotd;
  Char               path [PATH_MAX];
  Int4               idx;
  Char               lead_str[2];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  MemSet (&eotd, 0, sizeof (ExportOrgTableData));
  eotd.modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (eotd.modList == NULL) return;

  CountModifiers (eotd.modList, sep);

  eotd.organism_id_profile = NULL;
  GetOrganismIDProfile (sep, &eotd);

  if (SelectModifiersForExport (&eotd))
  {
    if (GetOutputFileName (path, sizeof (path), NULL))
    {
      eotd.fp = FileOpen (path, "w");
      if (eotd.fp == NULL) 
      {
        Message (MSG_ERROR, "Unable to open %s", path);
      }
	    else
	    {
	      lead_str [0] = 0;
		    lead_str [1] = 0;
        /* print a header line */
	      if (eotd.list_accession)
		    {
		      fprintf (eotd.fp, "Accession Number");
		      lead_str [0] = '\t';
		    }
		    if (eotd.list_local)
		    {
  	          fprintf (eotd.fp, "%sLocal ID", lead_str);
		      lead_str [0] = '\t';
		    }
        if (eotd.list_general)
		    {
		      fprintf (eotd.fp, "%sGeneral ID", lead_str);
		      lead_str [0] = '\t';
		    }
		    if (eotd.list_tax_name)
		    {
		      fprintf (eotd.fp, "%sTax Name", lead_str);
		      lead_str [0] = '\t';
		    }

        for (idx = 0; idx < NumDefLineModifiers (); idx++)
        {
          if (eotd.modList[idx].any_present)
          {
            fprintf (eotd.fp, "%s%s", lead_str, DefLineModifiers [idx].name);
            lead_str [0] = '\t';
          }
        }
        fprintf (eotd.fp, "\n");
        
        ExportOrganisms  (sep, &eotd);
        FileClose (eotd.fp);
#ifdef WIN_MSWIN
        if (eotd.launch_excel) {
          Nlm_MSWin_OpenApplication ("excel.exe", path);
        }
#endif
	    }
	  }
  }
  for (idx=0; idx < NumDefLineModifiers (); idx++)
  {
    ValNodeFree (eotd.modList[idx].values_seen);
  }
  eotd.modList = MemFree (eotd.modList);
  eotd.organism_id_profile = ValNodeFree (eotd.organism_id_profile);
}


typedef struct qualifierselect {
  LisT gene_quals;
  LisT cds_quals;
  LisT import_quals;
} QualifierSelectData, PNTR QualifierSelectPtr;

typedef struct featurequalloadformdata {
  FEATURE_FORM_BLOCK

  ValNodePtr             line_list;
  ValNodePtr             featlist;
  Int4                   num_feats;
  Uint2                  entityID;
  LisT                   feature_type;
  QualifierSelectData    match_qual;
  QualifierSelectData    apply_qual;
  PopuP                  sequence_id_type;
  PopuP                  sequence_column;
  PopuP                  feature_column;
  PopuP                  apply_column;
  GrouP                  seq_group;
  GrouP                  match_group;
  GrouP                  apply_group;
  Boolean                asked_about_replace;
  Boolean                do_replace;
  Boolean                use_semicolon;
  ButtoN                 replace_with_blank_btn;
  Boolean                replace_with_blank;
  ButtoN                 accept_button;
} FeatureQualLoadFormData, PNTR FeatureQualLoadFormPtr;
 
static void CleanupFeatureQualLoadForm (
  GraphiC g,
  VoidPtr data
)
{
  FeatureQualLoadFormPtr form_data;

  form_data = (FeatureQualLoadFormPtr)data;
  if (form_data == NULL) return;
  CleanUpTableData (form_data->line_list);
  StdCleanupFormProc (g, data);
} 

typedef struct qualifierlistdata {
  CharPtr item_name;
} QualifierListData, PNTR QualifierListPtr;

static Int4 GetSubtypeFromOffsetInFeatureList (Int4 list_offset,
                                               ValNodePtr list)
{
  ValNodePtr vnp;

  if (list == NULL || list_offset < 1) return -1;

  vnp = list;
  while (vnp != NULL && list_offset > 1)
  {
    list_offset --;
    vnp = vnp->next;
  }
  if (list_offset > 1) return -1;
  return vnp->choice;
}

static void ShowQualifierListBySubtype (QualifierSelectPtr qsp, Int4 feature_subtype)
{
  if (feature_subtype == -1) {
    Hide (qsp->gene_quals);
    Hide (qsp->cds_quals);
    Hide (qsp->import_quals);
  } else if (feature_subtype == FEATDEF_GENE) {
    Show (qsp->gene_quals);
    Hide (qsp->cds_quals);
    Hide (qsp->import_quals);
  } else if (feature_subtype == FEATDEF_CDS) {
    Hide (qsp->gene_quals);
    Show (qsp->cds_quals);
    Hide (qsp->import_quals);
  } else {
    Hide (qsp->gene_quals);
    Hide (qsp->cds_quals);
    Show (qsp->import_quals);
  }
}

static void ShowQualifierLists (FeatureQualLoadFormPtr form_data)
{
  Int4 feature_offset, feature_subtype;

  if (form_data == NULL) return;

  feature_offset = GetValue (form_data->feature_type);
  feature_subtype = GetSubtypeFromOffsetInFeatureList (feature_offset,
                                                       form_data->featlist);
  ShowQualifierListBySubtype (&(form_data->match_qual), feature_subtype);
  ShowQualifierListBySubtype (&(form_data->apply_qual), feature_subtype);
}

static QualifierListData GeneTableQualifiers [] = {
  { "allele"},
  { "locus" },
  { "locus_tag"},
  { "description"},
  { "product"}
};

#define NUM_GENE_TABLE_QUALIFIERS 5

static QualifierListData ProtTableQualifiers [] = {
  { "product"},
  { "description"},
  { "comment"}
};

#define NUM_PROT_TABLE_QUALIFIERS 3

/* check the appropriate list box for qualifier type,
 * keeping in mind that the first item is always "none"
 */
static Int4 GetQualSelectBySubtype (QualifierSelectPtr qsp, Int4 feature_subtype)
{
  Int4 qual_choice = -1;

  if (feature_subtype == -1) {
    qual_choice = -1;
  } else if (feature_subtype == FEATDEF_GENE) {
    qual_choice = GetValue (qsp->gene_quals);
    if (qual_choice < 2 
      || qual_choice > NUM_GENE_TABLE_QUALIFIERS + 2) {
      qual_choice = -1;
    } else {
      qual_choice = qual_choice - 1;
    }
  } else if (feature_subtype == FEATDEF_CDS) {
    qual_choice = GetValue (qsp->cds_quals);
    if (qual_choice < 2 
      || qual_choice > NUM_PROT_TABLE_QUALIFIERS + 2) {
      qual_choice = -1;
    } else {
      qual_choice = qual_choice - 1;
    }
  } else {
    qual_choice = GetValue (qsp->import_quals);
    if (qual_choice < 2
      || qual_choice > ParFlat_TOTAL_GBQUAL + 2) {
      qual_choice = -1;
    } else {
      qual_choice = qual_choice - 1;
    }
  }
  return qual_choice;

}

static void SetFeatureQualAcceptButton (Handle a)
{
  FeatureQualLoadFormPtr form_data;
  Boolean                have_feature_select_column;
  Boolean                have_feature_qual_select_column;
  Int4                   feature_select_choice;
  Int4                   apply_qual_choice;
  Int4                   match_qual_choice;
  Int4                   feature_subtype;
  Boolean                seqid_select_ok;
  Int4                   seqid_select_type;
  Int4                   seqid_select_column;

  form_data = GetObjectExtra (a);
  if (form_data == NULL) return;

  ShowQualifierLists (form_data);
  have_feature_qual_select_column = FALSE;
  have_feature_select_column = FALSE;
  seqid_select_ok = FALSE;

  feature_select_choice = GetValue (form_data->feature_type);
  if ( feature_select_choice > 0 && feature_select_choice < form_data->num_feats)
  {
    have_feature_select_column = TRUE;
  }
  
  feature_subtype = GetSubtypeFromOffsetInFeatureList (feature_select_choice,
                                                       form_data->featlist);
  
  apply_qual_choice = GetQualSelectBySubtype (&(form_data->apply_qual),
                                              feature_subtype);
  if (apply_qual_choice != -1)
  {
    have_feature_qual_select_column = TRUE;
  }

  /* if a column is selected for a sequence identifier, must select type */
  /* if a column is not selected for a sequence identifier, must not select type */
  seqid_select_type = GetValue (form_data->sequence_id_type);
  seqid_select_column = GetValue (form_data->sequence_column);
  if (seqid_select_type == 1 || seqid_select_type == 2) {
    if (seqid_select_column > 0) {
      seqid_select_ok = TRUE;
    }
  } else {
    seqid_select_ok = TRUE;
  }
 
  /* if the feature column is not specified, must not have match qualifier */
  /* and must have sequence ID */
  /* if the feature column is specified, must have match qualifier */
  match_qual_choice = GetQualSelectBySubtype (&(form_data->match_qual),
                                              feature_subtype);
  if (GetValue (form_data->feature_column) < 1)
  {
    if (match_qual_choice != -1)
    {
      have_feature_select_column = FALSE;
    }
    else if (seqid_select_column < 1)
    {
      have_feature_select_column = FALSE;
    }
  } else {
    if (match_qual_choice == -1)
    {
      have_feature_select_column = FALSE;
    }
  }

  if (GetValue (form_data->apply_column) < 1) {
    have_feature_qual_select_column = FALSE;
  }

  if ( have_feature_select_column && have_feature_qual_select_column
    && seqid_select_ok)
  {
    Enable (form_data->accept_button);
  }
  else
  {
    Disable (form_data->accept_button);
  }
}

static Boolean ProductNameMatches (
  SeqLocPtr slp,
  BioseqPtr  bsp,
  CharPtr   pString
)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

  if (slp == NULL || bsp == NULL || pString == NULL) return FALSE;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL)
  {
    if ( IsLocAInBonSameStrand (sfp->location, slp)
      && StringStr ( fcontext.label, pString) != NULL)
    {
      return TRUE;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }
  return FALSE;
}

static Boolean GeneHasThisQualifier (
  SeqFeatPtr sfp,
  Uint2      entityID,
  CharPtr    pString,
  Int2       popup_offset
)
{
  BioseqPtr  bsp;
  GeneRefPtr grp;

  if (sfp == NULL || pString == NULL) return FALSE;
  if (popup_offset < 1 || popup_offset > NUM_GENE_TABLE_QUALIFIERS) return FALSE;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return FALSE;

  if ( StringCmp (GeneTableQualifiers [popup_offset - 1].item_name, "allele") == 0) {
    if (StringCmp (pString, grp->allele) == 0) {
      return TRUE;
    } else {
      return FALSE;
    }
  } else if ( StringCmp (GeneTableQualifiers [popup_offset - 1].item_name, "locus") == 0) {
    if (StringCmp (pString, grp->locus) == 0) {
      return TRUE;
    } else {
      return FALSE;
    }
  } else if ( StringCmp (GeneTableQualifiers [popup_offset - 1].item_name, "locus_tag") == 0) {
    if (StringCmp (pString, grp->locus_tag) == 0) {
      return TRUE;
    } else {
      return FALSE;
    }
  } else if ( StringCmp (GeneTableQualifiers [popup_offset - 1].item_name, "description") == 0) {
    if (StringCmp (pString, grp->desc) == 0) {
      return TRUE;
    } else {
      return FALSE;
    }
  } else if ( StringCmp (GeneTableQualifiers [popup_offset - 1].item_name, "product") == 0)
  {
    bsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
    return ProductNameMatches (sfp->location, bsp, pString);
  }
  return FALSE;
}

static Boolean ProteinHasThisQualifier (
  SeqFeatPtr sfp,
  CharPtr    pString,
  Int2       popup_offset
)
{
  ProtRefPtr prp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || pString == NULL
    || popup_offset < 1 || popup_offset > NUM_PROT_TABLE_QUALIFIERS)
  {
    return FALSE;
  }

  prp = sfp->data.value.ptrvalue;
  if (prp == NULL) return FALSE;

  if (StringCmp (ProtTableQualifiers [ popup_offset - 1].item_name, "product") == 0)
  {
    if (prp->name != NULL
      && StringStr (prp->name->data.ptrvalue, pString) != NULL)
    {
      return TRUE;
    }
  } else if (StringCmp (ProtTableQualifiers [popup_offset - 1].item_name, "description") == 0) {
    if (StringStr (prp->desc, pString) != NULL)
    {
      return TRUE;
    }
  } else if (StringCmp (ProtTableQualifiers [popup_offset - 1].item_name, "comment") == 0) {
    if (StringStr (sfp->comment, pString) != NULL)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean CDSHasThisQualifier (
  SeqFeatPtr sfp,
  Uint2      entityID,
  CharPtr    pString,
  Int2       popup_offset
)
{
  SeqFeatPtr protein;

  if (sfp == NULL || pString == NULL) return FALSE;

  if (StringCmp (ProtTableQualifiers [popup_offset - 1].item_name, "comment") == 0) {
    if (StringStr (sfp->comment, pString) != NULL)
    {
      return TRUE;
    }
  }
  else
  {
    protein = FindBestProtein (entityID, sfp->product);
    if (protein == NULL) return FALSE;
    return ProteinHasThisQualifier (protein, pString, popup_offset);
  }
  return FALSE;
}
  
static Boolean FeatureHasThisQualifier (
  SeqFeatPtr sfp,
  Uint2      entityID,
  CharPtr    pString,
  Int2       popup_offset
)
{
  Int2      qualval;
  GBQualPtr gbqual;

  if (popup_offset < 1) return FALSE;
  if (sfp == NULL || pString == NULL) return FALSE;

  if (sfp->idx.subtype == FEATDEF_GENE)
  {
    return GeneHasThisQualifier (sfp, entityID, pString, popup_offset);
  }
  else if (sfp->idx.subtype == FEATDEF_CDS)
  {
    return CDSHasThisQualifier (sfp, entityID, pString, popup_offset);
  }
  else
  {
    for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
    {
      qualval = GBQualNameValid (gbqual->qual);
      if (qualval > -1 && qualval == popup_offset -1)
      {
        return TRUE;
      }
    }
  }
  return FALSE; 
}
 
static CharPtr GetDataForColumnOffset (ValNodePtr column_data, Int4 offset)
{
  ValNodePtr vnp;

  for (vnp = column_data; vnp != NULL && offset > 1; vnp = vnp->next)
  {
    offset --;
  }
  if (offset > 1 || vnp == NULL) {
    return NULL;
  } else {
    return vnp->data.ptrvalue;
  }
}

static Boolean ApplyQualToThisFeature (
  SeqFeatPtr sfp,
  ValNodePtr parts,
  FeatureQualLoadFormPtr form_data
)
{
  Int4       column_index;
  CharPtr    val;
  Int4       subtype;
  Int4       popup_offset;

  if (sfp == NULL || parts == NULL || form_data == NULL) return FALSE;

  column_index = GetValue (form_data->feature_column);
  subtype = GetSubtypeFromOffsetInFeatureList (
                     GetValue (form_data->feature_type), form_data->featlist);
  val = GetDataForColumnOffset (parts, column_index);

  if (val != NULL && sfp->idx.subtype == subtype)
  {
    popup_offset = GetQualSelectBySubtype (&(form_data->match_qual),
                                              subtype);
    if (popup_offset < 1
      || FeatureHasThisQualifier (sfp, form_data->entityID, val, popup_offset))
    {
      return TRUE;
    }
  }
  return FALSE;
}

typedef struct featurequaltabledata {
  ValNodePtr parts;
  FeatureQualLoadFormPtr form_data;
} FeatureQualTableData, PNTR FeatureQualTablePtr;

static void ApplyOneQualToProt (
  SeqFeatPtr sfp,
  Int2       popup_offset,
  CharPtr    pString,
  Boolean PNTR asked_about_replace,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon)
{
  ProtRefPtr prp;
  CharPtr    cp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || pString == NULL
    || popup_offset < 1 || popup_offset > NUM_PROT_TABLE_QUALIFIERS
    || asked_about_replace == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  prp = sfp->data.value.ptrvalue;
  if (prp == NULL) return;

  if (StringCmp (ProtTableQualifiers [ popup_offset - 1].item_name, "product") == 0)
  {
    if (prp->name == NULL)
    {
      ValNodeAddStr (&prp->name, 0, StringSave (pString));
    }
    else
    {
      cp = prp->name->data.ptrvalue;
      AppendOrReplaceString (&cp, pString,
                            asked_about_replace, do_replace, use_semicolon);
      prp->name->data.ptrvalue = cp;
    }
  } else if (StringCmp (ProtTableQualifiers [popup_offset - 1].item_name, "description") == 0) {
    AppendOrReplaceString (&(prp->desc), pString,
                            asked_about_replace, do_replace, use_semicolon);
  } else if (StringCmp (ProtTableQualifiers [popup_offset - 1].item_name, "comment") == 0) {
    AppendOrReplaceString (&(sfp->comment), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  }
}

static void ApplyOneQualToCDS (
  SeqFeatPtr sfp,
  Uint2      entityID,
  Int2       feat_qual_choice,
  CharPtr    pString,
  Boolean PNTR asked_about_replace,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon)
{
  SeqFeatPtr protein;

  if (sfp == NULL || pString == NULL
    || asked_about_replace == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  if (StringCmp (ProtTableQualifiers [feat_qual_choice - 1].item_name, "comment") == 0) {
    AppendOrReplaceString (&(sfp->comment), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  }
  else
  {
    protein = FindBestProtein (entityID, sfp->product);
    if (protein == NULL) return;
    ApplyOneQualToProt (protein, feat_qual_choice, pString, 
                            asked_about_replace, do_replace, use_semicolon);
  }
}

static void ApplyGeneProductName (
  SeqLocPtr slp,
  Uint2     entityID,
  CharPtr   pString,
  Boolean PNTR asked_about_replace,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  Int4              prot_popup_offset;
  BioseqPtr         bsp;

  if (slp == NULL || pString == NULL
    || asked_about_replace == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  bsp = GetBioseqGivenSeqLoc (slp, entityID);
  if (bsp == NULL) return;

  for (prot_popup_offset = 0;
       prot_popup_offset < NUM_PROT_TABLE_QUALIFIERS
         && StringCmp (ProtTableQualifiers [prot_popup_offset].item_name, "product") != 0;
       prot_popup_offset ++)
  {}
  if (prot_popup_offset >= NUM_PROT_TABLE_QUALIFIERS) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL)
  {
    if ( IsLocAInBonSameStrand (sfp->location, slp))
    {
      ApplyOneQualToCDS (sfp, entityID, prot_popup_offset + 1, pString, 
                            asked_about_replace, do_replace, use_semicolon);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }
}

static void ApplyOneQualToGene (
  SeqFeatPtr sfp,
  Uint2      entityID,
  Int2       feat_qual_choice,
  CharPtr    pString, 
  Boolean PNTR asked_about_replace,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon)
{
  GeneRefPtr grp;

  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_GENE
    || feat_qual_choice < 1
    || feat_qual_choice > NUM_GENE_TABLE_QUALIFIERS
    || pString == NULL
    || asked_about_replace == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;

  if (StringCmp ( GeneTableQualifiers [feat_qual_choice - 1].item_name,
                 "allele") == 0) {
    AppendOrReplaceString (&(grp->allele), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  } else if (StringCmp ( GeneTableQualifiers [feat_qual_choice - 1].item_name,
                        "locus") == 0) {
    AppendOrReplaceString (&(grp->locus), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  } else if (StringCmp ( GeneTableQualifiers [feat_qual_choice - 1].item_name,
                        "locus_tag") == 0) {
    AppendOrReplaceString (&(grp->locus_tag), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  } else if (StringCmp ( GeneTableQualifiers [feat_qual_choice - 1].item_name,
                        "description") == 0) {
    AppendOrReplaceString (&(grp->desc), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  } else if (StringCmp (GeneTableQualifiers [feat_qual_choice - 1].item_name,
                        "product") == 0) {
    ApplyGeneProductName (sfp->location, entityID, pString, 
                            asked_about_replace, do_replace, use_semicolon);
  }
}

static void ApplyOneQualToImportFeature (
  SeqFeatPtr sfp,
  Int2 feat_qual_choice,
  CharPtr pString,
  Boolean PNTR asked_about_replace,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon)
{
  GBQualPtr           gbqual;

  if (sfp == NULL
    || feat_qual_choice < 1
    || feat_qual_choice > ParFlat_TOTAL_GBQUAL
    || pString == NULL
    || asked_about_replace == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  gbqual = sfp->qual;
  while (gbqual != NULL
    && StringCmp (gbqual->qual,
                  ParFlat_GBQual_names [feat_qual_choice - 1].name) != 0)
  {
    gbqual = gbqual->next;
  }
  
  if (gbqual == NULL)
  {
    gbqual = GBQualNew ();
    if (gbqual == NULL) return;
    gbqual->qual = StringSave ( ParFlat_GBQual_names [feat_qual_choice - 1].name);
    gbqual->val = StringSave ( pString );
    gbqual->next = sfp->qual;
    sfp->qual = gbqual;
  }
  else
  {
    AppendOrReplaceString (&(gbqual->val), pString, 
                            asked_about_replace, do_replace, use_semicolon);
  }
}

static void ApplyOneQualToFeature (
  SeqFeatPtr sfp,
  Uint2      entityID,
  Int2       feat_qual_choice,
  CharPtr    pString,
  Boolean PNTR asked_about_replace,
  Boolean PNTR do_replace,
  Boolean PNTR use_semicolon)
{
  if (sfp == NULL || pString == NULL 
    || asked_about_replace == NULL
    || do_replace == NULL
    || use_semicolon == NULL)
  {
    return;
  }

  switch (sfp->idx.subtype)
  {
    case FEATDEF_GENE :
      ApplyOneQualToGene (sfp, entityID, feat_qual_choice, pString,
                            asked_about_replace, do_replace, use_semicolon);
      break;
    case FEATDEF_CDS :
      ApplyOneQualToCDS (sfp, entityID, feat_qual_choice, pString,
                            asked_about_replace, do_replace, use_semicolon);
      break;
    default :
      ApplyOneQualToImportFeature (sfp, feat_qual_choice, pString,
                            asked_about_replace, do_replace, use_semicolon);
      break;
  }
}

static CharPtr GetImpFeatKeyFromSubtype (Uint2 subtype)
{
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label;

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key == subtype)
    {
      return curr->typelabel;
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  return NULL;
}

static void ApplyQualsToFeaturesOnBsp (BioseqPtr bsp,
 FeatureQualLoadFormPtr form_data,
 ValNodePtr             parts)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  Boolean found_feature = FALSE;
  Int4                subtype;
  Int4                column_index;
  Int2                qualval;
  Int2                matchval;
  CharPtr             val;
  ImpFeatPtr          ifp;

  subtype = GetSubtypeFromOffsetInFeatureList (
                          GetValue (form_data->feature_type), 
                          form_data->featlist);
  column_index = GetValue (form_data->apply_column);
  val = GetDataForColumnOffset (parts, column_index);
  /* if the value to be applied is blank and we're not replacing with blanks, skip */
  if (! form_data->replace_with_blank && StringHasNoText (val)) {
    return;
  }
  qualval = GetQualSelectBySubtype (&(form_data->apply_qual), subtype);
  if (subtype == -1 || val == NULL || column_index < 1 || qualval < 1) {
    /* do nothing */
  }
  else
  {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (sfp != NULL)
    {
      if ( ApplyQualToThisFeature (sfp, parts, form_data))
      {
        found_feature = TRUE;
        ApplyOneQualToFeature (sfp, form_data->entityID, qualval, val,
                               &(form_data->asked_about_replace),
                               &(form_data->do_replace),
                               &(form_data->use_semicolon));
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
    }

    /* only want to create new features if no feature matching qualifier
     * was specified and no feature was found */
    matchval = GetQualSelectBySubtype (&(form_data->match_qual), subtype);

    if (! found_feature && matchval == -1)
    {
      if (subtype == FEATDEF_CDS)
      {
        sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_CDREGION, NULL);
        if (sfp == NULL) return;
        sfp->data.value.ptrvalue = CdRegionNew ();
      } else if (subtype == FEATDEF_GENE) {
        sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, NULL);
        if (sfp == NULL) return;
        sfp->data.value.ptrvalue = GeneRefNew ();
      } else {
        ifp = ImpFeatNew ();
        if (ifp == NULL) return;
        sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, NULL);
        if (sfp == NULL) return;
        ifp->key = StringSave ( GetImpFeatKeyFromSubtype (subtype));
        sfp->data.value.ptrvalue = ifp;
      }
      SeqMgrIndexFeatures (0, (Pointer) bsp);
      ApplyOneQualToFeature (sfp, form_data->entityID, qualval, val,
                               &(form_data->asked_about_replace),
                               &(form_data->do_replace),
                               &(form_data->use_semicolon));
    }
  }
}

static void ApplyTableQuals (BioseqPtr bsp, Pointer userdata)
{
  FeatureQualLoadFormPtr form_data;
  SeqDescrPtr            sdp;
  GBBlockPtr             gbp;
  Int4                   column_index;
  ValNodePtr             line;
  TableLinePtr           tlp;
  Boolean                use_local_id = FALSE;
  FeatureQualTableData   fqtd;
  Int2                   seq_match_choice;
  CharPtr                idval;
  
  form_data = (FeatureQualLoadFormPtr) userdata;
  if (form_data == NULL || bsp == NULL) return;

  gbp = NULL;
  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);
  if (sdp != NULL)
  {
    gbp = sdp->data.ptrvalue;
  }

  seq_match_choice = GetValue (form_data->sequence_id_type);
  if (seq_match_choice == 1) {
    use_local_id = FALSE;
  }

  fqtd.form_data = form_data;
  if (seq_match_choice == 1 || seq_match_choice == 2)
  {
    column_index = GetValue (form_data->sequence_column);

    for (line = form_data->line_list; line != NULL; line = line->next)
    {
      tlp = line->data.ptrvalue;
      if (tlp == NULL) continue;
      idval = GetDataForColumnOffset (tlp->parts, column_index);
      if (idval != NULL
        && ( IDListHasValue ( idval,
                              bsp->id, use_local_id, FALSE, FALSE)
          || (! use_local_id && 
              HasExtraAccession ( idval, gbp))))
      {
        fqtd.parts = tlp->parts;
        ApplyQualsToFeaturesOnBsp (bsp, form_data, tlp->parts);
      }
    }
  } else {
    for (line = form_data->line_list; line != NULL; line = line->next)
    {
      tlp = line->data.ptrvalue;
      if (tlp == NULL) continue;
      fqtd.parts = tlp->parts;
      ApplyQualsToFeaturesOnBsp (bsp, form_data, tlp->parts);
    }
  }
}

static void DoAcceptFeatureQuals (ButtoN b)
{
  FeatureQualLoadFormPtr form_data;
  SeqEntryPtr   sep;

  form_data = GetObjectExtra (b);
  if (form_data == NULL) return;

  form_data->replace_with_blank = GetStatus (form_data->replace_with_blank_btn);

  sep = GetTopSeqEntryForEntityID (form_data->entityID);
  if (sep == NULL) return;
  
  VisitBioseqsInSep (sep, form_data, ApplyTableQuals);
  Update ();
  ObjMgrSetDirtyFlag (form_data->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, form_data->entityID, 0, 0);
  Remove (form_data->form);
}

static PopuP BuildColumnSelector ( GrouP g,
                                  FeatureQualLoadFormPtr parent_form,
                                  ValNodePtr parts )
{
  PopuP      p;
  ValNodePtr this_part;

  if (g == NULL || parts == NULL) return NULL;

  p = PopupList (g, TRUE, (PupActnProc) SetFeatureQualAcceptButton);
  for (this_part = parts; this_part != NULL; this_part = this_part->next)
  {
    if (this_part->data.ptrvalue != NULL)
    {
      PopupItem (p, this_part->data.ptrvalue);
    }
  }
  SetObjectExtra (p, parent_form, NULL);
  return p;
}

static void BuildQualifierSelect ( GrouP parent, 
                                  FeatureQualLoadFormPtr parent_form,
                                  ValNodePtr parts,
                                  QualifierSelectPtr qsp )
{
  GrouP g;
  Int4  listheight = 4;
  Int4  listwidth  = 13;
  Int4  j;

  if (parent_form == NULL || parts == NULL || qsp == NULL)
  {
    return;
  }

  g = HiddenGroup (parent, 0, 0, NULL);

  qsp->gene_quals = SingleList (g, listwidth, listheight,
                                   (LstActnProc) SetFeatureQualAcceptButton);
  SetObjectExtra (qsp->gene_quals, parent_form, NULL);
  ListItem (qsp->gene_quals, "None");
  for (j=0; j < NUM_GENE_TABLE_QUALIFIERS; j++) {
    ListItem (qsp->gene_quals, GeneTableQualifiers [j].item_name);
  }
  SetValue (qsp->gene_quals, 1);

  qsp->cds_quals = SingleList (g, listwidth, listheight,
                                   (LstActnProc) SetFeatureQualAcceptButton);
  SetObjectExtra (qsp->cds_quals, parent_form, NULL);
  ListItem (qsp->cds_quals, "None");
  for (j=0; j < NUM_PROT_TABLE_QUALIFIERS; j++) {
    ListItem (qsp->cds_quals, ProtTableQualifiers [j].item_name);
  }
  SetValue (qsp->cds_quals, 1);

  qsp->import_quals = SingleList (g, listwidth, listheight,
                                   (LstActnProc) SetFeatureQualAcceptButton);
  SetObjectExtra (qsp->import_quals, parent_form, NULL);
  ListItem (qsp->import_quals, "None");
  for (j = 0; j < ParFlat_TOTAL_GBQUAL; j++) {
    ListItem (qsp->import_quals, ParFlat_GBQual_names [j].name);
  }
  SetValue (qsp->import_quals, 1);
}

static void BuildSequenceSelect ( GrouP parent,
                                  FeatureQualLoadFormPtr parent_form,
                                  ValNodePtr parts )
{
  if (parent == NULL || parent_form == NULL || parts == NULL) return;

  parent_form->seq_group = NormalGroup ( parent, 0, 2, "Sequence", programFont, NULL);
  StaticPrompt (parent_form->seq_group, "Identifier Type", 0, popupMenuHeight, programFont, 'l');
  parent_form->sequence_id_type = PopupList (parent_form->seq_group, TRUE, (PupActnProc) SetFeatureQualAcceptButton);
  PopupItem (parent_form->sequence_id_type, "Accession");
  PopupItem (parent_form->sequence_id_type, "Local ID");
  PopupItem (parent_form->sequence_id_type, "None");
  SetObjectExtra (parent_form->sequence_id_type, parent_form, NULL);
  SetValue (parent_form->sequence_id_type, 3);

  StaticPrompt (parent_form->seq_group, "Column", 0, popupMenuHeight, programFont, 'l');
  parent_form->sequence_column = BuildColumnSelector ( parent_form->seq_group, parent_form, parts);
}

static void BuildFeatureSelect ( GrouP parent,
                                 FeatureQualLoadFormPtr parent_form,
                                 ValNodePtr parts )
{
  Int4  listheight = 4;
  ValNodePtr vnp;

  if (parent == NULL || parent_form == NULL || parts == NULL) return;

  parent_form->match_group = NormalGroup ( parent, 0, 2,
                                    "Feature to Edit", programFont, NULL);
  StaticPrompt (parent_form->match_group, "Feature type", 0, 
                                   popupMenuHeight, programFont, 'l');
  parent_form->feature_type = SingleList (parent_form->match_group, 
                                   9, listheight,
                   (LstActnProc) SetFeatureQualAcceptButton);
  SetObjectExtra (parent_form->feature_type, parent_form, NULL);
  for (vnp = parent_form->featlist; vnp != NULL; vnp = vnp->next)
  {
    ListItem (parent_form->feature_type, (CharPtr) vnp->data.ptrvalue);
  }

  StaticPrompt (parent_form->match_group, "Qualifier to match", 0, popupMenuHeight, programFont, 'l');
  BuildQualifierSelect ( parent_form->match_group, parent_form, parts, &(parent_form->match_qual));

  StaticPrompt (parent_form->match_group, "Column", 0, popupMenuHeight, programFont, 'l');
  parent_form->feature_column = BuildColumnSelector ( parent_form->match_group, parent_form, parts);
  
}

static void BuildFeatureApply ( GrouP parent,
                               FeatureQualLoadFormPtr parent_form,
                               ValNodePtr parts )
{
  if (parent == NULL || parent_form == NULL || parts == NULL) return;

  parent_form->apply_group = NormalGroup ( parent, 0, 2, "Qualifier to Edit", programFont, NULL);
  
  StaticPrompt (parent_form->apply_group, "Qualifier Name", 0, popupMenuHeight, programFont, 'l');
  BuildQualifierSelect ( parent_form->apply_group, parent_form, parts, &(parent_form->apply_qual));

  StaticPrompt (parent_form->apply_group, "Column", 0, popupMenuHeight, programFont, 'l');
  parent_form->apply_column = BuildColumnSelector ( parent_form->apply_group, parent_form, parts);
  
}

extern void LoadFeatureQualifierTable (IteM i)
{
  BaseFormPtr   bfp;
  ValNodePtr    header_line;
  ValNodePtr    vnp;
  TableLinePtr  tlp;
  WindoW        w;
  GrouP         h, g, c;
  FeatureQualLoadFormPtr form_data;
  Int4          max_columns;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  header_line = ReadTableData ();
  if (header_line == NULL || header_line->data.ptrvalue == NULL) return;
  tlp = header_line->data.ptrvalue;
  max_columns = 0;
  for (vnp = tlp->parts; vnp != NULL; vnp = vnp->next)
  {
    max_columns ++;
  }
  
  form_data = MemNew (sizeof (FeatureQualLoadFormData));
  if (form_data == NULL) return;
  form_data->asked_about_replace = FALSE;
  form_data->do_replace = FALSE;
  form_data->use_semicolon = FALSE;
  form_data->entityID = bfp->input_entityID;
  form_data->line_list = header_line;
  form_data->featlist = BuildFeatureValNodeList (TRUE, NULL, 0, TRUE, FALSE);
  form_data->num_feats = 0;
  for (vnp = form_data->featlist; vnp != NULL; vnp = vnp->next)
  {
    form_data->num_feats ++;
  }

  /* now create a dialog to display values */
  w = FixedWindow (-50, -33, -10, -10, "Table Conversion", StdCloseWindowProc);
  SetObjectExtra (w, form_data, CleanupFeatureQualLoadForm);
  form_data->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  BuildSequenceSelect ( h, form_data, tlp->parts);
  BuildFeatureSelect (h, form_data, tlp->parts);
  BuildFeatureApply (h, form_data, tlp->parts);
  ShowQualifierLists (form_data);

  g = HiddenGroup (h, 1, 0, NULL);
  form_data->replace_with_blank_btn = CheckBox (g, "Erase current value when blank found in table", NULL);
  c = HiddenGroup (h, 4, 0, NULL);
  form_data->accept_button = DefaultButton (c, "Accept", DoAcceptFeatureQuals);
  SetObjectExtra (form_data->accept_button, form_data, NULL);
  Disable (form_data->accept_button);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER,
                (HANDLE) form_data->seq_group, 
                (HANDLE) form_data->match_group,
                (HANDLE) form_data->apply_group,
                (HANDLE) c, NULL);

  RealizeWindow (w);
  Show (w);
  Update ();
}


typedef struct qualloadformdata {
  FEATURE_FORM_BLOCK
  ValNodePtr             table;
  DialoG                 list_dlg;
  DialoG PNTR            column_list;
  Int4                   num_columns;
  ButtoN                 remove_quotes;
  ButtoN                 accept_button;
  ButtoN                 mail_button;
  ButtoN                 erase_all_when_blank;  /* J. Chen */
  PrompT                 err_msg;
} QualLoadFormData, PNTR QualLoadFormPtr;

static void CleanupQualLoadForm (
  GraphiC g,
  VoidPtr data
)
{
  QualLoadFormPtr form_data;

  form_data = (QualLoadFormPtr)data;
  if (form_data == NULL) return;
  form_data->table = FreeTabTable (form_data->table);
  StdCleanupFormProc (g, data);
} 


static void ChangeTabColumnChoice (Pointer data)
{
  QualLoadFormPtr form_data;
  Int4            i;
  Boolean         have_match = FALSE, have_apply = FALSE;
  TabColumnConfigPtr t;
  ValNodePtr         err_list = NULL;

  if (data == NULL) return;

  form_data = (QualLoadFormPtr) data;
  SetTitle (form_data->err_msg, "");
  /* must have one and only one match choice, must have at least one apply choice */

  if (form_data->list_dlg == NULL) {
    for (i = 0; i < form_data->num_columns; i++) {
      t = (TabColumnConfigPtr) DialogToPointer (form_data->column_list[i]);
      if (t != NULL) {
        if (t->match_type != NULL) {
          if (have_match) {
            Disable (form_data->accept_button);
            Disable (form_data->mail_button);
            return;
          } else {
            have_match = TRUE;
          }
        } else if (!IsFieldTypeEmpty(t->field)) {
          have_apply = TRUE;
        }
      }
      t = TabColumnConfigFree (t);
    }
    if (have_match && have_apply) {
      Enable (form_data->accept_button);
      Enable (form_data->mail_button);
    } else {
      Disable (form_data->accept_button);
      Disable (form_data->mail_button);
    }
  } else {
    err_list = TestDialog (form_data->list_dlg);
    if (err_list == NULL) {
      Enable (form_data->accept_button);
      Enable (form_data->mail_button);
    } else {
      Disable (form_data->accept_button);
      Disable (form_data->mail_button);
      SetTitle (form_data->err_msg, err_list->data.ptrvalue);
      err_list = ValNodeFree (err_list);
    }
  }
}


static Boolean ApplyTableValues (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list, vnp, obj_table, dup_dest_errs;
  LogInfoPtr lip;
  Int4       msg_len = 0;
  CharPtr    msg = NULL, msg_end = "Continue with errors?  Hit cancel to scroll through details";

  lip = OpenLog ("Table Problems");

  err_list = ValidateTabTableValues (table, columns);
  for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
    fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    lip->data_in_log = TRUE;
  }
  err_list = ValNodeFreeData (err_list);

  obj_table = GetObjectTableForTabTable (sep, table, columns, &err_list);

  dup_dest_errs = CheckObjTableForRowsThatApplyToTheSameDestination (obj_table);
  if (dup_dest_errs != NULL) {
    for (vnp = dup_dest_errs; vnp != NULL; vnp = vnp->next) {
      fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
      lip->data_in_log = TRUE;
    }
    CloseLog (lip);
    Message (MSG_ERROR, "For one or more columns, two or more rows in the table apply to the same object.  Cannot continue.");
    dup_dest_errs = ValNodeFreeData (dup_dest_errs);
    err_list = ValNodeFreeData (err_list);
    FreeLog (lip);
    obj_table = FreeObjectTableForTabTable (obj_table);
    DeleteMarkedObjects (SeqMgrGetEntityIDForSeqEntry (sep), 0, NULL);
    return FALSE;
  }

  ValNodeLink (&err_list, CheckObjTableForExistingText (sep, table, columns, obj_table));

  msg_len = StringLen (msg_end) + 1;
  /* cycle through errors twice - first time, just print the ones with choice 1 (and sum lengths) */
  for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) {
      fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
      lip->data_in_log = TRUE;
      msg_len += StringLen (vnp->data.ptrvalue) + 2;
    }
  }
  /* then cycle again, printing 0s */
  for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
      lip->data_in_log = TRUE;
    }
  }

  /* now produce error message */
  if (lip->data_in_log) {
    msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
    for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 1) {
        StringCat (msg, vnp->data.ptrvalue);
        StringCat (msg, "\n");
      }
    }
    StringCat (msg, msg_end);
  }
  
  err_list = ValNodeFreeData (err_list);

  CloseLog (lip);
  if (lip->data_in_log) {
    if (ANS_CANCEL == Message (MSG_OKC, msg)) {
      msg = MemFree (msg);
      FreeLog (lip);
      DeleteMarkedObjects (SeqMgrGetEntityIDForSeqEntry (sep), 0, NULL);
      return FALSE;
    }
  }
  msg = MemFree (msg);
  FreeLog (lip);

  err_list =  ApplyTableValuesToObjectTable (sep, table, columns, obj_table);

  lip = OpenLog ("Table Problems");
  for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
    fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    lip->data_in_log = TRUE;
  }
  err_list = ValNodeFreeData (err_list);
  CloseLog (lip);
  FreeLog (lip);
  obj_table = FreeObjectTableForTabTable (obj_table);
  DeleteMarkedObjects (SeqMgrGetEntityIDForSeqEntry (sep), 0, NULL);
  return TRUE;
}


static ValNodePtr CheckForDuplicateColumns (ValNodePtr columns)
{
  ValNodePtr vnp, vnp2;
  TabColumnConfigPtr t, t2;
  ValNodePtr duplicate_column_list = NULL;
  Boolean    this_dup;

  for (vnp = columns; vnp != NULL; vnp = vnp->next) {
    this_dup = FALSE;
    t = (TabColumnConfigPtr) vnp->data.ptrvalue;
    if (t != NULL && t->field != NULL) {
      for (vnp2 = vnp->next; vnp2 != NULL && !this_dup; vnp2 = vnp2->next) {
        t2 = (TabColumnConfigPtr) vnp2->data.ptrvalue;
        if (t2 != NULL && CompareFieldTypes(t->field, t2->field) == 0 && t2->existing_text != ExistingTextOption_add_qual) {
          ValNodeAddPointer (&duplicate_column_list, 0, SummarizeFieldType (t->field));
          this_dup = TRUE;
        }
      }
    }
  }
  duplicate_column_list = ValNodeSort (duplicate_column_list, SortVnpByString);
  ValNodeUnique (&duplicate_column_list, SortVnpByString, ValNodeFreeData);
  return duplicate_column_list;
}

  
static void FixDuplicateColumns (ValNodePtr columns)
{
  ValNodePtr vnp, vnp2;
  TabColumnConfigPtr t, t2;

  if (columns == NULL || columns->next == NULL) {
    return;
  }

  for (vnp = columns; vnp != NULL; vnp = vnp->next) {
    t = (TabColumnConfigPtr) vnp->data.ptrvalue;
    if (t->field != NULL) {
      for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
        t2 = (TabColumnConfigPtr) vnp2->data.ptrvalue;
        if (CompareFieldTypes(t->field, t2->field) == 0 && t2->existing_text != ExistingTextOption_add_qual) {
          t2->existing_text = ExistingTextOption_add_qual;
        }
      }
    }
  }
}


static Boolean HandleDuplicateColumns (QualLoadFormPtr form_data, ValNodePtr columns)
{
  ValNodePtr      dup_col_list;
  Int4            i;
  ValNodePtr      vnp;
  Boolean         rval = TRUE;

  if (form_data == NULL || columns == NULL) {
    return FALSE;
  }

  dup_col_list = CheckForDuplicateColumns (columns);
  if (dup_col_list != NULL) {
    if (ANS_CANCEL == Message (MSG_OKC, 
         "You have specified multiple columns for the same qualifier value.  This will probably not work as expected. "
         " Do you want to change the configuration to add them all as separate qualifiers?")) {
      rval = FALSE;
    } else {
      dup_col_list = ValNodeFreeData(dup_col_list);
      FixDuplicateColumns (columns);
      if (form_data->list_dlg == NULL) {
        for (i = 0, vnp = columns; i < form_data->num_columns && vnp != NULL; i++, vnp = vnp->next) {
          PointerToDialog (form_data->column_list[i], vnp->data.ptrvalue);
        }
      } else {
        PointerToDialog (form_data->list_dlg, columns);
      }
    }
  }
  dup_col_list = ValNodeFreeData(dup_col_list);

  return rval;
}


static void DoAcceptQuals (ButtoN b)
{
  QualLoadFormPtr form_data;
  ValNodePtr      columns = NULL;
  Int4            i;
  TabColumnConfigPtr t;
  SeqEntryPtr     sep;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;

  if (form_data->list_dlg == NULL) {
    for (i = 0; i < form_data->num_columns; i++) {
      t = DialogToPointer (form_data->column_list[i]);
      ValNodeAddPointer (&columns, 0, t);
    }
  } else {
    columns = DialogToPointer (form_data->list_dlg);
  }

  if (!HandleDuplicateColumns(form_data, columns)) {
    columns = TabColumnConfigListFree (columns);
    return;
  }

  if (GetStatus (form_data->remove_quotes)) {
    RemoveQuotesFromTabTable (form_data->table);
  }

  sep = GetTopSeqEntryForEntityID (form_data->input_entityID);
  if (ApplyTableValues (sep, form_data->table, columns)) {
    ObjMgrSetDirtyFlag (form_data->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, form_data->input_entityID, 0, 0);
    if (GetStatus (form_data->leave_dlg_up)) {
    } else {
      Remove (form_data->form);
    }
    Update();
  }
  columns = TabColumnConfigListFree (columns);

}


typedef struct mailreport {
  BioseqPtr    bsp;
  BioSourcePtr biop;
  CharPtr      old_name;
} MailReportData, PNTR MailReportPtr;


static MailReportPtr MailReportNew (BioseqPtr bsp, CharPtr old_name, BioSourcePtr biop)
{
  MailReportPtr m;

  m = (MailReportPtr) MemNew (sizeof (MailReportData));
  m->bsp = bsp;
  m->old_name = StringSave (old_name);
  m->biop = biop;
  return m;
}


static MailReportPtr MailReportFree (MailReportPtr m)
{
  if (m != NULL) {
    m->old_name = MemFree (m->old_name);
    m = MemFree (m);
  }
  return m;
}


static ValNodePtr MailReportListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = MailReportFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void MailReportCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;

  if (bsp == NULL || data == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &context)) {
    if ((biop = (BioSourcePtr) sdp->data.ptrvalue) != NULL
        && biop->org != NULL
        && !StringHasNoText (biop->org->taxname)) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, MailReportNew (bsp, biop->org->taxname, biop));
    }
  }
}


static void ReportOneMailReport (MailReportPtr m, BioSourcePtr biop, LogInfoPtr lip)
{
  Char              id_str[100];

  SeqIdWrite (SeqIdFindBest (m->bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
  fprintf (lip->fp, "%s\t%s\t%s\n", id_str, m->old_name == NULL ? "" : m->old_name, 
                 biop->org == NULL || biop->org->taxname == NULL ? "" : biop->org->taxname);
  lip->data_in_log = TRUE;
}


static void ReportOneMailReportNoTaxId (MailReportPtr m, LogInfoPtr lip)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;

  if (m == NULL || m->bsp == NULL) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (m->bsp, NULL, Seq_descr_source, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (m->bsp, sdp, Seq_descr_source, &context)) {
    if ((biop = (BioSourcePtr) sdp->data.ptrvalue) != NULL 
        && biop == m->biop
        &&!HasTaxonomyID(biop)) {
      ReportOneMailReport (m, biop, lip);
      lip->data_in_log = TRUE;
      break;
    }
  }
}


static void ReportOneMailReportSpToNoSp (MailReportPtr m, LogInfoPtr lip)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;

  if (m == NULL || m->bsp == NULL || StringStr (m->old_name, " sp.") == NULL) {
    return;
  }

  for (sdp = SeqMgrGetNextDescriptor (m->bsp, NULL, Seq_descr_source, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (m->bsp, sdp, Seq_descr_source, &context)) {
    if ((biop = (BioSourcePtr) sdp->data.ptrvalue) != NULL 
        && biop == m->biop
        && biop->org != NULL && biop->org->taxname != NULL
        && StringStr (biop->org->taxname, " sp.") == NULL) {
      ReportOneMailReport (m, biop, lip);
      lip->data_in_log = TRUE;
    }
  }
}


static void ReportUnpublishedTaxnames (ValNodePtr list, LogInfoPtr lip)
{
  ValNodePtr request_list = NULL, response_list;
  ValNodePtr vnp_m, vnp_r;
  MailReportPtr m;

  for (vnp_m = list; vnp_m != NULL; vnp_m = vnp_m->next) {
    m = (MailReportPtr) vnp_m->data.ptrvalue;
    ValNodeAddPointer (&request_list, 3, m->biop->org);
  }
  response_list = Taxon3GetOrgRefList (request_list);
  ValNodeFree (request_list);

  fprintf (lip->fp, "\n\nUnpublished Names\n");
  fprintf (lip->fp, "Accession\tOld Name\tNew Name\n");
  for (vnp_m = list, vnp_r = response_list;
       vnp_m != NULL && vnp_r != NULL;
       vnp_m = vnp_m->next, vnp_r = vnp_r->next) {
    if (vnp_r->choice & eReturnedOrgFlag_unpublished) {
      m = (MailReportPtr) vnp_m->data.ptrvalue;
      ReportOneMailReport (m, m->biop, lip);
    }
  }

  ValNodeFree (response_list);
}


static void WantBarcodeAdjustmentCallback (BioseqPtr bsp, Pointer data)
{
  if (bsp != NULL && data != NULL && HasBARCODETech(bsp)) {
    *((BoolPtr)data) = TRUE;
  }
}


static Boolean WantBarcodeAdjustment (SeqEntryPtr sep)
{
  Boolean want_adjustment = FALSE;
  Char    path [PATH_MAX];
  FILE   *fp;  

  if (!GetSequinAppParam ("BARCODE", "LIVELIST", NULL, path, sizeof (path))) {
    return FALSE;
  }
  if (StringHasNoText (path)) {
    return FALSE;
  }

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    return FALSE;
  }
  FileClose (fp);

  VisitBioseqsInSep (sep, &want_adjustment, WantBarcodeAdjustmentCallback);

  return want_adjustment;
}


static TwoStringHashPtr sBarcodeLiveListAccBin = NULL;
static TwoStringHashPtr sBarcodeLiveListAccId = NULL;
static TwoStringHashPtr sBarcodeLiveListIdBin = NULL;
static Boolean sAttemptedToReadBarcodeLiveList = FALSE;

static void FixBarcodeLiveListIds (ValNodePtr table)
{
  ValNodePtr row, col;
  CharPtr    cp;

  for (row = table; row != NULL; row = row->next) {
    col = row->data.ptrvalue;
    if (col != NULL && col->next != NULL && (cp = StringChr (col->next->data.ptrvalue, '.')) != NULL) {
      *cp = 0;
    }
  }
}


static Boolean ReadBarcodeLiveList (void)
{
  ValNodePtr table, col, vnp;
  Char    path [PATH_MAX];
  FILE   *fp; 
  Int4   len, i;

  if (sAttemptedToReadBarcodeLiveList) {
    if (sBarcodeLiveListAccBin == NULL) {
      return FALSE;
    } else {
      return TRUE;
    }
  }
  if (!GetSequinAppParam ("BARCODE", "LIVELIST", NULL, path, sizeof (path))) {
    return FALSE;
  }
  if (StringHasNoText (path)) {
    return FALSE;
  }

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    return FALSE;
  }
  table = ReadTabTableFromFile (fp);
  FileClose (fp);

  if (table == NULL) {
    return FALSE;
  }

  /* process IDs in second column are actually just the part before the first period */
  FixBarcodeLiveListIds (table);
  table = SortTableRowByAnyColumn (table, 1);
  len = ValNodeLen (table);
  sBarcodeLiveListAccBin = (TwoStringHashPtr) MemNew (sizeof (TwoStringHashData));
  sBarcodeLiveListAccBin->table = (CharPtr PNTR) MemNew (sizeof (CharPtr) * len * 2);
  sBarcodeLiveListAccBin->num_lines = 0;
  sBarcodeLiveListAccId = (TwoStringHashPtr) MemNew (sizeof (TwoStringHashData));
  sBarcodeLiveListAccId->table = (CharPtr PNTR) MemNew (sizeof (CharPtr) * len * 2);
  sBarcodeLiveListAccId->num_lines = 0;
  for (vnp = table; vnp != NULL; vnp = vnp->next) {
    col = vnp->data.ptrvalue;
    if (col != NULL && !StringHasNoText (col->data.ptrvalue)) {
      if (col->next != NULL && !StringHasNoText (col->next->data.ptrvalue)) {
        i = sBarcodeLiveListAccId->num_lines;
        sBarcodeLiveListAccId->table[2 * i] = StringSave (col->data.ptrvalue);
        sBarcodeLiveListAccId->table[2 * i + 1] = StringSave (col->next->data.ptrvalue);
        sBarcodeLiveListAccId->num_lines++;
      }
      if (col->next != NULL && col->next->next != NULL && col->next->next->next != NULL && !StringHasNoText (col->next->next->next->data.ptrvalue)) {
        i = sBarcodeLiveListAccBin->num_lines;
        sBarcodeLiveListAccBin->table[2 * i] = StringSave (col->data.ptrvalue);
        sBarcodeLiveListAccBin->table[2 * i + 1] = StringSave (col->next->next->next->data.ptrvalue);
        sBarcodeLiveListAccBin->num_lines++;
      }
    }
  }
  table = SortTableRowByAnyColumn (table, 2);
  sBarcodeLiveListIdBin = (TwoStringHashPtr) MemNew (sizeof (TwoStringHashData));
  sBarcodeLiveListIdBin->table = (CharPtr PNTR) MemNew (sizeof (CharPtr) * len * 2);
  sBarcodeLiveListIdBin->num_lines = 0;
  for (vnp = table; vnp != NULL; vnp = vnp->next) {
    col = vnp->data.ptrvalue;
    if (col != NULL && col->next != NULL && !StringHasNoText (col->next->data.ptrvalue)
        && col->next->next != NULL && col->next->next->next != NULL && !StringHasNoText (col->next->next->next->data.ptrvalue)) {
      i = sBarcodeLiveListIdBin->num_lines;
      sBarcodeLiveListIdBin->table[2 * i] = StringSave (col->next->data.ptrvalue);
      sBarcodeLiveListIdBin->table[2 * i + 1] = StringSave (col->next->next->next->data.ptrvalue);
      sBarcodeLiveListIdBin->num_lines++;
    }
  }

  table = FreeTabTable (table);

  sAttemptedToReadBarcodeLiveList = TRUE;
  return TRUE;
}


static const CharPtr sBarcodeSpecialWords[] = {
  " sp.", " aff.", " cf.", " group", " complex", " nr.", " near", NULL };


static Boolean IsJustOneWord (CharPtr organism)
{
  if (StringChr (organism, ' ') == NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean OrganismWantsBarcodeAdjustment (CharPtr organism)
{
  Boolean rval = FALSE;
  Int4    i;
  /* TODO - only if contains one of the special words, or is just one word */
  if (IsJustOneWord (organism)) {
    rval = TRUE;
  } else {
    for (i = 0; sBarcodeSpecialWords[i] != NULL && !rval; i++) {
      if (StringSearch (organism, sBarcodeSpecialWords[i]) != NULL) {
        rval = TRUE;
      }
    }
  }

  return rval;
}


static Boolean FieldTypeIsTaxname (ValNodePtr field)
{
  ValNodePtr srcq;
  if (field != NULL
      && field->choice == FieldType_source_qual
      && (srcq = field->data.ptrvalue) != NULL
      && srcq->choice == SourceQualChoice_textqual
      && srcq->data.intvalue == Source_qual_taxname) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static ValNodePtr AdjustTableForBarcodeLiveList (ValNodePtr table, QualLoadFormPtr form_data)
{
  ValNodePtr new_table;
  ValNodePtr         columns = NULL, vnp, org_vnp, acc_vnp;
  TabColumnConfigPtr t;
  Int4               accession_column = -1;
  Int4               organism_column = -1;
  Int4               i;
  CharPtr            uniquer, id, tmp_name;

  new_table = CopyTabTable (table);

  if (form_data->list_dlg == NULL) {
    for (i = 0; i < form_data->num_columns; i++) {
      t = DialogToPointer (form_data->column_list[i]);
      if (t != NULL && t->match_type != NULL && t->match_type->choice == eTableMatchNucID) {
        accession_column = i + 1;
      } else if (t != NULL && FieldTypeIsTaxname(t->field)) {
        organism_column = i + 1;
      }
      t = TabColumnConfigFree (t);
    }
  } else {
    columns = DialogToPointer (form_data->list_dlg);
    for (vnp = columns, i = 0; vnp != NULL; vnp = vnp->next, i++) {
      t = (TabColumnConfigPtr) vnp->data.ptrvalue;
      if (t != NULL && t->match_type != NULL && t->match_type->choice == eTableMatchNucID) {
        accession_column = i;
      } else if (t != NULL && FieldTypeIsTaxname(t->field)) {
        organism_column = i;
      }
    }
    columns = TabColumnConfigListFree (columns);
  }

  if (accession_column > -1 && organism_column > -1) {
    /* do replacements */
    for (vnp = new_table; vnp != NULL; vnp = vnp->next) {
      org_vnp = GetNthValNode(vnp->data.ptrvalue, organism_column);
      if (org_vnp != NULL && OrganismWantsBarcodeAdjustment (org_vnp->data.ptrvalue)) {
        acc_vnp = GetNthValNode(vnp->data.ptrvalue, accession_column);
        if (acc_vnp != NULL) {
          uniquer = GetValueFromTwoStringHash(acc_vnp->data.ptrvalue, sBarcodeLiveListAccBin);
          if (uniquer == NULL) {
            id = GetValueFromTwoStringHash(acc_vnp->data.ptrvalue, sBarcodeLiveListAccId);
            if (id != NULL) {
              uniquer = GetValueFromTwoStringHash(acc_vnp->data.ptrvalue, sBarcodeLiveListIdBin);
              if (uniquer == NULL) {
                uniquer = id;
              }
            }
          }
          if (uniquer != NULL) {
            /* make sure it's not already in the organism name */
            tmp_name = org_vnp->data.ptrvalue;
            if (StringSearch (tmp_name, uniquer) == NULL) {
              if (IsJustOneWord(tmp_name)) {
                SetStringValue (&tmp_name, "sp.", ExistingTextOption_append_space);
              }
              SetStringValue (&tmp_name, uniquer, ExistingTextOption_append_space);
              org_vnp->data.ptrvalue = tmp_name;
            }              
          }
        }
      }
    }
  }

  return new_table;
}


static void SpecialMailReport (ButtoN b)
{
  QualLoadFormPtr form_data;
  ValNodePtr      adjusted_table = NULL, tmp_table = NULL;
  ValNodePtr      list = NULL, vnp;
  SeqEntryPtr     sep;
  LogInfoPtr      lip;
  MailReportPtr   m;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;

  /* first, get list of accessions and current names */
  sep = GetTopSeqEntryForEntityID (form_data->input_entityID);
  VisitBioseqsInSep (sep, &list, MailReportCallback);

  /* adjust table if necessary for barcode livelist */
  tmp_table = form_data->table;
  if (WantBarcodeAdjustment (sep) && ReadBarcodeLiveList()) {
    adjusted_table = AdjustTableForBarcodeLiveList (form_data->table, form_data);
    form_data->table = adjusted_table;
  } 

  /* now apply table */
  DoAcceptQuals (b);

  /* put original table back */
  form_data->table = tmp_table;
  adjusted_table = FreeTabTable (adjusted_table);

  /* now do taxfix */
  Taxon3ReplaceOrgInSeqEntryEx (sep, FALSE, FALSE);
  ObjMgrSetDirtyFlag (form_data->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, form_data->input_entityID, 0, 0);
  Update();

  /* now create report for organisms without IDs */
  lip = OpenLog ("Organism Change Report");
  fprintf (lip->fp, "Failed Lookups\n");
  fprintf (lip->fp, "Accession\tOld Name\tNew Name\n");
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    m = (MailReportPtr) vnp->data.ptrvalue;
    ReportOneMailReportNoTaxId(m, lip);
  }
  fprintf (lip->fp, "\n\nSp. Replaced with Real\n");

  fprintf (lip->fp, "Accession\tOld Name\tNew Name\n");
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    m = (MailReportPtr) vnp->data.ptrvalue;
    ReportOneMailReportSpToNoSp(m, lip);
  }

  ReportUnpublishedTaxnames (list, lip);

  CloseLog (lip);
  lip = FreeLog (lip);
        
  MailReportListFree(list);

}

/* J. Chen */
static void ReMatchQuals(ButtoN b) 
{
  QualLoadFormPtr form_data;
  ValNodePtr      val, vnp;
  ValNodePtr      columns = NULL;
  Int4            i;
  TabColumnConfigPtr t;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;

  if (form_data->list_dlg == NULL) {
    for (i = 0; i < form_data->num_columns; i++) {
      t = DialogToPointer (form_data->column_list[i]);
      ValNodeAddPointer (&columns, 0, t);
    }
  } else {
    columns = DialogToPointer (form_data->list_dlg);
  }

  for (i = 0, val = form_data->table->data.ptrvalue, vnp = columns;
       val != NULL;
       val = val->next, vnp = vnp->next, i++) {
    if (vnp == NULL) {
      vnp = ValNodeNew (columns);
    }
    t = vnp->data.ptrvalue;
    if (t == NULL) {
      t = TabColumnConfigNew ();
      vnp->data.ptrvalue = t;
    }
    t->field = FieldTypeFromString (val->data.ptrvalue);
    if (t->field == NULL) {
	if (!strcmp(val->data.ptrvalue, "Accession")) { /* new MatchTypeData */
	    t->match_type = MatchTypeNew();   /* default: choice = eTableMatchNucID; */
	}
	else {
           t = TabColumnConfigFree (t);
           vnp->data.ptrvalue = NULL;
	}
	 
    }
    else  {  /* if match_type doesn't match the field, whatever, change "Match" to "Apply" */
	t->match_type = MatchTypeFree(t->match_type);	
    }
  }

  if (form_data->list_dlg == NULL) {
    for (i = 0, vnp = columns; i < form_data->num_columns && vnp != NULL; i++, vnp = vnp->next) {
      PointerToDialog (form_data->column_list[i], vnp->data.ptrvalue);
    }
  } else {
    PointerToDialog (form_data->list_dlg, columns);
  }
  columns = TabColumnConfigListFree (columns);
}  /* ReMatchQuals */


/* J. Chen */
static void EraseAllFieldsWhenBlank(ButtoN b)
{
  QualLoadFormPtr form_data;
  ValNodePtr      val, vnp;
  ValNodePtr      columns = NULL;
  Int4            i;
  TabColumnConfigPtr t;
  TabColumnConfigListDlgPtr dlg;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;


  if (form_data->list_dlg != NULL) {
     dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (form_data->list_dlg);
     for (i=0; i<  dlg->num_columns; i++) {
	if ( dlg->column_list[i] != NULL) {
	   if (dlg->blank_list[i] > 0) dlg->column_list[i]->skip_blank = FALSE;
	   else dlg->column_list[i]->skip_blank = TRUE;
	}
     }
     PopulateTabConfigListColumnListDoc(form_data->list_dlg);
  }

  return;

}  /* EraseAllFieldsWhenBlank */



static void AutoMatchQuals (ButtoN b)
{
  QualLoadFormPtr form_data;
  ValNodePtr      val, vnp;
  ValNodePtr      columns = NULL;
  Int4            i;
  TabColumnConfigPtr t;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;

  if (form_data->list_dlg == NULL) {
    for (i = 0; i < form_data->num_columns; i++) {
      t = DialogToPointer (form_data->column_list[i]);
      ValNodeAddPointer (&columns, 0, t);
    }
  } else {
    columns = DialogToPointer (form_data->list_dlg);
  }

  for (val = form_data->table->data.ptrvalue, vnp = columns;
       val != NULL;
       val = val->next, vnp = vnp->next) {
    if (vnp == NULL) {
      vnp = ValNodeNew (columns);
    }
    t = vnp->data.ptrvalue;
    if (t == NULL) {
      t = TabColumnConfigNew ();
      vnp->data.ptrvalue = t;
    }
    if (t->match_type == NULL && t->field == NULL) {
      t->field = FieldTypeFromString (val->data.ptrvalue);
      if (t->field == NULL) {
        t = TabColumnConfigFree (t);
        vnp->data.ptrvalue = NULL;
      }
    }
  }

  if (form_data->list_dlg == NULL) {
    for (i = 0, vnp = columns; i < form_data->num_columns && vnp != NULL; i++, vnp = vnp->next) {
      PointerToDialog (form_data->column_list[i], vnp->data.ptrvalue);
    }
  } else {
    PointerToDialog (form_data->list_dlg, columns);
  }
  columns = TabColumnConfigListFree (columns);
}


static void ResetColumns (ButtoN b)
{
  QualLoadFormPtr form_data;
  ValNodePtr      vnp;
  ValNodePtr      columns = NULL;
  Int4            i;
  TabColumnConfigPtr t;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;

  if (form_data->list_dlg == NULL) {
    for (i = 0; i < form_data->num_columns; i++) {
      t = DialogToPointer (form_data->column_list[i]);
      TabColumnConfigReset(t);
      PointerToDialog (form_data->column_list[i], t);
      t = TabColumnConfigFree(t);
    }
  } else {
    columns = DialogToPointer (form_data->list_dlg);
    for (vnp = columns;
         vnp != NULL;
         vnp = vnp->next) {
      vnp->data.ptrvalue = TabColumnConfigFree (vnp->data.ptrvalue);
    }
    PointerToDialog (form_data->list_dlg, columns);
    columns = TabColumnConfigListFree (columns);
  }  
}


static WindoW CreateTableReaderWindowWithTable (Uint2 entityID, ValNodePtr table);

static void LoadNewTable (ButtoN b)
{
  Char         path [PATH_MAX];
  FILE *fp;
  QualLoadFormPtr form_data, new_form_data;
  ValNodePtr      table, special_list, blank_list, columns = NULL, c_vnp;
  Int4            i;
  WindoW          w;
  TabColumnConfigPtr t;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL) return;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;
  
  fp = FileOpen (path, "r");

  table = ReadTabTableFromFile (fp);
  FileClose (fp);
  if (table == NULL) return;
  special_list = ScanTabTableForSpecialCharacters (table);
  if (special_list != NULL 
      && !FixSpecialCharactersForStringsInList (special_list, 
                                                "The table contains special characters\nand cannot be used until they are replaced.",
                                                FALSE)) {
    special_list = FreeContextList (special_list);
    return;
  }
  special_list = FreeContextList (special_list);
  if (form_data->list_dlg != NULL) {
    form_data->table = table;
    /* it's ok, can just repopulate dialogs */
    blank_list = CountTabTableBlanks (table);
    ChangeDataForTabColumnConfigListDialog (form_data->list_dlg, form_data->table->data.ptrvalue, blank_list);
    blank_list = ValNodeFree (blank_list);
  } else {
    /* build a new window */
    w = CreateTableReaderWindowWithTable (form_data->input_entityID, table);
    /* get old configurations */
    if (form_data->list_dlg == NULL) {
      for (i = 0; i < form_data->num_columns; i++) {
        t = DialogToPointer (form_data->column_list[i]);
        ValNodeAddPointer (&columns, 0, t);
      }
    } else {
      columns = DialogToPointer (form_data->list_dlg);
    }

    /* apply to new window */
    new_form_data = (QualLoadFormPtr) GetObjectExtra (w);
    if (new_form_data->list_dlg == NULL) {
      for (i = 0, c_vnp = columns; i < new_form_data->num_columns && c_vnp != NULL; i++, c_vnp = c_vnp->next) {
        PointerToDialog (new_form_data->column_list[i], c_vnp->data.ptrvalue);
      }
    } else {
      PointerToDialog (new_form_data->list_dlg, columns);
    }
    columns = ValNodeFree (columns);
    SetStatus (new_form_data->remove_quotes, GetStatus (form_data->remove_quotes));
    SetStatus (new_form_data->leave_dlg_up, GetStatus (form_data->leave_dlg_up));

    /* remove old window */
    Remove ((WindoW) form_data->form);
    
    RealizeWindow (w);
    Show (w);
    Update ();
  }


}

typedef void (*Nlm_ReparseTabTable) PROTO ((ValNodePtr));


static void ReDrawWindowAfterTableReload(QualLoadFormPtr form_data)
{
  QualLoadFormPtr new_form_data;
  ValNodePtr      blank_list, columns = NULL, c_vnp, table;
  Int4            i;
  WindoW          w;
  TabColumnConfigPtr t;

  if (form_data->list_dlg != NULL) {
    /* it's ok, can just repopulate dialogs */
    blank_list = CountTabTableBlanks (form_data->table);
    ChangeDataForTabColumnConfigListDialog (form_data->list_dlg, form_data->table->data.ptrvalue, blank_list);
    blank_list = ValNodeFree (blank_list);
  } else {
    /* build a new window */
    table = form_data->table;
    form_data->table = NULL;

    w = CreateTableReaderWindowWithTable (form_data->input_entityID, table);
    /* get old configurations */
    if (form_data->list_dlg == NULL) {
      for (i = 0; i < form_data->num_columns; i++) {
        t = DialogToPointer (form_data->column_list[i]);
        ValNodeAddPointer (&columns, 0, t);
      }
    } else {
      columns = DialogToPointer (form_data->list_dlg);
    }

    /* apply to new window */
    new_form_data = (QualLoadFormPtr) GetObjectExtra (w);
    if (new_form_data->list_dlg == NULL) {
      for (i = 0, c_vnp = columns; i < new_form_data->num_columns && c_vnp != NULL; i++, c_vnp = c_vnp->next) {
        PointerToDialog (new_form_data->column_list[i], c_vnp->data.ptrvalue);
      }
    } else {
      PointerToDialog (new_form_data->list_dlg, columns);
    }
    columns = ValNodeFree (columns);
    SetStatus (new_form_data->remove_quotes, GetStatus (form_data->remove_quotes));
    SetStatus (new_form_data->leave_dlg_up, GetStatus (form_data->leave_dlg_up));

    /* remove old window */
    Remove ((WindoW) form_data->form);
    
    RealizeWindow (w);
    Show (w);
    Update ();
  }
}


static void ReLoadTable (ButtoN b, Nlm_ReparseTabTable reparse_func)
{
  QualLoadFormPtr form_data;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL || reparse_func == NULL) return;

  reparse_func (form_data->table);

  ReDrawWindowAfterTableReload (form_data);
}


static void ReloadFirstSpace (ButtoN b)
{
  ReLoadTable (b, ReparseTabTableConvertFirstSpaceToTab);
}


static void ReloadMultiSpace (ButtoN b)
{
  ReLoadTable (b, ReparseTabTableConvertMultiSpaceToTab);
}


static void ReparseForFileId (ValNodePtr table)
{
  ValNodePtr column_list = NULL;
  ValNodeAddInt (&column_list, 0, 0);
  ValNodeAddInt (&column_list, 0, 1);
  CombineTabTableColumns (table, column_list, "/");
}


static void ReloadFileId (ButtoN b)
{
  ReLoadTable (b, ReparseForFileId);
}


static int LIBCALLBACK SortVnpByInt (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  if (vnp1->data.intvalue > vnp2->data.intvalue) {
    return 1;
  } else if (vnp1->data.intvalue < vnp2->data.intvalue) {
    return -1;
  }

  return 0;
}


static void CombineColumns (ButtoN b_in)
{
  WindoW w;
  GrouP  h, g, c;
  DialoG dlg;
  TexT   delimiter_txt;
  ValNodePtr options = NULL, vnp, combine_list = NULL;
  Int4       i;
  ButtoN b;
  ModalAcceptCancelData acd;
  QualLoadFormPtr form_data;
  CharPtr         label, delimiter;

  form_data = (QualLoadFormPtr) GetObjectExtra (b_in);
  if (form_data == NULL || form_data->table == NULL || form_data->table->data.ptrvalue == NULL) return;
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  acd.third_option = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  for (i = 0, vnp = form_data->table->data.ptrvalue;
       i < form_data->num_columns && vnp != NULL; 
       i++, vnp = vnp->next) 
  {
    label = (CharPtr) MemNew (sizeof (Char) * (StringLen (vnp->data.ptrvalue) + 15));
    sprintf (label, "%5d   %s", i + 1, (char *) vnp->data.ptrvalue);
    ValNodeAddPointer (&options, i, label);
  }

  dlg = ValNodeSelectionDialog (h, options, TALL_SELECTION_LIST, ValNodeStringName,
                                              ValNodeSimpleDataFree, ValNodeStringCopy,
                                              ValNodeStringMatch, "column", 
                                              NULL, NULL, TRUE);
  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Delimiter", 0, popupMenuHeight, programFont, 'l');
  delimiter_txt = DialogText (g, ",", 5, NULL);

  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Combine", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg, (HANDLE) g, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
    ProcessAnEvent ();
    if (acd.accepted)
    {
      /* combine here */
      options = DialogToPointer (dlg);
      if (options == NULL || options->next == NULL) {
        Message (MSG_ERROR, "Must choose at least two columns!");
        acd.accepted = FALSE;
      } else {
        for (vnp = options; vnp != NULL; vnp = vnp->next) {
          ValNodeAddInt (&combine_list, 0, vnp->choice);
        }
        combine_list = ValNodeSort (combine_list, SortVnpByInt);
        delimiter = JustSaveStringFromText (delimiter_txt);
        CombineTabTableColumns (form_data->table, combine_list, delimiter);
        delimiter = MemFree (delimiter);
        combine_list = ValNodeFree (combine_list);
        ReDrawWindowAfterTableReload (form_data);
      }
      options = ValNodeFree (options);
    }
  }
  Remove (w);
}


static void SeparateColumns (ButtoN b_in)
{
  WindoW w;
  GrouP  h, g, c;
  DialoG dlg;
  TexT   delimiter_txt;
  ButtoN just_first_delimiter_btn;
  ValNodePtr options = NULL, vnp;
  Int4       i;
  ButtoN b;
  ModalAcceptCancelData acd;
  QualLoadFormPtr form_data;
  CharPtr         label, delimiter;
  Boolean         just_first_delimiter;

  form_data = (QualLoadFormPtr) GetObjectExtra (b_in);
  if (form_data == NULL || form_data->table == NULL || form_data->table->data.ptrvalue == NULL) return;
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  acd.third_option = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  for (i = 0, vnp = form_data->table->data.ptrvalue;
       i < form_data->num_columns && vnp != NULL; 
       i++, vnp = vnp->next) 
  {
    label = (CharPtr) MemNew (sizeof (Char) * (StringLen (vnp->data.ptrvalue) + 15));
    sprintf (label, "%5d   %s", i + 1, (char *) vnp->data.ptrvalue);
    ValNodeAddPointer (&options, i, label);
  }

  dlg = ValNodeSelectionDialog (h, options, TALL_SELECTION_LIST, ValNodeStringName,
                                              ValNodeSimpleDataFree, ValNodeStringCopy,
                                              ValNodeStringMatch, "column", 
                                              NULL, NULL, TRUE);
  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Delimiter", 0, popupMenuHeight, programFont, 'l');
  delimiter_txt = DialogText (g, ",", 5, NULL);

  just_first_delimiter_btn = CheckBox (h, "Separate only at first delimiter", NULL);

  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Separate", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg, (HANDLE) g, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
    ProcessAnEvent ();
    if (acd.accepted)
    {
      /* combine here */
      options = DialogToPointer (dlg);
      if (options == NULL) {
        Message (MSG_ERROR, "Must choose at least one column!");
        acd.accepted = FALSE;
      } else {
        delimiter = JustSaveStringFromText (delimiter_txt);
        if (delimiter == NULL) {
          delimiter = StringSave (" ");
        }
        just_first_delimiter = GetStatus (just_first_delimiter_btn);
        for (vnp = options; vnp != NULL; vnp = vnp->next) {
          ReparseTabTableSeparateColumnAtDelimiter (form_data->table, *delimiter, vnp->choice, just_first_delimiter);
        }
        delimiter = MemFree (delimiter);
        ReDrawWindowAfterTableReload (form_data);
      }
      options = ValNodeFree (options);
    }
  }
  Remove (w);
}


static void ShowTableReaderTable (ButtoN b)
{
  QualLoadFormPtr form_data;
  Char            path [PATH_MAX];
  FILE           *fp;

  form_data = (QualLoadFormPtr) GetObjectExtra (b);
  if (form_data == NULL || form_data->table == NULL) return;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open temporary file %s", path);
    return;
  }
  WriteTabTableToFile (form_data->table, fp);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "WARNING - Editing in this window will not be reflected in the table reader!");
  FileRemove (path);
}


static WindoW CreateTableReaderWindowWithTable (Uint2 entityID, ValNodePtr table)
{
  ValNodePtr   blank_list;
  ValNodePtr   col_vnp, blank_vnp;
  Int4         index;
  WindoW        w;
  GrouP         h, g, g2, g3, g4, c;
  QualLoadFormPtr form_data;
  CharPtr         title;
  Int4            col_for_list = 3;
  Char            buf[15];
  ButtoN          b;
  
  form_data = MemNew (sizeof (QualLoadFormData));
  if (form_data == NULL) return NULL;
  form_data->input_entityID = entityID;
  form_data->table = table;
  blank_list = CountTabTableBlanks (form_data->table);
  form_data->num_columns = ValNodeLen (blank_list);
  form_data->column_list = (DialoG PNTR) MemNew (sizeof (DialoG) * form_data->num_columns);

  /* now create a dialog to display values */
  w = FixedWindow (-50, -33, -10, -10, "Apply Qualifiers", StdCloseWindowProc);
  SetObjectExtra (w, form_data, CleanupQualLoadForm);
  form_data->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  if (GetSequinAppParam ("SETTINGS", "tablecolumns", NULL, buf, sizeof (buf))) {
    col_for_list = atoi (buf);
  }


  if (form_data->num_columns >= col_for_list) {
    form_data->list_dlg = TabColumnConfigListDialog (h, form_data->table->data.ptrvalue, blank_list, ChangeTabColumnChoice, form_data);
    g = (GrouP) form_data->list_dlg;
  } else {
    g = HiddenGroup (h, 3, 0, NULL);
    col_vnp = form_data->table->data.ptrvalue;
    blank_vnp = blank_list;
    for (index = 0; index < form_data->num_columns; index++) {
      if (col_vnp == NULL || StringHasNoText (col_vnp->data.ptrvalue)) {
        title = StringSave ("First row value is blank");
      } else {
        title = StringSave (col_vnp->data.ptrvalue);
        if (StringLen (title) > 104) {
          StringCpy (title + 100, "...");
          *(title + 103) = 0;
        }
      }
      form_data->column_list[index] = TabColumnConfigDialog (g, title, blank_vnp->data.intvalue, ChangeTabColumnChoice, form_data);
      title = MemFree (title);
      if (col_vnp != NULL) {
        col_vnp = col_vnp->next;
      }
      blank_vnp = blank_vnp->next;
    }
  }

  g2 = HiddenGroup (h, 5, 0, NULL);
  b = PushButton (g2, "Automatch Qualifiers", AutoMatchQuals);
  SetObjectExtra (b, form_data, NULL);

  /* J. Chen */
  b = PushButton (g2, "Rematch Qualifiers", ReMatchQuals);
  SetObjectExtra (b, form_data, NULL);

  b = PushButton (g2, "Reset Columns", ResetColumns);
  SetObjectExtra (b, form_data, NULL);

  /* J. Chen */
  b = PushButton (g2, "Erase field when table cell is blank", EraseAllFieldsWhenBlank);
  SetObjectExtra(b, form_data, NULL); 

  form_data->remove_quotes = CheckBox (g2, "Remove quotes around values", NULL);
  SetStatus (form_data->remove_quotes, TRUE); /* TRUE for checked, FALSE for unchecked, J.C. */
  g4 = HiddenGroup (h, 8, 0, NULL);
  SetGroupSpacing (g4, 10, 10);
  b = PushButton (g4, "View Table (NO EDITING!)", ShowTableReaderTable);
  SetObjectExtra (b, form_data, NULL);
  g3 = HiddenGroup (h, 8, 0, NULL);
  b = PushButton (g3, "Load New Table", LoadNewTable);
  SetObjectExtra (b, form_data, NULL);
  b = PushButton (g3, "Reload, 1st Space=Tab", ReloadFirstSpace);
  SetObjectExtra (b, form_data, NULL);
  b = PushButton (g3, "Reload, Mult Space=Tab", ReloadMultiSpace);
  SetObjectExtra (b, form_data, NULL);
  b = PushButton (g3, "Reload, Combine First Two for ID", ReloadFileId);
  SetObjectExtra (b, form_data, NULL);
  b = PushButton (g3, "Combine Columns", CombineColumns);
  SetObjectExtra (b, form_data, NULL);
  b = PushButton (g3, "Separate Columns", SeparateColumns);
  SetObjectExtra (b, form_data, NULL);

  form_data->err_msg = StaticPrompt (h, "", 200, dialogTextHeight, systemFont, 'c');
  c = HiddenGroup (h, 7, 0, NULL);
  form_data->accept_button = DefaultButton (c, "Accept", DoAcceptQuals);
  SetObjectExtra (form_data->accept_button, form_data, NULL);
  Disable (form_data->accept_button);
  form_data->mail_button = DefaultButton (c, "Accept and Make Mail Report", SpecialMailReport);
  SetObjectExtra (form_data->mail_button, form_data, NULL);
  Disable (form_data->mail_button);
  
  PushButton (c, "Cancel", StdCancelButtonProc);
  form_data->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  AlignObjects (ALIGN_CENTER,
                (HANDLE) g, 
                (HANDLE) g2,
                (HANDLE) g3,
                (HANDLE) g4,
                (HANDLE) form_data->err_msg,
                (HANDLE) c, NULL);

  blank_list = ValNodeFree (blank_list);
  return w;
}


static WindoW CreateTableReaderWindow (Uint2 entityID)
{
  Char         path [PATH_MAX];
  FILE *fp;

  ValNodePtr   table;
  ValNodePtr      special_list;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return NULL;
  
  fp = FileOpen (path, "r");

  table = ReadTabTableFromFile (fp);
  FileClose (fp);
  if (table == NULL) return NULL;
  special_list = ScanTabTableForSpecialCharacters (table);
  if (special_list != NULL 
      && !FixSpecialCharactersForStringsInList (special_list, 
                                                "The table contains special characters\nand cannot be used until they are replaced.",
                                                FALSE)) {
    special_list = FreeContextList (special_list);
    return NULL;
  }
  special_list = FreeContextList (special_list);

  return CreateTableReaderWindowWithTable (entityID, table);
}


static void NewLoadFeatureQualifierTableEntityID (Uint2 entityID)
{
  WindoW        w;
  QualLoadFormPtr form_data;
  TabColumnConfigPtr t;
  ValNodePtr         column_list = NULL;


  w = CreateTableReaderWindow (entityID);

  if (w != NULL) {
    /* populate */
    form_data = (QualLoadFormPtr) GetObjectExtra (w);
    if (form_data != NULL) {
      t = TabColumnConfigNew ();
      t->match_type = MatchTypeNew ();
      t->match_type->match_location = eTableMatchNucID;
      if (form_data->list_dlg == NULL) {
        PointerToDialog (form_data->column_list[0], t);
      } else {
        ValNodeAddPointer (&column_list, 0, t);
        PointerToDialog (form_data->list_dlg, column_list);
        column_list = ValNodeFree (column_list);
      }
      t = TabColumnConfigFree (t);
    } 

    RealizeWindow (w);
    Show (w);
    Update ();
  }
}


extern void NewLoadFeatureQualifierTable (IteM i)
{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  NewLoadFeatureQualifierTableEntityID (bfp->input_entityID);
}


extern void NewLoadFeatureQualifierTableToolBtn (ButtoN b)
{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  NewLoadFeatureQualifierTableEntityID (bfp->input_entityID);
}


NLM_EXTERN void CreateTableReaderWindowWithStructuredComments (IteM i)
{
  BaseFormPtr  bfp;
  WindoW        w;
  QualLoadFormPtr form_data;
  TabColumnConfigPtr t;
  ValNodePtr         table, column_list = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  table = CreateStructuredCommentTableFromSeqEntry(GetTopSeqEntryForEntityID(bfp->input_entityID));
  w = CreateTableReaderWindowWithTable (bfp->input_entityID, table);

  if (w != NULL) {
    /* populate */
    form_data = (QualLoadFormPtr) GetObjectExtra (w);
    if (form_data != NULL) {
      t = TabColumnConfigNew ();
      t->match_type = MatchTypeNew ();
      t->match_type->match_location = eTableMatchNucID;
      if (form_data->list_dlg == NULL) {
        PointerToDialog (form_data->column_list[0], t);
      } else {
        ValNodeAddPointer (&column_list, 0, t);
        PointerToDialog (form_data->list_dlg, column_list);
        column_list = ValNodeFree (column_list);
      }
      t = TabColumnConfigFree (t);
    } 

    RealizeWindow (w);
    Show (w);
    Update ();
  }
}  
  

static void LoadTaxTableReaderBaseForm (BaseFormPtr bfp)
{
  WindoW        w;
  QualLoadFormPtr form_data;
  TabColumnConfigPtr t;
  ValNodePtr         column_list = NULL, s;

  if (bfp == NULL) return;

  w = CreateTableReaderWindow (bfp->input_entityID);

  if (w != NULL) {
    /* populate */
    form_data = (QualLoadFormPtr) GetObjectExtra (w);
    if (form_data != NULL) {
      if (form_data->list_dlg == NULL) {
        if (form_data->num_columns > 1) {
          t = TabColumnConfigNew ();
          t->match_type = MatchTypeNew ();
          t->match_type->choice = eTableMatchBioSource;
          PointerToDialog (form_data->column_list[0], t);
          t->match_type = MatchTypeFree(t->match_type);
          t->field = ValNodeNew(NULL);
          t->field->choice = FieldType_source_qual;
          s = ValNodeNew (NULL);
          t->field->data.ptrvalue = s;
          s->choice = SourceQualChoice_textqual;
          s->data.intvalue = Source_qual_taxname;
          PointerToDialog (form_data->column_list[1], t);
          t = TabColumnConfigFree(t);
        }
      } else {
        t = TabColumnConfigNew ();
        t->match_type = MatchTypeNew ();
        t->match_type->choice = eTableMatchBioSource;
        PointerToDialog (form_data->column_list[0], t);
        ValNodeAddPointer (&column_list, 0, t);
        t = TabColumnConfigNew ();
        t->field = ValNodeNew(NULL);
        t->field->choice = FieldType_source_qual;
        s = ValNodeNew (NULL);
        t->field->data.ptrvalue = s;
        s->choice = SourceQualChoice_textqual;
        s->data.intvalue = Source_qual_taxname;
        ValNodeAddPointer (&column_list, 0, t);
        PointerToDialog (form_data->list_dlg, column_list);
        column_list = TabColumnConfigListFree (column_list);
      }
    } 

    RealizeWindow (w);
    Show (w);
    Update ();
  }
}


extern void LoadTaxTableReader (IteM i)
{
  ButtonOrMenuItemTemplate (i, LoadTaxTableReaderBaseForm);
}


extern void LoadTaxTableReaderBtn (ButtoN b)
{
  ButtonOrMenuButtonTemplate (b, LoadTaxTableReaderBaseForm);
}


extern void NewLoadSourceQualifierTable (IteM i)
{
  BaseFormPtr  bfp;
  WindoW        w;
  QualLoadFormPtr form_data;
  TabColumnConfigPtr t1, t2;
  ValNodePtr         column_list = NULL, vnp;
  Int4               n;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  w = CreateTableReaderWindow (bfp->input_entityID);

  if (w != NULL) {
    /* populate */
    form_data = (QualLoadFormPtr) GetObjectExtra (w);
    if (form_data != NULL) {
      t1 = TabColumnConfigNew ();
      t1->match_type = MatchTypeNew ();
      t1->match_type->match_location = eTableMatchNucID;
      t2 = TabColumnConfigNew ();
      vnp = ValNodeNew (NULL);
      vnp->choice = SourceQualChoice_textqual;
      vnp->data.intvalue = Source_qual_acronym;
      t2->field = ValNodeNew (NULL);
      t2->field->choice = FieldType_source_qual;
      t2->field->data.ptrvalue = vnp;

      if (form_data->list_dlg == NULL) {
        PointerToDialog (form_data->column_list[0], t1);
        PointerToDialog (form_data->column_list[1], t2);
      } else {
        ValNodeAddPointer (&column_list, 0, t1);
        for (n = 1; n < form_data->num_columns; n++) {
          ValNodeAddPointer (&column_list, 0, t2);
        }
        PointerToDialog (form_data->list_dlg, column_list);
        column_list = ValNodeFree (column_list);
      }
      t1 = TabColumnConfigFree (t1);
      t2 = TabColumnConfigFree (t2);
    } 

    RealizeWindow (w);
    Show (w);
    Update ();
  }
}


extern void ExternalSourceQualifierTableReader (IteM i)
{
  Char         path [PATH_MAX];
  BaseFormPtr  bfp;
  ValNodePtr   table, columns = NULL, header_row, val;
  TabColumnConfigPtr t;
  LogInfoPtr         lip;
  Boolean            any_valid = FALSE;
  MsgAnswer          ans = ANS_OK;
  SeqEntryPtr        sep;
  FILE *fp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;
  
  fp = FileOpen (path, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to read from %s", path);
    return;
  }

  table = ReadTabTableFromFile (fp);
  FileClose (fp);
  if (table == NULL) {
    Message (MSG_ERROR, "Unable to read table from %s", path);
    return;
  }

  RemoveQuotesFromTabTable (table);

  /* first column is sequence ID */
  t = TabColumnConfigNew ();
  t->match_type = MatchTypeNew ();
  t->match_type->match_location = eTableMatchNucID;
  ValNodeAddPointer (&columns, 0, t);

  lip = OpenLog ("Table Problems");
  fprintf (lip->fp, "The following header values could not be matched to valid qualifier names.  They will be ignored.\n");

  /* automatch qualifiers */
  header_row = table->data.ptrvalue;
  if (header_row == NULL) {
    Message (MSG_ERROR, "First row of table must contain headers!");
  } else {
    for (val = header_row->next;
         val != NULL;
         val = val->next) {
      t = TabColumnConfigNew ();
      t->field = FieldTypeFromString (val->data.ptrvalue);
      if (t->field == NULL) {
        t = TabColumnConfigFree (t);
        fprintf (lip->fp, "%s\n", (CharPtr) val->data.ptrvalue);
        lip->data_in_log = TRUE;
      } else {
        any_valid = TRUE;
      }
      ValNodeAddPointer (&columns, 0, t);
    }
  }

  CloseLog (lip);
  if (lip->data_in_log) {
    ans = Message (MSG_OKC, "Some column headers were not recognized.  The values in this column will be ignored.  Continue anyway?");
  }
  lip = FreeLog (lip);

  if (ans == ANS_OK) {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (ApplyTableValues (sep, table->next, columns)) {
      ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
      Update();
    }
  }
  columns = TabColumnConfigListFree (columns);
  table = FreeTabTable (table);
}


typedef struct prefixformdata {
  FEATURE_FORM_BLOCK

  ModifierItemLocalPtr modList;
  PopuP PNTR           popup_list;
  ButtoN               add_org_name;
  ButtoN               use_labels;
} PrefixFormData, PNTR PrefixFormPtr;

static void CleanupPrefixForm (
  GraphiC g,
  VoidPtr data
)

{
  PrefixFormPtr  pfp;
  Int4           i;

  pfp = (PrefixFormPtr) data;
  if (pfp != NULL) {
    if (pfp->modList != NULL)
    {
      for (i=0; i < NumDefLineModifiers (); i++)
      {
        ValNodeFree (pfp->modList[i].values_seen);
      }
      MemFree (pfp->modList);
    }
    if (pfp->popup_list != NULL)
    {
      MemFree (pfp->popup_list);
    }
  }
  StdCleanupFormProc (g, data);
}

static Int4 GetDefLineModifierPopupIndex (
  Int4 popup_index, 
  ModifierItemLocalPtr modList
)
{
  Int4 index;
  Int4 want_index;
  
  want_index = 0;
  for (index = 0;
       index < NumDefLineModifiers () && want_index < popup_index;
       index++)
  {
    if (modList [index].any_present)
    {
      want_index ++;
    }
  }
  if (index >= NumDefLineModifiers () || index == 0
    || want_index < popup_index)
  {
    return -1;
  }
  return index - 1;
}

static void AddPrefixesToOneDefLine (
  PrefixFormPtr pfp,
  BioseqPtr bsp,
  SeqEntryPtr sep
)
{
  BioSourcePtr biop;
  SeqDescrPtr  sdp;
  ValNodePtr   strings;
  Int4         index;
  Int4         qual_index;
  Int4         popup_choice;
  CharPtr      new_defline;
  Char         taxName [196];
  Char         modifier_text [256];
  OrgRefPtr    orp;
  OrgNamePtr   onp;
  OrgModPtr    mod;
  SubSourcePtr ssp;
  CharPtr  org_desc;
  CharPtr  cp;
  SeqMgrDescContext  dcontext;
  Boolean      use_labels;
  Uint4        no_semicolon_len;
  Int4         prefix_len;

  if (bsp == NULL || pfp == NULL
    || pfp->popup_list == NULL
    || pfp->modList == NULL
    || (biop = GetBiopForBsp (bsp)) == NULL)
  {
    return;
  }

  if (biop->org == NULL) return;
  if (biop->org->taxname == NULL) return;
  StringNCpy (taxName, biop->org->taxname, sizeof (taxName) - 1);
  taxName [ sizeof (taxName) - 1] = 0;

  strings = NULL;
  CleanUpTaxName (taxName, TRUE);
  if (GetStatus (pfp->add_org_name))
  {
    ValNodeCopyStr( &strings, 0, taxName);
  }

  use_labels = GetStatus (pfp->use_labels);

  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  for (index = 0; index < NumDefLineModifiers (); index ++)
  {
    if (pfp->popup_list [ index] != NULL
      && (popup_choice = GetValue (pfp->popup_list [ index])) > 0)
    {
      qual_index = GetDefLineModifierPopupIndex (popup_choice, pfp->modList);
      if (qual_index < 0 || qual_index >= NumDefLineModifiers ()) continue;
      if (DefLineModifiers[qual_index].isOrgMod)
      {
        mod = onp->mod;
        while (mod != NULL
          && mod->subtype != DefLineModifiers[qual_index].subtype)
        {
          mod = mod->next;
        }
        if ( UseOrgModifier (mod, taxName, FALSE))
        {
          no_semicolon_len = StringCSpn (mod->subname, ";");
          if (mod->subtype == ORGMOD_nat_host)
          {
            sprintf (modifier_text, "from ");
            if (no_semicolon_len > sizeof (modifier_text) - 6) {
              prefix_len = sizeof (modifier_text) - 1;
              no_semicolon_len = sizeof (modifier_text) - 6;
            } else {
	            prefix_len = no_semicolon_len + 5;
            }
            StringNCpy (modifier_text + 5, mod->subname, no_semicolon_len);
            modifier_text[prefix_len] = 0;
          }
          else
          {
            AddModifierLabel (use_labels, TRUE, mod->subtype, modifier_text);
            if (modifier_text[0] != 0)
              StringCat (modifier_text, " ");
            if (no_semicolon_len > sizeof (modifier_text) - StringLen (modifier_text) - 1) {
              no_semicolon_len = sizeof (modifier_text) - StringLen (modifier_text) - 1;
	            prefix_len = sizeof (modifier_text) - 1;
            } else {
	            prefix_len = no_semicolon_len + StringLen (modifier_text);
            }

            StringNCat (modifier_text, mod->subname, no_semicolon_len);
            modifier_text[prefix_len] = 0;
          }
          ValNodeCopyStr( &strings, 0, modifier_text);
        }
      } else {
        ssp = biop->subtype;
        while (ssp != NULL
            && ssp->subtype != DefLineModifiers[qual_index].subtype)
        {
          ssp = ssp->next;
        }
        if (ssp != NULL)
        {
          no_semicolon_len = StringCSpn (ssp->name, ";");
          AddModifierLabel (use_labels, FALSE, ssp->subtype, modifier_text);
          if (ssp->subtype == SUBSRC_country)
          {
            sprintf (modifier_text, "from ");
            if (no_semicolon_len > sizeof (modifier_text) - 6) {
              no_semicolon_len = sizeof (modifier_text) - 6;
			  prefix_len = sizeof(modifier_text);
			} else {
			  prefix_len = StringLen (modifier_text) + no_semicolon_len;
			}
            StringNCpy (modifier_text + 5, ssp->name, no_semicolon_len);
            modifier_text[prefix_len] = 0;
            cp = StringChr (modifier_text, ':');
            if (cp != NULL) *cp = 0;
          }
          else if (ssp->name != NULL
            && (ssp->subtype != SUBSRC_plasmid_name
              || StringCmp (ssp->name, "unnamed") != 0))
          {
            if (modifier_text[0] != 0)
              StringCat (modifier_text, " ");
            if (no_semicolon_len > sizeof (modifier_text) - StringLen (modifier_text) - 1) {
              no_semicolon_len = sizeof (modifier_text) - StringLen (modifier_text) - 1;
			  prefix_len = sizeof (modifier_text);
			} else {
			  prefix_len = StringLen (modifier_text) + no_semicolon_len;
			}
 
            StringNCat (modifier_text, ssp->name, no_semicolon_len);
            modifier_text[prefix_len] = 0;
          }
          ValNodeCopyStr( &strings, 0, modifier_text);
        }
      }
    }
  }
  org_desc = MergeValNodeStrings (strings, FALSE);
  ValNodeFreeData (strings);

  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  if (sdp == NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
    if (sdp == NULL) return;
    cp = (CharPtr) sdp->data.ptrvalue;
    if (cp == NULL) return;
    sdp = SeqDescrAdd (&(bsp->descr));
    if (sdp == NULL) return;
    sdp->choice = Seq_descr_title;
    sdp->data.ptrvalue = StringSave (cp);
  }
  if (sdp == NULL) return;
  if (StringHasNoText (sdp->data.ptrvalue)) return;
  new_defline = MemNew ((StringLen (org_desc)
                        + StringLen (sdp->data.ptrvalue) + 4) * sizeof (Char));
  if (new_defline == NULL) return;
  StringCpy (new_defline, org_desc);
  StringCat (new_defline, " ");
  StringCat (new_defline, sdp->data.ptrvalue);
  MemFree (sdp->data.ptrvalue);
  sdp->data.ptrvalue = new_defline;
  ObjMgrSetDirtyFlag (pfp->input_entityID, TRUE);
}

static void AddPrefixesToDefLines (PrefixFormPtr pfp, SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  if (pfp == NULL || sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        AddPrefixesToDefLines (pfp, sep);
      }
      return;
    }
  }
  if (! IS_Bioseq (sep)) return;
  AddPrefixesToOneDefLine (pfp, sep->data.ptrvalue, sep);
}

static void DoPrefixDefLines (ButtoN b)
{
  PrefixFormPtr pfp;
  SeqEntryPtr   sep;

  pfp = GetObjectExtra (b);
  if (pfp == NULL) return;

  Hide (pfp->form);
  WatchCursor ();
  Update ();
  sep = GetTopSeqEntryForEntityID (pfp->input_entityID);
  if (sep == NULL) return;
  AddPrefixesToDefLines (pfp, sep);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (pfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, pfp->input_entityID, 0, 0);
}

extern void PrefixDefLines (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  ModifierItemLocalPtr modList;
  Int4         index;
  WindoW       w;
  GrouP        h, g, k;
  Int4         popup_index, item_index, listed_index;
  GrouP        c;
  ButtoN       b;
  PrefixFormPtr pfp;
  Char         label [256];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;
  CountModifiers (modList, sep);

  pfp = MemNew (sizeof (PrefixFormData));
  if (pfp == NULL) return;
  pfp->input_entityID = bfp->input_entityID;
  pfp->modList = modList;
  pfp->popup_list = MemNew (sizeof (PopuP) * NumDefLineModifiers ());
  if (pfp->popup_list == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, "Definition Line Prefixes",
                   StdCloseWindowProc);
  SetObjectExtra (w, pfp, CleanupPrefixForm);
  pfp->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 4, 0, NULL);
  for (index = 0; index < NumDefLineModifiers (); index++)
  {
    pfp->popup_list [index] = NULL;
  }
  popup_index = 0;
  for (index = 0; index < NumDefLineModifiers (); index++)
  {
    if (modList [ index].any_present)
    {
      k = HiddenGroup (g, 2, 0, NULL);
      sprintf (label, "%d", popup_index + 1);
      StaticPrompt (k, label, 0, popupMenuHeight, programFont, 'l');
      pfp->popup_list[popup_index] = PopupList (k, TRUE, NULL);
      listed_index = 0;
      for (item_index = 0; item_index < NumDefLineModifiers (); item_index ++)
      {
        if (modList [item_index].any_present)
        {
          PopupItem (pfp->popup_list[popup_index], 
                   DefLineModifiers[item_index].name);
          listed_index ++;
        }
      }
      PopupItem (pfp->popup_list[popup_index], "Ignore");
      listed_index++;
      SetValue (pfp->popup_list[popup_index], listed_index);
      popup_index ++;
    }
  }

  pfp->add_org_name = CheckBox (h, "Prefix with taxonomy name", NULL);
  SetStatus (pfp->add_org_name, TRUE);
  pfp->use_labels = CheckBox (h, "Use modifier labels", NULL);
  SetStatus (pfp->use_labels, TRUE);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoPrefixDefLines);
  SetObjectExtra (b, pfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) pfp->add_org_name,
                (HANDLE) pfp->use_labels,
                (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();

}

#ifndef WIN16
CharPtr objPrtMemStr = "PrintTemplateSet ::= {\n" \
"{ name \"StdSeqDesc\" ,\n" \
"format { asn1 \"Seqdesc\" , form block {\n" \
"components {\n" \
"{ asn1 \"Seqdesc.mol-type\" , label \"Molecule type\" , prefix \"\\n\" , form enum { } } ,\n" \
"{ asn1 \"Seqdesc.modif\" , label \"Modifiers\" , prefix \"\\n\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Seqdesc.modif.E\" , form enum { } } } } } ,\n" \
"{ asn1 \"Seqdesc.method\" , label \"Method\" , prefix \"\\n\" , form enum { } } ,\n" \
"{ asn1 \"Seqdesc.name\" , label \"Name\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.title\" , label \"Title\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.org\" , label \"Organism\" , prefix \"\\n\" , form use-template \"StdOrgRef\" } ,\n" \
"{ asn1 \"Seqdesc.comment\" , label \"Comment\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.num\" , label \"Numbering\" , prefix \"\\n\" , form use-template \"StdNumbering\" } ,\n" \
"{ asn1 \"Seqdesc.maploc\" , label \"Map location\" , prefix \"\\n\" , form use-template \"StdDbtag\" } ,\n" \
"{ asn1 \"Seqdesc.pir\" , label \"PIR block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.genbank\" , label \"GenBank block\" , prefix \"\\n\" , form use-template \"StdGBBlock\" } ,\n" \
"{ asn1 \"Seqdesc.pub\" , label \"Citation\" , prefix \"\\n\" , form use-template \"StdPubdesc\" } ,\n" \
"{ asn1 \"Seqdesc.region\" , label \"Region\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.user\" , label \"User Type\" , prefix \"\\n\" , form use-template \"StdUserObj\" } ,\n" \
"{ asn1 \"Seqdesc.sp\" , label \"SWISS-PROT block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.dbxref\" , label \"Cross reference\" , prefix \"\\n\" , form use-template \"StdDbtag\"  } ,\n" \
"{ asn1 \"Seqdesc.embl\" , label \"EMBL block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.create-date\" , label \"Create date\" , prefix \"\\n\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Seqdesc.update-date\" , label \"Update date\" , prefix \"\\n\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Seqdesc.prf\" , label \"PRF block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.pdb\" , label \"PDB block\" , prefix \"\\n\" , form null NULL } ,\n" \
"{ asn1 \"Seqdesc.het\" , label \"Heterogen\" , prefix \"\\n\" , form text { } } ,\n" \
"{ asn1 \"Seqdesc.source\" , label \"Biological Source\" , prefix \"\\n\" , form use-template \"StdBioSource\" } ,\n" \
"{ asn1 \"Seqdesc.molinfo\" , label \"Molecule Information\" , prefix \"\\n\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"MolInfo.biomol\" , form enum { } } ,\n" \
"{ asn1 \"MolInfo.tech\" , form enum { } } ,\n" \
"{ asn1 \"MolInfo.completeness\" , form enum { } } } } } } } } } ,\n" \
"{ name \"StdSeqFeatLocation\" ,\n" \
"format { asn1 \"Seq-feat.location\" , label \"Location\" , prefix \"\\t\" , form user { printfunc \"StdSeqLocPrint\" } } } ,\n" \
"{ name \"StdSeqFeatProduct\" ,\n" \
"format { asn1 \"Seq-feat.product\" , label \"Product\" , prefix \"\\t\" , form user { printfunc \"StdSeqLocPrint\" } } } ,\n" \
"{ name \"EntrySeqFeatData\" ,\n" \
"labelfrom \"Seq-feat.data\" ,\n" \
"format { asn1 \"Seq-feat.data\" , prefix \"\\t\" , form use-template \"StdSeqFeatData\" } } ,\n" \
"{ name \"StdSeqFeat\" ,\n" \
"labelfrom \"Seq-feat.data\" ,\n" \
"format { asn1 \"Seq-feat\" , prefix \"\\n\" , suffix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.data\" , form use-template \"StdSeqFeatData\" } ,\n" \
"{ asn1 \"Seq-feat\" , form use-template \"StdSeqFeatCommon\" } ,\n" \
"{ asn1 \"Seq-feat.product\" , label \"Product\" , prefix \" \" , form user { printfunc \"StdSeqLocPrint\" } } ,\n" \
"{ asn1 \"Seq-feat.location\" , label \"Location\" , prefix \" \" , form user { printfunc \"StdSeqLocPrint\" } } ,\n" \
"{ asn1 \"Seq-feat.cit\" , label \"Citations\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.cit.pub.E\" , form use-template \"StdPub\" } } } } ,\n" \
"{ asn1 \"Seq-feat.xref\" , label \"Cross-reference\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.xref.E\" , form use-template \"StdSeqFeatXref\" } } } } } } } } ,\n" \
"{ name \"StdSeqFeatData\" ,\n" \
"format { asn1 \"SeqFeatData\" , form block {\n" \
"components {\n" \
"{ asn1 \"SeqFeatData.gene\" , label \"Gene\" , form use-template \"StdGeneRef\" } ,\n" \
"{ asn1 \"SeqFeatData.org\" , label \"Organism\" , form use-template \"StdOrgRef\" } ,\n" \
"{ asn1 \"SeqFeatData.cdregion\" , label \"Coding Region\" , form use-template \"StdCdRegion\" } ,\n" \
"{ asn1 \"SeqFeatData.prot\" , label \"Protein\" , form use-template \"StdProtRef\" } ,\n" \
"{ asn1 \"SeqFeatData.rna\" , label \"RNA\" , form use-template \"StdRNARef\" } ,\n" \
"{ asn1 \"SeqFeatData.pub\" , label \"Citation\" , form use-template \"StdPubdesc\" } ,\n" \
"{ asn1 \"SeqFeatData.seq\" , label \"Sequence\" , form user { printfunc \"StdSeqLocPrint\" } } ,\n" \
"{ asn1 \"SeqFeatData.imp.key\" , label \"Import\" , form use-template \"StdImpFeat\" } ,\n" \
"{ asn1 \"SeqFeatData.region\" , label \"Region\" , form text { } } ,\n" \
"{ asn1 \"SeqFeatData.comment\" , label \"Comment\" , form null NULL } ,\n" \
"{ asn1 \"SeqFeatData.bond\" , label \"Bond\" , form enum { } } ,\n" \
"{ asn1 \"SeqFeatData.site\" , label \"Site\" , form enum { } } ,\n" \
"{ asn1 \"SeqFeatData.rsite\" , label \"Rest. Site\" , form use-template \"StdRsiteRef\" } ,\n" \
"{ asn1 \"SeqFeatData.user\" , label \"User Type\" , form use-template \"StdUserObj\" } ,\n" \
"{ asn1 \"SeqFeatData.txinit\" , label \"TxInit\" , form use-template \"StdTxInit\" } ,\n" \
"{ asn1 \"SeqFeatData.num\" , label \"Numbering\" , form use-template \"StdNumbering\" } ,\n" \
"{ asn1 \"SeqFeatData.psec-str\" , label \"Sec. Struct\" , form enum { } } ,\n" \
"{ asn1 \"SeqFeatData.non-std-residue\" , label \"NonStd Residue\" , form text { } } ,\n" \
"{ asn1 \"SeqFeatData.het\" , label \"Heterogen\" , form text { } } ,\n" \
"{ asn1 \"SeqFeatData.biosrc\" , label \"Biological Source\" , prefix \"\\n\" , form use-template \"StdBioSource\" } } } } } ,\n" \
"{ name \"StdGeneRef\" ,\n" \
"format { asn1 \"Gene-ref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Gene-ref\" , form block {\n" \
"components {\n" \
"{ asn1 \"Gene-ref.locus\" , form text { } } ,\n" \
"{ asn1 \"Gene-ref.allele\" , prefix \" \" , form text { } } } } } ,\n" \
"{ asn1 \"Gene-ref.desc\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Gene-ref.pseudo\" , form boolean {\n" \
"true \"This is a pseudogene.\" } } ,\n" \
"{ asn1 \"Gene-ref.syn\" , label \"Synonyms\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Gene-ref.syn.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Gene-ref.maploc\" , label \"Map Location\" , prefix \" \" , form text { } } ,\n" \
"{ asn1 \"Gene-ref.db\" , label \"Cross Reference\" , prefix \" \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Gene-ref.db.E\" , prefix \"(\" , suffix \")\" , form use-template \"StdDbtag\" } } } } } } } } ,\n" \
"{ name \"StdUserObj\" ,\n" \
"format { asn1 \"User-object\" , label \"User-object\" , prefix \"\\n\" , form block {\n" \
"separator \": \" ,\n" \
"components {\n" \
"{ asn1 \"User-object.class\" , form text { } } ,\n" \
"{ asn1 \"User-object.type\" , form use-template \"StdObjectId\" } } } } } ,\n" \
"{ name \"StdPubOnFeat\" ,\n" \
"format { asn1 \"Pub\" , label \"Citation\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pub\" , form use-template \"StdPub\" } } } } } ,\n" \
"{ name \"StdPub\" ,\n" \
"format { asn1 \"Pub\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pub.gen\" , form use-template \"StdCitGen\" } ,\n" \
"{ asn1 \"Pub.sub\" , form use-template \"StdCitSub\" } ,\n" \
"{ asn1 \"Pub.medline\" , form use-template \"StdMedlineEntry\" } ,\n" \
"{ asn1 \"Pub.muid\" , label \"MEDLINE uid: \" , form text { } } ,\n" \
"{ asn1 \"Pub.pmid\" , label \"PubMed id: \" , form text { } } ,\n" \
"{ asn1 \"Pub.article\" , form use-template \"StdCitArt\" } ,\n" \
"{ asn1 \"Pub.journal\" , form use-template \"StdCitJour\" } ,\n" \
"{ asn1 \"Pub.book\" , form use-template \"StdCitBook\" } ,\n" \
"{ asn1 \"Pub.proc\" , form use-template \"StdCitProc\" } ,\n" \
"{ asn1 \"Pub.patent\" , form use-template \"StdCitPat\" } ,\n" \
"{ asn1 \"Pub.pat-id\" , form use-template \"StdIdPat\" } ,\n" \
"{ asn1 \"Pub.man\" , form use-template \"StdCitLet\" } ,\n" \
"{ asn1 \"Pub.equiv\" , form use-template \"StdPubEquiv\" } } } } } ,\n" \
"{ name \"StdCitGen\" ,\n" \
"format { asn1 \"Cit-gen\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-gen.serial-number\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-gen.date\" , prefix \"(\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Cit-gen.title\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen.cit\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Cit-gen.journal\" , suffix \":\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-gen.issue\" , suffix \";\" , form text { } } ,\n" \
"{ asn1 \"Cit-gen.pages\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdCitSub\" ,\n" \
"format { asn1 \"Cit-sub\" , prefix \"Data Submission \" , form block {\n" \
"components {\n" \
"{ asn1 \"Cit-sub.medium\" , prefix \"on \" , suffix \" \" , form enum { } } ,\n" \
"{ asn1 \"Cit-sub.imp.date\" , prefix \"(\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Cit-sub.authors\" , prefix \"\\n\" , form use-template \"StdAuthList\" } } } } } ,\n" \
"{ name \"StdMedlineEntry\" ,\n" \
"format { asn1 \"Medline-entry\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry\" , form block {\n" \
"separator \"   \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.uid\" , label \"uid\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Medline-entry.em\" , label \"entry month\" , prefix \": \" , form user { printfunc \"StdDatePrint\" } } } } } ,\n" \
"{ asn1 \"Medline-entry.cit\" , form use-template \"StdCitArt\" } ,\n" \
"{ asn1 \"Medline-entry.abstract\" , label \"abstract\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Medline-entry.mesh\" , label \"Mesh Terms\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.mesh.E\" , form block {\n" \
"components {\n" \
"{ asn1 \"Medline-mesh.term\" , form text { } } ,\n" \
"{ asn1 \"Medline-mesh.mp\" , form boolean {\n" \
"true \" (Main Point)\" } } ,\n" \
"{ asn1 \"Medline-mesh.qual\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-mesh.qual.E\" , form block {\n" \
"components {\n" \
"{ asn1 \"Medline-qual.subh\" , form text { } } ,\n" \
"{ asn1 \"Medline-qual.mp\" , form boolean {\n" \
"true \" (Main Point)\" } } } } } } } } } } } } } } ,\n" \
"{ asn1 \"Medline-entry.substance\" , label \"Substance\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.substance.E\" , form block {\n" \
"components {\n" \
"{ asn1 \"Medline-rn.name\" , form text { } } ,\n" \
"{ asn1 \"Medline-rn.type\" , form enum {\n" \
"values {\n" \
"\"\" ,\n" \
"\" CAS: \" ,\n" \
"\"EC \" } } } ,\n" \
"{ asn1 \"Medline-rn.cit\" , form text { } } } } } } } } ,\n" \
"{ asn1 \"Medline-entry.xref\" , label \"Cross Reference\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.xref.E\" , form block {\n" \
"separator \": \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-si.type\" , form enum { } } ,\n" \
"{ asn1 \"Medline-si.cit\" , form text { } } } } } } } } ,\n" \
"{ asn1 \"Medline-entry.gene\" , label \"Possible Gene Symbols\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.gene.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Medline-entry.idnum\" , label \"Support\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Medline-entry.idnum.E\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdCitArt\" ,\n" \
"format { asn1 \"Cit-art\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-art.title\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-art.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-art.from.journal\" , form use-template \"StdCitJour\" } ,\n" \
"{ asn1 \"Cit-art.from.book\" , prefix \"(in) \" , form use-template \"StdCitBook\" } ,\n" \
"{ asn1 \"Cit-art.from.proc\" , prefix \"(in) \" , form use-template \"StdCitProc\" } } } } } ,\n" \
"{ name \"StdCitJour\" ,\n" \
"format { asn1 \"Cit-jour\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Cit-jour.title\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-jour.imp\" , form use-template \"StdImprint\" } } } } } ,\n" \
"{ name \"StdCitBook\" ,\n" \
"format { asn1 \"Cit-book\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-book.title\" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-book.coll\" , prefix \"Collection: \" , form use-template \"StdTitle\" } ,\n" \
"{ asn1 \"Cit-book.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-book.imp\" , form use-template \"StdImprint\" } } } } } ,\n" \
"{ name \"StdCitProc\" ,\n" \
"format { asn1 \"Cit-proc\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-proc.book\" , form use-template \"StdCitBook\" } ,\n" \
"{ asn1 \"Cit-proc.meet\" , label \"Meeting \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Meeting.number\" , form text { } } ,\n" \
"{ asn1 \"Meeting.date\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Meeting.place\" , form use-template \"StdAffil\" } } } } } } } } ,\n" \
"{ name \"StdCitPat\" ,\n" \
"format { asn1 \"Cit-pat\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-pat.title\" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.authors\" , form use-template \"StdAuthList\" } ,\n" \
"{ asn1 \"Cit-pat\" , form block {\n" \
"components {\n" \
"{ asn1 \"Cit-pat.country\" , suffix \" \" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.doc-type\" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.number\" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.date-issue\" , prefix \" (\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Cit-pat.app-number\" , prefix \" Appl: \" , form text { } } ,\n" \
"{ asn1 \"Cit-pat.app-date\" , prefix \" (\" , suffix \")\" , form user { printfunc \"StdDatePrint\" } } } } } } } } } ,\n" \
"{ name \"StdIdPat\" ,\n" \
"format { asn1 \"Id-pat\" , form block {\n" \
"components {\n" \
"{ asn1 \"Id-pat.country\" , suffix \" \" , form text { } } ,\n" \
"{ asn1 \"Id-pat.id.number\" , form text { } } ,\n" \
"{ asn1 \"Id-pat.id.app-number\" , prefix \"Appl: \" , form text { } } } } } } ,\n" \
"{ name \"StdCitLet\" ,\n" \
"format { asn1 \"Cit-let\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cit-let.type\" , prefix \"[\" , suffix \"]\" , form enum { } } ,\n" \
"{ asn1 \"Cit-let.man-id\" , form text { } } ,\n" \
"{ asn1 \"Cit-let.cit\" , form use-template \"StdCitBook\" } } } } } ,\n" \
"{ name \"StdPubEquiv\" ,\n" \
"format { asn1 \"Pub-equiv\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pub-equiv.E\" , form use-template \"StdPub\" } } } } } ,\n" \
"{ name \"StdTitle\" ,\n" \
"format { asn1 \"Title\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Title.E.trans\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Title.E.name\" , form text { } } ,\n" \
"{ asn1 \"Title.E.tsub\" , form text { } } ,\n" \
"{ asn1 \"Title.E.abr\" , form text { } } ,\n" \
"{ asn1 \"Title.E.iso-jta\" , form text { } } ,\n" \
"{ asn1 \"Title.E.ml-jta\" , label \"MEDLINE\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.jta\" , label \"jta\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.issn\" , label \"ISSN\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.coden\" , label \"CODEN\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Title.E.isbn\" , label \"ISBN\" , prefix \": \" , form text { } } } } } } ,\n" \
"{ name \"StdAuthList\" ,\n" \
"format { asn1 \"Auth-list\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Auth-list\" , form user { printfunc \"StdAuthListNamesPrint\" } } ,\n" \
"{ asn1 \"Auth-list.affil\" , form use-template \"StdAffil\" } } } } } ,\n" \
"{ name \"StdAffil\" ,\n" \
"format { asn1 \"Affil\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Affil.str\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.affil\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.div\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.street\" , form text { } } ,\n" \
"{ asn1 \"Affil.std\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Affil.std.city\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.sub\" , form text { } } ,\n" \
"{ asn1 \"Affil.std.country\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdImprint\" ,\n" \
"format { asn1 \"Imprint\" , form block {\n" \
"components {\n" \
"{ asn1 \"Imprint.date\" , prefix \"(\" , suffix \") \" , form user { printfunc \"StdDatePrint\" } } ,\n" \
"{ asn1 \"Imprint.volume\" , form text { } } ,\n" \
"{ asn1 \"Imprint.issue\" , prefix \" (\" , suffix \")\" , form text { } } ,\n" \
"{ asn1 \"Imprint.section\" , prefix \" (\" , suffix \")\" , form text { } } ,\n" \
"{ asn1 \"Imprint.part-sup\" , prefix \" (\" , suffix \")\" , form text { } } ,\n" \
"{ asn1 \"Imprint.pages\" , prefix \": \" , form text { } } ,\n" \
"{ asn1 \"Imprint.prepub\" , prefix \" (\" , suffix \")\" , form enum { } } ,\n" \
"{ asn1 \"Imprint.pub\" , label \"\nPublisher: \" , form use-template \"StdAffil\" } ,\n" \
"{ asn1 \"Imprint.cprt\" , label \" Copyright: \" , form user { printfunc \"StdDatePrint\" } } } } } } ,\n" \
"{ name \"StdSeqFeatXref\" ,\n" \
"format { asn1 \"SeqFeatXref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"SeqFeatXref.id\" , label \"Id=\" , form use-template \"StdFeatId\" } ,\n" \
"{ asn1 \"SeqFeatXref.data\" , form use-template \"StdSeqFeatData\" } } } } } ,\n" \
"{ name \"StdOrgRef\" ,\n" \
"format { asn1 \"Org-ref\" , label \"Org-ref\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref\" , form block {\n" \
"separator \" \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.taxname\" , form text { } } ,\n" \
"{ asn1 \"Org-ref.common\" , prefix \"(\" , suffix \")\" , form text { } } } } } ,\n" \
"{ asn1 \"Org-ref.mod\" , label \"Modifiers\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.mod.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Org-ref.db\" , label \"Cross Reference\" , prefix \" \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.db.E\" , prefix \"(\" , suffix \")\" , form use-template \"StdDbtag\" } } } } ,\n" \
"{ asn1 \"Org-ref.syn\" , label \"Synonyms\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Org-ref.syn.E\" , form text { } } } } } } } } } ,\n" \
"{ name \"StdBioSource\" ,\n" \
"format { asn1 \"BioSource\" , label \"BioSource\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"BioSource.genome\" , form enum { } } ,\n" \
"{ asn1 \"BioSource.org\" , label \"Organism\" , form use-template \"StdOrgRef\" } } } } } ,\n" \
"{ name \"StdCdRegion\" ,\n" \
"format { asn1 \"Cdregion\" , label \"Cdregion\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Cdregion.orf\" , form boolean {\n" \
"true \"Uncharacterized Open Reading Frame\" } } ,\n" \
"{ asn1 \"Cdregion.frame\" , label \"Reading Frame = \" , form enum { } } ,\n" \
"{ asn1 \"Cdregion.code\" , label \"Genetic Code: \" , suffix \";\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Genetic-code.E.name\" , form text { } } ,\n" \
"{ asn1 \"Genetic-code.E.id\" , label \"id= \" , form text { } } } } } ,\n" \
"{ asn1 \"Cdregion.conflict\" , form boolean {\n" \
"true \"Translation conflicts with protein sequence\" } } ,\n" \
"{ asn1 \"Cdregion.stops\" , prefix \"Translation contains \" , suffix \" stop codons\" , form text { } } ,\n" \
"{ asn1 \"Cdregion.gaps\" , prefix \"Translation contains \" , suffix \" gaps when aligned to protein\" , form text { } } ,\n" \
"{ asn1 \"Cdregion.mismatch\" , prefix \"Translation contains \" , suffix \" mismatches when aligned to protein\" , form text { } } } } } } ,\n" \
"{ name \"StdProtRef\" ,\n" \
"format { asn1 \"Prot-ref\" , label \"Prot-ref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.name\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.name.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Prot-ref.desc\" , prefix \"[\" , suffix \"]\" , form text { } } ,\n" \
"{ asn1 \"Prot-ref.processed\" , form enum { } } ,\n" \
"{ asn1 \"Prot-ref.ec\" , label \"ec\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.ec.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Prot-ref.activity\" , label \"activity\" , prefix \": \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.activity.E\" , form text { } } } } } ,\n" \
"{ asn1 \"Prot-ref.db\" , label \"Cross Reference\" , prefix \" \" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Prot-ref.db.E\" , prefix \"(\" , suffix \")\" , form use-template \"StdDbtag\" } } } } } } } } ,\n" \
"{ name \"StdRNARef\" ,\n" \
"format { asn1 \"RNA-ref\" , label \"RNA-ref\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"RNA-ref.type\" , form enum { } } ,\n" \
"{ asn1 \"RNA-ref.pseudo\" , form boolean {\n" \
"true \"This is an RNA pseudogene.\" } } ,\n" \
"{ asn1 \"RNA-ref.ext.name\" , form text { } } } } } } ,\n" \
"{ name \"StdPubdesc\" ,\n" \
"format { asn1 \"Pubdesc\" , label \"Pubdesc\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pubdesc.pub\" , form use-template \"StdPubEquiv\" } ,\n" \
"{ asn1 \"Pubdesc\" , prefix \"In this article:\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Pubdesc.name\" , label \"name=\" , form text { } } ,\n" \
"{ asn1 \"Pubdesc.fig\" , label \"figure=\" , form text { } } ,\n" \
"{ asn1 \"Pubdesc.poly-a\" , form boolean {\n" \
"true \"poly(A) shown\" } } ,\n" \
"{ asn1 \"Pubdesc.maploc\" , label \"map location=\" , form text { } } ,\n" \
"{ asn1 \"Pubdesc.num\" , form use-template \"StdNumbering\" } ,\n" \
"{ asn1 \"Pubdesc.numexc\" , form boolean {\n" \
"true \"numbering inconsistent\" } } } } } ,\n" \
"{ asn1 \"Pubdesc.comment\" , form text { } } } } } } ,\n" \
"{ name \"StdImpFeat\" ,\n" \
"format { asn1 \"Imp-feat.key\" , label \"Imp-feat\" , form text { } } } ,\n" \
"{ name \"StdRsiteRef\" ,\n" \
"format { asn1 \"Rsite-ref\" , label \"Rsite-ref\" , form block {\n" \
"components {\n" \
"{ asn1 \"Rsite-ref.str\" , form text { } } ,\n" \
"{ asn1 \"Rsite-ref.std\" , form use-template \"StdDbtag\" } } } } } ,\n" \
"{ name \"StdTxInit\" ,\n" \
"format { asn1 \"Txinit\" , label \"TxInit\" , form block {\n" \
"components {\n" \
"{ asn1 \"Txinit.name\" , form text { } } } } } } ,\n" \
"{ name \"StdNumbering\" ,\n" \
"format { asn1 \"Numbering\" , label \"Numbering\" , form null NULL } } ,\n" \
"{ name \"StdGBBlock\" ,\n" \
"format { asn1 \"GB-block\" , label \"GenBank-block\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"GB-block.extra-accessions\" , label \"Extra accessions\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"GB-block.extra-accessions.E\" , form text { } } } } } ,\n" \
"{ asn1 \"GB-block.keywords\" , label \"Keywords\" , prefix \" (\" , suffix \")\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"GB-block.keywords.E\" , form text { } } } } } ,\n" \
"{ asn1 \"GB-block.source\" , label \"Source: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.origin\" , label \"Origin: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.div\" , label \"Division: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.taxonomy\" , label \"Taxonomy: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.date\" , label \"Date: \" , form text { } } ,\n" \
"{ asn1 \"GB-block.entry-date\" , label \"Entry date: \" , form user { printfunc \"StdDatePrint\" } } } } } } ,\n" \
"{ name \"StdFeatId\" ,\n" \
"format { asn1 \"Feat-id\" , form block {\n" \
"components {\n" \
"{ asn1 \"Feat-id.gibb\" , label \"GenInfo Backbone: \" , form text { } } ,\n" \
"{ asn1 \"Feat-id.giim.id\" , label \"GenInfo Import Id: \" , form text { } } ,\n" \
"{ asn1 \"Feat-id.local\" , label \"Local: \" , form use-template \"StdObjectId\" } ,\n" \
"{ asn1 \"Feat-id.general\" , form use-template \"StdDbtag\" } } } } } ,\n" \
"{ name \"StdSeqFeatCommon\" ,\n" \
"format { asn1 \"Seq-feat\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.id\" , label \"Id=\" , form use-template \"StdFeatId\" } ,\n" \
"{ asn1 \"Seq-feat.title\" , form text { } } ,\n" \
"{ asn1 \"Seq-feat\" , suffix \";\" , form block {\n" \
"separator \", \" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.partial\" , form boolean {\n" \
"true \"Partial\" } } ,\n" \
"{ asn1 \"Seq-feat.except\" , form boolean {\n" \
"true \"Biological Exception\" } } ,\n" \
"{ asn1 \"Seq-feat.exp-ev\" , label \"Evidence\" , prefix \" is \" , form enum { } } } } } ,\n" \
"{ asn1 \"Seq-feat.comment\" , form text { } } ,\n" \
"{ asn1 \"Seq-feat.ext\" , form use-template \"StdUserObj\" } ,\n" \
"{ asn1 \"Seq-feat.qual\" , label \"Qualifiers\" , prefix \"\\n\" , form block {\n" \
"separator \"\\n\" ,\n" \
"components {\n" \
"{ asn1 \"Seq-feat.qual.E\" , prefix \"/\" , form block {\n" \
"separator \"= \" ,\n" \
"components {\n" \
"{ asn1 \"Gb-qual.qual\" , form text { } } ,\n" \
"{ asn1 \"Gb-qual.val\" , form text { } } } } } } } } } } } } ,\n" \
"{ name \"StdDbtag\" ,\n" \
"format { asn1 \"Dbtag\" , form block {\n" \
"components {\n" \
"{ asn1 \"Dbtag.db\" , suffix \": \" , form text { } } ,\n" \
"{ asn1 \"Dbtag.tag\" , form use-template \"StdObjectId\" } } } } } ,\n" \
"{ name \"StdObjectId\" ,\n" \
"format { asn1 \"Object-id\" , form block {\n" \
"components {\n" \
"{ asn1 \"Object-id.id\" , form text { } } ,\n" \
"{ asn1 \"Object-id.str\" , form text { } } } } } } };\n";
#else
CharPtr objPrtMemStr = "";
#endif

typedef Boolean (*MatchBioseqToTableData) PROTO ((BioseqPtr, ValNodePtr, Int4));
typedef Boolean (*BioseqAlreadyHasTableData) PROTO ((BioseqPtr, ValNodePtr, Int4));

static CharPtr GetMatchString (ValNodePtr columns, Int4 match_pos) 
{
  ValNodePtr vnp;

  vnp = columns;
  while (match_pos > 0 && vnp != NULL) {
    vnp = vnp->next;
    match_pos--;
  }
  if (vnp == NULL) {
    return NULL;
  } else {
    return vnp->data.ptrvalue;
  }
}

static Boolean MatchBioseqToTableDataByID (BioseqPtr bsp, ValNodePtr columns, Int4 match_pos)
{
  Char       seqid[100];
  Int4       seqid_len;
  SeqIdPtr   sip;
  CharPtr    match_string;
  Boolean    rval = FALSE;

  if (bsp == NULL || columns == NULL || match_pos < 0) return rval;

  match_string = GetMatchString (columns, match_pos);
  if (match_string == NULL) return rval;

  sprintf (seqid, "gb|");
  seqid_len = StringLen (match_string);
  if (seqid_len > 96) {
    seqid_len = 96;
  }
  StringNCpy (seqid + 3, match_string, seqid_len);
  seqid [seqid_len + 3] = 0;
  sip = MakeSeqID (seqid);
  if (SeqIdIn (sip, bsp->id)) {
    rval = TRUE;
  }
  return rval;
}

static Boolean MatchBioseqToTableDataByOrgname (BioseqPtr bsp, ValNodePtr columns, Int4 match_pos)
{
  CharPtr           match_string;
  Boolean           rval = FALSE;
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  BioSourcePtr      biop;

  if (bsp == NULL || columns == NULL || match_pos < 0) return rval;

  match_string = GetMatchString (columns, match_pos);
  if (match_string == NULL) return rval;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL && !rval) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL && biop->org != NULL && StringCmp (biop->org->taxname, match_string) == 0) {
      rval = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
  }

  return rval;
}


typedef struct bioseqtabledataassoc {
  BioseqPtr    bsp;
  TableLinePtr tlp;
} BioseqTableDataAssocData, PNTR BioseqTableDataAssocPtr;

typedef struct buildtabledataassoc {
  MatchBioseqToTableData match_func;
  Int4                   match_pos;
  ValNodePtr             table_data;
  ValNodePtr             assoc_list;
} BuildBioseqTableDataAssocData, PNTR BuildBioseqTableDataAssocPtr;

static void BuildTableDataAssociationCallback (BioseqPtr bsp, Pointer userdata)
{
  BuildBioseqTableDataAssocPtr bp;
  ValNodePtr                   vnp;
  TableLinePtr                 tlp;
  BioseqTableDataAssocPtr      bap;

  bp = (BuildBioseqTableDataAssocPtr) userdata;
  if (bsp == NULL || bp == NULL || bp->match_func == NULL || bp->table_data == NULL) return;
  for (vnp = bp->table_data; vnp != NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == NULL) continue;
    tlp = (TableLinePtr) vnp->data.ptrvalue;
    if (bp->match_func(bsp, tlp->parts, bp->match_pos)) {
      bap = (BioseqTableDataAssocPtr) MemNew (sizeof (BioseqTableDataAssocData));
      bap->bsp = bsp;
      bap->tlp = tlp;
      ValNodeAddPointer (&(bp->assoc_list), 0, bap);
    }
  }
}

static ValNodePtr BuildTableDataAssociation (SeqEntryPtr sep, ValNodePtr table_data, MatchBioseqToTableData match_func, Int4 match_pos)
{
  BuildBioseqTableDataAssocData bad;

  if (sep == NULL || table_data == NULL || match_func == NULL) return NULL;
  bad.assoc_list = NULL;
  bad.match_func = match_func;
  bad.match_pos = match_pos;
  bad.table_data = table_data;

  VisitBioseqsInSep (sep, &bad, BuildTableDataAssociationCallback);
  return bad.assoc_list;
}

static ValNodePtr 
CheckTableDataAssociation 
(ValNodePtr assoc_list,
 BaseFormPtr bfp,
 ValNodePtr table_data,
 Int4 match_pos,
 BioseqAlreadyHasTableData already_func,
 CharPtr already_fmt,
 CharPtr more_than_one_fmt,
 BoolPtr data_needed,
 Int4    num_columns)
{
  ValNodePtr vnp, vnp_match;
  BioseqTableDataAssocPtr bap1, bap2;
  ValNodePtr duplicated_bioseqs = NULL, already_has_list = NULL, err_list = NULL;
  ValNodePtr blanks = NULL;
  ClickableItemPtr cip;
  TableLinePtr     tlp;
  Boolean          found;
  CharPtr          no_match_fmt = "No match for %s";
  CharPtr          match_string;
  Int4             pos;

  if (assoc_list == NULL) {
    Message (MSG_ERROR, "No matches found!");
    return NULL;
  }
  for (vnp = assoc_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next) {
    bap1 = (BioseqTableDataAssocPtr) vnp->data.ptrvalue;
    for (vnp_match = vnp->next; vnp_match != NULL; vnp_match = vnp_match->next) {
      bap2 = (BioseqTableDataAssocPtr) vnp_match->data.ptrvalue;
      if (bap1->bsp == bap2->bsp) {
        for (vnp_match = duplicated_bioseqs;
             vnp_match != NULL && bap1->bsp != vnp_match->data.ptrvalue;
             vnp_match = vnp_match->next) {}
        if (vnp_match == NULL) {
          ValNodeAddPointer (&(duplicated_bioseqs), OBJ_BIOSEQ, bap1->bsp);
        }
        break;
      }
    }
  }
  for (vnp = assoc_list; vnp != NULL; vnp = vnp->next) {
    bap1 = (BioseqTableDataAssocPtr) vnp->data.ptrvalue;
    if (already_func != NULL && already_func (bap1->bsp, NULL, 0)) {
      ValNodeAddPointer (&already_has_list, OBJ_BIOSEQ, bap1->bsp);
    }
    if (data_needed != NULL) {
      for (vnp_match = bap1->tlp->parts, pos = 0;
           vnp_match != NULL && pos < num_columns;
           vnp_match = vnp_match->next, pos++) {
        if (data_needed[pos] && StringHasNoText (vnp_match->data.ptrvalue)) {
          break;
        }
      }
      if (pos < num_columns) {
        ValNodeAddPointer (&blanks, OBJ_BIOSEQ, bap1->bsp);
      }
    }
  }

  if (duplicated_bioseqs != NULL) {
    if (more_than_one_fmt == NULL) {
      cip = NewClickableItem (TABLE_DATA_MULTIPLE_VALUES, "%d sequences will receive data from more than one line.", duplicated_bioseqs);
    } else {
      cip = NewClickableItem (TABLE_DATA_MULTIPLE_VALUES, more_than_one_fmt, duplicated_bioseqs);
    }
    ValNodeAddPointer (&err_list, 0, cip);
  }

  if (already_has_list != NULL) {
    if (already_fmt == NULL) {
      cip = NewClickableItem (TABLE_DATA_ALREADY_HAS, "%d sequences already have data", already_has_list);
    } else {
      cip = NewClickableItem (TABLE_DATA_ALREADY_HAS, already_fmt, already_has_list);
    }
    ValNodeAddPointer (&err_list, 0, cip);
  }

  if (blanks != NULL) {
    cip = NewClickableItem (TABLE_DATA_CELL_BLANK, "%d sequences have blank values in the table", blanks);
    ValNodeAddPointer (&err_list, 0, cip);
  }

  for (vnp = table_data; vnp != NULL; vnp = vnp->next) {
    tlp = vnp->data.ptrvalue;
    found = FALSE;
    for (vnp_match = assoc_list; vnp_match != NULL && !found; vnp_match = vnp_match->next) {
      bap1 = (BioseqTableDataAssocPtr) vnp_match->data.ptrvalue;
      if (bap1->tlp == tlp) {
        found = TRUE;
      }
    }
    if (!found) {
      match_string = GetMatchString (tlp->parts, match_pos);
      if (match_string != NULL) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (cip, 0, sizeof (ClickableItemData));
        cip->clickable_item_type = TABLE_DATA_NOT_FOUND;
        cip->description = MemNew ((StringLen (no_match_fmt) + StringLen (match_string)) * sizeof (Char));
        sprintf (cip->description, no_match_fmt, match_string);
        ValNodeAddPointer (&err_list, 0, cip);
      }
    }
  }

  return err_list;
}

static Boolean BioseqHasGenomeProjectID (BioseqPtr bsp, ValNodePtr columns, Int4 match_pos)
{
  Boolean     rval = FALSE;
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
 
  if (bsp == NULL) return rval;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL && !rval) {
    rval = IsGenomeProjectIDDescriptor(sdp);
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }
  return rval;
}

static void ApplyGenomeProjectIDByTableDataAssociationList (ValNodePtr assoc_list, Int4 project_id_col, Boolean do_replace, Boolean blanks_erase, Uint2 entityID)
{
  BioseqTableDataAssocPtr bap;
  ValNodePtr              vnp;
  Int4                    projectID;
  SeqDescrPtr             sdp;
  UserObjectPtr           uop;
  ObjectIdPtr             oip;
  UserFieldPtr            ufp, ufp_last, ufp_next;
  CharPtr                 match_string;
  ObjValNodePtr           ovn;
   
  for (vnp = assoc_list; vnp != NULL; vnp = vnp->next) {
    bap = (BioseqTableDataAssocPtr) vnp->data.ptrvalue;
    sdp = NULL;
    if (BioseqHasGenomeProjectID (bap->bsp, NULL, 0)) {
      if (!do_replace) {
        continue;
      } else {
        sdp = GetGenomeProjectIDDescriptor(bap->bsp);
      }
    }
    match_string = GetMatchString (bap->tlp->parts, project_id_col);
    if (StringHasNoText (match_string)) {
      if (blanks_erase) {
        if (sdp != NULL && sdp->extended != 0) {
          ovn = (ObjValNodePtr) sdp;
          ovn->idx.deleteme = TRUE;
        }
      }
    } else {
      projectID = atoi (match_string);
      if (sdp == NULL) {
        sdp = SeqDescrNew (bap->bsp->descr);
        if (bap->bsp->descr == NULL) bap->bsp->descr = sdp;
        sdp->choice = Seq_descr_user;          
        uop = CreateGenomeProjectsDBUserObject ();
        AddIDsToGenomeProjectsDBUserObject (uop, projectID, 0);
        sdp->data.ptrvalue = uop;
      } else {
        uop = sdp->data.ptrvalue;
        if (uop == NULL) {
          uop = CreateGenomeProjectsDBUserObject ();
          sdp->data.ptrvalue = uop;
        }
        ufp_last = NULL;
        ufp = uop->data;
        while (ufp != NULL) {
          oip = ufp->label;
          if (oip != NULL && StringCmp (oip->str, "ProjectID") == 0) {
            break;
          } else {
            ufp_last = ufp;
            ufp = ufp->next;
          }
        }
        if (ufp == NULL) {
          ufp_next = NULL;
        } else {
          ufp_next = ufp->next;
        }

        if (ufp->choice != 2) {
          ufp_next = ufp->next;
          ufp = UserFieldFree(ufp);
        }
        if (ufp == NULL) {
          ufp = UserFieldNew ();
          oip = ObjectIdNew ();
          oip->str = StringSave ("ProjectID");
          ufp->label = oip;
          ufp->choice = 2; /* integer */
          if (ufp_last == NULL) {
            uop->data = ufp;
          } else {
            ufp_last->next = ufp;
          }
        }
        ufp->data.intvalue = projectID;
        ufp->next = ufp_next;
      }
    }
  }
  DeleteMarkedObjects (entityID, 0, NULL);
}

typedef struct genomeprojectid {
  FORM_MESSAGE_BLOCK
  PopuP id_column;
  PopuP match_column;
  PopuP gpid_column;

  ValNodePtr table_data;
  BaseFormPtr bfp;
} GenomeProjectIdData, PNTR GenomeProjectIdPtr;

static void CleanupGenomeProjectIDForm (GraphiC g, VoidPtr data)

{
  GenomeProjectIdPtr gpip;

  gpip = (GenomeProjectIdPtr) data;
  if (gpip != NULL)
  {
    CleanUpTableData (gpip->table_data);
  }
  StdCleanupFormProc (g, data);
}

static void ApplyGenomeProjectIDs (ButtoN b)
{
  GenomeProjectIdPtr gpip;
  ValNodePtr         assoc_list = NULL;
  SeqEntryPtr        sep;
  MatchBioseqToTableData match_func = NULL;
  Int4                   match_pos, project_id_col, col_num;
  Boolean                ok = FALSE;
  Boolean                skip_already_has = FALSE, blanks_erase = FALSE;
  ValNodePtr             err_list = NULL;
  SeqEntryPtr            oldscope;
  TableLinePtr           tlp;
  BoolPtr                data_needed = NULL;

  gpip = (GenomeProjectIdPtr) GetObjectExtra (b);
  if (gpip == NULL) return;

  match_pos = GetValue (gpip->match_column);
  if (match_pos < 1) return;
  match_pos--;
  project_id_col = GetValue (gpip->gpid_column);
  if (project_id_col < 1) return;
  project_id_col--;
  if (GetValue (gpip->id_column) == 1) {
    match_func = MatchBioseqToTableDataByID;
  } else {
    match_func = MatchBioseqToTableDataByOrgname;
  }

  tlp = gpip->table_data->data.ptrvalue;
  if (tlp != NULL) {
    data_needed = (BoolPtr) MemNew (sizeof (Boolean) * tlp->num_parts);
    for (col_num = 0; col_num < tlp->num_parts; col_num++) {
      if (col_num == project_id_col) {
        data_needed[col_num] = TRUE;
      } else {
        data_needed[col_num] = FALSE;
      }
    }
  }

  sep = GetTopSeqEntryForEntityID (gpip->input_entityID);
  oldscope = SeqEntrySetScope (sep);
  assoc_list = BuildTableDataAssociation (sep, gpip->table_data, match_func, match_pos);
  err_list = CheckTableDataAssociation (assoc_list, gpip->bfp, gpip->table_data, match_pos, 
                                 BioseqHasGenomeProjectID,
                                 "%d sequences already have a Genome Project ID.",
                                 "%d sequences will get more than one Genome Project ID from the table.",
                                 data_needed, tlp->num_parts);



  if (err_list == NULL || GetTableOptions (gpip->bfp, err_list, "Errors Matching Sequences from File", "", "",
                       "Skip sequences that already have Genome Project IDs",
                       "Erase Genome Project IDs when table cell is blank",
                       &skip_already_has, &blanks_erase)) {
    ApplyGenomeProjectIDByTableDataAssociationList (assoc_list, project_id_col, !skip_already_has, blanks_erase, gpip->input_entityID);
    ok = TRUE;
  }
  assoc_list = ValNodeFree (assoc_list);
  if (ok) {
    ObjMgrSetDirtyFlag (gpip->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, gpip->input_entityID, 0, 0);
    ArrowCursor ();
    Update ();
    Remove (gpip->form);
  }
  SeqEntrySetScope (oldscope);
}


extern void LoadGenomeProjectIDsFromFile (IteM i)
{
  BaseFormPtr   bfp;
  WindoW        w;
  GenomeProjectIdPtr gpip;
  GrouP              h, g, c;
  ButtoN             b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  gpip = (GenomeProjectIdPtr) MemNew (sizeof(GenomeProjectIdData));
  gpip->table_data = ReadTableData ();
  if (gpip->table_data == NULL) {
    gpip = MemFree (gpip);
    return;
  }
  gpip->bfp = bfp;

  w = FixedWindow (-50, -33, -10, -10, "Add Genome Project IDs", StdCloseWindowProc);
  gpip->form = (ForM) w;
  SetObjectExtra (w, gpip, CleanupGenomeProjectIDForm);
  gpip->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);

  g = HiddenGroup (h, 6, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Match", 0, dialogTextHeight, systemFont, 'c');
  gpip->match_column = PopupList (g, TRUE, NULL);
  FormatPopupWithTableDataColumns (gpip->match_column, gpip->table_data);
  SetValue (gpip->match_column, 1);

  StaticPrompt (g, "To", 0, dialogTextHeight, systemFont, 'c');
  gpip->id_column = PopupList (g, TRUE, NULL);
  PopupItem (gpip->id_column, "Accession");
  PopupItem (gpip->id_column, "Tax Name");
  SetValue (gpip->id_column, 1);
  StaticPrompt (g, "Use Project ID from", 0, dialogTextHeight, systemFont, 'c');
  gpip->gpid_column = PopupList (g, TRUE, NULL);
  FormatPopupWithTableDataColumns (gpip->gpid_column, gpip->table_data);
  SetValue (gpip->gpid_column, 1);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", ApplyGenomeProjectIDs);
  SetObjectExtra (b, gpip, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, gpip, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


static void SuppressGenesOnFeaturesInsideMobileElementsCallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr        mobile_element, sfp, gene;
  SeqMgrFeatContext fcontext_m, fcontext_s, fcontext_g;
  GeneRefPtr        grp;
  SeqFeatXrefPtr    xref;

  if (bsp == NULL) return;

  for (mobile_element = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_IMP, FEATDEF_repeat_region, &fcontext_m);
       mobile_element != NULL;
       mobile_element = SeqMgrGetNextFeature (bsp, mobile_element, SEQFEAT_IMP, FEATDEF_repeat_region, &fcontext_m)) {
    if (!IsMobileElement (mobile_element)) continue;
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext_s);
         sfp != NULL && fcontext_s.left <= fcontext_m.right;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext_s)) {
      if (SeqLocCompare (sfp->location, mobile_element->location) == SLC_A_IN_B
          && SeqMgrGetGeneXref (sfp) == NULL) {
        /* suppress gene if gene not contained in mobile_element */
        gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext_g);          
        if (gene != NULL && (fcontext_g.left <= fcontext_m.left || fcontext_g.right >= fcontext_m.right)) {
          grp = GeneRefNew ();
          if (grp != NULL) {
	          xref = SeqFeatXrefNew ();
	          xref->data.choice = SEQFEAT_GENE;
	          xref->data.value.ptrvalue = grp;
	          xref->next = sfp->xref;
	          sfp->xref = xref;
          }
        }
      }
    }
  }
}


extern void SuppressGenesOnFeaturesInsideMobileElements (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitBioseqsInSep (sep, NULL, SuppressGenesOnFeaturesInsideMobileElementsCallback);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}


/* Importing Protein ID Table */

typedef struct featurefieldcolumnchoice {
  Int4            match_type;
  ValNodePtr      field_choice;
  Boolean         change_mrna;
  Boolean         erase_when_blank;
  ExistingTextPtr etp;
} FeatureFieldColumnChoiceData, PNTR FeatureFieldColumnChoicePtr;

static FeatureFieldColumnChoicePtr FeatureFieldColumnChoiceFree (FeatureFieldColumnChoicePtr f)
{
  if (f != NULL) {
    f->field_choice = ValNodeFree (f->field_choice);
    f->etp = MemFree (f->etp);
    f = MemFree (f);
  }
  return f;
}


static ValNodePtr GetFeaturesForGene (SeqFeatPtr gene, Uint1 featdef)
{
  BioseqPtr bsp;
  SeqFeatPtr sfp;
  ValNodePtr feat_list = NULL;
  SeqMgrFeatContext fcontext;
  Int4              start, stop, swap;

  if (gene == NULL) return NULL;

  bsp = BioseqFindFromSeqLoc (gene->location);
  start = SeqLocStart (gene->location);
  stop = SeqLocStop (gene->location);
  if (stop < start) 
  {
    swap = start;
    start = stop;
    stop = swap;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
       sfp != NULL && fcontext.left < stop;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
  {
    if (fcontext.right >= start && gene == GetGeneForFeature (sfp))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}

static ValNodePtr GetFeatureListForProteinBioseq (Uint1 featdef, BioseqPtr bsp)
{
  ValNodePtr feat_list = NULL;
  SeqFeatPtr sfp, cds;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || !ISA_aa (bsp->mol)) 
  {
    return NULL;
  }

  if (featdef == FEATDEF_PROT || featdef == FEATDEF_mat_peptide_aa)
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  else
  {
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (cds != NULL) 
    {
      sfp = NULL;
      if (featdef == FEATDEF_CDS)
      {
        sfp = cds;
      }
      else if (featdef == FEATDEF_GENE)
      {
        sfp = GetGeneForFeature (cds);
      }
      else if (featdef == FEATDEF_mRNA)
      {
        sfp = SeqMgrGetOverlappingmRNA (cds->location, &fcontext);
      }
      if (sfp != NULL)
      {
        ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
      }
    }
  }
  return feat_list;
}


static ValNodePtr GetFeatureListForGene (Uint1 featdef, SeqFeatPtr gene)
{
  ValNodePtr feat_list = NULL, cds_list, vnp;
  SeqFeatPtr sfp, cds;
  SeqMgrFeatContext fcontext;
  BioseqPtr         protbsp;

  if (gene == NULL) 
  {
    return NULL;
  }

  if (featdef == FEATDEF_GENE)
  {
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, gene);
  }
  else if (featdef == FEATDEF_mRNA || featdef == FEATDEF_CDS)
  {
    feat_list = GetFeaturesForGene (gene, featdef);
  }
  else if (featdef == FEATDEF_PROT || featdef == FEATDEF_mat_peptide_aa)
  {
    cds_list = GetFeaturesForGene (gene, FEATDEF_CDS);
    for (vnp = cds_list; vnp != NULL; vnp = vnp->next) 
    {
      cds = vnp->data.ptrvalue;
      if (cds != NULL)
      {
        protbsp = BioseqFindFromSeqLoc (cds->product);
        for (sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, featdef, &fcontext);
             sfp != NULL;
             sfp = SeqMgrGetNextFeature (protbsp, sfp, 0, featdef, &fcontext))
        {
          ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
        }
      }
    }
    cds_list = ValNodeFree (cds_list);
  }

  return feat_list;
}


static SeqFeatPtr GetFeatureForFieldByMatchList (ValNodePtr field, ValNodePtr match_list, ValNodePtr PNTR errors)
{
  Uint1 featdef;
  ValNodePtr vnp;
  ValNodePtr tmp;
  CharPtr    msg;
  CharPtr    no_match_fmt = "No feature found for %s";
  CharPtr    too_many_fmt = "%d features found for %s (only one allowed!)";
  CharPtr    bad_comb_fmt = "Multiple match columns specify different features (%s, %s)!";
  SeqFeatPtr sfp = NULL;

  if (field == NULL || match_list == NULL || errors == NULL) return NULL;

  featdef = (Uint1)FeatDefTypeFromFieldList (field);
  if (featdef == FEATDEF_BAD || featdef == FEATDEF_ANY) 
  { 
    return NULL;
  }

  for (vnp = match_list; vnp != NULL; vnp = vnp->next)
  {
    tmp = NULL;
    if (vnp->choice == OBJ_BIOSEQ) 
    {
      tmp = GetFeatureListForProteinBioseq (featdef, vnp->data.ptrvalue);
    }
    else
    {
      tmp = GetFeatureListForGene (featdef, vnp->data.ptrvalue);
    }
    if (tmp == NULL)      
    {
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (vnp->data.ptrvalue)));
      sprintf (msg, no_match_fmt, vnp->data.ptrvalue);
      ValNodeAddPointer (errors, 0, msg);
    }
    else if (tmp->next != NULL)
    {
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (too_many_fmt) + StringLen (vnp->data.ptrvalue) + 15));
      sprintf (msg, too_many_fmt, ValNodeLen (tmp), vnp->data.ptrvalue);
      ValNodeAddPointer (errors, 0, msg);
    }
    else if (sfp != NULL && sfp != tmp->data.ptrvalue)
    {
      msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_comb_fmt) + StringLen (vnp->data.ptrvalue) + StringLen (match_list->data.ptrvalue)));
      sprintf (msg, too_many_fmt, match_list->data.ptrvalue, vnp->data.ptrvalue);
      ValNodeAddPointer (errors, 0, msg);
    }
    else
    {
      sfp = tmp->data.ptrvalue;
    }
    tmp = ValNodeFree (tmp);
  }
  return sfp;
}


static void 
GetExistingValueCounts 
(TableLinePtr                     tlp,
 FeatureFieldColumnChoicePtr PNTR f_list,
 Int4                             num_columns,
 ValNodePtr                       match_list,
 Int4Ptr                          existing_text,
 BoolPtr                          has_mrna,
 BoolPtr                          missing_mrna,
 ValNodePtr PNTR                  errors)
{
  Int4       i;
  ValNodePtr vnp;
  SeqFeatPtr sfp, mrna;
  CharPtr    val;

  if (tlp == NULL || f_list == NULL || existing_text == NULL) return;

  vnp = tlp->parts;
  for (i = 0; i < num_columns; i++)
  {
    if (f_list[i]->field_choice != NULL) 
    {
      sfp = GetFeatureForFieldByMatchList (f_list[i]->field_choice, match_list, errors);
      if (sfp != NULL) 
      {
        val = GetCDSGeneProtField (sfp, f_list[i]->field_choice, NULL);
        if (!StringHasNoText (val)) 
        {
          existing_text[i]++;
        }
        val = MemFree (val);
        if (IsCDSetProteinProductChoice (f_list[i]->field_choice) && f_list[i]->change_mrna)
        {
          mrna = GetmRNAForFeature (sfp);
          if (mrna == NULL)
          {
            *missing_mrna = TRUE;
          }
          else
          {
            *has_mrna = TRUE;
          }
        }
      }
    }
    if (vnp != NULL) 
    {
      vnp = vnp->next;
    }
  }
}


static Int4 
ApplyProteinIDTableLine 
(TableLinePtr                     tlp,
 FeatureFieldColumnChoicePtr PNTR f_list,
 Int4                             num_columns,
 ValNodePtr                       match_list)
{
  Int4       index;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  ApplyValueData avd;
  ValNodePtr     errors = NULL;
  Int4           fields_affected = 0;

  if (tlp == NULL || f_list == NULL || match_list == NULL) return 0;

  index = 0;
  vnp = tlp->parts;
  while (index < num_columns)
  {
    if (f_list[index]->field_choice != NULL) 
    {
      sfp = GetFeatureForFieldByMatchList (f_list[index]->field_choice, match_list, &errors);
      if (sfp != NULL)
      {
        if (vnp == NULL || StringHasNoText (vnp->data.ptrvalue)) 
        {
          if (f_list[index]->erase_when_blank)
          {
            RemoveCDSGeneProtField (sfp, f_list[index]->field_choice, NULL);
            fields_affected++;
          }
        }
        else
        {
          avd.etp = f_list[index]->etp;
          avd.field_list = f_list[index]->field_choice;
          avd.new_text = vnp->data.ptrvalue;
          avd.text_to_replace = NULL;
          avd.where_to_replace = EditApplyFindLocation_anywhere;
          SetCDSGeneProtField (sfp, f_list[index]->field_choice, &avd, NULL);
          fields_affected ++;
        }
        /* note - need special case to change mRNA product if CDS product changes */
        if (IsCDSetProteinProductChoice (f_list[index]->field_choice) && f_list[index]->change_mrna)
        {
          if (AdjustmRNAProductToMatchProteinProduct (sfp))
          {
            fields_affected++;
          }
        }
      }
    }
    index++;
    if (vnp != NULL) 
    {
      vnp = vnp->next;
    }
  }
  errors = ValNodeFreeData (errors);
  return fields_affected;
}


typedef struct featurefieldcolumnchoicedlg {
  DIALOG_MESSAGE_BLOCK
  GrouP                    match_or_field;
  PopuP                    match_type;
  DialoG                   field_choice;
  ButtoN                   change_mrna;
  /* for existing text options */
  GrouP                    apply_options;
  ButtoN                   erase_when_blank;
  GrouP                    existing_text_action_grp;
  GrouP                    existing_text_delim_grp;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} FeatureFieldColumnChoiceDlgData, PNTR FeatureFieldColumnChoiceDlgPtr;

static void FeatureFieldColumnChoicePopupChange (PopuP p)
{
  FeatureFieldColumnChoiceDlgPtr dlg;

  dlg = (FeatureFieldColumnChoiceDlgPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void FeatureFieldColumnChoiceFieldChange (Pointer data)
{
  FeatureFieldColumnChoiceDlgPtr dlg;
  ValNodePtr                     vnp;

  dlg = (FeatureFieldColumnChoiceDlgPtr) data;
  if (dlg != NULL) {
    vnp = DialogToPointer (dlg->field_choice);
    if (IsCDSetProteinProductChoice (vnp)) {
      Enable (dlg->change_mrna);
    } else {
      Disable (dlg->change_mrna);
    }
    vnp = ValNodeFree (vnp);
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static void FeatureFieldColumnChoiceGroupChange (GrouP p)
{
  FeatureFieldColumnChoiceDlgPtr dlg;

  dlg = (FeatureFieldColumnChoiceDlgPtr) GetObjectExtra (p);
  if (dlg != NULL) {
    if (GetValue (dlg->match_or_field) == 1) {
      Enable (dlg->match_type);
      Disable (dlg->field_choice);
      Disable (dlg->apply_options);
      Disable (dlg->change_mrna);
    } else {
      Disable (dlg->match_type);
      Enable (dlg->field_choice);
      Enable (dlg->apply_options);
      FeatureFieldColumnChoiceFieldChange (dlg);
    }
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static Pointer FeatureFieldColumnDialogToChoice (DialoG d)
{
  FeatureFieldColumnChoiceDlgPtr dlg;
  FeatureFieldColumnChoicePtr    f;
  Int4                           existing_text_action;
  Int4                           existing_text_delimiter;

  dlg = (FeatureFieldColumnChoiceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  f = (FeatureFieldColumnChoicePtr) MemNew (sizeof (FeatureFieldColumnChoiceData));
  f->change_mrna = FALSE;
  if (GetValue (dlg->match_or_field) == 1) {
    f->match_type = GetValue (dlg->match_type);
    f->etp = NULL;
    f->erase_when_blank = FALSE;
  } else {
    f->field_choice = DialogToPointer (dlg->field_choice);
    if (dlg->erase_when_blank == NULL) {
      f->erase_when_blank = FALSE;
    } else {
      f->erase_when_blank = GetStatus (dlg->erase_when_blank);
    }
    if (IsCDSetProteinProductChoice (f->field_choice)) {
      f->change_mrna = GetStatus (dlg->change_mrna);
    }
    f->etp = (ExistingTextPtr) MemNew (sizeof (ExistingTextData));
    existing_text_action = GetValue (dlg->existing_text_action_grp);
    if (existing_text_action == 1) {
      f->etp->existing_text_choice = eExistingTextChoiceReplaceOld;
    } else if (existing_text_action == 4) {
      f->etp->existing_text_choice = eExistingTextChoiceLeaveOld;
    } else {
      existing_text_delimiter = GetValue (dlg->existing_text_delim_grp);
      if (existing_text_action == 2) {
        f->etp->existing_text_choice = eExistingTextChoiceAppendSemi + existing_text_delimiter - 1;
      } else if (existing_text_action == 3) {
        f->etp->existing_text_choice = eExistingTextChoicePrefixSemi + existing_text_delimiter - 1;
      } else {
        f->etp->existing_text_choice = eExistingTextChoiceCancel;
      }
    }       
  }
  return f;

}


static void ChangeExistingTextActionChoice (GrouP g)
{
  FeatureFieldColumnChoiceDlgPtr dlg;
  Int4                           action_choice;

  dlg = (FeatureFieldColumnChoiceDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  
  action_choice = GetValue (dlg->existing_text_action_grp);
  if (action_choice == 2 || action_choice == 3) {
    Enable (dlg->existing_text_delim_grp);
  } else {
    Disable (dlg->existing_text_delim_grp);
  }
}


static DialoG FeatureFieldColumnChoiceDialog 
(GrouP                    h,
 CharPtr                  title,
 Int4                     num_blank,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  FeatureFieldColumnChoiceDlgPtr dlg;
  GrouP p, g1, g;
  ButtoN b1, b2;
  PrompT ppt;
  CharPtr real_title;
  CharPtr title_fmt = "(%d are blank)%s";

  dlg = (FeatureFieldColumnChoiceDlgPtr) MemNew (sizeof (FeatureFieldColumnChoiceDlgData));

  if (num_blank > 0) {
    real_title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (title) + 15));
    sprintf (real_title, title_fmt, num_blank, title);
  } else {
    real_title = title;
  }

  p = NormalGroup (h, -1, 0, real_title, programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  if (real_title != title) {
    real_title = MemFree (real_title);
  }
  dlg->dialog = (DialoG) p;
  dlg->fromdialog = FeatureFieldColumnDialogToChoice;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g1 = HiddenGroup (p, 2, 0, NULL);
  dlg->match_or_field = HiddenGroup (g1, 0, 2, FeatureFieldColumnChoiceGroupChange);
  SetObjectExtra (dlg->match_or_field, dlg, NULL);
  b1 = RadioButton (dlg->match_or_field, "Match to");
  b2 = RadioButton (dlg->match_or_field, "Apply to");
  SetValue (dlg->match_or_field, 1);

  g = HiddenGroup (g1, 0, 2, NULL);
  dlg->match_type = PopupList (g, TRUE, FeatureFieldColumnChoicePopupChange);
  SetObjectExtra (dlg->match_type, dlg, NULL);
  PopupItem (dlg->match_type, "None");
  PopupItem (dlg->match_type, "ID");
  PopupItem (dlg->match_type, "Gene locus tag");
  SetValue (dlg->match_type, 1);

  dlg->field_choice = CDSGeneProtFieldSelectionDialog (g, FALSE, FeatureFieldColumnChoiceFieldChange, dlg);
  Disable (dlg->field_choice);
  AlignObjects (ALIGN_MIDDLE, (HANDLE) b1, (HANDLE) dlg->match_type, NULL);
  AlignObjects (ALIGN_MIDDLE, (HANDLE) b2, (HANDLE) dlg->field_choice, NULL);

  dlg->change_mrna = CheckBox (p, "Also change mRNA product name", NULL);
  Disable (dlg->change_mrna);

  dlg->apply_options = HiddenGroup (p, -1, 0, NULL);
  if (num_blank > 0) {
    dlg->erase_when_blank = CheckBox (dlg->apply_options, "Erase field when table cell is blank", NULL);
  } else {
    dlg->erase_when_blank = NULL;
  }
  dlg->existing_text_action_grp = HiddenGroup (dlg->apply_options, 4, 0, ChangeExistingTextActionChoice);
  SetGroupSpacing (dlg->existing_text_action_grp, 10, 10);
  SetObjectExtra (dlg->existing_text_action_grp, dlg, NULL);
  RadioButton (dlg->existing_text_action_grp, "Overwrite existing text");
  RadioButton (dlg->existing_text_action_grp, "Append");
  RadioButton (dlg->existing_text_action_grp, "Prefix");
  RadioButton (dlg->existing_text_action_grp, "Ignore new text");
  SetValue (dlg->existing_text_action_grp, 1);

  ppt = StaticPrompt (dlg->apply_options, "Separate new text and old text with", 
                      0, dialogTextHeight, programFont, 'c');
  
  dlg->existing_text_delim_grp = HiddenGroup (dlg->apply_options, 4, 0, NULL);
  SetGroupSpacing (dlg->existing_text_delim_grp, 10, 10);
  RadioButton (dlg->existing_text_delim_grp, "Semicolon");
  RadioButton (dlg->existing_text_delim_grp, "Space");
  RadioButton (dlg->existing_text_delim_grp, "Colon");
  RadioButton (dlg->existing_text_delim_grp, "Do not separate");
  SetValue (dlg->existing_text_delim_grp, 1);
  Disable (dlg->existing_text_delim_grp);

  Disable (dlg->apply_options);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->existing_text_action_grp,
                              (HANDLE) ppt,
                              (HANDLE) dlg->existing_text_delim_grp,
                              (HANDLE) dlg->erase_when_blank, 
                              NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->change_mrna, (HANDLE) dlg->apply_options, NULL);

  return (DialoG) p;
}

typedef struct featurefieldtable {
  FEATURE_FORM_BLOCK
  ValNodePtr header_line;
  Int4       num_columns;
  DialoG PNTR columns;
  ButtoN      accept_button;
} FeatureFieldTableData, PNTR FeatureFieldTablePtr;


static void CleanupFeatureFieldTableForm (GraphiC g, VoidPtr data)
{
  FeatureFieldTablePtr  form;

  form = (FeatureFieldTablePtr) data;
  if (form != NULL) {
    CleanUpTableData (form->header_line);
    form->columns = MemFree (form->columns);
  }
  StdCleanupFormProc (g, data);
}


typedef struct findgenelocustag {
  CharPtr locus_tag;
  ValNodePtr gene_list;
} FindGeneLocusTagData, PNTR FindGeneLocusTagPtr;

static void FindGeneByLocusTagBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  FindGeneLocusTagPtr p;
  SeqFeatPtr          gene;
  SeqMgrFeatContext   fcontext;

  if (bsp == NULL || userdata == NULL || !ISA_na (bsp->mol)) {
    return;
  }

  p = (FindGeneLocusTagPtr) userdata;

  gene = SeqMgrGetGeneByLocusTag (bsp, p->locus_tag, &fcontext);
  if (gene != NULL) {
    ValNodeAddPointer (&p->gene_list, OBJ_SEQFEAT, gene);
  }
}


static ValNodePtr 
FindMatchListForRow 
(FeatureFieldColumnChoicePtr PNTR f_list,
 Int4                             num_columns,
 SeqEntryPtr                      sep,
 TableLinePtr                     tlp,
 ValNodePtr PNTR                  missing_ids,
 ValNodePtr PNTR                  missing_genes)
{
  Int4       index;
  ValNodePtr match_list = NULL, vnp_match;
  SeqIdPtr   sip;
  BioseqPtr  bsp;
  FindGeneLocusTagData fd;

  if (f_list == NULL || tlp == NULL || tlp->parts == NULL || sep == NULL) return NULL;
  for (index = 0, vnp_match = tlp->parts;
       index < num_columns && vnp_match != NULL;
       index++, vnp_match = vnp_match->next) 
  {
    if (f_list[index]->match_type == 2) 
    {
      sip = CreateSeqIdFromText (vnp_match->data.ptrvalue, sep);
      bsp = BioseqFind (sip);
      sip = SeqIdFree (sip);
      if (bsp == NULL) 
      {
        if (missing_ids != NULL) 
        {
          ValNodeAddPointer (missing_ids, 0, StringSave (vnp_match->data.ptrvalue));
        }
      }
      else
      {
        ValNodeAddPointer (&match_list, OBJ_BIOSEQ, bsp);
      }
    }        
    else if (f_list[index]->match_type == 3)
    {
      fd.locus_tag = vnp_match->data.ptrvalue;
      fd.gene_list = NULL;
      VisitBioseqsInSep (sep, &fd, FindGeneByLocusTagBioseqCallback);
      if (fd.gene_list == NULL) 
      {
        if (missing_genes != NULL) 
        {
          ValNodeAddPointer (missing_genes, 0, StringSave (vnp_match->data.ptrvalue));
        }
      }
      else
      {
        ValNodeLink (&match_list, fd.gene_list);
        fd.gene_list = NULL;
      }
    }
  }
  return match_list;
}


static void DoLoadFeatureFieldTable (ButtoN b)
{
  FeatureFieldTablePtr  form;
  SeqEntryPtr           sep, oldscope;
  ValNodePtr            vnp;
  FeatureFieldColumnChoicePtr PNTR f_list;  
  Int4                             index;
  Int4                             num_rows;
  ValNodePtr PNTR                  match_lists;
  ValNodePtr                       missing_ids = NULL, missing_genes = NULL;
  CharPtr                          msg;
  Int4Ptr                          existing_text;
  Boolean                          has_mrna = FALSE;
  Boolean                          missing_mrna = FALSE;
  ValNodePtr                       errors = NULL;
  LogInfoPtr                       lip;
  Int4                             fields_affected = 0;
  CharPtr                          label;

  form = (FeatureFieldTablePtr) GetObjectExtra (b);
  if (form == NULL) return;
  
  sep = GetTopSeqEntryForEntityID (form->input_entityID);

  f_list = (FeatureFieldColumnChoicePtr PNTR) MemNew (form->num_columns * sizeof (FeatureFieldColumnChoicePtr));
  for (index = 0; index < form->num_columns; index++) {
    f_list[index] = DialogToPointer (form->columns[index]);
  }

  /* find match items (either protein IDs or genes) for each row */
  /* make sure all IDs and gene locus_tags in table are in record */
  num_rows = ValNodeLen (form->header_line);
  match_lists = (ValNodePtr PNTR) MemNew (num_rows * sizeof (ValNodePtr));
  for (index = 0, vnp = form->header_line; index < num_rows && vnp != NULL; index++, vnp = vnp->next) {
    match_lists[index] = FindMatchListForRow (f_list, form->num_columns, sep, vnp->data.ptrvalue, &missing_ids, &missing_genes);
  }

  if (missing_ids != NULL) {
    msg = CreateListMessage ("ID", 
  	                         missing_ids->next == NULL 
  	                         ? " is not found in this record." 
  	                         : " not found in this record.",
  	                         missing_ids);
    ValNodeAddPointer (&errors, 0, msg);
    missing_ids = ValNodeFreeData (missing_ids);
  }
  if (missing_genes != NULL) {
    msg = CreateListMessage ("Gene", 
  	                         missing_genes->next == NULL 
  	                         ? " is not found in this record." 
  	                         : " not found in this record.",
  	                         missing_genes);
    ValNodeAddPointer (&errors, 0, msg);
    missing_genes = ValNodeFreeData (missing_genes);
  }

  existing_text = (Int4Ptr) MemNew (form->num_columns * sizeof (Int4));

  for (vnp = form->header_line, index = 0; vnp != NULL; vnp = vnp->next, index++) {
    if (vnp->data.ptrvalue == NULL) continue;
    GetExistingValueCounts (vnp->data.ptrvalue, f_list, form->num_columns, match_lists[index], existing_text, &has_mrna, &missing_mrna, &errors);
  }
  
  lip = OpenLog ("Table Problems");
  for (vnp = errors; vnp != NULL; vnp = vnp->next) {
    fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    lip->data_in_log = TRUE;
  }
  errors = ValNodeFreeData (errors);
  
  if (has_mrna && missing_mrna) {
    fprintf (lip->fp, "Some coding regions for which product names will be set have overlapping mRNA features, but some do not!\n");
    lip->data_in_log = TRUE;
  }
  for (index = 0; index < form->num_columns; index++) {
    if (f_list[index]->field_choice != NULL) {
      if (existing_text[index] > 0) {
        label = GetCDSGeneProtFieldName (f_list[index]->field_choice);
        fprintf (lip->fp, "%d features affected by column %d contain existing text in field %s.\n", existing_text[index], index + 1, label);
        label = MemFree (label);
        lip->data_in_log = TRUE;
      }
    }
  }
  CloseLog (lip);
  if (lip->data_in_log) {
    if (ANS_CANCEL == Message (MSG_OKC, "Continue with errors?")) {
      existing_text = MemFree (existing_text);
      /* free match lists */
      for (index = 0; index < num_rows; index++) {
        match_lists[index] = ValNodeFree (match_lists[index]);
      }
      /* free choices */
      for (index = 0; index < form->num_columns; index++) {
        f_list[index] = FeatureFieldColumnChoiceFree (f_list[index]);
      }
      return;
    }
  }

  oldscope = SeqEntrySetScope (sep);

  for (vnp = form->header_line, index = 0; vnp != NULL; vnp = vnp->next, index++) {
    if (vnp->data.ptrvalue == NULL) continue;
    fields_affected += ApplyProteinIDTableLine (vnp->data.ptrvalue, f_list, form->num_columns, match_lists[index]);
  }

  SeqEntrySetScope (oldscope);

  /* free match lists */
  for (index = 0; index < num_rows; index++) {
    match_lists[index] = ValNodeFree (match_lists[index]);
  }
  /* free choices */
  for (index = 0; index < form->num_columns; index++) {
    f_list[index] = FeatureFieldColumnChoiceFree (f_list[index]);
  }

  ObjMgrSetDirtyFlag (form->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, form->input_entityID, 0, 0);
  Remove (form->form);
  Message (MSG_OK, "%d fields were affected", fields_affected);
}

static void SetFeatureFieldTableAccept (Pointer data)
{
  FeatureFieldTablePtr form;
  Int4                 index;
  FeatureFieldColumnChoicePtr f;
  Boolean                     have_match = FALSE, have_apply = FALSE;

  form = (FeatureFieldTablePtr) data;
  if (form == NULL) return;

  for (index = 0; index < form->num_columns && (!have_match || !have_apply); index++)
  {
    f = DialogToPointer (form->columns[index]);
    if (f != NULL) 
    {
      if (f->match_type > 1)
      {
        have_match = TRUE;
      }
      else if (f->field_choice != NULL)
      {
        have_apply = TRUE;
      }
      f = FeatureFieldColumnChoiceFree (f);
    }
  }
  if (have_match && have_apply) 
  {
    Enable (form->accept_button);
  } 
  else
  {
    Disable (form->accept_button);
  }
}


extern void LoadFeatureFieldTable (IteM i)
{
  BaseFormPtr          bfp;
  SeqEntryPtr          sep;
  ValNodePtr           header_line;
  Int4                 index;
  FeatureFieldTablePtr form;
  WindoW               w;
  GrouP                h, g, c;
  TableLinePtr         tlp;
  ValNodePtr           vnp;
  Int4Ptr              blank_list = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);


  header_line = ReadTableData ();
  if (header_line == NULL || header_line->data.ptrvalue == NULL) {
    CleanUpTableData (header_line);
    return;
  }

  /* use form to pick columns for IDs, values */

  form = (FeatureFieldTablePtr) MemNew (sizeof (FeatureFieldTableData));
  if (form == NULL) return;
  form->input_entityID = bfp->input_entityID;
  form->header_line = header_line;
  tlp = (TableLinePtr) header_line->data.ptrvalue;
  form->num_columns = ValNodeLen (tlp->parts);
  form->columns = (DialoG PNTR) MemNew (form->num_columns * sizeof (DialoG));

  /* now create a dialog to display values */
  w = FixedWindow (-50, -33, -10, -10, "Feature Field Table", StdCloseWindowProc);
  SetObjectExtra (w, form, CleanupFeatureFieldTableForm);
  form->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  g = HiddenGroup (h, 3, 0, NULL);

  /* pre-analyze table for blanks, so we can display the information in the column dialogs */
  blank_list = GetColumnBlankCounts (header_line, form->num_columns);

  for (vnp = tlp->parts, index = 0; vnp != NULL; vnp = vnp->next, index++)
  {
    form->columns[index] = FeatureFieldColumnChoiceDialog (g, vnp->data.ptrvalue, 
                                                           blank_list == NULL ? 0 : blank_list[index], 
                                                           SetFeatureFieldTableAccept, form);
  }
  blank_list = MemFree (blank_list);

  c = HiddenGroup (h, 4, 0, NULL);
  form->accept_button = DefaultButton (c, "Accept", DoLoadFeatureFieldTable);
  SetObjectExtra (form->accept_button, form, NULL);
  Disable (form->accept_button);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

static Pointer GetRevSectionSequence (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BioseqPtr bsp;
  Char        id_str[45];

  if (data == NULL)
  {
    return NULL;
  } 
  else
  {
    bsp = (BioseqPtr) data;
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, 39);
    return StringSave (id_str);
  }
}


static Pointer GetRevSectionInterval (Uint1 data_choice, Pointer data, Pointer metadata)
{
  BioseqPtr bsp;
  DeltaSeqPtr dsp;
  Uint1       seg_num = 0;
  Int4        pos = 0;
  SeqLocPtr   loc;
  SeqLitPtr   slip;
  Char        buf[50];

  if ((bsp = (BioseqPtr)data) == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext == NULL)
  {
    return NULL;
  } 

  for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
    switch (dsp->choice)
    {
      case 1:      /* SeqLocPtr */
        if ((loc = (SeqLocPtr)dsp->data.ptrvalue) != NULL) {
          if (loc->choice != SEQLOC_NULL) {
            seg_num++;
            pos += SeqLocLen (loc);
          }
        }
        break;
      case 2:   /* SeqLitPtr */
        slip = (SeqLitPtr)(dsp->data.ptrvalue);
        if (slip != NULL) {
          if (slip->seq_data != NULL) {
            if (seg_num == data_choice) {
              sprintf (buf, "%d-%d", pos + 1, pos + 1 + slip->length);
              return StringSave (buf);
            } else {
              seg_num++;
            }
            pos += slip->length;
          }
        }
        break;
    }
  }
  return NULL;
}


/*
static BulkEdFieldData revsection_fields[] = {
  { "Sequence", NULL, NULL, GetRevSectionSequence, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { "Interval", NULL, NULL, GetRevSectionInterval, BulkDisplaySimpleText, BulkFreeSimpleText, NULL, BulkFormatSimpleText, NULL, NULL, BulkSimpleTextCopy },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};
*/


typedef struct deltaseqint {
  BioseqPtr bsp;
  Int4      start;
  Int4      stop;
  SeqLitPtr slip;
} DeltaSeqIntData, PNTR DeltaSeqIntPtr;


static DeltaSeqIntPtr DeltaSeqIntNew (BioseqPtr bsp, Int4 start, Int4 stop, SeqLitPtr slip)
{
  DeltaSeqIntPtr dsip;

  dsip = (DeltaSeqIntPtr) MemNew (sizeof (DeltaSeqIntData));
  dsip->bsp = bsp;
  dsip->start = start;
  dsip->stop = stop;
  dsip->slip = slip;
  return dsip;
}


static void ListDeltaSeqIntervalsCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR interval_list;
  DeltaSeqPtr dsp;
  Uint1       seg_num = 0;
  Int4        pos = 0;
  SeqLocPtr   loc;
  SeqLitPtr   slip;
  Char        id_str[45];
  ClickableItemPtr cip;
  
  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext == NULL || (interval_list = (ValNodePtr PNTR)userdata) == NULL) {
    return;
  }

  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);

  for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
    switch (dsp->choice)
    {
      case 1:      /* SeqLocPtr */
        if ((loc = (SeqLocPtr)dsp->data.ptrvalue) != NULL) {
          if (loc->choice != SEQLOC_NULL) {
            seg_num++;
            pos += SeqLocLen (loc);
          }
        }
        break;
      case 2:   /* SeqLitPtr */
        slip = (SeqLitPtr)(dsp->data.ptrvalue);
        if (slip != NULL) {
          if (slip->seq_data != NULL) {
            cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
            cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (id_str) + 65));
            sprintf (cip->description, "%s:%d-%d", id_str, pos + 1, pos + slip->length);
            cip->clickable_item_type = seg_num;
            ValNodeAddPointer (&cip->item_list, 0, DeltaSeqIntNew (bsp, pos, pos + slip->length - 1, slip));
            ValNodeAddPointer (interval_list, 0, cip);
            seg_num++;
          }
          pos += slip->length;
        }
        break;
    }
  }
}


static Int4 GetHighestInterval (ValNodePtr int_list)
{
  Int4 max = 0;
  ClickableItemPtr cip;

  while (int_list != NULL) {
    cip = (ClickableItemPtr) int_list->data.ptrvalue;
    if (cip != NULL && cip->clickable_item_type > max) {
      max = cip->clickable_item_type;
    }
    int_list = int_list->next;
  }
  return max;
}


typedef struct revseqintform {
  FORM_MESSAGE_BLOCK
  DialoG          clickable_list; 
  PopuP           interval_choice;
  ButtoN          also_reverse_feats;
  ValNodePtr      delta_list;
} RevSeqIntFormData, PNTR RevSeqIntFormPtr;

static void CleanupRevSeqIntForm (GraphiC g, VoidPtr data)
{
  RevSeqIntFormPtr f;
  ValNodePtr       vnp;
  ClickableItemPtr cip;

  f = (RevSeqIntFormPtr) data;
  if (f != NULL) {
    for (vnp = f->delta_list; vnp != NULL; vnp = vnp->next) {
      /* note - the rest of the ClickableItem will be freed by the clickableItemlist */
      cip = (ClickableItemPtr) vnp->data.ptrvalue;
      if (cip != NULL && cip->item_list != NULL) {
        cip->item_list->data.ptrvalue = MemFree (cip->item_list->data.ptrvalue);
      }
    }
  }
  StdCleanupFormProc (g, data);
}


static void CheckIntervals (ButtoN b)
{
  RevSeqIntFormPtr f;
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  Int4             val;

  f = (RevSeqIntFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  val = GetValue (f->interval_choice);
  for (vnp = f->delta_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL && cip->clickable_item_type == val - 1) {
      cip->chosen = TRUE;
    }
  }
  PointerToDialog (f->clickable_list, f->delta_list);
}


static Int4 RevCompCoordInInterval (Int4 coord, Int4 start, Int4 stop)
{
  Int4 offset;

  if (coord < start || coord > stop) return coord;

  offset = coord - start;
  coord = stop - offset;
  return coord;
}


static void RevCompIntFuzz (IntFuzzPtr ifp, Int4 start, Int4 stop)
{
  Int4 tmp;

  if (ifp == NULL) return;

	switch (ifp->choice)
	{
		case 1:      /* plus/minus - no changes */
		case 3:      /* percent - no changes */
			break;
		case 2:      /* range */
      tmp = RevCompCoordInInterval (ifp->a, start, stop);
      ifp->a = RevCompCoordInInterval (ifp->b, start, stop);
      ifp->b = tmp;
			break;
		case 4:     /* lim */
			switch (ifp->a)
			{
				case 1:    /* greater than */
					ifp->a = 2;
					break;
				case 2:    /* less than */
					ifp->a = 1;
					break;
				case 3:    /* to right of residue */
					ifp->a = 4;
					break;
				case 4:    /* to left of residue */
					ifp->a = 3;
					break;
				default:
					break;
			}
			break;
	}
}

static Uint1 RevStrand (Uint1 strand)
{
  if (strand == Seq_strand_minus) {
    strand = Seq_strand_plus;
  } else {
    strand = Seq_strand_minus;
  }
  return strand;
}


static void RevCompSeqPntForFlippedInterval (SeqPntPtr spp, BioseqPtr bsp, Int4 start, Int4 stop)
{

  if (spp == NULL || bsp == NULL) return;
  if (!SeqIdIn (spp->id, bsp->id)) return;
  if (spp->point < start || spp->point > stop) return;

  /* flip */
  spp->point = RevCompCoordInInterval (spp->point, start, stop);
  spp->strand = RevStrand (spp->strand);
  RevCompIntFuzz (spp->fuzz, start, stop);
}


static void RevCompLocationForFlippedInterval (SeqLocPtr head, BioseqPtr bsp, Int4 start, Int4 stop)
{
	SeqLocPtr slp;
	SeqIntPtr sip;
	SeqPntPtr spp;
	PackSeqPntPtr pspp, pspp2;
	SeqBondPtr sbp;
	SeqIdPtr oldids;
	Int4 numpnt, i, tpos, tmp;
	Boolean do_rev;
	IntFuzzPtr ifp;
	
	if ((head == NULL) || (bsp == NULL)) return;

	oldids = bsp->id;
	switch (head->choice)
	{    
    case SEQLOC_BOND:   /* bond -- 2 seqs */
			sbp = (SeqBondPtr)(head->data.ptrvalue);
      RevCompSeqPntForFlippedInterval (sbp->a, bsp, start, stop);
      RevCompSeqPntForFlippedInterval (sbp->b, bsp, start, stop);
			break;
    case SEQLOC_FEAT:   /* feat -- can't track yet */
    case SEQLOC_NULL:    /* NULL */
    case SEQLOC_EMPTY:    /* empty */
			break;
    case SEQLOC_WHOLE:    /* whole */
      /* nothing to reverse */
			break;
    case SEQLOC_EQUIV:    /* does it stay equiv? */
    case SEQLOC_MIX:    /* mix -- more than one seq */
    case SEQLOC_PACKED_INT:    /* packed int */
			for (slp = (SeqLocPtr)(head->data.ptrvalue); slp != NULL; slp = slp->next)
			{
        RevCompLocationForFlippedInterval (slp, bsp, start, stop);
      }
      break;
    case SEQLOC_INT:    /* int */
			sip = (SeqIntPtr)(head->data.ptrvalue);
			if (SeqIdIn(sip->id, oldids))
			{
        if (sip->from < start || sip->to > stop) {
          /* not contained in interval */
          break;
        }
        tmp = RevCompCoordInInterval (sip->from, start, stop);
        sip->from = RevCompCoordInInterval (sip->to, start, stop);
        sip->to = tmp;
        sip->strand = RevStrand (sip->strand);
        RevCompIntFuzz (sip->if_to, start, stop);
        RevCompIntFuzz (sip->if_from, start, stop);
        ifp = sip->if_to;
        sip->if_to = sip->if_from;
        sip->if_from = ifp;
      }
      break;
    case SEQLOC_PNT:    /* pnt */
			spp = (SeqPntPtr)(head->data.ptrvalue);
      RevCompSeqPntForFlippedInterval (spp, bsp, start, stop);
      break;
    case SEQLOC_PACKED_PNT:    /* packed pnt */
			pspp = (PackSeqPntPtr)(head->data.ptrvalue);
      do_rev = FALSE;
			if (SeqIdIn(pspp->id, oldids))
			{
        do_rev = TRUE;
				numpnt = PackSeqPntNum(pspp);
        /* can only reverse set if all are in interval */
				for (i = 0; i < numpnt && do_rev; i++)
				{
					tpos = PackSeqPntGet(pspp, i);
          if (tpos < start || tpos > stop) {
            do_rev = FALSE;
          }
        }
        if (do_rev) {
				  pspp2 = PackSeqPntNew();
  				pspp2->id = pspp->id;
          pspp2->fuzz = pspp->fuzz;
          pspp->fuzz = NULL;
          RevCompIntFuzz (pspp2->fuzz, start, stop);
				  for (i = 0; i < numpnt && do_rev; i++)
				  {
					  tpos = PackSeqPntGet(pspp, i);
            tpos = RevCompCoordInInterval (tpos, start, stop);
            PackSeqPntPut (pspp2, tpos);
          }
          pspp2->strand = RevStrand (pspp->strand);
          pspp = PackSeqPntFree (pspp);
          head->data.ptrvalue = pspp2;
			  }
      }
      break;
    default:
      break;
	}
}


static void RevCompOneFeatForBioseqInterval (SeqFeatPtr sfp, BioseqPtr bsp, Int4 start, Int4 stop)
{
  CodeBreakPtr cbp;
  CdRegionPtr  crp;
  RnaRefPtr    rrp;
  tRNAPtr      trp;

  if (sfp == NULL || bsp == NULL) return;

  RevCompLocationForFlippedInterval (sfp->location, bsp, start, stop);
  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          RevCompLocationForFlippedInterval (cbp->loc, bsp, start, stop);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          RevCompLocationForFlippedInterval (trp->anticodon, bsp, start, stop);
        }
      }
      break;
    default :
      break;
  }
}


static void FlipSelectedSequenceIntervals (ButtoN b)
{
  RevSeqIntFormPtr f;
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  Int4             num_disc;
  DeltaSeqIntPtr   dsip;
  Boolean          rev_feats;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;

  f = (RevSeqIntFormPtr) GetObjectExtra (b);
  if (f == NULL) return;
  num_disc = CountChosenDiscrepancies (f->delta_list, FALSE);
  if (num_disc == 0) {
    Message (MSG_ERROR, "No intervals selected");
    return;
  }

  rev_feats = GetStatus (f->also_reverse_feats);

  for (vnp = f->delta_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL && cip->chosen && cip->item_list != NULL) {
      dsip = (DeltaSeqIntPtr) cip->item_list->data.ptrvalue;
      if (dsip != NULL && dsip->slip != NULL) {
        ReverseSeqData (dsip->slip->seq_data_type, dsip->slip->length, dsip->slip->seq_data);
        ComplementSeqData (dsip->slip->seq_data_type, dsip->slip->length, dsip->slip->seq_data);
        if (rev_feats) {
          for (sfp = SeqMgrGetNextFeature (dsip->bsp, NULL, 0, 0, &fcontext);
               sfp != NULL && fcontext.left <= dsip->stop;
               sfp = SeqMgrGetNextFeature (dsip->bsp, sfp, 0, 0, &fcontext)) {
            if (fcontext.left >= dsip->start && fcontext.right <= dsip->stop) {
              RevCompOneFeatForBioseqInterval (sfp, dsip->bsp, dsip->start, dsip->stop);
            }
          }
        }
      }
    }
  }
  Update ();
  ObjMgrSetDirtyFlag (f->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, f->input_entityID, 0, 0);
  Remove (f->form);
}

NLM_EXTERN void FlipSequenceIntervals (IteM i)
{
  BaseFormPtr         bfp;
  RevSeqIntFormPtr    f;
  SeqEntryPtr         sep;
  WindoW              w;
  GrouP               h, g, c;
  ButtoN              b;
  ValNodePtr          delta_list = NULL;
  Int4                max, n;
  Char                buf[15];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitBioseqsInSep (sep, &delta_list, ListDeltaSeqIntervalsCallback);

  if (delta_list == NULL)
  {
    Message (MSG_OK, "No sequence intervals");
    return;
  }


  f = (RevSeqIntFormPtr) MemNew (sizeof (RevSeqIntFormData));
  if (f == NULL)
  {
    return;
  }
  f->delta_list = delta_list;
  
  f->input_entityID = bfp->input_entityID;
  w = FixedWindow (-50, -33, -10, -10, "Reverse Sequence Intervals", StdCloseWindowProc);
  SetObjectExtra (w, f, CleanupRevSeqIntForm);
  f->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  f->clickable_list = CreateClickableListDialog (h, "Intervals", "",
                                                 NULL, NULL, NULL, NULL);
  PointerToDialog (f->clickable_list, f->delta_list);

  max = GetHighestInterval (f->delta_list);
  
  g = HiddenGroup (h, 4, 0, NULL);
  StaticPrompt (g, "Check Interval ", 0, popupMenuHeight, programFont, 'l');
  f->interval_choice = PopupList (g, TRUE, NULL);
  for (n = 0; n <= max; n++) {
    sprintf (buf, "%d", n + 1);
    PopupItem (f->interval_choice, buf);
  }
  SetValue (f->interval_choice, 1);
  StaticPrompt (g, " for every sequence", 0, popupMenuHeight, programFont, 'l');
  b = PushButton (g, "Check", CheckIntervals);
  SetObjectExtra (b, f, NULL);

  f->also_reverse_feats = CheckBox (h, "Also reverse features", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);

  b = PushButton (c, "Flip Checked Intervals", FlipSelectedSequenceIntervals);
  SetObjectExtra (b, f, NULL);
    
  PushButton (c, "Dismiss", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) f->clickable_list, (HANDLE) g, (HANDLE) f->also_reverse_feats, (HANDLE) c, NULL);

  RealizeWindow (w);
  
  Show (w);
}


static void MacroAECRAction (IteM i, Boolean indexer_version, Uint1 action_type, Uint1 qual_type)
{
  BaseFormPtr        bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  SingleAECRMacroAction (bfp->input_entityID, indexer_version, action_type, qual_type);
}

extern void MacroApplyGBQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_apply, FieldType_feature_field);
}


extern void MacroApplySourceQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_apply, FieldType_source_qual);
}


extern void MacroApplyCDSGeneProt (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_apply, FieldType_cds_gene_prot);
}

extern void PublicMacroApplyCDSGeneProt (IteM i)
{
  MacroAECRAction (i, FALSE, ActionChoice_apply, FieldType_cds_gene_prot);
}


extern void MacroApplyRNAQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_apply, FieldType_rna_field);
}


extern void MacroRemoveGBQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_remove, FieldType_feature_field);
}


extern void MacroRemoveSourceQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_remove, FieldType_source_qual);
}


extern void MacroRemoveCDSGeneProt (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_remove, FieldType_cds_gene_prot);
}


extern void MacroRemoveRNAQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_remove, FieldType_rna_field);
}


extern void MacroConvertGBQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_convert, FieldType_feature_field);
}


extern void MacroConvertSourceQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_convert, FieldType_source_qual);
}


extern void MacroConvertCDSGeneProt (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_convert, FieldType_cds_gene_prot);
}


extern void MacroConvertRNAQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_convert, FieldType_rna_field);
}


extern void MacroSwapGBQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_swap, FieldType_feature_field);
}


extern void MacroSwapSourceQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_swap, FieldType_source_qual);
}


extern void MacroSwapCDSGeneProt (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_swap, FieldType_cds_gene_prot);
}


extern void MacroSwapRNAQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_swap, FieldType_rna_field);
}


extern void MacroEditGBQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_edit, FieldType_feature_field);
}


extern void MacroEditSourceQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_edit, FieldType_source_qual);
}


extern void MacroEditCDSGeneProt (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_edit, FieldType_cds_gene_prot);
}


extern void PublicMacroEditCDSGeneProt (IteM i)
{
  MacroAECRAction (i, FALSE, ActionChoice_edit, FieldType_cds_gene_prot);
}


extern void MacroEditRNAQual (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_edit, FieldType_rna_field);
}


extern void MacroApplyStructuredComment (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_apply, FieldType_struc_comment_field);
}

extern void PublicMacroApplyStructuredComment (IteM i)
{
  MacroAECRAction (i, FALSE, ActionChoice_apply, FieldType_struc_comment_field);
}

extern void MacroEditStructuredComment (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_edit, FieldType_struc_comment_field);
}

extern void PublicMacroEditStructuredComment (IteM i)
{
  MacroAECRAction (i, FALSE, ActionChoice_edit, FieldType_struc_comment_field);
}

extern void MacroRemoveStructuredComment (IteM i)
{
  MacroAECRAction (i, TRUE, ActionChoice_remove, FieldType_struc_comment_field);
}


NLM_EXTERN void RemoveAllSourceNoteBaseForm (BaseFormPtr bfp)
{
  SeqEntryPtr sep;
  RemoveActionPtr remove;
  AECRActionPtr aecr;
  ValNodePtr    aecr_a;
  ValNodePtr    macro;
  ValNodePtr    src_field;

  if (bfp == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) {
    return;
  }

  src_field = ValNodeNew (NULL);
  src_field->choice = SourceQualChoice_textqual;
  src_field->data.intvalue = Source_qual_all_notes;
  remove = RemoveActionNew ();
  remove->field = ValNodeNew (NULL);
  remove->field->choice = FieldType_source_qual;
  remove->field->data.ptrvalue = src_field;
  aecr_a = ValNodeNew (NULL);
  aecr_a->choice = ActionChoice_remove;
  aecr_a->data.ptrvalue = remove;
  aecr = AECRActionNew ();
  aecr->action = aecr_a;
  macro = ValNodeNew (NULL);
  macro->choice = MacroActionChoice_aecr;
  macro->data.ptrvalue = aecr;

  ApplyMacroToSeqEntryEx (sep, macro, NULL, Sequin_GlobalAlign2Seq);

  macro = MacroActionChoiceFree(macro);
}


NLM_EXTERN void ApplyIsolationSourceBaseForm (BaseFormPtr bfp)
{
  ForM           f;
  ValNodePtr     macro;
  ValNodePtr     src_field;
  ApplyActionPtr apply;
  AECRActionPtr  aecr;
  ValNodePtr     aecr_a;
  
  f = SingleMacroAction (bfp->input_entityID, TRUE);

  src_field = ValNodeNew (NULL);
  src_field->choice = SourceQualChoice_textqual;
  src_field->data.intvalue = Source_qual_isolation_source;

  apply = ApplyActionNew ();
  apply->field = ValNodeNew (NULL);
  apply->field->choice = FieldType_source_qual;
  apply->field->data.ptrvalue = src_field;

  aecr_a = ValNodeNew (NULL);
  aecr_a->choice = ActionChoice_apply;
  aecr_a->data.ptrvalue = apply;
  aecr = AECRActionNew ();
  aecr->action = aecr_a;

  macro = ValNodeNew (NULL);
  macro->choice = MacroActionChoice_aecr;
  macro->data.ptrvalue = aecr;

  PointerToForm (f, macro);
  macro = MacroActionChoiceFree (macro);
    
  Show (f);
  Select (f);
  SendMessageToForm (f, VIB_MSG_ENTER);
}


typedef struct convertboisourcedbxreftofeaturedbxref {
  FORM_MESSAGE_BLOCK
  DialoG feature_type;
} ConvertBioSourceDbxrefToFeatureDbxrefData, PNTR ConvertBioSourceDbxrefToFeatureDbxrefPtr;


static ValNodePtr DbxrefListFree (ValNodePtr vnp)

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


static ValNodePtr DbxrefListCopy (ValNodePtr vnp)

{
  ValNodePtr copy = NULL;
  
  while (vnp != NULL) {
    ValNodeAddPointer (&copy, 0, AsnIoMemCopy (vnp->data.ptrvalue, (AsnReadFunc) DbtagAsnRead, (AsnWriteFunc) DbtagAsnWrite));
    vnp = vnp->next;
  }
  return copy;
}


static ValNodePtr ExtractNonTaxonDbxrefs (ValNodePtr PNTR list)
{
  ValNodePtr prev = NULL, extracted = NULL, vnp, vnp_next;
  DbtagPtr dbtag;

  if (list == NULL || *list == NULL) return NULL;
  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    dbtag = (DbtagPtr) vnp->data.ptrvalue;
    if (dbtag != NULL && StringCmp (dbtag->db, "taxon") != 0) {
      if (prev == NULL) {
        *list = vnp_next;
      } else {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      ValNodeLink (&extracted, vnp);
    } else {
      prev = vnp;
    }
    vnp = vnp_next;
  }
  return extracted;
}


static void ConvertBioSourceDbxrefToFeatureDbxrefBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr        vnp;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  BioSourcePtr      biop;
  ValNodePtr        dbxref_list = NULL;
  Int4              featdef;
  Boolean           any_feat = FALSE;

  if (bsp == NULL || userdata == NULL) {
    return;
  }

  vnp = (ValNodePtr) userdata;
  featdef = GetFeatdefFromFeatureType (vnp->choice);

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp == NULL || sdp->data.ptrvalue == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop->org == NULL || biop->org->db == NULL) return;
  dbxref_list = ExtractNonTaxonDbxrefs (&(biop->org->db));
  if (dbxref_list == NULL) return;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext)) {
    ValNodeLink (&(sfp->dbxref), DbxrefListCopy(dbxref_list));
    any_feat = TRUE;
  }
  if (any_feat) {
    dbxref_list = DbxrefListFree(dbxref_list);
  } else {
    ValNodeLink (&(biop->org->db), dbxref_list);
  }
}


static void AcceptConvertBioSourceDbxrefToFeatureDbxref (ButtoN b)
{
  ConvertBioSourceDbxrefToFeatureDbxrefPtr dlg;
  ValNodePtr vnp;
  SeqEntryPtr sep;

  dlg = (ConvertBioSourceDbxrefToFeatureDbxrefPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  vnp = DialogToPointer (dlg->feature_type);
  if (vnp == NULL) {
    Message (MSG_ERROR, "Please choose feature type");
    return;
  }

  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);

  VisitBioseqsInSep (sep, vnp, ConvertBioSourceDbxrefToFeatureDbxrefBioseqCallback);
  vnp = ValNodeFree (vnp);

  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  Remove (dlg->form);
}
  
  
extern void ConvertBioSourceDbxrefToFeatureDbxref (IteM i)
{
  BaseFormPtr         bfp;
  ConvertBioSourceDbxrefToFeatureDbxrefPtr dlg;
  WindoW              w;
  GrouP               h, c;
  ButtoN              b;
  PrompT              ppt;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  dlg = (ConvertBioSourceDbxrefToFeatureDbxrefPtr) MemNew (sizeof (ConvertBioSourceDbxrefToFeatureDbxrefData));

  w = FixedWindow (-50, -33, -10, -10, "Convert BioSource Dbxrefs to Feature Dbxrefs", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupExtraProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt = StaticPrompt (h, "Choose feature type to receive BioSource dbxrefs", 0, dialogTextHeight, programFont, 'l');
  dlg->feature_type = FeatureTypeDialog (h, NULL, NULL);
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", AcceptConvertBioSourceDbxrefToFeatureDbxref);
  SetObjectExtra (b, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->feature_type, (HANDLE) c, NULL);
  Show (w);
  Select (w);
}


typedef struct structuredcommentsform {
  FORM_MESSAGE_BLOCK
  GrouP table_type;
  GrouP match_info;
  PopuP match_column;
  DialoG match_type;
  DialoG database_name;

  ValNodePtr table;
} StructuredCommentsFormData, PNTR StructuredCommentsFormPtr;

static void CleanupStructuredCommentsForm (GraphiC g, VoidPtr data)

{
  StructuredCommentsFormPtr  frm;

  frm = (StructuredCommentsFormPtr) data;
  if (frm != NULL) {
    frm->table = FreeTabTable(frm->table);
  }
  StdCleanupFormProc (g, data);
}


static UserObjectPtr UserObjectFromRow (ValNodePtr header, ValNodePtr line, Int4 col, CharPtr dbname, ValNodePtr PNTR err_list)
{
  ValNodePtr vnp_h, vnp_l;
  CharPtr    id_str = NULL;
  Int4       num;
  CharPtr    extra_data_fmt = "Too many fields in line for %s";
  CharPtr    msg;
  UserObjectPtr uop;
  CharPtr       prefix = NULL, suffix = NULL;
  CharPtr       prefix_fmt = "##%s-START##";
  CharPtr       suffix_fmt = "##%s-END##";
  CharPtr       tmp_dbname = NULL;
  Boolean       prefix_from_row = FALSE;

  if (header == NULL || line == NULL) {
    return NULL;
  }

  /* find id_str */
  vnp_l = line->data.ptrvalue;
  if (vnp_l == NULL) {
    return NULL;
  }
  num = 0;
  while (vnp_l != NULL && num < col) {
    num++;
    vnp_l = vnp_l->next;
  }
  if (vnp_l == NULL || StringHasNoText (vnp_l->data.ptrvalue)) {
    return NULL;
  } else {
    id_str = vnp_l->data.ptrvalue;
  }

  
  /* look for dbname in table if not supplied */
  if (StringHasNoText (dbname)) {
    for (vnp_h = header, vnp_l = line->data.ptrvalue;
         vnp_h != NULL && vnp_l != NULL && StringHasNoText (dbname);
         vnp_h = vnp_h->next, vnp_l = vnp_l->next) {
      if (StringICmp (vnp_h->data.ptrvalue, "StructuredCommentPrefix") == 0
          || StringICmp (vnp_h->data.ptrvalue, "StructuredCommentSuffix") == 0) {
        tmp_dbname = StructuredCommentDbnameFromString(vnp_l->data.ptrvalue);
        dbname = tmp_dbname;
        prefix_from_row = TRUE;
      }
    }
  } else {
    /* fix supplied dbname if necessary */
    tmp_dbname = StructuredCommentDbnameFromString(dbname);
    if (tmp_dbname != NULL) {
      dbname = tmp_dbname;
    }
  }

  /* build user object */
  vnp_h = header;
  vnp_l = line->data.ptrvalue;
  num = 0;

  if (!StringHasNoText (dbname)) {
    prefix = (CharPtr) MemNew (sizeof (Char) * (StringLen (prefix_fmt) + StringLen (dbname)));
    sprintf (prefix, prefix_fmt, dbname);
    suffix = (CharPtr) MemNew (sizeof (Char) * (StringLen (suffix_fmt) + StringLen (dbname)));
    sprintf (suffix, suffix_fmt, dbname);
  }

  uop = CreateStructuredCommentUserObject (prefix, suffix);

  prefix = MemFree (prefix);
  suffix = MemFree (suffix);

  while (vnp_h != NULL && vnp_l != NULL) {
    if (num != col && !StringHasNoText (vnp_l->data.ptrvalue)) {
      if (!prefix_from_row 
          || (StringCmp (vnp_h->data.ptrvalue, "StructuredCommentPrefix") != 0
              && StringCmp (vnp_h->data.ptrvalue, "StructuredCommentSuffix") != 0)) {
        AddItemStructuredCommentUserObject (uop, vnp_h->data.ptrvalue, vnp_l->data.ptrvalue);
      }
    }
    vnp_h = vnp_h->next;
    vnp_l = vnp_l->next;
    num++;
  }
  while (vnp_l != NULL && StringHasNoText (vnp_l->data.ptrvalue)) {
    vnp_l = vnp_l->next;
  }
  if (vnp_l != NULL) {
    msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (extra_data_fmt) + StringLen (id_str)));
    sprintf (msg, extra_data_fmt, id_str);
    ValNodeAddPointer (err_list, 0, msg);
  }
  tmp_dbname = MemFree (tmp_dbname);
  return uop;
}


static void ChangeStructuredCommentTableType(GrouP g)
{
  StructuredCommentsFormPtr frm;

  frm = (StructuredCommentsFormPtr)GetObjectExtra (g);
  if (frm == NULL) {
    return;
  }

  if (GetValue (frm->table_type) == 1) {
    Enable (frm->match_info);
  } else {
    Disable (frm->match_info);
  }
}


static void DoParseStructuredComments (ButtoN b)
{
  StructuredCommentsFormPtr frm;
  Int4                      col;
  MatchTypePtr              match_type;
  ValNodePtr                sequence_lists;
  SeqEntryPtr               sep;
  ValNodePtr                err_list = NULL, vnp, header, line, s_row, vnp_s, tmp_table, vnp_h;
  LogInfoPtr                lip;
  BioseqPtr                 bsp;
  UserObjectPtr             uop, uop_cpy;
  SeqDescrPtr               sdp;
  CharPtr                   database_name = NULL, spec_db, file_db;
  Boolean                   found_prefix;

  frm = (StructuredCommentsFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  database_name = DialogToPointer (frm->database_name);

  if (GetValue (frm->table_type) == 1) {
    col = GetValue (frm->match_column) - 1;
    match_type = DialogToPointer (frm->match_type);
    sequence_lists = GetSequenceListsForMatchTypeInTabTable (sep, frm->table->next, col, match_type, &err_list);

    if (err_list != NULL) {
      lip = OpenLog ("Table Problems");
      for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
        fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
      }
      lip->data_in_log = TRUE;
      CloseLog (lip);
      lip = FreeLog (lip);
      err_list = ValNodeFreeData (err_list);
      if (ANS_YES != Message (MSG_YN, "Continue with table problems")) {
        sequence_lists = FreeSequenceLists(sequence_lists);
        return;
      }
    }

    WatchCursor ();
    Update();

    header = frm->table->data.ptrvalue;
    line = frm->table->next;
    s_row = sequence_lists;


    while (line != NULL && s_row != NULL) {
      vnp_s = s_row->data.ptrvalue;
      if (vnp_s != NULL) {
        uop = UserObjectFromRow (header, line, col, database_name, &err_list);
        while (vnp_s != NULL) {
          bsp = vnp_s->data.ptrvalue;
          uop_cpy = (UserObjectPtr) AsnIoMemCopy (uop, (AsnReadFunc) UserObjectAsnRead, (AsnWriteFunc) UserObjectAsnWrite);
          sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_user);
          sdp->data.ptrvalue = uop_cpy;
          vnp_s = vnp_s->next;
        }
        uop = UserObjectFree (uop);
      }
      line = line->next;
      s_row = s_row->next;
    }


    sequence_lists = FreeSequenceLists(sequence_lists);

  } else {
    tmp_table = FlipTabTableAxes (frm->table);
    if (!StringHasNoText (database_name)) {
      found_prefix = FALSE;
      spec_db = StructuredCommentDbnameFromString(database_name);
      for (vnp_h = tmp_table->data.ptrvalue, vnp_s = tmp_table->next->data.ptrvalue; 
           vnp_h != NULL && vnp_s != NULL; 
           vnp_h = vnp_h->next, vnp_s = vnp_s->next) {
        if (StringICmp(vnp_h->data.ptrvalue, "StructuredCommentPrefix") == 0
            || StringICmp (vnp_h->data.ptrvalue, "StructuredCommentSuffix") == 0) {
          file_db = StructuredCommentDbnameFromString(vnp_s->data.ptrvalue);
          if (StringCmp (file_db, spec_db) != 0) {
            vnp_s->data.ptrvalue = MemFree (vnp_s->data.ptrvalue);
            vnp_s->data.ptrvalue = StringSave (spec_db);
          }
          file_db = MemFree (file_db);
          found_prefix = TRUE;
        }
      }
      if (!found_prefix) {
        vnp = ValNodeNew (NULL);
        vnp->data.ptrvalue = StringSave ("StructuredCommentPrefix");
        vnp->next = tmp_table->data.ptrvalue;
        tmp_table->data.ptrvalue = vnp;

        vnp = ValNodeNew (NULL);
        vnp->data.ptrvalue = StringSave (spec_db);
        vnp->next = tmp_table->next->data.ptrvalue;
        tmp_table->next->data.ptrvalue = vnp;
      }
      spec_db = MemFree (spec_db);
    }

    CreateStructuredCommentsForAllFromTable(sep, tmp_table->data.ptrvalue, tmp_table->next, &err_list);
    tmp_table = FreeTabTable(tmp_table);
  }

  /* show errors from creating structured comments */
  if (err_list != NULL) {
    lip = OpenLog ("Data Problems");
    for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
      fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }
    lip->data_in_log = TRUE;
    CloseLog (lip);
    lip = FreeLog (lip);
    err_list = ValNodeFreeData (err_list);
  }

  database_name = MemFree (database_name);

  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);

  Remove (frm->form);
  ArrowCursor ();
  Update ();
  
}


NLM_EXTERN void CreateStructuredCommentsBaseForm (BaseFormPtr bfp)
{
  SeqEntryPtr        sep;
  Char               path [PATH_MAX];
  ValNodePtr         vnp;
  FILE *fp;
  StructuredCommentsFormPtr frm;
  WindoW                    w;
  GrouP                     h, c;
  PrompT                    ppt1, ppt2, ppt3;
  CharPtr                   choice_fmt = "Column %d (%s)";
  CharPtr                   choice_str, val;
  Int4                      col;
  ButtoN                    b;

  if (bfp == NULL)
    return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;

  if (GetInputFileName (path, sizeof (path), "", "TEXT")) {
    fp = FileOpen (path, "r");
    if (fp == NULL) {
      Message (MSG_ERROR, "Unable to open %s");
    } else {
      frm = (StructuredCommentsFormPtr) MemNew (sizeof (StructuredCommentsFormData));
      if (frm == NULL) {
        FileClose (fp);
        return;
      }
      frm->table = ReadTabTableFromFile(fp);
      FileClose (fp);
      if (frm->table == NULL || frm->table->data.ptrvalue == NULL || frm->table->next == NULL) {
        Message (MSG_ERROR, "Unable to read table from file");
        frm = MemFree (frm);
        return;
      }

      w = FixedWindow (-50, -33, -10, -10, "Parse Structured Comment From File", StdCloseWindowProc);
      SetObjectExtra (w, frm, CleanupStructuredCommentsForm);
      frm->form = (ForM) w;

      frm->input_entityID = bfp->input_entityID;

      h = HiddenGroup (w, -1, 0, NULL);

      frm->table_type = HiddenGroup (h, 0, 2, ChangeStructuredCommentTableType);
      SetGroupSpacing (frm->table_type, 10, 10);
      SetObjectExtra (frm->table_type, frm, NULL);
      RadioButton (frm->table_type, "Table includes accessions for structured comments, one per line");
      RadioButton (frm->table_type, "Table is one structured comment, field names across the top");
      SetValue (frm->table_type, 1);

      frm->match_info = HiddenGroup (h, -1, 0, NULL);
      SetGroupSpacing (frm->match_info, 10, 10);
      ppt1 = StaticPrompt (frm->match_info, "Choose column that identifies sequence for structured comment", 0, dialogTextHeight,
		                       programFont, 'l');
      frm->match_column = PopupList (frm->match_info, TRUE, NULL);
      vnp = frm->table->data.ptrvalue;
      col = 1;
      while (vnp != NULL) {
        if (vnp->data.ptrvalue == NULL) {
          val = "blank";
        } else {
          val = vnp->data.ptrvalue;
        }
        choice_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (choice_fmt) + 15 + StringLen (val)));
        sprintf (choice_str, choice_fmt, col, val);
        PopupItem (frm->match_column, choice_str);
        choice_str = MemFree (choice_str);
        vnp = vnp->next;
        col++;
      }
      SetValue (frm->match_column, 1);

      ppt2 = StaticPrompt (frm->match_info, "Specify relationship of column to sequence", 0, dialogTextHeight,
		                       programFont, 'l');      
      frm->match_type = MatchTypeDialog (frm->match_info, NULL, NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) ppt1,
                                  (HANDLE) frm->match_column,
                                  (HANDLE) ppt2,
                                  (HANDLE) frm->match_type,
                                  NULL);

      ppt3 = StaticPrompt (h, "Choose database name", 0, dialogTextHeight,
		                  programFont, 'l');

      frm->database_name = StructuredCommentDatabaseNameDialog (h, NULL, NULL);

      c = HiddenGroup (h, 2, 0, NULL);
      b = PushButton (c, "Accept", DoParseStructuredComments);
      SetObjectExtra (b, frm, NULL);
      PushButton (c, "Cancel", StdCancelButtonProc);
      AlignObjects (ALIGN_CENTER, (HANDLE) frm->table_type,
                                  (HANDLE) frm->match_info,
                                  (HANDLE) ppt3,
                                  (HANDLE) frm->database_name,
                                  (HANDLE) c,
                                  NULL);
      RealizeWindow (w);
      Show (w);
      Update();
    }
  } 
}


NLM_EXTERN void CreateStructuredCommentsItem (IteM i)
{
  BaseFormPtr        bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  CreateStructuredCommentsBaseForm (bfp);
}


static Boolean AnySequencesInList (ValNodePtr list)
{
  Boolean rval = FALSE;

  while (list != NULL && !rval) {
    if (list->data.ptrvalue != NULL) {
      rval = TRUE;
    }
    list = list->next;
  }
  return rval;
}


static Boolean HeaderHasAtLeastTwoColumns (ValNodePtr header)
{
  if (header == NULL || header->next == NULL) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static Boolean HeaderLooksLikeSeqId (ValNodePtr header)
{
  if (header == NULL) {
    return FALSE;
  } else if (StringICmp (header->data.ptrvalue, "accession") == 0
    || StringICmp (header->data.ptrvalue, "seqid") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


NLM_EXTERN void SubmitterCreateStructuredComments (IteM i)
{
  BaseFormPtr               bfp;
  Char                      path [PATH_MAX];
  FILE *fp;
  Int4                      col;
  MatchTypeData             match_type;
  ValNodePtr                sequence_lists;
  SeqEntryPtr               sep;
  ValNodePtr                err_list = NULL, vnp, header, line, s_row, vnp_s, tmp;
  LogInfoPtr                lip;
  BioseqPtr                 bsp;
  UserObjectPtr             uop, uop_cpy;
  SeqDescrPtr               sdp;
  ValNodePtr                table;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL)
    return;

  if (!GetInputFileName (path, sizeof (path), "", "TEXT")) {
    return;
  }
  fp = FileOpen (path, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s");
    return;
  }
  table = ReadTabTableFromFile(fp);
  FileClose (fp);
  if (table == NULL || table->data.ptrvalue == NULL || table->next == NULL) {
    Message (MSG_ERROR, "Unable to read table from file");
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  col = 0;
  match_type.choice = eTableMatchNucID;
  match_type.match_location = String_location_equals;
  match_type.data = NULL;

  sequence_lists = GetSequenceListsForMatchTypeInTabTable (sep, table->next, col, &match_type, &err_list);

  if (err_list != NULL && !AnySequencesInList(sequence_lists) 
      && HeaderHasAtLeastTwoColumns(table->data.ptrvalue)
      && !HeaderLooksLikeSeqId(table->data.ptrvalue)) {
    err_list = ValNodeFreeData (err_list);
    tmp = FlipTabTableAxes (table);
    table = FreeTabTable (table);
    table = tmp;
    CreateStructuredCommentsForAllFromTable(sep, table->data.ptrvalue, table->next, &err_list);
  } else {
    if (err_list != NULL) {
      lip = OpenLog ("Table Problems");
      for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
        fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
      }
      lip->data_in_log = TRUE;
      CloseLog (lip);
      lip = FreeLog (lip);
      err_list = ValNodeFreeData (err_list);
      if (ANS_YES != Message (MSG_YN, "Continue with table problems")) {
        sequence_lists = FreeSequenceLists(sequence_lists);
        return;
      }
    }

    WatchCursor ();
    Update();

    header = table->data.ptrvalue;
    line = table->next;
    s_row = sequence_lists;

    while (line != NULL && s_row != NULL) {
      vnp_s = s_row->data.ptrvalue;
      if (vnp_s != NULL) {
        uop = UserObjectFromRow (header, line, col, NULL, &err_list);
        while (vnp_s != NULL) {
          bsp = vnp_s->data.ptrvalue;
          uop_cpy = (UserObjectPtr) AsnIoMemCopy (uop, (AsnReadFunc) UserObjectAsnRead, (AsnWriteFunc) UserObjectAsnWrite);
          sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_user);
          sdp->data.ptrvalue = uop_cpy;
          vnp_s = vnp_s->next;
        }
        uop = UserObjectFree (uop);
      }
      line = line->next;
      s_row = s_row->next;
    }
  }

  if (err_list != NULL) {
    lip = OpenLog ("Data Problems");
    for (vnp = err_list; vnp != NULL; vnp = vnp->next) {
      fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }
    lip->data_in_log = TRUE;
    CloseLog (lip);
    lip = FreeLog (lip);
    err_list = ValNodeFreeData (err_list);
  }

  table = FreeTabTable (table);
  sequence_lists = FreeSequenceLists(sequence_lists);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();

}


/*
 * Two BioSources are mergeable if the following conditions are met:
 * 1)	For the following items, either the value for one or both BioSources
 * is not set or the values are identical:
 *    BioSource.genome
 *    BioSource.origin
 *    BioSource.is-focus
 *    BioSource.org.taxname
 *    BioSource.org.common
 *    BioSource.org.orgname.name
 *    BioSource.org.orgname.attrib
 *    BioSource.org.orgname.lineage
 *    BioSource.org.orgname.gcode
 *    BioSource.org.orgname.mgcode
 *    BioSource.org.orgname.div
 * 2) For the following items, either one or both BioSources have empty lists
 * or the lists for the two BioSources are identical:
 *    BioSource.subtype
 *    BioSource.org.syn
 *    BioSource.org.orgname.mod
 *    BioSource.org.mod
 * 3)	No Dbtag in the BioSource.org.db list from one BioSource can have the same
 * db value as an item in the BioSource.org.db list from the other BioSource 
 * unless the two Dbtags have identical tag values.
*/


static ValNodePtr OrgModListsMatch (OrgModPtr mod1, OrgModPtr mod2)
{
  ValNodePtr rval = NULL;
  OrgModPtr cmp1, cmp2;
  Boolean    found;
  CharPtr    str, qual, mismatch_fmt = "Mismatched values for %s";
#if 0
  CharPtr    missing_fmt = "Only one has %s";
#endif

  /* look for missing and mismatch */
  for (cmp1 = mod1; cmp1 != NULL; cmp1 = cmp1->next) {
    found = FALSE;
    for (cmp2 = mod2; cmp2 != NULL; cmp2 = cmp2->next) {
      if (cmp1->subtype == cmp2->subtype) {
        if (StringCmp (cmp1->subname, cmp2->subname) != 0) {
          qual = GetSourceQualName (GetSrcQualFromSubSrcOrOrgMod (cmp1->subtype, TRUE));
          str = (CharPtr) MemNew (sizeof (Char) * (StringLen (mismatch_fmt) + StringLen (qual)));
          sprintf (str, mismatch_fmt, qual);
          ValNodeAddPointer (&rval, 0, str);
        }
        found = TRUE;
      }
    }
#if 0
    if (!found) {
      qual = GetSourceQualName (GetSrcQualFromSubSrcOrOrgMod (cmp1->subtype, TRUE));
      str = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual)));
      sprintf (str, missing_fmt, qual);
      ValNodeAddPointer (&rval, 0, str);
    }
#endif
  }
#if 0
  /* only look for missing this time */
  for (cmp2 = mod2; cmp2 != NULL; cmp2 = cmp2->next) {
    found = FALSE;
    for (cmp1 = mod1; cmp1 != NULL && !found; cmp1 = cmp1->next) {
      if (cmp1->subtype == cmp2->subtype) {
        found = TRUE;
      }
    }
    if (!found) {
      qual = GetSourceQualName (GetSrcQualFromSubSrcOrOrgMod (cmp2->subtype, TRUE));
      str = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual)));
      sprintf (str, missing_fmt, qual);
      ValNodeAddPointer (&rval, 0, str);
    }
  }
#endif

  return rval;
}


static ValNodePtr OkToMergeOrgNames (OrgNamePtr on1, OrgNamePtr on2)
{
  ValNodePtr rval = NULL;

  if (on1 == NULL || on2 == NULL) {
    return NULL;
  }

  /* attrib */
  if (!StringHasNoText (on1->attrib) && !StringHasNoText (on2->attrib)
      && StringCmp (on1->attrib, on2->attrib) != 0)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Orgname.attrib do not match"));
  }
  /* lineage */
  if (!StringHasNoText (on1->lineage) && !StringHasNoText (on2->lineage)
           && StringCmp (on1->lineage, on2->lineage) != 0) 
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Lineages do not match"));
  }
  /* div */
  if (!StringHasNoText (on1->div) && !StringHasNoText (on2->div)
           && StringCmp (on1->div, on2->div) != 0)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Divisions do not match"));
  }
  /* gcode */
  if (on1->gcode != 0 && on2->gcode != 0 && on1->gcode != on2->gcode)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Genetic codes do not match"));
  }
  /* mgcode */
  if (on1->mgcode != 0 && on2->mgcode != 0 && on1->mgcode != on2->mgcode)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Mitochondrial genetic codes do not match"));
  }
  /* OrgMods */
  if (on1->mod != NULL && on2->mod != NULL)
  {
    ValNodeLink (&rval, OrgModListsMatch (on1->mod, on2->mod));
  }
  /* name */
  if (on1->data != NULL && on2->data != NULL 
           && on1->choice != on2->choice) 
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Orgname data differs"));
  }

  return rval;
}


static Boolean ValNodeStringsMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  Boolean rval = TRUE;

  while (rval && vnp1 != NULL && vnp2 != NULL) {
    if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0) {
      rval = FALSE;
    }
    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 != NULL || vnp2 != NULL) {
    rval = FALSE;
  }
  return rval;
}


static Boolean OkToMergeDbtagLists (ValNodePtr vnp1, ValNodePtr vnp2)
{
  Boolean rval = TRUE;
  DbtagPtr dbt1, dbt2;
  ValNodePtr vnp_x;

  while (rval && vnp1 != NULL) {
    dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
    if (dbt1 != NULL) {
      vnp_x = vnp2;
      while (rval && vnp_x != NULL) {
        dbt2 = (DbtagPtr) vnp_x->data.ptrvalue;
        if (dbt2 != NULL && StringCmp (dbt1->db, dbt2->db) == 0
            && !DbtagMatch (dbt1, dbt2)) {
          rval = FALSE;
        }
        vnp_x = vnp_x->next;
      }
    }
    vnp1 = vnp1->next;
  }
  return rval;
}


static ValNodePtr OkToMergeOrgRefs (OrgRefPtr org1, OrgRefPtr org2)
{
  ValNodePtr rval = NULL;

  if (org1 == NULL || org2 == NULL) {
    return NULL;
  }

  /* taxname */
  if (!StringHasNoText (org1->taxname) && !StringHasNoText (org2->taxname)
      && StringCmp (org1->taxname, org2->taxname) != 0)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Taxnames do not match"));
  }
  /* common */
  if (!StringHasNoText (org1->common) && !StringHasNoText (org2->common)
           && StringCmp (org1->common, org2->common) != 0)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Common names do not match"));
  }
  /* synonyms */
  if (org1->syn != NULL && org2->syn != NULL 
           && !ValNodeStringsMatch (org1->syn, org2->syn)) 
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Synonyms do not match"));
  }
  /* mod */
  if (org1->mod != NULL && org2->mod != NULL
           && !ValNodeStringsMatch (org1->mod, org2->mod))
  {
    ValNodeAddPointer (&rval, 0, StringSave ("OrgRef Mods do not match"));
  }
  /* orgname */
  ValNodeLink (&rval, OkToMergeOrgNames (org1->orgname, org2->orgname));
  /* dbtags */
  if (!OkToMergeDbtagLists (org1->db, org2->db)) 
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Dbxrefs do not match"));
  }
  return rval;
}


static ValNodePtr SubtypeListsMatch (SubSourcePtr ssp1, SubSourcePtr ssp2)
{
  ValNodePtr rval = NULL;
  SubSourcePtr cmp1, cmp2;
  Boolean    found;
  CharPtr    str, qual, mismatch_fmt = "Mismatched values for %s";
#if 0
  CharPtr    missing_fmt = "Only one has %s";
#endif

  /* look for missing and mismatch */
  for (cmp1 = ssp1; cmp1 != NULL; cmp1 = cmp1->next) {
    found = FALSE;
    for (cmp2 = ssp2; cmp2 != NULL; cmp2 = cmp2->next) {
      if (cmp1->subtype == cmp2->subtype) {
        if (StringCmp (cmp1->name, cmp2->name) != 0) {
          qual = GetSourceQualName (GetSrcQualFromSubSrcOrOrgMod (cmp1->subtype, FALSE));
          str = (CharPtr) MemNew (sizeof (Char) * (StringLen (mismatch_fmt) + StringLen (qual)));
          sprintf (str, mismatch_fmt, qual);
          ValNodeAddPointer (&rval, 0, str);
        }
        found = TRUE;
      }
    }
#if 0
    if (!found) {
      qual = GetSourceQualName (GetSrcQualFromSubSrcOrOrgMod (cmp1->subtype, FALSE));
      str = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual)));
      sprintf (str, missing_fmt, qual);
      ValNodeAddPointer (&rval, 0, str);
    }
#endif
  }
#if 0
  /* only look for missing this time */
  for (cmp2 = ssp2; cmp2 != NULL; cmp2 = cmp2->next) {
    found = FALSE;
    for (cmp1 = ssp1; cmp1 != NULL && !found; cmp1 = cmp1->next) {
      if (cmp1->subtype == cmp2->subtype) {
        found = TRUE;
      }
    }
    if (!found) {
      qual = GetSourceQualName (GetSrcQualFromSubSrcOrOrgMod (cmp2->subtype, FALSE));
      str = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_fmt) + StringLen (qual)));
      sprintf (str, missing_fmt, qual);
      ValNodeAddPointer (&rval, 0, str);
    }
  }
#endif

  return rval;
}


static ValNodePtr OkToMergeBioSources (BioSourcePtr biop1, BioSourcePtr biop2)
{
  ValNodePtr rval = NULL;

  if (biop1 == NULL || biop2 == NULL) 
  {
    return NULL;
  }

  /* genome */
  if  (biop1->genome != GENOME_unknown && biop2->genome != GENOME_unknown
       && biop1->genome != biop2->genome)
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Locations do not match"));
  } 
  /* origin */
  if (biop1->origin != ORG_UNKNOWN && biop2->origin != ORG_UNKNOWN
             && biop1->origin != biop2->origin) 
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Origins do not match"));
  } 
  /* is-focus */
  if ((biop1->is_focus && !biop2->is_focus) || (!biop1->is_focus && biop2->is_focus))
  {
    ValNodeAddPointer (&rval, 0, StringSave ("Focus does not match"));
  }
  /* subtype */
  if (biop1->subtype != NULL && biop2->subtype != NULL) {
    ValNodeLink (&rval, SubtypeListsMatch(biop1->subtype, biop2->subtype));
  }
  /* orgref */
  ValNodeLink (&rval, OkToMergeOrgRefs(biop1->org, biop2->org));
  return rval;
}


static OrgModPtr OrgModListCopy (OrgModPtr orig)
{
  OrgModPtr list = NULL, prev = NULL, mod;

  while (orig != NULL) {
    mod = AsnIoMemCopy (orig, (AsnReadFunc)OrgModAsnRead, (AsnWriteFunc)OrgModAsnWrite);
    if (prev == NULL) {
      list = mod;
    } else {
      prev->next = mod;
    }
    prev = mod;
    orig = orig->next;
  }
  return list;
}


static void AddOrgNameToOrgName (OrgNamePtr on1, OrgNamePtr on2)
{
  OrgModPtr cpy, mod;

  if (on1 == NULL || on2 == NULL) {
    return;
  }
  /* OrgMods */
  if (on2->mod != NULL) {
    cpy = OrgModListCopy (on2->mod);
    if (on1->mod == NULL) {
      on1->mod = cpy;
    } else {
      mod = on1->mod;
      while (mod->next != NULL) {
        mod = mod->next;
      }
      mod->next = cpy;
    }
  }
  /* attrib */
  if (StringHasNoText (on1->attrib) && !StringHasNoText (on2->attrib)) {
    on1->attrib = MemFree (on1->attrib);
    on1->attrib = StringSave (on2->attrib);
  }
  /* lineage */
  if (StringHasNoText (on1->lineage) && !StringHasNoText (on2->lineage)) {
    on1->lineage = MemFree (on1->lineage);
    on1->lineage = StringSave (on2->lineage);
  }
  /* div */
  if (StringHasNoText (on1->div) && !StringHasNoText (on2->div)) {
    on1->div = MemFree (on1->div);
    on1->div = StringSave (on2->div);
  }
  /* gcode */
  if (on1->gcode == 0) {
    on1->gcode = on2->gcode;
  }
  /* mgcode */
  if (on1->mgcode == 0) {
    on1->mgcode = on2->mgcode;
  }
  
}


static void AddOrgRefToOrgRef (OrgRefPtr org1, OrgRefPtr org2)
{
  ValNodePtr vnp1, vnp2;
  DbtagPtr dbt1, dbt2;

  if (org1 == NULL || org2 == NULL) {
    return;
  }

  /* taxname */
  if (StringHasNoText (org1->taxname) && !StringHasNoText (org2->taxname)) {
    org1->taxname = MemFree (org1->taxname);
    org1->taxname = StringSave (org2->taxname);
  }
  /* common */
  if (StringHasNoText (org1->common) && !StringHasNoText (org2->common)) {
    org1->common = MemFree (org1->common);
    org1->common = StringSave (org1->common);
  }
  /* synonyms */
  if (org1->syn == NULL && org2->syn != NULL) {
    org1->syn = ValNodeDupStringList (org2->syn);
  }
  /* mod */
  if (org1->mod == NULL && org2->mod != NULL) {
    org1->mod = ValNodeDupStringList (org2->mod);
  } 
  /* orgname */
  if (org2->orgname != NULL) {
    if (org1->orgname == NULL) {
      org1->orgname = AsnIoMemCopy (org2->orgname, (AsnReadFunc) OrgNameAsnRead,
                                                   (AsnWriteFunc) OrgNameAsnWrite);
    } else {
      AddOrgNameToOrgName (org1->orgname, org2->orgname);
    }
  }
  /* dbtags */
  if (org2->db != NULL) {
    vnp2 = org2->db;
    while (vnp2 != NULL) {
      dbt2 = vnp2->data.ptrvalue;
      if (dbt2 != NULL) {
        vnp1 = org1->db;
        while (vnp1 != NULL) {
          dbt1 = vnp1->data.ptrvalue;
          if (dbt1 != NULL && StringCmp (dbt1->db, dbt2->db) == 0) {
            break;
          }
          vnp1 = vnp1->next;
        }
        if (vnp1 == NULL) {
          ValNodeAddPointer (&(org1->db), 0, DbtagDup (dbt2));
        }
      }
      vnp2 = vnp2->next;
    }
  }   
}


static SubSourcePtr SubSourceListCopy (SubSourcePtr orig)
{
  SubSourcePtr list = NULL, prev = NULL, mod;

  while (orig != NULL) {
    mod = AsnIoMemCopy (orig, (AsnReadFunc)SubSourceAsnRead, (AsnWriteFunc)SubSourceAsnWrite);
    if (prev == NULL) {
      list = mod;
    } else {
      prev->next = mod;
    }
    prev = mod;
    orig = orig->next;
  }
  return list;
}


static void AddBioSourceToBioSource (BioSourcePtr biop1, BioSourcePtr biop2)
{
  SubSourcePtr cpy, ssp;

  if (biop1 == NULL || biop2 == NULL) {
    return;
  }

  /* genome */
  if (biop1->genome == GENOME_unknown) {
    biop1->genome = biop2->genome;
  }
  /* origin */
  if (biop1->origin == ORG_UNKNOWN) {
    biop1->origin = biop2->origin;
  }
  /* subtype */
  if (biop2->subtype != NULL) {
    cpy = SubSourceListCopy (biop2->subtype);
    if (biop1->subtype == NULL) {
      biop1->subtype = cpy;
    } else {
      ssp = biop1->subtype;
      while (ssp->next != NULL) {
        ssp = ssp->next;
      }
      ssp->next = cpy;
    }
  }

  if (biop1->org == NULL) {
    biop1->org = AsnIoMemCopy (biop2->org, (AsnReadFunc) OrgRefAsnRead,
                                           (AsnWriteFunc) OrgRefAsnWrite);
  } else {
    AddOrgRefToOrgRef (biop1->org, biop2->org);  
  }
  
}


static void MergeBioSourcesInDescrs (SeqDescrPtr sdp, CharPtr id, LogInfoPtr lip)
{
  SeqDescrPtr sdp1, sdp2;
  ObjValNodePtr ovp;
  ValNodePtr errs, vnp;

  sdp1 = sdp;
  while (sdp1 != NULL && sdp1->choice != Seq_descr_source) {
    sdp1 = sdp1->next;
  }
  if (sdp1 != NULL) {
    sdp2 = sdp1->next;
    while (sdp2 != NULL) {
      if (sdp2->choice == Seq_descr_source && sdp2->extended == 1) {
        errs = OkToMergeBioSources (sdp1->data.ptrvalue, sdp2->data.ptrvalue);
        if (errs == NULL) {
          AddBioSourceToBioSource (sdp1->data.ptrvalue, sdp2->data.ptrvalue);
          ovp = (ObjValNodePtr) sdp2;
          ovp->idx.deleteme = TRUE;
          fprintf (lip->fp, "Successfully merged biosources on %s\n", id);
        } else {
          fprintf (lip->fp, "Unable to merge biosources on %s because:\n", id);
          for (vnp = errs; vnp != NULL; vnp = vnp->next) {
            fprintf (lip->fp, "\t%s\n", (CharPtr) vnp->data.ptrvalue);
          }
          errs = ValNodeFreeData (errs);
        }
        lip->data_in_log = TRUE;
      }
      sdp2 = sdp2->next;
    }
  }          
}


static void MergeBiosourcesCallback (BioseqPtr bsp, Pointer data)
{
  CharPtr label;

  if (bsp == NULL || bsp->descr == NULL) {
    return;
  }

  label = GetBioseqLabel (bsp);
  MergeBioSourcesInDescrs (bsp->descr, label, (LogInfoPtr) data);
  label = MemFree (label);
}


static void MergeBioSourcesOnSetCallback (BioseqSetPtr bssp, Pointer data)
{
  CharPtr label;

  if (bssp == NULL || bssp->descr == NULL) {
    return;
  }

  label = GetBioseqSetLabel (bssp);
  MergeBioSourcesInDescrs (bssp->descr, label, (LogInfoPtr) data);
  label = MemFree (label);
}


NLM_EXTERN void MergeBiosources (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr sep;
  LogInfoPtr lip;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  lip = OpenLog ("BioSource Merge Report");
  VisitBioseqsInSep (sep, lip, MergeBiosourcesCallback);
  VisitSetsInSep (sep, lip, MergeBioSourcesOnSetCallback);
  CloseLog (lip);
  if (!lip->data_in_log) {
    Message (MSG_ERROR, "No Biosource found for merging");
  }
  lip = FreeLog (lip);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


typedef struct exportqualifiersform {
  FORM_MESSAGE_BLOCK
  ButtoN source_qualifiers;
  DialoG src_qual_selection;
  GrouP  qual_choice;
  ButtoN launch_excel;
} ExportQualifiersFormData, PNTR ExportQualifiersFormPtr;

static void DoExportQualifiers (ButtoN b)
{
  ExportQualifiersFormPtr frm;
  Char               path [PATH_MAX];
  Uint1              field_type = 0;
  Int2               val;
  FILE *fp;
  ValNodePtr         field_list = NULL;

  frm = (ExportQualifiersFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  if (GetStatus (frm->source_qualifiers)) {
    field_list = DialogToPointer (frm->src_qual_selection);
  }
  val = GetValue (frm->qual_choice);
  if ((val < 2 || val > 6) && field_list == NULL) {
    Message (MSG_ERROR, "Must select qualifiers to export!");
    return;
  }
  switch (val) {
    case 1:
      field_type = 0;
      break;
    case 2:
      field_type = FieldType_feature_field;
      break;
    case 3:
      field_type = FieldType_cds_gene_prot;
      break;
    case 4:
      field_type = FieldType_rna_field;
      break;
    case 5:
      field_type = FieldType_misc;
      break;
    case 6:
      field_type = FieldType_pub;
      break;
  }

  if (GetOutputFileName (path, sizeof (path), NULL)) {
    fp = FileOpen (path, "w");
    if (fp == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
      return;
    }
    ExportFieldTable (field_type, field_list, GetTopSeqEntryForEntityID (frm->input_entityID), fp);
    FileClose (fp);
#ifdef WIN_MSWIN
    if (GetStatus (frm->launch_excel)) {
      Nlm_MSWin_OpenApplication ("excel.exe", path);
    }
#endif

  }
  field_list = FieldTypeListFree (field_list);
  Remove (frm->form);
}


static void ChooseSourceQualifiers (ButtoN b)
{
  ExportQualifiersFormPtr frm;

  frm = (ExportQualifiersFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }
  if (GetStatus (frm->source_qualifiers)) {
    Enable (frm->src_qual_selection);
  } else {
    Disable (frm->src_qual_selection);
  }
}


static void FieldTypeDataFree (ValNodePtr vnp)
{
  if (vnp != NULL) {
    vnp->data.ptrvalue = SourceQualChoiceFree (vnp->data.ptrvalue);
  }
}


static Boolean FieldTypeMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  return AsnIoMemComp (vnp1, vnp2, (AsnWriteFunc) FieldTypeAsnWrite);
}


NLM_EXTERN void ExportQualifiers (IteM i)
{
  BaseFormPtr        bfp;
  ExportQualifiersFormPtr frm;
  WindoW         w;
  GrouP          h, c;
  ButtoN         b;
  SeqEntryPtr    sep;
  ValNodePtr     fields;


#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;


  frm = (ExportQualifiersFormPtr) MemNew (sizeof (ExportQualifiersFormData));
  if (frm == NULL) return;
  frm->input_entityID = bfp->input_entityID;

  w = FixedWindow (-50, -33, -10, -10, "Export Qualifiers", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  frm->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  frm->source_qualifiers = CheckBox (h, "Source Qualifiers", ChooseSourceQualifiers);
  SetObjectExtra (frm->source_qualifiers, frm, NULL);
  SetStatus (frm->source_qualifiers, TRUE);

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  fields = GetFieldListForFieldType (FieldType_source_qual, sep);

  frm->src_qual_selection = ValNodeSelectionDialog (h, fields, TALL_SELECTION_LIST,
                                SummarizeFieldType,
                                FieldTypeDataFree, 
                                FieldTypeCopy,
                                FieldTypeMatch,
                                "qualifiers", 
                                NULL, NULL, TRUE);
  SendMessageToDialog (frm->src_qual_selection, NUM_VIB_MSG + 1);

  frm->qual_choice = NormalGroup (h, 0, 6, "Other Qualifiers", programFont, NULL);
  RadioButton (frm->qual_choice, "None");
  RadioButton (frm->qual_choice, "Feature Qualifiers");
  RadioButton (frm->qual_choice, "CDS-Gene-Prot Qualifiers");
  RadioButton (frm->qual_choice, "RNA Qualifiers");
  RadioButton (frm->qual_choice, "Definition Lines");
  RadioButton (frm->qual_choice, "Publication Titles");
  SetValue (frm->qual_choice, 1);

#ifdef WIN_MSWIN
  frm->launch_excel = CheckBox (h, "Launch Excel", NULL);
#endif

  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", DoExportQualifiers);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) frm->source_qualifiers, 
                              (HANDLE) frm->src_qual_selection, 
                              (HANDLE) frm->qual_choice, 
#ifdef WIN_MSWIN
                              (HANDLE) frm->launch_excel,
#endif
                              (HANDLE) c, 
                              NULL);
  Show (w);
}



static void ExportBankitComment (SeqDescrPtr sdp, Pointer data)
{
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  CharPtr            comment, pound;
  BioseqPtr          bsp;
  Char               id_text[100];
  FILE               *fp;
  CharPtr            buf = NULL;
  Int4               buf_len = 0, line_len;
  
  if (sdp == NULL || sdp->choice != Seq_descr_user || data == NULL) {
    return;
  }

  fp = (FILE *) data;
  
  /* Bankit Comments */
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop != NULL && StringCmp (uop->_class, "SMART_V1.0") != 0) {
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "Submission") == 0) {
      for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
        oip = ufp->label;
        if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {

          /* print out ID */
          bsp = GetSequenceForObject (OBJ_SEQDESC, sdp);
          if (bsp != NULL) {
            SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_text, PRINTID_FASTA_SHORT, sizeof (id_text) - 1);
            fprintf (fp, ">%s\n", id_text);
          }

          comment = (CharPtr) ufp->data.ptrvalue;
          pound = StringSearch (comment, "##");
          while (pound != NULL) {
            line_len = pound - comment + 1;
            if (buf_len < line_len) {
              buf = MemFree (buf);
              buf_len = line_len;
              buf = (CharPtr) MemNew (sizeof (Char) * buf_len);
            }
            StringNCpy (buf, comment, line_len - 1);
            buf[line_len - 1] = 0;
            fprintf (fp, "%s\n", buf);
            comment = pound + 2;
            pound = StringSearch (comment, "##");
          }
          fprintf (fp, "%s\n", comment);
          buf = MemFree (buf);
        }
      }
    }
  }
}



extern void ExportBankitComments (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr    sep;
  Char          path [PATH_MAX];
  FILE *fp;


#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  if (!GetOutputFileName (path, sizeof (path), NULL)
      || (fp = FileOpen (path, "w")) == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  VisitDescriptorsInSep (sep, fp, ExportBankitComment);
  FileClose (fp);

}


extern void TrimPrimerSeqJunk (IteM i)
{
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;


#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  WatchCursor ();
  Update();
  TrimPrimerSeqJunkInSeqEntry (sep, NULL);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


extern void FixiPCRPrimerSeqsMenuItem (IteM i)
{
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;


#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  WatchCursor ();
  Update();
  VisitBioSourcesInSep (sep, NULL, FixiPCRPrimerSeqsCallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void GetContigsInScaffoldsCallback (BioseqPtr bsp, Pointer data)
{
  ValNodePtr PNTR contig_list;
  DeltaSeqPtr dsp;
  SeqLocPtr   loc;
  SeqIdPtr    sip;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext == NULL
      || (contig_list = (ValNodePtr PNTR) data) == NULL) {
    return;
  }

  for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {
      loc = (SeqLocPtr)dsp->data.ptrvalue;
      sip = SeqLocId (loc);
      if (sip != NULL) {
        ValNodeAddPointer (contig_list, 0, SeqIdDup (sip));
      }
    }
  }
}


static void AdjustFeaturesForSequenceCut (BioseqPtr bsp, Int4 cut_start, Int4 cut_len)
{
  SeqFeatPtr  sfp;
  SeqMgrFeatContext fcontext;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
      sfp != NULL;
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext)) {
    if (sfp->idx.deleteme) {
      continue;
    }
    if (fcontext.right >= cut_start) {
      if (fcontext.left >= cut_start && fcontext.right < cut_start + cut_len) {
        /* feature is entirely in removal area */
        sfp->idx.deleteme = TRUE;
      } else {
        /* adjust positions of downstream features */
        SeqEdSeqFeatDelete (sfp, bsp, cut_start, cut_start + cut_len - 1, TRUE);
      }
    }
  }
}

/* need to:
 * remove loc from sequence
 * remove features that corresponded to that loc
 * adjust positions of downstream features
 */
static void RemoveContigFromScaffold (BioseqPtr bsp, SeqIdPtr remove_sip)
{
  DeltaSeqPtr dsp, dsp_prev = NULL, dsp_next;
  SeqLocPtr   loc;
  SeqLitPtr   slip;
  SeqIdPtr    sip;
  Int4        offset = 0, len, orig_offset = 0;

  if (bsp == NULL || remove_sip == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext == NULL) {
    return;
  }

  for (dsp = (DeltaSeqPtr)(bsp->seq_ext); dsp != NULL; dsp = dsp_next) {
    dsp_next = dsp->next;
    if (dsp->choice == 1) {
      loc = (SeqLocPtr)dsp->data.ptrvalue;
      sip = SeqLocId (loc);
      len = SeqLocLen (loc);
      if (SeqIdComp (sip, remove_sip) == SIC_YES) {
        AdjustFeaturesForSequenceCut (bsp, orig_offset, len);

        /* remove from sequence */
        bsp->length -= len;
        if (dsp_prev == NULL) {
          bsp->seq_ext = dsp->next;
        } else {
          dsp_prev->next = dsp->next;
        }
        dsp->next = NULL;
        dsp = DeltaSeqFree (dsp);
        orig_offset += len;
      } else {
        offset += len;
        orig_offset += len;
        dsp_prev = dsp;
      }
    } else if (dsp->choice == 2) {   /* SeqLitPtr */
      slip = (SeqLitPtr)(dsp->data.ptrvalue);
      if (slip != NULL) {
        offset += slip->length;
        orig_offset += slip->length;
      }
      dsp_prev = dsp;
    }
  }

  /* if last deltaseq is gap, remove */
  if (dsp_prev != NULL && dsp_prev->choice == 2 && (slip = dsp_prev->data.ptrvalue) != NULL
    && slip->seq_data == NULL && slip->fuzz == NULL) {
    if (slip->length == 0) {
      bsp->length -= 100;
      AdjustFeaturesForSequenceCut (bsp, orig_offset, 100);
    } else {
      bsp->length -= slip->length;
      AdjustFeaturesForSequenceCut (bsp, orig_offset, slip->length);
    }
    dsp = bsp->seq_ext;
    dsp_prev = NULL;
    while (dsp->next != NULL) {
      dsp_prev = dsp;
      dsp = dsp->next;
    }
    dsp_prev->next = NULL;
    dsp = DeltaSeqFree (dsp);
  }

  /* if first delta seq is gap, remove */
  if ((dsp = (DeltaSeqPtr) bsp->seq_ext) != NULL && dsp->choice == 2
      && (slip = dsp->data.ptrvalue) != NULL
      && slip->seq_data == NULL && slip->fuzz == NULL) {
    if (slip->length == 0) {
      bsp->length -= 100;
      AdjustFeaturesForSequenceCut (bsp, 0, 100);
    } else {
      bsp->length -= slip->length;
      AdjustFeaturesForSequenceCut (bsp, 0, slip->length);
    }
    bsp->seq_ext = dsp->next;
    dsp->next = NULL;
    dsp = DeltaSeqFree (dsp);
  }

  CombineAdjacentGapsOnBioseq (bsp, NULL);
}


static void RemoveContigFromScaffoldCallback (BioseqPtr bsp, Pointer data)
{
  RemoveContigFromScaffold (bsp, (SeqIdPtr)data);
}


typedef struct removecontigfromscaffoldform {
  FORM_MESSAGE_BLOCK

  DialoG list_dlg;

} RemoveContigFromScaffoldFormData, PNTR RemoveContigFromScaffoldFormPtr;


static void DoRemoveContigsFromScaffolds (ButtoN b)
{
  ValNodePtr contig_list, vnp;
  RemoveContigFromScaffoldFormPtr frm;
  SeqEntryPtr sep;

  frm = (RemoveContigFromScaffoldFormPtr)GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  contig_list = DialogToPointer (frm->list_dlg);
  if (contig_list == NULL) {
    Message (MSG_ERROR, "No sequences selected");
    return;
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  for (vnp = contig_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (sep, vnp->data.ptrvalue, RemoveContigFromScaffoldCallback);
  }

  DeleteMarkedObjects (frm->input_entityID, 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);   	
  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  Remove (frm->form);
  Update();
}

extern void RemoveContigFromScaffoldMenuItem (IteM i)
{
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;
  ValNodePtr     contig_list = NULL;
  RemoveContigFromScaffoldFormPtr frm;
  WindoW         w;
  GrouP          h, c;
  ButtoN         b;
  PrompT         p1;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  /* first, need to get list of contigs in scaffolds */
  VisitBioseqsInSep (sep, &contig_list, GetContigsInScaffoldsCallback);

  if (contig_list == NULL) {
    Message (MSG_ERROR, "No scaffolds found");
    return;
  }

  /* get list of contigs to remove */
  frm = (RemoveContigFromScaffoldFormPtr) MemNew (sizeof (RemoveContigFromScaffoldFormData));
  w = FixedWindow (-50, -33, -10, -10, "Remove Contigs from Scaffolds", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupFormProc);
  frm->form = (ForM) w;
  frm->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  p1 = StaticPrompt (h, "Choose contigs to remove", 0, dialogTextHeight, programFont, 'c');
  frm->list_dlg = ValNodeSelectionDialog (h, contig_list, TALL_SELECTION_LIST,
                                          ValNodeSeqIdName, ValNodeSeqIdFree,
                                          ValNodeSeqIdCopy, ValNodeSeqIdMatch,
                                          "seqid", NULL, NULL, TRUE);
  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoRemoveContigsFromScaffolds);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) p1, (HANDLE) frm->list_dlg, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


static void ConvertRptUnitRangeToLocationCallback (SeqFeatPtr sfp, Pointer data)
{
  SeqLocPtr new_slp;
  GBQualPtr gbqual, gbq_prev = NULL;
  LogInfoPtr lip;
  BioseqPtr bsp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) {
    return;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);

  lip = (LogInfoPtr) data;

  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
    if (StringICmp (gbqual->qual, "rpt_unit_range") == 0 && bsp != NULL) {
      new_slp = ParseSimpleSeqLoc (gbqual->val, bsp);
      if (new_slp == NULL) {
        if (lip != NULL) {
          fprintf (lip->fp, "Unable to parse location from '%s'\n", gbqual->val);
          lip->data_in_log = TRUE;
        }
        gbq_prev = gbqual;
      } else {
        sfp->location = SeqLocFree (sfp->location);
        sfp->location = new_slp;
        if (gbq_prev == NULL) {
          sfp->qual = gbqual->next;
        } else {
          gbq_prev->next = gbqual->next;
        }
        gbqual->next = NULL;
        gbqual = GBQualFree (gbqual);
        break;
      }
    } else {
      gbq_prev = gbqual;
    }
  }
}

extern void ConvertRptUnitRangeToLocation (IteM i)
{
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;
  LogInfoPtr     lip;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  lip = OpenLog ("Conversions");
  
  VisitFeaturesInSep (sep, lip, ConvertRptUnitRangeToLocationCallback);

  CloseLog (lip);
  FreeLog (lip);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}


typedef struct retranscribemrnaform {
  FORM_MESSAGE_BLOCK

  TexT locus_tag;

  ButtoN leave_dlg_up;
} RetranscribemRNAFormData, PNTR RetranscribemRNAFormPtr;

static void CleanupRetranscribemRNAForm (GraphiC g, VoidPtr data)

{
  RetranscribemRNAFormPtr frm;

  frm = (RetranscribemRNAFormPtr) data;
  if (frm != NULL) {
    ObjMgrFreeUserData (frm->input_entityID, frm->procid, frm->proctype, frm->userkey);
  }
  StdCleanupFormProc (g, data);
}


static void RetranscribemRNAFormMessage (ForM f, Int2 mssg)

{
  RetranscribemRNAFormPtr frm;

  frm = (RetranscribemRNAFormPtr) GetObjectExtra (f);
  if (frm != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      default :
        if (frm->appmessage != NULL) {
          frm->appmessage (f, mssg);
        }
        break;
    }
  }
}

static Int2 LIBCALLBACK RetranscribemRNAFormMsgFunc (OMMsgStructPtr ommsp)
{
  WindoW                   currentport,
                           temport;
  OMUserDataPtr            omudp;
  RetranscribemRNAFormPtr  frm = NULL;
  
  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  frm = (RetranscribemRNAFormPtr) omudp->userdata.ptrvalue;
  if (frm == NULL) return OM_MSG_RET_ERROR;

  currentport = ParentWindow (frm->form);
  temport = SavePort (currentport);
  UseWindow (currentport);
  Select (frm->form);
  switch (ommsp->message) 
  {
      case OM_MSG_UPDATE:
          break;
      case OM_MSG_DESELECT:
          break;

      case OM_MSG_SELECT: 
          break;
      case OM_MSG_DEL:
          Remove (frm->form);
          break;
      case OM_MSG_HIDE:
          break;
      case OM_MSG_SHOW:
          break;
      case OM_MSG_FLUSH:
          Remove (frm->form);	
          break;
      default:
          break;
  }
  RestorePort (temport);
  UseWindow (temport);
  return OM_MSG_RET_OK;
}


static SeqFeatPtr FindGeneInSeqEntryByLocusTag (SeqEntryPtr sep, CharPtr locus_tag)
{
  SeqMgrFeatContext fcontext;
  BioseqSetPtr bssp;
  SeqFeatPtr   gene = NULL;


  if (sep == NULL || sep->data.ptrvalue == NULL || StringHasNoText (locus_tag)) {
    return NULL;
  }

  if (IS_Bioseq (sep)) {
    gene = SeqMgrGetGeneByLocusTag (sep->data.ptrvalue, locus_tag, &fcontext);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL && gene == NULL; sep = sep->next) {
      gene = FindGeneInSeqEntryByLocusTag (sep, locus_tag);
    }
  }
  return gene;
}


static void FindLocusTag (ButtoN b)
{
  RetranscribemRNAFormPtr frm;
  CharPtr locus_tag;
  SeqEntryPtr sep;
  SeqFeatPtr gene, mrna;
  ValNode    vn;
  SeqMgrFeatContext fcontext;

  frm = (RetranscribemRNAFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  if (TextHasNoText (frm->locus_tag)) {
    return;
  }

  locus_tag = SaveStringFromText (frm->locus_tag);

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);
  gene = FindGeneInSeqEntryByLocusTag(sep, locus_tag);
  if (gene == NULL) {
    Message (MSG_ERROR, "Unable to find gene for locus tag %s", locus_tag);
  } else {
    mrna = SeqMgrGetOverlappingmRNA (gene->location, &fcontext);
    if (mrna == NULL) {
      vn.data.ptrvalue = gene;
    } else {
      vn.data.ptrvalue = mrna;
    }
    vn.choice = OBJ_SEQFEAT;
    vn.next = NULL;
    ScrollToDiscrepancyItem (&vn, NULL);
  }
  locus_tag = MemFree (locus_tag);
}


static void OpenViewerForProduct (BioseqPtr bsp)
{
  Uint2 entityID;
  SeqEntryPtr sep;
  Char  str[255];
  Int2          handled;

  if (bsp == NULL) {
    return;
  }

  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), str, PRINTID_REPORT, sizeof (str));

  if (bsp->idx.entityID == 0) {
    entityID = ObjMgrRegister (OBJ_BIOSEQ, bsp);
  } else {
    entityID = bsp->idx.entityID;
  }

  SeqEntrySetScope (NULL);
  WatchCursor ();
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) {
    sep = SeqEntryNew ();
    if (sep != NULL) {
      sep->choice = 1;
      sep->data.ptrvalue = bsp;
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
    }
    sep = GetTopSeqEntryForEntityID (entityID);
  }

  ArrowCursor ();
  Update ();

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
}


static void FeaturePropOneFeatureThruAlign (SeqFeatPtr sfp, SeqAlignPtr salp, Int4 new_len)
{
  Boolean   partial5, partial3;
  CdRegionPtr crp;
  CodeBreakPtr cbp;
  RnaRefPtr    rrp;
  tRNAPtr      trp;

  if (sfp == NULL || salp == NULL) {
    return;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  ReplaceComplexLocation (sfp->location, salp, new_len, 1, 2);

  SetSeqLocPartial (sfp->location, partial5, partial3);
  sfp->partial |= partial5 | partial3;

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          ReplaceComplexLocation (cbp->loc, salp, new_len, 1, 2);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          ReplaceComplexLocation (trp->anticodon, salp, new_len, 1, 2);
        }
      }
      break;
    default :
      break;
  }
}


static BioseqPtr RetranscribeProduct (SeqFeatPtr mrna, BioseqPtr product)
{
  Int4          len;
  ByteStorePtr  bs;
  CharPtr       rnaseq;
  SeqFeatPtr    sfp;
  SeqMgrFeatContext fcontext;
  BioseqPtr         new_product;
  SeqAlignPtr       salp;
  Boolean           revcomp;
  SeqLocPtr         slp;

  if (mrna == NULL) {
    return product;
  }

  rnaseq = GetSequenceByFeature (mrna);
  if (rnaseq == NULL) return product;
  len = (Int4) StringLen (rnaseq);

  bs = BSNew (len + 2);
  if (bs == NULL) {
    rnaseq = MemFree (rnaseq);
    return product;
  }
  BSWrite (bs, (VoidPtr) rnaseq, len);
  MemFree (rnaseq);

  new_product = BioseqNew ();

  /* have to create temporary ID or alignment will fail */
  new_product->id = MakeSeqID ("temporary_retranscription_id");
  new_product->repr = Seq_repr_raw;
  new_product->mol = Seq_mol_rna;
  new_product->seq_data_type = Seq_code_iupacna;
  new_product->seq_data = (SeqDataPtr) bs;
  new_product->length = BSLen (bs);
  BioseqPack (new_product);


  if (product == NULL) {
    product = new_product;
    /* free the ID, so that the calling function will know to use the product ID */
    product->id = SeqIdFree (product->id);
  } else {
    /* align new and old sequences */
    salp = Sequin_GlobalAlign2Seq (product, new_product, &revcomp);
    if (salp != NULL) {
      /* adjust feature locations */
      for (sfp = SeqMgrGetNextFeature (product, NULL, 0, 0, &fcontext);
           sfp != NULL;
           sfp = SeqMgrGetNextFeature (product, sfp, 0, 0, &fcontext)) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          /* extend gene to length of entire sequence */
          slp = SeqLocIntNew (0, new_product->length, SeqLocStrand (sfp->location), SeqLocId (sfp->location)); 
          sfp->location = SeqLocFree (sfp->location);
          sfp->location = slp;
        } else {
          FeaturePropOneFeatureThruAlign (sfp, salp, new_product->length);
        }
      }
    }

    product->seq_data = SeqDataFree(product->seq_data, product->seq_data_type);
    product->repr = Seq_repr_raw;
    product->mol = Seq_mol_rna;
    product->seq_data_type = new_product->seq_data_type;
    product->seq_data = new_product->seq_data;
    new_product->seq_data = NULL;
    product->length = new_product->length;
    new_product = BioseqFree (new_product);
  }
  return product;
}


static void DoRetranscription (ButtoN b)
{
  RetranscribemRNAFormPtr frm;
  SelStructPtr sel;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  BioseqPtr bsp;
  SeqIdPtr  sip;
  Char      id_buf[255];
  Boolean   found_any = FALSE, any_problems = FALSE, need_unlock;

  frm = (RetranscribemRNAFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  sel = ObjMgrGetSelected ();

  while (sel != NULL) {
    if (sel->itemtype == OBJ_SEQFEAT) {
      sfp = SeqMgrGetDesiredFeature (sel->entityID, NULL, sel->itemID, 0, NULL, &fcontext);
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_mRNA) {
        if (sfp->product == NULL) {
          Message (MSG_ERROR, "mRNA has no product!");
          any_problems = TRUE;
        } else {
          sip = SeqLocId (sfp->product);
          need_unlock = TRUE;
          bsp = BioseqLockById (sip);
          if (bsp == NULL) {
            need_unlock = FALSE;
            SeqIdWrite (sip, id_buf, PRINTID_REPORT, sizeof (id_buf));
            Message (MSG_ERROR, "Unable to fetch product sequence %s", id_buf);
          }
          /* change fetched sequence to match transcription of mRNA feature */
          bsp = RetranscribeProduct (sfp, bsp);
          if (bsp != NULL && bsp->id == NULL) {
            bsp->id = SeqIdDup (sip);
          }

          /* open a new viewer for the altered sequence (be sure to set dirty flag) */
          OpenViewerForProduct (bsp);
          if (need_unlock) {
            BioseqUnlock (bsp);
          }
          found_any = TRUE;
        }
      }
    }
    sel = sel->next;
  }

  if (found_any && !any_problems && !GetStatus (frm->leave_dlg_up)) {
    Remove (frm->form);
  }
  Update();
}


extern void RetranscribemRNA (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;
  RetranscribemRNAFormPtr frm;
  WindoW            w;
  GrouP             h, g, c;
  ButtoN            b;
  OMUserDataPtr            omudp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  frm = (RetranscribemRNAFormPtr) MemNew (sizeof (RetranscribemRNAFormData));
  if (frm == NULL)
  {
    return;
  }
  
  frm->input_entityID = bfp->input_entityID;
  w = FixedWindow (-50, -33, -10, -10, "Retranscribe Selected mRNA", StdCloseWindowProc);
  SetObjectExtra (w, frm, CleanupRetranscribemRNAForm);
  frm->form = (ForM) w;
  frm->formmessage = RetranscribemRNAFormMessage;

  /* register to receive update messages */
  frm->userkey = OMGetNextUserKey ();
  frm->procid = 0;
  frm->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (frm->input_entityID, frm->procid, frm->proctype, frm->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) frm;
    omudp->messagefunc = RetranscribemRNAFormMsgFunc;
  }


#ifndef WIN_MAC
  CreateStdValidatorFormMenus (w);
#endif
 
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 3, 0, NULL);
  StaticPrompt (g, "Navigate to Locus Tag", 0, popupMenuHeight, programFont, 'l');
  frm->locus_tag = DialogText (g, "", 10, NULL);
  b = PushButton (g, "Find", FindLocusTag);
  SetObjectExtra (b, frm, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Retranscribe", DoRetranscription);
  SetObjectExtra (b, frm, NULL);

  b = PushButton (c, "Close", StdCancelButtonProc);

  frm->leave_dlg_up = CheckBox (c, "Leave Dialog Up After Retranscription", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  Show (w);
}


static void ClearCDDCallback (SeqFeatPtr sfp, Pointer data)
{
  ObjectIdPtr    oip;
  UserObjectPtr  uop;
  ValNodePtr     vnp, vnp_n, vnp_p;
  DbtagPtr       dbtag;
  Int4           len;

  if (sfp == NULL) {
    return;
  }
  if ((uop = sfp->ext) != NULL && (oip = uop->type) != NULL
      && (len = StringLen (oip->str)) >= 12
      && StringNICmp (oip->str, "cdd", 3) == 0
      && StringICmp (oip->str + len - 9, "ScoreData") == 0) {
    sfp->idx.deleteme = TRUE;
  } else if (sfp->dbxref != NULL) {
    vnp_p = NULL;
    for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp_n) {
      vnp_n = vnp->next;
      dbtag = (DbtagPtr) vnp->data.ptrvalue;
      if (dbtag != NULL && StringICmp (dbtag->db, "CDD") == 0) {
        if (vnp_p == NULL) {
          sfp->dbxref = vnp_n;
        } else {
          vnp_p->next = vnp_n;
        }
        vnp->next = NULL;
        vnp->data.ptrvalue = DbtagFree (vnp->data.ptrvalue);
        vnp = ValNodeFree (vnp);
      } else {
        vnp_p = vnp;
      }
    }
  }
}


NLM_EXTERN void CleanupCDD (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  WatchCursor();
  Update();
  VisitFeaturesInSep (sep, NULL, ClearCDDCallback);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, sep);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct trimseqbyalignform {
  FORM_MESSAGE_BLOCK

  ButtoN trim_before_btn;
  TexT trim_before_coord_txt;
  ButtoN trim_after_btn;
  TexT trim_after_coord_txt;
  ButtoN accept_btn;

  SeqAlignPtr salp;
} TrimSeqByAlignFormData, PNTR TrimSeqByAlignFormPtr;


static void CleanupTrimSeqByAlignForm (GraphiC g, VoidPtr data)

{
  TrimSeqByAlignFormPtr frm;

  frm = (TrimSeqByAlignFormPtr) data;
  if (frm != NULL) {
    frm->salp = SeqAlignFree (frm->salp);
  }
  StdCleanupFormProc (g, data);
}



static void ChangeTrimSeqByAln(TrimSeqByAlignFormPtr frm)
{
  Boolean is_ok = FALSE;
  CharPtr str;
  Int4 val1 = 0, val2 = 0;

  if (frm == NULL) {
    return;
  }

  if (GetStatus (frm->trim_before_btn)) {
    Enable (frm->trim_before_coord_txt);
    str = SaveStringFromText (frm->trim_before_coord_txt);
    if (!StringHasNoText (str)) {
      val1 = atoi (str); 
      if (val1 > 0) {
        is_ok = TRUE;
      }
    }
    str = MemFree (str);
  }

  if (GetStatus (frm->trim_after_btn)) {
    Enable (frm->trim_after_coord_txt);
    str = SaveStringFromText (frm->trim_after_coord_txt);
    if (!StringHasNoText (str)) {
      val2 = atoi (str); 
      if (val2 > 0) {
        is_ok = TRUE;
      } else {
        is_ok = FALSE;
      }
    }
    str = MemFree (str);
  }

  if (val1 > 0 && val2 > 0 && val2 <= val1) {
    is_ok = FALSE;
  }

  if (is_ok) {
    Enable (frm->accept_btn);
  } else {
    Disable (frm->accept_btn);
  }

}


static void ChangeTrimSeqByAlnBtn(ButtoN b)
{
  TrimSeqByAlignFormPtr frm;
  frm = (TrimSeqByAlignFormPtr) GetObjectExtra (b);
  ChangeTrimSeqByAln (frm);
}
  

static void ChangeTrimSeqByAlnTxt(TexT t)
{
  TrimSeqByAlignFormPtr frm;
  frm = (TrimSeqByAlignFormPtr) GetObjectExtra (t);
  ChangeTrimSeqByAln (frm);
}


static Boolean AreAllSequencesInAlignmentAtPosition (SeqAlignPtr salp, Int4 pos)
{
  Int4 row, num_rows, val;

  if (salp == NULL || pos < 0) {
    return FALSE;
  }
  num_rows = AlnMgr2GetNumRows(salp);
  for (row = 1; row <= num_rows; row++) {
    val = AlnMgr2MapSeqAlignToBioseq (salp, pos, row);
    if (val == -2) {
      return FALSE;
    }
  }
  return TRUE;
}


static void DoTrimSeqsByAln (ButtoN b)
{
  TrimSeqByAlignFormPtr frm;
  CharPtr str;
  Int4 val_start = -1, val_stop = -1, num_rows, row, pos, val;
  SeqIdPtr sip;
  BioseqPtr bsp;
  SeqLocPtr slp;
  SeqEntryPtr sep;
  SeqAnnotPtr sap;
  LogInfoPtr  lip;

  frm = (TrimSeqByAlignFormPtr) GetObjectExtra (b);

  if (GetStatus (frm->trim_before_btn)) {
    str = SaveStringFromText (frm->trim_before_coord_txt);
    val_start = atoi (str) - 1;
    str = MemFree (str);
    if (!AreAllSequencesInAlignmentAtPosition(frm->salp, val_start)) {
      if (Message(MSG_OKC, "Not all sequences are in the alignment at position %d - are you sure you want to continue?", val_start) == ANS_CANCEL) {
        return;
      }
    }
  }
  if (GetStatus (frm->trim_after_btn)) {
    str = SaveStringFromText (frm->trim_after_coord_txt);
    val_stop = atoi (str) - 1;
    str = MemFree (str);
    if (!AreAllSequencesInAlignmentAtPosition(frm->salp, val_stop)) {
      if (Message(MSG_OKC, "Not all sequences are in the alignment at position %d - are you sure you want to continue?", val_stop) == ANS_CANCEL) {
        return;
      }
    }
  }
  lip = OpenLog ("Trimmed Locations");

  num_rows = AlnMgr2GetNumRows(frm->salp);

  /* for each sequence in alignment*/
  for (row = 1; row <= num_rows; row++) {
    sip = AlnMgr2GetNthSeqIdPtr (frm->salp, row);
    bsp = BioseqLockById (sip);
    if (bsp != NULL) {
      if (AlnMgr2GetNthStrand(frm->salp, row) == Seq_strand_minus) {
        /* trim 3' ends first */
        if (val_start > -1) {
          /* get sequence position for alignment coordinate */
          pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val_start, row);
          if (pos < 0) {
            val = val_start + 1;
            while (pos < 0 && val > -1) {
              pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val, row);
              val ++;
            }
          }

          if (pos > -1 && pos < bsp->length - 1) {
            /* remove sequence after position */
            slp = SeqLocIntNew (pos + 1, bsp->length - 1, Seq_strand_plus, sip);
            SeqDeleteByLoc (slp, TRUE, FALSE);
            LogTrimmedLocation (lip, slp);
            TrimQualityScores (bsp, bsp->length - 1 - pos, FALSE);
          }
        }
        /* now trim 5' end */
        if (val_stop > -1) {
          /* get sequence position for alignment coordinate */
          pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val_stop, row);
          if (pos < 0) {
            val = val_stop - 1;
            while (pos < 0 && val > -1) {
              pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val, row);
              val --;
            }
          }

          if (pos > 0) {
            /* remove sequence before position */
            slp = SeqLocIntNew (0, pos - 1, Seq_strand_plus, sip);
            SeqDeleteByLoc (slp, TRUE, FALSE);
            LogTrimmedLocation (lip, slp);
            TrimQualityScores (bsp, pos, TRUE);
          }
        }
      } else {
        /* trim 3' ends first */
        if (val_stop > -1) {
          /* get sequence position for alignment coordinate */
          pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val_stop, row);
          if (pos < 0) {
            val = val_stop - 1;
            while (pos < 0 && val > -1) {
              pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val, row);
              val --;
            }
          }

          if (pos > -1 && pos < bsp->length - 1) {
            /* remove sequence after position */
            slp = SeqLocIntNew (pos + 1, bsp->length - 1, Seq_strand_plus, sip);
            SeqDeleteByLoc (slp, TRUE, FALSE);
            LogTrimmedLocation (lip, slp);
            TrimQualityScores (bsp, bsp->length - 1 - pos, FALSE);
          }
        }
        /* now trim 5' end */
        if (val_start > -1) {
          /* get sequence position for alignment coordinate */
          pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val_start, row);
          if (pos < 0) {
            val = val_start + 1;
            while (pos < 0 && val > -1) {
              pos = AlnMgr2MapSeqAlignToBioseq (frm->salp, val, row);
              val ++;
            }
          }

          if (pos > -1) {
            /* remove sequence before position */
            slp = SeqLocIntNew (0, pos - 1, Seq_strand_plus, sip);
            SeqDeleteByLoc (slp, TRUE, FALSE);
            LogTrimmedLocation (lip, slp);
            TrimQualityScores (bsp, pos, TRUE);
          }
        }
      }
      BioseqUnlock (bsp);
    }
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);
  sap = (SeqAnnotPtr) FindSeqAlignInSeqEntry (sep, OBJ_SEQANNOT);
  if (sap != NULL) {
    CloseAlignmentEditor (sap->idx.entityID, sap->idx.itemID);

    sap->idx.deleteme = TRUE;
    DeleteMarkedObjects (frm->input_entityID, 0, NULL);
    fprintf (lip->fp, "Removed Alignment\n");
    lip->data_in_log = TRUE;
  }

  CloseLog (lip);
  lip = FreeLog (lip);

  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  Remove (frm->form);
}


NLM_EXTERN void TrimSequencesByAlignment (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;
  TrimSeqByAlignFormPtr frm;
  WindoW                w;
  ButtoN                b;
  GrouP                 h, g, c;
  SeqAlignPtr           salp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  salp  = (SeqAlignPtr)FindSeqAlignInSeqEntry (sep, OBJ_SEQALIGN);
  if (!salp) {
    Message (MSG_ERROR, "No alignments found!");
    return;
  }

  frm = (TrimSeqByAlignFormPtr) MemNew (sizeof (TrimSeqByAlignFormData));
  if (frm == NULL)
  {
    return;
  }

  frm->salp = SeqAlignDup (salp);
  AlnMgr2IndexSingleChildSeqAlign(frm->salp);

  frm->input_entityID = bfp->input_entityID;
  w = FixedWindow (-50, -33, -10, -10, "Trim Sequences by Alignment Coordinates", StdCloseWindowProc);
  SetObjectExtra (w, frm, CleanupTrimSeqByAlignForm);
  frm->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);
  frm->trim_before_btn = CheckBox (g, "Trim before alignment coordinate", ChangeTrimSeqByAlnBtn);
  SetObjectExtra (frm->trim_before_btn, frm, NULL);
  SetStatus (frm->trim_before_btn, FALSE);
  frm->trim_before_coord_txt = DialogText (g, "", 10, ChangeTrimSeqByAlnTxt);
  SetObjectExtra (frm->trim_before_coord_txt, frm, NULL);
  Disable (frm->trim_before_coord_txt);

  frm->trim_after_btn = CheckBox (g, "Trim after alignment coordinate", ChangeTrimSeqByAlnBtn);
  SetObjectExtra (frm->trim_after_btn, frm, NULL);
  SetStatus (frm->trim_after_btn, FALSE);
  frm->trim_after_coord_txt = DialogText (g, "", 10, ChangeTrimSeqByAlnTxt);
  SetObjectExtra (frm->trim_after_coord_txt, frm, NULL);
  Disable (frm->trim_after_coord_txt);

  c = HiddenGroup (h, 4, 0, NULL);
  frm->accept_btn = PushButton (c, "Accept", DoTrimSeqsByAln);
  SetObjectExtra (frm->accept_btn, frm, NULL);

  b = PushButton (c, "Close", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  Show (w);
}

/* add a repeat-region with satellite-type of microsatellite and repeat-type of tandem.  
* Make the feature complete over the entire span of the sequence.
*/
static void MakeCommonFeatureMicrosatelliteCallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr sfp;
  ImpFeatPtr ifp;
  GBQualPtr  gbq;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, NULL);
  ifp = ImpFeatNew ();
  ifp->key = StringSave ("repeat_region");
  sfp->data.value.ptrvalue = ifp;
  gbq = GBQualNew ();
  gbq->qual = StringSave ("rpt_type");
  gbq->val = StringSave ("tandem");
  sfp->qual = gbq;
  gbq = GBQualNew ();
  gbq->qual = StringSave ("satellite");
  gbq->val = StringSave ("microsatellite");
  gbq->next = sfp->qual;
  sfp->qual = gbq;
}


NLM_EXTERN void MakeCommonFeatureMicrosatellite (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, NULL, MakeCommonFeatureMicrosatelliteCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void MakeCommonFeatureControlRegionCallback (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr sfp;
  ImpFeatPtr ifp;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, NULL);
  ifp = ImpFeatNew ();
  ifp->key = StringSave ("misc_feature");
  sfp->data.value.ptrvalue = ifp;
  sfp->comment = StringSave ("control region");
  SetSeqLocPartial (sfp->location, TRUE, TRUE);
  sfp->partial = TRUE;
}


NLM_EXTERN void MakeCommonFeatureControlRegion (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, NULL, MakeCommonFeatureControlRegionCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void RemoveNonNIslandsCallback (BioseqPtr bsp, Pointer data)
{
  CharPtr     bases, cp;
  Int4        pos = 0, island_start = 0, island_len = 0, i;
  SeqDataPtr  seq_data;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) return;

  for (cp = bases; *cp != 0; cp++) {
    if (*cp == 'N') {
      if (island_len < 3 && island_len > 0) {
        for (i = island_start; i < pos; i++) {
          *(bases + i) = 'N';
        }
      }
      island_len = 0;
    } else {
      if (island_len == 0) {
        island_start = pos;
      }
      island_len++;
    }
    pos++;
  }
  if (island_len < 3 && island_len > 0) {
    for (i = island_start; i < pos; i++) {
      *(bases + i) = 'N';
    }
  }

  seq_data = (SeqDataPtr) BSNew (bsp->length);
  AddBasesToByteStore ((ByteStorePtr) seq_data, bases);
  bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data_type = Seq_code_iupacna;
  bsp->seq_data = seq_data;
  bases = MemFree (bases);
  BioseqPack (bsp);
}

NLM_EXTERN void RemoveNonNIslands (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, NULL, RemoveNonNIslandsCallback);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
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


static void RemoveUTRsFromCDSCallback (SeqFeatPtr sfp, Pointer data)
{
  BioseqPtr bsp;
  SeqFeatPtr utr;
  SeqMgrFeatContext fcontext;
  Int4              start, stop, right;
  Uint1             cds_strand;
  Boolean           changed, partial5, partial3, retranslate = FALSE;
  CdRegionPtr       crp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return;
  }

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) {
    return;
  }

  start = SeqLocStart (sfp->location);
  stop = SeqLocStop (sfp->location);
  if (start < stop) {
    right = stop;
  } else {
    right = start;
  }
  cds_strand = SeqLocStrand (sfp->location);

  for (utr = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_5UTR, &fcontext);
       utr != NULL;
       utr = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_5UTR, &fcontext)) {
    if (fcontext.left > right) {
      break;
    } else if (!StrandOk (cds_strand, SeqLocStrand (utr->location))) {
      continue;
    }
    if ((cds_strand == Seq_strand_minus && fcontext.right == start) 
        || (cds_strand != Seq_strand_minus && fcontext.left == start)) {
      sfp->location = SeqLocDelete (sfp->location, SeqLocId (sfp->location), fcontext.left, fcontext.right, FALSE, &changed);
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        crp->frame = 1;
      }
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      SetSeqLocPartial (sfp->location, FALSE, partial3);
      retranslate = TRUE;
      break;
    } 
  }

  for (utr = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_3UTR, &fcontext);
       utr != NULL;
       utr = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_3UTR, &fcontext)) {
    if (fcontext.left > right) {
      break;
    } else if (!StrandOk (cds_strand, SeqLocStrand (utr->location))) {
      continue;
    }
    if ((cds_strand == Seq_strand_minus && fcontext.left == stop) 
        || (cds_strand != Seq_strand_minus && fcontext.right == stop)) {
      sfp->location = SeqLocDelete (sfp->location, SeqLocId (sfp->location), fcontext.left, fcontext.right, FALSE, &changed);
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      SetSeqLocPartial (sfp->location, partial5, FALSE);
      retranslate = TRUE;
      break;
    } 
  }
  
  if (retranslate) {
    RetranslateOneCDS (sfp, sfp->idx.entityID, FALSE, FALSE);
  }
}


NLM_EXTERN void RemoveUTRsFromCDSs (IteM i)
{
  BaseFormPtr       bfp;
  SeqEntryPtr       sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitFeaturesInSep (sep, NULL, RemoveUTRsFromCDSCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef enum {
  eFeatureBetweenPolicy_NoFeatPartial = 1,
  eFeatureBetweenPolicy_NoFeatNoFeat,
  eFeatureBetweenPolicy_ExtendComplete,
  eFeatureBetweenPolicy_ExtendPartial
} EFeatureBetweenPolicy;

typedef struct featurebetweenendpolicy {
  Uint1 policy;
  ValNodePtr feat_type;
  ConstraintChoiceSetPtr constraint;
} FeatureBetweenEndPolicyData, PNTR FeatureBetweenEndPolicyPtr;


static FeatureBetweenEndPolicyPtr FeatureBetweenEndPolicyNew (void)
{
  FeatureBetweenEndPolicyPtr pol;

  pol = (FeatureBetweenEndPolicyPtr) MemNew (sizeof (FeatureBetweenEndPolicyData));
  pol->policy = eFeatureBetweenPolicy_ExtendPartial;
  return pol;
}


static FeatureBetweenEndPolicyPtr FeatureBetweenEndPolicyFree (FeatureBetweenEndPolicyPtr pol)
{
  if (pol != NULL) {
    pol->feat_type = ValNodeFree (pol->feat_type);
    pol = MemFree (pol);
  }
  return pol;
}


typedef struct addfeaturebetween {
  FeatureBetweenEndPolicyPtr left_policy;
  ValNodePtr feat_type_add;
  ApplyFeatureDetailsPtr details;
  FeatureBetweenEndPolicyPtr right_policy;
} AddFeatureBetweenData, PNTR AddFeatureBetweenPtr;

static AddFeatureBetweenPtr AddFeatureBetweenNew (void)
{
  AddFeatureBetweenPtr add;

  add = (AddFeatureBetweenPtr) MemNew (sizeof (AddFeatureBetweenData));
  return add;
}


static AddFeatureBetweenPtr AddFeatureBetweenFree (AddFeatureBetweenPtr add)
{
  if (add != NULL) {
    add->left_policy = FeatureBetweenEndPolicyFree (add->left_policy);
    add->feat_type_add = ValNodeFree (add->feat_type_add);
    add->details = ApplyFeatureDetailsFree (add->details);
    add->right_policy = FeatureBetweenEndPolicyFree (add->right_policy);
    add = MemFree (add);
  }
  return add;
}


static void 
AddOneFeatureBetween 
(BioseqPtr bsp,
 Int4 feat_left,
 Boolean partial_left,
 Int4 feat_right,
 Boolean partial_right,
 Uint1   strand,
 Uint1   featdef,
 ApplyFeatureDetailsPtr details)
{
  Int4       seqfeattype;
  SeqLocPtr  slp;

  seqfeattype = FindFeatFromFeatDefType (featdef);

  slp = SeqLocIntNew (feat_left, feat_right, strand, SeqIdFindWorst (bsp->id));
  SetSeqLocPartial (slp, partial_left, partial_right);
  
  ApplyOneFeatureToBioseq (bsp, featdef, slp, details == NULL ? NULL : details->fields, details == NULL ? NULL : details->src_fields, details->add_mrna);
}


static void AddFeaturesBetweenCallback (BioseqPtr bsp, Pointer data)
{
  AddFeatureBetweenPtr add;
  SeqFeatPtr           feat_before, feat_after;
  SeqMgrFeatContext    context_before, context_after;
  Boolean              found = FALSE, partial_left = TRUE, partial_right = TRUE;
  Uint1                featdef_before = FEATDEF_BAD, featdef_after = FEATDEF_BAD, featdef_new;
  Int4                 seqfeattype;
  Int4                 feat_left, feat_right;

  if (bsp == NULL || (add = (AddFeatureBetweenPtr) data) == NULL || add->left_policy == NULL || add->right_policy == NULL) {
    return;
  }
  if (add->feat_type_add == NULL) {
    return;
  }

  featdef_new = GetFeatdefFromFeatureType (add->feat_type_add->choice);
  seqfeattype = FindFeatFromFeatDefType (featdef_new);
  if (seqfeattype == SEQFEAT_PROT) {
    if (!ISA_aa (bsp->mol)) {
      return;
    }
  } else {
    if (ISA_aa (bsp->mol)) {
      return;
    }
  }

  if ((add->left_policy->policy == eFeatureBetweenPolicy_NoFeatPartial 
          || add->left_policy->policy == eFeatureBetweenPolicy_NoFeatNoFeat)
      && add->left_policy->feat_type != NULL) {
    featdef_before = GetFeatdefFromFeatureType (add->left_policy->feat_type->choice);
  }
  if ((add->right_policy->policy == eFeatureBetweenPolicy_NoFeatPartial 
          || add->right_policy->policy == eFeatureBetweenPolicy_NoFeatNoFeat)
      && add->right_policy->feat_type != NULL) {
    featdef_after = GetFeatdefFromFeatureType (add->right_policy->feat_type->choice);
  }

  if (featdef_before != FEATDEF_BAD) {
    for (feat_before = SeqMgrGetNextFeature (bsp, NULL, 0, featdef_before, &context_before);
         feat_before != NULL;
         feat_before = SeqMgrGetNextFeature (bsp, feat_before, 0, featdef_before, &context_before)) {
      if (DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, feat_before, add->left_policy->constraint)) {
        /* add feature starting at context_before.right + 1, ending with first feature matching right constraint */
        feat_left = context_before.right + 1;
        partial_left = FALSE;
        feat_right = bsp->length - 1;
        if (featdef_after != FEATDEF_BAD) {
          for (feat_after = SeqMgrGetNextFeature (bsp, NULL, 0, featdef_after, &context_after);
               feat_after != NULL && feat_right == bsp->length - 1;
               feat_after = SeqMgrGetNextFeature (bsp, feat_after, 0, featdef_after, &context_after)) {
            if (context_after.left > feat_left && DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, feat_after, add->right_policy->constraint)) {
              feat_right = context_after.left - 1;
              partial_right = FALSE;
            }
          }
        }
        if (feat_right > -1
            && feat_right > feat_left
            && (feat_right < bsp->length - 1 
                || add->right_policy->policy != eFeatureBetweenPolicy_NoFeatNoFeat)) {
          if (add->right_policy->policy == eFeatureBetweenPolicy_ExtendComplete) {
            partial_right = FALSE;
          }
          AddOneFeatureBetween (bsp, feat_left, partial_left, feat_right, partial_right, context_before.strand, featdef_new, add->details);
        }
        found = TRUE;
      }
    }
    if (!found && add->left_policy->policy == eFeatureBetweenPolicy_NoFeatPartial) {
      /* add feature starting at 0 and partial on this end, ending with first feature matching right constraint */
      feat_left = 0;
      feat_right = bsp->length - 1;
      if (featdef_after != FEATDEF_BAD) {
        for (feat_after = SeqMgrGetNextFeature (bsp, NULL, 0, featdef_after, &context_after);
              feat_after != NULL && feat_right == bsp->length - 1;
              feat_after = SeqMgrGetNextFeature (bsp, feat_after, 0, featdef_after, &context_after)) {
          if (context_after.left > feat_left && DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, feat_after, add->right_policy->constraint)) {
            feat_right = context_after.left - 1;
            partial_right = FALSE;
            AddOneFeatureBetween (bsp, feat_left, partial_left, feat_right, partial_right, context_after.strand, featdef_new, add->details);
            found = TRUE;
          }
        }
      }
      if (!found && add->right_policy->policy != eFeatureBetweenPolicy_NoFeatNoFeat) {
        if (add->right_policy->policy == eFeatureBetweenPolicy_ExtendComplete) {
          partial_right = FALSE;
        }
        if (feat_right > feat_left && feat_right > -1) {
          AddOneFeatureBetween (bsp, feat_left, partial_left, feat_right, partial_right, Seq_strand_plus, featdef_new, add->details);
        }
      }
    }
  } else if (featdef_after != FEATDEF_BAD) {
    feat_left = 0;
    if (add->left_policy->policy == eFeatureBetweenPolicy_ExtendComplete) {
      partial_left = FALSE;
    }
    for (feat_after = SeqMgrGetNextFeature (bsp, NULL, 0, featdef_after, &context_after);
         feat_after != NULL && !found;
         feat_after = SeqMgrGetNextFeature (bsp, feat_after, 0, featdef_after, &context_after)) {
      if (DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, feat_after, add->right_policy->constraint)) {
        /* add feature starting at 0 and partial on this end, ending with context_after.left - 1 */
        partial_right = FALSE;
        feat_right = context_after.left - 1;
        if (feat_right > feat_left && feat_right > -1) {
          AddOneFeatureBetween (bsp, feat_left, partial_left, feat_right, partial_right, context_after.strand, featdef_new, add->details);
        }
        found = TRUE;
      }
    }
    if (!found && add->right_policy->policy != eFeatureBetweenPolicy_NoFeatNoFeat) {
      AddOneFeatureBetween (bsp, feat_left, partial_left, bsp->length - 1, partial_right, Seq_strand_plus, featdef_new, add->details);      
    }
  } else {
    feat_left = 0;
    feat_right = bsp->length - 1;
    if (add->left_policy->policy == eFeatureBetweenPolicy_ExtendComplete) {
      partial_left = FALSE;
    }
    if (add->right_policy->policy == eFeatureBetweenPolicy_ExtendComplete) {
      partial_right = FALSE;
    }
    if (feat_right > feat_left && feat_right > -1) {
      AddOneFeatureBetween (bsp, feat_left, partial_left, feat_right, partial_right, Seq_strand_plus, featdef_new, add->details);          
    }
  }
}

typedef struct featurebetweenendpolicydlg {
  DIALOG_MESSAGE_BLOCK

  GrouP          missing_feat_policy;
  DialoG         feat_type_dlg;
  DialoG         constraint_dlg;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} FeatureBetweenEndPolicyDlgData, PNTR FeatureBetweenEndPolicyDlgPtr;


static void EnableDialogFeatureBetweenEndPolicy (FeatureBetweenEndPolicyDlgPtr dlg)
{
  Int4 val;

  if (dlg != NULL) {
    val = GetValue (dlg->missing_feat_policy);
    if (val == eFeatureBetweenPolicy_NoFeatPartial || val == eFeatureBetweenPolicy_NoFeatNoFeat) {
      Enable (dlg->feat_type_dlg);
      Enable (dlg->constraint_dlg);
    } else {
      Disable (dlg->feat_type_dlg);
      Disable (dlg->constraint_dlg);
    }
  }
}


static void FeatureBetweenEndPolicyToDialog (DialoG d, Pointer data)
{
  FeatureBetweenEndPolicyDlgPtr dlg;
  FeatureBetweenEndPolicyPtr    pol;
  ValNode                       vn;

  dlg = (FeatureBetweenEndPolicyDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  pol = (FeatureBetweenEndPolicyPtr) data;
  if (data == NULL) {
    SetValue (dlg->missing_feat_policy, eFeatureBetweenPolicy_NoFeatPartial);
    MemSet (&vn, 0, sizeof (ValNode));
    vn.choice = Macro_feature_type_cds;
    PointerToDialog (dlg->feat_type_dlg, &vn);
    PointerToDialog (dlg->constraint_dlg, NULL);
  } else {
    SetValue (dlg->missing_feat_policy, pol->policy);
    PointerToDialog (dlg->feat_type_dlg, pol->feat_type);
    PointerToDialog (dlg->constraint_dlg, pol->constraint);
  }
  EnableDialogFeatureBetweenEndPolicy (dlg);
}


static Pointer DialogToFeatureBetweenEndPolicy (DialoG d)
{
  FeatureBetweenEndPolicyDlgPtr dlg;
  FeatureBetweenEndPolicyPtr    pol;

  dlg = (FeatureBetweenEndPolicyDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  pol = FeatureBetweenEndPolicyNew();
  pol->policy = GetValue (dlg->missing_feat_policy);
  pol->feat_type = DialogToPointer (dlg->feat_type_dlg);
  pol->constraint = DialogToPointer (dlg->constraint_dlg);

  return pol;
}


static void ChangeFeatureBetweenPolicy (GrouP g)
{
  EnableDialogFeatureBetweenEndPolicy (GetObjectExtra (g));
}


static DialoG FeatureBetweenEndPolicyDialog (GrouP h, CharPtr label, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{  
  FeatureBetweenEndPolicyDlgPtr dlg;
  GrouP                    p;
  PrompT                   ppt;

  dlg = (FeatureBetweenEndPolicyDlgPtr) MemNew (sizeof (FeatureBetweenEndPolicyDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FeatureBetweenEndPolicyToDialog;
  dlg->fromdialog = DialogToFeatureBetweenEndPolicy;
  dlg->dialogmessage = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  ppt = StaticPrompt (p, label, 0, popupMenuHeight, programFont, 'l');
  dlg->missing_feat_policy = HiddenGroup (p, 0, 4, ChangeFeatureBetweenPolicy);
  SetObjectExtra (dlg->missing_feat_policy, dlg, NULL);
  RadioButton (dlg->missing_feat_policy, "If selected feature missing, extend to partial end");
  RadioButton (dlg->missing_feat_policy, "If selected feature missing, do not create feature");
  RadioButton (dlg->missing_feat_policy, "Extend to complete end");
  RadioButton (dlg->missing_feat_policy, "Extend to partial end");

  dlg->feat_type_dlg = FeatureTypeDialog (p, NULL, NULL);
  dlg->constraint_dlg = ConstraintSetDialog (p, NULL, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->missing_feat_policy, 
                              (HANDLE) dlg->feat_type_dlg, (HANDLE) dlg->constraint_dlg, NULL);

  return (DialoG) p;
}


typedef struct addfeaturebetweendlg {
  DIALOG_MESSAGE_BLOCK

  DialoG         left_policy_dlg;

  DialoG         feat_type_add_dlg;
  DialoG         details_dlg;

  DialoG         right_policy_dlg;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} AddFeatureBetweenDlgData, PNTR AddFeatureBetweenDlgPtr;


static Pointer DialogToAddFeatureBetween (DialoG d)
{
  AddFeatureBetweenDlgPtr dlg;
  AddFeatureBetweenPtr add;

  dlg = (AddFeatureBetweenDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  add = AddFeatureBetweenNew ();
  add->left_policy = DialogToPointer (dlg->left_policy_dlg);
  add->feat_type_add = DialogToPointer (dlg->feat_type_add_dlg);
  add->details = DialogToPointer (dlg->details_dlg);
  add->right_policy = DialogToPointer (dlg->right_policy_dlg);

  return add;
}


static void AddFeatureBetweenToDialog (DialoG d, Pointer data)
{
  AddFeatureBetweenDlgPtr dlg;
  AddFeatureBetweenPtr add;
  ValNode              vn;

  dlg = (AddFeatureBetweenDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  add = (AddFeatureBetweenPtr) data;
  if (add == NULL) {
    PointerToDialog (dlg->left_policy_dlg, NULL);
    MemSet (&vn, 0, sizeof (ValNode));
    vn.choice = Macro_feature_type_misc_feature;
    PointerToDialog (dlg->feat_type_add_dlg, &vn);
    PointerToDialog (dlg->details_dlg, NULL);
    PointerToDialog (dlg->right_policy_dlg, NULL);
  } else {
    PointerToDialog (dlg->left_policy_dlg, add->left_policy);
    PointerToDialog (dlg->feat_type_add_dlg, add->feat_type_add);
    PointerToDialog (dlg->details_dlg, add->details);
    PointerToDialog (dlg->right_policy_dlg, add->right_policy);
  }

}


static DialoG AddFeatureBetweenDialog (GrouP h, Uint1 add_type, ApplyFeatureDetailsPtr details, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  AddFeatureBetweenDlgPtr dlg;
  GrouP                    p, k2;
  PrompT                   ppt;

  dlg = (AddFeatureBetweenDlgPtr) MemNew (sizeof (AddFeatureBetweenDlgData));
  p = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = AddFeatureBetweenToDialog;
  dlg->fromdialog = DialogToAddFeatureBetween;
  dlg->dialogmessage = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->left_policy_dlg = FeatureBetweenEndPolicyDialog (p, "Left Feature", NULL, NULL);

  k2 = HiddenGroup (p, -1, 0, NULL);
  ppt = StaticPrompt (k2, "Feature to Add", 0, popupMenuHeight, programFont, 'l');
  dlg->feat_type_add_dlg = FeatureTypeDialog (k2, change_notify, change_userdata);
  dlg->details_dlg = ApplyFeatureDetailsDialog (k2, add_type, details, TRUE, NULL, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->feat_type_add_dlg, (HANDLE) dlg->details_dlg, NULL);

  dlg->right_policy_dlg = FeatureBetweenEndPolicyDialog (p, "Right Feature", NULL, NULL);

  return (DialoG) p;
}


typedef struct addfeaturebetweenform {
  FORM_MESSAGE_BLOCK

  DialoG add_between_dlg;
  ButtoN leave_dlg_up;

  Boolean changing;
} AddFeatureBetweenFormData, PNTR AddFeatureBetweenFormPtr;


static void ChangeAddFeatureFunc (Pointer data);

static void DoAddFeatureBetween (ButtoN b)
{
  AddFeatureBetweenFormPtr dlg;
  AddFeatureBetweenPtr     add;
  SeqEntryPtr              sep;

  dlg = (AddFeatureBetweenFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  WatchCursor();
  Update();
  add = DialogToPointer (dlg->add_between_dlg);
  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);

  VisitBioseqsInSep (sep, add, AddFeaturesBetweenCallback);

  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  
  if (!GetStatus (dlg->leave_dlg_up)) {
    Remove (dlg->form);
  }
  ArrowCursor();
  Update();
}


static WindoW AddFeatureBetweenForm (Uint2 entityID, AddFeatureBetweenPtr add)
{
  WindoW            w;
  GrouP             h, c;
  ButtoN            b;
  AddFeatureBetweenFormPtr dlg;

  dlg = (AddFeatureBetweenFormPtr) MemNew (sizeof (AddFeatureBetweenFormData));
  if (dlg == NULL) return NULL;
    
  w = FixedWindow (-50, -33, -10, -10, "Add Feature Between", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  dlg->changing = TRUE;
  dlg->add_between_dlg = AddFeatureBetweenDialog (h, (add == NULL || add->feat_type_add == NULL) ? FEATDEF_misc_feature : GetFeatdefFromFeatureType (add->feat_type_add->choice),
                                                  add == NULL ? NULL : add->details,
                                                  ChangeAddFeatureFunc, dlg);
  PointerToDialog (dlg->add_between_dlg, add);
  dlg->changing = FALSE;

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoAddFeatureBetween);
  SetObjectExtra (b, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  dlg->leave_dlg_up = CheckBox (c, "Leave dialog up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->add_between_dlg, (HANDLE) c, NULL);

  return w;
}


static void ChangeAddFeatureFunc (Pointer data)
{
  AddFeatureBetweenFormPtr dlg;
  AddFeatureBetweenPtr add;
  WindoW               w;

  dlg = (AddFeatureBetweenFormPtr) data;
  if (dlg == NULL || dlg->changing) {
    return;
  }

  add = DialogToPointer (dlg->add_between_dlg);
  Hide (dlg->form);
  /* build new form */
  w = AddFeatureBetweenForm (dlg->input_entityID, add);
  RealizeWindow (w);
  Show (w);
  Update ();

  /* remove old form */
  dlg->changing = TRUE;
  Remove (dlg->form);
  
  add = AddFeatureBetweenFree (add);

}


NLM_EXTERN void AddFeatureBetween (IteM i)
{
  BaseFormPtr       bfp;
  WindoW            w;
  AddFeatureBetweenFormPtr dlg;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  dlg = (AddFeatureBetweenFormPtr) MemNew (sizeof (AddFeatureBetweenFormData));
  if (dlg == NULL) return;

  w = AddFeatureBetweenForm (bfp->input_entityID, NULL);
    
  RealizeWindow (w);
  Show (w);
  Update ();
}


typedef struct changequalcaseform {
  FORM_MESSAGE_BLOCK

  DialoG qual_dlg;
  DialoG cap_dlg;
  ButtoN leave_dlg_up;

} ChangeQualCaseFormData, PNTR ChangeQualCaseFormPtr;


typedef struct changecase {
  ValNodePtr dst_field;
  Int4       change_cap;
  ValNodePtr orgnames;
} ChangeCaseData, PNTR ChangeCasePtr;

static void ChangeQualCaseCallbackDesc (SeqDescrPtr sdp, Pointer data)
{
  ChangeCasePtr ccp;
  ParseDstOrgPtr o;
  CharPtr        str;

  ccp = (ChangeCasePtr) data;
  if (sdp == NULL || ccp == NULL || ccp->dst_field == NULL || ccp->change_cap < 1) {
    return;
  }

  switch (ccp->dst_field->choice) {
    case ParseDest_defline:
      if (sdp->choice == Seq_descr_title) {
        FixCapitalizationInString ((CharPtr PNTR)&(sdp->data.ptrvalue), ccp->change_cap, ccp->orgnames);
      }
      break;
    case ParseDest_org:
      if ((o = (ParseDstOrgPtr) ccp->dst_field->data.ptrvalue) != NULL
          && o->type != Object_type_constraint_feature
          && sdp->choice == Seq_descr_source) {
        str = GetSourceQualFromBioSource (sdp->data.ptrvalue, o->field, NULL);
        FixCapitalizationInString(&str, ccp->change_cap, NULL);
        SetSourceQualInBioSource(sdp->data.ptrvalue, o->field, NULL, str, ExistingTextOption_replace_old);
        str = MemFree (str);
      }
      break;
    case ParseDest_comment_descriptor:
      if (sdp->choice == Seq_descr_comment) {
        FixCapitalizationInString ((CharPtr PNTR)&(sdp->data.ptrvalue), ccp->change_cap, ccp->orgnames);
      }
      break;
    case ParseDest_dbxref:
      if (sdp->choice == Seq_descr_source) {
        str = GetDBxrefFromBioSource (sdp->data.ptrvalue, ccp->dst_field->data.ptrvalue);
        if (str != NULL) {
          FixCapitalizationInString(&str, ccp->change_cap, ccp->orgnames);
          SetDBxrefForBioSource (sdp->data.ptrvalue, ccp->dst_field->data.ptrvalue, str, ExistingTextOption_replace_old);
          str = MemFree (str);
        }
      }
      break;
  }
}


static void ChangeQualCaseCallbackFeat (SeqFeatPtr sfp, Pointer data)
{
  ChangeCasePtr ccp;
  ParseDstOrgPtr o;
  CharPtr        str;
  FeatureFieldLegalPtr fl;
  FeatureField f;

  ccp = (ChangeCasePtr) data;
  if (sfp == NULL || ccp == NULL || ccp->dst_field == NULL || ccp->change_cap < 1) {
    return;
  }

  switch (ccp->dst_field->choice) {
    case ParseDest_org:
      if ((o = (ParseDstOrgPtr) ccp->dst_field->data.ptrvalue) != NULL
          && o->type != Object_type_constraint_descriptor
          && sfp->data.choice == SEQFEAT_BIOSRC) {
        str = GetSourceQualFromBioSource (sfp->data.value.ptrvalue, o->field, NULL);
        FixCapitalizationInString(&str, ccp->change_cap, NULL);
        SetSourceQualInBioSource(sfp->data.value.ptrvalue, o->field, NULL, str, ExistingTextOption_replace_old);
        str = MemFree (str);
      }
      break;
    case ParseDest_featqual:
      fl = (FeatureFieldLegalPtr) ccp->dst_field->data.ptrvalue;
      if (fl != NULL) {
        f.type = fl->type;
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = fl->field;        
        str = GetQualFromFeature (sfp, &f, NULL);
        if (str != NULL) {
          FixCapitalizationInString(&str, ccp->change_cap, ccp->orgnames);
          SetQualOnFeature (sfp, &f, NULL, str, ExistingTextOption_replace_old);
          f.field = ValNodeFree (f.field);
          str = MemFree (str);
        }
      }
      break;
    case ParseDest_dbxref:
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        str = GetDBxrefFromBioSource (sfp->data.value.ptrvalue, ccp->dst_field->data.ptrvalue);
        if (str != NULL) {
          FixCapitalizationInString(&str, ccp->change_cap, ccp->orgnames);
          SetDBxrefForBioSource (sfp->data.value.ptrvalue, ccp->dst_field->data.ptrvalue, str, ExistingTextOption_replace_old);
          str = MemFree (str);
        }
      }
      break;
  }
}


static void DoChangeQualCase(ButtoN b)
{
  ChangeQualCaseFormPtr dlg;
  ChangeCaseData        ccd;
  SeqEntryPtr           sep;
  ParseDstOrgPtr        o;

  dlg = (ChangeQualCaseFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  ccd.dst_field = DialogToPointer (dlg->qual_dlg);
  if (ccd.dst_field == NULL) {
    return;
  }

  ccd.change_cap = GetCapChangeDialogValue (dlg->cap_dlg);

  ccd.orgnames = NULL;
  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  VisitBioSourcesInSep (sep, &ccd.orgnames, GetOrgNamesInRecordCallback);


  switch (ccd.dst_field->choice) {
    case ParseDest_defline :
      VisitDescriptorsInSep (sep, &ccd, ChangeQualCaseCallbackDesc);
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) ccd.dst_field->data.ptrvalue;
      if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_descriptor) {
        VisitDescriptorsInSep (sep, &ccd, ChangeQualCaseCallbackDesc);
      }
      if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_feature) {
        VisitFeaturesInSep (sep, &ccd, ChangeQualCaseCallbackFeat);
      }
      break;
    case ParseDest_featqual :
      VisitFeaturesInSep (sep, &ccd, ChangeQualCaseCallbackFeat);
      break;
    case ParseDest_comment_descriptor :
      VisitDescriptorsInSep (sep, &ccd, ChangeQualCaseCallbackDesc);
      break;
    case ParseDest_dbxref :
      VisitDescriptorsInSep (sep, &ccd, ChangeQualCaseCallbackDesc);
      VisitFeaturesInSep (sep, &ccd, ChangeQualCaseCallbackFeat);
      break;
  }

  ccd.orgnames = ValNodeFree (ccd.orgnames);
  ccd.dst_field = ParseDestFree (ccd.dst_field);

  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);

  if (!GetStatus (dlg->leave_dlg_up)) {
    Remove (dlg->form);
  }
  ArrowCursor();
  Update();
}


NLM_EXTERN void ChangeQualifierCaseBaseForm (BaseFormPtr bfp)
{
  WindoW            w;
  GrouP             h, c;
  ButtoN            b;
  ChangeQualCaseFormPtr dlg;
  ValNode           vn;

  if (bfp == NULL || bfp->input_entityID == 0) return;

  dlg = (ChangeQualCaseFormPtr) MemNew (sizeof (ChangeQualCaseFormData));
  if (dlg == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, "Lowercase Qualifiers", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  dlg->qual_dlg = ParseDstDialog (h, NULL, NULL);
  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = ParseDest_org;
  PointerToDialog (dlg->qual_dlg, &vn);
  dlg->cap_dlg = CapChangeDialog (h, NULL, NULL);
  SetCapChangeDialogValue (dlg->cap_dlg, Cap_change_tolower);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoChangeQualCase);
  SetObjectExtra (b, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  dlg->leave_dlg_up = CheckBox (c, "Leave dialog up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->qual_dlg, (HANDLE) dlg->cap_dlg, (HANDLE) c, NULL);
    
  RealizeWindow (w);
  Show (w);
  Update ();
}

NLM_EXTERN void ChangeQualifierCase (IteM i)
{
  BaseFormPtr       bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  ChangeQualifierCaseBaseForm (bfp);
}


typedef struct applymoltypeform {
  FORM_MESSAGE_BLOCK

  DialoG mol_type_dlg;

  Uint1 mol_type;
  Uint1 mol_class;

} ApplyMolTypeFormData, PNTR ApplyMolTypeFormPtr;


static void ApplyMoleculeTypeCallback (BioseqPtr bsp, Pointer data)
{
  ApplyMolTypeFormPtr dlg;
  SeqDescrPtr         sdp;
  SeqMgrDescContext   context;
  MolInfoPtr          mip;

  if (bsp == NULL || ISA_aa(bsp->mol) || (dlg = (ApplyMolTypeFormPtr) data) == NULL) {
    return;
  }

  bsp->mol = dlg->mol_class;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp == NULL) {
    sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_molinfo);
  }

  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) {
    mip = MolInfoNew ();
    sdp->data.ptrvalue = mip;
  }
  mip->biomol = dlg->mol_type;
}


static void DoApplyMoleculeType(ButtoN b)
{
  ApplyMolTypeFormPtr dlg;
  ValNodePtr          vnp;
  SeqEntryPtr         sep;

  dlg = (ApplyMolTypeFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);

  vnp = DialogToPointer (dlg->mol_type_dlg);
  if (vnp == NULL) {
    return;
  }

  dlg->mol_type = BiomolFromMoleculeType (vnp->choice);
  if (StringSearch (vnp->data.ptrvalue, "RNA") == NULL) {
    dlg->mol_class = MOLECULE_CLASS_DNA;
  } else {
    dlg->mol_class = MOLECULE_CLASS_RNA;
  }

  VisitBioseqsInSep (sep, dlg, ApplyMoleculeTypeCallback);

  vnp = ValNodeFree (vnp);
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  Remove (dlg->form);
  ArrowCursor();
  Update();
}


NLM_EXTERN void ExternalApplyMoleculeType (IteM i)
{
  BaseFormPtr       bfp;
  WindoW            w;
  GrouP             h, c;
  ButtoN            b;
  ApplyMolTypeFormPtr dlg;
  ValNodePtr          moltype_list, vnp, vnp_new;


#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  dlg = (ApplyMolTypeFormPtr) MemNew (sizeof (ApplyMolTypeFormData));
  if (dlg == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, "Apply Molecule Type", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  moltype_list = GetMoleculeTypeList ();
  /* edit list - split genomic into genomic DNA and genomic RNA */
  for (vnp = moltype_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Molecule_type_unknown) {
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = StringSave ("unknown");
    } else if (vnp->choice == Molecule_type_genomic) {
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = StringSave ("genomic DNA");
      vnp_new = ValNodeNew (NULL);
      vnp_new->choice = Molecule_type_genomic;
      vnp_new->data.ptrvalue = StringSave ("genomic RNA");
      vnp_new->next = vnp->next;
      vnp->next = vnp_new;
      vnp = vnp_new;
    }
  }

  dlg->mol_type_dlg = ValNodeSelectionDialog (h, moltype_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                              ValNodeSimpleDataFree, ValNodeStringCopy,
                                              ValNodeStringMatch, "molecule type", 
                                              NULL, NULL, FALSE);
  vnp = ValNodeNew (NULL);
  vnp->choice = Molecule_type_unknown;
  vnp->data.ptrvalue = StringSave ("unknown");
  PointerToDialog (dlg->mol_type_dlg, vnp);
  vnp = ValNodeFree (vnp);
                                              
  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoApplyMoleculeType);
  SetObjectExtra (b, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->mol_type_dlg, (HANDLE) c, NULL);
    
  RealizeWindow (w);
  Show (w);
  Update ();
}


typedef struct formattedstruccomform {
  FORM_MESSAGE_BLOCK

  TexT * values;
  Int4   num_fields;

  CommentRulePtr    comment_rule;

} FormattedStrucComFormData, PNTR FormattedStrucComFormPtr;


static void CleanupFormattedStrucComFormProc (GraphiC g, VoidPtr data)

{
  FormattedStrucComFormPtr  dlg;

  if (data) {
    dlg = (FormattedStrucComFormPtr) data;
    dlg->values = MemFree (dlg->values);
  }
  StdCleanupFormProc (g, data);
}


static void DoApplyFormattedStructuredComment (ButtoN b)
{
  FormattedStrucComFormPtr dlg;
  UserObjectPtr            uop;
  CharPtr                  suffix, value;
  Int4                     len;
  SeqDescrPtr              sdp;
  SeqEntryPtr              sep;
  FieldRulePtr             field;
  Int4                     n;

  dlg = (FormattedStrucComFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  uop = UserObjectNew ();
  uop->type = ObjectIdNew ();
  uop->type->str = StringSave ("StructuredComment");

  AddItemStructuredCommentUserObject (uop, "StructuredCommentPrefix", dlg->comment_rule->prefix);

  for (field = dlg->comment_rule->fields, n = 0; field != NULL; field = field->next, n++) {
    if (!TextHasNoText (dlg->values[n])) {
      value = SaveStringFromText (dlg->values[n]);
      AddItemStructuredCommentUserObject (uop, field->field_name, value);
      value = MemFree (value);
    }
  }

  if ((len = StringLen (dlg->comment_rule->prefix)) > 7 
    && StringCmp (dlg->comment_rule->prefix + len - 7, "START##") == 0) {      
    suffix = StringSave (dlg->comment_rule->prefix);
    sprintf (suffix + len - 7, "%s", "END##");
    AddItemStructuredCommentUserObject (uop, "StructuredCommentSuffix", suffix);
    suffix = MemFree (suffix);
  }

  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  sdp = CreateNewDescriptor (sep, Seq_descr_user);
  sdp->data.ptrvalue = uop;

  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  Remove (dlg->form);
  Update ();
}

static void ApplyFormattedStructuredComment (IteM i, CharPtr prefix)
{
  BaseFormPtr       bfp;
  CharPtr           prefix_fmt = "##%s-Data-START##";
  CharPtr           title_fmt = "Apply %s Structured Comment";
  CharPtr           tmp;
  WindoW            w;
  GrouP             h, g, c;
  ButtoN            b;
  CommentRulePtr    comment_rule;
  FieldRulePtr      field;
  Int4              n;
  FormattedStrucComFormPtr dlg;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (prefix) + StringLen (prefix_fmt)));
  sprintf (tmp, prefix_fmt, prefix);

  comment_rule = GetCommentRuleFromRuleSet (tmp);
  tmp = MemFree (tmp);
  if (comment_rule == NULL) {
    return;
  }

  dlg = (FormattedStrucComFormPtr) MemNew (sizeof (FormattedStrucComFormData));
  if (dlg == NULL) return;

  dlg->comment_rule = comment_rule;
  n = 0;
  for (field = comment_rule->fields; field != NULL; field = field->next) {
    n++;
  }
  dlg->num_fields = n;
  dlg->values = (TexT *) MemNew (sizeof (TexT) * dlg->num_fields);

  tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen(prefix)));
  sprintf (tmp, title_fmt, prefix);
  w = FixedWindow (-50, -33, -10, -10, tmp, StdCloseWindowProc);
  tmp = MemFree (tmp);
  SetObjectExtra (w, dlg, CleanupFormattedStrucComFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);
  for (field = dlg->comment_rule->fields, n = 0; field != NULL; field = field->next, n++) {
    StaticPrompt (g, field->field_name, 0, dialogTextHeight, programFont, 'l');
    dlg->values[n] = DialogText (g, "", 10, NULL);
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoApplyFormattedStructuredComment);
  SetObjectExtra (b, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE)g, (HANDLE) c, NULL);
    
  RealizeWindow (w);
  Show (w);
  Update ();
}


NLM_EXTERN void ApplyGenomeAssemblyComment (IteM i)
{
  ApplyFormattedStructuredComment (i, "Genome-Assembly");
}


NLM_EXTERN void ApplyAssemblyComment (IteM i)
{
  ApplyFormattedStructuredComment (i, "Assembly");
}


typedef struct releasedateform {
  FORM_MESSAGE_BLOCK

  DialoG dlg;
} ReleaseDateFormData, PNTR ReleaseDateFormPtr;


static void DoApplyReleaseDate (ButtoN b)
{
  ReleaseDateFormPtr dlg;
  SeqEntryPtr        sep;
  SeqSubmitPtr       ssp;

  dlg = (ReleaseDateFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  if (sep == NULL) return;

  ssp = FindSeqSubmitForSeqEntry (sep);
  if (ssp == NULL) {
    return;
  }

  ssp->sub = SubmitBlockFree (ssp->sub);
  ssp->sub = DialogToPointer (dlg->dlg);

  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  Remove (dlg->form);
  Update ();
}


NLM_EXTERN void SetReleaseDate(IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  SeqSubmitPtr ssp;
  ReleaseDateFormPtr dlg;
  WindoW             w;
  GrouP              h, c;
  ButtoN             b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  ssp = FindSeqSubmitForSeqEntry (sep);
  if (ssp == NULL) {
    return;
  }

  dlg = (ReleaseDateFormPtr) MemNew (sizeof (ReleaseDateFormData));
  if (dlg == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, "Set Release Date", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  dlg->dlg = CreateSubmitDataDialog (h, "", ssp->sub == NULL, ssp->sub == NULL || ssp->sub->hup);
  PointerToDialog (dlg->dlg, ssp->sub);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoApplyReleaseDate);
  SetObjectExtra (b, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE)dlg->dlg, (HANDLE) c, NULL);
    
  RealizeWindow (w);
  Show (w);
  Update ();
}


typedef struct removetextform {
  FORM_MESSAGE_BLOCK
  DialoG text_portion;
  GrouP  not_found;
  DialoG field;
  DialoG constraint;
  ButtoN leave_dlg_up;
  ButtoN accept_btn;

} RemoveTextFormData, PNTR RemoveTextFormPtr;

typedef struct removetext {
  TextPortionPtr text_portion;
  Boolean    remove_if_not_found;
  ValNodePtr field;
  ValNodePtr constraint;
} RemoveTextData, PNTR RemoveTextPtr;


static RemoveTextPtr RemoveTextFree (RemoveTextPtr data)
{
  if (data != NULL) {
    data->text_portion = TextPortionFree (data->text_portion);
    data->field = ParseDestFree (data->field);
    data->constraint = ConstraintChoiceSetFree (data->constraint);
    data = MemFree (data);
  }
  return data;
}


static void RemoveTextCallbackDesc (SeqDescPtr sdp, Pointer data)
{
  RemoveTextPtr     rt;
  ParseDstOrgPtr    o;
  CharPtr           str;
  ObjValNodePtr     ovn;
  StringConstraint  sc;

  if (sdp == NULL || (rt = (RemoveTextPtr) data) == NULL
      || !DoesObjectMatchConstraintChoiceSet(OBJ_SEQDESC, sdp, rt->constraint)) {
    return;
  }
  switch (rt->field->choice) {
    case ParseDest_defline:
      if (sdp->choice == Seq_descr_title) {
        if (!RemoveTextPortionFromString ((CharPtr)(sdp->data.ptrvalue), rt->text_portion)
            && rt->remove_if_not_found) {
          if (sdp->extended) {
            ovn = (ObjValNodePtr) sdp;
            ovn->idx.deleteme = TRUE;
          }
        }
      }
      break;
    case ParseDest_org:
      if ((o = (ParseDstOrgPtr) rt->field->data.ptrvalue) != NULL
          && o->type != Object_type_constraint_feature
          && sdp->choice == Seq_descr_source) {
        str = GetSourceQualFromBioSource (sdp->data.ptrvalue, o->field, NULL);
        if (RemoveTextPortionFromString(str, rt->text_portion)) {
          SetSourceQualInBioSource(sdp->data.ptrvalue, o->field, NULL, str, ExistingTextOption_replace_old);
        } else if (rt->remove_if_not_found) {
          MemSet (&sc, 0, sizeof (StringConstraint));
          sc.case_sensitive = TRUE;
          sc.match_location = String_location_equals;
          sc.match_text = str;
          RemoveSourceQualFromBioSource (sdp->data.ptrvalue, o->field, &sc);
        }
        str = MemFree (str);
      }
      break;
    case ParseDest_comment_descriptor:
      if (sdp->choice == Seq_descr_comment) {
        if (!RemoveTextPortionFromString ((CharPtr)(sdp->data.ptrvalue), rt->text_portion)
          && rt->remove_if_not_found) {
          if (sdp->extended) {
            ovn = (ObjValNodePtr) sdp;
            ovn->idx.deleteme = TRUE;
          }
        }
      }
      break;
    case ParseDest_dbxref:
      if (sdp->choice == Seq_descr_source) {
        str = GetDBxrefFromBioSource (sdp->data.ptrvalue, rt->field->data.ptrvalue);
        if (str != NULL) {
          if (RemoveTextPortionFromString(str, rt->text_portion)) {
            SetDBxrefForBioSource (sdp->data.ptrvalue, rt->field->data.ptrvalue, str, ExistingTextOption_replace_old);
          } else if (rt->remove_if_not_found) {
            MemSet (&sc, 0, sizeof (StringConstraint));
            sc.case_sensitive = TRUE;
            sc.match_location = String_location_equals;
            sc.match_text = str;
            RemoveDBxrefForBioSource (sdp->data.ptrvalue, rt->field->data.ptrvalue, &sc);
          }
          str = MemFree (str);
        }
      }
      break;
  }
}


static void RemoveTextCallbackFeat (SeqFeatPtr sfp, Pointer data)
{
  RemoveTextPtr     rt;
  ParseDstOrgPtr    o;
  CharPtr           str;
  StringConstraint  sc;
  FeatureFieldLegalPtr fl;
  FeatureField         f;

  if (sfp == NULL || (rt = (RemoveTextPtr) data) == NULL
      || !DoesObjectMatchConstraintChoiceSet(OBJ_SEQFEAT, sfp, rt->constraint)) {
    return;
  }

  switch (rt->field->choice) {
    case ParseDest_org:
      if ((o = (ParseDstOrgPtr) rt->field->data.ptrvalue) != NULL
          && o->type != Object_type_constraint_descriptor
          && sfp->data.choice == SEQFEAT_BIOSRC) {
        str = GetSourceQualFromBioSource (sfp->data.value.ptrvalue, o->field, NULL);
        if (RemoveTextPortionFromString(str, rt->text_portion)) {
          SetSourceQualInBioSource(sfp->data.value.ptrvalue, o->field, NULL, str, ExistingTextOption_replace_old);
        } else if (rt->remove_if_not_found) {
          MemSet (&sc, 0, sizeof (StringConstraint));
          sc.case_sensitive = TRUE;
          sc.match_location = String_location_equals;
          sc.match_text = str;
          RemoveSourceQualFromBioSource (sfp->data.value.ptrvalue, o->field, &sc);
        }          
        str = MemFree (str);
      }
      break;
    case ParseDest_featqual:
      fl = (FeatureFieldLegalPtr) rt->field->data.ptrvalue;
      if (fl != NULL) {
        f.type = fl->type;
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = fl->field;        
        str = GetQualFromFeature (sfp, &f, NULL);
        if (str != NULL) {
          if (RemoveTextPortionFromString(str, rt->text_portion)) {
            SetQualOnFeature (sfp, &f, NULL, str, ExistingTextOption_replace_old);
          } else if (rt->remove_if_not_found) {
            MemSet (&sc, 0, sizeof (StringConstraint));
            sc.case_sensitive = TRUE;
            sc.match_location = String_location_equals;
            sc.match_text = str;
            RemoveQualFromFeature (sfp, &f, &sc);
          }
          f.field = ValNodeFree (f.field);
          str = MemFree (str);
        }
      }
      break;
    case ParseDest_dbxref:
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        str = GetDBxrefFromBioSource (sfp->data.value.ptrvalue, rt->field->data.ptrvalue);
        if (str != NULL) {
          if (RemoveTextPortionFromString (str, rt->text_portion)) {
            SetDBxrefForBioSource (sfp->data.value.ptrvalue, rt->field->data.ptrvalue, str, ExistingTextOption_replace_old);
          } else if (rt->remove_if_not_found) {
            MemSet (&sc, 0, sizeof (StringConstraint));
            sc.case_sensitive = TRUE;
            sc.match_location = String_location_equals;
            sc.match_text = str;
            RemoveDBxrefForBioSource (sfp->data.value.ptrvalue, rt->field->data.ptrvalue, &sc);
          }
          str = MemFree (str);
        }
      }
      break;
  }

}


static void DoRemoveText (ButtoN b)
{
  RemoveTextFormPtr dlg;
  RemoveTextPtr     data;
  SeqEntryPtr       sep;
  ParseDstOrgPtr    o;

  dlg = (RemoveTextFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  data = (RemoveTextPtr) MemNew (sizeof (RemoveTextData));

  data->text_portion = DialogToPointer (dlg->text_portion);
  data->field = DialogToPointer (dlg->field);
  data->remove_if_not_found = (GetValue (dlg->not_found) == 2);
  data->constraint = DialogToPointer (dlg->constraint);
  
  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);


  switch (data->field->choice) {
    case ParseDest_defline :
      VisitDescriptorsInSep (sep, data, RemoveTextCallbackDesc);
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) data->field->data.ptrvalue;
      if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_descriptor) {
        VisitDescriptorsInSep (sep, data, RemoveTextCallbackDesc);
      }
      if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_feature) {
        VisitFeaturesInSep (sep, data, RemoveTextCallbackFeat);
      }
      break;
    case ParseDest_featqual :
      VisitFeaturesInSep (sep, data, RemoveTextCallbackFeat);
      break;
    case ParseDest_comment_descriptor :
      VisitDescriptorsInSep (sep, data, RemoveTextCallbackDesc);
      break;
    case ParseDest_dbxref :
      VisitDescriptorsInSep (sep, data, RemoveTextCallbackDesc);
      VisitFeaturesInSep (sep, data, RemoveTextCallbackFeat);
      break;
  }

  data = RemoveTextFree (data);

  DeleteMarkedObjects (dlg->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);
  if (GetStatus (dlg->leave_dlg_up)) {
  } else {
    Remove (dlg->form);
  }
  Update();
}


static void EnableRemoveTextAccept (Pointer data)
{
  RemoveTextFormPtr dlg;
  ValNodePtr        field;

  dlg = (RemoveTextFormPtr) data;
  if (dlg == NULL) {
    return;
  }

  field = DialogToPointer (dlg->field);
  if (field == NULL) {
    Disable (dlg->accept_btn);
  } else {
    Enable (dlg->accept_btn);
  }
  field = ParseDestFree (field);
}


static void ClearRemoveText(ButtoN b)
{
  RemoveTextFormPtr dlg;

  dlg = (RemoveTextFormPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  PointerToDialog (dlg->text_portion, NULL);
  PointerToDialog (dlg->field, NULL);
  EnableRemoveTextAccept(dlg);
}


static void NewRemoveText (IteM i, Boolean inside)
{
  BaseFormPtr  bfp;
  RemoveTextFormPtr dlg;
  WindoW             w;
  GrouP              h, g1, c;
  PrompT             p1;
  ButtoN             clear_btn;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dlg = (RemoveTextFormPtr) MemNew (sizeof (RemoveTextFormData));
  if (dlg == NULL) return;

  w = FixedWindow (-50, -33, -10, -10, inside ? "Remove Text Inside String" : "Remove Text Outside String", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  dlg->text_portion = TextPortionDialog (h, inside, NULL, NULL);

  g1 = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g1, "If field doesn't contain searched for text:", 0, 0, programFont, 'l');

  dlg->not_found = HiddenGroup (g1, 0, 3, NULL);

  RadioButton (dlg->not_found, "     Do nothing to field");
  RadioButton (dlg->not_found, "     Remove entire field text");
  SetValue (dlg->not_found, 1);

  p1 = StaticPrompt (h, "Perform excision in", 0, popupMenuHeight, programFont, 'l');

  dlg->field = ParseDstDialog (h, EnableRemoveTextAccept, dlg);

  /* clear button */
  clear_btn = PushButton (h, "Clear", ClearRemoveText);
  SetObjectExtra (clear_btn, dlg, NULL);

  dlg->constraint = ConstraintSetDialog (h, NULL, NULL);


  c = HiddenGroup (h, 4, 0, NULL);
  dlg->accept_btn = DefaultButton (c, "Accept", DoRemoveText);
  SetObjectExtra (dlg->accept_btn, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  dlg->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE)dlg->text_portion, 
                              (HANDLE) g1, 
                              (HANDLE) p1, 
                              (HANDLE) dlg->field, 
                              (HANDLE) dlg->constraint,
                              (HANDLE) clear_btn, 
                              (HANDLE) c, NULL);
  EnableRemoveTextAccept (dlg);
  RealizeWindow (w);
  Show (w);
  Update ();
}


NLM_EXTERN void NewRemoveTextOutsideString (IteM i)
{
  NewRemoveText (i, FALSE);
}


static void PrintTPAHistErrors (LogInfoPtr lip, ValNodePtr errors)
{
  ValNodePtr vnp;

  if (lip == NULL || lip->fp == NULL || errors == NULL) return;

  for (vnp = errors; vnp != NULL; vnp = vnp->next)
  {
    fprintf (lip->fp, "%s\n", (CharPtr) vnp->data.ptrvalue);
    lip->data_in_log = TRUE;
  }
  fprintf (lip->fp, "\n\n");
}


NLM_EXTERN void ImportTpaAlignment (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  BioseqPtr    bsp, bsp_tmp;
  Char         path [PATH_MAX];
  AsnIoPtr     aip;
  SeqAlignPtr  align, last_align, next_align;
  SeqAnnotPtr  sap;
  SeqHistPtr   shp;
  ValNodePtr   err_list = NULL;
  LogInfoPtr   lip;
  UserObjectPtr uop;
  SeqDescPtr   sdp;
  Char         id_buf[255];
  SeqIdPtr     sip;
  Int4         primary_start, primary_stop;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) {
    return;
  }

  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) {
    Message (MSG_ERROR, "Need to select sequence to attach alignment to");
    return;
  }

  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;

  aip = AsnIoOpen (path, "rb");
  align = SeqAlignAsnRead (aip, NULL);

  if (align == NULL) {
    AsnIoClose (aip);
    aip = AsnIoOpen (path, "rb");
    sap = SeqAnnotAsnRead (aip, NULL);
    if (sap == NULL) {
      Message (MSG_ERROR, "Unable to read alignment");
      return;
    } else if (sap->type != 2 || sap->data == NULL) {
      Message (MSG_ERROR, "Unable to read alignment");
      sap = SeqAnnotFree (sap);
      return;
    } else {
      align = sap->data;
      sap->data = NULL;
      sap = SeqAnnotFree (sap);
    }
  } else {
    last_align = align;
    while ((next_align = SeqAlignAsnRead (aip, NULL)) != NULL) {
      last_align->next = next_align;
      last_align = next_align;
    }
  }
  AsnIoClose (aip);

  SortTPAAssembly (&align);

  uop = GetTPAAssembly(bsp);
  if (uop != NULL || (bsp->hist != NULL && bsp->hist->assembly != NULL)) {
    if (Message (MSG_OKC, "Sequence already has assembly - remove and replace with the new?") != ANS_OK) {
      while (align != NULL) {
        next_align = align->next;
        align->next = NULL;
        align = SeqAlignFree (align);
        align = next_align;
      }
      return;
    }
  }    

  if (bsp->hist != NULL) {
	  shp = bsp->hist;
	  if (shp->assembly != NULL)
	    SeqAlignSetFree(shp->assembly);
	  shp->assembly = align;
	}
  else
	{
	  shp = SeqHistNew();
	  shp->assembly = align;
	  bsp->hist = shp;
	}

  ValidateTPAHistAlign (bsp, &err_list);
  if (err_list != NULL) {
    lip = OpenLog ("TPA Problems");
    PrintTPAHistErrors (lip, err_list);
    lip->data_in_log = TRUE;
    err_list = ValNodeFreeData (err_list);
    CloseLog (lip);
    FreeLog (lip);
  }

  if (uop == NULL) 
  {
    sdp = CreateNewDescriptor (SeqMgrGetSeqEntryForData (bsp), Seq_descr_user);
    uop = CreateTpaAssemblyUserObject ();
    sdp->data.ptrvalue = uop;
  }

  RemoveOldTPAAccessions (uop);

  for (align = shp->assembly; align != NULL; align = align->next)
  {
    sip = AlnMgr2GetNthSeqIdPtr (align, 2);
    bsp_tmp = BioseqLockById (sip);
    if (bsp_tmp == NULL) {
      SeqIdWrite (sip, id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
    } else {
      SeqIdWrite (SeqIdFindBest (bsp_tmp->id, SEQID_GENBANK), id_buf, PRINTID_REPORT, sizeof (id_buf) - 1);
      BioseqUnlock (bsp_tmp);
    }
    sip = SeqIdFree (sip);
    AlnMgr2GetNthSeqRangeInSA (align, 1, &primary_start, &primary_stop);
    if (AlnMgr2GetNthStrand (align, 2) == Seq_strand_minus) {
      AddAccessionToTpaAssemblyUserObject (uop, id_buf, primary_stop, primary_start);
    } else {
      AddAccessionToTpaAssemblyUserObject (uop, id_buf, primary_start, primary_stop);
    }
  }


  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
  Update ();
}


NLM_EXTERN void AdjustForConsensusSplice (IteM i)
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
  if (sep == NULL) {
    return;
  }

  AdjustSeqEntryForConsensusSplice (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
  Update ();
}



static ValNodePtr BioseqListFromString(CharPtr str, SeqEntryPtr sep)
{
  ValNodePtr list = NULL;
  SeqIdPtr   sip;
  BioseqPtr  bsp;
  CharPtr    cp;
  Int4       len;
  Char       ch;
  ValNodePtr errors = NULL;

  if (str == NULL) {
    return NULL;
  }
  len = StringSpn (str, " \t\r\n,;");
  cp = str + len;

  while (*cp != 0) {
    len = StringCSpn (cp, " \t\r\n,;");
    ch = *(cp + len);
    *(cp + len) = 0;
    sip = CreateSeqIdFromText (cp, sep);
    if (sip == NULL) {
      ValNodeAddPointer (&errors, 0, StringSave (cp));
    } else {
      bsp = BioseqFind (sip);
      if (bsp == NULL) {
        ValNodeAddPointer (&errors, 0, StringSave (cp));
      } else {
        ValNodeAddPointer (&list, OBJ_BIOSEQ, bsp);
      }
      sip = SeqIdFree (sip);
    }
    *(cp + len) = ch;
    cp += len;
    len = StringSpn (cp, " \t\r\n,;");
    cp = cp + len;
  }
  
  if (errors != NULL) {
    cp = ValNodeMergeStrs (errors);
    Message (MSG_ERROR, "Unable to find Bioseqs for the following IDs: %s", cp);
    cp = MemFree (cp);
    list = ValNodeFree (list);
  } else if (list == NULL) {
    Message (MSG_ERROR, "No IDs listed!");
  }
  return list;
}


typedef struct withdrawsequencesfrm {
  FORM_MESSAGE_BLOCK

  TexT id_list;
  ButtoN leave_dlg_up;
} WithdrawSequencesFrmData, PNTR WithdrawSequencesFrmPtr;

static void DoWithdrawSequences (ButtoN b)
{
  WithdrawSequencesFrmPtr frm;
  SeqEntryPtr             sep;
  ValNodePtr              list = NULL;
  CharPtr                 list_str;
  PopSetRetroStatData     stat;

  frm = (WithdrawSequencesFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  /* get id list from text */
  list_str = SaveStringFromText (frm->id_list);
  list = BioseqListFromString (list_str, sep);
  list_str = MemFree (list_str);
  if (list == NULL) {
    return;
  }

  /* move listed sequences to wrapper */
  MoveSequencesFromSetToWrapper (list, frm->input_entityID);

  /* taxfix/cleanup */
  ForceCleanupEntityID (frm->input_entityID);

  /* redo set title */
  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));
  RemovePopsetTitles (sep);
  PopSetAutoDefRetro (sep, &stat);

  /* validate */
  ValSeqEntryForm (frm->form);

  if (!GetStatus (frm->leave_dlg_up)) {
    Remove (frm->form);
  }

}


NLM_EXTERN void WithdrawSequences (IteM i)
{
  BaseFormPtr  bfp;
  WithdrawSequencesFrmPtr frm;
  WindoW       w;
  GrouP        h, c;
  ButtoN       b;
  PrompT       ppt;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;


  frm = MemNew (sizeof (WithdrawSequencesFrmData));
  if (frm == NULL) return;
  frm->input_entityID = bfp->input_entityID;

  w = FixedWindow (-50, -33, -10, -10, "Withdraw Sequences", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupFormProc);
  frm->form = (ForM) w;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt = StaticPrompt (h, "List of Accessions", 0, popupMenuHeight, programFont, 'l');

  frm->id_list = ScrollText (h, 25, 10, programFont, TRUE, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoWithdrawSequences);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  frm->leave_dlg_up = CheckBox (c, "Leave dialog up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) frm->id_list,
                (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


extern void ParseTaxnameToSourceQual (IteM i)
{
  BaseFormPtr  bfp;
  ForM         f;
  ParseSrcOrgPtr org;
  ParseActionPtr parse;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  f = SingleParseAction (bfp->input_entityID);
  if (f != NULL) {
    /* create parse action with preselected fields */
    parse = ParseActionNew ();
    /* choose taxname after bionomial as source */
    org = ParseSrcOrgNew ();
    org->field = ValNodeNew (NULL);
    org->field->choice = ParseSrcOrgChoice_taxname_after_binomial;
    parse->src = ValNodeNew (NULL);
    parse->src->choice = ParseSrc_org;
    parse->src->data.ptrvalue = org;
    /* choose source qualifiers as dest */
    parse->dest = ValNodeNew (NULL);
    parse->dest->choice = ParseDest_org;
    PointerToForm (f, parse);
    parse = ParseActionFree (parse);
    Show (f);
  }
}


extern void ParseDefLineToSourceQual (IteM i)
{
  BaseFormPtr  bfp;
  ForM         f;
  ParseActionPtr parse;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  f = SingleParseAction (bfp->input_entityID);
  if (f != NULL) {
    /* create parse action with preselected fields */
    parse = ParseActionNew ();
    /* choose defline as source */
    parse->src = ValNodeNew (NULL);
    parse->src->choice = ParseSrc_defline;
    /* choose source qualifiers as dest */
    parse->dest = ValNodeNew (NULL);
    parse->dest->choice = ParseDest_org;
    PointerToForm (f, parse);
    parse = ParseActionFree (parse);
    Show (f);
  }
}

extern void ParseFlatfileToSourceQual (IteM i)
{
  BaseFormPtr  bfp;
  ForM         f;
  ParseActionPtr parse;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  f = SingleParseAction (bfp->input_entityID);
  if (f != NULL) {
    /* create parse action with preselected fields */
    parse = ParseActionNew ();
    /* choose defline as source */
    parse->src = ValNodeNew (NULL);
    parse->src->choice = ParseSrc_flatfile;
    /* choose source qualifiers as dest */
    parse->dest = ValNodeNew (NULL);
    parse->dest->choice = ParseDest_org;
    PointerToForm (f, parse);
    parse = ParseActionFree (parse);
    Show (f);
  }
}


extern void ParseLocalIDToSourceQual (ButtoN b)
{
  BaseFormPtr bfp;
  ForM         f;
  ParseActionPtr parse;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;

  f = SingleParseAction (bfp->input_entityID);
  if (f != NULL) {
    /* create parse action with preselected fields */
    parse = ParseActionNew ();
    /* choose localID as source */
    parse->src = ValNodeNew (NULL);
    parse->src->choice = ParseSrc_local_id;
    /* choose source qualifiers as dest */
    parse->dest = ValNodeNew (NULL);
    parse->dest->choice = ParseDest_org;
    PointerToForm (f, parse);
    parse = ParseActionFree (parse);
    Show (f);
  }
}


NLM_EXTERN DialoG SegregateIdDialog (GrouP g)
{
  DialoG dlg;
  StringConstraint   scd;

  dlg = StringConstraintDialog (g, "Mark sequences where sequence ID", TRUE, NULL, NULL);
  MemSet (&scd, 0, sizeof (StringConstraint));
  scd.match_location = eStringConstraintInList;
  PointerToDialog (dlg, &scd);
  return dlg;
}


static void ChooseCategoriesByStringConstraint (ValNodePtr value_list, StringConstraintPtr scp, Boolean do_choose)
{
  ClickableItemPtr cip;

  while (value_list != NULL) {
    cip = (ClickableItemPtr) value_list->data.ptrvalue;
      if (cip != NULL) {
      if (cip->item_list != NULL 
          && cip->item_list->choice == OBJ_BIOSEQ 
          && cip->item_list->data.ptrvalue != NULL
          && DoesSeqIDListMeetStringConstraint (((BioseqPtr) cip->item_list->data.ptrvalue)->id, scp)) {
        cip->chosen = do_choose;
      } else {
          ChooseCategoriesByStringConstraint (cip->subcategories, scp, do_choose);
      }
      }
      value_list = value_list->next;
  }
}


NLM_EXTERN void ChooseCategoriesForSegregateIdDialog(DialoG d, ValNodePtr list)
{
  StringConstraintPtr scp;

  scp = DialogToPointer (d);
  ChooseCategoriesByStringConstraint (list, scp, TRUE);
  scp = StringConstraintFree (scp);

}


NLM_EXTERN DialoG SegregateConstraintDialog (GrouP g)
{
  DialoG dlg;

  dlg = ComplexConstraintDialog (g, NULL, NULL);
  ChangeComplexConstraintFieldType (dlg, FieldType_molinfo_field, NULL, Macro_feature_type_any);
  return dlg;
}


static void ChooseCategoriesByConstraint (ValNodePtr value_list, ValNodePtr constraint, Boolean do_choose)
{
  ClickableItemPtr cip;

  while (value_list != NULL) {
    cip = (ClickableItemPtr) value_list->data.ptrvalue;
    if (cip != NULL) {
      if (cip->item_list != NULL 
          && DoesObjectMatchConstraintChoiceSet (cip->item_list->choice, cip->item_list->data.ptrvalue, constraint)) {
        cip->chosen = do_choose;
      } else {
          ChooseCategoriesByConstraint (cip->subcategories, constraint, do_choose);
      }
    }
    value_list = value_list->next;
  }
}


NLM_EXTERN void ChooseCategoriesForSegregateConstraintDialog(DialoG d, ValNodePtr list)
{
  ValNodePtr constraint;

  constraint = DialogToPointer (d);
  ChooseCategoriesByConstraint (list, constraint, TRUE);
  constraint = ConstraintChoiceSetFree (constraint);
}


NLM_EXTERN void UnchooseCategoriesForSegregateConstraintDialog(DialoG d, ValNodePtr list)
{
  ValNodePtr constraint;

  constraint = DialogToPointer (d);
  ChooseCategoriesByConstraint (list, constraint, FALSE);
  constraint = ConstraintChoiceSetFree (constraint);
}


static Boolean DomRNAAndProteinMatch (SeqFeatPtr mrna, ProtRefPtr prp)
{
  CharPtr mrna_product;
  Boolean rval = FALSE;

  if (mrna == NULL || prp == NULL) {
    return FALSE;
  }
  mrna_product = GetRNARefProductString (mrna->data.value.ptrvalue, NULL);
  if (StringHasNoText (mrna_product)) {
    if (prp->name == NULL || StringHasNoText (prp->name->data.ptrvalue)) {
      rval = TRUE;
    }
  } else {
    if (prp->name != NULL && StringCmp (mrna_product, prp->name->data.ptrvalue) == 0) {
      rval = TRUE;
    }
  }
  mrna_product = MemFree (mrna_product);
  return rval;
}


NLM_EXTERN void BulkCDSSetProtein (Pointer target, Pointer data)
{
  SeqFeatXrefPtr xref;
  ProtRefPtr     prp = NULL;
  BioseqPtr      prot_bsp;
  SeqFeatPtr     prot_feat;
  SeqMgrFeatContext fcontext;
  ValNodePtr        prot_name_list = NULL, vnp;
  SeqFeatPtr        sfp, mrna;

  if ((sfp = (SeqFeatPtr) target) == NULL) return;

  prot_name_list = (ValNodePtr) data;

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
      if (prot_feat == NULL && prot_name_list != NULL) {
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

  if (prp == NULL && prot_name_list != NULL) {
    xref = SeqFeatXrefNew ();
    prp = ProtRefNew ();
    xref->data.choice = SEQFEAT_PROT;
    xref->data.value.ptrvalue = prp;
    xref->next = sfp->xref;
    sfp->xref = xref;
  }

  mrna = SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
  if (mrna != NULL) {
    if (!DomRNAAndProteinMatch(mrna, prp)) {
      mrna = NULL;
    }
  }

  if (prp != NULL) {
    prp->name = ValNodeFreeData (prp->name);
    for (vnp = prot_name_list; vnp != NULL; vnp = vnp->next) {
      ValNodeAddPointer (&prp->name, vnp->choice, StringSave (vnp->data.ptrvalue));
    }
    if (mrna != NULL && prot_name_list != NULL) {
      SetRNARefProductString (mrna->data.value.ptrvalue, NULL, prot_name_list->data.ptrvalue, ExistingTextOption_replace_old);
    }
  }

}


typedef struct stringobjectpair {
  CharPtr str;
  Uint1   choice;
  Pointer data;
} StringObjectPairData, PNTR StringObjectPairPtr;


static StringObjectPairPtr StringObjectPairNew (CharPtr str, Uint1 choice, Pointer data)
{
  StringObjectPairPtr sop = (StringObjectPairPtr) MemNew (sizeof (StringObjectPairData));
  sop->str = str;
  sop->choice = choice;
  sop->data = data;
  return sop;
}


static StringObjectPairPtr StringObjectPairFree (StringObjectPairPtr sop)
{
  if (sop != NULL) {
    sop->str = MemFree (sop->str);
    sop = MemFree (sop);
  }
  return sop;
}


static ValNodePtr StringObjectPairListFree (ValNodePtr list)
{
  ValNodePtr vnp_next;

  while (list != NULL) {
    vnp_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = StringObjectPairFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = vnp_next;
  }
  return list;
}


static int LIBCALLBACK SortVnpByStringObjectPair (VoidPtr ptr1, VoidPtr ptr2)

{
  StringObjectPairPtr     sop1;
  StringObjectPairPtr     sop2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      sop1 = (StringObjectPairPtr) vnp1->data.ptrvalue;
      sop2 = (StringObjectPairPtr) vnp2->data.ptrvalue;
      if (sop1 != NULL && sop2 != NULL) {
        return StringCmp (sop1->str, sop2->str);
      }
    }
  }
  return 0;
}


typedef struct uniqstructuredcomment {
  StructuredCommentFieldPtr field;
  ValNodeBlock pair_list;
} UniqStructuredCommentData, PNTR UniqStructuredCommentPtr;


static void GetUniqueValuesForStructuredCommentFieldCallback (SeqDescrPtr sdp, Pointer data)
{
  UniqStructuredCommentPtr usc;
  UserObjectPtr uop;
  CharPtr str;

  if (sdp == NULL || sdp->choice != Seq_descr_user 
      || (uop = (UserObjectPtr) sdp->data.ptrvalue) == NULL
      || uop->type == NULL
      || StringCmp (uop->type->str, "StructuredComment") != 0
      || (usc = (UniqStructuredCommentPtr)data) == NULL) {
    return;
  }

  str = GetStructuredCommentFieldFromUserObject (uop, usc->field, NULL);
  if (str != NULL) {
    ValNodeAddPointerToEnd (&(usc->pair_list), 0, StringObjectPairNew(str, OBJ_SEQDESC, sdp));
  }
}


NLM_EXTERN ValNodePtr GetUniqueValuesForStructuredCommentField (SeqEntryPtr sep, CharPtr field_name)
{
  UniqStructuredCommentData usd;
  ValNodePtr list, vnp;
  ClickableItemPtr last_cip = NULL;
  StringObjectPairPtr sop;
  ValNodeBlock cip_list;

  usd.field = ValNodeNew (NULL);
  usd.field->choice = StructuredCommentField_named;
  usd.field->data.ptrvalue = field_name;
  InitValNodeBlock(&(usd.pair_list), NULL);

  VisitDescriptorsInSep (sep, &usd, GetUniqueValuesForStructuredCommentFieldCallback);
  usd.field = ValNodeFree (usd.field);
  list = usd.pair_list.head;
  list = ValNodeSort (list, SortVnpByStringObjectPair);

  InitValNodeBlock(&cip_list, NULL);

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    sop = (StringObjectPairPtr) vnp->data.ptrvalue;
    if (sop != NULL) {
      if (last_cip == NULL || StringCmp (last_cip->description, sop->str) != 0) {
        /* make new category */
        last_cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        MemSet (last_cip, 0, sizeof (ClickableItemData));
        last_cip->description = sop->str;
        sop->str = NULL;
        ValNodeAddPointerToEnd (&(cip_list), 0, last_cip);
      }
      /* add to category */
      ValNodeAddPointer (&(last_cip->item_list), sop->choice, sop->data);        
    }
  }
  list = StringObjectPairListFree(list);

  return cip_list.head;
}


static void GetStructuredCommentFieldNamesCallback (SeqDescrPtr sdp, Pointer data)
{
  UserObjectPtr uop;
  UserFieldPtr  ufp;

  if (sdp != NULL && data != NULL && sdp->choice == Seq_descr_user
      && (uop = sdp->data.ptrvalue) != NULL
      && IsUserObjectStructuredComment (uop)) {

    ufp = uop->data;
    while (ufp != NULL) {
      if (ufp->label != NULL && ufp->label->str != NULL
          && StringCmp (ufp->label->str, "StructuredCommentPrefix") != 0
          && StringCmp (ufp->label->str, "StructuredCommentSuffix") != 0) {
        ValNodeAddPointer ((ValNodePtr PNTR) data, 0, StringSave (ufp->label->str));
      }
      ufp = ufp->next;
    }
  }
}


NLM_EXTERN ValNodePtr GetStructuredCommentFieldNames (SeqEntryPtr sep)
{
  ValNodePtr list = NULL;

  VisitDescriptorsInSep (sep, &list, GetStructuredCommentFieldNamesCallback);
  list = ValNodeSort (list, SortVnpByString);
  ValNodeUnique (&list, SortVnpByString, ValNodeFreeData);
  return list;
}


typedef struct badqual {
  Boolean wrong;
  Boolean illegal;
  Boolean move_to_note;
} BadQualData, PNTR BadQualPtr;

static void RemoveBadQualsCommonCallback (SeqFeatPtr sfp, Pointer data)
{
  BadQualPtr bq;
  GBQualPtr  gbqual, gbqual_prev = NULL, gbqual_next;
  Int2       index;
  Int2       val;
  CharPtr    key;
  Boolean    do_remove;
  CharPtr    note_txt;
  CharPtr    note_fmt = "[%s=%s]";
  Int4       qvalid;

  if (sfp == NULL || (bq = (BadQualPtr) data) == NULL) {
    return;
  }
  key = GetGBFeatKeyForFeature (sfp);
  index = GBFeatKeyNameValid (&key, FALSE);
  key = MemFree (key);
  if (index < 0) {
    return;
  }

  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual_next) {
    gbqual_next = gbqual->next;
    do_remove = FALSE;
    qvalid = IsQualValidForFeature(gbqual, sfp);
    if (qvalid == 1) {
      val = GBQualNameValid (gbqual->qual);
      if (bq->illegal && !ShouldBeAGBQual (sfp->idx.subtype, val, FALSE)) {
        do_remove = TRUE;
      }
    } else {
      if (bq->illegal && (val = GBQualNameValid (gbqual->qual)) == -1) {
        do_remove = TRUE;
      } else if (bq->wrong && (qvalid == 0 || qvalid > 1)) {
        do_remove = TRUE;
      }
    }
    if (do_remove) {
      if (bq->move_to_note) {
        note_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (note_fmt) + StringLen (gbqual->qual) + StringLen (gbqual->val)));
        sprintf (note_txt, note_fmt, gbqual->qual == NULL ? "" : gbqual->qual, gbqual->val == NULL ? "" : gbqual->val);
        SetStringValue (&(sfp->comment), note_txt, ExistingTextOption_append_semi);
        note_txt = MemFree (note_txt);
      }
      if (gbqual_prev == NULL) {
        sfp->qual = gbqual_next;
      } else {
        gbqual_prev->next = gbqual_next;
      }
      gbqual->next = NULL;
      gbqual = GBQualFree (gbqual);
    } else {
      gbqual_prev = gbqual;
    }
  }
}


NLM_EXTERN void RemoveBadQualsCommon (Uint2 entityID, Boolean wrong, Boolean illegal, Boolean move_to_note)
{
  SeqEntryPtr  sep;
  BadQualData  bq;


  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) {
    return;
  }

  bq.wrong = wrong;
  bq.illegal = illegal;
  bq.move_to_note = move_to_note;

  VisitFeaturesInSep (sep, &bq, RemoveBadQualsCommonCallback);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  Update();
}


static void RemoveBadQualsCommonItem (IteM i, Boolean wrong, Boolean illegal, Boolean move_to_note)
{
  BaseFormPtr  bfp;
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;

  RemoveBadQualsCommon (bfp->input_entityID, wrong, illegal, move_to_note);
}


NLM_EXTERN void RemoveIllegalQuals (IteM i)
{
  RemoveBadQualsCommonItem (i, FALSE, TRUE, FALSE);
}


NLM_EXTERN void RemoveWrongQuals (IteM i)
{
  RemoveBadQualsCommonItem (i, TRUE, FALSE, FALSE);
}


NLM_EXTERN void MoveIllegalQuals (IteM i)
{
  RemoveBadQualsCommonItem (i, FALSE, TRUE, TRUE);
}


NLM_EXTERN void MoveWrongQuals (IteM i)
{
  RemoveBadQualsCommonItem (i, TRUE, FALSE, TRUE);
}


static SeqIdPtr GetProteinSeqIdFromBioseq (BioseqPtr bsp)
{
  SeqIdPtr sip;
  DbtagPtr dbtag;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL
        && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL
        && !IsSkippableDbtag(dbtag)) {
      return sip;
    }
  }
  return NULL;
}


static void FindProteinIDDBCallback (BioseqPtr bsp, Pointer data)
{
  CharPtr PNTR pDb;
  SeqIdPtr sip;
  DbtagPtr dbtag;

  if (bsp == NULL || !ISA_aa(bsp->mol) || (pDb = (CharPtr PNTR)data) == NULL || *pDb != NULL) {
    return;
  }

  sip = GetProteinSeqIdFromBioseq (bsp);
  if (sip != NULL 
      && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL
      && !IsSkippableDbtag(dbtag)) {
    *pDb = StringSave (dbtag->db);
  }
}



static void CreateOneProteinId (BioseqPtr bsp, CharPtr db, CharPtr locus_tag, LogInfoPtr lip)
{
  SeqIdPtr sip;
  DbtagPtr dbtag;
  SeqEntryPtr sep;
  SeqIdPtr sip_local, sip_last = NULL;
  Char     id1[255];
  Char     id2[255];

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL
        && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL
        && !IsSkippableDbtag(dbtag)) {
      return;
    }
    if (sip->choice == SEQID_LOCAL) {
      sip_local = sip;
    }
    sip_last = sip;
  }

  sip = ValNodeNew(NULL);
  sip->choice = SEQID_GENERAL;
  dbtag = DbtagNew ();
  dbtag->db = StringSave (db);
  dbtag->tag = ObjectIdNew();
  dbtag->tag->str = StringSave (locus_tag);
  sip->data.ptrvalue = dbtag;
  SeqIdWrite (sip, id2, PRINTID_REPORT, sizeof (id2) - 1);
  if (sip_local == NULL) {
    SeqIdWrite (SeqIdFindWorst (bsp->id), id1, PRINTID_REPORT, sizeof (id1) - 1);
    if (sip_last == NULL) {
      bsp->id = sip;
    } else {
      sip_last->next = sip;
    }
    /* replace in bioseq index */
    SeqMgrReplaceInBioseqIndex (bsp);
    if (lip != NULL) {
      fprintf (lip->fp, "Added protein ID %s to Bioseq %s\n", id2, id1);
      lip->data_in_log = TRUE;
    }
  } else {
    SeqIdWrite (sip_local, id1, PRINTID_REPORT, sizeof (id1) - 1);
    sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
    ReplaceSeqIdWithSeqId (sip_local, sip, sep);
    sip = SeqIdFree (sip);
    if (lip != NULL) {
      fprintf (lip->fp, "Replaced protein ID %s with %s\n", id1, id2);
      lip->data_in_log = TRUE;
    }
  }
}


static void CreateProteinIDsCallback (BioseqPtr bsp, Pointer data)
{
  CharPtr db;
  SeqIdPtr sip;
  DbtagPtr dbtag;
  SeqFeatPtr cds, gene = NULL;
  GeneRefPtr grp = NULL;
  SeqMgrFeatContext context;

  if (bsp == NULL || !ISA_aa(bsp->mol) || (db = (CharPtr)data) == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL
        && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL
        && !IsSkippableDbtag(dbtag)) {
      return;
    }
  }

  cds = SeqMgrGetCDSgivenProduct (bsp, &context);
  if (cds == NULL) {
    return;
  }
  GetGeneInfoForFeature (cds, &grp, &gene);
  if (grp != NULL && !StringHasNoText (grp->locus_tag)) {
    CreateOneProteinId (bsp, db, grp->locus_tag, NULL);
  }

}


typedef struct justtextform {
  FORM_MESSAGE_BLOCK
  TexT       db_txt;
} JustTextFormData, PNTR JustTextFormPtr;


static void DoCreateProteinIDs (ButtoN b)
{
  JustTextFormPtr frm;
  CharPtr db;
  SeqEntryPtr sep;

  frm = (JustTextFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  } else if (TextHasNoText (frm->db_txt)) {
    Message (MSG_ERROR, "Must provide database tag!");
    return;
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  WatchCursor();
  Update();
  db = SaveStringFromText (frm->db_txt);
  VisitBioseqsInSep (sep, db, CreateProteinIDsCallback);
  db = MemFree (db);
  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  Remove (frm->form);
  ArrowCursor();
  Update();
}


NLM_EXTERN void CreateProteinIDs (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  JustTextFormPtr frm;
  WindoW       w;
  GrouP        h, g, c;
  ButtoN       b;
  CharPtr      db = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, &db, FindProteinIDDBCallback);

  frm = (JustTextFormPtr) MemNew (sizeof (JustTextFormData));

  w = FixedWindow (-50, -33, -10, -10, "Create Protein IDs", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  frm->form = (ForM) w;
  frm->input_entityID = bfp->input_entityID;

  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Database", 0, dialogTextHeight, programFont, 'c');
  frm->db_txt = DialogText (g, db, 14, NULL);
  db = MemFree (db);

  c = HiddenGroup (h, 2, 0, NULL);
  b = DefaultButton (c, "Accept", DoCreateProteinIDs);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  Show (w);

}


typedef struct existinglocustag {
  CharPtr prefix;
  Int4    num_digits;
  Int4    max;
  CharPtr id_db;
  LogInfoPtr lip;
} ExistingLocusTagData, PNTR ExistingLocusTagPtr;


static Int4 GetLocusTagPrefixLen (CharPtr str)

{
  Int4    len;

  if (StringHasNoText (str)) return 0;

  len = StringLen (str) - 1;

  while (len > 0 && isdigit (*(str + len))) {
    len--;
  }
  return len;
}


static void PopulateExistingLocusTag (ExistingLocusTagPtr pLocusTag, CharPtr locus_tag)
{
  Int4 prefix_len;
  Int4 val;
  Int4 num_digits;

  if (pLocusTag == NULL || StringHasNoText (locus_tag)
      || (prefix_len = GetLocusTagPrefixLen(locus_tag)) < 1) {
    return;
  }

  if (pLocusTag->prefix == NULL) {
    pLocusTag->prefix = (CharPtr) MemNew (sizeof (Char) * (prefix_len + 1));
    StringNCpy (pLocusTag->prefix, locus_tag, prefix_len);
    pLocusTag->prefix[prefix_len] = 0;
  } 
  if (StringNCmp (locus_tag, pLocusTag->prefix, prefix_len) == 0) {
    num_digits = StringLen (locus_tag + prefix_len + 1);
    if (pLocusTag->num_digits < num_digits) {
      pLocusTag->num_digits = num_digits;
    }
    val = atol (locus_tag + prefix_len + 1);
    if (val > pLocusTag->max) {
      pLocusTag->max = val;
    }
  }
}


static void FindGenePrefix (SeqFeatPtr sfp, Pointer data)
{
  ExistingLocusTagPtr pLocusTag;
  GeneRefPtr grp;

  if (sfp == NULL 
      || sfp->data.choice != SEQFEAT_GENE
      || (pLocusTag = (ExistingLocusTagPtr)data) == NULL
      || (grp = (GeneRefPtr) sfp->data.value.ptrvalue) == NULL) {
    return;
  }

  PopulateExistingLocusTag (pLocusTag, grp->locus_tag);
}


static void FindGeneLocusTagFromProteinID (BioseqPtr bsp, Pointer data)
{
  ExistingLocusTagPtr pLocusTag;
  SeqIdPtr sip;
  DbtagPtr dbtag;

  if (bsp == NULL || !ISA_aa (bsp->mol)
      || (pLocusTag = (ExistingLocusTagPtr)data) == NULL) {
    return;
  }

  sip = GetProteinSeqIdFromBioseq (bsp);
  if (sip != NULL 
      && (dbtag = (DbtagPtr) sip->data.ptrvalue) != NULL
      && !IsSkippableDbtag(dbtag)) {
    PopulateExistingLocusTag (pLocusTag, dbtag->tag->str);
    if (pLocusTag->id_db == NULL) {
      pLocusTag->id_db = StringSave (dbtag->db);
    }
  }
}


static Boolean FeatureSubtypeWantsLocusTagGene (Int4 featdef)
{
  if (featdef == FEATDEF_CDS
      || featdef == FEATDEF_tRNA
      || featdef == FEATDEF_rRNA
      || featdef == FEATDEF_tmRNA
      || featdef == FEATDEF_ncRNA) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void CreateLocusTagGenesCallback (SeqFeatPtr sfp, Pointer data)
{
  ExistingLocusTagPtr locus_data;
  SeqFeatPtr gene = NULL, mrna;
  GeneRefPtr grp = NULL;
  Char buf[15];
  BioseqPtr bsp, pbsp = NULL;
  SeqLocPtr gene_loc;
  DbtagPtr  dbtag = NULL;
  SeqIdPtr  sip;
  Char      fmt_buf[30];
  ValNode   vn;
  CharPtr   txt;
  Boolean   created_gene = FALSE;

  /* create locus tag genes for features that don't have them */
  if (sfp == NULL || !FeatureSubtypeWantsLocusTagGene(sfp->idx.subtype) || (locus_data = (ExistingLocusTagPtr) data) == NULL) {
    return;
  }

  if (sfp->idx.subtype == FEATDEF_CDS) {
    pbsp = BioseqFindFromSeqLoc (sfp->product);
    sip = GetProteinSeqIdFromBioseq (pbsp);
    if (sip != NULL) { 
      dbtag = (DbtagPtr) sip->data.ptrvalue;
    }
  }
      
  GetGeneInfoForFeature (sfp, &grp, &gene);
  if (grp == NULL) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      mrna = GetmRNAforCDS(sfp);
      if (mrna == NULL) {
        gene_loc = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
      } else {
        gene_loc = SeqLocMerge (bsp, mrna->location, NULL, TRUE, FALSE, FALSE);
      }
      if (gene_loc != NULL) {
        gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, gene_loc);
        gene_loc = SeqLocFree (gene_loc);
        grp = GeneRefNew ();
        gene->data.value.ptrvalue = grp;
        created_gene = TRUE;
      }
    }
  }
  if (grp != NULL && StringHasNoText (grp->locus_tag)) {
    grp->locus_tag = MemFree (grp->locus_tag);

    if (dbtag != NULL && dbtag->tag != NULL && !StringHasNoText (dbtag->tag->str))
    {
      grp->locus_tag = StringSave (dbtag->tag->str);
    } else {
      sprintf (fmt_buf, "%%0%dd", locus_data->num_digits);
      sprintf (buf, fmt_buf, locus_data->max + 5);

      grp->locus_tag = (CharPtr) MemNew (sizeof (Char) * (StringLen (locus_data->prefix) + StringLen (buf) + 2));
      sprintf (grp->locus_tag, "%s_%s", locus_data->prefix, buf);
      locus_data->max += 5;
    }
    if (locus_data->lip != NULL) {
      MemSet (&vn, 0, sizeof (ValNode));
      vn.choice = OBJ_SEQFEAT;
      vn.data.ptrvalue = sfp;
      txt = GetDiscrepancyItemText (&vn);

      if (locus_data->lip->fp != NULL) {
        fprintf (locus_data->lip->fp, "Added locus tag %s%s for %s\n", 
                                      grp->locus_tag, 
                                      created_gene ? " (and created gene)":"", 
                                      txt);
      }
      locus_data->lip->data_in_log = TRUE;
      txt = MemFree (txt);
    }

    if (pbsp != NULL && !StringHasNoText (locus_data->id_db)) {
      CreateOneProteinId (pbsp, locus_data->id_db, grp->locus_tag, locus_data->lip);
    }
  }
}


typedef struct locustaggeneform {
  FORM_MESSAGE_BLOCK
  TexT       prefix_txt;
  ButtoN     also_protein_ids;
  TexT       db_txt;

  ExistingLocusTagData locus_data;
} LocusTagGeneFormData, PNTR LocusTagGeneFormPtr;


static void CleanupLocusTagGeneForm (
  GraphiC g,
  VoidPtr data
)

{
  LocusTagGeneFormPtr frm;

  frm = (LocusTagGeneFormPtr) data;
  if (frm != NULL) {
    frm->locus_data.prefix = MemFree (frm->locus_data.prefix);
    frm->locus_data.id_db = MemFree (frm->locus_data.id_db);
  }
  StdCleanupFormProc (g, data);
}


static void DoCreateLocusTagGenes (ButtoN b)
{
  LocusTagGeneFormPtr frm;
  SeqEntryPtr sep;

  frm = (LocusTagGeneFormPtr) GetObjectExtra (b);
  if (frm == NULL || TextHasNoText (frm->prefix_txt)) {
    return;
  }
  frm->locus_data.prefix = MemFree (frm->locus_data.prefix);
  frm->locus_data.prefix = SaveStringFromText (frm->prefix_txt);
  frm->locus_data.id_db = MemFree (frm->locus_data.id_db);
  if (GetStatus (frm->also_protein_ids)) {
    frm->locus_data.id_db = SaveStringFromText (frm->db_txt);
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);
  WatchCursor();
  Update();
  frm->locus_data.lip = OpenLog ("Locus Tags Created");
  VisitFeaturesInSep (sep, &(frm->locus_data), CreateLocusTagGenesCallback);
  CloseLog (frm->locus_data.lip);
  frm->locus_data.lip = FreeLog (frm->locus_data.lip);
  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  Remove (frm->form);
  ArrowCursor();
  Update();
}


static void AlsoProteinIdsBtn (ButtoN b)
{
  LocusTagGeneFormPtr frm;

  frm = (LocusTagGeneFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  if (GetStatus (frm->also_protein_ids)) {
    Enable (frm->db_txt);
  } else {
    Disable (frm->db_txt);
  }
}


NLM_EXTERN void AddLocusTagGenes (IteM i)
{
  BaseFormPtr         bfp;
  SeqEntryPtr         sep;
  LocusTagGeneFormPtr frm;
  WindoW              w;
  GrouP               h, g, g2, c;
  ButtoN              b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);

  frm = (LocusTagGeneFormPtr) MemNew (sizeof (LocusTagGeneFormData));
  frm->locus_data.prefix = NULL;
  frm->locus_data.max = 0;
  VisitFeaturesInSep (sep, &(frm->locus_data), FindGenePrefix);
  VisitBioseqsInSep (sep, &(frm->locus_data), FindGeneLocusTagFromProteinID);
  if (frm->locus_data.num_digits < 1) {
    frm->locus_data.num_digits = 4;
  }

  w = FixedWindow (-50, -33, -10, -10, "Create locus_tag genes", StdCloseWindowProc);
  SetObjectExtra (w, frm, CleanupLocusTagGeneForm);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  frm->form = (ForM) w;
  frm->input_entityID = bfp->input_entityID;

  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Prefix", 0, dialogTextHeight, programFont, 'c');
  frm->prefix_txt = DialogText (g, frm->locus_data.prefix, 14, NULL);

  g2 = HiddenGroup (h, 2, 0, NULL);
  frm->also_protein_ids = CheckBox (g2, "Also create protein IDs with database:", AlsoProteinIdsBtn);
  SetObjectExtra (frm->also_protein_ids, frm, NULL);
  frm->db_txt = DialogText (g2, frm->locus_data.id_db, 14, NULL);
  Disable (frm->db_txt);

  c = HiddenGroup (h, 2, 0, NULL);
  b = DefaultButton (c, "Accept", DoCreateLocusTagGenes);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  Show (w);

}
