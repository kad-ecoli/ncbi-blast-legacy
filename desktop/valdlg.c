/*   valdlg.c
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
* File Name:  valdlg.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   4/8/2009
*
* $Revision: 1.4 $
*
* File Description: Dialogs for editing structured comment rules
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <vibforms.h>
#include <document.h>
#include <objseq.h>
#include <sqnutils.h>

#define NLM_GENERATED_CODE_PROTO
#include <objvalid.h>
#include <valapi.h>
#include <valdlg.h>

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

static Uint2 commentruleedit_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_POPUP
};

static Uint2 commentruleedit_widths [] = {
  20, 20, 20
};

ENUM_ALIST(true_false)
  { "FALSE",  0 },
  { "TRUE",   1 },
END_ENUM_ALIST

static EnumFieldAssocPtr commentruleedit_alists [] = { NULL, NULL, true_false };

#define COMMENTRULEEDIT_FIELDNAME_COLUMN 0
#define COMMENTRULEEDIT_MATCH_COLUMN  1
#define COMMENTRULEEDIT_REQUIRED_COLUMN  2


static CharPtr TagStringFromFieldRule (FieldRulePtr rule)
{
  Int4 len;
  CharPtr fmt = "%s\t%s\t%d\n";
  CharPtr str;

  if (rule == NULL) {
    return NULL;
  }

  len = StringLen (fmt) + StringLen (rule->field_name) + StringLen (rule->match_expression);
  str = (CharPtr) MemNew (sizeof (Char) * len);
  sprintf (str, fmt, rule->field_name == NULL ? "" : rule->field_name,
           rule->match_expression == NULL ? "" : rule->match_expression,
           rule->required ? 1 : 0);
  return str;
}


static FieldRulePtr FieldRuleFromTagString (CharPtr str)
{
  FieldRulePtr rule;
  CharPtr tmp;

  if (StringHasNoText (str)) {
    return NULL;
  }

  rule = FieldRuleNew ();
  rule->field_name = ExtractTagListColumn (str, COMMENTRULEEDIT_FIELDNAME_COLUMN);
  rule->match_expression = ExtractTagListColumn (str, COMMENTRULEEDIT_MATCH_COLUMN);
  tmp = ExtractTagListColumn (str, COMMENTRULEEDIT_REQUIRED_COLUMN);
  if (tmp != NULL && atoi (tmp) > 0) {
    rule->required = TRUE;
  } else {
    rule->required = FALSE;
  }
  tmp = MemFree (tmp);

  return rule;
}


static void CommentFieldsToDialog (DialoG d, Pointer data)
{
  TagListPtr tlp;
  ValNodePtr fields, vnp;
  CharPtr    str;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) {
    return;
  }

  fields = (ValNodePtr) data;

  tlp->vnp = ValNodeFreeData (tlp->vnp);
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  for (vnp = fields; vnp != NULL; vnp = vnp->next) {
    str = TagStringFromFieldRule (vnp->data.ptrvalue);
    if (str != NULL) {
      ValNodeAddPointer (&(tlp->vnp), 0, str);
    }
  }

  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  tlp->max = MAX ((Int2) 0, (Int2) (ValNodeLen (tlp->vnp) - tlp->rows));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  if (tlp->max > 0) {
    SafeShow (tlp->bar);
  } else {
    SafeHide (tlp->bar);
  }
}


static Pointer CommentFieldsFromDialog (DialoG d)
{
  TagListPtr tlp;
  ValNodePtr fields = NULL, vnp;
  FieldRulePtr rule;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) {
    return NULL;
  }

  for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
    rule = FieldRuleFromTagString (vnp->data.ptrvalue);
    if (rule != NULL) {
      ValNodeAddPointer (&fields, 0, rule);
    }
  }
  return fields;
}


static DialoG CommentRuleFieldsDialog (GrouP h)
{
  DialoG dlg;

  dlg = CreateTagListDialogEx (h, 5, 3, 2,
                               commentruleedit_types, commentruleedit_widths,
                               commentruleedit_alists, TRUE, FALSE, 
                               CommentFieldsToDialog, CommentFieldsFromDialog);
  return dlg;
}


typedef struct commentruledlg
{
  DIALOG_MESSAGE_BLOCK
  TexT prefix;
  DialoG fields_dlg;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} CommentRuleDlgData, PNTR CommentRuleDlgPtr;


static void CommentRuleToDialog (DialoG d, Pointer data)
{
  CommentRuleDlgPtr dlg;
  CommentRulePtr    rule;

  dlg = (CommentRuleDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  rule = (CommentRulePtr) data;

  if (rule == NULL) {
    SetTitle (dlg->prefix, "");
    PointerToDialog (dlg->fields_dlg, NULL);
  } else {
    SetTitle (dlg->prefix, rule->prefix == NULL ? "" : rule->prefix);
    PointerToDialog (dlg->fields_dlg, rule->fields);
  }
}


static Pointer CommentRuleFromDialog (DialoG d)
{
  CommentRuleDlgPtr dlg;
  CommentRulePtr    rule;

  dlg = (CommentRuleDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  rule = CommentRuleNew ();
  rule->prefix = SaveStringFromText (dlg->prefix);
  rule->fields = DialogToPointer (dlg->fields_dlg);
  return rule;
}


static void ChangeCommentRuleDialogPrefix (TexT t)
{
  CommentRuleDlgPtr dlg;

  dlg = (CommentRuleDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


NLM_EXTERN DialoG CommentRuleDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  CommentRuleDlgPtr dlg;
  GrouP             p, g1, g2;
  
  dlg = (CommentRuleDlgPtr) MemNew (sizeof (CommentRuleDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = CommentRuleToDialog;
  dlg->fromdialog = CommentRuleFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "Prefix", 0, popupMenuHeight, programFont, 'r');
  dlg->prefix = DialogText (g1, "", 20, ChangeCommentRuleDialogPrefix);
  SetObjectExtra (dlg->prefix, dlg, NULL);

  g2 = HiddenGroup (p, 3, 0, NULL);
  StaticPrompt (g2, "Field Name", commentruleedit_widths[COMMENTRULEEDIT_FIELDNAME_COLUMN] * Nlm_stdCharWidth, popupMenuHeight, programFont, 'c');
  StaticPrompt (g2, "Match Expression", commentruleedit_widths[COMMENTRULEEDIT_MATCH_COLUMN] * Nlm_stdCharWidth, popupMenuHeight, programFont, 'c');
  StaticPrompt (g2, "Required", commentruleedit_widths[COMMENTRULEEDIT_REQUIRED_COLUMN] * Nlm_stdCharWidth, popupMenuHeight, programFont, 'c');

  dlg->fields_dlg = CommentRuleFieldsDialog (p);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, (HANDLE) dlg->fields_dlg, NULL);
  
  return (DialoG) p;
}


typedef struct commentsetdlg
{
  DIALOG_MESSAGE_BLOCK
  DoC rule_list_dlg;
  DialoG rule_editor;
  FonT   font;

  CommentRulePtr rule_list;
  Int4           curr_item;
} CommentSetDlgData, PNTR CommentSetDlgPtr;

static void ListRules (DoC doc, CommentRulePtr rule_list, FonT font)
{
  CommentRulePtr cr;
  RecT       r;


  Reset (doc);

  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  
  if (font == NULL) font = programFont;

  for (cr = rule_list; cr != NULL; cr = cr->next) {
    AppendText (doc, cr->prefix == NULL ? "" : cr->prefix, NULL, NULL, font);
  }
  AppendText (doc, "Click here to add new rule", NULL, NULL, font);

  UpdateDocument (doc, 0, 0);
}


static Int2 ReplaceNthCommentRule (CommentRulePtr PNTR list, CommentRulePtr new_cr, Int2 pos)
{
  Int2 match;
  CommentRulePtr replace_cr, last_cr = NULL;

  if (list == NULL) {
    return 0;
  }

  for (match = 1, replace_cr = *list;
       match < pos && replace_cr != NULL;
       match++, replace_cr = replace_cr->next) {
    last_cr = replace_cr;
  }

  if (replace_cr == NULL) {
    if (last_cr == NULL) {
      *list = new_cr;
    } else {
      last_cr->next = new_cr;
    }
  } else if (new_cr == NULL) {
    if (last_cr != NULL) {
      last_cr->next = replace_cr->next;
      replace_cr->next = NULL;
      replace_cr = CommentRuleFree (replace_cr);
    }
  } else {
    new_cr->next = replace_cr->next;
    replace_cr->next = NULL;
    replace_cr = CommentRuleFree (replace_cr);
    if (last_cr == NULL) {
      *list = new_cr;
    } else {
      last_cr->next = new_cr;
    }
  }
  return match;
}


static void ClickRulesDoc (DoC d, PoinT pt)
{
  Int2               item, row, col, match;
  RecT               rct;
  CommentSetDlgPtr   dlg;
  BaR                sb_vert = NULL;
  Int2               scroll_pos = 0;
  CommentRulePtr     cr, last_cr = NULL;
  
  dlg = (CommentSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->rule_list == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item == 0 && dlg->rule_list == NULL) {
    cr = CommentRuleNew ();
    dlg->rule_list = cr;
    dlg->curr_item = 1;
    PointerToDialog (dlg->rule_editor, cr);
    ListRules (dlg->rule_list_dlg, dlg->rule_list, dlg->font);
  } else if (item > 0 && row > 0) {
    if (item != dlg->curr_item) {
      cr = DialogToPointer (dlg->rule_editor);
      ReplaceNthCommentRule (&(dlg->rule_list), cr, dlg->curr_item);

      cr = dlg->rule_list;
      match = 1;
      while (match < item && cr != NULL) {
        match++;
        last_cr = cr;
        cr = cr->next;
      }
      if (cr == NULL) {
        cr = CommentRuleNew ();
        if (last_cr == NULL) {
          dlg->rule_list = cr;
        } else {
          last_cr->next = cr;
        }
      }

      PointerToDialog (dlg->rule_editor, cr);
      dlg->curr_item = item;
      ListRules (dlg->rule_list_dlg, dlg->rule_list, dlg->font);
    }
  }
}


static Boolean RuleHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  CommentSetDlgPtr dlg;
  
  dlg = (CommentSetDlgPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;
  
  if (item == dlg->curr_item) return TRUE;
  return FALSE;
}


static void SetupCommentSetDialogFont (CommentSetDlgPtr dlg)

{
  if (dlg == NULL) return;

#ifdef WIN_MAC
  dlg->font = ParseFont ("Times,12");
#endif
#ifdef WIN_MSWIN
  dlg->font = ParseFont ("Times New Roman,12");
#endif
#ifdef WIN_MOTIF
  dlg->font = ParseFont ("Times,12");
#endif
}


static void CommentSetToDialog (DialoG d, Pointer data)
{
  CommentSetDlgPtr dlg;

  dlg = (CommentSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  dlg->rule_list = CommentSetFree (dlg->rule_list);
  dlg->rule_list = (CommentRulePtr) data;
  ListRules (dlg->rule_list_dlg, dlg->rule_list, dlg->font);
  if (dlg->rule_list != NULL) {
    dlg->curr_item = 1;
    PointerToDialog (dlg->rule_editor, dlg->rule_list);
  }
}


static Pointer CommentSetFromDialog (DialoG d)
{
  CommentSetDlgPtr dlg;

  dlg = (CommentSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  return dlg->rule_list;
}


static void UpdateRuleList (Pointer data)
{
  CommentSetDlgPtr dlg;
  CommentRulePtr new_cr, last_cr = NULL;

  dlg = (CommentSetDlgPtr) data;
  if (dlg == NULL) {
    return;
  }

  new_cr = DialogToPointer (dlg->rule_editor);

  ReplaceNthCommentRule (&(dlg->rule_list), new_cr, dlg->curr_item);
  ListRules (dlg->rule_list_dlg, dlg->rule_list, dlg->font);  
}


static void CleanupCommentSetDialog (GraphiC g, VoidPtr data)

{
  CommentSetDlgPtr dlg;

  dlg = (CommentSetDlgPtr) data;
  if (dlg != NULL) {
    dlg->rule_list = CommentSetFree (dlg->rule_list);
  }
  StdCleanupFormProc (g, data);
}


NLM_EXTERN DialoG CommentSetDialog (GrouP h)
{
  CommentSetDlgPtr dlg;
  GrouP             p;
  
  dlg = (CommentSetDlgPtr) MemNew (sizeof (CommentSetDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupCommentSetDialog);

  dlg->dialog = (DialoG) p;
  dlg->todialog = CommentSetToDialog;
  dlg->fromdialog = CommentSetFromDialog;

  SetupCommentSetDialogFont (dlg);
  dlg->rule_list_dlg = DocumentPanel (p, stdCharWidth * 50, stdLineHeight * 12);
  SetObjectExtra (dlg->rule_list_dlg, dlg, NULL);
  SetDocProcs (dlg->rule_list_dlg, ClickRulesDoc, NULL, NULL, NULL);
  SetDocShade (dlg->rule_list_dlg, NULL, NULL, RuleHighlight, NULL);

  dlg->rule_editor = CommentRuleDialog (p, UpdateRuleList, dlg);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->rule_list_dlg, (HANDLE) dlg->rule_editor, NULL);
  
  return (DialoG) p;
}


static void SaveRuleChanges (ButtoN b)
{
  CommentSetDlgPtr    dlg;

  dlg = (CommentSetDlgPtr) GetObjectExtra (b);
}


NLM_EXTERN void LaunchCommentRulesEditor (IteM i)
{
  BaseFormPtr         bfp;
  WindoW              w;
  DialoG              d;
  CommentSetDlgPtr    dlg;
  GrouP               h, c;
  MenU                m;
  ButtoN              b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dlg = (CommentSetDlgPtr) MemNew (sizeof (CommentSetDlgData));
  if (dlg == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Comment Rules Editor", StdCloseWindowProc);

#if 0
  m = PulldownMenu (w, "File");
  FormCommandItem (m, "Open", (BaseFormPtr)f, VIB_MSG_OPEN);
  FormCommandItem (m, "Import", (BaseFormPtr)f, VIB_MSG_IMPORT);
  FormCommandItem (m, "Save", (BaseFormPtr)f, VIB_MSG_SAVE);
  FormCommandItem (m, "Save As", (BaseFormPtr)f, VIB_MSG_SAVE_AS);
  SeparatorItem (m);
  FormCommandItem (m, "Quit", (BaseFormPtr)f, VIB_MSG_QUIT);
  m = PulldownMenu (w, "Edit");
  FormCommandItem (m, "Copy All to Clipboard", (BaseFormPtr) f, VIB_MSG_COPY);
#endif
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  d = CommentSetDialog (h);
  dlg = GetObjectExtra (d);

  PointerToDialog (d, LoadCommentRuleSet());

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Accept", SaveRuleChanges);
  SetObjectExtra (b, d, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) d, (HANDLE) c, NULL);
  Show (w);
}


