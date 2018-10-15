/*   macrodlg.c
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
* File Name:  macrodlg.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   11/23/2007
*
* $Revision: 1.290 $
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

#include <objseq.h>
#include <objfdef.h>
#include <sqnutils.h>
#include <vibforms.h>
#include <document.h>
#include <dlogutil.h>
#include <cdrgn.h>
#include <seqpanel.h>
#include <biosrc.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>
#include <macrodlg.h>

/* macro editor dialog */
typedef struct macroeditorform {
  FORM_MESSAGE_BLOCK
  DoC    macro_summary;
  ButtoN run_btn;
  
  ValNodePtr macro_list;
  CharPtr    last_filename;
  Boolean    unsaved;
  Boolean    indexer_version;
  FonT       summary_font;

  MacroCloseCallback close_callback;
  Pointer            callback_data;
} MacroEditorFormData, PNTR MacroEditorFormPtr;

static Boolean SaveMacroFile (ForM f, CharPtr filename);

static Boolean BeforeCloseMacroEditor (MacroEditorFormPtr f)
{
  if (f != NULL) {
    if (f->unsaved) {
      if (Message (MSG_YN, "Do you want to save changes to the current file?") == ANS_YES) {
        if (!SaveMacroFile(f->form, f->last_filename)) {
          return FALSE;
        }
      }
    }
    if (f->close_callback != NULL) {
      (f->close_callback)(f->last_filename, f->callback_data);
    }
  }
  return TRUE;
}


extern void CloseMacroEditorWindowProc (WindoW w)

{
  MacroEditorFormPtr f;
  Boolean do_remove = FALSE;

  f = (MacroEditorFormPtr) GetObjectExtra (w);
  if (BeforeCloseMacroEditor(f)) {
    Remove (w);
  }
}


static void CleanupMacroEditorForm (GraphiC g, VoidPtr data)

{
  MacroEditorFormPtr f;

  f = (MacroEditorFormPtr) data;
  if (f != NULL) {
    f->macro_list = MacroActionListFree (f->macro_list);
    f->last_filename = MemFree (f->last_filename);
  }
  StdCleanupFormProc (g, data);
}


static void SetupMacroEditorFont (MacroEditorFormPtr f)

{
  if (f == NULL) return;

#ifdef WIN_MAC
  f->summary_font = ParseFont ("Times,12");
#endif
#ifdef WIN_MSWIN
  f->summary_font = ParseFont ("Times New Roman,12");
#endif
#ifdef WIN_MOTIF
  f->summary_font = ParseFont ("Times,12");
#endif
}


static void SetMacroEditorFileItems (MacroEditorFormPtr mefp);

static Boolean SaveMacroFile (ForM f, CharPtr filename)

{
  MacroEditorFormPtr mefp;
  Boolean            rval = FALSE;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;

  mefp = (MacroEditorFormPtr) GetObjectExtra (f);
  if (mefp != NULL) {
    path [0] = '\0';
    StringNCpy_0 (path, filename, sizeof (path));
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
      aip = AsnIoOpen (path, "w");
      if (aip == NULL) {
        Message (MSG_ERROR, "Unable to open %s", path);
      } else {
        MacroActionListAsnWrite (mefp->macro_list, aip, NULL);
        AsnIoClose (aip);
        mefp->last_filename = MemFree (mefp->last_filename);
        mefp->last_filename = StringSave (path);
        mefp->unsaved = FALSE;
        rval = TRUE;
      }
    }
  }
  return rval;
}



static CharPtr SummarizeMacroAction (ValNodePtr vnp);

static const Int4 kNumberColumnWidth = 40;

static void SummarizeMacro (DoC doc, ValNodePtr macro_list, FonT font)
{
  ValNodePtr vnp;
  CharPtr    str;
  CharPtr    tmp;
  RecT       r;
  Int4       pos;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {kNumberColumnWidth, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };


  Reset (doc);

  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  
  ColFmt[5].pixWidth = stdCharWidth;
  ColFmt[6].pixWidth = r.right - r.left - stdCharWidth - 48 - kNumberColumnWidth;

  if (font == NULL) font = programFont;

  if (macro_list == NULL) {
    AppendText (doc, "(Click here to start a new macro script)", NULL, NULL, font);
  } else {
    AppendText (doc, "(Click here to insert an action at the beginning of the script)", NULL, NULL, font);
    for (vnp = macro_list, pos = 1; vnp != NULL; vnp = vnp->next, pos++) {
      str = SummarizeMacroAction (vnp);
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 30));
      sprintf (tmp, "%d\t\t\t\t\tC\t%s\n", pos, str);
      str = MemFree (str);
      AppendText (doc, tmp, &ParFmt, ColFmt, font);
      tmp = MemFree (tmp);
    }
    AppendText (doc, "(Click here to insert an action at the end of the script)", NULL, NULL, font);
  }
  UpdateDocument (doc, 0, 0);
}


static Boolean OpenMacroFile (ForM f, CharPtr filename)

{
  MacroEditorFormPtr mefp;
  Boolean            rval = FALSE;
  CharPtr            extension;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;
  ValNodePtr         action_list;

  mefp = (MacroEditorFormPtr) GetObjectExtra (f);
  if (mefp != NULL) {
    if (mefp->unsaved) {
      if (Message (MSG_YN, "Do you want to save changes to the current file?") == ANS_YES) {
        if (!SaveMacroFile(f, mefp->last_filename)) {
           return FALSE;
        }
      }
    }
    path [0] = '\0';
    StringNCpy_0 (path, filename, sizeof (path));
    extension = NULL;
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip == NULL) {
        Message (MSG_ERROR, "Unable to open %s", path);
      } else {
        action_list = MacroActionListAsnRead (aip, NULL);
        if (action_list == NULL) {
          Message (MSG_ERROR, "Unable to read action list from %s.", path);
        } else {
          mefp->macro_list = MacroActionListFree (mefp->macro_list);
          mefp->macro_list = action_list;
          mefp->last_filename = MemFree (mefp->last_filename);
          mefp->last_filename = StringSave (path);
          mefp->unsaved = FALSE;
          rval = TRUE;
          SetMacroEditorFileItems (mefp);
          SummarizeMacro (mefp->macro_summary, mefp->macro_list, mefp->summary_font);
          SafeEnable (mefp->run_btn);
        }
        AsnIoClose (aip);
      }
    }    
  }
  return rval;
}


static Boolean ImportMacroFile (ForM f, CharPtr filename)

{
  MacroEditorFormPtr mefp;
  Boolean            rval = FALSE;
  CharPtr            extension;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;
  ValNodePtr         action_list;

  mefp = (MacroEditorFormPtr) GetObjectExtra (f);
  if (mefp == NULL) return FALSE;
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  extension = NULL;
  if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
    aip = AsnIoOpen (path, "r");
    if (aip == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
    } else {
      action_list = MacroActionListAsnRead (aip, NULL);
      if (action_list == NULL) {
        Message (MSG_ERROR, "Unable to read action list from %s.", path);
      } else {
        ValNodeLink (&(mefp->macro_list), action_list);
        mefp->unsaved = FALSE;
        rval = TRUE;
        SetMacroEditorFileItems (mefp);
        SummarizeMacro (mefp->macro_summary, mefp->macro_list, mefp->summary_font);
        SafeEnable (mefp->run_btn);
      }
      AsnIoClose (aip);
    }
  }
  return rval;
}


static void SetMacroEditorFileItems (MacroEditorFormPtr mefp)

{
  IteM  i;

  if (mefp != NULL) {
    i = FindFormMenuItem ((BaseFormPtr) mefp, VIB_MSG_OPEN);
    SafeSetTitle (i, "Load Macro File");
    SafeEnable (i);

    i = FindFormMenuItem ((BaseFormPtr) mefp, VIB_MSG_IMPORT);
    SafeSetTitle (i, "Add Macros from File to Script");
    SafeEnable (i);

    i = FindFormMenuItem ((BaseFormPtr) mefp, VIB_MSG_SAVE);
    SafeSetTitle (i, "Save Macro File");
    if (mefp->macro_list == NULL) {
      SafeDisable (i);
    } else {
      SafeEnable (i);
    }
    i = FindFormMenuItem ((BaseFormPtr) mefp, VIB_MSG_SAVE_AS);
    SafeSetTitle (i, "Save Macro File As");
    if (mefp->macro_list == NULL) {
      SafeDisable (i);
    } else {
      SafeEnable (i);
    }
  }
}


static void MacroEditorFormMessage (ForM f, Int2 mssg)

{
  MacroEditorFormPtr  mefp;
  ValNodePtr          vnp, string_list = NULL;
  Int4                len = 1;
  CharPtr             str;

  mefp = (MacroEditorFormPtr) GetObjectExtra (f);
  if (mefp != NULL) {
    switch (mssg) {
      case VIB_MSG_OPEN :
        OpenMacroFile (f, NULL);
        break;
      case VIB_MSG_IMPORT :
        ImportMacroFile (f, NULL);
        break;
      case VIB_MSG_SAVE :
        SaveMacroFile (f, mefp->last_filename);
        break;
      case VIB_MSG_SAVE_AS :
        SaveMacroFile (f, NULL);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        /* get length of entire macro list */
        for (vnp = mefp->macro_list; vnp != NULL; vnp = vnp->next) {
          str = SummarizeMacroAction (vnp);
          len += StringLen (str) + 2;
          ValNodeAddPointer (&string_list, 0, str);
        }
        str = (CharPtr) MemNew (sizeof (Char) * len);
        str[0] = 0;
        for (vnp = string_list; vnp != NULL; vnp = vnp->next) {
          StringCat (str, vnp->data.ptrvalue);
          StringCat (str, "\r\n");
        }
        string_list = ValNodeFreeData (string_list);
        StringToClipboard (str);
        str = MemFree (str);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        break;
      case VIB_MSG_CLOSE:
      case VIB_MSG_QUIT:
        CloseMacroEditorWindowProc ((WindoW)f);
        break;
      default :
        break;
    }
  }
}


static void RunMacro (ButtoN b)
{
  MacroEditorFormPtr  f;
  SeqEntryPtr         sep;
  ValNodePtr          sep_list, vnp;
  Uint2               entityID;
  LogInfoPtr          lip;

  f = (MacroEditorFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  if (f->macro_list == NULL) {
    Message (MSG_ERROR, "No macro loaded!");
  } else {
    sep_list = GetViewedSeqEntryList ();
    if (sep_list == NULL) {
      Message (MSG_ERROR, "No records open!");
    } else if (sep_list->next != NULL 
      && ANS_CANCEL == Message (MSG_OKC, "You have more than one record open - run macro for all open records?")) {
      /* do nothing */
    } else {
      WatchCursor();
      Update();
      lip = OpenLog ("Macro Actions");
      for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
        sep = vnp->data.ptrvalue;
        entityID = ObjMgrGetEntityIDForChoice(sep);
        lip->data_in_log |= ApplyMacroToSeqEntryEx (sep, f->macro_list, lip->fp, Sequin_GlobalAlign2Seq);
        ObjMgrSetDirtyFlag (entityID, TRUE);
        ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
      }
      sep_list = ValNodeFree (sep_list);
      ArrowCursor ();
      Update ();   
      if (!lip->data_in_log) {
        fprintf (lip->fp, "Macro had no effect\n");
        lip->data_in_log = TRUE;
      }
      CloseLog (lip);
      lip = FreeLog (lip);
    }
  }
}

static Boolean EditMacroAction (ValNodePtr action, Boolean indexer_version);

static Boolean EditMacroItem (MacroEditorFormPtr f, Int2 item)
{
  Boolean rval = FALSE;
  ValNodePtr vnp;

  if (f == NULL) return FALSE;

  for (vnp = f->macro_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {}
  if (vnp != NULL) {
    rval = EditMacroAction (vnp, f->indexer_version);
  }
  return rval;
}

static AECRActionPtr BuildDefaultAECRAction (Uint1 action_type, Uint1 qual_type)
{
  AECRActionPtr action = AECRActionNew();
  ApplyActionPtr apply;
  EditActionPtr  edit;
  ConvertActionPtr convert;
  CopyActionPtr copy;
  SwapActionPtr swap;
  AECRParseActionPtr parse;
  RemoveActionPtr remove;

  switch (action_type) {
    case ActionChoice_apply:
      apply = ApplyActionNew();
      apply->field = ValNodeNew (NULL);
      apply->field->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = apply;
      break;
    case ActionChoice_edit:
      edit = EditActionNew ();
      edit->field = ValNodeNew (NULL);
      edit->field->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = edit;
      break;
    case ActionChoice_convert:
      convert = ConvertActionNew();
      convert->fields = ValNodeNew (NULL);
      convert->fields->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = convert;
      break;
    case ActionChoice_copy:
      copy = CopyActionNew();
      copy->fields = ValNodeNew (NULL);
      copy->fields->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = copy;
      break;
    case ActionChoice_swap:
      swap = SwapActionNew();
      swap->fields = ValNodeNew (NULL);
      swap->fields->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = swap;
      break;
    case ActionChoice_remove:
      remove = RemoveActionNew();
      remove->field = ValNodeNew (NULL);
      remove->field->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = remove;
      break;
    case ActionChoice_parse:
      parse = AECRParseActionNew();
      parse->fields = ValNodeNew (NULL);
      parse->fields->choice = qual_type;
      action->action = ValNodeNew (NULL);
      action->action->choice = action_type;
      action->action->data.ptrvalue = parse;
      break;
  }

  return action;
}


static ApplyFeatureActionPtr BuildDefaultApplyFeatureAction (Uint2 feat_type)
{
  ApplyFeatureActionPtr action;
  FeatQualLegalValPtr   qual;

  action = ApplyFeatureActionNew ();
  action->type = feat_type;
  action->location = ValNodeNew (NULL);
  action->location->choice = LocationChoice_whole_sequence;
  action->seq_list = ValNodeNew (NULL);
  action->seq_list->choice = SequenceListChoice_all;
  qual = FeatQualLegalValNew ();
  qual->qual = Feat_qual_legal_codon_start;
  qual->val = StringSave ("best");
  action->fields = ValNodeNew (NULL);
  action->fields->choice = FeatQualLegalValChoice_qual;
  action->fields->data.ptrvalue = qual;
  return action;
}


static EditFeatureLocationActionPtr BuildDefaultEditFeatureLocationAction (Uint2 feat_type)
{
  EditFeatureLocationActionPtr a;
  EditLocationStrandPtr strand;

  a = EditFeatureLocationActionNew ();
  a->type = feat_type;
  a->action = ValNodeNew (NULL);
  a->action->choice = LocationEditType_strand;
  strand = EditLocationStrandNew ();
  strand->strand_to = Feature_location_strand_to_reverse;
  a->action->data.ptrvalue = strand;
  return a;
}


static RemoveFeatureActionPtr BuildDefaultRemoveFeatureAction (void)
{
  RemoveFeatureActionPtr remove;

  remove = RemoveFeatureActionNew ();
  return remove;
}


static ValNodePtr BuildDefaultNewMacroAction (void)
{
  ValNodePtr vnp;
  ApplyFeatureActionPtr apply_feat;

  apply_feat = BuildDefaultApplyFeatureAction (Macro_feature_type_cds);

  vnp = ValNodeNew (NULL);
  vnp->choice = MacroActionChoice_add_feature;
  vnp->data.ptrvalue = apply_feat;
  return vnp;
}


static void UpdateMacroSummary (MacroEditorFormPtr f, Int4 scroll_pos)
{
  Int4 scroll_max;
  BaR  sb_vert;

  if (f == NULL) return;

  f->unsaved = TRUE;
  SetMacroEditorFileItems (f);
  SummarizeMacro (f->macro_summary, f->macro_list, f->summary_font);
  if (f->macro_list == NULL) {
    SafeDisable (f->run_btn);
  } else {
    SafeEnable (f->run_btn);
  }
  if (scroll_pos > 0) {
    sb_vert = GetSlateVScrollBar ((SlatE) f->macro_summary);
    scroll_max = GetBarMax (sb_vert);
    if (scroll_pos > scroll_max) {
      scroll_pos = scroll_max;
    }
    CorrectBarValue (sb_vert, scroll_pos);
  }
}


static void AddMacroActions (MacroEditorFormPtr f, Int2 item);


static Boolean CloneMacroItem (MacroEditorFormPtr f, Int2 item)
{
  ValNodePtr            this_action = NULL, new_action;
  Int2                  pos;
  Int4                  scroll_pos;
  BaR                   sb_vert;

  if (f == NULL || item == 0 || f->macro_list == NULL) {
    return FALSE;
  }

  pos = 1;
  this_action = f->macro_list;
  while (pos < item && this_action != NULL) {
    pos++;
    this_action = this_action->next;
  }
  if (this_action == NULL) {
    return FALSE;
  }
  new_action = AsnIoMemCopy (this_action, (AsnReadFunc) MacroActionChoiceAsnRead, (AsnWriteFunc) MacroActionChoiceAsnWrite);
  if (new_action != NULL) {
    new_action->next = this_action->next;
    this_action->next = new_action;
  }

  /* get current scroll position */
  sb_vert = GetSlateVScrollBar ((SlatE) f->macro_summary);
  scroll_pos = GetBarValue (sb_vert);
  /* we will want to increase the scroll position after each addition
   * note that we need to get the scroll bar and check the initial position
   * each time - if there was no scroll bar after the last update, scroll_pos
   * needs to be zero to start.
   */
  scroll_pos++;

  /* update summary */
  UpdateMacroSummary (f, scroll_pos);
  return TRUE;
}


static void ClickMacroDoc (DoC d, PoinT pt)
{
  Int2               item, row, col;
  RecT               rct;
  MacroEditorFormPtr f;
  ValNodePtr         vnp, vnp_prev = NULL, two_prev = NULL, vnp_next;
  Boolean            changed_macro = FALSE;
  BaR                sb_vert = NULL;
  Int2               scroll_pos = 0;
  
  f = (MacroEditorFormPtr) GetObjectExtra (d);
  if (f == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item == 0 && row == 0 && f->macro_list == NULL) {
    AddMacroActions (f, 0);
  } else if (item > 0 && row > 0) {
    if (item == 1) {
      /* add to beginning of list */
      AddMacroActions (f, 0);
    } else if (item == ValNodeLen (f->macro_list) + 2) {
      /* add to end of list */
      AddMacroActions (f, item);
    } else {
      /* correct for explanatory line */
      item--;
      sb_vert = GetSlateVScrollBar ((SlatE) f->macro_summary);
      scroll_pos = GetBarValue (sb_vert);
      switch (col) {
        case 2:
          /* delete this item */
          for (vnp = f->macro_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
            vnp_prev = vnp;
          }
          if (vnp != NULL) {
            if (vnp_prev == NULL) {
              f->macro_list = vnp->next;
            } else {
              vnp_prev->next = vnp->next;
            }
            vnp->next = NULL;
            vnp = MacroActionListFree (vnp);
            changed_macro = TRUE;
          }
          break;
        case 3:
          /* move this item up */
          for (vnp = f->macro_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
            two_prev = vnp_prev;
            vnp_prev = vnp;
          }
          if (vnp != NULL && vnp_prev != NULL) {
            vnp_prev->next = vnp->next;
            vnp->next = vnp_prev;
            if (two_prev == NULL) {
              f->macro_list = vnp;
            } else {
              two_prev->next = vnp;
            }
            /* decrease the scroll position, so cursor will still be over the same item */
            scroll_pos--;
            changed_macro = TRUE;
          }
          break;
        case 4:
          /* move this item down */
          for (vnp = f->macro_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
            vnp_prev = vnp;
          }
          if (vnp != NULL && vnp->next != NULL) {
            vnp_next = vnp->next;
            vnp->next = vnp_next->next;
            vnp_next->next = vnp;
            if (vnp_prev == NULL) {
              f->macro_list = vnp_next;
            } else {
              vnp_prev->next = vnp_next;
            }
            /* increase the scroll position, so cursor will still be over the same item */
            scroll_pos++;
            changed_macro = TRUE;
          }
          break;
        case 5:
          /* insert item */
          if (pt.y >= rct.top && pt.y <= rct.top + 4) {
            /* insert macro before this one */
            AddMacroActions (f, item - 1);
          } else if (pt.y >= rct.bottom - 4 && pt.y <= rct.bottom) {
            /* insert macro after this one */
            AddMacroActions (f, item);
          }
          break;
        case 6:
          /* clone item */
          CloneMacroItem (f, item);
          changed_macro = TRUE;
          break;
        case 7:
          /* edit this item */
          changed_macro = EditMacroItem (f, item);
          break;
      }
    }
  }
  if (changed_macro) {
    UpdateMacroSummary (f, scroll_pos);
  }
}


static void DrawMacroDocControls (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  MacroEditorFormPtr dlg;
  RecT               rct;
  Int4               width;
  PoinT              pt1, pt2;

  dlg = (MacroEditorFormPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  /* don't draw controls for explanatory text */
  if (item == 1 || item >= ValNodeLen (dlg->macro_list) + 2) return;

  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;

    rct.left += kNumberColumnWidth;

    width = 10;
    /* draw X for deletion */
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1;
    pt2.x = pt1.x + width;
    pt2.y = pt1.y + width;
    DrawLine (pt1, pt2);
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1 + width;
    pt2.x = pt1.x + width;
    pt2.y = rct.top + 1;
    DrawLine (pt1, pt2);

    /* draw up arrow for moving step up */
    if (item > 2) {
      pt1.x = rct.left + width + 3;
      pt1.y = rct.top + 3;
      pt2.x = pt1.x + 5;
      pt2.y = rct.top + 1;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x + 5;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x;
      pt1.y = pt2.y + width;
      DrawLine (pt1, pt2);
    }
    /* draw up arrow for moving step up */
    if (item < ValNodeLen (dlg->macro_list) + 1) {
      pt1.x = rct.left + 2 * width + 5;
      pt1.y = rct.top + width - 2;
      pt2.x = pt1.x + 5;
      pt2.y = rct.top + width + 1;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x + 5;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x;
      pt1.y = rct.top + 1;
      DrawLine (pt1, pt2);
    }

    /* draw insertion controls */
    pt1.x = rct.left + 3 * width + 7;
    pt1.y = rct.top + 4;
    pt2.x = pt1.x + width;
    pt2.y = rct.top;
    if (item > 2) {
      DrawLine (pt1, pt2);
    }
    pt1.y = rct.bottom - 4;
    pt2.y = rct.bottom;
    if (item < ValNodeLen (dlg->macro_list) + 1) {
      DrawLine (pt1, pt2);
    }
  }
}


NLM_EXTERN void LaunchMacroEditor (Uint2 entityID, CharPtr filename, MacroCloseCallback close_callback, Pointer callback_data)
{
  WindoW              w;
  MacroEditorFormPtr  f;
  GrouP               h, c = NULL;
  MenU                m;
#ifdef TEST_MACRO_TEMPLATE_EDITOR
  ButtoN              b;
#endif

  f = (MacroEditorFormPtr) MemNew (sizeof (MacroEditorFormData));
  if (f == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Macro Editor", CloseMacroEditorWindowProc);
  SetObjectExtra (w, f, CleanupMacroEditorForm);
  f->form = (ForM) w;
  f->input_entityID = entityID;

  f->formmessage = MacroEditorFormMessage;

  f->macro_list = NULL;
  f->last_filename = NULL;
  f->unsaved = FALSE;
  f->indexer_version = TRUE;
  f->close_callback = close_callback;
  f->callback_data = callback_data;

  m = PulldownMenu (w, "File");
  FormCommandItem (m, "Open", (BaseFormPtr)f, VIB_MSG_OPEN);
  FormCommandItem (m, "Import", (BaseFormPtr)f, VIB_MSG_IMPORT);
  FormCommandItem (m, "Save", (BaseFormPtr)f, VIB_MSG_SAVE);
  FormCommandItem (m, "Save As", (BaseFormPtr)f, VIB_MSG_SAVE_AS);
  SeparatorItem (m);
  FormCommandItem (m, "Quit", (BaseFormPtr)f, VIB_MSG_QUIT);
  m = PulldownMenu (w, "Edit");
  FormCommandItem (m, "Copy All to Clipboard", (BaseFormPtr) f, VIB_MSG_COPY);

  SetupMacroEditorFont (f);
  SetMacroEditorFileItems (f);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  f->macro_summary = DocumentPanel (h, stdCharWidth * 50, stdLineHeight * 20);
  SetObjectExtra (f->macro_summary, f, NULL);
  SetDocProcs (f->macro_summary, ClickMacroDoc, NULL, NULL, NULL);
  SetDocShade (f->macro_summary, DrawMacroDocControls, NULL, NULL, NULL);

  if (f->input_entityID != 0) {
    c = HiddenGroup (h, 3, 0, NULL);
    f->run_btn = PushButton (c, "Run", RunMacro);
    SetObjectExtra (f->run_btn, f, NULL);
    Disable (f->run_btn);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) f->macro_summary, (HANDLE) c, NULL);
  SummarizeMacro (f->macro_summary, f->macro_list, f->summary_font);
  if (!StringHasNoText (filename)) {
    OpenMacroFile (f->form, filename);
  }
  Show (w);
}


NLM_EXTERN void LaunchMacroEditorMenuItem (IteM i)
{
  BaseFormPtr         bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  LaunchMacroEditor (bfp->input_entityID, NULL, NULL, NULL);
}


static void ClearDialogBtn (ButtoN b)
{
  DialoG d;
  
  d = (DialoG) GetObjectExtra (b);
  
  PointerToDialog (d, NULL);
}


typedef struct textmarkerdialog
{
  DIALOG_MESSAGE_BLOCK

  GrouP marker_choice;
  TexT  marker_text;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} TextMarkerDialogData, PNTR TextMarkerDialogPtr;


static void ResetTextMarkerDialog (DialoG d)
{
  TextMarkerDialogPtr tp;

  tp = (TextMarkerDialogPtr) GetObjectExtra (d);
  if (tp == NULL) {
    return;
  }
  SetValue (tp->marker_choice, 1);
  SetTitle (tp->marker_text, "");
  Enable (tp->marker_text);
}


static void TextMarkerToDialog(DialoG d, Pointer data)
{
  TextMarkerDialogPtr tp;
  TextMarkerPtr       tdata;

  tp = (TextMarkerDialogPtr) GetObjectExtra (d);
  if (tp == NULL) {
    return;
  }

  tdata = (TextMarkerPtr) data;

  if (tdata == NULL) {
    ResetTextMarkerDialog(d);
  } else {
    if (tdata->choice == TextMarker_free_text) {
      SetValue (tp->marker_choice, 1);
      SetTitle (tp->marker_text, tdata->data.ptrvalue == NULL ? "" : tdata->data.ptrvalue);
      Enable (tp->marker_text);
    } else if (tdata->choice == TextMarker_digits) {
      SetValue (tp->marker_choice, 3);
      SetTitle (tp->marker_text, "");
      Disable (tp->marker_text);
    } else if (tdata->choice == TextMarker_letters) {
      SetValue (tp->marker_choice, 4);
      SetTitle (tp->marker_text, "");
      Disable (tp->marker_text);
    } else {
      SetValue (tp->marker_choice, 1);
      SetTitle (tp->marker_text, "");
      Enable (tp->marker_text);
    }
  }
}


static Pointer DialogToTextMarker (DialoG d)
{
  TextMarkerDialogPtr tp;
  Int2                val;
  ValNodePtr          tdata = NULL;

  tp = (TextMarkerDialogPtr) GetObjectExtra (d);
  if (tp == NULL) {
    return NULL;
  }

  val = GetValue (tp->marker_choice);
  if (val == 1) {
    ValNodeAddPointer (&tdata, TextMarker_free_text, JustSaveStringFromText (tp->marker_text));
  } else if (val == 3) {
    tdata = ValNodeNew (NULL);
    tdata->choice = TextMarker_digits;
  } else if (val == 4) {
    tdata = ValNodeNew (NULL);
    tdata->choice = TextMarker_letters;
  }

  if (IsTextMarkerEmpty(tdata)) {
    tdata = TextMarkerFree (tdata);
  }
  return tdata;
}


static void TextMarkerMessage (DialoG d, Int2 mssg)

{
  TextMarkerDialogPtr tp;

  tp = (TextMarkerDialogPtr) GetObjectExtra (d);
  if (tp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetTextMarkerDialog (d);        
        break;
      case VIB_MSG_ENTER :
        Select (tp->marker_text);
        break;
      default :
        break;
    }
  }
}


static void ChangeTextMarkerChoice (GrouP g)
{
  TextMarkerDialogPtr tp;
  Int2 val;

  tp = (TextMarkerDialogPtr) GetObjectExtra (g);
  if (tp == NULL) {
    return;
  }

  val = GetValue (tp->marker_choice);
  if (val == 1) {
    Enable (tp->marker_text);
  } else {
    Disable (tp->marker_text);
  }

  if (tp->change_notify != NULL) {
    (tp->change_notify) (tp->change_userdata);
  }
}


static void ChangeTextMarkerText (TexT t)
{
  TextMarkerDialogPtr tp;

  tp = (TextMarkerDialogPtr) GetObjectExtra (t);
  if (tp == NULL) {
    return;
  }

  if (tp->change_notify != NULL) {
    (tp->change_notify) (tp->change_userdata);
  }
}


static DialoG TextMarkerDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TextMarkerDialogPtr tp;
  GrouP               p;
  
  tp = (TextMarkerDialogPtr) MemNew (sizeof (TextMarkerDialogData));
  if (tp == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, tp, StdCleanupExtraProc);

  tp->dialog = (DialoG) p;
  tp->todialog = TextMarkerToDialog;
  tp->fromdialog = DialogToTextMarker;
  tp->dialogmessage = TextMarkerMessage;
  tp->testdialog = NULL;

  tp->change_notify = change_notify;
  tp->change_userdata = change_userdata;

  tp->marker_choice = HiddenGroup (p, 4, 0, ChangeTextMarkerChoice);
  SetObjectExtra (tp->marker_choice, tp, NULL);
  SetGroupSpacing (tp->marker_choice, 10, 10);

  RadioButton (tp->marker_choice, "Text");
  tp->marker_text = DialogText (tp->marker_choice, "", 10, ChangeTextMarkerText);
  SetObjectExtra (tp->marker_text, tp, NULL);
  RadioButton (tp->marker_choice, "Digits");
  RadioButton (tp->marker_choice, "Letters");
  SetValue (tp->marker_choice, 1);

  return (DialoG) p;
}


typedef struct textportiondialog
{
  DIALOG_MESSAGE_BLOCK

  GrouP  start_choice;
  DialoG start_marker;
  GrouP  end_choice;
  DialoG end_marker;
  ButtoN rem_before;
  ButtoN also_rem_before;
  ButtoN rem_after;
  ButtoN also_rem_after;
  ButtoN insensitive;
  ButtoN whole_word;

  Boolean inside;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} TextPortionDialogData, PNTR TextPortionDialogPtr;


static void OutsideEnableDisable (TextPortionDialogPtr tp)
{
  if (tp == NULL || tp->inside)
  {
    return;
  }
  if (GetStatus (tp->rem_before)) 
  {
    Enable (tp->start_marker);
    Enable (tp->also_rem_before);
  }
  else
  {
    Disable (tp->start_marker);
    Disable (tp->also_rem_before);
  }
  if (GetStatus (tp->rem_after)) 
  {
    Enable (tp->end_marker);
    Enable (tp->also_rem_after);
  }
  else
  {
    Disable (tp->end_marker);
    Disable (tp->also_rem_after);
  }
}

static void ResetTextPortionDialog (TextPortionDialogPtr tp)
{
  if (tp == NULL)
  {
    return;
  }
  if (tp->inside)
  {
    SetValue (tp->start_choice, 1);
    SetValue (tp->end_choice, 1);
  }
  else
  {
    SetStatus (tp->rem_before, FALSE);
    SetStatus (tp->also_rem_before, FALSE);
    SetStatus (tp->rem_after, FALSE);
    SetStatus (tp->also_rem_after, FALSE);
  }
  PointerToDialog (tp->start_marker, NULL);
  PointerToDialog (tp->end_marker, NULL);
  SetStatus (tp->insensitive, FALSE);
  SetStatus (tp->whole_word, FALSE);
  OutsideEnableDisable (tp);
}

static void TextPortionToDialog (DialoG d, Pointer data)
{
  TextPortionDialogPtr tdlg;
  TextPortionPtr       tdata;
  
  tdlg = (TextPortionDialogPtr) GetObjectExtra (d);
  if (tdlg == NULL) {
    return;
  }
  tdata = (TextPortionPtr) data;
  ResetTextPortionDialog (tdlg);  
  if (tdata != NULL) {
    if (tdlg->inside) {
      if (tdata->include_left) {        
        SetValue (tdlg->start_choice, 2);
      } else {
        SetValue (tdlg->start_choice, 1);
      }
      if (tdata->include_right) {
        SetValue (tdlg->end_choice, 2);
      } else {
        SetValue (tdlg->end_choice, 1);
      }
    } else {
      SetStatus (tdlg->rem_before, !IsTextMarkerEmpty(tdata->left_marker));
      SetStatus (tdlg->rem_after, !IsTextMarkerEmpty(tdata->right_marker));
      SetStatus (tdlg->also_rem_before, tdata->include_left);
      SetStatus (tdlg->also_rem_after, tdata->include_right);
    }
    PointerToDialog (tdlg->start_marker, tdata->left_marker);
    PointerToDialog (tdlg->end_marker, tdata->right_marker);
    SetStatus (tdlg->insensitive, !tdata->case_sensitive);
    SetStatus (tdlg->whole_word, tdata->whole_word);
  }
  OutsideEnableDisable (tdlg);
}

static Pointer DialogToTextPortion (DialoG d)
{
  TextPortionDialogPtr tdlg;
  TextPortionPtr       tdata;

  tdlg = (TextPortionDialogPtr) GetObjectExtra (d);
  if (tdlg == NULL) {
    return NULL;
  }

  tdata = TextPortionNew();
  if (tdata != NULL) {
    if (tdlg->inside) {
      if (GetValue (tdlg->start_choice) == 1) {
        tdata->include_left = FALSE;
      } else {
        tdata->include_left = TRUE;
      }
      tdata->left_marker = DialogToPointer (tdlg->start_marker);

      if (GetValue (tdlg->end_choice) == 1) {
        tdata->include_right = FALSE;
      } else {
        tdata->include_right = TRUE;
      }
      tdata->right_marker = DialogToPointer (tdlg->end_marker);
    } else {
      if (GetStatus (tdlg->rem_before)) {
        tdata->left_marker = DialogToPointer (tdlg->start_marker);
      } else {
        tdata->left_marker = NULL;
      }
      if (GetStatus (tdlg->rem_after)) {
        tdata->right_marker = DialogToPointer (tdlg->end_marker);
      } else {
        tdata->right_marker = NULL;
      }

      tdata->include_left = GetStatus (tdlg->also_rem_before);
      tdata->include_right = GetStatus (tdlg->also_rem_after);
    }
        
    tdata->case_sensitive = !GetStatus (tdlg->insensitive);
    tdata->whole_word = GetStatus (tdlg->whole_word);
    tdata->inside = tdlg->inside;
  }
  return tdata;
}

static void TextPortionMessage (DialoG d, Int2 mssg)

{
  TextPortionDialogPtr tp;

  tp = (TextPortionDialogPtr) GetObjectExtra (d);
  if (tp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetTextPortionDialog (tp);        
        break;
      case VIB_MSG_ENTER :
        Select (tp->start_marker);
        break;
      default :
        break;
    }
  }
}

static void ChangeTextPortionGroup (GrouP g)
{
  TextPortionDialogPtr tp;

  tp = (TextPortionDialogPtr) GetObjectExtra (g);
  if (tp == NULL) return;

  if (tp->change_notify != NULL)
  {
    (tp->change_notify) (tp->change_userdata);
  }
}


static void ChangeTextPortionBtn (ButtoN b)
{
  TextPortionDialogPtr tp;

  tp = (TextPortionDialogPtr) GetObjectExtra (b);
  if (tp == NULL) return;
  OutsideEnableDisable (tp);

  if (tp->change_notify != NULL)
  {
    (tp->change_notify) (tp->change_userdata);
  }
}


static ValNodePtr TestTextPortionDialog (DialoG d)
{
  TextPortionDialogPtr tp;
  ValNodePtr err_list = NULL;

  tp = (TextPortionDialogPtr) GetObjectExtra (d);
  if (tp == NULL) return NULL;

  /* don't actually need to fill anything in.  Could want to copy the entire field. */

  return err_list;
}


static GrouP GetTextPortionStartChoiceGroup (DialoG d)
{
  TextPortionDialogPtr dlg;

  dlg = (TextPortionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  } else {
    return dlg->start_choice;
  }
}


static GrouP GetTextPortionEndChoiceGroup (DialoG d)
{
  TextPortionDialogPtr dlg;

  dlg = (TextPortionDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  } else {
    return dlg->end_choice;
  }
}


NLM_EXTERN DialoG TextPortionDialog (GrouP h, Boolean inside, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TextPortionDialogPtr tp;
  GrouP                p, g1, g2;
  
  tp = (TextPortionDialogPtr) MemNew (sizeof (TextPortionDialogData));
  if (tp == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, tp, StdCleanupExtraProc);

  tp->dialog = (DialoG) p;
  tp->todialog = TextPortionToDialog;
  tp->fromdialog = DialogToTextPortion;
  tp->dialogmessage = TextPortionMessage;
  tp->testdialog = TestTextPortionDialog;

  tp->change_notify = change_notify;
  tp->change_userdata = change_userdata;
  tp->inside = inside;

  g1 = HiddenGroup (p, 3, 0, NULL);
  SetGroupSpacing (g1, 10, 10);

  if (inside) 
  {
    StaticPrompt (g1, "Between", 0, popupMenuHeight, programFont, 'r');
    tp->start_choice = HiddenGroup (g1, 2, 0, ChangeTextPortionGroup);
    RadioButton (tp->start_choice, "just after");
    RadioButton (tp->start_choice, "starting at");
    SetValue (tp->start_choice, 1);
    SetObjectExtra (tp->start_choice, tp, NULL);

    tp->start_marker = TextMarkerDialog (g1, tp->change_notify, tp->change_userdata);
    
    StaticPrompt (g1, "And", 0, popupMenuHeight, programFont, 'r');
    tp->end_choice = HiddenGroup (g1, 2, 0, ChangeTextPortionGroup);
    RadioButton (tp->end_choice, "up to");
    RadioButton (tp->end_choice, "including");
    SetValue (tp->end_choice, 1);
    SetObjectExtra (tp->end_choice, tp, NULL);
      
    tp->end_marker = TextMarkerDialog (g1, tp->change_notify, tp->change_userdata);
  }
  else
  {
    tp->rem_before = CheckBox (g1, "Before", ChangeTextPortionBtn);
    SetObjectExtra (tp->rem_before, tp, NULL);
    tp->start_marker = TextMarkerDialog (g1, tp->change_notify, tp->change_userdata);
    tp->also_rem_before = CheckBox (g1, "Also Remove Entered Text", ChangeTextPortionBtn);

    tp->rem_after = CheckBox (g1, "After", ChangeTextPortionBtn);
    SetObjectExtra (tp->rem_after, tp, NULL);
    tp->end_marker = TextMarkerDialog (g1, tp->change_notify, tp->change_userdata);
    tp->also_rem_after = CheckBox (g1, "Also Remove Entered Text", ChangeTextPortionBtn);    
  }
  
  g2 = HiddenGroup (p, 2, 0, NULL);
  tp->insensitive = CheckBox (g2, "Case insensitive", ChangeTextPortionBtn);
  SetObjectExtra (tp->insensitive, tp, NULL);
  tp->whole_word = CheckBox (g2, "Whole word", ChangeTextPortionBtn);
  SetObjectExtra (tp->whole_word, tp, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);
  
  ResetTextPortionDialog (tp);
  return (DialoG) p;
}


typedef struct fieldeditdlg {
  DIALOG_MESSAGE_BLOCK
  TexT find_txt;
  TexT repl_txt;
  GrouP location;
  ButtoN case_insensitive;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} FieldEditDlgData, PNTR FieldEditDlgPtr;


static void ChangeFieldEditText (TexT t)
{
  FieldEditDlgPtr dlg;

  dlg = (FieldEditDlgPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL) {
    dlg->change_notify (dlg->change_userdata);
  }
}


static void ChangeFieldEditGroup (GrouP g)
{
  FieldEditDlgPtr dlg;

  dlg = (FieldEditDlgPtr) GetObjectExtra (g);
  if (dlg != NULL && dlg->change_notify != NULL) {
    dlg->change_notify (dlg->change_userdata);
  }
}


static void FieldEditToDialog (DialoG d, Pointer data)
{
  FieldEditDlgPtr dlg;
  FieldEditPtr edit;

  dlg = (FieldEditDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) return;

  edit = (FieldEditPtr) data;
  if (edit == NULL) {
    SetTitle (dlg->find_txt, "");
    SetTitle (dlg->repl_txt, "");
    SetValue (dlg->location, 1);
    SetStatus (dlg->case_insensitive, FALSE);
  } else {
    SetTitle (dlg->find_txt, edit->find_txt);
    SetTitle (dlg->repl_txt, edit->repl_txt);
    switch (edit->location) {
      case Field_edit_location_anywhere:
        SetValue (dlg->location, 1);
        break;
      case Field_edit_location_beginning:
        SetValue (dlg->location, 2);
        break;
      case Field_edit_location_end:
        SetValue (dlg->location, 3);
        break;
      default:
        SetValue (dlg->location, 1);
        break;
    }
    SetStatus (dlg->case_insensitive, edit->case_insensitive);
  }
  if (dlg->change_notify != NULL) {
    dlg->change_notify (dlg->change_userdata);
  }
}


static Pointer DialogToFieldEdit (DialoG d)
{
  FieldEditDlgPtr dlg;
  FieldEditPtr edit;

  dlg = (FieldEditDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) return NULL;

  edit = FieldEditNew ();
  edit->find_txt = JustSaveStringFromText (dlg->find_txt);
  edit->repl_txt = JustSaveStringFromText (dlg->repl_txt);
  if ((edit->find_txt == NULL || edit->find_txt[0] == 0) && (edit->repl_txt == NULL || edit->repl_txt[0] == 0)) {
    edit = FieldEditFree (edit);
  } else {
    switch (GetValue (dlg->location)) {
      case 1:
        edit->location = Field_edit_location_anywhere;
        break;
      case 2:
        edit->location = Field_edit_location_beginning;
        break;
      case 3:
        edit->location = Field_edit_location_end;
        break;
      default:
        edit->location = Field_edit_location_anywhere;
        break;
    }
    edit->case_insensitive = GetStatus (dlg->case_insensitive);
  }
  return edit;
}


static void FieldEditDlgCopy (ButtoN b)
{
  FieldEditDlgPtr dlg;
  CharPtr         str = NULL;

  dlg = (FieldEditDlgPtr) GetObjectExtra (b);
  if (dlg == NULL)
  {
    return;
  }
  str = JustSaveStringFromText (dlg->find_txt);
  SetTitle (dlg->repl_txt, str);
  str = MemFree (str);
}

static ValNodePtr TestFieldEditDialog (DialoG d)
{
  FieldEditDlgPtr dlg;
  ValNodePtr err_list = NULL;
  CharPtr tmp;

  dlg = (FieldEditDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  tmp = JustSaveStringFromText (dlg->find_txt);
  if (tmp == NULL || tmp[0] == 0) {
    ValNodeAddPointer (&err_list, 0, "no find text");
  }
  tmp = MemFree (tmp);
  return err_list;
}


static void SetFieldEditDialogText (DialoG d, CharPtr str_find, CharPtr str_repl)
{
  FieldEditDlgPtr dlg;

  dlg = (FieldEditDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (str_find != NULL) {
      SetTitle (dlg->find_txt, str_find);
    }
    if (str_repl != NULL) {
      SetTitle (dlg->repl_txt, str_repl);
    }
  }
}


static DialoG FieldEditDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  FieldEditDlgPtr dlg;
  GrouP p, g;
  ButtoN b;

  p = HiddenGroup (h, -1, 0, NULL);
  dlg = (FieldEditDlgPtr) MemNew (sizeof (FieldEditDlgData));
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->fromdialog = DialogToFieldEdit;
  dlg->todialog = FieldEditToDialog;
  dlg->testdialog = TestFieldEditDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  g = HiddenGroup (p, 3, 0, NULL);
  StaticPrompt (g, "Find", 0, dialogTextHeight, systemFont, 'r');
  dlg->find_txt = DialogText (g, "", 18, ChangeFieldEditText);
  SetObjectExtra (dlg->find_txt, dlg, NULL);
  b = PushButton (g, "Copy", FieldEditDlgCopy);
  SetObjectExtra (b, dlg, NULL);
  Hide (b);
  StaticPrompt (g, "Replace", 0, dialogTextHeight, systemFont, 'r');
  dlg->repl_txt = DialogText (g, "", 18, ChangeFieldEditText);
  SetObjectExtra (dlg->repl_txt, dlg, NULL);
  b = PushButton (g, "Copy", FieldEditDlgCopy);
  SetObjectExtra (b, dlg, NULL);

  dlg->location = HiddenGroup (p, 3, 0, ChangeFieldEditGroup);
  SetObjectExtra (dlg->location, dlg, NULL);
  RadioButton (dlg->location, "Anywhere in field");
  RadioButton (dlg->location, "At the beginning of the field");
  RadioButton (dlg->location, "At the end of the field");
  SetValue (dlg->location, 1);

  dlg->case_insensitive = CheckBox (p, "Case-insensitive", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlg->location, (HANDLE) dlg->case_insensitive, NULL);

  return (DialoG) p;
}


typedef struct existingtextdlg {
  DIALOG_MESSAGE_BLOCK
  GrouP                action_grp;
  ButtoN               append_btn;
  ButtoN               prefix_btn;
  ButtoN               add_field_btn;
  GrouP                delim_grp;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} ExistingTextDlgData, PNTR ExistingTextDlgPtr;

static void ChangeExistingTextActionChoice (GrouP g)
{
  ExistingTextDlgPtr dlg;
  Int4               action_choice;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  
  action_choice = GetValue (dlg->action_grp);
  if (action_choice == 2 || action_choice == 3) {
    Enable (dlg->delim_grp);
  } else {
    Disable (dlg->delim_grp);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SetExistingTextDialogValue (DialoG d, Uint2 existing_text)
{
  ExistingTextDlgPtr dlg;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  switch (existing_text) {
    case ExistingTextOption_append_semi :
      SetValue (dlg->action_grp, 2);
      SetValue (dlg->delim_grp, 1);
      break;
    case ExistingTextOption_append_space :
      SetValue (dlg->action_grp, 2);
      SetValue (dlg->delim_grp, 2);
      break;
    case ExistingTextOption_append_colon :
      SetValue (dlg->action_grp, 2);
      SetValue (dlg->delim_grp, 3);
      break;
    case ExistingTextOption_append_comma :
      SetValue (dlg->action_grp, 2);
      SetValue (dlg->delim_grp, 4);
      break;
    case ExistingTextOption_append_none :
      SetValue (dlg->action_grp, 2);
      SetValue (dlg->delim_grp, 5);
      break;
    case ExistingTextOption_prefix_semi :
      SetValue (dlg->action_grp, 3);
      SetValue (dlg->delim_grp, 1);
      break;
    case ExistingTextOption_prefix_space :
      SetValue (dlg->action_grp, 3);
      SetValue (dlg->delim_grp, 2);
      break;
    case ExistingTextOption_prefix_colon :
      SetValue (dlg->action_grp, 3);
      SetValue (dlg->delim_grp, 3);
      break;
    case ExistingTextOption_prefix_comma :
      SetValue (dlg->action_grp, 3);
      SetValue (dlg->delim_grp, 4);
      break;
    case ExistingTextOption_prefix_none :
      SetValue (dlg->action_grp, 3);
      SetValue (dlg->delim_grp, 5);
      break;
    case ExistingTextOption_leave_old :
      SetValue (dlg->action_grp, 4);
      break;
    case ExistingTextOption_add_qual :
      SetValue (dlg->action_grp, 5);
      break;
    case ExistingTextOption_replace_old :
      SetValue (dlg->action_grp, 1);
      break;
    default:
      SetValue (dlg->action_grp, 1);
      SetValue (dlg->delim_grp, 1);
      break;
  }
  ChangeExistingTextActionChoice (dlg->action_grp);
}


static void ChangeExistingTextDelimChoice (GrouP g)
{
  ExistingTextDlgPtr dlg;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Uint2 GetExistingTextDialogValue (DialoG d)
{
  ExistingTextDlgPtr dlg;
  Int4 action_choice, separator_choice;
  Uint2 existing_text = ExistingTextOption_replace_old;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return existing_text;
  action_choice = GetValue (dlg->action_grp);
  separator_choice = GetValue (dlg->delim_grp);
  switch (action_choice) {
    case 1:
      existing_text = ExistingTextOption_replace_old;
      break;
    case 2:
      switch (separator_choice) {
        case 1:
          existing_text = ExistingTextOption_append_semi;
          break;
        case 2:
          existing_text = ExistingTextOption_append_space;
          break;
        case 3:
          existing_text = ExistingTextOption_append_colon;
          break;
        case 4:
          existing_text = ExistingTextOption_append_comma;
          break;
        case 5:
          existing_text = ExistingTextOption_append_none;
          break;
      }
      break;
    case 3:
      switch (separator_choice) {
        case 1:
          existing_text = ExistingTextOption_prefix_semi;
          break;
        case 2:
          existing_text = ExistingTextOption_prefix_space;
          break;
        case 3:
          existing_text = ExistingTextOption_prefix_colon;
          break;
        case 4:
          existing_text = ExistingTextOption_prefix_comma;
          break;
        case 5:
          existing_text = ExistingTextOption_prefix_none;
          break;
      }
      break;
    case 4:
      existing_text = ExistingTextOption_leave_old;
      break;
    case 5:
      existing_text = ExistingTextOption_add_qual;
      break;
  }
  return existing_text;
}


static ValNodePtr TestExistingTextDialog (DialoG d)
{
  ExistingTextDlgPtr dlg;
  ValNodePtr err_list = NULL;
  Int2       val;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->action_grp);
  if ((val == 2 && !Enabled (dlg->append_btn)) || (val == 3 && !Enabled (dlg->prefix_btn))) {
    ValNodeAddPointer (&err_list, 0, "invalid existing text option for nontext value");
  } else if (val == 5 && !Enabled (dlg->add_field_btn)) {
    ValNodeAddPointer (&err_list, 0, "invalid existing text option for field type");
  }
  return err_list;
}


static void EnableNonTextOptions (DialoG d)
{
  ExistingTextDlgPtr dlg;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  Enable (dlg->append_btn);
  Enable (dlg->prefix_btn);
}


static void DisableNonTextOptions (DialoG d)
{
  ExistingTextDlgPtr dlg;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  Disable (dlg->append_btn);
  Disable (dlg->prefix_btn);
}


static void EnableMultiOptions (DialoG d)
{
  ExistingTextDlgPtr dlg;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  Enable (dlg->add_field_btn);
}


static void DisableMultiOptions (DialoG d)
{
  ExistingTextDlgPtr dlg;

  dlg = (ExistingTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  Disable (dlg->add_field_btn);
}


static DialoG ExistingTextDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ExistingTextDlgPtr dlg;
  GrouP              p;
  PrompT             ppt;

  dlg = (ExistingTextDlgPtr) MemNew (sizeof (ExistingTextDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->testdialog = TestExistingTextDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->action_grp = HiddenGroup (p, 5, 0, ChangeExistingTextActionChoice);
  SetGroupSpacing (dlg->action_grp, 10, 10);
  SetObjectExtra (dlg->action_grp, dlg, NULL);
  RadioButton (dlg->action_grp, "Overwrite existing text");
  dlg->append_btn = RadioButton (dlg->action_grp, "Append");
  dlg->prefix_btn = RadioButton (dlg->action_grp, "Prefix");
  RadioButton (dlg->action_grp, "Ignore new text");
  dlg->add_field_btn = RadioButton (dlg->action_grp, "Add new qual");
  Disable (dlg->add_field_btn);
  SetValue (dlg->action_grp, 1);

  ppt = StaticPrompt (p, "Separate new text and old text with", 
                      0, dialogTextHeight, programFont, 'c');
  
  dlg->delim_grp = HiddenGroup (p, 5, 0, ChangeExistingTextDelimChoice);
  SetObjectExtra (dlg->delim_grp, dlg, NULL);
  SetGroupSpacing (dlg->delim_grp, 10, 10);
  RadioButton (dlg->delim_grp, "Semicolon");
  RadioButton (dlg->delim_grp, "Space");
  RadioButton (dlg->delim_grp, "Colon");
  RadioButton (dlg->delim_grp, "Comma");
  RadioButton (dlg->delim_grp, "Do not separate");
  SetValue (dlg->delim_grp, 1);
  Disable (dlg->delim_grp);
  return (DialoG) p;
}


/* ExistingText handling dialog and structures */
typedef struct existingtexttwostepdlg 
{
  GrouP pre_app_grp;
  GrouP delim_grp;
} ExistingTextTwoStepDlgData, PNTR ExistingTextTwoStepDlgPtr;

static void ChangePreAppIgnoreChoice (GrouP g)
{
  ExistingTextTwoStepDlgPtr etdp;
  Int4               handle_choice;
  
  etdp = (ExistingTextTwoStepDlgPtr) GetObjectExtra (g);
  if (etdp == NULL)
  {
    return;
  }
  
  handle_choice = GetValue (etdp->pre_app_grp);
  if (handle_choice == 1 || handle_choice == 2)
  {
    Enable (etdp->delim_grp);
  }
  else
  {
    Disable (etdp->delim_grp);
  }
}


NLM_EXTERN Uint2 TwoStepExistingText (Int4 num_found, Boolean non_text, Boolean allow_multi)
{
  WindoW                w;
  GrouP                 h, c;
  Uint2                 existing_text = 0;
  ButtoN                b;
  ExistingTextTwoStepDlgData   etdd;
  ModalAcceptCancelData acd;
  Char                  txt [128];
  MsgAnswer             ans;
  PrompT                ppt;
  Int4                  handle_choice;

  if (num_found <= 0)
  {
    return ExistingTextOption_replace_old;
  }
  
  sprintf (txt, "%d affected fields already contain a value.  Do you wish to overwrite existing text?",
           num_found);
  ans = Message (MSG_YNC, txt, 0, dialogTextHeight, systemFont, 'l');
  if (ans == ANS_CANCEL)
  {
    return 0;
  }
  else if (ans == ANS_YES)
  {
    return ExistingTextOption_replace_old;
  }
    
  w = MovableModalWindow(-20, -13, -10, -10, "How to Add New Text", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  etdd.pre_app_grp = HiddenGroup (h, 0, 4, ChangePreAppIgnoreChoice);
  SetGroupSpacing (etdd.pre_app_grp, 10, 10);
  RadioButton (etdd.pre_app_grp, "Append");
  RadioButton (etdd.pre_app_grp, "Prefix");
  RadioButton (etdd.pre_app_grp, "Ignore new text");
  if (allow_multi) {
    RadioButton (etdd.pre_app_grp, "Add new qual");
  }
  SetValue (etdd.pre_app_grp, 1);
  SetObjectExtra (etdd.pre_app_grp, &etdd, NULL);
  
  ppt = StaticPrompt (h, "Separate new text and old text with", 
                      0, dialogTextHeight, programFont, 'c');
  etdd.delim_grp = HiddenGroup (h, 0, 5, NULL);
  SetGroupSpacing (etdd.delim_grp, 10, 10);
  RadioButton (etdd.delim_grp, "Semicolon");
  RadioButton (etdd.delim_grp, "Space");
  RadioButton (etdd.delim_grp, "Colon");
  RadioButton (etdd.delim_grp, "Comma");
  RadioButton (etdd.delim_grp, "Do not separate");
  SetValue (etdd.delim_grp, 1);
  
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) etdd.pre_app_grp,
                              (HANDLE) ppt, 
                              (HANDLE) etdd.delim_grp, 
                              (HANDLE) c, 
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
    existing_text = 0;
  }
  else
  {
    handle_choice = GetValue (etdd.pre_app_grp);
    if (handle_choice == 1)
    {
      switch (GetValue (etdd.delim_grp))
      {
        case 1:
          existing_text = ExistingTextOption_append_semi;
          break;
        case 2:
          existing_text = ExistingTextOption_append_space;
          break;
        case 3:
          existing_text = ExistingTextOption_append_colon;
          break;
        case 4:
          existing_text = ExistingTextOption_append_comma;
          break;
        case 5:
          existing_text = ExistingTextOption_append_none;
          break;
      }
    }
    else if (handle_choice == 2)
    {
      switch (GetValue (etdd.delim_grp))
      {
        case 1:
          existing_text = ExistingTextOption_prefix_semi;
          break;
        case 2:
          existing_text = ExistingTextOption_prefix_space;
          break;
        case 3:
          existing_text = ExistingTextOption_prefix_colon;
          break;
        case 4:
          existing_text = ExistingTextOption_prefix_comma;
          break;
        case 5:
          existing_text = ExistingTextOption_prefix_none;
          break;
      }
    }
    else if (handle_choice == 4) 
    {
      existing_text = ExistingTextOption_add_qual;
    }
    else
    {
      existing_text = ExistingTextOption_leave_old;
    }
  }
  Remove (w);
  return existing_text;
}


static Boolean WordSubstitutionIsEmpty (WordSubstitutionPtr word)
{
  ValNodePtr vnp;
  Boolean rval = FALSE;

  if (word == NULL) {
    rval = TRUE;
  } else if (word->word == NULL || word->word[0] == 0) {
    rval = TRUE;
    for (vnp = word->synonyms; vnp != NULL && rval; vnp = vnp->next) {
      if (vnp->data.ptrvalue != NULL || ((CharPtr)vnp->data.ptrvalue)[0] != 0) {
        rval = FALSE;
      }
    }
  } else {
    rval = FALSE;
  }
  return rval;
}


typedef struct editwordsubstitution {
  TexT  word;
  DialoG synonyms;
  ButtoN case_insensitive;
  ButtoN whole_word;

  ButtoN accept_btn;
} EditWordSubstitutionData, PNTR EditWordSubstitutionPtr;


static Boolean EditWordSubstitution (WordSubstitutionPtr word)
{
  ModalAcceptCancelData acd;
  EditWordSubstitutionData    ecd;
  Boolean               rval = FALSE;
  WindoW                w;
  ButtoN                b;
  GrouP                 h, g1, g2, c;
  
  if (word == NULL) return FALSE;

  w = MovableModalWindow(-20, -13, -10, -10, "Word Substitution", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g1 = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g1, "Pattern Word", 0, dialogTextHeight, programFont, 'c');
  ecd.word = DialogText (g1, "", 15, NULL);

  ecd.synonyms = CreateVisibleStringDialog (h, 3, -1, 25);

  g2 = HiddenGroup (h, 2, 0, NULL);
  ecd.case_insensitive = CheckBox (g2, "Case insensitive", NULL);
  ecd.whole_word = CheckBox (g2, "Whole Word", NULL);

  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  ecd.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (ecd.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ecd.word,
                              (HANDLE) g1,
                              (HANDLE) ecd.synonyms,
                              (HANDLE) g2,
                              (HANDLE) c, 
                              NULL);

  if (word->word != NULL) {
    SetTitle (ecd.word, word->word);
  }
  if (word->synonyms != NULL) {
    PointerToDialog (ecd.synonyms, word->synonyms);
  }
  SetStatus (ecd.case_insensitive, !word->case_sensitive);
  SetStatus (ecd.whole_word, word->whole_word);

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
  if (!acd.cancelled)
  {
    word->word = MemFree (word->word);
    word->word = JustSaveStringFromText (ecd.word);
    word->synonyms = ValNodeFreeData (word->synonyms);
    word->synonyms = DialogToPointer (ecd.synonyms);
    word->case_sensitive = !GetStatus (ecd.case_insensitive);
    word->whole_word = GetStatus (ecd.whole_word);

    rval = TRUE;
  }
  Remove (w);
  return rval;
}


typedef struct wordsubstitutionsetdlg {
  DIALOG_MESSAGE_BLOCK
  DoC                     word_doc;
  WordSubstitutionSetPtr  word_list;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} WordSubstitutionSetDlgData, PNTR WordSubstitutionSetDlgPtr;

static void PopulateWordSubstitutionDoc (DoC d, WordSubstitutionSetPtr word_list)
{
  WordSubstitutionPtr word;
  CharPtr    phrase, tmp;
  RecT       r;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };

  if (d == NULL) return;

  Reset (d);

  ObjectRect (d, &r);
  InsetRect (&r, 4, 4);
  
  ColFmt[1].pixWidth = r.right - r.left - 12;

  Reset (d);

  for (word = word_list; word != NULL; word = word->next) {
    phrase = SummarizeWordSubstitution (word);
    if (phrase == NULL) {
      AppendText (d, "\tUnable to summarize word substitution\n", &ParFmt, ColFmt, programFont);
    } else {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (phrase) + 3));
      sprintf (tmp, "\t%s\n", phrase);
      phrase = MemFree (phrase);
      AppendText (d, tmp, &ParFmt, ColFmt, programFont);
      tmp = MemFree (tmp);
    }
  }
  AppendText (d, "(Click here to add new word substitution)", NULL, NULL, programFont);
  UpdateDocument (d, 0, 0);
}


static Int4 WordSubstitutionSetLen (WordSubstitutionSetPtr word)
{
  Int4 num_words = 0;

  while (word != NULL) {
    num_words++;
    word = word->next;
  }
  return num_words;
}


static void DrawWordSubstitutionDocControls (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RecT                rct;
  Int4                width;
  PoinT               pt1, pt2;
  WordSubstitutionSetDlgPtr dlg;

  dlg = (WordSubstitutionSetDlgPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0 && item <= WordSubstitutionSetLen (dlg->word_list)) {
    rct = *r;
  
    /* draw X for deletion */
    width = 10;
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1;
    pt2.x = pt1.x + width;
    pt2.y = pt1.y + width;
    DrawLine (pt1, pt2);
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1 + width;
    pt2.x = pt1.x + width;
    pt2.y = rct.top + 1;
    DrawLine (pt1, pt2);
  }
}


static void ClickWordSubstitutionDoc (DoC d, PoinT pt)
{
  Int2      item, row, col;
  RecT      rct;
  BaR       sb_vert;
  Int4      scroll_pos = 0, scroll_max;
  WordSubstitutionSetDlgPtr f;
  WordSubstitutionPtr word, word_prev = NULL;
  Boolean             changed = FALSE;
  
  f = (WordSubstitutionSetDlgPtr) GetObjectExtra (d);
  if (f == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item == 0 && row == 0 && f->word_list == NULL) {
    /* create new constraint */
    word = WordSubstitutionNew();
    if (EditWordSubstitution (word) && !WordSubstitutionIsEmpty(word)) {
      f->word_list = word;
      changed = TRUE;
    } else {
      word = WordSubstitutionFree (word);
    }
  } else if (item > 0 && row > 0) {
    if (item == WordSubstitutionSetLen(f->word_list) + 1) {
      /* create new constraint */
      word = WordSubstitutionNew();
      if (EditWordSubstitution (word) && !WordSubstitutionIsEmpty(word)) {
        if (f->word_list == NULL) {
          f->word_list = word;
        } else {
          word_prev = f->word_list;
          while (word_prev->next != NULL) {
            word_prev = word_prev->next;
          }
          word_prev->next = word;
        }
        changed = TRUE;
      } else {
        word = WordSubstitutionFree (word);
      }
    } else {
      for (word = f->word_list; word != NULL && item > 1; word = word->next, item--) {
        word_prev = word;
      }
      if (word != NULL) {      
        sb_vert = GetSlateVScrollBar ((SlatE) f->word_doc);
        scroll_pos = GetBarValue (sb_vert);
        switch (col) {
          case 1:
            /* delete this item */
            if (word_prev == NULL) {
              f->word_list = word->next;
            } else {
              word_prev->next = word->next;
            }
            word->next = NULL;
            word = WordSubstitutionFree (word);
            changed = TRUE;
            break;
          case 2:
            /* edit */
            if (EditWordSubstitution (word)) {
              if (WordSubstitutionIsEmpty(word)) {
                if (word_prev == NULL) {
                  f->word_list = word->next;
                } else {
                  word_prev->next = word->next;
                }
                word->next = NULL;
                word = WordSubstitutionFree (word);
              }
              changed = TRUE;
            }
            break;
        }
      }
    }
  }
  if (changed) {
    PopulateWordSubstitutionDoc (f->word_doc, f->word_list);
    if (scroll_pos > 0) {
      sb_vert = GetSlateVScrollBar ((SlatE) f->word_doc);
      scroll_max = GetBarMax (sb_vert);
      if (scroll_pos > scroll_max) {
        scroll_pos = scroll_max;
      }
      CorrectBarValue (sb_vert, scroll_pos);
    }
    if (f->change_notify != NULL) {
      (f->change_notify) (f->change_userdata);
    }
  }
}

static void WordSubstitutionSetToDialog (DialoG d, Pointer data)
{
  WordSubstitutionSetDlgPtr dlg;

  dlg = (WordSubstitutionSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dlg->word_list = WordSubstitutionSetFree (dlg->word_list);
  if (data != NULL) {
    dlg->word_list = AsnIoMemCopy ((WordSubstitutionSetPtr) data,
                                         (AsnReadFunc) WordSubstitutionSetAsnRead,
                                         (AsnWriteFunc) WordSubstitutionSetAsnWrite);
  }
  PopulateWordSubstitutionDoc (dlg->word_doc, dlg->word_list);
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer WordSubstitutionSetFromDialog (DialoG d)
{
  WordSubstitutionSetDlgPtr dlg;
  WordSubstitutionSetPtr word_list = NULL;

  dlg = (WordSubstitutionSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->word_list != NULL) {
    word_list = AsnIoMemCopy ((WordSubstitutionSetPtr) dlg->word_list,
                                    (AsnReadFunc) WordSubstitutionSetAsnRead,
                                    (AsnWriteFunc) WordSubstitutionSetAsnWrite);
  }
  return (Pointer) word_list;
}

static void CleanupWordSubstitutionSetDialog (GraphiC g, VoidPtr data)

{
  WordSubstitutionSetDlgPtr dlg;

  dlg = (WordSubstitutionSetDlgPtr) data;
  if (dlg != NULL) {
    dlg->word_list = WordSubstitutionSetFree(dlg->word_list);
  }
  StdCleanupExtraProc (g, data);
}


static void AddWordSubstitutionBtn (ButtoN b)
{
  WordSubstitutionSetDlgPtr dlg;
  WordSubstitutionPtr       word, prev_word;
  BaR       sb_vert;
  Int4      scroll_pos = 0, scroll_max;

  dlg = (WordSubstitutionSetDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  sb_vert = GetSlateVScrollBar ((SlatE) dlg->word_doc);
  scroll_pos = GetBarValue (sb_vert);

  /* create new WordSubstitution */
  word = WordSubstitutionNew ();
  if (EditWordSubstitution (word) && !WordSubstitutionIsEmpty(word)) {
    if (dlg->word_list == NULL) {
      dlg->word_list = word;
    } else {
      prev_word = dlg->word_list;
      while (prev_word->next != NULL) {
        prev_word = prev_word->next;
      }
      prev_word->next = word;
    }
    PopulateWordSubstitutionDoc (dlg->word_doc, dlg->word_list);
    if (scroll_pos > 0) {
      sb_vert = GetSlateVScrollBar ((SlatE) dlg->word_doc);
      scroll_max = GetBarMax (sb_vert);
      if (scroll_pos > scroll_max) {
        scroll_pos = scroll_max;
      }
      CorrectBarValue (sb_vert, scroll_pos);
    }
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  } else {
    word = WordSubstitutionFree (word);
  }

}


NLM_EXTERN DialoG WordSubstitutionSetDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  WordSubstitutionSetDlgPtr dlg;
  GrouP               p, g;
  ButtoN              b;

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  dlg = (WordSubstitutionSetDlgPtr) MemNew (sizeof ( WordSubstitutionSetDlgData));
  SetObjectExtra (p, dlg, CleanupWordSubstitutionSetDialog);

  dlg->dialog = (DialoG) p;
  dlg->todialog = WordSubstitutionSetToDialog;
  dlg->fromdialog = WordSubstitutionSetFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  dlg->word_list = NULL;
  dlg->word_doc = DocumentPanel (p, stdCharWidth * 30, stdLineHeight * 3);
  SetObjectExtra (dlg->word_doc, dlg, NULL);
  SetDocProcs (dlg->word_doc, ClickWordSubstitutionDoc, NULL, NULL, NULL);
  SetDocShade (dlg->word_doc, DrawWordSubstitutionDocControls, NULL, NULL, NULL);
  PopulateWordSubstitutionDoc (dlg->word_doc, dlg->word_list);

  g = HiddenGroup (p, 2, 0, NULL);
  b = PushButton (g, "Add Word Substitution", AddWordSubstitutionBtn);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (g, "Clear Word Substitutions", ClearDialogBtn);
  SetObjectExtra (b, p, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->word_doc, (HANDLE) g, NULL);

  return (DialoG) p;
}


typedef struct stringconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  DialoG tbs;
  GrouP  tab_pages[2];
  PopuP  match_choice;
  TexT   match_text;
  ButtoN insensitive;
  ButtoN whole_word;  
  ButtoN ignore_space;
  ButtoN ignore_punct;
  ButtoN ignore_weasel;
  DialoG ignore_words;
  ButtoN is_all_caps;
  ButtoN is_all_lower;
  ButtoN is_all_punct;

  Int2   current_page;
  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} StringConstraintDialogData, PNTR StringConstraintDialogPtr;

static void ResetStringConstraintDialog (StringConstraintDialogPtr scdp)
{
  if (scdp == NULL) return;
  
  SetValue (scdp->match_choice, 1);
  SetTitle (scdp->match_text, "");
  SetStatus (scdp->insensitive, FALSE);
  SetStatus (scdp->whole_word, FALSE);
  SetStatus (scdp->ignore_space, FALSE);
  SetStatus (scdp->ignore_punct, FALSE);
  SetStatus (scdp->ignore_weasel, FALSE);
  SetStatus (scdp->is_all_caps, FALSE);
  SetStatus (scdp->is_all_lower, FALSE);
  SetStatus (scdp->is_all_punct, FALSE);
  PointerToDialog (scdp->ignore_words, NULL);
}


static void ClearStringConstraintDialogText (DialoG d)
{
  StringConstraintDialogPtr dlg;

  dlg = (StringConstraintDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    SetTitle (dlg->match_text, "");
  }
  PointerToDialog (dlg->ignore_words, NULL);
}


static Int4 GetPopupPosForStringConstraint (StringConstraintPtr scp)
{
  Int4 rval = 1;

  if (scp == NULL) return 1;

  switch (scp->match_location)
  {
    case String_location_contains:
      if (scp->not_present) 
      {
        rval = 2;
      } else {
        rval = 1;
      }
      break;
    case String_location_equals:
      if (scp->not_present) {
        rval = 4;
      } else {
        rval = 3;
      }
      break;
    case String_location_starts:
      if (scp->not_present) {
        rval = 9;
      } else {
        rval = 5;
      }
      break;
    case String_location_ends:
      if (scp->not_present) {
        rval = 10;
      } else {
        rval = 6;
      }
      break;
    case String_location_inlist:
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
  StringConstraintDialogPtr scdp;
  StringConstraintPtr       scp;

  scdp = (StringConstraintDialogPtr) GetObjectExtra (d);
  scp = (StringConstraintPtr) data;
  if (scdp == NULL)
  {
    return;
  }

  if (scp == NULL)
  {
    ResetStringConstraintDialog (scdp);
  }
  else
  {
    SetValue (scdp->match_choice, GetPopupPosForStringConstraint (scp));
    
    if (scp->match_text == NULL)
    {
      SetTitle (scdp->match_text, "");
    }
    else
    {
      SetTitle (scdp->match_text, scp->match_text);
    }
    
    SetStatus (scdp->insensitive, !scp->case_sensitive);
    SetStatus (scdp->whole_word, scp->whole_word);
    SetStatus (scdp->ignore_space, scp->ignore_space);
    SetStatus (scdp->ignore_punct, scp->ignore_punct);
    SetStatus (scdp->ignore_weasel, scp->ignore_weasel);
    SetStatus (scdp->is_all_caps, scp->is_all_caps);
    SetStatus (scdp->is_all_lower, scp->is_all_lower);
    SetStatus (scdp->is_all_punct, scp->is_all_punct);
    PointerToDialog (scdp->ignore_words, scp->ignore_words);
  }
  if (scdp->change_notify != NULL) {
    (scdp->change_notify) (scdp->change_userdata);
  }
}

static Pointer DialogToStringConstraint (DialoG d)

{
  StringConstraintDialogPtr scdp;
  StringConstraintPtr       scp;
  Int4                      match_choice;

  scdp = (StringConstraintDialogPtr) GetObjectExtra (d);
  if (scdp == NULL)
  {
    return NULL;
  }
  scp = StringConstraintNew();
  if (scp != NULL)
  {
    scp->match_text = JustSaveStringFromText (scdp->match_text);
    scp->case_sensitive = !GetStatus (scdp->insensitive);
    scp->whole_word = GetStatus (scdp->whole_word);
    scp->ignore_space = GetStatus (scdp->ignore_space);
    scp->ignore_punct = GetStatus (scdp->ignore_punct);
    scp->is_all_caps = GetStatus (scdp->is_all_caps);
    scp->is_all_lower = GetStatus (scdp->is_all_lower);
    scp->is_all_punct = GetStatus (scdp->is_all_punct);
    scp->ignore_weasel = GetStatus (scdp->ignore_weasel);
    scp->ignore_words = DialogToPointer (scdp->ignore_words);
    match_choice = GetValue (scdp->match_choice);
    switch (match_choice)
    {
      case 1:
        scp->match_location = String_location_contains;
        scp->not_present = FALSE;
        break;
      case 2:
        scp->match_location = String_location_contains;
        scp->not_present = TRUE;
        break;
      case 3:
        scp->match_location = String_location_equals;
        scp->not_present = FALSE;
        break;
      case 4:
        scp->match_location = String_location_equals;
        scp->not_present = TRUE;
        break;
      case 5:
        scp->match_location = String_location_starts;
        scp->not_present = FALSE;
        break;
      case 6:
        scp->match_location = String_location_ends;
        scp->not_present = FALSE;
        break;
      case 7:
        scp->match_location = String_location_inlist;
        scp->not_present = FALSE;
        break;
      case 8:
        scp->match_location = String_location_inlist;
        scp->not_present = TRUE;
        break;
      case 9:
        scp->match_location = String_location_starts;
        scp->not_present = TRUE;
        break;
      case 10:
        scp->match_location = String_location_ends;
        scp->not_present = TRUE;
        break;
      default:
        scp->match_location = String_location_contains;
        scp->not_present = FALSE;
        break;
    }
  }
  return scp;
}

static void StringConstraintMessage (DialoG d, Int2 mssg)

{
  StringConstraintDialogPtr scdp;

  scdp = (StringConstraintDialogPtr) GetObjectExtra (d);
  if (scdp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        ResetStringConstraintDialog (scdp);        
        break;
      case VIB_MSG_ENTER :
        Select (scdp->match_text);
        break;
      default :
        break;
    }
  }
}


static void ChangeStringConstraintDialogText (TexT t)
{
  StringConstraintDialogPtr dlg;

  dlg = (StringConstraintDialogPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeStringConstraintDialogBtn (ButtoN b)
{
  StringConstraintDialogPtr dlg;

  dlg = (StringConstraintDialogPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  if (b == dlg->is_all_caps) {
    if (GetStatus (dlg->is_all_caps)) {
      SetStatus (dlg->is_all_lower, FALSE);
      SetStatus (dlg->is_all_punct, FALSE);
    }
  }
  if (b == dlg->is_all_lower) {
    if (GetStatus (dlg->is_all_lower)) {
      SetStatus (dlg->is_all_caps, FALSE);
      SetStatus (dlg->is_all_punct, FALSE);
    }
  }
  if (b == dlg->is_all_punct) {
    if (GetStatus (dlg->is_all_punct)) {
      SetStatus (dlg->is_all_caps, FALSE);
      SetStatus (dlg->is_all_lower, FALSE);
    }
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static ValNodePtr TestStringConstraintDialog (DialoG d)

{
  ValNodePtr err_list = NULL;
  StringConstraintDialogPtr dlg;

  dlg = (StringConstraintDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
  }
  return err_list;
}


static CharPtr  sStringConstraintTabs [] = {
  "String", "Word Substitution List", NULL
};


static void ChangeStringConstraintPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  StringConstraintDialogPtr  dlg;

  dlg = (StringConstraintDialogPtr) data;
  if (dlg == NULL) {
    return;
  }
  dlg->current_page = newval;
  switch (dlg->current_page) {
    case 0:
      Show (dlg->tab_pages[0]);
      Hide (dlg->tab_pages[1]);
      break;
    case 1:
      Show (dlg->tab_pages[1]);
      Hide (dlg->tab_pages[0]);
      break;
  }

    Update ();
}


NLM_EXTERN DialoG StringConstraintDialog (GrouP h, CharPtr label, Boolean clear_btn, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  StringConstraintDialogPtr scdp;
  GrouP                     p, tab_grp, g, k, special_grp, special_grp2;
  ButtoN                    b = NULL;
  
  scdp = (StringConstraintDialogPtr) MemNew (sizeof (StringConstraintDialogData));
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
  scdp->testdialog = TestStringConstraintDialog;
  scdp->change_notify = change_notify;
  scdp->change_userdata = change_userdata;

  scdp->tbs = CreateFolderTabs (p, sStringConstraintTabs, 0, 
                                  0, 0, SYSTEM_FOLDER_TAB,
                                  ChangeStringConstraintPage, (Pointer) scdp);
  scdp->current_page = 0;
  tab_grp = HiddenGroup (p, 0, 0, NULL);
  scdp->tab_pages[0] = HiddenGroup (tab_grp, -1, 0, NULL);
  g = HiddenGroup (scdp->tab_pages[0], 3, 0, NULL);
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
  scdp->match_text = DialogText (g, "", 15, ChangeStringConstraintDialogText);
  SetObjectExtra (scdp->match_text, scdp, NULL);
  
  k = HiddenGroup (scdp->tab_pages[0], 4, 0, NULL);
  SetGroupSpacing (k, 10, 10);
  scdp->insensitive = CheckBox (k, "Case Insensitive", NULL);
  scdp->whole_word = CheckBox (k, "Whole Word", NULL);
  scdp->ignore_space = CheckBox (k, "Ignore Space", NULL);
  scdp->ignore_punct = CheckBox (k, "Ignore Punctuation", NULL);

  special_grp2 = HiddenGroup (scdp->tab_pages[0], 3, 0, NULL);
  SetGroupSpacing (special_grp2, 10, 10);
  scdp->ignore_weasel = CheckBox (special_grp2, "Ignore 'putative' synonyms", NULL);

  special_grp = HiddenGroup (scdp->tab_pages[0], 3, 0, NULL);
  SetGroupSpacing (special_grp, 10, 10);
  scdp->is_all_caps = CheckBox (special_grp, "All letters are uppercase", ChangeStringConstraintDialogBtn);
  SetObjectExtra (scdp->is_all_caps, scdp, NULL);
  scdp->is_all_lower = CheckBox (special_grp, "All letters are lowercase", ChangeStringConstraintDialogBtn);
  SetObjectExtra (scdp->is_all_lower, scdp, NULL);
  scdp->is_all_punct = CheckBox (special_grp, "All characters are punctuation", ChangeStringConstraintDialogBtn);
  SetObjectExtra (scdp->is_all_punct, scdp, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) special_grp2, (HANDLE) special_grp, NULL);


  scdp->tab_pages[1] = HiddenGroup (tab_grp, -1, 0, NULL);
  scdp->ignore_words = WordSubstitutionSetDialog (scdp->tab_pages[1], change_notify, change_userdata);
  Hide (scdp->tab_pages[1]);
  AlignObjects (ALIGN_CENTER, (HANDLE) scdp->tab_pages[0], scdp->tab_pages[1], NULL);


  if (clear_btn)
  {
    b = PushButton (p, "Clear Constraint", ClearDialogBtn);
    SetObjectExtra (b, p, NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) tab_grp, (HANDLE) b, NULL);
    
  return (DialoG) p;
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


typedef struct enddistancedialog
{
  DIALOG_MESSAGE_BLOCK
  PopuP dist_type;
  TexT  dist_val;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata; 
} EndDistanceDialogData, PNTR EndDistanceDialogPtr;


static void ChangeEndDistanceDialogPopup (PopuP p)
{
  EndDistanceDialogPtr dlg;
  Int2 val;

  dlg = (EndDistanceDialogPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }

  val = GetValue (dlg->dist_type);
  if (val == 2 || val == 3 || val == 4) {
    Show (dlg->dist_val);
  } else {
    Hide (dlg->dist_val);
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeEndDistanceDialogTexT (TexT t)
{
  EndDistanceDialogPtr dlg;

  dlg = (EndDistanceDialogPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void EndDistanceToDialog (DialoG d, Pointer data)
{
  EndDistanceDialogPtr dlg;
  ValNodePtr v;
  Char       buf[15];

  dlg = (EndDistanceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  v = (ValNodePtr) data;
  if (v == NULL) {
    SetValue (dlg->dist_type, 1);
    SetTitle (dlg->dist_val, "");
    Hide (dlg->dist_val);
  } else {
    switch (v->choice) {
      case LocationPosConstraint_dist_from_end:
        SetValue (dlg->dist_type, 2);
        sprintf (buf, "%d", v->data.intvalue);
        SetTitle (dlg->dist_val, buf);
        Show (dlg->dist_val);
        break;
      case LocationPosConstraint_max_dist_from_end:
        SetValue (dlg->dist_type, 3);
        sprintf (buf, "%d", v->data.intvalue);
        SetTitle (dlg->dist_val, buf);
        Show (dlg->dist_val);
        break;
      case LocationPosConstraint_min_dist_from_end:
        SetValue (dlg->dist_type, 4);
        sprintf (buf, "%d", v->data.intvalue);
        SetTitle (dlg->dist_val, buf);
        Show (dlg->dist_val);
        break;
      default:
        SetValue (dlg->dist_type, 1);
        SetTitle (dlg->dist_val, "");
        Hide (dlg->dist_val);
        break;
    }
  }
}


static Pointer DialogToEndDistance (DialoG d)
{
  EndDistanceDialogPtr dlg;
  Int2       val;
  ValNodePtr v = NULL;
  CharPtr    str;

  dlg = (EndDistanceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->dist_type);
  if ((val == 2 || val == 3 || val == 4)
       && !TextHasNoText(dlg->dist_val))
  {
    str = SaveStringFromText (dlg->dist_val);
    if (IsAllDigits (str)) {
      v = ValNodeNew (NULL);
      switch (val) {
        case 2:
          v->choice = LocationPosConstraint_dist_from_end;
          break;
        case 3:
          v->choice = LocationPosConstraint_max_dist_from_end;
          break;
        case 4:
          v->choice = LocationPosConstraint_min_dist_from_end;
          break;
      }
      str = SaveStringFromText (dlg->dist_val);
      v->data.intvalue = atoi (str);
    }
    str = MemFree (str);
  }
  return v;
}


static ValNodePtr TestEndDistanceDialog (DialoG d)
{
  EndDistanceDialogPtr dlg;
  Int2       val;
  CharPtr    str;
  ValNodePtr err_list = NULL;

  dlg = (EndDistanceDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    val = GetValue (dlg->dist_type);
    if (val == 2 || val == 3 || val == 4) {
      if (TextHasNoText (dlg->dist_val)) {
        ValNodeAddPointer (&err_list, 0, "No distance value");
      } else {
        str = SaveStringFromText (dlg->dist_val);
        if (!IsAllDigits (str)) {
          ValNodeAddPointer (&err_list, 0, "Non-numeric value");
        }
      }
    }
  }

  return err_list;
}


NLM_EXTERN DialoG EndDistanceDialog (GrouP h, CharPtr title, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  EndDistanceDialogPtr dlg;
  GrouP                p;
  
  dlg = (EndDistanceDialogPtr) MemNew (sizeof (EndDistanceDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 4, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = EndDistanceToDialog;
  dlg->fromdialog = DialogToEndDistance;
  dlg->testdialog = TestEndDistanceDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  StaticPrompt (p, "Where end of sequence is ", 0, dialogTextHeight, systemFont, 'l');
  dlg->dist_type = PopupList (p, TRUE, ChangeEndDistanceDialogPopup);
  SetObjectExtra (dlg->dist_type, dlg, NULL);
  PopupItem (dlg->dist_type, "Any distance");
  PopupItem (dlg->dist_type, "Exactly");
  PopupItem (dlg->dist_type, "No more than");
  PopupItem (dlg->dist_type, "No less than");
  SetValue (dlg->dist_type, 1);
  dlg->dist_val = DialogText (p, "", 5, ChangeEndDistanceDialogTexT);
  SetObjectExtra (dlg->dist_val, dlg, NULL);
  Hide (dlg->dist_val);
  StaticPrompt (p, title, 0, dialogTextHeight, systemFont, 'l');

  return (DialoG) p;
}


typedef struct locationconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  PopuP  strand;
  PopuP  sequence_type;
  PopuP  partial5;
  PopuP  partial3;
  PopuP  location_type;
  DialoG dist5;
  DialoG dist3;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata; 
} LocationConstraintDialogData, PNTR LocationConstraintDialogPtr;


static void ChangeLocationConstraintDialogPopup (PopuP p)
{
  LocationConstraintDialogPtr dlg;

  dlg = (LocationConstraintDialogPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void LocationConstraintToDialog (DialoG d, Pointer data)
{
  LocationConstraintDialogPtr dlg;
  LocationConstraintPtr l;

  dlg = (LocationConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  l = (LocationConstraintPtr) data;
  if (l == NULL) {
    SetValue (dlg->strand, 1);
    SetValue (dlg->sequence_type, 1);
    SetValue (dlg->partial5, 1);
    SetValue (dlg->partial3, 1);
    SetValue (dlg->location_type, 1);
    PointerToDialog (dlg->dist5, NULL);
    PointerToDialog (dlg->dist3, NULL);
  } else {
    switch (l->strand) {
      case Strand_constraint_any:
        SetValue (dlg->strand, 1);
        break;
      case Strand_constraint_plus:
        SetValue (dlg->strand, 2);
        break;
      case Strand_constraint_minus:
        SetValue (dlg->strand, 3);
        break;
      default:
        SetValue (dlg->strand, 1);
        break;
    }
    switch (l->seq_type) {
      case Seqtype_constraint_any:
        SetValue (dlg->sequence_type, 1);
        break;
      case Seqtype_constraint_nuc:
        SetValue (dlg->sequence_type, 2);
        break;
      case Seqtype_constraint_prot:
        SetValue (dlg->sequence_type, 3);
        break;
      default:
        SetValue (dlg->sequence_type, 1);
        break;
     }
    switch (l->partial5) {
      case Partial_constraint_either:
        SetValue (dlg->partial5, 1);
        break;
      case Partial_constraint_partial:
        SetValue (dlg->partial5, 2);
        break;
      case Partial_constraint_complete:
        SetValue (dlg->partial5, 3);
        break;
    }
    switch (l->partial3) {
      case Partial_constraint_either:
        SetValue (dlg->partial3, 1);
        break;
      case Partial_constraint_partial:
        SetValue (dlg->partial3, 2);
        break;
      case Partial_constraint_complete:
        SetValue (dlg->partial3, 3);
        break;
    }
    switch (l->location_type) {
      case Location_type_constraint_any:
        SetValue (dlg->location_type, 1);
        break;
      case Location_type_constraint_single_interval:
        SetValue (dlg->location_type, 2);
        break;
      case Location_type_constraint_joined:
        SetValue (dlg->location_type, 3);
        break;
      case Location_type_constraint_ordered:
        SetValue (dlg->location_type, 4);
        break;
    }
    PointerToDialog (dlg->dist5, l->end5);
    PointerToDialog (dlg->dist3, l->end3);
  }
  ChangeLocationConstraintDialogPopup (dlg->strand);
}


static Pointer DialogToLocationConstraint (DialoG d)
{
  LocationConstraintDialogPtr dlg;
  LocationConstraintPtr l;

  dlg = (LocationConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  l = LocationConstraintNew();
  if (l != NULL) {
    switch (GetValue (dlg->strand)) {
      case 1:
        l->strand = Strand_constraint_any;
        break;
      case 2:
        l->strand = Strand_constraint_plus;
        break;
      case 3:
        l->strand = Strand_constraint_minus;
        break;
      default:
        l->strand = Strand_constraint_any;
        break;
    }
    switch (GetValue (dlg->sequence_type)) {
      case 1:
        l->seq_type = Seqtype_constraint_any;
        break;
      case 2:
        l->seq_type = Seqtype_constraint_nuc;
        break;
      case 3:
        l->seq_type = Seqtype_constraint_prot;
        break;
      default:
        l->seq_type = Seqtype_constraint_any;
        break;
    }
    switch (GetValue (dlg->partial5)) {
      case 1:
        l->partial5 = Partial_constraint_either;
        break;
      case 2:
        l->partial5 = Partial_constraint_partial;
        break;
      case 3:
        l->partial5 = Partial_constraint_complete;
        break;
    }
    switch (GetValue (dlg->partial3)) {
      case 1:
        l->partial3 = Partial_constraint_either;
        break;
      case 2:
        l->partial3 = Partial_constraint_partial;
        break;
      case 3:
        l->partial3 = Partial_constraint_complete;
        break;
    }
    switch (GetValue (dlg->location_type)) {
      case 1:
        l->location_type = Location_type_constraint_any;
        break;
      case 2:
        l->location_type = Location_type_constraint_single_interval;
        break;
      case 3:
        l->location_type = Location_type_constraint_joined;
        break;
      case 4:
        l->location_type = Location_type_constraint_ordered;
        break;
    }
    l->end5 = DialogToPointer (dlg->dist5);
    l->end3 = DialogToPointer (dlg->dist3);
  }
  if (IsLocationConstraintEmpty (l)) {
    l = LocationConstraintFree (l);
  }
  return l;
}


static ValNodePtr TestLocationConstraintDialog (DialoG d)
{
  LocationConstraintDialogPtr dlg;
  ValNodePtr err_list = NULL;
  LocationConstraintPtr lcp;

  dlg = (LocationConstraintDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    lcp = DialogToPointer (d);
    if (IsLocationConstraintEmpty (lcp)) {
      ValNodeAddPointer (&err_list, 0, "location constraint empty");
    }
    lcp = LocationConstraintFree (lcp);
    ValNodeLink (&err_list, TestDialog (dlg->dist5));
    ValNodeLink (&err_list, TestDialog (dlg->dist3));
  }
  return err_list;
}


NLM_EXTERN DialoG LocationConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  LocationConstraintDialogPtr dlg;
  GrouP                     p, g1, g2;
  
  dlg = (LocationConstraintDialogPtr) MemNew (sizeof (LocationConstraintDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = LocationConstraintToDialog;
  dlg->fromdialog = DialogToLocationConstraint;
  dlg->testdialog = TestLocationConstraintDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g1 = HiddenGroup (p, 4, 0, NULL);
  SetGroupSpacing (g1, 10, 10);  
  StaticPrompt (g1, "on", 0, dialogTextHeight, systemFont, 'l');
  dlg->strand = PopupList (g1, TRUE, ChangeLocationConstraintDialogPopup);
  SetObjectExtra (dlg->strand, dlg, NULL);
  PopupItem (dlg->strand, "Any strand");
  PopupItem (dlg->strand, "Plus strand");
  PopupItem (dlg->strand, "Minus strand");
  SetValue (dlg->strand, 1);
  
  StaticPrompt (g1, "on", 0, dialogTextHeight, systemFont, 'l');
  dlg->sequence_type = PopupList (g1, TRUE, ChangeLocationConstraintDialogPopup);
  SetObjectExtra (dlg->sequence_type, dlg, NULL);
  PopupItem (dlg->sequence_type, "nucleotide and protein sequences");
  PopupItem (dlg->sequence_type, "nucleotide sequences only");
  PopupItem (dlg->sequence_type, "protein sequences only");
  SetValue (dlg->sequence_type, 1);

  StaticPrompt (g1, "5' end", 0, dialogTextHeight, systemFont, 'l');
  dlg->partial5 = PopupList (g1, TRUE, ChangeLocationConstraintDialogPopup);
  SetObjectExtra (dlg->partial5, dlg, NULL);
  PopupItem (dlg->partial5, "Partial or Complete");
  PopupItem (dlg->partial5, "Partial");
  PopupItem (dlg->partial5, "Complete");
  SetValue (dlg->partial5, 1);

  StaticPrompt (g1, "3' end", 0, dialogTextHeight, systemFont, 'l');
  dlg->partial3 = PopupList (g1, TRUE, ChangeLocationConstraintDialogPopup);
  SetObjectExtra (dlg->partial3, dlg, NULL);
  PopupItem (dlg->partial3, "Partial or Complete");
  PopupItem (dlg->partial3, "Partial");
  PopupItem (dlg->partial3, "Complete");
  SetValue (dlg->partial3, 1);

  StaticPrompt (g1, "location type", 0, dialogTextHeight, systemFont, 'l');
  dlg->location_type = PopupList (g1, TRUE, ChangeLocationConstraintDialogPopup);
  SetObjectExtra (dlg->location_type, dlg, NULL);
  PopupItem (dlg->location_type, "Any");
  PopupItem (dlg->location_type, "Single Interval");
  PopupItem (dlg->location_type, "Joined");
  PopupItem (dlg->location_type, "Ordered");
  SetValue (dlg->location_type, 1);

  g2 = HiddenGroup (p, 0, 2, NULL);
  SetGroupSpacing (g2, 10, 10);
  dlg->dist5 = EndDistanceDialog (g2, "from 5' end of feature", dlg->change_notify, dlg->change_userdata);
  dlg->dist3 = EndDistanceDialog (g2, "from 3' end of feature", dlg->change_notify, dlg->change_userdata);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);

  return (DialoG) p;
}


typedef struct sourcequalchoicedialog
{
  DIALOG_MESSAGE_BLOCK
  DialoG dlg;
  DialoG tax_dlg;
  DialoG loc_dlg;
  DialoG orig_dlg;
  GrouP  choice_type;
  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} SourceQualChoiceDialogData, PNTR SourceQualChoiceDialogPtr;


static void ChangeSourceQualChoice (GrouP g)
{
  SourceQualChoiceDialogPtr dlg;
  Int4 val;

  dlg = (SourceQualChoiceDialogPtr) GetObjectExtra (g);
  if (dlg == NULL || dlg->choice_type == NULL) return;

  val = GetValue (dlg->choice_type);
  switch (val) {
    case 1:
      Show (dlg->dlg);
      Hide (dlg->tax_dlg);
      Hide (dlg->loc_dlg);
      Hide (dlg->orig_dlg);
      break;
    case 2:
      Hide (dlg->dlg);
      Show (dlg->tax_dlg);
      Hide (dlg->loc_dlg);
      Hide (dlg->orig_dlg);
      break;
    case 3:
      Hide (dlg->dlg);
      Hide (dlg->tax_dlg);
      Show (dlg->loc_dlg);
      Hide (dlg->orig_dlg);
      break;
    case 4:
      Hide (dlg->dlg);
      Hide (dlg->tax_dlg);
      Hide (dlg->loc_dlg);
      Show (dlg->orig_dlg);
      break;
    default:
      Hide (dlg->dlg);
      Hide (dlg->tax_dlg);
      Hide (dlg->loc_dlg);
      Hide (dlg->orig_dlg);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Boolean IsSourceQualTaxQual (Int4 srcqual)
{
  if (srcqual == Source_qual_taxname
      || srcqual == Source_qual_lineage
      || srcqual == Source_qual_common_name
      || srcqual == Source_qual_division) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void SourceQualChoiceToDialog (DialoG d, Pointer data)
{
  SourceQualChoiceDialogPtr dlg;
  ValNodePtr                vnp;
  ValNode                   vn;

  dlg = (SourceQualChoiceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  vnp = (ValNodePtr) data;

  if (vnp == NULL) {
    SafeSetValue (dlg->choice_type, 1);
    PointerToDialog (dlg->dlg, NULL);
  } else if (vnp->choice == SourceQualChoice_textqual) {
    vn.choice = 0;
    vn.next = NULL;
    vn.data.ptrvalue = GetSourceQualName (vnp->data.intvalue);
    if (IsSourceQualTaxQual (vnp->data.intvalue)) {
      SafeSetValue (dlg->choice_type, 2);
      PointerToDialog (dlg->tax_dlg, &vn);
    } else {
      SafeSetValue (dlg->choice_type, 1);
      PointerToDialog (dlg->dlg, &vn);
    }
  } else if (vnp->choice == SourceQualChoice_location) {
    SafeSetValue (dlg->choice_type, 3);
    vn.choice = vnp->data.intvalue;
    vn.next = NULL;
    vn.data.ptrvalue = NULL;
    PointerToDialog (dlg->loc_dlg, &vn);
  } else if (vnp->choice == SourceQualChoice_origin) {
    SafeSetValue (dlg->choice_type, 4);
    vn.choice = vnp->data.intvalue;
    vn.next = NULL;
    vn.data.ptrvalue = NULL;
    PointerToDialog (dlg->orig_dlg, &vn);
  }

  ChangeSourceQualChoice (dlg->choice_type);
}


static Pointer DialogToSourceQualChoice (DialoG d)
{
  SourceQualChoiceDialogPtr dlg;
  ValNodePtr                vnp, sqc = NULL;
  Int4                      choice_type = SourceQualChoice_textqual;
  DialoG                    text_dlg = NULL;

  dlg = (SourceQualChoiceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->choice_type == NULL) {
    choice_type = SourceQualChoice_textqual;
    text_dlg = dlg->dlg;
  } else {
    switch (GetValue (dlg->choice_type)) {
      case 1:
        choice_type = SourceQualChoice_textqual;
        text_dlg = dlg->dlg;
        break;
      case 2:
        choice_type = SourceQualChoice_textqual;
        text_dlg = dlg->tax_dlg;
        break;
      case 3:
        choice_type = SourceQualChoice_location;
        break;
      case 4:
        choice_type = SourceQualChoice_origin;
        break;
    }
  }
  switch (choice_type) {
    case SourceQualChoice_textqual:
      vnp = DialogToPointer (text_dlg);
      if (vnp != NULL) { 
        sqc = ValNodeNew (NULL);
        sqc->choice = SourceQualChoice_textqual;
        sqc->data.intvalue = GetSourceQualTypeByName (vnp->data.ptrvalue);
      }
      vnp = ValNodeFree (vnp);
      break;
    case SourceQualChoice_location:
      vnp = DialogToPointer (dlg->loc_dlg);
      sqc = ValNodeNew (NULL);
      sqc->choice = SourceQualChoice_location;
      if (vnp == NULL) {
        sqc->data.intvalue = 0;
      } else {
        sqc->data.intvalue = vnp->choice;
      }
      vnp = ValNodeFree (vnp);
      break;
    case SourceQualChoice_origin:
      vnp = DialogToPointer (dlg->orig_dlg);
      sqc = ValNodeNew (NULL);
      sqc->choice = SourceQualChoice_origin;
      if (vnp == NULL) { 
        sqc->data.intvalue = 0;
      } else {
        sqc->data.intvalue = vnp->choice;
      }
      vnp = ValNodeFree (vnp);
      break;
  }
  return sqc;
}


static ValNodePtr TestSourceQualChoiceDialog (DialoG d)
{
  SourceQualChoiceDialogPtr dlg;
  ValNodePtr err_list = NULL;
  Int4 val;

  dlg = (SourceQualChoiceDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->choice_type == NULL) {
    err_list = TestDialog (dlg->dlg);
  } else {
    val = GetValue (dlg->choice_type);
    switch (val) {
      case 1:
        err_list = TestDialog (dlg->dlg);
        break;
      case 2:
        err_list = TestDialog (dlg->tax_dlg);
        break;
      case 3:
        err_list = TestDialog (dlg->loc_dlg);
        break;
      case 4:
        err_list = TestDialog (dlg->orig_dlg);
        break;
      default:
        ValNodeAddPointer (&err_list, 0, "no source qual choice");
        break;
    }
  }
  return err_list;
}


static void SourceQualChoiceDialogMessage (DialoG d, Int2 mssg)

{
  SourceQualChoiceDialogPtr dlg;
  ValNode vn;

  dlg = (SourceQualChoiceDialogPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_ENTER :
        switch (GetValue (dlg->choice_type)) {
          case 1:
            MemSet (&vn, 0, sizeof (ValNode));
            vn.choice = Source_qual_acronym;
            vn.data.ptrvalue = "acronym";
            PointerToDialog (dlg->dlg, &vn);
            SendMessageToDialog (dlg->dlg, VIB_MSG_ENTER);
            break;
          case 2:
            Select (dlg->tax_dlg);
            break;
          case 3:
            Select (dlg->loc_dlg);
            break;
          case 4:
            Select (dlg->orig_dlg);
            break;
        }
      break;
    }
  }
}


static DialoG SourceQualChoiceDialog (GrouP h, Boolean text_only, Boolean for_remove, Boolean tall, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  SourceQualChoiceDialogPtr dlg;
  GrouP                     p, g;
  ValNodePtr                qual_list, loc_list, orig_list, tax_list = NULL;
  ButtoN                    b1, b2, b3 = NULL, b4 = NULL;
  ValNode                   vn;
  Int2                      list_height;
  
  dlg = (SourceQualChoiceDialogPtr) MemNew (sizeof (SourceQualChoiceDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SourceQualChoiceToDialog;
  dlg->fromdialog = DialogToSourceQualChoice;
  dlg->testdialog = TestSourceQualChoiceDialog;
  dlg->dialogmessage = SourceQualChoiceDialogMessage;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  if (tall) {
    list_height = TALL_SELECTION_LIST;
  } else {
    list_height = SHORT_SELECTION_LIST;
  }

  dlg->choice_type = HiddenGroup (p, 0, 4, ChangeSourceQualChoice);
  SetObjectExtra (dlg->choice_type, dlg, NULL);
  b1 = RadioButton (dlg->choice_type, "Text Qualifier");
  b2 = RadioButton (dlg->choice_type, "Taxonomy");
  if (!text_only) {
    b3 = RadioButton (dlg->choice_type, "Location");
    b4 = RadioButton (dlg->choice_type, "Origin");
  }
  
  g = HiddenGroup (p, 0, 0, NULL);
  qual_list = GetSourceQualList (for_remove);
  dlg->dlg = ValNodeSelectionDialog (g, qual_list, list_height,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "source qual",
                                               change_notify, change_userdata, FALSE);

  ValNodeAddPointer (&tax_list, 0, StringSave (GetSourceQualName(Source_qual_taxname)));
  ValNodeAddPointer (&tax_list, 0, StringSave (GetSourceQualName(Source_qual_common_name)));
  ValNodeAddPointer (&tax_list, 0, StringSave (GetSourceQualName(Source_qual_division)));
  ValNodeAddPointer (&tax_list, 0, StringSave (GetSourceQualName(Source_qual_lineage)));
  dlg->tax_dlg = ValNodeSelectionDialog (g, tax_list, list_height,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "source qual",
                                               change_notify, change_userdata, FALSE);
  vn.choice = 0;
  vn.data.ptrvalue = GetSourceQualName(Source_qual_taxname);
  vn.next = NULL;
  PointerToDialog (dlg->tax_dlg, &vn);

  if (!text_only) {
    if (for_remove) {
      loc_list = GetLocationList (for_remove);
    } else {
      loc_list = GetLocListForBioSource (NULL);
    }
    dlg->loc_dlg = ValNodeSelectionDialog (g, loc_list, list_height, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
    if (for_remove) {
      vn.choice = Source_location_unknown;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->loc_dlg, &vn);
    }
    orig_list = GetOriginList (for_remove);
    dlg->orig_dlg = ValNodeSelectionDialog (g, orig_list, list_height, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "origin", 
                                           change_notify, change_userdata, FALSE);
    if (for_remove) {
      vn.choice = Source_origin_unknown;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->orig_dlg, &vn);
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->dlg, (HANDLE) dlg->tax_dlg, (HANDLE) dlg->loc_dlg, (HANDLE) dlg->orig_dlg, NULL);

  SetValue (dlg->choice_type, 1);
  ChangeSourceQualChoice (dlg->choice_type);

  return (DialoG) p;
}


typedef struct sourceconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  GrouP  constraint_type;

  DialoG qual_present;
  GrouP sc_group;
  DialoG qual_string;
  DialoG string_constraint;
  GrouP match_group;
  DialoG qualmatch1;
  DialoG qualmatch2;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} SourceConstraintDialogData, PNTR SourceConstraintDialogPtr;


static void ChangeSourceConstraintType (GrouP p)
{
  SourceConstraintDialogPtr dlg;
  Int4                      constraint_type;

  dlg = (SourceConstraintDialogPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  Hide (dlg->qual_present);
  Hide (dlg->sc_group);
  Hide (dlg->match_group);
  constraint_type = GetValue (dlg->constraint_type);
  switch (constraint_type) {
    case 1:
      Show (dlg->qual_present);
      break;
    case 2:
      Show (dlg->sc_group);
      break;
    case 3:
      Show (dlg->match_group);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void SourceConstraintToDialog (DialoG d, Pointer data)
{
  SourceConstraintDialogPtr dlg;
  SourceConstraintPtr s;

  dlg = (SourceConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  s = (SourceConstraintPtr) data;
  SetValue (dlg->constraint_type, 1);
  PointerToDialog (dlg->qual_present, NULL);
  PointerToDialog (dlg->qual_string, NULL);
  PointerToDialog (dlg->string_constraint, NULL);
  PointerToDialog (dlg->qualmatch1, NULL);
  PointerToDialog (dlg->qualmatch2, NULL);
  if (s != NULL) {
    if (s->field1 != NULL && s->field2 != NULL) {
      SetValue (dlg->constraint_type, 3);
      PointerToDialog (dlg->qualmatch1, s->field1);
      PointerToDialog (dlg->qualmatch2, s->field2);
    } else if ((s->field1 != NULL || s->field2 != NULL) && s->constraint != NULL) {
      SetValue (dlg->constraint_type, 2);
      if (s->field1 != NULL) {
        PointerToDialog (dlg->qual_string, s->field1);
      } else if (s->field2 != NULL) {
        PointerToDialog (dlg->qual_string, s->field2);
      }
      PointerToDialog (dlg->string_constraint, s->constraint);
    } else if (s->field1 != NULL || s->field2 != NULL) {
      SetValue (dlg->constraint_type, 1);
      if (s->field1 != NULL) {
        PointerToDialog (dlg->qual_present, s->field1);
      } else if (s->field2 != NULL) {
        PointerToDialog (dlg->qual_present, s->field2);
      }
    }
  }
  ChangeSourceConstraintType (dlg->constraint_type);
}


static Pointer DialogToSourceConstraint (DialoG d)
{
  SourceConstraintDialogPtr dlg;
  SourceConstraintPtr s = NULL;
  Int2                constraint_type;

  dlg = (SourceConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  
  constraint_type = GetValue (dlg->constraint_type);
  switch (constraint_type) {
    case 1:
      s = SourceConstraintNew();
      s->field1 = DialogToPointer (dlg->qual_present);
      break;
    case 2:
      s = SourceConstraintNew();
      s->field1 = DialogToPointer (dlg->qual_string);
      s->constraint = DialogToPointer (dlg->string_constraint);
      break;
    case 3:
      s = SourceConstraintNew();
      s->field1 = DialogToPointer (dlg->qualmatch1);
      s->field2 = DialogToPointer (dlg->qualmatch2);
      break;
  }
  return s;
}


static ValNodePtr TestSourceConstraintDialog (DialoG d)
{
  SourceConstraintPtr constraint;
  ValNodePtr err_list = NULL;

  constraint = DialogToPointer (d);
  if (IsSourceConstraintEmpty (constraint) || (constraint->field1 == NULL && constraint->field2 == NULL)) {
    ValNodeAddPointer (&err_list, 0, "empty source constraint");
  }
  constraint = SourceConstraintFree (constraint);
  return err_list;
}


static DialoG SourceConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  SourceConstraintDialogPtr dlg;
  GrouP                     p, g;
  ButtoN                    b1, b2, b3;
  
  dlg = (SourceConstraintDialogPtr) MemNew (sizeof (SourceConstraintDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SourceConstraintToDialog;
  dlg->fromdialog = DialogToSourceConstraint;
  dlg->testdialog = TestSourceConstraintDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->constraint_type = HiddenGroup (p, 3, 0, ChangeSourceConstraintType);
  SetObjectExtra (dlg->constraint_type, dlg, NULL);
  b1 = RadioButton (dlg->constraint_type, "When qualifier present");
  b2 = RadioButton (dlg->constraint_type, "String constraint");
  b3 = RadioButton (dlg->constraint_type, "When qualifiers match");

  g = HiddenGroup (p, 0, 0, NULL);
  
  dlg->qual_present = SourceQualChoiceDialog (g, TRUE, FALSE, TRUE, change_notify, change_userdata);

  dlg->sc_group = HiddenGroup (g, 2, 0, NULL);
  dlg->qual_string = SourceQualChoiceDialog (dlg->sc_group, TRUE, FALSE, TRUE, change_notify, change_userdata);
  dlg->string_constraint = StringConstraintDialog (dlg->sc_group, NULL, FALSE, change_notify, change_userdata);

  dlg->match_group = HiddenGroup (g, 3, 0, NULL);
  dlg->qualmatch1 = SourceQualChoiceDialog (dlg->match_group, TRUE, FALSE, TRUE, change_notify, change_userdata);
  StaticPrompt (dlg->match_group, "Equals", 0, dialogTextHeight, systemFont, 'l');
  dlg->qualmatch2 = SourceQualChoiceDialog (dlg->match_group, TRUE, FALSE, TRUE, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->qual_present, (HANDLE) dlg->sc_group, (HANDLE) dlg->match_group, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->constraint_type, (HANDLE) g, NULL);    
    
  SetValue (dlg->constraint_type, 1);
  ChangeSourceConstraintType (dlg->constraint_type); 
  return (DialoG) p;
}


static DialoG CDSGeneProtFieldDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);


typedef struct cdsgeneprotpseudoconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG feat;
  GrouP  pseudo;
} CDSGeneProtPseudoConstraintDlgData, PNTR CDSGeneProtPseudoConstraintDlgPtr;


static void CDSGeneProtPseudoConstraintToDialog (DialoG d, Pointer data)
{
  CDSGeneProtPseudoConstraintDlgPtr dlg;
  CDSGeneProtPseudoConstraintPtr pseudo;
  ValNode vn;

  dlg = (CDSGeneProtPseudoConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  pseudo = (CDSGeneProtPseudoConstraintPtr) data;

  if (pseudo == NULL) {
    PointerToDialog (dlg->feat, NULL);
    SetValue (dlg->pseudo, 1);
  } else {
    vn.choice = (Uint1) pseudo->feature;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->feat, &vn);
    if (pseudo->is_pseudo) {
      SetValue (dlg->pseudo, 1);
    } else {
      SetValue (dlg->pseudo, 2);
    }
  }
}


static Pointer DialogToCDSGeneProtPseudoConstraint (DialoG d)
{
  CDSGeneProtPseudoConstraintDlgPtr dlg;
  CDSGeneProtPseudoConstraintPtr constraint;
  Uint2 feat_type;
  ValNodePtr vnp;
  Int2       val;

  dlg = (CDSGeneProtPseudoConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  val = GetValue (dlg->pseudo);
  if (val < 1 || val > 2) {
    return NULL;
  }
  vnp = (ValNodePtr) DialogToPointer (dlg->feat);
  if (vnp == NULL) return NULL;
  feat_type = vnp->choice;
  vnp = ValNodeFree (vnp);

  constraint = CDSGeneProtPseudoConstraintNew ();
  constraint->feature = feat_type;
  constraint->is_pseudo = val == 1 ? TRUE : FALSE;
  
  return (Pointer) constraint;
}


static ValNodePtr TestCDSGeneProtPseudoConstraintDialog (DialoG d)
{
  CDSGeneProtPseudoConstraintDlgPtr dlg;
  ValNodePtr err_list = NULL;
  Int4       val;

  dlg = (CDSGeneProtPseudoConstraintDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    val = GetValue (dlg->pseudo);
    if (val < 1 || val > 2) {
      ValNodeAddPointer (&err_list, 0, "pseudo");
    }
    ValNodeLink (&err_list, TestDialog (dlg->feat));
  }
  return err_list;
}

static DialoG CDSGeneProtPseudoConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  CDSGeneProtPseudoConstraintDlgPtr dlg;
  ValNodePtr feat_list = NULL;
  GrouP p;
  
  p = HiddenGroup (h, -1, 0, NULL);
  dlg = (CDSGeneProtPseudoConstraintDlgPtr) MemNew (sizeof (CDSGeneProtPseudoConstraintDlgData));
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = CDSGeneProtPseudoConstraintToDialog;
  dlg->fromdialog = DialogToCDSGeneProtPseudoConstraint;
  dlg->testdialog = TestCDSGeneProtPseudoConstraintDialog;

  AddAllCDSGeneProtFeaturesToChoiceList (&feat_list);
  dlg->feat = ValNodeSelectionDialog (p, feat_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                change_notify, change_userdata, FALSE);

  dlg->pseudo = HiddenGroup (p, 2, 0, NULL);
  RadioButton (dlg->pseudo, "Is pseudo");
  RadioButton (dlg->pseudo, "Is not pseudo");
  SetValue (dlg->pseudo, 1);
  return (DialoG) p;
}


typedef struct cdsgeneprotqualconstraintdialog 
{
  DIALOG_MESSAGE_BLOCK
  GrouP  constraint_type;

  DialoG qual_present;
  GrouP sc_group;
  DialoG qual_string;
  DialoG string_constraint;
  GrouP match_group;
  DialoG qualmatch1;
  DialoG qualmatch2;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} CDSGeneProtQualConstraintDialogData, PNTR CDSGeneProtQualConstraintDialogPtr;


static void ChangeCDSGeneProtQualConstraintType (GrouP p)
{
  CDSGeneProtQualConstraintDialogPtr dlg;
  Int4                      constraint_type;

  dlg = (CDSGeneProtQualConstraintDialogPtr) GetObjectExtra (p);
  if (dlg == NULL)
  {
    return;
  }
  Hide (dlg->qual_present);
  Hide (dlg->sc_group);
  Hide (dlg->match_group);
  constraint_type = GetValue (dlg->constraint_type);
  switch (constraint_type) {
    case 1:
      Show (dlg->qual_present);
      break;
    case 2:
      Show (dlg->sc_group);
      break;
    case 3:
      Show (dlg->match_group);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void CDSGeneProtQualConstraintToDialog (DialoG d, Pointer data)
{
  CDSGeneProtQualConstraintDialogPtr dlg;
  CDSGeneProtQualConstraintPtr s;

  dlg = (CDSGeneProtQualConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  s = (CDSGeneProtQualConstraintPtr) data;
  SetValue (dlg->constraint_type, 1);
  PointerToDialog (dlg->qual_present, NULL);
  PointerToDialog (dlg->qual_string, NULL);
  PointerToDialog (dlg->string_constraint, NULL);
  PointerToDialog (dlg->qualmatch1, NULL);
  PointerToDialog (dlg->qualmatch2, NULL);
  if (s != NULL) {
    if (s->field1 != NULL && s->field2 != NULL) {
      SetValue (dlg->constraint_type, 3);
      PointerToDialog (dlg->qualmatch1, s->field1);
      PointerToDialog (dlg->qualmatch2, s->field2);
    } else if ((s->field1 != NULL || s->field2 != NULL) && s->constraint != NULL) {
      SetValue (dlg->constraint_type, 2);
      if (s->field1 != NULL) {
        PointerToDialog (dlg->qual_string, s->field1);
      } else if (s->field2 != NULL) {
        PointerToDialog (dlg->qual_string, s->field2);
      }
      PointerToDialog (dlg->string_constraint, s->constraint);
    } else if (s->field1 != NULL || s->field2 != NULL) {
      SetValue (dlg->constraint_type, 1);
      if (s->field1 != NULL) {
        PointerToDialog (dlg->qual_present, s->field1);
      } else if (s->field2 != NULL) {
        PointerToDialog (dlg->qual_present, s->field2);
      }
    }
  }
  ChangeCDSGeneProtQualConstraintType (dlg->constraint_type);
}


static Pointer DialogToCDSGeneProtQualConstraint (DialoG d)
{
  CDSGeneProtQualConstraintDialogPtr dlg;
  CDSGeneProtQualConstraintPtr s = NULL;
  Int2                constraint_type;
  ValNodePtr          qual;

  dlg = (CDSGeneProtQualConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  
  constraint_type = GetValue (dlg->constraint_type);
  switch (constraint_type) {
    case 1:
      s = CDSGeneProtQualConstraintNew();
      qual = DialogToPointer (dlg->qual_present);
      if (qual != NULL) {
        ValNodeAddInt (&(s->field1), CDSGeneProtConstraintField_field, qual->data.intvalue);
      }
      break;
    case 2:
      s = CDSGeneProtQualConstraintNew();
      qual = DialogToPointer (dlg->qual_string);
      if (qual != NULL) {
        ValNodeAddInt (&(s->field1), CDSGeneProtConstraintField_field, qual->data.intvalue);
      }
      s->constraint = DialogToPointer (dlg->string_constraint);
      break;
    case 3:
      s = CDSGeneProtQualConstraintNew();
      qual = DialogToPointer (dlg->qualmatch1);
      if (qual != NULL) {
        ValNodeAddInt (&(s->field1), CDSGeneProtConstraintField_field, qual->data.intvalue);
      }
      qual = DialogToPointer (dlg->qualmatch2);
      if (qual != NULL) {
        ValNodeAddInt (&(s->field2), CDSGeneProtConstraintField_field, qual->data.intvalue);
      }
      break;
  }
  return s;
}


static ValNodePtr TestCDSGeneProtQualConstraintDialog (DialoG d)
{
  CDSGeneProtQualConstraintPtr constraint;
  ValNodePtr err_list = NULL;

  constraint = DialogToPointer (d);
  if (IsCDSGeneProtQualConstraintEmpty (constraint)) {
    ValNodeAddPointer (&err_list, 0, "empty CDS-gene-prot qual constraint");
  }
  constraint = CDSGeneProtQualConstraintFree (constraint);
  return err_list;
}

  
static DialoG CDSGeneProtQualConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  CDSGeneProtQualConstraintDialogPtr dlg;
  GrouP                     p, g;
  ButtoN                    b1, b2, b3;
  
  dlg = (CDSGeneProtQualConstraintDialogPtr) MemNew (sizeof (CDSGeneProtQualConstraintDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = CDSGeneProtQualConstraintToDialog;
  dlg->fromdialog = DialogToCDSGeneProtQualConstraint;
  dlg->testdialog = TestCDSGeneProtQualConstraintDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->constraint_type = HiddenGroup (p, 3, 0, ChangeCDSGeneProtQualConstraintType);
  SetObjectExtra (dlg->constraint_type, dlg, NULL);
  b1 = RadioButton (dlg->constraint_type, "When qualifier present");
  b2 = RadioButton (dlg->constraint_type, "String constraint");
  b3 = RadioButton (dlg->constraint_type, "When qualifiers match");

  g = HiddenGroup (p, 0, 0, NULL);
  
  dlg->qual_present = CDSGeneProtFieldDialog (g, change_notify, change_userdata);

  dlg->sc_group = HiddenGroup (g, 2, 0, NULL);
  dlg->qual_string = CDSGeneProtFieldDialog (dlg->sc_group, change_notify, change_userdata);
  dlg->string_constraint = StringConstraintDialog (dlg->sc_group, NULL, FALSE, change_notify, change_userdata);

  dlg->match_group = HiddenGroup (g, 3, 0, NULL);
  dlg->qualmatch1 = CDSGeneProtFieldDialog (dlg->match_group, change_notify, change_userdata);
  StaticPrompt (dlg->match_group, "Equals", 0, dialogTextHeight, systemFont, 'l');
  dlg->qualmatch2 = CDSGeneProtFieldDialog (dlg->match_group, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->qual_present, (HANDLE) dlg->sc_group, (HANDLE) dlg->match_group, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->constraint_type, (HANDLE) g, NULL);    
  SetValue (dlg->constraint_type, 1);
  ChangeCDSGeneProtQualConstraintType (dlg->constraint_type); 
  return (DialoG) p;
}


typedef struct sequencequaldlg {
  DIALOG_MESSAGE_BLOCK
  PopuP field_type;
  DialoG molecule_dlg;
  DialoG technique_dlg;
  DialoG completedness_dlg;
  DialoG mol_class_dlg;
  DialoG topology_dlg;
  DialoG strand_dlg;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} SequenceQualDlgData, PNTR SequenceQualDlgPtr;


static void ChangeSequenceQualType (PopuP p)
{
  SequenceQualDlgPtr dlg;
  Int2               val;

  dlg = (SequenceQualDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  Hide (dlg->molecule_dlg);
  Hide (dlg->technique_dlg);
  Hide (dlg->completedness_dlg);
  Hide (dlg->mol_class_dlg);
  Hide (dlg->topology_dlg);
  Hide (dlg->strand_dlg);

  val = GetValue (dlg->field_type);
  switch (val) {
    case 1:
      Show (dlg->molecule_dlg);
      break;
    case 2:
      Show (dlg->technique_dlg);
      break;
    case 3:
      Show (dlg->completedness_dlg);
      break;
    case 4:
      Show (dlg->mol_class_dlg);
      break;
    case 5:
      Show (dlg->topology_dlg);
      break;
    case 6:
      Show (dlg->strand_dlg);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SequenceQualToDialog (DialoG d, Pointer data)
{
  SequenceQualDlgPtr dlg;
  ValNodePtr         field;
  ValNode            vn;

  dlg = (SequenceQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  field = (ValNodePtr) data;

  if (field == NULL) {
    SetValue (dlg->field_type, 1);
    PointerToDialog (dlg->molecule_dlg, NULL);
  } else {
    vn.next = NULL;
    vn.data.ptrvalue = NULL;

    switch (field->choice) {
      case MolinfoField_molecule:
        SetValue (dlg->field_type, 1);
        vn.choice = field->data.intvalue;
        PointerToDialog (dlg->molecule_dlg, &vn);
        break;
      case MolinfoField_technique:
        SetValue (dlg->field_type, 2);
        vn.choice = field->data.intvalue;
        PointerToDialog (dlg->technique_dlg, &vn);
        break;
      case MolinfoField_completedness:
        SetValue (dlg->field_type, 3);
        vn.choice = field->data.intvalue;
        PointerToDialog (dlg->completedness_dlg, &vn);
        break;
      case MolinfoField_mol_class:
        SetValue (dlg->field_type, 4);
        vn.choice = field->data.intvalue;
        PointerToDialog (dlg->mol_class_dlg, &vn);
        break;
      case MolinfoField_topology:
        SetValue (dlg->field_type, 5);
        vn.choice = field->data.intvalue;
        PointerToDialog (dlg->topology_dlg, &vn);
        break;
      case MolinfoField_strand:
        SetValue (dlg->field_type, 6);
        vn.choice = field->data.intvalue;
        PointerToDialog (dlg->strand_dlg, &vn);
        break;
      default:
        SetValue (dlg->field_type, 1);
        PointerToDialog (dlg->molecule_dlg, NULL);
        break;
    }
  }
  ChangeSequenceQualType (dlg->field_type);
}


static Pointer DialogToSequenceQual (DialoG d)
{
  SequenceQualDlgPtr dlg;
  ValNodePtr         field = NULL, vnp;
  Int2               val;

  dlg = (SequenceQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->field_type);
  switch (val) {
    case 1:
      vnp = DialogToPointer (dlg->molecule_dlg);
      if (vnp != NULL) {
        field = ValNodeNew (NULL);
        field->choice = MolinfoField_molecule;
        field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      break;
    case 2:
      vnp = DialogToPointer (dlg->technique_dlg);
      if (vnp != NULL) {
        field = ValNodeNew (NULL);
        field->choice = MolinfoField_technique;
        field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      break;
    case 3:
      vnp = DialogToPointer (dlg->completedness_dlg);
      if (vnp != NULL) {
        field = ValNodeNew (NULL);
        field->choice = MolinfoField_completedness;
        field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      break;
    case 4:
      vnp = DialogToPointer (dlg->mol_class_dlg);
      if (vnp != NULL) {
        field = ValNodeNew (NULL);
        field->choice = MolinfoField_mol_class;
        field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      break;
    case 5:
      vnp = DialogToPointer (dlg->topology_dlg);
      if (vnp != NULL) {
        field = ValNodeNew (NULL);
        field->choice = MolinfoField_topology;
        field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      break;
    case 6:
      vnp = DialogToPointer (dlg->strand_dlg);
      if (vnp != NULL) {
        field = ValNodeNew (NULL);
        field->choice = MolinfoField_strand;
        field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      break;
  }

  return field;    
  
}


static ValNodePtr TestSequenceQualDialog (DialoG d)
{
  SequenceQualDlgPtr dlg;
  Int2               val;
  ValNodePtr         err_list = NULL;

  dlg = (SequenceQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->field_type);
  if (val < 1 || val > 6) {
  }
  switch (val) {
    case 1:
      err_list = TestDialog (dlg->molecule_dlg);
      break;
    case 2:
      err_list = TestDialog (dlg->technique_dlg);
      break;
    case 3:
      err_list = TestDialog (dlg->completedness_dlg);
      break;
    case 4:
      err_list = TestDialog (dlg->mol_class_dlg);
      break;
    case 5:
      err_list = TestDialog (dlg->topology_dlg);
      break;
    case 6:
      err_list = TestDialog (dlg->strand_dlg);
      break;
    default:
      ValNodeAddPointer (&err_list, 0, "field type");
      break;
  }

  return err_list;
}


static DialoG SequenceQualDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  SequenceQualDlgPtr dlg;
  GrouP              p, g;
  ValNodePtr         val_list;
  
  dlg = (SequenceQualDlgPtr) MemNew (sizeof (SequenceQualDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SequenceQualToDialog;
  dlg->fromdialog = DialogToSequenceQual;
  dlg->testdialog = TestSequenceQualDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->field_type = PopupList (p, TRUE, ChangeSequenceQualType);
  SetObjectExtra (dlg->field_type, dlg, NULL);
  PopupItem (dlg->field_type, "molecule");
  PopupItem (dlg->field_type, "technique");
  PopupItem (dlg->field_type, "completedness");
  PopupItem (dlg->field_type, "class");
  PopupItem (dlg->field_type, "topology");
  PopupItem (dlg->field_type, "strand");
  SetValue (dlg->field_type, 1);

  g = HiddenGroup (p, 0, 0, NULL);
  
  val_list = GetMoleculeTypeList ();
  dlg->molecule_dlg = ValNodeSelectionDialog (g, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetTechniqueTypeList ();
  dlg->technique_dlg = ValNodeSelectionDialog (g, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetCompletednessTypeList ();
  dlg->completedness_dlg = ValNodeSelectionDialog (g, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetMoleculeClassTypeList ();
  dlg->mol_class_dlg = ValNodeSelectionDialog (g, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetTopologyTypeList ();
  dlg->topology_dlg = ValNodeSelectionDialog (g, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetStrandTypeList ();
  dlg->strand_dlg = ValNodeSelectionDialog (g, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  AlignObjects (ALIGN_MIDDLE, (HANDLE) dlg->molecule_dlg, 
                              (HANDLE) dlg->technique_dlg,
                              (HANDLE) dlg->completedness_dlg,
                              (HANDLE) dlg->mol_class_dlg,
                              (HANDLE) dlg->topology_dlg,
                              (HANDLE) dlg->strand_dlg,
                              NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->field_type, (HANDLE) g, NULL);    
  ChangeSequenceQualType (dlg->field_type); 
  return (DialoG) p;
}


typedef struct sequencequalpairdlg {
  DIALOG_MESSAGE_BLOCK
  PopuP field_type;
  GrouP molecule_grp;
  DialoG molecule_from_dlg;
  DialoG molecule_to_dlg;
  GrouP technique_grp;
  DialoG technique_from_dlg;
  DialoG technique_to_dlg;
  GrouP completedness_grp;
  DialoG completedness_from_dlg;
  DialoG completedness_to_dlg;
  GrouP mol_class_grp;
  DialoG mol_class_from_dlg;
  DialoG mol_class_to_dlg;
  GrouP topology_grp;
  DialoG topology_from_dlg;
  DialoG topology_to_dlg;
  GrouP strand_grp;
  DialoG strand_from_dlg;
  DialoG strand_to_dlg;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} SequenceQualPairDlgData, PNTR SequenceQualPairDlgPtr;


static void ChangeSequenceQualPairType (PopuP p)
{
  SequenceQualPairDlgPtr dlg;
  Int2                   val;

  dlg = (SequenceQualPairDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  Hide (dlg->molecule_grp);
  Hide (dlg->technique_grp);
  Hide (dlg->completedness_grp);
  Hide (dlg->mol_class_grp);
  Hide (dlg->topology_grp);
  Hide (dlg->strand_grp);

  val = GetValue (dlg->field_type);
  switch (val) {
    case 1:
      Show (dlg->molecule_grp);
      break;
    case 2:
      Show (dlg->technique_grp);
      break;
    case 3:
      Show (dlg->completedness_grp);
      break;
    case 4:
      Show (dlg->mol_class_grp);
      break;
    case 5:
      Show (dlg->topology_grp);
      break;
    case 6:
      Show (dlg->strand_grp);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SequenceQualPairToDialog (DialoG d, Pointer data)
{
  SequenceQualPairDlgPtr dlg;
  ValNodePtr m_fields;
  ValNode    vn;

  dlg = (SequenceQualPairDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  m_fields = (ValNodePtr) data;
  if (m_fields == NULL) {
    SetValue (dlg->field_type, 1);
    PointerToDialog (dlg->molecule_from_dlg, NULL);
    PointerToDialog (dlg->molecule_to_dlg, NULL);
  } else {
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    switch (m_fields->choice) {
      case MolinfoFieldPair_molecule:
        SetValue (dlg->field_type, 1);
        if (m_fields->data.ptrvalue == NULL) {
          PointerToDialog (dlg->molecule_from_dlg, NULL);
          PointerToDialog (dlg->molecule_to_dlg, NULL);
        } else {
          vn.choice = (Uint1)((MolinfoMoleculePairPtr) m_fields->data.ptrvalue)->from;
          PointerToDialog (dlg->molecule_from_dlg, &vn);
          vn.choice = (Uint1)((MolinfoMoleculePairPtr) m_fields->data.ptrvalue)->to;
          PointerToDialog (dlg->molecule_to_dlg, &vn);
        }
        break;
      case MolinfoFieldPair_technique:
        SetValue (dlg->field_type, 2);
        if (m_fields->data.ptrvalue == NULL) {
          PointerToDialog (dlg->technique_from_dlg, NULL);
          PointerToDialog (dlg->technique_to_dlg, NULL);
        } else {
          vn.choice = (Uint1)((MolinfoTechniquePairPtr) m_fields->data.ptrvalue)->from;
          PointerToDialog (dlg->technique_from_dlg, &vn);
          vn.choice = (Uint1)((MolinfoTechniquePairPtr) m_fields->data.ptrvalue)->to;
          PointerToDialog (dlg->technique_to_dlg, &vn);
        }
        break;
      case MolinfoFieldPair_completedness:
        SetValue (dlg->field_type, 3);
        if (m_fields->data.ptrvalue == NULL) {
          PointerToDialog (dlg->completedness_from_dlg, NULL);
          PointerToDialog (dlg->completedness_to_dlg, NULL);
        } else {
          vn.choice = (Uint1)((MolinfoCompletednessPairPtr) m_fields->data.ptrvalue)->from;
          PointerToDialog (dlg->completedness_from_dlg, &vn);
          vn.choice = (Uint1)((MolinfoCompletednessPairPtr) m_fields->data.ptrvalue)->to;
          PointerToDialog (dlg->completedness_to_dlg, &vn);
        }
        break;
      case MolinfoFieldPair_mol_class:
        SetValue (dlg->field_type, 4);
        if (m_fields->data.ptrvalue == NULL) {
          PointerToDialog (dlg->mol_class_from_dlg, NULL);
          PointerToDialog (dlg->mol_class_to_dlg, NULL);
        } else {
          vn.choice = (Uint1)((MolinfoMolClassPairPtr) m_fields->data.ptrvalue)->from;
          PointerToDialog (dlg->mol_class_from_dlg, &vn);
          vn.choice = (Uint1)((MolinfoMolClassPairPtr) m_fields->data.ptrvalue)->to;
          PointerToDialog (dlg->mol_class_to_dlg, &vn);
        }
        break;
      case MolinfoFieldPair_topology:
        SetValue (dlg->field_type, 5);
        if (m_fields->data.ptrvalue == NULL) {
          PointerToDialog (dlg->topology_from_dlg, NULL);
          PointerToDialog (dlg->topology_to_dlg, NULL);
        } else {
          vn.choice = (Uint1)((MolinfoTopologyPairPtr) m_fields->data.ptrvalue)->from;
          PointerToDialog (dlg->topology_from_dlg, &vn);
          vn.choice = (Uint1)((MolinfoTopologyPairPtr) m_fields->data.ptrvalue)->to;
          PointerToDialog (dlg->topology_to_dlg, &vn);
        }
        break;
      case MolinfoFieldPair_strand:
        SetValue (dlg->field_type, 6);
        if (m_fields->data.ptrvalue == NULL) {
          PointerToDialog (dlg->strand_from_dlg, NULL);
          PointerToDialog (dlg->strand_to_dlg, NULL);
        } else {
          vn.choice = (Uint1)((MolinfoStrandPairPtr) m_fields->data.ptrvalue)->from;
          PointerToDialog (dlg->strand_from_dlg, &vn);
          vn.choice = (Uint1)((MolinfoStrandPairPtr) m_fields->data.ptrvalue)->to;
          PointerToDialog (dlg->strand_to_dlg, &vn);
        }
        break;
      default:
        SetValue (dlg->field_type, 1);
        PointerToDialog (dlg->molecule_from_dlg, NULL);
        PointerToDialog (dlg->molecule_to_dlg, NULL);
        break;
    }
  }
}


static Pointer DialogToSequenceQualPair (DialoG d)
{
  SequenceQualPairDlgPtr dlg;
  Int2                   val;
  ValNodePtr             field = NULL, vnp_from, vnp_to;
  MolinfoMoleculePairPtr mol;
  MolinfoTechniquePairPtr tech;
  MolinfoCompletednessPairPtr comp;
  MolinfoMolClassPairPtr mol_class;
  MolinfoTopologyPairPtr topology;
  MolinfoStrandPairPtr   strand;

  dlg = (SequenceQualPairDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->field_type);
  switch (val) {
    case 1:
      field = ValNodeNew (NULL);
      field->choice = MolinfoFieldPair_molecule;
      vnp_from = DialogToPointer (dlg->molecule_from_dlg);
      vnp_to = DialogToPointer (dlg->molecule_to_dlg);
      mol = MolinfoMoleculePairNew ();
      field->data.ptrvalue = mol;
      if (vnp_from != NULL) {
        mol->from = vnp_from->choice;
        vnp_from = ValNodeFree (vnp_from);
      }
      if (vnp_to != NULL) {
        mol->to = vnp_to->choice;
        vnp_to = ValNodeFree (vnp_to);
      }
      break;
    case 2:
      field = ValNodeNew (NULL);
      field->choice = MolinfoFieldPair_technique;
      tech = MolinfoTechniquePairNew ();
      field->data.ptrvalue = tech;
      vnp_from = DialogToPointer (dlg->technique_from_dlg);
      vnp_to = DialogToPointer (dlg->technique_to_dlg);
      if (vnp_from != NULL) {
        tech->from = vnp_from->choice;
        vnp_from = ValNodeFree (vnp_from);
      }
      if (vnp_to != NULL) {
        tech->to = vnp_to->choice;
        vnp_to = ValNodeFree (vnp_to);
      }
      break;
    case 3:
      field = ValNodeNew (NULL);
      field->choice = MolinfoFieldPair_completedness;
      comp = MolinfoCompletednessPairNew ();
      field->data.ptrvalue = comp;
      vnp_from = DialogToPointer (dlg->completedness_from_dlg);
      vnp_to = DialogToPointer (dlg->completedness_to_dlg);
      if (vnp_from != NULL) {
        comp->from = vnp_from->choice;
        vnp_from = ValNodeFree (vnp_from);
      }
      if (vnp_to != NULL) {
        comp->to = vnp_to->choice;
        vnp_to = ValNodeFree (vnp_to);
      }
      break;
    case 4:
      field = ValNodeNew (NULL);
      field->choice = MolinfoFieldPair_technique;
      mol_class = MolinfoMolClassPairNew ();
      field->data.ptrvalue = mol_class;
      vnp_from = DialogToPointer (dlg->mol_class_from_dlg);
      vnp_to = DialogToPointer (dlg->mol_class_to_dlg);
      if (vnp_from != NULL) {
        mol_class->from = vnp_from->choice;
        vnp_from = ValNodeFree (vnp_from);
      }
      if (vnp_to != NULL) {
        mol_class->to = vnp_to->choice;
        vnp_to = ValNodeFree (vnp_to);
      }
      break;
    case 5:
      field = ValNodeNew (NULL);
      field->choice = MolinfoFieldPair_topology;
      topology = MolinfoTopologyPairNew ();
      field->data.ptrvalue = topology;
      vnp_from = DialogToPointer (dlg->topology_from_dlg);
      vnp_to = DialogToPointer (dlg->topology_to_dlg);
      if (vnp_from != NULL) {
        topology->from = vnp_from->choice;
        vnp_from = ValNodeFree (vnp_from);
      }
      if (vnp_to != NULL) {
        topology->to = vnp_to->choice;
        vnp_to = ValNodeFree (vnp_to);
      }
      break;
    case 6:
      field = ValNodeNew (NULL);
      field->choice = MolinfoFieldPair_strand;
      strand = MolinfoStrandPairNew ();
      field->data.ptrvalue = strand;
      vnp_from = DialogToPointer (dlg->strand_from_dlg);
      vnp_to = DialogToPointer (dlg->strand_to_dlg);
      if (vnp_from != NULL) {
        strand->from = vnp_from->choice;
        vnp_from = ValNodeFree (vnp_from);
      }
      if (vnp_to != NULL) {
        strand->to = vnp_to->choice;
        vnp_to = ValNodeFree (vnp_to);
      }
      break;
  }
  return field;
}


static ValNodePtr TestSequenceQualPairDialog (DialoG d)
{
  ValNodePtr err_list = NULL, field;
  MolinfoMoleculePairPtr mol;
  MolinfoTechniquePairPtr tech;
  MolinfoCompletednessPairPtr comp;
  MolinfoMolClassPairPtr mol_class;
  MolinfoTopologyPairPtr topology;
  MolinfoStrandPairPtr   strand;

  field = DialogToPointer (d);
  if (field == NULL) {
    ValNodeAddPointer (&err_list, 0, "field");
  } else {
    switch (field->choice) {
      case MolinfoFieldPair_molecule:
        mol = (MolinfoMoleculePairPtr) field->data.ptrvalue;
        if (mol == NULL || mol->from == mol->to) {
          ValNodeAddPointer (&err_list, 0, "value");
        } 
        break;
      case MolinfoFieldPair_technique:
        tech = (MolinfoTechniquePairPtr) field->data.ptrvalue;
        if (tech == NULL || tech->from == tech->to) {
          ValNodeAddPointer (&err_list, 0, "value");
        } 
        break;
      case MolinfoFieldPair_completedness:
        comp = (MolinfoCompletednessPairPtr) field->data.ptrvalue;
        if (comp == NULL || comp->from == comp->to) {
          ValNodeAddPointer (&err_list, 0, "value");
        } 
        break;
      case MolinfoFieldPair_mol_class:
        mol_class = (MolinfoMolClassPairPtr) field->data.ptrvalue;
        if (mol_class == NULL) {
          ValNodeAddPointer (&err_list, 0, "value");
        } 
        break;
      case MolinfoFieldPair_topology:
        topology = (MolinfoTopologyPairPtr) field->data.ptrvalue;
        if (topology == NULL || topology->from == topology->to) {
          ValNodeAddPointer (&err_list, 0, "value");
        } 
        break;
      case MolinfoFieldPair_strand:
        strand = (MolinfoStrandPairPtr) field->data.ptrvalue;
        if (strand == NULL || strand->from == strand->to) {
          ValNodeAddPointer (&err_list, 0, "value");
        } 
        break;
      default:
        ValNodeAddPointer (&err_list, 0, "type");
        break;
    }
    field = MolinfoFieldPairFree (field);
  }

  return err_list;
}


static DialoG SequenceQualPairDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{
  SequenceQualPairDlgPtr dlg;
  GrouP                  p, g;
  ValNodePtr             val_list;
  
  dlg = (SequenceQualPairDlgPtr) MemNew (sizeof (SequenceQualPairDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SequenceQualPairToDialog;
  dlg->fromdialog = DialogToSequenceQualPair;
  dlg->testdialog = TestSequenceQualPairDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->field_type = PopupList (p, TRUE, ChangeSequenceQualPairType);
  SetObjectExtra (dlg->field_type, dlg, NULL);
  PopupItem (dlg->field_type, "molecule");
  PopupItem (dlg->field_type, "technique");
  PopupItem (dlg->field_type, "completedness");
  PopupItem (dlg->field_type, "class");
  PopupItem (dlg->field_type, "topology");
  PopupItem (dlg->field_type, "strand");
  SetValue (dlg->field_type, 1);

  g = HiddenGroup (p, 0, 0, NULL);

  dlg->molecule_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->molecule_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->molecule_grp, "To", 0, dialogTextHeight, programFont, 'l');
  val_list = GetMoleculeTypeList ();
  dlg->molecule_from_dlg = ValNodeSelectionDialog (dlg->molecule_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetMoleculeTypeList ();
  dlg->molecule_to_dlg = ValNodeSelectionDialog (dlg->molecule_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  dlg->technique_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->technique_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->technique_grp, "To", 0, dialogTextHeight, programFont, 'l');
  val_list = GetTechniqueTypeList ();
  dlg->technique_from_dlg = ValNodeSelectionDialog (dlg->technique_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetTechniqueTypeList ();
  dlg->technique_to_dlg = ValNodeSelectionDialog (dlg->technique_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  dlg->completedness_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->completedness_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->completedness_grp, "To", 0, dialogTextHeight, programFont, 'l');
  val_list = GetCompletednessTypeList ();
  dlg->completedness_from_dlg = ValNodeSelectionDialog (dlg->completedness_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetCompletednessTypeList ();
  dlg->completedness_to_dlg = ValNodeSelectionDialog (dlg->completedness_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  dlg->mol_class_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->mol_class_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->mol_class_grp, "To", 0, dialogTextHeight, programFont, 'l');
  val_list = GetMoleculeClassTypeList ();
  dlg->mol_class_from_dlg = ValNodeSelectionDialog (dlg->mol_class_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetMoleculeClassTypeList ();
  dlg->mol_class_to_dlg = ValNodeSelectionDialog (dlg->mol_class_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  dlg->topology_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->topology_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->topology_grp, "To", 0, dialogTextHeight, programFont, 'l');
  val_list = GetTopologyTypeList ();
  dlg->topology_from_dlg = ValNodeSelectionDialog (dlg->topology_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetTopologyTypeList ();
  dlg->topology_to_dlg = ValNodeSelectionDialog (dlg->topology_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  dlg->strand_grp = HiddenGroup (g, 2, 0, NULL);
  StaticPrompt (dlg->strand_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->strand_grp, "To", 0, dialogTextHeight, programFont, 'l');

  val_list = GetStrandTypeList ();
  dlg->strand_from_dlg = ValNodeSelectionDialog (dlg->strand_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  val_list = GetStrandTypeList ();
  dlg->strand_to_dlg = ValNodeSelectionDialog (dlg->strand_grp, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  AlignObjects (ALIGN_MIDDLE, (HANDLE) dlg->molecule_grp, 
                              (HANDLE) dlg->technique_grp,
                              (HANDLE) dlg->completedness_grp,
                              (HANDLE) dlg->mol_class_grp,
                              (HANDLE) dlg->topology_grp,
                              (HANDLE) dlg->strand_grp,
                              NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->field_type, (HANDLE) g, NULL);    
  ChangeSequenceQualPairType (dlg->field_type); 
  return (DialoG) p;
}


static CharPtr  molinfo_block_labels [] = {
  "Molecule", "Technique", "Completedness",
  "Class", "Topology", "Strand", NULL
};

typedef struct molinfoblocklistdlg {
  DIALOG_MESSAGE_BLOCK

  DialoG molecule;
  DialoG technique;
  DialoG completedness;
  DialoG mol_class;
  DialoG topology;
  DialoG strand;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} MolInfoBlockListDlgData, PNTR MolInfoBlockListDlgPtr;


static void MolInfoBlockListToDialog (DialoG d, Pointer data)
{
  MolInfoBlockListDlgPtr dlg;
  ValNodePtr         mol_fields, vnp;
  ValNode            val;

  dlg = (MolInfoBlockListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  mol_fields = (ValNodePtr) data;

  vnp = ValNodeNew(NULL);
  vnp->choice = 255;
  vnp->data.ptrvalue = StringSave ("any");

  PointerToDialog (dlg->molecule, vnp);
  PointerToDialog (dlg->technique, vnp);
  PointerToDialog (dlg->completedness, vnp);
  PointerToDialog (dlg->mol_class, vnp);
  PointerToDialog (dlg->topology, vnp);
  PointerToDialog (dlg->strand, vnp);

  vnp = ValNodeFreeData (vnp);
  MemSet (&val, 0, sizeof (ValNode));

  for (vnp = mol_fields; vnp != NULL; vnp = vnp->next) {
    val.choice = vnp->data.intvalue;
    switch (vnp->choice) {
      case MolinfoField_molecule:
        PointerToDialog (dlg->molecule, &val);
        break;
      case MolinfoField_technique:
        PointerToDialog (dlg->technique, &val);
        break;
      case MolinfoField_completedness:
        PointerToDialog (dlg->completedness, &val);
        break;
      case MolinfoField_mol_class:
        PointerToDialog (dlg->mol_class, &val);
        break;
      case MolinfoField_topology:
        PointerToDialog (dlg->topology, &val);
        break;
      case MolinfoField_strand:
        PointerToDialog (dlg->strand, &val);
        break;
    }
  }
}


static Pointer MolInfoBlockListFromDialog (DialoG d)
{
  MolInfoBlockListDlgPtr dlg;
  ValNodePtr         mol_fields = NULL, vnp;

  dlg = (MolInfoBlockListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  vnp = DialogToPointer (dlg->molecule);
  if (vnp != NULL && vnp->choice != 255) {
    ValNodeAddInt (&mol_fields, MolinfoField_molecule, vnp->choice);
  }
  vnp = ValNodeFree (vnp);
 
  vnp = DialogToPointer (dlg->technique);
  if (vnp != NULL && vnp->choice != 255) {
    ValNodeAddInt (&mol_fields, MolinfoField_technique, vnp->choice);
  }
  vnp = ValNodeFree (vnp);
 
  vnp = DialogToPointer (dlg->completedness);
  if (vnp != NULL && vnp->choice != 255) {
    ValNodeAddInt (&mol_fields, MolinfoField_completedness, vnp->choice);
  }
  vnp = ValNodeFree (vnp);

  vnp = DialogToPointer (dlg->mol_class);
  if (vnp != NULL && vnp->choice != 255) {
    ValNodeAddInt (&mol_fields, MolinfoField_mol_class, vnp->choice);
  }
  vnp = ValNodeFree (vnp);
 
  vnp = DialogToPointer (dlg->topology);
  if (vnp != NULL && vnp->choice != 255) {
    ValNodeAddInt (&mol_fields, MolinfoField_topology, vnp->choice);
  }
  vnp = ValNodeFree (vnp);

  vnp = DialogToPointer (dlg->strand);
  if (vnp != NULL && vnp->choice != 255) {
    ValNodeAddInt (&mol_fields, MolinfoField_strand, vnp->choice);
  }
  vnp = ValNodeFree (vnp);

  return (Pointer) mol_fields;
}


static DialoG MolInfoBlockListDialog (GrouP h, CharPtr any_name, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  MolInfoBlockListDlgPtr dlg;
  GrouP                  p;
  Int2                   wid;
  ValNodePtr             val_list, vnp;
  
  dlg = (MolInfoBlockListDlgPtr) MemNew (sizeof (MolInfoBlockListDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = MolInfoBlockListToDialog;
  dlg->fromdialog = MolInfoBlockListFromDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  SelectFont (programFont);
  wid = MaxStringWidths (molinfo_block_labels);
  SelectFont (systemFont);

  /* molecule */
  StaticPrompt (p, molinfo_block_labels[0], wid, popupMenuHeight, programFont, 'l');
  val_list = GetMoleculeTypeList ();
  if (any_name != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = 255;
    vnp->data.ptrvalue = StringSave (any_name);
    vnp->next = val_list;
    val_list = vnp;
  }
  dlg->molecule = ValNodeSelectionDialog (p, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  /* technique */
  StaticPrompt (p, molinfo_block_labels[1], wid, popupMenuHeight, programFont, 'l');
  val_list = GetTechniqueTypeList ();
  if (any_name != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = 255;
    vnp->data.ptrvalue = StringSave (any_name);
    vnp->next = val_list;
    val_list = vnp;
  }
  dlg->technique = ValNodeSelectionDialogExEx (p, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE, FALSE, TRUE, NULL);

  /* completedness */
  StaticPrompt (p, molinfo_block_labels[2], wid, popupMenuHeight, programFont, 'l');
  val_list = GetCompletednessTypeList ();
  if (any_name != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = 255;
    vnp->data.ptrvalue = StringSave (any_name);
    vnp->next = val_list;
    val_list = vnp;
  }
  dlg->completedness = ValNodeSelectionDialog (p, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  
  /* mol_class */
  StaticPrompt (p, molinfo_block_labels[3], wid, popupMenuHeight, programFont, 'l');
  val_list = GetMoleculeClassTypeList ();
  if (any_name != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = 255;
    vnp->data.ptrvalue = StringSave (any_name);
    vnp->next = val_list;
    val_list = vnp;
  }
  dlg->mol_class = ValNodeSelectionDialog (p, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  
  /* topology */
  StaticPrompt (p, molinfo_block_labels[4], wid, popupMenuHeight, programFont, 'l');
  val_list = GetTopologyTypeList ();
  if (any_name != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = 255;
    vnp->data.ptrvalue = StringSave (any_name);
    vnp->next = val_list;
    val_list = vnp;
  }
  dlg->topology = ValNodeSelectionDialog (p, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);
  /* strand */
  StaticPrompt (p, molinfo_block_labels[5], wid, popupMenuHeight, programFont, 'l');
  val_list = GetStrandTypeList ();
  if (any_name != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = 255;
    vnp->data.ptrvalue = StringSave (any_name);
    vnp->next = val_list;
    val_list = vnp;
  }
  dlg->strand = ValNodeSelectionDialog (p, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "location", 
                                           change_notify, change_userdata, FALSE);

  return (DialoG) p;
}


typedef struct molinfoblockdlg
{
  DIALOG_MESSAGE_BLOCK

  DialoG to_list;
  DialoG from_list;
  DialoG constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} MolinfoBlockDlgData, PNTR MolinfoBlockDlgPtr;
 

static void MolInfoBlockToDialog (DialoG d, Pointer data)
{
  MolinfoBlockDlgPtr dlg;
  MolinfoBlockPtr    mib;

  dlg = (MolinfoBlockDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  mib = (MolinfoBlockPtr) data;

  if (mib == NULL) {
    PointerToDialog (dlg->from_list, NULL);
    PointerToDialog (dlg->to_list, NULL);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    PointerToDialog (dlg->from_list, mib->from_list);
    PointerToDialog (dlg->to_list, mib->to_list);
    PointerToDialog (dlg->constraint, mib->constraint);
  }
}


static Pointer MolInfoBlockFromDialog (DialoG d)
{
  MolinfoBlockDlgPtr dlg;
  MolinfoBlockPtr    mib;

  dlg = (MolinfoBlockDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  mib = MolinfoBlockNew ();
  mib->from_list = DialogToPointer (dlg->from_list);
  mib->to_list = DialogToPointer (dlg->to_list);
  mib->constraint = DialogToPointer (dlg->constraint);

  return (Pointer) mib;
}


NLM_EXTERN DialoG MolInfoBlockDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  MolinfoBlockDlgPtr dlg;
  GrouP              p, g;

  dlg = (MolinfoBlockDlgPtr) MemNew (sizeof (MolinfoBlockDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = MolInfoBlockToDialog;
  dlg->fromdialog = MolInfoBlockFromDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g = HiddenGroup (p, edit ? 2 : -1 , 0, NULL);

  if (edit) {
    StaticPrompt (g, "From", 0, popupMenuHeight, programFont, 'c');
    StaticPrompt (g, "To", 0, popupMenuHeight, programFont, 'c');
    dlg->from_list = MolInfoBlockListDialog (g, "any", change_notify, change_userdata);
  }
  dlg->to_list = MolInfoBlockListDialog (g, "no change", change_notify, change_userdata);

  dlg->constraint = ComplexConstraintDialog (p, change_notify, change_userdata);
  ChangeComplexConstraintFieldType (dlg->constraint, FieldType_molinfo_field, NULL, Macro_feature_type_any);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlg->constraint, NULL);

  return (DialoG) p;
}


typedef struct capsaction {
  DIALOG_MESSAGE_BLOCK
  Uint1 action_type;
} CapsActionDlgData, PNTR CapsActionDlgPtr;


static Pointer CapsActionFromDialog (DialoG d)
{
  ValNodePtr vnp;
  CapsActionDlgPtr dlg;

  dlg = (CapsActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  vnp = ValNodeNew (NULL);
  vnp->choice = dlg->action_type;
  return vnp;
}


static DialoG CapsActionDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata, Uint1 action_type)
{
  CapsActionDlgPtr dlg;
  GrouP              p;

  dlg = (CapsActionDlgPtr) MemNew (sizeof (CapsActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = CapsActionFromDialog;
  dlg->action_type = action_type;

  return (DialoG) p;
}


static DialoG MusMusculusActionDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return CapsActionDialog (h, edit, change_notify, change_userdata, FixCapsAction_mouse_strain);
}


NLM_EXTERN DialoG FixCapsSourceCountryDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return CapsActionDialog (h, edit, change_notify, change_userdata, FixCapsAction_src_country);
}

typedef struct srcqualcapsaction {
  DIALOG_MESSAGE_BLOCK
  DialoG field;
} SrcQualCapsActionDlgData, PNTR SrcQualCapsActionDlgPtr;


static Pointer SrcQualCapsActionFromDialog (DialoG d)
{
  ValNodePtr vnp;
  SrcQualCapsActionDlgPtr dlg;
  Int4 val;

  dlg = (SrcQualCapsActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  vnp = DialogToPointer (dlg->field);
  val = GetSourceQualTypeByName (vnp->data.ptrvalue);
  vnp = ValNodeFree (vnp);
  vnp = ValNodeNew (NULL);
  vnp->choice = FixCapsAction_src_qual;
  vnp->data.intvalue = val;
  return vnp;
}


static void SrcQualCapsActionToDialog (DialoG d, Pointer data)
{
  ValNodePtr vnp;
  ValNode vn;
  SrcQualCapsActionDlgPtr dlg;

  dlg = (SrcQualCapsActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  vnp = (ValNodePtr) data;

  MemSet (&vn, 0, sizeof (ValNode));
  if (vnp == NULL || vnp->data.intvalue == 0) {
    vn.data.ptrvalue = GetSourceQualName(Source_qual_sex);
  } else {
    vn.data.ptrvalue = GetSourceQualName (vnp->data.intvalue);
  }
  PointerToDialog (dlg->field, &vn);
}


static DialoG FixSrcQualCapsActionDialog (GrouP h,  Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  SrcQualCapsActionDlgPtr dlg;
  GrouP            p;
  ValNodePtr       list = NULL;

  dlg = (SrcQualCapsActionDlgPtr) MemNew (sizeof (SrcQualCapsActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SrcQualCapsActionToDialog;
  dlg->fromdialog = SrcQualCapsActionFromDialog;

  ValNodeAddPointer (&list, Source_qual_cell_type, StringSave (GetSourceQualName(Source_qual_cell_type)));
  ValNodeAddPointer (&list, Source_qual_dev_stage, StringSave (GetSourceQualName(Source_qual_dev_stage)));
  ValNodeAddPointer (&list, Source_qual_nat_host, StringSave (GetSourceQualName(Source_qual_nat_host)));
  ValNodeAddPointer (&list, Source_qual_isolation_source, StringSave (GetSourceQualName(Source_qual_isolation_source)));
  ValNodeAddPointer (&list, Source_qual_lab_host, StringSave (GetSourceQualName(Source_qual_lab_host)));
  ValNodeAddPointer (&list, Source_qual_sex, StringSave (GetSourceQualName(Source_qual_sex)));
  ValNodeAddPointer (&list, Source_qual_tissue_type, StringSave (GetSourceQualName(Source_qual_tissue_type)));

  dlg->field = ValNodeSelectionDialog (p, list, 4,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "source qual",
                                               change_notify, change_userdata, FALSE);

  return (DialoG) p;
}


typedef struct removexrefsaction {
  DIALOG_MESSAGE_BLOCK
  DialoG feature;
  PopuP  suppression;
  PopuP  necessary;
  DialoG constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} RemoveXrefsActionDlgData, PNTR RemoveXrefsActionDlgPtr;


static Pointer RemoveXrefsActionFromDialog (DialoG d)
{
  RemoveXrefsActionDlgPtr dlg;
  RemoveXrefsActionPtr    action;
  GeneXrefTypePtr         g;
  Int2                    val;

  dlg = (RemoveXrefsActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  g = GeneXrefTypeNew ();
  g->feature = GetFeatureTypeFromFeatureTypeDialog(dlg->feature);

  val = GetValue (dlg->necessary);
  switch (val) {
    case 1:
      g->necessary = Gene_xref_necessary_type_any;
      break;
    case 2:
      g->necessary = Gene_xref_necessary_type_necessary;
      break;
    case 3:
      g->necessary = Gene_xref_necessary_type_unnecessary;
      break;
    default:
      g->necessary = Gene_xref_necessary_type_any;
  }

  val = GetValue (dlg->suppression);
  switch (val) {
    case 1:
      g->suppression = Gene_xref_suppression_type_any;
      break;
    case 2:
      g->suppression = Gene_xref_suppression_type_suppressing;
      break;
    case 3:
      g->suppression = Gene_xref_suppression_type_non_suppressing;
      break;
    default:
      g->suppression = Gene_xref_suppression_type_any;
  }

  action = RemoveXrefsActionNew ();
  action->xref_type = ValNodeNew (NULL);
  action->xref_type->choice = XrefType_gene;
  action->xref_type->data.ptrvalue = g;
  action->constraint = DialogToPointer (dlg->constraint);

  return action;
}


static void RemoveXrefsActionToDialog (DialoG d, Pointer data)
{
  RemoveXrefsActionDlgPtr dlg;
  RemoveXrefsActionPtr    action;
  GeneXrefTypePtr         g;

  dlg = (RemoveXrefsActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  action = (RemoveXrefsActionPtr) data;
  if (action == NULL
      || action->xref_type == NULL || action->xref_type->choice != XrefType_gene
      || (g = (GeneXrefTypePtr) action->xref_type->data.ptrvalue) == NULL) {
    PointerToDialog (dlg->feature, NULL);
    SetValue (dlg->suppression, 1);
    SetValue (dlg->necessary, 1);
  } else {
    SetFeatureTypeInFeatureTypeDialog(dlg->feature, g->feature);
    switch (g->necessary) {
      case Gene_xref_necessary_type_any:
        SetValue (dlg->necessary, 1);
        break;
      case Gene_xref_necessary_type_necessary:
        SetValue (dlg->necessary, 2);
        break;
      case Gene_xref_necessary_type_unnecessary:
        SetValue (dlg->necessary, 3);
        break;
      default:
        SetValue (dlg->necessary, 1);
        break;
    }

    switch (g->suppression) {
      case Gene_xref_suppression_type_any:
        SetValue (dlg->suppression, 1);
        break;
      case Gene_xref_suppression_type_suppressing:
        SetValue (dlg->suppression, 2);
        break;
      case Gene_xref_suppression_type_non_suppressing:
        SetValue (dlg->suppression, 3);
        break;
      default:
        SetValue (dlg->suppression, 1);
        break;
    }
  }

  if (action == NULL) {
    PointerToDialog (dlg->constraint, NULL);
  } else {
    PointerToDialog (dlg->constraint, action->constraint);
  }
  
}


static void ChangeXrefsDialogPopup (PopuP p)
{
  RemoveXrefsActionDlgPtr dlg;

  dlg = (RemoveXrefsActionDlgPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->userdata);
  }
}


static DialoG RemoveXrefsDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  RemoveXrefsActionDlgPtr dlg;
  GrouP                   p, k;

  dlg = (RemoveXrefsActionDlgPtr) MemNew (sizeof (RemoveXrefsActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = RemoveXrefsActionFromDialog;
  dlg->todialog = RemoveXrefsActionToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->feature = FeatureTypeDialog(p, change_notify, change_userdata);

  k = NormalGroup (p, 0, 3, "Remove Gene Xrefs that are", programFont, NULL);
  dlg->suppression = PopupList (k, TRUE, ChangeXrefsDialogPopup);
  SetObjectExtra (dlg->suppression, dlg, NULL);
  PopupItem (dlg->suppression, "Suppressing or non-suppressing");
  PopupItem (dlg->suppression, "Suppressing");
  PopupItem (dlg->suppression, "Non-suppressing");
  SetValue (dlg->suppression, 1);

  dlg->necessary = PopupList (k, TRUE, ChangeXrefsDialogPopup);
  SetObjectExtra (dlg->necessary, dlg, NULL);
  PopupItem (dlg->necessary, "Necessary or unnecessary");
  PopupItem (dlg->necessary, "Necessary");
  PopupItem (dlg->necessary, "Unnecessary");
  SetValue (dlg->necessary, 1);

  dlg->constraint = ConstraintSetDialog (p, dlg->change_notify, dlg->change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature, (HANDLE) k, (HANDLE) dlg->constraint, NULL);

  return (DialoG) p;
}


typedef struct makegenexrefactiondialog {
  DIALOG_MESSAGE_BLOCK
  DialoG feature;
  DialoG constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} MakeGeneXrefActionDlgData, PNTR MakeGeneXrefActionDlgPtr;


static Pointer MakeGeneXrefActionFromDialog (DialoG d)
{
  MakeGeneXrefActionDlgPtr dlg;
  MakeGeneXrefActionPtr    action = NULL;
  ValNodePtr constraint;

  dlg = (MakeGeneXrefActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = MakeGeneXrefActionNew ();
  action->feature = GetFeatureTypeFromFeatureTypeDialog(dlg->feature);

  constraint = DialogToPointer (dlg->constraint);
  action->constraint = constraint;
  return action;
}


static void MakeGeneXrefActionToDialog (DialoG d, Pointer data)
{
  MakeGeneXrefActionDlgPtr dlg;
  MakeGeneXrefActionPtr    action;

  dlg = (MakeGeneXrefActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  
  if ((action = (MakeGeneXrefActionPtr) data) == NULL) {
    SetFeatureTypeInFeatureTypeDialog(dlg->feature, Macro_feature_type_any);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    SetFeatureTypeInFeatureTypeDialog(dlg->feature, action->feature);
    PointerToDialog (dlg->constraint, action->constraint);
  }
}


static DialoG MakeGeneXrefDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  MakeGeneXrefActionDlgPtr dlg;
  GrouP                    p;

  dlg = (MakeGeneXrefActionDlgPtr) MemNew (sizeof (MakeGeneXrefActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = MakeGeneXrefActionFromDialog;
  dlg->todialog = MakeGeneXrefActionToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->feature = FeatureTypeDialog(p, change_notify, change_userdata);
  dlg->constraint = ConstraintSetDialog (p, change_notify, change_userdata);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature, (HANDLE) dlg->constraint, NULL);
  return (DialoG) p;
}


typedef struct authorfixactiondialog {
  DIALOG_MESSAGE_BLOCK
  GrouP  fix_type;
  DialoG constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} AuthorFixActionDlgData, PNTR AuthorFixActionDlgPtr;


static Pointer AuthorFixActionFromDialog (DialoG d)
{
  AuthorFixActionDlgPtr dlg;
  AuthorFixActionPtr    action = NULL;
  ValNodePtr constraint;

  dlg = (AuthorFixActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = AuthorFixActionNew ();
  switch (GetValue (dlg->fix_type)) {
    case 1:
      action->fix_type = Author_fix_type_truncate_middle_initials;
      break;
    case 2:
      action->fix_type = Author_fix_type_strip_suffix;
      break;
    case 3:
      action->fix_type = Author_fix_type_move_middle_to_first;
      break;
    default:
      action->fix_type = Author_fix_type_truncate_middle_initials;
      break;
  }
  constraint = DialogToPointer (dlg->constraint);
  action->constraint = constraint;
  return action;
}


static void AuthorFixActionToDialog (DialoG d, Pointer data)
{
  AuthorFixActionDlgPtr dlg;
  AuthorFixActionPtr    action;

  dlg = (AuthorFixActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  
  if ((action = (AuthorFixActionPtr) data) == NULL) {
    SetValue (dlg->fix_type, 1);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    switch (action->fix_type) {
      case Author_fix_type_truncate_middle_initials:
        SetValue (dlg->fix_type, 1);
        break;
      case Author_fix_type_strip_suffix:
        SetValue (dlg->fix_type, 2);
        break;
      case Author_fix_type_move_middle_to_first:
        SetValue (dlg->fix_type, 3);
        break;
      default:
        SetValue (dlg->fix_type, 1);
        break;
    }
    PointerToDialog (dlg->constraint, action->constraint);
  }
}

static DialoG AuthorFixActionDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  AuthorFixActionDlgPtr dlg;
  GrouP                    p;

  dlg = (AuthorFixActionDlgPtr) MemNew (sizeof (AuthorFixActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = AuthorFixActionFromDialog;
  dlg->todialog = AuthorFixActionToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->fix_type = NormalGroup (p, 0, 3, "Author change", programFont, NULL);
  RadioButton (dlg->fix_type, "Truncate middle initials");
  RadioButton (dlg->fix_type, "Strip author suffix");
  RadioButton (dlg->fix_type, "Move middle name to first name");
  SetValue (dlg->fix_type, 1);
  dlg->constraint = ConstraintSetDialog (p, change_notify, change_userdata);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->fix_type, (HANDLE) dlg->constraint, NULL);
  return (DialoG) p;
}


typedef struct updatesequencesactiondialog {
  DIALOG_MESSAGE_BLOCK
  TexT filename;
  ButtoN add_cit_subs;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} UpdateSequencesActionDlgData, PNTR UpdateSequencesActionDlgPtr;


static Pointer UpdateSequencesActionFromDialog (DialoG d)
{
  UpdateSequencesActionDlgPtr dlg;
  UpdateSequencesActionPtr    action = NULL;

  dlg = (UpdateSequencesActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = UpdateSequencesActionNew ();
  action->filename = SaveStringFromText (dlg->filename);
  action->add_cit_subs = GetStatus (dlg->add_cit_subs);
  return action;
}


static void UpdateSequencesActionToDialog (DialoG d, Pointer data)
{
  UpdateSequencesActionDlgPtr dlg;
  UpdateSequencesActionPtr    action;

  dlg = (UpdateSequencesActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  
  if ((action = (UpdateSequencesActionPtr) data) == NULL) {
    SetTitle (dlg->filename, "");
    SetStatus (dlg->add_cit_subs, FALSE);
  } else {
    SetTitle (dlg->filename, action->filename);
    SetStatus (dlg->add_cit_subs, action->add_cit_subs);
  }
}


static void ChangeUpdateSequencesActionText (TexT t)
{
  UpdateSequencesActionDlgPtr dlg;

  dlg = (UpdateSequencesActionDlgPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ChangeUpdateSequencesActionButton (ButtoN b)
{
  UpdateSequencesActionDlgPtr dlg;

  dlg = (UpdateSequencesActionDlgPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static ValNodePtr TestUpdateSequencesActionDialog (DialoG d)
{
  UpdateSequencesActionDlgPtr dlg;
  ValNodePtr errs = NULL;

  dlg = (UpdateSequencesActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (TextHasNoText (dlg->filename)) {
      ValNodeAddPointer (&errs, 0, "Needs filename");
    }
  }
  return errs;
}


static DialoG UpdateSequencesActionDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  UpdateSequencesActionDlgPtr dlg;
  GrouP                    p, g;

  dlg = (UpdateSequencesActionDlgPtr) MemNew (sizeof (UpdateSequencesActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = UpdateSequencesActionFromDialog;
  dlg->todialog = UpdateSequencesActionToDialog;
  dlg->testdialog = TestUpdateSequencesActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "File Name", 0, dialogTextHeight, programFont, 'r');
  dlg->filename =  DialogText (g, "", 10, ChangeUpdateSequencesActionText);
  SetObjectExtra (dlg->filename, dlg, NULL);

  dlg->add_cit_subs = CheckBox (p, "Add CitSub", ChangeUpdateSequencesActionButton);
  SetObjectExtra (dlg->add_cit_subs, dlg, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlg->add_cit_subs, NULL);
  return (DialoG) p;
}


typedef struct formataction {
  DIALOG_MESSAGE_BLOCK
  Uint1 action_type;
} FormatActionDlgData, PNTR FormatActionDlgPtr;


static Pointer FormatActionFromDialog (DialoG d)
{
  ValNodePtr vnp;
  FormatActionDlgPtr dlg;

  dlg = (FormatActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  vnp = ValNodeNew (NULL);
  vnp->choice = dlg->action_type;
  return vnp;
}


static DialoG FormatActionDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata, Uint1 action_type)
{
  FormatActionDlgPtr dlg;
  GrouP              p;

  dlg = (FormatActionDlgPtr) MemNew (sizeof (FormatActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = NULL;
  dlg->fromdialog = FormatActionFromDialog;
  dlg->action_type = action_type;

  return (DialoG) p;
}


static DialoG FormatCollectionDateDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return FormatActionDialog (h, edit, change_notify, change_userdata, FixFormatAction_collection_date);
}


static DialoG FormatLatLonDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return FormatActionDialog (h, edit, change_notify, change_userdata, FixFormatAction_lat_lon);
}


static DialoG FormatPrimerDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return FormatActionDialog (h, edit, change_notify, change_userdata, FixFormatAction_primers);
}


static DialoG FormatProteinNameDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return FormatActionDialog (h, edit, change_notify, change_userdata, FixFormatAction_protein_name);
}



typedef struct removeduplicatefeataction {
  DIALOG_MESSAGE_BLOCK
  DialoG feature_type;
  ButtoN ignore_partials;
  ButtoN case_sensitive;
  ButtoN remove_proteins;
  DialoG constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} RemoveDuplicateFeatActionDlgData, PNTR RemoveDuplicateFeatActionDlgPtr;


static Pointer RemoveDuplicateFeatActionFromDialog (DialoG d)
{
  RemoveDuplicateFeatureActionPtr action;
  RemoveDuplicateFeatActionDlgPtr dlg;
  ValNodePtr vnp;

  dlg = (RemoveDuplicateFeatActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = RemoveDuplicateFeatureActionNew ();
  vnp = DialogToPointer (dlg->feature_type);
  if (vnp == NULL) {
    action->type = Macro_feature_type_any;
  } else {
    action->type = vnp->choice;
  }
  vnp = ValNodeFree (vnp);
  action->ignore_partials = GetStatus (dlg->ignore_partials);
  action->remove_proteins = GetStatus (dlg->remove_proteins);
  action->case_sensitive = GetStatus (dlg->case_sensitive);
  action->rd_constraint = DialogToPointer (dlg->constraint);
  return action;
}


static void RemoveDuplicateFeatActionToDialog (DialoG d, Pointer data)
{
  RemoveDuplicateFeatureActionPtr action;
  RemoveDuplicateFeatActionDlgPtr dlg;
  ValNode vn;

  dlg = (RemoveDuplicateFeatActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  action = (RemoveDuplicateFeatureActionPtr) data;
  if (data == NULL) {
    PointerToDialog (dlg->feature_type, NULL);
    SetStatus (dlg->ignore_partials, FALSE);
    SetStatus (dlg->remove_proteins, TRUE);
    SetStatus (dlg->case_sensitive, TRUE);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    MemSet (&vn, 0, sizeof (ValNode));
    vn.choice = (Uint1)action->type;
    PointerToDialog (dlg->feature_type, &vn);
    SetStatus (dlg->ignore_partials, action->ignore_partials);
    SetStatus (dlg->remove_proteins, action->remove_proteins);
    SetStatus (dlg->case_sensitive, action->case_sensitive);
    PointerToDialog (dlg->constraint, action->rd_constraint);
  }
}


static void ChangeRemoveDuplicateFeatActionButon (ButtoN b)
{
  RemoveDuplicateFeatActionDlgPtr dlg;

  dlg = (RemoveDuplicateFeatActionDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


NLM_EXTERN DialoG RemoveDuplicateFeatActionDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  RemoveDuplicateFeatActionDlgPtr dlg;
  GrouP              p;
  ValNodePtr         feature_list = NULL;

  dlg = (RemoveDuplicateFeatActionDlgPtr) MemNew (sizeof (RemoveDuplicateFeatActionDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = RemoveDuplicateFeatActionToDialog;
  dlg->fromdialog = RemoveDuplicateFeatActionFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  feature_list = ValNodeNew (NULL);
  feature_list->choice = Macro_feature_type_any;
  feature_list->data.ptrvalue = StringSave ("Any");
  AddAllFeaturesToChoiceList (&feature_list);
  dlg->feature_type = ValNodeSelectionDialog (p, feature_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                change_notify, change_userdata, FALSE);
  dlg->ignore_partials = CheckBox (p, "Ignore partials", ChangeRemoveDuplicateFeatActionButon);
  SetObjectExtra (dlg->ignore_partials, dlg, NULL);
  dlg->case_sensitive = CheckBox (p, "Case sensitive", ChangeRemoveDuplicateFeatActionButon);
  SetObjectExtra (dlg->case_sensitive, dlg, NULL);
  SetStatus (dlg->case_sensitive, TRUE);

  dlg->remove_proteins = CheckBox (p, "Remove proteins", ChangeRemoveDuplicateFeatActionButon);
  SetObjectExtra (dlg->remove_proteins, dlg, NULL);
  SetStatus (dlg->remove_proteins, TRUE);

  dlg->constraint = ComplexConstraintDialog (p, change_notify, change_userdata);
  ChangeComplexConstraintFieldType (dlg->constraint, FieldType_feature_field, NULL, Macro_feature_type_any);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_type, 
                              (HANDLE) dlg->ignore_partials, 
                              (HANDLE) dlg->case_sensitive, 
                              (HANDLE) dlg->remove_proteins, 
                              (HANDLE) dlg->constraint, 
                              NULL);

  return (DialoG) p;
}


static DialoG PubFieldDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)

{  
  return ValNodeSelectionDialog (h, GetPubFieldList(), TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "pub field", 
                                change_notify, change_userdata, FALSE);
}




typedef struct PublicationConstraintdlg {
  DIALOG_MESSAGE_BLOCK
  PopuP pub_type;
  DialoG field_type;
  DialoG string_constraint;

  DialoG special_field;
  PopuP  special_field_type;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} PublicationConstraintDialogData, PNTR PublicationConstraintDialogPtr;


static void ClearPublicationConstraintDialogText (DialoG d)
{
  PublicationConstraintDialogPtr dlg;

  dlg = (PublicationConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  ClearStringConstraintDialogText (dlg->string_constraint);
}


static Pointer DialogToPublicationConstraint (DialoG d) 
{
  PublicationConstraintDialogPtr dlg;
  StringConstraintPtr    string_constraint = NULL;
  PublicationConstraintPtr       constraint = NULL;
  PubFieldConstraintPtr  p;
  Int4                   val;
  ValNodePtr             vnp;

  dlg = (PublicationConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  constraint = PublicationConstraintNew ();
  val = GetValue (dlg->pub_type);
  switch (val) {
    case 2:
      constraint->type = Pub_type_published;
      break;
    case 3:
      constraint->type = Pub_type_unpublished;
      break;
    case 4:
      constraint->type = Pub_type_in_press;
      break;
    case 5:
      constraint->type = Pub_type_submitter_block;
      break;
    default:
      constraint->type = Pub_type_any;
      break;
  }
  string_constraint = DialogToPointer (dlg->string_constraint);
  vnp = DialogToPointer (dlg->field_type);
  if (vnp != NULL && !IsStringConstraintEmpty (string_constraint)) {
    p = PubFieldConstraintNew ();
    p->field = vnp->choice;
    p->constraint = string_constraint;
    string_constraint = NULL;
    constraint->field = p;
  }
  vnp = ValNodeFree (vnp);
  string_constraint = StringConstraintFree (string_constraint);

  vnp = DialogToPointer (dlg->special_field);
  if (vnp != NULL) {
    val = GetValue (dlg->special_field_type);
    if (val > 1 && val <= k_NumSpecialPubFieldWords + 1) {
      constraint->special_field = PubFieldSpecialConstraintNew ();
      constraint->special_field->field = vnp->choice;
      constraint->special_field->constraint = ValNodeNew (NULL);
      switch (val) {
        case 2:
          constraint->special_field->constraint->choice = PubFieldSpecialConstraintType_is_present;
          break;
        case 3:
          constraint->special_field->constraint->choice = PubFieldSpecialConstraintType_is_not_present;
          break;
        case 4:
          constraint->special_field->constraint->choice = PubFieldSpecialConstraintType_is_all_caps;
          break;
        case 5:
          constraint->special_field->constraint->choice = PubFieldSpecialConstraintType_is_all_lower;
          break;
        case 6:
          constraint->special_field->constraint->choice = PubFieldSpecialConstraintType_is_all_punct;
          break;
      }
    }
    vnp = ValNodeFree (vnp);
  }

  return (Pointer) constraint;
}


static void PublicationConstraintToDialog (DialoG d, Pointer data)
{
  PublicationConstraintDialogPtr dlg;
  ValNode vn; 
  PublicationConstraintPtr constraint;

  dlg = (PublicationConstraintDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  constraint = (PublicationConstraintPtr) data;
  if (constraint == NULL) {
    SetValue (dlg->pub_type, 1);
    PointerToDialog (dlg->field_type, NULL);
    PointerToDialog (dlg->string_constraint, NULL);
    PointerToDialog (dlg->special_field, NULL);
    SetValue (dlg->special_field_type, 1);
  } else {
    switch (constraint->type) {
      case Pub_type_published:
        SetValue (dlg->pub_type, 2);
        break;
      case Pub_type_unpublished:
        SetValue (dlg->pub_type, 3);
        break;
      case Pub_type_in_press:
        SetValue (dlg->pub_type, 4);
        break;
      case Pub_type_submitter_block:
        SetValue (dlg->pub_type, 5);
        break;
      case Pub_type_any:
      default:
        SetValue (dlg->pub_type, 1);
        break;
    }
    if (constraint->field == NULL) {
      PointerToDialog (dlg->field_type, NULL);
      PointerToDialog (dlg->string_constraint, NULL);
    } else {
      vn.choice = (Uint1)constraint->field->field;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->field_type, &vn);
      PointerToDialog (dlg->string_constraint, constraint->field->constraint);
    }
    if (constraint->special_field == NULL) {
      PointerToDialog (dlg->special_field, NULL);
      SetValue (dlg->special_field_type, 1);
    } else {
      vn.choice = (Uint1)(constraint->special_field->field);
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->special_field, &vn);
      if (constraint->special_field->constraint == NULL) {
        SetValue (dlg->special_field_type, 1);
      } else {
        switch (constraint->special_field->constraint->choice) {
          case PubFieldSpecialConstraintType_is_present:
            SetValue (dlg->special_field_type, 2);
            break;
          case PubFieldSpecialConstraintType_is_not_present:
            SetValue (dlg->special_field_type, 3);
            break;
          case PubFieldSpecialConstraintType_is_all_caps:
            SetValue (dlg->special_field_type, 4);
            break;
          case PubFieldSpecialConstraintType_is_all_lower:
            SetValue (dlg->special_field_type, 5);
            break;
          case PubFieldSpecialConstraintType_is_all_punct:
            SetValue (dlg->special_field_type, 6);
            break;
          default:
            SetValue (dlg->special_field_type, 1);
            break;
        }
      }
    }
  } 
}

static ValNodePtr TestPublicationConstraintDialog (DialoG d)
{
  PublicationConstraintPtr p;
  ValNodePtr err_list = NULL;

  p = DialogToPointer (d);
  if (IsPublicationConstraintEmpty (p)) {
    ValNodeAddPointer (&err_list, 0, "No pub constraint");
  }
  p = PublicationConstraintFree (p);
  return err_list;
}


static void ChangePublicationConstraintPopup (PopuP p)
{
  PublicationConstraintDialogPtr dlg;

  dlg = (PublicationConstraintDialogPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
} 


static DialoG PublicationConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  PublicationConstraintDialogPtr dlg;
  GrouP      p, g1, g2;
  Int4       i;
  
  dlg = (PublicationConstraintDialogPtr) MemNew (sizeof (PublicationConstraintDialogData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = PublicationConstraintToDialog;
  dlg->fromdialog = DialogToPublicationConstraint;
  dlg->testdialog = TestPublicationConstraintDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "Publication Status", 0, dialogTextHeight, programFont, 'r');
  dlg->pub_type = PopupList (g1, TRUE, ChangePublicationConstraintPopup);
  SetObjectExtra (dlg->pub_type, dlg, NULL);
  PopupItem (dlg->pub_type, "Any");
  PopupItem (dlg->pub_type, "Published");
  PopupItem (dlg->pub_type, "Unpublished");
  PopupItem (dlg->pub_type, "In-press");
  PopupItem (dlg->pub_type, "Submitter block");
  SetValue (dlg->pub_type, 1);

  g2 = HiddenGroup (p, 2, 0, NULL);
  dlg->field_type = PubFieldDialog (g2, change_notify, change_userdata);
  dlg->string_constraint = StringConstraintDialog (g2, NULL, FALSE, change_notify, change_userdata);

  dlg->special_field = PubFieldDialog (g2, change_notify, change_userdata);
  dlg->special_field_type = PopupList (g2, TRUE, ChangePublicationConstraintPopup);
  SetObjectExtra (dlg->special_field_type, dlg, NULL);
  PopupItem (dlg->special_field_type, "Any");
  for (i = 0; i < k_NumSpecialPubFieldWords; i++) {
    PopupItem (dlg->special_field_type, s_SpecialPubFieldWords[i]);
  }
  SetValue (dlg->special_field_type, 1);

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);

  return (DialoG) p;
}



typedef struct quantityconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  PopuP quantity_type;
  TexT  quantity;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} QuantityConstraintDlgDaga, PNTR QuantityConstraintDlgPtr;


static void QuantityConstraintToDialog (DialoG d, Pointer data)
{
  QuantityConstraintDlgPtr dlg;
  ValNodePtr               vnp;
  Char                     buf[15];

  dlg = (QuantityConstraintDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) {
    return;
  }

  vnp = (ValNodePtr) data;
  if (vnp == NULL) {
    SetValue (dlg->quantity_type, 1);
    Hide (dlg->quantity);
  } else if (vnp->choice == QuantityConstraint_equals) {
    SetValue (dlg->quantity_type, 2);
    sprintf (buf, "%d", vnp->data.intvalue);
    SetTitle (dlg->quantity, buf);
    Show (dlg->quantity);
  } else if (vnp->choice == QuantityConstraint_greater_than) {
    SetValue (dlg->quantity_type, 3);
    sprintf (buf, "%d", vnp->data.intvalue);
    SetTitle (dlg->quantity, buf);
    Show (dlg->quantity);
  } else if (vnp->choice == QuantityConstraint_less_than) {
    SetValue (dlg->quantity_type, 4);
    sprintf (buf, "%d", vnp->data.intvalue);
    SetTitle (dlg->quantity, buf);
    Show (dlg->quantity);
  } else {
    SetValue (dlg->quantity_type, 1);
    Hide (dlg->quantity);
  }
}


static Pointer DialogToQuantityConstraint (DialoG d)
{
  QuantityConstraintDlgPtr dlg;
  ValNodePtr               vnp = NULL;
  Int2                     val;
  CharPtr                  num_text;

  dlg = (QuantityConstraintDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) {
    return NULL;
  }

  if (!TextHasNoText (dlg->quantity)) {
    num_text = SaveStringFromText (dlg->quantity);
    if (IsAllDigits(num_text)) {
      val = GetValue (dlg->quantity_type);
      switch (val) {
        case 2:
          vnp = ValNodeNew (NULL);
          vnp->choice = QuantityConstraint_equals;
          vnp->data.intvalue = atoi (num_text);
          break;
        case 3:
          vnp = ValNodeNew (NULL);
          vnp->choice = QuantityConstraint_greater_than;
          vnp->data.intvalue = atoi (num_text);
          break;
        case 4:
          vnp = ValNodeNew (NULL);
          vnp->choice = QuantityConstraint_less_than;
          vnp->data.intvalue = atoi (num_text);
          break;
      }
    }
    num_text = MemFree (num_text);
  }
  return vnp;
}


static ValNodePtr TestQuantityConstraintDialog (DialoG d)
{
  QuantityConstraintDlgPtr dlg;
  Int2                     val;
  CharPtr                  num_text;
  ValNodePtr err_list = NULL;

  dlg = (QuantityConstraintDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    val = GetValue (dlg->quantity_type);
    if (val > 1 && val <= k_NumQuantityWords + 1) {
      if (TextHasNoText (dlg->quantity)) {
        ValNodeAddPointer (&err_list, 0, "missing value");
      } else {
        num_text = SaveStringFromText (dlg->quantity);
        if (!IsAllDigits (num_text)) {
          ValNodeAddPointer (&err_list, 0, "bad value");
        }
        num_text = MemFree (num_text);
      }
    }
  }
  return err_list;
}


static void ChangeQuantityConstraintQuantityType (PopuP p)
{
  QuantityConstraintDlgPtr dlg;
  Int2 val;

  dlg = (QuantityConstraintDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }
  val = GetValue (dlg->quantity_type);

  if (val < 2 || val > k_NumQuantityWords + 1) {
    Hide (dlg->quantity);
  } else {
    Show (dlg->quantity);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ChangeQuantityConstraintQuantity (TexT t)
{
  QuantityConstraintDlgPtr dlg;

  dlg = (QuantityConstraintDlgPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static DialoG QuantityConstraintDialog (GrouP h, CharPtr title, CharPtr default_val, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  QuantityConstraintDlgPtr dlg;
  GrouP                    p;
  Int4                     i;
  
  dlg = (QuantityConstraintDlgPtr) MemNew (sizeof (QuantityConstraintDlgDaga));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = QuantityConstraintToDialog;
  dlg->fromdialog = DialogToQuantityConstraint;
  dlg->testdialog = TestQuantityConstraintDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  StaticPrompt (p, title, 0, dialogTextHeight, programFont, 'r');
  dlg->quantity_type = PopupList (p, TRUE, ChangeQuantityConstraintQuantityType);
  SetObjectExtra (dlg->quantity_type, dlg, NULL);
  PopupItem (dlg->quantity_type, "Any");
  for (i = 0; i < k_NumQuantityWords; i++) {
    PopupItem (dlg->quantity_type, s_QuantityWords[i]);
  }
  SetValue (dlg->quantity_type, 1);
  dlg->quantity = DialogText (p, default_val, 5, ChangeQuantityConstraintQuantity);
  SetObjectExtra (dlg->quantity, dlg, NULL);
  Hide (dlg->quantity);

  return (DialoG) p;
}


typedef struct sequenceconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  GrouP seqtype;
  GrouP rna_subtype_grp;
  DialoG rna_subtype;
  ButtoN feat_present;
  DialoG feature_type;
  DialoG num_type_features;
  DialoG num_features;
  DialoG length;
  DialoG id;
  PopuP  strandedness;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} SequenceConstraintDlgData, PNTR SequenceConstraintDlgPtr;


static void ChangeSequenceConstraintButton (ButtoN b)
{
  SequenceConstraintDlgPtr dlg;

  dlg = (SequenceConstraintDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;

  if (GetStatus (dlg->feat_present)) {
    Enable (dlg->feature_type);
  } else {
    Disable (dlg->feature_type);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeSequenceConstraintGroup (GrouP g)
{
  SequenceConstraintDlgPtr dlg;

  dlg = (SequenceConstraintDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;

  if (GetValue (dlg->seqtype) == 4) {
    Show (dlg->rna_subtype_grp);
  } else {
    Hide (dlg->rna_subtype_grp);
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeSequenceConstraintPopup (PopuP p)
{
  SequenceConstraintDlgPtr dlg;

  dlg = (SequenceConstraintDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SequenceConstraintToDialog (DialoG d, Pointer data)
{
  SequenceConstraintDlgPtr dlg;
  SequenceConstraintPtr constraint;
  ValNode vn;

  dlg = (SequenceConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  constraint = (SequenceConstraintPtr) data;
  if (constraint == NULL) {
    SetValue (dlg->seqtype, 2);
    PointerToDialog (dlg->rna_subtype, NULL);
    Hide (dlg->rna_subtype_grp);
    SetStatus (dlg->feat_present, FALSE);
    PointerToDialog (dlg->feature_type, NULL);
    PointerToDialog (dlg->num_type_features, NULL);
    PointerToDialog (dlg->id, NULL);
    PointerToDialog (dlg->num_features, NULL);
    PointerToDialog (dlg->length, NULL);
    SetValue (dlg->strandedness, 1);
  } else {
    if (constraint->seqtype == NULL) {
      SetValue (dlg->seqtype, 2);
      PointerToDialog (dlg->rna_subtype, NULL);
      Hide (dlg->rna_subtype_grp);
    } else {
      switch (constraint->seqtype->choice) {
        case SequenceConstraintMolTypeConstraint_any :
          SetValue (dlg->seqtype, 1);
          PointerToDialog (dlg->rna_subtype, NULL);
          Hide (dlg->rna_subtype_grp);
          break;
        case SequenceConstraintMolTypeConstraint_nucleotide :
          SetValue (dlg->seqtype, 2);
          PointerToDialog (dlg->rna_subtype, NULL);
          Hide (dlg->rna_subtype_grp);
          break;
        case SequenceConstraintMolTypeConstraint_dna :
          SetValue (dlg->seqtype, 3);
          PointerToDialog (dlg->rna_subtype, NULL);
          Hide (dlg->rna_subtype_grp);
          break;
        case SequenceConstraintMolTypeConstraint_rna :
          SetValue (dlg->seqtype, 4);
          vn.choice = constraint->seqtype->data.intvalue;
          vn.data.ptrvalue = NULL;
          vn.next = NULL;
          PointerToDialog (dlg->rna_subtype, &vn);
          Show (dlg->rna_subtype_grp);
          break;
        case SequenceConstraintMolTypeConstraint_protein :
          SetValue (dlg->seqtype, 5);
          PointerToDialog (dlg->rna_subtype, NULL);
          Hide (dlg->rna_subtype_grp);
          break;
      }
    }
  
    if (constraint->feature == Macro_feature_type_any) {
      SetStatus (dlg->feat_present, FALSE);
      PointerToDialog (dlg->feature_type, NULL);
    } else {
      SetStatus (dlg->feat_present, TRUE);
      vn.choice = (Uint1)constraint->feature;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->feature_type, &vn);
    }
    PointerToDialog (dlg->num_type_features, constraint->num_type_features);

    PointerToDialog (dlg->id, constraint->id);
    PointerToDialog (dlg->num_features, constraint->num_features);
    PointerToDialog (dlg->length, constraint->length);

    switch (constraint->strandedness) {
      case Feature_strandedness_constraint_any:
        SetValue (dlg->strandedness, 1);
        break;
      case Feature_strandedness_constraint_minus_only:
        SetValue (dlg->strandedness, 2);
        break;
      case Feature_strandedness_constraint_plus_only:
        SetValue (dlg->strandedness, 3);
        break;
      case Feature_strandedness_constraint_at_least_one_minus:
        SetValue (dlg->strandedness, 4);
        break;
      case Feature_strandedness_constraint_at_least_one_plus:
        SetValue (dlg->strandedness, 5);
        break;
      case Feature_strandedness_constraint_no_minus:
        SetValue (dlg->strandedness, 6);
        break;
      case Feature_strandedness_constraint_no_plus:
        SetValue (dlg->strandedness, 7);
        break;
    }
  }
  ChangeSequenceConstraintButton (dlg->feat_present);
}


static Pointer SequenceConstraintFromDialog (DialoG d)
{
  SequenceConstraintDlgPtr dlg;
  SequenceConstraintPtr constraint;
  ValNodePtr vnp;
  Int4 val;

  dlg = (SequenceConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  constraint = SequenceConstraintNew();

  val = GetValue (dlg->seqtype);
  switch (val) {
    case 1:
      /* don't bother to fill in, it's optional */
      break;
    case 2:
      constraint->seqtype = ValNodeNew (NULL);
      constraint->seqtype->choice = SequenceConstraintMolTypeConstraint_nucleotide;
      break;
    case 3:
      constraint->seqtype = ValNodeNew (NULL);
      constraint->seqtype->choice = SequenceConstraintMolTypeConstraint_dna;
      break;
    case 4:
      constraint->seqtype = ValNodeNew (NULL);
      constraint->seqtype->choice = SequenceConstraintMolTypeConstraint_rna;
      vnp = DialogToPointer (dlg->rna_subtype);
      if (vnp == NULL) {
        constraint->seqtype->data.intvalue = Sequence_constraint_rnamol_any;
      } else {
        constraint->seqtype->data.intvalue = vnp->choice;
      }
      break;
    case 5:
      constraint->seqtype = ValNodeNew (NULL);
      constraint->seqtype->choice = SequenceConstraintMolTypeConstraint_protein;
      break;
  }

  if (GetStatus (dlg->feat_present)) {
    vnp = DialogToPointer (dlg->feature_type);
    if (vnp == NULL) {
      constraint->feature = Macro_feature_type_any;
    } else {
      constraint->feature = vnp->choice;
    }
    vnp = ValNodeFree (vnp);
    constraint->num_type_features = DialogToPointer (dlg->num_type_features);
  } else {
    constraint->feature = Macro_feature_type_any;
  }

  constraint->id = DialogToPointer (dlg->id);
  if (IsStringConstraintEmpty (constraint->id)) {
    constraint->id = StringConstraintFree (constraint->id);
  }
  constraint->num_features = DialogToPointer (dlg->num_features);
  constraint->length = DialogToPointer (dlg->length);
  val = GetValue (dlg->strandedness);
  switch (val) {
    case 1:
      constraint->strandedness = Feature_strandedness_constraint_any;
      break;
    case 2:
      constraint->strandedness = Feature_strandedness_constraint_minus_only;
      break;
    case 3:
      constraint->strandedness = Feature_strandedness_constraint_plus_only;
      break;
    case 4:
      constraint->strandedness = Feature_strandedness_constraint_at_least_one_minus;
      break;
    case 5:
      constraint->strandedness = Feature_strandedness_constraint_at_least_one_plus;
      break;
    case 6:
      constraint->strandedness = Feature_strandedness_constraint_no_minus;
      break;
    case 7:
      constraint->strandedness = Feature_strandedness_constraint_no_plus;
      break;
  }
  return constraint;
}


NLM_EXTERN DialoG FeatureTypeDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ValNodePtr feature_list = NULL;

  AddAllFeaturesToChoiceList (&feature_list);
  return ValNodeSelectionDialog (h, feature_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                change_notify, change_userdata, FALSE);
}


NLM_EXTERN Uint2 GetFeatureTypeFromFeatureTypeDialog (DialoG d)
{
  ValNodePtr vnp;
  Uint2 feature;

  vnp = DialogToPointer (d);
  if (vnp == NULL) {
    feature = Macro_feature_type_any;
  } else {
    feature = vnp->choice;
  }
  vnp = ValNodeFree (vnp);
  return feature;
}


NLM_EXTERN void SetFeatureTypeInFeatureTypeDialog (DialoG d, Uint2 feature)
{
  ValNode vn;

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = (Uint1)feature;
  PointerToDialog (d, &vn);
}


NLM_EXTERN DialoG FeatureTypeDialogMulti (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ValNodePtr feature_list = NULL;

  ValNodeAddPointer (&feature_list, Macro_feature_type_any, StringSave ("Any"));
  AddAllFeaturesToChoiceList (&feature_list);
  return ValNodeSelectionDialog (h, feature_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                change_notify, change_userdata, FALSE);
}


static ValNodePtr TestSequenceConstraintDialog (DialoG d)
{
  SequenceConstraintDlgPtr dlg;
  ValNodePtr err_list = NULL;
  SequenceConstraintPtr constraint;

  dlg = (SequenceConstraintDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    constraint = DialogToPointer (d);
    if (IsSequenceConstraintEmpty (constraint)) {
      ValNodeAddPointer (&err_list, 0, "empty constraint");
    }
    ValNodeLink (&err_list, TestDialog (dlg->num_features));
    ValNodeLink (&err_list, TestDialog (dlg->length));
    constraint = SequenceConstraintFree (constraint);
  }
  return err_list;
}
   

NLM_EXTERN DialoG SequenceConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  SequenceConstraintDlgPtr dlg;
  GrouP p, g, grp1, grp2;
  ValNodePtr rna_subtypes = NULL;

  dlg = (SequenceConstraintDlgPtr) MemNew (sizeof (SequenceConstraintDlgData));
  p = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = SequenceConstraintToDialog;
  dlg->fromdialog = SequenceConstraintFromDialog;
  dlg->testdialog = TestSequenceConstraintDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  grp1 = HiddenGroup (p, -1, 0, NULL);
  dlg->seqtype = HiddenGroup (grp1, 5, 0, ChangeSequenceConstraintGroup);
  SetObjectExtra (dlg->seqtype, dlg, NULL);
  RadioButton (dlg->seqtype, "Any sequence");
  RadioButton (dlg->seqtype, "Nucleotides");
  RadioButton (dlg->seqtype, "DNA");
  RadioButton (dlg->seqtype, "RNA");
  RadioButton (dlg->seqtype, "Proteins");
  SetValue (dlg->seqtype, 2);
  dlg->rna_subtype_grp = HiddenGroup (grp1, 2, 0, NULL);
  StaticPrompt (dlg->rna_subtype_grp, "RNA Type", 0, dialogTextHeight, programFont, 'r');
  AddAllRNASubtypesToChoiceList (&rna_subtypes);
  dlg->rna_subtype = ValNodeSelectionDialog (dlg->rna_subtype_grp, rna_subtypes, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "rna subtype",
                                change_notify, change_userdata, FALSE);

  Hide (dlg->rna_subtype_grp);

  g = HiddenGroup (grp1, 3, 0, NULL);
  dlg->feat_present = CheckBox (g, "Feature type", ChangeSequenceConstraintButton);
  SetObjectExtra (dlg->feat_present, dlg, NULL);
  dlg->feature_type = FeatureTypeDialog (g, change_notify, change_userdata);
  dlg->num_type_features = QuantityConstraintDialog (g, "Number present", "1", change_notify, change_userdata);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->seqtype, 
                              (HANDLE) dlg->rna_subtype_grp, 
                              (HANDLE) g, 
                               NULL);

  StaticPrompt (g, "Strandedness of features", 0, popupMenuHeight, programFont, 'r');
  dlg->strandedness = PopupList (g, TRUE, ChangeSequenceConstraintPopup);
  SetObjectExtra (dlg->strandedness, dlg, NULL);
  PopupItem (dlg->strandedness, "Any");
  PopupItem (dlg->strandedness, SummarizeFeatureStrandedness(Feature_strandedness_constraint_minus_only));
  PopupItem (dlg->strandedness, SummarizeFeatureStrandedness(Feature_strandedness_constraint_plus_only));
  PopupItem (dlg->strandedness, SummarizeFeatureStrandedness(Feature_strandedness_constraint_at_least_one_minus));
  PopupItem (dlg->strandedness, SummarizeFeatureStrandedness(Feature_strandedness_constraint_at_least_one_plus));
  PopupItem (dlg->strandedness, SummarizeFeatureStrandedness(Feature_strandedness_constraint_no_minus));
  PopupItem (dlg->strandedness, SummarizeFeatureStrandedness(Feature_strandedness_constraint_no_plus));

  grp2 = HiddenGroup (p, -1, 0, NULL);
  dlg->id = StringConstraintDialog (grp2, "Where sequence ID", FALSE, change_notify, change_userdata);
  dlg->num_features = QuantityConstraintDialog (grp2, "Number of features present (any type)", "1", change_notify, change_userdata);
  dlg->length = QuantityConstraintDialog (grp2, "Length of sequence", "200", change_notify, change_userdata);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->id, 
                              (HANDLE) dlg->num_features,
                              (HANDLE) dlg->length,
                               NULL);


  return (DialoG) p;
}


typedef struct molinfoconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG field;
  ButtoN is_not;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} MolinfoConstraintDlgDaga, PNTR MolinfoConstraintDlgPtr;


static void MolinfoConstraintToDialog (DialoG d, Pointer data)
{
  MolinfoConstraintDlgPtr   dlg;
  MolinfoFieldConstraintPtr constraint;

  dlg = (MolinfoConstraintDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) {
    return;
  }

  constraint = (MolinfoFieldConstraintPtr) data;
  if (constraint == NULL) {
    PointerToDialog (dlg->field, NULL);
    SetStatus (dlg->is_not, FALSE);
  } else {
    PointerToDialog (dlg->field, constraint->field);
    SetStatus (dlg->is_not, constraint->is_not);
  }
}


static Pointer DialogToMolinfoConstraint (DialoG d)
{
  MolinfoConstraintDlgPtr   dlg;
  MolinfoFieldConstraintPtr constraint = NULL;

  dlg = (MolinfoConstraintDlgPtr) GetObjectExtra (d);

  if (dlg == NULL) {
    return NULL;
  }

  constraint = MolinfoFieldConstraintNew();
  constraint->field = DialogToPointer (dlg->field);
  constraint->is_not = GetStatus (dlg->is_not);
  return constraint;
}


static ValNodePtr TestMolinfoConstraintDialog (DialoG d)
{
  MolinfoConstraintDlgPtr   dlg;
  ValNodePtr field;
  ValNodePtr err_list = NULL;

  dlg = (MolinfoConstraintDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    field = DialogToPointer (dlg->field);
    if (field == NULL) {
      ValNodeAddPointer (&err_list, 0, "missing field");
    }
    field = MolinfoFieldFree (field);
  }
  return err_list;
}


static void ChangeMolinfoConstraintIsNot (ButtoN b)
{
  MolinfoConstraintDlgPtr dlg;

  dlg = (MolinfoConstraintDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static DialoG MolinfoConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  MolinfoConstraintDlgPtr dlg;
  GrouP                    p;
  
  dlg = (MolinfoConstraintDlgPtr) MemNew (sizeof (MolinfoConstraintDlgDaga));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = MolinfoConstraintToDialog;
  dlg->fromdialog = DialogToMolinfoConstraint;
  dlg->testdialog = TestMolinfoConstraintDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->field = SequenceQualDialog (p, change_notify, change_userdata);
  dlg->is_not = CheckBox (p, "Is not", ChangeMolinfoConstraintIsNot);
  SetObjectExtra (dlg->is_not, dlg, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->field, (HANDLE) dlg->is_not, NULL);

  return (DialoG) p;
}


typedef struct translationconstraintdlg {
  DIALOG_MESSAGE_BLOCK

  DialoG actual_strings;
  DialoG transl_strings;
  PopuP  internal_stops;
  DialoG num_mismatches;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} TranslationConstraintDlgData, PNTR TranslationConstraintDlgPtr;


static void TranslationConstraintToDialog(DialoG d, Pointer data)
{
  TranslationConstraintDlgPtr dlg;
  TranslationConstraintPtr    constraint;

  dlg = (TranslationConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  constraint = (TranslationConstraintPtr) data;
  if (constraint == NULL) {
    PointerToDialog (dlg->actual_strings, NULL);
    PointerToDialog (dlg->transl_strings, NULL);
    PointerToDialog (dlg->num_mismatches, NULL);
    SetValue (dlg->internal_stops, 1);
  } else {
    PointerToDialog (dlg->actual_strings, constraint->actual_strings);
    PointerToDialog (dlg->transl_strings, constraint->transl_strings);
    PointerToDialog (dlg->num_mismatches, constraint->num_mismatches);
    switch (constraint->internal_stops) {
      case Match_type_constraint_dont_care:
        SetValue (dlg->internal_stops, 1);
        break;
      case Match_type_constraint_yes:
        SetValue (dlg->internal_stops, 2);
        break;
      case Match_type_constraint_no:
        SetValue (dlg->internal_stops, 3);
        break;
      default:
        SetValue (dlg->internal_stops, 1);
        break;
    }
  }
}


static Pointer DialogToTranslationConstraint (DialoG d)
{
  TranslationConstraintDlgPtr dlg;
  TranslationConstraintPtr    constraint;
  Int2 val;

  dlg = (TranslationConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  constraint = TranslationConstraintNew();
  constraint->actual_strings = DialogToPointer (dlg->actual_strings);
  constraint->transl_strings = DialogToPointer (dlg->transl_strings);
  constraint->num_mismatches = DialogToPointer (dlg->num_mismatches);
  val = GetValue (dlg->internal_stops);
  switch (val) {
    case 1:
      constraint->internal_stops = Match_type_constraint_dont_care;
      break;
    case 2:
      constraint->internal_stops = Match_type_constraint_yes;
      break;
    case 3:
      constraint->internal_stops = Match_type_constraint_no;
      break;
    default:
      constraint->internal_stops = Match_type_constraint_dont_care;
      break;
  }

  if (IsTranslationConstraintEmpty(constraint)) {
    constraint = TranslationConstraintFree (constraint);
  }
  return constraint;
}


static ValNodePtr TestTranslationConstraintDialog (DialoG d)
{
  TranslationConstraintPtr    constraint;
  ValNodePtr err_list = NULL;

  constraint = DialogToPointer (d);
  if (constraint == NULL) {
    ValNodeAddPointer (&err_list, 0, "constraint is empty");
  } else {
    constraint = TranslationConstraintFree (constraint);
  }
  return err_list;
}


static void ChangeTranslationConstraintPopup (PopuP p)
{
  TranslationConstraintDlgPtr dlg;

  dlg = (TranslationConstraintDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static DialoG TranslationConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TranslationConstraintDlgPtr dlg;
  GrouP                       p, g;
  
  dlg = (TranslationConstraintDlgPtr) MemNew (sizeof (TranslationConstraintDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = TranslationConstraintToDialog;
  dlg->fromdialog = DialogToTranslationConstraint;
  dlg->testdialog = TestTranslationConstraintDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->actual_strings = StringConstraintDialog (p, "Where protein sequence", TRUE, change_notify, change_userdata);
  dlg->transl_strings = StringConstraintDialog (p, "Where translation", TRUE, change_notify, change_userdata);

  g = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Internal Stops in Translation", 0, dialogTextHeight, programFont, 'r');
  dlg->internal_stops = PopupList (g, TRUE, ChangeTranslationConstraintPopup);
  SetObjectExtra (dlg->internal_stops, dlg, NULL);
  PopupItem (dlg->internal_stops, "Don't care");
  PopupItem (dlg->internal_stops, "Yes");
  PopupItem (dlg->internal_stops, "No");
  SetValue (dlg->internal_stops, 1);

  dlg->num_mismatches = QuantityConstraintDialog (p, "Number of mismatches between protein sequence and translation", "1", change_notify, change_userdata);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->actual_strings, (HANDLE) dlg->transl_strings, (HANDLE) g, (HANDLE) dlg->internal_stops, NULL);

  return (DialoG) p;
}


static DialoG RNAFieldConstraintDialog (GrouP h, CharPtr type_label, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
static DialoG FeatureFieldConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
static DialoG MiscFieldConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);


/* Notes:
 * There can be more than one constraint, and the user will need to be able to view, edit, add, and delete them.
 * Should have a DOC that lists the constraints, user can click on buttons to edit or delete, another button to add.
 * Need to have an Edit window to edit individual constraints.
 * Let user toggle between types of constraints only when creating a new one?
 */
typedef enum {
  eConstraintPopup_String = 1,
  eConstraintPopup_Location,
  eConstraintPopup_Source,
  eConstraintPopup_CDSGeneProt,
  eConstraintPopup_CDSGeneProtPseudo,
  eConstraintPopup_Sequence,
  eConstraintPopup_Pub,
  eConstraintPopup_RNA,
  eConstraintPopup_FeatureField,
  eConstraintPopup_MiscField,
  eConstraintPopup_MolinfoField,
  eConstraintPopup_MissingField,
  eConstraintPopup_Translation,
  eConstraintPopup_Max
} EConstraintPopup;

typedef struct constrainttable {
  CharPtr title;
  Uint1   constraint_choice;
} ConstraintTableData, PNTR ConstraintTablePtr;

static ConstraintTableData ConstraintTable[] = {
  { "String", ConstraintChoice_string },
  {  "Location", ConstraintChoice_location },
  {  "Source", ConstraintChoice_source },
  {  "CDS-Gene-Prot Qualifier", ConstraintChoice_cdsgeneprot_qual },
  {  "CDS-Gene-Prot Pseudo Feature", ConstraintChoice_cdsgeneprot_pseudo },
  {  "Sequence and Feature", ConstraintChoice_sequence },
  {  "Publication", ConstraintChoice_pub },
  {  "RNA Field", ConstraintChoice_field },
  {  "Feature Field", ConstraintChoice_field },
  {  "Misc Field", ConstraintChoice_field },
  {  "Molinfo Field", ConstraintChoice_molinfo },
  {  "Missing Field", ConstraintChoice_field_missing },
  {  "Translation", ConstraintChoice_translation }
};


static Int2 PopupFromConstraintChoice (Uint1 constraint_choice)
{
  Int2 val;
  for (val = 0; val < eConstraintPopup_Max - 1; val++) {
    if (ConstraintTable[val].constraint_choice == constraint_choice) {
      return val + 1;
    }
  }
  return 0;
}


typedef struct editconstraint {
  PopuP constraint_type;

  DialoG dlgs[eConstraintPopup_Max];

  ButtoN accept_btn;
} EditConstraintData, PNTR EditConstraintPtr;


static void EnableEditConstraintAccept (Pointer data)
{
  EditConstraintPtr ecp;
  Int2 val;
  ValNodePtr err_list = NULL;
  Boolean ok_to_accept = TRUE;
  ValNodePtr tmp;
  StringConstraintPtr scp;

  ecp = (EditConstraintPtr) data;
  if (ecp != NULL && ecp->accept_btn != NULL) {
    val = GetValue (ecp->constraint_type);
    switch (val) {
      case eConstraintPopup_String:
        scp = DialogToPointer (ecp->dlgs[val - 1]);
        if (IsStringConstraintEmpty (scp)) {
          ValNodeAddPointer (&err_list, 0, "empty string constraint");
        }
        scp = StringConstraintFree (scp);
        break;
      case eConstraintPopup_Location:
      case eConstraintPopup_Source:
      case eConstraintPopup_CDSGeneProt:
      case eConstraintPopup_CDSGeneProtPseudo:
      case eConstraintPopup_Sequence:
      case eConstraintPopup_Pub:
      case eConstraintPopup_MolinfoField:
      case eConstraintPopup_MissingField:
      case eConstraintPopup_Translation:
        err_list = TestDialog (ecp->dlgs[val - 1]);
        break;      
      case eConstraintPopup_RNA:
      case eConstraintPopup_FeatureField:
      case eConstraintPopup_MiscField:
        tmp = DialogToPointer (ecp->dlgs[val - 1]);
        if (tmp == NULL || IsFieldConstraintEmpty (tmp->data.ptrvalue)) {
          ValNodeAddPointer (&err_list, 0, "No constraint");
        }
        tmp = ConstraintChoiceFree (tmp);
        break;
      default:
        ok_to_accept = FALSE;
    }
    if (err_list != NULL) {
      ok_to_accept = FALSE;
      err_list = ValNodeFree (err_list);
    }
    if (ok_to_accept) {
      Enable (ecp->accept_btn);
    } else {
      Disable (ecp->accept_btn);
    }
  }
}


static void ChangeEditConstraintType (PopuP p)
{
  EditConstraintPtr ecp;
  Int2 val;

  ecp = (EditConstraintPtr) GetObjectExtra (p);
  if (ecp == NULL) return;

  for (val = 0; val < eConstraintPopup_Max - 1; val++) {
    Hide (ecp->dlgs[val]);
  }

  val = GetValue (ecp->constraint_type);
  if (val > 0 && val < eConstraintPopup_Max) {
    Show (ecp->dlgs[val - 1]);
  }
  EnableEditConstraintAccept (ecp);
}


static void FreeConstraintData (ValNodePtr constraint)
{
  if (constraint == NULL) return;

  switch (constraint->choice) {
    case ConstraintChoice_string:
      constraint->data.ptrvalue = StringConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_location:
      constraint->data.ptrvalue = LocationConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_source:
      constraint->data.ptrvalue = SourceConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_cdsgeneprot_qual:
      constraint->data.ptrvalue = CDSGeneProtQualConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_cdsgeneprot_pseudo:
      constraint->data.ptrvalue = CDSGeneProtPseudoConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_sequence:
      constraint->data.ptrvalue = SequenceConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_pub:
      constraint->data.ptrvalue = PublicationConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_field:
      constraint->data.ptrvalue = FieldConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_molinfo:
      constraint->data.ptrvalue = MolinfoFieldConstraintFree (constraint->data.ptrvalue);
      break;
    case ConstraintChoice_translation:
      constraint->data.ptrvalue = TranslationConstraintFree (constraint->data.ptrvalue);
      break;
  }
}

static DialoG SingleFieldTypeDialog (GrouP h, Boolean text_only, Boolean for_remove, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

static Boolean EditConstraint (ValNodePtr constraint)
{
  ModalAcceptCancelData acd;
  EditConstraintData    ecd;
  Boolean               rval = FALSE;
  WindoW                w;
  ButtoN                b;
  GrouP                 h, g, c;
  Int2                  val;
  FieldConstraintPtr    fcp;
  ValNodePtr            tmp;
  
  if (constraint == NULL) return FALSE;

  w = MovableModalWindow(-20, -13, -10, -10, "Constraint", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ecd.constraint_type = PopupList (h, TRUE, ChangeEditConstraintType);
  SetObjectExtra (ecd.constraint_type, &ecd, NULL);
  for (val = 1; val < eConstraintPopup_Max; val++) {
    PopupItem (ecd.constraint_type, ConstraintTable[val - 1].title);
  }

  g = HiddenGroup (h, 0, 0, NULL);
  ecd.accept_btn = NULL;
  ecd.dlgs[eConstraintPopup_String - 1] = StringConstraintDialog (g, "", TRUE, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_Location - 1] = LocationConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_Source - 1] = SourceConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_CDSGeneProt - 1] = CDSGeneProtQualConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_CDSGeneProtPseudo - 1] = CDSGeneProtPseudoConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_Sequence - 1] = SequenceConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_Pub - 1] = PublicationConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_RNA - 1] = RNAFieldConstraintDialog (g, "RNA Type", EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_FeatureField - 1] = FeatureFieldConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_MiscField - 1] = MiscFieldConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_MolinfoField - 1] = MolinfoConstraintDialog (g, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_MissingField - 1] = SingleFieldTypeDialog(g, TRUE, FALSE, EnableEditConstraintAccept, &ecd);
  ecd.dlgs[eConstraintPopup_Translation - 1] = TranslationConstraintDialog (g, EnableEditConstraintAccept, &ecd);

  AlignObjects (ALIGN_CENTER, (HANDLE) ecd.dlgs[eConstraintPopup_String - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_Location - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_Source - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_CDSGeneProt - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_CDSGeneProtPseudo - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_Sequence - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_Pub - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_RNA - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_FeatureField - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_MiscField - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_MolinfoField - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_MissingField - 1],
                              (HANDLE) ecd.dlgs[eConstraintPopup_Translation - 1],
                              NULL);

  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  ecd.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (ecd.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ecd.constraint_type,
                              (HANDLE) g,
                              (HANDLE) c, 
                              NULL);

  if (constraint->choice == ConstraintChoice_field) {
    fcp = (FieldConstraintPtr) constraint->data.ptrvalue;
    if (fcp == NULL || fcp->field == NULL || fcp->field->choice == FieldType_rna_field) {
      SetValue (ecd.constraint_type, eConstraintPopup_RNA);
      PointerToDialog (ecd.dlgs[eConstraintPopup_RNA - 1], constraint->data.ptrvalue);
    } else if (fcp->field->choice == FieldType_feature_field) {
      SetValue (ecd.constraint_type, eConstraintPopup_FeatureField);
      PointerToDialog (ecd.dlgs[eConstraintPopup_FeatureField - 1], constraint->data.ptrvalue);
    } else if (fcp->field->choice == FieldType_misc) {
      SetValue (ecd.constraint_type, eConstraintPopup_MiscField);
      PointerToDialog (ecd.dlgs[eConstraintPopup_MiscField - 1], constraint->data.ptrvalue);
    }
  } else {
    val = PopupFromConstraintChoice (constraint->choice);
    if (val > 0 && val < eConstraintPopup_Max) {
      SetValue (ecd.constraint_type, val);
      PointerToDialog (ecd.dlgs[val - 1], constraint->data.ptrvalue);
    } else {
      SetValue (ecd.constraint_type, eConstraintPopup_String);
    }
  } 

  ChangeEditConstraintType (ecd.constraint_type);
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
  if (!acd.cancelled)
  {
    val = GetValue (ecd.constraint_type);
    switch (val) {
      case eConstraintPopup_String:
      case eConstraintPopup_Location:
      case eConstraintPopup_Source:
      case eConstraintPopup_CDSGeneProt:
      case eConstraintPopup_CDSGeneProtPseudo:
      case eConstraintPopup_Sequence:
      case eConstraintPopup_Pub:
      case eConstraintPopup_MolinfoField:
      case eConstraintPopup_MissingField:
      case eConstraintPopup_Translation:
        FreeConstraintData (constraint);
        constraint->choice = ConstraintTable[val - 1].constraint_choice;
        constraint->data.ptrvalue = DialogToPointer (ecd.dlgs[val - 1]);
        rval = TRUE;
        break;
      case eConstraintPopup_RNA:
      case eConstraintPopup_FeatureField:
      case eConstraintPopup_MiscField:
        FreeConstraintData (constraint);
        constraint->choice = ConstraintChoice_field;
        tmp = DialogToPointer (ecd.dlgs[val - 1]);
        if (tmp != NULL) {
          constraint->data.ptrvalue = tmp->data.ptrvalue;
          tmp = ValNodeFree (tmp);
        }
        rval = TRUE;
        break;
    }
  }
  Remove (w);
  return rval;
}


typedef struct constraintsetdlg {
  DIALOG_MESSAGE_BLOCK
  DoC                  constraint_doc;
  ValNodePtr           constraint_list;
  ValNodePtr           default_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} ConstraintSetDlgData, PNTR ConstraintSetDlgPtr;

static void PopulateConstraintDoc (DoC d, ValNodePtr constraint_list)
{
  ValNodePtr vnp;
  CharPtr    phrase, tmp;
  RecT       r;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };

  if (d == NULL) return;

  Reset (d);

  ObjectRect (d, &r);
  InsetRect (&r, 4, 4);
  
  ColFmt[1].pixWidth = r.right - r.left - 12;

  Reset (d);

  for (vnp = constraint_list; vnp != NULL; vnp = vnp->next) {
    phrase = SummarizeConstraint (vnp);
    if (phrase == NULL) {
      AppendText (d, "\tUnable to summarize constraint\n", &ParFmt, ColFmt, programFont);
    } else {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (phrase) + 3));
      sprintf (tmp, "\t%s\n", phrase);
      phrase = MemFree (phrase);
      AppendText (d, tmp, &ParFmt, ColFmt, programFont);
      tmp = MemFree (tmp);
    }
  }
  AppendText (d, "(Click here to add new constraint)", NULL, NULL, programFont);
  UpdateDocument (d, 0, 0);
}

static void DrawConstraintDocControls (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RecT                rct;
  Int4                width;
  PoinT               pt1, pt2;
  ConstraintSetDlgPtr dlg;

  dlg = (ConstraintSetDlgPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0 && item <= ValNodeLen (dlg->constraint_list)) {
    rct = *r;
  
    /* draw X for deletion */
    width = 10;
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1;
    pt2.x = pt1.x + width;
    pt2.y = pt1.y + width;
    DrawLine (pt1, pt2);
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1 + width;
    pt2.x = pt1.x + width;
    pt2.y = rct.top + 1;
    DrawLine (pt1, pt2);
  }
}


static ValNodePtr DefaultFeatureFieldConstraint (Int2 feat_type);

static ValNodePtr NewConstraintFromDefault (ValNodePtr default_constraint)
{
  ValNodePtr vnp = NULL;
  FieldConstraintPtr f;
  FeatureFieldPtr ff;
  ErrSev     oldErrSev;

  if (default_constraint == NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = ConstraintChoice_string;
  } else {
    if (default_constraint->choice == ConstraintChoice_field
        && (f = default_constraint->data.ptrvalue) != NULL
        && f->field != NULL
        && f->field->choice == FieldType_feature_field
        && (ff = f->field->data.ptrvalue) != NULL
        && ff->field == NULL) {
      vnp = DefaultFeatureFieldConstraint (ff->type);
    } else if (default_constraint->data.ptrvalue != NULL) {
      oldErrSev = ErrSetMessageLevel (SEV_FATAL);
      vnp = AsnIoMemCopy (default_constraint, (AsnReadFunc) ConstraintChoiceAsnRead, (AsnWriteFunc) ConstraintChoiceAsnWrite);
      ErrSetMessageLevel (oldErrSev);
    }
    if (vnp == NULL) {
      vnp = ValNodeNew (NULL);
      vnp->choice = default_constraint->choice;
    }
  }
  return vnp;
}


static void ClickConstraintDoc (DoC d, PoinT pt)
{
  Int2      item, row, col;
  RecT      rct;
  BaR       sb_vert;
  Int4      scroll_pos = 0, scroll_max;
  ConstraintSetDlgPtr f;
  ValNodePtr         vnp, vnp_prev = NULL;
  Boolean            changed = FALSE;
  
  f = (ConstraintSetDlgPtr) GetObjectExtra (d);
  if (f == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item == 0 && row == 0 && f->constraint_list == NULL) {
    /* create new constraint */
    vnp = NewConstraintFromDefault (f->default_constraint);
    if (EditConstraint (vnp)) {
      f->constraint_list = vnp;
      changed = TRUE;
    } else {
      vnp = ConstraintChoiceFree (vnp);
    }
  } else if (item > 0 && row > 0) {
    if (item == ValNodeLen (f->constraint_list) + 1) {
      /* create new constraint */
      vnp = NewConstraintFromDefault (f->default_constraint);
      if (EditConstraint (vnp)) {
        ValNodeLink (&(f->constraint_list), vnp);
        changed = TRUE;
      } else {
        vnp = ConstraintChoiceFree (vnp);
      }
    } else {
      for (vnp = f->constraint_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
        vnp_prev = vnp;
      }
      if (vnp != NULL) {      
        sb_vert = GetSlateVScrollBar ((SlatE) f->constraint_doc);
        scroll_pos = GetBarValue (sb_vert);
        switch (col) {
          case 1:
            /* delete this item */
            if (vnp_prev == NULL) {
              f->constraint_list = vnp->next;
            } else {
              vnp_prev->next = vnp->next;
            }
            vnp->next = NULL;
            vnp = ConstraintChoiceFree (vnp);
            changed = TRUE;
            break;
          case 2:
            /* edit */
            changed = EditConstraint (vnp);
            break;
        }
      }
    }
  }
  if (changed) {
    PopulateConstraintDoc (f->constraint_doc, f->constraint_list);
    if (scroll_pos > 0) {
      sb_vert = GetSlateVScrollBar ((SlatE) f->constraint_doc);
      scroll_max = GetBarMax (sb_vert);
      if (scroll_pos > scroll_max) {
        scroll_pos = scroll_max;
      }
      CorrectBarValue (sb_vert, scroll_pos);
    }
    if (f->change_notify != NULL) {
      (f->change_notify) (f->change_userdata);
    }
  }
}

static void ConstraintSetToDialog (DialoG d, Pointer data)
{
  ConstraintSetDlgPtr dlg;

  dlg = (ConstraintSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dlg->constraint_list = ConstraintChoiceSetFree (dlg->constraint_list);
  if (data != NULL) {
    dlg->constraint_list = AsnIoMemCopy ((ConstraintChoiceSetPtr) data,
                                         (AsnReadFunc) ConstraintChoiceSetAsnRead,
                                         (AsnWriteFunc) ConstraintChoiceSetAsnWrite);
 }
 PopulateConstraintDoc (dlg->constraint_doc, dlg->constraint_list);
 if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer ConstraintSetFromDialog (DialoG d)
{
  ConstraintSetDlgPtr dlg;
  ValNodePtr constraint_list = NULL;

  dlg = (ConstraintSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->constraint_list != NULL) {
    constraint_list = AsnIoMemCopy ((ConstraintChoiceSetPtr) dlg->constraint_list,
                                    (AsnReadFunc) ConstraintChoiceSetAsnRead,
                                    (AsnWriteFunc) ConstraintChoiceSetAsnWrite);
  }
  return (Pointer) constraint_list;
}

static void CleanupConstraintSetDialog (GraphiC g, VoidPtr data)

{
  ConstraintSetDlgPtr dlg;

  dlg = (ConstraintSetDlgPtr) data;
  if (dlg != NULL) {
    dlg->constraint_list = ConstraintChoiceSetFree (dlg->constraint_list);
    dlg->default_constraint = ConstraintChoiceFree (dlg->default_constraint);
  }
  StdCleanupExtraProc (g, data);
}


static void SetConstraintSetDefaultConstraintType (DialoG d, Uint1 constraint_type)
{
  ConstraintSetDlgPtr dlg;

  dlg = (ConstraintSetDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    dlg->default_constraint = ConstraintChoiceFree (dlg->default_constraint);
    dlg->default_constraint = ValNodeNew (NULL);
    dlg->default_constraint->choice = constraint_type;
  }
}


static void SetConstraintSetDefaultConstraintTypeEx (DialoG d, ValNodePtr constraint)
{
  ConstraintSetDlgPtr dlg;

  dlg = (ConstraintSetDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    dlg->default_constraint = ConstraintChoiceFree (dlg->default_constraint);
    if (constraint != NULL) {
      dlg->default_constraint = NewConstraintFromDefault (constraint);
    }
  }
}


static void AddConstraintBtn (ButtoN b)
{
  ConstraintSetDlgPtr dlg;
  ValNodePtr          vnp;
  BaR       sb_vert;
  Int4      scroll_pos = 0, scroll_max;

  dlg = (ConstraintSetDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  sb_vert = GetSlateVScrollBar ((SlatE) dlg->constraint_doc);
  scroll_pos = GetBarValue (sb_vert);

  /* create new constraint */
  vnp = NewConstraintFromDefault (dlg->default_constraint);
  if (EditConstraint (vnp)) {
    ValNodeLink (&(dlg->constraint_list), vnp);
    PopulateConstraintDoc (dlg->constraint_doc, dlg->constraint_list);
    if (scroll_pos > 0) {
      sb_vert = GetSlateVScrollBar ((SlatE) dlg->constraint_doc);
      scroll_max = GetBarMax (sb_vert);
      if (scroll_pos > scroll_max) {
        scroll_pos = scroll_max;
      }
      CorrectBarValue (sb_vert, scroll_pos);
    }
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  } else {
    vnp = ConstraintChoiceFree (vnp);
  }

}


NLM_EXTERN DialoG ConstraintSetDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ConstraintSetDlgPtr dlg;
  GrouP               p, g;
  PrompT              ppt;
  ButtoN              b;

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  dlg = (ConstraintSetDlgPtr) MemNew (sizeof (ConstraintSetDlgData));
  SetObjectExtra (p, dlg, CleanupConstraintSetDialog);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ConstraintSetToDialog;
  dlg->fromdialog = ConstraintSetFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  dlg->constraint_list = NULL;
  ppt = StaticPrompt (p, "Constraint List", 0, dialogTextHeight, programFont, 'c');
  dlg->constraint_doc = DocumentPanel (p, stdCharWidth * 30, stdLineHeight * 3);
  SetObjectExtra (dlg->constraint_doc, dlg, NULL);
  SetDocProcs (dlg->constraint_doc, ClickConstraintDoc, NULL, NULL, NULL);
  SetDocShade (dlg->constraint_doc, DrawConstraintDocControls, NULL, NULL, NULL);
  PopulateConstraintDoc (dlg->constraint_doc, dlg->constraint_list);

  g = HiddenGroup (p, 2, 0, NULL);
  b = PushButton (g, "Add Contraint", AddConstraintBtn);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (g, "Clear Constraints", ClearDialogBtn);
  SetObjectExtra (b, p, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->constraint_doc, (HANDLE) g, NULL);

  dlg->default_constraint = ValNodeNew (NULL);
  dlg->default_constraint->choice = ConstraintChoice_string;
  return (DialoG) p;
}


typedef struct texttransformtable {
  CharPtr title;
  Uint1   transform_choice;
} TextTransformTableData, PNTR TextTransformTablePtr;

static TextTransformTableData TextTransformTable[] = {
  {  "Edit", TextTransform_edit },
  {  "Capitalization", TextTransform_caps },
  {  "Remove Text Outside", TextTransform_remove }
};


typedef enum {
  eTextTransformPopup_Edit = 1,
  eTextTransformPopup_Caps,
  eTextTransformPopup_RemoveOutside,
  eTextTransformPopup_Max
} ETextTransformPopup;


static Int2 PopupFromTextTransformChoice (Uint1 transform_choice)
{
  Int2 val;
  for (val = 0; val < eTextTransformPopup_Max - 1; val++) {
    if (TextTransformTable[val].transform_choice == transform_choice) {
      return val + 1;
    }
  }
  return 0;
}


typedef struct edittexttransform {
  PopuP transform_type;

  DialoG dlgs[eTextTransformPopup_Max];

  ButtoN accept_btn;
} EditTextTransformData, PNTR EditTextTransformPtr;


static void EnableEditTextTransformAccept (Pointer data)
{
  EditTextTransformPtr ecp;
  Int2 val;
  Boolean ok_to_accept = TRUE;
  ValNodePtr vnp;

  ecp = (EditTextTransformPtr) data;
  if (ecp != NULL && ecp->accept_btn != NULL) {
    val = GetValue (ecp->transform_type);
    switch (val) {
      case eTextTransformPopup_Edit:
      case eTextTransformPopup_RemoveOutside:
        vnp = ValNodeNew (NULL);
        vnp->choice = TextTransformTable[val - 1].transform_choice;
        vnp->data.ptrvalue = DialogToPointer (ecp->dlgs[val - 1]);
        if (IsTextTransformEmpty(vnp)) {
          ok_to_accept = FALSE;
        }
        vnp = TextTransformFree (vnp);
        break;
      case eTextTransformPopup_Caps:
        if (GetCapChangeDialogValue(ecp->dlgs[val - 1]) == Cap_change_none) {
          ok_to_accept = FALSE;
        }
        break;
      default:
        ok_to_accept = FALSE;
        break;
    }
    if (ok_to_accept) {
      Enable (ecp->accept_btn);
    } else {
      Disable (ecp->accept_btn);
    }
  }
}


static void ChangeEditTextTransformType (PopuP p)
{
  EditTextTransformPtr ecp;
  Int2 val;

  ecp = (EditTextTransformPtr) GetObjectExtra (p);
  if (ecp == NULL) return;

  for (val = 0; val < eTextTransformPopup_Max - 1; val++) {
    Hide (ecp->dlgs[val]);
  }

  val = GetValue (ecp->transform_type);
  if (val > 0 && val < eTextTransformPopup_Max) {
    Show (ecp->dlgs[val - 1]);
  }
  EnableEditTextTransformAccept (ecp);
}


static void FreeTextTransformData (ValNodePtr transform)
{
  if (transform == NULL) return;

  switch (transform->choice) {
    case TextTransform_edit:
      transform->data.ptrvalue = FieldEditFree (transform->data.ptrvalue);
      break;
  }
}


static Boolean EditTextTransform (ValNodePtr transform)
{
  ModalAcceptCancelData acd;
  EditTextTransformData    ecd;
  Boolean               rval = FALSE;
  WindoW                w;
  ButtoN                b;
  GrouP                 h, g, c;
  Int2                  val;
  
  if (transform == NULL) return FALSE;

  w = MovableModalWindow(-20, -13, -10, -10, "Text Transform", NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ecd.transform_type = PopupList (h, TRUE, ChangeEditTextTransformType);
  SetObjectExtra (ecd.transform_type, &ecd, NULL);
  for (val = 1; val < eTextTransformPopup_Max; val++) {
    PopupItem (ecd.transform_type, TextTransformTable[val - 1].title);
  }

  g = HiddenGroup (h, 0, 0, NULL);
  ecd.accept_btn = NULL;
  ecd.dlgs[eTextTransformPopup_Edit - 1] = FieldEditDialog (g, EnableEditTextTransformAccept, &ecd);
  ecd.dlgs[eTextTransformPopup_Caps - 1] = CapChangeDialog (g, EnableEditTextTransformAccept, &ecd);
  ecd.dlgs[eTextTransformPopup_RemoveOutside - 1] = TextPortionDialog (g, FALSE, EnableEditTextTransformAccept, &ecd);

  AlignObjects (ALIGN_CENTER, (HANDLE) ecd.dlgs[eTextTransformPopup_Edit - 1],
                              (HANDLE) ecd.dlgs[eTextTransformPopup_Caps - 1],
                              (HANDLE) ecd.dlgs[eTextTransformPopup_RemoveOutside - 1],
                              NULL);

  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  ecd.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (ecd.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ecd.transform_type,
                              (HANDLE) g,
                              (HANDLE) c, 
                              NULL);

  val = PopupFromTextTransformChoice (transform->choice);
  if (val > 0 && val < eTextTransformPopup_Max) {
    SetValue (ecd.transform_type, val);
    if (val == eTextTransformPopup_Caps) {
      SetCapChangeDialogValue(ecd.dlgs[val - 1], transform->data.intvalue);
    } else {
      PointerToDialog (ecd.dlgs[val - 1], transform->data.ptrvalue);
    }
  } else {
    SetValue (ecd.transform_type, eTextTransformPopup_Edit);
  }

  ChangeEditTextTransformType (ecd.transform_type);
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
  if (!acd.cancelled)
  {
    val = GetValue (ecd.transform_type);
    switch (val) {
      case eTextTransformPopup_Edit:
      case eTextTransformPopup_RemoveOutside:
        FreeTextTransformData (transform);
        transform->choice = TextTransformTable[val - 1].transform_choice;
        transform->data.ptrvalue = DialogToPointer (ecd.dlgs[val - 1]);
        rval = TRUE;
        break;
      case eTextTransformPopup_Caps:
        FreeTextTransformData (transform);
        transform->choice = TextTransform_caps;
        transform->data.intvalue = GetCapChangeDialogValue(ecd.dlgs[val - 1]);
        rval = TRUE;
        break;
    }
  }
  Remove (w);
  return rval;
}


typedef struct texttransformsetdlg {
  DIALOG_MESSAGE_BLOCK
  DoC                  transform_doc;
  ValNodePtr           transform_list;
  ValNodePtr           default_transform;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} TextTransformSetDlgData, PNTR TextTransformSetDlgPtr;

static void PopulateTextTransformDoc (DoC d, ValNodePtr transform_list)
{
  ValNodePtr vnp;
  CharPtr    phrase, tmp;
  RecT       r;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };

  if (d == NULL) return;

  Reset (d);

  ObjectRect (d, &r);
  InsetRect (&r, 4, 4);
  
  ColFmt[1].pixWidth = r.right - r.left - 12;

  Reset (d);

  for (vnp = transform_list; vnp != NULL; vnp = vnp->next) {
    phrase = SummarizeTextTransform (vnp);
    if (phrase == NULL) {
      AppendText (d, "\tUnable to summarize text transform\n", &ParFmt, ColFmt, programFont);
    } else {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (phrase) + 3));
      sprintf (tmp, "\t%s\n", phrase);
      phrase = MemFree (phrase);
      AppendText (d, tmp, &ParFmt, ColFmt, programFont);
      tmp = MemFree (tmp);
    }
  }
  AppendText (d, "(Click here to add new text transform)", NULL, NULL, programFont);
  UpdateDocument (d, 0, 0);
}

static void DrawTextTransformDocControls (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  RecT                   rct;
  Int4                   width;
  PoinT                  pt1, pt2;
  TextTransformSetDlgPtr dlg;

  dlg = (TextTransformSetDlgPtr) GetObjectExtra (d);
  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0 && item <= ValNodeLen (dlg->transform_list)) {
    rct = *r;
  
    /* draw X for deletion */
    width = 10;
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1;
    pt2.x = pt1.x + width;
    pt2.y = pt1.y + width;
    DrawLine (pt1, pt2);
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1 + width;
    pt2.x = pt1.x + width;
    pt2.y = rct.top + 1;
    DrawLine (pt1, pt2);
  }
}

static void ClickTextTransformDoc (DoC d, PoinT pt)
{
  Int2      item, row, col;
  RecT      rct;
  BaR       sb_vert;
  Int4      scroll_pos = 0, scroll_max;
  TextTransformSetDlgPtr f;
  ValNodePtr         vnp, vnp_prev = NULL;
  Boolean            changed = FALSE;
  
  f = (TextTransformSetDlgPtr) GetObjectExtra (d);
  if (f == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item == 0 && row == 0 && f->transform_list == NULL) {
    /* create new constraint */
    vnp = ValNodeNew (NULL);
    vnp->choice = TextTransform_edit;
    if (EditTextTransform (vnp)) {
      f->transform_list = vnp;
      changed = TRUE;
    } else {
      vnp = TextTransformFree (vnp);
    }
  } else if (item > 0 && row > 0) {
    if (item == ValNodeLen (f->transform_list) + 1) {
      /* create new constraint */
      vnp = ValNodeNew (NULL);
      vnp->choice = TextTransform_edit;
      if (EditTextTransform (vnp)) {
        ValNodeLink (&(f->transform_list), vnp);
        changed = TRUE;
      } else {
        vnp = TextTransformFree (vnp);
      }
    } else {
      for (vnp = f->transform_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
        vnp_prev = vnp;
      }
      if (vnp != NULL) {      
        sb_vert = GetSlateVScrollBar ((SlatE) f->transform_doc);
        scroll_pos = GetBarValue (sb_vert);
        switch (col) {
          case 1:
            /* delete this item */
            if (vnp_prev == NULL) {
              f->transform_list = vnp->next;
            } else {
              vnp_prev->next = vnp->next;
            }
            vnp->next = NULL;
            vnp = TextTransformFree (vnp);
            changed = TRUE;
            break;
          case 2:
            /* edit */
            changed = EditTextTransform (vnp);
            break;
        }
      }
    }
  }
  if (changed) {
    PopulateTextTransformDoc (f->transform_doc, f->transform_list);
    if (scroll_pos > 0) {
      sb_vert = GetSlateVScrollBar ((SlatE) f->transform_doc);
      scroll_max = GetBarMax (sb_vert);
      if (scroll_pos > scroll_max) {
        scroll_pos = scroll_max;
      }
      CorrectBarValue (sb_vert, scroll_pos);
    }
    if (f->change_notify != NULL) {
      (f->change_notify) (f->change_userdata);
    }
  }
}


static void TextTransformSetToDialog (DialoG d, Pointer data)
{
  TextTransformSetDlgPtr dlg;

  dlg = (TextTransformSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dlg->transform_list = TextTransformSetFree (dlg->transform_list);
  if (data != NULL) {
    dlg->transform_list = AsnIoMemCopy ((TextTransformSetPtr) data,
                                         (AsnReadFunc) TextTransformSetAsnRead,
                                         (AsnWriteFunc) TextTransformSetAsnWrite);
 }
 PopulateTextTransformDoc (dlg->transform_doc, dlg->transform_list);
 if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer TextTransformSetFromDialog (DialoG d)
{
  TextTransformSetDlgPtr dlg;
  ValNodePtr transform_list = NULL;

  dlg = (TextTransformSetDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->transform_list != NULL) {
    transform_list = AsnIoMemCopy ((TextTransformSetPtr) dlg->transform_list,
                                    (AsnReadFunc) TextTransformSetAsnRead,
                                    (AsnWriteFunc) TextTransformSetAsnWrite);
  }
  return (Pointer) transform_list;
}

static void CleanupTextTransformSetDialog (GraphiC g, VoidPtr data)

{
  TextTransformSetDlgPtr dlg;

  dlg = (TextTransformSetDlgPtr) data;
  if (dlg != NULL) {
    dlg->transform_list = TextTransformSetFree (dlg->transform_list);
    dlg->default_transform = TextTransformFree (dlg->default_transform);
  }
  StdCleanupExtraProc (g, data);
}


static void AddTextTransformBtn (ButtoN b)
{
  TextTransformSetDlgPtr dlg;
  ValNodePtr          vnp;
  BaR       sb_vert;
  Int4      scroll_pos = 0, scroll_max;

  dlg = (TextTransformSetDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  sb_vert = GetSlateVScrollBar ((SlatE) dlg->transform_doc);
  scroll_pos = GetBarValue (sb_vert);

  /* create new constraint */
  vnp = ValNodeNew (NULL);
  vnp->choice = TextTransform_edit;
  if (EditTextTransform (vnp)) {
    ValNodeLink (&(dlg->transform_list), vnp);
    PopulateTextTransformDoc (dlg->transform_doc, dlg->transform_list);
    if (scroll_pos > 0) {
      sb_vert = GetSlateVScrollBar ((SlatE) dlg->transform_doc);
      scroll_max = GetBarMax (sb_vert);
      if (scroll_pos > scroll_max) {
        scroll_pos = scroll_max;
      }
      CorrectBarValue (sb_vert, scroll_pos);
    }
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  } else {
    vnp = TextTransformFree (vnp);
  }

}


NLM_EXTERN DialoG TextTransformSetDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TextTransformSetDlgPtr dlg;
  GrouP               p, g;
  PrompT              ppt;
  ButtoN              b;

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  dlg = (TextTransformSetDlgPtr) MemNew (sizeof (TextTransformSetDlgData));
  SetObjectExtra (p, dlg, CleanupTextTransformSetDialog);

  dlg->dialog = (DialoG) p;
  dlg->todialog = TextTransformSetToDialog;
  dlg->fromdialog = TextTransformSetFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  dlg->transform_list = NULL;
  ppt = StaticPrompt (p, "Text Transform List", 0, dialogTextHeight, programFont, 'c');
  dlg->transform_doc = DocumentPanel (p, stdCharWidth * 30, stdLineHeight * 2);
  SetObjectExtra (dlg->transform_doc, dlg, NULL);
  SetDocProcs (dlg->transform_doc, ClickTextTransformDoc, NULL, NULL, NULL);
  SetDocShade (dlg->transform_doc, DrawTextTransformDocControls, NULL, NULL, NULL);
  PopulateTextTransformDoc (dlg->transform_doc, dlg->transform_list);

  g = HiddenGroup (p, 2, 0, NULL);
  b = PushButton (g, "Add Text Transform", AddTextTransformBtn);
  SetObjectExtra (b, dlg, NULL);
  b = PushButton (g, "Clear Text Transforms", ClearDialogBtn);
  SetObjectExtra (b, p, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->transform_doc, (HANDLE) g, NULL);

  dlg->default_transform = ValNodeNew (NULL);
  dlg->default_transform->choice = TextTransform_edit;
  return (DialoG) p;
}



typedef struct oldcdsgeneprotconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG string_qual_choice;
  DialoG string_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} OldCDSGeneProtConstraintDlgData, PNTR OldCDSGeneProtConstraintDlgPtr;


static void ConstraintToOldCDSGeneProtConstraintDlg (DialoG d, Pointer data)
{
  OldCDSGeneProtConstraintDlgPtr dlg;
  ValNodePtr                     constraint;
  CDSGeneProtQualConstraintPtr   cgp_constraint;

  dlg = (OldCDSGeneProtConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  constraint = (ValNodePtr) data;
  if (constraint == NULL || constraint->choice != ConstraintChoice_cdsgeneprot_qual) {
    PointerToDialog (dlg->string_qual_choice, NULL);
    PointerToDialog (dlg->string_constraint, NULL);
  } else {
    cgp_constraint = (CDSGeneProtQualConstraintPtr) constraint->data.ptrvalue;    
    if (cgp_constraint == NULL) {
      PointerToDialog (dlg->string_qual_choice, NULL);
      PointerToDialog (dlg->string_constraint, NULL);
    } else {
      PointerToDialog (dlg->string_qual_choice, cgp_constraint->field1);
      PointerToDialog (dlg->string_constraint, cgp_constraint->constraint);
    }
  }
 
  if (dlg->change_notify) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer OldCDSGeneProtConstraintDlgToConstraint (DialoG d)
{
  OldCDSGeneProtConstraintDlgPtr dlg;
  ValNodePtr constraint = NULL;
  CDSGeneProtQualConstraintPtr   cgp_constraint;
  ValNodePtr                     vnp;
  StringConstraintPtr            string_constraint;

  dlg = (OldCDSGeneProtConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  /* string constraint */
  string_constraint = DialogToPointer (dlg->string_constraint);
  if (IsStringConstraintEmpty (string_constraint)) {
    string_constraint = StringConstraintFree (string_constraint);
  } else {
    vnp = DialogToPointer (dlg->string_qual_choice);
    if (vnp == NULL) {
      string_constraint = StringConstraintFree (string_constraint);
    } else {
      vnp->choice = CDSGeneProtConstraintField_field;
      cgp_constraint = CDSGeneProtQualConstraintNew ();
      cgp_constraint->field1 = vnp;
      cgp_constraint->constraint = string_constraint;
      ValNodeAddPointer (&constraint, ConstraintChoice_cdsgeneprot_qual, cgp_constraint);
    }
  }
  return (Pointer) constraint;
}


static void ClearOldCDSGeneProtConstraintDialogText (DialoG d)
{
  OldCDSGeneProtConstraintDlgPtr dlg;

  dlg = (OldCDSGeneProtConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  ClearStringConstraintDialogText (dlg->string_constraint);
}


static ValNodePtr TestOldCDSGeneProtConstraintDlg (DialoG d)
{
  OldCDSGeneProtConstraintDlgPtr dlg;
  ValNodePtr                     err_list = NULL;
  StringConstraintPtr            scp;

  dlg = (OldCDSGeneProtConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  
  scp = DialogToPointer (dlg->string_constraint);
  if (!IsStringConstraintEmpty (scp)) {
    ValNodeLink (&err_list, TestDialog (dlg->string_qual_choice));
  }
  scp = StringConstraintFree (scp);
  return err_list;
}


static DialoG CDSGeneProtFieldDialogEx (GrouP h, Boolean tall, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

static DialoG OldCDSGeneProtConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  OldCDSGeneProtConstraintDlgPtr dlg;
  GrouP p, q;
  ButtoN b;

  dlg = (OldCDSGeneProtConstraintDlgPtr) MemNew (sizeof (OldCDSGeneProtConstraintDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ConstraintToOldCDSGeneProtConstraintDlg;
  dlg->fromdialog = OldCDSGeneProtConstraintDlgToConstraint;
  dlg->testdialog = TestOldCDSGeneProtConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  q = HiddenGroup (p, 2, 0, NULL);
  dlg->string_qual_choice = CDSGeneProtFieldDialogEx (q, FALSE, change_notify, change_userdata);
  dlg->string_constraint = StringConstraintDialog (q, "", FALSE, change_notify, change_userdata);
  
  
  b = PushButton (p, "Clear Constraint", ClearDialogBtn);
  SetObjectExtra (b, p, NULL);    
  
  AlignObjects (ALIGN_CENTER, (HANDLE) q, (HANDLE) b, NULL);
   
  return (DialoG) p;
}


typedef struct oldsrcconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG string_qual_choice;
  DialoG string_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} OldSrcConstraintDlgData, PNTR OldSrcConstraintDlgPtr;


static void ConstraintToOldSrcConstraintDlg (DialoG d, Pointer data)
{
  OldSrcConstraintDlgPtr dlg;
  SourceConstraintPtr    src_constraint;
  ValNodePtr             constraint;

  dlg = (OldSrcConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  constraint = (ValNodePtr) data;
  if (constraint == NULL || constraint->choice != ConstraintChoice_source) {
    PointerToDialog (dlg->string_qual_choice, NULL);
    PointerToDialog (dlg->string_constraint, NULL);
  } else {
    src_constraint = (SourceConstraintPtr) constraint->data.ptrvalue;
    if (src_constraint == NULL) {
      PointerToDialog (dlg->string_qual_choice, NULL);
      PointerToDialog (dlg->string_constraint, NULL);
    } else {
      PointerToDialog (dlg->string_qual_choice, src_constraint->field1);
      PointerToDialog (dlg->string_constraint, src_constraint->constraint);
    }
  }
 
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ClearOldSrcConstraintDialogText (DialoG d)
{
  OldSrcConstraintDlgPtr dlg;
  dlg = (OldSrcConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  ClearStringConstraintDialogText (dlg->string_constraint);
}


static Pointer OldSrcConstraintDlgToConstraint (DialoG d)
{
  OldSrcConstraintDlgPtr dlg;
  ValNodePtr             constraint = NULL;
  SourceConstraintPtr    src_constraint;
  StringConstraintPtr    string_constraint;
  ValNodePtr             vnp;

  dlg = (OldSrcConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  /* string constraint */
  string_constraint = DialogToPointer (dlg->string_constraint);
  if (IsStringConstraintEmpty (string_constraint)) {
    string_constraint = StringConstraintFree (string_constraint);
  } else {
    vnp = DialogToPointer (dlg->string_qual_choice);
    if (vnp == NULL) {
      string_constraint = StringConstraintFree (string_constraint);
    } else {
      src_constraint = SourceConstraintNew ();
      src_constraint->field1 = vnp;
      src_constraint->constraint = string_constraint;
      ValNodeAddPointer (&constraint, ConstraintChoice_source, src_constraint);
    }
  }
  return (Pointer) constraint;
}


static ValNodePtr TestOldSrcConstraintDlg (DialoG d)
{
  OldSrcConstraintDlgPtr dlg;
  ValNodePtr             err_list = NULL;
  StringConstraintPtr    scp;

  dlg = (OldSrcConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  scp = DialogToPointer (dlg->string_constraint);
  if (!IsStringConstraintEmpty (scp)) {
      ValNodeLink (&err_list, TestDialog (dlg->string_qual_choice));
  }
  scp = StringConstraintFree (scp);
  return err_list;
}


static DialoG OldSrcConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  OldSrcConstraintDlgPtr dlg;
  GrouP p, q;
  ButtoN b;

  dlg = (OldSrcConstraintDlgPtr) MemNew (sizeof (OldSrcConstraintDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ConstraintToOldSrcConstraintDlg;
  dlg->fromdialog = OldSrcConstraintDlgToConstraint;
  dlg->testdialog = TestOldSrcConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  q = HiddenGroup (p, 2, 0, NULL);
  dlg->string_qual_choice = SourceQualChoiceDialog (q, TRUE, FALSE, TRUE, change_notify, change_userdata);
  dlg->string_constraint = StringConstraintDialog (q, "", FALSE, change_notify, change_userdata);
    
  b = PushButton (p, "Clear Constraint", ClearDialogBtn);
  SetObjectExtra (b, p, NULL);    
  
  AlignObjects (ALIGN_CENTER, (HANDLE) q, (HANDLE) b, NULL);
   
  return (DialoG) p;
}


typedef struct simplesequenceconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  GrouP seqtype;
  GrouP rna_subtype_grp;
  DialoG rna_subtype;
  DialoG src_loc;
  DialoG id;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
  
} SimpleSequenceConstraintDlgData, PNTR SimpleSequenceConstraintDlgPtr;


static void ClearSimpleSequenceConstraintDialogText (DialoG d)
{
  SimpleSequenceConstraintDlgPtr dlg;

  dlg = (SimpleSequenceConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  ClearStringConstraintDialogText (dlg->id);
}


static void ChangeSimpleSequenceConstraintGroup (GrouP g)
{
  SimpleSequenceConstraintDlgPtr dlg;

  dlg = (SimpleSequenceConstraintDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;

  if (GetValue (dlg->seqtype) == 4) {
    Show (dlg->rna_subtype_grp);
  } else {
    Hide (dlg->rna_subtype_grp);
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SimpleSequenceConstraintToDialog (DialoG d, Pointer data)
{
  SimpleSequenceConstraintDlgPtr dlg;
  ConstraintChoiceSetPtr constraint;
  SequenceConstraintPtr  seq;
  SourceConstraintPtr    src;
  ValNode vn;
  Int4    loc_val;

  dlg = (SimpleSequenceConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  PointerToDialog (dlg->src_loc, NULL);
  SetValue (dlg->seqtype, 2);
  PointerToDialog (dlg->rna_subtype, NULL);
  Hide (dlg->rna_subtype_grp);
  constraint = (ConstraintChoiceSetPtr) data;
  PointerToDialog (dlg->id, NULL);
  
  while (constraint != NULL) {
    if (constraint->choice == ConstraintChoice_sequence
      && (seq = (SequenceConstraintPtr) constraint->data.ptrvalue) != NULL) {
      if (seq->seqtype != NULL) {
        switch (seq->seqtype->choice) {
          case SequenceConstraintMolTypeConstraint_any :
            SetValue (dlg->seqtype, 1);
            PointerToDialog (dlg->rna_subtype, NULL);
            Hide (dlg->rna_subtype_grp);
            break;
          case SequenceConstraintMolTypeConstraint_nucleotide :
            SetValue (dlg->seqtype, 2);
            PointerToDialog (dlg->rna_subtype, NULL);
            Hide (dlg->rna_subtype_grp);
            break;
          case SequenceConstraintMolTypeConstraint_dna :
            SetValue (dlg->seqtype, 3);
            PointerToDialog (dlg->rna_subtype, NULL);
            Hide (dlg->rna_subtype_grp);
            break;
          case SequenceConstraintMolTypeConstraint_rna :
            SetValue (dlg->seqtype, 5);
            vn.choice = seq->seqtype->data.intvalue;
            vn.data.ptrvalue = NULL;
            vn.next = NULL;
            PointerToDialog (dlg->rna_subtype, &vn);
            Show (dlg->rna_subtype_grp);
            break;
          case SequenceConstraintMolTypeConstraint_protein :
            SetValue (dlg->seqtype, 4);
            PointerToDialog (dlg->rna_subtype, NULL);
            Hide (dlg->rna_subtype_grp);
            break;
        }
      }
      PointerToDialog (dlg->id, seq->id);
    } else if (constraint->choice == ConstraintChoice_source
      && (src = (SourceConstraintPtr) constraint->data.ptrvalue) != NULL) {
      if (src->field1 != NULL && src->field1->choice == SourceQualChoice_location) {
        vn.choice = src->field1->data.intvalue;
        vn.data.ptrvalue = NULL;
        vn.next = NULL;
        PointerToDialog (dlg->src_loc, &vn);
      } else if (src->field2 != NULL && src->field2->choice == SourceQualChoice_location) {
        vn.choice = src->field2->data.intvalue;
        vn.data.ptrvalue = NULL;
        vn.next = NULL;
        PointerToDialog (dlg->src_loc, &vn);
      }
      if (src->constraint != NULL 
          && src->constraint->match_location == String_location_equals
          && (loc_val = GenomeFromLocName (src->constraint->match_text)) > 0) {
        vn.choice = SrcLocFromGenome (loc_val);
        vn.data.ptrvalue = NULL;
        vn.next = NULL;
        PointerToDialog (dlg->src_loc, &vn);
      }
    }
    constraint = constraint->next;
  }
}


static Pointer SimpleSequenceConstraintFromDialog (DialoG d)
{
  SimpleSequenceConstraintDlgPtr dlg;
  SequenceConstraintPtr seq = NULL;
  SourceConstraintPtr   src = NULL;
  StringConstraintPtr   id = NULL;
  ValNodePtr            constraint = NULL, vnp;
  Int4 val;

  dlg = (SimpleSequenceConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->seqtype);
  switch (val) {
    case 1:
      /* don't bother to fill in, it's optional */
      break;
    case 2:
      seq = SequenceConstraintNew ();
      seq->seqtype = ValNodeNew (NULL);
      seq->seqtype->choice = SequenceConstraintMolTypeConstraint_nucleotide;
      break;
    case 3:
      seq = SequenceConstraintNew ();
      seq->seqtype = ValNodeNew (NULL);
      seq->seqtype->choice = SequenceConstraintMolTypeConstraint_dna;
      break;
    case 4:
      seq = SequenceConstraintNew ();
      seq->seqtype = ValNodeNew (NULL);
      seq->seqtype->choice = SequenceConstraintMolTypeConstraint_protein;
      break;
    case 5:
      seq = SequenceConstraintNew ();
      seq->seqtype = ValNodeNew (NULL);
      seq->seqtype->choice = SequenceConstraintMolTypeConstraint_rna;
      vnp = DialogToPointer (dlg->rna_subtype);
      if (vnp == NULL) {
        seq->seqtype->data.intvalue = Sequence_constraint_rnamol_any;
      } else {
        seq->seqtype->data.intvalue = vnp->choice;
      }
      break;
  }

  id = DialogToPointer (dlg->id);
  if (IsStringConstraintEmpty (id)) {
    id = StringConstraintFree (id);
  } else {
    if (seq == NULL) {
      seq = SequenceConstraintNew ();
    }
    seq->id = id;
  }

  if (seq != NULL) {
    ValNodeAddPointer (&constraint, ConstraintChoice_sequence, seq);
  }

  vnp = DialogToPointer (dlg->src_loc);
  if (vnp != NULL && vnp->choice != Source_location_unknown) {
    src = SourceConstraintNew ();
    src->field1 = ValNodeNew (NULL);
    src->field1->choice = SourceQualChoice_location;
    src->field1->data.intvalue = vnp->choice;
    src->constraint = StringConstraintNew ();
    src->constraint->match_location = String_location_equals;
    src->constraint->match_text = StringSave (LocNameFromGenome (GenomeFromSrcLoc(vnp->choice)));
    ValNodeAddPointer (&constraint, ConstraintChoice_source, src);
  }
  vnp = ValNodeFree (vnp);

  return constraint;
}


static DialoG SimpleSequenceConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  SimpleSequenceConstraintDlgPtr dlg;
  GrouP p, g, g2;
  ValNodePtr rna_subtypes = NULL, loc_list;

  dlg = (SimpleSequenceConstraintDlgPtr) MemNew (sizeof (SimpleSequenceConstraintDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = SimpleSequenceConstraintToDialog;
  dlg->fromdialog = SimpleSequenceConstraintFromDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g2 = HiddenGroup (p, 2, 0, NULL);
  dlg->seqtype = HiddenGroup (g2, 5, 0, ChangeSimpleSequenceConstraintGroup);
  SetObjectExtra (dlg->seqtype, dlg, NULL);
  RadioButton (dlg->seqtype, "Any sequence");
  RadioButton (dlg->seqtype, "Nucleotides");
  RadioButton (dlg->seqtype, "DNA");
  RadioButton (dlg->seqtype, "Proteins");
  RadioButton (dlg->seqtype, "RNA");
  SetValue (dlg->seqtype, 2);
  dlg->rna_subtype_grp = HiddenGroup (g2, 2, 0, NULL);
  StaticPrompt (dlg->rna_subtype_grp, "RNA Type", 0, dialogTextHeight, programFont, 'r');
  AddAllRNASubtypesToChoiceList (&rna_subtypes);
  dlg->rna_subtype = ValNodeSelectionDialog (dlg->rna_subtype_grp, rna_subtypes, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "rna subtype",
                                change_notify, change_userdata, FALSE);

  Hide (dlg->rna_subtype_grp);

  g = HiddenGroup (p, 2, 0, NULL);
  loc_list = GetLocationList (TRUE);
  StaticPrompt (g, "Where source location is", 0, dialogTextHeight, programFont, 'r');
  dlg->src_loc = ValNodeSelectionDialog (g,loc_list, 2, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "src location",
                                change_notify, change_userdata, FALSE);

  dlg->id = StringConstraintDialog (p, "Where sequence ID", FALSE, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) g2, (HANDLE) g, (HANDLE) dlg->id, NULL);
  return (DialoG) p;
}


typedef struct locorstringconstraintdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG tbs;
  DialoG pages[2];

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
  Int4 current_page;
} LocOrStringConstraintDlgData, PNTR LocOrStringConstraintDlgPtr;


static void ChangeLocOrStringConstraintPage (VoidPtr data, Int2 newval, Int2 oldval)
{
  LocOrStringConstraintDlgPtr dlg;

  dlg = (LocOrStringConstraintDlgPtr) data;
  if (dlg == NULL) return;
  if (newval == 0) {
    Show (dlg->pages[0]);
    Hide (dlg->pages[1]);
  } else {
    Show (dlg->pages[1]);
    Hide (dlg->pages[0]);
  }
  dlg->current_page = newval;
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void DataToLocOrStringConstraintDlg (DialoG d, Pointer data)
{
  LocOrStringConstraintDlgPtr dlg;
  ValNodePtr                  constraint;

  dlg = (LocOrStringConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  constraint = (ValNodePtr) data;
  if (constraint == NULL) {
    PointerToDialog (dlg->pages[0], NULL);
    PointerToDialog (dlg->pages[1], NULL);
    SetValue (dlg->tbs, 0);
    dlg->current_page = 0;
  } else if (constraint->choice == ConstraintChoice_string) {
    PointerToDialog (dlg->pages[0], constraint->data.ptrvalue);
    SetValue (dlg->tbs, 0);
    dlg->current_page = 0;
  } else if (constraint->choice == ConstraintChoice_location) {
    PointerToDialog (dlg->pages[1], constraint->data.ptrvalue);
    SetValue (dlg->tbs, 1);
    dlg->current_page = 1;
  }
  ChangeLocOrStringConstraintPage (dlg, dlg->current_page, 0);
}


static Pointer DataFromLocOrStringConstraintDialog (DialoG d)
{
  LocOrStringConstraintDlgPtr dlg;
  ValNodePtr                  constraint = NULL;
  StringConstraintPtr         string_constraint;
  LocationConstraintPtr       location_constraint;

  dlg = (LocOrStringConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->current_page == 0) {
    string_constraint = DialogToPointer (dlg->pages[0]);
    if (IsStringConstraintEmpty (string_constraint)) {
      string_constraint = StringConstraintFree (string_constraint);
    } else {
      ValNodeAddPointer (&constraint, ConstraintChoice_string, string_constraint);
    }
  } else {
    location_constraint = DialogToPointer (dlg->pages[1]);
    if (location_constraint != NULL) {
      ValNodeAddPointer (&constraint, ConstraintChoice_location, location_constraint);
    }
  }
  return constraint;
}


static ValNodePtr TestLocOrStringConstraintDlg (DialoG d)
{
  LocOrStringConstraintDlgPtr dlg;
  ValNodePtr                  constraint = NULL;
  StringConstraintPtr         string_constraint;
  LocationConstraintPtr       location_constraint;
  ValNodePtr err_list = NULL;

  dlg = (LocOrStringConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->current_page == 0) {
    string_constraint = DialogToPointer (dlg->pages[0]);
    if (!IsStringConstraintEmpty (string_constraint)) {
      ValNodeLink (&err_list, TestDialog (dlg->pages[0]));
    }
    string_constraint = StringConstraintFree (string_constraint);
  } else {
    location_constraint = DialogToPointer (dlg->pages[1]);
    if (location_constraint != NULL) {
      ValNodeLink (&err_list, TestDialog (dlg->pages[1]));
    }
    location_constraint = LocationConstraintFree (location_constraint);
  }
  return constraint;

}


static DialoG LocOrStringConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  GrouP         p, k;
  LocOrStringConstraintDlgPtr dlg;
  Int4          num_pages = 0;
  CharPtr       filterTabs[5];
  Int4          i;
    
  dlg = (LocOrStringConstraintDlgPtr) MemNew (sizeof (LocOrStringConstraintDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = NormalGroup (h, -1, 0, NULL, programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToLocOrStringConstraintDlg;
  dlg->fromdialog = DataFromLocOrStringConstraintDialog;
  dlg->testdialog = TestLocOrStringConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  dlg->current_page = 0;

  filterTabs [num_pages++] = "String Constraint";
  filterTabs [num_pages++] = "Location Constraint";
  filterTabs [num_pages] = NULL;
  
  dlg->tbs = CreateFolderTabs (p, filterTabs, 0,
                                0, 0, PROGRAM_FOLDER_TAB,
                                ChangeLocOrStringConstraintPage, (Pointer) dlg);
  k = HiddenGroup (p, 0, 0, NULL);
  num_pages = 0;
  dlg->pages [num_pages++] = StringConstraintDialog (k, NULL, TRUE, dlg->change_notify, dlg->change_userdata);
  dlg->pages [num_pages++] = LocationConstraintDialog (k, dlg->change_notify, dlg->change_userdata);
  for (i = 1; i < num_pages; i++)
  {
    Hide (dlg->pages [i]);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->pages [0],
                              (HANDLE) dlg->pages [1],
                              NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->tbs, (HANDLE) k, NULL);

  return (DialoG) p;
}


typedef struct rnafieldconstraintdlg {
  DIALOG_MESSAGE_BLOCK

  DialoG RNA_field;
  DialoG string_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} RNAFieldConstraintDlgData, PNTR RNAFieldConstraintDlgPtr;


static void ClearRnaFieldConstraintDialogText (DialoG d)
{
  RNAFieldConstraintDlgPtr  dlg;

  dlg = (RNAFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  ClearStringConstraintDialogText (dlg->string_constraint);
}


static Pointer RNAFieldConstraintDlgToData (DialoG d)
{
  RNAFieldConstraintDlgPtr  dlg;
  ValNodePtr                constraint = NULL;
  FieldConstraintPtr        f = NULL;
  RnaQualPtr                q;

  dlg = (RNAFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  q = DialogToPointer (dlg->RNA_field);
  if (q != NULL) {
    f = FieldConstraintNew ();
    f->field = ValNodeNew (NULL);
    f->field->choice = FieldType_rna_field;
    f->field->data.ptrvalue = q;
    f->string_constraint = DialogToPointer (dlg->string_constraint);
    ValNodeAddPointer (&constraint, ConstraintChoice_field, f);
  }
  return (Pointer) constraint;
}


static void DataToRNAFieldConstraintDlg (DialoG d, Pointer data)
{
  RNAFieldConstraintDlgPtr dlg;
  FieldConstraintPtr       f;
  RnaQualPtr               rq;

  dlg = (RNAFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  f = (FieldConstraintPtr) data;
  if (f == NULL) {
    PointerToDialog (dlg->RNA_field, NULL);
    PointerToDialog (dlg->string_constraint, NULL);
  } else {
    PointerToDialog (dlg->string_constraint, f->string_constraint);
    if (f->field == NULL) {
      PointerToDialog (dlg->RNA_field, NULL);
    } else if (f->field->choice == FieldType_rna_field) {
      PointerToDialog (dlg->RNA_field, f->field->data.ptrvalue);
    } else if (f->field->choice == FieldType_feature_field) {
      rq = RnaQualFromFeatureField (f->field->data.ptrvalue);
      PointerToDialog (dlg->RNA_field, rq);
      rq = RnaQualFree (rq);
    }
  }
}


static ValNodePtr TestRNAFieldConstraintDlg (DialoG d)
{
  RNAFieldConstraintDlgPtr dlg;
  ValNodePtr               err_list = NULL;

  dlg = (RNAFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  return err_list;
}

static DialoG RnaQualDialog (GrouP h, CharPtr type_label, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

static DialoG RNAFieldConstraintDialog (GrouP h, CharPtr type_label, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  GrouP p;
  RNAFieldConstraintDlgPtr dlg;

  dlg = (RNAFieldConstraintDlgPtr) MemNew (sizeof (RNAFieldConstraintDlgData));

  if (type_label == NULL) {
    p = NormalGroup (h, 2, 0, "RNA Field Constraint", programFont, NULL);
  } else {
    p = NormalGroup (h, -1, 0, "RNA Field Constraint", programFont, NULL);
  }
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToRNAFieldConstraintDlg;
  dlg->fromdialog = RNAFieldConstraintDlgToData;
  dlg->testdialog = TestRNAFieldConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->RNA_field = RnaQualDialog (p, type_label, dlg->change_notify, dlg->change_userdata);
  dlg->string_constraint = StringConstraintDialog (p, NULL, TRUE, dlg->change_notify, dlg->change_userdata);
  if (type_label != NULL) {
    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->RNA_field, (HANDLE) dlg->string_constraint, NULL);
  }

  return (DialoG) p;
}


typedef struct featurefieldconstraintdlg {
  DIALOG_MESSAGE_BLOCK

  DialoG feature_type;
  DialoG feature_field;
  DialoG string_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} FeatureFieldConstraintDlgData, PNTR FeatureFieldConstraintDlgPtr;


static void ClearFeatureFieldConstraintDialogText (DialoG d)
{
  FeatureFieldConstraintDlgPtr  dlg;

  dlg = (FeatureFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  ClearStringConstraintDialogText (dlg->string_constraint);
}


static Pointer FeatureFieldConstraintDlgToData (DialoG d)
{
  FeatureFieldConstraintDlgPtr  dlg;
  ValNodePtr                constraint = NULL, vnp, field = NULL;
  FieldConstraintPtr        f = NULL;
  FeatureFieldPtr           ff;

  dlg = (FeatureFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  vnp = DialogToPointer (dlg->feature_type);
  if (vnp != NULL) {
    ff = FeatureFieldNew();
    ff->type = vnp->choice;
    ff->field = DialogToPointer (dlg->feature_field);
    field = ValNodeNew (NULL);
    field->choice = FieldType_feature_field;
    field->data.ptrvalue = ff;
    vnp = ValNodeFreeData (vnp);
  }

  if (field != NULL) {
    f = FieldConstraintNew ();
    f->field = field;
    f->string_constraint = DialogToPointer (dlg->string_constraint);
    ValNodeAddPointer (&constraint, ConstraintChoice_field, f);
  }
  return (Pointer) constraint;
}


static void DataToFeatureFieldConstraintDlg (DialoG d, Pointer data)
{
  FeatureFieldConstraintDlgPtr dlg;
  FieldConstraintPtr           f;
  FeatureFieldPtr              ff;
  ValNode                      vn;

  dlg = (FeatureFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  f = (FieldConstraintPtr) data;
  if (f == NULL) {
    PointerToDialog (dlg->feature_type, NULL);
    PointerToDialog (dlg->feature_field, NULL);
    PointerToDialog (dlg->string_constraint, NULL);
  } else {
    PointerToDialog (dlg->string_constraint, f->string_constraint);

    if (f->field == NULL || f->field->choice != FieldType_feature_field) {
      PointerToDialog (dlg->feature_type, NULL);
      PointerToDialog (dlg->feature_field, NULL);
    } else {
      ff = (FeatureFieldPtr) f->field->data.ptrvalue;
      if (ff == NULL) {
        PointerToDialog (dlg->feature_type, NULL);
        PointerToDialog (dlg->feature_field, NULL);
      } else {
        vn.choice = (Uint1)ff->type;
        vn.data.ptrvalue = NULL;
        vn.next = NULL;
        PointerToDialog (dlg->feature_type, &vn);
        PointerToDialog (dlg->feature_field, ff->field);
      }
    }
  }
}


static ValNodePtr TestFeatureFieldConstraintDlg (DialoG d)
{
  FeatureFieldConstraintDlgPtr dlg;
  ValNodePtr               err_list = NULL;
  StringConstraintPtr      string_constraint;

  dlg = (FeatureFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  string_constraint = DialogToPointer (dlg->string_constraint);
  if (!IsStringConstraintEmpty (string_constraint)) {
    err_list = TestDialog (dlg->feature_type);
    ValNodeLink (&err_list, TestDialog (dlg->feature_field));
  }
  string_constraint = StringConstraintFree (string_constraint);
  return err_list;
}


static DialoG FeatQualChoiceDialog (GrouP h, Boolean allow_illegal, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

static void ClearFeatureFieldConstraint (ButtoN b) 
{
  FeatureFieldConstraintDlgPtr dlg;

  dlg = (FeatureFieldConstraintDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  PointerToDialog (dlg->feature_type, NULL);
  PointerToDialog (dlg->feature_field, NULL);
  PointerToDialog (dlg->string_constraint, NULL);
}


static DialoG FeatureFieldConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  GrouP p, g;
  FeatureFieldConstraintDlgPtr dlg;
  ButtoN b;

  dlg = (FeatureFieldConstraintDlgPtr) MemNew (sizeof (FeatureFieldConstraintDlgData));

  p = NormalGroup (h, 3, 0, "Feature Field Constraint", programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToFeatureFieldConstraintDlg;
  dlg->fromdialog = FeatureFieldConstraintDlgToData;
  dlg->testdialog = TestFeatureFieldConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->feature_type = FeatureTypeDialog (p, change_notify, change_userdata);
  dlg->feature_field = FeatQualChoiceDialog (p, FALSE, change_notify, change_userdata);

  g = HiddenGroup (p, -1, 0, NULL);
  dlg->string_constraint = StringConstraintDialog (g, NULL, FALSE, dlg->change_notify, dlg->change_userdata);
  b = PushButton (g, "Clear Constraint", ClearFeatureFieldConstraint);
  SetObjectExtra (b, dlg, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->string_constraint, (HANDLE) b, NULL);

  return (DialoG) p;
}

typedef struct miscfieldconstraintdlg {
  DIALOG_MESSAGE_BLOCK

  DialoG field;
  DialoG string_constraint;

  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} MiscFieldConstraintDlgData, PNTR MiscFieldConstraintDlgPtr;


static Pointer MiscFieldConstraintDlgToData (DialoG d)
{
  MiscFieldConstraintDlgPtr  dlg;
  ValNodePtr                constraint = NULL, field = NULL;
  FieldConstraintPtr        f = NULL;

  dlg = (MiscFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  field = DialogToPointer (dlg->field);
  
  if (field != NULL) {
    f = FieldConstraintNew ();
    f->field = field;
    f->string_constraint = DialogToPointer (dlg->string_constraint);
    ValNodeAddPointer (&constraint, ConstraintChoice_field, f);
  }
  return (Pointer) constraint;
}


static void DataToMiscFieldConstraintDlg (DialoG d, Pointer data)
{
  MiscFieldConstraintDlgPtr dlg;
  FieldConstraintPtr        f;

  dlg = (MiscFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  f = (FieldConstraintPtr) data;
  if (f == NULL) {
    PointerToDialog (dlg->field, NULL);
    PointerToDialog (dlg->string_constraint, NULL);
  } else {
    PointerToDialog (dlg->string_constraint, f->string_constraint);

    if (f->field == NULL || f->field->choice != FieldType_misc) {
      PointerToDialog (dlg->field, NULL);
    } else {
      PointerToDialog (dlg->field, f->field);
    }
  }
}


static ValNodePtr TestMiscFieldConstraintDlg (DialoG d)
{
  MiscFieldConstraintDlgPtr dlg;
  ValNodePtr                err_list = NULL;

  dlg = (MiscFieldConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  err_list = TestDialog (dlg->field);
  return err_list;
}


static void ClearMiscFieldConstraint (ButtoN b) 
{
  MiscFieldConstraintDlgPtr dlg;

  dlg = (MiscFieldConstraintDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  PointerToDialog (dlg->field, NULL);
  PointerToDialog (dlg->string_constraint, NULL);
}


static DialoG MiscFieldDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);


static DialoG MiscFieldConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  GrouP p;
  MiscFieldConstraintDlgPtr dlg;
  ButtoN b;

  dlg = (MiscFieldConstraintDlgPtr) MemNew (sizeof (MiscFieldConstraintDlgData));

  p = NormalGroup (h, -1, 0, "Misc Field Constraint", programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToMiscFieldConstraintDlg;
  dlg->fromdialog = MiscFieldConstraintDlgToData;
  dlg->testdialog = TestMiscFieldConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->field = MiscFieldDialog (p, change_notify, change_userdata);
  dlg->string_constraint = StringConstraintDialog (p, NULL, FALSE, dlg->change_notify, dlg->change_userdata);
  b = PushButton (p, "Clear Constraint", ClearMiscFieldConstraint);
  SetObjectExtra (b, dlg, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->field, (HANDLE) dlg->string_constraint, (HANDLE) b, NULL);

  return (DialoG) p;
}


static EComplexConstraintType ComplexConstraintTypeFromFeatureFieldType (Uint2 qualtype)
{
  EComplexConstraintType rval = eComplexConstraintType_string;

  switch (qualtype) {
    case FieldType_source_qual:
      rval = eComplexConstraintType_source;
      break;
    case FieldType_cds_gene_prot:
      rval = eComplexConstraintType_cdsgeneprot;
      break;
    case FieldType_pub:
      rval = eComplexConstraintType_pub;
      break;
    case FieldType_feature_field:
      rval = eComplexConstraintType_feature_field;
      break;
    case FieldType_rna_field:
      rval = eComplexConstraintType_rna_field;
      break;
    case FieldType_molinfo_field:
      rval = eComplexConstraintType_molinfo_field;
      break;
    case FieldType_struc_comment_field:
      rval = eComplexConstraintType_seqid;
      break;
  }
  return rval;
}


static Uint2 FeatureFieldTypeFromComplexConstraintType (EComplexConstraintType constraint_type)
{
  Uint2 rval = 0;

  switch (constraint_type) {
    case eComplexConstraintType_source:
      rval = FieldType_source_qual;
      break;
    case eComplexConstraintType_cdsgeneprot:
      rval = FieldType_cds_gene_prot;
      break;
    case eComplexConstraintType_pub:
      rval = FieldType_pub;
      break;
    case eComplexConstraintType_feature_field:
      rval = FieldType_feature_field;
      break;
    case eComplexConstraintType_rna_field:
      rval = FieldType_rna_field;
      break;
    case eComplexConstraintType_molinfo_field:
      rval = FieldType_molinfo_field;
      break;
  }
  return rval;
}


typedef struct complexconstraintdlg {
  DIALOG_MESSAGE_BLOCK

  GrouP  constraint_choice;
  GrouP  old_constraint_grp;
  DialoG cgp_constraint;
  DialoG src_constraint;
  DialoG rna_constraint;
  DialoG pub_constraint;
  DialoG feat_field_constraint;
  DialoG seq_constraint;
  DialoG string_constraint;
  DialoG constraint_set_dlg;

  ValNodePtr rna_type;
  Int2       feat_type;

  EComplexConstraintType constraint_type;
  Nlm_ChangeNotifyProc change_notify;
  Pointer              change_userdata;
} ComplexConstraintDlgData, PNTR ComplexConstraintDlgPtr;


static void ShowSimpleConstraint (ComplexConstraintDlgPtr dlg)
{
  if (dlg == NULL) {
    return;
  }

  Hide (dlg->src_constraint);
  Hide (dlg->cgp_constraint);
  Hide (dlg->pub_constraint);
  Hide (dlg->feat_field_constraint);
  Hide (dlg->rna_constraint);
  Hide (dlg->seq_constraint);
  Hide (dlg->string_constraint);

  switch (dlg->constraint_type) {
    case eComplexConstraintType_source:
      Show (dlg->src_constraint);
      break;
    case eComplexConstraintType_cdsgeneprot:
      Show (dlg->cgp_constraint);
      break;
    case eComplexConstraintType_pub:
      Show (dlg->pub_constraint);
      break;
    case eComplexConstraintType_feature_field:
      Show (dlg->feat_field_constraint);
      break;
    case eComplexConstraintType_rna_field:
      Show (dlg->rna_constraint);
      break;
    case eComplexConstraintType_molinfo_field:
      Show (dlg->seq_constraint);
      break;
    case eComplexConstraintType_seqid:
      Show (dlg->seq_constraint);
      break;
    case eComplexConstraintType_string:
      Show (dlg->string_constraint);
      break;
    default:
      break;
  }
}


static ValNodePtr DefaultRNAFieldConstraint (ValNodePtr rna_type)
{
  ValNodePtr vnp;
  FieldConstraintPtr      fcp;
  RnaQualPtr rq;

  rq = RnaQualNew ();
  rq->field = Rna_field_product;
  rq->type = AsnIoMemCopy (rna_type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
  fcp = FieldConstraintNew ();
  fcp->field = ValNodeNew (NULL);
  fcp->field->choice = FieldType_rna_field;
  fcp->field->data.ptrvalue = rq;
  fcp->string_constraint = StringConstraintNew ();
  fcp->string_constraint->match_text = StringSave ("");
  vnp = ValNodeNew (NULL);
  vnp->choice = ConstraintChoice_field;
  vnp->data.ptrvalue = fcp;
  return vnp;
}


static ValNodePtr DefaultFeatureFieldConstraint (Int2 feat_type)
{
  ValNodePtr         vnp;
  FieldConstraintPtr fcp;
  FeatureFieldPtr    ff;

  ff = FeatureFieldNew ();
  ff->type = feat_type;
  fcp = FieldConstraintNew ();
  fcp->field = ValNodeNew (NULL);
  fcp->field->choice = FieldType_feature_field;
  fcp->field->data.ptrvalue = ff;
  fcp->string_constraint = StringConstraintNew ();
  fcp->string_constraint->match_text = StringSave ("");
  vnp = ValNodeNew (NULL);
  vnp->choice = ConstraintChoice_field;
  vnp->data.ptrvalue = fcp;
  return vnp;
}


NLM_EXTERN void ChangeComplexConstraintFieldType (DialoG d, Uint2 qual_type, ValNodePtr rna_type, Int2 feat_type)
{
  ComplexConstraintDlgPtr dlg;
  ValNodePtr              default_constraint;

  dlg = (ComplexConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  dlg->constraint_type = ComplexConstraintTypeFromFeatureFieldType (qual_type);
  dlg->feat_type = feat_type;
  dlg->rna_type = RnaFeatTypeFree (dlg->rna_type);
  dlg->rna_type = AsnIoMemCopy (rna_type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);

  /* set default constraint type for advanced constraints */
  switch (qual_type) {
    case FieldType_source_qual:
      SetConstraintSetDefaultConstraintType (dlg->constraint_set_dlg, ConstraintChoice_source);
      break;
    case FieldType_cds_gene_prot:
      SetConstraintSetDefaultConstraintType (dlg->constraint_set_dlg, ConstraintChoice_cdsgeneprot_qual);
      break;
    case FieldType_pub:
      SetConstraintSetDefaultConstraintType (dlg->constraint_set_dlg, ConstraintChoice_pub);
      break;
    case FieldType_molinfo_field:
      SetConstraintSetDefaultConstraintType (dlg->constraint_set_dlg, ConstraintChoice_sequence);
      break;
    case FieldType_feature_field:
      default_constraint = DefaultFeatureFieldConstraint (dlg->feat_type);
      SetConstraintSetDefaultConstraintTypeEx (dlg->constraint_set_dlg, default_constraint);
      PointerToDialog (dlg->feat_field_constraint, default_constraint->data.ptrvalue);
      default_constraint = ConstraintChoiceFree (default_constraint);
      break;
    case FieldType_rna_field: 
      default_constraint = DefaultRNAFieldConstraint (dlg->rna_type);
      SetConstraintSetDefaultConstraintTypeEx (dlg->constraint_set_dlg, default_constraint);
      default_constraint = ConstraintChoiceFree (default_constraint);
      break;
    case FieldType_struc_comment_field:
      SetConstraintSetDefaultConstraintType (dlg->constraint_set_dlg, ConstraintChoice_sequence);
      break;
  }    
  /* change which constraint dialog is visible, if using "simple" constraints */
  if (GetValue (dlg->constraint_choice) == 1) {
    ShowSimpleConstraint (dlg);
  }
}


NLM_EXTERN void SetComplexConstraintType (DialoG d, EComplexConstraintType constraint_type)
{
  Uint2 ftype;
  ComplexConstraintDlgPtr dlg;

  ftype = FeatureFieldTypeFromComplexConstraintType (constraint_type);
  if (ftype != 0) {
    ChangeComplexConstraintFieldType (d, ftype, NULL, Macro_feature_type_any);
  } else {
    dlg = (ComplexConstraintDlgPtr) GetObjectExtra (d);
    if (dlg != NULL) {
      dlg->constraint_type = constraint_type;
      /* change which constraint dialog is visible, if using "simple" constraints */
      if (GetValue (dlg->constraint_choice) == 1) {
        ShowSimpleConstraint (dlg);
      }
    }
  }
}


static void ChangeComplexConstraintChoice (GrouP g)
{
  ComplexConstraintDlgPtr dlg;

  dlg = (ComplexConstraintDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  if (GetValue (dlg->constraint_choice) == 1) {
    Hide (dlg->constraint_set_dlg);
    Show (dlg->old_constraint_grp);
    ShowSimpleConstraint (dlg);
  } else {
    Hide (dlg->old_constraint_grp);
    Show (dlg->constraint_set_dlg);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void DataToComplexConstraintDlg (DialoG d, Pointer data)
{
  ComplexConstraintDlgPtr dlg;
  ValNodePtr              constraint = NULL;
  FieldConstraintPtr      fcp;

  dlg = (ComplexConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  constraint = (ValNodePtr) data;

  if (GetValue (dlg->constraint_choice) == 2
      || (constraint != NULL && constraint->next != NULL)) {
    PointerToDialog (dlg->constraint_set_dlg, constraint);
    SetValue (dlg->constraint_choice, 2);
  } else if (constraint == NULL) {
    switch (dlg->constraint_type) {
      case eComplexConstraintType_source:
        PointerToDialog (dlg->src_constraint, NULL);
        break;
      case eComplexConstraintType_cdsgeneprot:
        PointerToDialog (dlg->cgp_constraint, NULL);
        break;
      case eComplexConstraintType_pub:
        PointerToDialog (dlg->pub_constraint, NULL);
        break;
      case eComplexConstraintType_feature_field:
        PointerToDialog (dlg->feat_field_constraint, NULL);
        break;
      case eComplexConstraintType_rna_field:
        PointerToDialog (dlg->rna_constraint, NULL);
        break;
      case eComplexConstraintType_molinfo_field:
        PointerToDialog (dlg->seq_constraint, NULL);
        break;
      case eComplexConstraintType_seqid:
        PointerToDialog (dlg->seq_constraint, NULL);
        break;
      case eComplexConstraintType_string:
        PointerToDialog (dlg->string_constraint, NULL);
        break;
    }
  } else {
    switch (constraint->choice) {
      case ConstraintChoice_field:
        fcp = constraint->data.ptrvalue;
        if (fcp->field == NULL || fcp->field->choice != FeatureFieldTypeFromComplexConstraintType(dlg->constraint_type)) {
          SetValue (dlg->constraint_choice, 2);
          PointerToDialog (dlg->constraint_set_dlg, constraint);
        } else if (dlg->constraint_type == eComplexConstraintType_feature_field) {
          PointerToDialog (dlg->feat_field_constraint, constraint);
        } else if (dlg->constraint_type == eComplexConstraintType_rna_field) {
          PointerToDialog (dlg->rna_constraint, constraint);
        }
        break;
      case ConstraintChoice_source:
        if (dlg->constraint_type == eComplexConstraintType_source) {
          PointerToDialog (dlg->src_constraint, constraint);
        } else {
          SetValue (dlg->constraint_choice, 2);
          PointerToDialog (dlg->constraint_set_dlg, constraint);
        }
        break;
      case ConstraintChoice_cdsgeneprot_qual:
      case ConstraintChoice_cdsgeneprot_pseudo:
        if (dlg->constraint_type == eComplexConstraintType_cdsgeneprot) {
          PointerToDialog (dlg->cgp_constraint, constraint);
        } else {
          SetValue (dlg->constraint_choice, 2);
          PointerToDialog (dlg->constraint_set_dlg, constraint);
        }
        break;
      case ConstraintChoice_sequence:
        SetValue (dlg->constraint_choice, 2);
        PointerToDialog (dlg->constraint_set_dlg, constraint);
        break;
      case ConstraintChoice_pub:
        if (dlg->constraint_type == eComplexConstraintType_pub) {
          PointerToDialog (dlg->pub_constraint, constraint->data.ptrvalue);
        } else {
          SetValue (dlg->constraint_choice, 2);
          PointerToDialog (dlg->constraint_set_dlg, constraint);
        }
        break;
    }
  }
    
  ChangeComplexConstraintChoice(dlg->constraint_choice);
}


static Pointer DataFromComplexConstraintDlg (DialoG d)
{
  ComplexConstraintDlgPtr  dlg;
  ValNodePtr               constraint = NULL;
  FieldConstraintPtr       fc;
  RnaQualPtr               rq;
  PublicationConstraintPtr pub;
  StringConstraintPtr      scp;

  dlg = (ComplexConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (GetValue (dlg->constraint_choice) == 2) {
    constraint = DialogToPointer (dlg->constraint_set_dlg);
  } else {
    switch (dlg->constraint_type) {
      case eComplexConstraintType_source:
        constraint = DialogToPointer (dlg->src_constraint);
        break;
      case eComplexConstraintType_cdsgeneprot:
        constraint = DialogToPointer (dlg->cgp_constraint);
        break;
      case eComplexConstraintType_pub:
        pub = DialogToPointer (dlg->pub_constraint);
        if (!IsPublicationConstraintEmpty (pub)) {
          ValNodeAddPointer (&constraint, ConstraintChoice_pub, pub);
        } else {
          pub = PublicationConstraintFree (pub);
        }
        break;
      case eComplexConstraintType_feature_field:
        constraint = DialogToPointer (dlg->feat_field_constraint);
        break;
      case eComplexConstraintType_rna_field:
        constraint = DialogToPointer (dlg->rna_constraint);
        if (constraint != NULL) {
          fc = (FieldConstraintPtr) constraint->data.ptrvalue;
          if (fc == NULL || fc->field == NULL) {
            constraint = ConstraintChoiceFree (constraint);
          } else if (fc->field->choice == FieldType_rna_field) {
            rq = (RnaQualPtr) fc->field->data.ptrvalue;
            if (rq != NULL) {
              rq->type = AsnIoMemCopy (dlg->rna_type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
            }
          }
        }
        break;
      case eComplexConstraintType_molinfo_field:
        constraint = DialogToPointer (dlg->seq_constraint);
        break;
      case eComplexConstraintType_seqid:
        constraint = DialogToPointer (dlg->seq_constraint);
        break;
      case eComplexConstraintType_string:
        scp = DialogToPointer (dlg->string_constraint);
        if (IsStringConstraintEmpty (scp)) {
          scp = StringConstraintFree(scp);
        } else {
          ValNodeAddPointer (&constraint, ConstraintChoice_string, scp);
        }
        break;
    }
  }
  return constraint;
}


static void ClearComplexConstraintDialogText (DialoG d)
{
  ComplexConstraintDlgPtr dlg;

  dlg = (ComplexConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  ClearOldCDSGeneProtConstraintDialogText (dlg->cgp_constraint);
  ClearOldSrcConstraintDialogText (dlg->src_constraint);
  ClearRnaFieldConstraintDialogText (dlg->rna_constraint);
  ClearPublicationConstraintDialogText (dlg->pub_constraint);
  ClearFeatureFieldConstraintDialogText (dlg->feat_field_constraint);
  ClearSimpleSequenceConstraintDialogText (dlg->seq_constraint);
  ClearStringConstraintDialogText (dlg->string_constraint);
}


static ValNodePtr TestComplexConstraintDlg (DialoG d)
{
  ValNodePtr err_list = NULL;
  ComplexConstraintDlgPtr dlg;
  PublicationConstraintPtr pub;

  dlg = (ComplexConstraintDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (GetValue (dlg->constraint_choice) == 2) {
    err_list = TestDialog (dlg->constraint_set_dlg);
  } else {
    switch (dlg->constraint_type) {
      case eComplexConstraintType_source:
        err_list = TestDialog (dlg->src_constraint);
        break;
      case eComplexConstraintType_cdsgeneprot:
        err_list = TestDialog (dlg->cgp_constraint);
        break;
      case eComplexConstraintType_pub:
        pub = DialogToPointer (dlg->pub_constraint);
        if (!IsPublicationConstraintEmpty (pub)) {
          err_list = TestDialog (dlg->pub_constraint);
        }
        pub = PublicationConstraintFree (pub);
        break;
      case eComplexConstraintType_feature_field:
        err_list = TestDialog (dlg->feat_field_constraint);
        break;
      case eComplexConstraintType_rna_field:
        err_list = TestDialog (dlg->rna_constraint);
        break;
      case eComplexConstraintType_molinfo_field:
        err_list = TestDialog (dlg->seq_constraint);
        break;
      case eComplexConstraintType_seqid:
        err_list = TestDialog (dlg->seq_constraint);
        break;
      case eComplexConstraintType_string:
        err_list = TestDialog (dlg->string_constraint);
        break;
    }
  }
  return err_list;
}


static void CleanupComplexConstraintDialog (GraphiC g, VoidPtr data)
{
  ComplexConstraintDlgPtr dlg;

  dlg = (ComplexConstraintDlgPtr) data;
  if (dlg != NULL) {
    dlg->rna_type = RnaFeatTypeFree (dlg->rna_type);
  }
  StdCleanupExtraProc (g, data);
}


NLM_EXTERN DialoG ComplexConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  GrouP p, constraint_grp;
  ComplexConstraintDlgPtr dlg;

  dlg = (ComplexConstraintDlgPtr) MemNew (sizeof (ComplexConstraintDlgData));

  p = NormalGroup (h, -1, 0, "Constraints", programFont, NULL);
  SetObjectExtra (p, dlg, CleanupComplexConstraintDialog);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToComplexConstraintDlg;
  dlg->fromdialog = DataFromComplexConstraintDlg;
  dlg->testdialog = TestComplexConstraintDlg;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  constraint_grp = HiddenGroup (p, 0, 0, NULL);
  dlg->old_constraint_grp = HiddenGroup (constraint_grp, 0, 0, NULL);
  dlg->cgp_constraint = OldCDSGeneProtConstraintDialog (dlg->old_constraint_grp, dlg->change_notify, dlg->change_userdata);
  dlg->src_constraint = OldSrcConstraintDialog (dlg->old_constraint_grp, dlg->change_notify, dlg->change_userdata);
  dlg->pub_constraint = PublicationConstraintDialog (dlg->old_constraint_grp, dlg->change_notify, dlg->change_userdata);
  dlg->feat_field_constraint = FeatureFieldConstraintDialog (dlg->old_constraint_grp, dlg->change_notify, dlg->change_userdata);
  dlg->rna_constraint = RNAFieldConstraintDialog (dlg->old_constraint_grp, NULL, dlg->change_notify, dlg->change_userdata);
  dlg->seq_constraint = SimpleSequenceConstraintDialog (dlg->old_constraint_grp, dlg->change_notify, dlg->change_userdata);
  dlg->string_constraint = StringConstraintDialog (dlg->old_constraint_grp, "Where object text", TRUE, dlg->change_notify, dlg->change_userdata);
  AlignObjects (ALIGN_CENTER,
                (HANDLE) dlg->cgp_constraint,
                (HANDLE) dlg->src_constraint, 
                (HANDLE) dlg->pub_constraint, 
                (HANDLE) dlg->feat_field_constraint, 
                (HANDLE) dlg->rna_constraint, 
                (HANDLE) dlg->seq_constraint,
                (HANDLE) dlg->string_constraint,
                NULL);

  dlg->constraint_set_dlg = ConstraintSetDialog (constraint_grp, dlg->change_notify, dlg->change_userdata);
  Hide (dlg->constraint_set_dlg);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->constraint_set_dlg, (HANDLE) dlg->old_constraint_grp, NULL);
  dlg->constraint_choice = HiddenGroup (p, 2, 0, ChangeComplexConstraintChoice);
  SetObjectExtra (dlg->constraint_choice, dlg, NULL);
  RadioButton (dlg->constraint_choice, "Simple");
  RadioButton (dlg->constraint_choice, "Advanced");
  SetValue (dlg->constraint_choice, 1);
  AlignObjects (ALIGN_CENTER, (HANDLE) constraint_grp, (HANDLE) dlg->constraint_choice, NULL);
  return (DialoG) p;
}


static DialoG LegalFeatQualChoiceDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ValNodePtr field_list = NULL;

  AddAllFeatureFieldsToChoiceList (&field_list);
  return ValNodeSelectionDialog (h, field_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type",
                                change_notify, change_userdata, FALSE);
}


typedef struct featqualchoicedlg {
  DIALOG_MESSAGE_BLOCK
  DialoG legal_qual;
  DialoG illegal_qual;
  Boolean is_legal;
} FeatQualChoiceDlgData, PNTR FeatQualChoiceDlgPtr;


static void FeatQualChoiceToDialog (DialoG d, Pointer data)
{
  FeatQualChoiceDlgPtr dlg;
  ValNodePtr vnp;
  ValNode vn;

  dlg = (FeatQualChoiceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  vnp = (ValNodePtr) data;
  if (vnp == NULL) {
    dlg->is_legal = TRUE;
    PointerToDialog (dlg->legal_qual, NULL);
    Show (dlg->legal_qual);
    SafeHide (dlg->illegal_qual);
  } else if (vnp->choice == FeatQualChoice_legal_qual) {
    dlg->is_legal = TRUE;
    vn.choice = vnp->data.intvalue;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->legal_qual, &vn);
    Show (dlg->legal_qual);
    SafeHide (dlg->illegal_qual);
  } else if (vnp->choice == FeatQualChoice_illegal_qual) {
    dlg->is_legal = FALSE;
    PointerToDialog (dlg->illegal_qual, vnp->data.ptrvalue);
    Hide (dlg->legal_qual);
    SafeShow (dlg->illegal_qual);
  } else {
    dlg->is_legal = TRUE;
    PointerToDialog (dlg->legal_qual, NULL);
    Show (dlg->legal_qual);
    SafeHide (dlg->illegal_qual);
  }
}


static Pointer DialogToFeatQualChoice (DialoG d)
{
  FeatQualChoiceDlgPtr dlg;
  ValNodePtr vnp = NULL, vnp2;
  StringConstraintPtr scp;

  dlg = (FeatQualChoiceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->is_legal) {
    vnp2 = (ValNodePtr) DialogToPointer (dlg->legal_qual);
    if (vnp2 != NULL) {
      vnp = ValNodeNew (NULL);
      vnp->choice = FeatQualChoice_legal_qual;
      vnp->data.intvalue = vnp2->choice;
      vnp2 = ValNodeFree (vnp2);
    }
  } else {
    scp = (StringConstraintPtr) DialogToPointer (dlg->illegal_qual);
    if (scp != NULL) {
      vnp = ValNodeNew (NULL);
      vnp->choice = FeatQualChoice_illegal_qual;
      vnp->data.ptrvalue = scp;
    }
  }
  return vnp;
}


static ValNodePtr TestFeatQualChoiceDialog (DialoG d)
{
  FeatQualChoiceDlgPtr dlg;
  ValNodePtr err_list = NULL;
  StringConstraintPtr scp;

  dlg = (FeatQualChoiceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  if (dlg->is_legal) {
    err_list = TestDialog (dlg->legal_qual);
  } else {
    scp = DialogToPointer (dlg->illegal_qual);
    if (scp == NULL || StringHasNoText (scp->match_text)) {
      ValNodeAddPointer (&err_list, 0, "match text");
    }
  }
  return err_list;
}

static DialoG FeatQualChoiceDialog (GrouP h, Boolean allow_illegal, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  FeatQualChoiceDlgPtr dlg;
  GrouP                 p;
  
  dlg = (FeatQualChoiceDlgPtr) MemNew (sizeof (FeatQualChoiceDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 0, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FeatQualChoiceToDialog;
  dlg->fromdialog = DialogToFeatQualChoice;
  dlg->testdialog = TestFeatQualChoiceDialog;
  dlg->is_legal = TRUE;

  dlg->legal_qual = LegalFeatQualChoiceDialog (p, change_notify, change_userdata);
  if (allow_illegal) {
    dlg->illegal_qual = StringConstraintDialog (p, NULL, FALSE, change_notify, change_userdata);
    Hide (dlg->illegal_qual);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->legal_qual, (HANDLE) dlg->illegal_qual, NULL);

  return (DialoG) p;
}


typedef struct cdsgeneprotfielddlg {
  DIALOG_MESSAGE_BLOCK
  DialoG dlg;
} CDSGeneProtFieldDlgData, PNTR CDSGeneProtFieldDlgPtr;


static void CDSGeneProtFieldToDialog (DialoG d, Pointer data)
{
  CDSGeneProtFieldDlgPtr dlg;
  ValNodePtr vnp;
  ValNode vn;
  
  dlg = (CDSGeneProtFieldDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    vnp = (ValNodePtr) data;
    if (vnp == NULL) {
      PointerToDialog (dlg->dlg, NULL);
    } else {
      vn.choice = vnp->data.intvalue;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->dlg, &vn);
    }
  }
}

static Pointer DialogToCDSGeneProtField (DialoG d)
{
  CDSGeneProtFieldDlgPtr dlg;
  ValNodePtr vnp, field = NULL;
  
  dlg = (CDSGeneProtFieldDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    vnp = (ValNodePtr) DialogToPointer (dlg->dlg);
    if (vnp != NULL) {
      field = ValNodeNew (NULL);
      field->choice = FieldType_cds_gene_prot;
      field->data.intvalue = vnp->choice;
    }
  }
  return field;
}


static ValNodePtr TestCDSGeneProtFieldDialog (DialoG d)
{
  CDSGeneProtFieldDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (CDSGeneProtFieldDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    err_list = TestDialog (dlg->dlg);
  }
  return err_list;
}


static DialoG CDSGeneProtFieldDialogEx (GrouP h, Boolean tall, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  CDSGeneProtFieldDlgPtr dlg;
  GrouP                  p;
  ValNodePtr             field_list = NULL;
  
  dlg = (CDSGeneProtFieldDlgPtr) MemNew (sizeof (CDSGeneProtFieldDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 0, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = CDSGeneProtFieldToDialog;
  dlg->fromdialog = DialogToCDSGeneProtField;
  dlg->testdialog = TestCDSGeneProtFieldDialog;

  AddAllCDSGeneProtFieldsToChoiceList (&field_list);
  dlg->dlg = ValNodeSelectionDialog (p, field_list, tall ? TALL_SELECTION_LIST : SHORT_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type",
                                change_notify, change_userdata, FALSE);

  return (DialoG) p;
}


static DialoG CDSGeneProtFieldDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return CDSGeneProtFieldDialogEx (h, TRUE, change_notify, change_userdata);
}


typedef struct rnafeaturetypedlg {
  DIALOG_MESSAGE_BLOCK
  DialoG rna_type_dlg;
  GrouP  ncrna_class_grp;
  DialoG ncrna_class_dlg;
  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} RnaFeatureTypeDlgData, PNTR RnaFeatureTypeDlgPtr;



static void ChangeRnaFeatureType (Pointer data)
{
  RnaFeatureTypeDlgPtr dlg;
  ValNodePtr           vnp;

  dlg = (RnaFeatureTypeDlgPtr) data;
  if (dlg == NULL) return;

  vnp = DialogToPointer (dlg->rna_type_dlg);
  if (vnp == NULL || vnp->choice != RnaFeatType_ncRNA) {
    Hide (dlg->ncrna_class_grp);
  } else {
    Show (dlg->ncrna_class_grp);
  }
  vnp = ValNodeFree (vnp);
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void RnaFeatureTypeToDialog (DialoG d, Pointer data)
{
  RnaFeatureTypeDlgPtr dlg;
  RnaFeatTypePtr r;
  ValNode        vn;

  dlg = (RnaFeatureTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  r = (RnaFeatTypePtr) data;

  if (r == NULL) {
    vn.choice = RnaFeatType_rRNA;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->rna_type_dlg, &vn);
    PointerToDialog (dlg->ncrna_class_dlg, NULL);
  } else {
    vn.choice = r->choice;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->rna_type_dlg, &vn);
    if (r->choice == RnaFeatType_ncRNA) {
      PointerToDialog (dlg->ncrna_class_dlg, r->data.ptrvalue);
    } else {
      PointerToDialog (dlg->ncrna_class_dlg, NULL);
    }
  }
  ChangeRnaFeatureType (dlg);
}


static Pointer DialogToRnaFeatureType (DialoG d)
{
  RnaFeatureTypeDlgPtr dlg;
  RnaFeatTypePtr r = NULL;
  ValNodePtr     vnp;

  dlg = (RnaFeatureTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  vnp = DialogToPointer (dlg->rna_type_dlg);
  if (vnp != NULL) {
    r = ValNodeNew (NULL);
    r->choice = vnp->choice;
    if (r->choice == RnaFeatType_ncRNA) {
      r->data.ptrvalue = DialogToPointer (dlg->ncrna_class_dlg);
    }
    vnp = ValNodeFree (vnp);
  }
  return r;
}


static ValNodePtr TestRnaFeatureTypeDialog (DialoG d)
{
  RnaFeatureTypeDlgPtr dlg;
  ValNodePtr     vnp, err_list = NULL;

  dlg = (RnaFeatureTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  vnp = DialogToPointer (dlg->rna_type_dlg);
  if (vnp == NULL) {
    ValNodeAddPointer (&err_list, 0, "no type");
  }
  vnp = ValNodeFree (vnp);
  return err_list;
}

  
static DialoG RnaFeatureTypeDialog (GrouP h, CharPtr type_label, Boolean allow_any, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  RnaFeatureTypeDlgPtr dlg;
  GrouP                  p;
  ValNodePtr             type_list = NULL;
  
  dlg = (RnaFeatureTypeDlgPtr) MemNew (sizeof (RnaFeatureTypeDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = NormalGroup (h, -1, 0, type_label, programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = RnaFeatureTypeToDialog;
  dlg->fromdialog = DialogToRnaFeatureType;
  dlg->testdialog = TestRnaFeatureTypeDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  if (allow_any) {
    ValNodeAddPointer (&type_list, RnaFeatType_any, StringSave ("any"));
  }
  ValNodeLink (&type_list, GetRNATypeList ());
  dlg->rna_type_dlg = ValNodeSelectionDialog (p, type_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type",
                                ChangeRnaFeatureType, dlg, FALSE);

  dlg->ncrna_class_grp = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (dlg->ncrna_class_grp, "ncRNA class", 0, dialogTextHeight, programFont, 'r');
  dlg->ncrna_class_dlg = CreatencRNAClassDialog (dlg->ncrna_class_grp, TRUE, change_notify, change_userdata);
  Hide (dlg->ncrna_class_grp);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->rna_type_dlg, (HANDLE) dlg->ncrna_class_grp, NULL);

  return (DialoG) p;

}


typedef struct rnaqualdlg {
  DIALOG_MESSAGE_BLOCK
  DialoG rnafeat_dlg;
  DialoG field_dlg;
  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} RnaQualDlgData, PNTR RnaQualDlgPtr;


static void RnaQualToDialog (DialoG d, Pointer data)
{
  RnaQualDlgPtr dlg;
  RnaQualPtr    q;
  ValNode       vn;

  dlg = (RnaQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  q = (RnaQualPtr) data;
  if (q == NULL) {
    PointerToDialog (dlg->rnafeat_dlg, NULL);
    PointerToDialog (dlg->field_dlg, NULL);
  } else {
    PointerToDialog (dlg->rnafeat_dlg, q->type);
    vn.choice = (Uint1)q->field;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->field_dlg, &vn);
  }
}


static Pointer DialogToRnaQual (DialoG d)
{
  RnaQualDlgPtr dlg;
  RnaQualPtr    q;
  ValNodePtr    vnp;

  dlg = (RnaQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  q = RnaQualNew ();
  q->type = DialogToPointer (dlg->rnafeat_dlg);
  vnp = DialogToPointer (dlg->field_dlg);
  if (vnp != NULL) {
    q->field = vnp->choice;
  }
  vnp = ValNodeFree (vnp);
  return (Pointer) q;
}


static ValNodePtr TestRnaQualDialog (DialoG d)
{
  ValNodePtr err_list = NULL;
  RnaQualDlgPtr dlg;

  dlg = (RnaQualDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  err_list = TestDialog (dlg->rnafeat_dlg);
  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  return err_list;
}


static DialoG RnaQualDialogEx (GrouP h, CharPtr type_label, Boolean allow_any, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  RnaQualDlgPtr dlg;
  GrouP          p;
  ValNodePtr     field_list = NULL;
  ValNode        vn;
  PrompT         ppt = NULL;
  
  dlg = (RnaQualDlgPtr) MemNew (sizeof (RnaQualDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = RnaQualToDialog;
  dlg->fromdialog = DialogToRnaQual;
  dlg->testdialog = TestRnaQualDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  if (type_label != NULL) {
    dlg->rnafeat_dlg = RnaFeatureTypeDialog (p, type_label, allow_any, change_notify, change_userdata);
  }
  ppt = StaticPrompt (p, "RNA Field", 0, dialogTextHeight, programFont, 'l');
  field_list = GetRnaFieldList ();
  dlg->field_dlg = ValNodeSelectionDialog (p, field_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type",
                                change_notify, change_userdata, FALSE);
  vn.choice = field_list->choice;
  vn.data.ptrvalue = NULL;
  vn.next = NULL;
  PointerToDialog (dlg->field_dlg, &vn);

  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->field_dlg, (HANDLE) dlg->rnafeat_dlg, NULL);
  return (DialoG) p;
}


static DialoG RnaQualDialog (GrouP h, CharPtr type_label, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  return RnaQualDialogEx (h, type_label, FALSE, change_notify, change_userdata);
}


typedef struct structuredcommentfielddlg {
  DIALOG_MESSAGE_BLOCK
  PopuP field_type;
  TexT  field_name;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} StructuredCommentFieldDlgData, PNTR StructuredCommentFieldDlgPtr;


static void ChangeStructureCommentFieldChoice (PopuP p)
{
  Int2 val;
  StructuredCommentFieldDlgPtr dlg;

  dlg = (StructuredCommentFieldDlgPtr) GetObjectExtra (p);

  if (dlg == NULL) {
    return;
  }

  val = GetValue (dlg->field_type);
  switch (val) {
    case 2:
      Hide (dlg->field_name);
      break;
    case 1:
      Show (dlg->field_name);
      break;
    case 3:
      Hide (dlg->field_name);
      break;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void StructuredCommentFieldToDialog (DialoG d, Pointer data)
{
  StructuredCommentFieldDlgPtr dlg;
  StructuredCommentFieldPtr    field;

  dlg = (StructuredCommentFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  field = (StructuredCommentFieldPtr) data;
  if (field == NULL) {
    SetValue (dlg->field_type, 1);
  } else {
    switch (field->choice) {
      case StructuredCommentField_database:
        SetValue (dlg->field_type, 2);
        break;
      case StructuredCommentField_named:
        SetValue (dlg->field_type, 1);
        SetTitle (dlg->field_name, field->data.ptrvalue);
        break;
      case StructuredCommentField_field_name:
        SetValue (dlg->field_type, 3);
        break;
    }
  }

  ChangeStructureCommentFieldChoice (dlg->field_type);
}


static Pointer StructuredCommentFieldFromDialog (DialoG d)
{
  StructuredCommentFieldDlgPtr dlg;
  StructuredCommentFieldPtr    field;
  Int2                         val;

  dlg = (StructuredCommentFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  field = ValNodeNew (NULL);
  val = GetValue (dlg->field_type);

  switch (val) {
    case 2:
      field->choice = StructuredCommentField_database;
      break;
    case 1:
      field->choice = StructuredCommentField_named;
      field->data.ptrvalue = SaveStringFromText (dlg->field_name);
      if (StringHasNoText (field->data.ptrvalue)) {
        field = ValNodeFree (field);
      }
      break;
    case 3:
      field->choice = StructuredCommentField_field_name;
      break;
    default:
      field = ValNodeFree (field);
      break;
  }
  return (Pointer) field;
}


static ValNodePtr TestStructuredCommentFieldDialog (DialoG d)
{
  ValNodePtr field, err_list = NULL;

  field = DialogToPointer (d);
  if (field == NULL) {
    ValNodeAddPointer (&err_list, 0, "bad field");
  } else {
    field = StructuredCommentFieldFree (field);
  }
  return err_list;
}


static void ChangeStructureCommentFieldDialogText (TexT t)
{
  StructuredCommentFieldDlgPtr dlg;

  dlg = (StructuredCommentFieldDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) { 
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG StructuredCommentFieldDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  StructuredCommentFieldDlgPtr dlg;
  GrouP                 p;
  
  dlg = (StructuredCommentFieldDlgPtr) MemNew (sizeof (StructuredCommentFieldDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = StructuredCommentFieldToDialog;
  dlg->fromdialog = StructuredCommentFieldFromDialog;
  dlg->testdialog = TestStructuredCommentFieldDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->field_type = PopupList (p, TRUE, ChangeStructureCommentFieldChoice);
  SetObjectExtra (dlg->field_type, dlg, NULL);
  PopupItem (dlg->field_type, "Field");
  PopupItem (dlg->field_type, "Database Name");
  PopupItem (dlg->field_type, "Field Name");
  SetValue (dlg->field_type, 1);

  dlg->field_name = DialogText (p, "", 20, ChangeStructureCommentFieldDialogText);
  SetObjectExtra (dlg->field_name, dlg, NULL);

  return (DialoG) p;
}


typedef struct miscfielddlg {
  DIALOG_MESSAGE_BLOCK
  PopuP field_type;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} MiscFieldDlgData, PNTR MiscFieldDlgPtr;


static void ChangeMiscFieldChoice (PopuP p)
{
  MiscFieldDlgPtr dlg;

  dlg = (MiscFieldDlgPtr) GetObjectExtra (p);

  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void MiscFieldToDialog (DialoG d, Pointer data)
{
  MiscFieldDlgPtr dlg;
  ValNodePtr    field;

  dlg = (MiscFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  field = (ValNodePtr) data;
  if (field == NULL) {
    SetValue (dlg->field_type, 1);
  } else {
    switch (field->data.intvalue) {
      case Misc_field_genome_project_id:
        SetValue (dlg->field_type, 1);
        break;
      case Misc_field_comment_descriptor:
        SetValue (dlg->field_type, 2);
        break;
      case Misc_field_defline:
        SetValue (dlg->field_type, 3);
        break;
      case Misc_field_keyword:
        SetValue (dlg->field_type, 4);
        break;
    }
  }

  ChangeMiscFieldChoice (dlg->field_type);
}


static Pointer MiscFieldFromDialog (DialoG d)
{
  MiscFieldDlgPtr dlg;
  ValNodePtr      field;
  Int2            val;

  dlg = (MiscFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  field = ValNodeNew (NULL);
  field->choice = FieldType_misc;
  val = GetValue (dlg->field_type);

  switch (val) {
    case 1:
      field->data.intvalue = Misc_field_genome_project_id;
      break;
    case 2:
      field->data.intvalue = Misc_field_comment_descriptor;
      break;
    case 3:
      field->data.intvalue = Misc_field_defline;
      break;
    case 4:
      field->data.intvalue = Misc_field_keyword;
      break;
    default:
      field = ValNodeFree (field);
      break;
  }
  return (Pointer) field;
}


static ValNodePtr TestMiscFieldDialog (DialoG d)
{
  ValNodePtr field, err_list = NULL;

  field = DialogToPointer (d);
  if (field == NULL) {
    ValNodeAddPointer (&err_list, 0, "bad field");
  } else {
    field = ValNodeFree (field);
  }
  return err_list;
}


static DialoG MiscFieldDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  MiscFieldDlgPtr dlg;
  GrouP           p;
  
  dlg = (MiscFieldDlgPtr) MemNew (sizeof (MiscFieldDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = MiscFieldToDialog;
  dlg->fromdialog = MiscFieldFromDialog;
  dlg->testdialog = TestMiscFieldDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->field_type = PopupList (p, TRUE, ChangeMiscFieldChoice);
  SetObjectExtra (dlg->field_type, dlg, NULL);
  PopupItem (dlg->field_type, "Genome Project ID");
  PopupItem (dlg->field_type, "Comment Descriptor");
  PopupItem (dlg->field_type, "Definition Line");
  PopupItem (dlg->field_type, "Keyword");
  SetValue (dlg->field_type, 1);

  return (DialoG) p;
}


typedef struct fieldtypedlg {
  DIALOG_MESSAGE_BLOCK
  Uint1 field_type;
  DialoG src_qual;
  GrouP feature_field_grp;
  DialoG feature_type;
  DialoG feature_field;
  DialoG cdsgeneprot;
  DialoG sequence_qual;
  DialoG pub_field;
  DialoG rna_field;
  DialoG structured_comment_field;
  DialoG misc_field;
  PopuP  dblink_field;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} FieldTypeDlgData, PNTR FieldTypeDlgPtr;


static void FieldTypeToDialog (DialoG d, Pointer data)
{
  FieldTypeDlgPtr dlg;
  FeatureFieldPtr ffp;
  ValNodePtr vnp;
  ValNode vn;

  dlg = (FieldTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  SafeHide (dlg->src_qual);
  SafeHide (dlg->feature_field_grp);
  SafeHide (dlg->cdsgeneprot);
  SafeHide (dlg->sequence_qual);
  SafeHide (dlg->pub_field);
  SafeHide (dlg->rna_field);
  SafeHide (dlg->structured_comment_field);
  SafeHide (dlg->misc_field);
  SafeHide (dlg->dblink_field);

  vnp = (ValNodePtr) data;
  if (vnp != NULL) {
    switch (vnp->choice) {
      case FieldType_source_qual:
        SafeShow (dlg->src_qual);
        PointerToDialog (dlg->src_qual, vnp->data.ptrvalue);
        dlg->field_type = FieldType_source_qual;
        break;
      case FieldType_feature_field:
        SafeShow (dlg->feature_field_grp);
        ffp = (FeatureFieldPtr) vnp->data.ptrvalue;
        if (ffp == NULL) {
          PointerToDialog (dlg->feature_type, NULL);
          PointerToDialog (dlg->feature_field, NULL);
        } else {
          vn.choice = (Uint1)ffp->type;
          vn.data.ptrvalue = NULL;
          vn.next = NULL;
          PointerToDialog (dlg->feature_type, &vn);
          PointerToDialog (dlg->feature_field, ffp->field);
        }
        dlg->field_type = FieldType_feature_field;
        break;
      case FieldType_cds_gene_prot:
        SafeShow (dlg->cdsgeneprot);
        PointerToDialog (dlg->cdsgeneprot, vnp);
        dlg->field_type = FieldType_cds_gene_prot;
        break;
      case FieldType_molinfo_field:
        SafeShow (dlg->sequence_qual);
        PointerToDialog (dlg->sequence_qual, vnp->data.ptrvalue);
        dlg->field_type = FieldType_molinfo_field;
        break;
      case FieldType_pub:
        SafeShow (dlg->pub_field);
        vn.choice = vnp->data.intvalue;
        vn.data.ptrvalue = NULL;
        vn.next = NULL;
        PointerToDialog (dlg->pub_field, &vn);
        dlg->field_type = FieldType_pub;
        break;
      case FieldType_rna_field:
        SafeShow (dlg->rna_field);
        PointerToDialog (dlg->rna_field, vnp->data.ptrvalue);
        dlg->field_type = FieldType_rna_field;
        break;
      case FieldType_struc_comment_field:
        SafeShow (dlg->structured_comment_field);
        PointerToDialog (dlg->structured_comment_field, vnp->data.ptrvalue);
        dlg->field_type = FieldType_struc_comment_field;
        break;
      case FieldType_misc:
        SafeShow (dlg->misc_field);
        PointerToDialog (dlg->misc_field, vnp);
        dlg->field_type = FieldType_misc;
        break;
      case FieldType_dblink:
        SafeShow (dlg->dblink_field);
        SetValue (dlg->dblink_field, vnp->data.intvalue);
        dlg->field_type = FieldType_dblink;
        break;
    }
  }
}


static ValNodePtr GetFieldOfTypeFromFieldType(DialoG d, Uint1 field_type)
{
  FieldTypeDlgPtr dlg;
  FeatureFieldPtr ffp;
  ValNodePtr vnp = NULL, vnp2;

  dlg = (FieldTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  switch (field_type) {
    case FieldType_source_qual:
      vnp = ValNodeNew(NULL);
      vnp->choice = FieldType_source_qual;
      vnp->data.ptrvalue = DialogToPointer (dlg->src_qual);
      break;
    case FieldType_feature_field:
      vnp2 = DialogToPointer (dlg->feature_type);
      if (vnp2 != NULL) {
        ffp = FeatureFieldNew();
        ffp->type = vnp2->choice;
        ffp->field = DialogToPointer (dlg->feature_field);
        vnp = ValNodeNew (NULL);
        vnp->choice = FieldType_feature_field;
        vnp->data.ptrvalue = ffp;
      }
      break;
    case FieldType_cds_gene_prot:
      vnp = DialogToPointer (dlg->cdsgeneprot);
      break;
    case FieldType_molinfo_field:
      vnp = ValNodeNew (NULL);
      vnp->choice = FieldType_molinfo_field;
      vnp->data.ptrvalue = DialogToPointer (dlg->sequence_qual);
      break;
    case FieldType_pub:
      vnp2 = DialogToPointer (dlg->pub_field);
      if (vnp2 != NULL) {
        vnp = ValNodeNew (NULL);
        vnp->choice = FieldType_pub;
        vnp->data.intvalue = vnp2->choice;
      }
      break;
    case FieldType_rna_field:
      vnp = ValNodeNew (NULL);
      vnp->choice = FieldType_rna_field;
      vnp->data.ptrvalue = DialogToPointer (dlg->rna_field);
      break;
    case FieldType_struc_comment_field:
      vnp = ValNodeNew (NULL);
      vnp->choice = FieldType_struc_comment_field;
      vnp->data.ptrvalue = DialogToPointer (dlg->structured_comment_field);
      break;
    case FieldType_misc:
      vnp = DialogToPointer (dlg->misc_field);
      break;
    case FieldType_dblink:
      vnp = ValNodeNew (NULL);
      vnp->choice = FieldType_dblink;
      vnp->data.intvalue = GetValue (dlg->dblink_field);
      break;
  }
  return vnp;
}


static Pointer DialogToFieldType (DialoG d)
{
  FieldTypeDlgPtr dlg;

  dlg = (FieldTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  return GetFieldOfTypeFromFieldType(d, dlg->field_type);
}


static ValNodePtr TestFieldTypeDialog (DialoG d)
{
  FieldTypeDlgPtr dlg;
  ValNodePtr      err_list = NULL;

  dlg = (FieldTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  switch (dlg->field_type) {
    case FieldType_source_qual:
      ValNodeLink (&err_list, TestDialog (dlg->src_qual));
      break;
    case FieldType_feature_field:
      ValNodeLink (&err_list, TestDialog (dlg->feature_type));
      ValNodeLink (&err_list, TestDialog (dlg->feature_field));
      break;
    case FieldType_cds_gene_prot:
      ValNodeLink (&err_list, TestDialog (dlg->cdsgeneprot));
      break;
    case FieldType_molinfo_field:
      ValNodeLink (&err_list, TestDialog (dlg->sequence_qual));
      break;
    case FieldType_pub:
      ValNodeLink (&err_list, TestDialog (dlg->pub_field));
      break;
    case FieldType_rna_field:
      ValNodeLink (&err_list, TestDialog (dlg->rna_field));
      break;
    case FieldType_struc_comment_field:
      ValNodeLink (&err_list, TestDialog (dlg->structured_comment_field));
      break;
    case FieldType_misc:
      ValNodeLink (&err_list, TestDialog (dlg->misc_field));
      break;
    case FieldType_dblink:
      break;
    default :
      ValNodeAddPointer (&err_list, 0, "No field type chosen");
      break;
  }
  return err_list;  
}


static void ChangeFieldTypeDialogPopup (PopuP p)
{
  FieldTypeDlgPtr dlg;

  dlg = (FieldTypeDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static PopuP MakeDbLinkFieldPopup (GrouP h, Nlm_PupActnProc actn, Pointer extradata)
{
  Int4            num, i;
  PopuP           p;

  p = PopupList (h, TRUE, actn);
  SetObjectExtra (p, extradata, NULL);
  num = GetNumDBLinkFields();
  for (i = 1; i <= num; i++) {
    PopupItem (p, GetDBLinkNameFromDBLinkFieldType(i));
  }
  SetValue (p, 1);
  return p;
}


static void FieldTypeDialogMessage (DialoG d, Int2 mssg)

{
  FieldTypeDlgPtr dlg;

  dlg = (FieldTypeDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case VIB_MSG_ENTER :
        switch (dlg->field_type) {
          case FieldType_source_qual:
            SendMessageToDialog (dlg->src_qual, VIB_MSG_ENTER);
            break;
          case FieldType_feature_field:
            Select (dlg->feature_field);
            break;
          case FieldType_cds_gene_prot:
            Select (dlg->cdsgeneprot);
            break;
          case FieldType_molinfo_field:
            Select (dlg->sequence_qual);
            break;
          case FieldType_pub:
            Select (dlg->pub_field);
            break;
          case FieldType_rna_field:
            Select (dlg->rna_field);
            break;
          case FieldType_struc_comment_field:
            Select (dlg->structured_comment_field);
            break;
          case FieldType_misc:
            Select (dlg->misc_field);
            break;
          case FieldType_dblink:
            Select (dlg->dblink_field);
            break;
        }
        break;
      default :
        break;
    }
  }
}


static DialoG FieldTypeDialog (GrouP h, Boolean text_only, Boolean for_remove, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  FieldTypeDlgPtr dlg;
  GrouP           p;
  
  dlg = (FieldTypeDlgPtr) MemNew (sizeof (FieldTypeDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 0, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FieldTypeToDialog;
  dlg->fromdialog = DialogToFieldType;
  dlg->testdialog = TestFieldTypeDialog;
  dlg->dialogmessage = FieldTypeDialogMessage;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->src_qual = SourceQualChoiceDialog (p, text_only, for_remove, TRUE, change_notify, change_userdata);
  dlg->feature_field_grp = HiddenGroup (p, 2, 0, NULL);
  dlg->feature_type = FeatureTypeDialog (dlg->feature_field_grp, change_notify, change_userdata);
  dlg->feature_field = FeatQualChoiceDialog (dlg->feature_field_grp, for_remove, change_notify, change_userdata);
  dlg->cdsgeneprot = CDSGeneProtFieldDialog (p, change_notify, change_userdata);
  dlg->sequence_qual = SequenceQualDialog (p, change_notify, change_userdata);
  dlg->pub_field = PubFieldDialog (p, change_notify, change_userdata);
  dlg->rna_field = RnaQualDialogEx (p, "RNA Type of Feature to be Edited", TRUE, change_notify, change_userdata);
  dlg->structured_comment_field = StructuredCommentFieldDialog (p, change_notify, change_userdata);
  dlg->misc_field = MiscFieldDialog (p, change_notify, change_userdata);
  dlg->dblink_field = MakeDbLinkFieldPopup (p, ChangeFieldTypeDialogPopup, dlg);

  SafeHide (dlg->src_qual);
  SafeHide (dlg->feature_field_grp);
  SafeHide (dlg->cdsgeneprot);
  SafeHide (dlg->sequence_qual);
  SafeHide (dlg->pub_field);
  SafeHide (dlg->rna_field);
  SafeHide (dlg->structured_comment_field);
  SafeHide (dlg->misc_field);
  SafeHide (dlg->dblink_field);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->src_qual,
                              (HANDLE) dlg->feature_field_grp,
                              (HANDLE) dlg->cdsgeneprot,
                              (HANDLE) dlg->sequence_qual,
                              (HANDLE) dlg->pub_field,
                              (HANDLE) dlg->rna_field,
                              (HANDLE) dlg->structured_comment_field,
                              (HANDLE) dlg->misc_field,
                              (HANDLE) dlg->dblink_field,
                              NULL);

  return (DialoG) p;
}


typedef struct fieldpairtypedlg {
  DIALOG_MESSAGE_BLOCK
  Uint1 field_type;
  GrouP src_qual_grp;
  DialoG src_qual_from;
  DialoG src_qual_to;
  GrouP feature_field_grp;
  DialoG feature_type;
  DialoG feature_field_from;
  DialoG feature_field_to;
  GrouP cdsgeneprot_grp;
  DialoG cdsgeneprot_from;
  DialoG cdsgeneprot_to;
  DialoG molinfo_dlg;
  GrouP  rna_grp;
  DialoG rna_feat_dlg;
  DialoG rna_field_from;
  DialoG rna_field_to;
  GrouP  dblink_grp;
  PopuP  dblink_field_from;
  PopuP  dblink_field_to;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} FieldPairTypeDlgData, PNTR FieldPairTypeDlgPtr;


static void FieldPairTypeToDialog (DialoG d, Pointer data)
{
  FieldPairTypeDlgPtr dlg;
  ValNodePtr vnp;
  ValNode      vn;
  FieldTypePtr from_field, to_field;  
  FeatureFieldPairPtr ffp;
  CDSGeneProtFieldPairPtr cgp;
  RnaQualPairPtr          rna_quals;
  DBLinkFieldPairPtr      dblink;

  dlg = (FieldPairTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  SafeHide (dlg->src_qual_grp);
  SafeHide (dlg->feature_field_grp);
  SafeHide (dlg->cdsgeneprot_grp);
  SafeHide (dlg->molinfo_dlg);
  SafeHide (dlg->rna_grp);
  SafeHide (dlg->dblink_grp);
  vnp = (ValNodePtr) data;
  if (vnp == NULL) {
    dlg->field_type = 0;
  } else {
    switch (vnp->choice) {
      case FieldPairType_source_qual:
        SafeShow (dlg->src_qual_grp);
        from_field = GetFromFieldFromFieldPair (vnp);
        if (from_field == NULL) {
          PointerToDialog (dlg->src_qual_from, NULL);
        } else {
          PointerToDialog (dlg->src_qual_from, from_field->data.ptrvalue);
        }
        from_field = FieldTypeFree (from_field);
        to_field = GetToFieldFromFieldPair (vnp);
        if (to_field == NULL) {
          PointerToDialog (dlg->src_qual_to, NULL);
        } else {
          PointerToDialog (dlg->src_qual_to, to_field->data.ptrvalue);
        }
        to_field = FieldTypeFree (to_field);
        dlg->field_type = FieldType_source_qual;
        break;
      case FieldPairType_feature_field:
        SafeShow (dlg->feature_field_grp);
        ffp = (FeatureFieldPairPtr) vnp->data.ptrvalue;
        if (ffp == NULL) {
          PointerToDialog (dlg->feature_type, NULL);
          PointerToDialog (dlg->feature_field_from, NULL);
          PointerToDialog (dlg->feature_field_to, NULL);
        } else {
          vn.choice = (Uint1)ffp->type;
          vn.data.ptrvalue = NULL;
          vn.next = NULL;
          PointerToDialog (dlg->feature_type, &vn);
          PointerToDialog (dlg->feature_field_from, ffp->field_from);
          PointerToDialog (dlg->feature_field_to, ffp->field_to);
        }
        dlg->field_type = FieldType_feature_field;
        break;
      case FieldPairType_cds_gene_prot:
        SafeShow (dlg->cdsgeneprot_grp);
        cgp = (CDSGeneProtFieldPairPtr) vnp->data.ptrvalue;
        if (cgp == NULL) {
          PointerToDialog (dlg->cdsgeneprot_from, NULL);
          PointerToDialog (dlg->cdsgeneprot_to, NULL);
        } else {
          from_field = GetFromFieldFromFieldPair (vnp);
          PointerToDialog (dlg->cdsgeneprot_from, from_field);
          from_field = FieldTypeFree (from_field);
          to_field = GetToFieldFromFieldPair (vnp);
          PointerToDialog (dlg->cdsgeneprot_to, to_field);
          to_field = FieldTypeFree (to_field);
        }
        dlg->field_type = FieldType_cds_gene_prot;
        break;
      case FieldPairType_molinfo_field:
        SafeShow (dlg->molinfo_dlg);
        PointerToDialog (dlg->molinfo_dlg, vnp->data.ptrvalue);
        dlg->field_type = FieldPairType_molinfo_field;
        break;
      case FieldPairType_rna_field:
        SafeShow (dlg->rna_grp);
        rna_quals = (RnaQualPairPtr) vnp->data.ptrvalue;
        if (rna_quals == NULL) {
          PointerToDialog (dlg->rna_feat_dlg, NULL);
          PointerToDialog (dlg->rna_field_from, NULL);
          PointerToDialog (dlg->rna_field_to, NULL);
        } else {
          PointerToDialog (dlg->rna_feat_dlg, rna_quals->type);
          vn.choice = (Uint1)rna_quals->field_from;
          vn.data.ptrvalue = NULL;
          vn.next = NULL;
          PointerToDialog (dlg->rna_field_from, &vn);
          vn.choice = (Uint1)rna_quals->field_to;
          PointerToDialog (dlg->rna_field_to, &vn);
        }
        dlg->field_type = FieldPairType_rna_field;
        break;
      case FieldPairType_dblink:
        SafeShow (dlg->dblink_grp);
        dblink = (DBLinkFieldPairPtr) vnp->data.ptrvalue;
        if (dblink == NULL) {
          SetValue (dlg->dblink_field_from, 1);
          SetValue (dlg->dblink_field_to, 1);
        } else {
          SetValue (dlg->dblink_field_from, dblink->from);
          SetValue (dlg->dblink_field_to, dblink->to);
        }
        dlg->field_type = FieldPairType_dblink;
        break;
    }
  }
}


static Pointer DialogToFieldPairType (DialoG d)
{
  FieldPairTypeDlgPtr dlg;
  ValNodePtr vnp = NULL, vnp1, vnp2;
  SourceQualPairPtr sqp;
  FeatureFieldPairPtr ffp;
  CDSGeneProtFieldPairPtr cgp;
  RnaQualPairPtr          rna_quals;
  DBLinkFieldPairPtr dblink;

  dlg = (FieldPairTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  switch (dlg->field_type) {
    case FieldPairType_source_qual:
      vnp1 = DialogToPointer (dlg->src_qual_from);
      vnp2 = DialogToPointer (dlg->src_qual_to);
      if (vnp1 != NULL || vnp2 != NULL) {
        sqp = SourceQualPairNew ();
        if (vnp1 != NULL) {
          sqp->field_from = vnp1->data.intvalue;
        }
        if (vnp2 != NULL) {
          sqp->field_to = vnp2->data.intvalue;
        }
        vnp = ValNodeNew (NULL);
        vnp->choice = FieldPairType_source_qual;
        vnp->data.ptrvalue = sqp;
      }
      vnp1 = SourceQualChoiceFree (vnp1);
      vnp2 = SourceQualChoiceFree (vnp2); 
      break;
    case FieldPairType_feature_field:
      vnp2 = (ValNodePtr) DialogToPointer (dlg->feature_type);
      if (vnp2 != NULL) {
        ffp = FeatureFieldPairNew ();
        ffp->type = vnp2->choice;
        vnp2 = ValNodeFree (vnp2);
        ffp->field_from = DialogToPointer (dlg->feature_field_from);
        ffp->field_to = DialogToPointer (dlg->feature_field_to);
        vnp = ValNodeNew (NULL);
        vnp->choice = FieldPairType_feature_field;
        vnp->data.ptrvalue = ffp;
      }
      break;
    case FieldPairType_cds_gene_prot:
      vnp1 = DialogToPointer (dlg->cdsgeneprot_from);
      vnp2 = DialogToPointer (dlg->cdsgeneprot_to);
      if (vnp1 != NULL || vnp2 != NULL) {
        cgp = CDSGeneProtFieldPairNew ();
        if (vnp1 != NULL) {
          cgp->field_from = vnp1->data.intvalue;
        }
        if (vnp2 != NULL) {
          cgp->field_to = vnp2->data.intvalue;
        }
        vnp = ValNodeNew (NULL);
        vnp->choice = FieldPairType_cds_gene_prot;
        vnp->data.ptrvalue = cgp;
      }
      vnp1 = ValNodeFree (vnp1);
      vnp2 = ValNodeFree (vnp2);
      break;
    case FieldPairType_molinfo_field:
      vnp1 = DialogToPointer (dlg->molinfo_dlg);
      if (vnp1 != NULL) {
        vnp = ValNodeNew (NULL);
        vnp->choice = FieldPairType_molinfo_field;
        vnp->data.ptrvalue = vnp1;
      }
      break;
    case FieldPairType_rna_field:
      rna_quals = RnaQualPairNew ();
      rna_quals->type = DialogToPointer (dlg->rna_feat_dlg);
      vnp1 = DialogToPointer (dlg->rna_field_from);
      if (vnp1 != NULL) {
        rna_quals->field_from = vnp1->choice;
        vnp1 = ValNodeFree (vnp1);
      }
      vnp1 = DialogToPointer (dlg->rna_field_to);
      if (vnp1 != NULL) {
        rna_quals->field_to = vnp1->choice;
        vnp1 = ValNodeFree (vnp1);
      }
      vnp = ValNodeNew (NULL);
      vnp->choice = FieldPairType_rna_field;
      vnp->data.ptrvalue = rna_quals;
      break;
    case FieldPairType_dblink:
      dblink = DBLinkFieldPairNew ();
      dblink->from = GetValue (dlg->dblink_field_from);
      dblink->to = GetValue (dlg->dblink_field_to);
      vnp = ValNodeNew (NULL);
      vnp->choice = FieldPairType_dblink;
      vnp->data.ptrvalue = dblink;
      break;
  }

  return vnp;
}


static ValNodePtr TestFieldPairTypeDialog (DialoG d)
{
  FieldPairTypeDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (FieldPairTypeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  switch (dlg->field_type) {
    case FieldPairType_source_qual:
      err_list = TestDialog (dlg->src_qual_from);
      ValNodeLink (&err_list, TestDialog (dlg->src_qual_to));
      break;
    case FieldPairType_feature_field:
      err_list = TestDialog (dlg->feature_type);
      ValNodeLink (&err_list, TestDialog (dlg->feature_field_from));
      ValNodeLink (&err_list, TestDialog (dlg->feature_field_to));
      break;
    case FieldPairType_cds_gene_prot:
      ValNodeLink (&err_list, TestDialog (dlg->cdsgeneprot_from));
      ValNodeLink (&err_list, TestDialog (dlg->cdsgeneprot_to));
      break;
    case FieldPairType_molinfo_field:
      ValNodeLink (&err_list, TestDialog (dlg->molinfo_dlg));
      break;
    case FieldPairType_rna_field:
      ValNodeLink (&err_list, TestDialog (dlg->rna_feat_dlg));
      ValNodeLink (&err_list, TestDialog (dlg->rna_field_from));
      ValNodeLink (&err_list, TestDialog (dlg->rna_field_to));
      break;
    case FieldPairType_dblink:
      break;
    default:
      ValNodeAddPointer (&err_list, 0, "No field type");
      break;
  }
  return err_list;
}


static void ChangeFieldPairTypeDialogPopup (PopuP p)
{
  FieldPairTypeDlgPtr dlg;

  dlg = (FieldPairTypeDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}



static DialoG FieldPairTypeDialog (GrouP h, Boolean for_convert, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  FieldPairTypeDlgPtr dlg;
  GrouP               p, g;
  ValNodePtr          field_list;
  
  dlg = (FieldPairTypeDlgPtr) MemNew (sizeof (FieldPairTypeDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, 0, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FieldPairTypeToDialog;
  dlg->fromdialog = DialogToFieldPairType;
  dlg->testdialog = TestFieldPairTypeDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->src_qual_grp = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (dlg->src_qual_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->src_qual_grp, "To", 0, dialogTextHeight, programFont, 'l');
  dlg->src_qual_from = SourceQualChoiceDialog (dlg->src_qual_grp, TRUE, FALSE, TRUE, change_notify, change_userdata);
  dlg->src_qual_to = SourceQualChoiceDialog (dlg->src_qual_grp, TRUE, FALSE, TRUE, change_notify, change_userdata);

  dlg->feature_field_grp = HiddenGroup (p, 3, 0, NULL);
  StaticPrompt (dlg->feature_field_grp, "Feature Type", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->feature_field_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->feature_field_grp, "To", 0, dialogTextHeight, programFont, 'l');
  dlg->feature_type = FeatureTypeDialogMulti (dlg->feature_field_grp, change_notify, change_userdata);
  dlg->feature_field_from = FeatQualChoiceDialog (dlg->feature_field_grp, FALSE, change_notify, change_userdata);
  dlg->feature_field_to = FeatQualChoiceDialog (dlg->feature_field_grp, FALSE, change_notify, change_userdata);
  
  dlg->cdsgeneprot_grp = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (dlg->cdsgeneprot_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->cdsgeneprot_grp, "To", 0, dialogTextHeight, programFont, 'l');
  dlg->cdsgeneprot_from = CDSGeneProtFieldDialog (dlg->cdsgeneprot_grp, change_notify, change_userdata);
  dlg->cdsgeneprot_to = CDSGeneProtFieldDialog (dlg->cdsgeneprot_grp, change_notify, change_userdata);

  if (for_convert) {
    dlg->molinfo_dlg = SequenceQualPairDialog (p, change_notify, change_userdata);
  }

  dlg->rna_grp = HiddenGroup (p, -1, 0, NULL);
  StaticPrompt (dlg->rna_grp, "RNA Type", 0, dialogTextHeight, programFont, 'l');
  dlg->rna_feat_dlg = RnaFeatureTypeDialog (dlg->rna_grp, "RNA Type of Feature to be Edited", TRUE, change_notify, change_userdata);
  g = HiddenGroup (dlg->rna_grp, 2, 0, NULL);
  StaticPrompt (g, "From RNA Field", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (g, "To RNA Field", 0, dialogTextHeight, programFont, 'l');
  field_list = GetRnaFieldList ();
  dlg->rna_field_from = ValNodeSelectionDialog (g, field_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type",
                                change_notify, change_userdata, FALSE);
  field_list = GetRnaFieldList ();
  dlg->rna_field_to = ValNodeSelectionDialog (g, field_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type",
                                change_notify, change_userdata, FALSE); 
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->rna_feat_dlg, (HANDLE) g, NULL); 

  dlg->dblink_grp = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (dlg->dblink_grp, "From", 0, dialogTextHeight, programFont, 'l');
  StaticPrompt (dlg->dblink_grp, "To", 0, dialogTextHeight, programFont, 'l');
  dlg->dblink_field_from = MakeDbLinkFieldPopup (dlg->dblink_grp, ChangeFieldPairTypeDialogPopup, dlg);
  dlg->dblink_field_to = MakeDbLinkFieldPopup (dlg->dblink_grp, ChangeFieldPairTypeDialogPopup, dlg);

  SafeHide (dlg->src_qual_grp);
  SafeHide (dlg->feature_field_grp);
  SafeHide (dlg->cdsgeneprot_grp);
  SafeHide (dlg->molinfo_dlg);
  SafeHide (dlg->rna_grp);
  SafeHide (dlg->dblink_grp);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->src_qual_grp,
                              (HANDLE) dlg->feature_field_grp,
                              (HANDLE) dlg->cdsgeneprot_grp, 
                              (HANDLE) dlg->molinfo_dlg,
                              (HANDLE) dlg->rna_grp,
                              (HANDLE) dlg->dblink_grp,
                              NULL);

  return (DialoG) p;
}


typedef struct singlefielddlg {
  DIALOG_MESSAGE_BLOCK

  DialoG field_type;
  DialoG field;

  Nlm_ChangeNotifyProc change_notify;
  Pointer change_userdata;
} SingleFieldDlgData, PNTR SingleFieldDlgPtr;


static void SingleFieldTypeToDialog (DialoG d, Pointer data)
{
  SingleFieldDlgPtr dlg;
  ValNodePtr field;
  ValNode vn;

  dlg = (SingleFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  MemSet (&vn, 0, sizeof (ValNode));
  field = (ValNodePtr) data;
  if (field == NULL) {
    vn.choice = FieldType_source_qual;
    PointerToDialog (dlg->field_type, &vn);
    PointerToDialog (dlg->field, &vn);
  } else {
    PointerToDialog (dlg->field_type, field);
    PointerToDialog (dlg->field, field);
  }

}


static Pointer DialogToSingleFieldType (DialoG d)
{
  SingleFieldDlgPtr dlg;

  dlg = (SingleFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  return DialogToPointer (dlg->field);
}


static ValNodePtr TestSingleFieldTypeDialog (DialoG d)
{
  SingleFieldDlgPtr dlg;

  dlg = (SingleFieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  return TestDialog (dlg->field);
}


static void ChangeSingleFieldType (Pointer data)
{
  SingleFieldDlgPtr dlg;
  ValNodePtr vnp;

  if ((dlg = (SingleFieldDlgPtr) data) == NULL) {
    return;
  }
  vnp = DialogToPointer (dlg->field_type);
  vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
  PointerToDialog (dlg->field, vnp);
  vnp = FieldTypeFree (vnp);

}


static DialoG SingleFieldTypeDialog (GrouP h, Boolean text_only, Boolean for_remove, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  SingleFieldDlgPtr dlg;
  GrouP             p;
  ValNodePtr        val_list = NULL;
  
  dlg = (SingleFieldDlgPtr) MemNew (sizeof (SingleFieldDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SingleFieldTypeToDialog;
  dlg->fromdialog = DialogToSingleFieldType;
  dlg->testdialog = TestSingleFieldTypeDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  
  ValNodeAddPointer (&val_list, FieldType_source_qual, StringSave ("Source Qual"));
  ValNodeAddPointer (&val_list, FieldType_feature_field, StringSave ("Feature Qual"));
  ValNodeAddPointer (&val_list, FieldType_cds_gene_prot, StringSave ("CDS-Gene-Prot Qual"));
  ValNodeAddPointer (&val_list, FieldType_rna_field, StringSave ("RNA Qual"));
  ValNodeAddPointer (&val_list, FieldType_molinfo_field, StringSave ("MolInfo Qual"));
  ValNodeAddPointer (&val_list, FieldType_pub, StringSave ("Pub Field"));
  ValNodeAddPointer (&val_list, FieldType_struc_comment_field, StringSave ("Structured Comment Field"));
  ValNodeAddPointer (&val_list, FieldType_misc, StringSave ("Misc"));
  ValNodeAddPointer (&val_list, FieldType_dblink, StringSave ("DBLink"));

  dlg->field_type = ValNodeSelectionDialog (p, val_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type", 
                                ChangeSingleFieldType, dlg, FALSE);
  val_list = NULL;
  dlg->field = FieldTypeDialog (p, FALSE, FALSE, change_notify, change_userdata);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->field_type, dlg->field, NULL);
  return (DialoG) p;
}


#define AECR_DLG_BLOCK       \
  DIALOG_MESSAGE_BLOCK \
  DialoG qual_type_dlg; \
  DialoG field_dlg; \
  ButtoN autopopulate_btn; \
  Nlm_ChangeNotifyProc     autopopulate; \
  Pointer                  autopopulate_data; \
  Nlm_ChangeNotifyProc     change_notify; \
  Pointer                  change_userdata; \
  Nlm_ChangeNotifyProc     redraw_notify; \
  Pointer                  redraw_userdata;

typedef struct aecractiondlg {
  AECR_DLG_BLOCK
} AECRActionDlgData, PNTR AECRActionDlgPtr;


static void SetAECRActionDlgFieldTypeDialogs (DialoG d, DialoG qual_type_dlg, DialoG field_dlg)
{
  AECRActionDlgPtr dlg;

  dlg = (AECRActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    dlg->qual_type_dlg = qual_type_dlg;
    dlg->field_dlg = field_dlg;
  }
}


static Uint1 FieldTypeChoiceFromAECRActionDlg (DialoG d)
{
  AECRActionDlgPtr dlg;
  ValNodePtr       vnp;
  Uint1            rval = FieldType_source_qual;
 
  dlg = (AECRActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    vnp = DialogToPointer (dlg->qual_type_dlg);
    if (vnp != NULL) {
      rval = vnp->choice;
      vnp = ValNodeFree (vnp);
    }
  }
  return rval;
}


static void SingleFieldToAECRActionDlg (DialoG d, ValNodePtr field)
{
  AECRActionDlgPtr dlg;
  ValNode          vn;

  dlg = (AECRActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (field == NULL) {
      PointerToDialog (dlg->qual_type_dlg, NULL);
    } else {
      vn.choice = field->choice;
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->qual_type_dlg, &vn);
    }
    PointerToDialog (dlg->field_dlg, field);
  }
}


static void FieldPairToAECRActionDlg (DialoG d, ValNodePtr fields)
{
  AECRActionDlgPtr dlg;
  ValNode          vn;

  dlg = (AECRActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    if (fields == NULL) {
      PointerToDialog (dlg->qual_type_dlg, NULL);
    } else {
      vn.choice = FieldTypeChoiceFromFieldPairTypeChoice (fields->choice);
      vn.data.ptrvalue = NULL;
      vn.next = NULL;
      PointerToDialog (dlg->qual_type_dlg, &vn);
    }
    PointerToDialog (dlg->field_dlg, fields);
  }
}


static Boolean GetAutopopulateStatus (DialoG d)
{
  AECRActionDlgPtr dlg;

  dlg = (AECRActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->autopopulate_btn == NULL) {
    return FALSE;
  } else {
    return GetStatus (dlg->autopopulate_btn);
  }
}


static void SetAutopopulateStatus (DialoG d, Boolean status)
{
  AECRActionDlgPtr dlg;

  dlg = (AECRActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL || dlg->autopopulate_btn == NULL) {
    return;
  } else {
    SetStatus (dlg->autopopulate_btn, status);
  }
}


static void ChangeAutopopulateStatus (ButtoN b)
{
  AECRActionDlgPtr dlg;

  dlg = (AECRActionDlgPtr) GetObjectExtra (b);
  if (dlg != NULL) {
    SetAppParam ("SEQUINCUSTOM", "BATCHDIALOG", "AUTOPOPULATE", GetStatus (b) ? "TRUE" : "FALSE");
    if (dlg->autopopulate) {
      (dlg->autopopulate)(dlg->autopopulate_data);
    }
  }
}


typedef struct applyactiondlg {
  AECR_DLG_BLOCK

  TexT   value_txt;
  DialoG existing_text;
} ApplyActionDlgData, PNTR ApplyActionDlgPtr;


static Pointer DialogToApplyAction (DialoG d)
{
  ApplyActionDlgPtr dlg;
  ApplyActionPtr     apply;

  dlg = (ApplyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  apply = ApplyActionNew();
  apply->field = DialogToPointer (dlg->field_dlg);
  if (apply->field == NULL) {
    ValNodeAddPointer (&apply->field, FieldTypeChoiceFromAECRActionDlg(d), NULL);
  }

  if (IsFieldTypeNonText (apply->field)) {
    apply->value = StringSave ("TRUE");
  } else {
    apply->value = SaveStringFromText (dlg->value_txt);
  }        
  apply->existing_text = GetExistingTextDialogValue (dlg->existing_text);
  return (Pointer) apply; 
}


static void ApplyActionToDialog (DialoG d, Pointer data)
{
  ApplyActionDlgPtr dlg;
  ApplyActionPtr    apply;

  dlg = (ApplyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  apply = (ApplyActionPtr) data;

  if (apply == NULL) {
    SingleFieldToAECRActionDlg (d, NULL);
    SetTitle (dlg->value_txt, "");
    SetExistingTextDialogValue(dlg->existing_text, 0);
  } else {
    SingleFieldToAECRActionDlg (d, apply->field);
    SetTitle (dlg->value_txt, apply->value);
    SetExistingTextDialogValue(dlg->existing_text, apply->existing_text);
  }
}


static ValNodePtr TestApplyActionDialog (DialoG d)
{
  ApplyActionDlgPtr dlg;
  Int2 field_type = 0;
  ValNodePtr err_list = NULL;
  ValNodePtr field, vnp;
  Boolean    field_is_nontext = FALSE;

  dlg = (ApplyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);
  vnp = DialogToPointer (dlg->qual_type_dlg);
  if (vnp != NULL) {
    field_type = vnp->choice;
    vnp = ValNodeFree (vnp);
  }

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  ValNodeLink (&err_list, TestDialog (dlg->existing_text));

  field = DialogToPointer (dlg->field_dlg);
  field_is_nontext = IsFieldTypeNonText (field);
  field = FieldTypeFree (field);

  if (TextHasNoText (dlg->value_txt) && !field_is_nontext) {
    ValNodeAddPointer (&err_list, 0, "no apply text");
  }
  return err_list;
}


static void ChangeDialogForApplyFieldChoice (DialoG d)
{
  ValNodePtr field;
  ApplyActionDlgPtr dlg;

  dlg = (ApplyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  field = DialogToPointer (dlg->field_dlg);
  if (IsFieldTypeNonText (field)) {
    DisableNonTextOptions (dlg->existing_text);
    Hide (dlg->value_txt);
  } else {
    EnableNonTextOptions (dlg->existing_text);
    Show (dlg->value_txt);
  }
  if (AllowFieldMulti (field)) {
    EnableMultiOptions (dlg->existing_text);
  } else {
    DisableMultiOptions (dlg->existing_text);
  }
  field = FieldTypeFree (field);
}


static void ChangeApplyActionDialogText (TexT t)
{
  ApplyActionDlgPtr dlg;

  dlg = (ApplyActionDlgPtr) GetObjectExtra (t);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}

static void SetApplyActionDialogText (DialoG d, CharPtr str)
{
  ApplyActionDlgPtr dlg;

  dlg = (ApplyActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    SetTitle (dlg->value_txt, str);
  }
}


static DialoG 
ApplyActionDialog 
(GrouP h,
 Boolean indexer_version,
 Boolean show_existing_text,
 Nlm_ChangeNotifyProc     autopopulate,
 Pointer                  autopopulate_data,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  ApplyActionDlgPtr     dlg;
  GrouP                 p;

  dlg = (ApplyActionDlgPtr) MemNew (sizeof (ApplyActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = ApplyActionToDialog;
  dlg->fromdialog = DialogToApplyAction;
  dlg->testdialog = TestApplyActionDialog;
  dlg->autopopulate = autopopulate;
  dlg->autopopulate_data = autopopulate_data;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  if (dlg->autopopulate != NULL) {
    dlg->autopopulate_btn = CheckBox (p, "Autopopulate", ChangeAutopopulateStatus);
    SetObjectExtra (dlg->autopopulate_btn, dlg, NULL);
  }
  dlg->value_txt = DialogText (p, "", 20, ChangeApplyActionDialogText);
  SetObjectExtra (dlg->value_txt, dlg, NULL);
  if (show_existing_text) {
    dlg->existing_text = ExistingTextDialog (p, change_notify, change_userdata);
  } else {
    dlg->existing_text = NULL;
  }
  if (dlg->autopopulate == NULL) {
    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->value_txt,
                                (HANDLE) dlg->existing_text,
                                NULL);
  } else {
    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->autopopulate_btn,
                                (HANDLE) dlg->value_txt,
                                (HANDLE) dlg->existing_text,
                                NULL);
  }
  ChangeDialogForApplyFieldChoice ((DialoG) p);
  return (DialoG) p;
}


typedef struct editactiondlg {
  AECR_DLG_BLOCK

  DialoG edit;
} EditActionDlgData, PNTR EditActionDlgPtr;


static void PointerToEditActionDlg (DialoG d, Pointer data)
{
  EditActionDlgPtr dlg;
  EditActionPtr    edit;

  dlg = (EditActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  edit = (EditActionPtr) data;

  if (edit == NULL) {
    SingleFieldToAECRActionDlg (d, NULL);
    PointerToDialog (dlg->edit, NULL);
  } else {
    SingleFieldToAECRActionDlg (d, edit->field);
    PointerToDialog (dlg->edit, edit->edit);
  }
}


static Pointer DialogToEditAction (DialoG d)
{
  EditActionDlgPtr dlg;
  EditActionPtr      edit;

  dlg = (EditActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  edit = EditActionNew ();
  edit->field = DialogToPointer (dlg->field_dlg);
  edit->edit = DialogToPointer (dlg->edit);
  return edit;
}


static ValNodePtr TestEditActionDialog (DialoG d)
{
  EditActionDlgPtr dlg;
  Int2 field_type = 0;
  ValNodePtr err_list = NULL;
  ValNodePtr field, vnp;
  Boolean    field_is_nontext = FALSE;

  dlg = (EditActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);
  vnp = DialogToPointer (dlg->qual_type_dlg);
  if (vnp != NULL) {
    field_type = vnp->choice;
    vnp = ValNodeFree (vnp);
  }

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));

  field = DialogToPointer (dlg->field_dlg);
  field_is_nontext = IsFieldTypeNonText (field);
  field = FieldTypeFree (field);

  if (field_is_nontext) {
    ValNodeAddPointer (&err_list, 0, "invalid action for field type");
  }
  ValNodeLink (&err_list, TestDialog (dlg->edit));
  return err_list;
}


static void SetEditActionDialogText (DialoG d, CharPtr str_find, CharPtr str_repl)
{
  EditActionDlgPtr dlg;

  dlg = (EditActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    SetFieldEditDialogText (dlg->edit, str_find, str_repl);
  }
}


static void EditActionDialogMessage (DialoG d, Int2 mssg)

{
  EditActionDlgPtr dlg;

  dlg = (EditActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    switch (mssg) {
      case NUM_VIB_MSG + 1 :
        SetFieldEditDialogText (dlg->edit, "", "");
        break;
      default :
        break;
    }
  }
}


static DialoG 
EditActionDialog 
(GrouP h,
 Boolean indexer_version,
 Nlm_ChangeNotifyProc     autopopulate,
 Pointer                  autopopulate_data,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  EditActionDlgPtr dlg;
  GrouP                 p;

  dlg = (EditActionDlgPtr) MemNew (sizeof (EditActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToEditActionDlg;
  dlg->fromdialog = DialogToEditAction;
  dlg->testdialog = TestEditActionDialog;
  dlg->dialogmessage = EditActionDialogMessage;
  dlg->autopopulate = autopopulate;
  dlg->autopopulate_data = autopopulate_data;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  if (dlg->autopopulate != NULL) {
    dlg->autopopulate_btn = CheckBox (p, "Autopopulate", ChangeAutopopulateStatus);
    SetObjectExtra (dlg->autopopulate_btn, dlg, NULL);
  }

  dlg->edit = FieldEditDialog (p, change_notify, change_userdata);
  if (dlg->autopopulate_btn != NULL) {
    AlignObjects (ALIGN_CENTER, (HANDLE) dlg->autopopulate_btn, (HANDLE) dlg->edit, NULL);
  }

  return (DialoG) p;
}


typedef struct convertactiondlg {
  AECR_DLG_BLOCK

  ButtoN strip_name;
  ButtoN keep_original;
  DialoG capitalization;
  DialoG existing_text;

} ConvertActionDlgData, PNTR ConvertActionDlgPtr;


static void PointerToConvertActionDlg (DialoG d, Pointer data)
{
  ConvertActionDlgPtr dlg;
  ConvertActionPtr    convert;

  dlg = (ConvertActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  convert = (ConvertActionPtr) data;

  if (convert == NULL) {
    FieldPairToAECRActionDlg (d, NULL);
    SetStatus (dlg->strip_name, FALSE);
    SetStatus (dlg->keep_original, FALSE);
    SetCapChangeDialogValue (dlg->capitalization, 0);
    SetExistingTextDialogValue(dlg->existing_text, 0);
  } else {    
    FieldPairToAECRActionDlg (d, convert->fields);
    SetStatus (dlg->strip_name, convert->strip_name);
    SetStatus (dlg->keep_original, convert->keep_original);
    SetCapChangeDialogValue (dlg->capitalization, convert->capitalization);
    SetExistingTextDialogValue(dlg->existing_text, convert->existing_text);
  }
}


static Pointer DialogToConvertAction (DialoG d)
{
  ConvertActionDlgPtr dlg;
  ConvertActionPtr   convert;

  dlg = (ConvertActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  convert = ConvertActionNew();
  convert->fields = DialogToPointer (dlg->field_dlg);
  convert->strip_name = GetStatus (dlg->strip_name);
  convert->keep_original = GetStatus (dlg->keep_original);
  convert->capitalization = GetCapChangeDialogValue (dlg->capitalization);
  convert->existing_text = GetExistingTextDialogValue (dlg->existing_text);
  return convert;
}


static ValNodePtr TestConvertActionDialog (DialoG d)
{
  ConvertActionDlgPtr dlg;
  Int2 field_type = 0;
  ValNodePtr err_list = NULL;
  ValNodePtr field_pair, field_from, field_to, vnp;
  Boolean    field_is_nontext = FALSE;

  dlg = (ConvertActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);
  vnp = DialogToPointer (dlg->qual_type_dlg);
  if (vnp != NULL) {
    field_type = vnp->choice;
    vnp = ValNodeFree (vnp);
  }

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  ValNodeLink (&err_list, TestDialog (dlg->capitalization));
  ValNodeLink (&err_list, TestDialog (dlg->existing_text));

  field_pair = DialogToPointer (dlg->field_dlg);
  field_from = GetFromFieldFromFieldPair (field_pair);
  field_to = GetFromFieldFromFieldPair (field_pair);
  field_pair = FieldPairTypeFree (field_pair);
  field_is_nontext = IsFieldTypeNonText (field_from) || IsFieldTypeNonText (field_to);
  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);

  if (field_is_nontext && field_type != 4) {
    ValNodeAddPointer (&err_list, 0, "invalid action for field type");
  }
  
  return err_list;
}


static DialoG 
ConvertActionDialog 
(GrouP h,
 Boolean indexer_version,
 Boolean show_existing_text,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  ConvertActionDlgPtr dlg;
  GrouP                 p, g;

  dlg = (ConvertActionDlgPtr) MemNew (sizeof (ConvertActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToConvertActionDlg;
  dlg->fromdialog = DialogToConvertAction;
  dlg->testdialog = TestConvertActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  g = HiddenGroup (p, 2, 0, NULL);
  dlg->keep_original = CheckBox (g, "Leave on original", NULL);
  dlg->strip_name = CheckBox (g, "Strip name from text", NULL);

  dlg->capitalization = CapChangeDialog(p, change_notify, change_userdata);

  if (show_existing_text) {
    dlg->existing_text = ExistingTextDialog (p, change_notify, change_userdata);
  } else {
    dlg->existing_text = NULL;
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlg->capitalization, (HANDLE) dlg->existing_text, NULL);
  
  return (DialoG) p;
}


typedef struct copyactiondlg {
  AECR_DLG_BLOCK

  DialoG existing_text;

} CopyActionDlgData, PNTR CopyActionDlgPtr;


static void PointerToCopyActionDlg (DialoG d, Pointer data)
{
  CopyActionDlgPtr dlg;
  CopyActionPtr    copy;

  dlg = (CopyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  copy = (CopyActionPtr) data;

  if (copy == NULL) {
    FieldPairToAECRActionDlg (d, NULL);
    SetExistingTextDialogValue(dlg->existing_text, 0);
  } else {    
    FieldPairToAECRActionDlg (d, copy->fields);
    SetExistingTextDialogValue(dlg->existing_text, copy->existing_text);
  }
}


static Pointer DialogToCopyAction (DialoG d)
{
  CopyActionDlgPtr dlg;
  CopyActionPtr    copy;

  dlg = (CopyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  copy = CopyActionNew();
  copy->fields = DialogToPointer (dlg->field_dlg);
  copy->existing_text = GetExistingTextDialogValue (dlg->existing_text);
  return copy;
}


static ValNodePtr TestCopyActionDialog (DialoG d)
{
  CopyActionDlgPtr dlg;
  Int2 field_type = 0;
  ValNodePtr err_list = NULL;
  ValNodePtr field, vnp;
  Boolean    field_is_nontext = FALSE;

  dlg = (CopyActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);
  vnp = DialogToPointer (dlg->qual_type_dlg);
  if (vnp != NULL) {
    field_type = vnp->choice;
    vnp = ValNodeFree (vnp);
  }

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  ValNodeLink (&err_list, TestDialog (dlg->existing_text));

  field = DialogToPointer (dlg->field_dlg);
  field_is_nontext = IsFieldTypeNonText (field);
  field = FieldTypeFree (field);

  if (field_is_nontext && field_type != 4) {
    ValNodeAddPointer (&err_list, 0, "invalid action for field type");
  }
  
  return err_list;
}


static DialoG 
CopyActionDialog 
(GrouP h,
 Boolean indexer_version,
 Boolean show_existing_text,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  CopyActionDlgPtr dlg;
  GrouP                 p;

  dlg = (CopyActionDlgPtr) MemNew (sizeof (CopyActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToCopyActionDlg;
  dlg->fromdialog = DialogToCopyAction;
  dlg->testdialog = TestCopyActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  if (show_existing_text) {
    dlg->existing_text = ExistingTextDialog (p, change_notify, change_userdata);
  } else {
    dlg->existing_text = NULL;
  }
  
  return (DialoG) p;
}


typedef struct removeactiondlg {
  AECR_DLG_BLOCK
} RemoveActionDlgData, PNTR RemoveActionDlgPtr;


static void RemoveActionToDialog (DialoG d, Pointer data)
{
  RemoveActionDlgPtr dlg;
  RemoveActionPtr    remove;

  dlg = (RemoveActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  remove = (RemoveActionPtr) data;

  if (remove == NULL) {
    SingleFieldToAECRActionDlg (d, NULL);
  } else {
    SingleFieldToAECRActionDlg (d, remove->field); 
    PointerToDialog (dlg->field_dlg, remove->field);
  }
}


static Pointer DialogToRemoveAction (DialoG d)
{
  RemoveActionDlgPtr dlg;
  RemoveActionPtr    remove;

  dlg = (RemoveActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  remove = RemoveActionNew();
  remove->field = DialogToPointer (dlg->field_dlg);
  return remove;

}


static ValNodePtr TestRemoveActionDialog (DialoG d)
{
  RemoveActionDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (RemoveActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  
  return err_list;
}


static DialoG 
RemoveActionDialog 
(GrouP h,
 Boolean indexer_version,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  RemoveActionDlgPtr dlg;
  GrouP                 p;

  dlg = (RemoveActionDlgPtr) MemNew (sizeof (RemoveActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = RemoveActionToDialog;
  dlg->fromdialog = DialogToRemoveAction;
  dlg->testdialog = TestRemoveActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  return (DialoG) p;
}


typedef struct swapactiondlg {
  AECR_DLG_BLOCK
} SwapActionDlgData, PNTR SwapActionDlgPtr;


static void SwapActionToDialog (DialoG d, Pointer data)
{
  SwapActionDlgPtr dlg;
  SwapActionPtr    swap;

  dlg = (SwapActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  swap = (SwapActionPtr) data;

  if (swap == NULL) {
    FieldPairToAECRActionDlg (d, NULL);
  } else {    
    FieldPairToAECRActionDlg (d, swap->fields);
  }

}


static Pointer DialogToSwapAction (DialoG d)
{
  SwapActionDlgPtr dlg;
  SwapActionPtr    swap;

  dlg = (SwapActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  swap = SwapActionNew();
  swap->fields = DialogToPointer (dlg->field_dlg);
  return swap;

}


static ValNodePtr TestSwapActionDialog (DialoG d)
{
  SwapActionDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (SwapActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  
  return err_list;
}


static DialoG 
SwapActionDialog 
(GrouP h,
 Boolean indexer_version,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  SwapActionDlgPtr dlg;
  GrouP                 p;

  dlg = (SwapActionDlgPtr) MemNew (sizeof (SwapActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = SwapActionToDialog;
  dlg->fromdialog = DialogToSwapAction;
  dlg->testdialog = TestSwapActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  return (DialoG) p;
}


typedef struct aecrparseactiondlg {
  AECR_DLG_BLOCK

  DialoG text_portion;
  DialoG transform;
  ButtoN remove_from_parsed;
  ButtoN remove_left;
  ButtoN remove_right;
  DialoG existing_text;

} AECRParseActionDlgData, PNTR AECRParseActionDlgPtr;


static void AECRParseActionToDialog (DialoG d, Pointer data)
{
  AECRParseActionDlgPtr dlg;
  AECRParseActionPtr    parse;

  dlg = (AECRParseActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  parse = (AECRParseActionPtr) data;

  if (parse == NULL) {
    FieldPairToAECRActionDlg (d, NULL);
    PointerToDialog (dlg->text_portion, NULL);
    PointerToDialog (dlg->transform, NULL);
    SetStatus (dlg->remove_from_parsed, FALSE);
    SetStatus (dlg->remove_left, FALSE);
    SetStatus (dlg->remove_right, FALSE);
    SetExistingTextDialogValue(dlg->existing_text, 0);
  } else {    
    FieldPairToAECRActionDlg (d, parse->fields);
    PointerToDialog (dlg->text_portion, parse->portion);
    PointerToDialog (dlg->transform, parse->transform);
    SetStatus (dlg->remove_from_parsed, parse->remove_from_parsed);
    SetStatus (dlg->remove_left, parse->remove_left);
    SetStatus (dlg->remove_right, parse->remove_right);
    SetExistingTextDialogValue(dlg->existing_text, parse->existing_text);
  }

}


static Pointer DialogToAECRParseAction (DialoG d)
{
  AECRParseActionDlgPtr dlg;
  AECRParseActionPtr    parse;

  dlg = (AECRParseActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  parse = AECRParseActionNew();
  parse->fields = DialogToPointer (dlg->field_dlg);
  parse->portion = DialogToPointer (dlg->text_portion);
  parse->transform = DialogToPointer (dlg->transform);
  parse->remove_from_parsed = GetStatus (dlg->remove_from_parsed);
  parse->remove_left = GetStatus (dlg->remove_left);
  parse->remove_right = GetStatus (dlg->remove_right);
  parse->existing_text = GetExistingTextDialogValue (dlg->existing_text);
  return parse;
}


static ValNodePtr TestAECRParseActionDialog (DialoG d)
{
  AECRParseActionDlgPtr dlg;
  ValNodePtr err_list = NULL;
  ValNodePtr field_pair, field;
  Boolean    field_is_nontext = FALSE;

  dlg = (AECRParseActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->qual_type_dlg);

  ValNodeLink (&err_list, TestDialog (dlg->field_dlg));
  ValNodeLink (&err_list, TestDialog (dlg->existing_text));

  field_pair = DialogToPointer (dlg->field_dlg);
  field = GetFromFieldFromFieldPair (field_pair);
  field_is_nontext = IsFieldTypeNonText (field);
  field = FieldTypeFree (field);
  if (!field_is_nontext) {
    field = GetToFieldFromFieldPair (field_pair);
    field_is_nontext = IsFieldTypeNonText (field);
    field = FieldTypeFree (field);
  }

  if (field_is_nontext) {
    ValNodeAddPointer (&err_list, 0, "invalid action for field type");
  }
  ValNodeLink (&err_list, TestDialog (dlg->text_portion));
  
  return err_list;
}


static void ChangeAECRParseActionTextPortion (Pointer data)
{
  AECRParseActionDlgPtr dlg;
  TextPortionPtr        tp;

  dlg = (AECRParseActionDlgPtr) data;
  if (dlg == NULL) {
    return;
  }

  if (GetStatus (dlg->remove_from_parsed)) {
    tp = DialogToPointer (dlg->text_portion);
    if (tp == NULL) {
      Hide (dlg->remove_left);
      Hide (dlg->remove_right);
    } else {
      if (!tp->include_left) {
        Show (dlg->remove_left);
      } else {
        Hide (dlg->remove_left);
      }
      if (!tp->include_right) {
        Show (dlg->remove_right);
      } else {
        Hide (dlg->remove_right);
      }
    }
  } else {
    Hide (dlg->remove_left);
    Hide (dlg->remove_right);
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeAECRParseActionRemove (ButtoN b)
{
  AECRParseActionDlgPtr dlg;

  dlg = (AECRParseActionDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  ChangeAECRParseActionTextPortion (dlg);
}


static DialoG 
AECRParseActionDialog 
(GrouP h,
 Boolean indexer_version,
 Boolean show_existing_text,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  AECRParseActionDlgPtr dlg;
  GrouP                 p, g1, g2;

  dlg = (AECRParseActionDlgPtr) MemNew (sizeof (AECRParseActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = AECRParseActionToDialog;
  dlg->fromdialog = DialogToAECRParseAction;
  dlg->testdialog = TestAECRParseActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  g1 = HiddenGroup (p, 2, 0, NULL);
  dlg->text_portion = TextPortionDialog (g1, TRUE, ChangeAECRParseActionTextPortion, dlg);
  g2 = HiddenGroup (g1, 0, 2, NULL);
  dlg->remove_left = CheckBox (g2, "Also remove entered text", NULL);
  dlg->remove_right = CheckBox (g2, "Also remove entered text", NULL);

  AlignObjects (ALIGN_VERTICAL, (HANDLE) GetTextPortionStartChoiceGroup (dlg->text_portion), (HANDLE) dlg->remove_left, NULL);
  AlignObjects (ALIGN_VERTICAL, (HANDLE) GetTextPortionEndChoiceGroup (dlg->text_portion), (HANDLE) dlg->remove_right, NULL);
  Hide (dlg->remove_left);
  Hide (dlg->remove_right);

  dlg->remove_from_parsed = CheckBox (p, "Remove from parsed field", ChangeAECRParseActionRemove);
  SetObjectExtra (dlg->remove_from_parsed, dlg, NULL);

  dlg->transform = TextTransformSetDialog (p, NULL, NULL);

  if (show_existing_text) {
    dlg->existing_text = ExistingTextDialog (p, change_notify, change_userdata);
  } else {
    dlg->existing_text = NULL;
  }

  AlignObjects (ALIGN_CENTER, 
                              (HANDLE) g1,
                              (HANDLE) dlg->transform,
                              (HANDLE) dlg->remove_from_parsed,
                              (HANDLE) dlg->qual_type_dlg,
                              (HANDLE) dlg->field_dlg,
                              (HANDLE) dlg->existing_text,
                              NULL);
  return (DialoG) p;
}


typedef struct editaecractiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG action_dlgs[7];
  DialoG constraint_dlg;
  ButtoN accept_btn;

  /* groups and dialogs for selecting field types */
  DialoG single_qual_type_dlg;
  DialoG single_field;
  DialoG single_field_remove;
  GrouP  single_field_grp;
  DialoG pair_qual_type_dlg;
  DialoG field_pair;
  DialoG field_pair_convert;
  GrouP  field_pair_grp;

  Uint1      action_type;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
  Nlm_ChangeNotifyProc     redraw_notify;
  Pointer                  redraw_userdata;
} EditAECRActionDlgData, PNTR EditAECRActionDlgPtr;



static Uint1 ActionTypeChoiceFromPopupValue (Int4 i)
{
  Uint1 action_type_choice = ActionChoice_apply;

  switch (i) {
    case 2:
      action_type_choice = ActionChoice_apply;
      break;
    case 3:
      action_type_choice = ActionChoice_edit;
      break;
    case 4:
      action_type_choice = ActionChoice_convert;
      break;
    case 5:
      action_type_choice = ActionChoice_copy;
      break;
    case 6:
      action_type_choice = ActionChoice_swap;
      break;
    case 8:
      action_type_choice = ActionChoice_remove;
      break;
    case 7:
      action_type_choice = ActionChoice_parse;
      break;
  }
  return action_type_choice;
}


static Uint1 EditAECRDlgNumFromAECRActionType (Uint1 action_type)
{
  Uint1 dlg_num = 0;

  switch (action_type) {
    case ActionChoice_apply:
      dlg_num = 0;
      break;
    case ActionChoice_edit:
      dlg_num = 1;
      break;
    case ActionChoice_convert:
      dlg_num = 2;
      break;
    case ActionChoice_copy:
      dlg_num = 3;
      break;
    case ActionChoice_swap:
      dlg_num = 4;
      break;
    case ActionChoice_parse:
      dlg_num = 5;
      break;
    case ActionChoice_remove:
      dlg_num = 6;
      break;
  }
  return dlg_num;
}


static Uint1 FieldTypeChoiceFromEditAECRActionDialog (DialoG d)
{
  EditAECRActionDlgPtr dlg;
  Uint1                rval = FieldType_source_qual;

  dlg = (EditAECRActionDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    rval = FieldTypeChoiceFromAECRActionDlg (dlg->action_dlgs[EditAECRDlgNumFromAECRActionType(dlg->action_type)]);
  }

  return rval;
}


static ValNodePtr TestAECRActionDialog (DialoG d)
{
  EditAECRActionDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (EditAECRActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  err_list = TestDialog (dlg->action_dlgs[EditAECRDlgNumFromAECRActionType(dlg->action_type)]);
  
  return err_list;
}


static Pointer AECRActionFromEditDialog (DialoG d)
{
  EditAECRActionDlgPtr dlg;
  AECRActionPtr      action = NULL;

  dlg = (EditAECRActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  action = AECRActionNew ();
  ValNodeAddPointer (&(action->action), dlg->action_type, DialogToPointer (dlg->action_dlgs[EditAECRDlgNumFromAECRActionType(dlg->action_type)]));
  action->constraint = DialogToPointer (dlg->constraint_dlg);

  return (Pointer) action;
}


static void AECRActionToEditDialog (DialoG d, Pointer data)
{
  EditAECRActionDlgPtr dlg;
  AECRActionPtr action;
  Uint1         dlg_num;

  dlg = (EditAECRActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  action = (AECRActionPtr) data;
  
  for (dlg_num = 0; dlg_num < 7; dlg_num++) {
    Hide (dlg->action_dlgs[dlg_num]);
  }
  if (action == NULL || action->action == NULL) {
    dlg->action_type = ActionChoice_apply;
    PointerToDialog (dlg->action_dlgs[0], NULL);
    Show (dlg->action_dlgs[0]);
  } else {
    dlg->action_type = action->action->choice;
    dlg_num = EditAECRDlgNumFromAECRActionType(dlg->action_type);
    PointerToDialog (dlg->action_dlgs[dlg_num], action->action->data.ptrvalue);
    Show (dlg->action_dlgs[dlg_num]);
  } 
  
  switch (dlg->action_type) {
    case ActionChoice_apply:
    case ActionChoice_edit:
      Show (dlg->single_field_grp);
      Hide (dlg->field_pair_grp);
      Show (dlg->single_field);
      Hide (dlg->single_field_remove);
      break;
    case ActionChoice_remove:
      Show (dlg->single_field_grp);
      Hide (dlg->field_pair_grp);
      Hide (dlg->single_field);
      Show (dlg->single_field_remove);
      break;
    case ActionChoice_convert:
      Show (dlg->field_pair_grp);
      Show (dlg->field_pair_convert);
      Hide (dlg->field_pair);
      Hide (dlg->single_field_grp);
      break;
    case ActionChoice_copy:
    case ActionChoice_swap:
    case ActionChoice_parse:
      Show (dlg->field_pair_grp);
      Hide (dlg->field_pair_convert);
      Show (dlg->field_pair);
      Hide (dlg->single_field_grp);
      break;
  }
  /* set constraint */
  if (action == NULL) {
    PointerToDialog (dlg->constraint_dlg, NULL);
  } else {
    PointerToDialog (dlg->constraint_dlg, action->constraint);
  }        
}


static void ChangeEditAECRActionQualType (Pointer data)
{
  EditAECRActionDlgPtr dlg;
  ValNodePtr           vnp;

  dlg = (EditAECRActionDlgPtr) data;
  if (dlg != NULL) {
    switch (dlg->action_type) {
      case ActionChoice_apply:
      case ActionChoice_edit:
        vnp = DialogToPointer (dlg->single_qual_type_dlg);
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        PointerToDialog (dlg->single_field, vnp);
        vnp = FieldTypeFree (vnp);
        break;
      case ActionChoice_remove:
        vnp = DialogToPointer (dlg->single_qual_type_dlg);
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        PointerToDialog (dlg->single_field_remove, vnp);
        vnp = FieldTypeFree (vnp);
        break;
      case ActionChoice_convert:
        vnp = DialogToPointer (dlg->pair_qual_type_dlg);
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        PointerToDialog (dlg->field_pair_convert, vnp);
        vnp = FieldTypeFree (vnp);
        break;
      case ActionChoice_copy:
      case ActionChoice_swap:
      case ActionChoice_parse:
        vnp = DialogToPointer (dlg->pair_qual_type_dlg);
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        PointerToDialog (dlg->field_pair, vnp);
        vnp = FieldTypeFree (vnp);
        break;
    }
  
    if (dlg->action_type == ActionChoice_apply) {
      ChangeDialogForApplyFieldChoice (dlg->action_dlgs[0]);
    }
    if (dlg->change_notify) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static void ChangeEditAECRActionField (Pointer data)
{
  EditAECRActionDlgPtr dlg;

  dlg = (EditAECRActionDlgPtr) data;
  if (dlg != NULL) {
    if (dlg->action_type == ActionChoice_apply) {
      ChangeDialogForApplyFieldChoice (dlg->action_dlgs[0]);
    }
    if (dlg->change_notify) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static DialoG 
EditAECRActionDialog 
(GrouP h,
 Boolean indexer_version,
 Boolean show_existing_text,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  EditAECRActionDlgPtr dlg;
  GrouP                 p, g, field_grp, k;
  ValNodePtr            val_list = NULL;

  dlg = (EditAECRActionDlgPtr) MemNew (sizeof (EditAECRActionDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  dlg->dialog = (DialoG) p;
  dlg->todialog = AECRActionToEditDialog;
  dlg->fromdialog = AECRActionFromEditDialog;
  dlg->testdialog = TestAECRActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->redraw_notify = redraw_notify;
  dlg->redraw_userdata = redraw_userdata;

  dlg->action_type = ActionChoice_apply;

  field_grp = HiddenGroup (p, 0, 0, NULL);
  dlg->single_field_grp = HiddenGroup (field_grp, -1, 0, NULL);
 
  ValNodeAddPointer (&val_list, FieldType_source_qual, StringSave ("Source Qual"));
  ValNodeAddPointer (&val_list, FieldType_feature_field, StringSave ("Feature Qual"));
  ValNodeAddPointer (&val_list, FieldType_cds_gene_prot, StringSave ("CDS-Gene-Prot Qual"));
  ValNodeAddPointer (&val_list, FieldType_rna_field, StringSave ("RNA Qual"));
  ValNodeAddPointer (&val_list, FieldType_molinfo_field, StringSave ("MolInfo Qual"));
  ValNodeAddPointer (&val_list, FieldType_pub, StringSave ("Pub Field"));
  ValNodeAddPointer (&val_list, FieldType_struc_comment_field, StringSave ("Structured Comment Field"));
  ValNodeAddPointer (&val_list, FieldType_misc, StringSave ("Misc"));
  ValNodeAddPointer (&val_list, FieldType_dblink, StringSave ("DBLink"));

  dlg->single_qual_type_dlg = ValNodeSelectionDialog (dlg->single_field_grp, val_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type", 
                                ChangeEditAECRActionQualType, dlg, FALSE);
  val_list = NULL;
  k = HiddenGroup (dlg->single_field_grp, 0, 0, NULL);
  dlg->single_field = FieldTypeDialog (k, FALSE, FALSE, ChangeEditAECRActionField, dlg);
  dlg->single_field_remove = FieldTypeDialog (k, FALSE, TRUE, ChangeEditAECRActionField, dlg);
  AlignObjects (ALIGN_CENTER, (HANDLE) (HANDLE) dlg->single_field, (HANDLE) dlg->single_field_remove, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->single_qual_type_dlg, k, NULL);

  dlg->field_pair_grp = HiddenGroup (field_grp, -1, 0, NULL);
  ValNodeAddPointer (&val_list, FieldType_source_qual, StringSave ("Source Qual"));
  ValNodeAddPointer (&val_list, FieldType_feature_field, StringSave ("Feature Qual"));
  ValNodeAddPointer (&val_list, FieldType_cds_gene_prot, StringSave ("CDS-Gene-Prot Qual"));
  ValNodeAddPointer (&val_list, FieldType_rna_field, StringSave ("RNA Qual"));
  ValNodeAddPointer (&val_list, FieldType_struc_comment_field, StringSave ("Structured Comment Field"));
  ValNodeAddPointer (&val_list, FieldType_molinfo_field, StringSave ("MolInfo Qual"));
  ValNodeAddPointer (&val_list, FieldType_dblink, StringSave ("DBLink"));
  dlg->pair_qual_type_dlg = ValNodeSelectionDialog (dlg->field_pair_grp, val_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type", 
                                ChangeEditAECRActionQualType, dlg, FALSE);
  val_list = NULL;
  
  k = HiddenGroup (dlg->field_pair_grp, 0, 0, NULL);
  dlg->field_pair = FieldPairTypeDialog (k, FALSE, ChangeEditAECRActionField, dlg);
  dlg->field_pair_convert = FieldPairTypeDialog (k, TRUE, ChangeEditAECRActionField, dlg);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->field_pair_convert, (HANDLE) dlg->field_pair, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->pair_qual_type_dlg, (HANDLE) k, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->single_field_grp, (HANDLE) dlg->field_pair_grp, NULL);

  Hide (dlg->single_field_grp);
  Hide (dlg->field_pair_grp);

  g = HiddenGroup (p, 0, 0, NULL);
  dlg->action_dlgs[0] = ApplyActionDialog (g, indexer_version, show_existing_text, NULL, NULL, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[0], dlg->single_qual_type_dlg, dlg->single_field);
  dlg->action_dlgs[1] = EditActionDialog (g, indexer_version, NULL, NULL, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[1], dlg->single_qual_type_dlg, dlg->single_field);
  dlg->action_dlgs[2] = ConvertActionDialog (g, indexer_version, show_existing_text, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[2], dlg->pair_qual_type_dlg, dlg->field_pair_convert);
  dlg->action_dlgs[3] = CopyActionDialog (g, indexer_version, show_existing_text, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[3], dlg->pair_qual_type_dlg, dlg->field_pair);
  dlg->action_dlgs[4] = SwapActionDialog (g, indexer_version, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[4], dlg->pair_qual_type_dlg, dlg->field_pair);
  dlg->action_dlgs[5] = AECRParseActionDialog (g, indexer_version, show_existing_text, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[5], dlg->pair_qual_type_dlg, dlg->field_pair);
  dlg->action_dlgs[6] = RemoveActionDialog (g, indexer_version, change_notify, change_userdata, redraw_notify, redraw_userdata);
  SetAECRActionDlgFieldTypeDialogs (dlg->action_dlgs[6], dlg->single_qual_type_dlg, dlg->single_field_remove);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->action_dlgs[0],
                              (HANDLE) dlg->action_dlgs[1],
                              (HANDLE) dlg->action_dlgs[2],
                              (HANDLE) dlg->action_dlgs[3],
                              (HANDLE) dlg->action_dlgs[4],
                              (HANDLE) dlg->action_dlgs[5],
                              (HANDLE) dlg->action_dlgs[6],
                              NULL);
  dlg->constraint_dlg = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) g,
                              (HANDLE) dlg->constraint_dlg,
                              NULL);
  return (DialoG) p;
}


static SeqFeatPtr CreateFeatureWithImportFeatureQuals (Uint1 featdef_type, ValNodePtr fields)
{
  SeqFeatPtr sfp;
  ValNodePtr vnp;
  FeatQualLegalValPtr q;
  GBQualPtr gbq;
  CharPtr   qualname, cp;

  sfp = SeqFeatNew();
  sfp->idx.subtype = featdef_type;
  sfp->data.choice = FindFeatFromFeatDefType (sfp->idx.subtype);

  for (vnp = fields; vnp != NULL; vnp = vnp->next) {
    q = (FeatQualLegalValPtr) vnp->data.ptrvalue;
    if (q != NULL && q->qual != Feat_qual_legal_gene && q->qual != Feat_qual_legal_gene_description && q->qual != Feat_qual_legal_note) {
      gbq = GBQualNew ();
      qualname = StringSave (GetFeatQualName(q->qual));
      cp = StringChr (qualname, '-');
      while (cp != NULL) {
        *cp = '_';
        cp = StringChr (cp + 1, '-');
      }
      gbq->qual = qualname;
      gbq->val = StringSave (q->val);
      gbq->next = sfp->qual;
      sfp->qual = gbq;
    }
  }
  return sfp;  
}


static CharPtr FindFeatQualValue (ValNodePtr vnp, Int4 featqual)
{
  FeatQualLegalValPtr q;
  CharPtr str = NULL;
  while (vnp != NULL && str == NULL) {
    q = (FeatQualLegalValPtr) vnp->data.ptrvalue;
    if (q != NULL && q->qual == featqual) {
      str = q->val;
    }
    vnp = vnp->next;
  }
  if (str == NULL) {
    str = "";
  }
  return str;
}


static void SetFeatQualValue (ValNodePtr PNTR fields, Int4 featqual, CharPtr val)
{
  FeatQualLegalValPtr q;

  q = FeatQualLegalValNew ();
  q->qual = featqual;
  q->val = StringSave (val);
  ValNodeAddPointer (fields, 1, q);
}


static void AddGBQualsToActionFields (GBQualPtr gbq, ValNodePtr PNTR fields)
{
  Int4 qualtype;

  if (fields == NULL) return;

  while (gbq != NULL) {
    qualtype = GetFeatQualByName (gbq->qual);
    if (qualtype > -1) {
      SetFeatQualValue (fields, qualtype, gbq->val);
    }
    gbq = gbq->next;
  }
}


NLM_EXTERN ApplyFeatureDetailsPtr ApplyFeatureDetailsNew (ApplyFeatureActionPtr action)
{
  ApplyFeatureDetailsPtr details;

  details = (ApplyFeatureDetailsPtr) MemNew (sizeof (ApplyFeatureDetailsData));

  if (action != NULL) {
    details->add_mrna = action->add_mrna;
    details->fields = AsnIoMemCopy (action->fields, (AsnReadFunc) FeatQualLegalSetAsnRead, (AsnWriteFunc) FeatQualLegalSetAsnWrite);
    details->src_fields = AsnIoMemCopy(action->src_fields, (AsnReadFunc) SourceQualValSetAsnRead, (AsnWriteFunc) SourceQualValSetAsnWrite);
  }
  return details;
}


NLM_EXTERN ApplyFeatureDetailsPtr ApplyFeatureDetailsFree (ApplyFeatureDetailsPtr details) 
{
  if (details != NULL) {
    details->fields = FeatQualLegalSetFree (details->fields);
    details->src_fields = SourceQualValSetFree (details->src_fields);
    details = MemFree (details);
  }
  return details;
}


typedef struct applyfeaturedetailsdlg {
  DIALOG_MESSAGE_BLOCK

  /* for CDS */
  ButtoN                   add_mrna;
  PopuP                    reading_frame;
  TexT                     protein_name;
  TexT                     protein_description;

  /* for RNA */
  DialoG                   ncrna_class;
  TexT                     rna_name;

  /* for import features */
  DialoG                   quals;
  
  /* for source features */
  DialoG                   src_quals;

  /* for all */
  TexT                     gene_locus;
  TexT                     gene_description;
  TexT                     comment;

  Uint2                    featdef_type;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} ApplyFeatureDetailsDlgData, PNTR ApplyFeatureDetailsDlgPtr;


static void ApplyFeatureDetailsToDialog (DialoG d, Pointer data)
{
  ApplyFeatureDetailsDlgPtr dlg;
  ApplyFeatureDetailsPtr    details;
  CharPtr                   txt;
  BioSourcePtr              biop;

  dlg = (ApplyFeatureDetailsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  details = (ApplyFeatureDetailsPtr) data;

  if (data == NULL) {
    /* for CDS */
    SafeSetStatus (dlg->add_mrna, FALSE);
    SafeSetValue (dlg->reading_frame, 4);
    SafeSetTitle (dlg->protein_name, "");
    SafeSetTitle (dlg->protein_description, "");

    /* for RNA */
    PointerToDialog (dlg->ncrna_class, NULL);
    SafeSetTitle (dlg->rna_name, "");

    /* for import features */
    PointerToDialog (dlg->quals, NULL);

    /* for source features */
    PointerToDialog (dlg->src_quals, NULL);
  
    /* for all */
    SafeSetTitle (dlg->gene_locus, "");
    SafeSetTitle (dlg->gene_description, "");
    SafeSetTitle (dlg->comment, "");
  } else {
    /* for CDS */
    SafeSetStatus (dlg->add_mrna, details->add_mrna);
    txt = FindFeatQualValue (details->fields, Feat_qual_legal_codon_start);
    if (StringICmp (txt, "best") == 0) {
      SafeSetValue (dlg->reading_frame, 4);
    } else {
      SafeSetValue (dlg->reading_frame, atoi (txt));
    }

    SafeSetTitle (dlg->protein_name, FindFeatQualValue (details->fields, Feat_qual_legal_product));
    SafeSetTitle (dlg->protein_description, FindFeatQualValue (details->fields, Feat_qual_legal_description));

    /* for source feature */
    biop = BioSourceFromSourceQualVals (details->src_fields);
    PointerToDialog (dlg->src_quals, biop);
    biop = BioSourceFree (biop);

    /* for RNA */
    PointerToDialog (dlg->ncrna_class, FindFeatQualValue (details->fields, Feat_qual_legal_ncRNA_class));
    SafeSetTitle (dlg->rna_name, FindFeatQualValue (details->fields, Feat_qual_legal_product));
  
    /* for all */
    SafeSetTitle (dlg->gene_locus, FindFeatQualValue (details->fields, Feat_qual_legal_gene));
    SafeSetTitle (dlg->gene_description, FindFeatQualValue (details->fields, Feat_qual_legal_gene_description));
    SafeSetTitle (dlg->comment, FindFeatQualValue (details->fields, Feat_qual_legal_note));
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static Pointer DialogToApplyFeatureDetails (DialoG d)
{
  ApplyFeatureDetailsDlgPtr dlg;
  ApplyFeatureDetailsPtr    details;
  CharPtr                  txt;
  Char                     num[15];
  Int4                     frame;
  GBQualPtr                gbq;
  BioSourcePtr             biop;

  dlg = (ApplyFeatureDetailsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  details = ApplyFeatureDetailsNew (NULL);

  /* for CDS */
  if (dlg->add_mrna != NULL && GetStatus (dlg->add_mrna)) {
    details->add_mrna = TRUE;
  } else {
    details->add_mrna = FALSE;
  }

  if (dlg->reading_frame != NULL) {
    frame = GetValue (dlg->reading_frame);
    if (frame > 0 && frame < 4) {
      sprintf (num, "%d", frame);
      SetFeatQualValue (&(details->fields), Feat_qual_legal_codon_start, num);
    } else {
      SetFeatQualValue (&(details->fields), Feat_qual_legal_codon_start, "best");
    }
  }

  if (dlg->protein_name != NULL && !TextHasNoText (dlg->protein_name)) {
    txt = SaveStringFromText (dlg->protein_name);
    SetFeatQualValue (&(details->fields), Feat_qual_legal_product, txt);
    txt = MemFree (txt);
  }
  if (dlg->protein_name != NULL && !TextHasNoText (dlg->protein_description)) {
    txt = SaveStringFromText (dlg->protein_description);
    SetFeatQualValue (&(details->fields), Feat_qual_legal_description, txt);
    txt = MemFree (txt);
  }

  /* for RNA */
  if (dlg->ncrna_class != NULL) {
    txt = DialogToPointer (dlg->ncrna_class);
    if (!StringHasNoText (txt)) {
      SetFeatQualValue (&(details->fields), Feat_qual_legal_ncRNA_class, txt);
    }
    txt = MemFree (txt);
  }
  if (dlg->rna_name != NULL && !TextHasNoText (dlg->rna_name)) {
    txt = SaveStringFromText (dlg->rna_name);
    SetFeatQualValue (&(details->fields), Feat_qual_legal_product, txt);
    txt = MemFree (txt);
  }

  /* for source features */
  if (dlg->src_quals != NULL) {
    biop = DialogToPointer (dlg->src_quals);
    details->src_fields = SourceQualValsFromBioSourcePtr (biop);
    biop = BioSourceFree (biop);
  }

  /* for import features */
  if (dlg->quals != NULL) {
    gbq = DialogToPointer (dlg->quals);
    AddGBQualsToActionFields (gbq, &(details->fields));
    gbq = GBQualFree (gbq);
  }
  
  /* for all */
  if (dlg->gene_locus != NULL && !TextHasNoText (dlg->gene_locus)) {
    txt = SaveStringFromText (dlg->gene_locus);
    SetFeatQualValue (&(details->fields), Feat_qual_legal_gene, txt);
    txt = MemFree (txt);
  }
  if (dlg->gene_description != NULL && !TextHasNoText (dlg->gene_description)) {
    txt = SaveStringFromText (dlg->gene_description);
    SetFeatQualValue (&(details->fields), Feat_qual_legal_gene_description, txt);
    txt = MemFree (txt);
  }
  if (dlg->comment != NULL && !TextHasNoText (dlg->comment)) {
    txt = SaveStringFromText (dlg->comment);
    SetFeatQualValue (&(details->fields), Feat_qual_legal_note, txt);
    txt = MemFree (txt);
  }

  return (Pointer) details;
}

static void AddTextToComment (ButtoN b, CharPtr text)
{
  ApplyFeatureDetailsDlgPtr dlg;
  CharPtr                  orig_comment;
  CharPtr                  new_comment;
  
  dlg = (ApplyFeatureDetailsDlgPtr) GetObjectExtra (b);
  if (dlg == NULL || StringHasNoText (text))
  {
    return;
  }
  
  orig_comment = SaveStringFromText (dlg->comment);
  if (StringHasNoText (orig_comment))
  {
    SetTitle (dlg->comment, text);
  }
  else
  {
    new_comment = (CharPtr) MemNew ((StringLen (orig_comment) + StringLen (text) + 3) * sizeof (Char));
    if (new_comment != NULL)
    {
      StringCpy (new_comment, orig_comment);
      StringCat (new_comment, "; ");
      StringCat (new_comment, text);
      SetTitle (dlg->comment, new_comment);
      new_comment = MemFree (new_comment);
    }
  } 
  orig_comment = MemFree (orig_comment);
}


static void Add18SITS28SToComment (ButtoN b) 
{
  AddTextToComment (b, "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA");
}


static void Add16SIGS23SToComment (ButtoN b)
{
  AddTextToComment (b, "contains 16S ribosomal RNA, 16S-23S ribosomal RNA intergenic spacer, and 23S ribosomal RNA");
}


static void ChangeApplyFeatureDetailsBtn (ButtoN b)
{
  ApplyFeatureDetailsDlgPtr dlg;

  dlg = (ApplyFeatureDetailsDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


NLM_EXTERN DialoG ApplyFeatureDetailsDialog (GrouP h, Uint1 featdef_type, ApplyFeatureDetailsPtr details, Boolean indexer_version, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ApplyFeatureDetailsDlgPtr dlg;
  GrouP                 p;
  GrouP                 r, frame_grp = NULL, comment_btns_grp = NULL;
  ButtoN                comment_btn;
  Int4                  seqfeattype;
  SeqFeatPtr            sfp;

  dlg = (ApplyFeatureDetailsDlgPtr) MemNew (sizeof (ApplyFeatureDetailsDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ApplyFeatureDetailsToDialog;
  dlg->fromdialog = DialogToApplyFeatureDetails;
  dlg->dialogmessage = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->featdef_type = featdef_type;

  seqfeattype = FindFeatFromFeatDefType (featdef_type);

  if (featdef_type == FEATDEF_CDS) {
    if (indexer_version) {
      dlg->add_mrna = CheckBox (p, "Also add mRNA", ChangeApplyFeatureDetailsBtn);
      SetObjectExtra (dlg->add_mrna, dlg, NULL);
    }

    frame_grp = HiddenGroup (p, 2, 0, NULL);
    StaticPrompt (frame_grp, "Reading Frame", 0, dialogTextHeight, programFont, 'l');
    dlg->reading_frame = PopupList (frame_grp, TRUE, NULL);
    PopupItem (dlg->reading_frame, "1");
    PopupItem (dlg->reading_frame, "2");
    PopupItem (dlg->reading_frame, "3");
    PopupItem (dlg->reading_frame, "Best");
    SetValue (dlg->reading_frame, 4);
  } else if (featdef_type == FEATDEF_source) {
    dlg->src_quals = BioSourceDialog (p);
  } 

  r = HiddenGroup (p, 2, 0, NULL);
  if (featdef_type == FEATDEF_CDS) {

    StaticPrompt (r, "Protein Name", 0, dialogTextHeight, programFont, 'l');
    dlg->protein_name = DialogText (r, "", 20, NULL);
    StaticPrompt (r, "Protein Description", 0, dialogTextHeight, programFont, 'l');
    dlg->protein_description = DialogText (r, "", 20, NULL);
  } else if (seqfeattype == SEQFEAT_RNA) {
    /* RNA name */
    StaticPrompt (r, "RNA Name", 0, dialogTextHeight, programFont, 'l');
    dlg->rna_name = DialogText (r, "", 20, NULL);
    /* ncRNA class */
    if (featdef_type == FEATDEF_ncRNA) {
      StaticPrompt (r, "ncRNA Class", 0, dialogTextHeight, programFont, 'l');
      dlg->ncrna_class = CreatencRNAClassDialog (r, FALSE, NULL, NULL);
    }
  } else if (featdef_type != FEATDEF_GENE && featdef_type != FEATDEF_source) {
    StaticPrompt (r, "Qualifiers", 0, dialogTextHeight, programFont, 'l');
    sfp = CreateFeatureWithImportFeatureQuals (featdef_type, details == NULL ? NULL : details->fields);
    dlg->quals = NewCreateImportFields (r, GetFeatureNameFromFeatureType (GetFeatureTypeFromFeatdef(featdef_type)), sfp, FALSE);
    sfp = SeqFeatFree (sfp);
  }

  /* gene qualifiers ( for all ) */
  StaticPrompt (r, "Gene Locus", 0, dialogTextHeight, programFont, 'l');
  dlg->gene_locus = DialogText (r, "", 20, NULL);
  StaticPrompt (r, "Gene Description", 0, dialogTextHeight, programFont, 'l');
  dlg->gene_description = DialogText (r, "", 20, NULL);
  StaticPrompt (r, "Comment", 0, dialogTextHeight, programFont, 'l');
  dlg->comment = DialogText (r, "", 20, NULL);  
  if ((featdef_type == FEATDEF_otherRNA || featdef_type == FEATDEF_misc_RNA) && indexer_version) {
    comment_btns_grp = HiddenGroup (p, 2, 0, NULL);
    comment_btn = PushButton (comment_btns_grp, "Add '18S-ITS-5.8S-ITS-28S' to comment", Add18SITS28SToComment);
    SetObjectExtra (comment_btn, dlg, NULL);
    comment_btn = PushButton (comment_btns_grp, "Add '16S-IGS-23S' to comment", Add16SIGS23SToComment);
    SetObjectExtra (comment_btn, dlg, NULL);
  }

  if (frame_grp != NULL) {
    AlignObjects (ALIGN_CENTER, (HANDLE) r,
                                (HANDLE) frame_grp,
                                (HANDLE) dlg->add_mrna,
                                NULL);
  } else if (comment_btns_grp != NULL) {
    AlignObjects (ALIGN_CENTER, (HANDLE) r,
                                (HANDLE) comment_btns_grp,
                                NULL);
  } else if (dlg->src_quals != NULL) {
    AlignObjects (ALIGN_CENTER, (HANDLE) r,
                                (HANDLE) dlg->src_quals,
                                NULL);
  } else {
    AlignObjects (ALIGN_CENTER, (HANDLE) r,
                                NULL);
  }
  return (DialoG) p;
}


typedef struct applyfeatureactiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG                   feature_type_dlg;
  ButtoN                   apply_to_parts;
  TexT                     only_this_part;
  ButtoN                   partial5;
  ButtoN                   partial3;
  GrouP                    strand_group;
  GrouP                    use_whole_interval;
  TexT                     left_end;
  TexT                     right_end;
  GrouP                    all_or_some_group;
  TexT                     accession_list_txt;
  ButtoN                   add_redundant;

  DialoG                   details;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;


} ApplyFeatureActionDlgData, PNTR ApplyFeatureActionDlgPtr;


static ValNodePtr TestApplyFeatureActionDialog (DialoG d)
{
  ApplyFeatureActionDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (ApplyFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  err_list = TestDialog (dlg->feature_type_dlg);

  return err_list;
}


static void ChangeApplyFeatureActionDlg (ApplyFeatureActionDlgPtr dlg)
{
  if (dlg == NULL) return;

  if (GetStatus (dlg->apply_to_parts)) {
    SafeEnable (dlg->only_this_part);
  } else {
    SafeDisable (dlg->only_this_part);
  }

  if (dlg->use_whole_interval != NULL) {
    if (GetValue (dlg->use_whole_interval) == 2) {
      Enable (dlg->left_end);
      Enable (dlg->right_end);
    } else {
      Disable (dlg->left_end);
      Disable (dlg->right_end);
    }
  }

  if (dlg->all_or_some_group != NULL) {
    if (GetValue (dlg->all_or_some_group) == 2) {
      Enable (dlg->accession_list_txt);
    } else {
      Disable (dlg->accession_list_txt);
    }
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static CharPtr SequenceListWrite (ValNodePtr seq_list)
{
  Int4 len = 0;
  ValNodePtr vnp;
  CharPtr    txt = NULL;

  if (seq_list == NULL) return NULL;
  for (vnp = seq_list; vnp != NULL; vnp = vnp->next) {
    len += StringLen (vnp->data.ptrvalue) + 2;
  }

  txt = (CharPtr) MemNew (sizeof (Char) * len);
  for (vnp = seq_list; vnp != NULL; vnp = vnp->next) {
    StringCat (txt, vnp->data.ptrvalue);
    if (vnp->next != NULL) {
      StringCat (txt, ";");
    }
  }
  return txt;
}


static ValNodePtr SequenceListCollect (CharPtr txt)
{
  ValNodePtr seq_list = NULL;
  CharPtr    cp;

  if (StringHasNoText (txt)) {
    return NULL;
  }
  cp = StringChr (txt, ';');
  while (cp != NULL && *txt != 0) {
    *cp = 0;
    if (!StringHasNoText (txt)) {
      ValNodeAddPointer (&seq_list, 0, StringSave (txt));
    }
    *cp = ';';
    txt = cp + 1;
    cp = StringChr (txt, ':');
  }
  if (!StringHasNoText (txt)) {
    ValNodeAddPointer (&seq_list, 0, StringSave (txt));
  }
  return seq_list;
}


static void ApplyFeatureActionToDialog (DialoG d, Pointer data)
{
  ApplyFeatureActionDlgPtr dlg;
  ApplyFeatureActionPtr    action;
  Char                     num[20];
  LocationIntervalPtr      lint;
  CharPtr                  txt;
  ValNode                  vn;
  ApplyFeatureDetailsPtr   details;

  dlg = (ApplyFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  action = (ApplyFeatureActionPtr) data;

  if (action == NULL) {
    PointerToDialog (dlg->feature_type_dlg, NULL);
    SetStatus (dlg->apply_to_parts, FALSE);
    SetTitle (dlg->only_this_part, "");
    SetStatus (dlg->partial5, FALSE);
    SetStatus (dlg->partial3, FALSE);
    SetValue (dlg->strand_group, 1);
    SetValue (dlg->use_whole_interval, 1);
    SetTitle (dlg->left_end, "");
    SetTitle (dlg->right_end, "");
    SetValue (dlg->all_or_some_group, 1);
    SetTitle (dlg->accession_list_txt, "");
    SetStatus (dlg->add_redundant, FALSE);

  } else {
    vn.choice = (Uint1) action->type;
    vn.next = NULL;
    vn.data.ptrvalue = NULL;
    PointerToDialog (dlg->feature_type_dlg, &vn);

    SetStatus (dlg->apply_to_parts, action->apply_to_parts);
    if (action->only_seg_num > -1) {
      sprintf (num, "%d", action->only_seg_num);  
      SetTitle (dlg->only_this_part, num);
    } else {
      SetTitle (dlg->only_this_part, "");
    }
    SetStatus (dlg->partial5, action->partial5);
    SetStatus (dlg->partial3, action->partial3);
    if (action->plus_strand) {
      SetValue (dlg->strand_group, 1);
    } else {
      SetValue (dlg->strand_group, 2);
    }
    if (action->location == NULL || action->location->choice == LocationChoice_whole_sequence) {    
      SetValue (dlg->use_whole_interval, 1);
      SetTitle (dlg->left_end, "");
      SetTitle (dlg->right_end, "");
    } else {
      SetValue (dlg->use_whole_interval, 2);
      if (action->location->choice == LocationChoice_interval) {
        lint = (LocationIntervalPtr) action->location->data.ptrvalue;
        if (lint == NULL) {
          SetTitle (dlg->left_end, "");
          SetTitle (dlg->right_end, "");
        } else {
          sprintf (num, "%d", lint->from);  
          SetTitle (dlg->left_end, num);
          sprintf (num, "%d", lint->to);  
          SetTitle (dlg->right_end, num);
        }
      } else if (action->location->choice == LocationChoice_point) {
        sprintf (num, "%d^", action->location->data.intvalue);
        SetTitle (dlg->left_end, num);
        sprintf (num, "%d", action->location->data.intvalue + 1);  
        SetTitle (dlg->right_end, num);
      } else {
        SetTitle (dlg->left_end, "");
        SetTitle (dlg->right_end, "");
      }
    }

    if (action->seq_list == NULL || action->seq_list->choice == SequenceListChoice_all) {
      SetValue (dlg->all_or_some_group, 1);
      SetTitle (dlg->accession_list_txt, "");
    } else {
      txt = SequenceListWrite (action->seq_list->data.ptrvalue);
      SetValue (dlg->all_or_some_group, 2);
      SetTitle (dlg->accession_list_txt, txt);
      txt = MemFree (txt);
    }

    SetStatus (dlg->add_redundant, action->add_redundant);
  }

  details = ApplyFeatureDetailsNew (action);
  PointerToDialog (dlg->details, details);
  details = ApplyFeatureDetailsFree (details);
  ChangeApplyFeatureActionDlg (dlg);
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer DialogToApplyFeatureAction (DialoG d)
{
  ApplyFeatureActionDlgPtr dlg;
  ApplyFeatureActionPtr    action;
  CharPtr                  num_txt;
  LocationIntervalPtr      lint;
  ValNodePtr               vnp;
  ApplyFeatureDetailsPtr   details;

  dlg = (ApplyFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  vnp = DialogToPointer (dlg->feature_type_dlg);
  if (vnp == NULL) {
    return NULL;
  }

  action = ApplyFeatureActionNew();

  action->type = vnp->choice;
  vnp = ValNodeFree (vnp);

  action->apply_to_parts = GetStatus (dlg->apply_to_parts);
  if (action->apply_to_parts && !TextHasNoText (dlg->only_this_part)) {
    num_txt = SaveStringFromText (dlg->only_this_part);
    action->only_seg_num = atoi (num_txt);
    num_txt = MemFree (num_txt);
  }

  action->partial5 = GetStatus (dlg->partial5);
  action->partial3 = GetStatus (dlg->partial3);
  if (GetValue (dlg->strand_group) == 1) {
    action->plus_strand = TRUE;
  } else {
    action->plus_strand = FALSE;
  }

  if (GetValue (dlg->use_whole_interval) == 1) {
    action->location = ValNodeNew(NULL);
    action->location->choice = LocationChoice_whole_sequence;
  } else {
    num_txt = SaveStringFromText (dlg->left_end);
    if (num_txt != NULL && num_txt[StringLen(num_txt) - 1] == '^') {
      action->location = ValNodeNew(NULL);
      action->location->choice = LocationChoice_point;
      action->location->data.intvalue = atoi (num_txt);
      num_txt = MemFree (num_txt);
    } else {
      num_txt = MemFree (num_txt);
      num_txt = SaveStringFromText (dlg->right_end);
      if (num_txt != NULL && num_txt[0] == '^') {
        action->location = ValNodeNew(NULL);
        action->location->choice = LocationChoice_point;
        action->location->data.intvalue = atoi (num_txt + 1) - 1;
        num_txt = MemFree (num_txt);
      } else {
        num_txt = MemFree (num_txt);
        lint = LocationIntervalNew ();
        num_txt = SaveStringFromText (dlg->left_end);
        lint->from = atoi (num_txt);
        num_txt = MemFree (num_txt);
        num_txt = SaveStringFromText (dlg->right_end);
        lint->to = atoi (num_txt);
        num_txt = MemFree (num_txt);
        action->location = ValNodeNew(NULL);
        action->location->choice = LocationChoice_interval;
        action->location->data.ptrvalue = lint;
      }
    }
  }
  if (GetValue (dlg->all_or_some_group) == 1) {
    action->seq_list = ValNodeNew(NULL);
    action->seq_list->choice = SequenceListChoice_all;
  } else {
    num_txt = SaveStringFromText (dlg->accession_list_txt);
    action->seq_list = ValNodeNew(NULL);
    action->seq_list->choice = SequenceListChoice_list;
    action->seq_list->data.ptrvalue = SequenceListCollect (num_txt);
    num_txt = MemFree (num_txt);
  }

  action->add_redundant = GetStatus (dlg->add_redundant);

  details = DialogToPointer (dlg->details);
  action->add_mrna = details->add_mrna;
  action->fields = details->fields;
  details->fields = NULL;
  action->src_fields = details->src_fields;
  details->src_fields = NULL;
  details = ApplyFeatureDetailsFree (details);

  return action;
}


static void ChangeApplyFeatureActionDlgBtn (ButtoN b)
{
  ApplyFeatureActionDlgPtr dlg;

  dlg = (ApplyFeatureActionDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) return;
  ChangeApplyFeatureActionDlg (dlg);
}


static void ChangeApplyFeatureActionDlgGrp (GrouP g)
{
  ApplyFeatureActionDlgPtr dlg;

  dlg = (ApplyFeatureActionDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  ChangeApplyFeatureActionDlg (dlg);
}


static DialoG 
ApplyFeatureActionDialog 
(GrouP h,
 ApplyFeatureActionPtr action,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 Nlm_ChangeNotifyProc     redraw_notify,
 Pointer                  redraw_userdata)
{
  ApplyFeatureActionDlgPtr dlg;
  GrouP                 p, g, parts_group, x, indexer_only_group;
  GrouP                 r2, r3, r4;
  ValNodePtr            feature_type_list = NULL;
  ApplyFeatureDetailsPtr details;

  dlg = (ApplyFeatureActionDlgPtr) MemNew (sizeof (ApplyFeatureActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ApplyFeatureActionToDialog;
  dlg->fromdialog = DialogToApplyFeatureAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestApplyFeatureActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  ValNodeAddPointer (&feature_type_list, Macro_feature_type_cds, StringSave ("CDS"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_gene, StringSave ("gene"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_preRNA, StringSave ("precursor RNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_mRNA, StringSave ("mRNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_tRNA, StringSave ("tRNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_rRNA, StringSave ("rRNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_snRNA, StringSave ("snRNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_scRNA, StringSave ("scRNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_otherRNA, StringSave ("misc RNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_ncRNA, StringSave ("ncRNA"));
  ValNodeAddPointer (&feature_type_list, Macro_feature_type_tmRNA, StringSave ("tmRNA"));
  AddImportFeaturesToChoiceList (&feature_type_list);

  /* note - the ValNodeSelectionDialog will free feature_type_list when done */                                            
  dlg->feature_type_dlg = ValNodeSelectionDialog (p, feature_type_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "feature type", 
                                redraw_notify, redraw_userdata, FALSE);

  parts_group = HiddenGroup (p, 1, 0, NULL);
  dlg->apply_to_parts = CheckBox (parts_group, "Apply to segmented parts, not segmented sequence", ChangeApplyFeatureActionDlgBtn);
  SetObjectExtra (dlg->apply_to_parts, dlg, NULL);
  x = HiddenGroup (parts_group, 2, 0, NULL);
  StaticPrompt (x, "Apply only to particular numbered segment", 0, dialogTextHeight, programFont, 'l');
  dlg->only_this_part = DialogText (x, "", 4, NULL);
  Disable (dlg->only_this_part);

  
  g = HiddenGroup (p, 2, 0, NULL);
  dlg->partial5 = CheckBox (g, "Incomplete at 5' end", NULL);
  dlg->partial3 = CheckBox (g, "Incomplete at 3' end", NULL);

  /* group for strand */
  dlg->strand_group = HiddenGroup (p, 2, 0, NULL);
  SetObjectExtra (dlg->strand_group, dlg, NULL);
  RadioButton (dlg->strand_group, "Plus Strand");
  RadioButton (dlg->strand_group, "Minus Strand");
  SetValue (dlg->strand_group, 1);

  /* coordinates */
  if (indexer_version)
  {
    indexer_only_group = HiddenGroup (p, -1, 0, NULL);
    r2 = HiddenGroup (indexer_only_group, 5, 0, NULL);
    dlg->use_whole_interval = HiddenGroup (r2, 0, 2, ChangeApplyFeatureActionDlgGrp);
    SetObjectExtra (dlg->use_whole_interval, dlg, NULL);
    RadioButton (dlg->use_whole_interval, "Use Whole Sequence Interval");
    RadioButton (dlg->use_whole_interval, "Use these coordinates:");
    r3 = HiddenGroup (r2, 0, 2, NULL);
    StaticPrompt (r3, "", 0, dialogTextHeight, programFont, 'l');
    r4 = HiddenGroup (r3, 4, 0, NULL);
    StaticPrompt (r4, "From", 0, dialogTextHeight, programFont, 'l');
    dlg->left_end = DialogText (r4, "1", 5, NULL);
    StaticPrompt (r4, "To", 0, dialogTextHeight, programFont, 'l');
    dlg->right_end = DialogText (r4, "1", 5, NULL);
    SetValue (dlg->use_whole_interval, 1);
    
    /* apply to some sequences or all sequences */
    dlg->all_or_some_group = HiddenGroup (indexer_only_group, 1, 0, ChangeApplyFeatureActionDlgGrp);
    SetObjectExtra (dlg->all_or_some_group, dlg, NULL);
    RadioButton (dlg->all_or_some_group, "Apply to all sequences");
    RadioButton (dlg->all_or_some_group, "Apply to sequences in this list");
    dlg->accession_list_txt = DialogText (dlg->all_or_some_group, "", 25, NULL);
    SetValue (dlg->all_or_some_group, 1);
    
    dlg->add_redundant = CheckBox (indexer_only_group, "Add even if feature of same type already present", ChangeApplyFeatureActionDlgBtn);
    SetObjectExtra (dlg->add_redundant, dlg, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) r2, (HANDLE) dlg->all_or_some_group, (HANDLE) dlg->add_redundant, NULL);
  }  

  details = ApplyFeatureDetailsNew (action);
  dlg->details = ApplyFeatureDetailsDialog (p, action == NULL ? FEATDEF_CDS : GetFeatdefFromFeatureType(action->type),
                                            details, indexer_version, change_notify, change_userdata);
  details = ApplyFeatureDetailsFree (details);
                                            
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_type_dlg,
                              (HANDLE) parts_group,
                              (HANDLE) g,
                              (HANDLE) dlg->strand_group,
                              (HANDLE) indexer_only_group,
                              (HANDLE) dlg->details,
                              NULL);

  return (DialoG) p;  
}


typedef struct removefeatureactiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG                   feature_type_dlg;
  DialoG                   constraint_dlg;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} RemoveFeatureActionDlgData, PNTR RemoveFeatureActionDlgPtr;


static void RemoveFeatureActionToDialog (DialoG d, Pointer data)
{
  RemoveFeatureActionDlgPtr dlg;
  RemoveFeatureActionPtr    action;
  ValNode                   vn;

  dlg = (RemoveFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  action = (RemoveFeatureActionPtr) data;
  if (action == NULL) {
    PointerToDialog (dlg->feature_type_dlg, NULL);
    PointerToDialog (dlg->constraint_dlg, NULL);
  } else {
    vn.choice = (Uint1)action->type;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->feature_type_dlg, &vn);
    PointerToDialog (dlg->constraint_dlg, action->constraint);
  }

}


static Pointer DialogToRemoveFeatureAction (DialoG d)
{
  RemoveFeatureActionDlgPtr dlg;
  RemoveFeatureActionPtr    action;
  ValNodePtr                vnp;

  dlg = (RemoveFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  vnp = DialogToPointer (dlg->feature_type_dlg);
  if (vnp == NULL) return NULL;

  action = RemoveFeatureActionNew ();
  action->type = vnp->choice;
  vnp = ValNodeFree (vnp);
  action->constraint = DialogToPointer (dlg->constraint_dlg);
  return action;
}


static ValNodePtr TestRemoveFeatureActionDialog (DialoG d)
{
  RemoveFeatureActionDlgPtr dlg;
  dlg = (RemoveFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  return TestDialog (dlg->feature_type_dlg);
}


static DialoG 
RemoveFeatureActionDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  RemoveFeatureActionDlgPtr dlg;
  GrouP                     p;

  dlg = (RemoveFeatureActionDlgPtr) MemNew (sizeof (RemoveFeatureActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = RemoveFeatureActionToDialog;
  dlg->fromdialog = DialogToRemoveFeatureAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestRemoveFeatureActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->feature_type_dlg = FeatureTypeDialogMulti (p, change_notify, change_userdata);
  dlg->constraint_dlg = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_type_dlg, (HANDLE) dlg->constraint_dlg, NULL);

  return (DialoG) p;
}


typedef struct convertfromcdsoptionsdlg {
  DIALOG_MESSAGE_BLOCK
  ButtoN remove_mRNA;
  ButtoN remove_gene;
  ButtoN remove_transcript_id;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} ConvertFromCDSOptionsDlgData, PNTR ConvertFromCDSOptionsDlgPtr;

static void PointerToConvertFromCDSOptionsDialog (DialoG d, Pointer data)
{
  ConvertFromCDSOptionsDlgPtr dlg;
  ConvertFromCDSOptionsPtr    options;

  dlg = (ConvertFromCDSOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  options = (ConvertFromCDSOptionsPtr) data;
  if (options == NULL) {
    SetStatus (dlg->remove_mRNA, FALSE);
    SetStatus (dlg->remove_gene, FALSE);
    SetStatus (dlg->remove_transcript_id, FALSE);
  } else {
    SetStatus (dlg->remove_mRNA, options->remove_mRNA);
    SetStatus (dlg->remove_gene, options->remove_gene);
    SetStatus (dlg->remove_transcript_id, options->remove_transcript_id);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}

static Pointer ConvertFromCDSOptionsDialogToPointer (DialoG d)
{
  ConvertFromCDSOptionsDlgPtr dlg;
  ConvertFromCDSOptionsPtr    options;

  dlg = (ConvertFromCDSOptionsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  options = ConvertFromCDSOptionsNew ();
  options->remove_mRNA = GetStatus (dlg->remove_mRNA);
  options->remove_gene = GetStatus (dlg->remove_gene);
  options->remove_transcript_id = GetStatus (dlg->remove_transcript_id);
  return (Pointer) options;
}


static void ChangeConvertFromCDSOptionsBtn (ButtoN b)
{
  ConvertFromCDSOptionsDlgPtr dlg;
  dlg = (ConvertFromCDSOptionsDlgPtr) GetObjectExtra (b);
  if (dlg != NULL) {
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static DialoG ConvertFromCDSOptionsDialog
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ConvertFromCDSOptionsDlgPtr dlg;
  GrouP                       p;

  dlg = (ConvertFromCDSOptionsDlgPtr) MemNew (sizeof (ConvertFromCDSOptionsDlgData));
  p = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToConvertFromCDSOptionsDialog;
  dlg->fromdialog = ConvertFromCDSOptionsDialogToPointer;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->remove_mRNA = CheckBox (p, "Remove overlapping mRNA", ChangeConvertFromCDSOptionsBtn);
  SetObjectExtra (dlg->remove_mRNA, dlg, NULL);
  dlg->remove_gene = CheckBox (p, "Remove overlapping gene", ChangeConvertFromCDSOptionsBtn);
  SetObjectExtra (dlg->remove_gene, dlg, NULL);
  dlg->remove_transcript_id = CheckBox (p, "Remove transcript ID", ChangeConvertFromCDSOptionsBtn);
  SetObjectExtra (dlg->remove_transcript_id, dlg, NULL);

  return (DialoG) p;
}


typedef struct convertfeatureactiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG                   feature_type_from_dlg;
  DialoG                   feature_type_to_dlg;
  PrompT                   supported;

  DialoG                   cds_options;
  DialoG                   bond_type;
  DialoG                   site_type;
  ButtoN                   create_prot_regions;
  DialoG                   ncrna_class_dlg;

  ButtoN                   leave_original;

  DialoG                   constraint_dlg;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} ConvertFeatureActionDlgData, PNTR ConvertFeatureActionDlgPtr;



static void ChangeConvertFeatureType (Pointer data)
{
  ConvertFeatureActionDlgPtr dlg;
  ValNodePtr                 vnp;
  Uint2                      type_from = 0, type_to = 0;

  dlg = (ConvertFeatureActionDlgPtr) data;
  if (dlg == NULL) return;

  /* if from choice is coding region, show coding region options */
  vnp = DialogToPointer (dlg->feature_type_from_dlg);
  if (vnp == NULL || vnp->choice != Macro_feature_type_cds) {
    Hide (dlg->cds_options);
  } else {
    Show (dlg->cds_options);
  }
  if (vnp != NULL) {
    type_from = vnp->choice;
  }
  vnp = ValNodeFree (vnp);

  /* dst options */
  /* if to choice is bond, show bond options */
  /* if to choice is site, show site options */
  /* if to choice is region, show region checkbox */
  /* if to choice is ncRNA, show ncRNA_class choice */
  vnp = DialogToPointer (dlg->feature_type_to_dlg);
  Hide (dlg->bond_type);
  Hide (dlg->site_type);
  Hide (dlg->create_prot_regions);
  Hide (dlg->ncrna_class_dlg);
  if (vnp == NULL) {
    /* do nothing */
  } else if (vnp->choice == Macro_feature_type_bond) {
    Show (dlg->bond_type);
  } else if (vnp->choice == Macro_feature_type_site) {
    Show (dlg->site_type);
  } else if (vnp->choice == Macro_feature_type_region) {
    Show (dlg->create_prot_regions);
  } else if (vnp->choice == Macro_feature_type_ncRNA) {
    Show (dlg->ncrna_class_dlg);
  }

  if (vnp != NULL) {
    type_to = vnp->choice;
  }
  vnp = ValNodeFree (vnp);

  if (type_to == 0 || type_from == 0 || IsConversionSupported (type_from, type_to)) {
    SetTitle (dlg->supported, "");
  } else {
    SetTitle (dlg->supported, "Conversion is not supported");
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ConvertFeatureActionToDialog (DialoG d, Pointer data)
{
  ConvertFeatureActionDlgPtr dlg;
  ConvertFeatureActionPtr action;
  ValNode                 vn;
  RegionTypePtr           r;

  dlg = (ConvertFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  action = (ConvertFeatureActionPtr) data; 
  if (action == NULL) {
    PointerToDialog (dlg->feature_type_from_dlg, NULL);
    PointerToDialog (dlg->feature_type_to_dlg, NULL);
    PointerToDialog (dlg->cds_options, NULL);
    SetStatus (dlg->leave_original, FALSE);
  } else {
    vn.next = NULL;
    vn.data.ptrvalue = NULL;
    vn.choice = (Uint1)action->type_from;
    PointerToDialog (dlg->feature_type_from_dlg, &vn);
    vn.choice = (Uint1)action->type_to;
    PointerToDialog (dlg->feature_type_to_dlg, &vn);
    /* set source options */
    if (action->src_options != NULL && action->src_options->choice == ConvertFeatureSrcOptions_cds) {
      PointerToDialog (dlg->cds_options, action->src_options->data.ptrvalue);
    } else {
      PointerToDialog (dlg->cds_options, NULL);
    }
    /* set dest options */
    if (action->dst_options == NULL) {
      PointerToDialog (dlg->bond_type, NULL);
      PointerToDialog (dlg->site_type, NULL);
      SetStatus (dlg->create_prot_regions, TRUE);
    } else {
      switch (action->dst_options->choice) {
        case ConvertFeatureDstOptions_bond:
          vn.choice = action->dst_options->data.intvalue;
          PointerToDialog (dlg->bond_type, &vn);
          break;
        case ConvertFeatureDstOptions_site:
          vn.choice = action->dst_options->data.intvalue;
          PointerToDialog (dlg->site_type, &vn);
          break;
        case ConvertFeatureDstOptions_region:
          r = (RegionTypePtr) action->dst_options->data.ptrvalue;
          if (r == NULL || !r->create_nucleotide) {
            SetStatus (dlg->create_prot_regions, TRUE);
          } else {
            SetStatus (dlg->create_prot_regions, FALSE);
          }
          break;
        case ConvertFeatureDstOptions_ncrna_class:
          PointerToDialog (dlg->ncrna_class_dlg, action->dst_options->data.ptrvalue);
          break;
      }
    }
    SetStatus (dlg->leave_original, action->leave_original);
    PointerToDialog (dlg->constraint_dlg, action->src_feat_constraint);
  }
  ChangeConvertFeatureType (dlg);
}


static Pointer DialogToConvertFeatureAction (DialoG d)
{
  ConvertFeatureActionDlgPtr dlg;
  ConvertFeatureActionPtr    action;
  ConvertFromCDSOptionsPtr   src_options;
  ValNodePtr                 vnp;
  RegionTypePtr              r;

  dlg = (ConvertFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  action = ConvertFeatureActionNew ();
  vnp = DialogToPointer (dlg->feature_type_from_dlg);
  if (vnp == NULL) {
    action = ConvertFeatureActionFree (action);
    return NULL;
  } else {
    action->type_from = vnp->choice;
  }
  vnp = ValNodeFree (vnp);

  vnp = DialogToPointer (dlg->feature_type_to_dlg);
  if (vnp == NULL) {
    action = ConvertFeatureActionFree (action);
    return NULL;
  } else {
    action->type_to = vnp->choice;
  }
  vnp = ValNodeFree (vnp);
  
  action->leave_original = GetStatus (dlg->leave_original);

  if (action->type_from == Macro_feature_type_cds) {
    src_options = (ConvertFromCDSOptionsPtr) DialogToPointer (dlg->cds_options);
    if (src_options != NULL) {
      ValNodeAddPointer (&(action->src_options), ConvertFeatureSrcOptions_cds, src_options);
    }
  }

  /* if type to is bond, site, or region, or ncRNA, get options */
  if (action->type_to == Macro_feature_type_bond) {
    vnp = DialogToPointer (dlg->bond_type);
    if (vnp == NULL) {
      action = ConvertFeatureActionFree (action);
      return NULL;
    } else {
      ValNodeAddInt (&(action->dst_options), ConvertFeatureDstOptions_bond, vnp->choice);
    }
    vnp = ValNodeFree (vnp);
  } else if (action->type_to == Macro_feature_type_site) {
    vnp = DialogToPointer (dlg->site_type);
    if (vnp == NULL) {
      action = ConvertFeatureActionFree (action);
      return NULL;
    } else {
      ValNodeAddInt (&(action->dst_options), ConvertFeatureDstOptions_site, vnp->choice);
    }
    vnp = ValNodeFree (vnp);
  } else if (action->type_to == Macro_feature_type_region) {
    r = RegionTypeNew ();
    r->create_nucleotide = !GetStatus (dlg->create_prot_regions);
    ValNodeAddPointer (&(action->dst_options), ConvertFeatureDstOptions_region, r);
  } else if (action->type_to == Macro_feature_type_ncRNA) {
    ValNodeAddPointer (&(action->dst_options), ConvertFeatureDstOptions_ncrna_class, DialogToPointer (dlg->ncrna_class_dlg));
  } else if (action->type_from == Macro_feature_type_cds && action->type_to == Macro_feature_type_mat_peptide_aa) {
    /* hack for converting from coding region to mat-peptide */
    action->dst_options = ValNodeNew (NULL);
    action->dst_options->choice = ConvertFeatureDstOptions_remove_original;
    action->dst_options->data.boolvalue = !action->leave_original;
  }


  action->src_feat_constraint = DialogToPointer (dlg->constraint_dlg);

  return action;
}


static ValNodePtr TestConvertFeatureActionDialog (DialoG d)
{
  ConvertFeatureActionDlgPtr dlg;
  ValNodePtr err_list = NULL, vnp;
  Uint2     type_from = 0, type_to = 0;

  dlg = (ConvertFeatureActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  ValNodeLink (&err_list, TestDialog (dlg->feature_type_from_dlg));
  ValNodeLink (&err_list, TestDialog (dlg->feature_type_to_dlg));

  vnp = DialogToPointer (dlg->feature_type_from_dlg);
  if (vnp != NULL) {
    type_from = vnp->choice;
    if (type_from == Macro_feature_type_cds) {
      ValNodeLink (&err_list, TestDialog (dlg->cds_options));
    }
  }
  vnp = ValNodeFree (vnp);

  vnp = DialogToPointer (dlg->feature_type_to_dlg);
  if (vnp == NULL) {
    /* do nothing */
  } else {
    type_to = vnp->choice;
    if (type_to == Macro_feature_type_bond) {
      ValNodeLink (&err_list, TestDialog (dlg->bond_type));
    } else if (type_to == Macro_feature_type_site) {
      ValNodeLink (&err_list, TestDialog (dlg->site_type));
    }
  }
  vnp = ValNodeFree (vnp);
  if (!IsConversionSupported (type_from, type_to)) {
    ValNodeAddPointer (&err_list, 0, "Unsupported conversion");
  }
  return err_list;  
}

static void ChangeConvertFeatureActionBtn (ButtoN b)
{
  ConvertFeatureActionDlgPtr dlg;

  dlg = (ConvertFeatureActionDlgPtr) GetObjectExtra (b);
  if (dlg != NULL) {
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}

  
static DialoG ConvertFeatureActionDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ConvertFeatureActionDlgPtr dlg;
  GrouP                      p, g1, g2;
  ValNodePtr                 val_list;

  dlg = (ConvertFeatureActionDlgPtr) MemNew (sizeof (ConvertFeatureActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ConvertFeatureActionToDialog;
  dlg->fromdialog = DialogToConvertFeatureAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestConvertFeatureActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g1 = HiddenGroup (p, 2, 0, NULL);
  dlg->feature_type_from_dlg = FeatureTypeDialog (g1, ChangeConvertFeatureType, dlg);
  dlg->feature_type_to_dlg = FeatureTypeDialog (g1, ChangeConvertFeatureType, dlg);

  dlg->supported = StaticPrompt (p, "Conversion is not supported", 0, dialogTextHeight, systemFont, 'l');

  dlg->cds_options = ConvertFromCDSOptionsDialog (p, change_notify, change_userdata);
  Hide (dlg->cds_options);

  g2 = HiddenGroup (p, 0, 0, NULL);
  val_list = GetBondTypeList ();
  dlg->bond_type = ValNodeSelectionDialog (g2, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "bond_type", 
                                           change_notify, change_userdata, FALSE);
  Hide (dlg->bond_type);
  val_list = GetSiteTypeList ();
  dlg->site_type = ValNodeSelectionDialog (g2, val_list, SHORT_SELECTION_LIST, ValNodeStringName,
                                           ValNodeSimpleDataFree, ValNodeStringCopy,
                                           ValNodeChoiceMatch, "site_type", 
                                           change_notify, change_userdata, FALSE);
  Hide (dlg->site_type);
  dlg->create_prot_regions = CheckBox (g2, "Create region features on protein sequence of overlapping coding region", ChangeConvertFeatureActionBtn);
  SetObjectExtra (dlg->create_prot_regions, dlg, NULL);
  Hide (dlg->create_prot_regions);

  dlg->ncrna_class_dlg = CreatencRNAClassDialog (g2, FALSE, change_notify, change_userdata);
  Hide (dlg->ncrna_class_dlg);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->bond_type,
                              (HANDLE) dlg->site_type,
                              (HANDLE) dlg->create_prot_regions,
                              (HANDLE) dlg->ncrna_class_dlg,
                              NULL);

  dlg->leave_original = CheckBox (p, "Leave original", ChangeConvertFeatureActionBtn);
  SetObjectExtra (dlg->leave_original, dlg, NULL);

  dlg->constraint_dlg = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->supported, (HANDLE) g1, (HANDLE) dlg->cds_options, (HANDLE) g2, (HANDLE) dlg->leave_original, (HANDLE) dlg->constraint_dlg, NULL);

  return (DialoG) p;
}


static DialoG DescriptorTypeDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ValNodePtr desc_list = NULL;

  AddAllDescriptorsToChoiceList (&desc_list);
  return ValNodeSelectionDialog (h, desc_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "descriptor type", 
                                change_notify, change_userdata, FALSE);
}


typedef struct removedescriptoractiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG                   descriptor_type_dlg;

  DialoG                   constraint_dlg;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} RemoveDescriptorActionDlgData, PNTR RemoveDescriptorActionDlgPtr;


static void RemoveDescriptorActionToDialog (DialoG d, Pointer data)
{
  RemoveDescriptorActionDlgPtr dlg;
  RemoveDescriptorActionPtr    action;
  ValNode vn;

  dlg = (RemoveDescriptorActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  vn.data.ptrvalue = NULL;
  vn.next = NULL;

  action = (RemoveDescriptorActionPtr) data;
  if (action == NULL) {
    vn.choice = Descriptor_type_all;
    PointerToDialog (dlg->descriptor_type_dlg, &vn);
    PointerToDialog (dlg->constraint_dlg, NULL);
  } else {
    vn.choice = (Uint1)action->type;
    PointerToDialog (dlg->descriptor_type_dlg, &vn);
    PointerToDialog (dlg->constraint_dlg, action->constraint);
  }
}


static Pointer DialogToRemoveDescriptorAction (DialoG d)
{
  RemoveDescriptorActionDlgPtr dlg;
  RemoveDescriptorActionPtr    action;
  ValNodePtr                   vnp;

  dlg = (RemoveDescriptorActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = RemoveDescriptorActionNew();

  vnp = DialogToPointer (dlg->descriptor_type_dlg);
  if (vnp == NULL) {
    action->type = Descriptor_type_all;
  } else {
    action->type = vnp->choice;
  }
  vnp = ValNodeFree (vnp);
  action->constraint = DialogToPointer (dlg->constraint_dlg);
  return action;
}


static ValNodePtr TestRemoveDescriptorActionDialog (DialoG d)
{
  RemoveDescriptorActionDlgPtr dlg;

  dlg = (RemoveDescriptorActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  return TestDialog (dlg->descriptor_type_dlg);
}


static DialoG RemoveDescriptorActionDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  RemoveDescriptorActionDlgPtr dlg;
  GrouP                      p;

  dlg = (RemoveDescriptorActionDlgPtr) MemNew (sizeof (RemoveDescriptorActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = RemoveDescriptorActionToDialog;
  dlg->fromdialog = DialogToRemoveDescriptorAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestRemoveDescriptorActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->descriptor_type_dlg = DescriptorTypeDialog (p, change_notify, change_userdata);

  dlg->constraint_dlg = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->descriptor_type_dlg, (HANDLE) dlg->constraint_dlg, NULL);

  return (DialoG) p;
}


typedef struct autodefactiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG                   modifiers_dlg;
  PopuP                    clause_list_type_popup;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} AutodefActionDlgData, PNTR AutodefActionDlgPtr;


static void AutodefActionToDialog(DialoG d, Pointer data)
{
  AutodefActionDlgPtr dlg;
  AutodefActionPtr    action;
  ValNodePtr          list = NULL, vnp;

  dlg = (AutodefActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  action = (AutodefActionPtr) data;
  if (action == NULL) {
    PointerToDialog (dlg->modifiers_dlg, NULL);
    SetValue (dlg->clause_list_type_popup, 2);
  } else {
    /* populate modifiers */
    for (vnp = action->modifiers; vnp != NULL; vnp = vnp->next) {
      ValNodeAddPointer (&list, 0, GetSourceQualName (vnp->data.intvalue));
    }
    PointerToDialog (dlg->modifiers_dlg, list);

    /* populate clause list type */
    switch (action->clause_list_type) {
      case Autodef_list_type_feature_list:
        SetValue (dlg->clause_list_type_popup, 1);
        break;
      case Autodef_list_type_complete_sequence:
        SetValue (dlg->clause_list_type_popup, 2);
        break;
      case Autodef_list_type_complete_genome:
        SetValue (dlg->clause_list_type_popup, 3);
        break;
      default:
        SetValue (dlg->clause_list_type_popup, 2);
        break;
    }
  }
}


static Pointer DialogToAutodefAction (DialoG d)
{
  AutodefActionDlgPtr dlg;
  AutodefActionPtr    action = NULL;
  Int2                v;
  ValNodePtr          list, vnp;

  dlg = (AutodefActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = AutodefActionNew();
  v = GetValue (dlg->clause_list_type_popup);
  switch (v) {
    case 1:
      action->clause_list_type = Autodef_list_type_feature_list;
      break;
    case 2:
      action->clause_list_type = Autodef_list_type_complete_sequence;
      break;
    case 3:
      action->clause_list_type = Autodef_list_type_complete_genome;
      break;
    default:
      action->clause_list_type = Autodef_list_type_feature_list;
      break;
  }

  list = DialogToPointer (dlg->modifiers_dlg);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    ValNodeAddInt (&(action->modifiers), SourceQualChoice_textqual, GetSourceQualTypeByName (vnp->data.ptrvalue));
  }
  list = ValNodeFree (list);

  return action;
}


static void ChangeAutodefPopup (PopuP p)
{
  AutodefActionDlgPtr dlg;

  dlg = (AutodefActionDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG AutodefActionDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  AutodefActionDlgPtr dlg;
  GrouP               p;
  PrompT              ppt1, ppt2;
  ValNodePtr          src_quals = NULL;
  ValNode             vn;

  dlg = (AutodefActionDlgPtr) MemNew (sizeof (AutodefActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = AutodefActionToDialog;
  dlg->fromdialog = DialogToAutodefAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = NULL;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  ppt1 = StaticPrompt (p, "Modifiers", 0, dialogTextHeight, systemFont, 'l');

  ValNodeAddPointer (&src_quals, 0, StringSave ("clone"));
  ValNodeAddPointer (&src_quals, 0, StringSave ("cultivar"));
  ValNodeAddPointer (&src_quals, 0, StringSave ("culture-collection"));
  ValNodeAddPointer (&src_quals, 0, StringSave ("haplogroup"));
  ValNodeAddPointer (&src_quals, 0, StringSave ("isolate"));
  ValNodeAddPointer (&src_quals, 0, StringSave ("strain"));
  ValNodeAddPointer (&src_quals, 0, StringSave ("specimen-voucher"));

  dlg->modifiers_dlg = ValNodeSelectionDialog (p, src_quals, SHORT_SELECTION_LIST,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "source qual",
                                               change_notify, change_userdata, TRUE);
  vn.choice = 0;
  vn.next = NULL;
  vn.data.ptrvalue = "haplogroup";
  PointerToDialog (dlg->modifiers_dlg, &vn);

  ppt2 = StaticPrompt (p, "Features or Complete", 0, dialogTextHeight, systemFont, 'l');

  dlg->clause_list_type_popup = PopupList (p, TRUE, ChangeAutodefPopup);
  SetObjectExtra (dlg->clause_list_type_popup, dlg, NULL);
  PopupItem (dlg->clause_list_type_popup, "List Features");
  PopupItem (dlg->clause_list_type_popup, "Complete Sequence");
  PopupItem (dlg->clause_list_type_popup, "Complete Genome");
  SetValue (dlg->clause_list_type_popup, 2);

  AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) dlg->modifiers_dlg, (HANDLE) ppt2, (HANDLE) dlg->clause_list_type_popup, NULL);

  return (DialoG) p;
}


typedef struct fixpubcapsdlg {
  DIALOG_MESSAGE_BLOCK

  ButtoN fix_title;
  ButtoN fix_author;
  ButtoN fix_affil;
  ButtoN fix_affil_country;
  ButtoN punct_only;
  DialoG constraint;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} FixPubCapsDlgData, PNTR FixPubCapsDlgPtr;


static void FixPubCapsActionToDialog (DialoG d, Pointer data)
{
  FixPubCapsDlgPtr dlg;
  FixPubCapsActionPtr action;

  dlg = (FixPubCapsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  action = (FixPubCapsActionPtr) data;
  if (action == NULL) {
    SetStatus (dlg->fix_title, FALSE);
    SetStatus (dlg->fix_author, FALSE);
    SetStatus (dlg->fix_affil, FALSE);
    SetStatus (dlg->fix_affil_country, FALSE);
    SetStatus (dlg->punct_only, FALSE);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    SetStatus (dlg->fix_title, action->title);
    SetStatus (dlg->fix_author, action->authors);
    SetStatus (dlg->fix_affil, action->affiliation);
    SetStatus (dlg->fix_affil_country, action->affil_country);
    SetStatus (dlg->punct_only, action->punct_only);
    PointerToDialog (dlg->constraint, action->constraint);
  }

}


static Pointer DialogToFixPubCapsAction (DialoG d)
{
  FixPubCapsDlgPtr dlg;
  FixPubCapsActionPtr action = NULL;

  dlg = (FixPubCapsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = FixPubCapsActionNew ();
  action->title = GetStatus (dlg->fix_title);
  action->authors = GetStatus (dlg->fix_author);
  action->affiliation = GetStatus (dlg->fix_affil);
  action->affil_country = GetStatus (dlg->fix_affil_country);
  action->punct_only = GetStatus (dlg->punct_only);
  if (IsFixPubCapsActionEmpty(action)) {
    action = FixPubCapsActionFree (action);
  } else {
    action->constraint = DialogToPointer (dlg->constraint);
  }
  return action;
}


static void ChangeFixPubCapsActionBtn (ButtoN b)
{
  FixPubCapsDlgPtr dlg;

  dlg = (FixPubCapsDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static ValNodePtr TestFixPubCapsActionDialog (DialoG d)
{
  FixPubCapsDlgPtr dlg;
  FixPubCapsActionPtr a;
  ValNodePtr err_list = NULL;

  dlg = (FixPubCapsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  a = (FixPubCapsActionPtr) DialogToPointer (d);
  if (a == NULL) {
    ValNodeAddPointer (&err_list, 0, "bad action");
  }
  return err_list;
}


static DialoG FixPubCapsDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  FixPubCapsDlgPtr dlg;
  GrouP            p, g1;

  dlg = (FixPubCapsDlgPtr) MemNew (sizeof (FixPubCapsDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FixPubCapsActionToDialog;
  dlg->fromdialog = DialogToFixPubCapsAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestFixPubCapsActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g1 = HiddenGroup (p, 0, 5, NULL);
  SetGroupSpacing (g1, 10, 10);
  dlg->fix_title = CheckBox (g1, "Fix title", ChangeFixPubCapsActionBtn);
  SetObjectExtra (dlg->fix_title, dlg, NULL);
  dlg->fix_author = CheckBox (g1, "Fix authors", ChangeFixPubCapsActionBtn);
  SetObjectExtra (dlg->fix_author, dlg, NULL);
  dlg->fix_affil = CheckBox (g1, "Fix affiliation", ChangeFixPubCapsActionBtn);
  SetObjectExtra (dlg->fix_affil, dlg, NULL);
  dlg->fix_affil_country = CheckBox (g1, "Fix affiliation country", ChangeFixPubCapsActionBtn);
  SetObjectExtra (dlg->fix_affil_country, dlg, NULL);
  dlg->punct_only = CheckBox (g1, "Punctuation only", ChangeFixPubCapsActionBtn);
  SetObjectExtra (dlg->punct_only, dlg, NULL);

  dlg->constraint = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) dlg->constraint, NULL);

  return (DialoG) p;
}


typedef struct sortfieldsdlg {
  DIALOG_MESSAGE_BLOCK

  /* Note - for now, the only field you can sort is protein name */
  GrouP  order;
  DialoG constraint;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} SortFieldsDlgData, PNTR SortFieldsDlgPtr;


static void SortFieldsActionToDialog (DialoG d, Pointer data)
{
  SortFieldsDlgPtr dlg;
  SortFieldsActionPtr action;

  dlg = (SortFieldsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  action = (SortFieldsActionPtr) data;
  if (action == NULL) {
    SetValue (dlg->order, 1);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    if (action->order > 0) {
      SetValue (dlg->order, action->order);
    } else {
      SetValue (dlg->order, 1);
    }
    PointerToDialog (dlg->constraint, action->constraint);
  }

}


static Pointer DialogToSortFieldsAction (DialoG d)
{
  SortFieldsDlgPtr dlg;
  SortFieldsActionPtr action = NULL;
  FeatureFieldPtr ffield;

  dlg = (SortFieldsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  action = SortFieldsActionNew ();
  action->order = GetValue (dlg->order);
  action->constraint = DialogToPointer (dlg->constraint);

  ffield = FeatureFieldNew ();
  ffield->type = Macro_feature_type_cds;
  ffield->field = ValNodeNew (NULL);
  ffield->field->choice = FeatQualChoice_legal_qual;
  ffield->field->data.intvalue = Feat_qual_legal_product;
  action->field = ValNodeNew (NULL);
  action->field->choice = FieldType_feature_field;
  action->field->data.ptrvalue = ffield;
  return action;
}


static void ChangeSortFieldsActionGrp (GrouP g)
{
  SortFieldsDlgPtr dlg;

  dlg = (SortFieldsDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static ValNodePtr TestSortFieldsActionDialog (DialoG d)
{
  SortFieldsDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (SortFieldsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  return err_list;
}


static DialoG SortFieldsDialog
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  SortFieldsDlgPtr dlg;
  GrouP            p;

  dlg = (SortFieldsDlgPtr) MemNew (sizeof (SortFieldsDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = SortFieldsActionToDialog;
  dlg->fromdialog = DialogToSortFieldsAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestSortFieldsActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->order = HiddenGroup (p, 0, 3, ChangeSortFieldsActionGrp);
  SetGroupSpacing (dlg->order, 10, 10);
  RadioButton (dlg->order, "By length, short to long");
  RadioButton (dlg->order, "By length, long to short");
  RadioButton (dlg->order, "Alphabetically");
  SetObjectExtra (dlg->order, dlg, NULL);
  SetValue (dlg->order, 1);

  dlg->constraint = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->order, (HANDLE) dlg->constraint, NULL);

  return (DialoG) p;
}


typedef struct editfeaturelocationactiondlg {
  DIALOG_MESSAGE_BLOCK
  DialoG                   feature_type_dlg;
  PopuP                    feature_edit_type;
  GrouP                    strand_edit_grp;
  PopuP                    strand_from;
  PopuP                    strand_to;                   
  GrouP                    set_5_partial_grp;
  PopuP                    set_5_type;
  ButtoN                   extend5;
  PopuP                    clear_5_type;
  GrouP                    set_3_partial_grp;
  PopuP                    set_3_type;
  ButtoN                   extend3;
  PopuP                    clear_3_type;
  GrouP                    set_both_partial_grp;
  PopuP                    set_both_type;
  ButtoN                   extendboth;
  PopuP                    clear_both_type;
  PopuP                    convert_loc;
  ButtoN                   retranslate_cds;
  DialoG                   constraint_dlg;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} EditFeatureLocationActionDlgData, PNTR EditFeatureLocationActionDlgPtr;


static void ChangeFeatureLocationEditType (PopuP p)
{
  EditFeatureLocationActionDlgPtr dlg;
  Int2 val;

  dlg = (EditFeatureLocationActionDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;
  
  Hide (dlg->strand_edit_grp);
  Hide (dlg->set_5_partial_grp);
  Hide (dlg->clear_5_type);
  Hide (dlg->set_3_partial_grp);
  Hide (dlg->clear_3_type);
  Hide (dlg->set_both_partial_grp);
  Hide (dlg->clear_both_type);
  Hide (dlg->convert_loc);

  val = GetValue (dlg->feature_edit_type);
  switch (val) {
    case 1:
      Show (dlg->strand_edit_grp);
      break;
    case 2:
      Show (dlg->set_5_partial_grp);
      break;
    case 3:
      Show (dlg->clear_5_type);
      break;
    case 4:
      Show (dlg->set_3_partial_grp);
      break;
    case 5:
      Show (dlg->clear_3_type);
      break;
    case 6:
      Show (dlg->set_both_partial_grp);
      break;
    case 7:
      Show (dlg->clear_both_type);
      break;
    case 8:
      Show (dlg->convert_loc);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ResetFeatureLocationActionEditType (EditFeatureLocationActionDlgPtr dlg)
{
  if (dlg != NULL) {
    SetValue (dlg->feature_edit_type, 1);
    SetValue (dlg->strand_from, 1);
    SetValue (dlg->strand_to, 1);
    SetValue (dlg->retranslate_cds, FALSE);
  }
}


static void EditFeatureLocationActionToDialog (DialoG d, Pointer data)
{
  EditFeatureLocationActionDlgPtr dlg;
  EditFeatureLocationActionPtr    action;
  ValNode                   vn;
  EditLocationStrandPtr    strand;
  Partial5SetActionPtr     set5;
  Partial3SetActionPtr     set3;
  PartialBothSetActionPtr  setboth;

  dlg = (EditFeatureLocationActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  action = (EditFeatureLocationActionPtr) data;
  if (action == NULL) {
    PointerToDialog (dlg->feature_type_dlg, NULL);
    ResetFeatureLocationActionEditType (dlg);
    PointerToDialog (dlg->constraint_dlg, NULL);
    SetStatus (dlg->retranslate_cds, FALSE);
  } else {
    vn.choice = (Uint1)action->type;
    vn.data.ptrvalue = NULL;
    vn.next = NULL;
    PointerToDialog (dlg->feature_type_dlg, &vn);
    PointerToDialog (dlg->constraint_dlg, action->constraint);
    SetStatus (dlg->retranslate_cds, action->retranslate_cds);
    if (action->action == NULL) {
      ResetFeatureLocationActionEditType (dlg);
    } else {
      switch (action->action->choice) {
        case LocationEditType_strand:
          SetValue (dlg->feature_edit_type, 1);
          strand = (EditLocationStrandPtr) action->action->data.ptrvalue;
          if (strand == NULL) {
            SetValue (dlg->strand_from, 1);
            SetValue (dlg->strand_to, 1);
          } else {
            switch (strand->strand_from) {
              case Feature_location_strand_from_any:
                SetValue (dlg->strand_from, 1);
                break;
              case Feature_location_strand_from_plus:
                SetValue (dlg->strand_from, 2);
                break;
              case Feature_location_strand_from_minus:
                SetValue (dlg->strand_from, 3);
                break;
              case Feature_location_strand_from_unknown:
                SetValue (dlg->strand_from, 4);
                break;
              case Feature_location_strand_from_both:
                SetValue (dlg->strand_from, 5);
                break;
              default:
                SetValue (dlg->strand_from, 1);
                break;
            }
            switch (strand->strand_to) {
              case Feature_location_strand_to_plus:
                SetValue (dlg->strand_to, 1);
                break;
              case Feature_location_strand_to_minus:
                SetValue (dlg->strand_to, 2);
                break;
              case Feature_location_strand_to_unknown:
                SetValue (dlg->strand_to, 3);
                break;
              case Feature_location_strand_to_both:
                SetValue (dlg->strand_to, 4);
                break;
              case Feature_location_strand_to_reverse:
                SetValue (dlg->strand_to, 5);
                break;
              default:
                SetValue (dlg->strand_to, 1);
                break;
            }
          }
          break;
        case LocationEditType_set_5_partial:
          SetValue (dlg->feature_edit_type, 2);
          set5 = (Partial5SetActionPtr) action->action->data.ptrvalue;
          if (set5 == NULL) {
            SetValue (dlg->set_5_type, 1);
            SetStatus (dlg->extend5, FALSE);
          } else {
            switch (set5->constraint) {
              case Partial_5_set_constraint_all:
                SetValue (dlg->set_5_type, 1);
                break;
              case Partial_5_set_constraint_at_end:
                SetValue (dlg->set_5_type, 2);
                break;
              case Partial_5_set_constraint_bad_start:
                SetValue (dlg->set_5_type, 3);
                break;
              case Partial_5_set_constraint_frame_not_one:
                SetValue (dlg->set_5_type, 4);
                break;
              default:
                SetValue (dlg->set_5_type, 1);
                break;
            }
            SetStatus (dlg->extend5, set5->extend);
          }
          break;   
        case LocationEditType_clear_5_partial:
          SetValue (dlg->feature_edit_type, 3);
          switch (action->action->data.intvalue) {
            case Partial_5_clear_constraint_all:
              SetValue (dlg->clear_5_type, 1);
              break;
            case Partial_5_clear_constraint_not_at_end:
              SetValue (dlg->clear_5_type, 2);
              break;
            case Partial_5_clear_constraint_good_start:
              SetValue (dlg->clear_5_type, 3);
              break;
            default:
              SetValue (dlg->clear_5_type, 1);
              break;
          }
          break;
        case LocationEditType_set_3_partial:
          SetValue (dlg->feature_edit_type, 4);
          set3 = (Partial3SetActionPtr) action->action->data.ptrvalue;
          if (set3 == NULL) {
            SetValue (dlg->set_3_type, 1);
            SetStatus (dlg->extend3, FALSE);
          } else {
            switch (set3->constraint) {
              case Partial_3_set_constraint_all:
                SetValue (dlg->set_3_type, 1);
                break;
              case Partial_3_set_constraint_at_end:
                SetValue (dlg->set_3_type, 2);
                break;
              case Partial_3_set_constraint_bad_end:
                SetValue (dlg->set_3_type, 3);
                break;
              default:
                SetValue (dlg->set_3_type, 1);
                break;
            }
            SetStatus (dlg->extend3, set3->extend);
          }
          break;   
        case LocationEditType_clear_3_partial:
          SetValue (dlg->feature_edit_type, 5);
          switch (action->action->data.intvalue) {
            case Partial_3_clear_constraint_all:
              SetValue (dlg->clear_3_type, 1);
              break;
            case Partial_3_clear_constraint_not_at_end:
              SetValue (dlg->clear_3_type, 2);
              break;
            case Partial_3_clear_constraint_good_end:
              SetValue (dlg->clear_3_type, 3);
              break;
            default:
              SetValue (dlg->clear_3_type, 1);
              break;
          }
          break;
        case LocationEditType_set_both_partial:
          SetValue (dlg->feature_edit_type, 6);
          setboth = (PartialBothSetActionPtr) action->action->data.ptrvalue;
          if (setboth == NULL) {
            SetValue (dlg->set_both_type, 1);
            SetStatus (dlg->extendboth, FALSE);
          } else {
            switch (setboth->constraint) {
              case Partial_both_set_constraint_all:
                SetValue (dlg->set_both_type, 1);
                break;
              case Partial_both_set_constraint_at_end:
                SetValue (dlg->set_both_type, 2);
                break;
              default:
                SetValue (dlg->set_both_type, 1);
                break;
            }
            SetStatus (dlg->extendboth, setboth->extend);
          }
          break;
        case LocationEditType_clear_both_partial:
          SetValue (dlg->feature_edit_type, 7);
          switch (action->action->data.intvalue) {
            case Partial_both_clear_constraint_all:
              SetValue (dlg->clear_both_type, 1);
              break;
            case Partial_both_clear_constraint_not_at_end:
              SetValue (dlg->clear_both_type, 2);
              break;
            default:
              SetValue (dlg->clear_both_type, 1);
              break;
          }
          break;
        case LocationEditType_convert:
          SetValue (dlg->feature_edit_type, 8);
          switch (action->action->data.intvalue) {
            case Convert_location_type_join:
              SetValue (dlg->convert_loc, 1);
              break;
            case Convert_location_type_order:
              SetValue (dlg->convert_loc, 2);
              break;
            case Convert_location_type_merge:
              SetValue (dlg->convert_loc, 3);
              break;
            default:
              SetValue (dlg->convert_loc, 1);
              break;
          }
          break;
        case LocationEditType_extend_5:
          SetValue (dlg->feature_edit_type, 9);
          break;
        case LocationEditType_extend_3:
          SetValue (dlg->feature_edit_type, 10);
          break;
        default:
          ResetFeatureLocationActionEditType (dlg);
          break;
      }
    }
  }
  ChangeFeatureLocationEditType (dlg->feature_edit_type);
}


static Pointer DialogToEditFeatureLocationAction (DialoG d)
{
  EditFeatureLocationActionDlgPtr dlg;
  EditFeatureLocationActionPtr    action;
  ValNodePtr                vnp;
  Int2                      val, val2;
  EditLocationStrandPtr    strand;
  Partial5SetActionPtr     set5;
  Partial3SetActionPtr     set3;
  PartialBothSetActionPtr  setboth;

  dlg = (EditFeatureLocationActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  vnp = DialogToPointer (dlg->feature_type_dlg);
  if (vnp == NULL) return NULL;

  action = EditFeatureLocationActionNew ();
  action->type = vnp->choice;
  vnp = ValNodeFree (vnp);
  action->retranslate_cds = GetStatus (dlg->retranslate_cds);
  action->constraint = DialogToPointer (dlg->constraint_dlg);

  val = GetValue (dlg->feature_edit_type);
  switch (val) {
    case 1:
      strand = EditLocationStrandNew ();
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_strand;
      action->action->data.ptrvalue = strand;
      val2 = GetValue (dlg->strand_from);
      switch (val2) {
        case 1:
          strand->strand_from = Feature_location_strand_from_any;
          break;
        case 2:
          strand->strand_from = Feature_location_strand_from_plus;
          break;
        case 3:
          strand->strand_from = Feature_location_strand_from_minus;
          break;
        case 4:
          strand->strand_from = Feature_location_strand_from_unknown;
          break;
        case 5:
          strand->strand_from = Feature_location_strand_from_both;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
      }
      val2 = GetValue (dlg->strand_to);
      switch (val2) {
        case 1:
          strand->strand_to = Feature_location_strand_to_plus;
          break;
        case 2:
          strand->strand_to = Feature_location_strand_to_minus;
          break;
        case 3:
          strand->strand_to = Feature_location_strand_to_unknown;
          break;
        case 4:
          strand->strand_to = Feature_location_strand_to_both;
          break;
        case 5:
          strand->strand_to = Feature_location_strand_to_reverse;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          break;
      }
      break;
    case 2:
      set5 = Partial5SetActionNew ();
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_set_5_partial;
      action->action->data.ptrvalue = set5;
      set5->extend = GetStatus (dlg->extend5);
      val2 = GetValue (dlg->set_5_type);
      switch (val2) {
        case 1:
          set5->constraint = Partial_5_set_constraint_all;
          break;
        case 2:
          set5->constraint = Partial_5_set_constraint_at_end;
          break;
        case 3:
          set5->constraint = Partial_5_set_constraint_bad_start;
          break;
        case 4:
          set5->constraint = Partial_5_set_constraint_frame_not_one;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 3:
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_clear_5_partial;
      val2 = GetValue (dlg->clear_5_type);
      switch (val2) {
        case 1:
          action->action->data.intvalue = Partial_5_clear_constraint_all;
          break;
        case 2:
          action->action->data.intvalue = Partial_5_clear_constraint_not_at_end;
          break;
        case 3:
          action->action->data.intvalue = Partial_5_clear_constraint_good_start;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 4:
      set3 = Partial3SetActionNew ();
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_set_3_partial;
      action->action->data.ptrvalue = set3;
      set3->extend = GetStatus (dlg->extend3);
      val2 = GetValue (dlg->set_3_type);
      switch (val2) {
        case 1:
          set3->constraint = Partial_3_set_constraint_all;
          break;
        case 2:
          set3->constraint = Partial_3_set_constraint_at_end;
          break;
        case 3:
          set3->constraint = Partial_3_set_constraint_bad_end;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 5:
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_clear_3_partial;
      val2 = GetValue (dlg->clear_3_type);
      switch (val2) {
        case 1:
          action->action->data.intvalue = Partial_3_clear_constraint_all;
          break;
        case 2:
          action->action->data.intvalue = Partial_3_clear_constraint_not_at_end;
          break;
        case 3:
          action->action->data.intvalue = Partial_3_clear_constraint_good_end;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 6:
      setboth = PartialBothSetActionNew ();
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_set_both_partial;
      action->action->data.ptrvalue = setboth;
      setboth->extend = GetStatus (dlg->extendboth);
      val2 = GetValue (dlg->set_both_type);
      switch (val2) {
        case 1:
          setboth->constraint = Partial_both_set_constraint_all;
          break;
        case 2:
          setboth->constraint = Partial_both_set_constraint_at_end;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 7:
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_clear_both_partial;
      val2 = GetValue (dlg->clear_both_type);
      switch (val2) {
        case 1:
          action->action->data.intvalue = Partial_both_clear_constraint_all;
          break;
        case 2:
          action->action->data.intvalue = Partial_both_clear_constraint_not_at_end;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 8:
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_convert;
      val2 = GetValue (dlg->convert_loc);
      switch (val2) {
        case 1:
          action->action->data.intvalue = Convert_location_type_join;
          break;
        case 2:
          action->action->data.intvalue = Convert_location_type_order;
          break;
        case 3:
          action->action->data.intvalue = Convert_location_type_merge;
          break;
        default:
          action = EditFeatureLocationActionFree (action);
          return NULL;
          break;
      }
      break;
    case 9:
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_extend_5;
      break;
    case 10:
      action->action = ValNodeNew (NULL);
      action->action->choice = LocationEditType_extend_3;
      break;
    default:
      action = EditFeatureLocationActionFree (action);
      break;
  }

  return action;
}


static ValNodePtr TestEditFeatureLocationActionDialog (DialoG d)
{
  EditFeatureLocationActionDlgPtr dlg;
  EditFeatureLocationActionPtr a;
  ValNodePtr err_list = NULL;

  dlg = (EditFeatureLocationActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  a = (EditFeatureLocationActionPtr) DialogToPointer (d);
  if (a == NULL) {
    ValNodeAddPointer (&err_list, 0, "bad action");
  } else {
    a = EditFeatureLocationActionFree (a);
    err_list = TestDialog (dlg->feature_type_dlg);
  }
  return err_list;
}


static void ChangeEditFeatureLocationPopup (PopuP p)
{
  EditFeatureLocationActionDlgPtr dlg;

  dlg = (EditFeatureLocationActionDlgPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeEditFeatureLocationButton (ButtoN b)
{
  EditFeatureLocationActionDlgPtr dlg;

  dlg = (EditFeatureLocationActionDlgPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG 
EditFeatureLocationActionDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  EditFeatureLocationActionDlgPtr dlg;
  GrouP                     p;
  GrouP                     g;

  dlg = (EditFeatureLocationActionDlgPtr) MemNew (sizeof (EditFeatureLocationActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = EditFeatureLocationActionToDialog;
  dlg->fromdialog = DialogToEditFeatureLocationAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestEditFeatureLocationActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->feature_type_dlg = FeatureTypeDialog (p, change_notify, change_userdata);

  dlg->feature_edit_type = PopupList (p, TRUE, ChangeFeatureLocationEditType);
  SetObjectExtra (dlg->feature_edit_type, dlg, NULL);
  PopupItem (dlg->feature_edit_type, "Edit Strand");
  PopupItem (dlg->feature_edit_type, "Set 5' Partial");
  PopupItem (dlg->feature_edit_type, "Clear 5' Partial");
  PopupItem (dlg->feature_edit_type, "Set 3' Partial");
  PopupItem (dlg->feature_edit_type, "Clear 3' Partial");
  PopupItem (dlg->feature_edit_type, "Set Both Ends Partial");
  PopupItem (dlg->feature_edit_type, "Clear Both Ends Partial");
  PopupItem (dlg->feature_edit_type, "Convert location");
  PopupItem (dlg->feature_edit_type, "Extend 5' end to end of sequence");
  PopupItem (dlg->feature_edit_type, "Extend 3' end to end of sequence");
  SetValue (dlg->feature_edit_type, 1);

  g = HiddenGroup (p, 0, 0, NULL);
  dlg->strand_edit_grp = HiddenGroup (g, 4, 0, NULL);

  /* strand */
  StaticPrompt (dlg->strand_edit_grp, "Convert location strand from", 0, dialogTextHeight, systemFont, 'l');
  dlg->strand_from = PopupList (dlg->strand_edit_grp, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->strand_from, dlg, NULL);
  PopupItem (dlg->strand_from, "Any");
  PopupItem (dlg->strand_from, "Plus");
  PopupItem (dlg->strand_from, "Minus");
  PopupItem (dlg->strand_from, "Unknown");
  PopupItem (dlg->strand_from, "Both");
  SetValue (dlg->strand_from, 1);
  dlg->strand_to = PopupList (dlg->strand_edit_grp, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->strand_to, dlg, NULL);
  PopupItem (dlg->strand_to, "Plus");
  PopupItem (dlg->strand_to, "Minus");
  PopupItem (dlg->strand_to, "Unknown");
  PopupItem (dlg->strand_to, "Both");
  PopupItem (dlg->strand_to, "Reverse");
  SetValue (dlg->strand_to, 5);

  /* partials */
  dlg->set_5_partial_grp = HiddenGroup (g, 2, 0, NULL);
  dlg->set_5_type = PopupList (dlg->set_5_partial_grp, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->set_5_type, dlg, NULL);
  PopupItem (dlg->set_5_type, "All");
  PopupItem (dlg->set_5_type, "If 5' end at end of sequence");
  PopupItem (dlg->set_5_type, "If bad start codon");
  PopupItem (dlg->set_5_type, "If CDS frame > 1");
  SetValue (dlg->set_5_type, 1);
  dlg->extend5 = CheckBox (dlg->set_5_partial_grp, "Extend to 5' end if partial is set", ChangeEditFeatureLocationButton);
  SetObjectExtra (dlg->extend5, dlg, NULL);

  dlg->clear_5_type = PopupList (g, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->clear_5_type, dlg, NULL);
  PopupItem (dlg->clear_5_type, "All");
  PopupItem (dlg->clear_5_type, "If 5' end not at end of sequence");
  PopupItem (dlg->clear_5_type, "If good start codon");
  SetValue (dlg->clear_5_type, 1);

  dlg->set_3_partial_grp = HiddenGroup (g, 2, 0, NULL);
  dlg->set_3_type = PopupList (dlg->set_3_partial_grp, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->set_3_type, dlg, NULL);
  PopupItem (dlg->set_3_type, "All");
  PopupItem (dlg->set_3_type, "If 3' end at end of sequence");
  PopupItem (dlg->set_3_type, "If bad stop codon");
  SetValue (dlg->set_3_type, 1);
  dlg->extend3 = CheckBox (dlg->set_3_partial_grp, "Extend to 3' end if partial is set", ChangeEditFeatureLocationButton);
  SetObjectExtra (dlg->extend3, dlg, NULL);

  dlg->clear_3_type = PopupList (g, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->clear_3_type, dlg, NULL);
  PopupItem (dlg->clear_3_type, "All");
  PopupItem (dlg->clear_3_type, "If 3' end not at end of sequence");
  PopupItem (dlg->clear_3_type, "If good stop codon");
  SetValue (dlg->clear_3_type, 1);

  dlg->set_both_partial_grp = HiddenGroup (g, 2, 0, NULL);
  dlg->set_both_type = PopupList (dlg->set_both_partial_grp, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->set_both_type, dlg, NULL);
  PopupItem (dlg->set_both_type, "All");
  PopupItem (dlg->set_both_type, "If both ends at end of sequence");
  SetValue (dlg->set_both_type, 1);
  dlg->extendboth = CheckBox (dlg->set_both_partial_grp, "Extend to ends of sequence if partials are set", ChangeEditFeatureLocationButton);
  SetObjectExtra (dlg->extendboth, dlg, NULL);

  dlg->clear_both_type = PopupList (g, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->clear_both_type, dlg, NULL);
  PopupItem (dlg->clear_both_type, "All");
  PopupItem (dlg->clear_both_type, "If both ends not at end of sequence");
  SetValue (dlg->clear_both_type, 1);

  /* convert location */
  dlg->convert_loc = PopupList (g, TRUE, ChangeEditFeatureLocationPopup);
  SetObjectExtra (dlg->convert_loc, dlg, NULL);
  PopupItem (dlg->convert_loc, "Join");
  PopupItem (dlg->convert_loc, "Order");
  PopupItem (dlg->convert_loc, "Single Interval");
  SetValue (dlg->convert_loc, 1);

  AlignObjects (ALIGN_MIDDLE, (HANDLE) dlg->strand_edit_grp,
                              (HANDLE) dlg->set_5_partial_grp,
                              (HANDLE) dlg->clear_5_type,
                              (HANDLE) dlg->set_3_partial_grp,
                              (HANDLE) dlg->clear_3_type,
                              (HANDLE) dlg->set_both_partial_grp,
                              (HANDLE) dlg->clear_both_type,
                              (HANDLE) dlg->convert_loc,
                              (HANDLE) NULL);

  /* retranslate option */
  dlg->retranslate_cds = CheckBox (p, "Retranslate coding region if location changes", ChangeEditFeatureLocationButton);
  SetObjectExtra (dlg->retranslate_cds, dlg, NULL);

  dlg->constraint_dlg = ConstraintSetDialog (p, change_notify, change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_type_dlg,
                              (HANDLE) dlg->feature_edit_type,
                              (HANDLE) g,
                              (HANDLE) dlg->retranslate_cds,
                              (HANDLE) dlg->constraint_dlg,
                              NULL);

  ChangeFeatureLocationEditType (dlg->feature_edit_type);

  return (DialoG) p;
}


typedef struct parsesrcgeneraliddlg {
 DIALOG_MESSAGE_BLOCK
 PopuP wanted;
 TexT  db_wanted;
 Nlm_ChangeNotifyProc change_notify;
 Pointer              change_userdata;
} ParseSrcGeneralIdDlgData, PNTR ParseSrcGeneralIdDlgPtr;


static void ChangeParseSrcGeneralIdDialogPopup (PopuP p)
{
  ParseSrcGeneralIdDlgPtr dlg;
  Int2           val;

  dlg = (ParseSrcGeneralIdDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  val = GetValue (dlg->wanted);
  if (val == 3) {
    Show (dlg->db_wanted);
  } else {
    Hide (dlg->db_wanted);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ParseSrcGeneralIdToDialog (DialoG d, Pointer data)
{
  ParseSrcGeneralIdDlgPtr dlg;
  ValNodePtr  id;

  dlg = (ParseSrcGeneralIdDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  Hide (dlg->db_wanted);
  if ((id = (ValNodePtr) data) == NULL) {
    SetValue (dlg->wanted, 1);
    SetTitle (dlg->db_wanted, "");
  } else {
    switch (id->choice) {
      case ParseSrcGeneralId_whole_text:
        SetValue (dlg->wanted, 1);
        SetTitle (dlg->db_wanted, "");
        break;
      case ParseSrcGeneralId_db:
        SetValue (dlg->wanted, 2);
        SetTitle (dlg->db_wanted, "");
        break;
      case ParseSrcGeneralId_tag:
        SetValue (dlg->wanted, 3);
        if (StringHasNoText (id->data.ptrvalue)) {
          SetTitle (dlg->db_wanted, "");
        } else {
          SetTitle (dlg->db_wanted, id->data.ptrvalue);
        }
        Show (dlg->db_wanted);
        break;
      default:
        SetValue (dlg->wanted, 1);
        SetTitle (dlg->db_wanted, "");
        break;
    }
  }
}


static Pointer ParseSrcGeneralIdFromDialog (DialoG d)
{
  ParseSrcGeneralIdDlgPtr dlg;
  ValNodePtr  id;
  Int2        val;

  dlg = (ParseSrcGeneralIdDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  id = ValNodeNew (NULL);
  val = GetValue (dlg->wanted);
  switch (val) {
    case 1:
      id->choice = ParseSrcGeneralId_whole_text;
      break;
    case 2:
      id->choice = ParseSrcGeneralId_db;
      break;
    case 3:
      id->choice = ParseSrcGeneralId_tag;
      if (!TextHasNoText (dlg->db_wanted)) {
        id->data.ptrvalue = SaveStringFromText (dlg->db_wanted);
      }
      break;
    default:
      id->choice = ParseSrcGeneralId_whole_text;
      break;
  }
  return id;
}


static void ChangeParseGeneralIdDlgText (TexT t)
{
  ParseSrcGeneralIdDlgPtr dlg;

  dlg = (ParseSrcGeneralIdDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG ParseGeneralIdDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ParseSrcGeneralIdDlgPtr dlg;
  GrouP          p;

  p = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dlg = (ParseSrcGeneralIdDlgPtr) MemNew (sizeof (ParseSrcGeneralIdDlgData));
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = ParseSrcGeneralIdFromDialog;
  dlg->todialog = ParseSrcGeneralIdToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->wanted = PopupList (p, TRUE, ChangeParseSrcGeneralIdDialogPopup);
  SetObjectExtra (dlg->wanted, dlg, NULL);
  PopupItem (dlg->wanted, "Entire General ID");
  PopupItem (dlg->wanted, "General ID db");
  PopupItem (dlg->wanted, "General ID tag");
  SetValue (dlg->wanted, 1);
  dlg->db_wanted = DialogText (p, "", 5, ChangeParseGeneralIdDlgText);
  SetObjectExtra (dlg->db_wanted, dlg, NULL);
  Hide (dlg->db_wanted);

  return (DialoG) p;
}


typedef struct parsesrcdlg {
 DIALOG_MESSAGE_BLOCK
 PopuP main_choice;
 GrouP src_grp;
 DialoG src_qual_choice;
 GrouP feat_or_desc;
 TexT  structured_comment_field;
 DialoG structured_comment_field_list;
 DialoG general_id_dlg;

 Nlm_ChangeNotifyProc change_notify;
 Pointer              change_userdata;
} ParseSrcDlgData, PNTR ParseSrcDlgPtr;


static void ChangeParseSrcDialogPopup (PopuP p)
{
  ParseSrcDlgPtr dlg;
  Int2           val;

  dlg = (ParseSrcDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  val = GetValue (dlg->main_choice);
  switch (val) {
    case 5:
      Show (dlg->src_grp);
      SafeHide (dlg->structured_comment_field);
      SafeHide (dlg->structured_comment_field_list);
      SafeHide (dlg->general_id_dlg);
      break;
    case 8:
      SafeShow (dlg->structured_comment_field);
      SafeShow (dlg->structured_comment_field_list);
      SafeHide (dlg->general_id_dlg);
      Hide (dlg->src_grp);
      break;
    case 11:
      SafeHide (dlg->structured_comment_field);
      SafeHide (dlg->structured_comment_field_list);
      SafeShow (dlg->general_id_dlg);
      Hide (dlg->src_grp);
      break;
    default:
      SafeHide (dlg->structured_comment_field);
      SafeHide (dlg->structured_comment_field_list);
      SafeHide (dlg->general_id_dlg);
      Hide (dlg->src_grp);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeParseSrcDialogText (TexT t)
{
  ParseSrcDlgPtr dlg;

  dlg = (ParseSrcDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) return;

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ParseSrcToDialog (DialoG d, Pointer data)
{
  ParseSrcDlgPtr dlg;
  ValNodePtr src;
  ParseSrcOrgPtr org;
  ValNode        vn;

  dlg = (ParseSrcDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  src = (ValNodePtr) data;
  if (src == NULL) {
    SetValue (dlg->main_choice, 4);
    vn.choice = 0;
    vn.data.ptrvalue = kTaxnameAfterBinomialString;
    vn.next = NULL;
    PointerToDialog (dlg->src_qual_choice, &vn);
    SetValue (dlg->feat_or_desc, 1);
  } else {
    switch (src->choice) {
      case ParseSrc_defline:
        SetValue (dlg->main_choice, 1);
        break;
      case ParseSrc_flatfile:
        SetValue (dlg->main_choice, 2);
        break;
      case ParseSrc_local_id:
        SetValue (dlg->main_choice, 3);
        break;
      case ParseSrc_org:
        org = (ParseSrcOrgPtr) src->data.ptrvalue;
        if (org == NULL) {
          SetValue (dlg->main_choice, 5);
          SetValue (dlg->feat_or_desc, 1);
          PointerToDialog (dlg->src_qual_choice, NULL);
        } else {
          switch (org->type) {
            case Object_type_constraint_any:
              SetValue (dlg->feat_or_desc, 1);
              break;
            case Object_type_constraint_feature:
              SetValue (dlg->feat_or_desc, 2);
              break;
            case Object_type_constraint_descriptor:
              SetValue (dlg->feat_or_desc, 3);
              break;
            default:
              SetValue (dlg->feat_or_desc, 1);
              break;
          }
          if (org->field == NULL) {
            PointerToDialog (dlg->src_qual_choice, NULL);
            SetValue (dlg->main_choice, 5);
          } else if (org->field->choice == ParseSrcOrgChoice_taxname_after_binomial) {
            SetValue (dlg->main_choice, 9);
          } else if (org->field->choice == ParseSrcOrgChoice_source_qual) {
            vn.choice = 0;
            vn.data.ptrvalue = GetSourceQualName (org->field->data.intvalue);
            vn.next = NULL;
            PointerToDialog (dlg->src_qual_choice, &vn);
            if (org->field->data.intvalue == Source_qual_taxname) {
              SetValue (dlg->main_choice, 4);
            } else {
              SetValue (dlg->main_choice, 5);
            }
          } else {
            PointerToDialog (dlg->src_qual_choice, NULL);
            SetValue (dlg->main_choice, 5);
          }
        }
        break;
      case ParseSrc_comment:
        SetValue (dlg->main_choice, 6);
        break;
      case ParseSrc_bankit_comment:
        SetValue (dlg->main_choice, 7);
        break;
      case ParseSrc_structured_comment:
        SetValue (dlg->main_choice, 8);
        if (dlg->structured_comment_field != NULL) {
          SetTitle (dlg->structured_comment_field, src->data.ptrvalue);
        } else {
          PointerToDialog (dlg->structured_comment_field_list, src->data.ptrvalue);
        }
        break;
      case ParseSrc_file_id:
        SetValue (dlg->main_choice, 10);
        break;
      case ParseSrc_general_id:
        SetValue (dlg->main_choice, 11);
        PointerToDialog (dlg->general_id_dlg, src->data.ptrvalue);
        break;
      default:
        SetValue (dlg->main_choice, 0);
        break;
    }
  }

  ChangeParseSrcDialogPopup (dlg->main_choice);
}


static Pointer DialogToParseSrc (DialoG d)
{
  ParseSrcDlgPtr dlg;
  ParseSrcOrgPtr org;
  ValNodePtr     src = NULL, vnp;
  Int2           val;

  dlg = (ParseSrcDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->main_choice);
  switch (val) {
    case 1:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_defline;
      break;
    case 2:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_flatfile;
      break;
    case 3:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_local_id;
      break;
    case 4:
      org = ParseSrcOrgNew ();
      switch (GetValue (dlg->feat_or_desc)) {
        case 1:
          org->type = Object_type_constraint_any;
          break;
        case 2:
          org->type = Object_type_constraint_feature;
          break;
        case 3:
          org->type = Object_type_constraint_descriptor;
          break;
        default:
          org->type = Object_type_constraint_any;
          break;
      }
      org->field = ValNodeNew (NULL);
      org->field->choice = ParseSrcOrgChoice_source_qual;
      org->field->data.intvalue = Source_qual_taxname;
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_org;
      src->data.ptrvalue = org;
      break;
    case 5:
      org = ParseSrcOrgNew ();
      switch (GetValue (dlg->feat_or_desc)) {
        case 1:
          org->type = Object_type_constraint_any;
          break;
        case 2:
          org->type = Object_type_constraint_feature;
          break;
        case 3:
          org->type = Object_type_constraint_descriptor;
          break;
        default:
          org->type = Object_type_constraint_any;
          break;
      }
      vnp = (ValNodePtr) DialogToPointer (dlg->src_qual_choice);
      if (vnp != NULL) {
        org->field = ValNodeNew (NULL);
        if (StringCmp (vnp->data.ptrvalue, kTaxnameAfterBinomialString) == 0) {
          org->field->choice = ParseSrcOrgChoice_taxname_after_binomial;
        } else {
          org->field->choice = ParseSrcOrgChoice_source_qual;
          org->field->data.intvalue = GetSourceQualTypeByName (vnp->data.ptrvalue);
        }
        vnp = ValNodeFree (vnp);
      }
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_org;
      src->data.ptrvalue = org;
      break;
    case 6:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_comment;
      break;
    case 7:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_bankit_comment;
      break;
    case 8:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_structured_comment;
      if (dlg->structured_comment_field != NULL) {
        src->data.ptrvalue = SaveStringFromText (dlg->structured_comment_field);
      } else {
        src->data.ptrvalue = DialogToPointer (dlg->structured_comment_field_list);
      }
      break;
    case 9:
      org = ParseSrcOrgNew ();
      switch (GetValue (dlg->feat_or_desc)) {
        case 1:
          org->type = Object_type_constraint_any;
          break;
        case 2:
          org->type = Object_type_constraint_feature;
          break;
        case 3:
          org->type = Object_type_constraint_descriptor;
          break;
        default:
          org->type = Object_type_constraint_any;
          break;
      }
      org->field = ValNodeNew (NULL);
      org->field->choice = ParseSrcOrgChoice_taxname_after_binomial;
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_org;
      src->data.ptrvalue = org;
      break;
    case 10:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_file_id;
      break;
    case 11:
      src = ValNodeNew (NULL);
      src->choice = ParseSrc_general_id;
      src->data.ptrvalue = DialogToPointer (dlg->general_id_dlg);
      break;
  } 
  return src;
}


static ValNodePtr TestParseSrcDialog (DialoG d)
{
  ParseSrcDlgPtr dlg;
  Int2 val;
  ValNodePtr err_list = NULL;

  dlg = (ParseSrcDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->main_choice);
  if (val < 1 || val > 11) {
    ValNodeAddPointer (&err_list, 0, "main choice");
  } else if (val == 5) {
    ValNodeLink (&err_list, TestDialog (dlg->src_qual_choice));
  } else if (val == 8 && 
             ((dlg->structured_comment_field != NULL && TextHasNoText (dlg->structured_comment_field)))) { 
    ValNodeAddPointer (&err_list, 0, "structured comment field");
  }
  return err_list;
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


static DialoG ParseSrcDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata, SeqEntryPtr sep)
{
  ParseSrcDlgPtr dlg;
  GrouP          p, g;
  ValNodePtr     qual_list;
  ValNodePtr     structured_comment_field_list = NULL;
  ValNode        vn;

  p = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dlg = (ParseSrcDlgPtr) MemNew (sizeof (ParseSrcDlgData));
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = DialogToParseSrc;
  dlg->todialog = ParseSrcToDialog;
  dlg->testdialog = TestParseSrcDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  VisitDescriptorsInSep (sep, &structured_comment_field_list, AddCommentFieldName);  

  dlg->main_choice = PopupList (p, TRUE, ChangeParseSrcDialogPopup);
  SetObjectExtra (dlg->main_choice, dlg, NULL);
  PopupItem (dlg->main_choice, "Definition Line");
  PopupItem (dlg->main_choice, "GenBank FlatFile");
  PopupItem (dlg->main_choice, "Local ID");
  PopupItem (dlg->main_choice, "Organism Name");
  PopupItem (dlg->main_choice, "Source Qualifier");
  PopupItem (dlg->main_choice, "Comment");
  PopupItem (dlg->main_choice, "BankIT Comment");
  PopupItem (dlg->main_choice, "Structured Comment");
  PopupItem (dlg->main_choice, "Taxname after Binomial");
  PopupItem (dlg->main_choice, "File ID");
  PopupItem (dlg->main_choice, "General ID");
  SetValue (dlg->main_choice, 4);

  g = HiddenGroup (p, 0, 0, NULL);
  dlg->src_grp = HiddenGroup (g, -1, 0, NULL);
  qual_list = ValNodeNew (NULL);
  qual_list->choice = 0;
  qual_list->data.ptrvalue = StringSave (kTaxnameAfterBinomialString);
  qual_list->next = GetSourceQualList (FALSE);  
  dlg->src_qual_choice = ValNodeSelectionDialog (dlg->src_grp, qual_list, SHORT_SELECTION_LIST,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "source qual",
                                               change_notify, change_userdata, FALSE);
  vn.choice = 0;
  vn.data.ptrvalue = "Taxname after binomial";
  vn.next = NULL;
  PointerToDialog (dlg->src_qual_choice, &vn);
  dlg->feat_or_desc = HiddenGroup (dlg->src_grp, 3, 0, NULL);
  RadioButton (dlg->feat_or_desc, "Descriptors and Features");
  RadioButton (dlg->feat_or_desc, "Descriptors Only");
  RadioButton (dlg->feat_or_desc, "Features Only");
  SetValue (dlg->feat_or_desc, 1);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->src_qual_choice, (HANDLE) dlg->feat_or_desc, NULL);
  Hide (dlg->src_grp);

  if (structured_comment_field_list == NULL) {
    dlg->structured_comment_field = DialogText (g, "", 20, ChangeParseSrcDialogText);
    SetObjectExtra (dlg->structured_comment_field, dlg, NULL);
    Hide (dlg->structured_comment_field);
  } else {
    dlg->structured_comment_field_list = StringComboDialog (g, structured_comment_field_list, SHORT_SELECTION_LIST, 20, change_notify, change_userdata);
    Hide (dlg->structured_comment_field_list);
  }

  dlg->general_id_dlg = ParseGeneralIdDialog(g, change_notify, change_userdata);
  Hide (dlg->general_id_dlg);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->src_grp, (HANDLE) dlg->structured_comment_field, (HANDLE) dlg->general_id_dlg, NULL);

  return (DialoG) p;
}


typedef struct parsedstdlg {
 DIALOG_MESSAGE_BLOCK
 PopuP main_choice;
 GrouP src_grp;
 PopuP biosrc_list;
 DialoG src_qual_choice;
 GrouP feat_or_desc;
 PopuP  gene_field;
 DialoG rna_field;
 PopuP protein_field;
 GrouP  feature_field_grp;
 DialoG feature_type_choice;
 DialoG feature_field_choice;
 DialoG feature_note;
 TexT  dbxref_db;
 Nlm_ChangeNotifyProc change_notify;
 Pointer              change_userdata;
} ParseDstDlgData, PNTR ParseDstDlgPtr;

typedef enum {
  eParseDstDialog_defline = 1,
  eParseDstDialog_biosrc,
  eParseDstDialog_srcqual,
  eParseDstDialog_genefield,
  eParseDstDialog_rnafield,
  eParseDstDialog_cdscomment,
  eParseDstDialog_proteinfield,
  eParseDstDialog_featqual,
  eParseDstDialog_featnote,
  eParseDstDialog_commentdesc,
  eParseDstDialog_dbxref 
} EParseDstDialog;

#define eParseDstDialog_last eParseDstDialog_dbxref

static void ChangeParseDstDialogPopup (PopuP p) 
{
  ParseDstDlgPtr dlg;
  Int2 val;

  dlg = (ParseDstDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  val = GetValue (dlg->main_choice);
  switch (val) {
    case eParseDstDialog_biosrc:
      Show (dlg->src_grp);
      Show (dlg->biosrc_list);
      Hide (dlg->src_qual_choice);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
    case eParseDstDialog_srcqual:
      Show (dlg->src_grp);
      Hide (dlg->biosrc_list);
      Show (dlg->src_qual_choice);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
    case eParseDstDialog_genefield:
      Show (dlg->gene_field);
      Hide (dlg->src_grp);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
    case eParseDstDialog_rnafield:
      Show (dlg->rna_field);
      Hide (dlg->src_grp);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
    case eParseDstDialog_proteinfield:
      Show (dlg->protein_field);
      Hide (dlg->rna_field);
      Hide (dlg->src_grp);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->gene_field);
      Hide (dlg->feature_note);
      break;
    case eParseDstDialog_featqual:
      Hide (dlg->src_grp);
      Show (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
    case eParseDstDialog_featnote:
      Show (dlg->feature_note);
      Hide (dlg->src_grp);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      break;
    case eParseDstDialog_dbxref:
      Hide (dlg->src_grp);
      Hide (dlg->feature_field_grp);
      Show (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
    default:
      Hide (dlg->src_grp);
      Hide (dlg->feature_field_grp);
      Hide (dlg->dbxref_db);
      Hide (dlg->rna_field);
      Hide (dlg->gene_field);
      Hide (dlg->protein_field);
      Hide (dlg->feature_note);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ChangeParseDstDialogText (TexT t) 
{
  ParseDstDlgPtr dlg;

  dlg = (ParseDstDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) return;

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ParseDstToDialog (DialoG d, Pointer data)
{
  ParseDstDlgPtr dlg;
  ValNodePtr     dst;
  ParseDstOrgPtr org;
  FeatureFieldLegalPtr ffp;
  ValNode        vn;

  dlg = (ParseDstDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  dst = (ValNodePtr) data;
  if (dst == NULL) {
    SetValue (dlg->main_choice, eParseDstDialog_defline);
  } else {
    switch (dst->choice) {
      case ParseDest_defline:
        SetValue (dlg->main_choice, eParseDstDialog_defline );
        break;
      case ParseDest_org:
        org = (ParseDstOrgPtr) dst->data.ptrvalue;
        if (org == NULL) {
          SetValue (dlg->main_choice, eParseDstDialog_srcqual);
          PointerToDialog (dlg->src_qual_choice, NULL);
          SetValue (dlg->feat_or_desc, 1);
        } else {
          switch (org->type) {
            case Object_type_constraint_any:
              SetValue (dlg->feat_or_desc, 1);
              break;
            case Object_type_constraint_feature:
              SetValue (dlg->feat_or_desc, 3);
              break;
            case Object_type_constraint_descriptor:
              SetValue (dlg->feat_or_desc, 2);
              break;
            default:
              SetValue (dlg->feat_or_desc, 1);
              break;
          }
          if (org->field == NULL) {
            PointerToDialog (dlg->src_qual_choice, NULL);
            SetValue (dlg->main_choice, eParseDstDialog_srcqual);
          } else {
            switch (org->field->data.intvalue) {
              case Source_qual_taxname:
                SetValue (dlg->biosrc_list, 1);
                SetValue (dlg->main_choice, eParseDstDialog_biosrc);
                break;
              case Source_qual_lineage:
                SetValue (dlg->biosrc_list, 2);
                SetValue (dlg->main_choice, eParseDstDialog_biosrc);
                break;
              case Source_qual_division:
                SetValue (dlg->biosrc_list, 3);
                SetValue (dlg->main_choice, eParseDstDialog_biosrc);
                break;
              default:
                vn.choice = 0;
                vn.data.ptrvalue = GetSourceQualName (org->field->data.intvalue);
                vn.next = NULL;
                PointerToDialog (dlg->src_qual_choice, &vn);
                SetValue (dlg->main_choice, eParseDstDialog_biosrc);
                break;
            }
          }
        }
        break;
      case ParseDest_featqual:
        /* TODO - look for RNA fields, gene fields, cds comment, feature note? */
        SetValue (dlg->main_choice, eParseDstDialog_featqual);
        ffp = (FeatureFieldLegalPtr) dst->data.ptrvalue;
        if (ffp == NULL) {
          PointerToDialog (dlg->feature_type_choice, NULL);
          PointerToDialog (dlg->feature_field_choice, NULL);
        } else {
          vn.choice = (Uint1)ffp->type;
          vn.data.ptrvalue = NULL;
          vn.next = NULL;
          PointerToDialog (dlg->feature_type_choice, &vn);
          vn.choice = (Uint1)ffp->field;
          PointerToDialog (dlg->feature_field_choice, &vn);
        }
        break;
      case ParseDest_comment_descriptor:
        SetValue (dlg->main_choice, eParseDstDialog_commentdesc);
        break;
      case ParseDest_dbxref:
        SetValue (dlg->main_choice, eParseDstDialog_dbxref);
        SetTitle (dlg->dbxref_db, dst->data.ptrvalue);
        break;
    }
  }

  ChangeParseDstDialogPopup (dlg->main_choice);
}


static Pointer DialogToParseDst (DialoG d)
{
  ParseDstDlgPtr dlg;
  ParseDstOrgPtr org;
  FeatureFieldLegalPtr ffp;
  FeatureFieldPtr ff;
  ValNodePtr     dst = NULL, vnp;
  Int2           val, subval;
  RnaQualPtr     rq;

  dlg = (ParseDstDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->main_choice);
  switch (val) {
    case eParseDstDialog_defline:
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_defline;
      break;
    case eParseDstDialog_biosrc:
      org = ParseDstOrgNew ();
      switch (GetValue (dlg->feat_or_desc)) {
        case 1:
          org->type = Object_type_constraint_any;
          break;
        case 3:
          org->type = Object_type_constraint_feature;
          break;
        case 2:
          org->type = Object_type_constraint_descriptor;
          break;
        default:
          org->type = Object_type_constraint_any;
          break;
      }
      subval = GetValue (dlg->biosrc_list);
      switch (subval) {
        case 1:
          org->field = ValNodeNew (NULL);
          org->field->choice = SourceQualChoice_textqual;
          org->field->data.intvalue = Source_qual_taxname;
          break;
        case 2:
          org->field = ValNodeNew (NULL);
          org->field->choice = SourceQualChoice_textqual;
          org->field->data.intvalue = Source_qual_lineage;
          break;
        case 3:
          org->field = ValNodeNew (NULL);
          org->field->choice = SourceQualChoice_textqual;
          org->field->data.intvalue = Source_qual_division;
          break;
      }
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_org;
      dst->data.ptrvalue = org;
      break;
    case eParseDstDialog_srcqual:
      org = ParseDstOrgNew ();
      switch (GetValue (dlg->feat_or_desc)) {
        case 1:
          org->type = Object_type_constraint_any;
          break;
        case 3:
          org->type = Object_type_constraint_feature;
          break;
        case 2:
          org->type = Object_type_constraint_descriptor;
          break;
        default:
          org->type = Object_type_constraint_any;
          break;
      }
      vnp = (ValNodePtr) DialogToPointer (dlg->src_qual_choice);
      if (vnp != NULL) {
        org->field = ValNodeNew (NULL);
        org->field->choice = SourceQualChoice_textqual;
        org->field->data.intvalue = GetSourceQualTypeByName (vnp->data.ptrvalue);
        vnp = ValNodeFree (vnp);
      }      
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_org;
      dst->data.ptrvalue = org;
      break;
    case eParseDstDialog_genefield:
      subval = GetValue (dlg->gene_field);
      switch (subval) {
        case 1: /* locus */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_gene;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 2: /* description */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_gene_description;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 3: /* comment */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_note;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 4: /* allele */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_allele;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 5: /* maploc */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_map;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 6: /* locus tag */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_locus_tag;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 7: /* synonym */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_synonym;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 8: /* old locus tag */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_gene;
          ffp->field = Feat_qual_legal_old_locus_tag;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
      }
      break;
    case eParseDstDialog_rnafield:
      /* get rna qual, translate to feature qual */
      rq = DialogToPointer (dlg->rna_field);
      ff = FeatureFieldFromRnaQual (rq);
      rq = RnaQualFree (rq);
      if (ff != NULL) {
        ffp = FeatureFieldLegalNew ();
        ffp->type = ff->type;
        if (ff->field != NULL) {
          ffp->field = ff->field->data.intvalue;
        }
        dst = ValNodeNew (NULL);
        dst->choice = ParseDest_featqual;
        dst->data.ptrvalue = ffp;
      }
      ff = FeatureFieldFree (ff);
      break;
    case eParseDstDialog_cdscomment:
      ffp = FeatureFieldLegalNew ();
      ffp->type = Macro_feature_type_cds;
      ffp->field = Feat_qual_legal_note;
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_featqual;
      dst->data.ptrvalue = ffp;
      break;
    case eParseDstDialog_proteinfield:
      subval = GetValue (dlg->protein_field);
      switch (subval) {
        case 1: /* name */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_cds;
          ffp->field = Feat_qual_legal_product;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 2: /* description */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_prot;
          ffp->field = Feat_qual_legal_description;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 3: /* e.c. number */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_prot;
          ffp->field = Feat_qual_legal_ec_number;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 4: /* activity */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_prot;
          ffp->field = Feat_qual_legal_activity;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 5: /* comment */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_prot;
          ffp->field = Feat_qual_legal_note;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 6: /* mat-peptide name */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_mat_peptide_aa;
          ffp->field = Feat_qual_legal_product;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 7: /* mat-peptide description */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_mat_peptide_aa;
          ffp->field = Feat_qual_legal_description;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
        case 8: /* mat-peptide comment */
          ffp = FeatureFieldLegalNew ();
          ffp->type = Macro_feature_type_mat_peptide_aa;
          ffp->field = Feat_qual_legal_note;
          dst = ValNodeNew (NULL);
          dst->choice = ParseDest_featqual;
          dst->data.ptrvalue = ffp;
          break;
      }
      break;
    case eParseDstDialog_featqual:
      ffp = FeatureFieldLegalNew ();
      vnp = (ValNodePtr) DialogToPointer (dlg->feature_type_choice);
      if (vnp != NULL) {
        ffp->type = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      vnp = (ValNodePtr) DialogToPointer (dlg->feature_field_choice);
      if (vnp != NULL) {
        ffp->field = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_featqual;
      dst->data.ptrvalue = ffp;
      break;
    case eParseDstDialog_featnote:
      ffp = FeatureFieldLegalNew ();
      vnp = (ValNodePtr) DialogToPointer (dlg->feature_note);
      if (vnp != NULL) {
        ffp->type = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      ffp->field = Feat_qual_legal_note;
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_featqual;
      dst->data.ptrvalue = ffp;
      break;
    case eParseDstDialog_commentdesc:
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_comment_descriptor;
      break;
    case eParseDstDialog_dbxref:
      dst = ValNodeNew (NULL);
      dst->choice = ParseDest_dbxref;
      dst->data.ptrvalue = SaveStringFromText (dlg->dbxref_db);
      break;
  }
  return dst;
}


static ValNodePtr TestParseDstDialog (DialoG d)
{
  ParseDstDlgPtr dlg;
  ValNodePtr     err_list = NULL;
  Int2           val;

  dlg = (ParseDstDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->main_choice);
  if (val < 1 || val > eParseDstDialog_last) {
    ValNodeAddPointer (&err_list, 0, "main choice");
  } else if (val == eParseDstDialog_srcqual) {
    ValNodeLink (&err_list, TestDialog (dlg->src_qual_choice));
  } else if (val == eParseDstDialog_rnafield) {
    ValNodeLink (&err_list, TestDialog (dlg->rna_field));
  } else if (val == eParseDstDialog_featqual) {
    ValNodeLink (&err_list, TestDialog (dlg->feature_type_choice));
    ValNodeLink (&err_list, TestDialog (dlg->feature_field_choice));
  } else if (val == eParseDstDialog_featnote) {
    ValNodeLink (&err_list, TestDialog (dlg->feature_note));
  } else if (val == eParseDstDialog_dbxref && TextHasNoText (dlg->dbxref_db)) {
    ValNodeAddPointer (&err_list, 0, "dbxref db");
  }
  return err_list;
}


NLM_EXTERN DialoG ParseDstDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  ParseDstDlgPtr dlg;
  GrouP          p, g, g2;
  ValNodePtr     qual_list;

  p = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  dlg = (ParseDstDlgPtr) MemNew (sizeof (ParseDstDlgData));
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ParseDstToDialog;
  dlg->fromdialog = DialogToParseDst;
  dlg->testdialog = TestParseDstDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->main_choice = PopupList (p, TRUE, ChangeParseDstDialogPopup);
  SetObjectExtra (dlg->main_choice, dlg, NULL);
  PopupItem (dlg->main_choice, "Definition Line");
  PopupItem (dlg->main_choice, "Biosource");
  PopupItem (dlg->main_choice, "Source Qualifier");
  PopupItem (dlg->main_choice, "Gene Field");
  PopupItem (dlg->main_choice, "RNA Field");
  PopupItem (dlg->main_choice, "CDS Comment");
  PopupItem (dlg->main_choice, "Protein Field");
  PopupItem (dlg->main_choice, "Feature Qualifier");
  PopupItem (dlg->main_choice, "Feature Note");
  PopupItem (dlg->main_choice, "Comment Descriptor");
  PopupItem (dlg->main_choice, "Dbxref");
  SetValue (dlg->main_choice, 1);

  g = HiddenGroup (p, 0, 0, NULL);
  dlg->src_grp = HiddenGroup (g, -1, 0, NULL);
  g2 = HiddenGroup (dlg->src_grp, 0, 0, NULL);
  qual_list = GetSourceQualList (FALSE);  
  dlg->src_qual_choice = ValNodeSelectionDialog (g2, qual_list, SHORT_SELECTION_LIST,
                                               ValNodeStringName,
                                               ValNodeSimpleDataFree,
                                               ValNodeStringCopy,
                                               ValNodeStringMatch,
                                               "source qual",
                                               change_notify, change_userdata, FALSE);
  dlg->biosrc_list = PopupList (g2, TRUE, ChangeParseDstDialogPopup);
  SetObjectExtra (dlg->biosrc_list, dlg, NULL);
  PopupItem (dlg->biosrc_list, "Organism Name");
  PopupItem (dlg->biosrc_list, "Lineage");
  PopupItem (dlg->biosrc_list, "Division");
  SetValue (dlg->biosrc_list, 1);

  dlg->feat_or_desc = HiddenGroup (dlg->src_grp, 3, 0, NULL);
  RadioButton (dlg->feat_or_desc, "Descriptors and Features");
  RadioButton (dlg->feat_or_desc, "Descriptors Only");
  RadioButton (dlg->feat_or_desc, "Features Only");
  SetValue (dlg->feat_or_desc, 1);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->src_qual_choice, (HANDLE) dlg->feat_or_desc, NULL);
  Hide (dlg->src_grp);

  dlg->gene_field = PopupList (g, TRUE, ChangeParseDstDialogPopup);
  SetObjectExtra (dlg->gene_field, dlg, NULL);
  PopupItem (dlg->gene_field, "locus");
  PopupItem (dlg->gene_field, "description");
  PopupItem (dlg->gene_field, "comment");
  PopupItem (dlg->gene_field, "allele");
  PopupItem (dlg->gene_field, "maploc");
  PopupItem (dlg->gene_field, "locus_tag");
  PopupItem (dlg->gene_field, "synonym");
  PopupItem (dlg->gene_field, "old_locus_tag");
  SetValue (dlg->gene_field, 1);
  Hide (dlg->gene_field);

  dlg->rna_field = RnaQualDialog (g, "RNA Type", change_notify, change_userdata);
  Hide (dlg->rna_field);

  dlg->protein_field = PopupList (g, TRUE, ChangeParseDstDialogPopup);
  SetObjectExtra (dlg->protein_field, dlg, NULL);
  PopupItem (dlg->protein_field, "name");
  PopupItem (dlg->protein_field, "description");
  PopupItem (dlg->protein_field, "E.C. number");
  PopupItem (dlg->protein_field, "activity");
  PopupItem (dlg->protein_field, "comment");
  PopupItem (dlg->protein_field, "mat_peptide name");
  PopupItem (dlg->protein_field, "mat_peptide description");
  PopupItem (dlg->protein_field, "mat_peptide comment");
  SetValue (dlg->protein_field, 1);
  Hide (dlg->protein_field);

  dlg->feature_field_grp = HiddenGroup (g, 2, 0, NULL);
  dlg->feature_type_choice = FeatureTypeDialog (dlg->feature_field_grp, change_notify, change_userdata);
  dlg->feature_field_choice = LegalFeatQualChoiceDialog (dlg->feature_field_grp, change_notify, change_userdata);
  Hide (dlg->feature_field_grp);

  dlg->feature_note = FeatureTypeDialog (g, change_notify, change_userdata);
  Hide (dlg->feature_note);

  dlg->dbxref_db = DialogText (g, "", 20, ChangeParseDstDialogText);
  SetObjectExtra (dlg->dbxref_db, dlg, NULL);
  Hide (dlg->dbxref_db);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->src_grp, (HANDLE) dlg->feature_field_grp, (HANDLE) dlg->dbxref_db, NULL);

  return (DialoG) p;
}


typedef struct capchangedlg {
 DIALOG_MESSAGE_BLOCK
 GrouP  cap_change;
 Nlm_ChangeNotifyProc change_notify;
 Pointer              change_userdata;
} CapChangeDlgData, PNTR CapChangeDlgPtr;


NLM_EXTERN void SetCapChangeDialogValue (DialoG d, Int4 cap_change)
{
  CapChangeDlgPtr dlg;

  dlg = (CapChangeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  switch (cap_change) {
    case Cap_change_none:
      SetValue (dlg->cap_change, 1);
      break;
    case Cap_change_tolower:
      SetValue (dlg->cap_change, 2);
      break;
    case Cap_change_toupper:
      SetValue (dlg->cap_change, 3);
      break;
    case Cap_change_firstcap:
      SetValue (dlg->cap_change, 4);
      break;
    case Cap_change_firstcaprestnochange:
      SetValue (dlg->cap_change, 5);
      break;
    default:
      SetValue (dlg->cap_change, 1);
      break;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


NLM_EXTERN Int2 GetCapChangeDialogValue (DialoG d)
{
  CapChangeDlgPtr dlg;
  Int2 cap_change = Cap_change_none;

  dlg = (CapChangeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return cap_change;

  switch (GetValue (dlg->cap_change)) {
    case 1:
      cap_change = Cap_change_none;
      break;
    case 2:
      cap_change = Cap_change_tolower;
      break;
    case 3:
      cap_change = Cap_change_toupper;
      break;
    case 4:
      cap_change = Cap_change_firstcap;
      break;
    case 5:
      cap_change = Cap_change_firstcaprestnochange;
      break;
    default:
      cap_change = Cap_change_none;
      break;
  }
  return cap_change;
}


static void ChangeCapChangeGrp (GrouP g)
{
  CapChangeDlgPtr dlg;
  dlg = (CapChangeDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) return;
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static ValNodePtr TestCapChangeDialog (DialoG d)
{
  CapChangeDlgPtr dlg;
  ValNodePtr err_list = NULL;
  Int2 val;

  dlg = (CapChangeDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  val = GetValue (dlg->cap_change);
  if (val < 1 || val > 5) {
    ValNodeAddPointer (&err_list, 0, "cap change");
  }
  return err_list;
}


NLM_EXTERN DialoG 
CapChangeDialog
(GrouP h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  CapChangeDlgPtr dlg;
  GrouP           p;

  dlg = (CapChangeDlgPtr) MemNew (sizeof (CapChangeDlgData));
  p = NormalGroup (h, 5, 0, "Capitalization", programFont, ChangeCapChangeGrp);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestCapChangeDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->cap_change = p;
  SetGroupSpacing (dlg->cap_change, 10, 10);

  RadioButton (dlg->cap_change, "No change");
  RadioButton (dlg->cap_change, "To lower");
  RadioButton (dlg->cap_change, "To upper");
  RadioButton (dlg->cap_change, "First cap, rest lower");
  RadioButton (dlg->cap_change, "First cap, rest no change");
  SetValue (dlg->cap_change, 1);

  return (DialoG) p;  
}




typedef struct parseactiondlg {
 DIALOG_MESSAGE_BLOCK
 DialoG text_portion;
 DialoG src;
 DialoG dst;
 DialoG cap_change;
 DialoG transform;
 ButtoN remove_from_parsed;
 DialoG existing_text;
 Nlm_ChangeNotifyProc change_notify;
 Pointer              change_userdata;
 Nlm_ChangeNotifyProc redraw_notify;
 Pointer              redraw_userdata;
} ParseActionDlgData, PNTR ParseActionDlgPtr;


static void ParseActionToDialog (DialoG d, Pointer data)
{
  ParseActionDlgPtr dlg;
  ParseActionPtr action;

  dlg = (ParseActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  action = (ParseActionPtr) data;
  if (action == NULL) {
    PointerToDialog (dlg->text_portion, NULL);
    PointerToDialog (dlg->src, NULL);
    PointerToDialog (dlg->dst, NULL);
    PointerToDialog (dlg->cap_change, NULL);
    PointerToDialog (dlg->transform, NULL);
    SetStatus (dlg->remove_from_parsed, FALSE);
    SetExistingTextDialogValue(dlg->existing_text, 0);
  } else {
    PointerToDialog (dlg->text_portion, action->portion);
    PointerToDialog (dlg->src, action->src);
    PointerToDialog (dlg->dst, action->dest);
    SetCapChangeDialogValue (dlg->cap_change, action->capitalization);
    PointerToDialog (dlg->transform, action->transform);
    SetStatus (dlg->remove_from_parsed, action->remove_from_parsed);
    SetExistingTextDialogValue(dlg->existing_text, action->existing_text);
  }    
}


static Pointer DialogToParseAction (DialoG d)
{
  ParseActionDlgPtr dlg;
  ParseActionPtr    action;

  dlg = (ParseActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  action = ParseActionNew();
  action->portion = DialogToPointer (dlg->text_portion);
  action->src = DialogToPointer (dlg->src);
  action->dest = DialogToPointer (dlg->dst);
  action->remove_from_parsed = GetStatus (dlg->remove_from_parsed);
  action->capitalization = GetCapChangeDialogValue (dlg->cap_change);
  action->transform = DialogToPointer (dlg->transform);
  action->existing_text = GetExistingTextDialogValue (dlg->existing_text);
  return action;
}


static ValNodePtr TestParseActionDialog (DialoG d)
{
  ParseActionDlgPtr dlg;
  ValNodePtr err_list = NULL;

  dlg = (ParseActionDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  ValNodeLink (&err_list, TestDialog (dlg->text_portion));
  ValNodeLink (&err_list, TestDialog (dlg->src));
  ValNodeLink (&err_list, TestDialog (dlg->dst));
  ValNodeLink (&err_list, TestDialog (dlg->existing_text));
  ValNodeLink (&err_list, TestDialog (dlg->cap_change));
  return err_list;
}


static DialoG 
ParseActionDialogEx 
(GrouP h,
 Boolean indexer_version, 
 Boolean use_existing_text,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata,
 SeqEntryPtr              sep)
{
  ParseActionDlgPtr dlg;
  GrouP                 p, g1, g2, g3;

  dlg = (ParseActionDlgPtr) MemNew (sizeof (ParseActionDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = ParseActionToDialog;
  dlg->fromdialog = DialogToParseAction;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestParseActionDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g3 =HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g3, "Select text", 0, dialogTextHeight, systemFont, 'l');
  dlg->text_portion = TextPortionDialog (g3, TRUE, change_notify, change_userdata);

  g1 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g1, "From", 0, dialogTextHeight, systemFont, 'r');
  dlg->src = ParseSrcDialog (g1, change_notify, change_userdata, sep);
  g2 = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g2, "And place in", 0, dialogTextHeight, systemFont, 'r');
  dlg->dst = ParseDstDialog (g2, change_notify, change_userdata);
  dlg->cap_change = CapChangeDialog (p, NULL, NULL);

  dlg->transform = TextTransformSetDialog (p, change_notify, change_userdata);

  dlg->remove_from_parsed = CheckBox (p, "Remove from parsed field", NULL);

  if (use_existing_text) {
    dlg->existing_text = ExistingTextDialog (p, change_notify, change_userdata);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->cap_change,
                              (HANDLE) dlg->transform,
                              (HANDLE) dlg->remove_from_parsed,
                              (HANDLE) g3,
                              (HANDLE) g1, (HANDLE) g2,
                              (HANDLE) dlg->existing_text,
                              NULL);

  return (DialoG) p;  
}


static DialoG 
ParseActionDialog 
(GrouP h,
 Boolean indexer_version, 
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  return ParseActionDialogEx (h, indexer_version, TRUE, change_notify, change_userdata, NULL);
}


typedef struct editmacroaction {
  PopuP action_type;

  DialoG action_dlg;

  ButtoN accept_btn;
  ButtoN test_btn;
  ButtoN undo_btn;
  ButtoN leave_dlg_up;
  WindoW w;
  Boolean rebuild_window;
  Char    undo_file [PATH_MAX];
  ValNodePtr action_copy;
  Boolean    indexer_version;
  Uint2      entityID;
  Uint1      last_qual_type;
} EditMacroActionData, PNTR EditMacroActionPtr;


/* rebuild dialog with new feature type */
static void ChangeMacroAction (Pointer data)
{
  EditMacroActionPtr e;

  e = (EditMacroActionPtr) data;
  if (e != NULL) {
    e->rebuild_window = TRUE;
  }
}


static void ChangeMacroActionPopup (PopuP p)
{
  ChangeMacroAction (GetObjectExtra (p));
}


static void EnableEditMacroActionAccept (Pointer data)
{
  EditMacroActionPtr e;
  ValNodePtr err_list = NULL;

  e = (EditMacroActionPtr) data;
  if (e != NULL) {
    err_list = TestDialog (e->action_dlg);
    if (err_list == NULL) {
      Enable (e->accept_btn);
      Enable (e->test_btn);
    } else {
      Disable (e->accept_btn);
      Disable (e->test_btn);
    }
    err_list = ValNodeFree (err_list);
  }
}  


typedef DialoG (*MacroActionChoiceDialog) PROTO((GrouP h, Boolean indexer_version, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata));

typedef struct macropopup {
  CharPtr popup_label;
  Int4 action_choice;
  MacroActionChoiceDialog dialog_func;
} MacroPopupData, PNTR MacroPopupPtr;

typedef enum {
  eMacroPopup_AddFeature = 1,
  eMacroPopup_ApplyQualifier,
  eMacroPopup_EditQualifier,
  eMacroPopup_ConvertQualifier,
  eMacroPopup_CopyQualifier,
  eMacroPopup_SwapQualifier,
  eMacroPopup_ParseQualifier,
  eMacroPopup_RemoveQualifier,
  eMacroPopup_ParseText,
  eMacroPopup_RemoveFeature,
  eMacroPopup_EditFeatureLocation,
  eMacroPopup_ConvertFeature,
  eMacroPopup_RemoveDescriptor,
  eMacroPopup_Autodef,
  eMacroPopup_RemoveDuplicateNestedSets,
  eMacroPopup_TrimJunkInPrimerSeqs,
  eMacroPopup_FixUSAandStateAbbreviationsInPublications,
  eMacroPopup_RemoveTrailingStopFromCodingRegions,
  eMacroPopup_SynchronizeCodingRegionPartials,
  eMacroPopup_AdjustCodingRegionsForConsensusSpliceSites,
  eMacroPopup_FixPublicationCapitalization,
  eMacroPopup_RemoveSegGaps,
  eMacroPopup_SortProteinNames,
  eMacroPopup_EditMolInfoFields,
  eMacroPopup_FixSourceCountryQualCapitalization,
  eMacroPopup_FixSourceQualCapitalization,
  eMacroPopup_FormatCollectionDate,
  eMacroPopup_FormatLatLon,
  eMacroPopup_FormatPrimers,
  eMacroPopup_FormatProteinName,
  eMacroPopup_FixSpell,
  eMacroPopup_RemoveDuplicateFeatures,
  eMacroPopup_RemoveLineageSourceNotes,
  eMacroPopup_FixCapsInCommonMusMusculusStrains,
  eMacroPopup_RemoveXrefs,
  eMacroPopup_MakeGeneXrefs,
  eMacroPopup_FixAuthorNames,
  eMacroPopup_UpdateSequences
} EMacroPopup;

static MacroPopupData s_MacroPopupList[] = {
  {"Apply Feature", MacroActionChoice_add_feature, NULL},
  {"Apply Qualifier", MacroActionChoice_aecr, NULL},
  {"Edit Qualifier", MacroActionChoice_aecr, NULL},
  {"Convert Qualifier", MacroActionChoice_aecr, NULL},
  {"Copy Qualifier", MacroActionChoice_aecr, NULL},
  {"Swap Qualifier", MacroActionChoice_aecr, NULL},
  {"Parse Qualifier", MacroActionChoice_aecr, NULL},
  {"Remove Qualifier", MacroActionChoice_aecr, NULL},
  {"Parse Text", MacroActionChoice_parse, ParseActionDialog },
  {"Remove Feature", MacroActionChoice_remove_feature, RemoveFeatureActionDialog },
  {"Edit Feature Location", MacroActionChoice_edit_location, EditFeatureLocationActionDialog },
  {"Convert Feature", MacroActionChoice_convert_feature, ConvertFeatureActionDialog },
  {"Remove Descriptor", MacroActionChoice_remove_descriptor, RemoveDescriptorActionDialog },
  {"Autodef", MacroActionChoice_autodef, AutodefActionDialog },
  {"Remove Duplicate Nested Sets", MacroActionChoice_removesets, NULL},
  {"Trim Junk in Primer Seqs", MacroActionChoice_trim_junk_from_primer_seq, NULL},
  {"Fix USA and state abbreviations in publications", MacroActionChoice_fix_usa_and_states, NULL},
  {"Remove trailing * from Coding Regions", MacroActionChoice_trim_stop_from_complete_cds, NULL},
  {"Synchronize Coding Region Partials", MacroActionChoice_synchronize_cds_partials, NULL},
  {"Adjust coding regions for consensus splice sites", MacroActionChoice_adjust_for_consensus_splice, NULL},
  {"Fix publication capitalization", MacroActionChoice_fix_pub_caps, FixPubCapsDialog },
  {"Remove Seg-gaps", MacroActionChoice_remove_seg_gaps, NULL},
  {"Sort Protein Names", MacroActionChoice_sort_fields, SortFieldsDialog },
  {"Edit MolInfo Fields", MacroActionChoice_apply_molinfo_block, MolInfoBlockDialog },
  {"Fix source country qual capitalization", MacroActionChoice_fix_caps, FixCapsSourceCountryDialog },
  {"Fix source qual capitalization", MacroActionChoice_fix_caps, FixSrcQualCapsActionDialog },
  {"Fix collection-date format", MacroActionChoice_fix_format, FormatCollectionDateDialog } ,
  {"Fix lat-lon format", MacroActionChoice_fix_format, FormatLatLonDialog } ,
  {"Fix primer format", MacroActionChoice_fix_format, FormatPrimerDialog } ,
  {"Fix protein name format", MacroActionChoice_fix_format, FormatProteinNameDialog } ,
  {"Fix spelling", MacroActionChoice_fix_spell, NULL } ,
  {"Remove Duplicate Features", MacroActionChoice_remove_duplicate_features, RemoveDuplicateFeatActionDialog } ,
  {"Remove Lineage Source Notes", MacroActionChoice_remove_lineage_notes, NULL } ,
  {"Fix capitalization in common Mus musculus strains", MacroActionChoice_fix_caps, MusMusculusActionDialog } ,
  {"Remove Xrefs", MacroActionChoice_remove_xrefs, RemoveXrefsDialog } ,
  {"Make gene xrefs from features", MacroActionChoice_make_gene_xrefs, MakeGeneXrefDialog } ,
  {"Make Barcode Xrefs", MacroActionChoice_make_bold_xrefs, NULL } ,
  {"Fix Author Names", MacroActionChoice_fix_author, AuthorFixActionDialog } ,
  {"Update Sequences", MacroActionChoice_update_sequences, UpdateSequencesActionDialog } ,
  {"Set trans-splicing exception in genes", MacroActionChoice_add_trans_splicing, NULL } ,
  {"Remove invalid EC_numbers", MacroActionChoice_remove_invalid_ecnumbers, NULL }
};

#define NUM_MacroPopup sizeof (s_MacroPopupList) / sizeof (MacroPopupData)


static Uint1 GetActionChoiceForPopupChoice (Int2 popup_choice)
{
  if (popup_choice > NUM_MacroPopup || popup_choice < 1) {
    return -1;
  } else {
    return s_MacroPopupList[popup_choice - 1].action_choice;
  }
}

static Int2 GetPopupChoiceForActionChoice (Uint1 action_choice)
{
  Int2 j;

  for (j = 0; j < NUM_MacroPopup; j++) {
    if (s_MacroPopupList[j].action_choice == action_choice) {
      return j + 1;
    }
  }
  return 1;
}


static void RunMacroInEditor (ButtoN b)
{
  EditMacroActionPtr e;
  ValNodePtr  action;
  ValNodePtr  sep_list, vnp;
  SeqEntryPtr sep;
  Uint2       entityID;
  LogInfoPtr  lip;

  e = (EditMacroActionPtr) GetObjectExtra (b);
  if (e == NULL) return;

  action = ValNodeNew (NULL);
  action->choice = GetActionChoiceForPopupChoice(GetValue (e->action_type));
  action->data.ptrvalue = DialogToPointer (e->action_dlg);

  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    Message (MSG_ERROR, "No records open!");
  } else if (sep_list->next != NULL 
    && ANS_CANCEL == Message (MSG_OKC, "You have more than one record open - run macro for all open records?")) {
    /* do nothing */
  } else {
    WatchCursor();
    Update();
    lip = OpenLog ("Macro Actions");
    for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
      sep = vnp->data.ptrvalue;
      entityID = ObjMgrGetEntityIDForChoice(sep);
      lip->data_in_log |= ApplyMacroToSeqEntryEx (sep, action, lip->fp, Sequin_GlobalAlign2Seq);
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
    sep_list = ValNodeFree (sep_list);
    ArrowCursor ();
    Update ();  
    if (!lip->data_in_log) {
      fprintf (lip->fp, "Macro had no effect\n");
      lip->data_in_log = TRUE;
    }
    CloseLog (lip);
    lip = FreeLog (lip);
  }  

  action = MacroActionChoiceFree (action);
  Enable (e->undo_btn);
}


static Boolean SaveBackupSeqEntryList (CharPtr path)
{
  ValNodePtr sep_list, vnp;
  AsnIoPtr   aip;

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) {
    Message (MSG_ERROR, "Unable to open file for backup");
    return FALSE;
  }

  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    AsnIoClose (aip);
    return FALSE;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    SeqEntryAsnWrite (vnp->data.ptrvalue, aip, NULL);
  }

  AsnIoClose (aip);
  sep_list = ValNodeFree (sep_list);
  return TRUE;  
}


static ValNodePtr ReadBackupSeqEntryList (CharPtr path)
{
  AsnIoPtr   aip;
  ValNodePtr sep_list = NULL;
  SeqEntryPtr sep;

  aip = AsnIoOpen (path, "r");
  if (aip == NULL) {
    Message (MSG_ERROR, "Unable to open %s", path);
  } else {
    while ((sep = SeqEntryAsnRead (aip, NULL)) != NULL) {
      ValNodeAddPointer (&sep_list, 0, sep);
    }
    AsnIoClose (aip);
  }
  return sep_list;  
}


static void RestoreList (ValNodePtr curr_sep_list, ValNodePtr backup_sep_list)
{
  ValNodePtr vnp_c, vnp_b;
  SeqEntryPtr sep_c, sep_b;
  Uint2       entityID;

  if (curr_sep_list == NULL) return;
  if (ValNodeLen (curr_sep_list) != ValNodeLen (backup_sep_list)) {
    Message (MSG_ERROR, "Backup list does not match current list.  Unable to undo.");
    return;
  }

  for (vnp_c = curr_sep_list, vnp_b = backup_sep_list; vnp_c != NULL && vnp_b != NULL; vnp_c = vnp_c->next, vnp_b = vnp_b->next) {
    sep_c = vnp_c->data.ptrvalue;
    sep_b = vnp_b->data.ptrvalue;
    SeqEntrySetScope (NULL);
    ReplaceSeqEntryWithSeqEntry (sep_c, sep_b, TRUE);
    entityID = ObjMgrGetEntityIDForChoice (sep_c);
    ObjMgrSetDirtyFlag (entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
  }
  Update();
}


static void UndoMacroActionTest (ButtoN b)
{
  EditMacroActionPtr e;
  ValNodePtr         curr_list;
  ValNodePtr         backup_list;
  
  e = (EditMacroActionPtr) GetObjectExtra (b);
  if (e == NULL) return;

  curr_list = GetViewedSeqEntryList ();
  if (curr_list == NULL) {
    Message (MSG_ERROR, "No records to undo");
    Disable (e->undo_btn);
    return;
  } else {
    backup_list = ReadBackupSeqEntryList (e->undo_file);
    RestoreList (curr_list, backup_list);
    Disable (e->undo_btn);
    curr_list = ValNodeFree (curr_list);
    backup_list = ValNodeFree (backup_list);
  }
  
}


static WindoW 
BuildEditMacroActionWindow 
(EditMacroActionPtr d,
 ModalAcceptCancelPtr acp,
 Boolean adding_new)
{
  AECRActionPtr         aecr = NULL;
  FixCapsActionPtr      fc;
  FixFormatActionPtr    fa;
  ApplyFeatureActionPtr apply_feat = NULL;
  WindoW                w;
  GrouP                 h, c;
  ButtoN                b;
  Int2                  action_type;

  w = MovableModalWindow(-20, -13, -10, -10, "Edit Action", NULL);
#ifndef WIN_MAC
  CreateStandardEditMenu (w);
#endif

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  d->action_type = PopupList (h, TRUE, ChangeMacroActionPopup);
  SetObjectExtra (d->action_type, d, NULL);
  for (action_type = 0; action_type < NUM_MacroPopup; action_type++) {
    PopupItem (d->action_type, s_MacroPopupList[action_type].popup_label);
  }

  if (d->action_copy == NULL) {
    SetValue (d->action_type, 1);
  } else {
    switch (d->action_copy->choice) {
      case MacroActionChoice_add_feature:
        SetValue (d->action_type, eMacroPopup_AddFeature);
        apply_feat = (ApplyFeatureActionPtr) d->action_copy->data.ptrvalue;
        break;
      case MacroActionChoice_aecr:
        aecr = (AECRActionPtr) d->action_copy->data.ptrvalue;
        if (aecr == NULL || aecr->action == NULL) {
          SetValue (d->action_type, eMacroPopup_ApplyQualifier);
        } else {
          switch (aecr->action->choice) {
            case ActionChoice_apply:
              SetValue (d->action_type, eMacroPopup_ApplyQualifier);
              break;
            case ActionChoice_edit:
              SetValue (d->action_type, eMacroPopup_EditQualifier);
              break;
            case ActionChoice_convert:
              SetValue (d->action_type, eMacroPopup_ConvertQualifier);
              break;
            case ActionChoice_copy:
              SetValue (d->action_type, eMacroPopup_CopyQualifier);
              break;
            case ActionChoice_swap:
              SetValue (d->action_type, eMacroPopup_SwapQualifier);
              break;
            case ActionChoice_parse:
              SetValue (d->action_type, eMacroPopup_ParseQualifier);
              break;
            case ActionChoice_remove:
              SetValue (d->action_type, eMacroPopup_RemoveQualifier);
              break;
            default:
              SetValue (d->action_type, eMacroPopup_ApplyQualifier);
              break;
          }
        }
        break;
      case MacroActionChoice_fix_caps:
        fc = (FixCapsActionPtr) d->action_copy->data.ptrvalue;
        if (fc == NULL) {
          SetValue (d->action_type, eMacroPopup_FixSourceCountryQualCapitalization);
        } else {
          switch (fc->choice) {
            case FixCapsAction_pub:
              SetValue (d->action_type, eMacroPopup_FixPublicationCapitalization);
              break;
            case FixCapsAction_src_country:
              SetValue (d->action_type, eMacroPopup_FixSourceCountryQualCapitalization);
              break;
            case FixCapsAction_mouse_strain:
              SetValue (d->action_type, eMacroPopup_FixCapsInCommonMusMusculusStrains);
              break;
            case FixCapsAction_src_qual:
              SetValue (d->action_type, eMacroPopup_FixSourceQualCapitalization);
              break;
            default:
              SetValue (d->action_type, eMacroPopup_FixSourceCountryQualCapitalization);
              break;
          }
        }
        break;
      case MacroActionChoice_fix_format:
        fa = (FixFormatActionPtr) d->action_copy->data.ptrvalue;
        if (fa == NULL) {
          SetValue (d->action_type, eMacroPopup_FormatCollectionDate);
        } else {
          switch (fa->choice) {
            case FixFormatAction_collection_date:
              SetValue (d->action_type, eMacroPopup_FormatCollectionDate);
              break;
            case FixFormatAction_lat_lon:
              SetValue (d->action_type, eMacroPopup_FormatLatLon);
              break;
            case FixFormatAction_primers:
              SetValue (d->action_type, eMacroPopup_FormatPrimers);
              break;
            case FixFormatAction_protein_name:
              SetValue (d->action_type, eMacroPopup_FormatProteinName);
              break;
            default:
              SetValue (d->action_type, eMacroPopup_FormatCollectionDate);
              break;
          }
        }
        break;
      default:
        SetValue (d->action_type, GetPopupChoiceForActionChoice(d->action_copy->choice));
        break;
    }
  } 
  d->action_dlg = NULL;

  action_type = GetValue (d->action_type);
  /* set up remaining controls */
  switch (action_type) {
    case eMacroPopup_AddFeature:  
      d->action_dlg = ApplyFeatureActionDialog (h, apply_feat, d->indexer_version, EnableEditMacroActionAccept, d, ChangeMacroAction, d);
      PointerToDialog (d->action_dlg, apply_feat);
      break;
    case eMacroPopup_ApplyQualifier:
    case eMacroPopup_EditQualifier:
    case eMacroPopup_ConvertQualifier:
    case eMacroPopup_CopyQualifier:
    case eMacroPopup_SwapQualifier:
    case eMacroPopup_ParseQualifier:
    case eMacroPopup_RemoveQualifier:
      d->action_dlg = EditAECRActionDialog (h, d->indexer_version, TRUE, EnableEditMacroActionAccept, d, EnableEditMacroActionAccept, d);
      PointerToDialog (d->action_dlg, aecr);
      break;
    default:
      if (s_MacroPopupList[action_type - 1].dialog_func != NULL) {
        d->action_dlg = s_MacroPopupList[action_type - 1].dialog_func (h, d->indexer_version, EnableEditMacroActionAccept, d);
        PointerToDialog (d->action_dlg, d->action_copy->data.ptrvalue);
      }
      break;
  }

  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  d->accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (d->accept_btn, acp, NULL);
  if (d->undo_file[0] != 0) {
    d->test_btn = PushButton (c, "Test", RunMacroInEditor);
    SetObjectExtra (d->test_btn, d, NULL);
    d->undo_btn = PushButton (c, "Undo", UndoMacroActionTest);
    SetObjectExtra (d->undo_btn, d, NULL);
    Disable (d->undo_btn);
  }
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, acp, NULL);

  if (adding_new) {
    if (d->undo_file[0] == 0) {
      d->leave_dlg_up = CheckBox (h, "Leave dialog up", NULL);
    } else {
      d->leave_dlg_up = CheckBox (h, "Leave dialog up (and add more macro actions)", NULL);
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) c, 
                              (HANDLE) d->action_type,
                              (HANDLE) d->action_dlg,
                              (HANDLE) d->leave_dlg_up,
                              NULL);
  EnableEditMacroActionAccept (d);
  return w;
}


static void AdjustMacroActionForRebuild (EditMacroActionPtr d)
{
  Int2                  action_type;
  Uint1 new_action_choice;
  ValNodePtr vnp;

  /* copy data from dialog if compatible */
  action_type = GetValue (d->action_type);
  new_action_choice = GetActionChoiceForPopupChoice (action_type);    
  if (d->action_copy != NULL && d->action_copy->choice == MacroActionChoice_aecr) {
    d->last_qual_type = FieldTypeChoiceFromEditAECRActionDialog (d->action_dlg);
  }
  if (d->action_copy != NULL && d->action_copy->choice == MacroActionChoice_aecr && new_action_choice == MacroActionChoice_aecr) {
    d->action_copy = MacroActionChoiceFree (d->action_copy);
    d->action_copy = ValNodeNew (NULL);
    d->action_copy->choice = MacroActionChoice_aecr;
    d->action_copy->data.ptrvalue = BuildDefaultAECRAction (ActionTypeChoiceFromPopupValue (action_type), d->last_qual_type);
  } else if (d->action_copy != NULL && new_action_choice == MacroActionChoice_fix_format) {
    d->action_copy = MacroActionChoiceFree (d->action_copy);
    d->action_copy = ValNodeNew (NULL);
    d->action_copy->choice = MacroActionChoice_fix_format;
    vnp = ValNodeNew (NULL);
    d->action_copy->data.ptrvalue = vnp;
    switch (action_type) {
      case eMacroPopup_FormatCollectionDate:
        vnp->choice = FixFormatAction_collection_date;
        break;
      case eMacroPopup_FormatLatLon:
        vnp->choice = FixFormatAction_lat_lon;
        break;
      case eMacroPopup_FormatPrimers:
        vnp->choice = FixFormatAction_primers;
        break;
      case eMacroPopup_FormatProteinName:
        vnp->choice = FixFormatAction_protein_name;
        break;
      default:
        vnp->choice = FixFormatAction_collection_date;
        break;
    }
  } else if (d->action_copy != NULL && d->action_copy->choice == MacroActionChoice_fix_caps) {
    d->action_copy = MacroActionChoiceFree (d->action_copy);
    d->action_copy = ValNodeNew (NULL);
    d->action_copy->choice = MacroActionChoice_fix_caps;
    vnp = ValNodeNew (NULL);
    d->action_copy->data.ptrvalue = vnp;
    switch (action_type) {
      case eMacroPopup_FixPublicationCapitalization:
        vnp->choice = FixCapsAction_pub;
        break;
      case eMacroPopup_FixSourceCountryQualCapitalization:
        vnp->choice = FixCapsAction_src_country;
        break;
      case eMacroPopup_FixCapsInCommonMusMusculusStrains:
        vnp->choice = FixCapsAction_mouse_strain;
        break;
      case eMacroPopup_FixSourceQualCapitalization:
        vnp->choice = FixCapsAction_src_qual;
        break;
      default:
        vnp->choice = FixCapsAction_src_country;
        break;
    }
  } else if (d->action_copy != NULL && d->action_copy->choice == new_action_choice) {
    d->action_copy = MacroActionChoiceFree (d->action_copy);
    d->action_copy = ValNodeNew (NULL);
    d->action_copy->choice = new_action_choice;
    d->action_copy->data.ptrvalue = DialogToPointer (d->action_dlg);
  } else {
    d->action_copy = MacroActionChoiceFree (d->action_copy);
    d->action_copy = ValNodeNew (NULL);
    d->action_copy->choice = new_action_choice;
    /* note - only need to set default actions where "NULL" isn't sufficient */
    switch (new_action_choice) {
      case MacroActionChoice_aecr:
        d->action_copy->data.ptrvalue = BuildDefaultAECRAction (ActionTypeChoiceFromPopupValue (action_type), d->last_qual_type);
        break;
      case MacroActionChoice_remove_feature:
        d->action_copy->data.ptrvalue = BuildDefaultRemoveFeatureAction ();
        break;
      case MacroActionChoice_edit_location:
        d->action_copy->data.ptrvalue = BuildDefaultEditFeatureLocationAction (Macro_feature_type_gene);
        break;
      case MacroActionChoice_fix_caps:
        vnp = ValNodeNew (NULL);
        d->action_copy->data.ptrvalue = vnp;
        switch (action_type) {
          case eMacroPopup_FixPublicationCapitalization:
            vnp->choice = FixCapsAction_pub;
            break;
          case eMacroPopup_FixSourceCountryQualCapitalization:
            vnp->choice = FixCapsAction_src_country;
            break;
          case eMacroPopup_FixCapsInCommonMusMusculusStrains:
            vnp->choice = FixCapsAction_mouse_strain;
            break;
          case eMacroPopup_FixSourceQualCapitalization:
            vnp->choice = FixCapsAction_src_qual;
            break;
          default:
            vnp->choice = FixCapsAction_src_country;
            break;
        }
        break;
    }
  }
}


static Boolean EditMacroAction (ValNodePtr action, Boolean indexer_version)
{
  EditMacroActionData   d;
  ModalAcceptCancelData acd;
  Boolean               rval = FALSE;
  WindoW                w, repl_w;
  ValNodePtr            tmp;
  
  if (action == NULL) return FALSE;

  d.indexer_version = indexer_version;
  d.accept_btn = NULL;
  d.test_btn = NULL;
  d.leave_dlg_up = NULL;

  TmpNam (d.undo_file);
  if (!SaveBackupSeqEntryList (d.undo_file)) {
    d.undo_file[0] = 0;
  }

  d.action_copy = AsnIoMemCopy (action, (AsnReadFunc) MacroActionChoiceAsnRead, (AsnWriteFunc) MacroActionChoiceAsnWrite);
  w = BuildEditMacroActionWindow (&d, &acd, FALSE);
  d.rebuild_window = FALSE;
  d.last_qual_type = FieldType_source_qual;
      
  Show (w);
  Select (w);
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    if (d.rebuild_window) {
      /* copy data from dialog if compatible */
      AdjustMacroActionForRebuild (&d);

      /* replace window */
      Remove (w);
      repl_w = BuildEditMacroActionWindow (&d, &acd, FALSE);
      w = repl_w;
      Show (w);
      Select (w);
      d.rebuild_window = FALSE;
    }        
    Update ();
  }
  ProcessAnEvent ();
  if (!acd.cancelled)
  {
    rval = TRUE;
    tmp = ValNodeNew (NULL);
    tmp->choice = action->choice;
    tmp->data.ptrvalue = action->data.ptrvalue;
    action->data.ptrvalue = NULL;
    tmp = MacroActionChoiceFree(tmp);
    action->choice = GetActionChoiceForPopupChoice(GetValue (d.action_type));
    action->data.ptrvalue = DialogToPointer (d.action_dlg);
  }
  Remove (w);
  FileRemove (d.undo_file);

  d.action_copy = MacroActionChoiceFree (d.action_copy);
  return rval;
}


static void AddMacroActions (MacroEditorFormPtr f, Int2 item)
{
  EditMacroActionData   d;
  ModalAcceptCancelData acd;
  Boolean               rval = FALSE;
  WindoW                w, repl_w;
  ValNodePtr            new_action, prev_action = NULL;
  Int2                  pos;
  Boolean               done = FALSE, leave_up;
  Int4                  scroll_pos;
  BaR                   sb_vert;

  if (f == NULL) return;

  if (item > 0 && f->macro_list != NULL) {
    pos = 1;
    prev_action = f->macro_list;
    while (pos < item && prev_action->next != NULL) {
      pos++;
      prev_action = prev_action->next;
    }
  }  
  
  d.action_copy = BuildDefaultNewMacroAction ();
  d.accept_btn = NULL;
  d.test_btn = NULL;
  d.leave_dlg_up = NULL;
  d.indexer_version = f->indexer_version;
  if (f->input_entityID != 0) {
    TmpNam (d.undo_file);
    if (!SaveBackupSeqEntryList (d.undo_file)) {
      d.undo_file[0] = 0;
    }
  }
  w = BuildEditMacroActionWindow (&d, &acd, TRUE);
  d.rebuild_window = FALSE;
  d.last_qual_type = FieldType_source_qual;
    
  Show (w);
  Select (w);
  while (!done) {
    acd.accepted = FALSE;
    acd.cancelled = FALSE;
    while (!acd.accepted && ! acd.cancelled)
    {
      ProcessExternalEvent ();
      if (d.rebuild_window) {
        /* copy data from dialog if compatible */
        AdjustMacroActionForRebuild (&d);

        /* replace window */
        leave_up = GetStatus (d.leave_dlg_up);
        Remove (w);
        repl_w = BuildEditMacroActionWindow (&d, &acd, TRUE);
        SetStatus (d.leave_dlg_up, leave_up);
        w = repl_w;
        Show (w);
        Select (w);
        d.rebuild_window = FALSE;
      }        
      Update ();
    }
    ProcessAnEvent ();
    if (acd.cancelled) {
      done = TRUE;
    } else {
      /* get current scroll position */
      sb_vert = GetSlateVScrollBar ((SlatE) f->macro_summary);
      scroll_pos = GetBarValue (sb_vert);
      /* we will want to increase the scroll position after each addition
       * note that we need to get the scroll bar and check the initial position
       * each time - if there was no scroll bar after the last update, scroll_pos
       * needs to be zero to start.
       */
      scroll_pos++;

      rval = TRUE;
      /* create new action based on contents of dialog */
      new_action = ValNodeNew (NULL);
      new_action->choice = GetActionChoiceForPopupChoice(GetValue (d.action_type));
      new_action->data.ptrvalue = DialogToPointer (d.action_dlg);
      /* add action to macro list */
      if (prev_action == NULL) {
        /* put at start of list */
        new_action->next = f->macro_list;
        f->macro_list = new_action;
      } else {
        /* insert in list */
        new_action->next = prev_action->next;
        prev_action->next = new_action;
      }
      /* next action will be inserted after this one */
      prev_action = new_action; 
      /* update summary */
      UpdateMacroSummary (f, scroll_pos);

      /* done if leave_dlg_up not checked */
      if (!GetStatus (d.leave_dlg_up)) {
        done = TRUE;
      }
    }
  }
  Remove (w);
  d.action_copy = MacroActionChoiceFree (d.action_copy);
}


/* for executing macro actions without creating a macro script */

typedef struct sampledialog {
  DIALOG_MESSAGE_BLOCK
  DoC                doc;
  Int2               start_item;
  Int2               start_row;
  Int2               start_col;
  Int2               end_item;
  Int2               end_row;
  Int2               end_col;

} SampleDialogData, PNTR SampleDialogPtr;


static void PointerToSampleDialog (DialoG d, Pointer data)
{
  SampleDialogPtr    dlg;
  ValNodePtr         sample_list;
  ValNodePtr         vnp;
  AECRSamplePtr      sample;
  CharPtr            line_txt;
  CharPtr            line_same_fmt = "%s\t(%d same)\t%s\n";
  CharPtr            line_mix_fmt = "%s\t(%d mixed)\t%s\n";
  CharPtr            fmt;
  RecT            r;
  ParData         ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData         ColFmt[] = 
  {
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };
  Int4            max_char_per_line;
  CharPtr         field_name;


  dlg = (SampleDialogPtr) GetObjectExtra (d);

  if (dlg == NULL) return;
  Reset (dlg->doc);

  sample_list = (ValNodePtr) data;
  if (sample_list == NULL) return;

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  ColFmt[0].pixWidth = (r.right - r.left) / 3;
  ColFmt[1].pixWidth = (r.right - r.left) / 4;
  ColFmt[2].pixWidth = (r.right - r.left)  - ColFmt[0].pixWidth - ColFmt[1].pixWidth;

  SelectFont (programFont);
  max_char_per_line = ColFmt[2].pixWidth / CharWidth ('0');

  for (vnp = sample_list; vnp != NULL; vnp = vnp->next) {
    sample = vnp->data.ptrvalue;
    if (sample != NULL) {
      field_name = SummarizeFieldType (sample->field);
      
      if (sample->all_same) {
        fmt = line_same_fmt;
      } else {
        fmt = line_mix_fmt;
      }

      line_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (field_name) + StringLen (sample->first_value) + 15));
      sprintf (line_txt, fmt, field_name, sample->num_found, sample->first_value == NULL ? "" : sample->first_value);
      AppendText (dlg->doc, line_txt, &ParFmt, ColFmt, programFont);
      MemFree (line_txt);
      field_name = MemFree (field_name);
    }
  }
  InvalDocument (dlg->doc);
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


static DialoG SampleDialog (GrouP h)
{
  SampleDialogPtr dlg;
  GrouP           p;
  PrompT          ppt;
  
  dlg = (SampleDialogPtr) MemNew (sizeof (SampleDialogData));
  if (dlg == NULL) 
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);

  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = PointerToSampleDialog;
  dlg->fromdialog = NULL;

  ppt = StaticPrompt (p, "Sample Values", 0, dialogTextHeight, programFont, 'l');

  dlg->doc = DocumentPanel (p, stdCharWidth * 27, stdLineHeight * 8);
  SetObjectExtra (dlg->doc, dlg, NULL);
  SetDocProcs (dlg->doc, SampleDialogOnClick, SampleDialogOnDrag, SampleDialogOnRelease, NULL);  
  SetDocShade (dlg->doc, NULL, NULL, SampleDialogInvert, NULL);
  SetDocAutoAdjust (dlg->doc, TRUE);
  SetSlateChar ((SlatE) dlg->doc, SampleDialogOnKey);

  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) dlg->doc, NULL);
    
  return (DialoG) p;  
}


static Int4 GetSampleDialogValueLen (DialoG d)
{
  SampleDialogPtr dlg;
  RecT            r;
  Int4            max_char_per_line;

  dlg = (SampleDialogPtr) GetObjectExtra (d);
  if (dlg == NULL) return 0;

  ObjectRect (dlg->doc, &r);
  InsetRect (&r, 4, 4);
  SelectFont (programFont);
  max_char_per_line = (r.right - r.left - ((r.right - r.left) / 3) - ((r.right - r.left) / 4)) / CharWidth ('0');
  return max_char_per_line;
}


typedef struct onemacroaction {
  FORM_MESSAGE_BLOCK

  DialoG sample_doc;
  PopuP  action_type;
  DialoG action_dlgs[7];
  DialoG constraint_dlg;

  /* groups and dialogs for selecting field types */
  DialoG single_qual_type_dlg;
  DialoG single_field;
  DialoG single_field_remove;
  GrouP  single_field_grp;
  DialoG pair_qual_type_dlg;
  DialoG field_pair;
  DialoG field_pair_convert;
  GrouP  field_pair_grp;

  ButtoN also_change_mrna;

  ButtoN accept_btn;
  ButtoN leave_dlg_up;

  ButtoN clear_on_change;

  WindoW  sample_window;
  Boolean indexer_version;
  Boolean no_callback;
  Int2   prev_action_type;
  Int2    win_num;
} OneMacroActionData, PNTR OneMacroActionPtr;

/* list for keeping track of macro action windows and their sample dialogs */
static ValNodePtr macro_action_window_list = NULL;

static void RemoveMacroActionWindow (WindoW w)
{
  ValNodePtr vnp, prev = NULL;
  OneMacroActionPtr frm;

  vnp = macro_action_window_list;
  while (vnp != NULL && vnp->data.ptrvalue != w) {
    prev = vnp;
    vnp = vnp->next;
  }
  if (vnp != NULL) {
    if (prev == NULL) {
      macro_action_window_list = vnp->next;
    } else {
      prev->next = vnp->next;
    }
    vnp->next = NULL;
    vnp = ValNodeFree (vnp);
  }
  frm = (OneMacroActionPtr) GetObjectExtra (w);
  if (frm != NULL && frm->sample_window != NULL) {
    Remove (frm->sample_window);
  }

  /* clean up userdata */
  ObjMgrFreeUserData (0, 0, frm->proctype, frm->userkey);
  Remove (w);
}


static Int2 AddMacroActionWindow (WindoW w)
{
  Int2 max = 0;
  ValNodePtr vnp;
  
  for (vnp = macro_action_window_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice > max) {
      max = vnp->choice;
    }
  }
  vnp = ValNodeNew (NULL);
  vnp->choice = max + 1;
  vnp->data.ptrvalue = w;
  vnp->next = macro_action_window_list;
  macro_action_window_list = vnp;
  return max + 1;
}


/* for creating, updating, and closing the sample window */
static void TruncateStringAfterTwoLines (CharPtr str, Uint4 char_per_line)
{
  CharPtr cp;
  Uint4    string_len;
  
  if (StringHasNoText (str)) {
    return;
  }
  
  string_len = StringLen (str);
  if (string_len < char_per_line) {
    return;
  }
  
  cp = str + char_per_line;
  while (! isspace ((Int4)(*cp)) && cp > str) {
    cp--;
  }
  if (cp == str) {
    if (string_len > 2 * char_per_line) {
      str [2 * char_per_line] = 0;
    }
  } else {
    if (StringLen (cp) > char_per_line) {
      cp[char_per_line] = 0;
    }
  }  
}


static void RefreshAECRSample (ButtoN b)
{
  OneMacroActionPtr d;
  Int2              val;
  Uint1             qual_type;
  ValNodePtr        sample_list = NULL, vnp;
  AECRSamplePtr     sample;
  Int4              max_char_per_line;
  
  d = (OneMacroActionPtr) GetObjectExtra (b);

  if (d == NULL || d->sample_doc == NULL) return;

  max_char_per_line = GetSampleDialogValueLen (d->sample_doc);
  val = GetValue (d->action_type);
  if (val >= 1 && val <= 7) {
    qual_type = FieldTypeChoiceFromAECRActionDlg (d->action_dlgs[val - 1]);
    sample_list = GetAECRSampleListForSeqEntry (qual_type, GetTopSeqEntryForEntityID (d->input_entityID));
  }

  for (vnp = sample_list; vnp != NULL; vnp = vnp->next) {
    sample = vnp->data.ptrvalue;
    if (sample != NULL) {
      TruncateStringAfterTwoLines (sample->first_value, max_char_per_line);
    }
  }

  PointerToDialog (d->sample_doc, sample_list);
  sample_list = AECRSampleListFree (sample_list);

}


static void ExportAECRSample (ButtoN b)
{
  OneMacroActionPtr d;
  FILE               *fp;
  Char               path [PATH_MAX];
  Int2               val;
  Uint1              qual_type;

  d = (OneMacroActionPtr) GetObjectExtra (b);

  if (d == NULL || d->sample_doc == NULL) return;
  val = GetValue (d->action_type);
  if (val >= 1 && val <= 7) {
    if (GetOutputFileName (path, sizeof (path), NULL)) {
      fp = FileOpen (path, "w");
      if (fp == NULL) {
        Message (MSG_ERROR, "Unable to open %s", path);
      } else {
        qual_type = FieldTypeChoiceFromAECRActionDlg (d->action_dlgs[val - 1]);
        GetAECRExistingTextList (qual_type, GetTopSeqEntryForEntityID (d->input_entityID), fp);
        FileClose (fp);
      }
    }
  }
}


static void CloseMacroSampleWindowProc (WindoW w)
{
  OneMacroActionPtr     frm;

  frm = (OneMacroActionPtr) GetObjectExtra (w);
  if (frm != NULL) {
    frm->sample_window = NULL;
    frm->sample_doc = NULL;
  }
  Remove (w);  
}


static void CloseSampleWindow (ButtoN b)
{
  OneMacroActionPtr     frm;

  frm = (OneMacroActionPtr) GetObjectExtra (b);
  if (frm != NULL) {
    CloseMacroSampleWindowProc (frm->sample_window);
  }
}


static void ShowMacroSampleWindow (ButtoN b)
{
  OneMacroActionPtr     frm;
  GrouP                 h, g;
  ButtoN                b2;
  CharPtr               title_fmt = "Sample Values (%d)";
  CharPtr               title;

  frm = (OneMacroActionPtr) GetObjectExtra (b);
  if (frm != NULL) {
    if (frm->sample_window == NULL) {
      title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + 15));
      sprintf (title, title_fmt, frm->win_num);
      frm->sample_window = FixedWindow (-20, -13, -10, -10, title, CloseMacroSampleWindowProc);
      title = MemFree (title);
      SetObjectExtra (frm->sample_window, frm, NULL);
      h = HiddenGroup (frm->sample_window, -1, 0, NULL);
      SetGroupSpacing (h, 10, 10);
      frm->sample_doc = SampleDialog (h);
      g = HiddenGroup (h, 3, 0, NULL);
      b2 = PushButton (g, "Refresh Sample", RefreshAECRSample);
      SetObjectExtra (b2, frm, NULL);
      b2 = PushButton (g, "Export Sample", ExportAECRSample);
      SetObjectExtra (b2, frm, NULL);
      b2 = PushButton (g, "Close", CloseSampleWindow);
      SetObjectExtra (b2, frm, NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) frm->sample_doc, (HANDLE) g, NULL);
    }
    RefreshAECRSample (frm->accept_btn);
    Show (frm->sample_window);
    Select (frm->sample_window);
  }
}


/* EnableSingleMacroActionAccept
 * This function should be called whenever anything on the form changes.
 */
static void EnableSingleMacroActionAccept (Pointer data)
{
  OneMacroActionPtr e;
  ValNodePtr err_list = NULL;
  Int2       val;
  IteM       i;

  e = (OneMacroActionPtr) data;
  if (e != NULL) {
    val = GetValue (e->action_type);
    if (val >= 1 && val <= 7) {
      err_list = TestDialog (e->action_dlgs[val - 1]);
      ValNodeLink (&err_list, TestDialog (e->constraint_dlg));
    }
    i = FindFormMenuItem ((BaseFormPtr) e, VIB_MSG_ACCEPT);
    if (val < 1 || val > 7) {
      Disable (e->accept_btn);
      SafeDisable (i);
    } else {
      if (err_list == NULL) {
        Enable (e->accept_btn);
        SafeEnable (i);
      } else {
        Disable (e->accept_btn);
        SafeDisable (i);
      }
    }
    err_list = ValNodeFree (err_list);
  }
}  


static CharPtr action_names[] = {
  "Apply Qualifier" ,
  "Edit Qualifier" ,
  "Convert Qualifier" ,
  "Copy Qualifier" ,
  "Swap Qualifier" ,
  "Parse Qualifier" ,
  "Remove Qualifier" };

static CharPtr GetQualTypeName (Uint1 qual_type)
{
  CharPtr txt = NULL;
  switch (qual_type) {
    case FieldType_source_qual:
      txt = "Source";
      break;
    case FieldType_feature_field:
      txt = "Feature";
      break;
    case FieldType_cds_gene_prot:
      txt = "CDS-Gene-Prot";
      break;
    case FieldType_molinfo_field:
      txt = "MolInfo";
      break;
    case FieldType_pub:
      txt = "Pub";
      break;
    case FieldType_rna_field:
      txt = "RNA";
      break;
    case FieldType_struc_comment_field:
      txt = "Structured Comment Field";
      break;
    case FieldType_misc:
      txt = "Misc";
      break;
    default:
      txt = "Unknown";
      break;
  }
  return txt;
}


/* AdjustSingleMacroActionTitle
 * This function should be called whenever the action or the qual type changes.
 */
static void AdjustSingleMacroActionTitle (OneMacroActionPtr d)
{
  Int2              val;
  Uint1             qual_type;
  CharPtr           qual_name, cp;
  CharPtr           new_title;
  CharPtr           title_fmt = " %s Qualifier (%d)";

  if (d == NULL) {
    return;
  }

  val = GetValue (d->action_type);
  if (val < 1 || val > 7) {
    return;
  }

  /* change window title */
  qual_type = FieldTypeChoiceFromAECRActionDlg (d->action_dlgs[val - 1]);
  qual_name = GetQualTypeName(qual_type);
  new_title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (action_names[val - 1]) + StringLen (qual_name) + 15));
  StringCpy (new_title, action_names[val - 1]);
  cp = StringChr (new_title, ' ');
  if (cp == NULL) {
    cp = new_title + StringLen (new_title);
  }
  sprintf (cp, title_fmt, qual_name, d->win_num);
  SetTitle (d->form, new_title);
  new_title = MemFree (new_title);
}



static Boolean IsSingleFieldActionVal (Int2 val)
{
  if (val == 1 || val == 2 || val == 7) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsConvertFieldActionVal (Int2 val)
{
  if (val == 3) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean IsDoubleFieldActionVal (Int2 val)
{
  if (val == 4 || val == 5 || val == 6) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean ActionAffectsCDSProduct (OneMacroActionPtr d)
{
  Int2 action_val;
  FieldTypePtr ft = NULL, tmp;
  Boolean show = FALSE;

  if (d == NULL) {
    return FALSE;
  }

  action_val = GetValue (d->action_type);
  switch (action_val) {
    case 1:
    case 2:
      ft = DialogToPointer (d->single_field);
      show = IsFieldTypeCDSProduct (ft);
      ft = FieldTypeFree (ft);
      break;
    case 7:
      ft = DialogToPointer (d->single_field_remove);
      show = IsFieldTypeCDSProduct (ft);
      ft = FieldTypeFree (ft);
      break;
    case 3:
      tmp = DialogToPointer (d->field_pair_convert);
      ft = GetToFieldFromFieldPair (tmp);
      show = IsFieldTypeCDSProduct (ft);
      ft = FieldTypeFree (ft);
      ft = GetFromFieldFromFieldPair (tmp);
      show |= IsFieldTypeCDSProduct (ft);
      ft = FieldTypeFree (ft);
      tmp = FieldPairTypeFree (tmp);
      break;
    case 4:
    case 5:
    case 6:
      tmp = DialogToPointer (d->field_pair);
      ft = GetToFieldFromFieldPair (tmp);
      show = IsFieldTypeCDSProduct (ft);
      ft = FieldTypeFree (ft);
      ft = GetFromFieldFromFieldPair (tmp);
      show |= IsFieldTypeCDSProduct (ft);
      ft = FieldTypeFree (ft);
      tmp = FieldPairTypeFree (tmp);
      break;
  }
  return show;
}


static void ShowOrHideAlsoChangeMrna (OneMacroActionPtr d)
{
  if (d == NULL) {
    return;
  }
  if (ActionAffectsCDSProduct (d)) {
    Show (d->also_change_mrna);
  } else {
    Hide (d->also_change_mrna);
  }
}


static ValNodePtr GetRNATypeFromSingleMacroActionForm (OneMacroActionPtr d, Int2 action_val)
{
  ValNodePtr   rna_type = NULL;
  RnaQualPtr   rq = NULL;
  FieldTypePtr ft = NULL;
  RnaQualPairPtr   rqp = NULL;
  FieldPairTypePtr fpt = NULL;
  
  if (d == NULL) {
    return NULL;
  }
  switch (action_val) {
    case 1:
    case 2:
      ft = DialogToPointer (d->single_field);
      if (ft != NULL && ft->choice == FieldType_rna_field) {
        rq = (RnaQualPtr) ft->data.ptrvalue;
        if (rq != NULL && rq->type != NULL) {
          rna_type = AsnIoMemCopy (rq->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
      }
      ft = FieldTypeFree (ft);
      break;
    case 7:
      ft = DialogToPointer (d->single_field_remove);
      if (ft != NULL && ft->choice == FieldType_rna_field) {
        rq = (RnaQualPtr) ft->data.ptrvalue;
        if (rq != NULL && rq->type != NULL) {
          rna_type = AsnIoMemCopy (rq->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
      }
      ft = FieldTypeFree (ft);
      break;
    case 3:
      fpt = DialogToPointer (d->field_pair_convert);
      if (fpt != NULL && fpt->choice == FieldPairType_rna_field) {
        rqp = (RnaQualPairPtr) fpt->data.ptrvalue;
        if (rqp != NULL && rqp->type != NULL) {
          rna_type = AsnIoMemCopy (rqp->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
      }
      fpt = FieldPairTypeFree (fpt);
      break;
    case 4:
    case 5:
    case 6:
      fpt = DialogToPointer (d->field_pair);
      if (fpt != NULL && fpt->choice == FieldPairType_rna_field) {
        rqp = (RnaQualPairPtr) fpt->data.ptrvalue;
        if (rqp != NULL && rqp->type != NULL) {
          rna_type = AsnIoMemCopy (rqp->type, (AsnReadFunc) RnaFeatTypeAsnRead, (AsnWriteFunc) RnaFeatTypeAsnWrite);
        }
      }
      fpt = FieldPairTypeFree (fpt);
      break;
  }
  return rna_type;
}


static void CopyQualTypeFromPreviousAction (Int2 action, Int2 prev_action, Boolean set_field, OneMacroActionPtr d)
{
  ValNodePtr prev_field = NULL;
  ValNodePtr new_field = NULL;

  if (d == NULL) {
    return;
  }

  /* get qual type (and field, if requested) from previous action */
  if (IsSingleFieldActionVal (prev_action)) {
    if (set_field) {
      if (prev_action == 7) {
        prev_field = DialogToPointer (d->single_field_remove);
      } else {
        prev_field = DialogToPointer (d->single_field);
      }
    }
    if (prev_field == NULL) {
      prev_field = DialogToPointer (d->single_qual_type_dlg);
      if (prev_field != NULL) {
        prev_field->data.ptrvalue = MemFree (prev_field->data.ptrvalue);
      }
    }
  } else {
    if (set_field) {
      if (IsConvertFieldActionVal (prev_action)) {
        prev_field = DialogToPointer (d->field_pair_convert);
      } else {
        prev_field = DialogToPointer (d->field_pair);
      }
    }
    if (prev_field == NULL) {
      prev_field = DialogToPointer (d->pair_qual_type_dlg);
      if (prev_field != NULL) {
        prev_field->data.ptrvalue = MemFree (prev_field->data.ptrvalue);
      }
    }
  }

  /* if necessary, convert from field pair to single field or vice versa */
  if (IsSingleFieldActionVal (prev_action) && !IsSingleFieldActionVal (action)) {
    new_field = BuildFieldPairFromFromField (prev_field);
  } else if (!IsSingleFieldActionVal (prev_action) && IsSingleFieldActionVal (action)) {
    new_field = GetFromFieldFromFieldPair (prev_field);
    if (new_field == NULL && prev_field != NULL) {
      new_field = ValNodeNew (NULL);
      new_field->choice = prev_field->choice;
    }
  } else {
    new_field = prev_field;
  }

  /* apply field to dialogs for new action */
  if (IsSingleFieldActionVal (action)) {
    PointerToDialog (d->single_qual_type_dlg, new_field);
    if (action == 7) {
      PointerToDialog (d->single_field_remove, new_field);
    } else {
      PointerToDialog (d->single_field, new_field);
    }
  } else {
    PointerToDialog (d->pair_qual_type_dlg, new_field);
    if (IsConvertFieldActionVal (action)) {
      PointerToDialog (d->field_pair_convert, new_field);
    } else {
      PointerToDialog (d->field_pair, new_field);
    }
  }

  /* free field */
  if (new_field != prev_field) {
    if (IsSingleFieldActionVal (action)) {
      new_field = FieldTypeFree (new_field);
    } else {
      new_field = FieldPairTypeFree (new_field);
    }
  }

  if (IsSingleFieldActionVal (prev_action)) {
    prev_field = FieldTypeFree (prev_field);
  } else {
    prev_field = FieldPairTypeFree (prev_field);
  }
  
  ShowOrHideAlsoChangeMrna (d);
}


static void MakeQualDialogsVisible (Int2 action, OneMacroActionPtr d)
{
  if (d == NULL) {
    return;
  }

  /* choose whether single or double fields are seen */
  if (IsSingleFieldActionVal (action)) {
    Show (d->single_field_grp);
    Hide (d->field_pair_grp);
    if (action == 7) {
      Show (d->single_field_remove);
      Hide (d->single_field);
    } else {
      Hide (d->single_field_remove);
      Show (d->single_field);
    }
  } else if (IsConvertFieldActionVal (action)) {
    Hide (d->single_field_grp);
    Show (d->field_pair_grp);
    Show (d->field_pair_convert);
    Hide (d->field_pair);
  } else if (IsDoubleFieldActionVal (action)) {
    Hide (d->single_field_grp);
    Show (d->field_pair_grp);
    Hide (d->field_pair_convert);
    Show (d->field_pair);
  }
}


static AECRActionPtr AECRActionFromSingleMacroDialog (OneMacroActionPtr d);

static void PrepopulateField (OneMacroActionPtr d)
{
  Int2 val;
  FieldTypePtr  single_field;
  Boolean       is_nontext;
  AECRActionPtr act;
  ValNodePtr    sample_list = NULL;
  AECRSamplePtr     sample;

  if (d == NULL) {
    return;
  }

  val = GetValue (d->action_type);
  if (val == 1 || val == 2) {
    if (GetAutopopulateStatus (d->action_dlgs[val - 1])) {
      single_field = DialogToPointer (d->single_field);
      if (single_field != NULL) {
        is_nontext = IsFieldTypeNonText (single_field);
        if (!is_nontext) {
          act = AECRActionFromSingleMacroDialog (d);
          sample_list = GetAECRSampleList (act, GetTopSeqEntryForEntityID (d->input_entityID));
          act = AECRActionFree (act);
          sample = GetFieldSampleFromList (sample_list, single_field);
          if (val == 1) {
            if (sample == NULL) {
              SetApplyActionDialogText (d->action_dlgs[0], "");
            } else {
              SetApplyActionDialogText (d->action_dlgs[0], sample->first_value);
            }
          } else {
            if (sample == NULL) {
              SetEditActionDialogText (d->action_dlgs[1], "", NULL);
            } else {
              SetEditActionDialogText (d->action_dlgs[1], sample->first_value, NULL);
            }
          }
          sample_list = AECRSampleListFree (sample_list);
        }
      }
      single_field = FieldTypeFree (single_field);
    }
  }
}


/* When the action changes, the following things must happen:
 * () change which field type and qual type dialogs are visible
 * () copy qual type from previously used qual type dialog to new one
 * () if we are not clearing on change, copy the field from the previous field type dialog
 * () make the correct action dialog visible
 * () adjust apply options (if the action is apply)
 * () prepopulate fields based on sample
 * () adjust the window title
 * () enable the accept button
 */
static void ChangeSingleMacroActionPopup (PopuP p)
{
  OneMacroActionPtr d;
  Int2              val;
  Int4              i;

  d = (OneMacroActionPtr) GetObjectExtra (p);
  if (d == NULL || d->no_callback) return;
  
  val = GetValue (d->action_type);
  if (val < 1 || val > 7) {
    return;
  }

  d->no_callback = TRUE;

  MakeQualDialogsVisible (val, d);
  CopyQualTypeFromPreviousAction (val, d->prev_action_type, !GetStatus (d->clear_on_change), d);

  /* change which dialog is displayed */
  for (i = 0; i < 7; i++) {
    Hide (d->action_dlgs[i]);
  }
  Show (d->action_dlgs[val - 1]);

  /* set action-specific options that depend on field type */
  if (val == 1) {
    ChangeDialogForApplyFieldChoice (d->action_dlgs[val - 1]);
  }

  PrepopulateField (d);

  AdjustSingleMacroActionTitle (d);

  d->prev_action_type = val;
  EnableSingleMacroActionAccept (d);  
  d->no_callback = FALSE;
}


static ValNodePtr GetFieldFromQualTypeDlg (DialoG d)
{
  ValNodePtr field;

  field = DialogToPointer (d);
  if (field != NULL) {
    field->data.ptrvalue = NULL;
  }
  return field;
}


static Int2 GetFeatureTypeFromSingleMacroActionForm (OneMacroActionPtr d)
{
  Int2 val, feat_type = Macro_feature_type_any;
  ValNodePtr field = NULL, field_pair;

  if (d == NULL) {
    return feat_type;
  }

  val = GetValue (d->action_type);
  if (IsSingleFieldActionVal (val)) {
    if (val == 7) {
      field = DialogToPointer (d->single_field_remove);
    } else {
      field = DialogToPointer (d->single_field);
    }
  } else {
    if (IsConvertFieldActionVal (val)) {
      field_pair = DialogToPointer (d->field_pair_convert);
    } else {
      field_pair = DialogToPointer (d->field_pair);
    }
    field = GetFromFieldFromFieldPair (field_pair);
    field_pair = FieldPairTypeFree (field_pair);
  }

  feat_type = FeatureTypeFromFieldType (field);
  field = FieldTypeFree (field);
  return feat_type;
}


static void SetFieldType (OneMacroActionPtr d)
{
  Int2              val;
  FieldTypePtr      single_field = NULL, tmp;
  FieldPairTypePtr  field_pair = NULL;
  Uint1             qual_type;
  ValNodePtr        rna_type = NULL;
  Int2              feat_type = Macro_feature_type_any;

  if (d == NULL) {
    return;
  }
  /* set field type */
  val = GetValue (d->action_type);
  if (IsSingleFieldActionVal (val)) {
    single_field = GetFieldFromQualTypeDlg (d->single_qual_type_dlg);
    if (single_field != NULL) {
      if (val == 7) {
        tmp = GetFieldOfTypeFromFieldType(d->single_field_remove, single_field->choice);
      } else {
        tmp = GetFieldOfTypeFromFieldType(d->single_field, single_field->choice);
      }
      if (tmp != NULL) {
        single_field = FieldTypeFree (single_field);
        single_field = tmp;
      }
    }
    if (val == 7) {
      PointerToDialog (d->single_field_remove, single_field);
    } else {
      PointerToDialog (d->single_field, single_field);
    }
  } else {
    single_field = GetFieldFromQualTypeDlg (d->pair_qual_type_dlg);
    field_pair = BuildFieldPairFromFromField (single_field);
    if (IsConvertFieldActionVal (val)) {
      PointerToDialog (d->field_pair_convert, field_pair);
    } else {
      PointerToDialog (d->field_pair, field_pair);
    }
  }

  qual_type = single_field->choice;
  single_field = FieldTypeFree (single_field);
  field_pair = FieldPairTypeFree (field_pair);

  /* change which constraint dialog is visible, if using "simple" constraints */
  if (qual_type == FieldType_rna_field) {
    rna_type = GetRNATypeFromSingleMacroActionForm (d, val);
  }

  feat_type = GetFeatureTypeFromSingleMacroActionForm (d);

  ChangeComplexConstraintFieldType (d->constraint_dlg, qual_type, rna_type, feat_type);
  rna_type = RnaFeatTypeFree (rna_type);

}

/* when qual type changes, must:
 * () set field type
 * () adjust constraints
 * () prepopulate field
 * () refresh sample
 * () adjust window title
 * () enable accept button
 */
static void ChangeQualType (Pointer data)
{
  OneMacroActionPtr d;

  d = (OneMacroActionPtr) data;
  if (d == NULL || d->no_callback) return;

  d->no_callback = TRUE;
  SetFieldType (d);

  /* set action-specific options that depend on field type */
  if (GetValue (d->action_type) == 1) {
    ChangeDialogForApplyFieldChoice (d->action_dlgs[0]);
  }

  PrepopulateField (d);

  /* refresh sample values */
  RefreshAECRSample (d->accept_btn);

  /* adjust title */
  AdjustSingleMacroActionTitle (d);

  /* enable accept button */
  EnableSingleMacroActionAccept (d);
  d->no_callback = FALSE;
}


/* when field changes, must:
 * () change constraint dialog (for RNA only)
 * () prepopulate field
 * () enable accept button
 */
static void ChangeSingleMacroActionFieldChoice (Pointer data)
{
  OneMacroActionPtr d;
  FieldTypePtr      single_field;
  Uint1             qual_type;
  ValNodePtr        rna_type = NULL;
  Int2              val, feat_type = Macro_feature_type_any;

  d = (OneMacroActionPtr) data;
  if (d == NULL || d->no_callback) return;

  val = GetValue (d->action_type);
  if (IsSingleFieldActionVal (val)) {
    single_field = GetFieldFromQualTypeDlg (d->single_qual_type_dlg);
  } else {
    single_field = GetFieldFromQualTypeDlg (d->pair_qual_type_dlg);
  }

  if (single_field == NULL) {
    return;
  }
  qual_type = single_field->choice;
  single_field = FieldTypeFree (single_field);

  d->no_callback = TRUE;

  if (qual_type == FieldType_rna_field) {
    rna_type = GetRNATypeFromSingleMacroActionForm (d, val);
  }
  feat_type = GetFeatureTypeFromSingleMacroActionForm (d);

  ChangeComplexConstraintFieldType (d->constraint_dlg, qual_type, rna_type, feat_type);
  rna_type = RnaFeatTypeFree (rna_type);

  /* set action-specific options that depend on field type */
  if (val == 1) {
    ChangeDialogForApplyFieldChoice (d->action_dlgs[0]);
  }

  ShowOrHideAlsoChangeMrna (d);

  PrepopulateField (d);

  /* enable accept button */
  EnableSingleMacroActionAccept (d);
  d->no_callback = FALSE;
}


static void AutopopulateSingleMacroDialog (OneMacroActionPtr frm)
{
  Int2 val;
  Boolean status;

  if (frm != NULL) {
    val = GetValue (frm->action_type);
    if (val == 1 || val == 2) {
      status = GetAutopopulateStatus (frm->action_dlgs[val - 1]);
      SetAutopopulateStatus (frm->action_dlgs[0], status);
      SetAutopopulateStatus (frm->action_dlgs[1], status);
      if (status) {
        frm->no_callback = TRUE;
        PrepopulateField (frm);
        EnableSingleMacroActionAccept (frm);
        frm->no_callback = FALSE;
      }
    }
  }
}


static AECRActionPtr AECRActionFromSingleMacroDialog (OneMacroActionPtr d)
{
  Int2 val;
  AECRActionPtr aecr = NULL;

  if (d == NULL) return NULL;

  /* create new action based on contents of dialog */
  val = GetValue (d->action_type);
  switch (val) {
    case 1:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_apply, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
    case 2:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_edit, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
    case 3:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_convert, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
    case 4:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_copy, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
    case 5:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_swap, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
    case 6:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_parse, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
    case 7:
      aecr = AECRActionNew ();
      ValNodeAddPointer (&(aecr->action), ActionChoice_remove, DialogToPointer (d->action_dlgs[val - 1]));
      aecr->constraint = DialogToPointer (d->constraint_dlg);
      break;
  }
  if (aecr != NULL && ActionAffectsCDSProduct (d)) {      
      aecr->also_change_mrna = GetStatus (d->also_change_mrna);
  }
  return aecr;
}


static void SingleMacroAcceptButton (ButtoN b)
{
  OneMacroActionPtr d;
  ValNodePtr        new_action;
  SeqEntryPtr       sep;
  Int4              num_fields, num_features;
  Int2              val;
  Uint2             existing_txt = ExistingTextOption_replace_old;
  AECRActionPtr     action;
  ApplyActionPtr    apply;
  ConvertActionPtr  convert;
  CopyActionPtr     copy;
  AECRParseActionPtr parse;
  ValNodePtr        sample_list = NULL;
  AECRSamplePtr     sample;
  FieldTypePtr      dest_field;

  d = (OneMacroActionPtr) GetObjectExtra (b);
  if (d == NULL) return;

  /* create new action based on contents of dialog */
  val = GetValue (d->action_type);

  new_action = ValNodeNew (NULL);
  switch (val) {
    case 1: /* apply */
      new_action->choice = MacroActionChoice_aecr;
      action = AECRActionFromSingleMacroDialog (d);
      new_action->data.ptrvalue = action;
      apply = action->action->data.ptrvalue;
      sample_list = GetAECRSampleList (action, GetTopSeqEntryForEntityID (d->input_entityID));
      sample = GetFieldSampleFromList (sample_list, apply->field);
      if (sample != NULL && sample->num_found > 0) {
        existing_txt = TwoStepExistingText (sample->num_found, IsFieldTypeNonText(apply->field), AllowFieldMulti (apply->field));
      }
      sample_list = AECRSampleListFree (sample_list);
      if (existing_txt == 0) {
        new_action = MacroActionChoiceFree (new_action);
      } else {
        apply->existing_text = existing_txt;
      }
      break;
    case 3: /* convert */
      new_action->choice = MacroActionChoice_aecr;
      action = AECRActionFromSingleMacroDialog (d);
      new_action->data.ptrvalue = action;
      convert = action->action->data.ptrvalue;
      sample_list = GetAECRSampleList (action, GetTopSeqEntryForEntityID (d->input_entityID));
      dest_field = GetToFieldFromFieldPair (convert->fields);
      sample = GetFieldSampleFromList (sample_list, dest_field);
      if (sample != NULL && sample->num_found > 0) {
        existing_txt = TwoStepExistingText (sample->num_found, IsFieldTypeNonText(dest_field), AllowFieldMulti (dest_field));
      }
      sample_list = AECRSampleListFree (sample_list);
      dest_field = FieldTypeFree (dest_field);
      if (existing_txt == 0) {
        new_action = MacroActionChoiceFree (new_action);
      } else {
        convert->existing_text = existing_txt;
      }
      break;
    case 4: /* copy */
      new_action->choice = MacroActionChoice_aecr;
      action = AECRActionFromSingleMacroDialog (d);
      new_action->data.ptrvalue = action;
      copy = action->action->data.ptrvalue;
      sample_list = GetAECRSampleList (action, GetTopSeqEntryForEntityID (d->input_entityID));
      dest_field = GetToFieldFromFieldPair (copy->fields);
      sample = GetFieldSampleFromList (sample_list, dest_field);
      if (sample != NULL && sample->num_found > 0) {
        existing_txt = TwoStepExistingText (sample->num_found, IsFieldTypeNonText(dest_field), AllowFieldMulti (dest_field));
      }
      sample_list = AECRSampleListFree (sample_list);
      dest_field = FieldTypeFree (dest_field);
      if (existing_txt == 0) {
        new_action = MacroActionChoiceFree (new_action);
      } else {
        copy->existing_text = existing_txt;
      }
      break;
    case 6: /* parse */
      new_action->choice = MacroActionChoice_aecr;
      action = AECRActionFromSingleMacroDialog (d);
      new_action->data.ptrvalue = action;
      parse = action->action->data.ptrvalue;
      sample_list = GetAECRSampleList (action, GetTopSeqEntryForEntityID (d->input_entityID));
      dest_field = GetToFieldFromFieldPair (parse->fields);
      sample = GetFieldSampleFromList (sample_list, dest_field);
      if (sample != NULL && sample->num_found > 0) {
        existing_txt = TwoStepExistingText (sample->num_found, IsFieldTypeNonText(dest_field), AllowFieldMulti (dest_field));
      }
      sample_list = AECRSampleListFree (sample_list);
      dest_field = FieldTypeFree (dest_field);
      if (existing_txt == 0) {
        new_action = MacroActionChoiceFree (new_action);
      } else {
        parse->existing_text = existing_txt;
      }
      break;
    case 2: /* edit */
    case 5: /* swap */
    case 7: /* remove */
      new_action->choice = MacroActionChoice_aecr;
      new_action->data.ptrvalue = AECRActionFromSingleMacroDialog (d);
      break;
  }

  if (new_action == NULL) {
    return;
  }

  /* execute action */
  WatchCursor();
  Update();

  sep = GetTopSeqEntryForEntityID (d->input_entityID);
  num_fields = 0;
  num_features = 0;
  ApplyMacroToSeqEntryEx (sep, new_action, NULL, Sequin_GlobalAlign2Seq);
  ObjMgrSetDirtyFlag (d->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, d->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();   
  new_action = MacroActionChoiceFree (new_action);

  /* done if leave_dlg_up not checked */
  if (!GetStatus (d->leave_dlg_up)) {
    RemoveMacroActionWindow ((WindoW)(d->form));
  }
}


static void SingleMacroCancelButton (ButtoN b)
{
  OneMacroActionPtr d;

  d = (OneMacroActionPtr) GetObjectExtra (b);
  if (d == NULL) return;

  RemoveMacroActionWindow ((WindoW)(d->form));
}


static void ActionToSingleMacroActionForm (ForM f, Pointer data)
{
  OneMacroActionPtr     frm;
  ValNodePtr            action;
  AECRActionPtr         aecr;
  Int2                  val, i;


  frm = (OneMacroActionPtr) GetObjectExtra (f);
  if (frm == NULL) return;

  frm->no_callback = TRUE;
  action = (ValNodePtr) data;

  if (action == NULL) {
    SetValue (frm->action_type, 1);
    PointerToDialog (frm->action_dlgs[0], NULL);
    val = 1;
  } else if (action->choice == MacroActionChoice_aecr) {
    aecr = (AECRActionPtr) action->data.ptrvalue;
    if (aecr == NULL || aecr->action == NULL) {
      SetValue (frm->action_type, 1);
      val = 1;
    } else {
      switch (aecr->action->choice) {
        case ActionChoice_apply:
          SetValue (frm->action_type, 1);
          break;
        case ActionChoice_edit:
          SetValue (frm->action_type, 2);
          break;
        case ActionChoice_convert:
          SetValue (frm->action_type, 3);
          break;
        case ActionChoice_copy:
          SetValue (frm->action_type, 4);
          break;
        case ActionChoice_swap:
          SetValue (frm->action_type, 5);
          break;
        case ActionChoice_parse:
          SetValue (frm->action_type, 6);
          break;
        case ActionChoice_remove:
          SetValue (frm->action_type, 7);
          break;
        default:
          SetValue (frm->action_type, 1);
          break;
      }
      val = GetValue (frm->action_type);
      PointerToDialog (frm->action_dlgs[val - 1], aecr->action->data.ptrvalue);
    }
  }

  /* change which dialog is displayed */
  for (i = 0; i < 7; i++) {
    Hide (frm->action_dlgs[i]);
  }
  Show (frm->action_dlgs[val - 1]);

  /* set field type */
  SetFieldType(frm);

  MakeQualDialogsVisible (val, frm);

  /* set action-specific options that depend on field type */
  if (val == 1) {
    ChangeDialogForApplyFieldChoice (frm->action_dlgs[0]);
  }

  PrepopulateField (frm);

  AdjustSingleMacroActionTitle (frm);

  RefreshAECRSample (frm->accept_btn);
  EnableSingleMacroActionAccept (frm);
  frm->no_callback = FALSE;
}


static void ClearSingleMacroAction (ButtoN b)
{
  OneMacroActionPtr frm;
  Int2              val;
  ValNodePtr        qual_type;
  DialoG            qual_dlg;

  frm = (OneMacroActionPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  val = GetValue (frm->action_type);

  /* get field type, so that we can put it back after the clear */
  if (IsSingleFieldActionVal (val)) {
    qual_dlg = frm->single_qual_type_dlg;
  } else {
    qual_dlg = frm->pair_qual_type_dlg;
  }
  qual_type = DialogToPointer (qual_dlg);

  /* clear action dialog */
  if (val > 0 && val < 8) {
    PointerToDialog (frm->action_dlgs[val - 1], NULL);
  }

  /* reapply field type */
  PointerToDialog (qual_dlg, qual_type);
  qual_type = ValNodeFreeData (qual_type);
  SetFieldType (frm);

  /* also clear constraints */
  PointerToDialog (frm->constraint_dlg, NULL);
}


static void ClearTextSingleMacroAction (ButtoN b)
{
  OneMacroActionPtr frm;
  Int2              val;

  frm = (OneMacroActionPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  val = GetValue (frm->action_type);

  switch (val) {
    case 1:
      SetApplyActionDialogText (frm->action_dlgs[0], "");
      break;
    case 2:
      SetEditActionDialogText (frm->action_dlgs[1], "", "");
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      break;
    case 7:
      break;
  }
  /* also clear text in simple constraint if shown */
  ClearComplexConstraintDialogText (frm->constraint_dlg);
  EnableSingleMacroActionAccept (frm);
}


static Int2 LIBCALLBACK SingleMacroActionMsgFunc (OMMsgStructPtr ommsp)

{
  ObjMgrDataPtr   omdp;
  OMUserDataPtr   omudp;
  OneMacroActionPtr frm;

  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  frm = (OneMacroActionPtr) omudp->userdata.ptrvalue;
  if (frm == NULL) return OM_MSG_RET_ERROR;
  switch (ommsp->message) {
    case OM_MSG_DEL:
      omdp = ObjMgrGetData (ommsp->entityID);
      if (omdp != NULL) {
        if (ObjMgrWholeEntity (omdp, ommsp->itemID, ommsp->itemtype)) {
          /* clear text here */
          SendMessageToDialog (frm->action_dlgs[1], NUM_VIB_MSG + 1);
          return OM_MSG_RET_OK;
        }
      }
      break;
    default :
      break;
  }
  return OM_MSG_RET_OK;
}


static void SingleMacroActionFormMessage (ForM f, Int2 mssg)
{
  OneMacroActionPtr frm;

  frm = (OneMacroActionPtr) GetObjectExtra (f);
  if (frm != NULL) {
    switch (mssg) {
      case VIB_MSG_ACCEPT :
        SingleMacroAcceptButton (frm->accept_btn);
        break;
      case VIB_MSG_CLOSE:
        RemoveMacroActionWindow ((WindoW)(frm->form));
        break;
      case VIB_MSG_ENTER :
        Select (f);
        Select (frm->single_field);
        break;
      default :
        break;
    }
  }
}


#ifndef WIN_MAC
static void CreateSingleMacroActionFormMenus (WindoW w)

{
  BaseFormPtr   bfp;
  MenU          m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    FormCommandItem (m, "Accept", bfp, VIB_MSG_ACCEPT);
    FormCommandItem (m, "Cancel", bfp, VIB_MSG_CLOSE);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif


/* for executing macro actions without creating a macro script */
NLM_EXTERN ForM 
SingleMacroAction (Uint2 entityID, Boolean indexer_version)
{
  WindoW                w;
  GrouP                 h, g5, c, field_grp, k, c1, not_constraint, split_grp;
  ButtoN                b, show_sample_btn;
  OneMacroActionPtr     frm;
  Int4                  i;
  ValNodePtr            val_list = NULL;
  Char                  tmpval[30];
  OMUserDataPtr         omudp;

  frm = (OneMacroActionPtr) MemNew (sizeof (OneMacroActionData));

  w = FixedWindow(-20, -13, -10, -10, "Edit Action", RemoveMacroActionWindow);
  SetObjectExtra (w, frm, NULL);

  frm->form = (ForM) w;
  frm->input_entityID = entityID;
  frm->toform = ActionToSingleMacroActionForm;
  frm->formmessage = SingleMacroActionFormMessage;
  frm->indexer_version = indexer_version;

  frm->win_num = AddMacroActionWindow (w);

  /* register to receive update messages */
  frm->userkey = OMGetNextUserKey ();
  frm->procid = 0;
  frm->proctype = OMPROC_EDIT;
  omudp = ObjMgrAddUserData (0, 0, frm->proctype, frm->userkey);
  if (omudp != NULL) {
    omudp->userdata.ptrvalue = (Pointer) frm;
    omudp->messagefunc = SingleMacroActionMsgFunc;
  }

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  frm->action_type = PopupList (h, TRUE, ChangeSingleMacroActionPopup);
  SetObjectExtra (frm->action_type, frm, NULL);
  for (i = 0; i < 7; i++) {
    PopupItem (frm->action_type, action_names[i]);
  }
  SetValue (frm->action_type, 1);
  frm->prev_action_type = 1;
  frm->clear_on_change = CheckBox (h, "Clear when changing actions", NULL);

  split_grp = h;

  not_constraint = HiddenGroup (split_grp, -1, 0, NULL);
  SetGroupSpacing (not_constraint, 10, 10);
  field_grp = HiddenGroup (not_constraint, 0, 0, NULL);
  frm->single_field_grp = HiddenGroup (field_grp, -1, 0, NULL);
 
  ValNodeAddPointer (&val_list, FieldType_source_qual, StringSave ("Source Qual"));
  ValNodeAddPointer (&val_list, FieldType_feature_field, StringSave ("Feature Qual"));
  ValNodeAddPointer (&val_list, FieldType_cds_gene_prot, StringSave ("CDS-Gene-Prot Qual"));
  ValNodeAddPointer (&val_list, FieldType_rna_field, StringSave ("RNA Qual"));
  ValNodeAddPointer (&val_list, FieldType_molinfo_field, StringSave ("MolInfo Qual"));
  ValNodeAddPointer (&val_list, FieldType_pub, StringSave ("Pub Field"));
  ValNodeAddPointer (&val_list, FieldType_struc_comment_field, StringSave ("Structured Comment Field"));
  ValNodeAddPointer (&val_list, FieldType_dblink, StringSave ("DBLink Field"));
  ValNodeAddPointer (&val_list, FieldType_misc, StringSave ("Misc"));

  frm->single_qual_type_dlg = ValNodeSelectionDialog (frm->single_field_grp, val_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type", 
                                ChangeQualType, frm, FALSE);
  val_list = NULL;
  k = HiddenGroup (frm->single_field_grp, 0, 0, NULL);
  frm->single_field = FieldTypeDialog (k, FALSE, FALSE, ChangeSingleMacroActionFieldChoice, frm);
  frm->single_field_remove = FieldTypeDialog (k, FALSE, TRUE, ChangeSingleMacroActionFieldChoice, frm);
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->single_field, (HANDLE) frm->single_field_remove, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->single_qual_type_dlg, (HANDLE) k, NULL);

  frm->field_pair_grp = HiddenGroup (field_grp, -1, 0, NULL);
  ValNodeAddPointer (&val_list, FieldType_source_qual, StringSave ("Source Qual"));
  ValNodeAddPointer (&val_list, FieldType_feature_field, StringSave ("Feature Qual"));
  ValNodeAddPointer (&val_list, FieldType_cds_gene_prot, StringSave ("CDS-Gene-Prot Qual"));
  ValNodeAddPointer (&val_list, FieldType_rna_field, StringSave ("RNA Qual"));
  ValNodeAddPointer (&val_list, FieldType_molinfo_field, StringSave ("MolInfo Qual"));
  ValNodeAddPointer (&val_list, FieldType_dblink, StringSave ("DBLink Field"));
  frm->pair_qual_type_dlg = ValNodeSelectionDialog (frm->field_pair_grp, val_list, TALL_SELECTION_LIST, ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "field type", 
                                ChangeQualType, frm, FALSE);
  val_list = NULL;
  
  k = HiddenGroup (frm->field_pair_grp, 0, 0, NULL);
  frm->field_pair = FieldPairTypeDialog (k, FALSE, ChangeSingleMacroActionFieldChoice, frm);
  frm->field_pair_convert = FieldPairTypeDialog (k, TRUE, ChangeSingleMacroActionFieldChoice, frm);
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->field_pair_convert, (HANDLE) frm->field_pair, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->pair_qual_type_dlg, (HANDLE) k, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) frm->single_field_grp, (HANDLE) frm->field_pair_grp, NULL);

  Hide (frm->single_field_grp);
  Hide (frm->field_pair_grp);

  g5 = HiddenGroup (not_constraint, 0, 0, NULL);
  /* temporarily shut off callbacks */
  frm->no_callback = TRUE;
  frm->action_dlgs[0] = ApplyActionDialog (g5, frm->indexer_version, FALSE, (Nlm_ChangeNotifyProc) AutopopulateSingleMacroDialog, frm, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[0], frm->single_qual_type_dlg, frm->single_field);
  frm->action_dlgs[1] = EditActionDialog (g5, frm->indexer_version, (Nlm_ChangeNotifyProc) AutopopulateSingleMacroDialog, frm, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[1], frm->single_qual_type_dlg, frm->single_field);
  frm->action_dlgs[2] = ConvertActionDialog (g5, frm->indexer_version, FALSE, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[2], frm->pair_qual_type_dlg, frm->field_pair_convert);
  frm->action_dlgs[3] = CopyActionDialog (g5, frm->indexer_version, FALSE, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[3], frm->pair_qual_type_dlg, frm->field_pair);
  frm->action_dlgs[4] = SwapActionDialog (g5, frm->indexer_version, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[4], frm->pair_qual_type_dlg, frm->field_pair);
  frm->action_dlgs[5] = AECRParseActionDialog (g5, frm->indexer_version, FALSE, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[5], frm->pair_qual_type_dlg, frm->field_pair);
  frm->action_dlgs[6] = RemoveActionDialog (g5, frm->indexer_version, EnableSingleMacroActionAccept, frm, EnableSingleMacroActionAccept, frm);
  SetAECRActionDlgFieldTypeDialogs (frm->action_dlgs[6], frm->single_qual_type_dlg, frm->single_field_remove);
  AlignObjects (ALIGN_CENTER, (HANDLE) frm->action_dlgs[0],
                              (HANDLE) frm->action_dlgs[1],
                              (HANDLE) frm->action_dlgs[2],
                              (HANDLE) frm->action_dlgs[3],
                              (HANDLE) frm->action_dlgs[4],
                              (HANDLE) frm->action_dlgs[5],
                              (HANDLE) frm->action_dlgs[6],
                              NULL);

  /* set default autopopulate value */
  if (GetAppParam ("SEQUINCUSTOM", "BATCHDIALOG", "AUTOPOPULATE", NULL, tmpval, sizeof (tmpval))
        && StringICmp (tmpval, "TRUE") == 0) {
    SetAutopopulateStatus (frm->action_dlgs[0], TRUE);
    SetAutopopulateStatus (frm->action_dlgs[1], TRUE);
  } else {
    SetAutopopulateStatus (frm->action_dlgs[0], FALSE);
    SetAutopopulateStatus (frm->action_dlgs[1], FALSE);
  }

  frm->also_change_mrna = CheckBox (not_constraint, "Make mRNA product match CDS protein name", NULL);
  Hide (frm->also_change_mrna);

  AlignObjects (ALIGN_CENTER, (HANDLE) field_grp, (HANDLE) g5, (HANDLE) frm->also_change_mrna, NULL);

  frm->constraint_dlg = ComplexConstraintDialog (split_grp, EnableSingleMacroActionAccept, frm);

  c1 = HiddenGroup (h, 4, 0, NULL);
  show_sample_btn = PushButton (c1, "Show Sample", ShowMacroSampleWindow);
  SetObjectExtra (show_sample_btn, frm, NULL);
  b = PushButton (c1, "Clear", ClearSingleMacroAction);
  SetObjectExtra (b, frm, NULL);
  b = PushButton (c1, "Clear Text", ClearTextSingleMacroAction);
  SetObjectExtra (b, frm, NULL);
  frm->leave_dlg_up = CheckBox (c1, "Leave dialog up", NULL);


  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  frm->accept_btn = PushButton (c, "Accept", SingleMacroAcceptButton);
  SetObjectExtra (frm->accept_btn, frm, NULL);
  b = PushButton (c, "Cancel", SingleMacroCancelButton);
  SetObjectExtra (b, frm, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) frm->action_type,
                              (HANDLE) frm->clear_on_change, 
                              (HANDLE) not_constraint,
                              (HANDLE) frm->constraint_dlg,
                              (HANDLE) c1,
                              (HANDLE) c,
                              NULL);

#ifndef WIN_MAC
  CreateSingleMacroActionFormMenus (w);
#endif

  /* turn callbacks back on */
  frm->no_callback = FALSE;
  ChangeSingleMacroActionPopup (frm->action_type);

  return (ForM) w;
}



NLM_EXTERN void SingleAECRMacroAction (Uint2 entityID, Boolean indexer_version, Uint1 AECR_action_type, Uint1 AECR_qual_type)
{
  ForM f;
  ValNodePtr action;
  
  f = SingleMacroAction (entityID, indexer_version);
  action = ValNodeNew (NULL);
  action->choice = MacroActionChoice_aecr;
  action->data.ptrvalue = BuildDefaultAECRAction (AECR_action_type, AECR_qual_type);
  PointerToForm (f, action);
  action = MacroActionChoiceFree (action);
    
  Show (f);
  Select (f);
  SendMessageToForm (f, VIB_MSG_ENTER);
}


NLM_EXTERN void MacroApplyKeyword (Uint2 entityID, Boolean indexer_version)
{
  ForM f;
  ValNodePtr action;
  AECRActionPtr aecr;
  ApplyActionPtr apply;
  
  f = SingleMacroAction (entityID, indexer_version);
  action = ValNodeNew (NULL);
  action->choice = MacroActionChoice_aecr;

  aecr = AECRActionNew ();
  aecr->action = ValNodeNew (NULL);
  aecr->action->choice = ActionChoice_apply;
  apply = ApplyActionNew ();
  apply->field = ValNodeNew (NULL);
  apply->field->choice = FieldType_misc;
  apply->field->data.intvalue = Misc_field_keyword;
  apply->value = StringSave ("");
  apply->existing_text = ExistingTextOption_add_qual;
  aecr->action->data.ptrvalue = apply;
  action->data.ptrvalue = aecr;
  PointerToForm (f, action);
  action = MacroActionChoiceFree (action);
    
  Show (f);
  Select (f);
}


typedef struct singleparseactionfrm {
  FORM_MESSAGE_BLOCK
  DialoG dlg;
  ButtoN leave_dlg_up;
  ButtoN accept_btn;
} SingleParseActionFrmData, PNTR SingleParseActionFrmPtr;


static void ChangeSingleParseAction (Pointer data)
{
  SingleParseActionFrmPtr frm;
  ValNodePtr err_list;

  frm = (SingleParseActionFrmPtr) data;
  if (frm == NULL) {
    return;
  }

  /* enable or disable accept button */
  err_list = TestDialog (frm->dlg);
  if (err_list == NULL) {
    Enable (frm->accept_btn);
  } else {
    Disable (frm->accept_btn);
  }
  err_list = ValNodeFree (err_list);
}


static void ClearSingleParseAction (ButtoN b)
{
  SingleParseActionFrmPtr frm;

  frm = (SingleParseActionFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  PointerToDialog (frm->dlg, NULL);
  ChangeSingleParseAction (frm);
}


static void ClearTextSingleParseAction (ButtoN b)
{
  SingleParseActionFrmPtr frm;
  ParseActionPtr parse;

  frm = (SingleParseActionFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  parse = DialogToPointer (frm->dlg);

  if (parse != NULL && parse->portion != NULL) {
    parse->portion->left_marker = TextMarkerFree (parse->portion->left_marker);
    parse->portion->right_marker = TextMarkerFree (parse->portion->right_marker);
    PointerToDialog (frm->dlg, parse);
    parse = ParseActionFree (parse);
    ChangeSingleParseAction (frm);
  }
}


static void ActionToSingleParseActionForm (ForM f, Pointer data)
{
  SingleParseActionFrmPtr frm;
  ParseActionPtr          parse;


  frm = (SingleParseActionFrmPtr) GetObjectExtra (f);
  if (frm == NULL) return;

  parse = (ParseActionPtr) data;

  PointerToDialog (frm->dlg, parse);
  ChangeSingleParseAction (frm);
}


static Pointer FormToParseAction (ForM f)
{
  SingleParseActionFrmPtr frm;
  ParseActionPtr          parse;


  frm = (SingleParseActionFrmPtr) GetObjectExtra (f);
  if (frm == NULL) return NULL;

  parse = DialogToPointer (frm->dlg);
  return parse;
}



static void SingleParseAcceptButton (ButtoN b)
{
  SingleParseActionFrmPtr frm;
  ParseActionPtr          parse;
  AECRSamplePtr           sample;
  SeqEntryPtr             sep;
  ValNodePtr              vnp;

  frm = (SingleParseActionFrmPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  parse = DialogToPointer (frm->dlg);
  if (parse == NULL) {
    return;
  }

  /* check for existing text, then set value in parse*/
  sample = GetExistingTextForParseAction (parse, sep);
  if (sample != NULL && sample->num_found > 0) {
    parse->existing_text = TwoStepExistingText (sample->num_found, FALSE, FALSE);
  } else {
    parse->existing_text = ExistingTextOption_replace_old;
  }
  sample = AECRSampleFree (sample);
  if (parse->existing_text == 0) {
    parse = ParseActionFree (parse);
    return;
  }

  vnp = ValNodeNew (NULL);
  vnp->choice = MacroActionChoice_parse;
  vnp->data.ptrvalue = parse;

  ApplyMacroToSeqEntryEx (sep, vnp, NULL, Sequin_GlobalAlign2Seq);

  vnp = MacroActionChoiceFree (vnp);

  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);
  if (!GetStatus (frm->leave_dlg_up)) {
    Remove (frm->form);
  }
  Update ();
}


/* todo - later allow different starting setups */
NLM_EXTERN ForM SingleParseAction (Uint2 entityID)
{
  WindoW                w;
  GrouP                 h, c, c1;
  ButtoN                b;
  SingleParseActionFrmPtr frm;
  SeqEntryPtr             sep;

  frm = (SingleParseActionFrmPtr) MemNew (sizeof (SingleParseActionFrmData));

  w = FixedWindow(-20, -13, -10, -10, "Parse Text", StdCloseWindowProc);
  SetObjectExtra (w, frm, NULL);

  frm->form = (ForM) w;
  frm->input_entityID = entityID;
  frm->toform = ActionToSingleParseActionForm;
  frm->fromform = FormToParseAction;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  sep = GetTopSeqEntryForEntityID (entityID);

  frm->dlg = ParseActionDialogEx (h, FALSE, FALSE, ChangeSingleParseAction, frm, sep);

  c1 = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c1, "Clear", ClearSingleParseAction);
  SetObjectExtra (b, frm, NULL);
  b = PushButton (c1, "Clear Text", ClearTextSingleParseAction);
  SetObjectExtra (b, frm, NULL);
  frm->leave_dlg_up = CheckBox (c1, "Leave dialog up", NULL);


  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  frm->accept_btn = PushButton (c, "Accept", SingleParseAcceptButton);
  SetObjectExtra (frm->accept_btn, frm, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) frm->dlg,
                              (HANDLE) c1,
                              (HANDLE) c,
                              NULL);

  return (ForM) w;

}


/* Functions for summarizing macro actions for display */

static CharPtr SummarizeApplyFeatureAction (ApplyFeatureActionPtr a)
{
  CharPtr    label = NULL;
  CharPtr    str;
  CharPtr    fmt = "Apply %s";

  if (a == NULL) {
    str = StringSave ("No action");
  } else {
    label = GetFeatureNameFromFeatureType (a->type);
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (label)));
    sprintf (str, fmt, label);
  }
  return str;
}


static CharPtr SummarizeRemoveFeatureAction (RemoveFeatureActionPtr a)
{
  CharPtr    label = NULL;
  CharPtr    constraint, str;
  CharPtr    fmt = "Remove %s";
  CharPtr    constraint_fmt = "Remove %s %s";

  if (a == NULL) {
    str = StringSave ("No action");
  } else {
    label = GetFeatureNameFromFeatureType (a->type);
    constraint = SummarizeConstraintSet (a->constraint);
    if (constraint == NULL) {
      str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (label)));
      sprintf (str, fmt, label);
    } else {
      str = (CharPtr) MemNew (sizeof (Char) * (StringLen (constraint_fmt) + StringLen (label) + StringLen (constraint)));
      sprintf (str, constraint_fmt, label, constraint);
      constraint = MemFree (constraint);
    }
  }

  return str;
}


static CharPtr SummarizeConvertSourceOptions (ValNodePtr vnp)
{
  ConvertFromCDSOptionsPtr options;
  CharPtr fmt = "(%sremove overlapping mRNA, %sremove overlapping gene, %sremove transcript ID)";
  CharPtr str;

  if (vnp == NULL || vnp->choice != ConvertFeatureSrcOptions_cds || vnp->data.ptrvalue == NULL) {
    return NULL;
  }

  options = (ConvertFromCDSOptionsPtr) vnp->data.ptrvalue;

  str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 21));
  sprintf (str, fmt, options->remove_mRNA ? "" : "do not ",
                     options->remove_gene ? "" : "do not ",
                     options->remove_transcript_id ? "" : "do not ");
  return str;
}


static CharPtr SummarizeConvertDestOptions (ValNodePtr vnp)
{
  RegionTypePtr r;
  CharPtr str = NULL;

  if (vnp == NULL) return NULL;

  switch (vnp->choice) {
    case ConvertFeatureDstOptions_bond:
      str = StringSave (GetMacroBondTypeName(vnp->data.intvalue));
      break;
    case ConvertFeatureDstOptions_site:
      str = StringSave (GetMacroSiteTypeName(vnp->data.intvalue));
      break;
    case ConvertFeatureDstOptions_region:
      r = (RegionTypePtr) vnp->data.ptrvalue;
      if (r != NULL) {
        if (r->create_nucleotide) {
          str = StringSave ("on nucleotide sequence");
        } else {
          str = StringSave ("on protein sequence");
        }
      }
      break;
  }
  return str;
}


static CharPtr SummarizeConvertFeatureAction (ConvertFeatureActionPtr a)
{
  CharPtr str = NULL, from_label, to_label, constraint, src_options, dst_options;
  CharPtr fmt = "Convert %s to %s";
  CharPtr keep_orig = ", keep original feature";
  CharPtr remove_orig = ", remove original feature";
  Int4    len;

  if (a == NULL) {
    str = StringSave ("No action");
  } else {
    from_label = GetFeatureNameFromFeatureType (a->type_from);
    to_label = GetFeatureNameFromFeatureType (a->type_to);
    src_options = SummarizeConvertSourceOptions (a->src_options);
    dst_options = SummarizeConvertDestOptions (a->dst_options);
    constraint = SummarizeConstraintSet (a->src_feat_constraint);
    len = StringLen (fmt) + StringLen (from_label) + StringLen (to_label);
    if (src_options != NULL) {
      len += StringLen (src_options) + 3;
    }
    if (dst_options != NULL) {
      len += StringLen (dst_options) + 1;
    }
    if (constraint != NULL) {
      len += StringLen (constraint) + 1;
    }
    if (a->leave_original) {
      len += StringLen (keep_orig);
    } else {
      len += StringLen (remove_orig);
    }
    str = (CharPtr) MemNew (sizeof (Char) * len);
    sprintf (str, fmt, from_label, to_label);
    if (dst_options != NULL) {
      StringCat (str, " ");
      StringCat (str, dst_options);
      dst_options = MemFree (dst_options);
    }
    if (src_options != NULL) {
      StringCat (str, ", ");
      StringCat (str, src_options);
      src_options = MemFree (src_options);
    }
    if (constraint != NULL) {
      StringCat (str, " ");
      StringCat (str, constraint);
      constraint = MemFree (constraint);
    }
    if (a->leave_original) {
      StringCat (str, keep_orig);
    } else {
      StringCat (str, remove_orig);
    }
  }
  return str;
}


static CharPtr SummarizeEditLocationStrand (EditLocationStrandPtr strand)
{
  CharPtr from_label = NULL, to_label = NULL;
  CharPtr fmt = "Convert %s strand to %s";
  CharPtr str = NULL;

  if (strand == NULL) return NULL;

  switch (strand->strand_from) {
    case Feature_location_strand_from_any:
      from_label = "any";
      break;
    case Feature_location_strand_from_plus:
      from_label = "plus";
      break;
    case Feature_location_strand_from_minus:
      from_label = "minus";
      break;
    case Feature_location_strand_from_unknown:
      from_label = "unknown";
      break;
    case Feature_location_strand_from_both:
      from_label = "both";
      break;
  }

  switch (strand->strand_to) {
    case Feature_location_strand_to_plus:
      to_label = "plus";
      break;
    case Feature_location_strand_to_minus:
      to_label = "minus";
      break;
    case Feature_location_strand_to_unknown:
      to_label = "unknown";
      break;
    case Feature_location_strand_to_both:
      to_label = "both";
      break;
    case Feature_location_strand_to_reverse:
      to_label = "reverse";
      break;
  }

  if (from_label != NULL && to_label != NULL) {
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (from_label) + StringLen (to_label)));
    sprintf (str, fmt, from_label, to_label);
  }
  return str;
}


static CharPtr SummarizePartial5SetAction (Partial5SetActionPtr a)
{
  CharPtr str = NULL;
  CharPtr constraint = NULL, extend = NULL;
  CharPtr fmt = "Set 5' partial%s%s";

  if (a == NULL) return NULL;

  switch (a->constraint) {
    case Partial_5_set_constraint_all:
      constraint = "";
      break;
    case Partial_5_set_constraint_at_end:
      constraint = " when 5' end of location is at end of sequence";
      break;
    case Partial_5_set_constraint_bad_start:
      constraint = " when coding region has no start codon";
      break;
    case Partial_5_set_constraint_frame_not_one:
      constraint = " when coding region frame > 1";
      break;
  }
  if (a->extend) {
    extend = ", extend 5' end of feature to end of sequence";
  } else {
    extend = "";
  }
  if (constraint != NULL) {
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt)
                                             + StringLen (constraint) 
                                             + StringLen (extend)));
    sprintf (str, fmt, constraint, extend);
  }
  return str;
}


static CharPtr SummarizePartial5ClearAction (Int4 a)
{
  CharPtr str = NULL;

  switch (a) {
    case Partial_5_clear_constraint_all:
      str = StringSave ("Clear 5' partial");
      break;
    case Partial_5_clear_constraint_not_at_end:
      str = StringSave ("Clear 5' partial when 5' end of feature is not at end of sequence");
      break;
    case Partial_5_clear_constraint_good_start:
      str = StringSave ("Clear 5' partial when coding region has start codon");
      break;
  }
  return str;
}


static CharPtr SummarizePartial3SetAction (Partial3SetActionPtr a)
{
  CharPtr str = NULL;
  CharPtr constraint = NULL, extend = NULL;
  CharPtr fmt = "Set 3' partial%s%s";

  if (a == NULL) return NULL;

  switch (a->constraint) {
    case Partial_3_set_constraint_all:
      constraint = "";
      break;
    case Partial_3_set_constraint_at_end:
      constraint = " when 3' end of location is at end of sequence";
      break;
    case Partial_3_set_constraint_bad_end:
      constraint = " when coding region has no stop codon";
      break;
  }
  if (a->extend) {
    extend = ", extend 3' end of feature to end of sequence";
  } else {
    extend = "";
  }
  if (constraint != NULL) {
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt)
                                             + StringLen (constraint) 
                                             + StringLen (extend)));
    sprintf (str, fmt, constraint, extend);
  }
  return str;
}


static CharPtr SummarizePartial3ClearAction (Int4 a)
{
  CharPtr str = NULL;

  switch (a) {
    case Partial_3_clear_constraint_all:
      str = StringSave ("Clear 3' partial");
      break;
    case Partial_3_clear_constraint_not_at_end:
      str = StringSave ("Clear 3' partial when 3' end of feature is not at end of sequence");
      break;
    case Partial_3_clear_constraint_good_end:
      str = StringSave ("Clear 3' partial when coding region has stop codon");
      break;
  }
  return str;
}


static CharPtr SummarizePartialBothSetAction (PartialBothSetActionPtr a)
{
  CharPtr str = NULL;
  CharPtr constraint = NULL, extend = NULL;
  CharPtr fmt = "Set both ends partial%s%s";

  if (a == NULL) return NULL;

  switch (a->constraint) {
    case Partial_5_set_constraint_all:
      constraint = "";
      break;
    case Partial_5_set_constraint_at_end:
      constraint = " when both ends of location are at end of sequence";
      break;
  }
  if (a->extend) {
    extend = ", extend both ends of feature to end of sequence";
  } else {
    extend = "";
  }
  if (constraint != NULL) {
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt)
                                             + StringLen (constraint) 
                                             + StringLen (extend)));
    sprintf (str, fmt, constraint, extend);
  }
  return str;
}


static CharPtr SummarizePartialBothClearAction (Int4 a)
{
  CharPtr str = NULL;

  switch (a) {
    case Partial_both_clear_constraint_all:
      str = StringSave ("Clear both ends partial");
      break;
    case Partial_3_clear_constraint_not_at_end:
      str = StringSave ("Clear both ends partial when both ends of feature are not at end of sequence");
      break;
  }
  return str;
}


static CharPtr SummarizeConvertLoc (Int4 a)
{
  CharPtr str = NULL;

  switch (a) {
    case Convert_location_type_join:
      str = StringSave ("Convert location to join");
      break;
    case Convert_location_type_order:
      str = StringSave ("Convert location to order");
      break;
    case Convert_location_type_merge:
      str = StringSave ("Convert location to single interval");
      break;
  }
  return str;
}


static CharPtr SummarizeEditFeatureLocationAction (EditFeatureLocationActionPtr a)
{
  CharPtr str = NULL, action_label = NULL, constraint, feature;
  CharPtr fmt = "%s for %s features";
  CharPtr constraint_fmt = "%s for %s features %s";

  if (a == NULL || a->action == NULL) {
    str = StringSave ("No action");
  } else {
    
    switch (a->action->choice) {
      case LocationEditType_strand:
        action_label = SummarizeEditLocationStrand (a->action->data.ptrvalue);
        break;
      case LocationEditType_set_5_partial:
        action_label = SummarizePartial5SetAction (a->action->data.ptrvalue);
        break;
      case LocationEditType_clear_5_partial:
        action_label = SummarizePartial5ClearAction (a->action->data.intvalue);
        break;
      case LocationEditType_set_3_partial:
        action_label = SummarizePartial3SetAction (a->action->data.ptrvalue);
        break;
      case LocationEditType_clear_3_partial:
        action_label = SummarizePartial3ClearAction (a->action->data.intvalue);
        break;
      case LocationEditType_set_both_partial:
        action_label = SummarizePartialBothSetAction (a->action->data.ptrvalue);
        break;
      case LocationEditType_clear_both_partial:
        action_label = SummarizePartialBothClearAction (a->action->data.intvalue);
        break;
      case LocationEditType_convert:
        action_label = SummarizeConvertLoc (a->action->data.intvalue);
        break;
      case LocationEditType_extend_5:
        action_label = StringSave ("Extend 5' end of feature to end of sequence");
        break;
      case LocationEditType_extend_3:
        action_label = StringSave ("Extend 3' end of feature to end of sequence");
        break;
    }
    if (action_label == NULL) {
      str = StringSave ("Invalid action");
    } else {
      feature = GetFeatureNameFromFeatureType (a->type);
      constraint = SummarizeConstraintSet (a->constraint);
      if (constraint == NULL) {
        str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (action_label) + StringLen (feature)));
        sprintf (str, fmt, action_label, feature);
      } else {
        str = (CharPtr) MemNew (sizeof (Char) * (StringLen (constraint_fmt) + StringLen (action_label) + StringLen (feature) + StringLen (constraint)));
        sprintf (str, constraint_fmt, action_label, feature, constraint);
        constraint = MemFree (constraint);
      }
    }
  } 
  return str;
}


static CharPtr s_Suppression[] = { NULL, "suppressing", "non-suppressing" };
static CharPtr s_Necessary[] = { NULL, "necessary", "unnecessary" };

static CharPtr SummarizeRemoveXref (RemoveXrefsActionPtr a)
{
  CharPtr str = NULL, label, constraint;
  GeneXrefTypePtr g;
  CharPtr fmt = "Remove %s%s%s%sgene xrefs from %s features";
  CharPtr suppression, necessary;
  Int4 len;

  if (a == NULL || a->xref_type == NULL) {
    str = StringSave ("No action");
  } else if (a->xref_type->choice != XrefType_gene
    || (g = (GeneXrefTypePtr) a->xref_type->data.ptrvalue) == NULL) {
    str = StringSave ("Invalid action");
  } else {
    label = GetFeatureNameFromFeatureType (g->feature);
    if (g->suppression < sizeof (s_Suppression) / sizeof (CharPtr)) {
      suppression = s_Suppression[g->suppression];
    } else {
      suppression = NULL;
    }
    if (g->necessary < sizeof (s_Necessary) / sizeof (CharPtr)) {
      necessary = s_Necessary[g->necessary];
    } else {
      necessary = NULL;
    }
    constraint = SummarizeConstraintSet (a->constraint);
    len = StringLen (label) + StringLen (fmt) + StringLen (suppression) + StringLen (necessary) + StringLen (constraint);
    str = (CharPtr) MemNew (sizeof (Char) * len);
    sprintf (str, fmt, 
             suppression == NULL ? "" : suppression, suppression == NULL ? "" : " ",
             necessary == NULL ? "" : necessary, necessary == NULL ? "" : " ",
             label);
    if (constraint != NULL) {
      StringCat (str, constraint);
      constraint = MemFree (constraint);
    }
  }
  return str;
}


static CharPtr SummarizeMakeGeneXrefs(MakeGeneXrefActionPtr a)
{
  CharPtr constraint, str, label;
  CharPtr fmt = "Make gene xrefs from overlapping gene features for %s features%s";

  if (a == NULL) {
    str = StringSave ("No action");
  } else {
    label = GetFeatureNameFromFeatureType (a->feature);
    constraint = SummarizeConstraintSet (a->constraint);
    str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (constraint) + StringLen (label)));
    sprintf (str, fmt, label, constraint == NULL ? "" : constraint);
    constraint = MemFree (constraint);
  }
  return str;
}


static CharPtr SummarizeMacroAction (ValNodePtr vnp)
{
  CharPtr str = NULL;

  if (vnp == NULL) {
    return StringSave ("No action");
  }
  switch (vnp->choice) {
    case MacroActionChoice_aecr:
      str = SummarizeAECRAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_parse:
      str = SummarizeParseAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_add_feature:
      str = SummarizeApplyFeatureAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_remove_feature:
      str = SummarizeRemoveFeatureAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_edit_location:
      str = SummarizeEditFeatureLocationAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_convert_feature:
      str = SummarizeConvertFeatureAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_remove_descriptor:
      str = SummarizeRemoveDescriptorAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_autodef:
      str = SummarizeAutodefAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_removesets:
      str = StringSave ("Remove duplicate nested sets");
      break;
    case MacroActionChoice_trim_junk_from_primer_seq:
      str = StringSave ("Trim junk from primer seqs");
      break;
    case MacroActionChoice_fix_usa_and_states:
      str = StringSave ("Fix USA and state abbreviations in publications");
      break;
    case MacroActionChoice_trim_stop_from_complete_cds:
      str = StringSave ("Remove trailing * from complete coding regions");
      break;
    case MacroActionChoice_synchronize_cds_partials:
      str = StringSave ("Synchronize coding region partials");
      break;
    case MacroActionChoice_adjust_for_consensus_splice:
      str = StringSave ("Adjust coding regions for consensus splice sites");
      break;
    case MacroActionChoice_fix_pub_caps:
      str = SummarizeFixPubCapsAction(vnp->data.ptrvalue);
      break;
    case MacroActionChoice_remove_seg_gaps:
      str = StringSave ("Remove seg-gaps");
      break;
    case MacroActionChoice_sort_fields:
      str = SummarizeSortFieldsAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_apply_molinfo_block:
      str = SummarizeMolinfoBlockAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_fix_caps:
      str = SummarizeFixCapsAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_fix_format:
      str = SummarizeFixFormatAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_fix_spell:
      str = StringSave ("Fix spelling");
      break;
    case MacroActionChoice_remove_duplicate_features:
      str = SummarizeRemoveDuplicateFeaturesAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_remove_lineage_notes:
      str = StringSave ("Remove lineage source notes");
      break;
    case MacroActionChoice_remove_xrefs:
      str = SummarizeRemoveXref(vnp->data.ptrvalue);
      break;
    case MacroActionChoice_make_gene_xrefs:
      str = SummarizeMakeGeneXrefs(vnp->data.ptrvalue);
      break;
    case MacroActionChoice_make_bold_xrefs:
      str = StringSave ("Make Barcode Xrefs");
      break;
    case MacroActionChoice_fix_author:
      str = SummarizeAuthorFixAction(vnp->data.ptrvalue);
      break;
    case MacroActionChoice_update_sequences:
      str = SummarizeUpdateSequencesAction (vnp->data.ptrvalue);
      break;
    case MacroActionChoice_add_trans_splicing:
      str = StringSave ("Set trans-splicing exception in genes");
      break;
    case MacroActionChoice_remove_invalid_ecnumbers:
      str = StringSave ("Remove invalid EC_numbers");
      break;
    default:
      str = StringSave ("Invalid action");
      break;
  }
  return str;
}


/* For applying selected actions from a macro script */
typedef struct macrorunform {
  FORM_MESSAGE_BLOCK
  DoC        macro_summary;
  FonT       summary_font;

  ValNodePtr macro_list;
  BoolPtr    ok_to_run;
  Int4       num_actions;
} MacroRunFormData, PNTR MacroRunFormPtr;


static void CleanupMacroRunForm (GraphiC g, VoidPtr data)

{
  MacroRunFormPtr f;

  f = (MacroRunFormPtr) data;
  if (f != NULL) {
    f->macro_list = MacroActionListFree (f->macro_list);
    f->ok_to_run = MemFree (f->ok_to_run);
  }
  StdCleanupFormProc (g, data);
}

static void ClickMacroRunDoc (DoC d, PoinT pt)
{
  Int2            item, row, col;
  RecT            rct;
  MacroRunFormPtr f;
  
  f = (MacroRunFormPtr) GetObjectExtra (d);
  if (f == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item > 0 && item <= f->num_actions && col == 1) {
    /* check or uncheck item */
    if (f->ok_to_run[item - 1]) {
      f->ok_to_run[item - 1] = FALSE;
    } else {
      f->ok_to_run[item - 1] = TRUE;
    }
    InvalDocument (f->macro_summary);
  }
}


static void DrawMacroRunControls (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  MacroRunFormPtr dlg;
  RecT            rct;
  Int4            width;
  PoinT           pt1, pt2;

  dlg = (MacroRunFormPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  /* don't draw controls for explanatory text */
  if (item == 0 || item > dlg->num_actions) return;

  if (dlg != NULL && r != NULL && item > 0 && item <= dlg->num_actions && firstLine == 0) {
    rct = *r;

    width = 10;
    /* draw box */
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1;
    pt2.x = pt1.x + width;
    pt2.y = pt1.y;
    DrawLine (pt1, pt2);
    pt1.x = pt2.x;
    pt1.y = pt2.y + width;
    DrawLine (pt1, pt2);
    pt2.x = pt1.x - width;
    pt2.y = pt1.y;
    DrawLine (pt1, pt2);
    pt1.x = pt2.x;
    pt1.y = pt2.y - width;
    DrawLine (pt1, pt2);

    if (dlg->ok_to_run[item - 1]) {
      /* draw X for check */
      pt1.x = rct.left + 1;
      pt1.y = rct.top + 1;
      pt2.x = pt1.x + width;
      pt2.y = pt1.y + width;
      DrawLine (pt1, pt2);
      pt1.x = rct.left + 1;
      pt1.y = rct.top + 1 + width;
      pt2.x = pt1.x + width;
      pt2.y = rct.top + 1;
      DrawLine (pt1, pt2);
    }
  }
}


static void RunSelectedMacro (ButtoN b)
{
  MacroRunFormPtr  f;
  ValNodePtr       vnp, tmp_list = NULL, last = NULL, newaction, tmp;
  Int4             i;
  ValNodePtr       sep_list;
  LogInfoPtr       lip;
  SeqEntryPtr      sep;
  Uint2            entityID;

  f = (MacroRunFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }

  /* make temporary macro list */
  for (vnp = f->macro_list, i = 0; vnp != NULL && i < f->num_actions; vnp = vnp->next, i++) {
    if (f->ok_to_run[i]) {
      tmp = vnp->next;
      vnp->next = NULL;
      newaction = AsnIoMemCopy (vnp, (AsnReadFunc) MacroActionListAsnRead, (AsnWriteFunc) MacroActionListAsnWrite);
      vnp->next = tmp;
      if (last == NULL) {
        tmp_list = newaction;
      } else {
        last->next = newaction;
      }
      last = newaction;
    }
  }

  if (tmp_list == NULL) {
    Message (MSG_ERROR, "No actions selected!");
    return;
  }

  /* run it */
  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    Message (MSG_ERROR, "No records open!");
  } else if (sep_list->next != NULL 
    && ANS_CANCEL == Message (MSG_OKC, "You have more than one record open - run macro for all open records?")) {
    /* do nothing */
  } else {
    WatchCursor();
    Update();
    lip = OpenLog ("Macro Actions");
    for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
      sep = vnp->data.ptrvalue;
      entityID = ObjMgrGetEntityIDForChoice(sep);
      lip->data_in_log |= ApplyMacroToSeqEntryEx (sep, tmp_list, lip->fp, Sequin_GlobalAlign2Seq);
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
    sep_list = ValNodeFree (sep_list);
    ArrowCursor ();
    Update ();  
    /* show log from actions */
    if (!lip->data_in_log) {
      Message (MSG_OK, "Macro had no effect\n");
    }
    CloseLog (lip);
    lip = FreeLog (lip);
  }
  /* free temporary macro list */
  tmp_list = MacroActionListFree (tmp_list);
}


static void PrintMacroList(ButtoN b)
{
  MacroRunFormPtr  f;
  ValNodePtr       vnp;
  CharPtr          str;
  Char             path [PATH_MAX];
  FILE * fp;

  f = (MacroRunFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }

  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp != NULL) {
    for (vnp = f->macro_list; vnp != NULL; vnp = vnp->next) {
      str = SummarizeMacroAction (vnp);
      fprintf (fp, "%s\n", str == NULL ? "Unknown action" : str);
      str = MemFree (str);
    }
    FileClose (fp);
    LaunchGeneralTextViewer (path, "Autofix");
  }
  FileRemove (path);
}


static void SummarizeRunMacro (DoC doc, ValNodePtr macro_list, FonT font)
{
  ValNodePtr vnp;
  CharPtr    str;
  CharPtr    tmp;
  RecT       r;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };


  Reset (doc);

  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  
  ColFmt[1].pixWidth = r.right - r.left - 12;

  if (font == NULL) font = programFont;

  for (vnp = macro_list; vnp != NULL; vnp = vnp->next) {
    str = SummarizeMacroAction (vnp);
    tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 6));
    sprintf (tmp, "\t%s\n", str);
    str = MemFree (str);
    AppendText (doc, tmp, &ParFmt, ColFmt, font);
    tmp = MemFree (tmp);
  }
  UpdateDocument (doc, 0, 0);
}


NLM_EXTERN void SelectiveMacroRun (ValNodePtr macro_list)
{
  WindoW           w;
  MacroRunFormPtr  f;
  GrouP            h, c;
  ButtoN           b;
  Int4             i;

  if (macro_list == NULL) {
    return;
  }
  f = (MacroRunFormPtr) MemNew (sizeof (MacroRunFormData));
  if (f == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Autofix", StdCloseWindowProc);
  SetObjectExtra (w, f, CleanupMacroRunForm);
  f->form = (ForM) w;

  f->macro_list = macro_list;
  f->num_actions = ValNodeLen (macro_list);
  f->ok_to_run = (BoolPtr) MemNew (sizeof (Boolean) * f->num_actions);
  /* by default, run everything */
  for (i = 0; i < f->num_actions; i++) {
    f->ok_to_run[i] = TRUE;
  }

  /* set up font */
#ifdef WIN_MAC
  f->summary_font = ParseFont ("Times,12");
#endif
#ifdef WIN_MSWIN
  f->summary_font = ParseFont ("Times New Roman,12");
#endif
#ifdef WIN_MOTIF
  f->summary_font = ParseFont ("Times,12");
#endif

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  f->macro_summary = DocumentPanel (h, stdCharWidth * 50, stdLineHeight * 30);
  SetObjectExtra (f->macro_summary, f, NULL);
  SetDocProcs (f->macro_summary, ClickMacroRunDoc, NULL, NULL, NULL);
  SetDocShade (f->macro_summary, DrawMacroRunControls, NULL, NULL, NULL);

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Run", RunSelectedMacro);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Print Macro List", PrintMacroList);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Close", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) f->macro_summary, (HANDLE) c, NULL);
  SummarizeRunMacro (f->macro_summary, f->macro_list, f->summary_font);
  Show (w);
}


typedef struct structuredcommentdatabasenamedialog {
  DIALOG_MESSAGE_BLOCK

  DialoG known_names;
  TexT   new_name;
  GrouP  known_or_new;

  ValNodePtr prefix_list;  /* note - prefix_list will be freed by the known_names dialg */
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} StructuredCommentDatabaseNameDlgData, PNTR StructuredCommentDatabaseNameDlgPtr;


static void ChangeStructuredCommentDatabaseNameNewOrKnown(GrouP g)
{
  StructuredCommentDatabaseNameDlgPtr dlg;

  dlg = (StructuredCommentDatabaseNameDlgPtr) GetObjectExtra (g);
  if (dlg == NULL) {
    return;
  }

  if (GetValue (dlg->known_or_new) == 1) {
    Enable (dlg->known_names);
    Disable (dlg->new_name);
  } else {
    Enable (dlg->known_names);
    Disable (dlg->new_name);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void StructuredCommentDatabaseNameToDialog (DialoG d, Pointer data)
{
  StructuredCommentDatabaseNameDlgPtr dlg;
  ValNodePtr vnp;
  CharPtr    str = (CharPtr) data;
  ValNode    vn;
  Boolean    is_known = FALSE;

  dlg = (StructuredCommentDatabaseNameDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  if (str == NULL) {
    PointerToDialog (dlg->known_names, NULL);
    SetTitle (dlg->new_name, "");
    SetValue (dlg->known_or_new, 2);
    ChangeStructuredCommentDatabaseNameNewOrKnown (dlg->known_or_new);
  } else {
    /* first, determine if known or unknown */
    for (vnp = dlg->prefix_list; vnp != NULL && !is_known; vnp = vnp->next) {
      if (StringCmp ((CharPtr) vnp->data.ptrvalue, str) == 0) {
        is_known = TRUE;
      }
    }
    if (is_known) {
      MemSet (&vn, 0, sizeof (ValNode));
      vn.data.ptrvalue = str;
      PointerToDialog (dlg->known_names, &vn);
      SetTitle (dlg->new_name, "");
      SetValue (dlg->known_or_new, 1);
    } else {
      PointerToDialog (dlg->known_names, NULL);
      SetTitle (dlg->new_name, str == NULL ? "" : str);
      SetValue (dlg->known_or_new, 2);
    }
    ChangeStructuredCommentDatabaseNameNewOrKnown (dlg->known_or_new);
  }
}


static Pointer StructuredCommentDatabaseNameFromDialog (DialoG d)
{
  StructuredCommentDatabaseNameDlgPtr dlg;
  ValNodePtr vnp;
  CharPtr    str = NULL;

  dlg = (StructuredCommentDatabaseNameDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  if (GetValue(dlg->known_or_new) == 1) {
    vnp = DialogToPointer (dlg->known_names);
    if (vnp != NULL) {
      str = StringSave (vnp->data.ptrvalue);
    }
  } else {
    str = SaveStringFromText (dlg->new_name);
  }
  return str;
}


static void ChangeStructuredCommentDatabaseNameText(TexT t)
{
  StructuredCommentDatabaseNameDlgPtr dlg;

  dlg = (StructuredCommentDatabaseNameDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


NLM_EXTERN DialoG StructuredCommentDatabaseNameDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  StructuredCommentDatabaseNameDlgPtr dlg;
  GrouP p;

  dlg = (StructuredCommentDatabaseNameDlgPtr) MemNew (sizeof (StructuredCommentDatabaseNameDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = StructuredCommentDatabaseNameFromDialog;
  dlg->todialog = StructuredCommentDatabaseNameToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->known_or_new = HiddenGroup (p, 0, 2, ChangeStructuredCommentDatabaseNameNewOrKnown);
  SetObjectExtra (dlg->known_or_new, dlg, NULL);
  SetGroupSpacing (dlg->known_or_new, 10, 10);
  RadioButton (dlg->known_or_new, "Known");
  RadioButton (dlg->known_or_new, "New");
  dlg->prefix_list = GetStructuredCommentPrefixList();
  dlg->known_names = ValNodeSelectionDialog (dlg->known_or_new, 
                                             dlg->prefix_list,
                                             TALL_SELECTION_LIST,
                                             ValNodeStringName,
                                             ValNodeSimpleDataFree,
                                             ValNodeStringCopy,
                                             ValNodeStringMatch,
                                             "dbname",
                                             change_notify, change_userdata, FALSE);
  dlg->new_name = DialogText (dlg->known_or_new, "", 10, ChangeStructuredCommentDatabaseNameText);
  SetObjectExtra (dlg->new_name, dlg, NULL);
  SetValue (dlg->known_or_new, 2);
  Disable (dlg->known_names);

  return (DialoG) p;
}


/* for editing macro "templates" where we only want to change the values */
/* any action that we don't have a template editor for, just print a string describing the action */

typedef struct macrotemplatestringdlg {
  DIALOG_MESSAGE_BLOCK
  DoC summ;
  ValNodePtr original_action;
} MacroTemplateStringDlgData, PNTR MacroTemplateStringDlgPtr;


static Pointer MacroFromMacroTemplateStringDialog (DialoG d)
{
  MacroTemplateStringDlgPtr dlg;

  dlg = (MacroTemplateStringDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  return AsnIoMemCopy (dlg->original_action, (AsnReadFunc) MacroActionChoiceAsnRead, (AsnWriteFunc) MacroActionChoiceAsnWrite);
}


static void MacroToMacroTemplateStringDialog (DialoG d, Pointer data)
{
  MacroTemplateStringDlgPtr dlg;
  CharPtr summ;

  dlg = (MacroTemplateStringDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  dlg->original_action = (ValNodePtr) data;

  summ = SummarizeMacroAction (dlg->original_action);
  Reset(dlg->summ);
  AppendText (dlg->summ, summ, NULL, NULL, programFont);
  summ = MemFree (summ);
}


NLM_EXTERN DialoG MacroTemplateStringDialog (GrouP h)
{
  MacroTemplateStringDlgPtr dlg;
  GrouP p;

  dlg = (MacroTemplateStringDlgPtr) MemNew (sizeof (MacroTemplateStringDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = MacroFromMacroTemplateStringDialog;
  dlg->todialog = MacroToMacroTemplateStringDialog;
  
  dlg->summ = DocumentPanel (p, stdCharWidth * 50, stdLineHeight * 1);


  return (DialoG) p;
}


typedef struct macrotemplateapplytextdlg {
  DIALOG_MESSAGE_BLOCK
  DoC summ;
  TexT text;
  ValNodePtr original_action;
} MacroTemplateApplyTextDlgData, PNTR MacroTemplateApplyTextDlgPtr;


static Pointer MacroFromMacroTemplateApplyTextDialog (DialoG d)
{
  MacroTemplateApplyTextDlgPtr dlg;
  ValNodePtr macro;
  AECRActionPtr action;
  ApplyActionPtr apply;

  dlg = (MacroTemplateApplyTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  macro = AsnIoMemCopy (dlg->original_action, (AsnReadFunc) MacroActionChoiceAsnRead, (AsnWriteFunc) MacroActionChoiceAsnWrite);
  if (macro != NULL && macro->choice == MacroActionChoice_aecr 
      && (action = (AECRActionPtr) macro->data.ptrvalue) != NULL
      && action->action != NULL
      && action->action->choice == ActionChoice_apply
      && (apply = (ApplyActionPtr) action->action->data.ptrvalue) != NULL) {
    apply->value = MemFree (apply->value);
    apply->value = JustSaveStringFromText (dlg->text);
  }
  return macro;
}


static void MacroToMacroTemplateApplyTextDialog (DialoG d, Pointer data)
{
  MacroTemplateApplyTextDlgPtr dlg;
  ValNodePtr macro;
  AECRActionPtr action;
  ApplyActionPtr apply;
  CharPtr field, existing, constraint;
  CharPtr fmt = "to %s (%s)";
  CharPtr summ;

  dlg = (MacroTemplateApplyTextDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  dlg->original_action = (ValNodePtr) data;
  macro = dlg->original_action;

  Reset(dlg->summ);
  if ((macro = (ValNodePtr)data) == NULL || macro->choice != MacroActionChoice_aecr 
      || (action = (AECRActionPtr) macro->data.ptrvalue) == NULL
      || action->action == NULL
      || action->action->choice != ActionChoice_apply
      || (apply = (ApplyActionPtr) action->action->data.ptrvalue) == NULL
      || IsNonTextFieldType(apply->field)) {
    SetTitle (dlg->text, "");
  } else {
    SetTitle (dlg->text, apply->value);
    field = SummarizeFieldType (apply->field);
    existing = SummarizeExistingText (apply->existing_text);
    constraint = SummarizeConstraintSet (action->constraint);
    summ = (CharPtr) MemNew (sizeof (Char) * 
          (StringLen (fmt) + StringLen (field) + StringLen (existing) + StringLen (constraint)));
    sprintf (summ, fmt, field, existing);
    if (constraint != NULL) {
      StringCat (summ, constraint);
    }
    field = MemFree (field);
    constraint = MemFree (constraint);
    AppendText (dlg->summ, summ, NULL, NULL, programFont);
    summ = MemFree (summ);
  }
}


NLM_EXTERN DialoG MacroTemplateApplyTextDialog (GrouP h)
{
  MacroTemplateApplyTextDlgPtr dlg;
  GrouP p;

  dlg = (MacroTemplateApplyTextDlgPtr) MemNew (sizeof (MacroTemplateApplyTextDlgData));

  p = HiddenGroup (h, 3, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = MacroFromMacroTemplateApplyTextDialog;
  dlg->todialog = MacroToMacroTemplateApplyTextDialog;
  
  StaticPrompt (p, "Apply", 0, popupMenuHeight, programFont, 'r');
  dlg->text = DialogText (p, "", 15, NULL);
  dlg->summ = DocumentPanel (p, stdCharWidth * 35, stdLineHeight * 1);

  return (DialoG) p;
}


static Boolean IsTextApply (ValNodePtr macro)
{
  AECRActionPtr action;
  ApplyActionPtr apply;

  if (macro == NULL || macro->choice != MacroActionChoice_aecr 
      || (action = (AECRActionPtr) macro->data.ptrvalue) == NULL
      || action->action == NULL
      || action->action->choice != ActionChoice_apply
      || (apply = (ApplyActionPtr) action->action->data.ptrvalue) == NULL
      || IsNonTextFieldType(apply->field)) {
    return FALSE;
  } else {
    return TRUE;
  }
}


typedef struct macrotemplateform {
  FORM_MESSAGE_BLOCK

  DialoG PNTR action_dlgs;
  ButtoN PNTR keep_btn;

  Int4 num_actions;
} MacroTemplateFormData, PNTR MacroTemplateFormPtr;

static void CleanupMacroTemplateForm (GraphiC g, VoidPtr data)

{
  MacroTemplateFormPtr f;

  f = (MacroTemplateFormPtr) data;
  if (f != NULL) {
    f->action_dlgs = MemFree (f->action_dlgs);
  }
  StdCleanupFormProc (g, data);
}


static void RunTemplateMacro (ButtoN b)
{
  ValNodePtr actions = NULL, last_action = NULL, new_action, vnp, sep_list;
  MacroTemplateFormPtr  f;
  SeqEntryPtr sep;
  LogInfoPtr lip;
  Uint2 entityID;
  Int4 i;

  f = (MacroTemplateFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }
  for (i = 0; i < f->num_actions; i++) {
    if (GetStatus (f->keep_btn[i])) {
      new_action = DialogToPointer (f->action_dlgs[i]);
      if (new_action != NULL) {
        if (last_action == NULL) {
          actions = new_action;
        } else {
          last_action->next = new_action;
        }
        last_action = new_action;
      }
    }
  }
  /* run it */
  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    Message (MSG_ERROR, "No records open!");
  } else if (sep_list->next != NULL 
    && ANS_CANCEL == Message (MSG_OKC, "You have more than one record open - run macro for all open records?")) {
    /* do nothing */
  } else {
    WatchCursor();
    Update();
    lip = OpenLog ("Macro Actions");
    for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
      sep = vnp->data.ptrvalue;
      entityID = ObjMgrGetEntityIDForChoice(sep);
      lip->data_in_log |= ApplyMacroToSeqEntryEx (sep, actions, lip->fp, Sequin_GlobalAlign2Seq);
      ObjMgrSetDirtyFlag (entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
    sep_list = ValNodeFree (sep_list);
    ArrowCursor ();
    Update ();  
    /* show log from actions */
    if (!lip->data_in_log) {
      Message (MSG_OK, "Macro had no effect\n");
    }
    CloseLog (lip);
    lip = FreeLog (lip);
  }
  /* free temporary macro list */
  actions = MacroActionListFree (actions);
  
}


static void SaveTemplateValues (ButtoN b)
{
  ValNodePtr actions = NULL, last_action = NULL, new_action;
  MacroTemplateFormPtr  f;
  Int4 i;
  Char             path [PATH_MAX];
  AsnIoPtr         aip;

  f = (MacroTemplateFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }
  for (i = 0; i < f->num_actions; i++) {
    new_action = DialogToPointer (f->action_dlgs[i]);
    if (new_action != NULL) {
      if (last_action == NULL) {
        actions = new_action;
      } else {
        last_action->next = new_action;
      }
      last_action = new_action;
    }
  }
  path [0] = '\0';
  if (GetOutputFileName (path, sizeof (path), NULL)) {
    aip = AsnIoOpen (path, "w");
    if (aip == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
    } else {
      MacroActionListAsnWrite (actions, aip, NULL);
      AsnIoClose (aip);
    }
  }
  actions = MacroActionListFree (actions);
}


NLM_EXTERN void EditMacroTemplate (void)
{
  WindoW           w;
  MacroTemplateFormPtr  f;
  GrouP            h, g, c;
  ButtoN           b;
  Int4             i;
  ValNodePtr       vnp;
  Char             path [PATH_MAX];
  AsnIoPtr         aip;
  ValNodePtr       macro_list;

  path [0] = '\0';
  if (!GetInputFileName (path, sizeof (path), NULL, "TEXT")) {
    return;
  } else if ((aip = AsnIoOpen (path, "r")) == NULL) {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  } else if ((macro_list = MacroActionListAsnRead (aip, NULL)) == NULL) {
    Message (MSG_ERROR, "Unable to read action list from %s.", path);
    aip = AsnIoClose (aip);
    return;
  }
  AsnIoClose (aip);

  f = (MacroTemplateFormPtr) MemNew (sizeof (MacroTemplateFormData));
  if (f == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Macro Template", StdCloseWindowProc);
  SetObjectExtra (w, f, CleanupMacroTemplateForm);
  f->form = (ForM) w;

  f->num_actions = ValNodeLen (macro_list);
  f->action_dlgs = (DialoG PNTR) MemNew (sizeof (DialoG) * f->num_actions);
  f->keep_btn = (ButtoN PNTR) MemNew (sizeof (ButtoN) * f->num_actions);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);
  for (i = 0, vnp = macro_list; i < f->num_actions && vnp != NULL; i++, vnp = vnp->next) {
    f->keep_btn[i] = CheckBox (g, "", NULL);
    SetStatus (f->keep_btn[i], TRUE);
    if (IsTextApply(vnp)) {
      f->action_dlgs[i] = MacroTemplateApplyTextDialog(g);
    } else {
      f->action_dlgs[i] = MacroTemplateStringDialog (g);
    }
    PointerToDialog (f->action_dlgs[i], vnp);
  }

  c = HiddenGroup (h, 3, 0, NULL);
  b = PushButton (c, "Run", RunTemplateMacro);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Save Template with New Values", SaveTemplateValues);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Close", StdCancelButtonProc); 

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  Show (w);
}


NLM_EXTERN void LaunchMacroTemplateEditor (IteM i)
{
  BaseFormPtr         bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  EditMacroTemplate ();
}


/* For using macro functions to apply tables */

typedef struct idmatchlocationdlg {
  DIALOG_MESSAGE_BLOCK
  PopuP  match_location;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} IdMatchLocationDlgData, PNTR IdMatchLocationDlgPtr;


typedef enum {
  eTabColumnConfigMatchLocation_equals = 1,
  eTabColumnConfigMatchLocation_contains,
  eTabColumnConfigMatchLocation_starts,
  eTabColumnConfigMatchLocation_ends,
  eTabColumnConfigMatchLocation_inlist
} ETabColumnConfigMatchLocation;


static Uint1 StringConstraintLocationFromTabColumnConfigMatchLocation (Int2 i)
{
  Uint1 val = String_location_equals;
  switch (i) {
    case eTabColumnConfigMatchLocation_equals:
      val = String_location_equals;
      break;
    case eTabColumnConfigMatchLocation_contains:
      val = String_location_contains;
      break;
    case eTabColumnConfigMatchLocation_starts:
      val = String_location_starts;
      break;
    case eTabColumnConfigMatchLocation_ends:
      val = String_location_ends;
      break;
    case eTabColumnConfigMatchLocation_inlist:
      val = String_location_inlist;
      break;
    default:
      val = String_location_equals;
      break;
  }
  return val;
} 


static Int2 TabColumnConfigMatchLocationFromStringConstraintLocation (Uint1 u)
{
  Int2 val = eTabColumnConfigMatchLocation_equals;

  switch (u) {
    case String_location_equals:
      val = eTabColumnConfigMatchLocation_equals;
      break;
    case String_location_contains:
      val = eTabColumnConfigMatchLocation_contains;
      break;
    case String_location_starts:
      val = eTabColumnConfigMatchLocation_starts;
      break;
    case String_location_ends:
      val = eTabColumnConfigMatchLocation_ends;
      break;
    case String_location_inlist:
      val = eTabColumnConfigMatchLocation_inlist;
      break;
    default:
      val = eTabColumnConfigMatchLocation_equals;
      break;
  }
  return val;
} 


static void ChangeIdMatchLocationPopup (PopuP p)
{
  IdMatchLocationDlgPtr dlg;

  dlg = (IdMatchLocationDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


NLM_EXTERN DialoG IdMatchLocationDlg (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  IdMatchLocationDlgPtr dlg;
  GrouP p;

  dlg = (IdMatchLocationDlgPtr) MemNew (sizeof (IdMatchLocationDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  dlg->match_location = PopupList (p, TRUE, ChangeIdMatchLocationPopup);
  SetObjectExtra (dlg->match_location, dlg, NULL);
  PopupItem (dlg->match_location, "Matches");
  PopupItem (dlg->match_location, "Is contained in");
  PopupItem (dlg->match_location, "Is start of");
  PopupItem (dlg->match_location, "Is end of");
  PopupItem (dlg->match_location, "List contains");
  SetValue (dlg->match_location, 1);

  return (DialoG) p;
}


NLM_EXTERN Uint1 GetMatchLocationFromIdMatchLocationDlg (DialoG d)
{
  IdMatchLocationDlgPtr dlg;

  dlg = (IdMatchLocationDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return String_location_equals;
  }

  return StringConstraintLocationFromTabColumnConfigMatchLocation (GetValue (dlg->match_location));
}


NLM_EXTERN void SetMatchLocationInIdMatchLocationDlg (DialoG d, Uint1 match_location)
{
  IdMatchLocationDlgPtr dlg;

  dlg = (IdMatchLocationDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  SetValue (dlg->match_location, TabColumnConfigMatchLocationFromStringConstraintLocation (match_location));
}


/* for configuring how a column in a table will be used */
typedef struct tabcolumnconfigdlg {
  DIALOG_MESSAGE_BLOCK
  PopuP                    column_action;
  GrouP                    match_grp;
  DialoG                   match_location;
  DialoG                   src_qual_match;
  GrouP                    qual_grp;
  DialoG                   src_qual;
  GrouP                    feature_field_grp;
  DialoG                   feature_type;
  DialoG                   feature_field;
  DialoG                   cdsgeneprot;
  DialoG                   pub_field;
  DialoG                   molinfo_field;
  DialoG                   struccomm_field;
  PopuP                    dblink_field;

  ButtoN                   change_mrna;
  GrouP                    apply_options;
  ButtoN                   erase_when_blank;
  DialoG                   existing_text;
  DialoG                   constraint;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} TabColumnConfigDlgData, PNTR TabColumnConfigDlgPtr;

typedef enum {
  eTabColumnConfig_Ignore = 1,
  eTabColumnConfig_Match_NucleotideID,
  eTabColumnConfig_Match_Taxname,
  eTabColumnConfig_Match_SourceQual,
  eTabColumnConfig_Match_FeatureID,
  eTabColumnConfig_Match_GeneLocusTag,
  eTabColumnConfig_Match_ProteinID,
  eTabColumnConfig_Match_ProteinName,  /* J. Chen */
  eTabColumnConfig_Match_Dbxref,
  eTabColumnConfig_Blank,
  eTabColumnConfig_Apply_Taxname,
  eTabColumnConfig_Apply_SourceQual,
  eTabColumnConfig_Apply_CDSGeneProtField,
  eTabColumnConfig_Apply_FeatureField,
  eTabColumnConfig_Apply_PubField,
  eTabColumnConfig_Apply_GenomeProjectId,
  eTabColumnConfig_Apply_CommentDescriptor,
  eTabColumnConfig_Apply_Defline,
  eTabColumnConfig_Apply_Keyword,
  eTabColumnConfig_Apply_Molinfo,
  eTabColumnConfig_Apply_StructuredCommentField,
  eTabColumnConfig_Apply_Dblink,
} ETabColumnConfig;

static void TabColumnConfigFieldChange (Pointer data)
{
  TabColumnConfigDlgPtr dlg;
  TabColumnConfigPtr    config;

  dlg = (TabColumnConfigDlgPtr) data;
  if (dlg != NULL) {
    config = (TabColumnConfigPtr) DialogToPointer (dlg->dialog);
    if (config == NULL) {
      Disable (dlg->change_mrna);
      DisableNonTextOptions (dlg->existing_text);
      DisableMultiOptions (dlg->existing_text);
    } else {
      if (!IsFieldTypeCDSProduct(config->field)) {
        Disable (dlg->change_mrna);
      } else {
        Enable (dlg->change_mrna);
      }
      if (AllowFieldMulti (config->field)) {
        EnableMultiOptions (dlg->existing_text);
      } else {
        DisableMultiOptions (dlg->existing_text);
      }
    }

    config = TabColumnConfigFree (config);
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static Pointer TabColumnConfigDialogToChoice (DialoG d)
{
  TabColumnConfigDlgPtr dlg;
  TabColumnConfigPtr    f;
  Int2                  val;
  FeatureFieldPtr       ff;
  ValNodePtr            vnp;

  dlg = (TabColumnConfigDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  f = TabColumnConfigNew ();
  f->match_mrna = FALSE;
  
  val = GetValue (dlg->column_action);
  if (val < 2 || val == eTabColumnConfig_Blank) {
    f = TabColumnConfigFree (f);
  } else if (val < eTabColumnConfig_Blank) {
    switch (val) {
      case eTabColumnConfig_Match_ProteinName:   /* J. Chen */
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchProteinName;
        break;
      case eTabColumnConfig_Match_FeatureID:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchFeatureID;
        break;
      case eTabColumnConfig_Match_GeneLocusTag:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchGeneLocusTag;
        break;
      case eTabColumnConfig_Match_ProteinID:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchProteinID;
        f->match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
      case eTabColumnConfig_Match_Dbxref:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchDbxref;
        f->match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
      case eTabColumnConfig_Match_NucleotideID:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchNucID;
        f->match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
      case eTabColumnConfig_Match_Taxname:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchBioSource;
        f->match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
      case eTabColumnConfig_Match_SourceQual:
        f->match_type = MatchTypeNew ();
        f->match_type->choice = eTableMatchSourceQual;
        f->match_type->data = DialogToPointer (dlg->src_qual_match);
        f->match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
    }
    f->existing_text = ExistingTextOption_replace_old;
    f->skip_blank = TRUE;
  } else {
    if (val == eTabColumnConfig_Apply_Taxname) {
      vnp = ValNodeNew (NULL);
      vnp->choice = SourceQualChoice_textqual;
      vnp->data.intvalue = Source_qual_taxname;
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_source_qual;
      f->field->data.ptrvalue = vnp;
    } else if (val == eTabColumnConfig_Apply_SourceQual) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_source_qual;
      f->field->data.ptrvalue = DialogToPointer (dlg->src_qual);
    } else if (val == eTabColumnConfig_Apply_FeatureField) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_feature_field;
      ff = FeatureFieldNew ();
      vnp = DialogToPointer (dlg->feature_type);
      if (vnp != NULL) {
        ff->type = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
      vnp = DialogToPointer (dlg->feature_field);
      if (vnp != NULL) {
        ValNodeAddInt (&(ff->field), FeatQualChoice_legal_qual, vnp->choice);
        ValNodeFree (vnp);
      }
      f->field->data.ptrvalue = ff;
    } else if (val == eTabColumnConfig_Apply_CDSGeneProtField) {
      f->field = DialogToPointer (dlg->cdsgeneprot);
    } else if (val == eTabColumnConfig_Apply_PubField) {
      vnp = DialogToPointer (dlg->pub_field);
      if (vnp != NULL) {
        f->field = ValNodeNew (NULL);
        f->field->choice = FieldType_pub;
        f->field->data.intvalue = vnp->choice;
        vnp = ValNodeFree (vnp);
      }
    } else if (val == eTabColumnConfig_Apply_GenomeProjectId) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_misc;
      f->field->data.intvalue = Misc_field_genome_project_id;
    } else if (val == eTabColumnConfig_Apply_CommentDescriptor) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_misc;
      f->field->data.intvalue = Misc_field_comment_descriptor;
    } else if (val == eTabColumnConfig_Apply_Defline) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_misc;
      f->field->data.intvalue = Misc_field_defline;
    } else if (val == eTabColumnConfig_Apply_Keyword) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_misc;
      f->field->data.intvalue = Misc_field_keyword;
    } else if (val == eTabColumnConfig_Apply_Molinfo) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_molinfo_field;
      f->field->data.ptrvalue = DialogToPointer (dlg->molinfo_field);
    } else if (val == eTabColumnConfig_Apply_StructuredCommentField) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_struc_comment_field;
      f->field->data.ptrvalue = DialogToPointer (dlg->struccomm_field);
    } else if (val == eTabColumnConfig_Apply_Dblink) {
      f->field = ValNodeNew (NULL);
      f->field->choice = FieldType_dblink;
      f->field->data.intvalue = GetValue (dlg->dblink_field);
    }


    if (dlg->erase_when_blank == NULL) {
      f->skip_blank = TRUE;
    } else {
      f->skip_blank = !GetStatus (dlg->erase_when_blank);
    }
    if (IsFieldTypeCDSProduct (f->field)) {
      f->match_mrna = GetStatus (dlg->change_mrna);
    }

    f->existing_text = GetExistingTextDialogValue (dlg->existing_text);
    f->constraint = DialogToPointer (dlg->constraint);
  }

  return f;
}


static void TabColumnActionChange (PopuP p)
{
  TabColumnConfigDlgPtr dlg;
  TabColumnConfigPtr    config;
  Int2 val;

  dlg = (TabColumnConfigDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  Hide (dlg->dblink_field);

  val = GetValue (p);

  switch (val) {
    case eTabColumnConfig_Ignore:
    case eTabColumnConfig_Blank:
    case eTabColumnConfig_Match_FeatureID:
    case eTabColumnConfig_Match_GeneLocusTag:
    case eTabColumnConfig_Match_ProteinID:
    case eTabColumnConfig_Match_Dbxref:
    case eTabColumnConfig_Match_NucleotideID:
    case eTabColumnConfig_Match_Taxname:
      /* matching, hide apply controls */
      Show (dlg->match_grp);
      Hide (dlg->src_qual_match);
      Hide (dlg->qual_grp);
      Disable (dlg->apply_options);
      break;
    case eTabColumnConfig_Match_SourceQual:
      /* matching to source qual, hide apply controls, show sourcequal match dialog */
      Hide (dlg->qual_grp);
      Disable (dlg->apply_options);
      Show (dlg->match_grp);
      Show (dlg->src_qual_match);
      break;
    case eTabColumnConfig_Apply_Taxname:
      /* choice is taxname, no need to show anything */
      Hide (dlg->qual_grp);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_SourceQual:
      /* show source qual, hide others */
      Show (dlg->qual_grp);
      Show (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      Disable (dlg->change_mrna);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_FeatureField:
      /* show feature qual, hide others */
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Show (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      config = DialogToPointer (dlg->dialog);
      if (config == NULL || !IsFieldTypeCDSProduct (config->field)) {
        SafeDisable (dlg->change_mrna);
      } else {
        SafeEnable (dlg->change_mrna);
      }
      config = TabColumnConfigFree (config);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_CDSGeneProtField:
      /* show cds-gene-prot qual, hide others */
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Show (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      config = DialogToPointer (dlg->dialog);
      if (config == NULL || !IsFieldTypeCDSProduct (config->field)) {
        SafeDisable (dlg->change_mrna);
      } else {
        SafeEnable (dlg->change_mrna);
      }
      config = TabColumnConfigFree (config);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_PubField:
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Show (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_GenomeProjectId:
    case eTabColumnConfig_Apply_CommentDescriptor:
    case eTabColumnConfig_Apply_Defline:
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      Disable (dlg->change_mrna);
      DisableMultiOptions (dlg->existing_text);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_Keyword:
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      Disable (dlg->change_mrna);
      EnableMultiOptions (dlg->existing_text);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_Molinfo:
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Show (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      Disable (dlg->change_mrna);
      DisableMultiOptions (dlg->existing_text);
      DisableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_StructuredCommentField:
      Show (dlg->qual_grp);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Show (dlg->struccomm_field);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      Disable (dlg->change_mrna);
      DisableMultiOptions (dlg->existing_text);
      EnableNonTextOptions (dlg->existing_text);
      break;
    case eTabColumnConfig_Apply_Dblink:
      Show (dlg->qual_grp);
      Show (dlg->dblink_field);
      Hide (dlg->src_qual);
      Hide (dlg->feature_field_grp);
      Hide (dlg->cdsgeneprot);
      Hide (dlg->pub_field);
      Hide (dlg->molinfo_field);
      Hide (dlg->struccomm_field);
      Enable (dlg->apply_options);
      Hide (dlg->match_grp);
      Disable (dlg->change_mrna);
      EnableMultiOptions (dlg->existing_text);
      DisableNonTextOptions (dlg->existing_text);
      break;
  }
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void TabColumnConfigToDialog (DialoG d, Pointer data)
{
  TabColumnConfigDlgPtr dlg;
  TabColumnConfigPtr    f;
  FeatureFieldPtr       ff;
  ValNode               vn;
  ValNodePtr            src_field;

  dlg = (TabColumnConfigDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  f = (TabColumnConfigPtr) data;
  if (f == NULL) {
    SetValue (dlg->column_action, 1);
    PointerToDialog (dlg->constraint, NULL);
  } else {
    f = TabColumnConfigCopy (f);
    if (f->match_type != NULL) {
      switch (f->match_type->choice) {
        case eTableMatchProteinName:    /* J. Chen */
          SetValue (dlg->column_action, eTabColumnConfig_Match_ProteinName);
          break;
        case eTableMatchFeatureID:
          SetValue (dlg->column_action, eTabColumnConfig_Match_FeatureID);
          break;
        case eTableMatchGeneLocusTag:
          SetValue (dlg->column_action, eTabColumnConfig_Match_GeneLocusTag);
          break;
        case eTableMatchProteinID:
          SetValue (dlg->column_action, eTabColumnConfig_Match_ProteinID);
          break;
        case eTableMatchDbxref:
          SetValue (dlg->column_action, eTabColumnConfig_Match_Dbxref);
          SetMatchLocationInIdMatchLocationDlg (dlg->match_location, f->match_type->match_location);
          break;
        case eTableMatchNucID:
          SetValue (dlg->column_action, eTabColumnConfig_Match_NucleotideID);
          break;
        case eTableMatchBioSource:
          SetValue (dlg->column_action, eTabColumnConfig_Match_Taxname);
          SetMatchLocationInIdMatchLocationDlg (dlg->match_location, f->match_type->match_location);
          break;
        case eTableMatchSourceQual:
          SetValue (dlg->column_action, eTabColumnConfig_Match_SourceQual);
          PointerToDialog (dlg->src_qual_match, f->match_type->data);
          SetMatchLocationInIdMatchLocationDlg (dlg->match_location, f->match_type->match_location);
          break;
      }
      PointerToDialog (dlg->constraint, NULL);
    } else if (f->field == NULL) {
      SetValue (dlg->column_action, eTabColumnConfig_Ignore);
      PointerToDialog (dlg->constraint, NULL);
    } else {
      if (f->field->choice == FieldType_source_qual) {
        src_field = f->field->data.ptrvalue;
        if (src_field != NULL 
            && src_field->choice == SourceQualChoice_textqual 
            && src_field->data.intvalue == Source_qual_taxname) {
          SetValue (dlg->column_action, eTabColumnConfig_Apply_Taxname);
        } else {
          SetValue (dlg->column_action, eTabColumnConfig_Apply_SourceQual);
          PointerToDialog (dlg->src_qual, f->field->data.ptrvalue);
        }
      } else if (f->field->choice == FieldType_feature_field) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_FeatureField);
        ff = (FeatureFieldPtr) f->field->data.ptrvalue;
        if (ff == NULL) {
          PointerToDialog (dlg->feature_type, NULL);
          PointerToDialog (dlg->feature_field, NULL);
        } else {
          vn.choice = (Uint1) ff->type;
          vn.data.ptrvalue = NULL;
          vn.next = NULL;
          PointerToDialog (dlg->feature_type, &vn);
          if (ff->field == NULL || ff->field->choice != FeatQualChoice_legal_qual) {
            PointerToDialog (dlg->feature_field, NULL);
          } else {
            vn.choice = ff->field->data.intvalue;
            vn.data.ptrvalue = NULL;
            vn.next = NULL;
            PointerToDialog (dlg->feature_field, &vn);
          }
        }          
      } else if (f->field->choice == FieldType_cds_gene_prot) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_CDSGeneProtField);
        PointerToDialog (dlg->cdsgeneprot, f->field);
      } else if (f->field->choice == FieldType_pub) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_PubField);
        vn.choice = f->field->data.intvalue;
        vn.data.ptrvalue = NULL;
        vn.next = NULL;      
        PointerToDialog (dlg->pub_field, &vn);
      } else if (f->field->choice == FieldType_misc && f->field->data.intvalue == Misc_field_genome_project_id) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_GenomeProjectId);
      } else if (f->field->choice == FieldType_misc && f->field->data.intvalue == Misc_field_comment_descriptor) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_CommentDescriptor);
      } else if (f->field->choice == FieldType_misc && f->field->data.intvalue == Misc_field_defline) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_Defline);
      } else if (f->field->choice == FieldType_misc && f->field->data.intvalue == Misc_field_keyword) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_Keyword);
      } else if (f->field->choice == FieldType_molinfo_field) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_Molinfo);
        PointerToDialog (dlg->molinfo_field, f->field->data.ptrvalue);
      } else if (f->field->choice == FieldType_struc_comment_field) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_StructuredCommentField);
        PointerToDialog (dlg->struccomm_field, f->field->data.ptrvalue);
      } else if (f->field->choice == FieldType_dblink) {
        SetValue (dlg->column_action, eTabColumnConfig_Apply_Dblink);
        SetValue (dlg->dblink_field, f->field->data.intvalue);
      }
        
      SafeSetStatus (dlg->change_mrna, f->match_mrna);
      SafeSetStatus (dlg->erase_when_blank, !f->skip_blank);
      SetExistingTextDialogValue(dlg->existing_text, f->existing_text);
      PointerToDialog (dlg->constraint, (Pointer) f->constraint);
    }
    f = TabColumnConfigFree (f);
  }
  TabColumnActionChange (dlg->column_action);
}


static void TabColumnConfigButtonChange (ButtoN b)
{
  TabColumnConfigDlgPtr dlg;

  dlg = (TabColumnConfigDlgPtr) GetObjectExtra (b);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void TabColumnConfigPopupChange (PopuP p)
{
  TabColumnConfigDlgPtr dlg;

  dlg = (TabColumnConfigDlgPtr) GetObjectExtra (p);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SetTabColumnConfigDialogTitle (DialoG d, CharPtr title, Int4 num_blank)
{
  CharPtr real_title;
  CharPtr title_fmt = "(%d are blank)%s";
  TabColumnConfigDlgPtr dlg;

  if (title == NULL) {
    real_title = "";
  } else if (num_blank > 0) {
    real_title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (title) + 15));
    sprintf (real_title, title_fmt, num_blank, title);
  } else {
    real_title = title;
  }

  SetTitle (d, real_title);

  if (real_title != title) {
    real_title = MemFree (real_title);
  }  

  dlg = (TabColumnConfigDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {  
    if (num_blank == 0) {
      SafeDisable (dlg->erase_when_blank);
    } else {
      SafeEnable (dlg->erase_when_blank);
    }
  }
}


NLM_EXTERN DialoG TabColumnConfigDialog 
(GrouP                    h,
 CharPtr                  title,
 Int4                     num_blank,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  TabColumnConfigDlgPtr dlg;
  GrouP p, k;
  CharPtr real_title;
  CharPtr title_fmt = "(%d are blank)%s";
  Boolean free_title = FALSE;

  dlg = (TabColumnConfigDlgPtr) MemNew (sizeof (TabColumnConfigDlgData));

  if (title == NULL) {
    real_title = "";
  } else if (num_blank > 0) {
    real_title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title_fmt) + StringLen (title) + 15));
    sprintf (real_title, title_fmt, num_blank, title);
    free_title = TRUE;
  } else {
    real_title = title;
  }

  p = NormalGroup (h, -1, 0, real_title, programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  if (free_title) {
    real_title = MemFree (real_title);
  }
  dlg->dialog = (DialoG) p;
  dlg->fromdialog = TabColumnConfigDialogToChoice;
  dlg->todialog = TabColumnConfigToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->column_action = PopupList (p, TRUE, TabColumnActionChange);
  SetObjectExtra (dlg->column_action, dlg, NULL);
  PopupItem (dlg->column_action, "Ignore column");
  PopupItem (dlg->column_action, "Match to Nucleotide ID");
  PopupItem (dlg->column_action, "Match to Taxname");
  PopupItem (dlg->column_action, "Match to Source Qual");
  PopupItem (dlg->column_action, "Match to Feature ID");
  PopupItem (dlg->column_action, "Match to Gene locus tag");
  PopupItem (dlg->column_action, "Match to Protein ID");
  PopupItem (dlg->column_action, "Match to Protein Name");   /* J. Chen */
  PopupItem (dlg->column_action, "Match to Dbxref");
  PopupItem (dlg->column_action, " ");
  PopupItem (dlg->column_action, "Apply to Taxname");
  PopupItem (dlg->column_action, "Apply to Source Qual");
  PopupItem (dlg->column_action, "Apply to CDS-Gene-Prot Field");
  PopupItem (dlg->column_action, "Apply to Feature Field");
  PopupItem (dlg->column_action, "Apply to Publication Field");
  PopupItem (dlg->column_action, "Apply to Genome Project ID");
  PopupItem (dlg->column_action, "Apply to Comment Descriptor");
  PopupItem (dlg->column_action, "Apply to Definition Line");
  PopupItem (dlg->column_action, "Apply to Keyword");
  PopupItem (dlg->column_action, "Apply to Molinfo");
  PopupItem (dlg->column_action, "Apply to Structured Comment Field");
  PopupItem (dlg->column_action, "Apply to Dblink");
  SetValue (dlg->column_action, 1);
  
  k = HiddenGroup (p, 0, 0, NULL);
  dlg->match_grp = HiddenGroup (k, 2, 0, NULL);
  dlg->match_location = IdMatchLocationDlg (dlg->match_grp, change_notify, change_userdata);
  dlg->src_qual_match = SourceQualChoiceDialog (dlg->match_grp, TRUE, FALSE, FALSE, change_notify, change_userdata);

  dlg->qual_grp = HiddenGroup (k, 0, 0, NULL);
  dlg->src_qual = SourceQualChoiceDialog (dlg->qual_grp, TRUE, FALSE, FALSE, TabColumnConfigFieldChange, dlg);
  dlg->feature_field_grp = HiddenGroup (dlg->qual_grp, 2, 0, NULL);
  dlg->feature_type = FeatureTypeDialog (dlg->feature_field_grp, TabColumnConfigFieldChange, dlg);
  dlg->feature_field = LegalFeatQualChoiceDialog (dlg->feature_field_grp, TabColumnConfigFieldChange, dlg);
  dlg->cdsgeneprot = CDSGeneProtFieldDialog (dlg->qual_grp, TabColumnConfigFieldChange, dlg);
  dlg->pub_field = PubFieldDialog (dlg->qual_grp, TabColumnConfigFieldChange, dlg);
  dlg->molinfo_field = MolinfoFieldChoiceDialog  (dlg->qual_grp, TabColumnConfigFieldChange, dlg);
  dlg->struccomm_field = StructuredCommentFieldDialog (dlg->qual_grp, TabColumnConfigFieldChange, dlg);
  dlg->dblink_field = MakeDbLinkFieldPopup (dlg->qual_grp, TabColumnConfigPopupChange, dlg);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->match_grp, (HANDLE) dlg->qual_grp, NULL);

  dlg->change_mrna = CheckBox (p, "Also change mRNA product name", TabColumnConfigButtonChange);
  SetObjectExtra (dlg->change_mrna, dlg, NULL);
  Disable (dlg->change_mrna);

  dlg->apply_options = HiddenGroup (p, -1, 0, NULL);
  if (num_blank > 0 || title == NULL) {
    dlg->erase_when_blank = CheckBox (dlg->apply_options, "Erase field when table cell is blank", TabColumnConfigButtonChange);
    SetObjectExtra (dlg->erase_when_blank, dlg, NULL);
    if (num_blank == 0) {
      Disable (dlg->erase_when_blank);
    }
  } else {
    dlg->erase_when_blank = NULL;
  }
  dlg->existing_text = ExistingTextDialog (dlg->apply_options, change_notify, change_userdata);
  dlg->constraint = ConstraintSetDialog (dlg->apply_options, change_notify, change_userdata);

  Disable (dlg->apply_options);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->existing_text,
                              (HANDLE) dlg->constraint,
                              (HANDLE) dlg->erase_when_blank, 
                              NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) dlg->apply_options, (HANDLE) dlg->change_mrna, NULL);

  TabColumnActionChange (dlg->column_action);

  return (DialoG) p;
}


static CharPtr SummarizeMatchType (MatchTypePtr match_type)
{
  CharPtr location_word = "Matches";
  CharPtr type_word = "feature ID";
  CharPtr match_fmt = "%s %s";
  CharPtr summ = NULL;
 
  if (match_type == NULL) {
    return NULL;
  }
  switch (match_type->match_location) {
    case String_location_contains :
      location_word = "Contained in";
      break;
    case String_location_equals :
      location_word = "Matches";
      break;
    case String_location_starts :
      location_word = "Is start of";
      break;
    case String_location_ends :
      location_word = "Is end of";
      break;
    case String_location_inlist :
      location_word = "List contains";
      break;
  }

  switch (match_type->choice) {
    case eTableMatchProteinName:     /* J. Chen */
      summ = StringSave ("Match to Protein Name");
      break;
    case eTableMatchFeatureID:
      summ = StringSave ("Match to feature ID");
      break;
    case eTableMatchGeneLocusTag:
      summ = StringSave ("Match to gene locus tag");
      break;
    case eTableMatchProteinID:
      summ = StringSave ("Match to protein ID");
      break;
    case eTableMatchDbxref:
      type_word = "feature dbxref";
      summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (match_fmt) + StringLen (type_word) + StringLen (location_word)));
      sprintf (summ, match_fmt, location_word, type_word);
      break;
    case eTableMatchNucID:
      summ = StringSave ("Match to nucleotide ID");
      break;
    case eTableMatchBioSource:
      type_word = "taxname";
      summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (match_fmt) + StringLen (type_word) + StringLen (location_word)));
      sprintf (summ, match_fmt, location_word, type_word);
      break;
    case eTableMatchSourceQual:
      if (match_type->data == NULL) {
        type_word = "unspecified source qualifier";
        summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (match_fmt) + StringLen (type_word) + StringLen (location_word)));
        sprintf (summ, match_fmt, location_word, type_word);
      } else {
        type_word = SummarizeSourceQual (match_type->data);
        summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (match_fmt) + StringLen (type_word) + StringLen (location_word)));
        sprintf (summ, match_fmt, location_word, type_word);
        type_word = MemFree (type_word);
      }
      break;
  }
  return summ;
}

static CharPtr SummarizeTabColumnConfig (TabColumnConfigPtr t)
{
  CharPtr summ = NULL;
  CharPtr apply_fmt = "Apply to %s %s%s(%s)";
  CharPtr field;
  CharPtr erase_when_blank = ", remove field when cell is blank";
  CharPtr also_mrna = ", make mRNA product match new CDS product";
  CharPtr existing_text;
  CharPtr constraint;
  Int4    summ_len = 0;

  if (t == NULL) return StringSave ("Ignore column");

  if (t->match_type == NULL) {
    field = SummarizeFieldType (t->field);
    if (field == NULL) {
      field = StringSave ("unspecified field");
    }
    existing_text = SummarizeExistingText (t->existing_text);
    constraint = SummarizeConstraintSet (t->constraint);
    summ_len = StringLen (apply_fmt) + StringLen (field) + StringLen (existing_text) + StringLen (constraint);
    if (!t->skip_blank) {
      summ_len += StringLen (erase_when_blank);
    }
    if (t->match_mrna) {
      summ_len += StringLen (also_mrna);
    }
    summ = (CharPtr) MemNew (sizeof (Char) * summ_len);
    sprintf (summ, apply_fmt, field, constraint == NULL ? "" : constraint, constraint == NULL ? "" : " ", existing_text);
    if (!t->skip_blank) {
      StringCat (summ, erase_when_blank);
    }
    if (t->match_mrna) {
      StringCat (summ, also_mrna);
    }
    field = MemFree (field);
  } else {
    summ = SummarizeMatchType (t->match_type);
  }
  return summ;
}


static void CleanupTabColumnConfigListDialog (GraphiC g, VoidPtr data)

{
  TabColumnConfigListDlgPtr dlg;
  Int4 i;

  dlg = (TabColumnConfigListDlgPtr) data;
  if (dlg != NULL) {
    for (i = 0; i < dlg->num_columns; i++) {
      dlg->column_list[i] = TabColumnConfigFree (dlg->column_list[i]);
      dlg->first_values[i] = MemFree (dlg->first_values[i]);
    }
    dlg->column_list = MemFree (dlg->column_list);
    dlg->first_values = MemFree (dlg->first_values);
    dlg->blank_list = MemFree (dlg->blank_list);
  }
  StdCleanupExtraProc (g, data);
}




static Pointer DialogToTabColumnConfigList (DialoG d)
{
  TabColumnConfigListDlgPtr dlg;
  Int4 i;
  ValNodePtr column_list = NULL;

  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;
  for (i = 0; i < dlg->num_columns; i++) {
    ValNodeAddPointer (&column_list, 0, TabColumnConfigCopy (dlg->column_list[i]));
  }
  return (Pointer) column_list;
}


static ValNodePtr TestTabColumnConfigListDialog (DialoG d)
{
  TabColumnConfigListDlgPtr dlg;
  Int4 i;
  ValNodePtr err_list = NULL;
  Int4 num_match = 0;
  Boolean have_apply = FALSE;

  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  for (i = 0; i < dlg->num_columns; i++) {
    if (dlg->column_list[i] != NULL) {
      if (dlg->column_list[i]->match_type != NULL) {
        if (dlg->column_list[i]->match_type->choice != eTableMatchSourceQual
            || dlg->column_list[i]->match_type->data != NULL) {
          num_match++;
        }
      } else if (!IsFieldTypeEmpty(dlg->column_list[i]->field)) {
        have_apply = TRUE;
      }
    }
  }
  if (num_match == 0) {
    ValNodeAddPointer (&err_list, 0, "No match column");
  } else if (num_match > 1) {
    ValNodeAddPointer (&err_list, 0, "Too many match columns");
  }
  if (!have_apply) {
    ValNodeAddPointer (&err_list, 0, "No apply column");
  }
  return err_list;
}


NLM_EXTERN void PopulateTabConfigListColumnListDoc (DialoG d)
{
  TabColumnConfigListDlgPtr dlg;
  Int4 i, len;
  CharPtr str, tmp;
  CharPtr    row_fmt = "%d\t%s\t%s\n";
  RecT       r;
  FonT       font = programFont;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {0, 0, 4, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE}, /* column number */
    {0, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE}, /* action */
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}   /* first value */
  };
  Int4 scroll_pos = 0;
  BaR  sb_vert;


  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  sb_vert = GetSlateVScrollBar ((SlatE) dlg->column_list_doc);
  if (sb_vert) {
    scroll_pos = GetBarValue (sb_vert);
  }

  Reset (dlg->column_list_doc);

  ObjectRect (dlg->column_list_doc, &r);
  InsetRect (&r, 4, 4);
  
  if (dlg->num_columns < 10) {    
    ColFmt[0].pixWidth = stdCharWidth * 2;
  } else {
    ColFmt[0].pixWidth = stdCharWidth * 3;
  }
  ColFmt[1].pixWidth = (r.right - r.left - ColFmt[0].pixWidth) / 2;
  ColFmt[2].pixWidth = ColFmt[1].pixWidth;
  for (i = 0; i < dlg->num_columns; i++) {
    str = SummarizeTabColumnConfig (dlg->column_list[i]);
    len = StringLen (str) + StringLen (dlg->first_values[i]) + 15 + StringLen (row_fmt);
    tmp = (CharPtr) MemNew (sizeof (Char) * len);
    sprintf (tmp, row_fmt, i + 1, str, dlg->first_values[i]);
    str = MemFree (str);
    AppendText (dlg->column_list_doc, tmp, &ParFmt, ColFmt, font);
    tmp = MemFree (tmp);
  }
  UpdateDocument (dlg->column_list_doc, 0, 0);
  if (scroll_pos > 0) {
    CorrectBarValue (sb_vert, scroll_pos);
  }
}


static void ClickColumnListDoc (DoC d, PoinT pt)
{
  Int2               item, row, col;
  RecT               rct;
  TabColumnConfigListDlgPtr dlg;
  TabColumnConfigPtr        data;
  
  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item > 0 && item <= dlg->num_columns) {
    dlg->current_column = item - 1;
    data = dlg->column_list[item - 1];
    PointerToDialog (dlg->edit_col_dlg, data);
    SetTabColumnConfigDialogTitle (dlg->edit_col_dlg, dlg->first_values[item - 1], dlg->blank_list[item - 1]);
  }
}

static Boolean ShowSelectedColumn (DoC doc, Int2 item, Int2 row, Int2 col)
{
  TabColumnConfigListDlgPtr dlg;
  
  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (doc);
  if (dlg == NULL) return FALSE;

  if (dlg->current_column == item - 1) {
    return TRUE;
  } else {
    return FALSE;
  }  
}


static void ChangeTabColumnConfig (Pointer userdata)
{
  TabColumnConfigListDlgPtr dlg;

  dlg = (TabColumnConfigListDlgPtr) userdata;
  if (dlg == NULL || dlg->column_list == NULL) return;

  dlg->column_list[dlg->current_column] = TabColumnConfigFree (dlg->column_list[dlg->current_column]);
  dlg->column_list[dlg->current_column] = DialogToPointer (dlg->edit_col_dlg);
  PopulateTabConfigListColumnListDoc (dlg->dialog);
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void TabColumnConfigListToDialog (DialoG d, Pointer data)
{
  TabColumnConfigListDlgPtr dlg;
  Int4 i;
  ValNodePtr vnp;

  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  vnp = (ValNodePtr) data;
  for (i = 0; i < dlg->num_columns; i++) {
    dlg->column_list[i] = TabColumnConfigFree (dlg->column_list[i]);
    if (vnp != NULL) {
      dlg->column_list[i] = TabColumnConfigCopy (vnp->data.ptrvalue);
      vnp = vnp->next;
    }
  }
  PopulateTabConfigListColumnListDoc (dlg->dialog);
  PointerToDialog (dlg->edit_col_dlg, dlg->column_list[dlg->current_column]);
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}




NLM_EXTERN DialoG TabColumnConfigListDialog (GrouP h, ValNodePtr first_values, ValNodePtr blank_list, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  TabColumnConfigListDlgPtr dlg;
  GrouP                p;
  ValNodePtr           vnp_v, vnp_b;
  Int4                 i;
  
  dlg = (TabColumnConfigListDlgPtr) MemNew (sizeof (TabColumnConfigListDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, CleanupTabColumnConfigListDialog);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = DialogToTabColumnConfigList;
  dlg->todialog = TabColumnConfigListToDialog;
  dlg->testdialog = TestTabColumnConfigListDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;


  dlg->column_list_doc = DocumentPanel (p, stdCharWidth * 50, stdLineHeight * 8);
  SetObjectExtra (dlg->column_list_doc, dlg, NULL);
  SetDocProcs (dlg->column_list_doc, ClickColumnListDoc, NULL, NULL, NULL);
  SetDocShade (dlg->column_list_doc, NULL, NULL, ShowSelectedColumn, NULL);

  dlg->edit_col_dlg = TabColumnConfigDialog (p, NULL, 0, ChangeTabColumnConfig, dlg);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->column_list_doc, (HANDLE) dlg->edit_col_dlg, NULL);

  /* populate column list */
  dlg->num_columns = ValNodeLen (blank_list);
  dlg->column_list = (TabColumnConfigPtr PNTR) MemNew (sizeof (TabColumnConfigPtr) * dlg->num_columns);
  MemSet (dlg->column_list, 0, sizeof (TabColumnConfigPtr) * dlg->num_columns);
  dlg->first_values = MemNew (sizeof (CharPtr) * dlg->num_columns);
  dlg->blank_list = (Int4Ptr) MemNew (sizeof (Int4) * dlg->num_columns);

  vnp_v = first_values;
  vnp_b = blank_list;
  i = 0;
  while (vnp_b != NULL) {
    dlg->blank_list[i] = vnp_b->data.intvalue;
    if (vnp_v == NULL || StringHasNoText (vnp_v->data.ptrvalue)) {
      dlg->first_values[i] = StringSave ("First row value is blank");
    } else {
      dlg->first_values[i] = StringSave (vnp_v->data.ptrvalue);
    }
    dlg->column_list[i] = NULL;
    vnp_b = vnp_b->next;
    if (vnp_v != NULL) {
      vnp_v = vnp_v->next;
    }
    i++;
  }
    
  PopulateTabConfigListColumnListDoc ((DialoG) p);

  dlg->current_column = 0;
  PointerToDialog (dlg->edit_col_dlg, dlg->column_list[0]);
  SetTabColumnConfigDialogTitle (dlg->edit_col_dlg, dlg->first_values[0], dlg->blank_list[0]);
    
  return (DialoG) p;
}


NLM_EXTERN void ChangeDataForTabColumnConfigListDialog (DialoG d, ValNodePtr first_values, ValNodePtr blank_list)
{
  TabColumnConfigListDlgPtr dlg;
  TabColumnConfigPtr PNTR column_list;
  Int4 num_columns, i;
  ValNodePtr           vnp_v, vnp_b;

  dlg = (TabColumnConfigListDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  num_columns = ValNodeLen (blank_list);
  if (num_columns < dlg->num_columns) {
    /* truncate existing list */
    for (i = num_columns; i < dlg->num_columns; i++) {
      dlg->column_list[i] = TabColumnConfigFree (dlg->column_list[i]);
      dlg->first_values[i] = MemFree (dlg->first_values[i]);
    }
    if (dlg->current_column >= num_columns) {
      dlg->current_column = num_columns - 1;
    }

    dlg->num_columns = num_columns;
  } else if (num_columns > dlg->num_columns) {
    /* need larger lists */
    for (i = 0; i < dlg->num_columns; i++) {
      dlg->first_values[i] = MemFree (dlg->first_values[i]);
    }
    dlg->first_values = MemFree (dlg->first_values);
    dlg->first_values = MemNew (sizeof (CharPtr) * num_columns);

    dlg->blank_list = MemFree (dlg->blank_list);
    dlg->blank_list = (Int4Ptr) MemNew (sizeof (Int4) * num_columns);

    column_list = (TabColumnConfigPtr PNTR) MemNew (sizeof (TabColumnConfigPtr) * num_columns);
    for (i = 0; i < dlg->num_columns; i++) {
      column_list[i] = dlg->column_list[i];
      dlg->column_list[i] = NULL;
    }
    dlg->column_list = MemFree (dlg->column_list);
    dlg->column_list = column_list;

    dlg->num_columns = num_columns;
  }

  /* populate column list */
  vnp_v = first_values;
  vnp_b = blank_list;
  i = 0;
  while (vnp_b != NULL) {
    dlg->blank_list[i] = vnp_b->data.intvalue;
    if (vnp_v == NULL || StringHasNoText (vnp_v->data.ptrvalue)) {
      dlg->first_values[i] = StringSave ("First row value is blank");
    } else {
      dlg->first_values[i] = StringSave (vnp_v->data.ptrvalue);
    }
    vnp_b = vnp_b->next;
    if (vnp_v != NULL) {
      vnp_v = vnp_v->next;
    }
    i++;
  }
    
  PopulateTabConfigListColumnListDoc (d);

  PointerToDialog (dlg->edit_col_dlg, dlg->column_list[dlg->current_column]);
  SetTabColumnConfigDialogTitle (dlg->edit_col_dlg, dlg->first_values[dlg->current_column], dlg->blank_list[dlg->current_column]);

}


typedef struct matchtypedlg {
  DIALOG_MESSAGE_BLOCK
  PopuP match_type;
  DialoG match_location;
  DialoG src_qual_match;
  
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} MatchTypeDlgData, PNTR MatchTypeDlgPtr;

typedef enum matchtypechoice {
  eMatchType_NucleotideID = 1,
  eMatchType_Taxname,
  eMatchType_SourceQual,
  eMatchType_FeatureID,
  eMatchType_GeneLocusTag,
  eMatchType_ProteinID,
  eMatchType_Dbxref,
  eMatchType_ProteinName
} EMatchTypeChoice;


static void ChangeMatchTypeChoice (PopuP p)
{
  MatchTypeDlgPtr dlg;
  Int2            val;

  dlg = (MatchTypeDlgPtr) GetObjectExtra (p);
  if (dlg != NULL) {
    val = GetValue (dlg->match_type);
    switch (val) {
      case eMatchType_Taxname:
      case eMatchType_Dbxref:
        Show (dlg->match_location);
        Hide (dlg->src_qual_match);
        break;
      case eMatchType_SourceQual:
        Show (dlg->match_location);
        Show (dlg->src_qual_match);
        break;
      case eMatchType_NucleotideID:
      case eMatchType_FeatureID:
      case eMatchType_GeneLocusTag:
      case eMatchType_ProteinID:
      default:
        Hide (dlg->match_location);
        Hide (dlg->src_qual_match);
        break;
    }
    if (dlg->change_notify != NULL) {
      (dlg->change_notify) (dlg->change_userdata);
    }
  }
}


static Pointer DialogToMatchType (DialoG d)
{
  MatchTypeDlgPtr dlg;
  MatchTypePtr match_type = NULL;
  Int2 val;

  dlg = (MatchTypeDlgPtr) GetObjectExtra (d);
  if (dlg != NULL) {
    val = GetValue (dlg->match_type);
    switch (val) {
      case eMatchType_ProteinName:    /* J. Chen */
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchProteinName;
        break;
      case eMatchType_FeatureID:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchFeatureID;
        break;
      case eMatchType_GeneLocusTag:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchGeneLocusTag;
        break;
      case eMatchType_ProteinID:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchProteinID;
        break;
      case eMatchType_Dbxref:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchDbxref;
        match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
      case eMatchType_NucleotideID:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchNucID;
        break;
      case eMatchType_Taxname:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchBioSource;
        match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
      case eMatchType_SourceQual:
        match_type = MatchTypeNew ();
        match_type->choice = eTableMatchSourceQual;
        match_type->data = DialogToPointer (dlg->src_qual_match);
        match_type->match_location = GetMatchLocationFromIdMatchLocationDlg (dlg->match_location);
        break;
    }
  }
  return match_type;
}


static ValNodePtr TestMatchTypeDialog (DialoG d)
{
  ValNodePtr err_list = NULL;
  MatchTypePtr match_type = NULL;

  match_type = DialogToPointer (d);
  if (match_type == NULL) {
    ValNodeAddPointer (&err_list, 0, "No match type");
  } 
  match_type = MatchTypeFree (match_type);
  return err_list;
}
  

NLM_EXTERN DialoG MatchTypeDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  MatchTypeDlgPtr dlg;
  GrouP           p, k;
  
  dlg = (MatchTypeDlgPtr) MemNew (sizeof (MatchTypeDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = DialogToMatchType;
  dlg->testdialog = TestMatchTypeDialog;

  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->match_type = PopupList (p, TRUE, ChangeMatchTypeChoice);
  SetObjectExtra (dlg->match_type, dlg, NULL);
  PopupItem (dlg->match_type, "Match to Nucleotide ID");
  PopupItem (dlg->match_type, "Match to Taxname");
  PopupItem (dlg->match_type, "Match to Source Qual");
  PopupItem (dlg->match_type, "Match to Feature ID");
  PopupItem (dlg->match_type, "Match to Gene locus tag");
  PopupItem (dlg->match_type, "Match to Protein ID");
  PopupItem (dlg->match_type, "Match to Protein Name");   /* J. Chen */
  PopupItem (dlg->match_type, "Match to Dbxref");
  SetValue (dlg->match_type, 1);
  
  k = HiddenGroup (p, 2, 0, NULL);
  dlg->match_location = IdMatchLocationDlg (k, change_notify, change_userdata);
  dlg->src_qual_match = SourceQualChoiceDialog (k, TRUE, FALSE, FALSE, change_notify, change_userdata);
  Hide (dlg->match_location);
  Hide (dlg->src_qual_match);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->match_type, (HANDLE) k, NULL);
  return (DialoG) p;
}

/* dialog for molinfo fields to be read from table */
typedef struct molinfofieldchoicedlg {
  DIALOG_MESSAGE_BLOCK
  PopuP molinfo_field;
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} MolinfoFieldChoiceDlgData, PNTR MolinfoFieldChoiceDlgPtr;


static Pointer MolinfoFieldChoiceDialogToChoice (DialoG d)
{
  MolinfoFieldChoiceDlgPtr dlg;
  ValNodePtr vnp;
  Int2       val;

  dlg = (MolinfoFieldChoiceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  vnp = ValNodeNew (NULL);
  val = GetValue (dlg->molinfo_field);
  switch (val) {
    case 1:
      vnp->choice = MolinfoField_molecule;
      break;
    case 2:
      vnp->choice = MolinfoField_technique;
      break;
    case 3:
      vnp->choice = MolinfoField_completedness;
      break;
    case 4:
      vnp->choice = MolinfoField_mol_class;
      break;
    case 5:
      vnp->choice = MolinfoField_topology;
      break;
    case 6:
      vnp->choice = MolinfoField_strand;
      break;
    default:
      vnp->choice = MolinfoField_molecule;
      break;
  }
  return (Pointer) vnp;
}


static void MolinfoFieldChoiceToDialog (DialoG d, Pointer data)
{
  MolinfoFieldChoiceDlgPtr dlg;
  ValNodePtr vnp;

  dlg = (MolinfoFieldChoiceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  vnp = (ValNodePtr) data;
  if (vnp == NULL) {
    SetValue (dlg->molinfo_field, 1);
  } else {
    switch (vnp->choice) {
      case MolinfoField_molecule:
        SetValue (dlg->molinfo_field, 1);
        break;
      case MolinfoField_technique:
        SetValue (dlg->molinfo_field, 2);
        break;
      case MolinfoField_completedness:
        SetValue (dlg->molinfo_field, 3);
        break;
      case MolinfoField_mol_class:
        SetValue (dlg->molinfo_field, 4);
        break;
      case MolinfoField_topology:
        SetValue (dlg->molinfo_field, 5);
        break;
      case MolinfoField_strand:
        SetValue (dlg->molinfo_field, 6);
        break;
      default:
        SetValue (dlg->molinfo_field, 1);
        break;
    }
  }
}


static void MolinfoFieldChoicePopupChange (PopuP p)
{
  MolinfoFieldChoiceDlgPtr dlg;

  dlg = (MolinfoFieldChoiceDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


NLM_EXTERN DialoG MolinfoFieldChoiceDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  MolinfoFieldChoiceDlgPtr dlg;
  GrouP p;

  dlg = (MolinfoFieldChoiceDlgPtr) MemNew (sizeof (MolinfoFieldChoiceDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = MolinfoFieldChoiceDialogToChoice;
  dlg->todialog = MolinfoFieldChoiceToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->molinfo_field = PopupList (p, TRUE, MolinfoFieldChoicePopupChange);
  SetObjectExtra (dlg->molinfo_field, dlg, NULL);
  PopupItem (dlg->molinfo_field, "molecule");
  PopupItem (dlg->molinfo_field, "technique");
  PopupItem (dlg->molinfo_field, "completedness");
  PopupItem (dlg->molinfo_field, "class");
  PopupItem (dlg->molinfo_field, "topology");
  PopupItem (dlg->molinfo_field, "strand");
  SetValue (dlg->molinfo_field, 1);

  return (DialoG) p;
}

/* suspect product rule editing*/

typedef Boolean (*IsObjectEmpty) PROTO ((Pointer));
typedef Pointer (*FreeObject) PROTO ((Pointer));
typedef DialoG  (*ObjEditDialog) PROTO ((GrouP, Nlm_ChangeNotifyProc, Pointer));

typedef struct editobject {
  DialoG dlg;
  ButtoN accept_btn;
  IsObjectEmpty empty_func;
  FreeObject free_func;
} EditObjectData, PNTR EditObjectPtr;


static void ChangeObject (Pointer data)
{
  EditObjectPtr e;
  Pointer obj;

  if ( (e = (EditObjectPtr) data) == NULL) {
    return;
  }

  /* if object is empty, disable accept button, otherwise enable it */
  obj = DialogToPointer (e->dlg);
  if (obj == NULL) {
    Disable (e->accept_btn);
  } else if (e->empty_func == NULL) {
    Enable (e->accept_btn);
  } else {
    if (e->empty_func(obj)) {
      Disable (e->accept_btn);
    } else {
      Enable (e->accept_btn);
    }
  }
  if (e->free_func == NULL) {
    obj = MemFree (obj);
  } else {
    obj = e->free_func(obj);
  }
}


static Pointer EditObject (Pointer orig, CharPtr title, IsObjectEmpty empty_func, FreeObject free_func, ObjEditDialog make_dlg)
{
  ModalAcceptCancelData acd;
  EditObjectData        esd;
  WindoW                w;
  ButtoN                b;
  GrouP                 h, c;
  Pointer               new_obj = NULL;
  
  w = MovableModalWindow(-20, -13, -10, -10, title, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  MemSet (&esd, 0, sizeof (EditObjectData));
  esd.empty_func = empty_func;
  esd.free_func = free_func;

  esd.dlg = make_dlg (h, ChangeObject, &esd);
  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  esd.accept_btn = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (esd.accept_btn, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) esd.dlg,
                              (HANDLE) c, 
                              NULL);

  PointerToDialog (esd.dlg, orig);

  ChangeObject (&esd);
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
  if (!acd.cancelled)
  {
    new_obj = DialogToPointer (esd.dlg);
  }
  Remove (w);
  return new_obj;
}


typedef struct searchfuncdlg {
  DIALOG_MESSAGE_BLOCK

  PopuP func_choice;
  DialoG string_constraint;
  TexT   text;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} SearchFuncDlgData, PNTR SearchFuncDlgPtr;


static void SearchFuncDlgPopupChange (PopuP p)
{
  SearchFuncDlgPtr dlg;
  Int2 val;

  dlg = (SearchFuncDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }

  val = GetValue (dlg->func_choice);
  switch (val) {
    case 1:
      Show (dlg->string_constraint);
      Hide (dlg->text);
      break;
    case 2:
    case 4:
    case 5:
    case 7:
    case 8:
      Hide (dlg->string_constraint);
      Hide (dlg->text);
      break;
    case 3:
    case 6:
    case 9:
    case 10:
      Hide (dlg->string_constraint);
      Show (dlg->text);
      break;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SearchFuncDlgTextChange (TexT t)
{
  SearchFuncDlgPtr dlg;

  dlg = (SearchFuncDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer SearchFuncFromDialog (DialoG d)
{
  SearchFuncDlgPtr dlg;
  Int2 val;
  CharPtr str;
  SearchFuncPtr func;

  dlg = (SearchFuncDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  func = ValNodeNew (NULL);
  val = GetValue (dlg->func_choice);
  switch (val) {
    case 1:
      func->choice = SearchFunc_string_constraint;
      func->data.ptrvalue = DialogToPointer (dlg->string_constraint);
      break;
    case 2:
      func->choice = SearchFunc_contains_plural;
      break;
    case 3:
      func->choice = SearchFunc_n_or_more_brackets_or_parentheses;
      str = SaveStringFromText (dlg->text);
      if (str == NULL) {
        func->data.intvalue = 0;
      } else {
        func->data.intvalue = atoi (str);
      }
      str = MemFree (str);
      break;
    case 4:
      func->choice = SearchFunc_three_numbers;
      break;
    case 5:
      func->choice = SearchFunc_underscore;
      break;
    case 6:
      func->choice = SearchFunc_prefix_and_numbers;
      func->data.ptrvalue = SaveStringFromText (dlg->text);
      break;
    case 7:
      func->choice = SearchFunc_all_caps;
      break;
    case 8:
      func->choice = SearchFunc_unbalanced_paren;
      break;
    case 9:
      func->choice = SearchFunc_too_long;
      str = SaveStringFromText (dlg->text);
      if (str == NULL) {
        func->data.intvalue = 0;
      } else {
        func->data.intvalue = atoi (str);
      }
      str = MemFree (str);
      break;
    case 10:
      func->choice = SearchFunc_has_term;
      func->data.ptrvalue = SaveStringFromText (dlg->text);
      break;
    default:
      func = ValNodeFree (func);
      break;
  }
  return (Pointer) func;
}


static void SearchFuncDialogChange (Pointer data)
{
  SearchFuncDlgPtr dlg;

  if ((dlg = (SearchFuncDlgPtr) data) != NULL 
    && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void SearchFuncToDialog (DialoG d, Pointer data)
{
  SearchFuncDlgPtr dlg;
  Char buf[15];
  SearchFuncPtr func;

  dlg = (SearchFuncDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  func = (SearchFuncPtr) data;
  if (func == NULL) {
    SetValue (dlg->func_choice, 1);
    PointerToDialog (dlg->string_constraint, NULL);
    return;
  }

  switch (func->choice) {
    case SearchFunc_string_constraint:
      SetValue (dlg->func_choice, 1);
      PointerToDialog (dlg->string_constraint, func->data.ptrvalue);
      break;
    case SearchFunc_contains_plural:
      SetValue (dlg->func_choice, 2);
      break;
    case SearchFunc_n_or_more_brackets_or_parentheses:
      SetValue (dlg->func_choice, 3);
      sprintf (buf, "%d", func->data.intvalue);
      SetTitle (dlg->text, buf);
      break;
    case SearchFunc_three_numbers:
      SetValue (dlg->func_choice, 4);
      break;
    case SearchFunc_underscore:
      SetValue (dlg->func_choice, 5);
      break;
    case SearchFunc_prefix_and_numbers:
      SetValue (dlg->func_choice, 6);
      SetTitle (dlg->text, func->data.ptrvalue == NULL ? "" : func->data.ptrvalue);
      break;
    case SearchFunc_all_caps:
      SetValue (dlg->func_choice, 7);
      break;
    case SearchFunc_unbalanced_paren:
      SetValue (dlg->func_choice, 8);
      break;
    case SearchFunc_too_long:
      SetValue (dlg->func_choice, 9);
      sprintf (buf, "%d", func->data.intvalue);
      SetTitle (dlg->text, buf);
      break;
    case SearchFunc_has_term:
      SetValue (dlg->func_choice, 10);
      SetTitle (dlg->text, func->data.ptrvalue == NULL ? "" : func->data.ptrvalue);
      break;
    default:
      SetValue (dlg->func_choice, 1);
      PointerToDialog (dlg->string_constraint, NULL);
      break;
  }
  SearchFuncDlgPopupChange(dlg->func_choice);
}



NLM_EXTERN DialoG SearchFuncDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  SearchFuncDlgPtr dlg;
  GrouP p, g;

  dlg = (SearchFuncDlgPtr) MemNew (sizeof (SearchFuncDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = SearchFuncFromDialog;
  dlg->todialog = SearchFuncToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->func_choice = PopupList (p, TRUE, SearchFuncDlgPopupChange);
  SetObjectExtra (dlg->func_choice, dlg, NULL);
  PopupItem (dlg->func_choice, "String Constraint");
  PopupItem (dlg->func_choice, "Contains plural");
  PopupItem (dlg->func_choice, "N or more brackets or parentheses");
  PopupItem (dlg->func_choice, "Three numbers");
  PopupItem (dlg->func_choice, "Contains underscore");
  PopupItem (dlg->func_choice, "Is prefix and numbers");
  PopupItem (dlg->func_choice, "Is all caps");
  PopupItem (dlg->func_choice, "Contains unbalanced parentheses");
  PopupItem (dlg->func_choice, "Is too long");
  PopupItem (dlg->func_choice, "Contains special term");
  SetValue (dlg->func_choice, 1);

  g = HiddenGroup (p, 0, 0, NULL);
  /* note  - need to also have string constraint list dialog */
  dlg->text = DialogText (g, "", 10, SearchFuncDlgTextChange);
  SetObjectExtra (dlg->text, dlg, NULL);
  Hide (dlg->text);
  dlg->string_constraint = StringConstraintDialog (g, NULL, FALSE, SearchFuncDialogChange, dlg);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->text, (HANDLE) dlg->string_constraint, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->func_choice, (HANDLE) g, NULL);

  return (DialoG) p;
}


typedef struct simplereplacedlg {
  DIALOG_MESSAGE_BLOCK
  TexT  replace;
  ButtoN whole_string;
  ButtoN weasel_to_putative;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} SimpleReplaceDlgData, PNTR SimpleReplaceDlgPtr;


static void ChangeSimpleReplaceText (TexT t)
{
  SimpleReplaceDlgPtr dlg;

  dlg = (SimpleReplaceDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void ChangeSimpleReplaceBtn (ButtoN b)
{
  SimpleReplaceDlgPtr dlg;

  dlg = (SimpleReplaceDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static Pointer SimpleReplaceFromDialog (DialoG d)
{
  SimpleReplaceDlgPtr dlg;
  SimpleReplacePtr    simple_replace;

  dlg = (SimpleReplaceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  simple_replace = SimpleReplaceNew();
  simple_replace->replace = JustSaveStringFromText (dlg->replace);
  simple_replace->whole_string = GetStatus (dlg->whole_string);
  simple_replace->weasel_to_putative = GetStatus (dlg->weasel_to_putative);

  return simple_replace;
}


static void SimpleReplaceToDialog (DialoG d, Pointer data)
{
  SimpleReplaceDlgPtr dlg;
  SimpleReplacePtr    simple_replace;

  dlg = (SimpleReplaceDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  simple_replace = (SimpleReplacePtr) data;
  if (simple_replace == NULL) {
    SetTitle (dlg->replace, "");
    SetStatus (dlg->whole_string, FALSE);
    SetStatus (dlg->weasel_to_putative, FALSE);
  } else {
    SetTitle (dlg->replace, simple_replace->replace == NULL ? "" : simple_replace->replace);
    SetStatus (dlg->whole_string, simple_replace->whole_string);
    SetStatus (dlg->weasel_to_putative, simple_replace->weasel_to_putative);
  }
}


NLM_EXTERN DialoG SimpleReplaceDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  SimpleReplaceDlgPtr dlg;
  GrouP p, g;

  dlg = (SimpleReplaceDlgPtr) MemNew (sizeof (SimpleReplaceDlgData));

  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = SimpleReplaceFromDialog;
  dlg->todialog = SimpleReplaceToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->replace = DialogText (p, "", 15, ChangeSimpleReplaceText);
  SetObjectExtra (dlg->replace, dlg, NULL);
  g = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  dlg->whole_string = CheckBox (g, "Replace entire string", ChangeSimpleReplaceBtn);
  SetObjectExtra (dlg->whole_string, dlg, NULL);
  dlg->weasel_to_putative = CheckBox (g, "Retain and normalize 'putative' synonym", ChangeSimpleReplaceBtn);
  SetObjectExtra (dlg->weasel_to_putative, dlg, NULL);

  return (DialoG) p;
}


typedef struct replacefuncdlg {
  DIALOG_MESSAGE_BLOCK

  PopuP func_type;
  DialoG simple_replace;
  TexT   text;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} ReplaceFuncDlgData, PNTR ReplaceFuncDlgPtr;


static void ChangeReplaceFuncPopup (PopuP p)
{
  ReplaceFuncDlgPtr dlg;
  Int2 val;

  dlg = (ReplaceFuncDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }
  
  val = GetValue (dlg->func_type);
  switch (val) {
    case 2:
      Show (dlg->simple_replace);
      Hide (dlg->text);
      break;
    case 3:
      Hide (dlg->simple_replace);
      Show (dlg->text);
      break;
    case 4:
      Hide (dlg->simple_replace);
      Hide (dlg->text);
      break;
    default:
      Hide (dlg->text);
      Hide (dlg->simple_replace);
      break;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ChangeReplaceFuncText (TexT t)
{
  ReplaceFuncDlgPtr dlg;

  dlg = (ReplaceFuncDlgPtr) GetObjectExtra (t);
  if (dlg == NULL) {
    return;
  }

  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static Pointer ReplaceFuncFromDialog (DialoG d)
{
  ReplaceFuncDlgPtr dlg;
  ReplaceFuncPtr    replace_func;
  SimpleReplacePtr  simple_replace;
  Int2              val;

  dlg = (ReplaceFuncDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  replace_func = ValNodeNew (NULL);
  val = GetValue (dlg->func_type);
  switch (val) {
    case 2:
      replace_func->choice = ReplaceFunc_simple_replace;
      replace_func->data.ptrvalue = DialogToPointer (dlg->simple_replace);
      break;
    case 3:
      replace_func->choice = ReplaceFunc_haem_replace;
      replace_func->data.ptrvalue = SaveStringFromText (dlg->text);
      break;
    case 4:
      replace_func->choice = ReplaceFunc_simple_replace;
      simple_replace = SimpleReplaceNew();
      simple_replace->replace = StringSave ("hypothetical protein");
      simple_replace->whole_string = TRUE;
      simple_replace->weasel_to_putative = FALSE;
      replace_func->data.ptrvalue = simple_replace;
      break;
    default:
      replace_func = ValNodeFree (replace_func);
      break;
  }
  /* note - need to return null if it's empty */
  return replace_func;
}


static Boolean IsHypotheticalSimpleReplace (SimpleReplacePtr simple_replace)
{
  if (simple_replace == NULL) {
    return FALSE;
  } else if (!simple_replace->whole_string || simple_replace->weasel_to_putative) {
    return FALSE;
  } else if (StringCmp (simple_replace->replace, "hypothetical protein") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void ReplaceFuncToDialog (DialoG d, Pointer data)
{
  ReplaceFuncDlgPtr dlg;
  ReplaceFuncPtr replace_func;

  dlg = (ReplaceFuncDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  replace_func = (ReplaceFuncPtr) data;
  if (replace_func == NULL) {
    SetValue (dlg->func_type, 1);
  } else {
    switch (replace_func->choice) {
      case ReplaceFunc_simple_replace:
        if (IsHypotheticalSimpleReplace(replace_func->data.ptrvalue)) {
          SetValue (dlg->func_type, 4);
        } else {
          SetValue (dlg->func_type, 2);
          PointerToDialog (dlg->simple_replace, replace_func->data.ptrvalue);
        }
        break;
      case ReplaceFunc_haem_replace:
        SetValue (dlg->func_type, 3);
        SetTitle (dlg->text, replace_func->data.ptrvalue == NULL ? "" : replace_func->data.ptrvalue);
        break;
      default:
        SetValue (dlg->func_type, 1);
        break;
    }
  }
  ChangeReplaceFuncPopup (dlg->func_type);
}


static DialoG ReplaceFuncDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ReplaceFuncDlgPtr dlg;
  GrouP p, g;

  dlg = (ReplaceFuncDlgPtr) MemNew (sizeof (ReplaceFuncDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = ReplaceFuncFromDialog;
  dlg->todialog = ReplaceFuncToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->func_type = PopupList (p, TRUE, ChangeReplaceFuncPopup);
  SetObjectExtra (dlg->func_type, dlg, NULL);
  PopupItem (dlg->func_type, "None");
  PopupItem (dlg->func_type, "Simple");
  PopupItem (dlg->func_type, "Haem");
  PopupItem (dlg->func_type, "Hypothetical");
  SetValue (dlg->func_type, 1);
 
  g = HiddenGroup (p, 0, 0, NULL);
  dlg->simple_replace = SimpleReplaceDialog (g, change_notify, change_userdata);
  Hide (dlg->simple_replace);
  dlg->text = DialogText (g, "", 15, ChangeReplaceFuncText);
  Hide (dlg->text);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->simple_replace, (HANDLE) dlg->text, NULL);

  return (DialoG) p;
}


typedef struct replaceruledlg {
  DIALOG_MESSAGE_BLOCK

  DialoG replace_func;
  ButtoN move_to_note;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} ReplaceRuleDlgData, PNTR ReplaceRuleDlgPtr;


static void ChangeReplaceRuleBtn (ButtoN b)
{
  ReplaceRuleDlgPtr dlg;

  dlg = (ReplaceRuleDlgPtr) GetObjectExtra (b);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static void ReplaceRuleToDialog (DialoG d, Pointer data)
{
  ReplaceRuleDlgPtr dlg;
  ReplaceRulePtr rule;

  dlg = (ReplaceRuleDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }
  rule = (ReplaceRulePtr) data;
  if (rule == NULL) {
    PointerToDialog (dlg->replace_func, NULL);
    SetStatus (dlg->move_to_note, FALSE);
  } else {
    PointerToDialog (dlg->replace_func, rule->replace_func);
    SetStatus (dlg->move_to_note, rule->move_to_note);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static Pointer ReplaceRuleFromDialog (DialoG d)
{
  ReplaceRuleDlgPtr dlg;
  ReplaceRulePtr rule;

  dlg = (ReplaceRuleDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }
  rule = ReplaceRuleNew();
  rule->replace_func =  DialogToPointer (dlg->replace_func);
  rule->move_to_note = GetStatus (dlg->move_to_note);
  if (rule->replace_func == NULL && !rule->move_to_note) {
    rule = ReplaceRuleFree (rule);
  }
  return rule;
}


static void ChangeReplaceRuleFunc (Pointer data)
{
  ReplaceRuleDlgPtr dlg;
  ValNodePtr replace_func;

  dlg = (ReplaceRuleDlgPtr) data;
  if (dlg == NULL) {
    return;
  }

  replace_func = DialogToPointer (dlg->replace_func);
  if (replace_func != NULL 
      && replace_func->choice == ReplaceFunc_simple_replace
      && IsHypotheticalSimpleReplace (replace_func->data.ptrvalue)) {
    SetStatus (dlg->move_to_note, TRUE);
  }
  
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static DialoG ReplaceRuleDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  ReplaceRuleDlgPtr dlg;
  GrouP p, g;

  dlg = (ReplaceRuleDlgPtr) MemNew (sizeof (ReplaceRuleDlgData));

  p = NormalGroup (h, -1, 0, "Replacement Action", programFont, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = ReplaceRuleFromDialog;
  dlg->todialog = ReplaceRuleToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  g = HiddenGroup (p, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Rule", 0, popupMenuHeight, programFont, 'r');
  dlg->replace_func = ReplaceFuncDialog (g, ChangeReplaceRuleFunc, dlg);
  dlg->move_to_note = CheckBox (p, "Move original to note", ChangeReplaceRuleBtn);
  SetObjectExtra (dlg->move_to_note, dlg, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) dlg->move_to_note, NULL);


  return (DialoG) p;
}


typedef struct suspectruledlg {
  DIALOG_MESSAGE_BLOCK

  PopuP rule_type;
  DialoG find;
  DialoG except;
  DialoG replace;
  DialoG feat_constraint;
  TexT   description;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} SuspectRuleDlgData, PNTR SuspectRuleDlgPtr;


static void SuspectRuleDlgPopupChange (PopuP p)
{
  SuspectRuleDlgPtr dlg;

  dlg = (SuspectRuleDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) {
    return;
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify)(dlg->change_userdata);
  }
}


static Pointer SuspectRuleFromDialog (DialoG d)
{
  SuspectRuleDlgPtr dlg;
  SuspectRulePtr new_rule;

  dlg = (SuspectRuleDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return NULL;
  }

  new_rule = SuspectRuleNew();
  new_rule->rule_type = GetValue (dlg->rule_type) - 1;
  new_rule->find = DialogToPointer (dlg->find);
  new_rule->except = DialogToPointer (dlg->except);
  if (IsSearchFuncEmpty(new_rule->except)) {
    new_rule->except = SearchFuncFree (new_rule->except);
  }
  new_rule->replace = DialogToPointer (dlg->replace);
  new_rule->feat_constraint = DialogToPointer (dlg->feat_constraint);
  new_rule->description = SaveStringFromText (dlg->description);
  return new_rule;
}


static void SuspectRuleToDialog (DialoG d, Pointer data)
{
  SuspectRuleDlgPtr dlg;
  SuspectRulePtr new_rule;
  SearchFuncPtr func;
  StringConstraintPtr scp;

  dlg = (SuspectRuleDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) {
    return;
  }

  new_rule = (SuspectRulePtr) data;
  if (new_rule == NULL) {
    SetValue (dlg->rule_type, 1);
    scp = StringConstraintNew();
    scp->case_sensitive = FALSE;
    scp->ignore_weasel = TRUE;
    func = ValNodeNew (NULL);
    func->choice = SearchFunc_string_constraint;
    func->data.ptrvalue = scp;
    PointerToDialog (dlg->find, func);
    func = SearchFuncFree (func);
    PointerToDialog (dlg->except, NULL);
    PointerToDialog (dlg->replace, NULL);
    PointerToDialog (dlg->feat_constraint, NULL);
    SetTitle (dlg->description, "");
  } else {
    SetValue (dlg->rule_type, new_rule->rule_type + 1);
    PointerToDialog (dlg->find, new_rule->find);
    PointerToDialog (dlg->except, new_rule->except);
    PointerToDialog (dlg->replace, new_rule->replace);
    PointerToDialog (dlg->feat_constraint, new_rule->feat_constraint);
    SetTitle (dlg->description, new_rule->description);
  }
}


NLM_EXTERN DialoG SuspectRuleDialog 
(GrouP                    h,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  SuspectRuleDlgPtr dlg;
  GrouP p;
  Int2 val;
  GrouP g1, g2;
  PrompT ppt;

  dlg = (SuspectRuleDlgPtr) MemNew (sizeof (SuspectRuleDlgData));

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->fromdialog = SuspectRuleFromDialog;
  dlg->todialog = SuspectRuleToDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  dlg->rule_type = PopupList (p, TRUE, SuspectRuleDlgPopupChange);
  SetObjectExtra (dlg->rule_type, dlg, NULL);
  for (val = 0; val <= Fix_type_gene; val++) {
    PopupItem (dlg->rule_type, SummarizeFixType(val));
  }
  SetValue (dlg->rule_type, 1);

  g1 = NormalGroup (p, 0, 0, "Where Product Matches", programFont, NULL);
  dlg->find = SearchFuncDialog (g1, change_notify, change_userdata);
  g2 = NormalGroup (p, 0, 0, "But Product does not Match", programFont, NULL);
  dlg->except = SearchFuncDialog (g2, change_notify, change_userdata);
  dlg->replace = ReplaceRuleDialog (p, change_notify, change_userdata);
  dlg->feat_constraint = ConstraintSetDialog (p, dlg->change_notify, dlg->change_userdata);
  ppt = StaticPrompt (p, "Description", 0, popupMenuHeight, programFont, 'c');
  dlg->description = ScrollText (p, 20, 3, programFont, TRUE, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->rule_type, (HANDLE) g1, (HANDLE) g2, (HANDLE) dlg->replace, 
                              (HANDLE) dlg->feat_constraint, 
                              (HANDLE) ppt,
                              (HANDLE) dlg->description,
                              NULL);

  return (DialoG) p;
}


static Pointer EditObjectFreeSuspectRule (Pointer rule)
{
  return (Pointer) SuspectRuleFree ((SuspectRulePtr)rule);
}


static SuspectRulePtr EditSuspectRule (SuspectRulePtr orig)
{
  return EditObject (orig, "Suspect Rule", (IsObjectEmpty) IsSuspectRuleEmpty, EditObjectFreeSuspectRule, SuspectRuleDialog);
}


/* suspect product rule editor */
typedef struct suspectproductruleeditorform {
  FORM_MESSAGE_BLOCK
  DoC    rule_summary;
  
  SuspectRuleSetPtr rule_list;
  CharPtr    last_filename;
  Boolean    unsaved;
  FonT       summary_font;
} SuspectProductRuleEditorFormData, PNTR SuspectProductRuleEditorFormPtr;

static void CleanupSuspectProductRuleEditorForm (GraphiC g, VoidPtr data)

{
  SuspectProductRuleEditorFormPtr f;

  f = (SuspectProductRuleEditorFormPtr) data;
  if (f != NULL) {
    f->rule_list = SuspectRuleSetFree (f->rule_list);
    f->last_filename = MemFree (f->last_filename);
  }
  StdCleanupFormProc (g, data);
}


static void SetupSuspectProductRuleEditorFont (SuspectProductRuleEditorFormPtr f)

{
  if (f == NULL) return;

#ifdef WIN_MAC
  f->summary_font = ParseFont ("Times,12");
#endif
#ifdef WIN_MSWIN
  f->summary_font = ParseFont ("Times New Roman,12");
#endif
#ifdef WIN_MOTIF
  f->summary_font = ParseFont ("Times,12");
#endif
}


static void SummarizeSuspectProductRule (DoC doc, SuspectRuleSetPtr rule_list, FonT font)
{
  SuspectRulePtr rule;
  CharPtr    str;
  CharPtr    tmp;
  Int4       pos = 1;
  RecT       r;
  ParData    ParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
  ColData    ColFmt[] = 
  {
    {kNumberColumnWidth, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {12, 0, 1, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, FALSE},
    {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
  };


  Reset (doc);

  ObjectRect (doc, &r);
  InsetRect (&r, 4, 4);
  
  ColFmt[5].pixWidth = stdCharWidth;
  ColFmt[6].pixWidth = r.right - r.left - stdCharWidth - 48 - kNumberColumnWidth;

  if (font == NULL) font = programFont;

  if (rule_list == NULL) {
    AppendText (doc, "(Click here to start a new rule list)", NULL, NULL, font);
  } else {
    AppendText (doc, "(Click here to insert a rule at the beginning of the list)", NULL, NULL, font);
    for (rule = rule_list; rule != NULL; rule = rule->next) {
      str = SummarizeSuspectRule (rule);
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 25));
      sprintf (tmp, "%d\t\t\t\t\tC\t%s\n", pos, str);
      str = MemFree (str);
      pos++;
      AppendText (doc, tmp, &ParFmt, ColFmt, font);
      tmp = MemFree (tmp);
    }
    AppendText (doc, "(Click here to insert a rule at the end of the list)", NULL, NULL, font);
  }
  UpdateDocument (doc, 0, 0);
}


static void SetSuspectProductRuleEditorFileItems (SuspectProductRuleEditorFormPtr frm)

{
  IteM  i;

  if (frm != NULL) {
    i = FindFormMenuItem ((BaseFormPtr) frm, VIB_MSG_OPEN);
    SafeSetTitle (i, "Load Rule List");
    SafeEnable (i);

    i = FindFormMenuItem ((BaseFormPtr) frm, VIB_MSG_IMPORT);
    SafeSetTitle (i, "Add Rules from File to List");
    SafeEnable (i);

    i = FindFormMenuItem ((BaseFormPtr) frm, VIB_MSG_SAVE);
    SafeSetTitle (i, "Save Rule List");
  }
}


static Boolean OpenSuspectProductRuleFile (ForM f, CharPtr filename)

{
  SuspectProductRuleEditorFormPtr mefp;
  Boolean            rval = FALSE;
  CharPtr            extension;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;
  SuspectRuleSetPtr  rule_list;

  mefp = (SuspectProductRuleEditorFormPtr) GetObjectExtra (f);
  if (mefp != NULL) {
    if (mefp->unsaved) {
      if (Message (MSG_YN, "Do you want to save changes to the current file?") == ANS_YES) {
        if (!SaveMacroFile(f, mefp->last_filename)) {
           return FALSE;
        }
      }
    }
    path [0] = '\0';
    StringNCpy_0 (path, filename, sizeof (path));
    extension = NULL;
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip == NULL) {
        Message (MSG_ERROR, "Unable to open %s", path);
      } else {
        rule_list = SuspectRuleSetAsnRead (aip, NULL);
        if (rule_list == NULL) {
          Message (MSG_ERROR, "Unable to read rule list from %s.", path);
        } else {
          mefp->rule_list = SuspectRuleSetFree (mefp->rule_list);
          mefp->rule_list = rule_list;
          mefp->last_filename = MemFree (mefp->last_filename);
          mefp->last_filename = StringSave (path);
          mefp->unsaved = FALSE;
          rval = TRUE;
          SetSuspectProductRuleEditorFileItems (mefp);
          SummarizeSuspectProductRule (mefp->rule_summary, mefp->rule_list, mefp->summary_font);
        }
        AsnIoClose (aip);
      }
    }    
  }
  return rval;
}


static Boolean SaveSuspectProductRuleFile (ForM f, CharPtr filename)

{
  SuspectProductRuleEditorFormPtr mefp;
  Boolean            rval = FALSE;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;

  mefp = (SuspectProductRuleEditorFormPtr) GetObjectExtra (f);
  if (mefp != NULL) {
    path [0] = '\0';
    StringNCpy_0 (path, filename, sizeof (path));
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
      aip = AsnIoOpen (path, "w");
      if (aip == NULL) {
        Message (MSG_ERROR, "Unable to open %s", path);
      } else {
        SuspectRuleSetAsnWrite (mefp->rule_list, aip, NULL);
        AsnIoClose (aip);
        mefp->last_filename = MemFree (mefp->last_filename);
        mefp->last_filename = StringSave (path);
        mefp->unsaved = FALSE;
        rval = TRUE;
      }
    }
  }
  return rval;
}


static Boolean ImportSuspectProductRuleFile (ForM f, CharPtr filename)

{
  SuspectProductRuleEditorFormPtr mefp;
  Boolean            rval = FALSE;
  CharPtr            extension;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;
  SuspectRuleSetPtr  rule_list, last;

  mefp = (SuspectProductRuleEditorFormPtr) GetObjectExtra (f);
  if (mefp == NULL) return FALSE;
  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  extension = NULL;
  if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
    aip = AsnIoOpen (path, "r");
    if (aip == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
    } else {
      rule_list = SuspectRuleSetAsnRead (aip, NULL);
      if (rule_list == NULL) {
        Message (MSG_ERROR, "Unable to read rule list from %s.", path);
      } else {
        if (mefp->rule_list == NULL) {
          mefp->rule_list = rule_list;
        } else {
          last = mefp->rule_list;
          while (last->next != NULL) {
            last = last->next;
          }
          last->next = rule_list;
        }
        mefp->unsaved = TRUE;
        rval = TRUE;
        SetSuspectProductRuleEditorFileItems (mefp);
        SummarizeSuspectProductRule (mefp->rule_summary, mefp->rule_list, mefp->summary_font);
      }
      AsnIoClose (aip);
    }
  }
  return rval;
}


static void SuspectProductRuleEditorFormMessage (ForM f, Int2 mssg)

{
  SuspectProductRuleEditorFormPtr  mefp;

  mefp = (SuspectProductRuleEditorFormPtr) GetObjectExtra (f);
  if (mefp != NULL) {
    switch (mssg) {
      case VIB_MSG_OPEN :
        OpenSuspectProductRuleFile (f, NULL);
        break;
      case VIB_MSG_IMPORT :
        ImportSuspectProductRuleFile (f, NULL);
        break;
      case VIB_MSG_SAVE :
        SaveSuspectProductRuleFile (f, mefp->last_filename);
        break;
      case VIB_MSG_SAVE_AS :
        SaveSuspectProductRuleFile (f, NULL);
        break;
      case VIB_MSG_CLOSE:
      case VIB_MSG_QUIT:
        if (mefp->unsaved) {
          if (Message (MSG_YN, "Do you want to save changes to the current file?") == ANS_YES) {
            if (SaveSuspectProductRuleFile(f, mefp->last_filename)) {
              Remove (mefp->form);
            }
          }
        }
        Remove (mefp->form);
        break;
      case NUM_VIB_MSG + 2:
        /* sort by find */
        SortSuspectRuleSetByFind(&(mefp->rule_list));
        SetSuspectProductRuleEditorFileItems (mefp);
        SummarizeSuspectProductRule (mefp->rule_summary, mefp->rule_list, mefp->summary_font);
        break;
      case NUM_VIB_MSG + 3:
        /* sort by fix-type, then find */
        SortSuspectRuleSetByFixTypeThenFind(&(mefp->rule_list));
        SetSuspectProductRuleEditorFileItems (mefp);
        SummarizeSuspectProductRule (mefp->rule_summary, mefp->rule_list, mefp->summary_font);
        break;
      default :
        break;
    }
  }
}


static void UpdateSuspectProductRuleSummary (SuspectProductRuleEditorFormPtr f, Int4 scroll_pos)
{
  Int4 scroll_max;
  BaR  sb_vert;

  if (f == NULL) return;

  f->unsaved = TRUE;
  SetSuspectProductRuleEditorFileItems (f);
  SummarizeSuspectProductRule (f->rule_summary, f->rule_list, f->summary_font);
  if (scroll_pos > 0) {
    sb_vert = GetSlateVScrollBar ((SlatE) f->rule_summary);
    scroll_max = GetBarMax (sb_vert);
    if (scroll_pos > scroll_max) {
      scroll_pos = scroll_max;
    }
    CorrectBarValue (sb_vert, scroll_pos);
  }
}


static Boolean CloneSuspectProductRuleItem (SuspectProductRuleEditorFormPtr f, Int2 item)
{
  SuspectRulePtr        this_rule = NULL, new_rule;
  Int2                  pos;
  Int4                  scroll_pos;
  BaR                   sb_vert;

  if (f == NULL || item == 0 || f->rule_list == NULL) {
    return FALSE;
  }

  pos = 1;
  this_rule = f->rule_list;
  while (pos < item && this_rule != NULL) {
    pos++;
    this_rule = this_rule->next;
  }
  if (this_rule == NULL) {
    return FALSE;
  }
  new_rule = AsnIoMemCopy (this_rule, (AsnReadFunc) SuspectRuleAsnRead, (AsnWriteFunc) SuspectRuleAsnWrite);
  if (new_rule != NULL) {
    new_rule->next = this_rule->next;
    this_rule->next = new_rule;
  }

  /* get current scroll position */
  sb_vert = GetSlateVScrollBar ((SlatE) f->rule_summary);
  scroll_pos = GetBarValue (sb_vert);
  /* we will want to increase the scroll position after each addition
   * note that we need to get the scroll bar and check the initial position
   * each time - if there was no scroll bar after the last update, scroll_pos
   * needs to be zero to start.
   */
  scroll_pos++;

  /* update summary */
  UpdateSuspectProductRuleSummary (f, scroll_pos);
  return TRUE;
}


static void AddSuspectProductRule (SuspectProductRuleEditorFormPtr f, Int4 item)
{
  SuspectRulePtr new_rule, prev_rule = NULL;
  Int2                  pos;
  Int4                  scroll_pos;
  BaR                   sb_vert;

  if (f == NULL) {
    return;
  }

  new_rule = EditSuspectRule(NULL);
  if (new_rule == NULL) {
    return;
  }

  if (f->rule_list == NULL || item == 0) {
    new_rule->next = f->rule_list;
    f->rule_list = new_rule;
  } else if (item > 0) {
    pos = 1;
    prev_rule = f->rule_list;
    while (pos < item && prev_rule->next != NULL) {
      pos++;
      prev_rule = prev_rule->next;
    }
    if (prev_rule == NULL) {
      new_rule->next = f->rule_list;
      f->rule_list = new_rule;
    } else {
      new_rule->next = prev_rule->next;
      prev_rule->next = new_rule;
    }
  }  

  /* get current scroll position */
  sb_vert = GetSlateVScrollBar ((SlatE) f->rule_summary);
  scroll_pos = GetBarValue (sb_vert);
  /* we will want to increase the scroll position after each addition
   * note that we need to get the scroll bar and check the initial position
   * each time - if there was no scroll bar after the last update, scroll_pos
   * needs to be zero to start.
   */
  scroll_pos++;

  UpdateSuspectProductRuleSummary (f, scroll_pos);
  f->unsaved = TRUE;
}


static Boolean EditSuspectProductRuleItem (SuspectProductRuleEditorFormPtr f, Int2 item)
{
  Boolean rval = FALSE;
  SuspectRulePtr vnp, vnp_prev = NULL, new_rule;

  if (f == NULL) return FALSE;

  for (vnp = f->rule_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
    vnp_prev = vnp;
  }
  if (vnp != NULL) {
    new_rule = EditSuspectRule (vnp);
    if (new_rule != NULL) {
      if (vnp_prev == NULL) {
        f->rule_list = new_rule;
      } else {
        vnp_prev->next = new_rule;
      }
      new_rule->next = vnp->next;
      vnp->next = NULL;
      vnp = SuspectRuleFree (vnp);
      rval = TRUE;
    }
  }
  return rval;
}


static void ClickSuspectProductRuleDoc (DoC d, PoinT pt)
{
  Int2               item, row, col;
  RecT               rct;
  SuspectProductRuleEditorFormPtr f;
  SuspectRuleSetPtr  vnp, vnp_prev = NULL, two_prev = NULL, vnp_next;
  Boolean            changed = FALSE;
  BaR                sb_vert = NULL;
  Int2               scroll_pos = 0;
  
  f = (SuspectProductRuleEditorFormPtr) GetObjectExtra (d);
  if (f == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, &col, &rct);
  if (item == 0 && row == 0 && f->rule_list == NULL) {
    AddSuspectProductRule (f, 0);
  } else if (item > 0 && row > 0) {
    if (item == 1) {
      /* add to beginning of list */
      AddSuspectProductRule (f, 0);
    } else if (item == CountSuspectRuleSet (f->rule_list) + 2) {
      /* add to end of list */
      AddSuspectProductRule (f, item);
    } else {
      /* correct for explanatory line */
      item--;
      sb_vert = GetSlateVScrollBar ((SlatE) f->rule_summary);
      scroll_pos = GetBarValue (sb_vert);
      switch (col) {
        case 2:
          /* delete this item */
          for (vnp = f->rule_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
            vnp_prev = vnp;
          }
          if (vnp != NULL) {
            if (vnp_prev == NULL) {
              f->rule_list = vnp->next;
            } else {
              vnp_prev->next = vnp->next;
            }
            vnp->next = NULL;
            vnp = SuspectRuleFree (vnp);
            changed = TRUE;
          }
          break;
        case 3:
          /* move this item up */
          for (vnp = f->rule_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
            two_prev = vnp_prev;
            vnp_prev = vnp;
          }
          if (vnp != NULL && vnp_prev != NULL) {
            vnp_prev->next = vnp->next;
            vnp->next = vnp_prev;
            if (two_prev == NULL) {
              f->rule_list = vnp;
            } else {
              two_prev->next = vnp;
            }
            /* decrease the scroll position, so cursor will still be over the same item */
            scroll_pos--;
            changed = TRUE;
          }
          break;
        case 4:
          /* move this item down */
          for (vnp = f->rule_list; vnp != NULL && item > 1; vnp = vnp->next, item--) {
            vnp_prev = vnp;
          }
          if (vnp != NULL && vnp->next != NULL) {
            vnp_next = vnp->next;
            vnp->next = vnp_next->next;
            vnp_next->next = vnp;
            if (vnp_prev == NULL) {
              f->rule_list = vnp_next;
            } else {
              vnp_prev->next = vnp_next;
            }
            /* increase the scroll position, so cursor will still be over the same item */
            scroll_pos++;
            changed = TRUE;
          }
          break;
        case 5:
          /* insert item */
          if (pt.y >= rct.top && pt.y <= rct.top + 4) {
            /* insert macro before this one */
            AddSuspectProductRule (f, item - 1);
          } else if (pt.y >= rct.bottom - 4 && pt.y <= rct.bottom) {
            /* insert macro after this one */
            AddSuspectProductRule (f, item);
          }
          break;
        case 6:
          /* clone item */
          CloneSuspectProductRuleItem (f, item);
          changed = TRUE;
          break;
        case 7:
          /* edit this item */
          changed = EditSuspectProductRuleItem (f, item);
          break;
      }
    }
  }
  if (changed) {
    f->unsaved = TRUE;
    UpdateSuspectProductRuleSummary (f, scroll_pos);
  }
}


static void DrawSuspectProductRuleDocControls (DoC d, RectPtr r, Int2 item, Int2 firstLine)

{
  SuspectProductRuleEditorFormPtr dlg;
  RecT               rct;
  Int4               width;
  PoinT              pt1, pt2;
  Int4               num_items;

  dlg = (SuspectProductRuleEditorFormPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  num_items = CountSuspectRuleSet (dlg->rule_list);

  /* don't draw controls for explanatory text */
  if (item == 1 || item >= num_items + 2) return;

  if (dlg != NULL && r != NULL && item > 0 && firstLine == 0) {
    rct = *r;

    rct.left += kNumberColumnWidth;

    width = 10;
    /* draw X for deletion */
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1;
    pt2.x = pt1.x + width;
    pt2.y = pt1.y + width;
    DrawLine (pt1, pt2);
    pt1.x = rct.left + 1;
    pt1.y = rct.top + 1 + width;
    pt2.x = pt1.x + width;
    pt2.y = rct.top + 1;
    DrawLine (pt1, pt2);

    /* draw up arrow for moving step up */
    if (item > 2) {
      pt1.x = rct.left + width + 3;
      pt1.y = rct.top + 3;
      pt2.x = pt1.x + 5;
      pt2.y = rct.top + 1;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x + 5;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x;
      pt1.y = pt2.y + width;
      DrawLine (pt1, pt2);
    }
    /* draw up arrow for moving step up */
    if (item < num_items + 1) {
      pt1.x = rct.left + 2 * width + 5;
      pt1.y = rct.top + width - 2;
      pt2.x = pt1.x + 5;
      pt2.y = rct.top + width + 1;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x + 5;
      DrawLine (pt1, pt2);
      pt1.x = pt2.x;
      pt1.y = rct.top + 1;
      DrawLine (pt1, pt2);
    }

    /* draw insertion controls */
    pt1.x = rct.left + 3 * width + 7;
    pt1.y = rct.top + 4;
    pt2.x = pt1.x + width;
    pt2.y = rct.top;
    if (item > 2) {
      DrawLine (pt1, pt2);
    }
    pt1.y = rct.bottom - 4;
    pt2.y = rct.bottom;
    if (item < num_items + 1) {
      DrawLine (pt1, pt2);
    }
  }
}


static void ListSuspectRuleMatches (ButtoN b)
{
  SuspectProductRuleEditorFormPtr  f;
  ValNodePtr sep_list;
  Char         path [PATH_MAX];
  FILE         *fp;

  f = (SuspectProductRuleEditorFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }

  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    return;
  }

  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp != NULL) {
    PrintSuspectRuleMatches (sep_list->data.ptrvalue, f->rule_list, fp);
    FileClose (fp);
    LaunchGeneralTextViewer (path, "Suspect Rule Matches");
    FileRemove (path);
  }

  sep_list = ValNodeFree (sep_list);
}


static void ApplySuspectRuleFixes (ButtoN b)
{
  SuspectProductRuleEditorFormPtr  f;
  ValNodePtr sep_list;
  Char         path [PATH_MAX];
  FILE         *fp;
  SeqEntryPtr  sep;

  f = (SuspectProductRuleEditorFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }

  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    return;
  }
  sep = sep_list->data.ptrvalue;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp != NULL) {
    ApplySuspectRuleFixesToSeqEntry (sep, f->rule_list, fp);

    FileClose (fp);
    LaunchGeneralTextViewer (path, "Suspect Rule Fixes");
    FileRemove (path);
  }

  sep_list = ValNodeFree (sep_list);
}


static void MakeProductUpdateTable (ButtoN b)
{
  SuspectProductRuleEditorFormPtr  f;
  ValNodePtr sep_list;
  Char         path [PATH_MAX];
  FILE         *fp;
  SeqEntryPtr  sep;

  f = (SuspectProductRuleEditorFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }

  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    return;
  }
  sep = sep_list->data.ptrvalue;

  if (GetOutputFileName (path, sizeof (path), NULL)) {
    fp = FileOpen (path, "w");
    if (fp == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
    } else {
      ExportProductUpdateTableWithPrecomputedSuggestions (fp, sep, f->rule_list);
      FileClose (fp);
#ifdef WIN_MSWIN
      Nlm_MSWin_OpenApplication ("excel.exe", path);
#endif
    }
  }
  sep_list = ValNodeFree (sep_list);
}


static void ExportRuleSetText (SuspectRuleSetPtr set, FILE *fp)
{
  CharPtr tmp;

  if (fp == NULL) {
    return;
  }
  while (set != NULL) {
    tmp = SummarizeSuspectRule(set);
    if (tmp != NULL) {
      fprintf (fp, "%s\n", tmp);
    }
    set = set->next;
  }
}


static void DisplayRuleText (SuspectRuleSetPtr set, CharPtr title)
{
  Char  path [PATH_MAX];
  FILE *fp;

  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp != NULL) {
    ExportRuleSetText(set, fp);
    FileClose (fp);
    LaunchGeneralTextViewer (path, title);
  }
  FileRemove (path);
}


static void DisplaySuspectProductRuleDescriptions (ButtoN b)
{
  SuspectProductRuleEditorFormPtr frm;

  frm = (SuspectProductRuleEditorFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }
  DisplayRuleText (frm->rule_list, "Suspect Rule Descriptions");
}


static void ShowSuspectRuleDiffs (ButtoN b)
{
  SuspectProductRuleEditorFormPtr  f;
  Char               path [PATH_MAX];
  AsnIoPtr           aip;
  SuspectRuleSetPtr  rule_list;
  SuspectRuleSetPtr  in1not2 = NULL, in2not1 = NULL;
  CharPtr            title;
  CharPtr            title1_fmt = "Found in current list but not %s";
  CharPtr            title2_fmt = "Found in %s but not in current list";

  f = (SuspectProductRuleEditorFormPtr) GetObjectExtra (b);
  if (f == NULL) {
    return;
  }

  path [0] = '\0';
  if (GetInputFileName (path, sizeof (path), NULL, "TEXT")) {
    aip = AsnIoOpen (path, "r");
    if (aip == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
    } else {
      rule_list = SuspectRuleSetAsnRead (aip, NULL);
      if (rule_list == NULL) {
        Message (MSG_ERROR, "Unable to read rule list from %s.", path);
      } else {
        FindDiffsBetweenRuleSets (f->rule_list, rule_list, &in1not2, &in2not1);
        title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title1_fmt) + StringLen (path)));
        sprintf (title, title1_fmt, path);
        DisplayRuleText (in1not2, title);
        title = MemFree (title);
        title = (CharPtr) MemNew (sizeof (Char) * (StringLen (title2_fmt) + StringLen (path)));
        sprintf (title, title2_fmt, path);
        DisplayRuleText (in2not1, title);
        title = MemFree (title);
        rule_list = SuspectRuleSetFree (rule_list);
        in1not2 = SuspectRuleSetFree (in1not2);
        in2not1 = SuspectRuleSetFree (in2not1);
      }
    }
  }
}


static void ApplyProductUpdateTableBtn (ButtoN b)
{
  ValNodePtr   sep_list, vnp, table;
  Char         path [PATH_MAX];
  FILE         *fp;
  SeqEntryPtr  sep;
  Uint2        entityID;
  LogInfoPtr   lip;

  sep_list = GetViewedSeqEntryList ();
  if (sep_list == NULL) {
    return;
  }

  path [0] = '\0';
  if (GetInputFileName (path, sizeof (path), "", "TEXT")) {
    fp = FileOpen (path, "r");
    if (fp == NULL) {
      Message (MSG_ERROR, "Unable to open %s", path);
    } else {
      table = ReadProductUpdateTable (fp);
      FileClose (fp);
      if (table == NULL) {
        Message (MSG_ERROR, "Unable to read table from %s", path);
      } else {
        lip = OpenLog("Product name changes");
        for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
          sep = vnp->data.ptrvalue;
          if (sep == NULL) continue;
          entityID = ObjMgrGetEntityIDForChoice(sep);
          lip->data_in_log |= ApplyProductUpdateTable (table, sep, lip->fp);
          ObjMgrSetDirtyFlag (entityID, TRUE);
          ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
        }
      }
      CloseLog (lip);
      FreeLog (lip);
      table = ProductUpdateTableFree (table);
    }
  }
  sep_list = ValNodeFree (sep_list);
}


NLM_EXTERN void LaunchSuspectProductRuleEditorBaseForm (BaseFormPtr bfp)
{
  WindoW              w;
  SuspectProductRuleEditorFormPtr  f;
  GrouP               h, c;
  MenU                m;
  ButtoN              b;

  if (bfp == NULL) return;

  f = (SuspectProductRuleEditorFormPtr) MemNew (sizeof (SuspectProductRuleEditorFormData));
  if (f == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Suspect Product Rule Editor", StdCloseWindowProc);
  SetObjectExtra (w, f, CleanupSuspectProductRuleEditorForm);
  f->form = (ForM) w;
  f->input_entityID = bfp->input_entityID;

  f->formmessage = SuspectProductRuleEditorFormMessage;

  f->rule_list = NULL;
  f->last_filename = NULL;
  f->unsaved = FALSE;

  m = PulldownMenu (w, "File");
  FormCommandItem (m, "Open", (BaseFormPtr)f, VIB_MSG_OPEN);
  FormCommandItem (m, "Add Rules from File to List", (BaseFormPtr)f, VIB_MSG_IMPORT);
  FormCommandItem (m, "Save", (BaseFormPtr)f, VIB_MSG_SAVE);
  FormCommandItem (m, "Save As", (BaseFormPtr)f, VIB_MSG_SAVE_AS);
  SeparatorItem (m);
  FormCommandItem (m, "Quit", (BaseFormPtr)f, VIB_MSG_QUIT);
  m = PulldownMenu (w, "Sort");
  FormCommandItem (m, "By Find", (BaseFormPtr) f, NUM_VIB_MSG + 2);
  FormCommandItem (m, "By Category, then Find", (BaseFormPtr) f, NUM_VIB_MSG + 3);


  SetupSuspectProductRuleEditorFont (f);
  SetSuspectProductRuleEditorFileItems (f);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  f->rule_summary = DocumentPanel (h, stdCharWidth * 50, stdLineHeight * 20);
  SetObjectExtra (f->rule_summary, f, NULL);
  SetDocProcs (f->rule_summary, ClickSuspectProductRuleDoc, NULL, NULL, NULL);
  SetDocShade (f->rule_summary, DrawSuspectProductRuleDocControls, NULL, NULL, NULL); 

  c = HiddenGroup (h, 6, 0, NULL);
  b = PushButton (c, "List current matches", ListSuspectRuleMatches);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Apply fixes", ApplySuspectRuleFixes);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Make Product Update Table", MakeProductUpdateTable);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Apply Product Update Table", ApplyProductUpdateTableBtn);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Show diffs with other file", ShowSuspectRuleDiffs);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Display Rule Text", DisplaySuspectProductRuleDescriptions);
  SetObjectExtra (b, f, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) f->rule_summary, (HANDLE) c, NULL);
  SummarizeSuspectProductRule (f->rule_summary, f->rule_list, f->summary_font);
  Show (w);
}


NLM_EXTERN void LaunchSuspectProductRuleEditor (IteM i)
{
  BaseFormPtr         bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  LaunchSuspectProductRuleEditorBaseForm (bfp);
}
