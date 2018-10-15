/*   sequin6.c
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
* File Name:  sequin6.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   11/12/97
*
* $Revision: 6.383 $
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
#include <ncbilang.h>
#include <gather.h>
#include <asn2gnbp.h>
#include <bspview.h>
#include <import.h>
#include <objsub.h>
#include <explore.h>
#include <subutil.h>
#include <gbftdef.h>
#include <edutil.h>
#include <salpanel.h>
#include <seqpanel.h>
#include <biosrc.h>
#include <vsm.h>
#include <actutils.h>
#include <findrepl.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macrodlg.h>
#include <macroapi.h>

#define NUMBER_OF_SUFFIXES    7

static ENUM_ALIST(name_suffix_alist)
  {" ",    0},
  {"Jr.",  1},
  {"Sr.",  2},
  {"II",   4},
  {"III",  5},
  {"IV",   6},
  {"V",    7},
  {"VI",   8},
END_ENUM_ALIST


#define PUBLICATION_PUBLISHED_FIELD 1
#define PUBLICATION_INPRESS_FIELD   2
#define PUBLICATION_UNPUB_FIELD     3

static ENUM_ALIST (publication_field_alist)
  {" ",                0},
  {"Published",        PUBLICATION_PUBLISHED_FIELD},
  {"In Press",         PUBLICATION_INPRESS_FIELD},
  {"Unpublished",      PUBLICATION_UNPUB_FIELD},
END_ENUM_ALIST


#define GENE_LOCUS_FIELD        1
#define GENE_DESCRIPTION_FIELD  2
#define GENE_ALLELE_FIELD       3
#define GENE_MAPLOC_FIELD       4
#define GENE_LOCUS_TAG_FIELD    5
#define GENE_SYNONYM_FIELD      6
#define GENE_COMMENT_FIELD      7

static ENUM_ALIST(gene_field_alist)
  {" ",                    0},
  {"Gene locus",           GENE_LOCUS_FIELD},
  {"Gene description",     GENE_DESCRIPTION_FIELD},
  {"Gene allele",          GENE_ALLELE_FIELD},
  {"Gene maploc",          GENE_MAPLOC_FIELD},
  {"Locus_tag",            GENE_LOCUS_TAG_FIELD},
  {"Gene synonym",         GENE_SYNONYM_FIELD},
  {"Gene comment",         GENE_COMMENT_FIELD},
END_ENUM_ALIST

#define CDS_COMMENT   1
#define CDS_GENE_XREF 2
#define CDS_DB_XREF   3

static ENUM_ALIST(cds_field_alist)
  {" ",                    0},
  {"CDS comment",          CDS_COMMENT},
  {"gene xref",            CDS_GENE_XREF},
  {"db_xref",              CDS_DB_XREF},
END_ENUM_ALIST

static ENUM_ALIST(cds_short_field_alist)
  {" ",                    0},
  {"CDS comment",          CDS_COMMENT},
END_ENUM_ALIST

#define PROT_NAME_FIELD        1
#define PROT_DESCRIPTION_FIELD 2
#define PROT_ECNUM_FIELD       3
#define PROT_ACTIVITY_FIELD    4
#define PROT_COMMENT_FIELD     5

static ENUM_ALIST(prot_field_alist)
  {" ",                    0},
  {"Protein name",         PROT_NAME_FIELD},
  {"Protein description",  PROT_DESCRIPTION_FIELD},
  {"Protein E.C. number",  PROT_ECNUM_FIELD},
  {"Protein activity",     PROT_ACTIVITY_FIELD},
  {"Protein comment",      PROT_COMMENT_FIELD},
END_ENUM_ALIST

#define RNA_NAME_FIELD       1
#define RNA_COMMENT_FIELD    2
#define RNA_GENE_XREF_FIELD  3

static ENUM_ALIST(rnax_field_alist)
  {" ",                    0},
  {"RNA Name",             RNA_NAME_FIELD},
  {"RNA Comment",          RNA_COMMENT_FIELD},
  {"gene xref",            RNA_GENE_XREF_FIELD},
END_ENUM_ALIST

static ENUM_ALIST(rnax_short_field_alist)
  {" ",                    0},
  {"RNA Name",             RNA_NAME_FIELD},
  {"RNA Comment",          RNA_COMMENT_FIELD},
  {"gene xref",            RNA_GENE_XREF_FIELD},
END_ENUM_ALIST

#define ORGREF_SCI_NAME_FIELD    1
#define ORGREF_COMMON_NAME_FIELD 2
#define ORGREF_LINEAGE_FIELD     3
#define ORGREF_DIVISION_FIELD    4

static ENUM_ALIST(orgref_field_alist)
  {" ",                    0},
  {"Scientific Name",      ORGREF_SCI_NAME_FIELD},
  {"Common Name",          ORGREF_COMMON_NAME_FIELD},
  {"Lineage",              ORGREF_LINEAGE_FIELD},
  {"Division",             ORGREF_DIVISION_FIELD},
END_ENUM_ALIST

#define IMPORT_GBQUAL_FIELD  1
#define IMPORT_COMMENT_FIELD 2

static ENUM_ALIST(impfeat_field_alist)
  {" ",        0},
  {"GBQual",   IMPORT_GBQUAL_FIELD},
  {"Comment",  IMPORT_COMMENT_FIELD},
END_ENUM_ALIST


#define NUM_SUBTARGET_POPUPS 12

static GbFeatName ParseQualifierList[] = {
 {"allele", Class_text}, {"anticodon", Class_pos_aa},
 {"bound_moiety", Class_text},
 {"chromosome", Class_text},
 {"citation", Class_bracket_int},
 {"codon", Class_seq_aa},
 {"codon_start", Class_int_or}, {"cons_splice", Class_site},
 {"db_xref", Class_text},
 {"direction", Class_L_R_B}, {"EC_number", Class_ecnum},
 {"evidence", Class_exper}, {"exception", Class_text},
 {"frequency", Class_text}, {"function", Class_text},
 {"gene", Class_text}, {"gdb_xref", Class_text},
 {"insertion_seq", Class_text},
 {"label", Class_token},
 {"map", Class_text},
 {"mod_base", Class_token}, {"note", Class_note},
 {"number", Class_number}, {"organism", Class_text},
 {"partial", Class_none}, {"PCR_conditions", Class_text},
 {"phenotype", Class_text},
 {"plasmid", Class_text}, {"product", Class_text},
 {"pseudo", Class_none},
 {"rearranged", Class_none}, { "replace", Class_text},
 {"rpt_family", Class_text}, {"rpt_type", Class_rpt},
 {"rpt_unit", Class_token},
 {"sequenced_mol", Class_text},
 {"standard_name", Class_text},
 {"translation", Class_text}, {"transl_except", Class_pos_aa},
 {"transl_table", Class_int}, {"transposon", Class_text},
 {"usedin", Class_token},
 {"focus", Class_none},
 {"protein_id", Class_text},
 {"organelle", Class_text}, {"transcript_id", Class_text},
 {"transgenic", Class_none}, {"environmental_sample", Class_none},
 {"locus_tag", Class_text}, {"mol_type", Class_text},
};

const Int4 NumParseQualifiers = sizeof (ParseQualifierList) / sizeof (GbFeatName);


static void PopulateParseQualPopup (PopuP p)
{
  Int4 i;
  for (i = 0; i < NumParseQualifiers; i++) {
    PopupItem (p, ParseQualifierList[i].name);
  }
}


static DialoG ImpFeatSelectDialog (GrouP g, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  EnumFieldAssocPtr list, ap;
  ValNodePtr   choice_list = NULL;

  list = import_featdef_alist (FALSE, FALSE, FALSE);
  SortEnumFieldAssocPtrArray (list, CompareImpFeatEnumFieldAssoc);
  for (ap = list; ap->name != NULL; ap++) {
    if (ap == list) {
      ValNodeAddPointer (&choice_list, 0, StringSave ("All import features"));
    } else {
      ValNodeAddPointer (&choice_list, 0, StringSave (ap->name));
    }
  }
  return ValNodeSelectionDialog (g, choice_list, TALL_SELECTION_LIST,
                                                ValNodeStringName,
                                ValNodeSimpleDataFree, ValNodeStringCopy,
                                ValNodeChoiceMatch, "Import Feature Type",
                                change_notify, change_userdata, FALSE);
}

typedef enum {
  eFieldTypeGene = 0,
  eFieldTypeCDS,
  eFieldTypeProtein,
  eFieldTypeRNA,
  eFieldTypeBioSource,
  eFieldTypeOrgModSubSource,
  eFieldTypeImport,
  eFieldTypeDefline,
  eFieldTypeCommentDescriptor,
  eFieldTypeFeatureNote,
  eFieldTypePublication,
  eFieldType_Max
} SegregateFieldType;

static CharPtr field_type_names[] = {
  "Gene",
  "CDS",
  "Protein",
  "RNA",
  "BioSource",
  "OrgMod and SubSource",
  "Import Feature",
  "DefLine",
  "Comment Descriptor",
  "Feature Note",
  "Publication" };


typedef struct fieldsubfield {
  Int4 field;
  Int4 subfield;
  CharPtr impfeat_key;
  ValNodePtr subfield_list;
} FieldSubfieldData, PNTR FieldSubfieldPtr;


static FieldSubfieldPtr FieldSubfieldFree (FieldSubfieldPtr f)
{
  if (f != NULL) {
    f->subfield_list = ValNodeFree (f->subfield_list);
    f->impfeat_key = MemFree (f->impfeat_key);
    f = MemFree (f);
  }
  return f;
}


typedef struct fieldsubfielddlg {
  DIALOG_MESSAGE_BLOCK
  PopuP field_list;
  DialoG subfield_dlg[eFieldType_Max];
  DialoG impfeat_type;
  Boolean allowed_fields[eFieldType_Max];
  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} FieldSubfieldDlgData, PNTR FieldSubfieldDlgPtr;


static Int4 FieldNumFromListVal (Int4 list_val, BoolPtr allowed_fields)
{
  Int4 i;
  Int4 field_num = -1;

  if (list_val < 1) return -1;

  for (i = 0; i < eFieldType_Max && field_num < 0; i++) {
    if (allowed_fields[i]) {
      if (list_val - 1 == i) {
        field_num = i;
      }
    } else {
      list_val++;
    }
  }
  return field_num;
}


static Int4 ListValFromFieldNum (Int4 field_num, BoolPtr allowed_fields)
{
  Int4 list_val = 0, i = 0;

  if (field_num < 0) return 0;

  while (i <= field_num) {
    if (allowed_fields[i]) {
      list_val++;
    }
    i++;
  }
  return list_val;
}


static void ChangeFieldType (PopuP p)
{
  FieldSubfieldDlgPtr dlg;
  Int4                list_val, field_num, i;

  dlg = (FieldSubfieldDlgPtr) GetObjectExtra (p);
  if (dlg == NULL) return;

  list_val = GetValue (dlg->field_list);
  field_num = FieldNumFromListVal (list_val, dlg->allowed_fields);

  for (i = 0; i < eFieldType_Max; i++) {
    if (dlg->subfield_dlg[i] == NULL) continue;
    if (i == field_num) {
      Show (dlg->subfield_dlg[i]);
    } else {
      Hide (dlg->subfield_dlg[i]);
    }
  }
  if (field_num == eFieldTypeImport) {
    Show (dlg->impfeat_type);
  } else {
    Hide (dlg->impfeat_type);
  }
  if (dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void FieldSubfieldToDialog (DialoG d, Pointer data)
{
  FieldSubfieldDlgPtr dlg;
  FieldSubfieldPtr    f;
  Int4                list_val;
  ValNode             vn;

  dlg = (FieldSubfieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  f = (FieldSubfieldPtr) data;
  if (f == NULL || f->field < 0) {
    SetValue (dlg->field_list, 0);
  } else {
    list_val = ListValFromFieldNum (f->field, dlg->allowed_fields);
    SetValue (dlg->field_list, list_val);
    if (dlg->subfield_dlg[f->field] != NULL) {    
      if (list_val == eFieldTypeFeatureNote) {
        PointerToDialog (dlg->subfield_dlg[eFieldTypeFeatureNote], f->subfield_list);
      } else {
        vn.choice = f->subfield;
        vn.data.ptrvalue = NULL;
        vn.next = NULL;
        PointerToDialog (dlg->subfield_dlg[f->field], &vn);
        if (f->field == eFieldTypeImport) {
          vn.choice = 0;
          vn.data.ptrvalue = f->impfeat_key;
          PointerToDialog (dlg->impfeat_type, &vn);
        }
      }
    }
  }
  ChangeFieldType (dlg->field_list);
}


static Pointer DialogToFieldSubfield (DialoG d)
{
  FieldSubfieldDlgPtr dlg;
  FieldSubfieldPtr f;
  Int4             val;
  ValNodePtr       vnp;
  SourceQualDescPtr sqdp;

  dlg = (FieldSubfieldDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  f = (FieldSubfieldPtr) MemNew (sizeof (FieldSubfieldData));
  f->field = -1;
  f->subfield = -1;
  f->subfield_list = NULL;
  f->impfeat_key = NULL;

  val = GetValue (dlg->field_list);
  if (val > 0) {
    f->field = FieldNumFromListVal (val, dlg->allowed_fields);    
    if (dlg->subfield_dlg[f->field] == NULL) {
      f->subfield = 0;
    } else {
      if (f->field == eFieldTypeFeatureNote) {
        f->subfield_list = DialogToPointer (dlg->subfield_dlg[f->field]);
      } else {
        vnp = DialogToPointer (dlg->subfield_dlg[f->field]);
        if (vnp != NULL) {
          if (f->field == eFieldTypeOrgModSubSource) {
            sqdp = (SourceQualDescPtr) vnp->data.ptrvalue;
            if (sqdp != NULL) {
              f->subfield = sqdp->subtype;
              if (!sqdp->isOrgMod) {
                f->subfield += 1000;
              }
            }
          } else {
            f->subfield = vnp->choice;
            vnp = ValNodeFree (vnp);
          }
        }
        if (f->field == eFieldTypeImport) {
          vnp = DialogToPointer (dlg->impfeat_type);
          if (vnp != NULL) {    
            f->impfeat_key = StringSave (vnp->data.ptrvalue);
            vnp = ValNodeFree (vnp);
          }
        }
      }
    }
  }

  return f;
}


static ValNodePtr TestFieldSubfieldDialog (DialoG d)
{
  ValNodePtr          err_list = NULL;
  FieldSubfieldPtr    f;
  
  f = DialogToPointer (d);
  if (f == NULL || f->field < 0 
      || (f->subfield < 0 && f->subfield_list == NULL)
      || (f->field == eFieldTypeImport && f->impfeat_key == NULL)) {
    ValNodeAddPointer (&err_list, 0, "Must choose field");
  }
  f = FieldSubfieldFree (f);
  return err_list;  
}


static DialoG 
CreateFieldSubfieldDlg 
(GrouP h,
 BoolPtr allowed_fields,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata)
{
  FieldSubfieldDlgPtr dlg;
  GrouP               p, g, g2;
  Int4                i;
  ValNode             vn;
  
  dlg = (FieldSubfieldDlgPtr) MemNew (sizeof (FieldSubfieldDlgData));
  if (dlg == NULL)
  {
    return NULL;
  }
  
  p = HiddenGroup (h, 2, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);

  dlg->dialog = (DialoG) p;
  dlg->todialog = FieldSubfieldToDialog;
  dlg->fromdialog = DialogToFieldSubfield;
  dlg->testdialog = TestFieldSubfieldDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;

  if (allowed_fields == NULL) {
    MemSet (dlg->allowed_fields, TRUE, sizeof (dlg->allowed_fields));
  } else {
    for (i = 0; i < eFieldType_Max; i++) {      
      dlg->allowed_fields[i] = allowed_fields[i];
    }
  }

  dlg->field_list = PopupList (p, TRUE, ChangeFieldType);
  SetObjectExtra (dlg->field_list, dlg, NULL);
  for (i = 0; i < eFieldType_Max; i++) {
    if (dlg->allowed_fields[i]) {
      PopupItem (dlg->field_list, field_type_names[i]);
    }
  }

  g = HiddenGroup (p, 0, 0, NULL);
  dlg->subfield_dlg[eFieldTypeGene] = EnumAssocSelectionDialog (g, gene_field_alist, 
                                                                NULL, FALSE, 
                                                                dlg->change_notify, 
                                                                dlg->change_userdata);
  dlg->subfield_dlg[eFieldTypeCDS] = EnumAssocSelectionDialog (g, cds_field_alist,
                                                                NULL, FALSE, 
                                                                dlg->change_notify, 
                                                                dlg->change_userdata);
  dlg->subfield_dlg[eFieldTypeProtein] = EnumAssocSelectionDialog (g, prot_field_alist,
                                                                NULL, FALSE, 
                                                                dlg->change_notify, 
                                                                dlg->change_userdata);
  dlg->subfield_dlg[eFieldTypeRNA] = EnumAssocSelectionDialog (g, rnax_field_alist,
                                                                NULL, FALSE, 
                                                                dlg->change_notify, 
                                                                dlg->change_userdata);
  dlg->subfield_dlg[eFieldTypeBioSource] = EnumAssocSelectionDialog (g, orgref_field_alist,
                                                                NULL, FALSE, 
                                                                dlg->change_notify, 
                                                                dlg->change_userdata);
  /* Set default BioSource field to scientific name */
  vn.choice = ORGREF_SCI_NAME_FIELD;
  vn.next = NULL;
  vn.data.ptrvalue = NULL;
  PointerToDialog (dlg->subfield_dlg[eFieldTypeBioSource], &vn);

  dlg->subfield_dlg[eFieldTypeOrgModSubSource] = SourceQualTypeSelectionDialog (g, FALSE,
                                                                                dlg->change_notify,
                                                                                dlg->change_userdata);
                                                                                  
  if (dlg->allowed_fields[eFieldTypeImport]) {
    g2 = HiddenGroup (g, 2, 0, NULL);
    dlg->impfeat_type = ImpFeatSelectDialog (g2, change_notify, change_userdata);
    Hide (dlg->impfeat_type);

    dlg->subfield_dlg[eFieldTypeImport] = EnumAssocSelectionDialog (g2, impfeat_field_alist,
                                                                NULL, FALSE, 
                                                                dlg->change_notify, 
                                                                dlg->change_userdata);
  }
  dlg->subfield_dlg[eFieldTypeFeatureNote] = FeatureSelectionDialog (g, TRUE, dlg->change_notify, dlg->change_userdata);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->subfield_dlg[eFieldTypeGene],
                              (HANDLE) dlg->subfield_dlg[eFieldTypeCDS],
                              (HANDLE) dlg->subfield_dlg[eFieldTypeProtein],
                              (HANDLE) dlg->subfield_dlg[eFieldTypeRNA],
                              (HANDLE) dlg->subfield_dlg[eFieldTypeBioSource],
                              (HANDLE) dlg->subfield_dlg[eFieldTypeOrgModSubSource],
                              (HANDLE) dlg->subfield_dlg[eFieldTypeFeatureNote],
                              (HANDLE) g2,
                              NULL);

  for (i = 0; i < eFieldType_Max; i++) {
    Hide (dlg->subfield_dlg[i]);
  }
  return (DialoG) p;
}


typedef struct convertformdata {
  FEATURE_FORM_BLOCK

  TexT               deleteText;
  CharPtr            deleteStr;
  GrouP              deleteLevel;
  Int2               deleteLevelInt;
  DialoG             target_dlg;
  Int2               impfeat_type;
  ButtoN             accept;
  ButtoN             leaveDlgUp;
  Int2               type;
  Int2               subtype;

  DialoG             textportion_dlg;
  TextPortionXPtr     textportion;
  Boolean            remove_inside;

  CharPtr            foundstr;
  GrouP              repeat_remove_grp;

  Boolean            isDirty;
  Boolean            repeat_remove;
  GrouP              ifNotFoundGroup;
  Int2               ifNotFound;
  BioseqSetPtr       target_set;
  Nlm_ChangeNotifyProc set_accept_proc;
  ValNodePtr         id_list;

} ConvertFormData, PNTR ConvertFormPtr;

static ConvertFormPtr ConvertFormNew (void)
{
  ConvertFormPtr cfp;

  cfp = (ConvertFormPtr) MemNew (sizeof (ConvertFormData));
  if (cfp == NULL) return NULL;
  MemSet (cfp, 0, sizeof (ConvertFormData));
  return cfp;
}


/* Values for ifNotFound field */

#define DO_NOTHING       2
#define REMOVE_ALL_TEXT  3

/* End values for string trimming */

#define TRIM_LEFT_END    1
#define TRIM_RIGHT_END   2

/*-------------------------------------------------------------------------*/
/*                                                                         */
/* TrimOffEndQuotes () -- Trim double-quotes off the ends of a string.     */
/*                                                                         */
/*-------------------------------------------------------------------------*/

static Boolean TrimOffEndQuotes (CharPtr trimString, 
				 Int2    whichEnd)
{
  Int4 strLen;
  Int4 i;

  strLen = StringLen(trimString);
  if (strLen == 0)
    return FALSE;

  /* If there is a quote at the end, remove it */

  if (TRIM_RIGHT_END == whichEnd) {
    if (trimString[strLen-1] == '"') {
      strLen--;
      trimString [strLen] = '\0';
      return TRUE;
    }
  }

  /* If there is a quote at the beginning, remove it */

  else {
    if (trimString[0] == '"') {
      for (i = 0; trimString[i] != '\0'; i++)
	trimString[i] = trimString [i+1];
      return TRUE;
    }
  }

  return FALSE;
}

/*-------------------------------------------------------------------------*/
/*                                                                         */
/* SaveStringFromTextNoStripSpaces () -                                    */
/*                                                                         */
/*-------------------------------------------------------------------------*/

static CharPtr SaveStringFromTextNoStripSpaces (TexT t)

{
  size_t   len;
  CharPtr  str;

  len = TextLength (t);
  if (len > 0) {
    str = (CharPtr) MemNew(len + 1);
    if (str != NULL) {
      GetTitle (t, str, len + 1);
      return str;
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}

static void ConvertFormTextCallback (TexT t)
{
  ConvertFormPtr cfp;

  cfp = (ConvertFormPtr) GetObjectExtra (t);
  if (cfp != NULL && cfp->set_accept_proc != NULL) {
    (cfp->set_accept_proc) (cfp);
  }
}


/*-------------------------------------------------------------------------*/
/*                                                                         */
/* SetDeleteAcceptButton () -- Enable/Disable the Accept button depending  */
/*                             on the condition of other window objects.   */
/*                                                                         */
/*-------------------------------------------------------------------------*/

static void SetDeleteAcceptButton (Pointer data)

{
  ConvertFormPtr  cfp;
  ValNodePtr      vnp;

  cfp = (ConvertFormPtr) data;
  if (cfp == NULL)
    return;

  /* Disable if a target or subtarget has not been selected */
  vnp = TestDialog (cfp->target_dlg);
  if (vnp != NULL) {
    vnp = ValNodeFree (vnp);
    SafeDisable (cfp->accept);
  } else if (TextHasNoText (cfp->deleteText)) {
    /* Disable if there is no delete string */
    SafeDisable (cfp->accept);
  } else {
    SafeEnable (cfp->accept);
  }
}


static void ChangeTargetFields (Pointer data)

{
  ConvertFormPtr  cfp;

  cfp = (ConvertFormPtr) data;
  if (cfp == NULL) return;

  if (cfp->set_accept_proc != NULL) {
    cfp->set_accept_proc (cfp);
  }
}


/*=========================================================================*/
/*                                                                         */
/* CheckForString () -- Searches for a given string with another string.   */
/*                                                                         */
/*=========================================================================*/

static Boolean CheckForString (CharPtr searchStr,
			       CharPtr sourceStr)
{
  if (NULL == SearchForString (sourceStr, searchStr, FALSE, FALSE))
    return FALSE;
  else
    return TRUE;
}

/*=========================================================================*/
/*                                                                         */
/* GeneRefPtr () -- Search a GeneRef feature for a given string.           */
/*                                                                         */
/*=========================================================================*/

static Boolean SearchGeneRef (GeneRefPtr     grp,
			      SeqFeatPtr     sfp,
			      ConvertFormPtr cfp)
{
  ValNodePtr vnp;

  /* Check parameters */

  if ((NULL == grp) || (NULL == sfp))
    return FALSE;

  /* Check each text field for the given string */

  switch (cfp->subtype) {
  case 1 :
    return CheckForString (cfp->deleteStr, grp->locus);
    break;
  case 2 :
    return CheckForString (cfp->deleteStr, grp->desc);
    break;
  case 3 :
    return CheckForString (cfp->deleteStr, grp->allele);
    break;
  case 4 :
    return CheckForString (cfp->deleteStr, grp->maploc);
    break;
  case 5 :
    return CheckForString (cfp->deleteStr, grp->locus_tag);
    break;
  case 6 :
    for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
      if (TRUE == CheckForString (cfp->deleteStr,
				  vnp->data.ptrvalue))
	return TRUE;
    }
    return FALSE;
    break;
  case 7 :
    return CheckForString (cfp->deleteStr, sfp->comment);
    break;
  default :
    break;
  }

  return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* SearchCDSFeat () -- Search a CDS feature for a given string.            */
/*                                                                         */
/*=========================================================================*/

static Boolean SearchCDSFeat (SeqFeatPtr     sfp,
			      ConvertFormPtr cfp)
{

  /* Check parameters */

  if (NULL == sfp)
    return FALSE;

  /* Check each text field for the given string */

  switch (cfp->subtype) {
  case 1 :
    return CheckForString (cfp->deleteStr, sfp->comment);
    break;
  case 2 :
  default :
    break;
  }

  /* If no match found, return FALSE */

  return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* SearchRnaRef () -- Search an RnaRef feature for a given string.         */
/*                                                                         */
/*=========================================================================*/

static Boolean SearchRnaRef (RnaRefPtr      rrp,
			     SeqFeatPtr     sfp,
			     ConvertFormPtr cfp)
{
  /* Check parameters */

  if ((NULL == rrp) || (NULL == sfp))
    return FALSE;

  /* Check each text field for the given string */

  switch (cfp->subtype) {
  case 1 :
    if ((0 == rrp->ext.choice) || (1 == rrp->ext.choice)) {
      return CheckForString (cfp->deleteStr, rrp->ext.value.ptrvalue);
    }
    break;
  case 2 :
    return CheckForString (cfp->deleteStr, sfp->comment);
    break;
  case 3 :
  default :
    break;
  }

  /* If no match found, return FALSE */

  return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* SearchProtRef () -- Search a ProtRef feature for a given string.        */
/*                                                                         */
/*=========================================================================*/

static Boolean SearchProtRef (ProtRefPtr     prp,
			      SeqFeatPtr     sfp,
			      ConvertFormPtr cfp)
{
  ValNodePtr vnp;

  /* Check parameters */

  if (NULL == prp)
    return FALSE;

  /* Check each text field for the given string */

  switch (cfp->subtype) {

    /* Search the name field */  

  case 1 :
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      if (TRUE == CheckForString (cfp->deleteStr,
				  (CharPtr) vnp->data.ptrvalue))
	return TRUE;
    }
    break;

    /* Search the desc field */

  case 2 :
    return CheckForString (cfp->deleteStr, prp->desc);
    break;

    /* Search the ec field */

  case 3 :
    for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
      if (TRUE == CheckForString (cfp->deleteStr, vnp->data.ptrvalue))
	return TRUE;
    }
    break;

    /* Search the activity field */

  case 4 :
    for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
      if (TRUE == CheckForString (cfp->deleteStr, vnp->data.ptrvalue))
	return TRUE;
    }
    break;

    /* Search the comment field */

  case 5 :
    return CheckForString (cfp->deleteStr, sfp->comment);
    break;

    /* Default check */

  default :
    break;
  }

  /* If we made it this far no match was found */

  return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* SearchImpFeat () -- Search an ImpFeat feature for a given string.       */
/*                                                                         */
/*=========================================================================*/

static Boolean SearchImpFeat (SeqFeatPtr     sfp,
			      ConvertFormPtr cfp)
{
  GBQualPtr       gbqp;

  /* Check parameters */

  if (NULL == sfp)
    return FALSE;

  switch (cfp->subtype) {

    /* Search the GB Quals */

  case IMPORT_GBQUAL_FIELD :
    gbqp = sfp->qual;
    while (NULL != gbqp) {
      if (NULL != gbqp->val)
	if (TRUE == CheckForString (cfp->deleteStr, gbqp->val))
	  return TRUE;
      gbqp = gbqp->next;
    }
    return FALSE;
    break;

    /* Search the comment field */

  case IMPORT_COMMENT_FIELD :
    return CheckForString (cfp->deleteStr, sfp->comment);
    break;
  default :
    break;
  }

  return FALSE;
}

/*=========================================================================*/
/*                                                                         */
/* MarkObjectsByText_Callback () - Called for each object in a SeqEntry    */
/*                                 this will mark the item for deletion    */
/*                                 if it matches the given criteria.       */
/*                                                                         */
/*=========================================================================*/

static void DeleteFeaturesByText_Callback (SeqEntryPtr sep,
					   Pointer     mydata,
					   Int4        index,
					   Int2        indent)
{
  ConvertFormPtr cfp;
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp = NULL;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  GeneRefPtr     grp;
  ProtRefPtr     prp;
  RnaRefPtr      rrp;
  Boolean        found = FALSE;

  /* Check parameters */

  if ((NULL == sep) || (NULL == sep->data.ptrvalue))
    return;

  cfp = (ConvertFormPtr) mydata;
  if (NULL == cfp)
    return;

  /* Get the list of annotations */

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  }
  else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  }
  else
    return;

  /* Search the requested item for the given string */

  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_GENE && cfp->type == eFieldTypeGene) {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
	  found = SearchGeneRef (grp, sfp, cfp);
        }
	else if (sfp->data.choice == SEQFEAT_CDREGION &&
		 cfp->type == eFieldTypeCDS) {
	  found = SearchCDSFeat (sfp, cfp);
        }
	else if (sfp->data.choice == SEQFEAT_PROT &&
		 cfp->type == eFieldTypeProtein) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
	  found = SearchProtRef (prp, sfp, cfp);
        }
	else if (sfp->data.choice == SEQFEAT_RNA && cfp->type == eFieldTypeRNA) {
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
	  found = SearchRnaRef (rrp, sfp, cfp);
        }
	else if (sfp->data.choice == SEQFEAT_IMP &&
		 cfp->type == eFieldTypeImport) {
	  found = SearchImpFeat (sfp, cfp);
        }
	if (TRUE == found)
	  break;
	else
	  sfp = sfp->next;
      }
    }
    if (TRUE == found)
      break;
    else
      sap = sap->next;
  }

  /* If we found the string, do the deletion */

  if (TRUE == found) {
    cfp->isDirty = TRUE;
    switch (cfp->deleteLevelInt) {
    case 1 :
      sfp->idx.deleteme = TRUE;
      break;
    case 2 :
      if (bsp != NULL)
	bsp->idx.deleteme = TRUE;
      break;
    case 3 :
      if (bssp != NULL)
	bssp->idx.deleteme = TRUE;
      else if (bsp != NULL) {
	if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
	  bssp = (BioseqSetPtr) bsp->idx.parentptr;
	  bssp->idx.deleteme = TRUE;
	}
      }
      break;
    default:
      break;
    }
  }

  /* Return successfully */

  return;
}

static Boolean DoesBioSourceContainText 
(BioSourcePtr biop,
 ConvertFormPtr cfp)
{
  OrgRefPtr     orp;
  Boolean       found;
  OrgNamePtr    onp;

  if (biop == NULL || cfp == NULL) return FALSE;

  found = FALSE;

  orp = biop->org;
  if (orp == NULL) return FALSE;
  switch (cfp->subtype) {
    case ORGREF_SCI_NAME_FIELD :
      found = CheckForString (cfp->deleteStr, orp->taxname);
      break;
    case ORGREF_COMMON_NAME_FIELD :
	  found = CheckForString (cfp->deleteStr, orp->common);
      break;
    case ORGREF_LINEAGE_FIELD :
      onp = orp->orgname;
      if (onp == NULL) {
        onp = OrgNameNew ();
        orp->orgname = onp;
      }
      if (onp != NULL)
	    found = CheckForString (cfp->deleteStr, onp->lineage);
	  else
	    found = FALSE;
      break;
    case ORGREF_DIVISION_FIELD :
      onp = orp->orgname;
      if (onp == NULL) {
        onp = OrgNameNew ();
        orp->orgname = onp;
      }
      if (onp != NULL)
        found = CheckForString (cfp->deleteStr, onp->div);
	  else
	    found = FALSE;
      break;
    default :
      break;
  }
  return found;
}

static Boolean DoSubSourcesContainText 
(BioSourcePtr biop,
 ConvertFormPtr cfp)
{
  OrgRefPtr     orp;
  OrgModPtr     mod;
  SubSourcePtr  ssp;
  OrgNamePtr    onp;
  Boolean found = FALSE;

  if (biop == NULL || cfp == NULL) return FALSE;
  if (cfp->subtype < 1000) {
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
	    mod = onp->mod;
	    while (mod != NULL && mod->subtype != cfp->subtype) {
	      mod = mod->next;
	    }
	    if (mod != NULL)
	      found = CheckForString (cfp->deleteStr, mod->subname);
	    else
	      found = FALSE;
	  }
	  else
	    found = FALSE;
	}
  } else {
	ssp = biop->subtype;
	while (ssp != NULL && ssp->subtype != (cfp->subtype - 1000)) {
	  ssp = ssp->next;
	}
	while (ssp != NULL) {
	  if (ssp->subtype == (cfp->subtype - 1000)) {
	    found = CheckForString (cfp->deleteStr, ssp->name);
	    if (found)
	      break;
	  }
	  ssp = ssp->next;
	}
	if (NULL == ssp)
	  found = FALSE;
  }
  return found;
}

/*=========================================================================*/
/*                                                                         */
/* DeleteSourceByText () - Given a text string, delete all items of a      */
/*                         given type that contain that string.            */
/*                                                                         */
/*=========================================================================*/

static void DeleteSourceByText (SeqDescrPtr    sdp,
				SeqEntryPtr    sep,
				BioseqPtr      bsp,
				ConvertFormPtr cfp)

{
  Boolean       found = FALSE;
  BioseqSetPtr  bssp;
  BioSourcePtr  biop;

  /* Check parameters */

  if (sdp == NULL || cfp == NULL)
    return;

  if (Seq_descr_source != sdp->choice)
    return;

  biop = sdp->data.ptrvalue;
  if (NULL == biop)
    return;

  /* Search the source for the given string */

  switch (cfp->type) {
    case eFieldTypeBioSource :
      found = DoesBioSourceContainText (biop, cfp);
      break;
    case eFieldTypeOrgModSubSource :
      found = DoSubSourcesContainText (biop, cfp);
    default:
      break;
  }

  /* If we found anything to delete then */
  /* delete it and set the dirty flag.   */

  if (TRUE == found) {
    cfp->isDirty = TRUE;
    switch (cfp->deleteLevelInt) {
    case 1 :
      break;
    case 2 :
      if (IS_Bioseq (sep))
	bsp->idx.deleteme = TRUE;
      break;
    case 3 :
      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
	bssp->idx.deleteme = TRUE;
      } else {
	if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
	  bssp = (BioseqSetPtr) bsp->idx.parentptr;
	  bssp->idx.deleteme = TRUE;
	}
      }
      break;
    default:
      break;
    }
  }

  /* Return successfully */

  return;

}


static void DeleteSimpleDescriptorsByString (BioseqPtr bsp, CharPtr deleteStr, Uint1 desc_choice, BoolPtr isDirty)
{
  SeqMgrDescContext dcontext;
  SeqDescrPtr        sdp;
  
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, desc_choice, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, desc_choice, &dcontext))
  {
    if (CheckForString (deleteStr, sdp->data.ptrvalue)) {
      if (sdp->extended != 0) {
        ((ObjValNodePtr)sdp)->idx.deleteme = TRUE;
        if (isDirty != NULL) {
          *isDirty = TRUE;
        }
      }
    }
  }
}


/*=========================================================================*/
/*                                                                         */
/* DeleteItemsByText () - Given a text string, delete all items of a given */
/*                        type that contain that string.                   */
/*                                                                         */
/*=========================================================================*/

static void DeleteItemsByText (SeqEntryPtr    sep,
			       ConvertFormPtr cfp)
{
  BioseqSetPtr      bssp;
  BioseqPtr         bsp;
  SeqMgrDescContext descContext;
  /*
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       parentSep;
  */
  SeqDescrPtr       sdp;

  /* If we have a Bioseq Set, then recurse until */
  /* we get down to an actual Bioseq.            */

  if (IS_Bioseq_set (sep)) {

    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (NULL == bssp)
      return;

    for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      DeleteItemsByText (sep, cfp);

    return;
  }

  /* If we made it this far, then we have a Bioseq */

  bsp = (BioseqPtr) sep->data.ptrvalue;

  switch (cfp->type) {
    case  eFieldTypeGene :
    case eFieldTypeCDS :
    case eFieldTypeProtein :
    case eFieldTypeRNA :
    case eFieldTypeImport :
      SeqEntryExplore (sep, (Pointer) cfp, DeleteFeaturesByText_Callback);
      break;
    case eFieldTypeBioSource :
    case eFieldTypeOrgModSubSource :
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &descContext);
      while (NULL != sdp) {
	    sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &descContext);
	    DeleteSourceByText (sdp, sep, bsp, cfp);
      }
      /*
      SeqEntryToBioSource (sep, NULL, NULL, 0, &biop);
      if (NULL == biop) {
	GetSeqEntryParent (sep, &parentptr, &parenttype);
	SeqEntryToBioSource (parentSep, NULL, NULL, 0, &biop);
      }
      DeleteSourceByText (sep, bsp, biop, cfp);
      */
      break;
    case eFieldTypeDefline :
      DeleteSimpleDescriptorsByString (bsp, cfp->deleteStr, Seq_descr_title, &(cfp->isDirty));
      break;
    case eFieldTypeCommentDescriptor :
      DeleteSimpleDescriptorsByString (bsp, cfp->deleteStr, Seq_descr_comment, &(cfp->isDirty));
      break;
    default:
      break;
  }

  /* Return succesfully */

  return;
}

/*=========================================================================*/
/*                                                                         */
/* DeleteByText_Callback () - Finds and deletes all items of a selected    */
/*                            type that contain a given text phrase.       */
/*                                                                         */
/*=========================================================================*/

static void DeleteByText_Callback (ButtoN b)
{
  ConvertFormPtr  cfp;
  SeqEntryPtr     sep;
  FieldSubfieldPtr f;

  /* Check the initial conditions and get the sequence */

  cfp = (ConvertFormPtr) GetObjectExtra (b);
  if (cfp == NULL || cfp->input_entityID == 0) {
    Remove (cfp->form);
    return;
  }

  sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
  if (sep == NULL) {
    Remove (cfp->form);
    return;
  }

  /* Get the string to search for */

  cfp->deleteStr = SaveStringFromTextNoStripSpaces (cfp->deleteText);
  if (StringHasNoText (cfp->deleteStr)){ 
    Remove (cfp->form);
    return;
  }

  /* Get the type of items to search */
  f = DialogToPointer (cfp->target_dlg);
  if (f == NULL || f->field < 0 || (f->subfield < 0 && f->subfield_list == NULL)) {
    Remove (cfp->form);
    f = FieldSubfieldFree (f);
    return;
  } else {
    cfp->type = f->field + 1;
    cfp->subtype = f->subfield;
  }

  /* Get the delete level */

  cfp->deleteLevelInt = GetValue (cfp->deleteLevel);

  /* Display the 'working' cursor */

  WatchCursor ();
  Update ();

  /* Do the search and mark and found objects for deletion */

  cfp->isDirty = FALSE;
  DeleteItemsByText (sep, cfp);

  /* Remove the window and update things */

  if (cfp->isDirty) {
    DeleteMarkedObjects (cfp->input_entityID, 0, NULL);
    ObjMgrSetDirtyFlag (cfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID, 0, 0);
  }

  ArrowCursor ();
  Update ();
  Remove (cfp->form);

  /* Return successfully */

  return;
}

static void CleanupDeleteByTextConvertForm (GraphiC g, VoidPtr data)
{
  ConvertFormPtr  cfp;

  cfp = (ConvertFormPtr) data;
  StdCleanupFormProc (g, data);
}

/*=========================================================================*/
/*                                                                         */
/* CreateDeleteByTextWindow () - Creates and then display the window for   */
/*                               getting delete by text info from the user.*/
/*                                                                         */
/*=========================================================================*/

extern Int2 LIBCALLBACK CreateDeleteByTextWindow (Pointer data)
{
  GrouP              c;
  ConvertFormPtr     cfp;
  GrouP              g;
  GrouP              h;
  GrouP              k;
  OMProcControlPtr   ompcp;
  GrouP              p;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  Boolean            allowed[eFieldType_Max];

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL)
    return OM_MSG_RET_ERROR;

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  cfp = ConvertFormNew ();
  if (cfp == NULL)
    return OM_MSG_RET_ERROR;
  cfp->set_accept_proc = SetDeleteAcceptButton;

  w = FixedWindow (-50, -33, -10, -10, "Delete By Text String",
		   StdCloseWindowProc);
  SetObjectExtra (w, cfp, CleanupDeleteByTextConvertForm);
  cfp->form = (ForM) w;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = ompcp->input_entityID;
  cfp->input_itemID = ompcp->input_itemID;
  cfp->input_itemtype = ompcp->input_itemtype;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  /* Add the popup lists */

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  k = HiddenGroup (h, 3, 0, NULL);

  StaticPrompt (k, "Delete at the level of", 0, dialogTextHeight,
		programFont, 'l');
  cfp->deleteLevel = HiddenGroup (k, 2, 0, NULL);
  RadioButton (cfp->deleteLevel, "Feature");
  RadioButton (cfp->deleteLevel, "Bioseq");
  RadioButton (cfp->deleteLevel, "Bioseq Set");
  SetValue (cfp->deleteLevel, 1);

  g = HiddenGroup (h, 3, 0, NULL);

  StaticPrompt (g, "Delete objects with the string", 0, dialogTextHeight,
		programFont, 'l');
  cfp->deleteText = DialogText (g, "", 10, ConvertFormTextCallback);
  SetObjectExtra (cfp->deleteText, cfp, NULL);

  p = HiddenGroup (h, 6, 0, NULL);
  StaticPrompt (p, "Find string in", 0, popupMenuHeight, programFont, 'l');
  MemSet (allowed, TRUE, sizeof (Boolean) * eFieldType_Max);
  allowed[eFieldTypeCommentDescriptor] = FALSE;
  allowed[eFieldTypeFeatureNote] = FALSE;
  allowed[eFieldTypePublication] = FALSE;

  cfp->target_dlg = CreateFieldSubfieldDlg (p, allowed, ChangeTargetFields, cfp);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 4, 0, NULL);
  cfp->accept = DefaultButton (c, "Accept", DeleteByText_Callback);
  SetObjectExtra (cfp->accept, cfp, NULL);
  Disable (cfp->accept);
  PushButton (c, "Cancel", StdCancelButtonProc);

  /* Line things up nicely */

  AlignObjects (ALIGN_LEFT, (HANDLE) p, (HANDLE) c, (HANDLE) h,
		(HANDLE) k, NULL);

  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (cfp->accept);
  Update ();
  return OM_MSG_RET_OK;
}

typedef Boolean (LIBCALLBACK *wantSegregateNucProtSetFunction) 
( BioseqSetPtr bssp,
  Pointer      userdata);

typedef Boolean (LIBCALLBACK *wantSegregateSequenceFunction) 
( BioseqPtr bsp,
  Pointer   userdata);
  
/*=========================================================================*/
/*                                                                         */
/* SegregateItemsRecursor () - Given a functions for determining which bioseqs     */
/*                     meet conditions, segregate bioseqs into separate    */
/*                     sets.                                               */
/*                                                                         */
/*=========================================================================*/

static void SegregateItemsRecursor 
(SeqEntryPtr                     seqlist,
 BioseqSetPtr                    set1,
 BioseqSetPtr                    set2,
 wantSegregateNucProtSetFunction do_np,
 wantSegregateSequenceFunction   do_seq,
 Pointer                         userdata
 )
{
  
  BioseqPtr         bsp;
  BioseqSetPtr      this_bssp;
  SeqEntryPtr       this_list;
  SeqEntryPtr       sep, next_sep;
  SeqEntryPtr       set1last, set2last;


  if (set1 == NULL || set2 == NULL || seqlist == NULL)
    return;

  set1last = set1->seq_set;
  while (set1last != NULL && set1last->next != NULL) {
    set1last = set1last->next;
  }
  set2last = set2->seq_set;
  while (set2last != NULL && set2last->next != NULL) {
    set2last = set2last->next;
  }

  sep = seqlist;
  while (sep != NULL) {
    next_sep = sep->next;
    if (IS_Bioseq_set (sep)) {
      this_bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (this_bssp->_class == BioseqseqSet_class_nuc_prot) {
        if (do_np != NULL && do_np (this_bssp, userdata)) {
          if (set2last == NULL) {
            set2->seq_set = sep;
          } else {
            set2last->next = sep;
          }
          set2last = sep;
        } else {
          if (set1last == NULL) {
            set1->seq_set = sep;
          } else {
            set1last->next = sep;
          }
          set1last = sep;
        }
        sep->next = NULL;
      } else {
        /* copy class types from this set if class types are not set */
        if (set1->_class == BioseqseqSet_class_genbank 
            && set2->_class == BioseqseqSet_class_genbank) {
          set1->_class = this_bssp->_class;
          set2->_class = this_bssp->_class;
        }
        /* copy descriptors from this set */
        if (this_bssp != set1) {
          ValNodeLink (&(set1->descr),
                       AsnIoMemCopy ((Pointer) this_bssp->descr,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        }
        if (this_bssp != set2) {
          ValNodeLink (&(set2->descr),
                       AsnIoMemCopy ((Pointer) this_bssp->descr,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        }
        if (this_bssp != set1 && this_bssp != set2) {
          this_bssp->descr = SeqDescrFree (this_bssp->descr);
        }
        
        this_list = this_bssp->seq_set;
        this_bssp->seq_set = NULL;
        SegregateItemsRecursor (this_list, set1, set2, do_np, do_seq, userdata);
      }
    } else if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (do_seq != NULL && do_seq (bsp, userdata)) {
        if (set2last == NULL) {
          set2->seq_set = sep;
        } else {
          set2last->next = sep;
        }
        set2last = sep;
      } else {
        if (set1last == NULL) {
          set1->seq_set = sep;
        } else {
          set1last->next = sep;
        }
        set1last = sep;
      }
      sep->next = NULL;
    }
    sep = next_sep;
  }
}


static Boolean DoFeaturesInAnnotContainText (SeqAnnotPtr sap, ConvertFormPtr cfp)
{
  SeqFeatPtr     sfp;
  GeneRefPtr     grp;
  ProtRefPtr     prp;
  RnaRefPtr      rrp;
  Boolean        found = FALSE;
  ValNodePtr     vnp;
  FieldSubfieldPtr f;
  
  if (sap == NULL || cfp == NULL) return FALSE;
  
  /* Search the requested item for the given string */

  while (sap != NULL && !found) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL && !found) {
        if (sfp->data.choice == SEQFEAT_GENE && cfp->type == eFieldTypeGene) {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
	      found = SearchGeneRef (grp, sfp, cfp);
        } else if (sfp->data.choice == SEQFEAT_CDREGION && cfp->type == eFieldTypeCDS) {
	      found = SearchCDSFeat (sfp, cfp);
        } else if (sfp->data.choice == SEQFEAT_PROT && cfp->type == eFieldTypeProtein) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
	      found = SearchProtRef (prp, sfp, cfp);
        } else if (sfp->data.choice == SEQFEAT_RNA && cfp->type == eFieldTypeRNA) {
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
	      found = SearchRnaRef (rrp, sfp, cfp);
        } else if (sfp->data.choice == SEQFEAT_IMP &&
		      cfp->type == eFieldTypeImport) {
	      found = SearchImpFeat (sfp, cfp);
        } else if (cfp->type == eFieldTypeFeatureNote) {
          f = DialogToPointer (cfp->target_dlg);
          if (f != NULL) {
            for (vnp = f->subfield_list; vnp != NULL; vnp = vnp->next) {
              if (vnp->choice == sfp->idx.subtype 
                  && CheckForString (cfp->deleteStr, sfp->comment)) {
                found = TRUE;
              }
            }
            f = FieldSubfieldFree (f);
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  return found;  
}

static Boolean DoFeaturesContainText_Callback 
(BioseqPtr      bsp,
 ConvertFormPtr cfp)
{
  
  SeqAnnotPtr    sap;

  /* Check parameters */

  if (bsp == NULL || cfp == NULL)
    return FALSE;

  /* Get the list of annotations */

  sap = bsp->annot;

  /* Search the requested item for the given string */

  return DoFeaturesInAnnotContainText (sap, cfp);
}

typedef struct objstringdata 
{
  CharPtr match;
  Boolean found;	
} ObjStringData, PNTR ObjStringPtr;

static void LIBCALLBACK AsnWriteRemoveForDCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr        pchFind;
  CharPtr        pchSource;
  ObjStringPtr   osp;

  osp = (ObjStringPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) {
	pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	pchFind = osp->match;
	if (StringSearch (pchSource, pchFind) != NULL) {
	  osp->found = TRUE;
	}
  }
}

static Boolean ObjectHasSubstring (ObjMgrTypePtr omtp, AsnIoPtr aip, Pointer ptr, ObjStringPtr osp)

{
  osp->found = FALSE;
  (omtp->asnwrite) (ptr, aip, NULL);
  return osp->found;
}

static Uint1 GetPubStatus (PubdescPtr pdp)
{
  ValNodePtr vnp;
  CitGenPtr  cgp;
  CitArtPtr  cap;
  CitJourPtr cjp;
  CitBookPtr cbp;
  CitSubPtr  csp;
  MedlineEntryPtr mlp;
  ImprintPtr ip = NULL;
  Uint1      status = 255; /* 255 is currently not a valid status */
  
  if (pdp == NULL) return status;
  
  for (vnp = pdp->pub; vnp != NULL && ip == NULL; vnp = vnp->next)
  {
  	switch (vnp->choice)
  	{
  	  case PUB_Gen:
        cgp = (CitGenPtr) vnp->data.ptrvalue;
  	    if (cgp != NULL && StringICmp (cgp->cit, "Unpublished"))
  	    {
  	  	  return PUB_STATUS_UNPUBLISHED;
  	    }
  	    break;
  	  case PUB_Article:
  	  case PUB_Medline:
  	    if (vnp->choice == PUB_Article)
  	    {
  	      cap = (CitArtPtr) vnp->data.ptrvalue;
  	    }
  	    else
  	    {
  	      cap = NULL;
  	      mlp = (MedlineEntryPtr) vnp->data.ptrvalue;
  	      if (mlp != NULL)
  	      {
  	        cap = mlp->cit;
  	      }
  	    }
  	    if (cap != NULL && cap->from == 1)
  	    {
  	      cjp = (CitJourPtr) cap->fromptr;
  	      if (cjp != NULL)
  	      {
  	      	ip = cjp->imp;
  	      }
  	    }
  	    break;
  	  case PUB_Man:
      case PUB_Book:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL)
        {
          ip = cbp->imp;
        }
        break;
  	  case PUB_Sub:
  	    csp = (CitSubPtr) vnp->data.ptrvalue;
  	    if (csp != NULL)
  	    {
  	      ip = csp->imp;
  	    }
  	    break; 
  	}
  }
  if (ip != NULL)
  {
  	status = ip->prepub;
  }
  return status;  
}

static Boolean DoesPubStatusMatch (PubdescPtr pdp, ConvertFormPtr cfp)
{
  Uint1 pub_status;
  
  if (pdp == NULL || cfp == NULL) return FALSE;
  if (cfp->subtype == 0) return TRUE;
  
  pub_status = GetPubStatus (pdp);
  
  if (cfp->subtype == PUBLICATION_PUBLISHED_FIELD 
      && pub_status == PUB_STATUS_PUBLISHED)
  {
  	return TRUE;
  }
  else if (cfp->subtype == PUBLICATION_INPRESS_FIELD
      && pub_status == PUB_STATUS_IN_PRESS)
  {
  	return TRUE;
  }
  else if (cfp->subtype == PUBLICATION_UNPUB_FIELD
      && pub_status == PUB_STATUS_UNPUBLISHED)
  {
  	return TRUE;
  }
  else
  {
  	return FALSE;
  }
}

static Boolean DoesSequenceHavePubWithText (BioseqPtr bsp, ConvertFormPtr cfp)
{
  AsnExpOptPtr      aeop;
  AsnIoPtr          aip;
  ObjStringData     osd;
  SeqMgrDescContext dcontext;
  SeqDescrPtr       sdp;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Boolean           rval = FALSE;
  ObjMgrPtr         omp;
  ObjMgrTypePtr     omtp;
  PubdescPtr        pdp;

  if (bsp == NULL || cfp == NULL) return FALSE;
  omp = ObjMgrGet ();
  if (omp == NULL) return FALSE;
  omtp = ObjMgrTypeFind (omp, OBJ_SEQDESC, NULL, NULL);
  if (omtp == NULL) return FALSE;
  
  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteRemoveForDCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &osd;
  }
  osd.match = cfp->deleteStr;

  /* look for publication descriptors */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp != NULL && !rval) {
    if (ObjectHasSubstring (omtp, aip, (Pointer) sdp, &osd)) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      if (DoesPubStatusMatch (pdp, cfp))
      {
        rval = TRUE;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }
  
  if (!rval)
  { 
    omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
    if (omtp != NULL) 
    {	
      /* look for publication features */
      sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_PUB, &fcontext);
      while (sfp != NULL && !rval) 
      {
        if (ObjectHasSubstring (omtp, aip, (Pointer) sfp, &osd))
        {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          if (DoesPubStatusMatch (pdp, cfp))
          {
            rval = TRUE;
          }
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_PUB, &fcontext);
      }
    }
  }
  
  AsnIoClose (aip); 
  return rval;
}

static Boolean DoesNucProtSetHavePubWithText (BioseqSetPtr bssp, ConvertFormPtr cfp)
{
  AsnExpOptPtr      aeop;
  AsnIoPtr          aip;
  ObjStringData     osd;
  SeqDescrPtr       sdp;
  Boolean           rval = FALSE;
  ObjMgrPtr         omp;
  ObjMgrTypePtr     omtp;
  PubdescPtr        pdp;
  
  if (bssp == NULL || cfp == NULL) return FALSE;  
  omp = ObjMgrGet ();
  if (omp == NULL) return FALSE;
  omtp = ObjMgrTypeFind (omp, OBJ_SEQDESC, NULL, NULL);
  if (omtp == NULL) return FALSE;

  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteRemoveForDCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) &osd;
  }
  osd.match = cfp->deleteStr;

  /* look for publication descriptors */
  sdp = bssp->descr;
  while (sdp != NULL && !rval) {
    if (sdp->choice == Seq_descr_pub && ObjectHasSubstring (omtp, aip, (Pointer) sdp, &osd)) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      if (DoesPubStatusMatch (pdp, cfp))
      {
        rval = TRUE;
      }
    }
    sdp = sdp->next;
  }
  
  AsnIoClose (aip); 
  return rval;
}

static Boolean DoesSimpleDescriptorForSequenceContainText (BioseqPtr bsp, CharPtr checkStr, Uint1 desc_choice)
{
  SeqMgrDescContext dcontext;
  SeqDescrPtr       sdp;
  Boolean           found = FALSE;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, desc_choice, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, desc_choice, &dcontext))
  {
    found = CheckForString (checkStr, sdp->data.ptrvalue);
  }
  return found;
}


static Boolean LIBCALLBACK DoesSequenceContainText (BioseqPtr bsp, Pointer userdata)
{
  BioSourcePtr biop;
  ConvertFormPtr cfp;
  Boolean      found = FALSE;
  SeqDescrPtr       sdp;
  SeqMgrDescContext descContext;

  cfp = (ConvertFormPtr) userdata;
  if (bsp == NULL || cfp == NULL) return FALSE;

  switch (cfp->type) {
    case eFieldTypeGene :
    case eFieldTypeCDS :
    case eFieldTypeProtein :
    case eFieldTypeRNA :
    case eFieldTypeImport :
    case eFieldTypeFeatureNote :
      found = DoFeaturesContainText_Callback (bsp, cfp);
      break;
    case eFieldTypeBioSource :
    case eFieldTypeOrgModSubSource :
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &descContext);
      while (NULL != sdp && ! found) {
        if (Seq_descr_source == sdp->choice 
            && (biop = sdp->data.ptrvalue) != NULL) {
          if (cfp->type == eFieldTypeBioSource) {
            found = DoesBioSourceContainText (biop, cfp);
          } else if (cfp->type == eFieldTypeOrgModSubSource) {
            found = DoSubSourcesContainText (biop, cfp);
          }
        }
        sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &descContext);
      }
      break;
    case eFieldTypeDefline :
      found = DoesSimpleDescriptorForSequenceContainText (bsp, cfp->deleteStr, Seq_descr_title);
      break;
    case eFieldTypeCommentDescriptor:
      found = DoesSimpleDescriptorForSequenceContainText (bsp, cfp->deleteStr, Seq_descr_comment);
      break;
  	case eFieldTypePublication :
  	  found = DoesSequenceHavePubWithText (bsp, cfp);
  	  break;
    default:
      break;
  }
  return found;
}

static Boolean LIBCALLBACK  DoesNucProtSetContainText (BioseqSetPtr bssp, Pointer userdata)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  ConvertFormPtr cfp;

  if (bssp == NULL) return FALSE;
  cfp = (ConvertFormPtr) userdata;
  if (cfp == NULL) return FALSE;
  if (cfp->type == eFieldTypePublication && DoesNucProtSetHavePubWithText (bssp, cfp))
  {
  	return TRUE;
  }
  if (DoFeaturesInAnnotContainText (bssp->annot, cfp))
  {
  	return TRUE;
  }

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (DoesSequenceContainText (bsp, (Pointer)cfp)) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void SegregateItemsByText 
(SeqEntryPtr                     seqlist,
 ConvertFormPtr                  cfp,
 BioseqSetPtr                    set1,
 BioseqSetPtr                    set2)
{
  SegregateItemsRecursor (seqlist, set1, set2, 
                  DoesNucProtSetContainText,
                  DoesSequenceContainText, (Pointer)cfp);
}

typedef void (LIBCALLBACK *segregateFunction) (
  SeqEntryPtr  seqlist,
  Pointer      userdata,
  BioseqSetPtr set1,
  BioseqSetPtr set2
);


static void CreateSetsForSegregate (BioseqSetPtr bssp, BioseqSetPtr PNTR pNewSet1, BioseqSetPtr PNTR pNewSet2)
{
  BioseqSetPtr parent_set;
  SeqEntryPtr  tmp1, tmp2, last_sep;
  BioseqSetPtr newset1, newset2;
  
  if (bssp == NULL || pNewSet1 == NULL || pNewSet2 == NULL) return;
  
  parent_set = (BioseqSetPtr)(bssp->idx.parentptr);

  if (parent_set == NULL || parent_set->seq_set == NULL) {
    /* this set has no parent, so make it the parent set, class GenBank,
     * and create two new sets using the original set class as members of this set
     */
    newset1 = BioseqSetNew ();
    if (newset1 == NULL) return;
    newset2 = BioseqSetNew ();
    if (newset2 == NULL) return;
    newset1->_class = bssp->_class;
    newset2->_class = bssp->_class;
    tmp1 = SeqEntryNew ();
    if (tmp1 == NULL) return;
    tmp1->choice = 2;
    tmp1->data.ptrvalue = (Pointer) newset1;
    tmp2 = SeqEntryNew ();
    if (tmp2 == NULL) return;
    tmp2->choice = 2;
    tmp2->data.ptrvalue = (Pointer) newset2;
    bssp->seq_set = tmp1;
    tmp1->next = tmp2;
    bssp->_class = BioseqseqSet_class_genbank;
    /* Propagate descriptors down */
    ValNodeLink (&(newset1->descr),
                 AsnIoMemCopy ((Pointer) bssp->descr,
                               (AsnReadFunc) SeqDescrAsnRead,
                               (AsnWriteFunc) SeqDescrAsnWrite));
    ValNodeLink (&(newset2->descr),
                 AsnIoMemCopy ((Pointer) bssp->descr,
                               (AsnReadFunc) SeqDescrAsnRead,
                               (AsnWriteFunc) SeqDescrAsnWrite));
    bssp->descr = SeqDescrFree (bssp->descr);
  } else {
    last_sep = parent_set->seq_set;
    newset1 = bssp;
    newset2 = BioseqSetNew ();
    if (newset2 == NULL) return;
    newset2->_class = newset1->_class;
    tmp1 = SeqEntryNew ();
    if (tmp1 == NULL) return;
    tmp1->choice = 2;
    tmp1->data.ptrvalue = (Pointer) newset2;
    while (last_sep != NULL && last_sep->next != NULL) {
      last_sep = last_sep->next;
    }
    if (last_sep == NULL) return;
    last_sep->next = tmp1;
    /* copy descriptors horizontally */
    ValNodeLink (&(newset2->descr),
                 AsnIoMemCopy ((Pointer) bssp->descr,
                               (AsnReadFunc) SeqDescrAsnRead,
                               (AsnWriteFunc) SeqDescrAsnWrite));
  }
  *pNewSet1 = newset1;
  *pNewSet2 = newset2;
}  
  

static void SegregateItemsGenericCallback
(SeqEntryPtr       sep,
 BioseqSetPtr      bssp,
 Uint2             entityID,
 Pointer           userdata,
 segregateFunction seg_func)
{
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  BioseqSetPtr      newset1;
  BioseqSetPtr      newset2;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       seqlist;

  if (sep == NULL || seg_func == NULL) return;
  /* Display the 'working' cursor */

  WatchCursor ();
  Update ();


  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  seqlist = bssp->seq_set;
  bssp->seq_set = NULL;

  CreateSetsForSegregate (bssp, &newset1, &newset2);

  /* Do the search and move sequences */
  (*seg_func)(seqlist, userdata, newset1, newset2);

  /* Remove the window and update things */
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);

  ArrowCursor ();
  Update ();
  
  /* Return successfully */
  return;
  
}
 

static Boolean OkToSegregate (OMProcControlPtr ompcp)
{
  if (ompcp == NULL
      || ompcp->input_itemtype != OBJ_BIOSEQSET 
      || ompcp->input_data == NULL) {
    Message (MSG_ERROR, "You must select a set to segregate!");
    return FALSE;
  } else {
    return TRUE;
  }
}


/*=========================================================================*/
/*                                                                         */
/* SegregateByText_Callback () - Finds and deletes all items of a selected */
/*                            type that contain a given text phrase.       */
/*                                                                         */
/*=========================================================================*/

static void SegregateByText_Callback (ButtoN b)
{
  ConvertFormPtr    cfp;
  SeqEntryPtr       sep;
  BioseqSetPtr      bssp;
  SeqEntryPtr       seqlist;
  BioseqSetPtr      newset1, newset2;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  Boolean           leaveDlgUp = FALSE;
  FieldSubfieldPtr  f;

  /* Check the initial conditions and get the sequence */
  cfp = (ConvertFormPtr) GetObjectExtra (b);
  if (cfp == NULL || cfp->input_entityID == 0 || cfp->target_set == NULL) {
    Remove (cfp->form);
    return;
  }

  sep = GetTopSeqEntryForEntityID (cfp->input_entityID); 
  if (sep == NULL) {
    Remove (cfp->form);
    return;
  }
  
  if (GetStatus (cfp->leaveDlgUp))
  {
  	leaveDlgUp = TRUE;
  }

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  bssp = cfp->target_set;
  seqlist = bssp->seq_set;
  bssp->seq_set = NULL;

  CreateSetsForSegregate (bssp, &newset1, &newset2);

  /* Get the string to search for */

  cfp->deleteStr = SaveStringFromTextNoStripSpaces (cfp->deleteText);
  if (StringHasNoText (cfp->deleteStr)){ 
    if (!leaveDlgUp)
    {
      Remove (cfp->form);
    }
    return;
  }

  /* Get the type of items to search */
  f = DialogToPointer (cfp->target_dlg);
  if (f == NULL || f->field < 0 || (f->subfield < 0 && f->subfield_list == NULL)) {
    f = FieldSubfieldFree (f);
    if (!leaveDlgUp)
    {
      Remove (cfp->form);
    }
    return;
  }
  
  cfp->type = f->field;
  cfp->subtype = f->subfield;
  f = FieldSubfieldFree (f);

  /* Display the 'working' cursor */

  WatchCursor ();
  Update ();

  /* Do the search and move sequences */
  SegregateItemsByText (seqlist, cfp, newset1, newset2);

  /* Remove the window and update things */
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  ObjMgrSetDirtyFlag (cfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();
  if (!leaveDlgUp)
  {
    Remove (cfp->form);
  }

  /* Return successfully */
  return;
}




/*-------------------------------------------------------------------------*/
/*                                                                         */
/* SetSegregateAcceptButton () -- Enable/Disable the Accept button depending  */
/*                             on the condition of other window objects.   */
/*                                                                         */
/*-------------------------------------------------------------------------*/

static void SetSegregateAcceptButton (Pointer data)

{
  ConvertFormPtr  cfp;
  FieldSubfieldPtr f;

  cfp = (ConvertFormPtr) data;
  if (cfp == NULL)
    return;

  /* Disable if a target has not been selected */
  f = DialogToPointer (cfp->target_dlg);
  if (f == NULL || f->field < 0 || (f->subfield < 0 && f->subfield_list == NULL)) {
    SafeDisable (cfp->accept);
    f = FieldSubfieldFree (f);
    return;
  }
  f = FieldSubfieldFree (f);

  if (cfp->deleteText != NULL) {
    /* Disable if there is no delete string */ 
    if (TextHasNoText (cfp->deleteText)) {
      SafeDisable (cfp->accept);
      return;
    }
  }

  /* If we made it to here, then we passed all */
  /* test and can enable the accept button.    */

  SafeEnable (cfp->accept);
}

static void CleanupSegregatePage (GraphiC g, VoidPtr data)

{
  ConvertFormPtr  cfp;

  cfp = (ConvertFormPtr) data;
  if (cfp != NULL) {
  }
  StdCleanupFormProc (g, data);
}

/*=========================================================================*/
/*                                                                         */
/* CreateSegregateByTextWindow () - Creates and then displays the window   */
/*                        for getting segregate by text info from the user.*/
/*                                                                         */
/*=========================================================================*/

extern Int2 LIBCALLBACK CreateSegregateByTextWindow (Pointer data)
{
  Boolean            allowed[eFieldType_Max];
  GrouP              c;
  ConvertFormPtr     cfp;
  GrouP              g;
  GrouP              h;
  OMProcControlPtr   ompcp;
  GrouP              p;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (!OkToSegregate(ompcp)) {
    return OM_MSG_RET_ERROR;
  }

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  cfp = ConvertFormNew ();
  if (cfp == NULL)
    return OM_MSG_RET_ERROR;
  cfp->target_set = (BioseqSetPtr)ompcp->input_data;
  cfp->set_accept_proc = SetSegregateAcceptButton;

  w = FixedWindow (-50, -33, -10, -10, "Segregate By Text String",
		   StdCloseWindowProc);
  SetObjectExtra (w, cfp, CleanupSegregatePage);
  cfp->form = (ForM) w;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = ompcp->input_entityID;
  cfp->input_itemID = ompcp->input_itemID;
  cfp->input_itemtype = ompcp->input_itemtype;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  /* Add the popup lists */

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 3, 0, NULL);

  StaticPrompt (g, "Segregate sequences with the string", 0, dialogTextHeight,
		programFont, 'l');
  cfp->deleteText = DialogText (g, "", 10, ConvertFormTextCallback);
  SetObjectExtra (cfp->deleteText, cfp, NULL);

  p = HiddenGroup (h, 6, 0, NULL);
  StaticPrompt (p, "Find string in", 0, popupMenuHeight, programFont, 'l');
  MemSet (allowed, TRUE, sizeof (Boolean) * eFieldType_Max);
  allowed[eFieldTypeCommentDescriptor] = TRUE;

  cfp->target_dlg = CreateFieldSubfieldDlg (p, allowed, ChangeTargetFields, cfp);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 4, 0, NULL);
  cfp->accept = DefaultButton (c, "Accept", SegregateByText_Callback);
  SetObjectExtra (cfp->accept, cfp, NULL);
  Disable (cfp->accept);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cfp->leaveDlgUp = CheckBox (c, "Leave Dialog Up", NULL);

  /* Line things up nicely */

  AlignObjects (ALIGN_LEFT, (HANDLE) p, (HANDLE) c, (HANDLE) h, NULL);

  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (cfp->accept);
  Update ();
  return OM_MSG_RET_OK;
}


typedef struct segregatefeatdata {
  FEATURE_FORM_BLOCK

  PopuP        type_popup;
  ValNodePtr   type_list;
  ButtoN       accept;
  
  BioseqSetPtr target_set;
  Uint2        segregate_type;
  Boolean      is_feat;
} SegregateFeatData, PNTR SegregateFeatPtr;

static Boolean LIBCALLBACK DoesSequenceContainFeatureType (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        feat;
  SegregateFeatPtr  sfp;
  
  sfp = (SegregateFeatPtr) userdata;
  if (sfp == NULL || bsp == NULL) return FALSE;
  
  feat = NULL;
  while ((feat = SeqMgrGetNextFeature (bsp, feat, 0, 0, &context)) != NULL)
  {
  	if (feat->idx.subtype == sfp->segregate_type)
  	{
  	  return TRUE;
  	}
  }
  return FALSE;
}

static Boolean LIBCALLBACK DoesNucProtSetContainFeatureType (BioseqSetPtr bssp, Pointer userdata)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  SegregateFeatPtr  sfp;
  
  sfp = (SegregateFeatPtr) userdata;
  if (sfp == NULL || bssp == NULL) return FALSE;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (DoesSequenceContainFeatureType (bsp, sfp)) {
        return TRUE;
      }
    }
  }
  return FALSE;	
}

static Boolean LIBCALLBACK DoesSequenceContainDescriptorType (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrDescContext context;
  SeqDescPtr        desc;
  SegregateFeatPtr  sfp;
  
  sfp = (SegregateFeatPtr) userdata;
  if (sfp == NULL || bsp == NULL) return FALSE;
  
  if((desc = SeqMgrGetNextDescriptor (bsp, NULL, sfp->segregate_type, &context)) != NULL)
  {
	return TRUE;
  }
  return FALSE;
}

typedef struct checkdescdata {
  Uint2   segregate_type;
  Boolean found;
} CheckDescData, PNTR CheckDescPtr;

static void DoesSetContainDescriptorType_Callback (SeqDescPtr sdp, Pointer userdata)
{
  CheckDescPtr p;
  
  if (sdp == NULL || userdata == NULL) return;
  p = (CheckDescPtr) userdata;
  if (p->found) return;
  if (sdp->choice == p->segregate_type) p->found = TRUE;
}

static Boolean LIBCALLBACK DoesNucProtSetContainDescriptorType (BioseqSetPtr bssp, Pointer userdata)
{
  CheckDescData d;
  SegregateFeatPtr  sfp;
  
  sfp = (SegregateFeatPtr) userdata;
  if (sfp == NULL || bssp == NULL) return FALSE;

  d.found = FALSE;
  d.segregate_type = sfp->segregate_type;
  VisitDescriptorsInSet (bssp, &d, DoesSetContainDescriptorType_Callback);
  return d.found;	
}

/*=========================================================================*/
/*                                                                         */
/* SegregateItemsByFeature () - Given a feature type, move bioseqs         */
/*                        containing those features to a new popset.       */
/*                                                                         */
/*=========================================================================*/

static void SegregateItemsByFeatureOrDescriptor 
(SeqEntryPtr      seqlist,
 SegregateFeatPtr sfp,
 BioseqSetPtr     set1,
 BioseqSetPtr     set2)
{
  

  if (sfp == NULL || set1 == NULL || set2 == NULL || seqlist == NULL)
    return;
  
  if (sfp->is_feat)
  {
  	SegregateItemsRecursor (seqlist, set1, set2, 
  	                               DoesNucProtSetContainFeatureType,
  	                               DoesSequenceContainFeatureType,
  	                               (Pointer) sfp);
  } 
  else 
  {
  	SegregateItemsRecursor (seqlist, set1, set2, 
  	                               DoesNucProtSetContainDescriptorType,
  	                               DoesSequenceContainDescriptorType,
  	                               (Pointer) sfp);
  }
  
}


/*=========================================================================*/
/*                                                                         */
/* SegregateByFeatureOrDescriptor_Callback () - Segregates sequences that  */
/*                            contain a selected feature.                  */
/*                                                                         */
/*=========================================================================*/

static void SegregateByFeatureOrDescriptor_Callback (ButtoN b)
{
  SegregateFeatPtr  sfp;
  SeqEntryPtr       sep;
  UIEnum            val;
  BioseqSetPtr      bssp;
  SeqEntryPtr       seqlist;
  BioseqSetPtr      newset1, newset2;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  ValNodePtr        vnp;

  /* Check the initial conditions and get the sequence */
  sfp = (SegregateFeatPtr) GetObjectExtra (b);
  if (sfp == NULL || sfp->input_entityID == 0 || sfp->target_set == NULL) {
    Remove (sfp->form);
    return;
  }

  sep = GetTopSeqEntryForEntityID (sfp->input_entityID); 
  if (sep == NULL) {
    Remove (sfp->form);
    return;
  }

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  bssp = sfp->target_set;
  seqlist = bssp->seq_set;
  bssp->seq_set = NULL;

  CreateSetsForSegregate (bssp, &newset1, &newset2);

  /* Get the feature to look for */
  val = GetValue (sfp->type_popup);
  for (vnp = sfp->type_list; vnp != NULL && val > 1; vnp = vnp->next, val--)
  {  	
  }
  if (vnp == NULL || val != 1)
  {
    Remove (sfp->form);
    return;
  }
  sfp->segregate_type = vnp->choice;
  
  /* Display the 'working' cursor */

  WatchCursor ();
  Update ();

  /* Do the search and move sequences */
  SegregateItemsByFeatureOrDescriptor (seqlist, sfp, newset1, newset2);

  /* Remove the window and update things */
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  ObjMgrSetDirtyFlag (sfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, sfp->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();
  Remove (sfp->form);

  /* Return successfully */
  return;
}


/*=========================================================================*/
/*                                                                         */
/* CreateSegregateByFeatureWindow () - Creates and then displays the window*/
/*                        for getting segregate by text info from the user.*/
/*                                                                         */
/*=========================================================================*/

static Int2 LIBCALLBACK CreateSegregateByFeatureOrDescriptorWindow (Pointer data, Boolean is_feat)
{
  GrouP              c;
  SegregateFeatPtr   sfp;
  GrouP              g;
  GrouP              h;
  OMProcControlPtr   ompcp;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  ValNodePtr         vnp;

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (!OkToSegregate(ompcp)) {
    return OM_MSG_RET_ERROR;
  }

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  sfp = (SegregateFeatPtr) MemNew (sizeof (SegregateFeatData));
  if (sfp == NULL)
    return OM_MSG_RET_ERROR;
  sfp->is_feat = is_feat;

  if (sfp->is_feat)
  {
    w = FixedWindow (-50, -33, -10, -10, "Segregate By Feature",
	                 StdCloseWindowProc);
  }
  else
  {
    w = FixedWindow (-50, -33, -10, -10, "Segregate By Descriptor",
	                 StdCloseWindowProc);
  }
		   
  SetObjectExtra (w, sfp, StdCleanupFormProc);
  sfp->form = (ForM) w;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    sfp->appmessage = sepp->handleMessages;
  }

  sfp->input_entityID = ompcp->input_entityID;
  sfp->input_itemID = ompcp->input_itemID;
  sfp->input_itemtype = ompcp->input_itemtype;
  sfp->target_set = (BioseqSetPtr)ompcp->input_data;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    sfp->appmessage = sepp->handleMessages;
  }

  /* Add the popup lists */

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 3, 0, NULL);

  if (sfp->is_feat)
  {
    StaticPrompt (g, "Segregate sequences with the feature", 0, dialogTextHeight,
	              programFont, 'l');
    sfp->type_list = BuildFeatureValNodeList (TRUE, NULL, 0, TRUE, FALSE);
  }
  else
  {
    StaticPrompt (g, "Segregate sequences with the descriptor", 0, dialogTextHeight,
	              programFont, 'l');
	sfp->type_list = BuildDescriptorValNodeList ();
  }

  sfp->type_popup = PopupList (g, TRUE, NULL);
  SetObjectExtra (sfp->type_popup, sfp, NULL);
  for (vnp = sfp->type_list; vnp != NULL; vnp = vnp->next)
  {
    PopupItem (sfp->type_popup, (CharPtr) vnp->data.ptrvalue);
  }
  SetValue (sfp->type_popup, 1);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 4, 0, NULL);
  sfp->accept = DefaultButton (c, "Accept", SegregateByFeatureOrDescriptor_Callback);
  SetObjectExtra (sfp->accept, sfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  /* Line things up nicely */

  AlignObjects (ALIGN_LEFT, (HANDLE) g, (HANDLE) c, (HANDLE) h, NULL);

  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (sfp->accept);
  Update ();
  return OM_MSG_RET_OK;
}

extern Int2 LIBCALLBACK CreateSegregateByFeatureWindow (Pointer data)
{
  return CreateSegregateByFeatureOrDescriptorWindow (data, TRUE);
}

extern Int2 LIBCALLBACK CreateSegregateByDescriptorWindow (Pointer data)
{
  return CreateSegregateByFeatureOrDescriptorWindow (data, FALSE);
}

typedef struct segregatemoltypedata
{
  FEATURE_FORM_BLOCK

  PopuP        type_popup;
  ButtoN       use_mol_type;
  PopuP        class_popup;
  ButtoN       use_mol_class;
  ButtoN       accept;
  Uint1        moltype;
  Uint1        molclass;
  BioseqSetPtr target_set;
} SegregateMolTypeData, PNTR SegregateMolTypePtr;

static ENUM_ALIST(molinfo_biomol_alist)
  {"Genomic DNA or RNA",     1},
  {"Precursor RNA",          2},
  {"mRNA [cDNA]",            3},
  {"Ribosomal RNA",          4},
  {"Transfer RNA",           5},
  {"Peptide",                8},
  {"Other-Genetic",          9},
  {"Genomic-mRNA",          10},
  {"cRNA",                  11},
  {"Transcribed RNA",       13},
  {"Other",                255},
END_ENUM_ALIST

static ENUM_ALIST(mol_class_alist)
{"DNA",             Seq_mol_dna},
{"RNA",             Seq_mol_rna},
{"Protein",         Seq_mol_aa},
{"Nucleotide",      Seq_mol_na},
{"Other",           Seq_mol_other},
END_ENUM_ALIST

static Boolean LIBCALLBACK DoesSequenceHaveMoleculeType (BioseqPtr bsp, Pointer userdata)
{
  SegregateMolTypePtr  smp;
  ValNodePtr           sdp;
  MolInfoPtr           mip;
  
  smp = (SegregateMolTypePtr) userdata;	
  if (bsp == NULL || smp == NULL) return FALSE;
  
  if (GetStatus (smp->use_mol_class) && bsp->mol != smp->molclass)
  {
    return FALSE;
  }

  if (GetStatus (smp->use_mol_type))
  {
    sdp = bsp->descr;
    while (sdp != NULL)
    {
      if (sdp->choice == Seq_descr_molinfo && sdp->data.ptrvalue != NULL) 
      {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip->biomol == smp->moltype)
        {
      	  return TRUE;
        }
      }
      sdp = sdp->next;
    }
    return FALSE;
  }
  return TRUE;
}

typedef struct lookformoltype
{
  Boolean found;
  Uint1   moltype;  
  Uint1   molclass;
} LookForMolTypeData, PNTR LookForMolTypePtr;

static void DoesNucProtSetHaveMoleculeTypeCallback (SeqDescPtr sdp, Pointer userdata)
{
  LookForMolTypePtr l;
  MolInfoPtr           mip;
  
  if (sdp == NULL || userdata == NULL) return;
  l = (LookForMolTypePtr) userdata;
  
  if (sdp->choice == Seq_descr_molinfo && sdp->data.ptrvalue != NULL)
  {
  	mip = (MolInfoPtr) sdp->data.ptrvalue;
  	if (mip->biomol == l->moltype)
  	{
  	  l->found = TRUE;
  	}
  }
}

typedef struct lookformolclass
{
  Boolean found;
  Uint1   molclass;
} LookForMolClassData, PNTR LookForMolClassPtr;

static void FindMolClassCallback (BioseqPtr bsp, Pointer userdata)
{
  LookForMolClassPtr lcp;
  
  if (bsp == NULL || userdata == NULL) return;
  lcp = (LookForMolClassPtr) userdata;
  
  if (bsp->mol == lcp->molclass)
  {
    lcp->found = TRUE;
  }
}

static Boolean LIBCALLBACK DoesNucProtSetHaveMoleculeType (BioseqSetPtr bssp, Pointer userdata)
{
  LookForMolTypeData   lm;
  LookForMolClassData  lc;
  SegregateMolTypePtr  smp;

  if (bssp == NULL || userdata == NULL) return FALSE;
  smp = (SegregateMolTypePtr) userdata;
  
  if (GetStatus (smp->use_mol_class))
  {
    lc.found = FALSE;
    lc.molclass = smp->molclass;
    VisitBioseqsInSet (bssp, &lc, FindMolClassCallback);
    if (!lc.found) return FALSE;
  }
  
  if (GetStatus (smp->use_mol_type))
  {
    lm.moltype = smp->moltype;
    lm.found = FALSE;
    VisitDescriptorsInSet (bssp, &lm, DoesNucProtSetHaveMoleculeTypeCallback);
    if (!lm.found) return FALSE;
  }
  return TRUE;
}

static void LIBCALLBACK SegregateByMoleculeType 
(SeqEntryPtr  seqlist,
 Pointer      userdata,
 BioseqSetPtr set1,
 BioseqSetPtr set2)
{
  SegregateItemsRecursor (seqlist, set1, set2, 
                          DoesNucProtSetHaveMoleculeType,
                          DoesSequenceHaveMoleculeType, userdata);
}


/*=========================================================================*/
/*                                                                         */
/* SegregateByMoleculeType_Callback () - Segregates sequences that  */
/*                            have a selected molecule type.                  */
/*                                                                         */
/*=========================================================================*/

static void SegregateByMoleculeType_Callback (ButtoN b)
{
  SegregateMolTypePtr  smp;
  SeqEntryPtr          sep;
  UIEnum               val;

  /* Check the initial conditions and get the sequence */
  smp = (SegregateMolTypePtr) GetObjectExtra (b);
  if (smp == NULL || smp->input_entityID == 0 || smp->target_set == NULL) {
    Remove (smp->form);
    return;
  }

  sep = GetTopSeqEntryForEntityID (smp->input_entityID); 
  if (sep == NULL) {
    Remove (smp->form);
    return;
  }
  
  if (!GetEnumPopup (smp->type_popup, molinfo_biomol_alist, &val)) 
  {
    return;
  }
  smp->moltype = val;
  if (!GetEnumPopup (smp->class_popup, mol_class_alist, &val)) 
  {
    return;
  }
  smp->molclass = val;

  SegregateItemsGenericCallback (sep, smp->target_set, smp->input_entityID,
                                 (Pointer) smp, SegregateByMoleculeType);

  Remove (smp->form);

  /* Return successfully */
  return;
}

static void EnableMolInfoPopups (ButtoN b)
{
  SegregateMolTypePtr smp;
  Boolean             ok_to_accept = FALSE;
  
  smp = (SegregateMolTypePtr) GetObjectExtra (b);
  if (smp == NULL) return;
  
  if (GetStatus (smp->use_mol_type))
  {
    Enable (smp->type_popup);
    ok_to_accept = TRUE;
  }
  else
  {
    Disable (smp->type_popup);
  }
  if (GetStatus (smp->use_mol_class))
  {
    Enable (smp->class_popup);
    ok_to_accept = TRUE;
  }
  else
  {
    Disable (smp->class_popup);
  }  
  if (ok_to_accept)
  {
    Enable (smp->accept);
  }
  else
  {
    Disable (smp->accept);
  }
}

/*=========================================================================*/
/*                                                                         */
/* CreateSegregateByFeatureWindow () - Creates and then displays the window*/
/*                        for getting segregate by text info from the user.*/
/*                                                                         */
/*=========================================================================*/

extern Int2 LIBCALLBACK CreateSegregateByMoleculeTypeWindow (Pointer data)
{
  GrouP               c;
  SegregateMolTypePtr smp;
  GrouP               g, k;
  GrouP               h;
  OMProcControlPtr    ompcp;
  StdEditorProcsPtr   sepp;
  WindoW              w;

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (!OkToSegregate(ompcp)) {
    return OM_MSG_RET_ERROR;
  }

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  smp = (SegregateMolTypePtr) MemNew (sizeof (SegregateMolTypeData));
  if (smp == NULL)
    return OM_MSG_RET_ERROR;
  
  w = FixedWindow (-50, -33, -10, -10, "Segregate By Molecule Type",
	               StdCloseWindowProc);
		   
  SetObjectExtra (w, smp, StdCleanupFormProc);
  smp->form = (ForM) w;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    smp->appmessage = sepp->handleMessages;
  }

  smp->input_entityID = ompcp->input_entityID;
  smp->input_itemID = ompcp->input_itemID;
  smp->input_itemtype = ompcp->input_itemtype;
  smp->target_set = (BioseqSetPtr)ompcp->input_data;

  /* Add the popup lists */

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 1, 0, NULL);

  StaticPrompt (g, "Segregate sequences with:", 0, dialogTextHeight,
	              programFont, 'l');
	k = HiddenGroup (g, 2, 0, NULL);
	smp->use_mol_type = CheckBox (k, "Molecule Type", EnableMolInfoPopups);
  SetObjectExtra (smp->use_mol_type, smp, NULL);
  smp->type_popup = PopupList (k, TRUE, NULL);
  InitEnumPopup (smp->type_popup, molinfo_biomol_alist, NULL);
  SetValue (smp->type_popup, 1);
  SetStatus (smp->use_mol_type, FALSE);
  Disable (smp->type_popup);
	k = HiddenGroup (g, 2, 0, NULL);
  smp->use_mol_class = CheckBox (k, "Molecule Class", EnableMolInfoPopups);
  SetObjectExtra (smp->use_mol_class, smp, NULL);
  smp->class_popup = PopupList (k, TRUE, NULL);
  InitEnumPopup (smp->class_popup, mol_class_alist, NULL);
  SetValue (smp->class_popup, 1);
  SetStatus (smp->use_mol_class, FALSE);
  Disable (smp->class_popup);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 4, 0, NULL);
  smp->accept = DefaultButton (c, "Accept", SegregateByMoleculeType_Callback);
  SetObjectExtra (smp->accept, smp, NULL);
  Disable (smp->accept);
  PushButton (c, "Cancel", StdCancelButtonProc);

  /* Line things up nicely */

  AlignObjects (ALIGN_LEFT, (HANDLE) g, (HANDLE) c, (HANDLE) h, NULL);

  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (smp->accept);
  Update ();
  return OM_MSG_RET_OK;
}

static Boolean LIBCALLBACK DoesSequenceHaveID (BioseqPtr bsp, Pointer userdata)
{
  ConvertFormPtr cfp;
  Boolean      found = FALSE;
  ValNodePtr     vnp;

  cfp = (ConvertFormPtr) userdata;
  if (bsp == NULL || cfp == NULL) return FALSE;

  for (vnp = cfp->id_list; vnp != NULL && !found; vnp = vnp->next) {
    if (SeqIdIn (vnp->data.ptrvalue, bsp->id)) {
      found = TRUE;
    }
  }
  return found;
}

static Boolean LIBCALLBACK  DoesNucProtSetContainID (BioseqSetPtr bssp, Pointer userdata)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  ConvertFormPtr cfp;

  if (bssp == NULL) return FALSE;
  cfp = (ConvertFormPtr) userdata;
  if (cfp == NULL) return FALSE;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (DoesSequenceHaveID (bsp, (Pointer)cfp)) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void SegregateItemsById 
(SeqEntryPtr                     seqlist,
 ConvertFormPtr                  cfp,
 BioseqSetPtr                    set1,
 BioseqSetPtr                    set2)
{
  SegregateItemsRecursor (seqlist, set1, set2, 
                  DoesNucProtSetContainID,
                  DoesSequenceHaveID, (Pointer)cfp);
}

static void SetSegregateByIDAcceptButton (Pointer data)

{
  ConvertFormPtr  cfp;
  CharPtr         id_str;

  cfp = (ConvertFormPtr) data;
  if (cfp == NULL)
    return;

  id_str = SaveStringFromText (cfp->deleteText);
  if (StringHasNoText (id_str)) {
    SafeDisable (cfp->accept);
  } else {
    SafeEnable (cfp->accept);
  }
}

/*=========================================================================*/
/*                                                                         */
/* SegregateByID_Callback () - Finds and deletes all items of a selected */
/*                            type that contain a given text phrase.       */
/*                                                                         */
/*=========================================================================*/

static void SegregateById_Callback (ButtoN b)
{
  ConvertFormPtr    cfp;
  SeqEntryPtr       sep;
  BioseqSetPtr      bssp;
  SeqEntryPtr       seqlist;
  BioseqSetPtr      newset1, newset2;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  Boolean           leaveDlgUp = FALSE;
  CharPtr           id_str;

  /* Check the initial conditions and get the sequence */
  cfp = (ConvertFormPtr) GetObjectExtra (b);
  if (cfp == NULL || cfp->input_entityID == 0 || cfp->target_set == NULL) {
    Remove (cfp->form);
    return;
  }

  sep = GetTopSeqEntryForEntityID (cfp->input_entityID); 
  if (sep == NULL) {
    Remove (cfp->form);
    return;
  }
  
  if (GetStatus (cfp->leaveDlgUp))
  {
  	leaveDlgUp = TRUE;
  }

  /* Get the string to search for */
  id_str = SaveStringFromText(cfp->deleteText);
  if (StringHasNoText (id_str)){ 
    if (!leaveDlgUp)
    {
      Remove (cfp->form);
    }
    id_str = MemFree (id_str);
    return;
  }
  
  cfp->id_list = ParseAccessionNumberListFromString (id_str, SeqMgrGetSeqEntryForData (cfp->target_set));
  id_str = MemFree (id_str);
  if (cfp->id_list == NULL) {
    Message (MSG_ERROR, "No IDs specified!");
    return;
  } 

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  bssp = cfp->target_set;
  seqlist = bssp->seq_set;
  bssp->seq_set = NULL;

  CreateSetsForSegregate (bssp, &newset1, &newset2);

  /* Display the 'working' cursor */

  WatchCursor ();
  Update ();

  /* Do the search and move sequences */
  SegregateItemsById (seqlist, cfp, newset1, newset2);

  cfp->id_list = FreeSeqIdList (cfp->id_list);
  /* Remove the window and update things */
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  ObjMgrSetDirtyFlag (cfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();
  if (!leaveDlgUp)
  {
    Remove (cfp->form);
  }

  /* Return successfully */
  return;
}

/*=========================================================================*/
/*                                                                         */
/* CreateSegregateByIdWindow () - Creates and then displays the window   */
/*                        for getting segregate by ID info from the user.*/
/*                                                                         */
/*=========================================================================*/

extern Int2 LIBCALLBACK CreateSegregateByIdWindow (Pointer data)
{
  GrouP              c;
  ConvertFormPtr     cfp;
  GrouP              g;
  GrouP              h;
  OMProcControlPtr   ompcp;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  /* Check parameters and get a pointer to the current data */

  ompcp = (OMProcControlPtr) data;
  if (!OkToSegregate(ompcp)) {
    return OM_MSG_RET_ERROR;
  }

  /* Create a new window, and a struct */
  /* to pass around the data in.       */

  cfp = ConvertFormNew ();
  if (cfp == NULL)
    return OM_MSG_RET_ERROR;
  cfp->target_set = (BioseqSetPtr)ompcp->input_data;
  cfp->set_accept_proc = SetSegregateByIDAcceptButton;

  w = FixedWindow (-50, -33, -10, -10, "Segregate By ID",
		   StdCloseWindowProc);
  SetObjectExtra (w, cfp, CleanupSegregatePage);
  cfp->form = (ForM) w;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = ompcp->input_entityID;
  cfp->input_itemID = ompcp->input_itemID;
  cfp->input_itemtype = ompcp->input_itemtype;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  /* Add the popup lists */

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 3, 0, NULL);

  StaticPrompt (g, "Segregate sequences with the following IDs", 0, dialogTextHeight,
		programFont, 'l');
  cfp->deleteText = DialogText (g, "", 10, ConvertFormTextCallback);
  SetObjectExtra (cfp->deleteText, cfp, NULL);

  /* Add Accept and Cancel buttons */

  c = HiddenGroup (h, 4, 0, NULL);
  cfp->accept = DefaultButton (c, "Accept", SegregateById_Callback);
  SetObjectExtra (cfp->accept, cfp, NULL);
  Disable (cfp->accept);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cfp->leaveDlgUp = CheckBox (c, "Leave Dialog Up", NULL);

  /* Line things up nicely */

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);

  /* Display the window now */

  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (cfp->accept);
  Update ();
  return OM_MSG_RET_OK;
}


static void SearchAndExciseTextInside (CharPtr PNTR strptr, ConvertFormPtr cfp)

{
  CharPtr  lft = NULL;
  CharPtr  rgt = NULL;
  CharPtr  string;
  CharPtr  next_delete = NULL;
  Int4     strLen;

  if (strptr == NULL || cfp == NULL || cfp->textportion == NULL) return;
  string = *strptr;
  if (string == NULL) return;

  FindTextPortionXInString (string, cfp->textportion, &lft, &strLen);
  if (lft == NULL) return;
  rgt = lft + strLen;
  next_delete = lft;
  if (lft < rgt) {  /* No text to delete */
    while (*rgt != 0) {
      *lft = *rgt;
      lft++;
      rgt++;
    }
    *lft = '\0';
  }
  if (cfp->repeat_remove && next_delete != NULL && *next_delete != 0) {
    SearchAndExciseTextInside (&next_delete, cfp);
  }
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* RemOutTxt_SearchAndExcise () -- Removes all text from a string that */
/*                                 does not match a specified substring*/
/*                                                                     */
/*---------------------------------------------------------------------*/

static void SearchAndExciseTextOutside (CharPtr sourceStr, ConvertFormPtr cfp)
{
  CharPtr  leftEnd = NULL;
  CharPtr  rightEnd = NULL;
  Int4     strLen;
  Int4     i;

  /* Check parameters */

  if (cfp == NULL || cfp->textportion == NULL || sourceStr == NULL)
  {
    return;
  }

  FindTextPortionXInString (sourceStr, cfp->textportion, &leftEnd, &strLen);

  if (leftEnd != NULL) 
  {
    rightEnd = leftEnd + strLen;
  }

  if (NULL == leftEnd)  /* String not found */
  {
    if (DO_NOTHING == cfp->ifNotFound)
    {
    	return;
    }
    else
    {
      sourceStr [0] = '\0';
      return;
    }
  }
  else
  {
    /* End the original string at rightEnd and  */
    /* then move everything starting at leftEnd */
    /* to the beginning of the original.        */

    rightEnd[0] = '\0';
    if (leftEnd == sourceStr)
      return;

    strLen = StringLen (leftEnd);

    for (i = 0; i <= strLen; i++)
      sourceStr[i] = leftEnd[i];

    sourceStr[i] = '\0';
  }
}


static void SearchAndExciseText (CharPtr str, ConvertFormPtr cfp)
{
  if (str == NULL || cfp == NULL) return;
  if (cfp->remove_inside) 
  {
    SearchAndExciseTextInside (&str, cfp);
  } 
  else
  {
    SearchAndExciseTextOutside (str, cfp);
  }
}


static void RemoveAFeatureText (SeqFeatPtr sfp, Pointer mydata)
{
  ConvertFormPtr  cfp;
  GBQualPtr       gbqp;
  GeneRefPtr      grp;
  ProtRefPtr      prp;
  RnaRefPtr       rrp;
  CharPtr         str;
  ValNodePtr      vnp;

  if (sfp == NULL) return;
  cfp = (ConvertFormPtr) mydata;
  if (cfp == NULL) return;

  if (sfp->data.choice == SEQFEAT_GENE && cfp->type == eFieldTypeGene) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      switch (cfp->subtype) {
        case 1 :
          SearchAndExciseText (grp->locus, cfp);
          break;
        case 2 :
          SearchAndExciseText (grp->desc, cfp);
          break;
        case 3 :
          SearchAndExciseText (grp->allele, cfp);
          break;
        case 4 :
          SearchAndExciseText (grp->maploc, cfp);
          break;
        case 5 :
          SearchAndExciseText (grp->locus_tag, cfp);
          break;
        case 6 :
          for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
            str = (CharPtr) vnp->data.ptrvalue;
            SearchAndExciseText (str, cfp);
          }
          break;
        case 7 :
          SearchAndExciseText (sfp->comment, cfp);
          break;
        default :
          break;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_CDREGION && cfp->type == eFieldTypeCDS) {
    switch (cfp->subtype) {
      case CDS_COMMENT :
        SearchAndExciseText (sfp->comment, cfp);
        break;
      case CDS_GENE_XREF :
        break;
      default :
        break;
    }
  } else if (sfp->data.choice == SEQFEAT_PROT && cfp->type == eFieldTypeProtein) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp != NULL) {
      switch (cfp->subtype) {
        case 1 :
          for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
            str = (CharPtr) vnp->data.ptrvalue;
            SearchAndExciseText (str, cfp);
          }
          break;
        case 2 :
          SearchAndExciseText (prp->desc, cfp);
          break;
        case 3 :
          for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
            str = (CharPtr) vnp->data.ptrvalue;
            SearchAndExciseText (str, cfp);
          }
          break;
        case 4 :
          for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
            str = (CharPtr) vnp->data.ptrvalue;
            SearchAndExciseText (str, cfp);
          }
          break;
        case 5 :
          SearchAndExciseText (sfp->comment, cfp);
          break;
        default :
          break;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_RNA && cfp->type == eFieldTypeRNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL) {
      switch (cfp->subtype) {
        case 1 :
          if (rrp->ext.choice == 0 || rrp->ext.choice == 1) {
            rrp->ext.choice = 1;
            str = (CharPtr) rrp->ext.value.ptrvalue;
            SearchAndExciseText (str, cfp);
          }
          break;
        case 2 :
          SearchAndExciseText (sfp->comment, cfp);
          break;
        case 3 :
        default :
          break;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_IMP && cfp->type == eFieldTypeImport) {
    switch (cfp->subtype) {
      case IMPORT_GBQUAL_FIELD :
        gbqp = sfp->qual;
        while (NULL != gbqp)
        {
          if (NULL != gbqp->val)
            SearchAndExciseText (gbqp->val, cfp);
          gbqp = gbqp->next;
        }
        break;
      case IMPORT_COMMENT_FIELD :
        SearchAndExciseText (sfp->comment, cfp);
        break;
      default :
        break;
    }
	}
}


static void RemoveOrgModText (BioSourcePtr biop, ConvertFormPtr cfp)
{
  OrgModPtr mod, prev_mod = NULL, next_mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL || cfp == NULL) return;

  mod = biop->org->orgname->mod;
  while (mod != NULL) {
    next_mod = mod->next;
    if (mod->subtype == cfp->subtype) {
	    SearchAndExciseText (mod->subname, cfp);
      if (StringHasNoText (mod->subname) && !IsNonTextModifier (GetOrgModQualName (mod->subtype))) {
        if (prev_mod == NULL) {
          biop->org->orgname->mod = mod->next;
        } else {
          prev_mod->next = mod->next;
        }
        mod->next = NULL;
        mod = OrgModFree (mod);
      } else {
        prev_mod = mod;
      }
    } else {
      prev_mod = mod;
    }
    mod = mod->next;
  }
}


static void RemoveSubSourceText (BioSourcePtr biop, ConvertFormPtr cfp)
{
  SubSourcePtr ssp, next_ssp, prev_ssp = NULL;

  if (biop == NULL || cfp == NULL) return;

  ssp = biop->subtype;
  while (ssp != NULL) {
    next_ssp = ssp->next;
    if (ssp->subtype == (cfp->subtype - 1000)) {
	    SearchAndExciseText (ssp->name, cfp);
      if (StringHasNoText (ssp->name) && !IsNonTextModifier (GetSubsourceQualName (ssp->subtype))) {
        if (prev_ssp == NULL) {
          biop->subtype = ssp->next;
        } else {
          prev_ssp->next = ssp->next;
        }
        ssp->next = NULL;
        ssp = SubSourceFree (ssp);
      } else {
        prev_ssp = ssp;
      }
    } else {
      prev_ssp = ssp;
    }
	  ssp = next_ssp;
  }
}

static void RemoveASourceText (BioSourcePtr biop, Pointer data)
{
  ConvertFormPtr cfp;
  OrgRefPtr      orp;
  CharPtr        tmp_str;

  cfp = (ConvertFormPtr) data;
  if (biop == NULL || cfp == NULL) return;
  switch (cfp->type) {
    case eFieldTypeBioSource :
      orp = biop->org;
      if (orp == NULL) return;
      switch (cfp->subtype) {
        case ORGREF_SCI_NAME_FIELD :
          tmp_str = StringSave (orp->taxname);
          SearchAndExciseText (tmp_str, cfp);
          if (StringCmp (tmp_str, orp->taxname) == 0) {
            /* no change, no need to remove taxref */
            tmp_str = MemFree (tmp_str);
          } else {
            SetTaxNameAndRemoveTaxRef (orp, tmp_str);
          }
          break;
        case ORGREF_COMMON_NAME_FIELD :
          SearchAndExciseText (orp->common, cfp);
          break;
        case ORGREF_LINEAGE_FIELD :
          if (orp->orgname != NULL) {
            SearchAndExciseText (orp->orgname->lineage, cfp);
          }
          break;
        case ORGREF_DIVISION_FIELD :
          if (orp->orgname != NULL) {
            SearchAndExciseText (orp->orgname->div, cfp);
          }
          break;
        default:
          break;
      }
      break;
    case eFieldTypeOrgModSubSource :
      if (cfp->subtype < 1000) {  /* Orgmod Note */
        RemoveOrgModText (biop, cfp);
      } else {  /* subsource note */
        RemoveSubSourceText (biop, cfp);
      }
      break;
    default:
      break;
  }
}


static void RemoveDescriptorTextCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ConvertFormPtr    cfp;
  CharPtr           title;
  ObjValNodePtr     ovp;

  cfp = (ConvertFormPtr) userdata;

  if (sdp == NULL || cfp == NULL) return;

  if (cfp->type == eFieldTypeDefline) {
    if (sdp->choice != Seq_descr_title) {
      return;
    }
  } else if (cfp->type == eFieldTypeCommentDescriptor) {
    if (sdp->choice != Seq_descr_comment) {
      return;
    }
  } else {
    return;
  }

  title = (CharPtr) sdp->data.ptrvalue;
  SearchAndExciseText (title, cfp);
  
  if (StringHasNoText (title) && (sdp->extended > 0)) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.deleteme = TRUE;
    cfp->isDirty = TRUE;
  }
}

  
static void DoRemoveText (ConvertFormPtr cfp)

{
  SeqEntryPtr     sep;
  FieldSubfieldPtr f;

  /* Get the current sequence and associated data */

  if (cfp == NULL || cfp->input_entityID == 0)
    return;

  sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
  if (sep == NULL)
    return;

  /* Get the feature and subfeature types */
  f = DialogToPointer (cfp->target_dlg);
  if (f == NULL || f->field < 0 || (f->subfield < 0 && f->subfield_list == NULL)) {
    f = FieldSubfieldFree (f);
    return;
  }

  cfp->type = f->field;
  cfp->subtype = f->subfield;

  f = FieldSubfieldFree (f);

  /* Hide the window and set the 'working' cursor */

  Hide (cfp->form);
  WatchCursor ();
  Update ();


  /* get information about text to remove */
  cfp->textportion = (TextPortionXPtr) DialogToPointer (cfp->textportion_dlg);
  
  if (GetValue (cfp->repeat_remove_grp) == 2) {
    cfp->repeat_remove = TRUE;
  } else {
    cfp->repeat_remove = FALSE;
  }

  /* Do actual work of removing text */
  cfp->isDirty = FALSE;
  if (cfp->type == eFieldTypeDefline || cfp->type == eFieldTypeCommentDescriptor) {
    VisitDescriptorsInSep (sep, cfp, RemoveDescriptorTextCallback);
  }
  else if (cfp->type == eFieldTypeBioSource || cfp->type == eFieldTypeOrgModSubSource) {
    VisitBioSourcesInSep (sep, cfp, RemoveASourceText);
  }
  else {
    VisitFeaturesInSep (sep, cfp, RemoveAFeatureText);
  }
  
  /* Clean up and exit */

  cfp->textportion = TextPortionXFree (cfp->textportion);

  if (cfp->isDirty) {
    DeleteMarkedObjects (cfp->input_entityID, 0, NULL);
  }

  ObjMgrSetDirtyFlag (cfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID, 0, 0);
  if (GetStatus (cfp->leaveDlgUp))
  {
  	Show (cfp->form);
  }
  else
  {
    Remove (cfp->form);  	
  }
  ArrowCursor ();
  Update ();
}

static void ConvertMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}


/* gather information specific ro remove inside text, then call common removetext function */
static void DoRemoveInsideText (ButtoN b)
{
  ConvertFormPtr cfp;

  cfp = (ConvertFormPtr) GetObjectExtra (b);
  if (cfp == NULL) return;

  if (GetValue (cfp->repeat_remove_grp) == 2) {
    cfp->repeat_remove = TRUE;
  } else {
    cfp->repeat_remove = FALSE;
  }

  DoRemoveText (cfp);
}


/* gather information specific ro remove inside text, then call common removetext function */
static void DoRemoveOutsideText (ButtoN b)
{
  ConvertFormPtr cfp;

  cfp = (ConvertFormPtr) GetObjectExtra (b);
  if (cfp == NULL) return;

  if (GetValue (cfp->ifNotFoundGroup) == DO_NOTHING)
    cfp->ifNotFound = DO_NOTHING;
  else
    cfp->ifNotFound = REMOVE_ALL_TEXT;


  DoRemoveText (cfp);
}


/*-------------------------------------------------------------------------*/
/*                                                                         */
/* SetRemoveTextAcceptButton () -- Enable/Disable the Accept button depending */
/*                              on the condition of other window objects.  */
/*                                                                         */
/*-------------------------------------------------------------------------*/

static void SetRemoveTextAcceptButton (Pointer data)

{
  ConvertFormPtr   cfp;
  FieldSubfieldPtr f;
  Boolean          ok_to_accept = FALSE;
  TextPortionXPtr   tp;

  cfp = (ConvertFormPtr) data;
  if (cfp == NULL) return;

  f = DialogToPointer (cfp->target_dlg);
  if (f != NULL && f->field >= 0 && (f->subfield >= 0 || f->subfield_list != NULL)) {
    if (f->field == eFieldTypeFeatureNote) {
      if (f->subfield_list != NULL) {
        ok_to_accept = TRUE;
      }
    } else { 
      ok_to_accept = TRUE;
    }
  }
  f = FieldSubfieldFree (f);
  if (ok_to_accept) {
    tp = (TextPortionXPtr) DialogToPointer (cfp->textportion_dlg);
    if (tp == NULL || (StringHasNoText (tp->start_text) && StringHasNoText (tp->end_text))) {
      ok_to_accept = FALSE;        
    }
    tp = TextPortionXFree (tp);
  }
  if (ok_to_accept) {
    SafeEnable (cfp->accept);
  } else {
    SafeDisable (cfp->accept);
  }
  f = FieldSubfieldFree (f);
}


static void ClearRemoveOutsideText (ButtoN b)
{
  ConvertFormPtr cfp;

  cfp = (ConvertFormPtr) GetObjectExtra (b);
  if (cfp == NULL) return;

  PointerToDialog (cfp->target_dlg, NULL);
  PointerToDialog (cfp->textportion_dlg, NULL);

  SetRemoveTextAcceptButton (cfp);
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* RemoveTextOutsideString ()                                          */
/*                                                                     */
/*---------------------------------------------------------------------*/

extern void RemoveTextOutsideString (IteM i)
{
  Boolean            allowed[eFieldType_Max];
  BaseFormPtr        bfp;
  GrouP              buttonGroup;
  ConvertFormPtr     cfp;
  GrouP              mainGroup;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  PrompT             p, p2;
  ButtoN             b;

  /* Get current sequence */

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

  /* Create a form for passing user data on to callbacks */

  cfp = ConvertFormNew ();
  if (cfp == NULL)
    return;

  /* Create a window for getting user input */

  w = FixedWindow (-50, -33, -10, -10, "Remove Text Outside String",
		   StdCloseWindowProc);
  SetObjectExtra (w, cfp, StdCleanupFormProc);
  cfp->form = (ForM) w;
  cfp->formmessage = ConvertMessageProc;

  /* Attach some basic data and callbacks to the data form */

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = bfp->input_entityID;
  cfp->input_itemID   = bfp->input_itemID;
  cfp->input_itemtype = bfp->input_itemtype;

  cfp->set_accept_proc = SetRemoveTextAcceptButton ;
  cfp->remove_inside = FALSE;

  /* Main object group for holding all the others */

  mainGroup = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (mainGroup, 10, 10);

  /* Set up Text group and headers */

  p = StaticPrompt (mainGroup, "Remove Text", 0, 0, programFont, 'l');

  cfp->textportion_dlg = TextPortionXDialogEx (mainGroup, FALSE, SetRemoveTextAcceptButton, cfp);

  /* The "If not found" group */

  cfp->ifNotFoundGroup = HiddenGroup (mainGroup, 0, 3, NULL);
  StaticPrompt (cfp->ifNotFoundGroup, "If field doesn't contain searched" 
		" for text:", 0, 0, programFont, 'l');

  RadioButton (cfp->ifNotFoundGroup, "     Do nothing to field");
  RadioButton (cfp->ifNotFoundGroup, "     Remove entire field text");
  SetValue (cfp->ifNotFoundGroup, DO_NOTHING);

  p2 = StaticPrompt (mainGroup, "Perform excision in", 0, popupMenuHeight, programFont, 'l');

  /* selecting the target field */
  MemSet (allowed, TRUE, sizeof (Boolean) * eFieldType_Max);
  allowed[eFieldTypeFeatureNote] = FALSE;
  allowed[eFieldTypePublication] = FALSE;
  cfp->target_dlg = CreateFieldSubfieldDlg (mainGroup, allowed, ChangeTargetFields, cfp);

  /* clear button */
  b = PushButton (mainGroup, "Clear", ClearRemoveOutsideText);
  SetObjectExtra (b, cfp, NULL);

  /* Accept and Cancel buttons */

  buttonGroup = HiddenGroup (mainGroup, 4, 0, NULL);
  cfp->accept = DefaultButton (buttonGroup, "Accept", DoRemoveOutsideText);
  SetObjectExtra (cfp->accept, cfp, NULL);
  Disable (cfp->accept);
  PushButton (buttonGroup, "Cancel", StdCancelButtonProc);
  cfp->leaveDlgUp = CheckBox (buttonGroup, "Leave Dialog Up", NULL);

  /* Line things up and display the window */

  AlignObjects (ALIGN_CENTER,
                (HANDLE) p,
                (HANDLE) cfp->textportion_dlg,
                (HANDLE) cfp->ifNotFoundGroup,
                (HANDLE) p2,
		            (HANDLE) cfp->target_dlg,
                (HANDLE) b,
		            (HANDLE) buttonGroup,
		            NULL);

  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}

extern void RemoveTextInsideString (IteM i)
{
  BaseFormPtr        bfp;
  GrouP              c;
  ConvertFormPtr     cfp;
  GrouP              h;
  PrompT             ppt, ppt2;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  Boolean            allowed[eFieldType_Max];

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  cfp = ConvertFormNew ();
  if (cfp == NULL) return;

  cfp->set_accept_proc = SetRemoveTextAcceptButton;
  cfp->remove_inside = TRUE;

  w = FixedWindow (-50, -33, -10, -10, "Remove Text Inside String", StdCloseWindowProc);
  SetObjectExtra (w, cfp, StdCleanupFormProc);
  cfp->form = (ForM) w;
  cfp->formmessage = ConvertMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = bfp->input_entityID;
  cfp->input_itemID = bfp->input_itemID;
  cfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt = StaticPrompt (h, "Remove Text", 0, popupMenuHeight, programFont, 'c');
  cfp->textportion_dlg = TextPortionXDialogEx (h, TRUE, SetRemoveTextAcceptButton, cfp);
  
  cfp->repeat_remove_grp = HiddenGroup (h, 2, 0, NULL);
  RadioButton (cfp->repeat_remove_grp, "Remove first instance in each string");
  RadioButton (cfp->repeat_remove_grp, "Remove all instances");
  SetValue (cfp->repeat_remove_grp, 1);

  ppt2 = StaticPrompt (h, "Perform excision in", 0, popupMenuHeight, programFont, 'c');
  
  MemSet (allowed, TRUE, sizeof (Boolean) * eFieldType_Max);
  allowed[eFieldTypeCommentDescriptor] = TRUE;
  allowed[eFieldTypeFeatureNote] = FALSE;
  allowed[eFieldTypePublication] = FALSE;
  cfp->target_dlg = CreateFieldSubfieldDlg (h, allowed, ChangeTargetFields, cfp);

  c = HiddenGroup (h, 4, 0, NULL);
  cfp->accept = DefaultButton (c, "Accept", DoRemoveInsideText);
  SetObjectExtra (cfp->accept, cfp, NULL);
  Disable (cfp->accept);
  PushButton (c, "Cancel", StdCancelButtonProc);
  cfp->leaveDlgUp = CheckBox (c, "Leave Dialog Up", NULL);


  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) cfp->textportion_dlg, (HANDLE) cfp->repeat_remove_grp, (HANDLE) ppt2, (HANDLE) cfp->target_dlg, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (cfp->textportion_dlg);
  Update ();
}

/* AAForCodon is extern in seqport.c */
/* NLM_EXTERN Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes); */

/*
static Boolean CorrectStartCodonCallback (GatherContextPtr gcp)

{
  Uint1           aa;
  Boolean         bad_base;
  CodeBreakPtr    cbp;
  Uint1           codon [3];
  CharPtr         codes;
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  Int2            i;
  Uint1           residue;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  SeqLocPtr       slp;
  SeqPntPtr       spntp;
  SeqPortPtr      spp;
  ValNodePtr      vnp;

  if (gcp == NULL) return TRUE;
  sep = (SeqEntryPtr) gcp->userdata;
  if (sep == NULL ) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return TRUE;

  if (crp->code_break != NULL) return TRUE;

	gc = NULL;
	if (crp->genetic_code != NULL)
	{
		vnp = (ValNodePtr)(crp->genetic_code->data.ptrvalue);
		while ((vnp != NULL) && (gcp == NULL))
		{
			switch (vnp->choice)
			{
			case 1:
				gc = GeneticCodeFind(0, (CharPtr)vnp->data.ptrvalue);
				break;
			case 2:
				gc = GeneticCodeFind(vnp->data.intvalue, NULL);
				break;
			case 3:
			case 6:
			case 4:
			case 5:
			case 7:
			case 8:
			default:
				break;
			}
			vnp = vnp->next;
		}
	}
	if (gc == NULL)
		gc = GeneticCodeFind(1, NULL);
	if (gc == NULL) return TRUE;

	codes = NULL;
	for (vnp = (ValNodePtr)gc->data.ptrvalue; vnp != NULL; vnp = vnp->next)
	{
		if (vnp->choice == 3)
			codes = (CharPtr)vnp->data.ptrvalue;
	}
	if (codes == NULL) return TRUE;

  if (crp->frame == 2 || crp->frame == 3) return TRUE;

  spp = SeqPortNewByLoc (sfp->location, Seq_code_ncbi4na);
  if (spp == NULL) return TRUE;
  bad_base = FALSE;
  for (i = 0; i < 3; i++) {
    residue = SeqPortGetResidue (spp);
    if (residue == SEQPORT_EOF)
      break;
    if (residue == INVALID_RESIDUE)
      bad_base = TRUE;
    codon[i] = residue;
  }
  SeqPortFree (spp);
  if (i != 3 || bad_base) return TRUE;
  aa = AAForCodon (codon, codes);

  spp = SeqPortNewByLoc (sfp->product, Seq_code_ncbieaa);
  if (spp == NULL) return TRUE;
  residue = SeqPortGetResidue (spp);
  SeqPortFree (spp);
  if (residue == SEQPORT_EOF || residue == INVALID_RESIDUE) return TRUE;

  if (residue != aa) {
    cbp = CodeBreakNew ();
    if (cbp != NULL) {
      spntp = SeqPntNew ();
      slp = ValNodeNew (NULL);
      slp->choice = SEQLOC_PNT;
      slp->data.ptrvalue = (Pointer) spntp;
      spntp->point = (Int4) 0;
      spntp->id = SeqIdStripLocus (SeqIdDup (SeqLocId (sfp->product)));
      cbp->loc = slp;
      slp = aaLoc_to_dnaLoc (sfp, cbp->loc);
      cbp->loc = SeqLocFree (cbp->loc);
      cbp->loc = slp;
      cbp->aa.value.intvalue = aa;
      cbp->aa.choice = 1;
      cbp->next = crp->code_break;
      crp->code_break = cbp;
    }
  }

  return TRUE;
}

static void CorrectStartCodon (SeqEntryPtr sep, Uint2 entityID)

{
  BioseqSetPtr  bssp;
  GatherScope   gs;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        CorrectStartCodon (sep, entityID);
      }
      return;
    }
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) sep, CorrectStartCodonCallback, &gs);
}

extern void CorrectCDSStartCodon (IteM i);
extern void CorrectCDSStartCodon (IteM i)

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
  CorrectStartCodon (sep, bfp->input_entityID);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}
*/


typedef struct applyformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  Int2           errcount;
  ValNodePtr     ambigList;
  Boolean        noLeft;
  Boolean        noRight;
  ButtoN         applyToParts;
  ButtoN         partial5;
  ButtoN         partial3;
  TexT           onlyThisPart;
  GrouP          all_or_some_grp;
  TexT           accession_list_txt;
  
  DialoG         feature_details_dlg;
  
  BatchApplyFeatureDetailsPtr feature_details_data;
  
  GrouP          strand_group;
  GrouP          use_whole_interval;
  TexT           left_end;
  TexT           right_end;
  ButtoN         add_to_seq_with_like_feature;
  ButtoN         also_add_mRNA_btn;
  ButtoN         accept;
  ButtoN         leaveDlgUp;
  
  GetSamplePtr   gsp;
  ExistingTextPtr etp;
} ApplyFormData, PNTR ApplyFormPtr;

typedef struct alreadyhas {
  Boolean        rsult;
  Uint1          featchoice;
  Uint1          descchoice;
  Uint1          rnatype;
} AlreadyHas, PNTR AlreadyHasPtr;

static Boolean SeeIfAlreadyHasGatherFunc (GatherContextPtr gcp)

{
  AlreadyHasPtr  ahp;
  RnaRefPtr      rrp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;

  ahp = (AlreadyHasPtr) gcp->userdata;
  if (ahp == NULL ) return TRUE;

  if (gcp->thistype == OBJ_SEQFEAT && ahp->featchoice != 0) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == ahp->featchoice && sfp->data.value.ptrvalue != NULL) {
      if (sfp->data.choice == SEQFEAT_RNA) {
        rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
        if (rrp->type != ahp->rnatype) return TRUE;
      }
      ahp->rsult = TRUE;
      return FALSE;
    }
  } else if (gcp->thistype == OBJ_SEQDESC && ahp->descchoice != 0) {
    sdp = (ValNodePtr) gcp->thisitem;
    if (sdp != NULL && sdp->choice == ahp->descchoice && sdp->data.ptrvalue != NULL) {
      ahp->rsult = TRUE;
      return FALSE;
    }
  }
  return TRUE;
}

static Boolean AlreadyHasFeatOrDesc (SeqEntryPtr sep, Uint1 featchoice, Uint1 descchoice, Uint1 rnatype)

{
  AlreadyHas   ah;
  BioseqPtr    bsp;
  GatherScope  gs;
  SeqEntryPtr  nsep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  ah.rsult = FALSE;
  ah.featchoice = featchoice;
  ah.descchoice = descchoice;
  ah.rnatype = rnatype;
  if (sep == NULL) return FALSE;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  if (descchoice != 0) {
    nsep = FindNucSeqEntry (sep);
    if (nsep != NULL && IS_Bioseq (nsep)) {
      bsp = (BioseqPtr) nsep->data.ptrvalue;
      if (bsp != NULL) {
        slp = ValNodeNew (NULL);
        slp->choice = SEQLOC_WHOLE;
        sip = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
        slp->data.ptrvalue = sip;
        gs.target = slp;
      }
    }
  }
  GatherSeqEntry (sep, (Pointer) (&ah), SeeIfAlreadyHasGatherFunc, &gs);
  gs.target = SeqLocFree (gs.target);
  return ah.rsult;
}

static void setSeqFeatStrand(SeqFeatPtr sfp, Uint1 strand)
{
  SeqIntPtr          sqip;

  if (sfp == NULL) return;
  if (sfp->location == NULL) return;
  if (sfp->location->choice != SEQLOC_INT) return;
  sqip = (SeqIntPtr) sfp->location->data.ptrvalue;
  if (sqip != NULL)
  {
    sqip->strand = Seq_strand_minus;
  }
}

/* must be called after partials and strand are set*/
static void AdjustSeqLocForApply (SeqFeatPtr sfp, ApplyFormPtr afp)
{
  BioseqPtr bsp;
  SeqLocPtr slp;
  SeqIntPtr sip;
  CharPtr   num_str;
  Int4      from;
  Int4      to;
  Int4      tmp;
  Boolean   partial3, partial5;
  Boolean   is_caret = FALSE;
  Uint1     strand;
  
  if (sfp == NULL || afp == NULL) return;
  
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  
  num_str = SaveStringFromText (afp->left_end);
  if (num_str == NULL) return;
  from = atoi (num_str);
  tmp = StringLen (num_str);
  if (tmp > 1 && num_str[tmp - 1] == '^')
  {
    is_caret = TRUE;
  } 
  num_str = MemFree (num_str);
  num_str = SaveStringFromText (afp->right_end);
  if (num_str == NULL) return;
  if (num_str[0] == '^') {
    is_caret = TRUE;
    to = atoi (num_str + 1);
  } else {
    to = atoi (num_str);
  }
  num_str = MemFree (num_str);
  if (from > to)
  {
  	tmp = from;
  	from = to;
  	to = tmp;
  }
  if (from < 1)
  {
  	from = 1;
  }
  if (to < 1)
  {
    if (to == 0 && is_caret) {
      to = from + 1;
    } else {
  	  to = 1;
  	}
  }
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  strand = SeqLocStrand (sfp->location);
  if (to > bsp->length)
  {
  	to = bsp->length;
  	partial3 = TRUE;
  }

  if (is_caret && to != from + 1) {
    is_caret = FALSE;
  }

  slp = NULL;
  if (is_caret) {
    AddSeqLocPoint (&slp, SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0))),
                    from, FALSE, TRUE, strand);
  } else {
    slp = ValNodeNew (NULL);
    if (slp != NULL) {
      sip = SeqIntNew ();
      if (sip != NULL) {
        sip->from = from - 1;
        sip->to = to - 1;
        sip->strand = strand;
        sip->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
        slp->choice = SEQLOC_INT;
        slp->data.ptrvalue = (Pointer) sip;
        SetSeqLocPartial (slp, partial5, partial3);
      }
    }
  }
  
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  
  if (partial5 || partial3)
  {
  	sfp->partial = TRUE;
  }
  
}

static void SetApplyFeatureLocation (SeqFeatPtr sfp, ApplyFormPtr afp)
{
  if (sfp == NULL || afp == NULL)
  {
    return;
  }
  
  if (GetValue (afp->strand_group) == 2) 
  {
    /* reverse strand direction - strand direction is plus by default */
    setSeqFeatStrand (sfp, Seq_strand_minus);
  }

  SetSeqLocPartial (sfp->location, afp->noLeft, afp->noRight);

  if (afp->use_whole_interval != NULL
      && GetValue (afp->use_whole_interval) != 1)
  {
    /* adjust location to match coordinates from user */
    AdjustSeqLocForApply (sfp, afp);
  }      
    
  sfp->partial = (afp->noLeft || afp->noRight);
}

static void AddGeneXrefToFeat (SeqFeatPtr sfp, CharPtr str)
{
  SeqFeatXrefPtr    xref;
  GeneRefPtr        grp;
  
  if (sfp == NULL || StringHasNoText (str)) return;
  
  /* add gene xref to feature */
  xref = SeqFeatXrefNew ();
  if (xref != NULL)
  {
    grp = CreateNewGeneRef (str, NULL, NULL, FALSE);
    if (grp != NULL) 
    {
      xref->data.choice = SEQFEAT_GENE;
      xref->data.value.ptrvalue = grp;
      xref->next = sfp->xref;
      sfp->xref = xref;
    }
  }
}

static SeqFeatPtr ApplyGene (CharPtr str, ApplyFormPtr afp, SeqEntryPtr gene_sep, SeqFeatPtr sfp)
{
  GeneRefPtr        grp;
  SeqFeatPtr        gene_sfp;
  SeqFeatXrefPtr    xref;
  SeqMgrFeatContext fcontext;
  BioseqPtr         bsp = NULL;
  SeqFeatPtr        other_feat;
  SeqFeatPtr        overlap_gene;
  Boolean           added_xrefs = FALSE;
  SeqFeatPtr        misc_feat = NULL;
  SeqLocPtr         overlap_loc;
  CharPtr           gene_desc = NULL;

  if (afp == NULL || gene_sep == NULL 
	  || (StringHasNoText (str) 
	      && (afp->feature_details_data == NULL 
		      || StringHasNoText (afp->feature_details_data->geneDesc))))
  {
    return NULL;
  }

  if (afp != NULL && afp->feature_details_data != NULL)
  {
    gene_desc = afp->feature_details_data->geneDesc;
  }

  /* we need a location to use when we're checking for feature-stealing genes */
  if (sfp != NULL)
  {
    overlap_loc = sfp->location;
  }
  else
  {
    misc_feat = CreateNewFeature (gene_sep, NULL, SEQFEAT_COMMENT, NULL);
    if (NULL == misc_feat)
    return NULL;
    
    SetApplyFeatureLocation (misc_feat, afp);

    overlap_loc = misc_feat->location;
  }
  
  /* first, add gene xrefs to all features on bioseq that are contained in the location */
  /* maintain list of features that had xrefs before, should not remove them later */
  if (IS_Bioseq (gene_sep))
  {
    bsp = (BioseqPtr) gene_sep->data.ptrvalue;
  }
  else if (sfp != NULL)
  {
    bsp = BioseqFindFromSeqLoc (sfp->location);
  }
  if (bsp != NULL)
  {
    other_feat = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (other_feat != NULL)
    {
      if (other_feat != sfp && other_feat->data.choice != SEQFEAT_GENE)
      {
        for (xref = other_feat->xref;
             xref != NULL && xref->data.choice != SEQFEAT_GENE;
             xref = xref->next)
        {}
        if (xref == NULL
            && SeqLocCompare (other_feat->location, overlap_loc) == SLC_A_EQ_B)
        {
          overlap_gene = SeqMgrGetOverlappingGene (other_feat->location, &fcontext);
          if (overlap_gene != NULL)
          {
            AddGeneXrefToFeat (other_feat, fcontext.label);
            added_xrefs = TRUE;
          }
        }
      }
      other_feat = SeqMgrGetNextFeature (bsp, other_feat, 0, 0, &fcontext);
    }   
  }
  
  if (misc_feat != NULL)
  {
    misc_feat->idx.deleteme = TRUE;
    DeleteMarkedObjects (0, OBJ_SEQENTRY, gene_sep);
  }
  
  grp = CreateNewGeneRef (str, NULL, gene_desc, FALSE);
  if (NULL == grp)
    return NULL;

  gene_sfp = CreateNewFeature (gene_sep, NULL, SEQFEAT_GENE, NULL);
  if (NULL == gene_sfp)
    return NULL;

  gene_sfp->data.value.ptrvalue = (Pointer) grp;
  
  SetApplyFeatureLocation (gene_sfp, afp);

  if (added_xrefs && sfp != NULL)
  {
    /* add gene xref to feature */
    AddGeneXrefToFeat (sfp, str);
  }
  
  return gene_sfp;
}


typedef struct adjustfeatforgapdialog {
  DIALOG_MESSAGE_BLOCK
  DialoG feature_select;
  ButtoN unknown_gaps;
  ButtoN known_gaps;
  GrouP  partial_grp;
  ButtoN trim_ends;
  ButtoN split_internal;
  ButtoN split_in_intron;

  Nlm_ChangeNotifyProc     change_notify;
  Pointer                  change_userdata;
} AdjustFeatForGapDialogData, PNTR AdjustFeatForGapDialogPtr;


static void AdjustFeaturesForGapToDialog (DialoG d, Pointer udata)
{
  AdjustFeatForGapDialogPtr dlg;
  AdjustFeatForGapPtr       data;

  dlg = (AdjustFeatForGapDialogPtr) GetObjectExtra (d);
  data = (AdjustFeatForGapPtr) udata;

  if (dlg == NULL) return;
 
  if (data == NULL) {
    PointerToDialog (dlg->feature_select, NULL);
    SetStatus (dlg->unknown_gaps, FALSE);
    SetStatus (dlg->known_gaps, FALSE);
    SetStatus (dlg->trim_ends, FALSE);
    SetStatus (dlg->split_internal, FALSE);
    SetStatus (dlg->split_in_intron, FALSE);
    Disable (dlg->split_in_intron);
    SetValue (dlg->partial_grp, 2);
  } else {
    PointerToDialog (dlg->feature_select, data->feature_list);
    SetStatus (dlg->unknown_gaps, data->options & eAdjustFeatForGap_unknown_gaps);
    SetStatus (dlg->known_gaps, data->options & eAdjustFeatForGap_known_gaps);
    SetStatus (dlg->split_internal, data->options & eAdjustFeatForGap_split_internal);
    SetStatus (dlg->split_in_intron, data->options & eAdjustFeatForGap_split_in_intron);
    if (!(data->options & eAdjustFeatForGap_split_internal)) {
      Disable (dlg->split_in_intron);
    }
    SetStatus (dlg->trim_ends, data->options & eAdjustFeatForGap_trim_ends);

    if ((data->options & eAdjustFeatForGap_make_partial) && (data->options & eAdjustFeatForGap_partial_for_pseudo)) {
      SetValue (dlg->partial_grp, 1);
    } else if ((data->options & eAdjustFeatForGap_make_partial) && !(data->options & eAdjustFeatForGap_partial_for_pseudo)) {
      SetValue (dlg->partial_grp, 2);
    } else {
      SetValue (dlg->partial_grp, 3);
    }
  }    
}


static Pointer AdjustFeaturesForGapToPointer (DialoG d)
{
  AdjustFeatForGapDialogPtr dlg;
  AdjustFeatForGapPtr       data;
  Int4                      grp_val;

  dlg = (AdjustFeatForGapDialogPtr) GetObjectExtra (d);

  if (dlg == NULL) return NULL;

  data = (AdjustFeatForGapPtr) MemNew (sizeof (AdjustFeatForGapData));

  /* get list of feature types to ask on */
  data->feature_list = (ValNodePtr) DialogToPointer (dlg->feature_select);
  data->features_in_gap = NULL;

  if (GetStatus (dlg->unknown_gaps)) {
    data->options |= eAdjustFeatForGap_unknown_gaps;
  }
  if (GetStatus (dlg->known_gaps)) {
    data->options |= eAdjustFeatForGap_known_gaps;
  }

  grp_val = GetValue (dlg->partial_grp);
  switch (grp_val) {
    case 1:
      data->options |= eAdjustFeatForGap_make_partial | eAdjustFeatForGap_partial_for_pseudo;
      break;
    case 2:
      data->options |= eAdjustFeatForGap_make_partial;
      break;
    case 3:
      /* both false */
      break;
  }
  
  if (GetStatus (dlg->split_internal)) {
    data->options |= eAdjustFeatForGap_split_internal;
    if (GetStatus (dlg->split_in_intron)) {
      data->options |= eAdjustFeatForGap_split_in_intron;
    }
  }
  if (GetStatus (dlg->trim_ends)) {
    data->options |= eAdjustFeatForGap_trim_ends;
  }

  return (Pointer) data;
}


static ValNodePtr TestAdjustForGapsDialog (DialoG d)
{
  AdjustFeatForGapDialogPtr dlg;
  AdjustFeatForGapPtr       data;
  ValNodePtr                err_list = NULL;
  
  dlg = (AdjustFeatForGapDialogPtr) GetObjectExtra (d);
  if (dlg == NULL)
  {
    ValNodeAddPointer (&err_list, 0, "No dialog");
    return err_list;
  }

  data = (AdjustFeatForGapPtr) DialogToPointer (d);
  if (data == NULL) 
  {
    ValNodeAddPointer (&err_list, 0, "No data");
    return err_list;
  }
    
  if (data->feature_list == NULL) 
  {
    ValNodeAddPointer (&err_list, 0, "No features");
  }

  if (!(data->options & eAdjustFeatForGap_unknown_gaps) && !(data->options & eAdjustFeatForGap_known_gaps))
  {
    ValNodeAddPointer (&err_list, 0, "No gaps");
  }

  if (!(data->options & eAdjustFeatForGap_split_internal) && !(data->options & eAdjustFeatForGap_trim_ends))
  {
    ValNodeAddPointer (&err_list, 0, "No action");
  }

  data = AdjustFeatForGapFree (data);

  return err_list;    
}


static void AdjustFeaturesForGapsChangeNotifyButton (ButtoN b)
{
  AdjustFeatForGapDialogPtr dlg;

  dlg = (AdjustFeatForGapDialogPtr) GetObjectExtra (b);
  if (GetStatus (dlg->split_internal)) {
    Enable (dlg->split_in_intron);
  } else {
    Disable (dlg->split_in_intron);
  }
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static void AdjustFeaturesForGapsChangeNotifyGroup (GrouP g)
{
  AdjustFeatForGapDialogPtr dlg;

  dlg = (AdjustFeatForGapDialogPtr) GetObjectExtra (g);
  if (dlg != NULL && dlg->change_notify != NULL) {
    (dlg->change_notify) (dlg->change_userdata);
  }
}


static DialoG AdjustFeaturesForGapDialog (GrouP h, Uint2 entityID, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata)
{
  AdjustFeatForGapDialogPtr dlg;
  GrouP                     p, g, g2;
  SeqEntryPtr               sep;

  dlg = (AdjustFeatForGapDialogPtr) MemNew (sizeof (AdjustFeatForGapDialogData));
  if (dlg == NULL) return NULL;

  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = AdjustFeaturesForGapToDialog;
  dlg->fromdialog = AdjustFeaturesForGapToPointer;
  dlg->dialogmessage = NULL;
  dlg->testdialog = TestAdjustForGapsDialog;
  dlg->change_notify = change_notify;
  dlg->change_userdata = change_userdata;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  dlg->feature_select = FeatureSelectionDialogEx (p, TRUE, sep,
                                                  change_notify, 
                                                  change_userdata);

  g = HiddenGroup (p, 2, 0, NULL);
  dlg->unknown_gaps = CheckBox (g, "Unknown length gaps", AdjustFeaturesForGapsChangeNotifyButton);
  SetObjectExtra (dlg->unknown_gaps, dlg, NULL);
  dlg->known_gaps = CheckBox (g, "Known length gaps", AdjustFeaturesForGapsChangeNotifyButton);
  SetObjectExtra (dlg->known_gaps, dlg, NULL);

  dlg->partial_grp = NormalGroup (p, 3, 0, "Make truncated ends partial", programFont, AdjustFeaturesForGapsChangeNotifyGroup);
  SetObjectExtra (dlg->partial_grp, dlg, NULL);
  RadioButton (dlg->partial_grp, "Always");
  RadioButton (dlg->partial_grp, "Unless pseudo");
  RadioButton (dlg->partial_grp, "Never");
  SetValue (dlg->partial_grp, 2);

  g2 = HiddenGroup (p, 3, 0, NULL);
  dlg->trim_ends = CheckBox (g2, "Trim ends in gaps", AdjustFeaturesForGapsChangeNotifyButton);
  SetObjectExtra (dlg->trim_ends, dlg, NULL);
  dlg->split_internal = CheckBox (g2, "Split for internal gaps", AdjustFeaturesForGapsChangeNotifyButton);
  SetObjectExtra (dlg->split_internal, dlg, NULL);
  dlg->split_in_intron = CheckBox (g2, "(Even when gaps are in introns)", AdjustFeaturesForGapsChangeNotifyButton);
  SetObjectExtra (dlg->split_in_intron, dlg, NULL);
  Disable (dlg->split_in_intron);

  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->feature_select, (HANDLE) g, (HANDLE) dlg->partial_grp, (HANDLE) g2, NULL);
  return (DialoG) p;
}

typedef struct adjustfeatforgapform {
  FORM_MESSAGE_BLOCK

  DialoG dlg;
  DialoG clickable_list;
  ButtoN accept_btn;
  
  ValNodePtr feat_list;
} AdjustFeatforGapFormData, PNTR AdjustFeatForGapFormPtr;


typedef struct featforadjust {
  AdjustFeatForGapPtr afgp;
  ValNodePtr          features_to_adjust;
  ValNodePtr          features_contain_gap;
  ValNodePtr          features_in_gap;  
} FeatForAdjustData, PNTR FeatForAdjustPtr;


static void FindFeaturesToBeAdjustedForGapsCallback (SeqFeatPtr sfp, Pointer data)
{
  FeatForAdjustPtr faap;
  BioseqPtr   gapped_bioseq;
  Boolean     terminal_gaps = FALSE;
  Boolean     entirely_in_gap = FALSE;
  Boolean     internal_gaps = FALSE;

  if (sfp == NULL || data == NULL) return;

  faap = (FeatForAdjustPtr) data;
  if (faap->afgp == NULL) return;

  if (!FeatureOkForFeatureList(sfp, faap->afgp->feature_list)) return;

  gapped_bioseq = BioseqFind (SeqLocId (sfp->location));

  LocationContainsGaps (sfp->location, gapped_bioseq, faap->afgp->options, &terminal_gaps, &internal_gaps, &entirely_in_gap);
  if (entirely_in_gap)
  {
    ValNodeAddPointer (&(faap->features_in_gap), OBJ_SEQFEAT, sfp);
  }
  
  if (internal_gaps)
  {
    ValNodeAddPointer (&(faap->features_contain_gap), OBJ_SEQFEAT, sfp);
  }

  if (((faap->afgp->options & eAdjustFeatForGap_split_internal) && internal_gaps) 
      || ((faap->afgp->options & eAdjustFeatForGap_trim_ends) && terminal_gaps)) 
  {
    ValNodeAddPointer (&(faap->features_to_adjust), OBJ_SEQFEAT, sfp);
  }
}


static void AdjustFeaturesForGapsChangeNotify (Pointer data)
{
  AdjustFeatForGapFormPtr dlg;
  FeatForAdjustData       ffad;
  ValNodePtr              err_list;
  SeqEntryPtr             sep;
  ClickableItemPtr        cip;

  dlg = (AdjustFeatForGapFormPtr) data;

  if (dlg == NULL) return;

  err_list = TestDialog (dlg->dlg);
  if (err_list == NULL) 
  {
    Enable (dlg->accept_btn);
  } 
  else
  {
    err_list = ValNodeFree (err_list);
    Disable (dlg->accept_btn);
  }

  /* clear clickable list display */
  PointerToDialog (dlg->clickable_list, NULL);
  /* re-create list of features */
  dlg->feat_list = FreeClickableList (dlg->feat_list);

  ffad.afgp = DialogToPointer (dlg->dlg);
  ffad.features_contain_gap = NULL;
  ffad.features_in_gap = NULL;
  ffad.features_to_adjust = NULL;

  if (ffad.afgp->feature_list != NULL)
  {
    sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
    VisitFeaturesInSep (sep, &(ffad), FindFeaturesToBeAdjustedForGapsCallback);
  }

  if (ffad.features_to_adjust != NULL) {
    cip = NewClickableItem (0, "%d features will be adjusted", ffad.features_to_adjust);
    ValNodeAddPointer (&(dlg->feat_list), 0, cip);
  }
  if (ffad.features_in_gap != NULL) {
    cip = NewClickableItem (0, "%d features are completely contained in gaps and will be deleted", ffad.features_in_gap);
    ValNodeAddPointer (&(dlg->feat_list), 0, cip);
  }
  if (ffad.features_contain_gap != NULL) {
    cip = NewClickableItem (0, "%d features contain internal gaps", ffad.features_contain_gap);
    ValNodeAddPointer (&(dlg->feat_list), 0, cip);
  }

  /* pass list of features to clickable list display */
  PointerToDialog (dlg->clickable_list, dlg->feat_list);

}


static void DoAdjustFeaturesForGaps (ButtoN b)
{
  AdjustFeatForGapFormPtr dlg;
  AdjustFeatForGapPtr     data;
  SeqEntryPtr             sep;

  dlg = (AdjustFeatForGapFormPtr) GetObjectExtra (b);

  if (dlg == NULL) return;

  data = (AdjustFeatForGapPtr) DialogToPointer (dlg->dlg);

  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  VisitFeaturesInSep (sep, data, AdjustFeatureForGapsCallback);

  /* remove features entirely in the gap */
  MarkFeaturesInGapsForDeletion (data);
  DeleteMarkedObjects (dlg->input_entityID, 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);

  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);

  data = AdjustFeatForGapFree(data);
  
  Remove (dlg->form);  
}


static void MakeGapFeatureReport (ButtoN b)
{
  FILE                    *fp;
  Char                    path [PATH_MAX];
  AdjustFeatForGapFormPtr dlg;
  Int4                    num_disc = 0;
  Boolean                 show_all = FALSE;


  dlg = (AdjustFeatForGapFormPtr) GetObjectExtra (b);

  if (dlg == NULL) return;
  
  num_disc = CountChosenDiscrepancies (dlg->feat_list, FALSE);

  if (num_disc == 0) 
  {
    if (ANS_CANCEL == Message (MSG_OKC, "No items selected!  Export all?"))
    {
      return;
    }
    else
    {
      show_all = TRUE;
    }
  }
  
  path [0] = '\0';
  if (GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      FileClose (fp);
    } else {
      FileCreate (path, "TEXT", "ttxt");
    }
#endif
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      WriteClickableListReport (fp, dlg->feat_list, show_all, FALSE);
      FileClose (fp);
      return;
    }
  }
}


extern void AdjustFeaturesForGaps (IteM i)
{
  BaseFormPtr             bfp;
  AdjustFeatForGapFormPtr dlg;
  WindoW                  w;
  GrouP                   h, c;
  ButtoN                  b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  dlg = (AdjustFeatForGapFormPtr) MemNew (sizeof (AdjustFeatforGapFormData));
  if (dlg == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Adjust Features for Gaps", StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;
  dlg->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  dlg->clickable_list =  CreateClickableListDialogEx (h, "", "Features", NULL, NULL,
                                                      ScrollToDiscrepancyItem, EditDiscrepancyItem, NULL,
                                                      GetDiscrepancyItemText,
                                                      stdCharWidth * 15,
                                                      stdCharWidth * 30,
                                                      TRUE, TRUE);

  dlg->dlg = AdjustFeaturesForGapDialog (h, bfp->input_entityID, AdjustFeaturesForGapsChangeNotify, dlg);

  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  dlg->accept_btn = PushButton (c, "Accept", DoAdjustFeaturesForGaps);
  SetObjectExtra (dlg->accept_btn, dlg, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  b = PushButton (c, "Make Report", MakeGapFeatureReport);
  SetObjectExtra (b, dlg, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) dlg->clickable_list,
                              (HANDLE) dlg->dlg,
                              (HANDLE) c,
                              NULL);
  RealizeWindow (w);
  AdjustFeaturesForGapsChangeNotify (dlg);
  Show (w);
  Select (w);
}


static ValNodePtr GapLocationsFromNs (BioseqPtr bsp, Int4Ptr gap_sizes)

{
  CharPtr     bases, txt;
  Char        ch;
  Int4        len;
  ValNodePtr  seq_ext;
  Boolean     unknown_greater_than_or_equal = FALSE;
  Boolean     known_greater_than_or_equal = FALSE;
  Int4        unknown_gap_size = 0;
  Int4        known_gap_size = 0;
  Int4        gap_len;
  Boolean     make_unknown_size;
  GapLocInfoPtr glip;
  ValNodePtr  result_list = NULL;
  Int4        start_pos = 0;
  

  if (bsp == NULL || bsp->repr != Seq_repr_raw || ISA_aa (bsp->mol)) 
  {
    return NULL;
  }
  
  if (gap_sizes == NULL)
  {
    known_greater_than_or_equal = TRUE;
  }
  else
  {
    unknown_gap_size = gap_sizes[0];
    known_gap_size = gap_sizes[1];
    if (unknown_gap_size < 0)
    {
      unknown_greater_than_or_equal = TRUE;
      unknown_gap_size = 0 - unknown_gap_size;
    }
    if (known_gap_size < 0)
    {
      known_greater_than_or_equal = TRUE;
      known_gap_size = 0 - known_gap_size;
    }
  }

  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) return NULL;

  for (txt = bases, ch = *txt; ch != '\0'; txt++, ch = *txt) {
    if (ch == 'N') break;
  }
  if (ch != 'N') {
    MemFree (bases);
    return NULL;
  }

  seq_ext = NULL;
  len = 0;

  txt = bases;
  ch = *txt;

  gap_len = 0;
  while (ch != '\0') {

    if (ch == 'N') {
      gap_len = StringSpn (txt, "N");
      if (gap_len == unknown_gap_size
          || (gap_len > unknown_gap_size && unknown_greater_than_or_equal)
          || gap_len == known_gap_size
          || (gap_len > known_gap_size && known_greater_than_or_equal))
      {
        make_unknown_size = FALSE;
        if (gap_len == 0)
        {
          make_unknown_size = FALSE;
        }
        else if (gap_len == unknown_gap_size)
        {
          make_unknown_size = TRUE;
        }
        else if (gap_len == known_gap_size)
        {
          make_unknown_size = FALSE;
        }
        else if (gap_len > unknown_gap_size && unknown_greater_than_or_equal)
        {
          if (!known_greater_than_or_equal)
          {
          	make_unknown_size = TRUE;
          }
          else if (unknown_gap_size > known_gap_size)
          {
          	make_unknown_size = TRUE;
          }
          else if (gap_len < known_gap_size)
          {
          	make_unknown_size = TRUE;
          }
        }
        
        /* Add Location to List */
        glip = (GapLocInfoPtr) MemNew (sizeof (GapLocInfoData));
        if (glip != NULL)
        {
          glip->start_pos = start_pos;
          glip->is_known = ! make_unknown_size;
          glip->length = gap_len;
          glip->replace = TRUE;
          ValNodeAddPointer (&result_list, 0, glip);
        }
      }
      txt += gap_len;
      start_pos += gap_len;
      ch = *txt;
      gap_len = 0;
    } else {
      gap_len = StringCSpn (txt, "N");
      txt+=gap_len;
      start_pos += gap_len;
      ch = *txt;
      gap_len = 0;
    }
  }
  
  MemFree (bases);
  return result_list;
}

static SeqLocPtr 
RemoveGapLocationsFromSeqLoc 
(SeqLocPtr  slp, 
 ValNodePtr gap_locs, 
 BioseqPtr  bsp,
 BoolPtr    loc_changed)
{
  ValNodePtr    gap_vnp;
  GapLocInfoPtr glip;
  SeqLocPtr   loc_list = NULL, prev_loc = NULL;
  SeqIdPtr    sip, before_sip, after_sip;
  Boolean     changed, partial5, partial3;
  SeqLocPtr   before = NULL, after = NULL;

  if (loc_changed != NULL)
  {
    *loc_changed = FALSE;
  }
  
  if (slp == NULL || gap_locs == NULL || bsp == NULL)
  {
    return slp;
  }
  
  sip = SeqLocId (slp);
  if (sip == NULL)
  {
    return slp;
  }
  
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  
  before = SeqLocCopy (slp);
  loc_list = before;
  
  for (gap_vnp = gap_locs;
       gap_vnp != NULL;
       gap_vnp = gap_vnp->next)
  {
    glip = (GapLocInfoPtr) gap_vnp->data.ptrvalue;
    if (glip == NULL || glip->is_known)
    {
      continue;
    }
    if (GapInLocation (glip->start_pos, glip->length, before))
    {
      if (loc_changed != NULL)
      {
        *loc_changed = TRUE;
      }
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
                              0, glip->start_pos + glip->length - 1,
                              FALSE, &changed, &partial5, &partial3);
                              
      /* in the "before" location, we free everything after the
       * beginning of the gap.
       */
      before = SeqLocDeleteEx (before, before_sip, 
                               glip->start_pos, bsp->length,
                               FALSE, &changed, &partial5, &partial3);
        
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

  slp = SeqLocFree (slp);
  return loc_list;  
}


static void AdjustCodingRegionLocationsForGapLocations (BioseqPtr bsp, ValNodePtr gap_list)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  BioseqPtr         protbsp = NULL, new_protbsp = NULL;
  SeqFeatPtr        new_sfp;
  CdRegionPtr       crp;
  Boolean           partial5, partial3;
  Uint2             entityID;
  Boolean           loc_changed;
  
  if (bsp == NULL || gap_list == NULL)
  {
    return;
  }
  
  entityID = bsp->idx.entityID;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext))
  {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    loc_changed = FALSE;
    sfp->location = RemoveGapLocationsFromSeqLoc (sfp->location, gap_list, bsp, &loc_changed);
    if (!loc_changed) {
        continue;
    }
    if (loc_changed)
    {
      AddCDSGapComment (sfp);
    }
   
    while (sfp->location->next != NULL)
    {
      new_sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_CDREGION, NULL);
		  if (new_sfp != NULL)
		  {
		    /* create copy of coding region data */
		    crp = (CdRegionPtr) AsnIoMemCopy ((CdRegionPtr) sfp->data.value.ptrvalue,
		                                       (AsnReadFunc) CdRegionAsnRead,
		                                       (AsnWriteFunc) CdRegionAsnWrite);
        new_sfp->data.value.ptrvalue = crp;
        new_sfp->location = sfp->location->next;
        new_sfp->comment = StringSave (sfp->comment);

        protbsp = BioseqFindFromSeqLoc (sfp->product);        
        if (protbsp != NULL)
        {
          new_protbsp = AddProteinSequenceCopy (protbsp, bsp, new_sfp, entityID);                                                  
        }
        
        /* still to do:
         * adjust frame to match prior protein
         */
      
        sfp->location->next = NULL;
        
        /* fix partials */
        SetSeqLocPartial (sfp->location, partial5, TRUE);
        sfp->partial = TRUE;
        
        /* adjust frame */
        AdjustFrame (sfp, protbsp);

        /* retranslate coding region */
        SeqEdTranslateOneCDS (sfp, bsp, entityID, Sequin_GlobalAlign2Seq);
        
        /* set partials on product */
        if (protbsp == NULL)
        {
          protbsp = BioseqFindFromSeqLoc (sfp->product);
        }
        SetProductSequencePartials (protbsp, partial5, partial3);

        partial5 = TRUE;
        sfp = new_sfp;
        protbsp = new_protbsp;
		  }
		  else
		  {
		    return; /* bail */
		  }
    }
    /* fix partials for last feature */
    SetSeqLocPartial (sfp->location, partial5, partial3);
    sfp->partial = partial5 || partial3;
    
    /* adjust frame */
    AdjustFrame (sfp, protbsp);

    /* retranslate coding region */
    SeqEdTranslateOneCDS (sfp, bsp, entityID, Sequin_GlobalAlign2Seq);
    /* set partials on product */
    if (protbsp == NULL)
    {
      protbsp = BioseqFindFromSeqLoc (sfp->product);
    }
    SetProductSequencePartials (protbsp, partial5, partial3);
  }
}

extern void 
PrepareCodingRegionLocationsForDeltaConversionCallback
(BioseqPtr bsp, Pointer userdata)
{
  Int4Ptr gap_sizes;
  ValNodePtr gap_locations;
  
  if (bsp == NULL)
  {
    return;
  }
  
  gap_sizes = (Int4Ptr) userdata;
  
  gap_locations = GapLocationsFromNs (bsp, gap_sizes);
  
  AdjustCodingRegionLocationsForGapLocations (bsp, gap_locations);
  
  gap_locations = ValNodeFreeData (gap_locations);
}

static Int4 
CalculateGapCoverage 
(BioseqPtr         bsp, 
 SeqLocPtr         slp, 
 SeqMgrFeatContext context,
 Boolean           include_terminal_gaps)
{
  DeltaSeqPtr       dsp;
  SeqLocPtr         loc;
  SeqLitPtr         slip;
  Int4              seq_offset = 0;
  Int4              covered = 0;
  Int4              k, start, stop;

  if (slp == NULL
      || bsp == NULL || !ISA_na (bsp->mol) 
      || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4)
  {
    return 0;
  }
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL && seq_offset < context.right)
  {
    if (dsp->choice == 1)
    {
      loc = (SeqLocPtr) dsp->data.ptrvalue;
      if (loc != NULL)
      {
        seq_offset += SeqLocLen (loc);
      }
    }
    else if (dsp->choice == 2)
    {
      slip = (SeqLitPtr) dsp->data.ptrvalue;
      if (slip != NULL)
      {
        if (IsDeltaSeqKnownGap (dsp) && !DoesDeltaSeqHaveGapTypeOrLinkage(dsp))
        {        
          if (seq_offset <= context.left && seq_offset + slip->length >= context.right)
          {
            /* gap covers entire location */
            return SeqLocLen (slp);
          }
          else if (include_terminal_gaps || 
                   (seq_offset > context.left 
                    && seq_offset + slip->length < context.right))
          {
            /* we only count internal gaps */
            for (k = 0; k < 2 * context.numivals; k += 2) 
            {
              start = context.ivals [k];
              stop = context.ivals [k + 1];
              if (seq_offset <= start && seq_offset + slip->length > stop)
              {
                /* gap covers entire interval */
                covered += stop - start;
              }
              else if (seq_offset > start && seq_offset + slip->length < stop)
              {
                /* interval covers entire gap */
                covered += slip->length;
              }
              else if (seq_offset < start && seq_offset + slip->length > start)
              {
                /* gap covers left end of interval */
                covered += seq_offset + slip->length - start;
              }
              else if (seq_offset < stop && seq_offset + slip->length > stop)
              {
                /* gap covers right end of interval */
                covered += stop - seq_offset;
              }
            }
          }
        }
        seq_offset += slip->length;
      }
    }
    dsp = dsp->next;
  }

  return covered;
  
}

static SeqLocPtr 
RemoveTerminalGapsFromLocation
(BioseqPtr         bsp, 
 SeqLocPtr         slp,
 SeqMgrFeatContext context)
{
  DeltaSeqPtr       dsp;
  SeqLocPtr         loc, slp_copy;
  SeqLitPtr         slip;
  Int4              seq_offset = 0;
  SeqIdPtr          sip;
  Boolean           changed = FALSE, partial5, partial3;

  if (slp == NULL)
  {
    return NULL;
  }
  
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  
  slp_copy = (SeqLocPtr) AsnIoMemCopy (slp, 
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
  if (bsp == NULL || !ISA_na (bsp->mol) 
      || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4)
  {
    return slp_copy;
  }
  
  dsp = (DeltaSeqPtr) bsp->seq_ext;
  while (dsp != NULL && seq_offset < context.right && slp_copy != NULL)
  {
    if (dsp->choice == 1)
    {
      loc = (SeqLocPtr) dsp->data.ptrvalue;
      if (loc != NULL)
      {
        seq_offset += SeqLocLen (loc);
      }
    }
    else if (dsp->choice == 2)
    {
      slip = (SeqLitPtr) dsp->data.ptrvalue;
      if (slip != NULL)
      {
        if (IsDeltaSeqKnownGap (dsp) && !DoesDeltaSeqHaveGapTypeOrLinkage(dsp))
        {        
          if (seq_offset <= context.left && seq_offset + slip->length >= context.right)
          {
            slp_copy = SeqLocFree (slp_copy);
          }
          else if (seq_offset <= context.left && seq_offset + slip->length >= context.left)
          {
            sip = SeqIdDup (SeqLocId (slp_copy));
            slp_copy = SeqLocDeleteEx (slp_copy, sip,
                              context.left, seq_offset + slip->length - 1,
                              FALSE, &changed, &partial5, &partial3);
            sip = SeqIdFree (sip);

          }
          
          if (seq_offset <= context.right && seq_offset + slip->length >= context.right)
          {
            sip = SeqIdDup (SeqLocId (slp_copy));
            slp_copy = SeqLocDeleteEx (slp_copy, sip,
                              seq_offset, context.right,
                              FALSE, &changed, &partial5, &partial3);
            sip = SeqIdFree (sip);
          }
        }
        seq_offset += slip->length;
      }
    }
    dsp = dsp->next;
  }

  if (slp_copy != NULL)  
  {
    SetSeqLocPartial (slp_copy, partial5, partial3);
  }

  return slp_copy;
}


static Int4 TrailingNLength (CharPtr seq)
{
  Int4 seq_len = StringLen (seq);
  Int4 trailing_ns = 0;
  CharPtr cp;

  if (seq == NULL) {
    return 0;
  }
  cp = seq + seq_len - 1;
  while (cp > seq && *cp == 'N') {
    cp--;
    trailing_ns++;
  }
  if (*cp == 'N') {
    trailing_ns++;
  }
  return trailing_ns;
}


static SeqLocPtr 
RemoveTerminalNsFromLocation
(SeqLocPtr slp)
{
  SeqLocPtr         loc = NULL, slp_copy;
  Int4              seq_offset = 0, start, stop;
  Boolean           first = TRUE;
  BioseqPtr         bsp;
  SeqIdPtr          sip;
  Boolean           changed = FALSE, partial5, partial3, tmp5, tmp3;
  CharPtr           buffer = NULL;
  Int4              buflen = 0;
  Int4              leading_ns;
  Int4              trailing_n_start = -1, this_ns;
  Int4              loc_len;
  Uint1             strand;

  if (slp == NULL)
  {
    return NULL;
  }
  
  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp == NULL) {
    return NULL;
  }
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  
  slp_copy = (SeqLocPtr) AsnIoMemCopy (slp, 
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
  if (!ISA_na (bsp->mol))
  {
    return slp_copy;
  }

  while ((loc = SeqLocFindNext (slp, loc)) != NULL) {
    start = SeqLocStart (loc);
    stop = SeqLocStop (loc);
    strand = SeqLocStrand (loc);
    loc_len = stop - start + 1;
    if (buflen < loc_len + 1) {
      buffer = MemFree (buffer);
      buflen = loc_len + 1;
      buffer = (CharPtr) MemNew (sizeof (Char) * buflen);
    }
    SeqPortStreamInt (bsp, start, stop, strand, EXPAND_GAPS_TO_DASHES, (Pointer) (buffer), NULL);
    leading_ns = 0;
    if (first) {
      leading_ns = StringSpn (buffer, "N");
      if (leading_ns > 0) {
        sip = SeqIdDup (SeqLocId (slp_copy));
        tmp5 = FALSE;
        tmp3 = FALSE;
        if (strand == Seq_strand_minus) {
          slp_copy = SeqLocDeleteEx (slp_copy, sip,
                                     stop - leading_ns + 1, stop,
                                     FALSE, &changed, &tmp5, &tmp3);
        } else {
          slp_copy = SeqLocDeleteEx (slp_copy, sip,
                            start, start + leading_ns - 1,
                            FALSE, &changed, &tmp5, &tmp3);
        }
        partial5 |= tmp5;
        partial3 |= tmp3;
        sip = SeqIdFree (sip);
      }
    }

    if (!first || leading_ns != loc_len) {
      this_ns = TrailingNLength (buffer);
      if (this_ns > 0) {
        if (this_ns == loc_len && trailing_n_start > -1) {
          /* add to existing trailing Ns */
        } else {
          if (strand == Seq_strand_minus) {
            trailing_n_start = start + this_ns - 1;
          } else {
            trailing_n_start = stop - this_ns + 1;
          }
        }
      }
    } else {
      trailing_n_start = -1;
    }

    if (first && leading_ns != loc_len) {
      first = FALSE;
    }
  }

  if (trailing_n_start > -1) {
    sip = SeqIdDup (SeqLocId (slp_copy));
    tmp5 = FALSE;
    tmp3 = FALSE;
    if (strand == Seq_strand_minus) {
      slp_copy = SeqLocDeleteEx (slp_copy, sip,
                              start, trailing_n_start,
                              FALSE, &changed, &tmp5, &tmp3);
    } else {
      slp_copy = SeqLocDeleteEx (slp_copy, sip,
                              trailing_n_start, stop,
                              FALSE, &changed, &tmp5, &tmp3);
    }
    partial5 |= tmp5;
    partial3 |= tmp3;
    sip = SeqIdFree (sip);
  }

  if (slp_copy != NULL)  
  {
    SetSeqLocPartial (slp_copy, partial5, partial3);
  }

  buffer = MemFree (buffer);

  return slp_copy;
}


static Boolean AdjustCodingRegionAfterLocAdjustment
(SeqFeatPtr sfp,
 BioseqPtr  bsp,
 SeqLocPtr  adjusted_loc)
{
  Int4              loc_len, adjusted_len;
  Boolean           partial5, partial3;
  BioseqPtr         protbsp;
  SeqFeatPtr        gene_sfp;
  SeqMgrFeatContext gene_context;
  Boolean           rval = FALSE;

  loc_len = SeqLocLen (sfp->location);

  adjusted_len = SeqLocLen (adjusted_loc);
  if (adjusted_loc != NULL && adjusted_len < loc_len)
  {
    gene_sfp = SeqMgrGetOverlappingGene (sfp->location, &gene_context);
    if (gene_sfp != NULL 
        && SeqLocCompare (gene_sfp->location, sfp->location) == SLC_A_EQ_B)
    {
      gene_sfp->location = SeqLocFree (gene_sfp->location);
      gene_sfp->location = (SeqLocPtr) AsnIoMemCopy (adjusted_loc,
                                                     (AsnReadFunc) SeqLocAsnRead,
                                                     (AsnWriteFunc) SeqLocAsnWrite);
    }
  
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = adjusted_loc;
    rval = TRUE;
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    sfp->partial = partial5 || partial3;
    
    protbsp = BioseqFindFromSeqLoc (sfp->product);

    if (sfp->data.choice == SEQFEAT_CDREGION)
    {
      /* adjust frame */
      AdjustFrame (sfp, protbsp);

      /* retranslate coding region */
      SeqEdTranslateOneCDS (sfp, bsp, sfp->idx.entityID, Sequin_GlobalAlign2Seq);
              
      /* set partials on product */
      SetProductSequencePartials (protbsp, partial5, partial3);
    }
  }
  return rval;
}


static void 
AdjustOneCodingRegionWithTerminalGapsOnBioseq 
(BioseqPtr         bsp,
 SeqFeatPtr        sfp,
 SeqMgrFeatContext context)
{
  SeqLocPtr         adjusted_loc;

  if (bsp == NULL || !ISA_na (bsp->mol) 
      || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4
      || sfp == NULL
      || (sfp->data.choice != SEQFEAT_CDREGION && sfp->idx.subtype != FEATDEF_mRNA))
  {
    return;
  }
    
  adjusted_loc = RemoveTerminalGapsFromLocation (bsp, sfp->location, context);

  if (AdjustCodingRegionAfterLocAdjustment (sfp, bsp, adjusted_loc)) 
  {
    if (sfp->idx.subtype == SEQFEAT_CDREGION) {
      AddCDSGapComment (sfp);
    }
  }
  else
  {
    adjusted_loc = SeqLocFree (adjusted_loc);
  }
  
}

static void AdjustCodingRegionsEndingInGapOnBioseq (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        sfp;
  
  if (bsp == NULL || !ISA_na (bsp->mol) 
      || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4)
  {
    return;
  }
    
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context))
  {
    AdjustOneCodingRegionWithTerminalGapsOnBioseq (bsp, sfp, context);
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_mRNA, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_mRNA, &context))
  {
    AdjustOneCodingRegionWithTerminalGapsOnBioseq (bsp, sfp, context);
  }
}

extern void AdjustCodingRegionsEndingInGap (IteM i)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, AdjustCodingRegionsEndingInGapOnBioseq);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void AdjustCodingRegionsWithTerminalNsBioseqCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        sfp;
  SeqLocPtr         adjusted_loc;
  
  if (bsp == NULL || !ISA_na (bsp->mol))
  {
    return;
  }
    
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context))
  {
    adjusted_loc = RemoveTerminalNsFromLocation (sfp->location);
    if (AdjustCodingRegionAfterLocAdjustment (sfp, bsp, adjusted_loc)) 
    {
      AddCDSGapComment (sfp);
    }
    else
    {
      adjusted_loc = SeqLocFree (adjusted_loc);
    }
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_mRNA, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_mRNA, &context))
  {
    adjusted_loc = RemoveTerminalNsFromLocation (sfp->location);
    if (!AdjustCodingRegionAfterLocAdjustment (sfp, bsp, adjusted_loc)) 
    {
      adjusted_loc = SeqLocFree (adjusted_loc);
    }
  }
}


extern void AdjustCodingRegionsWithTerminalNs (IteM i)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, NULL, AdjustCodingRegionsWithTerminalNsBioseqCallback);
  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


typedef struct cdstomiscfeatform
{
  FEATURE_FORM_BLOCK
  
  GrouP  all_or_percent_grp;
  DialoG string_constraint_dlg;
  
  Boolean             convert_all;
  FilterSetPtr        fsp;
  StringConstraintXPtr scp;
  BioseqPtr           bsp;
} CDSToMiscFeatFormData, PNTR CDSToMiscFeatFormPtr;

static void ConvertCDSToMiscFeatFeatureCallback (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  CDSToMiscFeatFormPtr cmffp;
  Int4                 orig_len, covered_len;
  SeqMgrFeatContext    context;
  SeqEntryPtr          oldscope;
  CDStoMiscFeatData    cmfd;

  cmffp = (CDSToMiscFeatFormPtr) userdata;
  
  if (sfp == NULL || cmffp == NULL || cmffp->bsp == NULL)
  {
    return;
  }
  
  sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, cmffp->bsp, sfp->idx.itemID, 0, sfp, &context);
  if (sfp == NULL)
  {
    return;
  }
  
  oldscope = SeqEntrySetScope (NULL);

  orig_len = SeqLocLen (sfp->location);
  covered_len = CalculateGapCoverage (cmffp->bsp, sfp->location, context, cmffp->convert_all);
  if (covered_len >= orig_len / 2 || (cmffp->convert_all && covered_len > 0))
  {
    cmfd.must_have_stops = FALSE;
    cmfd.viral = FALSE;
    cmfd.opts = NULL;
    ConvertCDSToMiscFeat (sfp, &cmfd);
  }
  
  SeqEntrySetScope (oldscope);
}

static void ConvertGappedCodingRegionsToMiscFeatBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  CDSToMiscFeatFormPtr cmffp;
  SeqEntryPtr          sep;
  
  if (bsp == NULL || !ISA_na (bsp->mol) 
      || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4
      || userdata == NULL)
  {
    return;
  }
  
  sep = SeqMgrGetSeqEntryForData (bsp);
  
  cmffp = (CDSToMiscFeatFormPtr) userdata;
  cmffp->bsp = bsp;
  
  OperateOnSeqEntryConstrainedObjects (sep, cmffp->fsp, 
                                       ConvertCDSToMiscFeatFeatureCallback,
                                       NULL, SEQFEAT_CDREGION, FEATDEF_CDS, 0, cmffp);

}

static void AcceptCDSToMiscFeat (ButtoN b)
{
  CDSToMiscFeatFormPtr cmffp;
  SeqEntryPtr          sep;
  
  cmffp = (CDSToMiscFeatFormPtr) GetObjectExtra (b);
  if (cmffp == NULL)
  {
    return;
  }
  
  if (GetValue (cmffp->all_or_percent_grp) == 2)
  {
    cmffp->convert_all = TRUE;
  }
  else
  {
    cmffp->convert_all = FALSE;
  }
  
  cmffp->fsp = FilterSetNew();
  cmffp->fsp->scp = DialogToPointer (cmffp->string_constraint_dlg);

  sep = GetTopSeqEntryForEntityID (cmffp->input_entityID);
  if (sep == NULL) return;

  VisitBioseqsInSep (sep, cmffp, ConvertGappedCodingRegionsToMiscFeatBioseqCallback);
  DeleteMarkedObjects (cmffp->input_entityID, 0, NULL);
  RenormalizeNucProtSets (sep, TRUE);
  ObjMgrSetDirtyFlag (cmffp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, cmffp->input_entityID, 0, 0);
  Remove (cmffp->form);  
}

extern void ConvertCodingRegionsWithInternalKnownGapToMiscFeat (IteM i)
{
  BaseFormPtr          bfp;
  CDSToMiscFeatFormPtr cmffp;
  WindoW               w;
  GrouP                h, c;
  ButtoN               b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  cmffp = (CDSToMiscFeatFormPtr) MemNew (sizeof (CDSToMiscFeatFormData));
  if (cmffp == NULL) return;
    
  w = FixedWindow (-50, -33, -10, -10, "Convert Coding Regions With Gaps to Misc_Feat", StdCloseWindowProc);
  SetObjectExtra (w, cmffp, StdCleanupFormProc);
  cmffp->form = (ForM) w;
  cmffp->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  cmffp->all_or_percent_grp = HiddenGroup (h, 0, 2, NULL);
  SetGroupSpacing (cmffp->all_or_percent_grp, 10, 10);
  RadioButton (cmffp->all_or_percent_grp, "Convert only when internal gap covers 50% or more of the coding region");
  RadioButton (cmffp->all_or_percent_grp, "Convert all coding regions with gaps (both terminal and internal)");
  SetValue (cmffp->all_or_percent_grp, 1);
  
  cmffp->string_constraint_dlg = StringConstraintDialogX (h, "Where feature text", FALSE);

  c = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", AcceptCDSToMiscFeat);
  SetObjectExtra (b, cmffp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) cmffp->all_or_percent_grp,
                              (HANDLE) cmffp->string_constraint_dlg,
                              (HANDLE) c,
                              NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* Apply_AddCDS () --                                                  */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void Apply_AddCDS (Uint2        entityID,
			  SeqEntryPtr  sep,
			  SeqEntryPtr  nsep,
			  ApplyFormPtr afp,
			  Boolean      suppressDups)
{
  ByteStorePtr       bs;
  BioseqPtr          bsp;
  Char               ch;
  CdRegionPtr        crp;
  RnaRefPtr          rrp;
  ValNodePtr         descr;
  Int2               genCode;
  Int2               i;
  MolInfoPtr         mip;
  SeqEntryPtr        old;
  CharPtr            prot;
  ProtRefPtr         prp;
  SeqEntryPtr        psep;
  CharPtr            ptr;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  Char               str [128];
  ValNodePtr         vnp;
  SeqEntryPtr        parent_sep;
  SeqEntryPtr        gene_sep;
  SeqFeatPtr         prot_sfp, mRNA_sfp;

  /* If necessary then check for duplication before adding */

  if (suppressDups &&
      entityID > 0 &&
      AlreadyHasFeatOrDesc (sep, SEQFEAT_CDREGION, 0, 0))
    return;

  /* determine the parent of this sequence (for use when segmented) */
  parent_sep = NULL;
  if (IS_Bioseq (sep))
  {
    parent_sep = GetBestTopParentForData (entityID, sep->data.ptrvalue);
  }
  if (parent_sep == NULL)
  {
    parent_sep = sep;
  }

  /*Create a new CDS feature */

  genCode = SeqEntryToGeneticCode (parent_sep, NULL, NULL, 0);
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (NULL == crp)
    return;
  
  sfp = CreateNewFeature (nsep, NULL, SEQFEAT_CDREGION, NULL); 	

  if (NULL == sfp)
    return;
  
  sfp->data.value.ptrvalue = (Pointer) crp;
  
  /* set the comment */
  if (afp->feature_details_data->featcomment != NULL) {
      sfp->comment = StringSave (afp->feature_details_data->featcomment);
  }

  /* adjust the location of the new feature */
  SetApplyFeatureLocation (sfp, afp);
  
  /* Fill in the fields of the new CDS feature */
  if (afp->feature_details_data->reading_frame < 1 
      || afp->feature_details_data->reading_frame > 3)
  {
    if (!SetBestFrameByLocation (sfp)) {
      str [0] = '\0';
      if (IS_Bioseq (nsep)) {
        bsp = (BioseqPtr) nsep->data.ptrvalue;
        if (bsp != NULL) {
          sip = SeqIdFindBest (bsp->id, 0);
          SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str) - 1);
        }
      }
      (afp->errcount)++;
      ValNodeCopyStr (&(afp->ambigList), 0, str);
    }
  }
  else
  {
    crp->frame = afp->feature_details_data->reading_frame;
  }

  /* Create corresponding protein sequence data for the CDS */

  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (NULL == bs)
    return;

  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (NULL == prot)
    return;

  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = (Int2) StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  }

  /* Create the product protein Bioseq */
  
  bsp = BioseqNew ();
  if (NULL == bsp)
    return;
  
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  old = SeqEntrySetScope (sep);
  bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  SeqEntrySetScope (old);
  
  /* Create a new SeqEntry for the Prot Bioseq */
  
  psep = SeqEntryNew ();
  if (NULL == psep)
    return;
  
  psep->choice = 1;
  psep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
  
  /* Add a descriptor to the protein Bioseq */
  
  mip = MolInfoNew ();
  if (NULL == mip)
    return;
  
  mip->biomol = 8;
  mip->tech = 8;
  if (afp->noLeft && afp->noRight) {
    mip->completeness = 5;
  } else if (afp->noLeft) {
    mip->completeness = 3;
  } else if (afp->noRight) {
    mip->completeness = 4;
  }
  vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
  if (NULL == vnp)
    return;
  
  vnp->data.ptrvalue = (Pointer) mip;
  
  /**/
  
  descr = ExtractBioSourceAndPubs (parent_sep);

  AddSeqEntryToSeqEntry (parent_sep, psep, TRUE);
  nsep = FindNucSeqEntry (parent_sep);
  ReplaceBioSourceAndPubs (parent_sep, descr);
  SetSeqFeatProduct (sfp, bsp);
  
  /* create a full-length protein feature for the new protein sequence */
  if (! StringHasNoText (afp->feature_details_data->protName) 
      && ! StringHasNoText (afp->feature_details_data->protDesc))
  {
    prp = CreateNewProtRef (afp->feature_details_data->protName, 
                            afp->feature_details_data->protDesc, 
                            NULL, NULL);
  }
  else if (!StringHasNoText (afp->feature_details_data->protName))
  {
    prp = CreateNewProtRef (afp->feature_details_data->protName, NULL, NULL, NULL);
  }
  else if (!StringHasNoText (afp->feature_details_data->protDesc))
  {
    prp = CreateNewProtRef (NULL, afp->feature_details_data->protDesc, NULL, NULL);
  }
  else
  { 
    prp = ProtRefNew ();
  }
  
  if (prp != NULL) {
    prot_sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (prot_sfp != NULL) {
      prot_sfp->data.value.ptrvalue = (Pointer) prp;
      SetSeqLocPartial (prot_sfp->location, afp->noLeft, afp->noRight);
      prot_sfp->partial = (afp->noLeft || afp->noRight);
    }
  }
  
  /* after the feature has been created, then adjust it for gaps */
  /* Note - this step may result in multiple coding regions being created. */
  AdjustCDSLocationsForUnknownGapsCallback (sfp, NULL);

  /* if requested, also add mRNA feature */
  if (afp->also_add_mRNA_btn != NULL && GetStatus(afp->also_add_mRNA_btn)) {
    mRNA_sfp = CreateNewFeature (nsep, NULL, SEQFEAT_RNA, NULL); 
    rrp = RnaRefNew ();
    if (rrp != NULL) {
      rrp->type = 2;
      if (!StringHasNoText(afp->feature_details_data->protName)) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave (afp->feature_details_data->protName);
      }    
    }
	mRNA_sfp->data.value.ptrvalue = rrp;
    SetApplyFeatureLocation (mRNA_sfp, afp);
    AdjustCDSLocationsForUnknownGapsCallback (mRNA_sfp, NULL);
  }

  /* Create a Gene ref feature on the nuc seq or segment */
  /* we can only create a feature where the sep->choice is 1 */
  if (sep->choice == 1)
  {
    gene_sep = sep;
  }
  else
  {
    gene_sep = nsep;
  }

  if (! StringHasNoText (afp->feature_details_data->geneName)) {
    if (entityID > 0 
        && suppressDups
        && AlreadyHasFeatOrDesc (gene_sep, SEQFEAT_GENE, 0, 0))
    {
      return;
    }

    ApplyGene (afp->feature_details_data->geneName, afp, gene_sep, sfp);
  }
}

static void CheckTitle (SeqEntryPtr sep, GetSamplePtr gsp)
{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp = NULL, part_set;
  SeqDescrPtr       sdp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL
      || gsp == NULL)
  {
    return;
  }
  
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL 
        && bssp->_class == BioseqseqSet_class_nuc_prot 
        && bssp->seq_set != NULL) 
    {
      sep = bssp->seq_set;
    }
  }
  
  if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL 
        && bssp->_class == BioseqseqSet_class_segset 
        && bssp->seq_set != NULL) 
    {
      /* first do parts */
      if (bssp->seq_set->next != NULL 
          && IS_Bioseq_set (bssp->seq_set->next)
          && bssp->seq_set->next->data.ptrvalue != NULL)
      {
        part_set = (BioseqSetPtr) bssp->seq_set->next->data.ptrvalue;
        if (part_set->_class == BioseqseqSet_class_parts)
        {
          for (sep = part_set->seq_set; sep != NULL; sep = sep->next)
          {
            CheckTitle (sep, gsp);
          }
        }
      }
      
      /* now do master */
      sep = bssp->seq_set;
    }
  }  

  if (!IS_Bioseq (sep))
  {
    return;
  }
    
  bsp = (BioseqPtr) sep->data.ptrvalue;
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp != NULL)
  {
    gsp->num_found ++;
    if (gsp->sample_text == NULL)
    {
      gsp->sample_text = StringSave (sdp->data.ptrvalue);
    }
    else if (StringCmp (gsp->sample_text, sdp->data.ptrvalue) != 0)
    {
      gsp->all_same = FALSE;
    }
  }
  
  
}

static void ApplyOneTitle (SeqEntryPtr sep, ExistingTextPtr etp, CharPtr defline)
{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp, part_set;
  SeqDescrPtr       sdp;
  
  if (sep == NULL || sep->data.ptrvalue == NULL
      || StringHasNoText (defline))
  {
    return;
  }
  
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL 
        && bssp->_class == BioseqseqSet_class_nuc_prot 
        && bssp->seq_set != NULL) 
    {
      sep = bssp->seq_set;
    }
  }
  
  if (IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL 
        && bssp->_class == BioseqseqSet_class_segset 
        && bssp->seq_set != NULL) 
    {
      /* first do parts */
      if (bssp->seq_set->next != NULL 
          && IS_Bioseq_set (bssp->seq_set->next)
          && bssp->seq_set->next->data.ptrvalue != NULL)
      {
        part_set = (BioseqSetPtr) bssp->seq_set->next->data.ptrvalue;
        if (part_set->_class == BioseqseqSet_class_parts)
        {
          for (sep = part_set->seq_set; sep != NULL; sep = sep->next)
          {
            ApplyOneTitle (sep, etp, defline);
          }
        }
      }
      
      /* now do master */
      sep = bssp->seq_set;
    }
  }    

  if (!IS_Bioseq (sep))
  {
    return;
  }
    
  bsp = (BioseqPtr) sep->data.ptrvalue;
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }

  if (sdp == NULL)
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_title);
  }
  
  if (sdp != NULL) {
    sdp->data.ptrvalue = HandleExistingText (sdp->data.ptrvalue,
                                             StringSave (defline), etp);
  }
  
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* RealApplyBioFeatToAll () --                                         */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void RealApplyBioFeatToAll (Uint2        entityID,
				                           SeqEntryPtr  sep,
				                           SeqEntryPtr  nsep,
                                   ApplyFormPtr afp,
				                           Boolean      suppressDups)

{
  ImpFeatPtr         ifp;
  SeqFeatPtr         sfp = NULL;
  SeqFeatPtr         gene_sfp;
  Boolean            put_comment_on_gene = FALSE;

  /* Check parameters */

  if (sep == NULL || nsep == NULL || afp == NULL)
    return;

  /* Add a title feature */

  if (afp->type == CHECK_TITLE)
  {
    CheckTitle (sep, afp->gsp);
    return;
  }
  else if (afp->type == ADD_TITLE) 
  {
    if (entityID == 0 && SeqEntryGetTitle (sep) != NULL)
    {
      return;
    }
    ApplyOneTitle (sep, afp->etp, afp->feature_details_data->defline);
    return;
  }

  /* Add a CDS feature */

  else if (afp->type == ADD_CDS)
    Apply_AddCDS (entityID, sep, nsep, afp, suppressDups);

  /* Add an rRNA feature */

  else if (afp->type == ADD_RRNA) {
    if (suppressDups && entityID > 0 
        && AlreadyHasRNA (sep, afp->feature_details_data->rnaType)) {
      return;
    }
   
    sfp = CreateNewFeature (nsep, NULL, SEQFEAT_RNA, NULL);
    ApplyRnaTypeToSeqFeat (sfp, afp->feature_details_data->rnaType);
    if (! StringHasNoText (afp->feature_details_data->rnaName)) {
      ApplyProductToRNA (sfp, afp->feature_details_data->rnaName);
    }
    
    SetApplyFeatureLocation (sfp, afp);      
    AddToComment (sfp, afp->feature_details_data->featcomment);

    if (! StringHasNoText (afp->feature_details_data->geneName)) {
      if (entityID > 0 
        && suppressDups
	      && AlreadyHasFeatOrDesc (sep, SEQFEAT_GENE, 0, 0))
      {
	       return;        
      }
      ApplyGene (afp->feature_details_data->geneName, afp, nsep, sfp);

    }
  }

  /* Add an Import feature */

  else if (afp->type == ADD_IMP) {
    if (afp->feature_details_data->featdef_choice == FEATDEF_GENE)
    {
      put_comment_on_gene = TRUE;
    }
    else 
    {
      ifp = ImpFeatNew ();
      if (ifp != NULL) {
        ifp->key = StringSave (afp->feature_details_data->featdef_name);
        sfp = CreateNewFeature (nsep, NULL, SEQFEAT_IMP, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = (Pointer) ifp;
          SetApplyFeatureLocation (sfp, afp);
          if (! StringHasNoText (afp->feature_details_data->featcomment)) 
          {
            sfp->comment = StringSave (afp->feature_details_data->featcomment);
          }
          sfp->qual = DialogToPointer (afp->gbquals);
        }
      }
    }
    if (! StringHasNoText (afp->feature_details_data->geneName) 
		|| ! StringHasNoText (afp->feature_details_data->geneDesc)) 
    {
      if (entityID > 0 
          && suppressDups
	        && AlreadyHasFeatOrDesc (sep, SEQFEAT_GENE, 0, 0)) 
      {
	      return;        
      }
      gene_sfp = ApplyGene (afp->feature_details_data->geneName, afp, nsep, sfp);

      if (gene_sfp != NULL && ! StringHasNoText (afp->feature_details_data->featcomment) && put_comment_on_gene)
      {
        gene_sfp->comment = StringSave (afp->feature_details_data->featcomment);
      }
    }
  }
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyBioFeatToRaw () --                                             */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void ApplyBioFeatToRaw (Uint2        entityID,
			       SeqEntryPtr  parentSep,
			       SeqEntryPtr  sep,
			       ApplyFormPtr afp,
			       Int2         onlythis)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int2          count;
  SeqEntryPtr   nucSep;

  /* Check parameters */

  if (sep == NULL || afp == NULL)
    return;

  /* If it is a set then recurse until we get to the raw Bioseq */

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      if (onlythis != 0 && bssp->_class == BioseqseqSet_class_parts) {
        for (sep = bssp->seq_set, count = 1;
             sep != NULL && count != onlythis;
             sep = sep->next, count++)
	  continue;
        if (sep != NULL) {
          ApplyBioFeatToRaw (entityID, parentSep, sep, afp, onlythis);
        }
      } else {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          ApplyBioFeatToRaw (entityID, parentSep, sep, afp, onlythis);
        }
      }
      return;
    }
  }

  /* Get the nucleotide Bioseq */

  nucSep = FindNucSeqEntry (sep);
  if (nucSep == NULL)
    return;

  bsp = (BioseqPtr) nucSep->data.ptrvalue;
  if (bsp == NULL)
    return;
  
  /* If we've got a raw Bioseq then do the apply */

  if (bsp->repr == Seq_repr_raw) {
    RealApplyBioFeatToAll (entityID, nucSep, nucSep, afp, FALSE);
    /*
    RealApplyBioFeatToAll (entityID, parentSep, nucSep, afp, FALSE);
    */
  }
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyBioFeatToAll () --                                             */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void ApplyBioFeatToAll (Uint2        entityID,
			       SeqEntryPtr  sep,
			       ApplyFormPtr afp)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   nsep;
  Int2          onlythis;
  Char          str [32];
  Boolean       suppressDups = FALSE;

  /* Check parameters */

  if (sep == NULL || afp == NULL)
    return;
  
  if (afp->add_to_seq_with_like_feature != NULL)
  {
    suppressDups = ! GetStatus (afp->add_to_seq_with_like_feature);
  }

  /* If it is a set then recurse until we get to a Bioseq */

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      if (bssp->_class == BioseqseqSet_class_genbank || IsPopPhyEtcSet (bssp->_class)) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          ApplyBioFeatToAll (entityID, sep, afp);
        }
        return;
      } else if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
        /* just the first item */
        ApplyBioFeatToAll (entityID, bssp->seq_set, afp);
        return;
      }
    }
  }

  /* Get the nucleotide Bioseq */

  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL)
    return;

  /* Apply the feature */

  if (afp->applyToParts != NULL && GetStatus (afp->applyToParts)) {
    GetTitle (afp->onlyThisPart, str, sizeof (str));
    if (! StrToInt (str, &onlythis)) {
      onlythis = 0;
    }
    ApplyBioFeatToRaw (entityID, sep, sep, afp, onlythis);
  } else {
    RealApplyBioFeatToAll (entityID, sep, nsep, afp, suppressDups);
  }
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* ApplyBioFeatToAll () --                                             */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void ApplyBioFeatToList (Uint2        entityID,
			       SeqEntryPtr  sep,
			       ApplyFormPtr afp,
			       ValNodePtr   id_list)

{
  BioseqPtr     bsp;
  SeqEntryPtr   nsep, oldscope;
  Int2          onlythis;
  Char          str [32];
  Boolean       suppressDups = FALSE;
  ValNodePtr    vnp;

  /* Check parameters */

  if (sep == NULL || afp == NULL || id_list == NULL)
    return;
  
  if (afp->add_to_seq_with_like_feature != NULL)
  {
    suppressDups = ! GetStatus (afp->add_to_seq_with_like_feature);
  }

  oldscope = SeqEntrySetScope (sep);
  for (vnp = id_list; vnp != NULL; vnp = vnp->next) {
     bsp = BioseqFind (vnp->data.ptrvalue);
     if (bsp != NULL) {
       sep = SeqMgrGetSeqEntryForData (bsp);
       
       /* Get the nucleotide Bioseq */
       nsep = FindNucSeqEntry (sep);
       if (nsep == NULL)
         continue;

       /* Apply the feature */

       if (afp->applyToParts != NULL && GetStatus (afp->applyToParts)) {
         GetTitle (afp->onlyThisPart, str, sizeof (str));
         if (! StrToInt (str, &onlythis)) {
           onlythis = 0;
         }
         ApplyBioFeatToRaw (entityID, sep, sep, afp, onlythis);
       } else {
         RealApplyBioFeatToAll (entityID, sep, nsep, afp, suppressDups);
       }
    }
  }
}


Int2 ApplyAnnotationToAll (Int2 type, SeqEntryPtr sep,
                           ButtoN partialLft, ButtoN partialRgt,
                           TexT geneName, TexT protName, 
                           TexT protDesc, TexT rnaName,
                           TexT featcomment, TexT defline)

{
  ApplyFormData  afd;
  RnaTypeData    rtd;

  MemSet ((Pointer) (&afd), 0, sizeof (ApplyFormData));
  afd.type = type;
  afd.errcount = 0;
  afd.ambigList = NULL;
  afd.partial5 = partialLft;
  afd.partial3 = partialRgt;
  afd.noLeft = GetStatus (afd.partial5);
  afd.noRight = GetStatus (afd.partial3);
  afd.feature_details_data = BatchApplyFeatureDetailsNew ();
  if (afd.feature_details_data == NULL)
  {
    return 1;
  }
  
  if (ADD_RRNA == type)
  {
    rtd.ncrna_class = NULL;
    rtd.rna_featdef = FEATDEF_rRNA;
    afd.feature_details_data->rnaType = &rtd;
  }

  if (geneName != NULL && ! TextHasNoText (geneName))
    afd.feature_details_data->geneName = SaveStringFromText (geneName);
  if (protName != NULL && ! TextHasNoText (protName))
    afd.feature_details_data->protName = SaveStringFromText (protName);
  if (protDesc != NULL && ! TextHasNoText (protDesc))
    afd.feature_details_data->protDesc = SaveStringFromText (protDesc);
  if (rnaName != NULL && ! TextHasNoText (rnaName))
    afd.feature_details_data->rnaName = SaveStringFromText (rnaName);
  if (featcomment != NULL && ! TextHasNoText (featcomment))
    afd.feature_details_data->featcomment = SaveStringFromText (featcomment);
  if (defline != NULL && ! TextHasNoText (defline))
    afd.feature_details_data->defline = SaveStringFromText (defline);
  ApplyBioFeatToAll (0, sep, &afd);
  if (afd.type == ADD_CDS) {
    return afd.errcount;
  }
  return 0;
}

static void CommonApplyToAllProcBfpInfo (Uint2 entityID,
                                         Uint4 itemID,
                                         Uint2 itemtype,
                                         Int2 type);

static void NowReadyToApplyToAll (ApplyFormPtr afp, DialoG gbquals)

{
  Uint2        parenttype;
  Pointer      parentptr;
  CharPtr      plural;
  SeqEntryPtr  sep;
  CharPtr      tmp;
  SeqEntryPtr  top;
  ValNodePtr   vnp;
  Char         path [PATH_MAX];
  FILE         *fp;
  ValNodePtr   id_list = NULL;

  if (afp == NULL) return;
  afp->gbquals = gbquals;
  sep = GetTopSeqEntryForEntityID (afp->input_entityID);
  if (sep == NULL) return;
  Hide (afp->form);
  WatchCursor ();
  Update ();
  afp->noLeft = GetStatus (afp->partial5);
  afp->noRight = GetStatus (afp->partial3);
  top = sep;
  if (afp->type == ADD_CDS) {
    GetSeqEntryParent (top, &parentptr, &parenttype);
  }

  if (afp->all_or_some_grp != NULL && GetValue (afp->all_or_some_grp) == 2) {  
    tmp = SaveStringFromText (afp->accession_list_txt);
    id_list = ParseAccessionNumberListFromString(tmp, sep);
    tmp = MemFree (tmp);
    if (id_list == NULL) {
      ArrowCursor();
      Update();
      Show (afp->form);
      return;
    }
    ApplyBioFeatToList (afp->input_entityID, sep, afp, id_list);
    id_list = FreeSeqIdList (id_list);
    tmp = MemFree (tmp);
  } else {
    ApplyBioFeatToAll (afp->input_entityID, sep, afp);
  }
  ArrowCursor ();
  Update ();
  if (afp->errcount > 0 && afp->type == ADD_CDS) {
    TmpNam (path);
    fp = FileOpen (path, "w");
    if (fp != NULL) {
      if (afp->errcount > 1) {
        plural = "records";
      } else {
        plural = "record";
      }
      fprintf (fp, "Possible ambiguous frames detected in %d %s\n",
               (int) afp->errcount, plural);
      for (vnp = afp->ambigList; vnp != NULL; vnp = vnp->next) {
        tmp = (CharPtr) vnp->data.ptrvalue;
        fprintf (fp, "%s\n", tmp);
      }
      FileClose (fp);
      LaunchGeneralTextViewer (path, "Ambiguous Frames");
      FileRemove (path);
    }
  }
  ObjMgrSetDirtyFlag (afp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, afp->input_entityID, 0, 0);
  if (GetStatus (afp->leaveDlgUp))
  {
    afp->ambigList = ValNodeFreeData (afp->ambigList);
    Show (afp->form);
  }
  else
  {
    Remove (afp->form);
  }
}

typedef struct qualsform {
  FEATURE_FORM_BLOCK

  ApplyFormPtr  afp;
} QualsForm, PNTR QualsFormPtr;

static void CallNowReady (ButtoN b)

{
  QualsFormPtr  qfp;

  qfp = (QualsFormPtr) GetObjectExtra (b);
  if (qfp == NULL) return;
  Hide (qfp->form);
  NowReadyToApplyToAll (qfp->afp, qfp->gbquals);
  Remove (qfp->form);
}


static void DoTheApplyToAllProc (ButtoN b)

{
  ApplyFormPtr       afp;
  CharPtr            name;
  QualsFormPtr       qfp;
  WindoW             w;

  afp = GetObjectExtra (b);
  if (afp == NULL) {
    Remove (ParentWindow (b));
    return;
  }
  
  afp->feature_details_data = DialogToPointer (afp->feature_details_dlg);
  
  if (afp->type == ADD_IMP) {
    /* if gene, do not collect quals */
    if (afp->feature_details_data->featdef_choice != FEATDEF_GENE)
    {      
      qfp = (QualsFormPtr) MemNew (sizeof (QualsForm));
      if (qfp != NULL) {
        Hide (afp->form);
        Update ();
        name = afp->feature_details_data->featdef_name;
        qfp->afp = afp;
        w = FixedWindow (-50, -33, -10, -10, "Qualifiers", StdCloseWindowProc);
        SetObjectExtra (w, qfp, StdCleanupFormProc);
        qfp->form = (ForM) w;
        CreateStandardEditMenu (w);
        qfp->gbquals = NewCreateImportFields (w, name, NULL, FALSE);
        b = PushButton (w, "Okay", CallNowReady);
        SetObjectExtra (b, qfp, NULL);
        AlignObjects (ALIGN_CENTER, (HANDLE) qfp->gbquals, (HANDLE) b, NULL);
        RealizeWindow (w);
        Show (w);
        Select (w);
      }
    } else {
      NowReadyToApplyToAll (afp, NULL);
    }
  }
  else if (afp->type == ADD_TITLE)
  {
    afp->gsp = GetSampleNew ();
    afp->type = CHECK_TITLE;
    NowReadyToApplyToAll (afp, NULL);
    afp->etp = GetExistingTextHandlerInfo (afp->gsp == NULL ? 0 : afp->gsp->num_found, FALSE);
    afp->gsp = GetSampleFree (afp->gsp);
    afp->type = ADD_TITLE;
    if (afp->etp == NULL
        || afp->etp->existing_text_choice != eExistingTextChoiceCancel)
    {
      NowReadyToApplyToAll (afp, NULL);
    }
    afp->etp = MemFree (afp->etp);
  } else {
    NowReadyToApplyToAll (afp, NULL);
  }
}

static void ApplyToPartsProc (ButtoN b)

{
  ApplyFormPtr  afp;

  afp = (ApplyFormPtr) GetObjectExtra (b);
  if (afp == NULL) return;
  if (GetStatus (b)) {
    SafeEnable (afp->onlyThisPart);
  } else {
    SafeDisable (afp->onlyThisPart);
  }
}

static void LookForParts (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;
  BoolPtr       rsult;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
      rsult = (BoolPtr) mydata;
      if (rsult != NULL) {
        *rsult = TRUE;
      }
    }
  }
}

extern Boolean HasPartsSet (SeqEntryPtr sep);
extern Boolean HasPartsSet (SeqEntryPtr sep)

{
  Boolean  rsult = FALSE;

  SeqEntryExplore (sep, (Pointer) (&rsult), LookForParts);
  return rsult;
}

static void EnableApplyCdsCoords (GrouP g)
{
  ApplyFormPtr afp;
  
  afp = (ApplyFormPtr) GetObjectExtra (g);
  if (afp == NULL) return;
  if (GetValue (g) == 1)
  {
  	Disable (afp->left_end);
  	Disable (afp->right_end);
  }
  else
  {
  	Enable (afp->left_end);
  	Enable (afp->right_end);
  }
}

static void ApplyMessageProc (ForM f, Int2 mssg)

{
  ApplyFormPtr  afp;

  afp = (ApplyFormPtr) GetObjectExtra (f);
  if (afp != NULL) {
    if (afp->appmessage != NULL) {
      afp->appmessage (f, mssg);
    }
  }
}

static void CleanupApplyToAllForm (GraphiC g, VoidPtr data)

{
  ApplyFormPtr  afp;

  afp = (ApplyFormPtr) data;
  if (afp != NULL) {
    afp->feature_details_data = BatchApplyFeatureDetailsFree (afp->feature_details_data);
    ValNodeFreeData (afp->ambigList);
  }
  StdCleanupFormProc (g, data);
}

static void ChangeAllOrSome (GrouP g)
{
  ApplyFormPtr       afp;

  afp = (ApplyFormPtr) GetObjectExtra (g);
  if (afp != NULL) {
    if (GetValue(g) == 1) {
      Disable (afp->accession_list_txt);
    } else {
      Enable (afp->accession_list_txt);
    }
  }
}


static void EnableApplyFeatureAccept (Pointer data)
{
  ApplyFormPtr afp;

  afp = (ApplyFormPtr) data;
  if (afp == NULL) return;

  if (OkToAcceptBatchApplyFeatureDetails (afp->feature_details_dlg))
  {
    Enable (afp->accept);
  } 
  else
  {
    Disable (afp->accept);
  }
}


static void CommonApplyToAllProcBfpInfo (Uint2 entityID,
                                         Uint4 itemID,
                                         Uint2 itemtype,
                                         Int2 type)
{
  ApplyFormPtr       afp;
  GrouP              c;
  GrouP              g = NULL;
  GrouP              h;
  GrouP              r2, r3, r4;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  GrouP              x;
  GrouP              parts_group = NULL;
  GrouP              feature_details = NULL;
  GrouP              indexer_only_group = NULL;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  afp = (ApplyFormPtr) MemNew (sizeof (ApplyFormData));
  if (afp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Automatic Processing", StdCloseWindowProc);
  SetObjectExtra (w, afp, CleanupApplyToAllForm);
  afp->form = (ForM) w;
  afp->formmessage = ApplyMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    afp->appmessage = sepp->handleMessages;
  }

  afp->input_entityID = entityID;
  afp->input_itemID = itemID;
  afp->input_itemtype = itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  afp->type = type;
  afp->errcount = 0;

  if (HasPartsSet (sep)) {
    parts_group = HiddenGroup (h, 1, 0, NULL);
    afp->applyToParts = CheckBox (parts_group, "Apply to segmented parts, not segmented sequence", ApplyToPartsProc);
    SetObjectExtra (afp->applyToParts, afp, NULL);
    x = HiddenGroup (parts_group, 2, 0, NULL);
    StaticPrompt (x, "Apply only to particular numbered segment", 0, dialogTextHeight, programFont, 'l');
    afp->onlyThisPart = DialogText (x, "", 4, NULL);
    Disable (afp->onlyThisPart);
  }

  afp->strand_group = NULL;
  afp->add_to_seq_with_like_feature = NULL;
  afp->also_add_mRNA_btn = NULL;

  if (type == ADD_CDS || type == ADD_RRNA || type == ADD_IMP) {
    /* create group to hold feature details */
    feature_details = HiddenGroup (h, -1, 0, NULL);
    
    /* group for feature completeness */
    g = HiddenGroup (feature_details, 2, 0, NULL);
    afp->partial5 = CheckBox (g, "Incomplete at 5' end", NULL);
    afp->partial3 = CheckBox (g, "Incomplete at 3' end", NULL);
    
    /* group for strand */
    afp->strand_group = HiddenGroup (feature_details, 2, 0, NULL);
    SetObjectExtra (afp->strand_group, afp, NULL);
    RadioButton (afp->strand_group, "Plus Strand");
    RadioButton (afp->strand_group, "Minus Strand");
    SetValue (afp->strand_group, 1);
    
    /* coordinates */
    if (indexerVersion)
    {
      indexer_only_group = HiddenGroup (feature_details, -1, 0, NULL);
      r2 = HiddenGroup (indexer_only_group, 5, 0, NULL);
      afp->use_whole_interval = HiddenGroup (r2, 0, 2, EnableApplyCdsCoords);
      SetObjectExtra (afp->use_whole_interval, afp, NULL);
      RadioButton (afp->use_whole_interval, "Use Whole Sequence Interval");
      RadioButton (afp->use_whole_interval, "Use these coordinates:");
      r3 = HiddenGroup (r2, 0, 2, NULL);
      StaticPrompt (r3, "", 0, dialogTextHeight, programFont, 'l');
      r4 = HiddenGroup (r3, 4, 0, NULL);
      StaticPrompt (r4, "From", 0, dialogTextHeight, programFont, 'l');
      afp->left_end = DialogText (r4, "1", 5, NULL);
      StaticPrompt (r4, "To", 0, dialogTextHeight, programFont, 'l');
      afp->right_end = DialogText (r4, "1", 5, NULL);
      SetValue (afp->use_whole_interval, 1);
      Disable (afp->left_end);
      Disable (afp->right_end);
      
      /* apply to some sequences or all sequences */
      afp->all_or_some_grp = HiddenGroup (indexer_only_group, 1, 0, ChangeAllOrSome);
      SetObjectExtra (afp->all_or_some_grp, afp, NULL);
      RadioButton (afp->all_or_some_grp, "Apply to all sequences");
      RadioButton (afp->all_or_some_grp, "Apply to sequences in this list");
      afp->accession_list_txt = DialogText (afp->all_or_some_grp, "", 25, NULL);
      SetValue (afp->all_or_some_grp, 1);
      ChangeAllOrSome(afp->all_or_some_grp);
      
      /* add to features that already have same */
      if (type == ADD_CDS)
      {
        afp->add_to_seq_with_like_feature = CheckBox (indexer_only_group, "Also add to sequences that already have a CDS", NULL);
        afp->also_add_mRNA_btn = CheckBox (indexer_only_group, "Also add mRNA", NULL);
        SetStatus (afp->also_add_mRNA_btn, FALSE);
      }
      else if (type == ADD_RRNA)
      {
        afp->add_to_seq_with_like_feature = CheckBox (indexer_only_group, "Also add to sequences that already have an rRNA", NULL);
      }
      else if (type == ADD_IMP)
      {
        afp->add_to_seq_with_like_feature = CheckBox (indexer_only_group, "Also add to sequences that already have this feature", NULL);
      }
      SafeSetStatus (afp->add_to_seq_with_like_feature, TRUE);
      AlignObjects (ALIGN_CENTER, (HANDLE) r2,
                                   (HANDLE) afp->add_to_seq_with_like_feature,
                                   (HANDLE) afp->also_add_mRNA_btn,
                                   NULL);
    }
    else
    {
      afp->use_whole_interval = NULL;
      afp->all_or_some_grp = NULL;
      afp->accession_list_txt = NULL;
    }
    
    AlignObjects (ALIGN_CENTER, (HANDLE) g,
                                (HANDLE) afp->strand_group,
                                (HANDLE) indexer_only_group,
                                NULL);
  }
  afp->feature_details_dlg = BatchApplyFeatureDetailsDialog (h, type, EnableApplyFeatureAccept, afp);
        

  afp->gbquals = NULL;

  c = HiddenGroup (h, 4, 0, NULL);
  afp->accept = DefaultButton (c, "Accept", DoTheApplyToAllProc);
  SetObjectExtra (afp->accept, afp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  afp->leaveDlgUp = CheckBox (c, "Leave Dialog Up", NULL);

  if (parts_group == NULL)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) c,
                                (HANDLE) afp->feature_details_dlg,    
                                (HANDLE) feature_details, 
                                NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) parts_group,
                                (HANDLE) c, 
                                (HANDLE) afp->feature_details_dlg,
                                (HANDLE) feature_details, 
                                NULL);
  }
  
  EnableApplyFeatureAccept (afp);
  RealizeWindow (w);
  Show (w);
  SendMessageToDialog (afp->feature_details_dlg, VIB_MSG_ENTER);
  Update ();
}

extern void CommonApplyToAllProc (BaseFormPtr bfp, Int2 type)
{
  if (bfp == NULL) return;
  CommonApplyToAllProcBfpInfo (bfp->input_entityID,
                               bfp->input_itemID,
                               bfp->input_itemtype, type);
}

static void CommonApplyToAllProcMenuItem (IteM i, Int2 type)
{
  BaseFormPtr bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  CommonApplyToAllProc (bfp, type);
}


extern void ApplyCDSBtn (ButtoN b)
{
  BaseFormPtr bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) {
    return;
  }
  CommonApplyToAllProc (bfp, ADD_CDS);
}


extern void ApplyTitle (IteM i)

{
  CommonApplyToAllProcMenuItem (i, ADD_TITLE);
}

extern void ApplyCDS (IteM i)

{
  CommonApplyToAllProcMenuItem (i, ADD_CDS);
}

extern void ApplyRRNA (IteM i)

{
  CommonApplyToAllProcMenuItem (i, ADD_RRNA);
}

extern void ApplyImpFeat (IteM i)

{
  CommonApplyToAllProcMenuItem (i, ADD_IMP);
}

#define SUBMISSION_PAGE   0
#define CONTACT_PAGE      1
#define AUTHOR_PAGE       2
#define AFFILIATION_PAGE  3

typedef struct submitform {
  FORM_MESSAGE_BLOCK
  GrouP           pages [4];
  Int2            currentPage;
  DialoG          tbs;

  GrouP           hup;
  GrouP           dateGrp;

  TexT            title;
  DialoG          reldate;
  TexT            firstname;
  TexT            middleinit;
  TexT            lastname;
  PopuP           suffix;
  DialoG          phonefaxemail;
  DialoG          authors;
  ButtoN          use_consortium;
  TexT            consortium;
  DialoG          affil;

  Boolean         visitedContact;
  Boolean         visitedAuthor;

  ButtoN          nextBtn;
  ButtoN          prevBtn;
  BtnActnProc     goToNext;
  BtnActnProc     goToPrev;
} SubmitForm, PNTR SubmitFormPtr;

static AuthListPtr AddConsortiumToAuthList (AuthListPtr alp, TexT consortium)

{
  AuthorPtr    ap;
  ValNodePtr   names;
  PersonIdPtr  pid;

  if (TextHasNoText (consortium)) return alp;
  if (alp == NULL) {
    alp = AuthListNew ();
    alp->choice = 1;
  }
  pid = PersonIdNew ();
  if (pid == NULL) return NULL;
  pid->choice = 5;
  pid->data = SaveStringFromText (consortium);
  ap = AuthorNew ();
  if (ap == NULL) return NULL;
  ap->name = pid;
  names = ValNodeAdd (&(alp->names));
  names->choice = 1;
  names->data.ptrvalue = ap;
  return alp;
}

static Boolean AuthListToConsortium (AuthListPtr alp, TexT consortium)

{
  AuthorPtr    ap;
  ValNodePtr   names;
  PersonIdPtr  pid;
  CharPtr      str;
  Boolean      rval = FALSE;

  if (alp == NULL || consortium == NULL) return FALSE;
  if (alp->choice != 1) return FALSE;
  for (names = alp->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap == NULL) continue;
    pid = ap->name;
    if (pid == NULL || pid->choice != 5 || StringHasNoText (pid->data)) continue;
    str = (CharPtr) pid->data;
    SafeSetTitle (consortium, str);
    rval = TRUE;
  }
  return rval;
}

static void SequinBlockPtrToSubmitForm (ForM f, Pointer data)

{
  AuthorPtr       ap;
  DatePtr         dp;
  NameStdPtr      nsp;
  PersonIdPtr     pid;
  SubmitFormPtr   sbfp;
  SequinBlockPtr  sbp;
  CharPtr         str;
  CharPtr         txt;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  sbp = (SequinBlockPtr) data;
  if (sbfp != NULL) {
    if (sbp != NULL) {
      SafeSetTitle (sbfp->title, sbp->citsubtitle);
      PointerToDialog (sbfp->reldate, (Pointer) sbp->releasedate);
      ap = sbp->contactperson;
      if (ap != NULL) {
        pid = ap->name;
        if (pid != NULL && pid->choice == 2) {
          nsp = pid->data;
          if (nsp != NULL) {
            str = NameStdPtrToAuthorSpreadsheetString (nsp);
            if (str != NULL) {
              txt = ExtractTagListColumn (str, 0);
              SafeSetTitle (sbfp->firstname, txt);
              MemFree (txt);
              txt = ExtractTagListColumn (str, 1);
              SafeSetTitle (sbfp->middleinit, txt);
              MemFree (txt);
              txt = ExtractTagListColumn (str, 2);
              SafeSetTitle (sbfp->lastname, txt);
              MemFree (txt);
              txt = ExtractTagListColumn (str, 3);
              if (! StringHasNoText (txt)) {
                SetEnumPopupByName (sbfp->suffix, name_suffix_alist, txt);
              }
              /*
              SafeSetTitle (sbfp->suffix, txt);
              */
              MemFree (txt);
              MemFree (str);
            }
          }
        }
        PointerToDialog (sbfp->phonefaxemail, (Pointer) ap->affil);
      }
      PointerToDialog (sbfp->authors, (Pointer) sbp->citsubauthors);
      if (AuthListToConsortium (sbp->citsubauthors, sbfp->consortium)) {
        SetStatus (sbfp->use_consortium, TRUE);
        Enable (sbfp->consortium);
      } else {
        SetStatus (sbfp->use_consortium, FALSE);
        Disable (sbfp->consortium);
      }

      PointerToDialog (sbfp->affil, (Pointer) sbp->citsubaffil);
      if (sbp->holduntilpublished) {
        SafeSetValue (sbfp->hup, 2);
        SafeShow (sbfp->dateGrp);
      } else {
        SafeSetValue (sbfp->hup, 1);
        SafeHide (sbfp->dateGrp);
      }
      sbfp->visitedAuthor = TRUE;
      SetValue (sbfp->tbs, 0);
    } else {
      SafeSetTitle (sbfp->title, NULL);
      dp = DateCurr ();
      if (dp != NULL) {
        dp->data [3] = 0; /* force to end of month */
        dp = DateAdvance (dp, 12);
        /*
        (dp->data [1])++;
        if (dp->data [2] == 2 && dp->data [3] > 28) {
          dp->data [3] = 28;
        }
        */
        PointerToDialog (sbfp->reldate, (Pointer) dp);
      } else {
        PointerToDialog (sbfp->reldate, NULL);
      }
      DateFree (dp);
      SafeSetTitle (sbfp->firstname, "");
      SafeSetTitle (sbfp->middleinit, "");
      SafeSetTitle (sbfp->lastname, "");
      PointerToDialog (sbfp->phonefaxemail, NULL);
      PointerToDialog (sbfp->authors, NULL);
      SafeSetTitle (sbfp->consortium, "");
      PointerToDialog (sbfp->affil, NULL);
      SafeSetValue (sbfp->hup, 1);
      SafeHide (sbfp->dateGrp);
      SetValue (sbfp->tbs, 0);
    }
  }
}


static CharPtr ContactNameFromSubmitForm (SubmitFormPtr sbfp)
{
  Char            str [128];
  CharPtr         txt;
  Uint2           suffixVal;
  Char            sfx [32];

  if (sbfp == NULL) {
    return NULL;
  }

  str [0] = '\0';
  txt = SaveStringFromText (sbfp->firstname);
  StringCat (str, txt);
  StringCat (str, "\t");
  MemFree (txt);
  txt = SaveStringFromText (sbfp->middleinit);
  StringCat (str, txt);
  StringCat (str, "\t");
  MemFree (txt);
  txt = SaveStringFromText (sbfp->lastname);
  StringCat (str, txt);
  StringCat (str, "\t");
  MemFree (txt);
  suffixVal = GetValue (sbfp->suffix);
  sprintf (sfx, "%d", (int) (suffixVal - 1));
  StringCat (str, sfx);
  StringCat (str, "\n");
  txt = StringSave (str);
  return txt;
}


static void SetContactNameFromString (SubmitFormPtr sbfp, CharPtr txt)
{
  CharPtr tmp;

  if (sbfp == NULL) {
    return;
  }
  tmp = ExtractTagListColumn (txt, 0);
  SetTitle (sbfp->firstname, tmp);
  tmp = MemFree (tmp);
  tmp = ExtractTagListColumn (txt, 1);
  SetTitle (sbfp->middleinit, tmp);
  tmp = MemFree (tmp);
  tmp = ExtractTagListColumn (txt, 2);
  SetTitle (sbfp->lastname, tmp);
  tmp = MemFree (tmp);
}


static void FixSubmitFormSpecialCharacters (SubmitFormPtr sbfp)
{
  ValNodePtr list = NULL;
  CharPtr         citsubtitle;
  AffilPtr        affil;
  CharPtr         firstname;
  CharPtr         middleinit;
  CharPtr         lastname;

  if (sbfp == NULL) {
    return;
  }

  citsubtitle = SaveStringFromText (sbfp->title);
  SpecialCharFindWithContext (&citsubtitle, &list, NULL, NULL);

  firstname = SaveStringFromText (sbfp->firstname);
  middleinit = SaveStringFromText (sbfp->middleinit);
  lastname = SaveStringFromText (sbfp->lastname);

  SpecialCharFindWithContext (&firstname, &list, NULL, NULL);
  SpecialCharFindWithContext (&middleinit, &list, NULL, NULL);
  SpecialCharFindWithContext (&lastname, &list, NULL, NULL);

  affil = (AffilPtr) DialogToPointer (sbfp->phonefaxemail);
  if (affil != NULL) {
    SpecialCharFindWithContext (&(affil->fax), &list, NULL, NULL);
    SpecialCharFindWithContext (&(affil->phone), &list, NULL, NULL);
    SpecialCharFindWithContext (&(affil->email), &list, NULL, NULL);
  }

  if (list != NULL) {
    FixSpecialCharactersForStringsInList (list, "You must replace all non-ASCII characters.", TRUE);
    SetTitle (sbfp->title, citsubtitle);
    SetTitle (sbfp->firstname, firstname);
    SetTitle (sbfp->middleinit, middleinit);
    SetTitle (sbfp->lastname, lastname);
    PointerToDialog (sbfp->phonefaxemail, affil);
    FreeContextList (list);
  }

  citsubtitle = MemFree (citsubtitle);
  firstname = MemFree (firstname);
  middleinit = MemFree (middleinit);
  lastname = MemFree (lastname);
  affil = AffilFree (affil);
}


static Pointer SubmitFormToSequinBlockPtr (ForM f)

{
  AffilPtr        affil;
  AuthorPtr       ap;
  NameStdPtr      nsp;
  PersonIdPtr     pid;
  SubmitFormPtr   sbfp;
  SequinBlockPtr  sbp;
  CharPtr         txt;

  sbp = NULL;
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    FixSubmitFormSpecialCharacters(sbfp);
    sbp = (SequinBlockPtr) MemNew (sizeof (SequinBlock));
    if (sbp != NULL) {
      sbp->citsubtitle = SaveStringFromTextAndStripNewlines (sbfp->title);
      ap = AuthorNew ();
      if (ap != NULL) {
        pid = PersonIdNew ();
        ap->name = pid;
        if (pid != NULL) {
          pid->choice = 2;
          txt = ContactNameFromSubmitForm (sbfp);
          nsp = AuthorSpreadsheetStringToNameStdPtr (txt);
          MemFree (txt);
          pid->data = nsp;
          if (nsp != NULL) {
            if (StringHasNoText (nsp->names [0])) {
              ap = AuthorFree (ap);
            }
          }
        }
        affil = (AffilPtr) DialogToPointer (sbfp->phonefaxemail);
        if (affil != NULL) {
          if (affil->choice == 2) {
            affil->affil = MemFree (affil->affil);
            affil->div = MemFree (affil->div);
            affil->city = MemFree (affil->city);
            affil->sub = MemFree (affil->sub);
            affil->country = MemFree (affil->country);
            affil->street = MemFree (affil->street);
            affil->postal_code = MemFree (affil->postal_code);
            if (affil->phone == NULL && affil->fax == NULL &&
                affil->email == NULL) {
              affil = AffilFree (affil);
            }
          } else {
            affil = AffilFree (affil);
          }
        }
        if (affil != NULL) {
          if (ap == NULL) {
            ap = AuthorNew();
          }
          ap->affil = affil;
        }
      }
      sbp->contactperson = ap;
      sbp->citsubauthors = (AuthListPtr) DialogToPointer (sbfp->authors);
      if (GetStatus (sbfp->use_consortium)) {
        sbp->citsubauthors = AddConsortiumToAuthList (sbp->citsubauthors, sbfp->consortium);
      }
      sbp->citsubaffil = (AffilPtr) DialogToPointer (sbfp->affil);
      if (GetValue (sbfp->hup) == 2) {
        sbp->holduntilpublished = TRUE;
        sbp->releasedate = (DatePtr) DialogToPointer (sbfp->reldate);
      }
      if (sbp->contactperson == NULL &&
          sbp->citsubauthors == NULL &&
          sbp->citsubaffil == NULL) {
        sbp = SequinBlockFree (sbp);
      }
    }
  }
  return (Pointer) sbp;
}

static void SubmitBlockPtrToSubmitForm (ForM f, Pointer data)

{
  AuthorPtr       ap;
  AuthListPtr     authors;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  DatePtr         dp;
  NameStdPtr      nsp;
  PersonIdPtr     pid;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;
  CharPtr         str;
  CharPtr         txt;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  sbp = (SubmitBlockPtr) data;
  if (sbfp != NULL) {
    if (sbp != NULL) {
      PointerToDialog (sbfp->reldate, (Pointer) sbp->reldate);
      cip = sbp->contact;
      if (cip != NULL) {
        ap = cip->contact;
        if (ap != NULL) {
          pid = ap->name;
          if (pid != NULL && pid->choice == 2) {
            nsp = pid->data;
            if (nsp != NULL) {
              str = NameStdPtrToAuthorSpreadsheetString (nsp);
              if (str != NULL) {
                txt = ExtractTagListColumn (str, 0);
                SafeSetTitle (sbfp->firstname, txt);
                MemFree (txt);
                txt = ExtractTagListColumn (str, 1);
                SafeSetTitle (sbfp->middleinit, txt);
                MemFree (txt);
                txt = ExtractTagListColumn (str, 2);
                SafeSetTitle (sbfp->lastname, txt);
                MemFree (txt);
                txt = ExtractTagListColumn (str, 3);
                if (! StringHasNoText (txt)) {
                  SetEnumPopupByName (sbfp->suffix, name_suffix_alist, txt);
                }
                /*
                SafeSetTitle (sbfp->suffix, txt);
                */
                MemFree (txt);
                MemFree (str);
              }
            }
          }
          PointerToDialog (sbfp->phonefaxemail, (Pointer) ap->affil);
        }
      }
      csp = sbp->cit;
      if (csp != NULL) {
        authors = csp->authors;
        if (authors != NULL) {
          PointerToDialog (sbfp->authors, (Pointer) authors);
          if (AuthListToConsortium (authors, sbfp->consortium)) {
            SetStatus (sbfp->use_consortium, TRUE);
            Enable (sbfp->consortium);
          } else {
            SetStatus (sbfp->use_consortium, FALSE);
            Disable (sbfp->consortium);
          }
          PointerToDialog (sbfp->affil, (Pointer) authors->affil);
        }
      }
      if (sbp->hup) {
        SafeSetValue (sbfp->hup, 2);
        SafeShow (sbfp->dateGrp);
      } else {
        SafeSetValue (sbfp->hup, 1);
        SafeHide (sbfp->dateGrp);
      }
      sbfp->visitedAuthor = TRUE;
      SetValue (sbfp->tbs, 0);
    } else {
      SafeSetTitle (sbfp->title, NULL);
      dp = DateCurr ();
      if (dp != NULL) {
        dp->data [3] = 0; /* force to end of month */
        dp = DateAdvance (dp, 12);
        PointerToDialog (sbfp->reldate, (Pointer) dp);
      } else {
        PointerToDialog (sbfp->reldate, NULL);
      }
      DateFree (dp);
      SafeSetTitle (sbfp->firstname, "");
      SafeSetTitle (sbfp->middleinit, "");
      SafeSetTitle (sbfp->lastname, "");
      SafeSetValue (sbfp->suffix, 0);
      PointerToDialog (sbfp->phonefaxemail, NULL);
      PointerToDialog (sbfp->authors, NULL);
      SafeSetTitle (sbfp->consortium, "");
      PointerToDialog (sbfp->affil, NULL);
      SafeSetValue (sbfp->hup, 1);
      SafeHide (sbfp->dateGrp);
      SetValue (sbfp->tbs, 0);
    }
  }
}

static void TitleToSubmitBlockForm(ForM f, Pointer data)

{
  CharPtr         man_title;
  SubmitFormPtr   sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  man_title = (CharPtr) data;
  if (sbfp != NULL) 
  {
    if (man_title == NULL) {
      SetTitle (sbfp->title, "");
    } else {
      SetTitle (sbfp->title, man_title);
    }
  }
}

static ValNodePtr TestSubmitForm (ForM f)

{
  AffilPtr       affil;
  AuthListPtr    authors;
  DatePtr        dp;
  ValNodePtr     head;
  SubmitFormPtr  sbfp;

  head = NULL;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {

    if (TextHasNoText (sbfp->firstname)) {
      head = AddStringToValNodeChain (head, "first name", 1);
    }
    if (TextHasNoText (sbfp->lastname)) {
      head = AddStringToValNodeChain (head, "last name", 1);
    }
    affil = DialogToPointer (sbfp->phonefaxemail);
    if (affil != NULL) {
      if (StringHasNoText (affil->phone)) {
        head = AddStringToValNodeChain (head, "telephone number", 0);
      }
      if (StringHasNoText (affil->fax)) {
        head = AddStringToValNodeChain (head, "fax number", 0);
      }
      if (StringHasNoText (affil->email)) {
        head = AddStringToValNodeChain (head, "e-mail address", 0);
      }
    } else {
      head = AddStringToValNodeChain (head, "telephone number", 0);
      head = AddStringToValNodeChain (head, "fax number", 0);
      head = AddStringToValNodeChain (head, "e-mail address", 0);
    }
    affil = AffilFree (affil);

    if (GetValue (sbfp->hup) == 2) {
      dp = DialogToPointer (sbfp->reldate);
      if (dp == NULL) {
        head = AddStringToValNodeChain (head, "release date", 1);
      }
      dp = DateFree (dp);
    }
    if (TextHasNoText (sbfp->title)) {
      head = AddStringToValNodeChain (head, "manuscript title", 1);
    }

    authors = DialogToPointer (sbfp->authors);
    if (GetStatus (sbfp->use_consortium)) {
      authors = AddConsortiumToAuthList (authors, sbfp->consortium);
    }
    if (authors == NULL) {
      head = AddStringToValNodeChain (head, "author names", 1);
    }
    authors = AuthListFree (authors);
    affil = DialogToPointer (sbfp->affil);
    if (affil == NULL) {
      head = AddStringToValNodeChain (head, "affiliation", 1);
    }
    affil = AffilFree (affil);

  }

  return head;
}

/*
static Boolean ReadSubmitBlock (ForM f, CharPtr filename)

{
  AsnIoPtr        aip;
  Char            path [PATH_MAX];
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        sbp = SubmitBlockAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (sbp != NULL) {
          SubmitBlockPtrToSubmitForm (f, (Pointer) sbp);
          sbp = SubmitBlockFree (sbp);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}
*/


static Pointer ReadNextASNObject (FILE *fp, Uint2Ptr datatypeptr, Uint2Ptr entityIDptr)
{
  AsnIoPtr       aip;
  Pointer        dataptr = NULL;
  Int4           pos;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp = NULL;
  FileCache      fc;
  CharPtr        str, tag, pEnd;
  Char           line [4096];
    
  if (fp == NULL) {
      return NULL;
  }
  
  FileCacheSetup (&fc, fp);

  pos = FileCacheTell (&fc);
  
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL && StringHasNoText (str)) {
    str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }

  if (str == NULL) return NULL; /* already at end of file */

  if (StringStr (str, "::=") == NULL) return NULL;

  /* first skip past empty space at start of line */
  tag = str;
  while (*tag != '\0' && IS_WHITESP (*tag)) {
    tag++;
  }
  pEnd = tag;
  while (*pEnd != '\0' && !IS_WHITESP (*pEnd)) {
    pEnd++;
  }
  *pEnd = 0;

  omp = ObjMgrReadLock ();
  omtp = ObjMgrTypeFind (omp, 0, tag, NULL);
  ObjMgrUnlock ();

  if (omtp == NULL) {
      return NULL;
  }
  FileCacheFree (&fc, FALSE);
  fseek (fp, pos, SEEK_SET);

  aip = AsnIoNew (ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  aip->scan_for_start = TRUE;
  dataptr = (*(omtp->asnread)) (aip, NULL);
  pos = AsnIoTell (aip);
  AsnIoFree (aip, FALSE);
  fseek (fp, pos, SEEK_SET);

  if (dataptr == NULL) {
    ErrPostEx (SEV_ERROR, 0, 0, "Couldn't read type [%s]", omtp->asnname);
  } else {
    if (datatypeptr != NULL) {
      *datatypeptr = omtp->datatype;
    }
    if (entityIDptr != NULL) {
      *entityIDptr = ObjMgrRegister (omtp->datatype, dataptr);
    }
  }
  return dataptr;
}


static Boolean ReadSubmitBlock (ForM f, CharPtr filename)

{
  Pointer         dataptr;
  Uint2           datatype;
  Uint2           entityID;
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;
  SeqSubmitPtr    ssp;
  Char            path [PATH_MAX];
  CharPtr         man_title;
  CitGenPtr       cgp;
  FILE           * fp;
  PubdescPtr      pdp;
  ValNodePtr      sdp;
  

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
	  fp = FileOpen(path, "r");
      dataptr = ReadNextASNObject (fp, &datatype, &entityID);
      sbp = NULL;
      man_title = NULL;
      while (dataptr != NULL) {
        if (entityID > 0) {
          switch (datatype) {
            case OBJ_SUBMIT_BLOCK :
              if (sbp == NULL) {
                sbp = (SubmitBlockPtr) AsnIoMemCopy (dataptr,
                                                     (AsnReadFunc) SubmitBlockAsnRead,
                                                     (AsnWriteFunc) SubmitBlockAsnWrite);
              }
              break;
            case OBJ_SEQSUB :
              if (sbp == NULL) {
                ssp = (SeqSubmitPtr) dataptr;
                if (ssp != NULL) {
                  sbp = (SubmitBlockPtr) AsnIoMemCopy (ssp->sub,
                                                       (AsnReadFunc) SubmitBlockAsnRead,
                                                       (AsnWriteFunc) SubmitBlockAsnWrite);
                }
              }
              break;
            case OBJ_SEQDESC:
              sdp = (ValNodePtr) dataptr;
              if (sdp->choice == Seq_descr_pub && sdp->data.ptrvalue != NULL && man_title == NULL) {
                pdp = (PubdescPtr) sdp->data.ptrvalue;
                sdp = pdp->pub;
                while (sdp != NULL && sdp->choice != PUB_Gen) {
                  sdp = sdp->next;
                }
                if (sdp != NULL && sdp->data.ptrvalue != NULL) {
                  cgp = (CitGenPtr) sdp->data.ptrvalue;
                  if (StringICmp (cgp->cit, "unpublished") == 0 
                      && !StringHasNoText (cgp->title)) {
                    man_title = StringSave (cgp->title);
                  }
                }
              }
              break;
            default :
              break;
          }
        }
        ObjMgrDelete (datatype, dataptr);
        dataptr = ReadNextASNObject (fp, &datatype, &entityID);
      }
      FileClose (fp);
      if (sbp != NULL || man_title != NULL) {
        if (sbp != NULL) {
          SubmitBlockPtrToSubmitForm (f, sbp);
          sbp = SubmitBlockFree (sbp);
        }
        if (man_title != NULL) {
          TitleToSubmitBlockForm (f, man_title);
          man_title = MemFree (man_title);
        }
        Update ();
        return TRUE;
      }
    }
  }
  return FALSE;
}

static SeqDescrPtr MakeUnpubPub (CharPtr title, AuthListPtr alp, AffilPtr affil)

{
  CitGenPtr    cgp;
  PubdescPtr   pdp;
  ValNodePtr   pep;
  SeqDescrPtr  vnp;

  if (StringHasNoText (title) || alp == NULL) return NULL;
  pdp = PubdescNew ();
  if (pdp == NULL) return NULL;
  vnp = SeqDescrNew (NULL);
  if (vnp == NULL) return NULL;
  vnp->choice = Seq_descr_pub;
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
      cgp->authors = alp;
      if (alp != NULL) {
        alp->affil = affil;
        if (affil != NULL) {
          affil->phone = MemFree (affil->phone);
          affil->fax = MemFree (affil->fax);
          affil->email = MemFree (affil->email);
        }
      }
      cgp->title = StringSave (title);
    }
  }
  return vnp;
}

extern SubmitBlockPtr ConvertSequinBlockToSubmitBlock (SequinBlockPtr sqp);

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));

static Boolean WriteSubmitBlock (ForM f, CharPtr filename)

{
  AffilPtr        affil;
  AsnIoPtr        aip;
  AuthListPtr     alp;
  Char            path [PATH_MAX];
  SubmitFormPtr   sbfp;
  SubmitBlockPtr  sbp;
  SeqDescrPtr     sdp;
  SequinBlockPtr  sqp;
  CharPtr         title;
#ifdef WIN_MAC
  FILE            *fp;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      fp = FileOpen (path, "r");
      if (fp != NULL) {
        FileClose (fp);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      sqp = (SequinBlockPtr) SubmitFormToSequinBlockPtr (f);
      if (sqp != NULL) {
        title = StringSaveNoNull (sqp->citsubtitle);
        alp = AsnIoMemCopy ((Pointer) sqp->citsubauthors,
                            (AsnReadFunc) AuthListAsnRead,
                            (AsnWriteFunc) AuthListAsnWrite);
        affil = AsnIoMemCopy ((Pointer) sqp->citsubaffil,
                              (AsnReadFunc) AffilAsnRead,
                              (AsnWriteFunc) AffilAsnWrite);
        sbp = ConvertSequinBlockToSubmitBlock (sqp);
        if (sbp != NULL) {
          aip = AsnIoOpen (path, "w");
          if (aip != NULL) {
            sbp->tool = MemFree (sbp->tool);
            SubmitBlockAsnWrite (sbp, aip, NULL);
            AsnPrintNewLine (aip);
            sdp = MakeUnpubPub (title, alp, affil);
            if (sdp != NULL) {
              AsnIoReset (aip);
              SeqDescAsnWrite (sdp, aip, NULL);
              AsnPrintNewLine (aip);
              SeqDescFree (sdp);
            }
            AsnIoClose (aip);
            sbp = SubmitBlockFree (sbp);
            return TRUE;
          }
        }
      }
    }
  }
  return FALSE;
}

static Boolean ReadContactPage (ForM f, CharPtr filename)

{
  AsnIoPtr        aip;
  AuthorPtr       ap;
  ContactInfoPtr  cip;
  NameStdPtr      nsp;
  Char            path [PATH_MAX];
  PersonIdPtr     pid;
  SubmitFormPtr   sbfp;
  CharPtr         str;
  CharPtr         txt;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        cip = ContactInfoAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (cip != NULL) {
          ap = cip->contact;
          if (ap != NULL) {
            pid = ap->name;
            if (pid != NULL && pid->choice == 2) {
              nsp = pid->data;
              if (nsp != NULL) {
                str = NameStdPtrToAuthorSpreadsheetString (nsp);
                if (str != NULL) {
                  txt = ExtractTagListColumn (str, 0);
                  SafeSetTitle (sbfp->firstname, txt);
                  MemFree (txt);
                  txt = ExtractTagListColumn (str, 1);
                  SafeSetTitle (sbfp->middleinit, txt);
                  MemFree (txt);
                  txt = ExtractTagListColumn (str, 2);
                  SafeSetTitle (sbfp->lastname, txt);
                  MemFree (txt);
                  txt = ExtractTagListColumn (str, 3);
                  if (! StringHasNoText (txt)) {
                    SetEnumPopupByName (sbfp->suffix, name_suffix_alist, txt);
                  }
                  MemFree (txt);
                  MemFree (str);
                }
              }
            }
            PointerToDialog (sbfp->phonefaxemail, (Pointer) ap->affil);
          }
          cip = ContactInfoFree (cip);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteContactPage (ForM f, CharPtr filename)

{
  AffilPtr        affil;
  AsnIoPtr        aip;
  AuthorPtr       ap;
  ContactInfoPtr  cip;
  NameStdPtr      nsp;
  Char            path [PATH_MAX];
  PersonIdPtr     pid;
  SubmitFormPtr   sbfp;
  Char            sfx [32];
  Char            str [128];
  Uint2           suffixVal;
  CharPtr         txt;
#ifdef WIN_MAC
  FILE            *fp;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      fp = FileOpen (path, "r");
      if (fp != NULL) {
        FileClose (fp);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        cip = ContactInfoNew ();
        if (cip != NULL) {
          ap = AuthorNew ();
          if (ap != NULL) {
            pid = PersonIdNew ();
            ap->name = pid;
            if (pid != NULL) {
              pid->choice = 2;
              str [0] = '\0';
              txt = SaveStringFromText (sbfp->firstname);
              StringCat (str, txt);
              StringCat (str, "\t");
              MemFree (txt);
              txt = SaveStringFromText (sbfp->middleinit);
              StringCat (str, txt);
              StringCat (str, "\t");
              MemFree (txt);
              txt = SaveStringFromText (sbfp->lastname);
              StringCat (str, txt);
              StringCat (str, "\t");
              MemFree (txt);
              suffixVal = GetValue (sbfp->suffix);
              sprintf (sfx, "%d", (int) (suffixVal - 1));
              StringCat (str, sfx);
              StringCat (str, "\n");
              txt = StringSave (str);
              nsp = AuthorSpreadsheetStringToNameStdPtr (txt);
              MemFree (txt);
              pid->data = nsp;
              if (nsp != NULL) {
                if (StringHasNoText (nsp->names [0])) {
                  ap = AuthorFree (ap);
                }
              }
            }
            affil = (AffilPtr) DialogToPointer (sbfp->phonefaxemail);
            if (affil != NULL) {
              if (affil->choice == 2) {
                affil->affil = MemFree (affil->affil);
                affil->div = MemFree (affil->div);
                affil->city = MemFree (affil->city);
                affil->sub = MemFree (affil->sub);
                affil->country = MemFree (affil->country);
                affil->street = MemFree (affil->street);
                affil->postal_code = MemFree (affil->postal_code);
                if (affil->phone == NULL && affil->fax == NULL &&
                    affil->email == NULL) {
                  affil = AffilFree (affil);
                }
              } else {
                affil = AffilFree (affil);
              }
            }
            ap->affil = affil;
          }
          cip->contact = ap;
          ContactInfoAsnWrite (cip, aip, NULL);
          AsnIoClose (aip);
          cip = ContactInfoFree (cip);
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean ReadAuthListPage (ForM f, CharPtr filename)

{
  AsnIoPtr       aip;
  AuthListPtr    alp;
  Char           path [PATH_MAX];
  SubmitFormPtr  sbfp;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        alp = AuthListAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (alp != NULL) {
          PointerToDialog (sbfp->authors, alp);
          if (AuthListToConsortium (alp, sbfp->consortium)) {
            SetStatus (sbfp->use_consortium, TRUE);
            Enable (sbfp->consortium);
          } else {
            SetStatus (sbfp->use_consortium, FALSE);
            Disable (sbfp->consortium);
          }
          alp = AuthListFree (alp);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteAuthListPage (ForM f, CharPtr filename)

{
  AsnIoPtr       aip;
  AuthListPtr    alp;
  Char           path [PATH_MAX];
  SubmitFormPtr  sbfp;
#ifdef WIN_MAC
  FILE           *fp;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      fp = FileOpen (path, "r");
      if (fp != NULL) {
        FileClose (fp);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        alp = DialogToPointer (sbfp->authors);
        if (GetStatus (sbfp->use_consortium)) {
          alp = AddConsortiumToAuthList (alp, sbfp->consortium);
        }
        AuthListAsnWrite (alp, aip, NULL);
        AsnIoClose (aip);
        alp = AuthListFree (alp);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean ReadAffilPage (ForM f, CharPtr filename)

{
  AffilPtr       affil;
  AsnIoPtr       aip;
  Char           path [PATH_MAX];
  SubmitFormPtr  sbfp;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
      aip = AsnIoOpen (path, "r");
      if (aip != NULL) {
        affil = AffilAsnRead (aip, NULL);
        AsnIoClose (aip);
        if (affil != NULL) {
          affil->phone = MemFree (affil->phone);
          affil->fax = MemFree (affil->fax);
          affil->email = MemFree (affil->email);
          PointerToDialog (sbfp->affil, affil);
          affil = AffilFree (affil);
          Update ();
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean WriteAffilPage (ForM f, CharPtr filename)

{
  AffilPtr       affil;
  AsnIoPtr       aip;
  Char           path [PATH_MAX];
  SubmitFormPtr  sbfp;
#ifdef WIN_MAC
  FILE           *fp;
#endif

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      fp = FileOpen (path, "r");
      if (fp != NULL) {
        FileClose (fp);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        affil = DialogToPointer (sbfp->affil);
        affil->phone = MemFree (affil->phone);
        affil->fax = MemFree (affil->fax);
        affil->email = MemFree (affil->email);
        AffilAsnWrite (affil, aip, NULL);
        AsnIoClose (aip);
        affil = AffilFree (affil);
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean ImportSubmitForm (ForM f, CharPtr filename)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    switch (sbfp->currentPage) {
      case SUBMISSION_PAGE :
        return ReadSubmitBlock (f, filename);
      case CONTACT_PAGE :
        return ReadContactPage (f, filename);
      case AUTHOR_PAGE :
        return ReadAuthListPage (f, filename);
      case AFFILIATION_PAGE :
        return ReadAffilPage (f, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportSubmitForm (ForM f, CharPtr filename)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    switch (sbfp->currentPage) {
      case SUBMISSION_PAGE :
        return WriteSubmitBlock (f, filename);
      case CONTACT_PAGE :
        return WriteContactPage (f, filename);
      case AUTHOR_PAGE :
        return WriteAuthListPage (f, filename);
      case AFFILIATION_PAGE :
        return WriteAffilPage (f, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static void CopyContactToAuthors (SubmitFormPtr sbfp)

{
  AuthListPtr  alp;
  AuthorPtr    ap;
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      txt;
  ValNodePtr   list = NULL;

  if (sbfp == NULL) return;
  ap = NULL;
  alp = AuthListNew ();
  if (alp != NULL) {
    alp->choice = 1;
    names = ValNodeNew (NULL);
    alp->choice = 1;
    alp->names = names;
    if (names != NULL) {
      ap = AuthorNew ();
      if (ap != NULL) {
        pid = PersonIdNew ();
        ap->name = pid;
        if (pid != NULL) {
          pid->choice = 2;
          FixSubmitFormSpecialCharacters (sbfp);
          txt = ContactNameFromSubmitForm (sbfp);
          nsp = AuthorSpreadsheetStringToNameStdPtr (txt);
          MemFree (txt);
          pid->data = nsp;
          if (nsp != NULL) {
            if (StringHasNoText (nsp->names [0])) {
              ap = AuthorFree (ap);
            }
          }
        }
      }
      names->choice = 1;
      names->data.ptrvalue = ap;
    }
    if (ap == NULL) {
      alp = AuthListFree (alp);
    }
    if (alp != NULL) {
      PointerToDialog (sbfp->authors, (Pointer) alp);
      if (AuthListToConsortium (alp, sbfp->consortium)) {
        SetStatus (sbfp->use_consortium, TRUE);
        Enable (sbfp->consortium);
      } else {
        SetStatus (sbfp->use_consortium, FALSE);
        Disable (sbfp->consortium);
      }
    }
  }
  alp = AuthListFree (alp);
}

static void SetSubmitterImportExportItems (SubmitFormPtr sbfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (sbfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) sbfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) sbfp, VIB_MSG_EXPORT);
    switch (sbfp->currentPage) {
      case SUBMISSION_PAGE :
        SafeSetTitle (importItm, "Import Submitter Info...");
        SafeSetTitle (exportItm, "Export Submitter Info...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      case CONTACT_PAGE :
        SafeSetTitle (importItm, "Import Contact...");
        SafeSetTitle (exportItm, "Export Contact...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      case AUTHOR_PAGE :
        SafeSetTitle (importItm, "Import Authors...");
        SafeSetTitle (exportItm, "Export Authors...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      case AFFILIATION_PAGE :
        SafeSetTitle (importItm, "Import Affiliation...");
        SafeSetTitle (exportItm, "Export Affiliation...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void EnterSubmitPage (SubmitFormPtr sbfp, Int2 page)

{
  AuthListPtr  alp;

  if (sbfp != NULL) {
    switch (page) {
      case SUBMISSION_PAGE :
        Select (sbfp->title);
        SafeSetTitle (sbfp->prevBtn, "<< Prev Form");
        SafeSetTitle (sbfp->nextBtn, "Next Page >>");
        break;
      case CONTACT_PAGE :
        Select (sbfp->firstname);
        sbfp->visitedContact = TRUE;
        SafeSetTitle (sbfp->prevBtn, "<< Prev Page");
        SafeSetTitle (sbfp->nextBtn, "Next Page >>");
        break;
      case AUTHOR_PAGE :
        alp = (AuthListPtr) DialogToPointer (sbfp->authors);
        if (GetStatus (sbfp->use_consortium)) {
          alp = AddConsortiumToAuthList (alp, sbfp->consortium);
        }
        if (sbfp->visitedContact && alp == NULL) {
          CopyContactToAuthors (sbfp);
        }
        AuthListFree (alp);
        /*
        if (sbfp->visitedContact && (! sbfp->visitedAuthor)) {
          CopyContactToAuthors (sbfp);
        }
        */
        SendMessageToDialog (sbfp->authors, VIB_MSG_ENTER);
        sbfp->visitedAuthor = TRUE;
        SafeSetTitle (sbfp->prevBtn, "<< Prev Page");
        SafeSetTitle (sbfp->nextBtn, "Next Page >>");
        break;
      case AFFILIATION_PAGE :
        SendMessageToDialog (sbfp->affil, VIB_MSG_ENTER);
        SafeSetTitle (sbfp->prevBtn, "<< Prev Page");
        SafeSetTitle (sbfp->nextBtn, "Next Form >>");
        break;
      default :
        break;
    }
  }
}

static void ChangeSubmitFormPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) data;
  if (sbfp != NULL) {
    sbfp->currentPage = newval;
    SafeHide (sbfp->pages [oldval]);
    Update ();
#ifdef WIN_MAC
    EnterSubmitPage (sbfp, newval);
#endif
    SetSubmitterImportExportItems (sbfp);
    SafeShow (sbfp->pages [newval]);
#ifndef WIN_MAC
    EnterSubmitPage (sbfp, newval);
#endif
    Update ();
    switch (newval) {
      case SUBMISSION_PAGE :
        SendHelpScrollMessage (helpForm, "Submitting Authors Form", "Submission Page");
        break;
      case CONTACT_PAGE :
        SendHelpScrollMessage (helpForm, "Submitting Authors Form", "Contact Page");
        break;
      case AUTHOR_PAGE :
        SendHelpScrollMessage (helpForm, "Submitting Authors Form", "Authors Page");
        break;
      case AFFILIATION_PAGE :
        SendHelpScrollMessage (helpForm, "Submitting Authors Form", "Affiliation Page");
        break;
      default :
        break;
    }
  }
}

static void NextSubmitFormBtn (ButtoN b)

{
  SubmitFormPtr  sbfp;
  CharPtr        txt;
  Boolean        ok = TRUE;

  sbfp = (SubmitFormPtr) GetObjectExtra (b);
  if (sbfp != NULL) {
    if (sbfp->currentPage == 0)
    {
      txt = SaveStringFromText (sbfp->title);
      if (StringHasNoText (txt))
      {
      	Message (MSG_ERROR, "You must supply a tentative title for your manuscript.");
      	MemFree (txt);
      	return;
      }
      MemFree (txt);
    }
    else if (sbfp->currentPage == 1)
    {
      txt = SaveStringFromText (sbfp->firstname);
      if (StringHasNoText (txt))
      {
      	ok = FALSE;
      }
      MemFree (txt);
      txt = SaveStringFromText (sbfp->lastname);
      if (StringHasNoText (txt))
      {
      	ok = FALSE;
      }
      MemFree (txt);
      if (!ok)
      {
      	Message (MSG_ERROR, "You must supply a first and last name for the contact.");
      	return;
      }
    }
    if (sbfp->currentPage < 3) {
      SetValue (sbfp->tbs, sbfp->currentPage + 1);
    } else if (sbfp->goToNext != NULL) {
      (sbfp->goToNext) (b);
    }
  }
}

static void PrevSubmitFormBtn (ButtoN b)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (b);
  if (sbfp != NULL) {
    if (sbfp->currentPage > 0) {
      SetValue (sbfp->tbs, sbfp->currentPage - 1);
    } else if (sbfp->goToPrev != NULL) {
      (sbfp->goToPrev) (b);
    }
  }
}

static void ChangeHup (GrouP g)

{
  Boolean        hup;
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (g);
  if (sbfp != NULL) {
    hup = (Boolean) (GetValue (sbfp->hup) == 2);
    if (hup) {
      SafeShow (sbfp->dateGrp);
    } else {
      SafeHide (sbfp->dateGrp);
    }
  }
}

static CharPtr  submitFormTabs [] = {
  "Submission", "Contact", "Authors", "Affiliation", NULL
};

static void SubmitFormMessage (ForM f, Int2 mssg)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (f);
  if (sbfp != NULL) {
    switch (mssg) {
      case VIB_MSG_IMPORT :
        ImportSubmitForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportSubmitForm (f, NULL);
        break;
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
        if (sbfp->appmessage != NULL) {
          sbfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void InitSubmitterFormActivate (WindoW w)

{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (w);
  if (sbfp != NULL) {
    if (sbfp->activate != NULL) {
      sbfp->activate (w);
    }
    SetSubmitterImportExportItems (sbfp);
  }
}

static void EnableConsortium (ButtoN b)
{
  SubmitFormPtr  sbfp;

  sbfp = (SubmitFormPtr) GetObjectExtra (b);
  if (sbfp != NULL) {
    if (GetStatus (sbfp->use_consortium)) {
      Enable (sbfp->consortium);
    } else {
      Disable (sbfp->consortium);
    }
  }
}

extern ForM CreateInitSubmitterForm (Int2 left, Int2 top, CharPtr title,
                                     BtnActnProc goToNext,
                                     BtnActnProc goBack,
                                     WndActnProc activateForm)

{
  GrouP              c;
  GrouP              g;
  GrouP              h;
  GrouP              j;
  GrouP              m;
  GrouP              n;
  GrouP              q;
  SubmitFormPtr      sbfp;
  StdEditorProcsPtr  sepp;
  GrouP              t;
  WindoW             w;
  GrouP              x;
  GrouP              z;
  GrouP              g1, g2;
  GrouP              p;
  DatePtr            dp;

  w = NULL;
  sbfp = MemNew (sizeof (SubmitForm));
  if (sbfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, sbfp, StdCleanupFormProc);
    sbfp->form = (ForM) w;
    sbfp->toform = SequinBlockPtrToSubmitForm;
    sbfp->fromform = SubmitFormToSequinBlockPtr;
    sbfp->testform = TestSubmitForm;
    sbfp->importform = ImportSubmitForm;
    sbfp->exportform = ExportSubmitForm;
    sbfp->formmessage = SubmitFormMessage;

#ifndef WIN_MAC
    CreateSqnInitialFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      sbfp->appmessage = sepp->handleMessages;
    }

    SetGroupSpacing (w, 10, 10);

    j = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (j, 10, 10);

    sbfp->tbs = CreateFolderTabs (j, submitFormTabs, 0, 0, 0,
                                  SYSTEM_FOLDER_TAB,
                                  ChangeSubmitFormPage, (Pointer) sbfp);
    sbfp->currentPage = SUBMISSION_PAGE;

    sbfp->visitedContact = FALSE;
    sbfp->visitedAuthor = FALSE;

    h = HiddenGroup (w, 0, 0, NULL);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    m = HiddenGroup (q, -1, 0, NULL);
    g = HiddenGroup (m, 0, -4, NULL);
    SetGroupSpacing (g, 3, 10);
    StaticPrompt (g, "When may we release your sequence record?",
                  0, stdLineHeight, programFont, 'l');
    sbfp->hup = HiddenGroup (g, 0, -2, ChangeHup);
    SetObjectExtra (sbfp->hup, sbfp, NULL);
    RadioButton (sbfp->hup, "Immediately After Processing");
    RadioButton (sbfp->hup, "Release Date:");
    SetValue (sbfp->hup, 1);
    sbfp->dateGrp = HiddenGroup (m, -1, 0, NULL);
    /* StaticPrompt (sbfp->dateGrp, "Release Date: ", 0, popupMenuHeight, programFont, 'l'); */

    dp = DateCurr ();
    if (dp != NULL && dp->data[0] == 1) {
      sbfp->reldate = CreateDateDialogEx (sbfp->dateGrp, NULL, dp->data[1] + 1900, 10);
    } else {
      sbfp->reldate = CreateDateDialog (sbfp->dateGrp, NULL);
    }
    p = MultiLinePrompt (sbfp->dateGrp, 
                         "NOTE: Sequences must be released when the accession number or any portion of the sequence is published.",
                         25 * stdCharWidth, programFont);
    AlignObjects (ALIGN_CENTER, (HANDLE) sbfp->reldate, (HANDLE) p, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) sbfp->dateGrp, NULL);
    Hide (sbfp->dateGrp);
    t = HiddenGroup (q, 1, 0, NULL);
    StaticPrompt (t, "Tentative title for manuscript (required)", 0, 0, programFont, 'c');
    sbfp->title = ScrollText (t, 25, 4, programFont, TRUE, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) sbfp->reldate,
                  (HANDLE) sbfp->title, (HANDLE) t, NULL);
    sbfp->pages [SUBMISSION_PAGE] = q;
    Hide (sbfp->pages [SUBMISSION_PAGE]);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    n = HiddenGroup (q, 4, 0, NULL);
    SetGroupSpacing (n, -1, 2);
    StaticPrompt (n, "First Name", 0, 0, programFont, 'c');
    StaticPrompt (n, "M.I.", 0, 0, programFont, 'c');
    StaticPrompt (n, "Last Name", 0, 0, programFont, 'c');
    StaticPrompt (n, "Sfx", 0, 0, programFont, 'c');
    sbfp->firstname = DialogText (n, "", 8, NULL);
    sbfp->middleinit = DialogText (n, "", 4, NULL);
    sbfp->lastname = DialogText (n, "", 9, NULL);
    /*
    sbfp->suffix = DialogText (n, "", 3, NULL);
    sbfp->suffix = PopupList (n, TRUE, SuffixPopup_Callback);
    */
    sbfp->suffix = PopupList (n, TRUE, NULL);
    SetObjectExtra (sbfp->suffix, sbfp, NULL);
    InitEnumPopup (sbfp->suffix, name_suffix_alist, NULL);
    SetEnumPopup (sbfp->suffix, name_suffix_alist, 0);

    sbfp->phonefaxemail = CreateExtAffilDialog (q, NULL, NULL, &x);
    Show (x);
    Show (sbfp->phonefaxemail);
    AlignObjects (ALIGN_CENTER, (HANDLE) n, (HANDLE) sbfp->phonefaxemail, NULL);
    sbfp->pages [CONTACT_PAGE] = q;
    Hide (sbfp->pages [CONTACT_PAGE]);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    sbfp->authors = CreateAuthorDialog (q, 3, -1);
    z = HiddenGroup (q, 0, 2, NULL);
    g1 = HiddenGroup (z, 2, 0, NULL);
    sbfp->use_consortium = CheckBox (g1, "Consortium", EnableConsortium);
    SetObjectExtra (sbfp->use_consortium, sbfp, NULL);
    sbfp->consortium = DialogText (g1, "", 16, NULL);
    Disable (sbfp->consortium);
    g2 = HiddenGroup (z, 1, 0, NULL);
    MultiLinePrompt (g2, "The consortium field should be used when a "
                         "consortium is responsible for the sequencing or "
                         "publication of the data.  Individual authors may "
                         "be listed along with a consortium name.",
                          25 * stdCharWidth, programFont);
    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) sbfp->authors, (HANDLE) z, NULL);
    sbfp->pages [AUTHOR_PAGE] = q;
    Hide (sbfp->pages [AUTHOR_PAGE]);

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    sbfp->affil = CreateExtAffilDialog (q, NULL, &g, NULL);
    Show (g);
    Show (sbfp->affil);
    sbfp->pages [AFFILIATION_PAGE] = q;
    Hide (sbfp->pages [AFFILIATION_PAGE]);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    sbfp->goToPrev = goBack;
    sbfp->prevBtn = PushButton (c, " << Prev Form ", PrevSubmitFormBtn);
    SetObjectExtra (sbfp->prevBtn, sbfp, NULL);
    sbfp->goToNext = goToNext;
    sbfp->nextBtn = PushButton (c, " Next Page >> ", NextSubmitFormBtn);
    SetObjectExtra (sbfp->nextBtn, sbfp, NULL);

    AlignObjects (ALIGN_CENTER,
                  (HANDLE) sbfp->pages [SUBMISSION_PAGE],
                  (HANDLE) sbfp->pages [CONTACT_PAGE],
                  (HANDLE) sbfp->pages [AUTHOR_PAGE],
                  (HANDLE) sbfp->pages [AFFILIATION_PAGE],
                  (HANDLE) sbfp->tbs, (HANDLE) c, NULL);

    RealizeWindow (w);

    SafeSetTitle (sbfp->prevBtn, "<< Prev Form");
    SafeSetTitle (sbfp->nextBtn, "Next Page >>");

    sbfp->activate = activateForm;
    SetActivate (w, InitSubmitterFormActivate);

    Show (sbfp->pages [sbfp->currentPage]);
    EnterSubmitPage (sbfp, sbfp->currentPage);
  }
  return (ForM) w;
}

extern void AddCodonListTotRNA (tRNAPtr trna, ValNodePtr codons)
{
  ValNodePtr vnp;
  Int2       j, k, q;
  Char       str [8];
  Char       ch;
  Uint1      codon [4]; 
  Uint1      code;

  if (trna == NULL) return;

  for (j = 0; j < 6; j++) {
    trna->codon [j] = 255;
  }
  
  if (codons == NULL)
  {
    return;
  }
  
  for (vnp = codons, j = 0; vnp != NULL && j < 6; vnp = vnp->next) {
    str [0] = '\0';
    StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
    if (str [0] != '\0') {
      k = 0;
      q = 0;
      ch = str [k];
      while (ch != '\0' && q < 3) {
        ch = TO_UPPER (ch);
        if (StringChr ("ACGTU", ch) != NULL) {
          if (ch == 'U') {
            ch = 'T';
          }
          codon [q] = (Uint1) ch;
          q++;
        }
        k++;
        ch = str [k];
      }
      codon [q] = 0;
      if (q == 3) {
        code = IndexForCodon (codon, Seq_code_iupacna);
        if (code != INVALID_RESIDUE) {
          trna->codon [j] = code;
        }
      }
      j++;
    }
  }
}


static void ParseCodonsFromtRNACommentProc (SeqFeatPtr sfp, Pointer userdata)
{
  RnaRefPtr rrp;
  tRNAPtr   trna;
  CharPtr   cp, sep;
  CharPtr   codon_name;
  Char      codon[4];
  Int4      k, q;
  Char      ch;
  Uint1     code;

  if (sfp == NULL 
      || sfp->comment == NULL
      || sfp->data.choice != SEQFEAT_RNA
      || (rrp = (RnaRefPtr)sfp->data.value.ptrvalue) == NULL
      || rrp->type != 3 || rrp->ext.choice != 2) {
    return;
  }
  trna = rrp->ext.value.ptrvalue;
  if (trna == NULL) return;

  cp = StringStr (sfp->comment, "recognized codon");
  if (cp == NULL) {
    cp = StringStr (sfp->comment, "codon recognized");
  }
  if (cp == NULL) return;
  sep = StringStr (cp, "=");
  if (sep == NULL) {
    sep = StringStr (cp, ":");
  }
  cp = sep;
  if (cp == NULL) {
    Message (MSG_ERROR, "Unable to read codon from %s", sfp->comment);
    return;
  }

  codon_name = cp + 1;

  q = 0;
  for (k = 0; q < 3 && codon_name[k] != 0; k++) {
    ch = TO_UPPER (codon_name [k]);
    if (StringChr ("ACGTU", ch) != NULL) {
      if (ch == 'U') {
        ch = 'T';
      }
      codon [q] = (Uint1) ch;
      q++;
    }
  }
  codon [q] = 0;
  if (q == 3) {
    code = IndexForCodon ((Uint1Ptr) codon, Seq_code_iupacna);
    if (code != INVALID_RESIDUE) {
      trna->codon[0] = code;
      for (k = 1; k < 6; k++) {
        trna->codon[k] = 255;
      }
    } else {
      Message (MSG_ERROR, "Invalid codon in %s", sfp->comment);
    }
  } else {
    Message (MSG_ERROR, "Invalid codon in %s", sfp->comment);
  }
}

extern void ParseCodonsFromtRNAComment (IteM i)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, ParseCodonsFromtRNACommentProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean GetAntiCodonPositionFromComment (CharPtr comment, Int4Ptr start, Int4Ptr stop)
{
  CharPtr cp, cpend;
  Int4    val;
  Char    ch;
  
  if (comment == NULL || start == NULL || stop == NULL) return FALSE;
  cp = StringStr (comment, "anticodon");
  if (cp == NULL) return FALSE;
  cp += 8;
  while (*cp != 0 && !isdigit ((Int4)(*cp)))
  {
  	cp++;
  }
  if (*cp == 0)
  {
  	Message (MSG_ERROR, "Found 'anticodon' but no numbers for position in %s", comment);
  	return FALSE;
  }
  
  cpend = cp + 1;
  while (*cpend != 0 && isdigit ((Int4)(*cpend)))
  {
  	cpend ++;
  }
  if (*cpend == 0)
  {
  	Message (MSG_ERROR, "Found 'anticodon' but no numbers for position in %s", comment);
  	return FALSE;
  }
  ch = *cpend;
  *cpend = 0;
  val = atoi (cp);
  *cpend = ch;
  cp = cpend + 1;
  
  while (*cp != 0 && !isdigit ((Int4)(*cp)))
  {
  	cp++;
  }
  if (*cp == 0)
  {
  	Message (MSG_ERROR, "Found 'anticodon' but no numbers for position in %s", comment);
  	return FALSE;
  }
  cpend = cp + 1;
  while (*cpend != 0 && isdigit ((Int4)(*cpend)))
  {
  	cpend ++;
  }
  ch = *cpend;
  *cpend = 0;
  *start = val;
  *stop = atoi (cp);
  *cpend = ch;
  return TRUE;
}

static void ParseAntiCodonsFromtRNACommentProc (SeqFeatPtr sfp, Pointer userdata)
{
  RnaRefPtr rrp;
  tRNAPtr   trna;
  Int4      start;
  Int4      stop;
  SeqLocPtr slp;
  SeqIntPtr sip;

  if (sfp == NULL 
      || sfp->comment == NULL
      || sfp->data.choice != SEQFEAT_RNA
      || (rrp = (RnaRefPtr)sfp->data.value.ptrvalue) == NULL
      || rrp->type != 3 || rrp->ext.choice != 2) {
    return;
  }
  trna = rrp->ext.value.ptrvalue;
  if (trna == NULL) return;

  if (! GetAntiCodonPositionFromComment (sfp->comment, &start, &stop)) return;
  sip = SeqIntNew ();
  if (sip == NULL) return;
  if (start <= stop) {
    sip->from = start - 1;
    sip->to = stop - 1;
  } else {
    sip->from = stop - 1;
    sip->to = start - 1;
  }
  sip->strand = SeqLocStrand (sfp->location);
  sip->id = SeqIdDup (SeqLocId (sfp->location));
  slp = ValNodeNew (NULL);
  if (slp == NULL) return;
  slp->choice = 4;
  slp->data.ptrvalue = sip;

  if (trna->anticodon != NULL) {
    SeqLocFree (trna->anticodon);
  }
  trna->anticodon = slp;
}

extern void ParseAntiCodonsFromtRNAComment (IteM i)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, ParseAntiCodonsFromtRNACommentProc);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct orgnamechangedata {
  CharPtr oldname;
  CharPtr newname;
} OrgNameChangeData, PNTR OrgNameChangePtr;

static void ReplaceOldOrgNameInTitle (SeqDescrPtr sdp, OrgNameChangePtr oncp)
{
  CharPtr oldtitle, newtitle, orgnamestart;
  Int4    oldorglen;

  if (sdp == NULL || oncp == NULL
    || oncp->oldname == NULL || oncp->newname == NULL)
  {
    return;
  }

  oldtitle = sdp->data.ptrvalue;

  if (oldtitle == NULL) return;

  orgnamestart = StringRChr (oldtitle, '[');
  if (orgnamestart == NULL) return;
  oldorglen = StringLen (oncp->oldname);
  if (StringNCmp (orgnamestart + 1, oncp->oldname, oldorglen) == 0
    && orgnamestart [ oldorglen + 1] == ']')
  {
    newtitle = MemNew (StringLen (oldtitle) - StringLen (oncp->oldname) + StringLen (oncp->newname) + 2);   
    if (newtitle == NULL) return;
    StringNCpy (newtitle, oldtitle, orgnamestart - oldtitle + 1);
    StringCat (newtitle, oncp->newname);
    StringCat (newtitle, "]");
    MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = newtitle;
  }
}

static void FixProteinTitleAfterOrganismNameChange (BioseqPtr bsp, Pointer userdata)
{
  OrgNameChangePtr oncp;
  SeqDescrPtr sdp;
  BioSourcePtr biop;

  if (bsp == NULL
    || ! ISA_aa (bsp->mol)
    || (oncp = (OrgNameChangePtr)userdata) == NULL
    || oncp->oldname == NULL
    || oncp->newname == NULL)
  {
    return;
  }

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title)
    {
      ReplaceOldOrgNameInTitle (sdp, oncp);
    }
    else if (sdp->choice == Seq_descr_source)
    {
      biop = sdp->data.ptrvalue;
      if (biop != NULL && biop->org != NULL 
        && StringCmp (biop->org->taxname, oncp->oldname) == 0)
      {
        SetTaxNameAndRemoveTaxRef (biop->org, StringSave (oncp->newname));
      }   
    }   
  }
}

typedef struct parseformdata {
  FEATURE_FORM_BLOCK

  TexT           atleft;
  TexT           atright;
  DialoG         orgmod;
  DialoG         subsource;
  PopuP          taxname;
  Boolean        parsedef;
  Char           path [PATH_MAX];
  FILE           *fp;
  Boolean        replaceOldAsked;
  Boolean        doReplaceAll;
  Boolean        use_semicolon;
  ButtoN         leaveDlgUp;
} ParseFormData, PNTR ParseFormPtr;

typedef struct localidfinddata {
  CharPtr str;
  Int4    str_size;
} LocalIDFindData, PNTR LocalIDFindPtr;
  
static void LookForLocalID (BioseqPtr bsp, Pointer userdata)
{
  LocalIDFindPtr lidfp;
  SeqIdPtr      sip;

  if (bsp == NULL || (lidfp = (LocalIDFindPtr) userdata) == NULL) {
    return;
  }
  if (lidfp->str [0] != 0) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      SeqIdWrite (sip, lidfp->str, PRINTID_REPORT, lidfp->str_size);
      if (StringNICmp (lidfp->str, "tmpseq_", 7) == 0 ||
          StringNICmp (lidfp->str, "segseq_", 7) == 0 ||
          StringNICmp (lidfp->str, "SEG_dna", 7) == 0)
      {
        /* don't want to use this one */
        lidfp->str [0] = 0;
      } else {
        return;
      }
    }
  }
}

static void FindLocalIDForOrganism (
  SeqEntryPtr  sep,
  CharPtr      str,
  Int4         str_size
)
{
  LocalIDFindData lidfd;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  lidfd.str = str;
  lidfd.str_size = str_size;
  VisitBioseqsInSep (sep, (Pointer) &lidfd, LookForLocalID);
}

static void GetDataForOrganism (
  SeqEntryPtr  sep,
  CharPtr      str,
  Int4         str_size,
  CharPtr      useme,
  ParseFormPtr pfp
)
{
  SeqEntryPtr nsep;
  ValNodePtr    ttl;
  Char          txt [128];
  CharPtr       ptr;

  str [0] = 0;
  if (sep == NULL || pfp == NULL) return;
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL || nsep->data.ptrvalue == NULL) return;

  if (useme != NULL) {
    StringNCpy_0 (str, useme, str_size);
  } else if (pfp->parsedef) {
    ttl = SeqEntryGetSeqDescr (nsep, Seq_descr_title, NULL);
    if (ttl == NULL && sep != nsep) {
      ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
    }
    if (ttl == NULL || ttl->data.ptrvalue == NULL) return;
    StringNCpy_0 (str, (CharPtr) ttl->data.ptrvalue, str_size);
    GetTitle (pfp->atleft, txt, sizeof (txt));
    if (! StringHasNoText (txt)) {
      ptr = StringISearch (str, txt);
      if (ptr != NULL) {
        StringNCpy_0 (str, ptr + StringLen (txt), str_size);
      }
    }
    GetTitle (pfp->atright, txt, sizeof (txt));
    if (txt [0] != '\0' /* ! StringHasNoText (txt) */) {
      ptr = StringISearch (str, txt);
      if (ptr != NULL) {
        *ptr = '\0';
      }
    }
  } else {
    FindLocalIDForOrganism (nsep, str, str_size);
    if (str [0] == 0) {
      FindLocalIDForOrganism (sep, str, str_size);
    }
  }
  return;
}

static void ApplyDataToBioSourcePtr (
  ParseFormPtr pfp,
  BioSourcePtr biop,
  BioSourcePtr topbiop,
  SeqEntryPtr  sep,
  CharPtr      data
)
{
  Int2          taxnameval;
  ValNodePtr    vnp;
  OrgRefPtr     orp;
  OrgModPtr     mod;
  Uint1         orgmodtype;
  Uint1         subsourcetype;
  Uint1         subtype;
  OrgNamePtr    onp;
  OrgModPtr     tmpmod;
  SubSourcePtr  ssp;
  SubSourcePtr  tmpssp;

  taxnameval = GetValue (pfp->taxname);
  if (taxnameval == 2) {
    if (biop == NULL) {
      biop = BioSourceNew ();
      if (biop != NULL) {
        orp = OrgRefNew ();
        biop->org = orp;
        vnp = CreateNewDescriptor (sep, Seq_descr_source);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (Pointer) biop;
        }
      }
    }
    if (biop == NULL) return;
    orp = biop->org;
    if (orp == NULL) return;
    SetTaxNameAndRemoveTaxRef (orp, StringSave (data));
  }
  if (biop == NULL && topbiop != NULL) {
    biop = (BioSourcePtr) AsnIoMemCopy ((Pointer) topbiop,
                                        (AsnReadFunc) BioSourceAsnRead,
                                        (AsnWriteFunc) BioSourceAsnWrite);
    if (biop != NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_source);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) biop;
      }
    }
  }
  if (biop == NULL) return;
  orgmodtype = 0;
  subsourcetype = 0;
  if (taxnameval == 3) {
    orgmodtype = 255;
  } else if (taxnameval == 4) {
    subsourcetype = 255;
  } else {
    vnp = DialogToPointer (pfp->orgmod);
    if (vnp != NULL) {
      orgmodtype = vnp->choice;
      vnp = ValNodeFreeData (vnp);
    }
    if (orgmodtype == 0) {
      vnp = DialogToPointer (pfp->subsource);
      if (vnp != NULL) {
        subsourcetype = vnp->choice;
        vnp = ValNodeFreeData (vnp);
      }
    }
  }
  if (! StringHasNoText (data)) {
    if (orgmodtype > 0) {
      subtype = orgmodtype;
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
          mod = onp->mod;
          while (mod != NULL && mod->subtype != subtype) {
            mod = mod->next;
          }
          if (mod == NULL) {
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
          }
          if (mod != NULL) {
            mod->subtype = subtype;
            AppendOrReplaceString (&(mod->subname), data,
                                   &(pfp->replaceOldAsked),
                                   &(pfp->doReplaceAll),
                                   &(pfp->use_semicolon));
          }
        }
      }
    } else if (subsourcetype > 0) {
      subtype = subsourcetype;
      ssp = biop->subtype;
      while (ssp != NULL && ssp->subtype != subtype) {
        ssp = ssp->next;
      }
      if (ssp == NULL) {
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
      }
      if (ssp != NULL) {
        ssp->subtype = subtype;
        AppendOrReplaceString (&(ssp->name), data,
                               &(pfp->replaceOldAsked),
                               &(pfp->doReplaceAll),
                               &(pfp->use_semicolon));
      }
    }
  }
}

static void DoOneParseItem (Uint2 entityID, SeqEntryPtr sep, ParseFormPtr pfp, BioSourcePtr topbiop, CharPtr useme, SeqEntryPtr nps)

{
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  Char          str [256];
  SeqEntryPtr   tmp;
  OrgNameChangeData oncd;

  if (sep == NULL || pfp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      /* if nucprot set, look for biosource here */
      biop = NULL;
      SeqEntryToBioSource (sep, NULL, NULL, 0, &biop);
      if (biop != NULL && topbiop != NULL) {
        topbiop = biop;
        if (bssp->_class == BioseqseqSet_class_nuc_prot)
        {
          GetDataForOrganism ( sep, str, sizeof (str), useme, pfp);
          if (topbiop->org != NULL && topbiop->org->taxname != NULL)
          {
            oncd.oldname = topbiop->org->taxname;
            oncd.newname = str;
            VisitBioseqsInSep (sep, &oncd, FixProteinTitleAfterOrganismNameChange);
          }
          ApplyDataToBioSourcePtr (pfp, topbiop, NULL, sep, str);
          return;
        }
      }
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        DoOneParseItem (entityID, tmp, pfp, topbiop, useme, sep);
      }
      return;
    }
  }
  GetDataForOrganism ( sep, str, sizeof (str), useme, pfp);
  if (StringHasNoText (str)) return;
  biop = NULL;
  SeqEntryToBioSource (sep, NULL, NULL, 0, &biop);
  ApplyDataToBioSourcePtr (pfp, biop, topbiop, sep, str);
}

static void DoOneParse (Uint2 entityID, SeqEntryPtr sep, ParseFormPtr pfp, BioSourcePtr topbiop)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Char          gb [256];
  CharPtr       ptr;
  SeqIdPtr      sip;
  Char          str [250];
  SeqEntryPtr   tmp;

  if (sep == NULL || pfp == NULL) return;
  if (pfp->fp != NULL) {
    ReadLine (pfp->fp, str, sizeof (str)); /* FileGets on Mac sometimes sees '\r' instead of '\n' */
    while (Nlm_fileDone || str [0] != 0) {
      if (! StringHasNoText (str)) {
        ptr = StringChr (str, '\t');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
          TrimSpacesAroundString (str);
          TrimSpacesAroundString (ptr);
          sip = MakeSeqID (str);
          bsp = BioseqFind (sip);
          SeqIdFree (sip);
          if (bsp == NULL) {
            sprintf (gb, "gb|%s|", str);
            sip = MakeSeqID (gb);
            bsp = BioseqFind (sip);
            SeqIdFree (sip);
          }
          if (bsp != NULL) {
            tmp = GetBestTopParentForData (entityID, bsp);
             DoOneParseItem (entityID, tmp, pfp, topbiop, ptr, NULL);
          }
        }
      }
      str [0] = 0;
      ReadLine (pfp->fp, str, sizeof (str));
    }
    return;
  }
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        DoOneParse (entityID, tmp, pfp, topbiop);
      }
      return;
    }
  }
  DoOneParseItem (entityID, sep, pfp, topbiop, NULL, NULL);
}

static void DoParseDeflineProc (ButtoN b)

{
  BioSourcePtr  topbiop;
  ParseFormPtr  pfp;
  SeqEntryPtr   sep;

  pfp = (ParseFormPtr) GetObjectExtra (b);
  if (pfp == NULL || pfp->input_entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (pfp->input_entityID);
  if (sep == NULL) return;
  Hide (pfp->form);
  WatchCursor ();
  Update ();
  /* if popset, look for top biosource */
  SeqEntryToBioSource (sep, NULL, NULL, 0, &topbiop);
  if (! StringHasNoText (pfp->path)) {
    pfp->fp = FileOpen (pfp->path, "r");
  }
  DoOneParse (pfp->input_entityID, sep, pfp, topbiop);
  FileClose (pfp->fp);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (pfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, pfp->input_entityID, 0, 0);
  if (GetStatus (pfp->leaveDlgUp)) {
    Show (pfp->form);
  } else {
    Remove (pfp->form);
  }
}

static void ParseDeflineMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

static void ParseDefOrLocalIDToSource (BaseFormPtr bfp, Boolean parsedef, CharPtr path)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  GrouP              p;
  ParseFormPtr       pfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  pfp = (ParseFormPtr) MemNew (sizeof (ParseFormData));
  if (pfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Parse file to source", StdCloseWindowProc);
  SetObjectExtra (w, pfp, StdCleanupFormProc);
  pfp->form = (ForM) w;
  pfp->formmessage = ParseDeflineMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    pfp->appmessage = sepp->handleMessages;
  }

  pfp->input_entityID = bfp->input_entityID;
  pfp->input_itemID = bfp->input_itemID;
  pfp->input_itemtype = bfp->input_itemtype;

  pfp->parsedef = parsedef;
  StringNCpy_0 (pfp->path, path, sizeof (pfp->path));
  pfp->fp = NULL;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = NULL;
  if (parsedef) {
    g = HiddenGroup (h, 4, 0, NULL);
    StaticPrompt (g, "Get text between", 0, dialogTextHeight, programFont, 'l');
    pfp->atleft = DialogText (g, "", 10, NULL);
    StaticPrompt (g, "and", 0, dialogTextHeight, programFont, 'l');
    pfp->atright = DialogText (g, "", 10, NULL);
  }

  p = HiddenGroup (h, 6, 0, NULL);
  if (parsedef || path != NULL) {
    StaticPrompt (p, "Place in", 0, popupMenuHeight, programFont, 'l');
  } else {
    StaticPrompt (p, "Place localID in", 0, popupMenuHeight, programFont, 'l');
  }
  pfp->taxname = PopupList (p, TRUE, NULL);
  SetObjectExtra (pfp->taxname, pfp, NULL);
  PopupItem (pfp->taxname, " ");
  PopupItem (pfp->taxname, "Taxname");
  PopupItem (pfp->taxname, "OrgMod Note");
  PopupItem (pfp->taxname, "SubSource Note");
  SetValue (pfp->taxname, 1);
  StaticPrompt (p, "or", 0, popupMenuHeight, programFont, 'l');
  pfp->orgmod = OrgModTypeDialog (p, SHORT_SELECTION_LIST, NULL, NULL, FALSE, FALSE, TRUE);
  StaticPrompt (p, "or", 0, popupMenuHeight, programFont, 'l');
  pfp->subsource = SubSourceTypeDialog (p, SHORT_SELECTION_LIST, NULL, NULL, FALSE, FALSE, TRUE);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoParseDeflineProc);
  SetObjectExtra (b, pfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  pfp->leaveDlgUp = CheckBox (c, "Leave Dialog Up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) c, (HANDLE) g, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  Select (pfp->atleft);
  Update ();
}

static void ParseDefOrLocalIDToSourceMenuItem (IteM i, Boolean parsedef, CharPtr path)
{
  BaseFormPtr bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;

  ParseDefOrLocalIDToSource (bfp, parsedef, path);
}

extern void ParseFileToSource (IteM i)

{
  Char  path [PATH_MAX];

  if (GetInputFileName (path, sizeof (path), "", "TEXT")) {
    ParseDefOrLocalIDToSourceMenuItem (i, FALSE, path);
  }
}

static void TrimOrgNameCallback (BioSourcePtr biop, Pointer userdata)
{
  OrgRefPtr     orp;
  CharPtr       tmp, cp;
  CharPtr       word_break;
  Int4          space_len;

  if (biop == NULL) return;

  orp = biop->org;
  if (orp == NULL || StringHasNoText (orp->taxname)) return;

  tmp = StringSave (orp->taxname);
  cp = tmp;
  if (StringNICmp (tmp, "uncultured ", 11) == 0) {
    cp = tmp + 11;
  }

  word_break = StringStr (cp, " ");
  if (word_break == NULL) return;
  
  space_len = StringSpn (word_break, " ");

  if (StringNCmp (word_break + space_len, "sp.", 3) == 0
      || StringNCmp (word_break + space_len, "cf.", 3) == 0
      || StringNCmp (word_break + space_len, "aff.", 4) == 0)
  {
    word_break = StringStr (word_break + space_len, " ");
    if (word_break == NULL) return;
  }

  word_break = StringStr (word_break + space_len, " ");
  if (word_break == NULL) return;
  *word_break = 0;
  
  SetTaxNameAndRemoveTaxRef (orp, tmp);
}

extern void TrimOrganismName (IteM i)
{
  BaseFormPtr bfp;
  SeqEntryPtr sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  VisitBioSourcesInSep (sep, NULL, TrimOrgNameCallback);
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct addmodinfo {
  Boolean isOrgMod;
  Uint1 subtype;
  Boolean only_sp;
  Boolean only_cf;
  Boolean only_aff;
  Boolean no_taxid;
  CharPtr abbreviation;
  Boolean use_abbreviation;
} AddModInfo, PNTR AddModInfoPtr;

typedef struct addmodlistitem {
  Boolean isOrgMod;
  Uint1 subtype;
  CharPtr name;
  CharPtr abbreviation;
} AddModListItem, PNTR AddModListItemPtr;

static AddModListItem mods_list[] = {
{ TRUE, ORGMOD_strain, "Strain", NULL },
{ TRUE, ORGMOD_isolate, "Isolate", NULL },
{ TRUE, ORGMOD_authority, "Authority", NULL },
{ TRUE, ORGMOD_bio_material, "Bio-material", NULL },
{ TRUE, ORGMOD_biovar, "Biovar", "bv." },
{ FALSE, SUBSRC_clone, "Clone", NULL },
{ FALSE, SUBSRC_collection_date, "Collection-date", NULL },
{ FALSE, SUBSRC_country, "Country", NULL },
{ TRUE, ORGMOD_culture_collection, "Culture-collection", NULL },
{ TRUE, ORGMOD_forma, "Forma", "f." },
{ TRUE, ORGMOD_forma_specialis, "Forma-specialis", "f. sp." },
{ FALSE, SUBSRC_genotype, "Genotype", NULL },
{ FALSE, SUBSRC_haplotype, "Haplotype", NULL },
{ TRUE, ORGMOD_nat_host, "Host", NULL },
{ TRUE, ORGMOD_isolate, "Isolate", NULL },
{ TRUE, ORGMOD_pathovar, "Pathovar", "pv."},
{ TRUE, ORGMOD_serotype, "Serotype", NULL },
{ TRUE, ORGMOD_serovar, "Serovar", "serovar" },
{ TRUE, ORGMOD_specimen_voucher, "Specimen voucher", NULL },
{ TRUE, ORGMOD_strain, "Strain", NULL },
{ TRUE, ORGMOD_sub_species, "Sub-species", "subsp." },
{ TRUE, ORGMOD_variety, "Variety", "var." }
};


static void AddModToOrgProc (BioSourcePtr biop, Pointer userdata)
{
  OrgRefPtr     orp;
  OrgModPtr     mod;
  OrgNamePtr    onp;
  Boolean       influenza;
  CharPtr	str_to_add;
  size_t        len;
  CharPtr	ptr;
  CharPtr	str, tmp_name;
  AddModInfoPtr	ami;
  SubSourcePtr  ssp;
  Boolean       ok_to_add = FALSE;

  if (biop == NULL) return;

  ami = (AddModInfoPtr) userdata;

  orp = biop->org;
  if (orp == NULL) return;

  if (! ami->only_sp && ! ami->only_cf && ! ami->only_aff)
  {
    ok_to_add = TRUE;
  }
  if (ami->only_sp && StringStr (orp->taxname, " sp.") != NULL)
  {
    ok_to_add = TRUE;
  }
  if (ami->only_cf && StringStr (orp->taxname, " cf.") != NULL)
  {
    ok_to_add = TRUE;
  }
  if (ami->only_aff && StringStr (orp->taxname, " aff.") != NULL)
  {
    ok_to_add = TRUE;
  }

  if (ami->no_taxid && HasTaxonomyID (biop)) {
    ok_to_add = FALSE;
  }

  if (!ok_to_add)
  {
    return;
  }
  
  onp = orp->orgname;
  if (onp == NULL && ami->isOrgMod) 
  {
    return;
  }

  str_to_add = NULL;

  if (ami->isOrgMod)
  {
    mod = onp->mod;
    str_to_add = NULL;
    while (mod != NULL) {
      if (mod->subtype == ami->subtype) {
        str_to_add = StringSave (mod->subname);
      }
      mod = mod->next;
    }
  } else {
    if (biop->subtype != NULL) {
      ssp = biop->subtype;
      while (ssp != NULL) {
        if (ssp->subtype == ami->subtype) {
          str_to_add = StringSave (ssp->name);
        }
        ssp = ssp->next;
      }
    }
  }
  if (str_to_add != NULL) {
    ptr = StringStr (orp->taxname, " (");
    if (ptr != NULL) {
      *ptr = '\0';
    } else {
      ptr = StringStr (orp->taxname, " (");
      if (ptr != NULL) {
        *ptr = '\0';
      }
    }

    /* strip colons for structured modifiers */
    if (ami->isOrgMod
        && (ami->subtype == ORGMOD_bio_material
            || ami->subtype == ORGMOD_culture_collection
            || ami->subtype == ORGMOD_specimen_voucher)) {
      FindReplaceString (&str_to_add, ":", " ", FALSE, FALSE);
    }

    influenza = FALSE;
    if (StringICmp (orp->taxname, "Influenza A virus") == 0 ||
        StringICmp (orp->taxname, "Influenza B virus") == 0) {
      influenza = TRUE;
    }
    len = StringLen (orp->taxname) + StringLen (str_to_add) + 6;
    if (influenza) {
      len += 2;
    }
    if (ami->abbreviation != NULL && ami->use_abbreviation)
    {
      len += StringLen (ami->abbreviation) + 1;
    }
    str = MemNew (len);
    if (str != NULL) {
      StringCpy (str, orp->taxname);
      StringCat (str, " ");
      if (ami->abbreviation != NULL && ami->use_abbreviation)
      {
        StringCat (str, ami->abbreviation);
        StringCat (str, " ");
      }
      if (influenza) {
        StringCat (str, "(");
      }
      StringCat (str, str_to_add);
      if (influenza) {
        StringCat (str, ")");
      }
      
      if (influenza)
      {
        tmp_name = FixInfluenzaVirusName (str);
        if (tmp_name != NULL)
        {
          str = MemFree (str);
          str = tmp_name;
          tmp_name = NULL;
        }
      }
      
      SetTaxNameAndRemoveTaxRef (biop->org, str);
      RemoveOldName (biop->org);
    }
    str_to_add = MemFree (str_to_add);
  }
}

static void AddModToOrgFeat (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_ORG 
      || sfp->data.value.ptrvalue == NULL || userdata == NULL)
  {
    return;
  }
  else
  {
    AddModToOrgProc (sfp->data.value.ptrvalue, userdata);
  }
}

static void AddModToOrgDesc (SeqDescrPtr sdp, Pointer userdata, FilterSetPtr fsp)
{
  if (sdp == NULL || sdp->choice != Seq_descr_source 
      || sdp->data.ptrvalue == NULL || userdata == NULL)
  {
    return;
  }
  else
  {
    AddModToOrgProc (sdp->data.ptrvalue, userdata);
  }
}


typedef struct addmodformptr {
  FEATURE_FORM_BLOCK

  PopuP  modifier_to_add;
  ButtoN only_sp;
  ButtoN only_cf;
  ButtoN only_aff;
  ButtoN no_taxid;
  ButtoN use_abbreviation;
  DialoG constraint;

} AddModFormData, PNTR AddModFormPtr;

static void DoAddModToOrg (ButtoN b)
{
  AddModFormPtr	amfp;
  AddModInfoPtr ami;
  SeqEntryPtr   sep;
  Int4		      mod_index;
  FilterSetPtr  fsp;

  amfp = GetObjectExtra (b);
  if (amfp == NULL) return;
  Hide (amfp->form);

  sep = GetTopSeqEntryForEntityID (amfp->input_entityID);
  if (sep == NULL) return;
  ami = (AddModInfoPtr) MemNew (sizeof (AddModInfo));
  mod_index = GetValue (amfp->modifier_to_add);
  if (mod_index < 1
    || mod_index > sizeof (mods_list) / sizeof (AddModListItem))
    return;
  mod_index = mod_index - 1;

  ami->isOrgMod = mods_list[mod_index].isOrgMod;
  ami->subtype = mods_list[mod_index].subtype;
  ami->abbreviation = mods_list[mod_index].abbreviation;
  ami->only_sp = GetStatus (amfp->only_sp);
  ami->only_cf = GetStatus (amfp->only_cf);
  ami->only_aff = GetStatus (amfp->only_aff);
  ami->no_taxid = GetStatus (amfp->no_taxid);
  ami->use_abbreviation = GetStatus (amfp->use_abbreviation);
  /* always use abbreviation for serovar */
  if (ami->isOrgMod && ami->subtype == ORGMOD_serovar)
  {
  	ami->use_abbreviation = TRUE;
  }
  fsp = (FilterSetPtr) DialogToPointer (amfp->constraint);
  OperateOnSeqEntryConstrainedObjects (sep, fsp, 
                                       AddModToOrgFeat,
                                       AddModToOrgDesc,
                                       SEQFEAT_ORG,
                                       FEATDEF_ORG,
                                       Seq_descr_source,
                                       ami);
  fsp = FilterSetFree (fsp);

  Update ();
  ObjMgrSetDirtyFlag (amfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, amfp->input_entityID, 0, 0);

  if (GetStatus (amfp->leave_dlg_up)) 
  {
    Show (amfp->form);
  }
  else
  {
    Remove (amfp->form);
  }
}

static void ChangeModPopup (PopuP p)
{
  AddModFormPtr amfp;
  Int4          mod_index;
  
  amfp = (AddModFormPtr) GetObjectExtra (p);
  if (amfp == NULL) return;

  mod_index = GetValue (amfp->modifier_to_add);
  if (mod_index < 1
    || mod_index > sizeof (mods_list) / sizeof (AddModListItem))
    return;
  mod_index = mod_index - 1;
  if ( mods_list [mod_index].abbreviation != NULL)
  {
    Enable (amfp->use_abbreviation);
  }
  else
  {
    Disable (amfp->use_abbreviation);
  }
}


extern void AddModToOrgBaseForm (BaseFormPtr bfp)
{
  AddModFormPtr amfp;
  WindoW w;
  GrouP	g;
  GrouP	h;
  GrouP c;
  Int2 index;
  ButtoN b;

  if (bfp == NULL) return;

  amfp = (AddModFormPtr) MemNew (sizeof (AddModFormData));
  if (amfp == NULL) return;
  w = FixedWindow (-50, -33, -20, -10, "Append to Organism Name",
	StdCloseWindowProc);
  SetObjectExtra (w, amfp, StdCleanupFormProc);
  amfp->form = (ForM) w;

  amfp->input_entityID = bfp->input_entityID;
  amfp->input_itemID = bfp->input_itemID;

  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);

  amfp->modifier_to_add = PopupList (g, TRUE, ChangeModPopup);
  for (index=0; index< sizeof (mods_list) / sizeof (AddModListItem); index++)
  {
    PopupItem (amfp->modifier_to_add, mods_list[index].name);
  }
  SetObjectExtra (amfp->modifier_to_add, amfp, NULL);
  SetValue (amfp->modifier_to_add, 1);

  h = HiddenGroup (g, 1, 0, NULL);
  amfp->only_sp = CheckBox (h, "Only append to sp. organisms", NULL);
  amfp->only_cf = CheckBox (h, "Only append to cf. organisms", NULL);
  amfp->only_aff = CheckBox (h, "Only append to aff. organisms", NULL);
  amfp->no_taxid = CheckBox (h, "Only to organisms with no taxonomy ID", NULL);
  amfp->use_abbreviation = CheckBox (h, "Use abbreviation (for example, pv., subsp., etc.)", NULL);
  SetStatus (amfp->use_abbreviation, TRUE);
  Disable (amfp->use_abbreviation);
  
  amfp->constraint = FilterGroup (g, FALSE, TRUE, FALSE, FALSE, FALSE, NULL);

  c = HiddenGroup (g, 3, 0, NULL);
  b = DefaultButton(c, "Accept", DoAddModToOrg);
  SetObjectExtra(b, amfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  amfp->leave_dlg_up = CheckBox (c, "Leave Dialog Up", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) amfp->modifier_to_add, 
                              (HANDLE) h, 
                              (HANDLE) amfp->constraint,
                              (HANDLE) c, NULL);
  RealizeWindow(w);
  Show(w);
  Update();
}


extern void AddModToOrg (IteM i)
{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  AddModToOrgBaseForm (bfp);
}

typedef struct descformdata {
  FEATURE_FORM_BLOCK

  Uint2          oldEntityID;
  Uint4          oldItemID;
  Uint2          oldItemtype;
  Uint2          lookfor;
  Boolean        found;
  ObjMgrProcPtr  ompp;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqEntryPtr    sep;
  Handle         target;
  Boolean        usePopupForTarget;
  Boolean        hasMutPopPhySet;
  Int2           targetScratchSpace;
  Int2           nucProtCount;
  Int2           segSetCount;
  Uint2          descsubtype;
  ButtoN         createNewBtn;
  CharPtr        user_object_tag;
} DescFormData, PNTR DescFormPtr;


static Boolean IsUserObjectCorrectTag (UserObjectPtr uop, CharPtr user_object_tag)
{
  if (uop == NULL) {
    return FALSE;
  } else if (user_object_tag == NULL) {
    return TRUE;
  } else if (StringCmp (uop->type->str, user_object_tag) == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean FindDescrFunc (GatherContextPtr gcp)

{
  DescFormPtr  dfp;
  ValNodePtr   vnp;

  if (gcp == NULL) return TRUE;
  dfp = (DescFormPtr) gcp->userdata;
  if (dfp == NULL) return TRUE;
  if (gcp->thistype == OBJ_SEQDESC) {
    vnp = (ValNodePtr) gcp->thisitem;
    if (vnp != NULL && vnp->choice == dfp->lookfor 
        && (vnp->choice != Seq_descr_user || IsUserObjectCorrectTag(vnp->data.ptrvalue, dfp->user_object_tag))) {
      dfp->oldEntityID = gcp->entityID;
      dfp->oldItemID = gcp->itemID;
      dfp->oldItemtype = gcp->thistype;
      dfp->found = TRUE;
      return FALSE;
    }
  }
  return TRUE;
}

static CharPtr segClassList [] = {
  " ", "Nucleotide-Protein Set", "Segmented Nucleotide Set",
  "conset", "parts", "gibb", "gi", "genbank", "pir", "pub-set",
  "equiv", "swissprot", "pdb-entry", "MUTATION SET",
  "POPULATION SET", "PHYLOGENETIC SET", "ENVIRONMENTAL SAMPLES",
  "genomic product set", "other", NULL
};

static Int4 AllButPartsList (SeqEntryPtr sep, Pointer mydata,
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
  if (Bioseq_set_class (sep) != 4) {
    if (mycallback != NULL)
      (*mycallback) (sep, mydata, index, indent);
    index++;
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  sep = bssp->seq_set;
  indent++;
  while (sep != NULL) {
    index = AllButPartsList (sep, mydata, mycallback, index, indent);
    sep = sep->next;
  }
  return index;
}

#define AllButPartsCount( a )  AllButPartsList( a ,NULL,NULL,0,0);
#define AllButPartsExplore(a,b,c) AllButPartsList(a, b, c, 0L, 0);

static void PopTargetProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  DescFormPtr   dfp;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Uint1         _class;
  CharPtr       ptr;
  SeqIdPtr      sip;
  Char          str [128];

  dfp = (DescFormPtr) mydata;
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (sep->choice == 1) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      ptr = StringChr (str, '|');
      if (ptr == NULL) {
        ptr = str;
      } else {
        ptr++;
      }
      if (dfp->usePopupForTarget) {
        PopupItem (dfp->target, ptr);
      } else {
        ListItem (dfp->target, ptr);
      }
      if (bsp == dfp->bsp) {
        dfp->targetScratchSpace = index + 1;
      }
    } else if (sep->choice == 2) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      _class = bssp->_class;
      if (_class > 17) {
        _class = 18;
      }
      if (_class == 7 || (IsPopPhyEtcSet (bssp->_class))) {
        dfp->hasMutPopPhySet = TRUE;
      }
      if (! dfp->hasMutPopPhySet) {
        sprintf (str, "[%s]", segClassList [_class]);
      } else if (_class == 1) {
        (dfp->nucProtCount)++;
        sprintf (str, "[%s %d]", segClassList [_class], (int) dfp->nucProtCount);
        (dfp->segSetCount)++;
      } else if (_class == 2) {
        sprintf (str, "[%s %d]", segClassList [_class], (int) dfp->segSetCount);
      } else {
        sprintf (str, "[%s]", segClassList [_class]);
      }
      if (dfp->usePopupForTarget) {
        PopupItem (dfp->target, str);
      } else {
        ListItem (dfp->target, str);
      }
    }
  }
}

static Int2 PopulateTarget (DescFormPtr dfp, BioseqPtr bsp, Uint2 entityID)

{
  SeqEntryPtr  sep;
  Int2         val;

  val = 0;
  if (dfp != NULL) {
    if (bsp != NULL) {
      sep = SeqMgrGetSeqEntryForData (bsp);
      entityID = ObjMgrGetEntityIDForPointer (sep);
    }
    sep = GetTopSeqEntryForEntityID (entityID);
    dfp->sep = sep;
    if (sep != NULL) {
      dfp->hasMutPopPhySet = FALSE;
      dfp->targetScratchSpace = 1;
      dfp->nucProtCount = 0;
      dfp->segSetCount = 0;
      AllButPartsExplore (sep, (Pointer) dfp, PopTargetProc);
      val = dfp->targetScratchSpace;
    }
  }
  return val;
}

static void FindNewTargetProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  DescFormPtr   dfp;

  dfp = (DescFormPtr) mydata;
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (index + 1 == dfp->targetScratchSpace) {
      if (sep->choice == 1) {
        dfp->bsp = (BioseqPtr) sep->data.ptrvalue;
      } else if (sep->choice == 2) {
        dfp->bssp = (BioseqSetPtr) sep->data.ptrvalue;
      }
    }
  }
}

static Boolean ChangeTargetItemID (GatherContextPtr gcp)

{
  DescFormPtr  dfp;

  if (gcp == NULL) return TRUE;
  dfp = (DescFormPtr) gcp->userdata;
  if (dfp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (dfp->bsp == (BioseqPtr) gcp->thisitem) {
      dfp->input_entityID = gcp->entityID;
      dfp->input_itemID = gcp->itemID;
      dfp->input_itemtype = gcp->thistype;
      return FALSE;
    }
  } else if (gcp->thistype == OBJ_BIOSEQSET) {
    if (dfp->bssp == (BioseqSetPtr) gcp->thisitem) {
      dfp->input_entityID = gcp->entityID;
      dfp->input_itemID = gcp->itemID;
      dfp->input_itemtype = gcp->thistype;
      return FALSE;
    }
  }
  return TRUE;
}

static void PreventDupTitleCreation (DescFormPtr dfp)

{
  ValNodePtr  sdp;

  if (dfp == NULL) return;
  if (dfp->bsp != NULL) {
    sdp = dfp->bsp->descr;
  } else if (dfp->bssp != NULL) {
    sdp = dfp->bssp->descr;
  } else return;
  if (dfp->descsubtype != Seq_descr_title) return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_title) {
      SafeDisable (dfp->createNewBtn);
      return;
    }
    sdp = sdp->next;
  }
  SafeEnable (dfp->createNewBtn);
}

static void ChangeNewDescTarget (Handle obj)

{
  DescFormPtr  dfp;
  GatherScope  gs;

  dfp = (DescFormPtr) GetObjectExtra (obj);
  if (dfp != NULL) {
    dfp->bsp = NULL;
    dfp->bssp = NULL;
    dfp->targetScratchSpace = GetValue (dfp->target);
    AllButPartsExplore (dfp->sep, (Pointer) dfp, FindNewTargetProc);
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 0;
    MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQSET] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    GatherEntity (dfp->input_entityID, (Pointer) dfp, ChangeTargetItemID, &gs);
    PreventDupTitleCreation (dfp);
  }
}



static void CreateNewDescProc (ButtoN b)

{
  DescFormPtr    dfp;
  OMProcControl  ompc;
  ObjMgrProcPtr  ompp;
  Int2           retval;

  dfp = (DescFormPtr) GetObjectExtra (b);
  if (dfp != NULL) {
    Hide (dfp->form);
    ompp = dfp->ompp;
    MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
    ompc.input_entityID = dfp->input_entityID;
    ompc.input_itemID = dfp->input_itemID;
    ompc.input_itemtype = dfp->input_itemtype;
    GatherDataForProc (&ompc, FALSE);
    ompc.proc = ompp;

    retval = (*(ompp->func)) (&ompc);
    if (retval == OM_MSG_RET_ERROR) {
      ErrShow ();
    }
    Update ();
    Remove (dfp->form);
  }
}

static void EditOldDescProc (ButtoN b)

{
  DescFormPtr  dfp;
  Int2         handled;

  dfp = (DescFormPtr) GetObjectExtra (b);
  if (dfp != NULL) {
    Hide (dfp->form);
    if (dfp->oldEntityID > 0 && dfp->oldItemID > 0 && dfp->oldItemtype > 0) {
      WatchCursor ();
      handled = GatherProcLaunch (OMPROC_EDIT, FALSE, dfp->oldEntityID,
                                  dfp->oldItemID, dfp->oldItemtype,
                                  0, 0, dfp->oldItemtype, 0);
      ArrowCursor ();
      if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
        Message (MSG_ERROR, "Unable to edit existing descriptor.");
      }
    }
    Update ();
    Remove (dfp->form);
  }
}

static CharPtr  newDescMsg = "\
Descriptors may apply to a single sequence or to a set of \
sequences.  Please set the target control to the desired \
set or sequence.";

static CharPtr  editOldMsg = "\
A descriptor of this type already exists at or above the \
displayed target.  You may want to edit it instead of \
creating a new descriptor.";

static void NewDescriptorMessageProc (ForM f, Int2 mssg)

{
  DescFormPtr  dfp;

  dfp = (DescFormPtr) GetObjectExtra (f);
  if (dfp != NULL) {
    if (dfp->appmessage != NULL) {
      dfp->appmessage (f, mssg);
    }
  }
}

extern void NewDescriptorMenuFuncEx (ObjMgrProcPtr ompp, BaseFormPtr bfp, Uint2 descsubtype, CharPtr user_object_tag)

{
  ButtoN             b;
  BioseqPtr          bsp;
  GrouP              c;
  Int4               count;
  DescFormPtr        dfp;
  Uint2              entityID;
  GrouP              g;
  GatherScope        gs;
  GrouP              h;
  GrouP              p1, p2;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  Int2               val;
  WindoW             w;
  GrouP              target_grp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#endif
  if (bfp == NULL) return;
  if (ompp == NULL || ompp->func == NULL) return;
  if (ompp->inputtype != OBJ_SEQDESC) return;
  dfp = (DescFormPtr) MemNew (sizeof (DescFormData));
  if (dfp == NULL) return;
  dfp->ompp = ompp;
  dfp->lookfor = ompp->subinputtype;
  if (dfp->lookfor == Seq_descr_user) {
    dfp->user_object_tag = user_object_tag;
  }
  dfp->found = FALSE;
  dfp->descsubtype = descsubtype;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 0;
  gs.get_feats_location = TRUE;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQSET] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  if (bsp != NULL) {
    slp = ValNodeNew (NULL);
    slp->choice = SEQLOC_WHOLE;
    sip = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
    slp->data.ptrvalue = sip;
    gs.target = slp;
  }
  dfp->bsp = bsp;
  GatherEntity (bfp->input_entityID, (Pointer) dfp, FindDescrFunc, &gs);
  SeqLocFree (gs.target);
  w = FixedWindow (-50, -33, -10, -10, "Descriptor Target Control", StdCloseWindowProc);
  SetObjectExtra (w, dfp, StdCleanupFormProc);
  dfp->form = (ForM) w;
  dfp->formmessage = NewDescriptorMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    dfp->appmessage = sepp->handleMessages;
  }

  dfp->input_entityID = bfp->input_entityID;
  dfp->input_itemID = bfp->input_itemID;
  dfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  target_grp = HiddenGroup (h, -1, 0, NULL);
  
  p1 = MultiLinePrompt (target_grp, newDescMsg, 25 * stdCharWidth, programFont);
    
  entityID = bfp->input_entityID;
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
  }
  sep = GetTopSeqEntryForEntityID (entityID);
  count = AllButPartsCount (sep);
  if (count < 32) {
    g = HiddenGroup (target_grp, 2, 0, NULL);
    StaticPrompt (g, "Target", 0, popupMenuHeight, programFont, 'l');
    dfp->usePopupForTarget = TRUE;
    dfp->target = PopupList (g, TRUE, (PupActnProc) ChangeNewDescTarget);
  } else {
    g = HiddenGroup (target_grp, 0, 2, NULL);
    StaticPrompt (g, "Target", 0, 0, programFont, 'c');
    dfp->usePopupForTarget = FALSE;
    dfp->target = SingleList (g, 14, 3, (LstActnProc) ChangeNewDescTarget);
  }
  SetObjectExtra (dfp->target, dfp, NULL);
  val = PopulateTarget (dfp, bsp, bfp->input_entityID);
  SetValue (dfp->target, val);
  ChangeNewDescTarget ((Handle) dfp->target);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p1, (HANDLE) g, NULL);
  
  p2 = NULL;
  if (dfp->found) {
    p2 = MultiLinePrompt (h, editOldMsg, 25 * stdCharWidth, programFont);
  }
  
  c = HiddenGroup (h, 4, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  dfp->createNewBtn = PushButton (c, "Create New", CreateNewDescProc);
  SetObjectExtra (dfp->createNewBtn, dfp, NULL);
  if (dfp->found) {
    b = DefaultButton (c, "Edit Old", EditOldDescProc);
    SetObjectExtra (b, dfp, NULL);
  }
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) target_grp, (HANDLE) c, (HANDLE) p2, NULL);
  RealizeWindow (w);
  PreventDupTitleCreation (dfp);
  Show (w);
  Update ();
}

extern void NewDescriptorMenuFunc (ObjMgrProcPtr ompp, BaseFormPtr bfp, Uint2 descsubtype)
{
  NewDescriptorMenuFuncEx (ompp, bfp, descsubtype, NULL);
}


static void NewDescriptorMenuProc (IteM i)

{
  NewObjectPtr  nop;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
  NewDescriptorMenuFunc (nop->ompp, nop->bfp, nop->descsubtype);
}


static void NewDescriptorMenuProcExtra (IteM i)

{
  NewObjectPtr  nop;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
  if (nop->blurb == NULL) {
    NewDescriptorMenuFunc (nop->ompp, nop->bfp, nop->descsubtype);
  } else {
    NewDescriptorMenuFuncEx (nop->ompp, nop->bfp, nop->descsubtype, nop->blurb);
  }
}


extern void SetupNewDescriptorsMenu (MenU m, BaseFormPtr bfp)

{
  Boolean        allowgenbank, allow_structured_comment = FALSE;
  IteM           i;
  NewObjectPtr   nop;
  ObjMgrPtr      omp;
  ObjMgrProcPtr  ompp;
  ObjMgrTypePtr  omtp;

  if (m == NULL) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, 0, 0, NULL);
  if (ompp == NULL) return;
  omtp = NULL;
  allowgenbank = FALSE;
/*#ifdef EXTRA_SERVICES*/
  if (extraServices) {
    allowgenbank = TRUE;
  }
/*#endif*/
  if (indexerVersion) {
    allow_structured_comment = TRUE;
  }

  while ((omtp = ObjMgrTypeFindNext (omp, omtp)) != NULL) {
    ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, omtp->datatype, 0, NULL);
    if (ompp != NULL) {
      switch (omtp->datatype) {
        case OBJ_SEQDESC :
          ompp = NULL;
          while ((ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT,
				  omtp->datatype, 0, ompp)) != NULL) {
            if (ompp->subinputtype != Seq_descr_pub) {
              if (!allowgenbank && ompp->subinputtype == Seq_descr_genbank) {
                /* skip */
              } else if (!allow_structured_comment && ompp->subinputtype == Seq_descr_user && StringCmp (ompp->proclabel, "Structured Comment") == 0) {
                /* skip */
              } else if (ompp->subinputtype == Seq_descr_user) {
                i = CommandItem (m, ompp->proclabel, NewDescriptorMenuProcExtra);
                nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
                if (nop != NULL) {
                  nop->ompp = ompp;
                  nop->bfp = bfp;
                  nop->descsubtype = ompp->subinputtype;
                  if (StringCmp (ompp->proclabel, "DBLink") == 0) {
                    nop->blurb = "DBLink"; 
                  } else if (StringCmp (ompp->proclabel, "GenomeProjectsDB") == 0) {
                    nop->blurb = "GenomeProjectsDB";
                  } else if (StringCmp (ompp->proclabel, "Structured Comment") == 0) {
                    nop->blurb = "StructuredComment";
                  } else if (StringCmp (ompp->proclabel, "TPA Assembly") == 0) {
                    nop->blurb = "TpaAssembly";
                  } else if (StringCmp (ompp->proclabel, "RefGene Tracking") == 0) {
                    nop->blurb = "RefGeneTracking";
                  }
                }
                SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);               
              } else {
                i = CommandItem (m, ompp->proclabel, NewDescriptorMenuProc);
                nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
                if (nop != NULL) {
                  nop->ompp = ompp;
                  nop->bfp = bfp;
                  nop->descsubtype = ompp->subinputtype;
                }
                SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);
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

static CharPtr  editOldDescMsg = "\
You may really want to edit an existing BioSource descriptor instead.\n\
Proceed anyway?";

static void NewFeatureMenuProc (IteM i)

{
  MsgAnswer      ans;
  BaseFormPtr    bfp;
  NewObjectPtr   nop;
  OMProcControl  ompc;
  ObjMgrProcPtr  ompp;
  Int2           retval;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  ompp = nop->ompp;
  if (ompp == NULL || ompp->func == NULL) return;
  if (ompp->subinputtype == FEATDEF_BIOSRC) {
    ans = Message (MSG_YN, editOldDescMsg);
    if (ans == ANS_NO) return;
  }
  MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
  ompc.input_entityID = bfp->input_entityID;
  ompc.input_itemID = bfp->input_itemID;
  ompc.input_itemtype = bfp->input_itemtype;
  GatherDataForProc (&ompc, FALSE);
  ompc.proc = ompp;
  retval = (*(ompp->func)) (&ompc);
  if (retval == OM_MSG_RET_ERROR) {
    ErrShow ();
  }
}

#ifdef WIN_MAC
VoidPtr macUserDataPtr = NULL;
#endif

extern void EnableFeaturesPerTarget (BaseFormPtr bfp)

{
  BioseqPtr     bsp;
  Uint1         mol;
  NewObjectPtr  nop;

  if (bfp == NULL) return;
#ifdef WIN_MAC
  nop = (NewObjectPtr) macUserDataPtr;
#else
  nop = (NewObjectPtr) bfp->userDataPtr;
#endif
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  mol = 0;
  if (bsp != NULL) {
    mol = bsp->mol;
  }
  while (nop != NULL) {
    if (nop->kind == 2) {
      /* analysis menu item, ignore it */
    } else if (mol == 0) {
      SafeDisable (nop->item);
    } else if (ISA_na (mol) && (nop->molgroup == 2 || nop->molgroup == 3)) {
      SafeEnable (nop->item);
    } else if (ISA_aa (mol) && (nop->molgroup == 1 || nop->molgroup == 3)) {
      SafeEnable (nop->item);
    } else {
      SafeDisable (nop->item);
    }
    nop = nop->next;
  }
}

static VoidPtr LinkNewObjectLists (NewObjectPtr list1, NewObjectPtr list2)

{
  NewObjectPtr  nop;

  if (list1 == NULL) return list2;
  nop = list1;
  while (nop->next != NULL) {
    nop = nop->next;
  }
  nop->next = list2;
  return list1;
}


extern Boolean IsUnwantedFeatureType (Uint1 key)
{
#if 1
  return FALSE;
#else
  if (key == FEATDEF_satellite || key == FEATDEF_repeat_unit) {
    return TRUE;
  } else {
    return FALSE;
  }
#endif
}

extern void SetupNewFeaturesMenu (MenU m, BaseFormPtr bfp)

{
  FeatDispGroupPtr  fdgp;
  FeatDefPtr        fdp;
  NewObjectPtr      first;
  IteM              i;
  Uint1             key;
  CharPtr           label;
  NewObjectPtr      last;
  NewObjectPtr      nop;
  ObjMgrPtr         omp;
  ObjMgrProcPtr     ompp;
  ObjMgrTypePtr     omtp;
  MenU              sub;
  Uint2             subtype;

  if (m == NULL) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, 0, 0, NULL);
  if (ompp == NULL) return;
  omtp = NULL;
  first = NULL;
  last = NULL;
  while ((omtp = ObjMgrTypeFindNext (omp, omtp)) != NULL) {
    ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, omtp->datatype, 0, NULL);
    if (ompp != NULL) {
      switch (omtp->datatype) {
        case OBJ_SEQFEAT :
          fdgp = NULL;
          while ((fdgp = DispGroupFindNext (fdgp, &key, &label)) != NULL) {
            if (fdgp->groupkey != 0) {
              sub = SubMenu (m, fdgp->groupname);
              fdp = NULL;
              label = NULL;
              while ((fdp = FeatDefFindNext (fdp, &key, &label,
                     fdgp->groupkey, FALSE)) != NULL) {
                if (key != FEATDEF_BAD) {
                  ompp = NULL;
                  while ((ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT,
				          omtp->datatype, 0, ompp)) != NULL) {
                    if (ompp->subinputtype == fdp->featdef_key &&
                        ompp->subinputtype != FEATDEF_PUB) {
                      i = CommandItem (sub, ompp->proclabel, NewFeatureMenuProc);
                      nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
                      if (nop != NULL) {
                        nop->kind = 1; /* feature creation item */
                        nop->ompp = ompp;
                        nop->bfp = bfp;
                        nop->item = i;
                        nop->molgroup = fdp->molgroup;
                      }
                      if (first == NULL) {
                        first = nop;
                      }
                      if (last != NULL) {
                        last->next = nop;
                      }
                      last = nop;
                      SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);
                    }
                  }
                }
              }
            }
          }
          /* if (indexerVersion) { */
            sub = SubMenu (m, "Remaining Features");
            fdp = NULL;
            label = NULL;
            while ((fdp = FeatDefFindNext (fdp, &key, &label, 0, FALSE)) != NULL) {
              if (key != FEATDEF_BAD && !IsUnwantedFeatureType(key)) {
                ompp = NULL;
                while ((ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT,
				                                           omtp->datatype, 0, ompp)) != NULL) {
				          subtype = ompp->subinputtype;
                  if (subtype == fdp->featdef_key && OkToListFeatDefInRemainingFeatures (subtype)) {
                    i = CommandItem (sub, ompp->proclabel, NewFeatureMenuProc);
                    nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
                    if (nop != NULL) {
                      nop->kind = 1; /* feature creation item */
                      nop->ompp = ompp;
                      nop->bfp = bfp;
                      nop->item = i;
                      nop->molgroup = fdp->molgroup;
                    }
                    if (first == NULL) {
                      first = nop;
                    }
                    if (last != NULL) {
                      last->next = nop;
                    }
                    last = nop;
                    SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);
                  }
                }
              }
            }
          /* } */
          break;
        default :
          break;
      }
    }
  }
#ifdef WIN_MAC
  macUserDataPtr = LinkNewObjectLists (macUserDataPtr, first);
#else
  bfp->userDataPtr = LinkNewObjectLists (bfp->userDataPtr, first);
#endif
}

extern void SetupNewPublicationsMenu (MenU m, BaseFormPtr bfp)

{
  NewObjectPtr   first;
  IteM           i;
  NewObjectPtr   nop;
  ObjMgrPtr      omp;
  ObjMgrProcPtr  ompp;
  ObjMgrTypePtr  omtp;

  if (m == NULL) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, 0, 0, NULL);
  if (ompp == NULL) return;
  omtp = NULL;
  while ((omtp = ObjMgrTypeFindNext (omp, omtp)) != NULL) {
    ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, omtp->datatype, 0, NULL);
    if (ompp != NULL) {
      switch (omtp->datatype) {
        case OBJ_SEQFEAT :
            ompp = NULL;
            while ((ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT,
				    omtp->datatype, 0, ompp)) != NULL) {
              if (ompp->subinputtype == FEATDEF_PUB) {
                i = CommandItem (m, "Publication Feature", NewFeatureMenuProc);
                nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
                if (nop != NULL) {
                  nop->kind = 1; /* feature creation item */
                  nop->ompp = ompp;
                  nop->bfp = bfp;
                  nop->item = i;
                  nop->molgroup = 3;
#ifdef WIN_MAC
                  first = (NewObjectPtr) macUserDataPtr;
#else
                  first = (NewObjectPtr) bfp->userDataPtr;
#endif
                  if (first != NULL) {
                    while (first->next != NULL) {
                      first = first->next;
                    }
                    first->next = nop;
                  }
                }
                SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);
              }
            }
          break;
        case OBJ_SEQDESC :
          ompp = NULL;
          while ((ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT,
				  omtp->datatype, 0, ompp)) != NULL) {
            if (ompp->subinputtype == Seq_descr_pub) {
              i = CommandItem (m, "Publication Descriptor", NewDescriptorMenuProc);
              nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
              if (nop != NULL) {
                nop->ompp = ompp;
                nop->bfp = bfp;
              }
              SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);
            }
          }
          break;
        default :
          break;
      }
    }
  }
}

#ifdef WIN_MAC
extern IteM  addSecondaryItem;
#endif

extern void SetupEditSecondary (MenU m, BaseFormPtr bfp)

{
  IteM           i;
  NewObjectPtr   nop;
  ObjMgrPtr      omp;
  ObjMgrProcPtr  ompp;
  ObjMgrTypePtr  omtp;

  if (m == NULL) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, 0, 0, NULL);
  if (ompp == NULL) return;
  omtp = NULL;
  while ((omtp = ObjMgrTypeFindNext (omp, omtp)) != NULL) {
    ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT, omtp->datatype, 0, NULL);
    if (ompp != NULL) {
      switch (omtp->datatype) {
        case OBJ_SEQDESC :
          ompp = NULL;
          while ((ompp = ObjMgrProcFindNext (omp, OMPROC_EDIT,
				  omtp->datatype, 0, ompp)) != NULL) {
            if (ompp->subinputtype != Seq_descr_pub) {
              if (ompp->subinputtype == Seq_descr_genbank) {
                i = CommandItem (m, "Add Secondary", NewDescriptorMenuProc);
                nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
                if (nop != NULL) {
                  nop->ompp = ompp;
                  nop->bfp = bfp;
                }
                SetObjectExtra (i, (Pointer) nop, StdCleanupExtraProc);
#ifdef WIN_MAC
                if (addSecondaryItem == NULL) {
                  addSecondaryItem = i;
                }
#endif
                return;
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

typedef struct gbformdata {
  FEATURE_FORM_BLOCK

  ButtoN         clearExtraAccns;
  ButtoN         clearSource;
  ButtoN         clearKeywords;
  ButtoN         clearOrigin;
  ButtoN         clearOldDate;
  ButtoN         clearEntryDate;
  ButtoN         clearDivision;
  ButtoN         clearTaxonomy;
} GbeditFormData, PNTR GbeditFormPtr;

static void EditGenbankCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Boolean        empty;
  GBBlockPtr     gbp;
  GbeditFormPtr  gfp;
  ValNodePtr     nextsdp;
  Pointer PNTR   prevsdp;
  ValNodePtr     sdp;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  gfp = (GbeditFormPtr) mydata;
  if (gfp == NULL) return;
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
      if (GetStatus (gfp->clearExtraAccns)) {
        gbp->extra_accessions = ValNodeFreeData (gbp->extra_accessions);
      }
      if (GetStatus (gfp->clearSource)) {
        gbp->source = MemFree (gbp->source);
      }
      if (GetStatus (gfp->clearKeywords)) {
        gbp->keywords = ValNodeFreeData (gbp->keywords);
      }
      if (GetStatus (gfp->clearOrigin)) {
        gbp->origin = MemFree (gbp->origin);
      }
      if (GetStatus (gfp->clearOldDate)) {
        gbp->date = MemFree (gbp->date);
      }
      if (GetStatus (gfp->clearEntryDate)) {
        gbp->entry_date = DateFree (gbp->entry_date);
      }
      if (GetStatus (gfp->clearDivision)) {
        gbp->div = MemFree (gbp->div);
      }
      if (GetStatus (gfp->clearTaxonomy)) {
        gbp->taxonomy = MemFree (gbp->taxonomy);
      }
      if (gbp->extra_accessions == NULL && gbp->source == NULL &&
          gbp->keywords == NULL && gbp->origin == NULL &&
          gbp->date == NULL && gbp->entry_date == NULL &&
          gbp->div == NULL && gbp->taxonomy == NULL) {
        empty = TRUE;
        ObjMgrDeSelect (0, 0, 0, 0, NULL);
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

static void DoEditGenbank (ButtoN b)

{
  GbeditFormPtr  gfp;
  SeqEntryPtr    sep;

  gfp = GetObjectExtra (b);
  if (gfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (gfp->input_entityID);
  if (sep == NULL) return;
  Hide (gfp->form);
  WatchCursor ();
  Update ();
  SeqEntryExplore (sep, (Pointer) gfp, EditGenbankCallback);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (gfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, gfp->input_entityID, 0, 0);
  Remove (gfp->form);
}

static void EditGenbankMessageProc (ForM f, Int2 mssg)

{
  GbeditFormPtr  gfp;

  gfp = (GbeditFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    if (gfp->appmessage != NULL) {
      gfp->appmessage (f, mssg);
    }
  }
}

extern void EditGenbankElements (Handle i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GbeditFormPtr      gfp;
  GrouP              h;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  gfp = (GbeditFormPtr) MemNew (sizeof (GbeditFormData));
  if (gfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "GenBank Block Removal", StdCloseWindowProc);
  SetObjectExtra (w, gfp, StdCleanupFormProc);
  gfp->form = (ForM) w;
  gfp->formmessage = EditGenbankMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    gfp->appmessage = sepp->handleMessages;
  }

  gfp->input_entityID = bfp->input_entityID;
  gfp->input_itemID = bfp->input_itemID;
  gfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, -1, 0, NULL);

  gfp->clearExtraAccns = CheckBox (g, "Clear Secondary Accessions", NULL);
  gfp->clearSource = CheckBox (g, "Clear Source Line", NULL);
  gfp->clearKeywords = CheckBox (g, "Clear Keywords", NULL);
  gfp->clearOrigin = CheckBox (g, "Clear Origin", NULL);
  gfp->clearOldDate = CheckBox (g, "Clear Old Date", NULL);
  gfp->clearEntryDate = CheckBox (g, "Clear Entry Date", NULL);
  gfp->clearDivision = CheckBox (g, "Clear Division", NULL);
  gfp->clearTaxonomy = CheckBox (g, "Clear Lineage", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoEditGenbank);
  SetObjectExtra (b, gfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

typedef struct helpindex {
  Int2           item;
  CharPtr        heading;
  CharPtr        section;
  CharPtr        combined;
} HelpIndex, PNTR HelpIndexPtr;

typedef struct helpform {
  FORM_MESSAGE_BLOCK
  DoC            doc;
  DoC            list;
  TexT           findTxt;
  ButtoN         findBtn;
  GrouP          dismissGrp;
  ValNodePtr     mainStrings;
  ValNodePtr     indexStrings;
  ValNodePtr     index;
  Char           file [PATH_MAX];
} HelpForm, PNTR HelpFormPtr;

#define TBL_FMT 1
#define TXT_FMT 2
#define HDG_FMT 3
#define SUB_FMT 4  /* and beyond, encoding indent */

static ParData lstParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData lstColFmt = {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

static ParData hdgParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData hdgColFmt = {0, 0, 80, 0, NULL, 'c', TRUE, FALSE, FALSE, FALSE, TRUE};

static ParData subParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData subColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static ParData tblParFmt = {TRUE, FALSE, FALSE, FALSE, TRUE, 0, 0};
static ColData tblColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static ParData txtParFmt = {TRUE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData txtColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static void ResetHelpLists (HelpFormPtr hfp)

{
  HelpIndexPtr  hip;
  ValNodePtr    vnp;

  if (hfp != NULL) {
    hfp->mainStrings = ValNodeFreeData (hfp->mainStrings);
    hfp->indexStrings = ValNodeFreeData (hfp->indexStrings);
    vnp = hfp->index;
    while (vnp != NULL) {
      hip = (HelpIndexPtr) vnp->data.ptrvalue;
      if (hip != NULL) {
        hip->heading = MemFree (hip->heading);
        hip->section = MemFree (hip->section);
        hip->combined = MemFree (hip->combined);
      }
      vnp = vnp->next;
    }
    hfp->index = ValNodeFreeData (hfp->index);
  }
}

static void TrimTrailingSpaces (CharPtr str)

{
  size_t  len;

  if (str != NULL && str [0] != '\0') {
    len = StringLen (str);
    while (len > 0 && str [len - 1] == ' ') {
      len--;
    }
    str [len] = '\0';
  }
}

static Boolean ParseHelpFile (HelpFormPtr hfp, Boolean printPath)

{
  Char          ch;
  Uint1         choice;
  FileCache     fc;
  FILE          *fp;
  Boolean       goOn;
  Char          heading [512];
  HelpIndexPtr  hip;
  Boolean       inTable;
  Int2          level;
  ValNodePtr    list;
  Int2          numItems;
  Char          path [PATH_MAX];
  CharPtr       ptr;
  Char          section [512];
  Char          str [512];
  ValNodePtr    vnp;

  if (hfp == NULL || hfp->doc == NULL) return FALSE;
  Hide (hfp->list);
  Hide (hfp->doc);
  Update ();
  Reset (hfp->list);
  Reset (hfp->doc);
  ResetHelpLists (hfp);
  if (hfp->file == NULL || hfp->file [0] == '\0') return FALSE;
  numItems = 0;
  ProgramPath (path, sizeof (path));
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    *ptr = '\0';
  }
  FileBuildPath (path, NULL, hfp->file);
  fp = FileOpen (path, "r");
  if (fp == NULL) {
    if (GetAppParam ("NCBI", "ErrorProcessing", "MsgPath", NULL, path, sizeof (path))) {
      FileBuildPath (path, NULL, hfp->file);
      fp = FileOpen (path, "r");
    }
  }
  if (fp == NULL) {
    if (GetAppParam ("NCBI", "NCBI", "DATA", NULL, path, sizeof (path))) {
      FileBuildPath (path, NULL, hfp->file);
      fp = FileOpen (path, "r");
    }
  }
  if (fp != NULL) {
    if (printPath) {
      SetTitle (hfp->form, path);
    }
    list = NULL;
    heading [0] = '\0';
    section [0] = '\0';
    inTable = FALSE;
    if (! FileCacheSetup (&fc, fp)) return FALSE;
    goOn = (FileCacheGetString (&fc, str, sizeof (str)) != NULL);
    while (goOn) {
      ptr = str;
      ch = *ptr;
      while (ch != '\n' && ch != '\r' && ch != '\0') {
        ptr++;
        ch = *ptr;
      }
      *ptr = '\0';
      if (inTable) {
        TrimTrailingSpaces (str);
      } else {
        TrimSpacesAroundString (str);
      }
      ch = str [0];
      if (ch == '>' || ch == '*' || ch == '#' || ch == '!') {
        if (list != NULL) {
          ptr = MergeValNodeStrings (list, inTable);
          if (inTable) {
            choice = TBL_FMT;
          } else {
            choice = TXT_FMT;
          }
          numItems++;
          vnp = ValNodeAdd (&(hfp->mainStrings));
          if (vnp != NULL) {
            vnp->choice = choice;
            vnp->data.ptrvalue = ptr;
          }
          /* ptr = MemFree (ptr); */
          list = ValNodeFreeData (list);
        }
        ch = str [0];
        if (ch == '>') {
          StringNCpy_0 (heading, str + 1, sizeof (heading));
          numItems++;
          vnp = ValNodeAdd (&(hfp->mainStrings));
          if (vnp != NULL) {
            vnp->choice = HDG_FMT;
            vnp->data.ptrvalue = StringSave (heading);
          }
          vnp = ValNodeAdd (&(hfp->indexStrings));
          if (vnp != NULL) {
            vnp->choice = 0;
            vnp->data.ptrvalue = StringSave (heading);
          }
          vnp = ValNodeAdd (&(hfp->index));
          if (vnp != NULL) {
            hip = (HelpIndexPtr) MemNew (sizeof (HelpIndex));
            vnp->data.ptrvalue = (Pointer) hip;
            if (hip != NULL) {
              hip->item = numItems;
              hip->heading = StringSave (heading);
            }
          }
        } else if (ch == '*') {
          level = 1;
          ch = str [level];
          while (ch == '*') {
            level++;
            ch = str [level];
          }
          StringNCpy_0 (section, str + level, sizeof (section));
          numItems++;
          vnp = ValNodeAdd (&(hfp->mainStrings));
          if (vnp != NULL) {
            if (level < 2) {
              vnp->choice = SUB_FMT;
            } else {
              vnp->choice = 5 * (level - 1);
            }
            vnp->data.ptrvalue = StringSave (section);
          }
          vnp = ValNodeAdd (&(hfp->indexStrings));
          if (vnp != NULL) {
            vnp->choice = 5 * level;
            vnp->data.ptrvalue = StringSave (section);
          }
          vnp = ValNodeAdd (&(hfp->index));
          if (vnp != NULL) {
            hip = (HelpIndexPtr) MemNew (sizeof (HelpIndex));
            vnp->data.ptrvalue = (Pointer) hip;
            if (hip != NULL) {
              hip->item = numItems;
              hip->heading = StringSave (heading);
              hip->section = StringSave (section);
              sprintf (str, "%s|%s", heading, section);
              hip->combined = StringSave (str);
            }
          }
        } else if (ch == '#' || ch == '!') {
          inTable = (Boolean) (ch == '!');
          vnp = ValNodeAdd (&list);
          if (vnp != NULL) {
            if (! StringHasNoText (str + 1)) {
              vnp->data.ptrvalue = StringSave (str + 1);
            }
          }
        }
      } else if (ch == '<') {
      } else if (ch != '\0') {
        vnp = ValNodeAdd (&list);
        if (vnp != NULL) {
          if (! StringHasNoText (str)) {
            vnp->data.ptrvalue = StringSave (str);
          }
        }
      }
      goOn = (Boolean) (goOn && (FileCacheGetString (&fc, str, sizeof (str)) != NULL));
    }
    if (list != NULL) {
      ptr = MergeValNodeStrings (list, inTable);
      if (inTable) {
        choice = TBL_FMT;
      } else {
        choice = TXT_FMT;
      }
      numItems++;
      vnp = ValNodeAdd (&(hfp->mainStrings));
      if (vnp != NULL) {
        vnp->choice = choice;
        vnp->data.ptrvalue = ptr;
      }
      /* ptr = MemFree (ptr); */
      list = ValNodeFreeData (list);
    }
    FileClose (fp);
    return TRUE;
  } else {
    return FALSE;
  }
}

static FonT GetHelpFontFromConfig (CharPtr param, FonT dfault)

{
  FonT  f;
  Char  str [128];

  f = dfault;
  if (GetSequinAppParam ("SCREEN", param, NULL, str, sizeof (str))) {
    f = Nlm_ParseFontEx (str, NULL);
  }
  if (f == NULL) {
    f = dfault;
  }
  return f;
}

static void SetupHelpFonts (void)

{
  if (IsJapanese()) {
#ifdef WIN_MAC
    /* Osaka is common font for Japanese on Macintosh. */
    hdgColFmt.font = Nlm_ParseFontEx ("Osaka,14,b", NULL);
    subColFmt.font = Nlm_ParseFontEx ("Osaka,10,b", NULL);
    txtColFmt.font = Nlm_ParseFontEx ("Osaka,10", NULL);
    lstColFmt.font = Nlm_ParseFontEx ("Osaka\x81\x7c\x93\x99\x95\x9d,12", NULL);
    tblColFmt.font = Nlm_ParseFontEx ("Osaka\x81\x7c\x93\x99\x95\x9d,12", NULL);
/*    hdgColFmt.font = GetResidentFont (Nlm_Minchou(14, STYLE_BOLD));	*/
/*    subColFmt.font = GetResidentFont (Nlm_Minchou(10, STYLE_BOLD));	*/
/*    txtColFmt.font = GetResidentFont (Nlm_Gothic(10, STYLE_REGULAR));	*/
/*    lstColFmt.font = GetResidentFont (Nlm_MinchouFixed(12, STYLE_REGULAR));	*/
/*    tblColFmt.font = GetResidentFont (Nlm_MinchouFixed(12, STYLE_REGULAR));	*/
#endif
#ifdef WIN_MSWIN
    hdgColFmt.font = Nlm_ParseFontEx ("\x82\x6c\x82\x72\x20\x82\x6f\x96\xbe\x92\xa9,14,b,Kanji", NULL);
    subColFmt.font = Nlm_ParseFontEx ("\x82\x6c\x82\x72\x20\x82\x6f\x96\xbe\x92\xa9,11,b,Kanji", NULL);
    txtColFmt.font = Nlm_ParseFontEx ("\x82\x6c\x82\x72\x20\x82\x6f\x83\x53\x83\x56\x83\x62\x83\x4e,11,,Kanji", NULL);
    lstColFmt.font = Nlm_ParseFontEx ("\x82\x6c\x82\x72\x20\x96\xbe\x92\xa9,12,f,Kanji", NULL);
    tblColFmt.font = Nlm_ParseFontEx ("\x82\x6c\x82\x72\x20\x96\xbe\x92\xa9,12,f,Kanji", NULL);
/*    hdgColFmt.font = GetResidentFont (Nlm_Minchou(14, STYLE_BOLD));	*/
/*    subColFmt.font = GetResidentFont (Nlm_Minchou(11, STYLE_BOLD));	*/
/*    txtColFmt.font = GetResidentFont (Nlm_Gothic(11, STYLE_REGULAR));	*/
/*    lstColFmt.font = GetResidentFont (Nlm_MinchouFixed(12, STYLE_REGULAR));	*/
/*    tblColFmt.font = GetResidentFont (Nlm_MinchouFixed(12, STYLE_REGULAR));	*/
#endif
#ifdef WIN_MOTIF
    hdgColFmt.font = ParseFont ("Times,18,b");
    subColFmt.font = ParseFont ("Helvetica,12,b");
    txtColFmt.font = ParseFont ("fixed,13");
    lstColFmt.font = programFont;
    tblColFmt.font = programFont;
#endif
  } else if (IsEnglish() || IsFrench() || IsGerman() || IsItalian() || IsSystemLang()) {
#ifdef WIN_MAC
    hdgColFmt.font = ParseFont ("Times,14,b");
    subColFmt.font = ParseFont ("Geneva,10,b");
    txtColFmt.font = ParseFont ("Geneva,10");
#endif /*WIN_MAC*/
#ifdef WIN_MSWIN
    hdgColFmt.font = ParseFont ("Times New Roman,14,b");
    subColFmt.font = ParseFont ("Arial,11,b");
    txtColFmt.font = ParseFont ("Times New Roman,11");
#endif /*WIN_MSWIN*/
#ifdef WIN_MOTIF
    hdgColFmt.font = ParseFont ("Times,18,b");
    subColFmt.font = ParseFont ("Helvetica,12,b");
    txtColFmt.font = ParseFont ("Times,14");
#endif /*WIN_MOTIF*/
    lstColFmt.font = programFont;
    tblColFmt.font = programFont;
  } else {
  /* above call to IsSystemLang should override this section */
    /* default system font is often native character set */
    /* because native character sets have a set of ascii letters,
       english help file can also be drawn as well as native letters. */
    hdgColFmt.font = systemFont;
    subColFmt.font = systemFont;
    txtColFmt.font = systemFont;
    lstColFmt.font = systemFont;
    tblColFmt.font = systemFont;
  }

  /* now allow override by sequin config file */
  hdgColFmt.font = GetHelpFontFromConfig ("HEADING", hdgColFmt.font);
  subColFmt.font = GetHelpFontFromConfig ("SUBHEAD", subColFmt.font);
  txtColFmt.font = GetHelpFontFromConfig ("TEXT", txtColFmt.font);
  lstColFmt.font = GetHelpFontFromConfig ("LIST", lstColFmt.font);
  tblColFmt.font = GetHelpFontFromConfig ("TABLE", tblColFmt.font);
}

static Boolean PopulateHelpForm (HelpFormPtr hfp)

{
  Int2        firstLine;
  Int2        firstShown;
  RecT        r;
  BaR         sb;
  Int4        startsAt;
  CharPtr     text;
  ValNodePtr  vnp;

  if (hfp == NULL || hfp->doc == NULL) return FALSE;
  Hide (hfp->list);
  Hide (hfp->doc);
  Update ();
  if (! GetScrlParams4 (hfp->doc, NULL, &firstShown, &firstLine)) {
    firstShown = 0;
    firstLine = 0;
  }
  sb = GetSlateVScrollBar ((SlatE) hfp->doc);
  Reset (hfp->list);
  Reset (hfp->doc);
  if (hfp->file == NULL || hfp->file [0] == '\0') return FALSE;
  if (hfp->mainStrings == NULL || hfp->indexStrings == NULL) return FALSE;
  ObjectRect (hfp->doc, &r);
  InsetRect (&r, 4, 4);
  lstColFmt.pixWidth = r.right - r.left;
  hdgColFmt.pixWidth = r.right - r.left;
  subColFmt.pixWidth = r.right - r.left;
  /*
  tblColFmt.pixWidth = screenRect.right - screenRect.left;
  */
  tblColFmt.pixWidth = r.right - r.left;
  tblColFmt.pixInset = 10;
  txtColFmt.pixWidth = r.right - r.left;
  txtColFmt.pixInset = 10;

  for (vnp = hfp->indexStrings; vnp != NULL; vnp = vnp->next) {
    lstColFmt.pixInset = vnp->choice;
    AppendText (hfp->list, vnp->data.ptrvalue, &lstParFmt, &lstColFmt, programFont);
  }

  for (vnp = hfp->mainStrings; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case TBL_FMT :
        AppendText (hfp->doc, vnp->data.ptrvalue, &tblParFmt, &tblColFmt, programFont);
        break;
      case TXT_FMT :
        AppendText (hfp->doc, vnp->data.ptrvalue, &txtParFmt, &txtColFmt, programFont);
        break;
      case HDG_FMT :
        AppendText (hfp->doc, vnp->data.ptrvalue, &hdgParFmt, &hdgColFmt, systemFont);
        break;
      case SUB_FMT :
        subColFmt.pixInset = 0;
        AppendText (hfp->doc, vnp->data.ptrvalue, &subParFmt, &subColFmt, programFont);
        break;
      default :
        subColFmt.pixInset = vnp->choice;
        AppendText (hfp->doc, vnp->data.ptrvalue, &subParFmt, &subColFmt, programFont);
        break;
    }
  }

  UpdateDocument (hfp->list, 0, 0);
  /*
  UpdateDocument (hfp->doc, 0, 0);
  */
  text = GetDocText (hfp->doc, firstShown, 0, 0);
  MemFree (text);
  AdjustDocScroll (hfp->doc);
  GetItemParams4 (hfp->doc, firstShown, &startsAt, NULL, NULL, NULL, NULL);
  CorrectBarValue (sb, startsAt + firstLine);
  Show (hfp->list);
  Show (hfp->doc);
  Update ();
  return TRUE;
}

static void RefreshHelpForm (ButtoN b)

{
  HelpFormPtr  hfp;

  hfp = (HelpFormPtr) GetObjectExtra (b);
  if (hfp == NULL) return;
  ParseHelpFile (hfp, TRUE);
  PopulateHelpForm (hfp);
  Update ();
}

static void CleanupHelpForm (GraphiC g, VoidPtr data)

{
  ResetHelpLists ((HelpFormPtr) data);
  StdCleanupFormProc (g, data);
}

static void HelpListNotify (DoC d, Int2 item, Int2 row, Int2 col, Boolean dblclick)

{
  HelpFormPtr   hfp;
  HelpIndexPtr  hip;
  BaR           sb;
  Int2          startsAt;
  ValNodePtr    vnp;

  hfp = (HelpFormPtr) GetObjectExtra (d);
  if (hfp == NULL || hfp->doc == NULL) return;
  if (item == 0 || row == 0 || col == 0) return;
  vnp = hfp->index;
  while (vnp != NULL && item > 1) {
    item--;
    vnp = vnp->next;
  }
  if (vnp != NULL) {
    hip = (HelpIndexPtr) vnp->data.ptrvalue;
    if (hip != NULL) {
      GetItemParams (hfp->doc, hip->item, &startsAt, NULL, NULL, NULL, NULL);
      ResetClip ();
      sb = GetSlateVScrollBar ((SlatE) hfp->doc);
      SetValue (sb, startsAt);
      Update ();
    }
  }
}

static void ResizeHelpForm (WindoW w)

{
  Int2         delta;
  Int2         diff;
  Int2         gap;
  Int2         height;
  HelpFormPtr  hfp;
  RecT         r;
  RecT         s;
  RecT         t;
  Int2         width;

  hfp = (HelpFormPtr) GetObjectExtra (w);
  if (hfp != NULL) {
    ObjectRect (w, &r);
    width = r.right - r.left;
    height = r.bottom - r.top;
    GetPosition (hfp->doc, &s);
    GetPosition (hfp->dismissGrp, &t);
    diff = t.bottom - t.top;
    gap = t.top - s.bottom;
    t.bottom = height - s.left;
    t.top = t.bottom - diff;
    delta = (width - t.right - t.left) / 2;
    t.left += delta;
    t.right += delta;
    s.right = width - s.left;
    /*
    s.bottom = height - s.left;
    */
    s.bottom = t.top - gap;
    SetPosition (hfp->dismissGrp, &t);
    SetPosition (hfp->doc, &s);
    AdjustPrnt (hfp->doc, &s, FALSE);
    PopulateHelpForm (hfp);
    Update ();
  }
}

static Boolean ScrollToTextInDoc (DoC d, Int2 item, CharPtr text)

{
  Pointer  data;
  Boolean  found;
  RecT     rct;
  BaR      sb;
  Int4     startsAt;
  CharPtr  str;
  WindoW   tempPort;

  found = FALSE;
  GetItemParams4 (d, item, &startsAt, NULL, NULL, NULL, &data);
  str = (CharPtr) data;
  if (StringISearch (str, text) != NULL) {
    found = TRUE;
  }
  if (found) {
    tempPort = SavePort (d);
    Select (d);
    sb = GetSlateVScrollBar ((SlatE) d);
    CorrectBarValue (sb, startsAt);
    ObjectRect (d, &rct);
    InsetRect (&rct, 4, 4);
    InsetRect (&rct, -1, -1);
    InvalRect (&rct);
    RestorePort (tempPort);
    Update ();
  }
  return found;
}

static void FindHelpBtnProc (ButtoN b)

{
  Int2         firstShown;
  HelpFormPtr  hfp;
  Int2         i;
  Int2         numItems;
  Char         str [256];

  hfp = (HelpFormPtr) GetObjectExtra (b);
  if (hfp == NULL) return;
  GetTitle (hfp->findTxt, str, sizeof (str) - 1);
  GetDocParams (hfp->doc, &numItems, NULL);
  if (GetScrlParams (hfp->doc, NULL, &firstShown, NULL)) {
    for (i = firstShown + 1; i <= numItems; i++) {
      if (ScrollToTextInDoc (hfp->doc, i, str)) {
        return;
      }
    }
    for (i = 1; i < firstShown; i++) {
      if (ScrollToTextInDoc (hfp->doc, i, str)) {
        return;
      }
    }
  }
}

static void FindHelpTextProc (TexT t)

{
  HelpFormPtr  hfp;

  hfp = (HelpFormPtr) GetObjectExtra (t);
  if (hfp == NULL) return;
  if (TextLength (t) > 0) {
    SafeEnable (hfp->findBtn);
  } else {
    SafeDisable (hfp->findBtn);
  }
}

static void HelpFormMessage (ForM f, Int2 mssg)

{
  FILE         *fp;
  HelpFormPtr  hfp;
  Char         path [PATH_MAX];

  hfp = (HelpFormPtr) GetObjectExtra (f);
  if (hfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Hide (f);
        break;
      case VIB_MSG_PRINT :
        PrintDocument (hfp->doc);
        break;
      case VIB_MSG_EXPORT :
        if (GetOutputFileName (path, sizeof (path), "help.txt")) {
          WatchCursor ();
#ifdef WIN_MAC
          fp = FileOpen (path, "r");
          if (fp != NULL) {
            FileClose (fp);
          } else {
            FileCreate (path, "TEXT", "ttxt");
          }
#endif
          fp = FileOpen (path, "w");
          if (fp != NULL) {
            SaveDocument (hfp->doc, fp);
            FileClose (fp);
          }
          ArrowCursor ();
        }
        break;
      default :
        if (hfp->appmessage != NULL) {
          hfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void HelpFormActivate (WindoW w)

{
  IteM         exportItm;
  HelpFormPtr  hfp;

  hfp = (HelpFormPtr) GetObjectExtra (w);
  if (hfp != NULL) {
    if (hfp->activate != NULL) {
      hfp->activate (w);
    }
    exportItm = FindFormMenuItem ((BaseFormPtr) hfp, VIB_MSG_EXPORT);
    SafeSetTitle (exportItm, "Export Help...");
  }
}

extern ForM CreateHelpForm (Int2 left, Int2 top, CharPtr title,
                            CharPtr file, BtnActnProc closeForm,
                            WndActnProc activateForm)

{
  ButtoN             b;
  GrouP              c;
  GrouP              h;
  Int2               height;
  HelpFormPtr        hfp;
  StdEditorProcsPtr  sepp;
  WindoW             w;
#ifndef WIN_MAC
  MenU               m;
#endif

  w = NULL;
  hfp = MemNew (sizeof (HelpForm));
  if (hfp != NULL) {
    w = DocumentWindow (left, top, -10, -10, title,
                        StdSendCloseWindowMessageProc, ResizeHelpForm);
    SetObjectExtra (w, hfp, CleanupHelpForm);
    hfp->form = (ForM) w;
    hfp->formmessage = HelpFormMessage;

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      hfp->appmessage = sepp->handleMessages;
    }

    StringNCpy_0 (hfp->file, file, sizeof (hfp->file));

#ifndef WIN_MAC
    m = PulldownMenu (w, "File");
    FormCommandItem (m, "Export", (BaseFormPtr) hfp, VIB_MSG_EXPORT);
    SeparatorItem (m);
#ifdef WIN_MSWIN
    FormCommandItem (m, "Print", (BaseFormPtr) hfp, VIB_MSG_PRINT);
    SeparatorItem (m);
#endif
    FormCommandItem (m, "Close", (BaseFormPtr) hfp, VIB_MSG_CLOSE);
#endif

    c = HiddenGroup (w, 4, 0, NULL);
    StaticPrompt (c, "Find", 0, dialogTextHeight, programFont, 'l');
    hfp->findTxt = DialogText (c, "", 20, FindHelpTextProc);
    SetObjectExtra (hfp->findTxt, hfp, NULL);
    hfp->findBtn = PushButton (c, "Find", FindHelpBtnProc);
    SetObjectExtra (hfp->findBtn, hfp, NULL);
    Disable (hfp->findBtn);

    h = HiddenGroup (w, -1, 0, NULL);
    SelectFont (programFont);
    height = LineHeight ();
    SelectFont (systemFont);
    hfp->list = DocumentPanel (h, stdCharWidth * 28, height * 5);
    SetDocAutoAdjust (hfp->list, FALSE);
    SetObjectExtra (hfp->list, hfp, NULL);
    SetDocNotify (hfp->list, HelpListNotify);

    hfp->doc = DocumentPanel (h, stdCharWidth * 33, stdLineHeight * 20);
    SetDocAutoAdjust (hfp->doc, FALSE);
    SetObjectExtra (hfp->doc, hfp, NULL);

    hfp->dismissGrp = HiddenGroup (w, 2, 0, NULL);
    PushButton (hfp->dismissGrp, "Dismiss", closeForm);
    if (indexerVersion) {
      b = PushButton (hfp->dismissGrp, "Refresh", RefreshHelpForm);
      SetObjectExtra (b, hfp, NULL);
    }

    AlignObjects (ALIGN_CENTER, (HANDLE) hfp->doc, (HANDLE) hfp->dismissGrp, NULL);

    RealizeWindow (w);

    SetupHelpFonts ();

    if (activateForm != NULL) {
      hfp->activate = activateForm;
    } else {
      if (sepp != NULL) {
        hfp->activate = sepp->activateForm;
      }
    }
    SetActivate (w, HelpFormActivate);
    HelpFormActivate ((WindoW) hfp->form);

    if (ParseHelpFile (hfp, FALSE)) {
      if (! PopulateHelpForm (hfp)) {
        w = Remove (w);
      }
    } else {
      w = Remove (w);
    }
  }
  return (ForM) w;
}

extern void SendHelpScrollMessage (ForM f, CharPtr heading, CharPtr section)

{
  HelpFormPtr   hfp;
  HelpIndexPtr  hip;
  BaR           sb;
  Int2          startsAt;
  CharPtr       str;
  Char          txt [256];
  Boolean       useBoth;
  Boolean       useHeading;
  Boolean       useSection;
  ValNodePtr    vnp;

  hfp = (HelpFormPtr) GetObjectExtra (f);
  if (hfp != NULL) {
    vnp = hfp->index;
    txt [0] = '\0';
    useBoth = FALSE;
    useHeading = FALSE;
    useSection = FALSE;
    if (heading != NULL && *heading != '\0' && section != NULL && *section != '\0') {
      useBoth = TRUE;
      if (StringLen (heading) + StringLen (section) < sizeof (txt) - 2) {
        StringCpy (txt, heading);
        StringCat (txt, "|");
        StringCat (txt, section);
      }
    } else if (heading != NULL && *heading != '\0') {
      useHeading = TRUE;
      StringNCpy_0 (txt, heading, sizeof (txt));
    } else if (section != NULL && *section != '\0') {
      useSection = TRUE;
      StringNCpy_0 (txt, section, sizeof (txt));
    }
    while (vnp != NULL) {
      hip = (HelpIndexPtr) vnp->data.ptrvalue;
      if (hip != NULL) {
        if (useBoth) {
          str = hip->combined;
        } else if (useHeading) {
          str = hip->heading;
        } else if (useSection) {
          str = hip->section;
        } else {
          str = NULL;
        }
        if (str != NULL && StringICmp (txt, str) == 0) {
          GetItemParams (hfp->doc, hip->item, &startsAt, NULL, NULL, NULL, NULL);
          sb = GetSlateVScrollBar ((SlatE) hfp->doc);
          ResetClip ();
          SetValue (sb, startsAt);
          Update ();
          return;
        }
      }
      vnp = vnp->next;
    }
  }
}

typedef struct applycdsframe 
{
  FEATURE_FORM_BLOCK
  DialoG constraint_dlg;
  PopuP  current_frame_popup;
  PopuP  new_frame_popup;
  ButtoN retranslate_btn;
  
  Int4 current_frame_flag;
  Int4 new_frame;
  Boolean retranslate_flag;
  ConstraintChoiceSetPtr constraint;
  LogInfoPtr lip;
} ApplyCDSFrameData, PNTR ApplyCDSFramePtr;

static void ApplyCDSFrameCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CdRegionPtr       crp;
  ApplyCDSFramePtr acfp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION 
      || (acfp = (ApplyCDSFramePtr) userdata) == NULL
      || !DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, sfp, acfp->constraint))
  {
    return;
  }
  
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL)
  {
    crp = CdRegionNew ();
    sfp->data.value.ptrvalue = crp;
  }
  
  if (acfp->current_frame_flag == 4 || crp->frame == acfp->current_frame_flag)
  {
    if (acfp->new_frame == 4) {
      if (!SetBestFrameByLocation (sfp)) {
        LogCDSAmbiguousFrame (acfp->lip, sfp);
      }
    } else {
      crp->frame = acfp->new_frame;
    }
    if (acfp->retranslate_flag)
    {
      RetranslateOneCDS (sfp, acfp->input_entityID, TRUE, FALSE);
    }
  }
}

static void DoApplyCDSFrame (ButtoN b)
{
  ApplyCDSFramePtr acfp;
  SeqEntryPtr      sep;
  
  acfp = (ApplyCDSFramePtr) GetObjectExtra (b);
  if (acfp == NULL)
  {
    return;
  }
  
  sep = GetTopSeqEntryForEntityID (acfp->input_entityID);
  if (sep == NULL)
  {
    return;
  }
  
  WatchCursor ();
  Update ();
  
  acfp->lip = OpenLog ("Ambiguous frames");
  acfp->constraint = DialogToPointer (acfp->constraint_dlg);
  
  acfp->current_frame_flag = GetValue (acfp->current_frame_popup);
  acfp->new_frame = GetValue (acfp->new_frame_popup);
  acfp->retranslate_flag = GetStatus (acfp->retranslate_btn);

  VisitFeaturesInSep (sep, acfp, ApplyCDSFrameCallback);

  acfp->constraint = ConstraintChoiceSetFree (acfp->constraint);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (acfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, acfp->input_entityID, 0, 0);  
  Remove (acfp->form);
  CloseLog (acfp->lip);
  acfp->lip = FreeLog (acfp->lip);

}

extern void ApplyCDSFrame (IteM i)
{
  BaseFormPtr        bfp;
  GrouP              g, h, c;
  ApplyCDSFramePtr   acfp;
  WindoW             w;
  ButtoN             b;

  /* Get current sequence */

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;
  
  acfp = (ApplyCDSFramePtr) MemNew (sizeof (ApplyCDSFrameData));
  if (acfp == NULL)
  {
    return;
  }
  acfp->lip = NULL;

  w = FixedWindow (-50, -33, -20, -10, "Apply CDS Frame",	StdCloseWindowProc);
  SetObjectExtra (w, acfp, StdCleanupFormProc);
  acfp->form = (ForM) w;

  acfp->input_entityID = bfp->input_entityID;
  acfp->input_itemID = bfp->input_itemID;

  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);

  h = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  StaticPrompt (h, "Set Coding Region frame to:", 0, dialogTextHeight,
                programFont, 'l');
  acfp->new_frame_popup = PopupList (h, TRUE, NULL);
  PopupItem (acfp->new_frame_popup, "1");
  PopupItem (acfp->new_frame_popup, "2");
  PopupItem (acfp->new_frame_popup, "3");
  PopupItem (acfp->new_frame_popup, "Best");
  SetValue (acfp->new_frame_popup, 1);

  StaticPrompt (h, "Where Coding Region frame is:", 0, dialogTextHeight,
                programFont, 'l');
  acfp->current_frame_popup = PopupList (h, TRUE, NULL);
  PopupItem (acfp->current_frame_popup, "1");
  PopupItem (acfp->current_frame_popup, "2");
  PopupItem (acfp->current_frame_popup, "3");
  PopupItem (acfp->current_frame_popup, "Any frame");
  
  SetValue (acfp->current_frame_popup, 4);
  
  acfp->constraint_dlg = ComplexConstraintDialog (g, NULL, NULL);
  ChangeComplexConstraintFieldType (acfp->constraint_dlg, FieldType_cds_gene_prot, NULL, Macro_feature_type_cds);

  acfp->retranslate_btn = CheckBox (g, "Retranslate adjusted coding regions", NULL);
  SetStatus (acfp->retranslate_btn, TRUE);
  
  c = HiddenGroup (g, 2, 0, NULL);
  b = DefaultButton(c, "Accept", DoApplyCDSFrame);
  SetObjectExtra(b, acfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) h, 
                              (HANDLE) acfp->constraint,
                              (HANDLE) acfp->retranslate_btn,
                              (HANDLE) c, NULL);
  RealizeWindow(w);
  Show(w);
  Update();
  
}

extern ValNodePtr FreeSeqIdList (ValNodePtr id_list)
{
  SeqIdPtr sip;
  ValNodePtr vnp;
  
  for (vnp = id_list; vnp != NULL; vnp = vnp->next) {
       sip = vnp->data.ptrvalue;
       sip = SeqIdFree (sip);
     vnp->data.ptrvalue = NULL;
  }
  id_list = ValNodeFree (id_list);
  return NULL;
}


static ValNodePtr SplitByCommasAndSpaces (CharPtr list)
{
  CharPtr cp;
  Char    ch;
  ValNodePtr str_list = NULL;
  Int4       len;

  if (StringHasNoText (list)) {
    return NULL;
  }

  /* skip any leading blanks */
  list += StringSpn (list, " ,\t");
  len = StringCSpn (list, " ,\t");
  while (len != 0) {
    cp = list + len;
    ch = *cp;
    *cp = 0;
    ValNodeAddPointer (&str_list, 0, StringSave (list));
    *cp = ch;
    list = cp + StringSpn (cp, " ,\t");
    len = StringCSpn (list, " ,\t");
  }
  return str_list;
}


static Boolean ParseIDStr (CharPtr id_str, CharPtr PNTR prefix, Int4Ptr num, Int4Ptr num_len)
{
  CharPtr start_number;
  Char    ch;

  if (StringHasNoText (id_str) || prefix == NULL || num == NULL || num_len == NULL) {
    return FALSE;
  }

  start_number = id_str + StringLen (id_str);
  while (start_number > id_str && isdigit (*(start_number - 1))) {
    start_number--;
  }
  *num_len = StringLen (start_number);
  if (num_len == 0) {
    return FALSE;
  }
  *num = atoi (start_number);
  if (start_number == id_str) {
    *prefix = NULL;
  } else {
    ch = *start_number;
    *start_number = 0;
    *prefix = StringSave (id_str);
    *start_number = ch;
  }
  return TRUE;
}

static ValNodePtr ExpandAccessionRanges (CharPtr list_str, SeqEntryPtr sep)
{
  ValNodePtr id_list = NULL;
  CharPtr    cp;
  Char       ch;
  SeqIdPtr   sip;
  CharPtr    last_prefix = NULL, this_prefix = NULL, tmpstr;
  Int4       num_len = 0, num_len2;
  Int4       range_start, range_end, sw;
  
  if (StringHasNoText (list_str)) {
    return NULL;
  } else if ((cp = StringChr (list_str, '-')) == NULL) {
    sip = CreateSeqIdFromText(list_str, sep);      
    if (sip == NULL) {/* Bad SeqId string */
      Message (MSG_ERROR, "Unable to parse %s as ID\n", list_str);
      return NULL;
    } else {
      ValNodeAddPointer (&id_list, 0, sip);
      return id_list;
    }
  } else {
    ch = *cp;
    *cp = 0;
    if (!ParseIDStr (list_str, &this_prefix, &range_start, &num_len)) {
      *cp = ch;
      Message (MSG_ERROR, "Unable to parse %s", list_str);
      return NULL;
    }
    *cp = ch;
    if (!ParseIDStr (cp + 1, &last_prefix, &range_end, &num_len2)) {
      Message (MSG_ERROR, "Unable to parse %s", list_str);
      this_prefix = MemFree (this_prefix);
      return NULL;
    }
    if (last_prefix != NULL && StringCmp (this_prefix, last_prefix) != 0) {
      Message (MSG_ERROR, "Unable to parse range for %s", list_str);
      this_prefix = MemFree (this_prefix);
      last_prefix = MemFree (last_prefix);
      return NULL;
    }

    if (range_start > range_end) {
      sw = range_start;
      range_start = range_end;
      range_end = sw;
    }
    tmpstr = (CharPtr) MemNew (sizeof (Char) * StringLen (this_prefix) + 15);
    for (sw = range_start; sw <= range_end; sw++) {
      sprintf (tmpstr, "%s%0*d", this_prefix == NULL ? "" : this_prefix, num_len, sw);
      sip = CreateSeqIdFromText(tmpstr, sep);
      if (sip == NULL) {
        if (sw == range_start || sw == range_end) {
          Message (MSG_ERROR, "Unable to parse range for %s", list_str);
          id_list = FreeSeqIdList (id_list);
          break;
        }
      } else {
        ValNodeAddPointer (&id_list, 0, sip);
      }
    }
    tmpstr = MemFree (tmpstr);
  }

  this_prefix = MemFree (this_prefix);  
  last_prefix = MemFree (last_prefix);
  
  return id_list;  
}


extern ValNodePtr ParseAccessionNumberListFromString (CharPtr list_str, SeqEntryPtr sep)
{
  ValNodePtr token_list, vnp, id_list = NULL, tmp;

  if (StringHasNoText (list_str)) {
      Message (MSG_ERROR, "No accession numbers listed!");
      return NULL;
  }

  token_list = SplitByCommasAndSpaces (list_str);
  for (vnp = token_list; vnp != NULL; vnp = vnp->next) {
    tmp = ExpandAccessionRanges (vnp->data.ptrvalue, sep);
    if (tmp == NULL) {
      id_list = FreeSeqIdList (id_list);
      break;
    } else {
      ValNodeLink (&id_list, tmp);
    }
  }
  token_list = ValNodeFreeData (token_list);
  return id_list;
}


extern Int2 LIBCALLBACK CopyDescriptorToList (Pointer data)

{
  OMProcControlPtr      ompcp;
  SeqEntryPtr           sep = NULL, oldscope;
  SeqDescrPtr           sdp, last_sdp, new_sdp, next_sdp = NULL;
  PrompT                p;
  TexT                  accession_list_txt;
  GrouP                 c;
  GrouP                 g;
  WindoW                w;
  ModalAcceptCancelData acd;
  ButtoN                b;
  ValNodePtr            id_list, vnp;
  CharPtr               str;
  BioseqPtr             bsp;
  BioseqSetPtr          bssp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC:
        sdp = ompcp->input_data;
        next_sdp = sdp->next;
        break;
    case OBJ_BIOSEQ:
        bsp = ompcp->input_data;
        sdp = bsp->descr;
        break;
    case OBJ_BIOSEQSET:
        bssp = ompcp->input_data;
        sdp = bssp->descr;
        break;
    default:
        Message (MSG_ERROR, "Must select a descriptor, a sequence, or a set!");
        return OM_MSG_RET_ERROR;
        break;
  }
  
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  
  w = MovableModalWindow (-50, -33, -10, -10, "Copy Descriptor to Sequences", NULL);
  SetGroupSpacing (w, 10, 10);

  g = HiddenGroup (w, -1, 0, NULL);
  p = StaticPrompt (g, "Enter comma-delimited Sequence ID list", 0, stdLineHeight, programFont, 'l');
  accession_list_txt = DialogText (g, "", 25, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) accession_list_txt, NULL);
  
  c = HiddenGroup (w, 4, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  b = DefaultButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);

  Select (accession_list_txt);
  Show (w);
  Select (w);
  Update ();

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
      str = SaveStringFromText (accession_list_txt);
      id_list = ParseAccessionNumberListFromString (str, sep);
      if (id_list == NULL) {
        acd.accepted = FALSE;
      } else {
        Hide (w);
        WatchCursor ();
        Update ();
        oldscope = SeqEntrySetScope (sep);
        if (next_sdp != NULL) {
          sdp->next = NULL;
        }
        for (vnp = id_list; vnp != NULL; vnp = vnp->next)
        {
          bsp = BioseqFind (vnp->data.ptrvalue);
          if (bsp != NULL) {
            new_sdp = (SeqDescrPtr) AsnIoMemCopy (sdp, (AsnReadFunc) SeqDescrAsnRead,
                                                  (AsnWriteFunc) SeqDescrAsnWrite);
            last_sdp = bsp->descr;
            while (last_sdp != NULL && last_sdp->next != NULL) {
              last_sdp = last_sdp->next;
            }
            if (last_sdp == NULL) {
              bsp->descr = new_sdp;
            } else {
              last_sdp->next = new_sdp;
            }
          }
        }
        if (next_sdp != NULL) {
          sdp->next = next_sdp;
        }
        
        id_list = FreeSeqIdList (id_list);
        SeqEntrySetScope (oldscope);
        ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
        ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
        ArrowCursor ();
        Update ();
      }
    }
  }
  Remove (w);
  return OM_MSG_RET_DONE;
}


static Boolean IsSeqAnnotPtrtRNAScan (SeqAnnotPtr sap)
{
  AnnotDescrPtr desc;

  if (sap == NULL || sap->type != 1 
      || sap->desc == NULL)
  {
    return FALSE;
  }

  for (desc = sap->desc; desc != NULL; desc = desc->next)
  {
    if (desc->choice == Annot_descr_name
        && StringCmp (desc->data.ptrvalue, "tRNAscan-SE") == 0)
    {
      return TRUE;
    }
  }
  
  return FALSE;
}  

static void FindtRNAScanAnnot (SeqAnnotPtr sap, Pointer userdata)
{
  if (userdata != NULL && IsSeqAnnotPtrtRNAScan(sap))
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQANNOT, sap);
  }
}


static void FindOthertRNAAnnot (SeqAnnotPtr sap, Pointer userdata)
{
  SeqFeatPtr sfp;
  Boolean    has_trna = FALSE;

  if (sap == NULL || userdata == NULL || sap->type != 1 
      || IsSeqAnnotPtrtRNAScan (sap))
  {
    return;
  }
  for (sfp = sap->data; sfp != NULL && !has_trna; sfp = sfp->next)
  {
    if (sfp->idx.subtype == FEATDEF_tRNA)
    {
      has_trna = TRUE;
    }
  }
  if (has_trna)
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_SEQANNOT, sap);
  }
}

typedef struct trnaupdatepair {
  SeqFeatPtr sfp_trna;
  SeqFeatPtr sfp_other;
} tRNAUpdatePairData, PNTR tRNAUpdatePairPtr;


static Boolean EndpointsOkFortRNAUpdate (Int4 p1, Int4 p2)
{
  if ((p1 - p2 < 4 && p1 - p2 >= 0) 
      || (p2 - p1 < 4 && p2 - p1 >= 0))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static ValNodePtr FindtRNAUpdatePairs (SeqAnnotPtr tRNAScan, SeqAnnotPtr otherAnnot, LogInfoPtr lip)
{
  ValNodePtr pair_list = NULL;
  tRNAUpdatePairPtr t;
  SeqFeatPtr        sfp_trna, sfp_other;
  RnaRefPtr         scan_rrp, other_rrp;
  tRNAPtr           scan_trp, other_trp;
  Boolean           any_match;
  Int4              start_scan, stop_scan, start_other, stop_other;

  if (!IsSeqAnnotPtrtRNAScan(tRNAScan)
      || otherAnnot == NULL 
      || otherAnnot->type != 1 
      || IsSeqAnnotPtrtRNAScan(otherAnnot))
  {
    return NULL;
  }

  for (sfp_trna = tRNAScan->data; sfp_trna != NULL; sfp_trna = sfp_trna->next)
  {
    any_match = FALSE;
    if (sfp_trna->idx.subtype != FEATDEF_tRNA) continue;
    scan_rrp = sfp_trna->data.value.ptrvalue;
    if (scan_rrp == NULL) continue;
    scan_trp = scan_rrp->ext.value.ptrvalue;
    if (scan_trp == NULL) continue;
    start_scan = SeqLocStart (sfp_trna->location);
    stop_scan = SeqLocStop (sfp_trna->location);
    for (sfp_other = otherAnnot->data; sfp_other != NULL; sfp_other = sfp_other->next)
    {
      if (sfp_other->idx.subtype != FEATDEF_tRNA) continue;
      other_rrp = sfp_other->data.value.ptrvalue;
      if (other_rrp == NULL) continue;
      other_trp = other_rrp->ext.value.ptrvalue;
      if (other_trp == NULL) continue;
      /* make sure amino acids match */
      if (scan_trp->aatype != other_trp->aatype
          || scan_trp->aa != other_trp->aa) continue;
      start_other = SeqLocStart (sfp_other->location);
      stop_other = SeqLocStop (sfp_other->location);
      if (EndpointsOkFortRNAUpdate (start_scan, start_other) && EndpointsOkFortRNAUpdate (stop_scan, stop_other))
      {
        t = (tRNAUpdatePairPtr) MemNew (sizeof (tRNAUpdatePairData));
        t->sfp_trna = sfp_trna;
        t->sfp_other = sfp_other;
        ValNodeAddPointer (&pair_list, 0, t);
        any_match = TRUE;
      }
    }
    if (!any_match && lip != NULL && lip->fp != NULL)
    {
      fprintf (lip->fp, "No match for %c %d-%d\n", scan_trp->aa, start_scan, stop_scan);
      lip->data_in_log = TRUE;
    }
  }
  return pair_list;
}


static void UpdateOnetRNAScanPair (tRNAUpdatePairPtr t, LogInfoPtr lip)
{
  RnaRefPtr rrp_scan, rrp_other;
  tRNAPtr   trp_scan, trp_other;
  SeqMgrFeatContext context;
  CharPtr           location;
  SeqFeatPtr        sfp;  

  if (t == NULL || t->sfp_trna == NULL || t->sfp_other == NULL
      || t->sfp_trna->idx.subtype != FEATDEF_tRNA
      || t->sfp_other->idx.subtype != FEATDEF_tRNA
      || t->sfp_trna->data.value.ptrvalue == NULL
      || t->sfp_other->data.value.ptrvalue == NULL)
  {
    return;
  }

  rrp_scan = (RnaRefPtr) t->sfp_trna->data.value.ptrvalue;
  rrp_other = (RnaRefPtr) t->sfp_other->data.value.ptrvalue;
  if (rrp_scan == NULL || rrp_other == NULL) return;

  trp_scan = rrp_scan->ext.value.ptrvalue;
  trp_other = rrp_other->ext.value.ptrvalue;
  if (trp_scan == NULL || trp_other == NULL) return;

  /* update anticodon location */
  if (trp_scan->anticodon != NULL && trp_other->anticodon == NULL)
  {
    trp_other->anticodon = SeqLocCopy (trp_scan->anticodon);
    if (lip != NULL && lip->fp != NULL)
    {
      sfp = SeqMgrGetDesiredFeature (t->sfp_other->idx.entityID, NULL, 0, 0, t->sfp_other, &context);
      location = SeqLocPrintUseBestID (sfp->location);
      fprintf (lip->fp, "Updated anticodon for %s to %s\n", context.label, location);
      location = MemFree (location);
      lip->data_in_log = TRUE;
    }
    /* remove tRNA scan */
    t->sfp_trna->idx.deleteme = TRUE;
  }
  else if (SeqLocCompare (trp_scan->anticodon, trp_other->anticodon) == SLC_A_EQ_B)
  {
    t->sfp_trna->idx.deleteme = TRUE;
  }
  else if (t->sfp_trna->pseudo)
  {
    /* delete because scan is pseudo */
    t->sfp_trna->idx.deleteme = TRUE;
  }
  else
  {
    sfp = SeqMgrGetDesiredFeature (t->sfp_other->idx.entityID, NULL, 0, 0, t->sfp_other, &context);
    fprintf (lip->fp, "Mismatched anticodon for %s\n", context.label);
    lip->data_in_log = TRUE;
  }
}

extern void tRNAScanUpdate (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  ValNodePtr         scan_annot = NULL, other_annot = NULL;
  ValNodePtr         scan_vnp, other_vnp, pair_vnp;
  ValNodePtr         pair_list = NULL;
  LogInfoPtr         lip;

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

  VisitAnnotsInSep (sep, &scan_annot, FindtRNAScanAnnot);

  if (scan_annot == NULL)
  {
    Message (MSG_ERROR, "No tRNAscan feature tables found!");
    return;
  }

  VisitAnnotsInSep (sep, &other_annot, FindOthertRNAAnnot);
  if (other_annot == NULL)
  {
    Message (MSG_ERROR, "No feature tables with tRNAs found!");
    scan_annot = ValNodeFree (scan_annot);
    return;
  }

  lip = OpenLog ("Updated tRNAs");  
  for (scan_vnp = scan_annot; scan_vnp != NULL; scan_vnp = scan_vnp->next)
  {
    for (other_vnp = other_annot; other_vnp != NULL; other_vnp = other_vnp->next)
    {
      ValNodeLink (&pair_list, FindtRNAUpdatePairs (scan_vnp->data.ptrvalue, other_vnp->data.ptrvalue, lip));
    }
  }

  for (pair_vnp = pair_list; pair_vnp != NULL; pair_vnp = pair_vnp->next)
  {
    if (pair_vnp->choice == 0)
    {
      UpdateOnetRNAScanPair (pair_vnp->data.ptrvalue, lip);
    }
  }  
  CloseLog (lip);
  lip = FreeLog (lip);
  scan_annot = ValNodeFree (scan_annot);
  other_annot = ValNodeFree (other_annot);
  pair_list = ValNodeFreeData (pair_list);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

typedef struct resolvefeatureoverlaps {
  Uint1 trim_type;
  ValNodePtr trim_constraint;
  Uint1 intersect_type;
  ValNodePtr intersect_constraint;
  LogInfoPtr lip;
} ResolveFeatureOverlapsData, PNTR ResolveFeatureOverlapsPtr;

static Boolean TrimLocationForIntersectingFeatures (SeqLocPtr PNTR p_slp, ValNodePtr feat_list)
{
  SeqLocPtr  slp;
  SeqFeatPtr sfp;
  ValNodePtr vnp;
  Int4       orig_left, new_left, feat_left;
  Int4       orig_right, new_right, feat_right;
  Int4       tmp;
  SeqIdPtr   sip;
  Boolean    changed = FALSE, end_changed;

  if (p_slp == NULL || *p_slp == NULL || feat_list == NULL) return FALSE;

  slp = *p_slp;
  orig_left = SeqLocStart (slp);
  orig_right = SeqLocStop (slp);
  if (orig_left > orig_right)
  {
    tmp = orig_left;
    orig_left = orig_right;
    orig_right = tmp;
  }

  new_left = orig_left;
  new_right = orig_right;

  for (vnp = feat_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    feat_left = SeqLocStart (sfp->location);
    feat_right = SeqLocStop (sfp->location);
    if (feat_left > feat_right)
    {
      tmp = feat_left;
      feat_left = feat_right;
      feat_right = tmp;
    }
    if (feat_left <= orig_left && feat_right >= orig_left && feat_right <= orig_right)
    {
      /* trim on left */
      if (new_left < feat_right + 1)
      {
        new_left = feat_right + 1;
      }
    }
    else if (feat_left <= orig_right && feat_right >= orig_right && feat_left >= orig_left)
    {
      /* trim on right */
      if (new_right > feat_left - 1)
      {
        new_right = feat_left - 1;
      }
    }
  }

  if (new_left >= new_right)
  {
    *p_slp = SeqLocFree (*p_slp);
  }
  else
  {
    if (new_left > orig_left)
    {
      sip = SeqLocId (slp);
      end_changed = FALSE;
      *p_slp = SeqLocDelete (*p_slp, sip, orig_left, new_left - 1, FALSE, &end_changed);
      changed |= end_changed;
    }
    if (new_right < orig_right)
    {
      sip = SeqLocId (slp);
      end_changed = FALSE;
      *p_slp = SeqLocDelete (*p_slp, sip, new_right + 1, orig_right, FALSE, &end_changed);
      changed |= end_changed;
    }
  }
  return changed;
}


static void ResolveFeatureOverlapsBioseqCallback (BioseqPtr bsp, Pointer data)
{
  ResolveFeatureOverlapsPtr r;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  ValNodePtr        overlap_list;
  CharPtr           feat_desc;
  ValNode           vn;

  r = (ResolveFeatureOverlapsPtr) data;
  if (bsp == NULL || r == NULL) return;

  MemSet (&vn, 0, sizeof (ValNode));
  vn.next = NULL;
  vn.choice = OBJ_SEQFEAT;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, r->trim_type, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, r->trim_type, &fcontext))
  {
    vn.data.ptrvalue = sfp;
    if (!DoesObjectMatchConstraintChoiceSet(OBJ_SEQFEAT, sfp, r->trim_constraint)) {
      continue;
    }
    feat_desc = GetDiscrepancyItemText (&vn);
    overlap_list = ListFeaturesOverlappingLocationEx (bsp, sfp->location, 0, r->intersect_type, r->intersect_constraint);
    if (TrimLocationForIntersectingFeatures (&(sfp->location), overlap_list))
    {
      if (sfp->location == NULL)
      {
        fprintf (r->lip->fp, "Feature completely overlapped and removed:\n%s\n", feat_desc);
        sfp->idx.deleteme = TRUE;
        r->lip->data_in_log = TRUE;
      }
      else
      {
        fprintf (r->lip->fp, "Trimmed Feature: %sNew location:", feat_desc);
        LogTrimmedLocation (r->lip, sfp->location);
        fprintf (r->lip->fp, "\n");
      }
    }
    feat_desc = MemFree (feat_desc);
    overlap_list = ValNodeFree (overlap_list);
  }  
}


typedef struct resolvefeatureoverlapsform {
  FORM_MESSAGE_BLOCK
  DialoG trim_type;
  DialoG trim_constraint;
  DialoG intersect_type;
  DialoG intersect_constraint;
  DialoG accept_cancel;

} ResolveFeatureOverlapsFormData, PNTR ResolveFeatureOverlapsFormPtr;


static void SetResolveFeatureOverlapsFormAccept (Pointer data)
{
  ResolveFeatureOverlapsFormPtr dlg;
  ValNodePtr                    vnp1, vnp2;

  dlg = (ResolveFeatureOverlapsFormPtr) data;
  if (dlg == NULL) return;

  vnp1 = DialogToPointer (dlg->trim_type);
  vnp2 = DialogToPointer (dlg->intersect_type);
  if (vnp1 == NULL || vnp2 == NULL)
  {
    DisableAcceptCancelDialogAccept (dlg->accept_cancel);
  }
  else
  {
    EnableAcceptCancelDialogAccept (dlg->accept_cancel);
  }

  vnp1 = ValNodeFree (vnp1);
  vnp2 = ValNodeFree (vnp2);
}


static void ChangeTrimFeatureType (Pointer data)
{
  ResolveFeatureOverlapsFormPtr dlg;
  ValNodePtr                    vnp;

  dlg = (ResolveFeatureOverlapsFormPtr) data;
  if (dlg == NULL) return;

  vnp = DialogToPointer (dlg->trim_type);
  if (vnp != NULL) {
    ChangeComplexConstraintFieldType (dlg->trim_constraint, FieldType_feature_field, NULL, vnp->choice);
  }
  vnp = ValNodeFree (vnp);
  SetResolveFeatureOverlapsFormAccept (data);
}


static void ChangeIntersectFeatureType (Pointer data)
{
  ResolveFeatureOverlapsFormPtr dlg;
  ValNodePtr                    vnp;

  dlg = (ResolveFeatureOverlapsFormPtr) data;
  if (dlg == NULL) return;

  vnp = DialogToPointer (dlg->intersect_type);
  if (vnp != NULL) {
    ChangeComplexConstraintFieldType (dlg->intersect_constraint, FieldType_feature_field, NULL, vnp->choice);
  }
  vnp = ValNodeFree (vnp);
  SetResolveFeatureOverlapsFormAccept (data);
}


static Boolean ResolveFeatureOverlapAction (Pointer data)
{
  ResolveFeatureOverlapsFormPtr dlg;
  ValNodePtr                    vnp1, vnp2;
  ResolveFeatureOverlapsData    r;
  SeqEntryPtr                   sep;

  dlg = (ResolveFeatureOverlapsFormPtr) data;
  if (dlg == NULL) return FALSE;

  vnp1 = DialogToPointer (dlg->trim_type);
  vnp2 = DialogToPointer (dlg->intersect_type);
 
  WatchCursor();
  Update();
  r.trim_type = vnp1->choice;
  r.trim_constraint = DialogToPointer (dlg->trim_constraint);
  r.intersect_type = vnp2->choice;
  r.intersect_constraint = DialogToPointer (dlg->intersect_constraint);
  r.lip = OpenLog ("Trimmed Features");

  sep = GetTopSeqEntryForEntityID (dlg->input_entityID);
  VisitBioseqsInSep (sep, &r, ResolveFeatureOverlapsBioseqCallback);
  
  CloseLog (r.lip);
  FreeLog (r.lip);

  r.trim_constraint = ConstraintChoiceSetFree (r.trim_constraint);
  r.intersect_constraint = ConstraintChoiceSetFree (r.intersect_constraint);

  DeleteMarkedObjects (dlg->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (dlg->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dlg->input_entityID, 0, 0);

  ArrowCursor ();
  Update ();
  return TRUE;
}


extern void ResolveFeatureOverlaps (IteM i)
{
  BaseFormPtr        bfp;
  GrouP              h, g;
  SeqEntryPtr        sep;
  WindoW             w;
  ResolveFeatureOverlapsFormPtr dlg;

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

  dlg = (ResolveFeatureOverlapsFormPtr) MemNew (sizeof (ResolveFeatureOverlapsFormData));
  w = FixedWindow (-50, -33, -10, -10, "Resolve Feature Overlaps",
		   StdCloseWindowProc);
  SetObjectExtra (w, dlg, StdCleanupFormProc);
  dlg->form = (ForM) w;

  dlg->input_entityID = bfp->input_entityID;
  dlg->input_itemID   = bfp->input_itemID;
  dlg->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 3, 0, NULL);
  StaticPrompt (g, "Trim features of type", 0, popupMenuHeight, programFont, 'c');
  dlg->trim_type = FeatureSelectionDialogEx (g, FALSE, sep, ChangeTrimFeatureType, dlg);
  dlg->trim_constraint = ComplexConstraintDialog (g, NULL, NULL);
  ChangeComplexConstraintFieldType (dlg->trim_constraint, FieldType_feature_field, NULL, Macro_feature_type_any);

  StaticPrompt (g, "Where they overlap features of type", 0, popupMenuHeight, programFont, 'c');
  dlg->intersect_type = FeatureSelectionDialogEx (g, FALSE, sep, ChangeIntersectFeatureType, dlg);
  dlg->intersect_constraint = ComplexConstraintDialog (g, NULL, NULL);
  ChangeComplexConstraintFieldType (dlg->intersect_constraint, FieldType_feature_field, NULL, Macro_feature_type_any);

  /* Accept and Cancel buttons */

  dlg->accept_cancel = AcceptCancelDialog (h, ResolveFeatureOverlapAction, NULL, 
                                          NULL, NULL, (Pointer) dlg, w);

  /* Line things up and display the window */

  AlignObjects (ALIGN_CENTER,
                (HANDLE) g,
		            (HANDLE) dlg->accept_cancel,
		            NULL);

  RealizeWindow (w);
  Show (w);
  Select (w);
  Update ();
}


static Boolean IsExternalGeneralID (SeqIdPtr sip)
{
  DbtagPtr dbtag;

  if (sip == NULL || sip->data.ptrvalue == NULL || sip->choice != SEQID_GENERAL) return FALSE;
  dbtag = (DbtagPtr) sip->data.ptrvalue;
  if (StringCmp (dbtag->db, "NCBIFILE") != 0
      && StringCmp (dbtag->db, "TMSMART") != 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void ListBioseqsWithExternalGeneralID (BioseqPtr bsp, Pointer userdata)
{
  SeqIdPtr sip;

  if (bsp != NULL && userdata != NULL)
  {
    for (sip = bsp->id; sip != NULL && !IsExternalGeneralID (sip); sip = sip->next)
    {
    }
    if (sip != NULL) 
    {
      ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_BIOSEQ, bsp);
    }
  }      
}


static SeqIdPtr LocalIdFromGeneralId (SeqIdPtr sip)
{
  DbtagPtr dbtag;
  ObjectIdPtr oip;
  SeqIdPtr sip_new;
  Char     label[15];

  if (!IsExternalGeneralID (sip)) return NULL;

  dbtag = (DbtagPtr) sip->data.ptrvalue;
  if (dbtag->tag == NULL) return NULL;
  
  sip_new = ValNodeNew (NULL);
  sip_new->choice = SEQID_LOCAL;
  oip = ObjectIdNew ();
  if (dbtag->tag->id > 0)
  {
    sprintf (label, "%d", dbtag->tag->id);
    oip->str = StringSave (label);
  }
  else
  {
    oip->str = StringSave (dbtag->tag->str);
  }
  sip_new->data.ptrvalue = oip;
  return sip_new;
}

extern void ConvertGeneralIdToLocalID (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  ValNodePtr         bsp_list = NULL, vnp;
  BioseqPtr          bsp;
  SeqIdPtr           sip_new, sip_next, sip_prev, sip_list, sip;

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

  VisitBioseqsInSep (sep, &bsp_list, ListBioseqsWithExternalGeneralID);

  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp == NULL) continue;
    sip = bsp->id;
    sip_prev = NULL;
    while (sip != NULL)
    {
      sip_next = sip->next;
      if (IsExternalGeneralID (sip) && (sip_new = LocalIdFromGeneralId (bsp->id)) != NULL)
      {
        sip_new->next = sip->next;
        sip->next = NULL;
        if (sip_prev == NULL)
        {
          bsp->id = sip_new;
        }
        else
        {
          sip_prev->next = sip_new;
        }
        sip_list = bsp->id;
        /* now use just old ID for replace */
        bsp->id = sip;
        BioseqReplaceID (bsp, sip_new);
        /* put list back */
        bsp->id = SeqIdFree (bsp->id);
        bsp->id = sip_list;
        SeqMgrReplaceInBioseqIndex (bsp);
      }
      else
      {
        sip_prev = sip;
      }
      sip = sip_next;
    }
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


extern void AbbreviateCitSubAffilStates (IteM i)
{
  BaseFormPtr        bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL)
    return;

  FixUsaAndStateAbbreviations(bfp->input_entityID, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}


static void CreateRefSeqProteinIDForCodingRegion (SeqFeatPtr sfp, SeqIdPtr master) 
{
  BioseqPtr prot_bsp;
  SeqIdPtr  sip, sip_new;
  ObjectIdPtr oip;
  Int4        start, stop;
  CharPtr     new_id_txt;
  CharPtr     id_fmt = "%s:%d-%d";
  Char        id_buf[255];

  if (sfp == NULL || master == NULL) return;

  start = SeqLocStart (sfp->location);
  stop = SeqLocStop (sfp->location);

  SeqIdWrite (master, id_buf, PRINTID_TEXTID_ACC_ONLY, sizeof (id_buf) - 1);
  new_id_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (id_fmt) + StringLen (id_buf) + 30));
  sprintf (new_id_txt, id_fmt, id_buf, start + 1, stop + 1);

  prot_bsp = BioseqFindFromSeqLoc (sfp->product);
  if (prot_bsp == NULL) return;

  /* find local ID */
  sip = prot_bsp->id;
  while (sip != NULL && sip->choice != SEQID_LOCAL) {
    sip = sip->next;
  }
  if (sip == NULL || sip->data.ptrvalue == NULL) return;

  oip = (ObjectIdPtr) sip->data.ptrvalue;
  if (StringCmp (oip->str, new_id_txt) != 0) {
    sip_new = ValNodeNew (NULL);
    sip_new->choice = SEQID_LOCAL;
    oip = ObjectIdNew ();
    oip->str = new_id_txt;
    new_id_txt = NULL;
    sip_new->data.ptrvalue = oip;
    BioseqReplaceID (prot_bsp, sip_new);
    sfp->product = SeqLocFree (sfp->product);
    sfp->product = ValNodeNew (NULL);
    sfp->product->choice = SEQLOC_WHOLE;
    sfp->product->data.ptrvalue = (Pointer) sip_new;
  } 
  new_id_txt = MemFree (new_id_txt);
}


static void CreateRefSeqProteinIDsForBioseq (BioseqPtr bsp, Pointer data)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  SeqIdPtr sip;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }

  sip = SeqIdFindBest (bsp->id, SEQID_OTHER);
  if (sip == NULL || sip->choice != SEQID_OTHER) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &context);
  while (sfp != NULL) {
    CreateRefSeqProteinIDForCodingRegion (sfp, sip);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &context);
  }
}


NLM_EXTERN void CreateRefSeqProteinIDs (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;

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

  VisitBioseqsInSep (sep, NULL, CreateRefSeqProteinIDsForBioseq);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


typedef struct flucommentlist {
  BioseqPtr bsp;
  BioSourcePtr biop;
} FluCommentListData, PNTR FluCommentListPtr;


static FluCommentListPtr FluCommentListNew (BioseqPtr bsp, BioSourcePtr biop)
{
  FluCommentListPtr f;

  f = (FluCommentListPtr) MemNew (sizeof (FluCommentListData));
  f->bsp = bsp;
  f->biop = biop;
  return f;
}


static int LIBCALLBACK SortFluCommentList (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  FluCommentListPtr f1, f2;
  Int4              rval = 0;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      f1 = (FluCommentListPtr) vnp1->data.ptrvalue;
      f2 = (FluCommentListPtr) vnp2->data.ptrvalue;
      if (f1 != NULL && f2 != NULL) {
        if (f1->biop == NULL && f2->biop == NULL) {
          rval = 0;
        } else if (f1->biop == NULL) {
          rval = -1;
        } else if (f2->biop == NULL) {
          rval = 1;
        } else if (f1->biop->org == NULL && f2->biop->org == NULL) {
          rval = 0;
        } else if (f1->biop->org == NULL) {
          rval = -1;
        } else if (f2->biop->org == NULL) {
          rval = 1;
        } else {
          rval = StringCmp (f1->biop->org->taxname, f2->biop->org->taxname);
        }
      }
    }
  }
  return rval;
}

static void GetBioSourceAndSeqIdPairs (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, 0, FluCommentListNew (bsp, sdp->data.ptrvalue));
  }
}


static CharPtr MakeFluCommentIdList (ValNodePtr vnp_start, ValNodePtr vnp_end)
{
  FluCommentListPtr  f;
  ValNodePtr         id_list = NULL;
  ValNodePtr         vnp, vnp_a;
  Int4               prefix_len, count, last_number, new_number, id_len;
  Char               id_buf[255];
  CharPtr            id_start, id_next, last_id;
  SeqIdPtr           sip;
  Boolean            mismatch = FALSE;
  Int4               list_len = 0;
  CharPtr            list = NULL;

  if (vnp_start == NULL || vnp_start->next == vnp_end) {
    return NULL;
  }

  /* make temporary list of id strings */
  for (vnp_a = vnp_start; vnp_a != vnp_end; vnp_a = vnp_a->next) {
    f = (FluCommentListPtr) vnp_a->data.ptrvalue;
    sip = SeqIdFindBest (f->bsp->id, SEQID_GENBANK);
    SeqIdWrite (sip, id_buf, PRINTID_TEXTID_ACC_ONLY, sizeof (id_buf) - 1);
    ValNodeAddPointer (&id_list, 0, StringSave (id_buf));
  }

  /* sort id_list */
  id_list = ValNodeSort (id_list, SortVnpByString);

  /* calculate length of id list */
  vnp_a = id_list;
  while (vnp_a != NULL) {
    id_start = vnp_a->data.ptrvalue;
    list_len += StringLen (vnp_a->data.ptrvalue) + 2;
    prefix_len = StringCSpn (id_start, "1234567890");
    if (prefix_len == 0) {
      vnp_a = vnp_a->next;
    } else if (StringSpn (id_start + prefix_len, "1234567890") != StringLen (id_start + prefix_len)) {
      vnp_a = vnp_a->next;
    } else {
      last_number = atoi (id_start + prefix_len);
      count = 1;
      vnp = vnp_a->next;
      mismatch = FALSE;
      while (vnp != NULL && !mismatch) {
        id_next = vnp->data.ptrvalue;
        if (StringNICmp (id_next, id_start, prefix_len) != 0) {
          mismatch = TRUE;
        } else if (StringSpn (id_next + prefix_len, "1234567890") != StringLen (id_next + prefix_len)) {
          mismatch = TRUE;
        } else if ((new_number = atoi (id_next + prefix_len)) != last_number + 1) {
          mismatch = TRUE;
        } else {
          count++;
          last_number = new_number;
          id_len = StringLen (id_next);
          vnp = vnp->next;
        }
      }
      if (count > 1) {
        list_len += id_len + 2;
      }
      vnp_a = vnp;
    }
  }

  /* allocate memory */
  list = (CharPtr) MemNew (sizeof (Char) * list_len);
  list[0] = 0;

  /* make list */
  vnp_a = id_list;
  while (vnp_a != NULL) {
    if (list[0] != 0) {
      StringCat (list, ", ");
    }
    id_start = vnp_a->data.ptrvalue;
    StringCat (list, id_start);
    prefix_len = StringCSpn (id_start, "1234567890");
    if (prefix_len == 0) {
      vnp_a = vnp_a->next;
    } else if (StringSpn (id_start + prefix_len, "1234567890") != StringLen (id_start + prefix_len)) {
      vnp_a = vnp_a->next;
    } else {
      last_number = atoi (id_start + prefix_len);
      count = 1;
      vnp = vnp_a->next;
      mismatch = FALSE;
      while (vnp != NULL && !mismatch) {
        id_next = vnp->data.ptrvalue;
        if (StringNICmp (id_next, id_start, prefix_len) != 0) {
          mismatch = TRUE;
        } else if (StringSpn (id_next + prefix_len, "1234567890") != StringLen (id_next + prefix_len)) {
          mismatch = TRUE;
        } else if ((new_number = atoi (id_next + prefix_len)) != last_number + 1) {
          mismatch = TRUE;
        } else {
          count++;
          last_number = new_number;
          last_id = id_next;
          id_len = StringLen (id_next);
          vnp = vnp->next;
        }
      }
      if (count > 1) {
        if (count == 2) {
          StringCat (list, ", ");
          StringCat (list, last_id);
        } else {
          StringCat (list, "-");
          StringCat (list, last_id);
        }
      }
      vnp_a = vnp;
    }
  }

  id_list = ValNodeFreeData (id_list);

  return list;
}


static void AddCommentsToList (ValNodePtr vnp_start, ValNodePtr vnp_end, CharPtr last_taxname)
{
  CharPtr            id_list_txt;
  CharPtr            comment_fmt = "GenBank Accession Numbers %s represent sequences from the %d segments of %s";
  CharPtr            comment;
  Int4               comment_len = 0, num_segments = 0;
  SeqDescrPtr        sdp;
  FluCommentListPtr  f;
  ValNodePtr         vnp_a;

  if (vnp_start == NULL || StringHasNoText (last_taxname)) {
    return;
  } else if (vnp_start->next == vnp_end) {
    /* only one in the group, don't bother with comment */
  } else {
    /* get number of segments */
    for (vnp_a = vnp_start; vnp_a != vnp_end; vnp_a = vnp_a->next) {
      num_segments++;
    }

    id_list_txt = MakeFluCommentIdList (vnp_start, vnp_end);

    /* make comment */
    comment = (CharPtr) MemNew (sizeof (Char) * (StringLen (comment_fmt) + 15 + StringLen (id_list_txt) + StringLen (last_taxname) ));
    sprintf (comment, comment_fmt, id_list_txt, num_segments, last_taxname);
    id_list_txt = MemFree (id_list_txt);
    /* add comments */
    for (vnp_a = vnp_start; vnp_a != vnp_end; vnp_a = vnp_a->next) {
      f = (FluCommentListPtr) vnp_a->data.ptrvalue;
      if (f != NULL && f->bsp != NULL) {
        sdp = CreateNewDescriptorOnBioseq (f->bsp, Seq_descr_comment);
        sdp->data.ptrvalue = StringSave (comment);
      }
    }
    comment = MemFree (comment);
  }
}


NLM_EXTERN void AddFluComments (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  ValNodePtr         list = NULL, vnp, vnp_start;
  CharPtr            last_taxname = NULL;
  FluCommentListPtr  f;

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

  /* get list of biosources and sequence IDs */
  VisitBioseqsInSep (sep, &list, GetBioSourceAndSeqIdPairs);
  /* sort so that identical organisms are together */
  list = ValNodeSort (list, SortFluCommentList);
  /* create comments for all sequences with the same organism name */
  vnp_start = NULL;
  vnp = list;
  while (vnp != NULL) {
    f = (FluCommentListPtr) vnp->data.ptrvalue;
    
    if (f == NULL || f->biop == NULL || f->biop->org == NULL) {
      /* skip this one */
    } else if (last_taxname == NULL) {
      last_taxname = f->biop->org->taxname;
      vnp_start = vnp;
    } else if (StringCmp (last_taxname, f->biop->org->taxname) != 0) {      
      if (vnp_start->next == vnp) {
        /* only one in the group, don't bother with comment */
      } else {
        AddCommentsToList (vnp_start, vnp, last_taxname);
      }
      vnp_start = vnp;
      f = (FluCommentListPtr) vnp->data.ptrvalue;       
      last_taxname = f->biop->org->taxname;
    }
    vnp = vnp->next;
  }
  AddCommentsToList (vnp_start, NULL, last_taxname);

  list = ValNodeFreeData (list);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


static CharPtr RNA_words[] = {
  "ITS1",
  "ITS2",
  "5.8S",
  "16S",
  "18S",
  "28S"
};

static const Int4 k_num_RNA_words = sizeof (RNA_words) / sizeof (CharPtr);


static ValNodePtr RNAWordsFromString (CharPtr str, TextFsaPtr tags)
{
  Int4        state;
  Char        ch;
  CharPtr     ptr;
  ValNodePtr  matches;
  CharPtr     last_hit = NULL, last_pos = NULL;
  Int4        match_len;
  ValNodePtr  tokens = NULL;

  if (StringHasNoText (str) || tags == NULL) {
    return NULL;
  }

  state = 0;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    matches = NULL;
    state = TextFsaNext (tags, state, ch, &matches);
    if (matches != NULL && (isspace (*(ptr + 1)) || ispunct (*(ptr + 1)) || *(ptr + 1) == 0) && (match_len = StringLen (matches->data.ptrvalue)) > 0
        && (ptr - match_len + 1 == str || isspace (*(ptr - match_len)) || ispunct (*(ptr - match_len)))) {
      ValNodeAddPointer (&tokens, 0, ptr - match_len + 1);
      last_pos = ptr;
      last_hit = (CharPtr) matches->data.ptrvalue;
    }
    ptr++;
    ch = *ptr;
  }
  return tokens;
}


static void ExpandFeatureForIntervals (SeqFeatPtr sfp)
{
  SeqLocPtr   slp, slp_next, orig;
  SeqFeatPtr  sfp_new;
  Boolean     part5, part3;

  orig = sfp->location;
  sfp->location = orig->data.ptrvalue;
  orig->data.ptrvalue = NULL;
  orig = SeqLocFree (orig);
  slp = sfp->location->next;
  sfp->location->next = NULL;
  sfp->partial = CheckSeqLocForPartial (sfp->location, &part5, &part3);
  while (slp != NULL) {
    sfp_new = SeqFeatCopy (sfp);
    sfp_new->next = sfp->next;
    sfp->next = sfp_new;
    sfp_new->location = SeqLocFree (sfp_new->location);
    slp_next = slp->next;
    slp->next = NULL;
    sfp_new->location = slp;
    sfp_new->partial = CheckSeqLocForPartial (sfp_new->location, &part5, &part3);
    slp = slp_next;
    sfp = sfp_new;
  }
} 
    

/* TODO - for each string, remove trailing delimiters and spaces, replace abbreviations */

static void FixExplodedRNAProduct (CharPtr PNTR pProduct)
{
  CharPtr cp;
  Int4    len, i;
  Char    repl_buf[200];
  CharPtr repl_fmt = "%s ribosomal RNA";

  if (pProduct == NULL || *pProduct == NULL) {
    return;
  }

  FindReplaceString (pProduct, " and ", "", FALSE, FALSE);

  len = StringLen (*pProduct);
  cp = *pProduct + len - 1;

  while (ispunct (*cp) || isspace (*cp)) {
    *cp = 0;
    cp--;
  }
   
  if (StringStr (*pProduct, " ribosomal RNA") != NULL) {
    return;
  }
  for (i = 2; i < k_num_RNA_words; i++) {
    sprintf (repl_buf, repl_fmt, RNA_words[i]);
    FindReplaceString (pProduct, RNA_words[i], repl_buf, FALSE, TRUE);
  }
}


static void ExplodeRNACallback (SeqFeatPtr sfp, Pointer data)
{
  TextFsaPtr  tags;
  Int4        num_interval = 0, len;
  SeqLocPtr   tmp;
  ValNodePtr  tokens = NULL, vnp;
  RnaRefPtr   rrp, rrp_new;
  RNAGenPtr   rgp, rgp_new;
  SeqFeatPtr  sfp_new;
  CharPtr     new_product;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || (tags = (TextFsaPtr) data) == NULL
      || sfp->location == NULL
      || (sfp->location->choice != SEQLOC_MIX && sfp->location->choice != SEQLOC_PACKED_INT)) {
    return;
  }

  for (tmp = sfp->location->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
    num_interval++;
  }

  if (!StringHasNoText (sfp->comment)
      && (tokens = RNAWordsFromString (sfp->comment, tags)) != NULL
      && (num_interval == ValNodeLen (tokens))) {
    ExpandFeatureForIntervals (sfp);
    for (vnp = tokens->next, sfp_new = sfp->next; vnp != NULL && sfp_new != NULL; vnp = vnp->next, sfp_new = sfp_new->next) {
      sfp_new->comment = MemFree (sfp_new->comment);
      sfp_new->comment = StringSave (vnp->data.ptrvalue);
      if (vnp->next != NULL) {
        len = (CharPtr) vnp->next->data.ptrvalue - (CharPtr)vnp->data.ptrvalue;
        sfp_new->comment[len] = 0;
      }
      FixExplodedRNAProduct(&(sfp_new->comment));
    }
    /* now truncate original */
    *((CharPtr) tokens->next->data.ptrvalue) = 0;
    FixExplodedRNAProduct(&(sfp->comment));
    tokens = ValNodeFree (tokens);
  } else if ((rrp = (RnaRefPtr) sfp->data.value.ptrvalue) != NULL) {
    if (rrp->ext.choice == 1) {
      tokens = RNAWordsFromString (rrp->ext.value.ptrvalue, tags);
      if (num_interval == ValNodeLen (tokens)) {
        ExpandFeatureForIntervals (sfp);
        for (vnp = tokens->next, sfp_new = sfp->next; vnp != NULL && sfp_new != NULL; vnp = vnp->next, sfp_new = sfp_new->next) {
          rrp_new = sfp_new->data.value.ptrvalue;
          rrp_new->ext.value.ptrvalue = MemFree (rrp_new->ext.value.ptrvalue);
          new_product = StringSave (vnp->data.ptrvalue);
          if (vnp->next != NULL) {
            len = (CharPtr) vnp->next->data.ptrvalue - (CharPtr)vnp->data.ptrvalue;
            new_product[len] = 0;
          }
          FixExplodedRNAProduct (&new_product);
          rrp_new->ext.value.ptrvalue = new_product;
        }  
        *((CharPtr) tokens->next->data.ptrvalue) = 0;
        new_product = rrp->ext.value.ptrvalue;
        FixExplodedRNAProduct (&new_product);
        rrp->ext.value.ptrvalue = new_product;
      }
      tokens = ValNodeFree (tokens);
    } else if (rrp->ext.choice == 3 && (rgp = (RNAGenPtr) rrp->ext.value.ptrvalue) != NULL) {
      tokens = RNAWordsFromString (rgp->product, tags);
      if (num_interval == ValNodeLen (tokens)) {
        ExpandFeatureForIntervals (sfp);
        for (vnp = tokens->next, sfp_new = sfp->next; vnp != NULL && sfp_new != NULL; vnp = vnp->next, sfp_new = sfp_new->next) {
          rrp_new = sfp_new->data.value.ptrvalue;
          rgp_new = rrp_new->ext.value.ptrvalue;
          rgp_new->product = MemFree (rgp_new->product);
          rgp_new->product = StringSave (vnp->data.ptrvalue);
          if (vnp->next != NULL) {
            len = (CharPtr) vnp->next->data.ptrvalue - (CharPtr)vnp->data.ptrvalue;
            rgp_new->product[len] = 0;
          }
          FixExplodedRNAProduct(&(rgp_new->product));
        }      
        *((CharPtr) tokens->next->data.ptrvalue) = 0;
        FixExplodedRNAProduct (&(rgp->product));
      }
      tokens = ValNodeFree (tokens);
    }
  }
}


NLM_EXTERN void ExplodeRNA (IteM i)
{
  BaseFormPtr        bfp;
  SeqEntryPtr        sep;
  TextFsaPtr         tags;
  Int4               j;

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

  
  tags = TextFsaNew();

  for (j = 0; j < k_num_RNA_words; j++) {
    TextFsaAdd (tags, RNA_words[j]);
  }

  VisitFeaturesInSep (sep, tags, ExplodeRNACallback);

  tags = TextFsaFree (tags);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


static void RemoveFeaturesLikeGapsCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr           sfp, gap_feat;
  SeqMgrFeatContext    context;
  DeltaSeqPtr          dsp;
  SeqLocPtr            slp;
  Int4                 pos = 0, len;
  CharPtr              new_note;
  ImpFeatPtr           imp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) {
    return;
  }

  dsp = (DeltaSeqPtr)(bsp->seq_ext);
  sfp = SeqMgrGetNextFeature(bsp, NULL, 0, 0, &context);
  while (dsp != NULL && sfp != NULL) {
    while (dsp != NULL && !IsDeltaSeqGap(dsp)) {
      if (dsp->choice == 1) {
        pos += SeqLocLen (dsp->data.ptrvalue);
      } else if (dsp->choice == 2) {
        pos += ((SeqLitPtr)dsp->data.ptrvalue)->length;
      }
      dsp = dsp->next;
    }
    if (dsp != NULL) {
      len = ((SeqLitPtr)dsp->data.ptrvalue)->length;
      while (context.left < pos && sfp != NULL) {
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
      }
      if (sfp != NULL) {
        gap_feat = NULL;
        new_note = NULL;
        while (context.left == pos && sfp != NULL) {
          if (context.right == pos + len - 1) {
            if (sfp->idx.subtype == FEATDEF_gap) {
              gap_feat = sfp;
            } else {
              if (sfp->comment != NULL) {
                SetStringValue(&new_note, sfp->comment, ExistingTextOption_append_semi);
              }
              sfp->idx.deleteme = TRUE;
            }
          }
          sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
        }
        /* only need to instantiate gap if there is a note to add */
        if (new_note != NULL) {
          if (gap_feat == NULL) {
            slp = SeqLocIntNew (pos, pos + len - 1, Seq_strand_plus, SeqIdFindWorst (bsp->id));
            gap_feat = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, slp);
            imp = ImpFeatNew ();
            imp->key = StringSave ("gap");
            gap_feat->data.value.ptrvalue = imp;
          }
          if (gap_feat->comment == NULL) {
            gap_feat->comment = new_note;
          } else {
            SetStringValue (&(gap_feat->comment), new_note, ExistingTextOption_append_semi);
            new_note = MemFree (new_note);
          }
        }
      }
      pos += len;
      dsp = dsp->next;
    }
  }
}


NLM_EXTERN void RemoveFeaturesLikeGaps (IteM i)
{
  BaseFormPtr          bfp;
  SeqEntryPtr          sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitBioseqsInSep (sep, NULL, RemoveFeaturesLikeGapsCallback);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


NLM_EXTERN void SplitPCRPrimersByPositionMenuItem (IteM i)
{
  BaseFormPtr          bfp;
  SeqEntryPtr          sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  SplitPCRPrimersByPosition (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


typedef struct splitpcrprimersbyconstraintform {
  FORM_MESSAGE_BLOCK

  DialoG fwd_constraint;
  DialoG rev_constraint;
} SplitPCRPrimersByConstraintFormData, PNTR SplitPCRPrimersByConstraintFormPtr;


static void SplitPCRPrimersByConstraintAccept (ButtoN b)
{
  SplitPCRPrimersByConstraintFormPtr frm;
  StringConstraintPtr fwd_scp, rev_scp;
  SeqEntryPtr sep;

  frm = (SplitPCRPrimersByConstraintFormPtr) GetObjectExtra (b);
  if (frm == NULL) {
    return;
  }

  fwd_scp = DialogToPointer (frm->fwd_constraint);
  rev_scp = DialogToPointer (frm->rev_constraint);
  sep = GetTopSeqEntryForEntityID (frm->input_entityID);

  SplitPCRPrimersByConstraints (sep, fwd_scp, rev_scp);
  fwd_scp = StringConstraintFree (fwd_scp);
  rev_scp = StringConstraintFree (rev_scp);
  ObjMgrSetDirtyFlag (frm->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, frm->input_entityID, 0, 0);  
  Remove (frm->form);
}

  
NLM_EXTERN void SplitPCRPrimersByConstraintsMenuItem (IteM i)
{
  BaseFormPtr          bfp;
  SplitPCRPrimersByConstraintFormPtr frm;
  WindoW              w;
  GrouP               h, g, c;
  ButtoN              b;
  PrompT              ppt;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    return;
  }

  frm = (SplitPCRPrimersByConstraintFormPtr) MemNew (sizeof (SplitPCRPrimersByConstraintFormData));

  w = FixedWindow (-50, -33, -10, -10, "Convert BioSource Dbxrefs to Feature Dbxrefs", StdCloseWindowProc);
  SetObjectExtra (w, frm, StdCleanupExtraProc);
  frm->form = (ForM) w;
  frm->input_entityID = bfp->input_entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt = StaticPrompt (h, "Choose string constraints for forward and reverse primer names for new sets", 0, dialogTextHeight, programFont, 'l');
  g = HiddenGroup (h, 0, 2, NULL);
  SetGroupSpacing (g, 10, 10);
  frm->fwd_constraint = StringConstraintDialog (g, "Forward", FALSE, NULL, NULL);
  frm->rev_constraint = StringConstraintDialog (g, "Reverse", FALSE, NULL, NULL);

  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton (c, "Accept", SplitPCRPrimersByConstraintAccept);
  SetObjectExtra (b, frm, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) g, (HANDLE) c, NULL);
  Show (w);
  Select (w);
}


NLM_EXTERN void MergePCRPrimerSetsMenuItem (IteM i)
{
  BaseFormPtr          bfp;
  SeqEntryPtr          sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  MergePCRPrimers (sep);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


static void ConvertToDelayedGenProdSetQualifiersCallback (SeqFeatPtr sfp, Pointer data)
{
  GBQualPtr gbq;

  if (sfp == NULL) {
    return;
  } else if (sfp->idx.subtype != FEATDEF_CDS && sfp->idx.subtype != FEATDEF_mRNA) {
    return;
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "protein_id") == 0) {
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("orig_protein_id");
    } else if (StringICmp (gbq->qual, "transcript_id") == 0) {
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("orig_transcript_id");
    }
  }
}


NLM_EXTERN void ConvertToDelayedGenProdSetQualifiers (IteM i)
{
  BaseFormPtr          bfp;
  SeqEntryPtr          sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) {
    return;
  }

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  VisitFeaturesInSep (sep, NULL, ConvertToDelayedGenProdSetQualifiersCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);  
}


#define NUM_FEAT_ED_PAGES 2

#define NUM_FEAT_ED_ACTIONS 8

typedef struct feated
{
  FORM_MESSAGE_BLOCK
  DialoG  feature_type;
  DialoG  constraint;
  DialoG  accept_cancel;

  /* for evidence */  
  Int4    evidence_val;
      
  /* for exceptions */
  PopuP   explanation_choice;
  ButtoN  move_to_comment;
  GrouP   product_choice;
  TexT    other_explanation;
  CharPtr explanation_text;
  Boolean do_move_to_comment;
  Int4    product_val;
  
  /* for partials */
  PopuP  convert_choice;
  PopuP  partial5_choice;
  ButtoN extend5_btn;
  PopuP  partial3_choice;
  ButtoN extend3_btn;
  ButtoN retranslate_btn;
  ButtoN adjust_gene_btn;
  Int4   orderjoinpolicy;
  Int4   leftpolicy;
  Int4   rightpolicy;
  Boolean extend5;
  Boolean extend3;
  Boolean retranslate;
  Boolean adjust_gene;
  
  /* for strand */
  PopuP  from_strand;
  PopuP  to_strand;
  Uint1  from_strand_val;
  Uint1  to_strand_val;

  GrouP  action_grps [NUM_FEAT_ED_ACTIONS];
  ButtoN action_btns [NUM_FEAT_ED_ACTIONS];
  GrouP  action_choice_grp;
  
  /* for citations */
  DialoG     citation_list;
  ButtoN     explanation_constraint;
  PopuP      explanation_constraint_choice;
  ButtoN     citation_translation;
  Boolean    use_explanation_constraint;
  Int4       explanation_constraint_choice_val;
  Boolean    use_product_match;

  /* for experiment */
  GrouP      experiment_action;
  TexT       experiment_text;
  
  /* for pseudo */
  GrouP      pseudo_action;
  Boolean    pseudo_val;
  
  /* for inferences */
  DialoG     inference_dlg;
} FeatEdData, PNTR FeatEdPtr;


static void FeatEdChangeNotify (Pointer userdata)
{
  FeatEdPtr mp;
  ValNodePtr      err_list = NULL;
  CharPtr         str;

  mp = (FeatEdPtr) userdata;
  if (mp == NULL) return;
  
  if (mp->explanation_choice != NULL)
  {
    if (GetValue (mp->explanation_choice) == 9 
        && ((mp->action_btns[1] == NULL && mp->action_choice_grp == NULL)
            || (mp->action_btns[1] != NULL && GetStatus (mp->action_btns[1]))
            || (mp->action_choice_grp != NULL && GetValue (mp->action_choice_grp) == 2)))
    {
      Enable (mp->other_explanation);
      str = SaveStringFromText (mp->other_explanation);
      if (StringHasNoText (str))
      {
        DisableAcceptCancelDialogAccept (mp->accept_cancel);
        str = MemFree (str);
        return;
      }
      str = MemFree (str);
    }
    else
    {
      Disable (mp->other_explanation);
    }
  }
  
  if (GetValue (mp->action_choice_grp) == FEAT_ED_INFERENCE)
  {
    ValNodeLink (&err_list, TestDialog (mp->inference_dlg));
  }
  
  if (err_list == NULL)
  {
    EnableAcceptCancelDialogAccept (mp->accept_cancel);
  }
  else
  {
    DisableAcceptCancelDialogAccept (mp->accept_cancel);
  }
  ValNodeFree (err_list);  
}


static void FeatEdClear (Pointer data)
{
  FeatEdPtr mp;

  mp = (FeatEdPtr) data;
  if (mp == NULL) return;
 
  SetFeatureTypeInFeatureTypeDialog(mp->feature_type, Macro_feature_type_any);
  PointerToDialog (mp->constraint, NULL);

  if (mp->explanation_choice != NULL)
  {
    SetValue (mp->explanation_choice, 1);
  }
  if (mp->move_to_comment != NULL)
  {
    SetStatus (mp->move_to_comment, FALSE);
  }
  if (mp->product_choice != NULL)
  {
    SetValue (mp->product_choice, 3);
  }
  if (mp->other_explanation != NULL)
  {
    SetTitle (mp->other_explanation, "");
    Disable (mp->other_explanation);
  }
  
  if (mp->convert_choice != NULL)
  {
    SetValue (mp->convert_choice, 1);
  }
  if (mp->partial5_choice != NULL)
  {
    SetValue (mp->partial5_choice, 8);
  }
  if (mp->extend5_btn != NULL)
  {
    SetStatus (mp->extend5_btn, FALSE);
  }
  if (mp->partial3_choice != NULL)
  {
    SetValue (mp->partial3_choice, 7);
  }
  if (mp->extend3_btn != NULL)
  {
    SetStatus (mp->extend3_btn, FALSE);
  }
  if (mp->retranslate_btn != NULL)
  {
    SetStatus (mp->retranslate_btn, FALSE);
  }
  if (mp->adjust_gene_btn != NULL) 
  {
    SetStatus (mp->adjust_gene_btn, FALSE);
  }

  /* for strand */
  if (mp->from_strand != NULL)
  {
    SetValue (mp->from_strand, 1);
  }
  if (mp->to_strand != NULL)
  {
    SetValue (mp->to_strand, 1);
  }

  /* for experiment */
  if (mp->experiment_text != NULL)
  {
    SetTitle (mp->experiment_text, "");
  }
  
  /* for pseudo */
  if (mp->pseudo_action != NULL)
  {
    SetValue (mp->pseudo_action, 1);
  }

  FeatEdChangeNotify (mp);
}


static void FeatEdClearText (Pointer data)
{
  FeatEdPtr     mp;

  mp = (FeatEdPtr) data;
  if (mp == NULL) return;
   
  /* for explanation */
  if (mp->other_explanation != NULL)
  {
    SetTitle (mp->other_explanation, "");
  }

  /* for experiment */
  if (mp->experiment_text != NULL)
  {
    SetTitle (mp->experiment_text, "");
  }

  FeatEdChangeNotify (mp);
}

static GrouP FeaturePseudoGroup (GrouP h, FeatEdPtr mp)
{
  GrouP g;
  
  if (mp == NULL)
  {
    return NULL;
  }
  
  g = HiddenGroup (h, 0, 0, NULL);

  mp->pseudo_action = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (mp->pseudo_action, 10, 10);
  RadioButton (mp->pseudo_action, "Set Pseudo");
  RadioButton (mp->pseudo_action, "Clear Pseudo");
  SetValue (mp->pseudo_action, 1);
  
  return g;  
}

static GrouP FeatureEvidenceGroup (GrouP h, FeatEdPtr mp)
{
  GrouP g;
  
  if (mp == NULL)
  {
    return NULL;
  }

  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  
  StaticPrompt (g, "Click Accept to clear feature evidence.", 0, dialogTextHeight, systemFont, 'l');
  
  return g;
}


static void FeatExplanationChange (PopuP p)
{
  FeatEdPtr mp;

  mp = (FeatEdPtr) GetObjectExtra (p);
  if (mp == NULL) return;
  
  if (GetValue (mp->explanation_choice) == 9)
  {
    Enable (mp->other_explanation);
  }
  else
  {
    Disable (mp->other_explanation);
  }
  
  FeatEdChangeNotify (mp);
}

static void FeatOtherExplanationChange (TexT t)
{
  FeatEdPtr mp;
 
  mp = (FeatEdPtr) GetObjectExtra (t);
  FeatEdChangeNotify (mp);
}

static CharPtr exception_explanations[] = 
{
  "RNA Editing", "reasons given in citation", "ribosomal slippage",
  "trans splicing", "artificial frameshift", "nonconsensus splice site",
  "rearrangement required", "other->"
};

static Int4 num_exception_explanations = sizeof (exception_explanations) / sizeof (CharPtr);

static GrouP FeatureExceptionGroup (GrouP h, FeatEdPtr mp)
{
  GrouP  g, k1, k2;
  Int4   i;
  
  if (mp == NULL)
  {
    return NULL;
  }
  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  
  k1 = HiddenGroup (g, 3, 0, NULL);
  SetGroupSpacing (k1, 10, 10);
  StaticPrompt (k1, "Set Explanation to ", 0, dialogTextHeight, systemFont, 'l');
  mp->explanation_choice = PopupList (k1, TRUE, FeatExplanationChange);
  SetObjectExtra (mp->explanation_choice, mp, NULL);
  PopupItem (mp->explanation_choice, " ");
  for (i = 0; i < num_exception_explanations; i++)
  {
    PopupItem (mp->explanation_choice, exception_explanations [i]);
  }
  SetValue (mp->explanation_choice, 1);
  mp->other_explanation = DialogText (k1, "", 10, FeatOtherExplanationChange);
  Disable (mp->other_explanation);
  SetObjectExtra (mp->other_explanation, mp, NULL);

  k2 = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (k2, 10, 10);
  StaticPrompt (k2, "Where feature product is", 0, dialogTextHeight, systemFont, 'l');
  mp->product_choice = HiddenGroup (k2, 3, 0, NULL);
  RadioButton (mp->product_choice, "Present");
  RadioButton (mp->product_choice, "Absent");
  RadioButton (mp->product_choice, "Either");
  SetValue (mp->product_choice, 3);
  
  mp->move_to_comment = CheckBox (g, "Move explanation to comment", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) k1,
                              (HANDLE) k2,
                              (HANDLE) mp->move_to_comment,
                              NULL);
  
  return g;
}

static void ChangeExplanationConstraint (ButtoN b)
{
  FeatEdPtr mp;
  
  mp = (FeatEdPtr) GetObjectExtra (b);
  if (mp == NULL)
  {
    return;
  }
  if (GetStatus (mp->explanation_constraint))
  {
    Enable (mp->explanation_constraint_choice);
  }
  else
  {
    Disable (mp->explanation_constraint_choice);
  }  
}

static GrouP FeatureEditorCitationGroup (GrouP h, FeatEdPtr mp)
{
  GrouP  g, k1, k2;
  Int4   i;
  
  if (mp == NULL)
  {
    return NULL;
  }
  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  
  /* need citation selector */
  mp->citation_list = FeatCitEditDialog (g, mp->input_entityID);
  PointerToDialog (mp->citation_list, NULL);

  k1 = HiddenGroup (g, 3, 0, NULL);
  SetGroupSpacing (k1, 10, 10);
  mp->explanation_constraint = CheckBox (k1, "Where explanation is", ChangeExplanationConstraint);
  SetObjectExtra (mp->explanation_constraint, mp, NULL);
  SetStatus (mp->explanation_constraint, TRUE);
  mp->explanation_constraint_choice = PopupList (k1, TRUE, NULL);
  PopupItem (mp->explanation_constraint_choice, "Any");
  PopupItem (mp->explanation_constraint_choice, " ");
  for (i = 0; i < num_exception_explanations; i++)
  {
    PopupItem (mp->explanation_constraint_choice, exception_explanations [i]);
  }
  SetValue (mp->explanation_constraint_choice, 4);

  k2 = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (k2, 10, 10);
  mp->citation_translation = CheckBox (k2, "Where translation does not match product", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) mp->citation_list,
                              (HANDLE) k1,
                              (HANDLE) k2,
                              (HANDLE) mp->citation_translation,
                              NULL);
  
  return g;
}


static GrouP FeatureEditorPartialGroup (GrouP h, FeatEdPtr mp)
{
  GrouP g, g1, g2,k1, k2, top;
  
  if (mp == NULL)
  {
    return NULL;
  }
  
  top = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (top, 10, 10);

  g = HiddenGroup (top, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  
  StaticPrompt (g, "Partials:", 0, dialogTextHeight, systemFont, 'l');
  g1 = NormalGroup (g, 1, 0, "", programFont, NULL);
  k1 = HiddenGroup (g1, 2, 0, NULL);

  SetGroupSpacing (k1, 10, 10);
  StaticPrompt (k1, "5' partial:", 0, dialogTextHeight, systemFont, 'l');
  mp->partial5_choice = PopupList (k1, TRUE, NULL);
  PopupItem (mp->partial5_choice, "Set");
  PopupItem (mp->partial5_choice, "Set only if at 5' end");
  PopupItem (mp->partial5_choice, "Set if bad start codon");
  PopupItem (mp->partial5_choice, "Set if CDS frame > 1");
  PopupItem (mp->partial5_choice, "Clear");
  PopupItem (mp->partial5_choice, "Clear if not at 5' end");
  PopupItem (mp->partial5_choice, "Clear if good start codon");
  PopupItem (mp->partial5_choice, "Do not change");
  SetValue (mp->partial5_choice, 8);
  
  mp->extend5_btn = CheckBox (g1, "Extend to 5' end of sequence if setting 5' partial", NULL);

  StaticPrompt (g1, "", 0, dialogTextHeight, systemFont, 'l');

  k2 = HiddenGroup (g1, 2, 0, NULL);
  SetGroupSpacing (k2, 10, 10);
  StaticPrompt (k2, "3' partial:", 0, dialogTextHeight, systemFont, 'l');
  mp->partial3_choice = PopupList (k2, TRUE, NULL);
  PopupItem (mp->partial3_choice, "Set");
  PopupItem (mp->partial3_choice, "Set only if at 3' end");
  PopupItem (mp->partial3_choice, "Set if bad stop codon");
  PopupItem (mp->partial3_choice, "Clear");
  PopupItem (mp->partial3_choice, "Clear if not at 3' end");
  PopupItem (mp->partial3_choice, "Clear if good stop codon");
  PopupItem (mp->partial3_choice, "Do not change");
  SetValue (mp->partial3_choice, 7);
  
  mp->extend3_btn = CheckBox (g1, "Extend to 3' end of sequence if setting 3' partial", NULL);

  StaticPrompt (g, "Join/Order/Merge:", 0, dialogTextHeight, systemFont, 'l');
  g2 = NormalGroup (g, 1, 0, "", programFont, NULL);
  mp->convert_choice = PopupList (g2, TRUE, NULL);
  PopupItem (mp->convert_choice, "Do not convert location");
  PopupItem (mp->convert_choice, "Convert location to join");
  PopupItem (mp->convert_choice, "Convert location to order");
  PopupItem (mp->convert_choice, "Convert location to single interval");
  SetValue (mp->convert_choice, 1);

  mp->retranslate_btn = CheckBox (top, "Retranslate coding regions", NULL);
  mp->adjust_gene_btn = CheckBox (top, "Adjust gene location", NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) mp->retranslate_btn, (HANDLE) mp->adjust_gene_btn, NULL);
  
  return top;
}

#define FEATED_STRAND_ANY_REVERSE 1
#define FEATED_STRAND_UNKNOWN     2
#define FEATED_STRAND_PLUS        3
#define FEATED_STRAND_MINUS       4
#define FEATED_STRAND_BOTH        5

static GrouP FeatureEditorStrandGroup (GrouP h, FeatEdPtr mp)
{
  GrouP  g;
  
  if (mp == NULL)
  {
    return NULL;
  }
  g = HiddenGroup (h, 5, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  StaticPrompt (g, "Convert location strand", 0, dialogTextHeight, systemFont, 'l');
  StaticPrompt (g, "From", 0, dialogTextHeight, systemFont, 'l');
  mp->from_strand = PopupList (g, TRUE, NULL);
  PopupItem (mp->from_strand, "Any");
  PopupItem (mp->from_strand, "Unknown");
  PopupItem (mp->from_strand, "Plus");
  PopupItem (mp->from_strand, "Minus");
  PopupItem (mp->from_strand, "Both");
  SetValue (mp->from_strand, FEATED_STRAND_PLUS); 
  StaticPrompt (g, "To", 0, dialogTextHeight, systemFont, 'l');
  mp->to_strand = PopupList (g, TRUE, NULL);
  PopupItem (mp->to_strand, "Reverse");
  PopupItem (mp->to_strand, "Unknown");
  PopupItem (mp->to_strand, "Plus");
  PopupItem (mp->to_strand, "Minus");
  PopupItem (mp->to_strand, "Both");
  SetValue (mp->to_strand, FEATED_STRAND_MINUS); 
  return g;
}


static void ChangeExperimentAction (GrouP g)
{
  FeatEdPtr mp;
  mp = (FeatEdPtr) GetObjectExtra (g);
  if (mp != NULL)
  {
    if (GetValue (mp->experiment_action) == 1)
    {
      Enable (mp->experiment_text);
    }
    else
    {
      Disable (mp->experiment_text);
    }
  }
}

static GrouP FeatureEditorExperimentGroup (GrouP h, FeatEdPtr mp)
{
  if (mp == NULL)
  {
    return NULL;
  }
  mp->experiment_action = HiddenGroup (h, 2, 0, ChangeExperimentAction);
  SetObjectExtra (mp->experiment_action, mp, NULL);
  SetGroupSpacing (mp->experiment_action, 10, 10);
  RadioButton (mp->experiment_action, "Set Experiment");
  mp->experiment_text = DialogText (mp->experiment_action, "", 20, NULL);
  RadioButton (mp->experiment_action, "Remove Experiment");
  SetValue (mp->experiment_action, 1);
  return mp->experiment_action;
}

static GrouP FeatureEditorInferenceGroup (GrouP h, FeatEdPtr mp)
{
  GrouP g;

  g = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (g, 10, 10);

  mp->inference_dlg = CreateInferenceEditDialog (g, FeatEdChangeNotify, mp);
  return g;
}

static void ChangeFeatureEditorActionGroup (GrouP g)
{
  FeatEdPtr mp;
  Int4            i;
  
  mp = (FeatEdPtr) GetObjectExtra (g);
  if (mp == NULL)
  {
    return;
  }
  
  for (i = 0; i < NUM_FEAT_ED_ACTIONS; i++)
  {
    Hide (mp->action_grps[i]);
  }
  
  i = GetValue (mp->action_choice_grp);
  if (i > 0 && i < NUM_FEAT_ED_ACTIONS + 1)
  {
    Show (mp->action_grps [i-1]);
  }

  FeatEdChangeNotify (mp);
  
}

static GrouP FeatureEditorActionGroup (GrouP h, FeatEdPtr mp, Int4 first_action)
{
  GrouP g, k, n = NULL;
  PrompT p1;
  SeqEntryPtr sep;
  
  if (mp == NULL)
  {
    return NULL;
  }
  g = HiddenGroup (h, -1, 0, NULL);

  sep = GetTopSeqEntryForEntityID(mp->input_entityID);
  p1 = StaticPrompt (g, "For", 0, dialogTextHeight, systemFont, 'l');
  mp->feature_type = FeatureTypeDialogMulti (g, FeatEdChangeNotify, mp);
  
  k = HiddenGroup (g, -1, 0, NULL);
  mp->action_choice_grp = HiddenGroup (g, NUM_FEAT_ED_ACTIONS, 0, ChangeFeatureEditorActionGroup);
  SetObjectExtra (mp->action_choice_grp, mp, NULL);
  RadioButton (mp->action_choice_grp, "Evidence");
  RadioButton (mp->action_choice_grp, "Exceptions");
  RadioButton (mp->action_choice_grp, "Location");
  RadioButton (mp->action_choice_grp, "Strands");
  RadioButton (mp->action_choice_grp, "Citations");
  RadioButton (mp->action_choice_grp, "Experiment");
  RadioButton (mp->action_choice_grp, "Inference");
  RadioButton (mp->action_choice_grp, "Pseudo");
  if (first_action >= 1 && first_action <= NUM_FEAT_ED_ACTIONS)
  {
    SetValue (mp->action_choice_grp, first_action);
  }
  else
  {
    SetValue (mp->action_choice_grp, FEAT_ED_PSEUDO);
  }
    
  n = HiddenGroup (g, 0, 0, NULL);
  mp->action_grps [0] = FeatureEvidenceGroup (n, mp);
  mp->action_grps[1] = FeatureExceptionGroup (n, mp);
  mp->action_grps[2] = FeatureEditorPartialGroup (n, mp);
  mp->action_grps[3] = FeatureEditorStrandGroup (n, mp);
  mp->action_grps[4] = FeatureEditorCitationGroup (n, mp);
  mp->action_grps[5] = FeatureEditorExperimentGroup (n, mp);
  mp->action_grps[6] = FeatureEditorInferenceGroup (n, mp);
  mp->action_grps[7] = FeaturePseudoGroup (n, mp);
  AlignObjects (ALIGN_CENTER, (HANDLE) mp->action_grps [0],
                              (HANDLE) mp->action_grps [1],
                              (HANDLE) mp->action_grps [2],
                              (HANDLE) mp->action_grps [3],
                              (HANDLE) mp->action_grps [4], 
                                            (HANDLE) mp->action_grps [5],
                                            (HANDLE) mp->action_grps [6],
                                            (HANDLE) mp->action_grps [7],
                              NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) mp->action_choice_grp, (HANDLE) n, NULL);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) p1,
                              (HANDLE) mp->feature_type,
                              (HANDLE) k, (HANDLE) n, NULL);
      
  return g;
}

static void DoPseudoFeatureProc (SeqFeatPtr sfp, Pointer userdata, FilterSetPtr fsp)
{
  FeatEdPtr mp;

  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  mp = (FeatEdPtr) userdata;
  
  sfp->pseudo = mp->pseudo_val;  
}


static Boolean DoesInferenceMatchCategory (CharPtr inf_cat, CharPtr match_cat)
{
  if (StringHasNoText (inf_cat)) return FALSE;

  if (StringHasNoText (match_cat) 
      || StringICmp (match_cat, "ANY") == 0
      || StringICmp (inf_cat, match_cat) == 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }

}


static void DoInferenceFeatureProc (SeqFeatPtr sfp, Pointer userdata)
{
  InferenceEditPtr iep;
  GBQualPtr        gbqual, gbqual_last = NULL, gbqual_next;
  CharPtr          new_val;
  InferenceParsePtr ipp;
  Boolean           changed;
  ApplyValueData    avd;

  iep = (InferenceEditPtr) userdata;

  if (iep == NULL || sfp == NULL) return;

  /* apply requested changes */
  gbqual = sfp->qual;
  while (gbqual != NULL)
  {
    gbqual_next = gbqual->next;
    if (StringCmp (gbqual->qual, "inference") == 0)
    {
      if (iep->action == eInferenceRemove)
      {
        if (gbqual_last == NULL)
        {
          sfp->qual = gbqual->next;
        }
        else
        {
          gbqual_last->next = gbqual->next;
        }
        gbqual->next = NULL;
        GBQualFree (gbqual);
      }
      else
      {
        /* perform editing here */
        ipp = ParseInferenceText (gbqual->val);
        changed = FALSE;
        if (ipp != NULL)
        {
          switch (iep->action)
          {
            case eInferenceEditCategory:
              if (DoesInferenceMatchCategory(ipp->category, iep->category_from))
              {
                ipp->category = MemFree (ipp->category);
                ipp->category = StringSave (iep->category_to);
                changed = TRUE;
              }
              break;
            case eInferenceEditType:
              if (DoesInferenceMatchCategory(ipp->type, iep->type_from))
              {
                ipp->type = MemFree (ipp->type);
                ipp->type = StringSave (iep->type_to);
                changed = TRUE;
              }
              break;
            case eInferenceApplyTypeFields:
            case eInferenceEditTypeFields:
              if (iep->field_edit != NULL 
                  && iep->field_edit->edit_apply != NULL
                  && DoesInferenceMatchCategory (ipp->type, iep->field_edit->field_type))
              {
                avd.where_to_replace = EditApplyFindLocation_anywhere;
                avd.field_list = NULL;
                avd.etp = NULL;
                AddEditApplyDataToApplyValue (iep->action == eInferenceApplyTypeFields ? 1 : 2,
                                              iep->field_edit->edit_apply, &avd);
                if (iep->field_edit->field_choice == 0)
                {
                  ipp->first_field = HandleApplyValue (ipp->first_field, &avd);
                }
                else
                {
                  ipp->second_field = HandleApplyValue (ipp->second_field, &avd);
                }
                avd.text_to_replace = MemFree (avd.text_to_replace);
                avd.new_text = MemFree (avd.new_text);
                changed = TRUE;
              }
              break;
            default:
              break;
          }
          if (changed)
          {
            new_val = InferenceTextFromStruct (ipp);
            gbqual->val = MemFree (gbqual->val);
            gbqual->val = new_val;
          }
        }
        gbqual_last = gbqual;
      }
    }
    else
    {
      gbqual_last = gbqual;
    }
    gbqual = gbqual_next;
  }
}


static void DoEvidenceFeatureProc (SeqFeatPtr sfp, Pointer userdata)
{
  FeatEdPtr mp;

  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  mp = (FeatEdPtr) userdata;
  
  sfp->exp_ev = (Uint1) mp->evidence_val;
}


static void DoExceptionFeatureProc (SeqFeatPtr sfp, Pointer userdata)
{
  FeatEdPtr  mp;
  ExistingTextData etd;
  CharPtr          old_explanation = NULL;

  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  mp = (FeatEdPtr) userdata;
  
  if ((mp->product_val == 1 && sfp->product != NULL)
      || (mp->product_val == 2 && sfp->product == NULL)
      || mp->product_val == 3)
  {
    /* first, take care of any pre-existing exception explanation */
    old_explanation = sfp->except_text;
    sfp->except_text = NULL;
    if (mp->do_move_to_comment)
    {
      etd.existing_text_choice = eExistingTextChoiceAppendSemi;
      sfp->comment = HandleExistingText (sfp->comment, old_explanation, &etd);
      old_explanation = NULL;
    }
    else
    {
      old_explanation = MemFree (old_explanation);
    }
    
    if (StringHasNoText (mp->explanation_text))
    {
      sfp->excpt = FALSE;
    }
    else
    {
      sfp->excpt = TRUE;
      sfp->except_text = StringSave (mp->explanation_text);
    }
  }
}

#define DO_NOT_CONVERT             1
#define CONVERT_TO_JOIN            2
#define CONVERT_TO_ORDER           3
#define CONVERT_TO_SINGLE_INTERVAL 4

static void DoPartialFeatureProc (SeqFeatPtr sfp, Pointer userdata)
{
  FeatEdPtr mp;
  Boolean      partial5, partial3;
  Boolean      dash_at_end = FALSE;
  Boolean      star_at_end = FALSE;
  ByteStorePtr bs;
  CharPtr      prot;
  CharPtr      ptr;
  Char         first_char = '\0', ch;
  Boolean      hasNulls;
  SeqLocPtr    slp, firstSlp = NULL, lastSlp = NULL;
  Uint1        strand;
  Boolean      atEnd, moved5 = FALSE;
  BioseqPtr    bsp;
  CdRegionPtr  crp = NULL;
  SeqFeatPtr   gene = NULL;
  
  if (sfp == NULL || sfp->location == NULL || userdata == NULL)
  {
    return;
  }
  
  mp = (FeatEdPtr) userdata;
  
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  bsp = BioseqFindFromSeqLoc (sfp->location);

  if (mp->adjust_gene && sfp->data.choice != SEQFEAT_GENE) {
    gene = GetGeneForFeature (sfp);
  }

  if (sfp->idx.subtype == FEATDEF_CDS) 
  {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    if (bs != NULL) 
    {
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot != NULL) 
      {
        first_char = *prot;
        ptr = prot;
        ch = *ptr;
        if (ch == '-') 
        {
          dash_at_end = TRUE;
        }
        while (ch != '\0') 
        {
          /* *ptr = TO_UPPER (ch); */
          if (ch == '*') 
          {
            star_at_end = TRUE;
          } else {
            star_at_end = FALSE;
          }
          ptr++;
          ch = *ptr;
        }
      }
      prot = MemFree (prot);
    }
  }

    /* If requested, convert joins to orders */
    /* or orders to joins.                   */

    hasNulls = LocationHasNullsBetween (sfp->location);
    switch (mp->orderjoinpolicy) 
    {
      case CONVERT_TO_JOIN :
        if (hasNulls) 
        {
          slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, FALSE);
            sfp->location = SeqLocFree (sfp->location);
            sfp->location = slp;
            if (bsp->repr == Seq_repr_seg) 
            {
              slp = SegLocToPartsEx (bsp, sfp->location, FALSE);
              sfp->location = SeqLocFree (sfp->location);
              sfp->location = slp;
              hasNulls = LocationHasNullsBetween (sfp->location);
              sfp->partial = (sfp->partial || hasNulls);
            }
            FreeAllFuzz (sfp->location);
            SetSeqLocPartial (sfp->location, partial5, partial3);
        }
        break;
      case CONVERT_TO_ORDER :
        if (!hasNulls) 
        {
            slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, TRUE);
        sfp->location = SeqLocFree (sfp->location);
            sfp->location = slp;
            if (bsp->repr == Seq_repr_seg) 
            {
              slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
              sfp->location = SeqLocFree (sfp->location);
              sfp->location = slp;
              hasNulls = LocationHasNullsBetween (sfp->location);
              sfp->partial = (sfp->partial || hasNulls);
            }
            FreeAllFuzz (sfp->location);
            SetSeqLocPartial (sfp->location, partial5, partial3);
        }
        break;
      case CONVERT_TO_SINGLE_INTERVAL :
        slp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
        sfp->location = SeqLocFree (sfp->location);
        sfp->location = slp;
          SetSeqLocPartial (sfp->location, partial5, partial3);
      default:
        break;
    }

  /* find first and last location parts in a multi-part location */
  slp = SeqLocFindNext (sfp->location, NULL);
  while (slp != NULL) 
  {
    if (firstSlp == NULL) 
    {
      firstSlp = slp;
    }
    lastSlp = slp;
    slp = SeqLocFindNext (sfp->location, slp);
  }

  /* Set the 5' partial */

  if (firstSlp != NULL) 
  {
    strand = SeqLocStrand (firstSlp);
    atEnd = FALSE;
    if (strand == Seq_strand_minus) 
    {
      if (GetOffsetInBioseq (firstSlp, bsp, SEQLOC_START) == bsp->length - 1)
          {
        atEnd = TRUE;
      }
    } else {
      if (GetOffsetInBioseq (firstSlp, bsp, SEQLOC_START) == 0) 
      {
        atEnd = TRUE;
      }
    }
    switch (mp->leftpolicy) 
    {
      case 1 :
        partial5 = TRUE;
        if (! atEnd && mp->extend5)
        {
          moved5 = TRUE;
          ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
        }
        break;
      case 2 :
        if (atEnd) 
        {
          partial5 = TRUE;
        }
        break;
      case 3 :
        if (crp != NULL && (crp->frame > 1 || first_char != 'M')) 
        {
          partial5 = TRUE;
          if (! atEnd && mp->extend5)
          {
            ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
            moved5 = TRUE;
          }
        }
        break;
      case 4:
        if (crp != NULL && crp->frame > 1)
        {
          partial5 = TRUE;
          if (! atEnd && mp->extend5)
          {
            ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
            moved5 = TRUE;
          }
        }
        break;
      case 5 :
        partial5 = FALSE;
        break;
      case 6 :
        if (! atEnd) 
        {
          partial5 = FALSE;
        }
        break;
      case 7 :
        if (crp != NULL && crp->frame == 1 && first_char == 'M')
        {
          partial5 = FALSE;
        }
        break;
      default :
        break;
    }
  }

  /* Set the 3' partial */

  if (lastSlp != NULL) 
  {
    strand = SeqLocStrand (firstSlp);
    atEnd = FALSE;
    if (strand == Seq_strand_minus) 
    {
      if (GetOffsetInBioseq (lastSlp, bsp, SEQLOC_STOP) == 0) 
      {
        atEnd = TRUE;
      }
    } 
    else 
    {
      if (GetOffsetInBioseq (lastSlp, bsp, SEQLOC_STOP) == bsp->length - 1) 
      {
        atEnd = TRUE;
      }
    }
    switch (mp->rightpolicy) 
    {
      case 1 :
        partial3 = TRUE;
        if (! atEnd && mp->extend3)
        {
          ExtendSeqLocToEnd (sfp->location, bsp, FALSE);
        }
        break;
      case 2 :
        if (atEnd) 
        {
          partial3 = TRUE;
        }
        break;
      case 3 :
        if ((! star_at_end) && crp != NULL) 
        {
          partial3 = TRUE;
          if (! atEnd && mp->extend3)
          {
            ExtendSeqLocToEnd (sfp->location, bsp, FALSE);
          }
        }
        break;
      case 4 :
        partial3 = FALSE;
        break;
      case 5 :
        if (! atEnd) 
        {
          partial3 = FALSE;
        }
        break;
      case 6 :
        if (star_at_end && crp != NULL) 
        {
          partial3 = FALSE;
        }
        break;
      default :
        break;
    }
  }

  SetSeqLocPartial (sfp->location, partial5, partial3);
  sfp->partial = (partial5 || partial3 ||
                  LocationHasNullsBetween (sfp->location));
  if (partial5 && moved5)
  {
       SetBestFrame (sfp);
  }

  /* set protein partials for coding regions */
  if (sfp->idx.subtype == FEATDEF_CDS) 
  {
    ResynchCDSPartials (sfp, NULL);
  }

  if (mp->retranslate) {
    RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, TRUE);
  }
  if (gene != NULL) {
    gene->location = SeqLocFree (gene->location);
    gene->location = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
    SetSeqLocPartial (gene->location, partial5, partial3);
    gene->partial = partial5 || partial3;
  }
}


static Boolean DoesStrandMatch (Int4 strand_choice, Uint1 strand_val)
{
  Boolean rval = FALSE;
  
  switch (strand_choice)
  {
    case FEATED_STRAND_ANY_REVERSE:
      rval = TRUE;
      break;
    case FEATED_STRAND_UNKNOWN:
      if (strand_val == Seq_strand_unknown)
      {
        rval = TRUE;
      }
      break;
    case FEATED_STRAND_PLUS:
      if (strand_val == Seq_strand_plus)
      {
        rval = TRUE;
      }
      break;
    case FEATED_STRAND_MINUS:
      if (strand_val == Seq_strand_minus)
      {
        rval = TRUE;
      }
      break;
    case FEATED_STRAND_BOTH:
      if (strand_val == Seq_strand_both)
      {
        rval = TRUE;
      }
      break;
  }
  return rval;
}

static Uint1 GetNewStrandValue (Int4 strand_choice, Uint1 strand_val)
{
  Uint1 rval = Seq_strand_unknown;
  
  switch (strand_choice)
  {
    case FEATED_STRAND_ANY_REVERSE:
      switch (strand_val)
      {
        case Seq_strand_plus:
          rval = Seq_strand_minus;
          break;
        case Seq_strand_minus:
          rval = Seq_strand_plus;
          break;
        default:
          rval = strand_val;
          break;
      }
      break;
    case FEATED_STRAND_UNKNOWN:
      rval = Seq_strand_unknown;
      break;
    case FEATED_STRAND_PLUS:
      rval = Seq_strand_plus;
      break;
    case FEATED_STRAND_MINUS:
      rval = Seq_strand_minus;
      break;
    case FEATED_STRAND_BOTH:
      rval = Seq_strand_both;
      break;
  }  
  return rval;
}

static void ConvertLocationStrand (SeqLocPtr slp, Int4 fromStrand, Int4 toStrand)
{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL && DoesStrandMatch (fromStrand, sinp->strand)) 
        {
          sinp->strand = GetNewStrandValue (toStrand, sinp->strand);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand))
        {
          spp->strand = GetNewStrandValue (toStrand, spp->strand);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL && DoesStrandMatch (fromStrand, psp->strand)) 
        {
          psp->strand = GetNewStrandValue (toStrand, psp->strand);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          ConvertLocationStrand (loc, fromStrand, toStrand);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand)) 
          {
            spp->strand = GetNewStrandValue (toStrand, spp->strand);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand)) 
          {
            spp->strand = GetNewStrandValue (toStrand, spp->strand);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
}

static void DoStrandFeatureProc (SeqFeatPtr sfp, Pointer userdata)
{
  FeatEdPtr mp;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  mp = (FeatEdPtr) userdata;
  
  ConvertLocationStrand (sfp->location, mp->from_strand_val, mp->to_strand_val);
  
}

static void DoCitationFeatureProc (SeqFeatPtr sfp, Pointer userdata)
{
  FeatEdPtr mp;
  ByteStorePtr    bs;
  BioseqPtr       protbsp;
  CharPtr         seq1 = NULL, seq2 = NULL;
  Boolean         add_citations = TRUE;
  
  if (sfp == NULL || userdata == NULL)
  {
    return;
  }
  
  mp = (FeatEdPtr) userdata;
  
  if (mp->use_explanation_constraint && mp->explanation_constraint_choice_val != 1)
  {
    /* do not continue unless feature matches explanation constraint */
    if (StringHasNoText (sfp->except_text))
    {
      if (mp->explanation_constraint_choice_val != 2)
      {
        add_citations = FALSE;
      }
    }
    else if (mp->explanation_constraint_choice_val > num_exception_explanations + 2)
    {
      add_citations = FALSE;
    }
    else if (StringICmp (sfp->except_text, 
             exception_explanations [mp->explanation_constraint_choice_val - 3]) != 0)
    {
      add_citations = FALSE;
    }
  }

  if (mp->use_product_match)
  {
    /* do not continue unless translation does not match product */
    if (sfp->data.choice != SEQFEAT_CDREGION)
    {
      add_citations = FALSE;
    }
    else if (sfp->product != NULL)
    {
      protbsp = BioseqFindFromSeqLoc (sfp->product);
      seq1 = GetSequenceByBsp (protbsp);
      bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
      if (bs != NULL)
      {
        seq2 = BSMerge (bs, NULL);
      }
      bs = BSFree (bs);
      if (StringICmp (seq1, seq2) == 0)
      {
        add_citations = FALSE;
      }
      seq1 = MemFree (seq1);
      seq2 = MemFree (seq2);
    }
  }

  /* add citations to feature */
  if (add_citations)
  {
    sfp->cit = PubSetFree (sfp->cit);
    sfp->cit = DialogToPointer (mp->citation_list);
  }
  
}


static ValNodePtr MakeRemoveFeatureFieldAction (Uint2 ftype, Int4 legalqual, DialoG constraint_dlg)
{
  ValNodePtr      single_macro = NULL;
  RemoveActionPtr remove;
  AECRActionPtr   act;

  remove = RemoveActionNew();
  remove->field = MakeFeatureFieldField (ftype, legalqual);
  act = AECRActionNew();
  act->constraint = DialogToPointer (constraint_dlg);
  act->action = ValNodeNew (NULL);
  act->action->choice = ActionChoice_remove;
  act->action->data.ptrvalue = remove;
  single_macro = ValNodeNew (NULL);
  single_macro->choice = MacroActionChoice_aecr;
  single_macro->data.ptrvalue = act;
  return single_macro;
}


static ValNodePtr 
MakeApplyFeatureFieldAction 
(Uint2 ftype, Int4 legalqual, CharPtr val, DialoG constraint_dlg, SeqEntryPtr sep)
{
  ValNodePtr      single_macro = NULL;
  ApplyActionPtr  apply;
  AECRActionPtr   act;
  ValNodePtr      sample_list;
  AECRSamplePtr   sample;

  apply = ApplyActionNew();
  apply->field = MakeFeatureFieldField (ftype, legalqual);
  apply->value = StringSave (val);
  act = AECRActionNew();
  act->constraint = DialogToPointer (constraint_dlg);
  act->action = ValNodeNew (NULL);
  act->action->choice = ActionChoice_apply;
  act->action->data.ptrvalue = apply;
  apply->existing_text = ExistingTextOption_replace_old;
  if (sep != NULL) {
    sample_list = GetAECRSampleList (act, sep);
    sample = GetFieldSampleFromList (sample_list, apply->field);
    sample_list = AECRSampleListFree (sample_list);
    if (sample != NULL && sample->num_found > 0) {
      apply->existing_text = TwoStepExistingText (sample->num_found, IsFieldTypeNonText(apply->field), AllowFieldMulti (apply->field));
    }
    sample = AECRSampleFree (sample);
  }
  if (apply->existing_text == 0) {
    act = AECRActionFree (act);
  } else {
    single_macro = ValNodeNew (NULL);
    single_macro->choice = MacroActionChoice_aecr;
    single_macro->data.ptrvalue = act;
  }
  return single_macro;
}


static void 
FeatureEditorDoOneAction 
(SeqEntryPtr     sep, 
 Uint2           ftype, 
 FeatEdPtr mp,
 Int4            action_choice)
{
  Int4            exception_choice;
  ValNodePtr      vnp;
  MsgAnswer       ans;
  Boolean         also_remove_evidence = FALSE;
  InferenceEditPtr iep;
  AECRActionPtr   act = NULL;
  ValNodePtr      sample_list;
  ValNodePtr      single_field;
  AECRSamplePtr   sample;
  ValNodePtr      macro_list = NULL;
  ValNodePtr      object_list = NULL;
  CharPtr         val;
  VisitFeaturesFunc callback_func = NULL;

  switch (action_choice) 
  {
    case FEAT_ED_EXPERIMENT:
      if (TextHasNoText (mp->experiment_text)) {
        macro_list = MakeRemoveFeatureFieldAction (ftype, Feat_qual_legal_experiment, mp->constraint);
      } else {
        val = SaveStringFromText (mp->experiment_text);
        macro_list = MakeApplyFeatureFieldAction (ftype, Feat_qual_legal_experiment, val, mp->constraint, sep);
        val = MemFree (val);
      }
      break;
    case FEAT_ED_INFERENCE:
      iep = DialogToPointer (mp->inference_dlg);
      if (iep->action == eInferenceRemove)
      {
        macro_list = MakeRemoveFeatureFieldAction (ftype, Feat_qual_legal_inference, mp->constraint);
        vnp = MakeApplyFeatureFieldAction (ftype, Feat_qual_legal_evidence, "X", mp->constraint, NULL);
        sample_list = GetAECRSampleList (vnp->data.ptrvalue, sep);
        single_field = MakeFeatureFieldField (ftype, Feat_qual_legal_evidence);
        sample = GetFieldSampleFromList (sample_list, single_field);
        single_field = FieldTypeFree (single_field);
        sample_list = AECRSampleListFree (sample_list);
        if (sample != NULL && sample->num_found > 0) {
          ans = Message (MSG_YNC, "There are evidence qualifiers on your features, which are displayed as inference qualifiers in the flat file.  Would you like to remove them as well?");
          if (ans == ANS_CANCEL)
          {
            macro_list = MacroActionListFree (macro_list);
          } 
          else if (ans == ANS_YES) 
          {
            also_remove_evidence = TRUE;
            mp->evidence_val = 0;
          }
        }
        sample = AECRSampleFree (sample);
      }
      else 
      {
        macro_list = MakeApplyFeatureFieldAction (ftype, Feat_qual_legal_inference, "X", mp->constraint, NULL);
      }
      if (macro_list != NULL)
      {
        object_list = GetObjectListForAECRAction (sep, macro_list->data.ptrvalue);
        macro_list = MacroActionListFree (macro_list);
        for (vnp = object_list; vnp != NULL; vnp = vnp->next) 
        {
          DoInferenceFeatureProc (vnp->data.ptrvalue, iep);
          if (also_remove_evidence) 
          {
            DoEvidenceFeatureProc (vnp->data.ptrvalue, mp);
          }
        }
        object_list = FreeObjectList (object_list);

      }
      iep = InferenceEditFree (iep);

      break;
    case FEAT_ED_PSEUDO:
      if (GetValue (mp->pseudo_action) == 1)
      {
        macro_list = MakeApplyFeatureFieldAction (ftype, Feat_qual_legal_pseudo, "TRUE", mp->constraint, NULL);
      } 
      else
      {
        macro_list = MakeRemoveFeatureFieldAction (ftype, Feat_qual_legal_pseudo, mp->constraint);
      }
      break;
    case FEAT_ED_EVIDENCE:
      macro_list = MakeRemoveFeatureFieldAction (ftype, Feat_qual_legal_evidence, mp->constraint);
      break;
    case FEAT_ED_EXCEPTION:
      mp->do_move_to_comment = GetStatus (mp->move_to_comment);
      mp->product_val = GetValue (mp->product_choice);
      exception_choice = GetValue (mp->explanation_choice);
      if (exception_choice == 1)
      {
        mp->explanation_text = NULL;
      }
      else if (exception_choice > 1 && exception_choice < num_exception_explanations + 1)
      {
        mp->explanation_text = StringSave (exception_explanations [exception_choice - 2]);
      }
      else if (exception_choice == num_exception_explanations + 1)
      {
        mp->explanation_text = SaveStringFromText (mp->other_explanation);
      }
      callback_func = DoExceptionFeatureProc;
      break;
    case FEAT_ED_PARTIAL:
      mp->orderjoinpolicy = GetValue (mp->convert_choice);
      mp->leftpolicy = GetValue (mp->partial5_choice);
      mp->rightpolicy = GetValue (mp->partial3_choice);
      mp->extend5 = GetStatus (mp->extend5_btn);
      mp->extend3 = GetStatus (mp->extend3_btn);
      mp->retranslate = GetStatus (mp->retranslate_btn);
      mp->adjust_gene = GetStatus (mp->adjust_gene_btn);
      callback_func = DoPartialFeatureProc;
      break;
    case FEAT_ED_STRAND:
      mp->from_strand_val = (Uint1) GetValue (mp->from_strand);
      mp->to_strand_val = (Uint1) GetValue (mp->to_strand);
      callback_func = DoStrandFeatureProc;
      break;
    case FEAT_ED_CITATION:
      mp->use_explanation_constraint = GetStatus (mp->explanation_constraint);
      mp->explanation_constraint_choice_val = GetValue (mp->explanation_constraint_choice);
      mp->use_product_match = GetStatus (mp->citation_translation);
      callback_func = DoCitationFeatureProc;
      break;
  }

  if (callback_func != NULL) 
  {
    macro_list = MakeApplyFeatureFieldAction (ftype, Feat_qual_legal_note, "X", mp->constraint, NULL);
    object_list = GetObjectListForAECRAction (sep, macro_list->data.ptrvalue);
    macro_list = MacroActionListFree (macro_list);
    for (vnp = object_list; vnp != NULL; vnp = vnp->next) 
    {
      callback_func (vnp->data.ptrvalue, mp);
    }
    object_list = FreeObjectList (object_list);
  }
  
  if (macro_list != NULL) 
  {
    ApplyMacroToSeqEntryEx (sep, macro_list, NULL, Sequin_GlobalAlign2Seq);
    macro_list = MacroActionListFree (macro_list);
  }

  mp->explanation_text = MemFree (mp->explanation_text);
}

static Boolean FeatureEditorAction (Pointer userdata)
{
  FeatEdPtr mp;
  Uint2           ftype;
  Int4            action_choice;
  SeqEntryPtr     sep;

  mp = (FeatEdPtr) userdata;
  if (mp == NULL)
  {
    return FALSE;
  }

  sep = GetTopSeqEntryForEntityID (mp->input_entityID);
  if (sep == NULL)
  {
    return FALSE;
  }

  WatchCursor ();
  Update ();
  
  ftype = GetFeatureTypeFromFeatureTypeDialog (mp->feature_type);
    
  action_choice = GetValue (mp->action_choice_grp);
  FeatureEditorDoOneAction (sep, ftype, mp, action_choice);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (mp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, mp->input_entityID, 0, 0);  
  return TRUE;
}

NLM_EXTERN void FeatureEditorBaseForm (BaseFormPtr bfp, Int4 first_action)
{
  FeatEdPtr mp;
  WindoW          w;
  GrouP           h, g;

  if (bfp == NULL) return;
  mp = (FeatEdPtr) MemNew (sizeof (FeatEdData));
  if (mp == NULL) return;
  
  w = FixedWindow (-50, -33, -10, -10, "Feature Editor", StdCloseWindowProc);
  SetObjectExtra (w, mp, StdCleanupFormProc);
  mp->form = (ForM) w;
  mp->input_entityID = bfp->input_entityID;
  
  h = HiddenGroup (w, -1, 0, NULL);
  
  g = FeatureEditorActionGroup (h, mp, first_action);
  
  mp->constraint = ComplexConstraintDialog (h, NULL, NULL);
  SetComplexConstraintType (mp->constraint, eComplexConstraintType_string);
  mp->accept_cancel = AcceptCancelDialog (h, FeatureEditorAction, NULL, 
                                          FeatEdClear, 
                                          FeatEdClearText, 
                                          (Pointer)mp, w);
  AlignObjects (ALIGN_CENTER, (HANDLE) g,
                              (HANDLE) mp->constraint,
                              (HANDLE) mp->accept_cancel, NULL);
  
  ChangeFeatureEditorActionGroup (mp->action_choice_grp);
  Show (w);
}
