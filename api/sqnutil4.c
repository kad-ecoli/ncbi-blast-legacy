/*   sqnutil4.c
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
* File Name:  sqnutil4.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   12/27/2007
*
* $Revision: 1.161 $
*
* File Description: 
* This file contains functions for automatically generating definition lines.
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
#include <ncbilang.h>
#include <objfdef.h>
#include <gather.h>
#include <explore.h>
#include <edutil.h>
#include <salutil.h>
#include <tofasta.h>
#include <gbftdef.h>
#include <gbfeat.h>
#include <findrepl.h>
#include <salpacc.h>
#include <salpedit.h>
#include <alignmgr.h>
#include <alignmgr2.h>
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>
#include <tax3api.h>

/* This is a list of the modifiers that are of interest */
/* Note that if you modify the DefLineModifiers array, */ 
/* you should make the corresponding change to the DefLinePos enum. */

ModifierItemGlobalData DefLineModifiers[] = {
  { "Acronym"              , TRUE , ORGMOD_acronym              },
  { "Anamorph"             , TRUE , ORGMOD_anamorph             },
  { "Authority"            , TRUE , ORGMOD_authority            },
  { "Bio-material"         , TRUE,  ORGMOD_bio_material         },
  { "Biotype"              , TRUE , ORGMOD_biotype              },
  { "Biovar"               , TRUE , ORGMOD_biovar               },
  { "Breed"                , TRUE , ORGMOD_breed                },
  { "Cell-line"            , FALSE, SUBSRC_cell_line            },
  { "Cell-type"            , FALSE, SUBSRC_cell_type            },
  { "Chemovar"             , TRUE , ORGMOD_chemovar             },
  { "Chromosome"           , FALSE, SUBSRC_chromosome           },
  { "Clone"                , FALSE, SUBSRC_clone                },
  { "Clone-lib"            , FALSE, SUBSRC_clone_lib            },
  { "Collected-by"         , FALSE, SUBSRC_collected_by         },
  { "Collection-date"      , FALSE, SUBSRC_collection_date      },
  { "Common"               , TRUE , ORGMOD_common               },
  { "Country"              , FALSE, SUBSRC_country              },
  { "Cultivar"             , TRUE , ORGMOD_cultivar             },
  { "Culture-collection"   , TRUE , ORGMOD_culture_collection   },
  { "Dev-stage"            , FALSE, SUBSRC_dev_stage            },
  { "Ecotype"              , TRUE , ORGMOD_ecotype              },
  { "Endogenous-virus-name", FALSE, SUBSRC_endogenous_virus_name},
  { "Environmental-sample" , FALSE, SUBSRC_environmental_sample },
  { "Forma"                , TRUE , ORGMOD_forma                },
  { "Forma-specialis"      , TRUE , ORGMOD_forma_specialis      },
  { "Frequency"            , FALSE, SUBSRC_frequency            },
  { "Genotype"             , FALSE, SUBSRC_genotype             },
  { "Germline"             , FALSE, SUBSRC_germline             },
  { "Group"                , TRUE , ORGMOD_group                },
  { "Haplogroup"           , FALSE, SUBSRC_haplogroup           },
  { "Haplotype"            , FALSE, SUBSRC_haplotype            },
  { "Host"                 , TRUE , ORGMOD_nat_host             },
  { "Identified-by"        , FALSE, SUBSRC_identified_by        },
  { "Isolate"              , TRUE , ORGMOD_isolate              },
  { "Isolation-source"     , FALSE, SUBSRC_isolation_source     },
  { "Lab-host"             , FALSE, SUBSRC_lab_host             },
  { "Lat-lon"              , FALSE, SUBSRC_lat_lon              },
  { "Linkage-group"        , FALSE, SUBSRC_linkage_group        },
  { "Map"                  , FALSE, SUBSRC_map                  },
  { "Mating-type"          , FALSE, SUBSRC_mating_type          },
  { "Metagenomic"          , FALSE, SUBSRC_metagenomic          },
  { "Note-OrgMod"          , TRUE,  ORGMOD_other                },
  { "Note-SubSrc"          , FALSE, SUBSRC_other                },
  { "Pathovar"             , TRUE , ORGMOD_pathovar             },
  { "Plasmid-name"         , FALSE, SUBSRC_plasmid_name         },
  { "Plastid-name"         , FALSE, SUBSRC_plastid_name         },
  { "Pop-variant"          , FALSE, SUBSRC_pop_variant          },
  { "Rearranged"           , FALSE, SUBSRC_rearranged           },
  { "Segment"              , FALSE, SUBSRC_segment              },
  { "Serogroup"            , TRUE , ORGMOD_serogroup            },
  { "Serotype"             , TRUE , ORGMOD_serotype             },
  { "Serovar"              , TRUE , ORGMOD_serovar              },
  { "Sex"                  , FALSE, SUBSRC_sex                  },
  { "Specimen voucher"     , TRUE , ORGMOD_specimen_voucher     },
  { "Strain"               , TRUE , ORGMOD_strain               },
  { "Subclone"             , FALSE, SUBSRC_subclone             },
  { "Subgroup"             , TRUE , ORGMOD_subgroup             },
  { "Sub-species"          , TRUE , ORGMOD_sub_species          },
  { "Substrain"            , TRUE , ORGMOD_substrain            },
  { "Subtype"              , TRUE , ORGMOD_subtype              },
  { "Synonym"              , TRUE , ORGMOD_synonym              },
  { "Teleomorph"           , TRUE , ORGMOD_teleomorph           },
  { "Tissue-lib"           , FALSE, SUBSRC_tissue_lib           },
  { "Tissue-type"          , FALSE, SUBSRC_tissue_type          },
  { "Transgenic"           , FALSE, SUBSRC_transgenic           },
  { "Type"                 , TRUE , ORGMOD_type                 },
  { "Variety"              , TRUE , ORGMOD_variety              }
};

#define numDefLineModifiers (sizeof (DefLineModifiers) / sizeof (ModifierItemGlobalData))

NLM_EXTERN size_t NumDefLineModifiers (void)

{
  return numDefLineModifiers;
}

NLM_EXTERN CharPtr MergeValNodeStrings (ValNodePtr list, Boolean useReturn)

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
        tmp = NULL;
        while (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            if (tmp == NULL) {
              tmp = ptr;
            } else if (useReturn) {
              tmp = StringMove (tmp, "\n");
            } else if (IsJapanese () && (tmp - ptr > 2) &&
                    IsMBLetter (tmp - 2) && IsMBLetter (str)) {
              /* no space required between two Japanese letters. */
              tmp = tmp;
            } else if (str [0] != ',' && str [0] != ';' && str [0] != ':') {
              tmp = StringMove (tmp, " ");
            } else {
              tmp = StringMove (tmp, " ");
            }
            tmp = StringMove (tmp, str);
          }
          vnp = vnp->next;
        }
      }
    }
  }
  return ptr;
}


/* The matchFunction functions are used to identify features that meet
 * specific requirements, usually that the feature is of a particular type.
 * This function is used instead of simply using the subtype for the feature
 * because some features are identified based on the contents or presence of
 * certain modifiers.
 * Functions of this type should always return FALSE if handed a NULL argument.
 */
typedef Boolean (LIBCALLBACK *matchFunction) (
  SeqFeatPtr sfp
);

static void ListClauses (
  ValNodePtr clauselist,
  ValNodePtr PNTR strings,
  Boolean    allow_semicolons,
  Boolean    suppress_final_and,
  Boolean    suppress_allele
);

static void LabelClauses 
( ValNodePtr clause_list,
  Uint1      biomol,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp);

static CharPtr GetProductName 
( SeqFeatPtr cds,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp);

#define DEFLINE_FEATLIST    1
#define DEFLINE_CLAUSEPLUS  2
#define DEFLINE_REMOVEFEAT  3

typedef struct featurelabeldata {
  Boolean pluralizable;
  Boolean is_typeword_first;
  CharPtr typeword;
  CharPtr description;
  CharPtr productname;
} FeatureLabelData, PNTR FeatureLabelPtr;


typedef struct featureclause {
  ValNodePtr       featlist;
  FeatureLabelData feature_label_data;
  CharPtr          allelename;
  CharPtr          interval;
  Boolean          is_alt_spliced;
  Boolean          has_mrna;
  SeqLocPtr        slp;
  GeneRefPtr       grp;
  Boolean          clause_info_only;
  Boolean          is_unknown;
  Boolean          make_plural;
  Boolean          delete_me;
  /* this information used only for segments */
  Int2             numivals;
  Int4Ptr          ivals;
} FeatureClauseData, PNTR FeatureClausePtr;

FeatureClausePtr NewFeatureClause (
  SeqFeatPtr sfp,
  BioseqPtr bsp,
  DeflineFeatureRequestListPtr rp);

static void PluralizeConsolidatedClauseDescription (
  FeatureClausePtr fcp
);

typedef Boolean (LIBCALLBACK *ShouldRemoveFunction) (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr        bsp,
  Boolean          isLonely,
  Boolean          isRequested,
  Boolean          isSegment,
  DeflineFeatureRequestListPtr rp
);

/* This section of the code contains some functions for dealing with
 * linked lists of strings */

/* This function finds the first occurrence of "search" in one of the
 * strings in list "strings".
 * "search" could be part of the string or could be the entire string.
 */
static ValNodePtr FindStringInStrings (
  ValNodePtr strings,
  CharPtr search
)
{
  while (strings != NULL)
  {
    if (DoesStringContainPhrase (strings->data.ptrvalue, search, TRUE, TRUE))
    {
      return strings;
    }
    strings = strings->next;
  }
  return NULL;
}


NLM_EXTERN ValNodePtr FindExactStringListMatch (
  ValNodePtr list,
  CharPtr value
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringCmp (str, value) == 0) return vnp;
  }

  return NULL;
}

/* This function creates a new linked list of strings with copies of
 * contents of orig.
 */
static ValNodePtr CopyStrings (
  ValNodePtr orig
)
{
  ValNodePtr new_string_start = NULL;

  while (orig != NULL)
  {
    ValNodeAddStr (&new_string_start, 0,
      StringSave (orig->data.ptrvalue));
    orig = orig->next;
  }
  return new_string_start;
}

/*
 * This section of the code contains functions and structures for obtaining a
 * description of the organism in the record, including functions for finding
 * the combination of modifiers that will make each organism description 
 * unique.
 * The method used for determining the best combination of modifiers involves
 * creating a list of required modifiers, and then creating a list of
 * combinations of modifiers by adding modifiers one at a time
 * to see if the additional modifiers provide any more differentiation in
 * the list.
 * In order to do this, I start with a list of required modifiers, and 
 * then create copies of this list.  For each copy I add one of the modifiers
 * that are present in the bio sources and not already on the list.
 * If adding the modifier increases the differentiation, I add that copy to
 * the list of possible combinations, otherwise I discard it.
 * I then make copies of all of the new items I added to the list and
 * add another modifier to each list, keeping the combinations that increase
 * the differentiation and discarding the rest.
 * This process continues until I have a combination that produces completely
 * differentiated bio sources, or I run out of possible combinations.
 * If I run out of possible combinations, I select the best combination from
 * the list.
 * This search process occurs in FindBestCombo.  The majority of the functions
 * in this section are here to support FindBestCombo, specifically to create,
 * copy, and grow lists of combinations. 
 */

/* BioSrcDescData is used to calculate the best possible combination of
 * source and organism modifiers for uniqueness.
 * biop contains the BioSourcePtr from a sequence in the record.
 * strings contains a list of string representations of the modifiers
 * for this combination for this organism.
 */
typedef struct biosrcdescdata {
  BioSourcePtr  biop;
  ValNodePtr    strings;
  Pointer       next;
} BioSrcDescData, PNTR BioSrcDescPtr;

/* OrgGroupData is used to calculate the best possible combination of
 * source and organism modifiers for uniqueness.
 * org_list is a list of all organisms that have identical descriptions
 * using the current set of modifiers.
 * num_organisms contains the number of organisms with identical descriptions.
 */
typedef struct orggroupdata {
  BioSrcDescPtr org_list;
  Int4          num_organisms;
  Pointer       next;
} OrgGroupData, PNTR OrgGroupPtr;

/* ModifierCombinationData is used to calculate the best possible combination
 * of source and organism modifiers for uniqueness.
 * num_groups is the number of groups of organisms with identical descriptions
 *           using the modifiers specified in modifier_indices.
 * num_mods is the number of modifiers specified in modifier_indices.
 * max_orgs_in_group is the maximum number of organisms in any one group.
 * num_unique_orgs is the number of organisms that are alone in their groups
 *           i.e., their description is unique.
 * modifier_indices is the list of modifier indices for this combination.
 * group_list is the list of groups of organisms with identical descriptions
 *           using the modifiers specified in modifier_indices.
 */
typedef struct modifiercombination {
  Int4         num_groups;
  Int4         num_mods;
  Int4         max_orgs_in_group;
  Int4         num_unique_orgs;
  ValNodePtr   modifier_indices;
  OrgGroupPtr  group_list;
  Pointer      next;
} ModifierCombinationData, PNTR ModifierCombinationPtr;

static Boolean IsDeflineModifierRequiredByDefault (Boolean is_orgmod, Int2 index)
{
  if (!is_orgmod
      && (index == SUBSRC_endogenous_virus_name
          || index == SUBSRC_plasmid_name
          || index == SUBSRC_transgenic)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void AddOneSubtypeField (ValNodePtr PNTR sq_list, SourceQualDescPtr orig, CharPtr str, Uint1 subfield)
{
  SourceQualDescPtr sqdp_cpy;

  if (sq_list == NULL || orig == NULL) {
    return;
  }
  sqdp_cpy = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
  MemCpy (sqdp_cpy, orig, sizeof (SourceQualDescData));

  sqdp_cpy->name = str;
  sqdp_cpy->subfield = subfield;

  ValNodeAddPointer (sq_list, 0, sqdp_cpy);  
}


static void AddSubtypeFields (ValNodePtr PNTR sq_list, SourceQualDescPtr orig)
{
  if (sq_list == NULL || orig == NULL) return;

  if (orig->isOrgMod) {
    switch (orig->subtype) {
      case ORGMOD_specimen_voucher:
        AddOneSubtypeField (sq_list, orig, "specimen-voucher INST", 1);
        AddOneSubtypeField (sq_list, orig, "specimen-voucher COLL", 2);
        AddOneSubtypeField (sq_list, orig, "specimen-voucher SpecID", 3);
        break;  
      case ORGMOD_culture_collection:
        AddOneSubtypeField (sq_list, orig, "culture-collection INST", 1);
        AddOneSubtypeField (sq_list, orig, "culture-collection COLL", 2);
        AddOneSubtypeField (sq_list, orig, "culture-collection SpecID", 3);
        break;  
      case ORGMOD_bio_material:
        AddOneSubtypeField (sq_list, orig, "bio-material INST", 1);
        AddOneSubtypeField (sq_list, orig, "bio-material COLL", 2);
        AddOneSubtypeField (sq_list, orig, "bio-material SpecID", 3);
        break; 
    }
  } 
}


static void AddQualList (ValNodePtr PNTR list, Nlm_QualNameAssocPtr qual_list, Boolean is_orgmod, Boolean use_alternate_note_name, Boolean get_subfields)
{
  Int4              k;
  SourceQualDescPtr sqdp;
  
  for (k = 0; qual_list[k].name != NULL; k++) {
    if (StringHasNoText (qual_list[k].name)) {
      continue;
    } 
    sqdp = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
    if (sqdp != NULL)
    {
      if (use_alternate_note_name 
          && ((is_orgmod && qual_list[k].value == ORGMOD_other)
              || (!is_orgmod && qual_list[k].value == SUBSRC_other))) 
      {
        if (is_orgmod) {
          sqdp->name = "Note -- OrgMod";
        } else {
          sqdp->name = "Note -- SubSource";
        }
      } else {
        sqdp->name = qual_list[k].name;
      }
      sqdp->isOrgMod = is_orgmod;
      sqdp->subtype = qual_list[k].value;
      sqdp->subfield = 0;
      ValNodeAddPointer (list, 0, sqdp);
    }
    if (get_subfields) {
      AddSubtypeFields (list, sqdp);
    }
  }
}


static void AddNoteQual (ValNodePtr PNTR list, Boolean is_orgmod, Boolean use_alternate_note_name)
{
  SourceQualDescPtr sqdp;

  if (list == NULL) return;

  sqdp = (SourceQualDescPtr) MemNew (sizeof (SourceQualDescData));
  if (sqdp != NULL)
  {
    if (use_alternate_note_name) 
    {
      if (is_orgmod) {
        sqdp->name = "Note -- OrgMod";
      } else {
        sqdp->name = "Note -- SubSource";
      }
    } else {
      sqdp->name = "Note";
    }
    sqdp->isOrgMod = is_orgmod;
    if (is_orgmod) {
      sqdp->subtype = ORGMOD_other;
    } else {
      sqdp->subtype = SUBSRC_other;
    }
    sqdp->subfield = 0;
    ValNodeAddPointer (list, 0, sqdp);
  }
}


NLM_EXTERN int LIBCALLBACK SortVnpBySourceQualDesc (VoidPtr ptr1, VoidPtr ptr2)

{
  SourceQualDescPtr     str1;
  SourceQualDescPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (SourceQualDescPtr) vnp1->data.ptrvalue;
      str2 = (SourceQualDescPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL
          && str1->name != NULL && str2->name != NULL) {
        return StringICmp (str1->name, str2->name);
      }
    }
  }
  return 0;
}


extern ValNodePtr GetSourceQualDescListEx (Boolean get_subsrc, Boolean get_orgmod, Boolean get_discouraged, Boolean get_discontinued, Boolean get_subfields)
{
  ValNodePtr        source_qual_list = NULL;

  if (get_orgmod) {
    AddQualList (&source_qual_list, current_orgmod_subtype_alist, TRUE, get_subsrc, get_subfields);
    if (get_discouraged) {
      AddQualList (&source_qual_list, discouraged_orgmod_subtype_alist, TRUE, get_subsrc, get_subfields);
    }
    if (get_discontinued) {
      AddQualList (&source_qual_list, discontinued_orgmod_subtype_alist, TRUE, get_subsrc, get_subfields);
    }
    AddNoteQual (&source_qual_list, TRUE, get_subsrc);
  }
  if (get_subsrc) {
    AddQualList (&source_qual_list, current_subsource_subtype_alist, FALSE, get_orgmod, get_subfields);
    if (get_discouraged) {
      AddQualList (&source_qual_list, discouraged_subsource_subtype_alist, FALSE, get_orgmod, get_subfields);
    }
    if (get_discontinued) {
      AddQualList (&source_qual_list, discontinued_subsource_subtype_alist, FALSE, get_orgmod, get_subfields);
    }
    AddNoteQual (&source_qual_list, FALSE, get_orgmod);
  }

  source_qual_list = ValNodeSort (source_qual_list, SortVnpBySourceQualDesc);
  return source_qual_list;
}

extern ValNodePtr GetSourceQualDescList (Boolean get_subsrc, Boolean get_orgmod, Boolean get_discouraged, Boolean get_discontinued)
{
  return GetSourceQualDescListEx (get_subsrc, get_orgmod, get_discouraged, get_discontinued, TRUE);
}

/*
 * The CountModifiersProc is used as the callback function for
 * VisitBioSourcesInSep when we are getting a list of all the modifiers
 * that appear in the sources.  We also obtain, for each modifier class,
 * the first value seen, whether or not each value seen is unique for
 * for the modifier, and whether or not the modifier is present for all
 * sources.
 */
static void CountModifiersProc (
  BioSourcePtr biop,
  Pointer userdata
)
{
  ModifierItemLocalPtr ItemList;
  OrgModPtr     mod;
  SubSourcePtr  ssp;
  Int2 i;
  Boolean       found_this_modifier;

  if (biop == NULL) return;
  ItemList = (ModifierItemLocalPtr) userdata;

  for (i=0; i < numDefLineModifiers; i++)
  {
    found_this_modifier = FALSE;
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
        if (mod != NULL && mod->subname != NULL)
        {
          found_this_modifier = TRUE;
          if (ItemList[i].first_value_seen != NULL)
          {
            if (StringCmp (ItemList[i].first_value_seen, mod->subname) != 0)
            {
              ItemList[i].is_unique = FALSE;
            }
          }
          else
          {
            ItemList[i].first_value_seen = mod->subname;
          }
          if ( FindExactStringListMatch (ItemList[i].values_seen, mod->subname)
            == NULL)
          {
            ValNodeAddStr (&ItemList[i].values_seen, 0, mod->subname);
          }
          else
          {
            ItemList[i].all_unique = FALSE;
          }
        }
      }
    } else {
      ssp = biop->subtype;
      while (ssp != NULL && ssp->subtype != DefLineModifiers[i].subtype)
      {
        ssp = ssp->next;
      }
      if (ssp != NULL && ssp->name != NULL)
      {
        found_this_modifier = TRUE;
        if (ItemList[i].first_value_seen != NULL)
        {
          if (StringCmp (ItemList[i].first_value_seen, ssp->name) != 0)
          {
            ItemList[i].is_unique = FALSE;
          }
        }
        else
        {
          ItemList[i].first_value_seen = ssp->name;
        }
        if ( FindExactStringListMatch (ItemList[i].values_seen, ssp->name)
          == NULL)
        {
          ValNodeAddStr (&ItemList[i].values_seen, 0, ssp->name);
        }
        else
        {
          ItemList[i].all_unique = FALSE;
        }
      }
    }
    if (found_this_modifier)
    {
      ItemList[i].any_present = TRUE;
    } else {
      ItemList[i].all_present = FALSE;
    }
  }
} 

/* The CountModifiers function visits all of the bio sources, determining
 * which modifiers are present, which modifiers have only one value,
 * which modifiers have all different values, and which modifiers are
 * present in all sources.
 * After this survey is complete, the function prepares a short summary
 * of the above information for each modifier, which is used in the 
 * definition line options dialog.
 */ 
NLM_EXTERN void CountModifiers (
  ModifierItemLocalPtr ItemList,
  SeqEntryPtr sep
)
{
  Int2 i;

  for (i=0; i < numDefLineModifiers; i++)
  {
    ItemList[i].all_present = TRUE;
    ItemList[i].is_unique = TRUE;
    ItemList[i].first_value_seen = NULL;
    ItemList[i].values_seen = NULL;
    ItemList[i].all_unique = TRUE;
  }

  VisitBioSourcesInSep (sep, ItemList, CountModifiersProc);
 
  for (i=0; i < numDefLineModifiers; i++)
  {
    if (ItemList[i].all_present && ItemList[i].all_unique)
    {
      ItemList[i].status = "All present, all unique";
    }
    else if (ItemList[i].all_present && ItemList[i].is_unique)
    {
      ItemList[i].status = "All present, one unique";
    }
    else if (ItemList[i].all_present && ! ItemList[i].is_unique)
    {
      ItemList[i].status = "All present, mixed";
    }
    else if (! ItemList[i].all_present && ItemList[i].all_unique)
    {
      ItemList[i].status = "Some missing, all unique";
    }
    else if (! ItemList[i].all_present && ItemList[i].is_unique)
    {
      ItemList[i].status = "Some missing, one unique";
    }
    else if (! ItemList[i].all_present && ! ItemList[i].is_unique)
    {
      ItemList[i].status = "Some missing, mixed";
    }
  }
}

/* The BioSrcDescData structure is used to hold a BioSourcePtr, a list
 * of strings used to describe the biosource, including the taxonomy name
 * and the values of all of the modifiers selected so far for this bio
 * source, and a pointer to the next BioSrcDescData structure in the list.
 */

/* The CopyBioSrcDescPtr function creates a copy of the linked list of
 * BioSrcDescData structures.
 */
static BioSrcDescPtr CopyBioSrcDescPtr (
  BioSrcDescPtr orig
)
{
  BioSrcDescPtr new_bsdp_start;

  if (orig == NULL) return NULL;

  new_bsdp_start = (BioSrcDescPtr) MemNew (sizeof (BioSrcDescData));
  if (new_bsdp_start == NULL) return NULL;

  new_bsdp_start->biop = orig->biop;
  new_bsdp_start->strings = CopyStrings (orig->strings);
  new_bsdp_start->next = CopyBioSrcDescPtr (orig->next);
  return new_bsdp_start;
}

/* The FreeBioSrcDescPtr function frees the memory associated with a
 * linked list of BioSrcDescData structures.
 */
static void FreeBioSrcDescPtr (
  BioSrcDescPtr bsdp
)
{
  if (bsdp == NULL) return;
  FreeBioSrcDescPtr (bsdp->next);
  bsdp->biop = NULL;
  ValNodeFreeData (bsdp->strings);
  MemFree (bsdp);
}

/* The AddQualToBioSrcDescPtr function finds the qualifier at the
 * feature_index position in the DefLineModifiers array in the
 * BioSourcePtr and adds the value for that modifier to the array
 * of strings describing the bio source.
 */
static void AddQualToBioSrcDescPtr (
  BioSrcDescPtr bsdp,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  OrgModPtr          mod;
  SubSourcePtr       ssp;
  CharPtr            tmp;

  if (bsdp == NULL) return;
  if (bsdp->biop == NULL) return;

  if (DefLineModifiers[feature_index].isOrgMod)
  {
    if (bsdp->biop->org == NULL || bsdp->biop->org->orgname == NULL) return;
    mod = bsdp->biop->org->orgname->mod;
    while (mod != NULL
        && mod->subtype != DefLineModifiers[feature_index].subtype)
    {
      mod = mod->next;
    }
    if (mod != NULL && mod->subname != NULL)
    {
      if (mod->subtype == ORGMOD_specimen_voucher && StringNICmp (mod->subname, "personal:", 9) == 0)
      {
        tmp = mod->subname + 9;
        while (isspace (*tmp)) 
        {
          tmp++;
        }
      }
      else
      {
        tmp = mod->subname;
      }
      ValNodeCopyStr( &(bsdp->strings), 0, tmp);
    } 
  } else {
    ssp = bsdp->biop->subtype;
    while (ssp != NULL
        && ssp->subtype != DefLineModifiers[feature_index].subtype)
    {
      ssp = ssp->next;
    }
    if (ssp != NULL)
    {
      if (ssp->subtype == SUBSRC_transgenic)
      {
        ValNodeCopyStr( &(bsdp->strings), 0, "transgenic");
      }
      else if (ssp->name != NULL)
      {
        ValNodeCopyStr( &(bsdp->strings), 0, ssp->name);
      }
    }
  }
}
 
/* The CompareOrganismDescriptors function compares the contents of the
 * lists of strings for each BioSrcDesc item.
 * The function returns:
 *     -1 if org1 < org2
 *      0 if org1 = org2
 *      1 if org1 > org2
 */
static int CompareOrganismDescriptors (
  BioSrcDescPtr org1,
  BioSrcDescPtr org2
)
{
  ValNodePtr vnp1, vnp2;
  int cmpval;

  vnp1 = org1->strings;
  vnp2 = org2->strings;

  while (vnp1 != NULL && vnp2 != NULL)
  {
    cmpval = StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    if (cmpval != 0) return cmpval;

    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return 0;
  }
  else if (vnp1 != NULL && vnp2 == NULL)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

/* The OrgGroupData structure contains a list of BioSrcDescData items
 * for which the contents of the descriptive strings list are identical,
 * i.e., all the organisms in the group would have the same description
 * if you used the modifiers used to generate this list of strings.
 * The structure also contains the number of organisms in the list
 * so that it will be easy to tell that the OrgGroup now contains a
 * single organism with a unique description.
 */

/* The CopyOrgGroupList function creates a copy of the list of OrgGroups */
static OrgGroupPtr CopyOrgGroupList (
  OrgGroupPtr orig
)
{
  OrgGroupPtr new_ogp_start = NULL, new_ogp;
 
  if (orig == NULL) return NULL;

  new_ogp_start = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
  if (new_ogp_start == NULL) return NULL;

  new_ogp_start->num_organisms = orig->num_organisms;
  new_ogp_start->org_list = CopyBioSrcDescPtr (orig->org_list);
  new_ogp_start->next = NULL;
  orig = orig->next;
  new_ogp = new_ogp_start;
  while (orig != NULL) {
    new_ogp->next = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
    new_ogp = new_ogp->next;
    new_ogp->num_organisms = orig->num_organisms;
    new_ogp->org_list = CopyBioSrcDescPtr (orig->org_list);
    new_ogp->next = NULL;
    orig = orig->next;
  }

  return new_ogp_start;
}

/* The FreeOrgGroupPtr function frees the memory associated with a
 * list of OrgGroups */
static void FreeOrgGroupPtr (
  OrgGroupPtr ogp
)
{
  OrgGroupPtr ogp_next;

  while (ogp != NULL) {
    ogp_next = ogp->next;
    FreeBioSrcDescPtr (ogp->org_list);
    ogp = MemFree (ogp);
    ogp = ogp_next;
  }
  return;
}

/* The ReorderGroupOrgs function sorts the OrgGroup list based on the results
 * of the CompareOrganismDescriptors function.
 */
static void ReorderGroupOrgs (
  OrgGroupPtr this_group
)
{
  BioSrcDescPtr bsdp;
  BioSrcDescPtr nextBsdp;
  BioSrcDescPtr prevBsdp;
  Boolean swap_needed = TRUE;

  if (this_group->org_list == NULL) return;
  if (this_group->org_list->next == NULL) return;

  while (swap_needed)
  {
    swap_needed = FALSE;
    bsdp = this_group->org_list;
    prevBsdp = NULL;
    while (bsdp->next != NULL)
    {
      nextBsdp = bsdp->next;
      if (CompareOrganismDescriptors (bsdp, nextBsdp) > 0)
      {
        swap_needed = TRUE;
        bsdp->next = nextBsdp->next;
        nextBsdp->next = bsdp;
        if (prevBsdp == NULL)
        {
          this_group->org_list = nextBsdp;
        }
        else
        {
          prevBsdp->next = nextBsdp;
        }
        prevBsdp = nextBsdp;
      }
      else
      {
        prevBsdp = bsdp;
        bsdp = bsdp->next;
      }
    }
  }
}

/* The ReGroupOrgs function operates on a single OrgGroup item.
 * If any of the BioSrcDesc items in the group now have different
 * descriptions, the function breaks it up into smaller, homogenous OrgGroups.
 */
static void ReGroupOrgs (
  OrgGroupPtr this_group
)
{
  BioSrcDescPtr bsdp;
  OrgGroupPtr new_group;
  int num_organisms;

  if (this_group == NULL) return;
  bsdp = this_group->org_list;
  if (bsdp == NULL) return;
  num_organisms = 0;
  while (bsdp->next != NULL)
  {
    num_organisms ++;
    if (CompareOrganismDescriptors (bsdp, bsdp->next) != 0)
    {
      /* create new group to hold next set of organisms */
      new_group = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
      if (new_group == NULL) return;
      new_group->org_list = bsdp->next;
      new_group->num_organisms = this_group->num_organisms - num_organisms;
      new_group->next = this_group->next;
      this_group->next = new_group;
      this_group->num_organisms = num_organisms;
      bsdp->next = NULL;
      ReGroupOrgs (new_group);
    }
    else
    {
      bsdp = bsdp->next;
    }
  }
}

/* The AddQualToGroup function operates on a single OrgGroup item.
 * The function adds a qualifier to each BioSrcDesc item in the OrgGroup,
 * breaks the group into multiple groups if the group is no longer
 * homogenous, and sorts the new list.
 */
static void AddQualToGroup (
  OrgGroupPtr this_group,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  BioSrcDescPtr bsdp;

  if (this_group == NULL) return;

  bsdp = this_group->org_list;
  while (bsdp != NULL)
  {
    AddQualToBioSrcDescPtr (bsdp, qual, feature_index);
    bsdp= bsdp->next;
  }

  /* now reorder organisms and break up group */
  ReorderGroupOrgs (this_group);

  ReGroupOrgs (this_group);
}

/* The AddQualToGroupList function operates on a list of OrgGroup items.
 * It calls AddQualToGroup for each item in the list.
 */
static void AddQualToGroupList (
  OrgGroupPtr group_list,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  OrgGroupPtr ogp;

  ogp = group_list;
  while (ogp != NULL)
  {
    AddQualToGroup (ogp, qual, feature_index);
    ogp = ogp->next;
  }
}

/* The CopyModifierIndices function creates a new ValNode list with the
 * same data.intvalue values for each node as the original modifier_indices
 * ValNode list.
 */
static ValNodePtr CopyModifierIndices (
  ValNodePtr modifier_indices
)
{
  ValNodePtr new_indices;

  if (modifier_indices == NULL) return NULL;
  new_indices = ValNodeNew (NULL);
  if (new_indices == NULL) return NULL;
  new_indices->choice = modifier_indices->choice;
  new_indices->data.intvalue = modifier_indices->data.intvalue;
  new_indices->next = CopyModifierIndices (modifier_indices->next);
  return new_indices;
}
 
/* The CopyModifierCombo creates a copy of a ModificationCombination item.
 * This includes creating a copy of the number and list of modifiers
 * and a copy of the number and list of OrgGroups, as well as copying the
 * maximum number of organisms in any one group and the number of unique
 * organism descriptions produced by this combination of modifiers.
 */
static ModifierCombinationPtr CopyModifierCombo (
  ModifierCombinationPtr m
)
{
  ModifierCombinationPtr newm;
  ValNodePtr  vnp;
  ValNodePtr  newval;

  newm = (ModifierCombinationPtr) MemNew (sizeof (ModifierCombinationData));
  if (newm == NULL) return NULL;

  newm->next = NULL;

  /* copy list of modifier indices */
  newm->num_mods = m->num_mods;
  newm->modifier_indices = NULL;
  vnp = m->modifier_indices;
  if (vnp != NULL)
  {
    newm->modifier_indices = ValNodeNew (NULL);
    if (newm->modifier_indices == NULL) return NULL;
    newm->modifier_indices->data.intvalue = vnp->data.intvalue;
    vnp = vnp->next;
    while (vnp != NULL)
    {
      newval = ValNodeNew (newm->modifier_indices);
      if (newval == NULL) return NULL;
      newval->data.intvalue = vnp->data.intvalue;
      vnp = vnp->next;
    }
  }
  
  /* copy groups */
  newm->num_groups = m->num_groups;
  newm->group_list = CopyOrgGroupList (m->group_list);
 
  return newm; 
}

/* This function creates a new ModifierCombination item using the supplied
 * OrgGroup list.  It calculates the number of groups, maximum number of 
 * organisms in any one group, and number of unique organisms.
 * Initially there are no modifiers.
 */
static ModifierCombinationPtr NewModifierCombo (
  OrgGroupPtr group_list
)
{
  ModifierCombinationPtr newm;
  OrgGroupPtr  ogp;

  newm = (ModifierCombinationPtr) MemNew (sizeof (ModifierCombinationData));
  if (newm == NULL) return NULL;

  newm->num_mods = 0;
  newm->modifier_indices = NULL;
  newm->num_unique_orgs = 0;
  
  /* copy groups */
  newm->group_list = CopyOrgGroupList (group_list);
 
  ogp = newm->group_list;
  newm->max_orgs_in_group = 0;
  newm->num_groups = 0;
  while (ogp != NULL)
  {
    if (newm->max_orgs_in_group < ogp->num_organisms)
      newm->max_orgs_in_group = ogp->num_organisms;
    if (ogp->num_organisms == 1)
      newm->num_unique_orgs ++;
    newm->num_groups ++;
    ogp = ogp->next;
  }

  newm->next = NULL;
  return newm; 
}

/* This function frees the memory associated with a list of
 * ModifierCombination items.
 */
static void FreeModifierCombo (
  ModifierCombinationPtr m
)
{
  if (m == NULL) return;
  FreeModifierCombo (m->next);
  ValNodeFree (m->modifier_indices);
  FreeOrgGroupPtr (m->group_list);
  MemFree (m);
}

/* This function adds the qualifier at the feature_index position in the
 * DefLineModifiers array to each OrgGroup in the list and recalculates
 * the maximum number of organisms in any one group and the number of
 * unique organism descriptions generated by this new combination of 
 * modifiers.
 */
static void AddQualToModifierCombo (
  ModifierCombinationPtr m,
  ModifierItemLocalPtr qual,
  Int2 feature_index
)
{
  OrgGroupPtr ogp;
  ValNodePtr vnp;

  if (m == NULL) return;

  /* now try adding the modifier, see if the number of groups goes up */
  /* if the number of organisms in each group is one, we can stop */
  vnp = ValNodeNew (m->modifier_indices);
  if (vnp == NULL) return;
  if (m->modifier_indices == NULL)
  {
    m->modifier_indices = vnp;
  }
  vnp->data.intvalue = feature_index;
  m->num_mods ++;
  AddQualToGroupList (m->group_list, qual, feature_index);
  ogp = m->group_list;
  m->max_orgs_in_group = 0;
  m->num_unique_orgs = 0;
  m->num_groups = 0;
  while (ogp != NULL)
  {
    if (m->max_orgs_in_group < ogp->num_organisms)
      m->max_orgs_in_group = ogp->num_organisms;
    if (ogp->num_organisms == 1)
      m->num_unique_orgs ++;
    m->num_groups ++;
    ogp = ogp->next;
  }
}

/* This function creates the initial OrgGroup list that is copied for every
 * ModifierCombination item.
 */
static void BuildTaxOrgGroupList (
  BioSourcePtr biop,
  Pointer userdata
)
{
  OrgGroupPtr   ogp;
  OrgGroupPtr   prevOgp;
  OrgGroupPtr PNTR pogp;
  BioSrcDescPtr newBsdp;
  OrgRefPtr     orp;
  int cmpval;

  pogp = (OrgGroupPtr PNTR) userdata;
  ogp = *pogp;

  newBsdp = (BioSrcDescPtr) MemNew (sizeof (BioSrcDescData));
  if (newBsdp == NULL) return;
  newBsdp->biop = biop;
  newBsdp->next = NULL;
  newBsdp->strings = NULL;

  /* add tax name as first string */
  /* later, move this into a separate function and add special handling */
  orp = biop->org;
  if (orp != NULL && orp->taxname != NULL)
  {
    ValNodeCopyStr (&(newBsdp->strings), 0, orp->taxname);
  }

  prevOgp = NULL;
  cmpval = -1;
  while (ogp != NULL && cmpval < 0)
  {
    if (ogp->org_list != NULL)
    {
      cmpval = CompareOrganismDescriptors (ogp->org_list, newBsdp);
      if (cmpval == 0)
      {
        newBsdp->next = ogp->org_list;
        ogp->org_list = newBsdp;
        ogp->num_organisms ++;
      }
    }
    if (cmpval < 0)
    {
      prevOgp = ogp;
      ogp = ogp->next;
    }
  }
  if (cmpval != 0)
  {
    /* create new group */
    ogp = (OrgGroupPtr) MemNew (sizeof (OrgGroupData));
    if (ogp == NULL) return;
    ogp->org_list = newBsdp;
    ogp->num_organisms = 1;
    ogp->next = NULL;
    if (prevOgp == NULL)
    {
      ogp->next = *pogp;
      *pogp = ogp;
    }
    else
    {
      ogp->next = prevOgp->next;
      prevOgp->next = ogp;
    }
  }
}

typedef struct bestsortdata {
  Int4    feature_index;
  Boolean all_unique;
  Boolean all_present;
  Boolean is_unique;
} BestSortData, PNTR BestSortPtr;

static Boolean Index1FoundBeforeIndex2 (
  Int4 index1,
  Int4 index2,
  ValNodePtr list
)
{
  ValNodePtr  vnp;
  BestSortPtr bsp;
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if ((bsp = vnp->data.ptrvalue) == NULL)
    {
      continue;
    }
    if (bsp->feature_index == index1) return TRUE;
    if (bsp->feature_index == index2) return FALSE;
  }
  return FALSE;
}

/* This function determines whether or not we should try adding this modifier
 * to our combination.  If we've already tried it and not added it to the list,
 * there's no reason to try adding it again.
 */
static Boolean OkToTryAddingQual (
  ModifierCombinationPtr m,
  ModifierItemLocalPtr ItemList,
  ValNodePtr           available_modifiers_list,
  Int2 feature_index
)
{
  ValNodePtr vnp;

  /* if feature_index indicates a value we don't use for best combos, skip */
  if (feature_index == DEFLINE_POS_Map || feature_index == DEFLINE_POS_Specific_host)
  {
    return FALSE;
  }

  if (m == NULL) return TRUE;

  /* if feature_index is lower than anything else on list (other than */
  /* a required value, this is a repeat combination, so skip it */
  vnp = m->modifier_indices;
  while (vnp != NULL)
  {
    if (feature_index == m->modifier_indices->data.intvalue)
      return FALSE;
    if (! ItemList[m->modifier_indices->data.intvalue].required &&
      Index1FoundBeforeIndex2 (feature_index,
                               m->modifier_indices->data.intvalue, 
                               available_modifiers_list))
    {
      return FALSE;
    }
    vnp = vnp->next;
  }
  return TRUE;
}

static ValNodePtr GetListOfAvailableModifiers ( ModifierItemLocalPtr ItemList)
{
  ValNodePtr  vnp, head;
  Int2        feature_index;
  BestSortPtr bsp;

  head = NULL;
  for (feature_index = 0; feature_index < numDefLineModifiers; feature_index++)
  {
    if ( ItemList[feature_index].any_present)
    {
      bsp = (BestSortPtr) MemNew (sizeof (BestSortData));
      if (bsp == NULL) return NULL;
      bsp->feature_index = feature_index;
      bsp->all_unique = ItemList[feature_index].all_unique;
      bsp->all_present = ItemList[feature_index].all_present;
      bsp->is_unique = ItemList[feature_index].is_unique;
      vnp = ValNodeNew (head);
      if (vnp == NULL) return NULL;
      vnp->data.ptrvalue = bsp;
      if (head == NULL) head = vnp;
    }
  }
  return head;
}

static Int4 DefLineQualSortOrder [] = {
  DEFLINE_POS_Transgenic,
  DEFLINE_POS_Plasmid_name,
  DEFLINE_POS_Endogenous_virus_name,
  DEFLINE_POS_Strain,
  DEFLINE_POS_Clone,
  DEFLINE_POS_Isolate,
  DEFLINE_POS_Haplotype,
  DEFLINE_POS_Cultivar,
  DEFLINE_POS_Specimen_voucher,
  DEFLINE_POS_Ecotype,
  DEFLINE_POS_Type,
  DEFLINE_POS_Serotype,
  DEFLINE_POS_Authority,
  DEFLINE_POS_Breed
};

static int LIBCALLBACK SortByImportanceAndPresence (
  VoidPtr ptr1,
  VoidPtr ptr2
)
{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  BestSortPtr bsp1, bsp2; 
  Int4       num_defline_qual_sort_order, index;

  if (ptr1 == NULL && ptr2 == NULL) return 0;
  
  if (ptr1 == NULL && ptr2 != NULL) return -1;
  if (ptr1 != NULL && ptr2 == NULL) return 1;
 
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;

  bsp1 = vnp1->data.ptrvalue;
  bsp2 = vnp2->data.ptrvalue;
  if (bsp1->feature_index == bsp2->feature_index) return 0;

  if (bsp1->all_present && bsp1->all_unique
    && (! bsp2->all_present || ! bsp2->all_unique))
  {
    return -1;
  }
  if (bsp2->all_present && bsp2->all_unique
    && (! bsp1->all_present || ! bsp1->all_unique))
  {
    return 1;
  }

  if ( ! bsp1->is_unique && bsp2->is_unique) return -1;
  if ( ! bsp2->is_unique && bsp1->is_unique) return 1;

  num_defline_qual_sort_order = sizeof (DefLineQualSortOrder) / sizeof (Int4);
  for (index = 0; index < num_defline_qual_sort_order; index++)
  {
    if (bsp1->feature_index == DefLineQualSortOrder [ index ]) return -1;
    if (bsp2->feature_index == DefLineQualSortOrder [ index ]) return 1;
  }

  if (bsp1->feature_index > bsp2->feature_index) return 1;
  if (bsp1->feature_index < bsp2->feature_index) return -1;
  return 0;
}


/* The function FindBestCombo tries to find the best combination of modifiers
 * to create unique organism descriptions.  This is accomplished by
 * creating a list of required modifiers, and then creating a list of
 * combinations of modifiers by adding modifiers one at a time
 * to see if the additional modifiers provide any more differentiation in
 * the list.
 * In order to do this, I start with a list of required modifiers, and 
 * then create copies of this list.  For each copy I add one of the modifiers
 * that are present in the bio sources and not already on the list.
 * If adding the modifier increases the differentiation, I add that copy to
 * the list of possible combinations, otherwise I discard it.
 * The function then makes copies of all of the new items added to the list,
 * starting with the item pointed to by start_of_expand, and adds another
 * modifier to each combination, keeping the combinations that increase
 * the differentiation and discarding the rest.
 * This process continues until I have a combination that produces completely
 * differentiated bio sources, or I run out of possible combinations.
 * If the list of possible combinations is exhausted before each organism
 * has a unique description, the function selects the combination from the
 * list with the largest number of unique organism descriptions.  If more 
 * than one combination produces the largest number of unique organisms,
 * the combination with the largest number of unique organisms and the
 * largest number of groups will be selected.
 */
static ModifierCombinationPtr FindBestCombo(
  SeqEntryPtr sep,
  ModifierItemLocalPtr ItemList
)
{
  OrgGroupPtr group_list;
  ModifierCombinationPtr mc_list, start_of_expand, best_found, end_of_list;
  ModifierCombinationPtr next_start_of_expand, m, newm;
  Int4 num_to_expand, next_num_to_expand;
  Int2 i;
  ValNodePtr available_modifier_list, vnp;
  BestSortPtr bsp;

  best_found = NULL;

  /* first, get list of organisms */
  group_list = NULL;
  VisitBioSourcesInSep (sep, &group_list, BuildTaxOrgGroupList);

  /* create combo with just the org groups */
  mc_list = NewModifierCombo (group_list);
  if (mc_list == NULL) return NULL;

  available_modifier_list = GetListOfAvailableModifiers (ItemList);

  /* next, add in any required qualifiers */
  for (vnp = available_modifier_list; vnp != NULL; vnp = vnp->next)
  {
    bsp = vnp->data.ptrvalue;
    if (bsp == NULL) return NULL;
    if (ItemList[bsp->feature_index].required)
    {
      AddQualToModifierCombo (mc_list, ItemList + bsp->feature_index,
                                       bsp->feature_index);
    }
  }
  if (mc_list->max_orgs_in_group == 1)
  {
    /* we're done - they're all unique */
    best_found = mc_list;
    return best_found;
  }

  available_modifier_list = ValNodeSort (available_modifier_list,
                                         SortByImportanceAndPresence);
  start_of_expand = mc_list;
  end_of_list = mc_list;
  num_to_expand = 1;
  while (best_found == NULL && start_of_expand != NULL)
  {
    next_num_to_expand = 0;
    next_start_of_expand = NULL;
    for (i=0; i < num_to_expand && start_of_expand != NULL; i++)
    {
      /* try adding qualifiers */ 
      for (vnp = available_modifier_list;
           vnp != NULL && best_found == NULL;
           vnp = vnp->next)
      {
        bsp = vnp->data.ptrvalue;
        if (bsp == NULL) return NULL;
        if (OkToTryAddingQual (start_of_expand, ItemList,
                               available_modifier_list,
                               bsp->feature_index))
        {
          newm = CopyModifierCombo (start_of_expand);
          AddQualToModifierCombo (newm, ItemList + bsp->feature_index,
                                  bsp->feature_index);
          if (start_of_expand->num_groups >= newm->num_groups)
          {
            /* situation didn't get better, don't bother to add this one */
            FreeModifierCombo (newm);
            newm = NULL;
          }
          else if (newm->max_orgs_in_group == 1)
          {
            best_found = newm;
          }
          else
          {
            end_of_list->next = newm;
            end_of_list = end_of_list->next;
            if (next_start_of_expand == NULL)
              next_start_of_expand = newm;
            next_num_to_expand++;
          }
        }
      }
      if (start_of_expand != NULL)
      {
        start_of_expand = start_of_expand->next;
      }
    }
    num_to_expand = next_num_to_expand;
    if (start_of_expand != NULL)
    {
      start_of_expand = start_of_expand->next;
    }
  }

  if (best_found != NULL)
  {
    FreeModifierCombo (mc_list);
    return best_found;
  }
  
  /* we want to find the one with the highest number of unique organisms */
  best_found = mc_list;
  m = mc_list->next;
  while (m!= NULL)
  {
    if (m->num_unique_orgs > best_found->num_unique_orgs)
    {
      best_found = m;
    }
    else if (m->num_unique_orgs == best_found->num_unique_orgs
           && m->num_groups > best_found->num_groups)
    {
      best_found = m;
    }
    else if (m->num_unique_orgs == best_found->num_unique_orgs
           && m->num_groups == best_found->num_groups
           && m->num_mods < best_found->num_mods)
    {
      best_found = m;
    }
    m = m->next;
  }

  m = mc_list;
  while (m != NULL)
  {
    if (m != best_found)
    {
      newm = m->next;
      m->next = NULL;
      FreeModifierCombo (m);
      m = newm;
    }
    else
    {
      FreeModifierCombo (m->next);
      m->next = NULL;
      m = NULL;
    }
  }
  return best_found;
}


/* create combo with the specified modifiers */
NLM_EXTERN ValNodePtr GetModifierIndicesFromModList (
  ModifierItemLocalPtr modList
)
{
  Int4       feature_index;
  ValNodePtr modifier_indices = NULL;

  if (modList == NULL) return NULL;
  for (feature_index = 0; feature_index < numDefLineModifiers; feature_index++)
  {
    if (modList[feature_index].any_present && modList [feature_index].required)
    {
      ValNodeAddInt (&modifier_indices, 0, feature_index);
    }
  }
  return modifier_indices;
}


/* This is the callback function for sorting the modifier list.  It 
 * implements an order specified by the indexers.
 */
static Int4 DefLineQualPresentationOrder [] = {
  DEFLINE_POS_Transgenic,
  DEFLINE_POS_Strain,
  DEFLINE_POS_Isolate,
  DEFLINE_POS_Cultivar,
  DEFLINE_POS_Specimen_voucher,
  DEFLINE_POS_Ecotype,
  DEFLINE_POS_Type,
  DEFLINE_POS_Serotype,
  DEFLINE_POS_Authority,
  DEFLINE_POS_Breed
};

static int LIBCALLBACK SortByImportance (
  VoidPtr ptr1,
  VoidPtr ptr2
)
{
  ValNodePtr vnp1;
  ValNodePtr vnp2;
  Int4       num_defline_qual_sort_order, index;

  if (ptr1 == NULL && ptr2 == NULL) return 0;
  
  if (ptr1 == NULL && ptr2 != NULL) return -1;
  if (ptr1 != NULL && ptr2 == NULL) return 1;
 
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.intvalue == vnp2->data.intvalue) return 0;

  num_defline_qual_sort_order = sizeof (DefLineQualPresentationOrder) / sizeof (Int4);
  for (index = 0; index < num_defline_qual_sort_order; index++)
  {
    if (vnp1->data.intvalue == DefLineQualPresentationOrder [ index ]) return -1;
    if (vnp2->data.intvalue == DefLineQualPresentationOrder [ index ]) return 1;
  }

  if ((vnp1->data.intvalue < 0 || vnp1->data.intvalue > numDefLineModifiers)
    && (vnp2->data.intvalue < 0 || vnp2->data.intvalue > numDefLineModifiers))
  {
    return 0;
  }
  if (vnp1->data.intvalue < 0 || vnp1->data.intvalue > numDefLineModifiers)
  {
    return 1;
  }
  if (vnp2->data.intvalue < 0 || vnp2->data.intvalue > numDefLineModifiers)
  {
    return -1;
  }

  if (DefLineModifiers [ vnp1->data.intvalue].isOrgMod
    && (! DefLineModifiers [ vnp2->data.intvalue].isOrgMod
      || vnp2->data.intvalue == DEFLINE_POS_Plasmid_name
      || vnp2->data.intvalue == DEFLINE_POS_Endogenous_virus_name))
  {
    return -1;
  }
  if (DefLineModifiers [ vnp2->data.intvalue].isOrgMod
    && (! DefLineModifiers [ vnp1->data.intvalue].isOrgMod
      || vnp1->data.intvalue == DEFLINE_POS_Plasmid_name
      || vnp1->data.intvalue == DEFLINE_POS_Endogenous_virus_name))
  {
    return 1;
  }
  
  if (vnp1->data.intvalue == DEFLINE_POS_Plasmid_name)
  {
    return -1;
  }
  if (vnp2->data.intvalue == DEFLINE_POS_Plasmid_name)
  {
    return 1;
  }
    
  if (vnp1->data.intvalue == DEFLINE_POS_Endogenous_virus_name)
  {
    return -1;
  }
  if (vnp2->data.intvalue == DEFLINE_POS_Endogenous_virus_name)
  {
    return 1;
  }

  if (! DefLineModifiers [ vnp1->data.intvalue].isOrgMod
     && vnp2->data.intvalue == DEFLINE_POS_Clone)
  {
    return 1;
  }
  if (! DefLineModifiers [ vnp2->data.intvalue].isOrgMod
     && vnp1->data.intvalue == DEFLINE_POS_Clone)
  {
    return -1;
  }

  if (! DefLineModifiers [ vnp1->data.intvalue].isOrgMod
     && vnp2->data.intvalue == DEFLINE_POS_Haplotype)
  {
    return 1;
  }
  if (! DefLineModifiers [ vnp2->data.intvalue].isOrgMod
     && vnp1->data.intvalue == DEFLINE_POS_Haplotype)
  {
    return -1;
  }

  if (vnp1->data.intvalue > vnp2->data.intvalue) return 1;
  if (vnp1->data.intvalue < vnp2->data.intvalue) return -1;
  return 0;
}

static Boolean RecordHasModifier (
  BioSourcePtr biop,
  Int4         modifier_index
)
{
  OrgModPtr     mod;
  OrgNamePtr    onp;
  SubSourcePtr  ssp;

  if (biop == NULL
    || modifier_index < 0
    || modifier_index >= numDefLineModifiers)
  {
    return FALSE;
  }
  if (DefLineModifiers[modifier_index].isOrgMod)
  {
    if (biop->org == NULL || (onp = biop->org->orgname) == NULL)
    {
      return FALSE;
    }
    mod = onp->mod;
    while (mod != NULL
        && mod->subtype != DefLineModifiers[modifier_index].subtype)
    {
      mod = mod->next;
    }
    if (mod != NULL && mod->subname != NULL)
    {
      return TRUE;
    }
  } else {
    ssp = biop->subtype;
    while (ssp != NULL && ssp->subtype != DefLineModifiers[modifier_index].subtype)
    {
      ssp = ssp->next;
    }
    if (ssp != NULL && ssp->name != NULL)
    {
      return TRUE;
    }
  }
  return FALSE;
}
  
/* This function adds in required modifiers for HIV sequences */
static void AddHIVModifierIndices (
  ValNodePtr PNTR modifier_indices,
  BioSourcePtr biop,
  ModifierItemLocalPtr modList,
  CharPtr taxName,
  Int4    clone_isolate_HIV_rule_num
)
{
  ValNodePtr  vnp;
  Boolean have_country_in_list;
  Boolean have_isolate_in_list;
  Boolean have_clone_in_list;
  Boolean have_country_mod;
  Boolean have_isolate_mod;
  Boolean have_clone_mod;

  /* special handling for HIV */
  if (StringNICmp (taxName, "HIV-1", 5) != 0
    && StringNICmp (taxName, "HIV-2", 5) != 0)
  {
    return;
  }

  have_country_in_list = FALSE;
  have_isolate_in_list = FALSE;
  have_clone_in_list = FALSE;
  have_country_mod = RecordHasModifier (biop, DEFLINE_POS_Country);
  have_isolate_mod = RecordHasModifier (biop, DEFLINE_POS_Isolate);
  have_clone_mod = RecordHasModifier (biop, DEFLINE_POS_Clone);

  if (modifier_indices != NULL)
  {
    for (vnp = *modifier_indices;
         vnp != NULL
           && (! have_country_in_list
             || ! have_isolate_in_list
             || ! have_clone_in_list);
         vnp = vnp->next)
    {
      if (vnp->data.intvalue == DEFLINE_POS_Country)
      {
        have_country_in_list = TRUE;
      }
      else if (vnp->data.intvalue == DEFLINE_POS_Isolate)
      {
        have_isolate_in_list = TRUE;
      }
      else if (vnp->data.intvalue == DEFLINE_POS_Clone)
      {
        have_clone_in_list = TRUE;
      }
    }
  }

  if ( ! have_country_in_list && have_country_mod)
  {
    vnp = ValNodeNew (*modifier_indices);
    vnp->data.intvalue = DEFLINE_POS_Country;
    if (*modifier_indices == NULL) *modifier_indices = vnp;
  }

  if ((have_clone_in_list && have_clone_mod)
      || (have_isolate_in_list && have_isolate_mod))
  {
    /* don't need HIV rule */
  }
  else
  {
    if ( ! have_isolate_in_list
        && have_isolate_mod
        && ( clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_prefer_isolate
          || clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_want_both
          || ! have_clone_mod))
    {
      vnp = ValNodeNew (*modifier_indices);
      vnp->data.intvalue = DEFLINE_POS_Isolate;
      if (*modifier_indices == NULL) *modifier_indices = vnp;
    }
    
    if ( ! have_clone_in_list
        && have_clone_mod
        && ( clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_prefer_clone
          || clone_isolate_HIV_rule_num == clone_isolate_HIV_rule_want_both
          || ! have_isolate_mod))
    {
      vnp = ValNodeNew (*modifier_indices);
      vnp->data.intvalue = DEFLINE_POS_Clone;
      if (*modifier_indices == NULL) *modifier_indices = vnp;
    }
  }
}

/* This function looks for an OrgMod note that contains the phrase
 * "type strain of".  This function is used to determine whether
 * strain is a required modifier for the defline for this source.
 */
static Boolean HasTypeStrainComment (BioSourcePtr biop)
{
  OrgModPtr mod;
  
  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL)
  {
    return FALSE;
  }
  
  mod = biop->org->orgname->mod;
  while (mod != NULL && mod->subtype != ORGMOD_strain)
  {
    mod = mod->next;
  }
  
  if (mod == NULL)
  {
    return FALSE;
  }
  
  if (!UseOrgModifier (mod, biop->org->taxname, FALSE))
  {
    return FALSE;
  }
  
  mod = biop->org->orgname->mod;
  while (mod != NULL)
  {
    if (mod->subtype == 255
        && StringISearch (mod->subname, "type strain of") != NULL)
    {
      return TRUE;
    }
    mod = mod->next;
  }
  return FALSE;
}


/* This function checks to see if there is a type strain comment on
 * the bio source.  If there is one, it checks to see whether strain
 * is already in the list of modifiers for the definition line.
 * If strain is not already in the list, it is added.
 */
static void 
AddTypeStrainModifierIndices 
(ValNodePtr PNTR modifier_indices,
 BioSourcePtr    biop)
{
  ValNodePtr vnp;
  
  if (modifier_indices == NULL || biop == NULL || ! HasTypeStrainComment (biop))
  {
    return;
  }

  for (vnp = *modifier_indices;
       vnp != NULL && vnp->data.intvalue != DEFLINE_POS_Strain;
       vnp = vnp->next)
  {
  }
  
  if (vnp == NULL)
  {
    ValNodeAddInt (modifier_indices, 0, DEFLINE_POS_Strain);
  }
}

static Boolean SpecialHandlingForSpecialTechniques (
  BioseqPtr bsp
);

/* This function checks to see if the Bioseq has a WGS technique.
 * If so, and if the strain text is not present in the taxname,
 * and strain is not already in the list of modifiers for the 
 * definition line, add strain.
 */
static void 
AddWGSModifierIndices 
(ValNodePtr PNTR modifier_indices,
 BioSourcePtr    biop,
 BioseqPtr       bsp)
{
  ValNodePtr vnp;
  OrgModPtr  omp;
  
  if (modifier_indices == NULL || biop == NULL 
      || biop->org == NULL
      || biop->org->orgname == NULL
      || biop->org->orgname->mod == NULL
      || ! SpecialHandlingForSpecialTechniques (bsp))
  {
    return;
  }

  for (vnp = *modifier_indices;
       vnp != NULL && vnp->data.intvalue != DEFLINE_POS_Strain;
       vnp = vnp->next)
  {
  }
  
  if (vnp == NULL)
  {
    omp = biop->org->orgname->mod;
    while (omp != NULL && omp->subtype != ORGMOD_strain) 
    {
      omp = omp->next;
    }
    if (omp != NULL) 
    {
      if (StringStr (biop->org->taxname, omp->subname) != NULL) 
      {
        /* don't add, present already */
      } else {
        /* add strain modifier */
        ValNodeAddInt (modifier_indices, 0, DEFLINE_POS_Strain);
      }
    }
  }
}

/* This function provides a label to be used in the definition line for
 * each modifier that requires one.  Most modifiers use a label that is 
 * similar to the name of the modifier displayed in the definition line
 * options dialog.
 */
NLM_EXTERN void AddModifierLabel (
  Boolean use_labels,
  Boolean is_orgmod,
  Uint1   subtype,
  CharPtr modifier_text
)
{
  CharPtr cp;
  if (!is_orgmod && subtype == SUBSRC_endogenous_virus_name)
  {
    StringCpy (modifier_text, "endogenous virus");
  }
  else if (is_orgmod && subtype == ORGMOD_specimen_voucher)
  {
    if (use_labels)
    {
      StringCpy (modifier_text, "voucher");
    } 
    else
    {
      modifier_text [0] = 0;
    }
  }
  else if (use_labels 
           || (!is_orgmod 
               && (subtype == SUBSRC_transgenic
                   || subtype == SUBSRC_plasmid_name)))
  {
    if (is_orgmod) 
    {
      StringCpy (modifier_text, GetOrgModQualName (subtype));
    } else {
      StringCpy (modifier_text, GetSubsourceQualName (subtype));
    }
    modifier_text[0] = tolower(modifier_text[0]);
    cp = StringStr (modifier_text, "-name");
    if (cp != NULL) *cp = 0;
  }
  else
  {
    modifier_text[0] = 0;
  }
}

typedef struct orgmodabbrevdata {
  Int2    subtype;
  CharPtr abbrev;
} OrgModAbbrevData, PNTR OrgModAbbrevPtr;

static OrgModAbbrevData orgmod_abbrevs[] = {
  { ORGMOD_variety, "var." },
  { ORGMOD_forma, "f." },
  { ORGMOD_forma_specialis, "f. sp." },
  { ORGMOD_pathovar, "pv." }
};

#define NUM_orgmod_abbrevs sizeof (orgmod_abbrevs) / sizeof (OrgModAbbrevData)



static Boolean FindModifierTextInTaxname (CharPtr search_text, Int2 subtype, CharPtr taxName, Boolean allow_at_end)
{
  CharPtr value_found, abbrev_start;
  Int4    value_len, i;
  Boolean other_abbrev_found;

  value_found = StringStr (taxName, search_text);
  value_len = StringLen (search_text);
  while (value_found != NULL)
  {
    if (value_found == taxName)
    {
      value_found = StringStr (value_found + 1, search_text);
      continue;
    }
    if (*(value_found - 1) != ' ' && *(value_found - 1) != '(')
    {
      value_found = StringStr (value_found + 1, search_text);
      continue;
    }
    if (*(value_found - 1) == ')' && *(value_found + value_len) != ')')
    {
      value_found = StringStr (value_found + 1, search_text);
      continue;
    }
    if (*(value_found + value_len) != ' ' && *(value_found + value_len) != 0)
    {
      value_found = StringStr (value_found + 1, search_text);
      continue;
    }
    if (allow_at_end && value_found != NULL && StringCmp (value_found, search_text) == 0) {
      return FALSE;
    }
    other_abbrev_found = FALSE;
    for (i = 0; i < NUM_orgmod_abbrevs; i++)
    {
      abbrev_start = value_found - StringLen (orgmod_abbrevs[i].abbrev) - 1;
      if (abbrev_start > taxName
        && StringNCmp (abbrev_start,
                        orgmod_abbrevs[i].abbrev,
                        StringLen (orgmod_abbrevs[i].abbrev)) == 0)
      {
        if (subtype == orgmod_abbrevs[i].subtype)
        {
          return TRUE;
        }
        else
        {
          other_abbrev_found = TRUE;
        }
      }
    }
    if ( ! other_abbrev_found 
      && ( subtype == ORGMOD_strain
        || subtype == ORGMOD_sub_species
        || subtype == ORGMOD_specimen_voucher
        || subtype == ORGMOD_isolate
        || subtype == ORGMOD_cultivar))
    {
      return TRUE;
    }
    value_found = StringStr (value_found + 1, search_text);
  }
  return FALSE;
}


/* The UseOrgModifier function looks for the values of certain kinds of 
 * modifiers in the taxonomy name, so that they will not be added to the
 * definition line as modifiers if they are already present in the
 * taxonomy name.
 */
NLM_EXTERN Boolean UseOrgModifier (
  OrgModPtr mod,
  CharPtr   taxName,
  Boolean   allow_at_end
)
{
  Boolean value_found = FALSE;
  CharPtr search_text;
  CharPtr cp;

  if (mod == NULL || mod->subname == NULL) return FALSE;

  /* If selected modifiers already appear in the tax Name, */
  /* don't use them in the organism description again */
  if (mod->subtype == ORGMOD_strain
    || mod->subtype == ORGMOD_variety
    || mod->subtype == ORGMOD_sub_species
    || mod->subtype == ORGMOD_forma
    || mod->subtype == ORGMOD_forma_specialis
    || mod->subtype == ORGMOD_pathovar
    || mod->subtype == ORGMOD_specimen_voucher
    || mod->subtype == ORGMOD_isolate
    || mod->subtype == ORGMOD_cultivar)
  {
    if (FindModifierTextInTaxname (mod->subname, mod->subtype, taxName, allow_at_end)) {
      value_found = TRUE;
    } else if (mod->subtype == ORGMOD_specimen_voucher && (cp = StringChr (mod->subname, ':')) != NULL) {
      search_text = StringSave (mod->subname);
      search_text[cp - mod->subname] = ' ';
      value_found = FindModifierTextInTaxname (search_text, mod->subtype, taxName, allow_at_end);
      search_text = MemFree (search_text);
    }
  }
  return !value_found;
}

/* The SetRequiredModifiers function copies the default required values from
 * the global DefLineModifiers array into the local list of modifier
 * information.
 */
NLM_EXTERN void SetRequiredModifiers (
  ModifierItemLocalPtr modList
)
{
  Int4  item_index;

  for (item_index = 0; item_index < numDefLineModifiers; item_index++)
  {
    modList[item_index].required = IsDeflineModifierRequiredByDefault(DefLineModifiers[item_index].isOrgMod,
                                                                      DefLineModifiers[item_index].subtype);
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

/* This function generates the modifiers for "AutoDefID" */
NLM_EXTERN void SetAutoDefIDModifiers (ModifierItemLocalPtr modList)
{
  Int4 index;
  Boolean added_required = FALSE;

  if (modList == NULL) return;
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
}


/* This function fixes HIV abbreviations, removes items in parentheses,
 * and trims spaces around the taxonomy name.
 */
NLM_EXTERN void CleanUpTaxName (
  CharPtr taxName,
  Boolean keep_in_paren
)
{
  CharPtr ptr;

  if (StringICmp (taxName, "Human immunodeficiency virus type 1") == 0
    || StringICmp (taxName, "Human immunodeficiency virus 1") == 0)
  {
    StringCpy (taxName, "HIV-1");
  }
  else if (StringICmp (taxName, "Human immunodeficiency virus type 2") == 0
    || StringICmp (taxName, "Human immunodeficiency virus 2") == 0)
  {
    StringCpy (taxName, "HIV-2");
  }
  else
  {
    if (! keep_in_paren)
    {
      ptr = StringStr (taxName, "(");
      if (ptr != NULL)
        *ptr = '\0';
    }
    TrimSpacesAroundString (taxName);
  }
}

/* This function gets the BioSource descriptor for the BioSeq. */
NLM_EXTERN BioSourcePtr GetBiopForBsp (
  BioseqPtr bsp
)
{
  SeqMgrDescContext  dcontext;
  SeqDescrPtr    sdp;
  BioSourcePtr    biop;

  if (bsp == NULL) return NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    return biop;
  }
  
  return NULL;
}


NLM_EXTERN Boolean IsSpName (CharPtr taxName)
{
  CharPtr cp;

  cp = StringStr (taxName, " sp.");
  /* check to make sure not "f. sp." */
  if (cp != NULL && cp[4] == ' '
      && (cp - taxName < 2 || *(cp - 2) != 'f' || *(cp - 1) != '.'))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static ValNodePtr ValNodeIntCopy (ValNodePtr orig)
{
  ValNodePtr cpy = NULL, last = NULL, vnp;

  while (orig != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = orig->choice;
    vnp->data.intvalue = orig->data.intvalue;
    if (last == NULL) {
      cpy = vnp;
    } else {
      last->next = vnp;
    }
    last = vnp;
    orig = orig->next;
  }
  return cpy;
}


NLM_EXTERN Boolean IsTSA (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;
  MolInfoPtr        mip;
  Boolean           rval = FALSE;

  if (bsp == NULL) return FALSE;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
       sdp != NULL && !rval;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_molinfo, &context)) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->tech == MI_TECH_tsa) {
      rval = TRUE;
    }
  }
  return rval;
}


NLM_EXTERN Boolean IsGenomeProjectIDDescriptor (SeqDescrPtr sdp) 
{
  UserObjectPtr        uop;
  ObjectIdPtr          oip;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return FALSE;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop != NULL) {
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "GenomeProjectsDB") == 0) {
      return TRUE;
    }
  }
  return FALSE;
}


NLM_EXTERN SeqDescrPtr GetGenomeProjectIDDescriptor (BioseqPtr bsp)
{
  SeqDescrPtr sdp;

  if (bsp == NULL) return NULL;
  sdp = bsp->descr;
  while (sdp != NULL) {
    if (IsGenomeProjectIDDescriptor(sdp)) {
      return sdp;
    }
    sdp = sdp->next;
  }
  return NULL;
}



NLM_EXTERN Int4 GetGenomeProjectID (BioseqPtr bsp)
{
  SeqMgrDescContext context;
  SeqDescrPtr       sdp;
  UserObjectPtr     uop;
  UserFieldPtr      ufp;
  Int4              gpid = 0;

  if (bsp == NULL) return 0;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL && gpid == 0) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL && uop->type != NULL && StringCmp (uop->type->str, "GenomeProjectsDB") == 0)
    {
      ufp = uop->data;
      while (ufp != NULL && gpid == 0) {
        if (ufp->label != NULL 
            && StringCmp (ufp->label->str, "ProjectID") == 0
            && ufp->choice == 2) {
          gpid = ufp->data.intvalue;
        }
        ufp = ufp->next;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }

  return gpid;
}


static void AddSpTaxnameToList (SeqDescrPtr sdp, Pointer userdata)
{
  BioSourcePtr biop;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL || biop->org == NULL || !IsSpName (biop->org->taxname)) return;

  ValNodeAddPointer ((ValNodePtr PNTR) userdata, 0, biop->org->taxname);
}


NLM_EXTERN Boolean ShouldExcludeSp (SeqEntryPtr sep)
{
  ValNodePtr name_list = NULL, vnp1, vnp2;
  Boolean    all_diff = TRUE;

  if (sep == NULL) return TRUE;
  VisitDescriptorsInSep (sep, &name_list, AddSpTaxnameToList);

  name_list = ValNodeSort (name_list, SortVnpByString);

  if (name_list != NULL && name_list->next != NULL)
  {
    for (vnp1 = name_list; vnp1 != NULL && vnp1->next != NULL && all_diff; vnp1 = vnp1->next)
    {
      for (vnp2 = vnp1->next; vnp2 != NULL && all_diff; vnp2 = vnp2->next)
      {
        if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) == 0)
        {
          all_diff = FALSE;
        }
      }
    }
  }
  name_list = ValNodeFree (name_list);
  return all_diff;
}


/* This function sets the default values for the organism description settings */
NLM_EXTERN void InitOrganismDescriptionModifiers(OrganismDescriptionModifiersPtr odmp, SeqEntryPtr sep)
{
  if (odmp == NULL) {
    return;
  }
  MemSet (odmp, 0, sizeof (OrganismDescriptionModifiers));
  odmp->use_labels = TRUE;
  odmp->max_mods = -99;
  odmp->keep_paren = TRUE;
  odmp->exclude_sp = ShouldExcludeSp (sep);
  odmp->exclude_cf = FALSE;
  odmp->exclude_aff = FALSE;
  odmp->exclude_nr = FALSE;
  odmp->include_country_extra = FALSE;
  odmp->clone_isolate_HIV_rule_num = clone_isolate_HIV_rule_want_both;
  odmp->use_modifiers = FALSE;
  odmp->allow_semicolon_in_modifier = FALSE;
  odmp->allow_mod_at_end_of_taxname = FALSE;

}


static ValNodePtr AddPlasmid (SubSourcePtr ssp_list, Boolean use_labels, Boolean allow_semicolon_in_modifier)
{
  Char         modifier_text [256];
  ValNodePtr   strings = NULL;
  Uint4        no_semicolon_len, label_len;

  while (ssp_list != NULL) {
    if (ssp_list->name != NULL && ssp_list->name[0] != 0
        && ssp_list->subtype == SUBSRC_plasmid_name
        && StringCmp (ssp_list->name, "unnamed") != 0)
    {
      AddModifierLabel (use_labels, FALSE, ssp_list->subtype, modifier_text);
      if (modifier_text[0] != 0)
      {
        StringCat (modifier_text, " ");
      }
      label_len = StringLen (modifier_text);
      if (allow_semicolon_in_modifier)
      {
        no_semicolon_len = StringLen (ssp_list->name);
      }
      else
      {
        no_semicolon_len = StringCSpn (ssp_list->name, ";");
      }
      if (no_semicolon_len > sizeof (modifier_text) - 1 - label_len)
      {
        no_semicolon_len = sizeof (modifier_text) - 1 - label_len;
      }
      StringNCat (modifier_text, ssp_list->name, no_semicolon_len);
      modifier_text [ no_semicolon_len + label_len ] = 0;
      ValNodeCopyStr( &strings, 0, modifier_text);
    }
    ssp_list = ssp_list->next;
  }
  return strings;
}


/* This function generates a string describing the organism based on the
 * modifiers selected and other organism description options.
 */
static CharPtr GetOrganismDescription (
  BioseqPtr bsp, 
  ModifierItemLocalPtr modList,
  ValNodePtr   modifier_indices,
  OrganismDescriptionModifiersPtr odmp
)
{
  Char         taxName [196];
  Char         modifier_text [256];
  ValNodePtr   strings = NULL;
  BioSourcePtr biop;
  OrgModPtr    mod;
  SubSourcePtr ssp;
  ValNodePtr   vnp;
  Int2         feature_index;
  CharPtr      org_desc;
  CharPtr      cp;
  Uint4        no_semicolon_len, label_len;
  CharPtr      tmp;
  Char         id[255];
  Int4         gpid;
  SeqIdPtr     sip;
  DbtagPtr     dbtag;

  taxName [0] = '\0';

  biop = GetBiopForBsp (bsp);
  if (biop == NULL) return NULL;
  if (biop->org == NULL) return NULL;
  if (biop->org->taxname == NULL) return NULL;
  StringNCpy (taxName, biop->org->taxname, sizeof (taxName) - 1);
  taxName [ sizeof (taxName) - 1] = 0;

  CleanUpTaxName (taxName, odmp->keep_paren);

  if (biop->origin == ORG_MUT)
  {
    ValNodeAddStr (&strings, 0, StringSave ("Mutant"));
  }

  ValNodeAddStr (&strings, 0, StringSave (taxName));

  if (odmp->exclude_sp && IsSpName(taxName))
  {
    ValNodeLink (&strings, AddPlasmid (biop->subtype, odmp->use_labels, odmp->allow_semicolon_in_modifier));
    org_desc = MergeValNodeStrings (strings, FALSE);
    ValNodeFreeData (strings);
    return org_desc;
  }

  if (odmp->exclude_cf)
  {
    cp = StringStr (taxName, " cf.");
    if (cp != NULL)
    {
      org_desc = MergeValNodeStrings (strings, FALSE);
      ValNodeFreeData (strings);
      return org_desc;
    }
  }
  
  if (odmp->exclude_aff)
  {
    cp = StringStr (taxName, " aff.");
    if (cp != NULL)
    {
      org_desc = MergeValNodeStrings (strings, FALSE);
      ValNodeFreeData (strings);
      return org_desc;
    }
  }
  if (odmp->exclude_nr)
  {
    cp = StringStr (taxName, " nr.");
    if (cp != NULL)
    {
      org_desc = MergeValNodeStrings (strings, FALSE);
      ValNodeFreeData (strings);
      return org_desc;
    }
  }
  

  if (HasTypeStrainComment (biop))
  {
    
  } 

  /* copy modifier indices list */
  modifier_indices = ValNodeIntCopy (modifier_indices);
  AddHIVModifierIndices (&modifier_indices, biop, modList, taxName,
                         odmp->clone_isolate_HIV_rule_num);
  AddTypeStrainModifierIndices (&modifier_indices, biop);
  AddWGSModifierIndices (&modifier_indices, biop, bsp);
  
  modifier_indices = ValNodeSort (modifier_indices, SortByImportance);
  for (vnp = modifier_indices;
       vnp != NULL && (odmp->max_mods == -99 || odmp->max_mods > 0);
       vnp = vnp->next)
  {
    feature_index = vnp->data.intvalue;
    if (! odmp->use_modifiers && !IsDeflineModifierRequiredByDefault(DefLineModifiers[feature_index].isOrgMod,
                                                                     DefLineModifiers[feature_index].subtype))
    {
      /* do nothing */
    }
    else if (DefLineModifiers[feature_index].isOrgMod)
    {
      if (biop->org == NULL || biop->org->orgname == NULL) continue;
      mod = biop->org->orgname->mod;
      while (mod != NULL
        && mod->subtype != DefLineModifiers[feature_index].subtype)
      {
        mod = mod->next;
      }
      if ( UseOrgModifier (mod, taxName, odmp->allow_mod_at_end_of_taxname))
      {
        if (odmp->allow_semicolon_in_modifier) {
          no_semicolon_len = StringLen (mod->subname);
        } else {
          no_semicolon_len = StringCSpn (mod->subname, ";");
        }

        if (mod->subtype == ORGMOD_nat_host)
        {
          sprintf (modifier_text, "from ");
          if (no_semicolon_len > sizeof (modifier_text) - 6)
          {
            no_semicolon_len = sizeof (modifier_text) - 6;
          }
          StringNCpy (modifier_text + 5, mod->subname,
                      no_semicolon_len);
          modifier_text[no_semicolon_len + 5] = 0;
        }
        else
        {
          AddModifierLabel (odmp->use_labels, TRUE, mod->subtype, modifier_text);
          if (modifier_text[0] != 0)
            StringCat (modifier_text, " ");
          label_len = StringLen (modifier_text);
          if (no_semicolon_len > (Int4) sizeof (modifier_text) - label_len - 1)
          {
            no_semicolon_len = (Int4) sizeof (modifier_text) - label_len - 1;
          } 
          if (mod->subtype == ORGMOD_specimen_voucher && StringNICmp (mod->subname, "personal:", 9) == 0)
          {
            tmp = mod->subname + 9;
            while (isspace (*tmp)) 
            {
              tmp++;
            }
            if (odmp->allow_semicolon_in_modifier) {
              no_semicolon_len = StringLen (tmp);
            } else {
              no_semicolon_len = StringCSpn (tmp, ";");
            }
          }
          else
          {
            tmp = mod->subname;
          }

          StringNCat (modifier_text, tmp,
                      no_semicolon_len);
          modifier_text [ no_semicolon_len + label_len] = 0;
        }
        ValNodeCopyStr( &strings, 0, modifier_text);
        if (odmp->max_mods != -99)
          odmp->max_mods --;
      }
    } else {
      ssp = biop->subtype;
      while (ssp != NULL
          && ssp->subtype != DefLineModifiers[feature_index].subtype)
      {
        ssp = ssp->next;
      }
      if (ssp != NULL)
      {
        if (odmp->include_country_extra || odmp->allow_semicolon_in_modifier)
        {
          no_semicolon_len = StringLen (ssp->name);
        }
        else
        {
          no_semicolon_len = StringCSpn (ssp->name, ";");
        }
        AddModifierLabel (odmp->use_labels, FALSE, ssp->subtype, modifier_text);
        if (ssp->subtype == SUBSRC_transgenic)
        {
          /* do nothing, transgenic already captured from label */
        }
        else if (ssp->subtype == SUBSRC_country)
        {
          sprintf (modifier_text, "from ");
          if (no_semicolon_len > sizeof (modifier_text) - 6)
          {
            no_semicolon_len = sizeof (modifier_text) - 6;
          }
          StringNCpy (modifier_text + 5, ssp->name, no_semicolon_len);
          modifier_text[5 + no_semicolon_len] = 0;
          if (!odmp->include_country_extra)
          {
            cp = StringChr (modifier_text, ':');
            if (cp != NULL) *cp = 0;
          }
        }
        else if (ssp->name != NULL && ssp->name[0] != 0
          && (ssp->subtype != SUBSRC_plasmid_name
            || StringCmp (ssp->name, "unnamed") != 0))
        {
          if (modifier_text[0] != 0)
            StringCat (modifier_text, " ");
          label_len = StringLen (modifier_text);
          if (no_semicolon_len > sizeof (modifier_text) - 1 - label_len)
          {
            no_semicolon_len = sizeof (modifier_text) - 1 - label_len;
          }
          StringNCat (modifier_text, ssp->name, no_semicolon_len);
          modifier_text [ no_semicolon_len + label_len ] = 0;
        }

        ValNodeCopyStr( &strings, 0, modifier_text);
        if (odmp->max_mods != -99)
          odmp->max_mods --;
      }
    }
  }

  /* add TSA project ID if necessary */
  if (IsTSA (bsp)) {
    gpid = GetGenomeProjectID (bsp);
    if (gpid > 0) {
      sprintf (id, "gpid:%d", gpid);
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GENERAL && sip->data.ptrvalue != NULL) {
          dbtag = (DbtagPtr) sip->data.ptrvalue;
          if (StringCmp (dbtag->db, id) == 0 && dbtag->tag != NULL) {
            if (dbtag->tag->str != NULL) {
              ValNodeAddPointer (&strings, 0, StringSave (dbtag->tag->str));
            } else {
              sprintf (id, "%d", dbtag->tag->id);
              ValNodeAddPointer (&strings, 0, StringSave (id));
            }
            break;
          }
        }
      }
    }
  }
  
  org_desc = MergeValNodeStrings (strings, FALSE);
  ValNodeFreeData (strings);
  modifier_indices = ValNodeFree (modifier_indices);
  return org_desc;

}

/* end of organism description section */

/* This section of code contains functions which are useful for dealing
 * with locations of features (SeqLocPtr objects).
 */

/* This function determines whether location A is on the same strand as
 * location B
 */
static Boolean AreAAndBOnSameStrand (
  SeqLocPtr slp1,
  SeqLocPtr slp2
)
{
  Uint1 strand1;
  Uint2 strand2;

  strand1 = SeqLocStrand (slp1);
  strand2 = SeqLocStrand (slp2);
  if (strand1 == Seq_strand_minus && strand2 != Seq_strand_minus)
    return FALSE;
  else if (strand1 != Seq_strand_minus && strand2 == Seq_strand_minus)
    return FALSE;
  else
    return TRUE;
}

/* This function determines whether location A is contained in or equal to
 * location B and on the same strand as location B.
 */
NLM_EXTERN Boolean IsLocAInBonSameStrand (
  SeqLocPtr slp1,
  SeqLocPtr slp2
)
{
  if (! AreAAndBOnSameStrand ( slp1, slp2))
  {
    return FALSE;
  }
  else if ( SeqLocAinB (slp1, slp2) < 0)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

/* This function calculates the intersection between two locations.
 */
static SeqLocPtr SeqLocIntersection (
  SeqLocPtr slp1,
  SeqLocPtr slp2,
  BioseqPtr bsp
)
{
  SeqLocPtr diff1, diff2, result;

  diff1 = SeqLocMerge ( bsp, slp1, NULL, FALSE, TRUE, FALSE);
  diff1 = SeqLocSubtract (diff1, slp2);
  diff2 = SeqLocMerge ( bsp, slp2, NULL, FALSE, TRUE, FALSE);
  diff2 = SeqLocSubtract (diff2, slp1);
  result = SeqLocMerge ( bsp, slp1, slp2, FALSE, TRUE, FALSE);
  
  if (diff1 != NULL)
  {
    result = SeqLocSubtract (result, diff1);
    SeqLocFree (diff1);
    if (result == NULL) return NULL;
  }
  if (diff2 != NULL)
  {
    result = SeqLocSubtract (result, diff2);
    SeqLocFree (diff2);
    if (result == NULL) return NULL;
  }
  return result;
}

#define ADJACENT_TYPE_ANY        0
#define ADJACENT_TYPE_UPSTREAM   1
#define ADJACENT_TYPE_DOWNSTREAM 2

/* This function determines whether A is "next to" B and upstream or downstream
 * from B.  A cannot overlap B.  If allow_interval is TRUE, there can be
 * space between A and B.
 */
static Boolean IsAAdjacentToB (
  SeqLocPtr a,
  SeqLocPtr b,
  BioseqPtr bsp,
  Int2      adjacent_type,
  Boolean   allow_interval
)
{
  Int4      a_end, b_end;
  Uint2     strand;
 
  if (adjacent_type != ADJACENT_TYPE_ANY
    && adjacent_type != ADJACENT_TYPE_UPSTREAM
    && adjacent_type != ADJACENT_TYPE_DOWNSTREAM)
  {
    return FALSE;
  }
  
  if ( ! AreAAndBOnSameStrand (a, b))
  {
    return FALSE;
  }

  strand = SeqLocStrand (a);
  if ( adjacent_type == ADJACENT_TYPE_ANY)
  {
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_RIGHT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_LEFT_END);
    if ((allow_interval && b_end < a_end)
      || b_end == a_end + 1)
    {
      return TRUE;
    }
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_LEFT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_RIGHT_END);
    if ((allow_interval && b_end > a_end)
      || a_end == b_end + 1)
    {
      return TRUE;
    }
    return FALSE;
  }
  else if ( (strand == Seq_strand_minus
      && adjacent_type == ADJACENT_TYPE_UPSTREAM)
    || (strand != Seq_strand_minus
      && adjacent_type == ADJACENT_TYPE_DOWNSTREAM))
  {
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_RIGHT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_LEFT_END);
    if ((allow_interval && b_end < a_end)
      || b_end == a_end + 1)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  else
  {
    a_end = GetOffsetInBioseq (a, bsp, SEQLOC_LEFT_END);
    b_end = GetOffsetInBioseq (b, bsp, SEQLOC_RIGHT_END);
    if ((allow_interval && b_end > a_end)
      || a_end == b_end + 1)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
}

static Boolean IsAEmptyIntervalOfB (SeqLocPtr a, SeqLocPtr b, BioseqPtr bsp)
{
  Int4 a_right, a_left, b_right, b_left, prev_right, prev_left;
  SeqLocPtr slp;
  Uint1 a_strand, b_strand;
  
  if (a == NULL || b == NULL || bsp == NULL) return FALSE;
  
  a_strand = SeqLocStrand (a);
  b_strand = SeqLocStrand (b);
  if ((a_strand == Seq_strand_minus && b_strand != Seq_strand_minus)
      || (a_strand != Seq_strand_minus && b_strand == Seq_strand_minus)) {
      return FALSE;
  }
  
  a_right = GetOffsetInBioseq (a, bsp, SEQLOC_RIGHT_END);
  a_left = GetOffsetInBioseq (a, bsp, SEQLOC_LEFT_END);
  
  slp = SeqLocFindNext (b, NULL);
  prev_right = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
  prev_left = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
  slp = SeqLocFindNext (b, slp);
  while (slp != NULL) {
    b_right = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
    b_left = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
    if (a_left == prev_right + 1 && a_right == b_left - 1) {
      return TRUE;
    } else if (a_left == b_right + 1 && a_right == prev_left - 1) {
      return TRUE;
    } else {
      prev_right = b_right;
      prev_left = b_left;
      slp = SeqLocFindNext (b, slp);
    }
  }
  return FALSE;
}


static Boolean LocAContainsIntervalOfB (SeqLocPtr a, SeqLocPtr b)
{
  SeqLocPtr interval;
  Boolean   rval = FALSE;

  if (a == NULL || b == NULL) return FALSE;

  interval = SeqLocFindNext (b, NULL);
  while (interval != NULL && !rval) {
    if (IsLocAInBonSameStrand (interval, a)) {
      rval = TRUE;
    } else {
      interval = SeqLocFindNext (b, interval);
    }
  }
  return rval;
}


/* This section of code deals with identifying and labeling features
 * for the definition line.
 * The features currently handled are:
 *     genes
 *     exons
 *     introns
 *     LTRs
 *     3' UTRs
 *     5' UTRs
 *     CDSs
 *     rRNA
 *     mRNA
 *     misc RNA
 *     snRNA
 *     snoRNA
 *     insertion sequences
 *     integrons
 *     D-loops
 *     mRNA
 *     tRNA
 *     control regions
 *     misc feature listed as intergenic spacer in comment
 *     satellite sequences
 *     promoter regions
 *     endogenous virus source features
 *     transposons
 */

static Boolean LIBCALLBACK IsGene (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return FALSE;
  return TRUE;
}

static CharPtr GetGeneName (GeneRefPtr grp, Boolean suppress_locus_tag)
{
  ValNodePtr syn;

  if (grp == NULL) return NULL;
  if (SeqMgrGeneIsSuppressed (grp)) return NULL;
  if (StringDoesHaveText (grp->locus)) return grp->locus;
  if (! suppress_locus_tag && StringDoesHaveText (grp->locus_tag)) 
      return grp->locus_tag;
  if (StringDoesHaveText (grp->desc)) return grp->desc;
  for (syn = grp->syn; syn != NULL; syn = syn->next)
  {
    if (syn != NULL && syn->data.ptrvalue != NULL)
      return syn->data.ptrvalue;
  }
  return NULL;
}

static CharPtr GetAlleleName (GeneRefPtr grp, Boolean suppress_locus_tag)
{
  size_t  lenallele;
  size_t  lengenename;
  CharPtr  gene_name;
  CharPtr  buffer;

  if (grp == NULL) return NULL;
  if (StringHasNoText (grp->allele)) return NULL;
  gene_name = GetGeneName (grp, suppress_locus_tag);
  if (StringHasNoText (gene_name)) return NULL;
  lenallele = StringLen (grp->allele);
  lengenename = StringLen (gene_name);
  
  if (lenallele > lengenename
    && StringNICmp (gene_name, grp->allele, lengenename) == 0)
  {
    return StringSave (grp->allele);
  }
  else if (grp->allele[0] == '-')
  {
    buffer = MemNew (lenallele + lengenename + 1);
    if (buffer == NULL) return NULL;
    StringCpy (buffer, gene_name);
    StringCat (buffer, grp->allele);
  }
  else
  {
    buffer = MemNew (lenallele + lengenename + 2);
    if (buffer == NULL) return NULL;
    StringCpy (buffer, gene_name);
    StringCat (buffer, "-");
    StringCat (buffer, grp->allele);
  }

  return buffer;
}

/* This function compares the gene names and allele names of the gene
 * to see if they match.
 */
static Boolean DoGenesMatch 
(GeneRefPtr grp1,
 GeneRefPtr grp2,
 Boolean suppress_locus_tag)
{
  CharPtr name1;
  CharPtr name2;

  name1 = GetGeneName (grp1, suppress_locus_tag);
  name2 = GetGeneName (grp2, suppress_locus_tag);
  if (StringCmp (name1, name2) != 0) return FALSE;
  
  name1 = GetAlleleName (grp1, suppress_locus_tag);
  name2 = GetAlleleName (grp2, suppress_locus_tag);
  if ((name1 == NULL && name2 != NULL)
    || (name1 != NULL && name2 == NULL))
  {
    if (name1 != NULL) MemFree (name1);
    if (name2 != NULL) MemFree (name2);
    return FALSE;
  }

  if ((name1 == NULL && name2 == NULL)
           || (StringCmp (name1, name2) == 0))
  {
    if (name1 != NULL) MemFree (name1);
    if (name2 != NULL) MemFree (name2);
    return TRUE;
  }
  
  if (name1 != NULL) MemFree (name1);
  if (name2 != NULL) MemFree (name2);
  return  FALSE;
}

/* This function looks at the pseudo flag on the object itself as well as
 * the pseudo flag on the gene reference for the object (if one is present).
 */
NLM_EXTERN Boolean IsPseudo (
  SeqFeatPtr sfp
)
{
  GeneRefPtr grp;
  SeqMgrFeatContext context;

  if (sfp == NULL) return FALSE;
  if (sfp->pseudo) return TRUE;
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    grp = sfp->data.value.ptrvalue;
  }
  else
  { 
    grp = SeqMgrGetGeneXref (sfp);
  }
  if (grp == NULL) 
  {
    if (sfp->data.choice != SEQFEAT_GENE) {
      sfp = SeqMgrGetOverlappingGene(sfp->location, &context);
      return IsPseudo(sfp);
    } else {
      return FALSE;
    }
  } else {
    return grp->pseudo;
  }
}

static Boolean LIBCALLBACK IsExon (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_exon) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsIntron (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_intron) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsExonOrIntron (SeqFeatPtr sfp)
{
  return IsExon(sfp) || IsIntron(sfp);
}

static Boolean LIBCALLBACK IsLTR (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_LTR) return FALSE;
  return TRUE;
}

static CharPtr GetLTRDescription (
  SeqFeatPtr sfp
)
{
  CharPtr description;
  size_t comment_len;
  if (sfp == NULL) return NULL;
  if (sfp->comment == NULL) return NULL;
  comment_len = StringLen (sfp->comment);
  if (comment_len > 3 && StringCmp (sfp->comment + comment_len - 3, "LTR") == 0)
  {
    description = (CharPtr) MemNew (comment_len - 3);
    if (description == NULL) return NULL;
    StringNCpy (description, sfp->comment, comment_len - 4);
    description[comment_len - 4] = 0;
  }
  else
  {
    description = StringSave (sfp->comment);
  }
  return description;
}

static Boolean LIBCALLBACK Is3UTR (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_3UTR) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK Is5UTR (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_5UTR) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsCDS (SeqFeatPtr sfp)
{
  if (sfp == NULL) return FALSE;
  if (sfp->data.choice == SEQFEAT_CDREGION)
    return TRUE;
  return FALSE;
}


static Boolean LIBCALLBACK IsuORF (SeqFeatPtr sfp)
{
  SeqMgrFeatContext context;
  CharPtr           cp;
  Int4              len;

  if (sfp == NULL) return FALSE;
  if (sfp->data.choice == SEQFEAT_CDREGION 
      && sfp->product != NULL
      && SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &context) == sfp) {
    cp = StringStr (context.label, "uORF");
    if (cp != NULL && (cp == context.label || *(cp - 1) == ' ')
        && (*(cp + 4) == 0 || isspace (*(cp + 4)) || isdigit (*(cp + 4)))) {
      return TRUE;
    }
    if ((len = StringLen (context.label)) >= 14
        && StringCmp (context.label + len - 14, "leader peptide") == 0) {
      return TRUE;
    }
  }
  return FALSE;
}


static Boolean LIBCALLBACK IsrRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_rRNA) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsMiscRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL 
    || (sfp->idx.subtype != FEATDEF_misc_RNA
      && sfp->idx.subtype != FEATDEF_otherRNA)) 
  {
    return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK IsncRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return FALSE;
  if (sfp->idx.subtype == FEATDEF_scRNA 
      || sfp->idx.subtype == FEATDEF_snRNA 
      || sfp->idx.subtype == FEATDEF_snoRNA
      || sfp->idx.subtype == FEATDEF_ncRNA
      || sfp->idx.subtype == FEATDEF_tmRNA
      || sfp->idx.subtype == FEATDEF_misc_RNA
      || sfp->idx.subtype == FEATDEF_otherRNA)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static CharPtr GetncRNAProduct (SeqFeatPtr sfp, Boolean use_ncrna_note)
{
  GBQualPtr gbq = NULL;
  CharPtr product = NULL, q_class = NULL, q_product = NULL;
  CharPtr tmp_class = NULL, cp;
  RnaRefPtr rrp;
  RNAGenPtr rgp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA)
  {
    return NULL;
  }
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;

  if (rrp->ext.choice == 3) {
    rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
    if (rgp != NULL) {
      q_class = rgp->_class;
      q_product = rgp->product;
    }
  } else {  
    gbq = sfp->qual;
    while (gbq != NULL && q_product == NULL) {
      if (StringICmp (gbq->qual, "ncRNA_class") == 0
          && !StringHasNoText (gbq->val)) {
        q_class = gbq->val;
      } else if (StringICmp (gbq->qual, "product") == 0
          && !StringHasNoText (gbq->val)) {
        q_product = gbq->val;
      }
      gbq = gbq->next;
    }
  }
  
  if (q_class != NULL) {
    tmp_class = StringSave (q_class);
    cp = tmp_class;
    while (*cp != 0) {
      if (*cp == '_') {
        *cp = ' ';
      }
      cp++;
    }
  }
  if (q_product != NULL) {
    if (tmp_class == NULL 
        || StringStr (q_product, tmp_class) != NULL 
        || StringCmp (tmp_class, "other") == 0
        || StringCmp (tmp_class, "RNase P RNA") == 0) {
      product = StringSave (q_product);
    } else {
      product = (CharPtr) MemNew (sizeof (Char) * (StringLen (q_product) + StringLen (tmp_class) + 2));
      sprintf (product, "%s %s", q_product, tmp_class);
    }
  } else if (q_class != NULL) {
    if (use_ncrna_note && !StringHasNoText (sfp->comment)) {
      product = StringSave (sfp->comment);
    } else if (StringCmp (tmp_class, "other") == 0) {
      product = StringSave ("non-coding RNA");
    } else {
      product = StringSave (tmp_class);
    }
  } else if ((use_ncrna_note || 
              (rrp != NULL && rrp->type == 10))
             && !StringHasNoText (sfp->comment)) {
    product = StringSave (sfp->comment);
  } else {
    product = StringSave ("non-coding RNA");
  }
  tmp_class = MemFree (tmp_class);
  cp = StringChr (product, ';');
  if (cp != NULL) {
    *cp = 0;
  }
  return product;
}


static Boolean LIBCALLBACK IsPrecursorRNA (SeqFeatPtr sfp)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_preRNA) return FALSE;
  return TRUE;
}


static CharPtr mobile_element_keywords [] = {
  "insertion sequence",
  "retrotransposon",
  "non-LTR retrotransposon",
  "transposon",
  "integron",
  "other",
  "SINE",
  "MITE",
  "LINE"
};

enum mobile_element_keyword_nums 
{
  eMobileElementInsertionSequence = 0,
  eMobileElementRetrotransposon,
  eMobileElementNonLTRRetrotransposon,
  eMobileElementTransposon,
  eMobileElementIntegron,
  eMobileElementOther,
  eMobileElementSINE,
  eMobileElementMITE,
  eMobileElementLINE
};

static Int4 StartsWithMobileElementKeyword (CharPtr txt)
{
  Int4 i, keyword_len;
    
  for (i=0; i < sizeof (mobile_element_keywords) / sizeof (CharPtr); i++) {
    keyword_len = StringLen (mobile_element_keywords[i]);
    if (StringNCmp (txt, mobile_element_keywords[i], keyword_len) == 0
        && (*(txt + keyword_len) == ':' || *(txt + keyword_len) == 0)) {      
      return i;
    }
  }
  return -1;
}

static Int4 IsMobileElementGBQual (GBQualPtr gbqual)
{
  Int4 keyword_idx;
  if (gbqual == NULL || gbqual->qual == NULL || gbqual->val == NULL) return -1;
  if (StringCmp (gbqual->qual, "mobile_element") != 0 && StringCmp (gbqual->qual, "mobile_element_type") != 0) return -1;
  keyword_idx = StartsWithMobileElementKeyword (gbqual->val);
  if (keyword_idx < 0) return -1;
  if (keyword_idx == eMobileElementOther
      && StringStr (gbqual->val, "transposable element") == NULL
      && StringStr (gbqual->val, "P element") == NULL) {
    return -1;
  } else {
    return keyword_idx;
  }
}


static Boolean FeatureDoesNotGetPartialComplete (SeqFeatPtr sfp)
{
  GBQualPtr gbqual;
  Int4 keyword_idx;
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    keyword_idx = IsMobileElementGBQual(gbqual);
    if (keyword_idx == eMobileElementSINE
        || keyword_idx == eMobileElementLINE) {
      return TRUE;
    }
  }
  return FALSE;
}


NLM_EXTERN Boolean LIBCALLBACK IsMobileElement (SeqFeatPtr sfp)
{
  GBQualPtr gbqual;
  if (sfp == NULL || (sfp->idx.subtype != FEATDEF_repeat_region && sfp->idx.subtype != FEATDEF_mobile_element)) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    if (IsMobileElementGBQual(gbqual) > -1) {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean LIBCALLBACK IsRemovableMobileElement (SeqFeatPtr sfp)
{
  GBQualPtr gbqual;
  Int4 keyword_idx;
  if (sfp == NULL || (sfp->idx.subtype != FEATDEF_repeat_region && sfp->idx.subtype != FEATDEF_mobile_element)) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    keyword_idx = IsMobileElementGBQual(gbqual);
    if (keyword_idx >= eMobileElementSINE
        && keyword_idx <= eMobileElementLINE) {
      return TRUE;
    }
  }
  return FALSE;
}

static CharPtr GetMobileElementTypeword (CharPtr desc_start, Int4 keyword_idx)
{
  if (keyword_idx < 0) return NULL;
  if (StringHasNoText (desc_start)) {
    return mobile_element_keywords[keyword_idx];
  }
  switch (keyword_idx) {
    case eMobileElementTransposon:
      if (StringStr (desc_start, "P-element") != NULL) {
        return "P-element";
      } else if (StringStr (desc_start, "MITE") != NULL) {
        return "MITE";
      } else {
        return mobile_element_keywords[keyword_idx];
      }
      break;
    case eMobileElementOther:
      return "transposable element";
      break;
    case eMobileElementIntegron:
      if (StringStr (desc_start, "superintegron") != NULL) {
        return "superintegron";
      } else {
        return mobile_element_keywords[keyword_idx];
      }
      break;
    default:
      return mobile_element_keywords[keyword_idx];
      break;
  }      
}


static void LIBCALLBACK GetMobileElementFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  GBQualPtr  gbqual;
  Int4    keyword_idx = -1;
  Int4    keyword_len;
  Int4    val_len;
  SeqFeatPtr      sfp;
  CharPtr         desc_start = NULL, typeword, cp;

  flp->pluralizable = TRUE;
  flp->is_typeword_first = FALSE;
  flp->typeword = NULL;
  flp->description = NULL;

  if (featlist == NULL) return;
  sfp = featlist->data.ptrvalue;
  if (sfp == NULL) return;

  gbqual = sfp->qual;
  while (gbqual != NULL 
         && (keyword_idx = IsMobileElementGBQual(gbqual)) < 0)
  {
    gbqual = gbqual->next;
  }
  if (gbqual == NULL) return;
    
  keyword_len = StringLen (mobile_element_keywords[keyword_idx]);
  desc_start = gbqual->val + keyword_len;
  while (isspace (*desc_start) || *desc_start == ':') {
    desc_start++;
  }

  /* find alternate typewords */
  typeword = GetMobileElementTypeword(desc_start, keyword_idx);
  if (typeword == NULL) return;
  keyword_len = StringLen (typeword);

  flp->typeword = StringSave (typeword);
  val_len = StringLen (desc_start);
  
  if (StringHasNoText (desc_start))
  {
    flp->is_typeword_first = FALSE;
    flp->description = NULL;
  } else if (StringCmp (desc_start, typeword) == 0) {
    /* just the keyword */
    flp->is_typeword_first = FALSE;
    flp->description = NULL;
    return;
  } else if (StringNCmp (desc_start, typeword, keyword_len) == 0) {
    /* starts with keyword */
    /* if keyword is hyphenated portion of name, no pluralization */
    if (desc_start[keyword_len] == '-') {
      flp->description = StringSave (desc_start);
      flp->typeword = MemFree (flp->typeword);
      flp->typeword = StringSave ("");
      flp->pluralizable = FALSE;
    } else {
      flp->is_typeword_first = TRUE;
      flp->description = StringSave (desc_start + keyword_len + 1);
    }
    return;
  } else if (val_len > 8 && StringCmp (desc_start + val_len - keyword_len, typeword) == 0
             && val_len - keyword_len - 1 >= 0
             && isspace (*(desc_start + val_len - keyword_len - 1))) {
    /* ends with keyword */
    flp->is_typeword_first = FALSE;
    flp->description = MemNew (val_len - keyword_len);
    if (flp->description == NULL) return;
    StringNCpy (flp->description, desc_start, val_len - keyword_len - 1);
    flp->description[val_len - keyword_len -1] = 0;
  } else if ((cp = StringStr (desc_start, typeword)) != NULL
             && cp != desc_start
             && isspace (*(cp -1))) {
    /* keyword in the middle */
    flp->description = StringSave (desc_start);
    flp->typeword = MemFree (flp->typeword);
    flp->typeword = StringSave ("");
    flp->pluralizable = FALSE;
  } else {
    /* keyword not in description */
    if (StringICmp (flp->typeword, "integron") == 0) {
      flp->is_typeword_first = FALSE;
    } else {
      flp->is_typeword_first = TRUE;
    }
    flp->description = StringSave (desc_start);
    if (StringCmp (flp->description, "") == 0) {
      flp->is_typeword_first = FALSE;
    }
  }
  if (StringCmp (flp->description, "unnamed") == 0) {
    flp->description = MemFree (flp->description);
  }
}


static Boolean LIBCALLBACK IsEndogenousVirusSequence (
  SeqFeatPtr sfp
)
{
  GBQualPtr gbqual;
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_repeat_region) return FALSE;
  
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    if (StringCmp (gbqual->qual, "endogenous_virus") == 0)
      return TRUE;
  }
  return FALSE;
}

static CharPtr GetEndogenousVirusSequenceDescription (
  SeqFeatPtr sfp
)
{
  GBQualPtr gbqual;

  if (sfp == NULL) return NULL;
  
  gbqual = sfp->qual;
  while (gbqual != NULL && StringCmp (gbqual->qual, "endogenous_virus") != 0)
  {
    gbqual = gbqual->next;
  }
  if (gbqual != NULL)
  {
    if (StringDoesHaveText (gbqual->val)
      && StringCmp (gbqual->val, "unnamed") != 0)
    {
      return StringSave (gbqual->val);
    }
  }
  return NULL;
}

static Boolean LIBCALLBACK IsDloop (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_D_loop) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsmRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mRNA) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IstRNA (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_tRNA) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsControlRegion (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || StringNCmp (sfp->comment, "control region", StringLen ("control region")) != 0)
  {
    return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK IsGeneCluster (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || (StringStr (sfp->comment, "gene cluster") == NULL
        && StringStr (sfp->comment, "gene locus") == NULL))
  {
    return FALSE;
  }
  return TRUE;
}


static void LIBCALLBACK GetGeneClusterFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  CharPtr    cp;
  Int4       datalen;

  if (featlist == NULL || featlist->data.ptrvalue == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (StringHasNoText (main_feat->comment)) return;
  cp = StringStr (main_feat->comment, "gene cluster");
  if (cp == NULL)
  {
    cp = StringStr (main_feat->comment, "gene locus");
    if (cp == NULL) return;
    flp->typeword = StringSave ("gene locus");
  }
  else
  {
    flp->typeword = StringSave ("gene cluster");
  }
  flp->pluralizable = FALSE;
  flp->is_typeword_first = FALSE;
  datalen = cp - main_feat->comment;
  if (datalen > 0)
  {
    flp->description = (CharPtr) MemNew ((datalen + 1) * sizeof (Char));
    StringNCpy (flp->description, main_feat->comment, datalen);
    flp->description [datalen] = 0;
    TrimSpacesAroundString (flp->description);
  }
  else
  {
    flp->description = NULL;
  }
}


static Boolean LIBCALLBACK IsIntergenicSpacer (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || StringStr (sfp->comment, "intergenic spacer") == NULL)
  {
    return FALSE;
  }
  return TRUE;
}

static ValNodePtr GettRNAGenesAndSpacers (CharPtr str);
static ValNodePtr FreetRNAGenesAndSpacersList (ValNodePtr list);

static Boolean LIBCALLBACK IsParsableList (
  SeqFeatPtr sfp
)
{
  ValNodePtr list;

  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL)
  {
    return FALSE;
  }
  
  list = GettRNAGenesAndSpacers (sfp->comment);
  if (list == NULL) 
  {
    return FALSE;
  } 
  else
  {
    FreetRNAGenesAndSpacersList (list);
    return TRUE;
  }
}


/* This function produces the default definition line label for a misc_feature 
 * that has the word "intergenic spacer" in the comment.  If the comment starts
 * with the word "contains", "contains" is ignored.  If "intergenic spacer"
 * appears first in the comment (or first after the word "contains", the text
 * after the words "intergenic spacer" but before the first semicolon (if any)
 * appear after the words "intergenic spacer" in the definition line.  If there
 * are words after "contains" or at the beginning of the comment before the words
 * "intergenic spacer", this text will appear in the definition line before the words
 * "intergenic spacer".
 */
static void LIBCALLBACK GetIntergenicSpacerFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  CharPtr    cp, buffer;
  Int4       datalen, offset;

  if (featlist == NULL || featlist->data.ptrvalue == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (StringHasNoText (main_feat->comment)) return;
  if (StringNCmp (main_feat->comment, "contains ", 9) == 0)
  {
    buffer = main_feat->comment + 9;
  }
  else if (StringNCmp (main_feat->comment, "may contain ", 12) == 0)
  {
    buffer = main_feat->comment + 12;
  }
  else
  {
    buffer = main_feat->comment;
  }
  cp = StringStr (buffer, "intergenic spacer");
  if (cp == NULL) return;
  flp->typeword = StringSave ("intergenic spacer");
  flp->pluralizable = FALSE;
  if (cp == buffer)
  {
    flp->is_typeword_first = TRUE;
    offset = StringLen ("intergenic spacer") + 1;
    if (StringNCmp (cp + offset, "and ", 4) == 0
        || *(cp + StringLen("intergenic spacer")) == 0)
    {
      flp->description = NULL;
    }
    else
    {
      flp->description = StringSave (cp + StringLen ("intergenic spacer") + 1);
      cp = StringChr (flp->description, ';');
      if (cp != NULL)
      {
        *cp = 0;
      }
    }
  }
  else
  {
    flp->is_typeword_first = FALSE;
    datalen = cp - buffer;
    flp->description = MemNew ( datalen + 1);
    if (flp->description == NULL) return;
    StringNCpy (flp->description, buffer, datalen);
    flp->description [datalen] = 0;
    TrimSpacesAroundString (flp->description);
  }
}

/* These structures are used for parsing tRNA and intergenic spacer information
 * from misc_feature comments. 
 */
typedef struct commentfeat
{
  CharPtr product_name;
  CharPtr gene_name;
} CommentFeatData, PNTR CommentFeatPtr;


static CommentFeatPtr CommentFeatFree (CommentFeatPtr cfp)
{
  if (cfp != NULL) {
    cfp->product_name = MemFree (cfp->product_name);
    cfp->gene_name = MemFree (cfp->gene_name);
    cfp = MemFree (cfp);
  }
  return cfp;
}


typedef struct intergenicspacerdef
{
  CharPtr first_gene;
  CharPtr second_gene;
} IntergenicSpacerDefData, PNTR IntergenicSpacerDefPtr;

static IntergenicSpacerDefPtr IntergenicSpacerDefFree (IntergenicSpacerDefPtr ip)
{
  if (ip != NULL) {
    ip->first_gene = MemFree (ip->first_gene);
    ip->second_gene = MemFree (ip->second_gene);
    ip = MemFree (ip);
  }
  return ip;
}


CharPtr s_tRNAGeneFromProduct (CharPtr product)
{
    CharPtr gene = NULL;

    if (StringNCmp (product, "tRNA-", 5) != 0) {
      return NULL;
    }
    product += 5;

    if (StringICmp (product, "Ala") == 0) {
        gene = "trnA";
    } else if (StringICmp (product, "Asx") == 0) {
        gene = "trnB";
    } else if (StringICmp (product, "Cys") == 0) {
        gene = "trnC";
    } else if (StringICmp (product, "Asp") == 0) {
        gene = "trnD";
    } else if (StringICmp (product, "Glu") == 0) {
        gene = "trnE";
    } else if (StringICmp (product, "Phe") == 0) {
        gene = "trnF";
    } else if (StringICmp (product, "Gly") == 0) {
        gene = "trnG";
    } else if (StringICmp (product, "His") == 0) {
        gene = "trnH";
    } else if (StringICmp (product, "Ile") == 0) {
        gene = "trnI";
    } else if (StringICmp (product, "Xle") == 0) {
        gene = "trnJ";
    } else if (StringICmp (product, "Lys") == 0) {
        gene = "trnK";
    } else if (StringICmp (product, "Leu") == 0) {
        gene = "trnL";
    } else if (StringICmp (product, "Met") == 0) {
        gene = "trnM";
    } else if (StringICmp (product, "Asn") == 0) {
        gene = "trnN";
    } else if (StringICmp (product, "Pyl") == 0) {
        gene = "trnO";
    } else if (StringICmp (product, "Pro") == 0) {
        gene = "trnP";
    } else if (StringICmp (product, "Gln") == 0) {
        gene = "trnQ";
    } else if (StringICmp (product, "Arg") == 0) {
        gene = "trnR";
    } else if (StringICmp (product, "Ser") == 0) {
        gene = "trnS";
    } else if (StringICmp (product, "Thr") == 0) {
        gene = "trnT";
    } else if (StringICmp (product, "Sec") == 0) {
        gene = "trnU";
    } else if (StringICmp (product, "Val") == 0) {
        gene = "trnV";
    } else if (StringICmp (product, "Trp") == 0) {
        gene = "trnW";
    } else if (StringICmp (product, "OTHER") == 0) {
        gene = "trnX";
    } else if (StringICmp (product, "Tyr") == 0) {
        gene = "trnY";
    } else if (StringICmp (product, "Glx") == 0) {
        gene = "trnZ";
    }
    return gene;
}


static CommentFeatPtr ParseGeneFromNoteForDefLine (CharPtr PNTR comment)
{
  CommentFeatPtr tdp;
  CharPtr    product_start, product_end, gene_start, gene_end, cp;
  Int4       product_len, gene_len;
  
  if (comment == NULL || *comment == NULL)
  {
    return NULL;
  }
  /* spacers are not genes */
  if (StringNICmp (*comment, "intergenic", 10) == 0 || StringNICmp (*comment, "spacer", 6) == 0) 
  {
    return NULL;
  }
  
  /* tRNA name must start with "tRNA-" and be followed by one uppercase letter and
   * two lowercase letters.
   */
  product_start = *comment;
  gene_start = product_start;
  while (*gene_start != 0 && !isspace (*gene_start) && *gene_start != ',') {
    gene_start++;
  }
  if (gene_start == product_start) {
    return NULL;
  }
  product_end = gene_start;

  /* if tRNA, don't require gene name, but parse if present */
  while (isspace (*gene_start)) {
    gene_start++;
  }
  if (*gene_start == '(') {
    /* parse in gene name */
    gene_start++;
    gene_end = gene_start;
    while (*gene_end != 0 && *gene_end != ')') {
      gene_end++;
    }
    if (*gene_end == 0) {
      return NULL;
    }
    cp = gene_end + 1;
    while (*cp != 0 && isspace (*cp)) {
      cp++;
    }
  } else if (StringNICmp (gene_start, "intergenic", 10) == 0 || StringNICmp (gene_start, "spacer", 6) == 0) {
    /* spacers are not genes */
    return NULL;
  } else if (StringNCmp (product_start, "tRNA", 4) != 0) {
    /* only tRNAs can leave out gene names */
    return NULL;
  } else {
    cp = gene_start;
    gene_start = NULL;

  }

  /* skip over gene or genes if present */
  if (StringNCmp (cp, "genes", 5) == 0) {
    cp +=5;
    while (*cp != 0 && isspace (*cp)) {
      cp++;
    }
  } else if (StringNCmp (cp, "gene", 4) == 0) {
    cp += 4;
    while (*cp != 0 && isspace (*cp)) {
      cp++;
    }
  }

  tdp = (CommentFeatPtr) MemNew (sizeof (CommentFeatData));
  if (tdp == NULL)
  {
    return NULL;
  }
  product_len = product_end - product_start;
  tdp->product_name = (CharPtr) MemNew (sizeof (Char) * (1 + product_len));
  StringNCpy (tdp->product_name, product_start, product_len);
  tdp->product_name[product_len] = 0;
  
  if (gene_start != NULL) {
    gene_len = gene_end - gene_start;
    tdp->gene_name = (CharPtr) MemNew (sizeof (Char) * (1 + gene_len));
    StringNCpy (tdp->gene_name, gene_start, gene_len);
    tdp->gene_name[gene_len] = 0;
  }
  
  *comment = cp;
  return tdp;
}


static IntergenicSpacerDefPtr ParseIntergenicSpacerFromNoteForDef (CharPtr PNTR comment)
{
  IntergenicSpacerDefPtr idp;
  CharPtr                first_gene_start, dash, second_gene_start, second_gene_end, cp;
  Int4                   first_gene_len, second_gene_len;
  
  if (comment == NULL || *comment == NULL)
  {
    return NULL;
  }
  
  /* description must start with "trn" and be followed by one uppercase letter, followed
   * by a dash, followed by "trn", followed by one uppercase letter, followed by whitespace,
   * followed by the phrase "intergenic spacer".
   */
  first_gene_start = *comment;
  dash = first_gene_start;
  while (*dash != 0 && isalpha (*dash)) {
    dash++;
  }
  if (*dash != '-') {
    return NULL;
  }
  second_gene_start = dash + 1;
  second_gene_end = second_gene_start;
  while (*second_gene_end != 0 && isalpha (*second_gene_end)) {
    second_gene_end ++;
  }
  if (!isspace (*second_gene_end)) {
    return NULL;
  }
  cp = second_gene_end;
  while (isspace (*cp)) {
    cp++;
  }
  if (StringNCmp (cp, "intergenic spacer", 17) != 0) {
    return NULL;
  }
  
  idp = (IntergenicSpacerDefPtr) MemNew (sizeof (IntergenicSpacerDefData));
  if (idp == NULL)
  {
    return NULL;
  }
  
  first_gene_len = dash - first_gene_start;
  idp->first_gene = (CharPtr) MemNew (sizeof (Char) * (1 + first_gene_len));
  StringNCpy (idp->first_gene, first_gene_start, first_gene_len);
  idp->first_gene [first_gene_len] = 0;

  second_gene_len = second_gene_end - second_gene_start;
  idp->second_gene = (CharPtr) MemNew (sizeof (Char) * (1 + second_gene_len));
  StringNCpy (idp->second_gene, second_gene_start, second_gene_len);
  idp->second_gene [second_gene_len] = 0;

  *comment = cp + 17;
  return idp;
}

/* This creates a feature clause from a tRNADef structure. */
static FeatureClausePtr 
FeatureClauseFromParsedComment 
(CommentFeatPtr tdp,
 SeqFeatPtr misc_feat,
 Boolean    is_partial,
 BioseqPtr  bsp,
 DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  CharPtr gene_fmt = "%s (%s)";
  
  if (tdp == NULL)
  {
    return NULL;
  }
  
  fcp = NewFeatureClause ( misc_feat, bsp, rp);
  if (fcp != NULL)
  {
    fcp->feature_label_data.is_typeword_first = FALSE;
    fcp->feature_label_data.typeword = StringSave ("gene");
    if (tdp->gene_name == NULL) {
      fcp->feature_label_data.description = StringSave (tdp->product_name);
    } else {
      fcp->feature_label_data.description = (CharPtr) MemNew (sizeof (Char) * (StringLen (gene_fmt)
                                                                               + StringLen (tdp->gene_name)
                                                                               + StringLen (tdp->product_name)));
      if (fcp->feature_label_data.description != NULL)
      {
        sprintf (fcp->feature_label_data.description, gene_fmt,
                tdp->product_name, tdp->gene_name);
      }
    }
    if (is_partial)
    {
      fcp->interval = StringSave ("partial sequence");
    }
    else
    {
      fcp->interval = StringSave ("complete sequence");
    }
  }
  return fcp;
}


static CharPtr AdvancePastSeparators (CharPtr cp)
{
  if (cp == NULL || *cp == '0') return cp;

  if (*cp == ',')
  {
    cp++;
  }
  while (isspace (*cp))
  {
    cp++;
  }
  if (StringNCmp (cp, "and", 3) == 0)
  {
    cp += 3;
  }
  while (isspace (*cp))
  {
    cp++;
  }
  return cp;
}

#define MISCFEAT_TRNA_GENE 1
#define MISCFEAT_TRNA_SPACER 2

static ValNodePtr FreetRNAGenesAndSpacersList (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    if (list->choice == MISCFEAT_TRNA_GENE) {
      list->data.ptrvalue = CommentFeatFree (list->data.ptrvalue);
    } else if (list->choice == MISCFEAT_TRNA_SPACER) {
      list->data.ptrvalue = IntergenicSpacerDefFree (list->data.ptrvalue);
    }
    list->next = NULL;
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static ValNodePtr GettRNAGenesAndSpacers (CharPtr str)
{
  ValNodePtr list = NULL;
  CharPtr    cp;
  CommentFeatPtr gene, last_gene = NULL;
  IntergenicSpacerDefPtr spacer, last_spacer = NULL;
  Boolean                none_left = FALSE, names_correct = TRUE, alternating = TRUE;
  Int4                   last_item_type = 0;

  if (StringNICmp (str, "contains ", 9) == 0) {
    cp = str + 9;
  } else if (StringNICmp (str, "may contain ", 12) == 0) {
    cp = str + 12;
  } else {
    return NULL;
  }

  while (isspace (*cp))
  {
    cp ++;
  }
  
  while (!none_left && *cp != 0 && *cp != ';' && alternating && names_correct) {
    none_left = TRUE;
    gene = ParseGeneFromNoteForDefLine (&cp);
    if (gene != NULL) {
      /* if previous item was spacer, spacer names and gene names must agree */
      if (last_item_type == MISCFEAT_TRNA_SPACER && last_spacer != NULL) {
        if (gene->gene_name == NULL) {
          if (StringCmp (last_spacer->second_gene, s_tRNAGeneFromProduct (gene->product_name)) != 0) {
            names_correct = FALSE;
          }
        } else if (StringCmp (last_spacer->second_gene, gene->gene_name) != 0) {
          names_correct = FALSE;
        }
      }
      ValNodeAddPointer (&list, MISCFEAT_TRNA_GENE, gene);
      cp = AdvancePastSeparators (cp);
      none_left = FALSE;
      last_item_type = MISCFEAT_TRNA_GENE;
      last_gene = gene;
    }
 
    spacer = ParseIntergenicSpacerFromNoteForDef (&cp);
    if (spacer != NULL) {
      /* must alternate between genes and spacers */
      if (last_item_type == MISCFEAT_TRNA_SPACER) {
        alternating = FALSE;
      }
      /* spacer names and gene names must agree */
      if (last_gene != NULL) {
        if (last_gene->gene_name == NULL) {
          if (StringCmp (s_tRNAGeneFromProduct (last_gene->product_name), spacer->first_gene) != 0) {
            names_correct = FALSE;
          }
        } else if (StringCmp (last_gene->gene_name, spacer->first_gene) != 0) {
          names_correct = FALSE;
        }
      }
      ValNodeAddPointer (&list, MISCFEAT_TRNA_SPACER, spacer);
      cp = AdvancePastSeparators (cp);
      none_left = FALSE;
      last_item_type = MISCFEAT_TRNA_SPACER;
      last_spacer = spacer;
    }      
  }
  if ((*cp != 0 && *cp != ';') || !alternating || !names_correct) {
    list = FreetRNAGenesAndSpacersList (list);
  }
  return list;
}


/* This function produces a feature clause list that should replace the original
 * single clause for a misc_feat that contains a note with one or more tRNAs and
 * an intergenic spacer.
 */
static ValNodePtr 
ParsetRNAIntergenicSpacerElements 
(SeqFeatPtr misc_feat,
BioseqPtr   bsp,
DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  ValNodePtr head = NULL, vnp, list;
  Boolean partial5, partial3;
  IntergenicSpacerDefPtr spacer = NULL;
  Boolean                current_is_partial;

  if (misc_feat == NULL 
      || StringHasNoText (misc_feat->comment)) 
  {
    return NULL;
  }

  list = GettRNAGenesAndSpacers (misc_feat->comment);
  if (list != NULL) {
    if (StringNICmp (misc_feat->comment, "may contain ", 12) == 0) {
      fcp = NewFeatureClause (misc_feat, bsp, rp);
      fcp->feature_label_data.description = StringSave (misc_feat->comment + 12);
      fcp->interval = StringSave ("region");
      ValNodeAddPointer (&head, DEFLINE_CLAUSEPLUS, fcp);
    } else {
      CheckSeqLocForPartial (misc_feat->location, &partial5, &partial3);
      for (vnp = list; vnp != NULL; vnp = vnp->next) {
        current_is_partial = (partial5 && vnp == list) || (partial3 && vnp->next == NULL);
        if (vnp->data.ptrvalue == NULL) continue;
        if (vnp->choice == MISCFEAT_TRNA_GENE) {
          fcp = FeatureClauseFromParsedComment (vnp->data.ptrvalue, misc_feat, current_is_partial, bsp, rp);
          if (fcp != NULL) {
            ValNodeAddPointer (&head, DEFLINE_CLAUSEPLUS, fcp);
          }
        } else if (vnp->choice == MISCFEAT_TRNA_SPACER) {
          spacer = (IntergenicSpacerDefPtr) vnp->data.ptrvalue;
          fcp = NewFeatureClause ( misc_feat, bsp, rp);
          if (fcp != NULL)
          {
            fcp->feature_label_data.is_typeword_first = FALSE;
            fcp->feature_label_data.typeword = StringSave ("intergenic spacer");
            fcp->feature_label_data.description = (CharPtr) MemNew (10 * sizeof (Char));
            if (fcp->feature_label_data.description != NULL)
            {
              sprintf (fcp->feature_label_data.description, "%s-%s",
                      spacer->first_gene, spacer->second_gene);
            }
            if (current_is_partial)
            {
              fcp->interval = StringSave ("partial sequence");
            }
            else
            {
              fcp->interval = StringSave ("complete sequence");
            }
            ValNodeAddPointer (&head, DEFLINE_CLAUSEPLUS, fcp);
          }
        }
      }
    }
    list = FreetRNAGenesAndSpacersList (list);
  }
  return head;
}

static Boolean LIBCALLBACK IsSatelliteSequence (
  SeqFeatPtr sfp
)
{
  GBQualPtr gbqual;
  if (sfp == NULL
    || sfp->idx.subtype != FEATDEF_repeat_region)
  {
    return FALSE;
  }
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next)
  {
    if (StringCmp (gbqual->qual, "satellite") == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean LIBCALLBACK IsPromoter (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_promoter) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK IsEndogenousVirusSourceFeature (
  SeqFeatPtr sfp
)
{
  BioSourcePtr biop;
  SubSourcePtr  ssp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_BIOSRC) return FALSE;
  if ((biop = sfp->data.value.ptrvalue) == NULL) return FALSE;
  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != SUBSRC_endogenous_virus_name)
  {
    ssp = ssp->next;
  }
  if (ssp != NULL) return TRUE;
  return FALSE;
}

static CharPtr GetEndogenousVirusSourceFeatureDescription (
  SeqFeatPtr sfp
)
{
  BioSourcePtr biop;
  SubSourcePtr  ssp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_BIOSRC) return NULL;
  if ((biop = sfp->data.value.ptrvalue) == NULL) return NULL;
  ssp = biop->subtype;
  while (ssp != NULL && ssp->subtype != SUBSRC_endogenous_virus_name)
  {
    ssp = ssp->next;
  }
  if (ssp != NULL && ssp->name != NULL)
  {
    return StringSave (ssp->name);
  }
  return NULL;
}


static CharPtr noncoding_feature_keywords[] = {
  "similar to ",
  "contains "
};

static CharPtr find_noncoding_feature_keyword (
  CharPtr comment
)
{
  Int4 i, num_noncoding_feature_keywords, keywordlen;
  CharPtr cp, buffer;

  if (comment == NULL) return NULL;
  num_noncoding_feature_keywords = sizeof (noncoding_feature_keywords) / sizeof (CharPtr);
  for (i=0; i < num_noncoding_feature_keywords; i++)
  {
    keywordlen = StringLen (noncoding_feature_keywords [i]);
    buffer = comment;
    while ((cp = StringStr (buffer, noncoding_feature_keywords [i])) != NULL)
    {
      if ( StringNCmp (cp + keywordlen,
                       "GenBank Accession Number",
                       StringLen ("GenBank Accession Number")) != 0)
      {
        return cp + keywordlen;
      }
      else
      {
        buffer = cp + 1;
      }
    }
  }
  return NULL;
}

static Boolean LIBCALLBACK IsNoncodingProductFeat (
  SeqFeatPtr sfp
)
{
  if ( sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL
    || StringStr (sfp->comment, "intergenic") != NULL
    || IsParsableList (sfp)
    || (find_noncoding_feature_keyword (sfp->comment) == NULL
      && (StringStr (sfp->comment, "nonfunctional ") == NULL
        || StringStr (sfp->comment, " due to ") == NULL)))
  {
    return FALSE;
  }

  return TRUE;
}

static CharPtr GetNoncodingProductFeatProduct (
  SeqFeatPtr sfp
)
{
  CharPtr productname;
  Int4    namelen, compare_len;
  CharPtr name_start, sep;

  if (sfp == NULL || sfp->comment == NULL) return NULL;

  if ((name_start = StringStr (sfp->comment, "nonfunctional ")) != NULL
    && (sep = StringStr (sfp->comment, " due to ")) != NULL
    && sep > name_start)
  {
    productname = StringSave (name_start);
    productname [ sep - name_start] = 0;
    return productname;
  }

  name_start = find_noncoding_feature_keyword (sfp->comment);
  if (name_start == NULL) return NULL;

  sep = StringStr (name_start, ";");
  if (sep == NULL)
  {
    namelen = StringLen (name_start);
  }
  else
  {
    namelen = sep - name_start;
  }

  productname = MemNew (namelen + 6);
  if (productname == NULL) return NULL;

  StringNCpy (productname, name_start, namelen);
  productname [namelen] = 0;
  
  /* remove sequence from end of name if present */
  compare_len = StringLen (" sequence");
  if (StringCmp (productname + namelen - compare_len, " sequence") == 0)
  {
    productname [ namelen - compare_len] = 0;
    namelen = StringLen (productname);
  }
  /* add "-like" if not present */
  compare_len = StringLen ("-like");
  if (StringCmp (productname + namelen - compare_len, "-like") != 0)
  {
    StringCat (productname, "-like");
    namelen = StringLen (productname);
  }
  return productname;
}

static Boolean LIBCALLBACK IsMiscFeat (
  SeqFeatPtr sfp
)
{
  if ( sfp == NULL
    || sfp->idx.subtype != FEATDEF_misc_feature
    || sfp->comment == NULL)
  {
    return FALSE;
  }

  return TRUE;
}

static Boolean LIBCALLBACK IsOperon (
  SeqFeatPtr sfp
)
{
  if (sfp == NULL 
    || sfp->idx.subtype != FEATDEF_operon)
  {
    return FALSE;
  }

  return TRUE;
}

static Boolean IsRecognizedFeature (
  SeqFeatPtr sfp
)
{
  if (IsGene (sfp)
    || IsCDS (sfp)
    || IsExon (sfp)
    || IsIntron (sfp)
    || IsLTR (sfp)
    || IsrRNA (sfp)
    || IstRNA (sfp)
    || IsmRNA (sfp)
    || IsMiscRNA (sfp)
    || IsncRNA (sfp)
    || IsPrecursorRNA (sfp)
    || Is3UTR (sfp)
    || Is5UTR (sfp)
    || IsMobileElement (sfp)
    || IsEndogenousVirusSequence (sfp)
    || IsEndogenousVirusSourceFeature (sfp)
    || IsDloop (sfp)
    || IsSatelliteSequence (sfp)
    || IsControlRegion (sfp)
    || IsIntergenicSpacer (sfp)
    || IsGeneCluster (sfp)
    || IsNoncodingProductFeat (sfp)
    || IsPromoter (sfp)
    || IsMiscFeat (sfp)
    || IsOperon (sfp))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/* The following section of code contains functions for dealing with lists of
 * clauses.
 */

/* The functions for freeing the memory associated with lists of clauses 
 * are recursive.
 */
static void FreeListElement (ValNodePtr element);

/* This function simply frees the ValNodePtr, since there is no extra
 * memory associated with a DEFLINE_FEATLIST item - the sfp that is
 * pointed to by data.ptrvalue came from the sequence indexing functions
 * and should under no circumstances be freed.
 */
static void FreeFeatlist (
  ValNodePtr featlist
)
{

  if (featlist == NULL) return;
  ValNodeFree (featlist);
}
 
/* This function frees the memory associated with a FeatureClause, including
 * the memory associated with any subclauses.
 */
static void FreeClausePlusData (
  FeatureClausePtr fcp
)
{
  if (fcp->interval != NULL)
  {
    MemFree (fcp->interval);
    fcp->interval = NULL;
  }
  if (fcp->allelename != NULL)
  {
    MemFree (fcp->allelename);
    fcp->allelename = NULL;
  }
  if (fcp->feature_label_data.typeword != NULL)
  {
    MemFree (fcp->feature_label_data.typeword);
    fcp->feature_label_data.typeword = NULL;
  }
  if (fcp->feature_label_data.description != NULL)
  {
    MemFree (fcp->feature_label_data.description);
    fcp->feature_label_data.description = NULL;
  }
  if (fcp->feature_label_data.productname != NULL)
  {
    MemFree (fcp->feature_label_data.productname);
    fcp->feature_label_data.productname = NULL;
  }
  if (fcp->featlist != NULL)
  {
    FreeListElement (fcp->featlist);
    fcp->featlist = NULL;
  }
  if (fcp->slp != NULL)
  {
    SeqLocFree (fcp->slp);
  }
}
    
/* This function frees the data associated with the FeatureClause
 * and then frees the ValNode.
 */
static void FreeClausePlus (
  ValNodePtr clauseplus
)
{
  FeatureClausePtr data_struct;

  if (clauseplus == NULL) return;
  data_struct = (FeatureClausePtr) clauseplus->data.ptrvalue;
  if (data_struct != NULL)
  {
    FreeClausePlusData (data_struct); 
    MemFree (data_struct);
    clauseplus->data.ptrvalue = NULL;
  }
  ValNodeFree (clauseplus);
}
 
/* This function frees a list of DEFLINE_FEATLIST, DEFLINE_REMOVEFEAT,
 * and DEFLINE_CLAUSEPLUS items, starting with the last item in the list.
 * It recursively frees memory associated with subclauses.
 */ 
static void FreeListElement (
  ValNodePtr element
)
{
  if (element == NULL) return;

  FreeListElement (element->next);
  element->next = NULL;
  if (element->choice == DEFLINE_FEATLIST 
    || element->choice == DEFLINE_REMOVEFEAT)
  {
    FreeFeatlist (element);
  }
  else if (element->choice == DEFLINE_CLAUSEPLUS)
  {
    FreeClausePlus (element);
  }
}

/* This function excises from the list pointed to by head all of the clauses
 * with the delete_me flag set to TRUE and all of the ValNodes with a choice
 * of DEFLINE_REMOVEFEAT.
 */
static void DeleteFeatureClauses (
  ValNodePtr PNTR head
)
{
  ValNodePtr vnp, prev;
  FeatureClausePtr fcp;
  Boolean          delete_this_one;

  if (head == NULL) return;

  prev = NULL;
  vnp = *head;
  while (vnp != NULL)
  {
    delete_this_one = FALSE;

    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp == NULL || fcp->delete_me || fcp->featlist == NULL) 
      {
        delete_this_one = TRUE;
      }
      else
      {
        DeleteFeatureClauses (&fcp->featlist);
        if (fcp->featlist == NULL) delete_this_one = TRUE;
      }
    }
    else if (vnp->choice == DEFLINE_REMOVEFEAT)
    {
      delete_this_one = TRUE;
    }

    if (delete_this_one)
    {
      if (prev == NULL)
      {
        *head = vnp->next;
        vnp->next = NULL;
        FreeListElement (vnp);
        if (*head == NULL) return;
        vnp = *head;
      }
      else
      {
        prev->next = vnp->next;
        vnp->next = NULL;
        FreeListElement (vnp);
        vnp = prev->next;
      }
    }
    else
    {
      prev = vnp;
      vnp = vnp->next;
    }
  }
}

/* This function counts the number of features in the feature list that
 * satisfy the itemmatch function (or all of them, if itemmatch is NULL).
 * If recurse_past_found_item, the function will not count features in
 * subclauses of features that satisfy itemmatch.
 */
static Int4 CountFeatures (
  ValNodePtr clause_list,
  matchFunction  itemmatch,
  Boolean    recurse_past_found_item
)
{
  ValNodePtr       vnp;
  Int4             num_features;
  FeatureClausePtr fcp;

  num_features = 0;
  for (vnp = clause_list;
       vnp != NULL;
       vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST 
      && (itemmatch == NULL || itemmatch (vnp->data.ptrvalue)))
    {
      num_features++;
      if (! recurse_past_found_item)
      {
        return num_features;
      }
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL)
    {
      num_features += CountFeatures (fcp->featlist,
                                     itemmatch,
                                     recurse_past_found_item);   
    }
  }
  return num_features;
}

/* The following section of code contains functions for grouping features. */

typedef struct matchruledata {
  matchFunction is_item;
  Int4          num_match_rules;
  matchFunction *match_rules;
} MatchRuleData, PNTR MatchRulePtr;

static void InitRuleForTopLevelClauses (MatchRulePtr mrp)
{
  if (mrp == NULL)
  {
    return;
  }
  mrp->num_match_rules = 4;
  mrp->match_rules = MemNew (mrp->num_match_rules
                                    * sizeof (matchFunction));
  if (mrp->match_rules == NULL) return;
  mrp->match_rules[0] = IsMobileElement;
  mrp->match_rules[1] = IsEndogenousVirusSourceFeature;
  mrp->match_rules[2] = IsOperon;
  mrp->match_rules[3] = IsGeneCluster;
}

static void InitRuleForBottomLevelClauses (MatchRulePtr mrp)
{
  if (mrp == NULL)
  {
    return;
  }
  mrp->num_match_rules = 6;
  mrp->match_rules = MemNew (mrp->num_match_rules
                                    * sizeof (matchFunction));
  if (mrp->match_rules == NULL) return;
  mrp->match_rules[0] = IsCDS;
  mrp->match_rules[1] = IsmRNA;
  mrp->match_rules[2] = IsGene;
  mrp->match_rules[3] = IsEndogenousVirusSourceFeature;
  mrp->match_rules[4] = IsOperon;
  mrp->match_rules[5] = IsGeneCluster;
}

/* NumGroupingRules is the number of features for which there is a list of
 * features to group under.
 * When grouping features, each feature in the list is examined sequentially.
 * If there is a rule set that applies to that feature, the entire feature 
 * list is searched for each feature type that this feature might group
 * beneath.  This preserves the biological order that was generated by the
 * original listing of features by sequence indexing.
 */
#define  NumGroupingRules 13
static MatchRulePtr InitializeGroupingRules()
{
  MatchRulePtr grouping_rules;

  grouping_rules = MemNew (NumGroupingRules * sizeof (MatchRuleData));
  if (grouping_rules == NULL) return NULL;

  grouping_rules[0].is_item = IsExon;
  grouping_rules[0].num_match_rules = 8;
  grouping_rules[0].match_rules = MemNew (grouping_rules[0].num_match_rules
                                    * sizeof (matchFunction));
  if (grouping_rules[0].match_rules == NULL) return NULL;
  grouping_rules[0].match_rules[0] = IsCDS;
  grouping_rules[0].match_rules[1] = IsNoncodingProductFeat;
  grouping_rules[0].match_rules[2] = IsDloop;
  grouping_rules[0].match_rules[3] = IsmRNA;
  grouping_rules[0].match_rules[4] = IsGene;
  grouping_rules[0].match_rules[5] = IsEndogenousVirusSourceFeature;
  grouping_rules[0].match_rules[6] = IsOperon;
  grouping_rules[0].match_rules[7] = IsGeneCluster;

  grouping_rules[1].is_item = IsIntron;
  grouping_rules[1].num_match_rules = 8;
  grouping_rules[1].match_rules = MemNew (grouping_rules[1].num_match_rules
                                    * sizeof (matchFunction));
  if (grouping_rules[1].match_rules == NULL) return NULL;
  grouping_rules[1].match_rules[0] = IsCDS;
  grouping_rules[1].match_rules[1] = IsNoncodingProductFeat;
  grouping_rules[1].match_rules[2] = IstRNA;
  grouping_rules[1].match_rules[3] = IsDloop;
  grouping_rules[1].match_rules[4] = IsGene;
  grouping_rules[1].match_rules[5] = IsEndogenousVirusSourceFeature;
  grouping_rules[1].match_rules[6] = IsOperon;
  grouping_rules[1].match_rules[7] = IsGeneCluster;
  
  grouping_rules[2].is_item = IsPromoter;
  InitRuleForBottomLevelClauses (grouping_rules + 2);
 
  grouping_rules[3].is_item = IsCDS;
  grouping_rules[3].num_match_rules = 5;
  grouping_rules[3].match_rules = MemNew (grouping_rules[3].num_match_rules
                                    * sizeof (matchFunction));
  if (grouping_rules[3].match_rules == NULL) return NULL;
  grouping_rules[3].match_rules[0] = IsmRNA;
  grouping_rules[3].match_rules[1] = IsMobileElement;
  grouping_rules[3].match_rules[2] = IsEndogenousVirusSourceFeature;
  grouping_rules[3].match_rules[3] = IsOperon;
  grouping_rules[3].match_rules[4] = IsGeneCluster;
  
  grouping_rules[4].is_item = IsMobileElement;
  InitRuleForTopLevelClauses (grouping_rules + 4);
 
  grouping_rules[5].is_item = Is3UTR;
  InitRuleForBottomLevelClauses (grouping_rules + 5);
 
  grouping_rules[6].is_item = Is5UTR;
  InitRuleForBottomLevelClauses (grouping_rules + 6);
 
  grouping_rules[7].is_item = IsLTR;
  InitRuleForBottomLevelClauses (grouping_rules + 7);
 
  grouping_rules[8].is_item = IsGene;
  InitRuleForTopLevelClauses (grouping_rules + 8);
  
  grouping_rules[9].is_item = IsIntergenicSpacer;
  InitRuleForTopLevelClauses (grouping_rules + 9);

  grouping_rules[10].is_item = IsNoncodingProductFeat;
  InitRuleForTopLevelClauses (grouping_rules + 10);
 
  grouping_rules[11].is_item = IsOperon;
  InitRuleForTopLevelClauses (grouping_rules + 11);
 
  grouping_rules[12].is_item = IsGeneCluster;
  InitRuleForTopLevelClauses (grouping_rules + 12);

  return grouping_rules; 
}

static void FreeGroupingRules(
  MatchRulePtr grouping_rules
)
{
  Int4 i;

  if (grouping_rules == NULL) return;

  for (i = 0; i < NumGroupingRules; i++)
  {
    if (grouping_rules[i].match_rules != NULL)
    MemFree (grouping_rules[i].match_rules);
    grouping_rules[i].match_rules = NULL;
  }
  MemFree (grouping_rules);
}

static Boolean IsmRNASequence (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  MolInfoPtr  mip;
  
  if (bsp == NULL || bsp->mol != Seq_mol_rna || bsp->descr == NULL)
  {
    return FALSE;
  }
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_molinfo)
  {
    sdp = sdp->next;
  }
  if (sdp == NULL || sdp->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  
  if (mip->biomol == 3)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

typedef struct matchcandidate 
{
  ValNodePtr matched_clause;
  SeqLocPtr  slp;
} MatchCandidateData, PNTR MatchCandidatePtr;

/* This function searches the search_list for features that satisfy the
 * match function and satisfy locational requirements relative to the
 * clause.
 * If more than one clause meets the match requirements, the smallest one
 * is chosen.
 */
static void FindBestMatchCandidate 
(FeatureClausePtr  clause,
 ValNodePtr        search_list,
 FeatureClausePtr  search_parent,
 matchFunction     match,
 Boolean           gene_cluster_opp_strand,
 BioseqPtr         bsp,
 MatchCandidatePtr current_candidate)
{
  ValNodePtr       search_clause;
  SeqFeatPtr       addsfp, clause_sfp;
  FeatureClausePtr searchfcp;
  SeqLocPtr        slp;

  if (clause == NULL || clause->slp == NULL || current_candidate == NULL) return;
  
  clause_sfp = (SeqFeatPtr) (clause->featlist->data.ptrvalue);

  for (search_clause = search_list;
       search_clause != NULL;
       search_clause = search_clause->next)
  {
    if (search_clause->data.ptrvalue == clause)
      continue;
    if (search_clause->choice == DEFLINE_FEATLIST
      && search_clause->data.ptrvalue != NULL)
    {
      addsfp = search_clause->data.ptrvalue;
      /* slp is the location of the feature we are trying to
       * group this feature with
       */ 
      if (search_parent != NULL)
      {
        slp = search_parent->slp;
      }
      else
      {
        slp = addsfp->location;
      }
      if (match (search_clause->data.ptrvalue))
      { 
        /* Transposons, insertion sequences, integrons, and endogenous virii
         * take subfeatures regardless of whether the subfeature is
         * on the same strand.
         * Gene Clusters can optionally take subfeatures on either
         * strand (gene_cluster_opp_strand is flag).
         * Promoters will match up to features that are adjacent.
         * Introns will match up to coding regions whose intervals
         * are adjacent to the endpoints of the intron, or to other
         * features if the intron location is inside the other feature.
         * All other feature matches must be that the feature to
         * go into the clause must fit inside the location of the
         * other clause.
         */
        if (((match == IsMobileElement
             || match == IsEndogenousVirusSourceFeature
             || (match == IsGeneCluster && gene_cluster_opp_strand))
            && SeqLocAinB (clause->slp, slp) > -1)
          || IsLocAInBonSameStrand (clause->slp, slp)
          || ( IsPromoter (clause_sfp)
            && IsAAdjacentToB (clause->slp, search_parent->slp, bsp,
                               ADJACENT_TYPE_UPSTREAM, TRUE))
          || (IsmRNASequence (bsp) 
              && match != IsMobileElement 
              && match != IsEndogenousVirusSourceFeature
              && match != IsGeneCluster)
          || (match == IsCDS 
              && IsIntron (clause_sfp)
              && IsAEmptyIntervalOfB (clause->slp, slp, bsp))
          || (match == IsCDS
              && IsExon (clause_sfp)
              && LocAContainsIntervalOfB (clause->slp, slp)))
        {
          /* if we don't already have a candidate, or if this
           * candidate's location is inside the current candidate,
           * take this candidate.
           */
          if (current_candidate->matched_clause == NULL
              || SeqLocAinB (slp, current_candidate->slp) > 0)
          {
            current_candidate->matched_clause = search_clause;
            current_candidate->slp = slp;
          }
        }
      }
    }
    else if (search_clause->choice == DEFLINE_CLAUSEPLUS 
      && search_clause->data.ptrvalue != NULL)
    {
      searchfcp = search_clause->data.ptrvalue;
      FindBestMatchCandidate (clause, searchfcp->featlist, searchfcp,
                              match, gene_cluster_opp_strand, bsp,
                              current_candidate);
    }
  }
}


/* This function iterates through the matches in the specified grouping rule.
 * If more than one match is found, the clause with the smallest location is
 * used.
 * If a match is found, the clause is added to the list of clauses for that
 * feature's parent clause.
 */
static Boolean GroupClauseByRule (
  FeatureClausePtr clause,
  ValNodePtr       search_list,
  MatchRulePtr     grouping_rule,
  Boolean          gene_cluster_opp_strand,
  BioseqPtr        bsp
)
{
  Int4               rule_index;
  MatchCandidateData mcd;  
  Boolean            rval = FALSE;
  ValNodePtr         newfeat;
  
  mcd.slp = NULL;
  mcd.matched_clause = NULL;

  for (rule_index = 0;
       rule_index < grouping_rule->num_match_rules;
       rule_index ++)
  {
  
    FindBestMatchCandidate (clause, search_list, NULL, 
                            grouping_rule->match_rules[rule_index],
                            gene_cluster_opp_strand, bsp, &mcd);
  }
  if (mcd.matched_clause != NULL)
  {
    newfeat = ValNodeNew (mcd.matched_clause);
    if (newfeat == NULL) return FALSE;
    newfeat->choice = DEFLINE_CLAUSEPLUS;
    newfeat->data.ptrvalue = clause;
    rval = TRUE;
  }
  return rval;
}


/* This function determines whether a subclause contains just a 3' UTR feature
 * and no other details.
 */
static Boolean Is3UTRClause (FeatureClausePtr clause)
{  
  if (clause == NULL
      || clause->featlist == NULL 
      || clause->featlist->choice != DEFLINE_FEATLIST
      || clause->featlist->data.ptrvalue == NULL
      || clause->featlist->next != NULL)
  {
    return FALSE;
  }
  return Is3UTR(clause->featlist->data.ptrvalue);
}


/* This function will move 3' UTRs to the end of any subfeat lists
 * so that they can be listed after the partial/complete CDS.
 */
static void Move3UTRToEndOfSubFeatList (ValNodePtr clause_list)
{
  ValNodePtr vnp, prev, last_vnp;
  FeatureClausePtr clause;
  
  if (clause_list == NULL || clause_list->next == NULL) 
  {
    return;
  }
  prev = clause_list;
  for (vnp = clause_list->next; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL && vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      clause = vnp->data.ptrvalue;
      if (Is3UTRClause (clause))
      {
        if (vnp->next != NULL)
        {
          /* move to end of clause list */
          last_vnp = vnp->next;
          while (last_vnp->next != NULL)
          {
            last_vnp = last_vnp->next;
          }
          prev->next = vnp->next;
          last_vnp->next = vnp;
          vnp->next = NULL;
        }
      }
      else
      {
        prev = vnp;
        Move3UTRToEndOfSubFeatList (clause->featlist);
      }
    }
    else
    {
      prev = vnp;
    }
  }  
}

/* This function iterates over the list of features, attempting to find and
 * apply grouping rules for each feature.
 */
static void GroupAllClauses (
  ValNodePtr PNTR clause_list,
  Boolean         gene_cluster_opp_strand,
  BioseqPtr       bsp
)
{
  MatchRulePtr     grouping_rules;
  ValNodePtr       vnp, prev;
  FeatureClausePtr clause;
  SeqFeatPtr       main_feat;
  Int4             rule_index;

  grouping_rules = InitializeGroupingRules();
  if (grouping_rules == NULL) return;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      clause = vnp->data.ptrvalue;
      if (clause->featlist != NULL
        && clause->featlist->choice == DEFLINE_FEATLIST
        && clause->featlist->data.ptrvalue != NULL)
      {
        main_feat = clause->featlist->data.ptrvalue;
        for (rule_index = 0;
             rule_index < NumGroupingRules 
               && ! grouping_rules[rule_index].is_item (main_feat);
             rule_index++)
        {
        }
        if (rule_index < NumGroupingRules)
        {
          if ( GroupClauseByRule (clause, *clause_list,
                                  grouping_rules + rule_index, 
                                  gene_cluster_opp_strand,
                                  bsp))
          {
            vnp->data.ptrvalue = NULL;
          }
        }
      }
    }
  }
  FreeGroupingRules(grouping_rules);

  vnp = *clause_list;
  prev = NULL;
  while (vnp != NULL)
  {
    if (vnp->data.ptrvalue == NULL)
    {
      if (prev == NULL)
      {
        *clause_list = vnp->next;
        vnp->next = NULL;
        ValNodeFree (vnp);
        vnp = *clause_list;
      }
      else
      {
        prev->next = vnp->next;
        vnp->next = NULL;
        ValNodeFree (vnp);
        vnp = prev->next;
      }
    }
    else
    {
      prev = vnp;
      vnp = vnp->next;
    }
  }
  
  Move3UTRToEndOfSubFeatList (*clause_list);
}

/* This function exists to handle the special case where two or more exons
 * are alternatively spliced, but there are no CDSs to represent some of the
 * alternatively spliced forms.  In order to make sure that all of the exons
 * that are alternatively spliced together appear with the CDS, they are
 * temporarily consolidated into a single clause with a location that
 * is the intersection of the exons' locations.  The clause will be
 * re-expanded after grouping by the ExpandAltSplicedExons function.
 */
static void GroupAltSplicedExons (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  Boolean         delete_now
)
{
  ValNodePtr       clause, search_clause, vnp;
  FeatureClausePtr fcp, search_fcp;
  SeqFeatPtr       sfp, search_sfp;
  SeqLocPtr        new_slp;

  if (clause_list == NULL) return;
 
  for (clause = *clause_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || clause->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = clause->data.ptrvalue;
    if ( ! fcp->is_alt_spliced
      || fcp->delete_me
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      continue;
    }
    sfp = fcp->featlist->data.ptrvalue;
    if ( ! IsExon (sfp))
    {
      continue;
    }

    for ( search_clause = clause->next; 
          search_clause != NULL
            && search_clause->choice == DEFLINE_CLAUSEPLUS
            && search_clause->data.ptrvalue != NULL
            && (search_fcp = search_clause->data.ptrvalue) != NULL
            && ! search_fcp->delete_me
            && search_fcp->is_alt_spliced
            && search_fcp->featlist != NULL
            && search_fcp->featlist->choice == DEFLINE_FEATLIST
            && (search_sfp = search_fcp->featlist->data.ptrvalue) != NULL
            && IsExon (search_sfp)
            && TestFeatOverlap (sfp, search_sfp, SIMPLE_OVERLAP) != -1;
          search_clause = search_clause->next)
    {
      vnp = ValNodeNew (fcp->featlist);
      if (vnp == NULL) return;
      vnp->choice = DEFLINE_FEATLIST;
      vnp->data.ptrvalue = search_sfp;
      search_fcp->delete_me = TRUE;
      new_slp = SeqLocIntersection (fcp->slp, search_fcp->slp, bsp);
      SeqLocFree (fcp->slp);
      fcp->slp = new_slp;
    }
  }
  if (delete_now)
  {
    DeleteFeatureClauses (clause_list);
  }
}

/* This function expands a clause filled with alternatively-spliced exons
 * that was created in the GroupAltSplicedExons function.
 */
static void ExpandAltSplicedExons (
  ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr clause, rest_of_list, featlist, new_clause = NULL;
  FeatureClausePtr fcp, new_fcp;
  SeqFeatPtr sfp;

  for (clause = clause_list;
       clause != NULL;
       clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || (fcp = clause->data.ptrvalue) == NULL
      || fcp->featlist == NULL)
    {
      continue;
    }
    if ( fcp->featlist->choice == DEFLINE_FEATLIST
      && (sfp = fcp->featlist->data.ptrvalue) != NULL
      && IsExon (sfp)
      && fcp->featlist->next != NULL
      && fcp->featlist->next->choice == DEFLINE_FEATLIST
      && IsExon (fcp->featlist->next->data.ptrvalue))
    {
      rest_of_list = clause->next;
      clause->next = NULL;
      for (featlist = fcp->featlist->next;
           featlist != NULL
             && featlist->choice == DEFLINE_FEATLIST
             && IsExon (featlist->data.ptrvalue);
           featlist = featlist->next)
      {
        new_clause = ValNodeNew (clause);
        if (new_clause == NULL) return;
        new_fcp = NewFeatureClause (featlist->data.ptrvalue, bsp, rp);
        if (new_fcp == NULL) return;
        new_fcp->grp = fcp->grp;
        new_fcp->is_alt_spliced = fcp->is_alt_spliced;
        new_fcp->make_plural = fcp->make_plural;
        new_clause->choice = DEFLINE_CLAUSEPLUS;
        new_clause->data.ptrvalue = new_fcp;
      }
      ValNodeFree (fcp->featlist->next);
      fcp->featlist->next = NULL;
      new_clause->next = rest_of_list;

      /* put back location for first exon - was reduced to union of 
       * all exon intervals in GroupAltSplicedExons
       */
      SeqLocFree (fcp->slp);
      sfp = fcp->featlist->data.ptrvalue;
      fcp->slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
    }
    else
    {
      ExpandAltSplicedExons (fcp->featlist, bsp, rp);
    }
  }
}



static Boolean DoFeaturesShareGene (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  Boolean share_gene = FALSE;
  SeqFeatPtr found_gene1, found_gene2;
  
  if (sfp1 != NULL && sfp2 != NULL
      && !SeqMgrGeneIsSuppressed (SeqMgrGetGeneXref(sfp1))
      && !SeqMgrGeneIsSuppressed (SeqMgrGetGeneXref(sfp2)))
  {
    found_gene1 = SeqMgrGetOverlappingGene (sfp1->location, NULL);
    found_gene2 = SeqMgrGetOverlappingGene (sfp2->location, NULL);
    if (found_gene1 == found_gene2 && found_gene1 != NULL)
    {
      share_gene = TRUE;
    }
  }
  return share_gene;
}

/* This function determines whether two features share the same product name */
static Boolean 
DoProductNamesMatch 
(SeqFeatPtr sfp1,
 SeqFeatPtr sfp2, 
 BioseqPtr  bsp,
 DeflineFeatureRequestListPtr rp)
{
  CharPtr productname1;
  CharPtr productname2;
  Boolean names_match = FALSE;
   
  productname1 = GetProductName (sfp1, bsp, rp);
  productname2 = GetProductName (sfp2, bsp, rp);
  if (StringHasNoText (productname1) && StringHasNoText (productname2))
  {
    names_match = TRUE;
  }
  else if (StringCmp (productname1, productname2) == 0)
  {
    names_match = TRUE;
  }
  
  productname1 = MemFree (productname1);
  productname2 = MemFree (productname2);
  
  return names_match;  
}

/* This function should combine CDSs that do not have a joined location
 * but are part of the same gene and have the same protein name.
 */
static void GroupSegmentedCDSs (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  Boolean         delete_now,
  DeflineFeatureRequestListPtr rp
)
{
  ValNodePtr       clause, search_clause, vnp;
  FeatureClausePtr fcp, search_fcp;
  SeqFeatPtr       sfp, search_sfp;
  SeqLocPtr        new_slp;

  if (clause_list == NULL) return;
 
  for (clause = *clause_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || clause->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = clause->data.ptrvalue;
    if (fcp->delete_me
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      continue;
    }
    sfp = fcp->featlist->data.ptrvalue;
    if ( ! IsCDS (sfp))
    {
      continue;
    }

    for ( search_clause = clause->next; 
          search_clause != NULL;
          search_clause = search_clause->next)
    {
      if (search_clause->choice != DEFLINE_CLAUSEPLUS
          || search_clause->data.ptrvalue == NULL
          || (search_fcp = search_clause->data.ptrvalue) == NULL
          || search_fcp->delete_me
          || search_fcp->featlist == NULL
          || search_fcp->featlist->choice != DEFLINE_FEATLIST
          || (search_sfp = search_fcp->featlist->data.ptrvalue) == NULL
          || ! IsCDS (search_sfp)
          || ! DoFeaturesShareGene (sfp, search_sfp)
          || ! DoProductNamesMatch (sfp, search_sfp, bsp, rp))
      {
        continue;
      }
      vnp = ValNodeNew (fcp->featlist);
      if (vnp == NULL) return;
      vnp->choice = DEFLINE_FEATLIST;
      vnp->data.ptrvalue = search_sfp;
      search_fcp->delete_me = TRUE;
      new_slp = SeqLocMerge (bsp, fcp->slp, search_fcp->slp,
                             FALSE, TRUE, FALSE);

      SeqLocFree (fcp->slp);
      fcp->slp = new_slp;
    }
  }
  if (delete_now)
  {
    DeleteFeatureClauses (clause_list);
  }
}


/* This function searches this list for clauses to which this gene should
 * apply.  This is not taken care of by the GroupAllClauses function
 * because genes are added to clauses as a GeneRefPtr instead of as an
 * additional feature in the list, and because a gene can apply to more
 * than one clause, while other features should really only belong to
 * one clause.
 */
static Boolean AddGeneToClauses 
( SeqFeatPtr gene,
  CharPtr    gene_productname,
  ValNodePtr clause_list,
  Boolean    suppress_locus_tag)
{
  ValNodePtr    clause;
  FeatureClausePtr fcp;
  SeqFeatPtr    sfp, found_gene;
  GeneRefPtr    grp;
  Boolean    used_gene;
  
  if (gene == NULL || gene->data.value.ptrvalue == NULL) return FALSE;
  if (clause_list == NULL) return FALSE;

  used_gene = FALSE;
  grp = gene->data.value.ptrvalue;
  for (clause = clause_list; clause != NULL; clause = clause->next)
  {
    fcp = clause->data.ptrvalue;
    if (fcp == NULL || fcp->featlist == NULL) return FALSE;
    sfp = fcp->featlist->data.ptrvalue;
    if (sfp != NULL && !SeqMgrGeneIsSuppressed (SeqMgrGetGeneXref(sfp))
        && (IsCDS (sfp)
            || IsrRNA (sfp)
            || IstRNA (sfp)
            || IsmRNA (sfp)
            || IsMiscRNA (sfp)
            || IsncRNA (sfp)
            || IsPrecursorRNA (sfp)
            || IsNoncodingProductFeat (sfp)))
    {
      if (fcp->grp == NULL)
      {
        found_gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
        if (found_gene != NULL)
        {
          fcp->grp = (GeneRefPtr) found_gene->data.value.ptrvalue;
        }
      }

      if (fcp->grp != NULL && DoGenesMatch (fcp->grp, grp, suppress_locus_tag))
      {
        used_gene = TRUE;
        if (gene_productname != NULL
          && fcp->feature_label_data.productname == NULL
          && IsCDS (sfp))
        {
          fcp->feature_label_data.productname =
              StringSave (gene_productname);
        }
      }
      else if (fcp->grp == NULL
        && IsLocAInBonSameStrand (sfp->location, gene->location))
      {
        fcp->grp = grp;
        used_gene = TRUE;
        if (gene_productname != NULL
          && fcp->feature_label_data.productname == NULL
          && IsCDS (sfp))
        {
          fcp->feature_label_data.productname =
              StringSave (gene_productname);
        }
      }
    }
  }
  return used_gene;
}

/* This function iterates through the list of features and calls
 * AddGeneToClauses for each gene feature it finds.
 */
static void GroupGenes (ValNodePtr PNTR clause_list, Boolean suppress_locus_tag)
{
  ValNodePtr  vnp;
  ValNodePtr  featlist;
  FeatureClausePtr fcp;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS) return;
    fcp = (FeatureClausePtr) vnp->data.ptrvalue;
    if (fcp == NULL) return;

    featlist = fcp->featlist;
    if (featlist != NULL
      && featlist->choice == DEFLINE_FEATLIST
      && IsGene (featlist->data.ptrvalue))
    {
      AddGeneToClauses (featlist->data.ptrvalue,
                        fcp->feature_label_data.productname,
                        vnp->next, suppress_locus_tag);
    }
  } 
}

/* This function searches this list for clauses to which this mRNA should
 * apply.  This is not taken care of by the GroupAllClauses function
 * because when an mRNA is added to a CDS, the product for the clause is
 * replaced and the location for the clause is expanded, rather than simply
 * adding the mRNA as an additional feature in the list, and because an 
 * mRNA can apply to more than one clause, while other features should 
 * really only belong to one clause.
 */
static Boolean AddmRNAToClauses 
( SeqFeatPtr mRNA,
  ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr    clause;
  FeatureClausePtr fcp;
  SeqFeatPtr    sfp;
  Boolean    used_mRNA;
  CharPtr       productname;
  SeqLocPtr     new_slp;
  
  if (mRNA == NULL || mRNA->data.value.ptrvalue == NULL) return FALSE;
  if (clause_list == NULL) return FALSE;

  used_mRNA = FALSE;
  productname = GetProductName (mRNA, bsp, rp);
  if (productname == NULL) return TRUE;

  for (clause = clause_list; clause != NULL; clause = clause->next)
  {
    fcp = clause->data.ptrvalue;
    if (fcp == NULL || fcp->featlist == NULL) return FALSE;
    sfp = fcp->featlist->data.ptrvalue;
    if (sfp == NULL)
    {
    }
    else if (IsCDS (sfp)
      && fcp->feature_label_data.productname != NULL
      && StringCmp (fcp->feature_label_data.productname, productname) == 0)
    {
      used_mRNA = TRUE;
      fcp->has_mrna = TRUE;
      if (IsLocAInBonSameStrand (sfp->location, mRNA->location))
      {
        new_slp = SeqLocMerge (bsp, fcp->slp, mRNA->location,
                                 FALSE, TRUE, FALSE);
        if (new_slp == NULL) return FALSE;
        if (fcp->slp != NULL)
        {
          SeqLocFree (fcp->slp);
        }
        fcp->slp = new_slp;
      }
    }
    else if (fcp->feature_label_data.productname == NULL
      && (IsCDS (sfp) || IsGene (sfp))
      && (IsLocAInBonSameStrand (sfp->location, mRNA->location)
        || IsLocAInBonSameStrand (mRNA->location, sfp->location)))
    {
      fcp->feature_label_data.productname = StringSave (productname);
      used_mRNA = TRUE;
      fcp->has_mrna = TRUE;
      new_slp = SeqLocMerge (bsp, fcp->slp, mRNA->location,
                                 FALSE, TRUE, FALSE);
      if (new_slp == NULL) return FALSE;
      if (fcp->slp != NULL)
      {
        SeqLocFree (fcp->slp);
      }
      fcp->slp = new_slp;
    }
  }
  return used_mRNA;
}

/* This function iterates through the list of features and calls
 * AddmRNAToClauses for each mRNA feature it finds.
 */
static void GroupmRNAs (
  ValNodePtr PNTR clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp
)
{
  ValNodePtr  vnp;
  ValNodePtr  featlist;
  FeatureClausePtr fcp;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS) return;
    fcp = (FeatureClausePtr) vnp->data.ptrvalue;
    if (fcp == NULL) return;

    featlist = fcp->featlist;
    if (featlist != NULL
      && featlist->choice == DEFLINE_FEATLIST
      && IsmRNA (featlist->data.ptrvalue))
    {
      if (AddmRNAToClauses (featlist->data.ptrvalue, *clause_list, bsp, rp))
      {
        fcp->delete_me = TRUE;
      }
    }
  } 
  DeleteFeatureClauses (clause_list);
}

/* This section of code contains functions for generating labels for
 * clauses for the definition lines.
 */

/* This function examines the specified typeword and determines whether it
 * should appear before or after the description of the feature in the
 * definition line.
 */
static Boolean IsTypeWordFirst (
  CharPtr typeword
)
{
  Int4 i;
  if (typeword == NULL) return FALSE;
  if (StringCmp (typeword, "exon") == 0
    || StringCmp (typeword, "intron") == 0
    || StringCmp (typeword, "endogenous virus") == 0)
  {
    return TRUE;
  }
  else
  {
    i = StartsWithMobileElementKeyword (typeword);
    if (i >= 0 && i != eMobileElementIntegron) {
      return TRUE;
    }
    return FALSE;
  }
}

/* This function determines the word to use to indicate what type of feature
 * is being described in the definition line.  For certain feature types,
 * the word to use in the definition line varies based on the type of
 * molecule in the record.
 */
static CharPtr GetFeatureTypeWord (
  Uint1 biomol,
  SeqFeatPtr sfp
)
{
  if (sfp == NULL) return NULL;
  if ( IsExon (sfp))
  {
    return StringSave ("exon");
  } 
  else if(IsIntron (sfp))
  {
    return StringSave ("intron");
  }
  else if (IsEndogenousVirusSequence (sfp))
  {
    return StringSave ("endogenous virus");
  }
  else if (IsControlRegion (sfp))
  {
    return StringSave ("control region");
  }
  else if (IsEndogenousVirusSourceFeature (sfp))
  {
    return StringSave ("endogenous virus");
  }
  else if (IsDloop (sfp))
  {
    return StringSave ("D-loop");
  }
  else if (IsLTR (sfp))
  {
    return StringSave ("LTR");
  }
  else if (Is3UTR (sfp))
  {
    return StringSave ("3' UTR");
  }
  else if (Is5UTR (sfp))
  {
    return StringSave ("5' UTR");
  }
  else if (IsOperon (sfp))
  {
    return StringSave ("operon");
  }
  else if (biomol == MOLECULE_TYPE_GENOMIC || biomol == MOLECULE_TYPE_CRNA)
  {
    if (IsPseudo (sfp))
    {
      return StringSave ("pseudogene");
    }
    else
    {
      return StringSave ("gene");
    }
  }
  else if ( IsrRNA (sfp) || IsncRNA (sfp))
  {
    return NULL;
  }
  else if (IsPrecursorRNA (sfp))
  {
    return StringSave ("precursor RNA");
  }
  else if (biomol == MOLECULE_TYPE_MRNA)
  {
    if (IsPseudo (sfp))
    {
      return StringSave ("pseudogene mRNA");
    }
    else
    {
      return StringSave ("mRNA");
    }
  }
  else if (biomol == MOLECULE_TYPE_PRE_MRNA)
  {
    if (IsPseudo (sfp))
    {
      return StringSave ("pseudogene precursor RNA");
    }
    else
    {
      return StringSave ("precursor RNA");
    }
  }
  else if (biomol == MOLECULE_TYPE_OTHER_GENETIC_MATERIAL)
  {
    return StringSave ("gene");
  }
  return StringSave ("");
}

/* Frequently the product associated with a feature is listed as part of the
 * description of the feature in the definition line.  This function determines
 * the name of the product associated with this specific feature.  Some
 * features will be listed with the product of a feature that is associated
 * with the feature being described - this function does not look at other
 * features to determine a product name.
 * If the feature is a misc_feat with particular keywords in the comment,
 * the product will be determined based on the contents of the comment.
 * If the feature is a CDS and is marked as pseudo, the product will be
 * determined based on the contents of the comment.
 * If the feature is a gene and has different strings in the description than
 * in the locus or locus tag, the description will be used as the product for
 * the gene.
 * If none of the above conditions apply, the sequence indexing context label
 * will be used to obtain the product name for the feature.
 */
static CharPtr GetProductName 
( SeqFeatPtr cds,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr protein_name;
  CharPtr semicolon;
  size_t len_to_copy;
  SeqMgrFeatContext  context;
  GeneRefPtr  grp;
  CharPtr gene_name;
  RnaRefPtr rrp;
  Boolean suppress_locus_tag = FALSE;

  if (cds == NULL) return NULL;
  protein_name = NULL;
  if (rp != NULL)
  {
    suppress_locus_tag = rp->suppress_locus_tags;
  }
  if (IsNoncodingProductFeat (cds))
  {
    return GetNoncodingProductFeatProduct (cds);
  }
  else if (cds->data.choice == SEQFEAT_CDREGION && cds->pseudo)
  {
    if (cds->comment != NULL)
    {
      semicolon = StringChr (cds->comment, ';');
      if (semicolon != NULL)
      {
        len_to_copy = semicolon - cds->comment;
      }
      else
      {
        len_to_copy = StringLen (cds->comment);
      }
      protein_name = MemNew (len_to_copy + 1);
      if (protein_name == NULL) return NULL;
      StringNCpy (protein_name, cds->comment, len_to_copy);
      protein_name[len_to_copy] = 0;
    }
    return protein_name;
  }
  else if (cds->data.choice == SEQFEAT_GENE)
  {
    grp = (GeneRefPtr) cds->data.value.ptrvalue;
    if (grp == NULL) return NULL;
    gene_name = GetGeneName (grp, suppress_locus_tag);
    if (grp->desc != NULL
      && StringCmp (grp->desc, gene_name) != 0)
    {
      return StringSave (grp->desc);
    }
#if 0
    /* removed by request from Linda Yankie */    
    if (grp->locus_tag != NULL && ! suppress_locus_tag
      && StringCmp (grp->locus_tag, gene_name) != 0)
    {
      return StringSave (grp->locus_tag);
    }
#endif    
  }
  else if (IsncRNA (cds))
  {
    return GetncRNAProduct(cds, rp == NULL ? FALSE : rp->use_ncrna_note);
  }
  else if (IstRNA (cds) 
           && SeqMgrGetDesiredFeature (0, bsp, 0, 0, cds, &context) == cds
           && context.label != NULL)
  {
    if (StringCmp (context.label, "Xxx") == 0) {
      protein_name = StringSave ("tRNA-OTHER");
    } else {
      protein_name = MemNew ( StringLen (context.label) + 6);
      if ( protein_name == NULL) return NULL;
      sprintf (protein_name, "tRNA-%s", context.label);
    }
    return protein_name;
  }
  else if (cds->data.choice == SEQFEAT_RNA)
  {
    
    rrp = (RnaRefPtr) cds->data.value.ptrvalue;
    if (rrp != NULL && rrp->ext.choice == 1 && !StringHasNoText (rrp->ext.value.ptrvalue))
    {
      return StringSave (rrp->ext.value.ptrvalue);
    }
  }
  else if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, cds, &context) == cds
           && context.label != NULL)
  {
    if ((IsCDS(cds) && StringCmp (context.label, "CDS") != 0)
        || (IsmRNA(cds) && StringCmp (context.label, "mRNA") != 0)
        || (! IsCDS(cds) && ! IsmRNA(cds)))
    {
      protein_name = StringSave (context.label);
      return protein_name;
    }
  }
  return NULL;
}

/* This function searches a list of features recursively for a
 * feature that satisfies the itemmatch condition and is associated with
 * the same gene as the fcp clause passed to the function.
 * This is used to obtain a product for a feature that may share a gene with
 * a product-producing feature but may not be contained in the interval of
 * the product-producing feature.
 */
static FeatureClausePtr FindProductInFeatureList (
  FeatureClausePtr fcp,
  ValNodePtr       clause_list,
  matchFunction    itemmatch,
  Boolean          suppress_locus_tag)
{
  ValNodePtr       vnp;
  FeatureClausePtr vnp_fcp;
  
  for (vnp = clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      vnp_fcp = vnp->data.ptrvalue;
      if (DoGenesMatch (vnp_fcp->grp, fcp->grp, suppress_locus_tag)
        && vnp_fcp->featlist != NULL
        && vnp_fcp->featlist->choice == DEFLINE_FEATLIST
        && itemmatch (vnp_fcp->featlist->data.ptrvalue))
      {
        return vnp_fcp;
      }
      else
      {
        vnp_fcp = FindProductInFeatureList (fcp, vnp_fcp->featlist,
                                            itemmatch, suppress_locus_tag);
        if (vnp_fcp != NULL) return vnp_fcp;
      }
    }
  }
  return NULL;
}

/* This function uses the available information in the clause to generate
 * a description from the name of the gene (if any) and the name of the
 * product for the feature (if any).
 * If there is only a gene, the description will be the name of the gene.
 * If there is only a product, the description will be the name of the product.
 * If there is a gene and a product, the description will be the name of
 * the product followed by the name of the gene in parentheses.
 */
static CharPtr GetGeneProtDescription 
( FeatureClausePtr fcp,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr    sfp;
  CharPtr    protein_name;
  CharPtr    gene_name;
  size_t    description_length;
  CharPtr    description;

  if (fcp == NULL
    || fcp->featlist == NULL
    || fcp->featlist->data.ptrvalue == NULL)
  {
    return NULL;
  }
  sfp = fcp->featlist->data.ptrvalue;

  description_length = 0;

  if (fcp->feature_label_data.productname != NULL)
  {
    protein_name = StringSave (fcp->feature_label_data.productname);
  }
  else
  {
    protein_name = GetProductName (sfp, bsp, rp);
    if (protein_name == NULL && IsGene (sfp))
    {
      
    }
  }
  if (protein_name != NULL)
  {
    description_length += StringLen (protein_name);
  }
     
  gene_name = GetGeneName (fcp->grp, rp == NULL ? FALSE : rp->suppress_locus_tags);
  if (gene_name != NULL)
  {
    description_length += StringLen (gene_name);
    if (protein_name != NULL)
    {
      description_length += 3;
    }
  }
  description = (CharPtr) MemNew (description_length + 1);
  if (description == NULL) return NULL;
  if (protein_name != NULL)
  {
    if (gene_name != NULL)
    {
      sprintf (description, "%s (%s)", protein_name, gene_name);
    }
    else
    {
      sprintf (description, "%s", protein_name);
    }
  }
  else
  {
    if (gene_name != NULL)
      sprintf (description, "%s", gene_name);
  }
  if (protein_name != NULL) MemFree (protein_name);
  if (StringHasNoText (description)) {
    description = MemFree (description);
  }
  return description;
}

/* This array of match functions is used to identify, in order of preference,
 * the features that might be used to generate a product for a gene-protein
 * description if the feature has not already been grouped with a product
 * feature.
 */
static matchFunction productfeatures[] = {
  IsCDS, IsmRNA, IstRNA
}; 

/* This function finds gene features without products and looks for 
 * features that might provide products for them.
 */
static void FindGeneProducts 
( ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr       vnp;
  FeatureClausePtr fcp, productfcp;
  Int4             i, NumProductFeatureTypes;
  Boolean          suppress_locus_tag = (rp == NULL ? FALSE : rp->suppress_locus_tags);

  NumProductFeatureTypes = sizeof (productfeatures) / sizeof (matchFunction);

  for (vnp = clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL
      && fcp->featlist != NULL)
    {
      if (fcp->featlist->choice == DEFLINE_FEATLIST
        && IsGene (fcp->featlist->data.ptrvalue)
        && fcp->feature_label_data.productname == NULL)
      {
        productfcp = NULL;
        for (i=0; i < NumProductFeatureTypes && productfcp == NULL; i++)
        {
          productfcp = FindProductInFeatureList (fcp, clause_list,
                                                 productfeatures[i],
                                                 suppress_locus_tag);
        }
        if (productfcp != NULL)
        {
          fcp->is_alt_spliced = productfcp->is_alt_spliced;
          if (productfcp->feature_label_data.productname != NULL)
          {
            fcp->feature_label_data.productname =
                  StringSave (productfcp->feature_label_data.productname);
          }
          else
          {
            fcp->feature_label_data.productname
                  = GetProductName (productfcp->featlist->data.ptrvalue, 
                                    bsp, rp);
          }
          if (fcp->feature_label_data.description != NULL)
          {
            MemFree (fcp->feature_label_data.description);
            fcp->feature_label_data.description = NULL;
          }
          fcp->feature_label_data.description =
            GetGeneProtDescription (fcp, bsp, rp);
        }
      }
      else
      {
        FindGeneProducts (fcp->featlist, bsp, rp);
      }
    }
  }
}

static Boolean ShowInterval (
  SeqFeatPtr sfp
)
{
  if (IsSatelliteSequence (sfp) || IsExon (sfp) || IsIntron (sfp)
    || IsPromoter (sfp) || Is3UTR (sfp) || Is5UTR (sfp))
    return FALSE;
  return TRUE;
}

static CharPtr GetExonDescription (
  BioseqPtr bsp,
  SeqFeatPtr sfp
)
{
  SeqMgrFeatContext  context;
  SeqFeatPtr new_sfp;
  CharPtr    label;

  if ((new_sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, bsp, 0, 0, sfp, &context)) != sfp
      || context.label == NULL)
  {
    if ((new_sfp = SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &context)) != sfp
      || context.label == NULL)
    {
      return NULL;
    }
  }
  if ((IsExon (sfp) && StringCmp (context.label, "exon") == 0)
    || (IsIntron (sfp) && StringCmp (context.label, "intron") == 0))
  {
    return NULL;
  }
  
  label = StringSave (context.label);
  return label;
}

static CharPtr GetFeatureDescription 
( FeatureClausePtr fcp,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr    sfp;

  if ( fcp == NULL
    || fcp->featlist == NULL
    || fcp->featlist->data.ptrvalue == NULL)
  {
    return NULL;
  }
  sfp = fcp->featlist->data.ptrvalue;
  if (sfp == NULL) return NULL;

  if (IsExon (sfp) || IsIntron (sfp))
  {
    return GetExonDescription (bsp, sfp);
  }
  else if (IsEndogenousVirusSequence (sfp))
  {
    return GetEndogenousVirusSequenceDescription (sfp);
  }
  else if (IsEndogenousVirusSourceFeature (sfp))
  {
    return GetEndogenousVirusSourceFeatureDescription (sfp);
  }
  else if (IsControlRegion (sfp))
  {
    return NULL;
  }
  else if (IsDloop (sfp))
  {
    return NULL;
  }
  else if (Is3UTR (sfp))
  {
    return NULL;
  }
  else if (Is5UTR (sfp))
  {
    return NULL;
  }
  else if (IsLTR (sfp))
  {
    return GetLTRDescription (sfp);
  }
  else
  {
    return GetGeneProtDescription (fcp, bsp, rp);
  }
}

static void LIBCALLBACK GetSatelliteFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  CharPtr    semicolon, colon;
  GBQualPtr  qual;
  Boolean    found = FALSE;

  flp->description = NULL; 
  flp->typeword = StringSave ("sequence");
  flp->pluralizable = FALSE;
  flp->is_typeword_first = FALSE;
 
  if (featlist == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (main_feat == NULL) return;
  for (qual = main_feat->qual; qual != NULL && !found; qual = qual->next)
  {
    if (StringCmp (qual->qual, "satellite") == 0) 
    {
      flp->description = StringSave (qual->val);
      if ((semicolon = StringStr (flp->description, ";")) != NULL)
      {
        *semicolon = 0;
      }
      if ((colon = StringChr (flp->description, ':')) != NULL)
      {
        *colon = ' ';
      }
    }
  }
}

static void LIBCALLBACK GetPromoterFeatureLabel (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
)
{
  SeqFeatPtr main_feat;
  
  flp->description = NULL;
  flp->typeword = StringSave ("promoter region");

  if (featlist == NULL) return;
  main_feat = featlist->data.ptrvalue;
  if (main_feat == NULL) return;

  flp->description =  NULL;
  flp->pluralizable = FALSE;
  flp->is_typeword_first = FALSE;

}

/* This function temporarily removes a 3' UTR clause from the end of
 * a clause list so that it will not be included in the list of subfeatures
 * before a CDS in the definition line.
 * The 3' UTR clause should be put back if it was not the only clause in the
 * list.
 */
static ValNodePtr Remove3UTRFromEndOfFeatList (ValNodePtr PNTR featlist)
{
  ValNodePtr vnp, prev = NULL;
  
  if (featlist == NULL || *featlist == NULL) return NULL;
  
  for (vnp = *featlist; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    prev = vnp;
  }
  if (vnp->choice == DEFLINE_CLAUSEPLUS && Is3UTRClause (vnp->data.ptrvalue))
  {
    if (prev == NULL)
    {
      *featlist = NULL;
    }
    else
    {
      prev->next = NULL;        
    }
  }
  else
  {
    vnp = NULL;
  }
  return vnp;
}

static Uint1 GetMoleculeType (BioseqPtr bsp, Uint2     entityID);
static void ConsolidateClauses (
  ValNodePtr PNTR list,
  BioseqPtr  bsp,
  Uint1      biomol,
  Boolean    delete_now,
  DeflineFeatureRequestListPtr rp);


/* This function calculates the "interval" for a clause in the definition
 * line.  The interval could be an empty string, it could indicate whether
 * the location of the feature is partial or complete and whether or not
 * the feature is a CDS, the interval could be a description of the
 * subfeatures of the clause, or the interval could be a combination of the
 * last two items if the feature is a CDS.
 */
static CharPtr GetGenericInterval 
( FeatureClausePtr fcp,
  Uint1            biomol,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr    interval;
  Boolean    partial5, partial3;
  SeqFeatPtr sfp;
  ValNodePtr featlist, strings, prev_feat;
  CharPtr    subfeatlist;
  Int4       len;
  Boolean    suppress_final_and;
  ValNodePtr utr3vnp = NULL;
  ValNodePtr last_feat;
  Uint1      molecule_type;

  if ( fcp == NULL || fcp->featlist == NULL) return NULL;
  if (fcp->is_unknown) return NULL;
  featlist = fcp->featlist;
  sfp = featlist->data.ptrvalue;
  if (sfp == NULL) return NULL;
  if ( IsExon (sfp) && fcp->is_alt_spliced)
  {
    return StringSave ("alternatively spliced");
  }
  if ( ! ShowInterval (sfp)) return NULL;

  if (IsIntergenicSpacer (sfp) && StringNCmp (sfp->comment, "may contain ", 12) == 0) {
    return StringSave ("region");
  }

  subfeatlist = NULL;
  len = 50;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  strings = NULL;
  prev_feat = NULL;
  while (featlist != NULL && featlist->choice != DEFLINE_CLAUSEPLUS)
  {
    prev_feat = featlist;
    featlist = featlist->next;
  }
  if (IsCDS (sfp))
  {
    utr3vnp = Remove3UTRFromEndOfFeatList (&featlist);
  }
  if (featlist != NULL)
  {
    suppress_final_and = FALSE;
    if (( IsCDS (sfp) && ! fcp->clause_info_only)
        || utr3vnp != NULL)
    {
      suppress_final_and = TRUE;
    }
    LabelClauses (featlist, biomol, bsp, rp);
    
    molecule_type = GetMoleculeType (bsp, bsp->idx.entityID);
    /* consolidate genes/proteins with the same names (usually hypothetical proteins) */
    ConsolidateClauses (&featlist, bsp, molecule_type, TRUE,
                        rp);

    /* make sure featlist is still intact - may have consolidated it */
    if (prev_feat == NULL)
    {
      fcp->featlist = featlist;
    }
    else
    {
      prev_feat->next = featlist;
    }

    ListClauses (featlist, &strings, FALSE, suppress_final_and, rp->suppress_allele);
    subfeatlist = MergeValNodeStrings (strings, FALSE);
      ValNodeFreeData (strings);
    len += StringLen (subfeatlist) + 7;
    
    if (utr3vnp != NULL)
    {
      len += 14;
    }
  }

  interval = (CharPtr) MemNew (len * sizeof (Char));
  if (interval == NULL) return NULL;
  interval[0] = 0;

  if (StringDoesHaveText (subfeatlist))
  {
    StringCat (interval, subfeatlist);
    if ( ! IsCDS (sfp) || fcp->clause_info_only)
    {
      if (utr3vnp != NULL)
      {
        if (featlist->next != NULL)
        {
          StringCat (interval, ",");
        }
        StringCat (interval, " and 3' UTR");
        /* put 3' UTR back at end of featlist */
        if (featlist != NULL)
        {
          last_feat = featlist;
          while (last_feat != NULL && last_feat->next != NULL)
          {
            last_feat = last_feat->next;
          }
          last_feat->next = utr3vnp;
        }
      }
      if (subfeatlist != NULL) MemFree (subfeatlist);
      return interval;
    }
    if (utr3vnp == NULL)
    {
      StringCat (interval, " and ");
    }
    else
    {
      StringCat (interval, ", ");
    }
  }

  if (FeatureDoesNotGetPartialComplete (sfp)) 
  {
    /* don't add partial or complete */
  } 
  else if (partial5 || partial3)
  {
    StringCat (interval, "partial ");
  }
  else
  {
    StringCat (interval, "complete ");
  }
  if (IsCDS (sfp) && ! IsPseudo (sfp))
  {
    StringCat (interval, "cds");
    if (fcp->is_alt_spliced)
      StringCat (interval, ", alternatively spliced");
  }
  else
  {
    StringCat (interval, "sequence");
    if (IsNoncodingProductFeat (sfp) && fcp->is_alt_spliced)
    {
      StringCat (interval, ", alternatively spliced");
    }
  }
  
  if (utr3vnp != NULL)
  {
    /* tack UTR3 on at end of clause */
    if (StringDoesHaveText (subfeatlist))
    {
      StringCat (interval, ",");
    }
    StringCat (interval, " and 3' UTR");
    
    /* put 3' UTR back at end of featlist */
    if (featlist != NULL)
    {
      last_feat = featlist;
      while (last_feat != NULL && last_feat->next != NULL)
      {
        last_feat = last_feat->next;
      }
      last_feat->next = utr3vnp;
    }
  }
  
  if (subfeatlist != NULL) MemFree (subfeatlist);    
  
  return interval;
} 


/* This function is used to generate feature label information for
 * a feature clause.  It is called by the LabelFeature function if
 * a "GetFeatureLabel" function is not found for the specific feature
 * type.
 * In the future it may be advisable to create "GetFeatureLabel" functions
 * for more of the specific feature types, to reduce the number of times
 * that the feature must be identified as being a certain type.
 */ 
static void LIBCALLBACK GetGenericFeatureLabel 
( FeatureClausePtr fcp,
  BioseqPtr        bsp,
  Uint1            biomol,
  FeatureLabelPtr  flp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr main_feat;
  
  if (fcp == NULL
    || fcp->featlist == NULL
    || fcp->featlist->data.ptrvalue == NULL)
  {
    return;
  }
  main_feat = fcp->featlist->data.ptrvalue;
  if (main_feat == NULL) return;

  if (flp->typeword == NULL)
  {
    flp->typeword = GetFeatureTypeWord (biomol, main_feat);
    flp->is_typeword_first = IsTypeWordFirst (flp->typeword);
    flp->pluralizable = TRUE;
  }
  if (flp->productname == NULL)
  {
    flp->productname = GetProductName (main_feat, bsp, rp);
  }
  if (flp->description == NULL
    && (! IsMiscRNA (main_feat)
      || StringStr (flp->productname, "spacer") == NULL ))
  {
    flp->description = GetFeatureDescription (fcp, bsp, rp);
  }

}

typedef void (LIBCALLBACK *GetFeatureLabelFunction) (
  ValNodePtr      featlist,
  BioseqPtr       bsp,
  Uint1           biomol,
  FeatureLabelPtr flp
);

typedef struct matchlabelfunction {
  matchFunction           itemmatch;
  GetFeatureLabelFunction labelfunction;
} MatchLabelFunctionData, PNTR MatchLabelFunctionPtr;

static MatchLabelFunctionData label_functions[] = {
 { IsSatelliteSequence, GetSatelliteFeatureLabel         },
 { IsMobileElement,     GetMobileElementFeatureLabel        },
 { IsPromoter,          GetPromoterFeatureLabel          },
 { IsIntergenicSpacer,  GetIntergenicSpacerFeatureLabel  },
 { IsGeneCluster,       GetGeneClusterFeatureLabel       }
};

typedef enum {
 DEFLINE_FEATLABEL_Satellite = 0,
 DEFLINE_FEATLABEL_Transposon,
 DEFLINE_FEATLABEL_Promoter,
 DEFLINE_FEATLABEL_IntergenicSpacer,
 DEFLINE_FEATLABEL_GeneCluster,
 NumDefLineFeatLabels
} DefLineFeatLabel;

static void LabelFeature 
( BioseqPtr        bsp,
  Uint1            biomol,
  FeatureClausePtr new_clauseplus,
  DeflineFeatureRequestListPtr rp)
{
  Int4             i;
  SeqFeatPtr       main_feat;

  if (new_clauseplus == NULL || new_clauseplus->featlist == NULL) return;

  if (new_clauseplus->featlist->choice == DEFLINE_FEATLIST)
  {
    main_feat = (SeqFeatPtr) new_clauseplus->featlist->data.ptrvalue;
    
    new_clauseplus->allelename = GetAlleleName (new_clauseplus->grp,
                                                rp == NULL ? FALSE : rp->suppress_locus_tags);
    if (new_clauseplus->interval == NULL)
    {
      new_clauseplus->interval =
                  GetGenericInterval (new_clauseplus, biomol, bsp, rp);
    }

    for (i=0; i < NumDefLineFeatLabels; i++)
    {
      if (label_functions [i].itemmatch (main_feat))
      {
        label_functions [i].labelfunction ( new_clauseplus->featlist,
                                          bsp, biomol,
                                          &new_clauseplus->feature_label_data);
        return;
      }
    }

    GetGenericFeatureLabel ( new_clauseplus, bsp, biomol, 
                           &new_clauseplus->feature_label_data, rp);
    return;
  }
}

/* This function is used to calculate the parts of a product name that
 * are "the same" for use as the name of an alternatively spliced product.
 * The common portion of the string must end at a recognized separator,
 * such as a space, comma, or dash instead of in the middle of a word.
 * The matching portions of the string could occur at the beginning or end
 * of the string, or even occasionally at the beginning and end of a
 * string, but not as the center of the string with a different beginning
 * and ending.
 */
static CharPtr FindStringIntersection (
  CharPtr str1,
  CharPtr str2,
  Boolean str1_previously_stripped
)
{
  Int4 matchleftlen = 0;
  Int4 matchlefttoken = 0;
  Int4 matchrightidx1 = 0;
  Int4 matchrightidx2 = 0;
  Int4 matchrighttoken = 0;
  CharPtr match_string;
  Int4 len1;
  Int4 len2;
  Int4 match_len;

  if (str1 == NULL || str2 == NULL) return NULL;
  if (StringCmp (str1, str2) == 0) return StringSave (str1);
  len1 = StringLen (str1);
  len2 = StringLen (str2);

  while (str1[matchleftlen] != 0 && str2[matchleftlen] != 0
         && str1[matchleftlen] == str2[matchleftlen])
  {
    if (str1 [matchleftlen] == ','
      || str1 [matchleftlen] == '-')
    {
      matchlefttoken = matchleftlen;
    }
    else if (str1 [matchleftlen] == ' '
      && matchlefttoken != matchleftlen - 1)
    {
      matchlefttoken = matchleftlen;
    }
    matchleftlen++;
  }
  if (matchleftlen == len1 && str1_previously_stripped) 
  {
    matchlefttoken = matchleftlen;
  }
  else
  {
    matchleftlen = matchlefttoken;
  }

  matchrightidx1 = len1;
  matchrightidx2 = len2;

  while (matchrightidx1 > -1 && matchrightidx2 > -1
         && str1[matchrightidx1] == str2[matchrightidx2])
  {
    if (str1 [matchrightidx1] == ' '
      || str1[matchrightidx1] == ','
      || str1[matchrightidx1] == '-')
    {
      matchrighttoken = matchrightidx1;
    }
    matchrightidx1--;
    matchrightidx2--;
  }
  if (matchrightidx1 == -1)
  {
    matchrighttoken = matchrightidx1;
  }
  else if (matchrighttoken > 0) 
  {
    matchrightidx1 = matchrighttoken;
  } 
  else if (str1_previously_stripped && matchrightidx1 < len1 - 1)
  {
    /* matchrightidx1 = matchrighttoken; */
    /* do nothing, leave right index where it is */
  }
  else
  {
    matchrightidx1 = len1;
  }

  match_len = matchleftlen;
  if (matchrightidx1 < len1 - 1)
  {
    match_len += len1 - matchrightidx1 - 1;
  }

  if (match_len <= 0) return NULL;

  match_string = MemNew (match_len + 2);
  if (match_string == NULL) return NULL;
  if (matchleftlen != 0)
  {
    StringNCpy (match_string, str1, matchleftlen);
    match_string[matchleftlen] = 0;
  }
  else
  {
    match_string[0] = 0;
  }
  if (matchrightidx1 < len1)
  {
    if (match_string[0] != 0) StringCat (match_string, " ");
    StringCat (match_string, str1 + matchrightidx1 + 1);
  }
  return match_string;
}

/* These are the words that are used to introduced the part of the protein
 * name that differs in alt-spliced products - they should not be part of
 * the alt-spliced product name.
 * Note that splice variant is listed before "variant" so that it will be
 * found first and "variant" will not be removed from "splice variant", leaving
 * splice as an orphan.
 */

static CharPtr UnwantedWords [] = {
 "splice variant",
 "splice product",
 "variant",
 "isoform"
};

static void TrimUnwantedWordsFromAltSpliceProductName (
  CharPtr productname
)
{
  Int4    num_unwanted_words, i;
  size_t  unwanted_word_len;
  CharPtr cp, tmp;
  
  num_unwanted_words = sizeof (UnwantedWords) / sizeof (CharPtr);
  for (i = 0; i < num_unwanted_words; i++)
  {
    unwanted_word_len = StringLen (UnwantedWords [i]);
    cp = StringStr (productname, UnwantedWords [i]);
    if (cp != NULL)
    {
      if (cp == productname)
      {
        /* word occurs in beginning of phrase */
        tmp = StringSave (productname + unwanted_word_len);
        StringCpy (productname, tmp);
        MemFree (tmp);
      }
      else if (cp - productname < StringLen (productname) - unwanted_word_len)
      {
        /* word occurs in middle of phrase */
        tmp = StringSave (cp + unwanted_word_len);
        StringCpy (cp - 1, tmp);
        MemFree (tmp);
      }
      else
      {
        /* word occurs at end of phrase */
        *cp = 0;
      }
    }
  }
}


static Boolean PreviouslyStripped (SeqFeatPtr cds, BioseqPtr bsp, CharPtr productname)
{
  CharPtr expected_product_name;
  Boolean rval = FALSE;
  
  if (cds == NULL || StringHasNoText (productname)) return FALSE;
  expected_product_name = GetProductName (cds, bsp, FALSE);
  if (StringCmp (productname, expected_product_name) != 0) {
    rval = TRUE;
  }
  expected_product_name = MemFree (expected_product_name);
  return rval;
}

/* This function determines whether two CDSs meet the conditions for
 * alternative splicing, and if so, it returns the name of the alternatively
 * spliced product.  In order to be alternatively spliced, the two CDSs 
 * must have the same gene, must share a complete interval, and must have
 * similarly named products.
 */
static CharPtr MeetAltSpliceRules 
( FeatureClausePtr cdsfcp1,
  FeatureClausePtr cdsfcp2,
  BioseqPtr        bsp,
  DeflineFeatureRequestListPtr rp)
{
  SeqFeatPtr cds1, cds2;
  CharPtr match_string;
  Int4    res;
  
  if (cdsfcp1 == NULL || cdsfcp2 == NULL
    || cdsfcp1->featlist == NULL || cdsfcp2->featlist == NULL)
  {
    return NULL;
  }

  cds1 = cdsfcp1->featlist->data.ptrvalue;
  cds2 = cdsfcp2->featlist->data.ptrvalue;
  if (! DoGenesMatch (cdsfcp1->grp, cdsfcp2->grp, rp == NULL ? FALSE : rp->suppress_locus_tags))
    return NULL;

  if ( (res = TestFeatOverlap (cds1, cds2, COMMON_INTERVAL)) != -1)
  {
    match_string = FindStringIntersection (
                     cdsfcp1->feature_label_data.productname,
                     cdsfcp2->feature_label_data.productname,
                     PreviouslyStripped(cds1, bsp, cdsfcp1->feature_label_data.productname));
    return match_string;
  }
  return NULL;
}

/* This function is used by the FindAltSplices function to locate the
 * next CDS in a list of feature clauses.
 */
static ValNodePtr FindNextCDSClause (ValNodePtr vnp)
{
  FeatureClausePtr fcp;

  while (vnp != NULL)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp != NULL && !fcp->delete_me && fcp->featlist != NULL
        && IsCDS (fcp->featlist->data.ptrvalue))
      {
        return vnp;
      }
    }
    vnp = vnp->next;
  }
  return NULL;
}

/* This function is used by the FindAltSplices function to move the features
 * and subclauses from the second CDS in an alternatively spliced pair of 
 * CDSs to the feature clause for the first CDS, so that the subfeatures
 * can be properly listed.
 */
static void MoveSubclauses (
  FeatureClausePtr dstfcp,
  FeatureClausePtr srcfcp
)
{
  ValNodePtr dst_last_feat, dst_first_clause, dst_last_clause;
  ValNodePtr src_last_feat, src_first_clause;

  if (dstfcp == NULL || srcfcp == NULL || srcfcp->featlist == NULL) return;

  dst_first_clause = NULL;
  dst_last_clause = NULL;
  src_first_clause = NULL;

  dst_last_feat = dstfcp->featlist;
  while (dst_last_feat != NULL 
      && dst_last_feat->next != NULL
      && dst_last_feat->next->choice == DEFLINE_FEATLIST)
  {
    dst_last_feat = dst_last_feat->next;
  }
  if (dst_last_feat != NULL)
  {
    dst_first_clause = dst_last_feat->next;
  }
  dst_last_clause = dst_first_clause;
  while (dst_last_clause != NULL && dst_last_clause->next != NULL)
  {
    dst_last_clause = dst_last_clause->next;
  }

  src_last_feat = srcfcp->featlist;
  while (src_last_feat != NULL 
      && src_last_feat->next != NULL
      && src_last_feat->next->choice == DEFLINE_FEATLIST)
  {
    src_last_feat = src_last_feat->next;
  }
  if (src_last_feat != NULL)
  {
    src_first_clause = src_last_feat->next;
  }

  /* insert features before clauses */
  if (dst_last_feat == NULL)
  {
    dstfcp->featlist = srcfcp->featlist;
    dst_last_feat = src_last_feat;
  }
  else
  {
    dst_last_feat->next = srcfcp->featlist;
  }
  /* insert clauses after feats */
  if (dst_first_clause != NULL)
  {
    src_last_feat->next = dst_first_clause;
    dst_last_clause->next = src_first_clause;
  }
  srcfcp->featlist = NULL;
}
 
/* we want to look through the list for CDS features */
/* if we find two CDSs that are alternatively spliced, */
/* we replace the first alternatively spliced CDS feature */
/* with a new CDS feature that has the new protein name as */
/* a comment and a data.choice value that indicates alt splicing */
/* we remove the second alternatively spliced CDS feature from the list */

static void FindAltSplices 
( ValNodePtr clause_list,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr  fcp1, fcp2;
  ValNodePtr cdsclause1, cdsclause2;
  ValNodePtr searchclause;
  CharPtr  combined_protein_name;
  Boolean    partial3_1, partial5_1, partial3_2, partial5_2;
  Int4       left1, left2, right1, right2;

  if (clause_list == NULL) return;
  
  cdsclause1 = FindNextCDSClause (clause_list);
  while (cdsclause1 != NULL)
  {
    fcp1 = (FeatureClausePtr) cdsclause1->data.ptrvalue;
    if (fcp1->feature_label_data.productname == NULL)
    {
      fcp1->feature_label_data.productname = 
           GetProductName (fcp1->featlist->data.ptrvalue, bsp, rp);
    }
    searchclause = cdsclause1->next;
    cdsclause2 = FindNextCDSClause (searchclause);
    while (cdsclause2 != NULL) 
    {
      fcp2 = (FeatureClausePtr) cdsclause2->data.ptrvalue;
      if (fcp2->feature_label_data.productname == NULL)
      {
        fcp2->feature_label_data.productname =
           GetProductName (fcp2->featlist->data.ptrvalue, bsp, rp);
      }
      combined_protein_name = MeetAltSpliceRules (fcp1, fcp2, bsp, rp);
      if (combined_protein_name != NULL)
      {
        /* get rid of variant, splice variant, splice product, isoform, etc.*/
        TrimUnwantedWordsFromAltSpliceProductName (combined_protein_name);

        /* get rid of trailing spaces in protein name */
        TrimSpacesAroundString (combined_protein_name);

        /* copy new protein name into first clause */
        MemFree (fcp1->feature_label_data.productname);
        fcp1->feature_label_data.productname = combined_protein_name;
        CheckSeqLocForPartial (fcp1->slp, &partial5_1, &partial3_1);
        left1 = GetOffsetInBioseq (fcp1->slp, bsp, SEQLOC_LEFT_END);
        right1 = GetOffsetInBioseq (fcp1->slp, bsp, SEQLOC_RIGHT_END);
        CheckSeqLocForPartial (fcp2->slp, &partial5_2, &partial3_2);
        left2 = GetOffsetInBioseq (fcp2->slp, bsp, SEQLOC_LEFT_END);
        right2 = GetOffsetInBioseq (fcp2->slp, bsp, SEQLOC_RIGHT_END);
        fcp1->slp = SeqLocMerge (bsp, fcp1->slp, fcp2->slp,
                                 FALSE, TRUE, FALSE);
        if (left1 == left2)
        {
          partial5_1 |= partial5_2;
        }
        else
        {
          partial5_1 = left1 < left2 ? partial5_1 : partial5_2;
        }
        if (right1 == right2)
        {
          partial3_1 |= partial3_2;
        }
        else
        {
          partial3_1 = right1 > right2 ? partial3_1 : partial3_2;
        }
        SetSeqLocPartial (fcp1->slp, partial5_1, partial3_1);
        fcp1->is_alt_spliced = TRUE;

        /* copy over fcp2 subclauses */
        MoveSubclauses (fcp1, fcp2);

        /* remove second clause */
        fcp2->delete_me = TRUE;
      }
      searchclause = cdsclause2->next;
      cdsclause2 = FindNextCDSClause (searchclause);
    }
    cdsclause1 = FindNextCDSClause (cdsclause1->next);
  } 
  DeleteFeatureClauses (&clause_list);
}

static void LabelClauses 
( ValNodePtr clause_list,
  Uint1      biomol,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr clause;
 
  clause = clause_list;
  while (clause != NULL)
  { 
    LabelFeature ( bsp, biomol, clause->data.ptrvalue, rp);
    clause = clause->next;
  }
}

static CharPtr misc_words [] = {
  "internal transcribed spacer",
  "external transcribed spacer",
  "ribosomal RNA intergenic spacer",
  "ribosomal RNA",
  "intergenic spacer"
};

typedef enum {
  MISC_RNA_WORD_INTERNAL_SPACER = 0,
  MISC_RNA_WORD_EXTERNAL_SPACER,
  MISC_RNA_WORD_RNA_INTERGENIC_SPACER,
  MISC_RNA_WORD_RNA,
  MISC_RNA_WORD_INTERGENIC_SPACER,
  NUM_MISC_RNA_WORDS
} MiscWord;

/* note - must put substrings of other separators after the longer version */
static CharPtr separators [] = {
  ", and ",
  " and ",
  ", ",
  "; "
};

#define num_separators 3


static ValNodePtr TokenListFromMiscRNAString (CharPtr str)
{
  ValNodePtr token_list = NULL;
  CharPtr cansep [num_separators];
  CharPtr token_start, next_sep, token;
  Int4    i, sep_len, datalen;
  Uint1   word_i;
  Boolean found_unparseable = FALSE;

  if ( StringStr (str, "spacer") == NULL) {
    return NULL;
  }

  token_start = str;
  for (i = 0; i < num_separators; i++) {
    cansep[i] = StringStr (token_start, separators[i]);
  }

  while (*token_start != 0 && !found_unparseable) {
    next_sep = NULL;
    sep_len = 0;
    for (i = 0; i < num_separators; i++) {
      if (cansep[i] != NULL) {
        if (cansep[i] < token_start) {
          cansep[i] = StringStr (token_start, separators[i]);
        }
      }
      if (cansep[i] != NULL && (next_sep == NULL || next_sep > cansep[i])) {
        next_sep = cansep[i];
        sep_len = StringLen (separators[i]);
      }
    }
    if (next_sep == NULL) {
      token = StringSave (token_start);
      datalen = StringLen (token);
    } else {
      datalen = next_sep - token_start;
      token = (CharPtr) MemNew (sizeof (Char) * (datalen + 1));
      StringNCpy (token, token_start, datalen);
      token[datalen] = 0;
    }
    /* determine which word is part of the token */
    for (word_i=0;
         word_i < NUM_MISC_RNA_WORDS
           && StringStr (token, misc_words [word_i]) == NULL;
         word_i++) {}
    if (word_i < NUM_MISC_RNA_WORDS) {
      ValNodeAddPointer (&token_list, word_i, token);
    } else {
      found_unparseable = TRUE;
    }      
    token_start += datalen + sep_len;
  }
  if (found_unparseable) {
    token_list = ValNodeFreeData (token_list);
  }
  return token_list;
}


static ValNodePtr 
GetFeatureClausesFromMiscRNATokens 
( ValNodePtr token_list,
  SeqFeatPtr misc_rna,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp) 
{
  ValNodePtr clause_list = NULL;
  ValNodePtr vnp;
  Boolean    partial5, partial3, unparseable = FALSE;
  CharPtr    word_loc;
  FeatureClausePtr fcp;

  if (token_list == NULL || misc_rna == NULL) {
    return NULL;
  }

  CheckSeqLocForPartial (misc_rna->location, &partial5, &partial3);
  
  for (vnp = token_list; vnp != NULL && !unparseable; vnp = vnp->next) {
    word_loc = StringStr (vnp->data.ptrvalue, misc_words [vnp->choice]);
    if (word_loc == NULL) {
      unparseable = TRUE;
    } else {
      fcp = NewFeatureClause ( misc_rna, bsp, rp);
      if (fcp == NULL) {
        unparseable = TRUE;
      } else {
        if (vnp->choice == MISC_RNA_WORD_INTERNAL_SPACER
            || vnp->choice == MISC_RNA_WORD_EXTERNAL_SPACER
            || vnp->choice == MISC_RNA_WORD_RNA_INTERGENIC_SPACER
            || vnp->choice == MISC_RNA_WORD_INTERGENIC_SPACER) {
          if (word_loc == vnp->data.ptrvalue) {
            fcp->feature_label_data.is_typeword_first = TRUE;
            fcp->feature_label_data.typeword = StringSave (misc_words [vnp->choice]);
            if (StringLen (misc_words [vnp->choice]) + 1 < StringLen (vnp->data.ptrvalue)) {
              fcp->feature_label_data.description =
                    StringSave ( ((CharPtr)vnp->data.ptrvalue) + StringLen (misc_words [vnp->choice]) + 1);
            }
          } else {
            fcp->feature_label_data.is_typeword_first = FALSE;
            fcp->feature_label_data.typeword = StringSave (misc_words [vnp->choice]);
            if (StringLen (misc_words [vnp->choice]) + 1 < StringLen (vnp->data.ptrvalue)) {
              fcp->feature_label_data.description = StringSave ( vnp->data.ptrvalue);
              fcp->feature_label_data.description [word_loc - ((CharPtr) vnp->data.ptrvalue) - 1] = 0;
            }
          }
        } else if (vnp->choice == MISC_RNA_WORD_RNA) {
          fcp->feature_label_data.description = StringSave (vnp->data.ptrvalue);
        }
        if ((vnp == token_list && partial5) || (vnp->next == NULL && partial3)) {
          fcp->interval = StringSave ("partial sequence");
        } else {
          fcp->interval = StringSave ("complete sequence");
        }
        ValNodeAddPointer (&clause_list, DEFLINE_CLAUSEPLUS, fcp);
      }
    }
  }
  if (unparseable) {
    DefLineFeatClauseListFree (clause_list);
    clause_list = NULL;
  }
  return clause_list;
}


static Boolean s_ItemNeedsGene (CharPtr item_str)
{
  Int4 item_len;

  item_len = StringLen (item_str);
  if (StringSearch (item_str, "RNA") != NULL 
      && !((item_len > 4 && StringCmp (item_str + item_len - 4, "gene") == 0)
      || (item_len > 5 && StringCmp (item_str + item_len - 5, "genes") == 0))) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static CharPtr GetRegionDescription 
(ValNodePtr token_list,
 SeqFeatPtr misc_rna,
 BioseqPtr bsp,
 DeflineFeatureRequestListPtr rp)
{
  ValNodePtr vnp;
  CharPtr desc, item_str;
  Int4    len = 5, item_len;

  if (token_list == NULL) {
    return NULL;
  } 

  for (vnp = token_list; vnp != NULL; vnp = vnp->next) {
    item_str = (CharPtr) vnp->data.ptrvalue;
    item_len = StringLen (item_str);
    len += item_len + 2;
    if (s_ItemNeedsGene(item_str)) {
      len += 5;
    }
  }

  desc = (CharPtr) MemNew (sizeof (Char) * (len));
  desc[0] = 0;
  for (vnp = token_list; vnp != NULL; vnp = vnp->next) {
    StringCat (desc, vnp->data.ptrvalue);
    if (s_ItemNeedsGene(vnp->data.ptrvalue)) {
      StringCat (desc, " gene");
    }
    if (vnp->next != NULL) {
      StringCat (desc, ", ");
      if (vnp->next->next == NULL) {
        StringCat (desc, "and ");
      }
    }
  }

  return desc;
}


static ValNodePtr GetMiscRNAelements 
( SeqFeatPtr misc_rna,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr buffer;
  ValNodePtr token_list, clause_list = NULL;
  FeatureClausePtr fcp;

  if (misc_rna == NULL) return NULL;
  buffer = GetProductName (misc_rna, bsp, rp);
  if (buffer == NULL) 
  {
    buffer = StringSave (misc_rna->comment);
  }
  else if (StringNCmp (buffer, misc_rna->comment, StringLen (buffer) -1) == 0
    && buffer [ StringLen (buffer) - 1] == '>')
  {
    MemFree (buffer);
    buffer = StringSave (misc_rna->comment);
  }

  if (StringNCmp (buffer, "contains ", 9) == 0) {
    token_list = TokenListFromMiscRNAString (buffer + 9);
    clause_list = GetFeatureClausesFromMiscRNATokens (token_list, misc_rna, bsp, rp);
    token_list = ValNodeFreeData (token_list);
  } else if (StringNCmp (buffer, "may contain ", 12) == 0) {
    token_list = TokenListFromMiscRNAString (buffer + 12);
    if (token_list != NULL) {
      fcp = NewFeatureClause ( misc_rna, bsp, rp);
      fcp->feature_label_data.description = GetRegionDescription (token_list, misc_rna, bsp, rp);
      fcp->feature_label_data.typeword = StringSave ("");
      fcp->interval = StringSave ("region");
      ValNodeAddPointer (&clause_list, DEFLINE_CLAUSEPLUS, fcp);
    }
    token_list = ValNodeFreeData (token_list);
  } else {
    token_list = TokenListFromMiscRNAString (buffer);
    clause_list = GetFeatureClausesFromMiscRNATokens (token_list, misc_rna, bsp, rp);
    token_list = ValNodeFreeData (token_list);
  }

  buffer = MemFree (buffer);
  return clause_list;
}

/* Some misc_RNA clauses have a comment that actually lists multiple
 * features.  This function creates a clause for each element in the
 * comment and inserts the list of new clauses into the feature list
 * at the point where the single previous clause was.
 */
static void ReplaceRNAClauses (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  SeqFeatPtr main_feat;
  ValNodePtr clause, replacement_clauses, nextclause, vnp;

  if (clause_list == NULL || *clause_list == NULL) return;
  clause = *clause_list;
  while (clause != NULL)
  {
    nextclause = clause->next;
    fcp = (clause->data.ptrvalue);
    if (fcp == NULL
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      return;
    }
    main_feat = (SeqFeatPtr) fcp->featlist->data.ptrvalue;
  
    if (IsrRNA (main_feat) || IsMiscRNA (main_feat))
    {
      replacement_clauses = GetMiscRNAelements ( main_feat, bsp, rp );
      if (replacement_clauses != NULL)
      {
        for (vnp = replacement_clauses; vnp->next != NULL; vnp = vnp->next) {}
        vnp->next = clause->next;
        clause->next = replacement_clauses;
        fcp->delete_me = TRUE;
      }
    }
    clause = nextclause;
  }
  DeleteFeatureClauses (clause_list);
}

/* Some misc_feat clauses have a comment that lists one or more tRNAs and
 * an intergenic spacer.  This function creates a clause for each element 
 * in the comment and inserts the list of new clauses into the feature list
 * at the point where the single previous clause was.
 */
static void ReplaceIntergenicSpacerClauses (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  SeqFeatPtr main_feat;
  ValNodePtr clause, replacement_clauses, nextclause, vnp;

  if (clause_list == NULL || *clause_list == NULL) return;
  clause = *clause_list;
  while (clause != NULL)
  {
    nextclause = clause->next;
    fcp = (clause->data.ptrvalue);
    if (fcp == NULL
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      return;
    }
    main_feat = (SeqFeatPtr) fcp->featlist->data.ptrvalue;
  
    if (IsParsableList (main_feat)) 
    {
      if ((replacement_clauses = ParsetRNAIntergenicSpacerElements ( main_feat, bsp, rp)) != NULL)
      {
        for (vnp = replacement_clauses; vnp->next != NULL; vnp = vnp->next) {}
        vnp->next = clause->next;
        clause->next = replacement_clauses;
        fcp->delete_me = TRUE;
      }
      else
      {
        fcp->delete_me = TRUE;
      }
    }
    clause = nextclause;
  }
  DeleteFeatureClauses (clause_list);
}

/* If we are applying a different rule for misc_feats, we need to recalculate
 * their descriptions.
 */
static void RenameMiscFeats (ValNodePtr clause_list, Uint1 biomol)
{
  ValNodePtr       vnp, featlist;
  FeatureClausePtr fcp, featlistclause;
  SeqFeatPtr       sfp;
  Int4             name_len;

  for (vnp = clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS || vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = vnp->data.ptrvalue;
    for (featlist = fcp->featlist; featlist != NULL; featlist = featlist->next)
    {
      if ( featlist->data.ptrvalue == NULL)
      {
        continue;
      }
      if (featlist->choice == DEFLINE_CLAUSEPLUS)
      {
        featlistclause = featlist->data.ptrvalue;
        RenameMiscFeats (featlistclause->featlist, biomol);
        continue;
      }
      if (featlist->choice != DEFLINE_FEATLIST)
      {
        continue;
      }
      sfp = featlist->data.ptrvalue;
      if (sfp->idx.subtype != FEATDEF_misc_feature
        || sfp->comment == NULL
        || IsIntergenicSpacer (sfp)
        || IsGeneCluster (sfp)
        || IsControlRegion (sfp)) 
      {
        continue;
      }
      if (fcp->feature_label_data.description != NULL) 
      {
        fcp->feature_label_data.description
                   = MemFree (fcp->feature_label_data.description);
      }
      name_len = StringCSpn (sfp->comment, ";");
      /* make sure we have space for terminating NULL */
      fcp->feature_label_data.description = MemNew ((name_len + 1) * sizeof (Char));
      if (fcp->feature_label_data.description == NULL) return;
      StringNCpy (fcp->feature_label_data.description, sfp->comment, name_len);
      fcp->feature_label_data.description [ name_len ] = 0;
      fcp->feature_label_data.typeword =
            MemFree (fcp->feature_label_data.typeword);
      if (StringCmp (fcp->feature_label_data.description + name_len - 9, " sequence") == 0) 
      {
        fcp->feature_label_data.description[name_len - 9] = 0;
        fcp->feature_label_data.typeword = StringSave ("sequence");
      }
      else if (biomol == MOLECULE_TYPE_GENOMIC)
      {
        fcp->feature_label_data.typeword = StringSave ("genomic sequence");
      }
      else if (biomol == MOLECULE_TYPE_MRNA)
      {
        fcp->feature_label_data.typeword = StringSave ("mRNA sequence");
      }
      else
      {
        fcp->feature_label_data.typeword = StringSave ("sequence");
      }
     
      fcp->interval = MemFree (fcp->interval);
      fcp->interval = StringSave ("");
    }
  }
}

static void RemoveUnwantedMiscFeats (
  ValNodePtr PNTR clause_list,
  Boolean delete_now
)
{
  ValNodePtr       vnp, featlist;
  FeatureClausePtr fcp, featlistclause;
  SeqFeatPtr       sfp;

  for (vnp = *clause_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS || vnp->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = vnp->data.ptrvalue;
    for (featlist = fcp->featlist; featlist != NULL; featlist = featlist->next)
    {
      if ( featlist->data.ptrvalue == NULL)
      {
        continue;
      }
      if (featlist->choice == DEFLINE_CLAUSEPLUS)
      {
        featlistclause = featlist->data.ptrvalue;
        RemoveUnwantedMiscFeats (&(featlistclause->featlist), FALSE);
        continue;
      }
      if (featlist->choice != DEFLINE_FEATLIST)
      {
        continue;
      }
      sfp = featlist->data.ptrvalue;
      if ( sfp->idx.subtype == FEATDEF_misc_feature
        && ! IsNoncodingProductFeat (sfp)
        && ! IsControlRegion (sfp)
        && ! IsIntergenicSpacer (sfp)
        && ! IsGeneCluster (sfp)
        && ! IsParsableList (sfp))
      {
        fcp->delete_me = TRUE;
      }
    }
  }
  DeleteFeatureClauses (clause_list);
}

/* When a feature is on the minus strand, the clauses are listed by 
 * sequence indexing in reverse biological order - we reverse the subclauses
 * for the feature in order to have them listed in the definition line
 * in biological order.
 * This is most noticeable when the main feature is a CDS with multiple
 * exons numbered sequentially.  If the exons are on the minus strand and
 * appear as 9, 8, 7, 6, we want to display them in the definition line as
 * 6, 7, 8, 9.
 */
static void ReverseClauses (
  ValNodePtr PNTR clause_list,
  matchFunction   itemmatch
)
{
  ValNodePtr vnp, last_feat, first_feat, next_item, new_list;
  FeatureClausePtr fcp;

  if (clause_list == NULL || *clause_list == NULL) return;

  last_feat = NULL;
  first_feat = NULL;
  new_list = NULL;
  vnp = *clause_list;
  while (vnp != NULL)
  {
    next_item = vnp->next;
    fcp = NULL;
    if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL
      && fcp->slp != NULL
      && SeqLocStrand (fcp->slp) == Seq_strand_minus
      && fcp->featlist != NULL
      && fcp->featlist->choice == DEFLINE_FEATLIST
      && itemmatch (fcp->featlist->data.ptrvalue))
    {
      vnp->next = new_list;
      new_list = vnp;
    }
    else
    {
      if (first_feat == NULL)
      {
        first_feat = vnp;
        last_feat = vnp;
      }
      else
      {
        last_feat->next = vnp;
        last_feat = vnp;
        last_feat->next = NULL;
      }
    }
    if (fcp != NULL)
    {
      ReverseClauses (&(fcp->featlist), itemmatch);
    }
    vnp = next_item;
  }
  if (first_feat == NULL)
  {
    *clause_list = new_list;
  }
  else
  {
    last_feat->next = new_list;
    *clause_list = first_feat;
  }
}

/* This function is used to determine whether two features are both exons
 * and whether they are numerically sequential - i.e., exon 7 and exon 8
 * are a pair of consecutive exons, exon 7 and exon 9 are not, and exon 7
 * and intron 9 are not.
 */
static Boolean ClausePairIsTwoConsecutiveExons (
  ValNodePtr vnp1,
  ValNodePtr vnp2,
  BioseqPtr  bsp
)
{
  FeatureClausePtr fcp1, fcp2;
  SeqFeatPtr       exon1, exon2;
  Int4 num1, num2;
  CharPtr          exdesc1, exdesc2;

  if (vnp1 == NULL || vnp2 == NULL 
    || vnp1->choice != DEFLINE_CLAUSEPLUS
    || vnp2->choice != DEFLINE_CLAUSEPLUS
    || vnp1->data.ptrvalue == NULL
    || vnp2->data.ptrvalue == NULL)
  {
    return FALSE;
  }
  fcp1 = vnp1->data.ptrvalue;
  fcp2 = vnp2->data.ptrvalue;
  if ( fcp1->featlist == NULL
    || fcp1->featlist->data.ptrvalue == NULL
    || fcp2->featlist == NULL
    || fcp2->featlist->data.ptrvalue == NULL
    || fcp1->featlist->choice != DEFLINE_FEATLIST
    || fcp2->featlist->choice != DEFLINE_FEATLIST
    || ! IsExon (fcp1->featlist->data.ptrvalue)
    || ! IsExon (fcp2->featlist->data.ptrvalue)
    || (fcp1->is_alt_spliced && ! fcp2->is_alt_spliced)
    || (! fcp1->is_alt_spliced && fcp2->is_alt_spliced))
  {
    return FALSE;
  }

  exon1 = (SeqFeatPtr)(fcp1->featlist->data.ptrvalue);
  exon2 = (SeqFeatPtr)(fcp2->featlist->data.ptrvalue);

  exdesc1 = GetExonDescription (bsp, exon1);
  exdesc2 = GetExonDescription (bsp, exon2);
  if (exdesc1 == NULL || exdesc2 == NULL)
  {
    if (exdesc1 != NULL) MemFree (exdesc1);
    if (exdesc2 != NULL) MemFree (exdesc2);
    return FALSE;
  }
  
  num1 = atoi (exdesc1);
  num2 = atoi (exdesc2);
  MemFree (exdesc1);
  MemFree (exdesc2);

  if (abs (num1 - num2) == 1)
  {
    return TRUE;
  }

  return FALSE; 
}

/* This function counts the number of consecutive exons in a list.
 */
static Int4 GetNumberOfConsecutiveExons (
  ValNodePtr list,
  BioseqPtr  bsp
)
{
  ValNodePtr check;
  Int4       num_exons;
 
  num_exons = 0;
  check = list->next;
  if ( ! ClausePairIsTwoConsecutiveExons (list, check, bsp)) return 0;
  
  num_exons = 2;
  while ( check != NULL
    && ClausePairIsTwoConsecutiveExons (check, check->next, bsp))
  {
    num_exons++;
    check = check->next;
  }
  return num_exons;
}

/* This function replaces a list of three or more consecutive exon clauses
 * with a single "summary" clause that gives the range of exons present -
 * i.e., if you have exons 1, 2, 3, and 4, a clause will be created that
 * contains all four of those features and has a description of "1 through 4".
 */
static void ReplaceExonClauseList (
  FeatureClausePtr fcp,
  ValNodePtr       clause,
  Int4             num_exons,
  BioseqPtr        bsp
)
{
  ValNodePtr       lastfeat, tmpclause;
  FeatureClausePtr tmpfcp;
  Int4             i;
  CharPtr          new_description;
  Int4             new_description_len;
  CharPtr          exdesc1, exdesc2;

  if (fcp == NULL || clause == NULL) return;

  lastfeat = fcp->featlist;
  while (lastfeat != NULL && lastfeat->next != NULL)
  {
    lastfeat = lastfeat->next;
  }
  tmpclause = clause->next;
  for (i=0; i < num_exons - 1 && tmpclause != NULL; i++)
  {
    tmpfcp = tmpclause->data.ptrvalue;
    tmpfcp->delete_me = TRUE;
    if (lastfeat == NULL)
    {
      fcp->featlist = tmpfcp->featlist;
    }
    else
    {
      lastfeat->next = tmpfcp->featlist;
    }
    tmpfcp->featlist = NULL;
    while (lastfeat != NULL && lastfeat->next != NULL)
    {
      lastfeat = lastfeat->next;
    }
          
    tmpclause = tmpclause->next;
  }

  exdesc1 = GetExonDescription (bsp, fcp->featlist->data.ptrvalue);
  exdesc2 = GetExonDescription (bsp, lastfeat->data.ptrvalue);
  if (exdesc1 == NULL || exdesc2 == NULL)
  {
    if (exdesc1 != NULL) MemFree (exdesc1);
    if (exdesc2 != NULL) MemFree (exdesc2);
    return;
  }
  new_description_len =
        StringLen (exdesc1)
      + StringLen (exdesc2)
      + StringLen (" through ")
      + 1;
  new_description = MemNew (new_description_len * sizeof (Char));
  if (new_description == NULL) return;
  sprintf (new_description, "%s through %s", exdesc1, exdesc2);
  MemFree (exdesc1);
  MemFree (exdesc2);
  if (fcp->feature_label_data.description != NULL)
  {
    MemFree (fcp->feature_label_data.description);
  }
  fcp->feature_label_data.description = new_description;
}

/* This function recursively searches for lists of consecutive exons
 * and calls ReplaceExonClauseList to consolidate the exons into a list
 * clause.
 */
static void RenameExonSequences (
  ValNodePtr PNTR list,
  BioseqPtr       bsp,
  Boolean         delete_now
)
{
  ValNodePtr       clause;
  Int4             num_exons;
  FeatureClausePtr fcp;
 
  if (list == NULL) return; 
  clause = *list;
  while (clause != NULL)
  {
    if (clause->choice == DEFLINE_CLAUSEPLUS
      && clause->data.ptrvalue != NULL)
    {
      fcp = clause->data.ptrvalue;
      if ( ! fcp->delete_me)
      {
        num_exons = GetNumberOfConsecutiveExons (clause, bsp);
        if (num_exons > 2)
        {
          ReplaceExonClauseList (fcp, clause, num_exons, bsp);
        }
        else
        {
          RenameExonSequences (&fcp->featlist, bsp, FALSE);
        }
      }
    }
    clause = clause->next;
  }
  if (delete_now) DeleteFeatureClauses (list);
}

static CharPtr organelleByGenome [] = {
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
  "",
  "",
  "",
  "apicoplast",
  "leucoplast",
  "proplastid",
  "",
  "hydrogenosome",
  "",
  "chromatophore",
  NULL,
};

static CharPtr organelleByPopup [] = {
  NULL,
  "mitochondrial",
  "chloroplast",
  "kinetoplast",
  "plastid",
  "chromoplast",
  "cyanelle",
  "apicoplast",
  "leucoplast",
  "proplastid",
  NULL
};

static void 
AddProductEnding 
(CharPtr    str, 
 BioseqPtr  bsp,
 Int2       mitochloroflag,
 ValNodePtr strings)
{
  Char orgnelle [80];
  BioSourcePtr  biop;
  ValNodePtr last_string;
  Int4 num_genes;
  SubSourcePtr  ssp;

  num_genes = 0;
  biop = GetBiopForBsp (bsp);

  if (biop != NULL) {
    if (FindStringInStrings (strings, "genes"))
    {
      num_genes = 2;
    }
    else if ((last_string = FindStringInStrings (strings, "gene")) != NULL
      && last_string->next != NULL
      && (last_string = FindStringInStrings (last_string->next, "gene")) != NULL)
    {
      num_genes = 2;
    }
    else
    {
      num_genes = 1;
    }

    orgnelle [0] = '\0';
  
    switch (biop->genome) {
    case GENOME_macronuclear :
      StringCat (str, "; macronuclear");
      break;
    case GENOME_nucleomorph :
      StringCat (str, "; nucleomorph");
      break;
    case GENOME_apicoplast :
    case GENOME_chloroplast :
    case GENOME_chromoplast :
    case GENOME_kinetoplast :
    case GENOME_mitochondrion :
    case GENOME_plastid :
    case GENOME_cyanelle :
    case GENOME_leucoplast :
    case GENOME_proplastid :
    case GENOME_hydrogenosome :
    case GENOME_chromatophore :
      sprintf (orgnelle, "; %s", organelleByGenome [biop->genome]);
      StringCat (str, orgnelle);
      break;
    default :
      ssp = biop->subtype;
      while (ssp != NULL && ssp->subtype != 255)
      {
        ssp = ssp->next;
      }
      if (ssp != NULL
        && ssp->name != NULL
        && StringStr (ssp->name, "micronuclear"))
      {
        StringCat (str, "; micronuclear");
      }
      else if (mitochloroflag > 0) {
        if (mitochloroflag > 9) {
          /* beyond list */
        }
        else {
          if (num_genes > 1)
          {
            sprintf (orgnelle, "; nuclear genes for %s products",
                     organelleByPopup [mitochloroflag]);
          }
          else 
          {
            sprintf (orgnelle, "; nuclear gene for %s product",
                     organelleByPopup [mitochloroflag]);
          }
          StringCat (str, orgnelle);
        }
      }
      break;
    }
  }  
}

/*---------------------------------------------------------------------*/
/*                                                                     */
/* AutoDef_AddEnding () -- Add an ending on to the definition line     */
/*                         after the last feature.                     */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void AutoDef_AddEnding (
  ValNodePtr   clause_list,
  ValNodePtr PNTR strings,
  BioseqPtr    bsp,
  Int2         mitochloroflag,
  Boolean      alternate_splice_flag
)
{
  Char str [200];
  ValNodePtr last_string;
  Int4 new_data_len;
  CharPtr new_data;

  str[0] = 0;
  AddProductEnding (str, bsp, mitochloroflag, *strings);
  if (alternate_splice_flag) {
    StringCat (str, ", alternatively spliced");
  }

  StringCat (str, ".");

  last_string = *strings;
  if (last_string == NULL)
  {
    ValNodeAddStr (strings, 0, StringSave ( str));
  }
  else
  {
    while (last_string->next != NULL) last_string = last_string->next;
    new_data_len = StringLen (last_string->data.ptrvalue) + StringLen (str) + 1;
    new_data = (CharPtr) MemNew (new_data_len);
    if (new_data == NULL) return;
    StringCpy (new_data, last_string->data.ptrvalue);
    StringCat (new_data, str);
    MemFree (last_string->data.ptrvalue);
    last_string->data.ptrvalue = new_data;
  }
}

static Boolean LastIntervalChangeBeforeEnd (
  FeatureClausePtr onebefore,
  FeatureClausePtr thisclause,
  ValNodePtr rest_of_list
)
{
  ValNodePtr       vnp;
  FeatureClausePtr fcp;

  if (onebefore == NULL || rest_of_list == NULL) return FALSE;
 
  if (StringCmp (onebefore->interval, thisclause->interval) == 0) return FALSE;
  
  for (vnp = rest_of_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      fcp = vnp->data.ptrvalue;
      if (StringCmp (thisclause->interval, fcp->interval) != 0) return FALSE;
    }
  }
  return TRUE;
    
}

static void PluralizeClauseIntervals (
  FeatureClausePtr fcp
)
{
  CharPtr new_interval, cp;

  if (fcp->interval != NULL
    && (cp = StringStr (fcp->interval, "gene, ")) != NULL)
  {
    new_interval = MemNew (StringLen (fcp->interval) + 2);
    if (new_interval == NULL) return;
    StringCpy (new_interval, fcp->interval);
    new_interval [ cp - fcp->interval + 4] = 's';
    StringCpy (new_interval + (cp - fcp->interval) + 5,
               cp + 4);
    MemFree (fcp->interval);
    fcp->interval = new_interval;
  }
}  

static Boolean DisplayAlleleName (FeatureClausePtr thisclause)
{
  if (thisclause == NULL) return FALSE;
  if (StringCmp (thisclause->feature_label_data.typeword, "gene") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "pseudogene") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "mRNA") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "pseudogene mRNA") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "precursor RNA") == 0
    || StringCmp (thisclause->feature_label_data.typeword, "pseudogene precursor RNA") == 0)
  {
    return TRUE;
  }
  return FALSE;
}

static void ListClauses (
  ValNodePtr clauselist,
  ValNodePtr PNTR strings,
  Boolean    allow_semicolons,
  Boolean    suppress_final_and,
  Boolean    suppress_allele
)
{
  FeatureClausePtr thisclause, onebefore, twobefore, oneafter, twoafter;
  Boolean print_typeword;
  Boolean print_and;
  Boolean print_comma;
  Boolean print_semicolon;
  Boolean print_comma_between_description_and_typeword;
  Boolean typeword_is_plural;
  size_t clause_len;
  CharPtr clause_string;
  Boolean oneafter_has_detail_change;
  Boolean oneafter_has_interval_change;
  Boolean oneafter_has_typeword_change;
  Boolean onebefore_has_detail_change;
  Boolean onebefore_has_interval_change;
  Boolean onebefore_has_typeword_change;
  SeqFeatPtr main_feat;
  CharPtr new_interval;
  ValNodePtr voneafter, vtwoafter;

  while (clauselist != NULL && clauselist->choice != DEFLINE_CLAUSEPLUS)
  {
    clauselist = clauselist->next;
  }
  if (clauselist == NULL) return;

  thisclause = clauselist->data.ptrvalue;
  onebefore = NULL;
  twobefore = NULL;
  
  while (thisclause != NULL)
  {
    oneafter_has_detail_change = FALSE;
    oneafter_has_interval_change = FALSE;
    oneafter_has_typeword_change = FALSE;
    onebefore_has_detail_change = FALSE;
    onebefore_has_interval_change = FALSE;
    onebefore_has_typeword_change = FALSE;
    if (onebefore != NULL)
    {
      if (StringCmp (onebefore->interval, thisclause->interval) != 0)
        onebefore_has_interval_change = TRUE;
      if (StringCmp (onebefore->feature_label_data.typeword,
                     thisclause->feature_label_data.typeword) != 0)
      {
        onebefore_has_typeword_change = TRUE;
      }
      if (onebefore_has_typeword_change || onebefore_has_interval_change
          || (!suppress_allele && DisplayAlleleName (onebefore) && StringLen (onebefore->allelename) != 0)
          || (!suppress_allele && DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) != 0))
     {
        onebefore_has_detail_change = TRUE;  
      }
    }
    voneafter = clauselist->next;
    while (voneafter != NULL && voneafter->choice != DEFLINE_CLAUSEPLUS)
    {
      voneafter = voneafter->next;
    }
    if (voneafter == NULL)
    {
      vtwoafter = NULL;
    }
    else
    {
      vtwoafter = voneafter->next;
      while (vtwoafter != NULL && vtwoafter->choice != DEFLINE_CLAUSEPLUS)
      {
        vtwoafter = vtwoafter->next;
      }
    }

    if (voneafter != NULL)
    {
      oneafter = voneafter->data.ptrvalue;
      if (StringCmp (oneafter->interval, thisclause->interval) != 0)
        oneafter_has_interval_change = TRUE;
      if (StringCmp (oneafter->feature_label_data.typeword,
                     thisclause->feature_label_data.typeword) != 0)
      {
        oneafter_has_typeword_change = TRUE;
      }
      if (oneafter_has_typeword_change  || oneafter_has_interval_change
          || (!suppress_allele && DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) != 0)
          || (!suppress_allele && DisplayAlleleName (oneafter) && StringLen (oneafter->allelename) != 0))
      {
        oneafter_has_detail_change = TRUE;
      }
      if (vtwoafter != NULL)
      {
        twoafter = vtwoafter->data.ptrvalue;
      }
      else
      {
        twoafter = NULL;
      }
    }
    else
    {
      oneafter = NULL;
      twoafter = NULL;
    }
    print_typeword = FALSE;
    typeword_is_plural = FALSE;
    print_and = FALSE;
    print_comma = FALSE;
    print_semicolon = FALSE;

    if (thisclause->feature_label_data.is_typeword_first)
    {
      if (onebefore == NULL || onebefore_has_detail_change)
      {
        print_typeword = TRUE;
        if (oneafter != NULL && ! oneafter_has_detail_change)
        {
          typeword_is_plural = TRUE;
        }
        else if (StringStr (thisclause->feature_label_data.description, " through ") != NULL
          && StringCmp (thisclause->feature_label_data.typeword, "exon") == 0)
        {
          typeword_is_plural = TRUE;
        }
      }
    }
    else
    {
      if (oneafter == NULL || oneafter_has_detail_change)
      {
        print_typeword = TRUE;
        if (onebefore != NULL && ! onebefore_has_detail_change)
        {
          typeword_is_plural = TRUE;
        }
      }
    }

    /* when to print and before this section */
    if ( onebefore != NULL
         && ! onebefore_has_detail_change
         && (oneafter == NULL || oneafter_has_detail_change))
    {
      print_and = TRUE;
    }
    else if (oneafter == NULL && onebefore != NULL)
    {
      print_and = TRUE;
    }
    else if (onebefore != NULL
         && ! onebefore_has_interval_change
         && oneafter_has_interval_change)
    {
      print_and = TRUE;
    }
    else if ( LastIntervalChangeBeforeEnd ( onebefore, 
                                            thisclause,
                                            clauselist->next))
    {
      print_and = TRUE;
    }

    if (suppress_final_and && oneafter == NULL)
    {
      print_and = FALSE;
    }
    if (suppress_final_and && oneafter != NULL && twoafter == NULL)
    {
      print_comma = TRUE;
    }
    
    /* when to print semicolon after this section */
    /* after every interval change except when exons change "interval" */
    /* exons changing interval are going from alt-spliced to not */
    /* or vice versa, in either case we don't want a semicolon or comma */
    if (oneafter != NULL && oneafter_has_interval_change
      && (StringCmp (thisclause->feature_label_data.typeword, "exon") != 0
         || StringCmp (oneafter->feature_label_data.typeword, "exon") != 0))
    {
      print_semicolon = TRUE;
    }

    /* when to print comma after this section */
    if (onebefore != NULL && oneafter != NULL
      && ! onebefore_has_detail_change
      && ! oneafter_has_detail_change )
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && onebefore != NULL
      && ! onebefore_has_interval_change && ! oneafter_has_interval_change
      &&  onebefore_has_typeword_change &&  oneafter_has_typeword_change)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL
      && ! oneafter_has_detail_change
      && StringCmp (twoafter->feature_label_data.typeword,
                    thisclause->feature_label_data.typeword) == 0
      && StringCmp (twoafter->interval,
                    thisclause->interval) == 0)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL  && twoafter != NULL
      && oneafter_has_typeword_change
      && StringCmp (twoafter->feature_label_data.typeword,
                    oneafter->feature_label_data.typeword) == 0
      && StringCmp (twoafter->interval,
                    oneafter->interval) == 0
      && ! print_and)
    {
      print_comma = TRUE;
    }
    else if (((oneafter_has_interval_change || oneafter == NULL)
      && StringDoesHaveText (thisclause->interval))
      || (oneafter_has_interval_change && oneafter != NULL && ! print_semicolon))
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL
      && !oneafter_has_interval_change
      && StringCmp (thisclause->interval, twoafter->interval) == 0
      && oneafter_has_typeword_change
      && StringCmp (thisclause->feature_label_data.typeword,
                    twoafter->feature_label_data.typeword) != 0)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && onebefore != NULL && twoafter != NULL
      && ! oneafter_has_interval_change && ! onebefore_has_interval_change
      && StringCmp (thisclause->interval, twoafter->interval) == 0
      && oneafter_has_typeword_change)
    {
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL
      && oneafter_has_typeword_change
      && StringCmp (oneafter->feature_label_data.typeword,
                    twoafter->feature_label_data.typeword) != 0
      && ! oneafter_has_interval_change
      && StringCmp (oneafter->interval, twoafter->interval) == 0)
    {
      /* spacer 1, foo RNA gene, and spacer2, complete sequence */
      /*         ^ */
      print_comma = TRUE;
    }
    else if (oneafter != NULL && twoafter != NULL 
      && ! oneafter_has_interval_change && StringCmp (thisclause->interval, twoafter->interval) == 0
      && ((!suppress_allele && DisplayAlleleName (oneafter) && StringLen (oneafter->allelename) > 0)
        || (!suppress_allele && DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) > 0)))
    {
      print_comma = TRUE;      
    }
    else if (oneafter != NULL && onebefore != NULL
      && ! oneafter_has_interval_change && ! onebefore_has_interval_change
      && ((!suppress_allele && DisplayAlleleName (oneafter) && StringLen (oneafter->allelename) > 0)
        || (!suppress_allele && DisplayAlleleName (thisclause) && StringLen (thisclause->allelename) > 0)))
    {
      print_comma = TRUE;      
    }

    if (thisclause->featlist != NULL
      && thisclause->featlist->data.ptrvalue != NULL
      && StringDoesHaveText (thisclause->interval)
      && StringNCmp (thisclause->interval, "partial", 7) != 0
      && StringNCmp (thisclause->interval, "complete", 8) != 0)
    {
      main_feat = thisclause->featlist->data.ptrvalue;
      if (IsMobileElement (main_feat)
        || IsEndogenousVirusSourceFeature (main_feat) )
      {
        print_comma = FALSE;
      }
    }

    if (onebefore != NULL
      && ! onebefore_has_interval_change
      && (oneafter_has_interval_change || oneafter == NULL))
    {
      PluralizeClauseIntervals (thisclause);
    }

    if ( thisclause->make_plural )
    {
      if ((onebefore != NULL && ! onebefore_has_detail_change)
        || (oneafter != NULL && !oneafter_has_detail_change))
      {
        PluralizeConsolidatedClauseDescription (thisclause);
      }
      else
      {
        typeword_is_plural = TRUE;
      }
    }

    clause_len = StringLen (thisclause->feature_label_data.description) + 1;
    
    /* add one in case we need to add the semicolon to this clause (when
     * the interval has changed because this clause has no interval and
     * the next one does).
     */
    clause_len++;

    /* we need to place a comma between the description and the type word 
     * when the description ends with "precursor" or when the type word
     * starts with "precursor"
     */
    if ( thisclause->feature_label_data.description != NULL
      && ! thisclause->feature_label_data.is_typeword_first
      && print_typeword
      && ! StringHasNoText (thisclause->feature_label_data.typeword)
      && ((StringNCmp (thisclause->feature_label_data.typeword, "precursor", 9) == 0
            && thisclause->feature_label_data.description [StringLen (thisclause->feature_label_data.description) - 1] != ')')
          || (clause_len > StringLen ("precursor")
              && StringCmp ( thisclause->feature_label_data.description
                     + clause_len - StringLen ("precursor") - 2,
                     "precursor") == 0)))
    {
      print_comma_between_description_and_typeword = TRUE;
      clause_len += 1;
    }
    else
    {
      print_comma_between_description_and_typeword = FALSE;
    }

    if (print_typeword)
      clause_len += StringLen (thisclause->feature_label_data.typeword) + 1;
    if (typeword_is_plural)
      clause_len += 1;
    if (print_and)
      clause_len += 4;
    if (print_comma)
      clause_len += 2;
    if (!suppress_allele && DisplayAlleleName (thisclause))
    {
      clause_len += StringLen (thisclause->allelename) + 10;
      if (StringLen (thisclause->allelename) > 0) 
      {
        clause_len += StringLen (thisclause->allelename) + StringLen ("allele ");
      }
    }
    
    clause_string = (CharPtr) MemNew (clause_len);
    if (clause_string == NULL)
      return;
    clause_string[0] = 0;
    if (print_and)
      StringCat (clause_string, "and ");
    if (thisclause->feature_label_data.is_typeword_first && print_typeword
      && thisclause->feature_label_data.typeword != NULL
      && ! StringHasNoText (thisclause->feature_label_data.typeword))
    {
      StringCat (clause_string, thisclause->feature_label_data.typeword);
      if (typeword_is_plural)
        StringCat (clause_string, "s");
      if (thisclause->feature_label_data.description != NULL)
        StringCat (clause_string, " ");
    }
    if (thisclause->feature_label_data.description != NULL)
    {
      StringCat (clause_string, thisclause->feature_label_data.description);
      if (print_comma_between_description_and_typeword)
      {
        StringCat (clause_string, ",");
      }
    }
    if (! thisclause->feature_label_data.is_typeword_first && print_typeword
      && thisclause->feature_label_data.typeword != NULL
      && ! StringHasNoText (thisclause->feature_label_data.typeword))
    {
      if (!StringHasNoText (thisclause->feature_label_data.description))
        StringCat (clause_string, " ");
      StringCat (clause_string, thisclause->feature_label_data.typeword);
      if (typeword_is_plural)
        StringCat (clause_string, "s");
      if (!suppress_allele && DisplayAlleleName (thisclause)
        && thisclause->allelename != NULL)
      {
        StringCat (clause_string, ", ");
        StringCat (clause_string, thisclause->allelename);
        StringCat (clause_string, " allele");
      }
    }
    if (StringLen (clause_string) > 0 ) 
    {
      if (print_comma)
        StringCat (clause_string, ",");
      ValNodeAddStr (strings, 0, clause_string);
    }
    else 
    {
        MemFree (clause_string);
        clause_string = NULL;
    }
 
    if (oneafter == NULL || oneafter_has_interval_change)
    {
      if (print_semicolon) {
        if (thisclause->interval == NULL
          || StringHasNoText(thisclause->interval)) {
          if (clause_string != NULL) {
            StringCat (clause_string, ";");
          }
        } else if (thisclause->interval[StringLen (thisclause->interval)] != ';') {
          new_interval = MemNew (StringLen (thisclause->interval) + 2);
          if (new_interval == NULL) return;
          StringCpy (new_interval, thisclause->interval);
          if (allow_semicolons) 
          {
            StringCat (new_interval, ";");
          }
          else
          {
            StringCat (new_interval, ",");
          }
          MemFree (thisclause->interval);
          thisclause->interval = new_interval;
        }
      }
      if (thisclause->interval != NULL
        && !StringHasNoText (thisclause->interval))
      {
        ValNodeAddStr (strings, 0, StringSave (thisclause->interval));
      }
    }
    twobefore = onebefore;
    onebefore = thisclause;   
    thisclause = oneafter;
    clauselist = voneafter;
  }
}

static Uint1 GetMoleculeType 
(BioseqPtr bsp,
 Uint2     entityID)
{
  SeqDescPtr         sdp;
  MolInfoPtr         mip;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL) return MOLECULE_TYPE_GENOMIC;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp == NULL) return MOLECULE_TYPE_GENOMIC;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return MOLECULE_TYPE_GENOMIC;
  return mip->biomol;
}

static Boolean SpecialHandlingForSpecialTechniques (
  BioseqPtr bsp
)
{
  SeqDescPtr sdp;
  MolInfoPtr mip;

  if (bsp == NULL) return MOLECULE_TYPE_GENOMIC;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, NULL);
  if (sdp == NULL)
  {
    for (sdp = bsp->descr;
         sdp != NULL && sdp->choice != Seq_descr_molinfo;
         sdp = sdp->next)
    {}
  }
  if (sdp == NULL) return FALSE;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return FALSE;
  if (mip->tech == MI_TECH_htgs_0 ||
      mip->tech == MI_TECH_htgs_1 ||
      mip->tech == MI_TECH_htgs_2 ||
      mip->tech == MI_TECH_est ||
      mip->tech == MI_TECH_sts ||
      mip->tech == MI_TECH_survey ||
      mip->tech == MI_TECH_wgs) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, NULL);
    if (sdp != NULL) {
      return TRUE;
    }
  }

  return FALSE;
}

static Boolean LIBCALLBACK ShouldRemoveExon (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  Boolean partial3, partial5;
  SeqFeatPtr main_feat;

  if (isSegment || isLonely || isRequested) return FALSE;
  if (parent_fcp == NULL
    || parent_fcp->featlist == NULL
    || parent_fcp->featlist->data.ptrvalue == NULL)
  {
    return TRUE;
  }
  
  main_feat = parent_fcp->featlist->data.ptrvalue;
  if ( IsCDS (main_feat))
  {
    CheckSeqLocForPartial (main_feat->location, &partial5, &partial3);
    if (partial5 || partial3) return FALSE;
  }
  else if (IsmRNA (main_feat) || parent_fcp->has_mrna)
  {
    return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK ShouldRemoveCDS (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr description;
  Boolean retval = FALSE;

  description = GetGeneProtDescription (this_fcp, bsp, rp);
  if (StringHasNoText (description))
  {
    retval = TRUE;
  }
  if (description != NULL) MemFree (description);
  return retval;
}

static Boolean LIBCALLBACK ShouldRemoveNoncodingProductFeat (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp, Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isRequested) return FALSE;
  return TRUE;
}

static Boolean LIBCALLBACK ShouldRemovePromoter (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp, Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  /* remove a promoter if it is in an mRNA or gene clause */
  if (isRequested)
  {
    return FALSE;
  }
  else if (parent_fcp != NULL 
      && (parent_fcp->has_mrna 
        || (parent_fcp->featlist != NULL
           && parent_fcp->featlist->choice == DEFLINE_FEATLIST
           && parent_fcp->featlist->data.ptrvalue != NULL
           && IsmRNA (parent_fcp->featlist->data.ptrvalue))))
  {
    return TRUE;
  }
  else if (isLonely)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static Boolean LIBCALLBACK ShouldRemoveLTR (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isRequested)
  {
    return FALSE;
  }
  else if (parent_fcp != NULL)
  {
    return TRUE;
  }
  else if (isLonely)
    return FALSE;
  else
    return TRUE;
}

static Boolean LIBCALLBACK ShouldRemove3UTR (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{ 
  if (isLonely || isRequested)
    return FALSE;
  else
    return TRUE;
}

static Boolean LIBCALLBACK ShouldRemove5UTR (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isLonely || isRequested)
    return FALSE;
  else
    return TRUE;
}


static Boolean LIBCALLBACK ShouldRemoveuORF (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isLonely || isRequested)
    return FALSE;
  else
    return TRUE;
}


static Boolean LIBCALLBACK ShouldRemoveIntron (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp, Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp
)
{
  if (isRequested)
  {
    return FALSE;
  }
  else if (parent_fcp != NULL 
      && (parent_fcp->has_mrna 
        || (parent_fcp->featlist != NULL
           && parent_fcp->featlist->choice == DEFLINE_FEATLIST
           && parent_fcp->featlist->data.ptrvalue != NULL
           && IsmRNA (parent_fcp->featlist->data.ptrvalue))))
  {
    return TRUE;
  }
  else if (isLonely)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static Boolean LIBCALLBACK ShouldRemoveMobileElement
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  return (!isLonely && !isRequested);
}

static Boolean LIBCALLBACK ShouldRemovencRNA
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  Boolean rval = FALSE;
  SeqMgrFeatContext context;
  SeqFeatPtr precursor_rna;
    Int2 cmpval;

  if (isRequested) {
    return FALSE;
  }
  for (precursor_rna = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_preRNA, &context);
       precursor_rna != NULL && !rval;
       precursor_rna = SeqMgrGetNextFeature (bsp, precursor_rna, 0, FEATDEF_preRNA, &context)) {
    cmpval = SeqLocCompare (sfp->location, precursor_rna->location);
    if (cmpval != SLC_NO_MATCH) {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean LIBCALLBACK ShouldRemoveGeneric 
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  CharPtr productname;
  Boolean rval;

  rval = FALSE;
  if (IsMiscRNA (sfp) && ( productname = GetProductName (sfp, bsp, rp)) != NULL)
  {
    if (StringStr (productname, "trans-spliced leader") != NULL)
    {
      rval = TRUE;
    }
    MemFree (productname);
  }
    
  return rval;
}


static Boolean IsBioseqPrecursorRNA (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  MolInfoPtr        mol;
  
  if (bsp == NULL) return FALSE;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
  if (sdp != NULL && sdp->data.ptrvalue != NULL)
  {
        mol = (MolInfoPtr) sdp->data.ptrvalue;
    if (mol->biomol == 2)
    {
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean LIBCALLBACK ShouldRemovePrecursorRNA
( SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isRequested,
  Boolean isSegment,
  DeflineFeatureRequestListPtr rp)
{
  if (!isLonely && IsBioseqPrecursorRNA(bsp))
  {
    return TRUE;
  }
  else
  {
    return ShouldRemoveGeneric (sfp, parent_fcp, this_fcp, bsp, isLonely, 
                                isRequested, isSegment, rp);
  }
}


typedef struct removableitemglobal {
  matchFunction  itemmatch;
  ShouldRemoveFunction ShouldRemove;
  CharPtr  group_name;
} RemovableItemGlobalData, PNTR RemovableItemGlobalPtr;

typedef struct removableitemlocal {
/*  ButtoN  keep_request; */
  Boolean  keep;
} RemovableItemLocalData, PNTR RemovableItemLocalPtr;

static RemovableItemGlobalData remove_items[] = {
  { IsExon, ShouldRemoveExon, "Exons" },
  { IsIntron, ShouldRemoveIntron, "Introns" },
  { Is5UTR, ShouldRemove5UTR, "5' UTRs" },
  { Is3UTR, ShouldRemove3UTR, "3' UTRs" },
  { IsuORF, ShouldRemoveuORF, "uORFs"},
  { IsCDS,  ShouldRemoveCDS, "CDSs" },
  { IsPromoter, ShouldRemovePromoter, "Promoters:" },
  { IsLTR, ShouldRemoveLTR, "LTRs" },
  { IsNoncodingProductFeat,  ShouldRemoveNoncodingProductFeat, "Misc feats with comments:" },
  { IsRemovableMobileElement, ShouldRemoveMobileElement, "Optional Mobile Element" },
  { IsPrecursorRNA, ShouldRemovePrecursorRNA, "Precursor RNAs" },
  { IsncRNA, ShouldRemovencRNA, "ncRNAs that overlap precursor RNAs"}
};


NLM_EXTERN CharPtr GetRemovableItemName (Int4 i)
{
  if (i < 0 || i >= NumRemovableItems) {
    return NULL;
  } else {
    return remove_items[i].group_name;
  }
}

NLM_EXTERN void InitFeatureRequests (
  DeflineFeatureRequestListPtr feature_requests
)
{
  Int4 i;
  for (i=0; i < NumRemovableItems; i++)
  {
    feature_requests->keep_items[i] = FALSE;
  }
  feature_requests->add_fake_promoters = TRUE;
  feature_requests->suppress_alt_splice_phrase = FALSE;
  feature_requests->remove_subfeatures = FALSE;
  feature_requests->feature_list_type = DEFLINE_USE_FEATURES;
  feature_requests->misc_feat_parse_rule = 2;
  feature_requests->suppress_locus_tags = FALSE;
  feature_requests->suppressed_feature_list = NULL;
  feature_requests->use_ncrna_note = FALSE;
  feature_requests->suppress_allele = FALSE;
}


NLM_EXTERN DeflineFeatureRequestListPtr FreeDeflineFeatureRequestList (DeflineFeatureRequestListPtr feature_requests)
{
  if (feature_requests != NULL) {
    feature_requests->suppressed_feature_list = ValNodeFree (feature_requests->suppressed_feature_list);
    feature_requests = MemFree (feature_requests);
  }
  return feature_requests;
}


static Boolean RemoveCondition (
  SeqFeatPtr sfp,
  FeatureClausePtr parent_fcp,
  FeatureClausePtr this_fcp,
  BioseqPtr bsp,
  Boolean isLonely,
  Boolean isSegment,
  DeflineFeatureRequestList *feature_requests
)
{
  Int4 i;
  if (sfp == NULL) return TRUE;
  for (i=0; i < NumRemovableItems; i++)
  {
    if (remove_items[i].itemmatch (sfp))
      return remove_items[i].ShouldRemove (sfp, parent_fcp, this_fcp, bsp,
                                           isLonely, feature_requests->keep_items[i],
                                           isSegment,
                                           feature_requests);
  }
  return ShouldRemoveGeneric(sfp, parent_fcp, this_fcp, bsp, isLonely, FALSE,
                             isSegment, feature_requests);
}

static Boolean FindOtherGeneClause 
( ValNodePtr feature_list,
  ValNodePtr me,
  GeneRefPtr grp,
  Boolean    suppress_locus_tag)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;

  if (grp == NULL) return FALSE;

  for (vnp = feature_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp == me) continue;
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp->delete_me) continue;
      if ( fcp->grp == grp
        || (fcp->grp != NULL && DoGenesMatch (fcp->grp, grp, suppress_locus_tag)))
      {
        return TRUE;
      }
      if ( FindOtherGeneClause (fcp->featlist, me, grp, suppress_locus_tag))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}
 
static void RemoveGenesMentionedElsewhere 
( ValNodePtr PNTR feature_list,
  ValNodePtr      search_list,
  Boolean         delete_now,
  Boolean         suppress_locus_tag)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;

  for (vnp = *feature_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS && vnp->data.ptrvalue != NULL)
    {
      fcp = vnp->data.ptrvalue;
      if (fcp->featlist == NULL)
      {
        continue;
      }
      if ( IsGene (fcp->featlist->data.ptrvalue)
        && fcp->featlist->next == NULL
        && FindOtherGeneClause ( search_list, vnp, fcp->grp, suppress_locus_tag))
      {
        fcp->delete_me = TRUE;
      }
      else
      {
        RemoveGenesMentionedElsewhere ( &(fcp->featlist), search_list, FALSE, suppress_locus_tag);
      }
    }
  }
  if (delete_now)
  {
    DeleteFeatureClauses (feature_list);
  }
}

static void MarkUnwantedFeatureClauseForRemoval (
  ValNodePtr clause,
  BioseqPtr  bsp,
  Boolean    isLonely,
  FeatureClausePtr parent_fcp,
  Boolean    isSegment,
  DeflineFeatureRequestList PNTR feature_requests
)
{
  FeatureClausePtr fcp;
  ValNodePtr       featlist;
  ValNodePtr       firstfeat;
  Int4             clause_count;
  SeqFeatPtr       sfp;

  if (clause == NULL
    || clause->choice != DEFLINE_CLAUSEPLUS
    || clause->data.ptrvalue == NULL)
  {
    return;
  }

  fcp = clause->data.ptrvalue;
  firstfeat = fcp->featlist;
  clause_count = 0;
  for (featlist = firstfeat;
       featlist != NULL && isLonely;
       featlist = featlist->next)
  {
    if (featlist->choice == DEFLINE_CLAUSEPLUS)
    {
      clause_count ++;
      if (clause_count > 1)
      {
        isLonely = FALSE;
      }
    }
  }
    
  featlist = firstfeat;
  while (featlist != NULL)
  {  
    if (featlist->choice == DEFLINE_FEATLIST
      && featlist->data.ptrvalue != NULL)
    {
      sfp = (SeqFeatPtr) featlist->data.ptrvalue;
      if (RemoveCondition (featlist->data.ptrvalue, parent_fcp, fcp, bsp, 
                          isLonely, isSegment, feature_requests))
      {
        fcp->delete_me = TRUE;
      }
      else if (! IsGene (sfp) && ! IsmRNA (sfp))
      {
        isLonely = FALSE;
      }
    }
    else if (featlist->choice == DEFLINE_CLAUSEPLUS 
      && featlist->data.ptrvalue != NULL)
    {
      MarkUnwantedFeatureClauseForRemoval (featlist, bsp, isLonely, fcp,
                                           isSegment,
                                           feature_requests);
    }
    featlist = featlist->next;
  }
}
  
static void RemoveUnwantedFeatures (
  ValNodePtr PNTR list,
  BioseqPtr bsp,
  Boolean   isSegment,
  DeflineFeatureRequestList PNTR feature_requests
)
{
  ValNodePtr       vnp;
  Boolean          isLonely;

  if (list == NULL) return;

  isLonely = TRUE;

  for (vnp = *list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->next != NULL) isLonely = FALSE;
    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      MarkUnwantedFeatureClauseForRemoval (vnp, bsp, isLonely, NULL,
                                           isSegment, feature_requests);
    }
  }
  DeleteFeatureClauses (list);
}

static Boolean IsFeatureInSelectionList (SeqFeatPtr sfp, ValNodePtr feat_list)
{
  ValNodePtr       vnp;

  if (sfp == NULL || feat_list == NULL)
  {
    return FALSE;
  }
  
  for (vnp = feat_list; vnp != NULL && sfp->idx.subtype != vnp->choice; vnp = vnp->next)
  {
  }
  if (vnp == NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}

static void MarkSuppressedFeatureClauseForRemoval (
  ValNodePtr clause,
  ValNodePtr suppressed_feature_list
)
{
  FeatureClausePtr fcp;
  ValNodePtr       featlist;
  ValNodePtr       firstfeat;
  SeqFeatPtr       sfp;

  if (clause == NULL
    || clause->choice != DEFLINE_CLAUSEPLUS
    || clause->data.ptrvalue == NULL)
  {
    return;
  }

  fcp = clause->data.ptrvalue;
  firstfeat = fcp->featlist;
    
  featlist = firstfeat;
  while (featlist != NULL)
  {  
    if (featlist->choice == DEFLINE_FEATLIST
      && featlist->data.ptrvalue != NULL)
    {
      sfp = (SeqFeatPtr) featlist->data.ptrvalue;
      if (IsFeatureInSelectionList (sfp, suppressed_feature_list))
      {
        fcp->delete_me = TRUE;
      }
    }
    else if (featlist->choice == DEFLINE_CLAUSEPLUS 
      && featlist->data.ptrvalue != NULL)
    {
      MarkSuppressedFeatureClauseForRemoval (featlist, suppressed_feature_list);
    }
    featlist = featlist->next;
  }
}
  
static void RemoveSuppressedFeatures (ValNodePtr PNTR list,
                                      ValNodePtr suppressed_feature_list)
{
  ValNodePtr vnp;
  
  if (list == NULL || *list == NULL || suppressed_feature_list == NULL)
  {
    return;
  }
  
  for (vnp = *list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      MarkSuppressedFeatureClauseForRemoval (vnp, suppressed_feature_list);
    }
  }
  DeleteFeatureClauses (list);  
}

static Boolean LIBCALLBACK IsMasterClause (
  SeqFeatPtr sfp
)
{
  if ( IsMobileElement (sfp)) return TRUE;
  return FALSE;
}

static void DeleteSubfeatures (
  ValNodePtr PNTR feature_list,
  Boolean         delete_now
)
{
  ValNodePtr       clause, featlist;
  FeatureClausePtr clause_fcp, fcp;

  if (feature_list == NULL) return;
  for (clause = *feature_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || (clause_fcp = clause->data.ptrvalue) == NULL
      || clause_fcp->featlist == NULL)
    {
      continue;
    }
    if (clause_fcp->featlist->choice == DEFLINE_FEATLIST
      && IsMasterClause (clause_fcp->featlist->data.ptrvalue))
    {
      for (featlist = clause_fcp->featlist->next;
           featlist != NULL;
           featlist = featlist->next)
      {
        if (featlist->choice == DEFLINE_CLAUSEPLUS
          && (fcp = featlist->data.ptrvalue) != NULL)
        {
          fcp->delete_me = TRUE;
        }
      }
    }
    else
    {
      DeleteSubfeatures ( &(clause_fcp->featlist), FALSE);
    }
  }
  if (delete_now) 
  {
    DeleteFeatureClauses (feature_list);
  }
}

static void DeleteOperonAndGeneClusterSubfeatures (
  ValNodePtr PNTR feature_list,
  Boolean         delete_now
)
{
  ValNodePtr       clause, featlist;
  FeatureClausePtr clause_fcp, fcp;

  if (feature_list == NULL) return;
  for (clause = *feature_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS
      || (clause_fcp = clause->data.ptrvalue) == NULL
      || clause_fcp->featlist == NULL)
    {
      continue;
    }
    if (clause_fcp->featlist->choice == DEFLINE_FEATLIST
      && (IsOperon (clause_fcp->featlist->data.ptrvalue) 
          || IsGeneCluster (clause_fcp->featlist->data.ptrvalue)))
    {
      for (featlist = clause_fcp->featlist->next;
           featlist != NULL;
           featlist = featlist->next)
      {
        if (featlist->choice == DEFLINE_CLAUSEPLUS
          && (fcp = featlist->data.ptrvalue) != NULL)
        {
          fcp->delete_me = TRUE;
        }
      }
    }
    else
    {
      DeleteOperonAndGeneClusterSubfeatures ( &(clause_fcp->featlist), FALSE);
    }
  }
  if (delete_now) 
  {
    DeleteFeatureClauses (feature_list);
  }
}

static void RemoveFeats (
  ValNodePtr    list,
  matchFunction itemmatch
)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;
  
  if (list == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST
      && itemmatch (vnp->data.ptrvalue))
    {
      vnp->choice = DEFLINE_REMOVEFEAT;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL)
    {
      RemoveFeats (fcp->featlist, itemmatch);
    }
  }
}

/* A clause is "tall" if it has only one clause at any level */
static Boolean IsClauseTall (
  FeatureClausePtr fcp
)
{
  ValNodePtr featlist;
  Int4       num_clauses;
  FeatureClausePtr subclause;

  num_clauses = 0;
  if (fcp == NULL) return FALSE;
  subclause = NULL;
  if (fcp->featlist == NULL) return FALSE;
  for (featlist = fcp->featlist;
       featlist != NULL;
       featlist = featlist->next)
  {
    if (featlist->choice == DEFLINE_CLAUSEPLUS)
    {
      subclause = featlist->data.ptrvalue;
      if (subclause == NULL || ! IsClauseTall (subclause))
      {
        return FALSE;
      }
      num_clauses ++;
      if (num_clauses > 1) return FALSE;
    }
  }
  if (subclause == NULL || ! subclause->feature_label_data.is_typeword_first)
  {
    return TRUE;
  }
  return FALSE;
}

static void SmashOneTallClause (
  FeatureClausePtr fcp
)
{
  FeatureClausePtr subclause;
  ValNodePtr       featlist;
  ValNodePtr       subclause_featlist;
  ValNodePtr       subclause_firstclause;
  CharPtr          new_description;
  Int4             new_description_len;
  SeqFeatPtr       main_feat;

  if (fcp == NULL || fcp->featlist == NULL) return;

  /* move features up */
  featlist = fcp->featlist;
  if (featlist->choice == DEFLINE_FEATLIST)
  {
    main_feat = fcp->featlist->data.ptrvalue;
  }
  else
  {
    main_feat = NULL;
  }
  
  while (featlist != NULL && featlist->choice != DEFLINE_CLAUSEPLUS)
  {
    featlist = featlist->next;
  }
  if (featlist == NULL) return;
  subclause = featlist->data.ptrvalue;
  if (subclause == NULL) return;
 
  /* move subclause feats to top of list */
  if (subclause->featlist != NULL
    && subclause->featlist->choice == DEFLINE_FEATLIST)
  {
    subclause_featlist = subclause->featlist;
    while (subclause->featlist != NULL
           && subclause->featlist->next != NULL
           && subclause->featlist->next->choice == DEFLINE_FEATLIST)
    {
      subclause->featlist = subclause->featlist->next;
    }
    if (subclause->featlist != NULL)
    {
      subclause_firstclause = subclause->featlist->next;
      subclause->featlist->next = fcp->featlist;
      fcp->featlist = subclause->featlist;
      subclause->featlist = subclause_firstclause;
    }
  }

  /* create new description */
  new_description_len = StringLen (subclause->feature_label_data.description)
                   + StringLen (fcp->feature_label_data.description)
                   + StringLen (fcp->feature_label_data.typeword)
                   + 4;
  new_description = (CharPtr) MemNew (new_description_len);
  if (new_description == NULL) return;
  new_description [0] = 0;
  if ( fcp->feature_label_data.is_typeword_first)
  {
    StringCat (new_description, fcp->feature_label_data.typeword);
    StringCat (new_description, " ");
  }
  StringCat (new_description, fcp->feature_label_data.description);
  if ( ! fcp->feature_label_data.is_typeword_first)
  {
    StringCat (new_description, fcp->feature_label_data.typeword);
  }

  if ( ! IsMobileElement (main_feat)
    && ! IsEndogenousVirusSourceFeature (main_feat))
  {
    StringCat (new_description, ",");
  }
  StringCat (new_description, " ");
  StringCat (new_description, subclause->feature_label_data.description);

  if (fcp->feature_label_data.description != NULL)
  {
    MemFree (fcp->feature_label_data.description);
  }
  fcp->feature_label_data.description = new_description;
 
  /* move interval up */
  if (fcp->interval != NULL)
  {
    MemFree (fcp->interval);
  }
  fcp->interval = subclause->interval;
  subclause->interval = NULL;

  /* move typeword up */
  fcp->feature_label_data.typeword = subclause->feature_label_data.typeword;
  fcp->feature_label_data.is_typeword_first = 
               subclause->feature_label_data.is_typeword_first;
  subclause->feature_label_data.typeword = NULL;
  subclause->delete_me = TRUE;

}


static void SmashTallClauses (
  ValNodePtr PNTR clause_list,
  Boolean         delete_now
)
{
  ValNodePtr clause;
  FeatureClausePtr fcp;

  if (clause_list == NULL) return;
  for (clause = *clause_list; clause != NULL; clause = clause->next)
  {
    if (clause->choice != DEFLINE_CLAUSEPLUS || clause->data.ptrvalue == NULL)
    {
      continue;
    }
    fcp = clause->data.ptrvalue;
    if ( IsClauseTall (fcp))
    {
      SmashOneTallClause (fcp);
    }
    else
    {
      SmashTallClauses (& (fcp->featlist), FALSE);
    }
  }
  if (delete_now) 
  {
    DeleteFeatureClauses (clause_list);
  }
}

static ValNodePtr RemoveAllButLastCDS (
  ValNodePtr list,
  ValNodePtr last_cds
)
{
  ValNodePtr vnp;
  FeatureClausePtr fcp;

  /* now remove all CDSs except the last one */
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST
      && IsCDS (vnp->data.ptrvalue))
    {
      if (last_cds != NULL)
      {
        last_cds->choice = DEFLINE_REMOVEFEAT;
      }
      last_cds = vnp;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL)
    {
      last_cds = RemoveAllButLastCDS (fcp->featlist, last_cds);
    }
  }
  return last_cds;
}

static Boolean OkToConsolidate (
  CharPtr last_desc,
  CharPtr new_desc,
  Boolean last_partial,
  Boolean new_partial,
  FeatureClausePtr last_fcp,
  FeatureClausePtr fcp
)
{
  if (StringCmp (last_desc, new_desc) != 0) return FALSE;
  if (new_partial != last_partial) return FALSE;
  if ( ( fcp->is_alt_spliced && ! last_fcp->is_alt_spliced)
      || (! fcp->is_alt_spliced && last_fcp->is_alt_spliced))
  {
    return FALSE;
  }
  if (fcp->featlist == NULL || last_fcp->featlist == NULL) return FALSE;
  if ( fcp->featlist->choice != DEFLINE_FEATLIST) return FALSE;
  if ( last_fcp->featlist->choice != DEFLINE_FEATLIST) return FALSE;
  if ( (IsCDS (fcp->featlist->data.ptrvalue)
        && ! IsCDS (last_fcp->featlist->data.ptrvalue)
        && ! IsGene (last_fcp->featlist->data.ptrvalue))
      || (! IsCDS (fcp->featlist->data.ptrvalue)
        && ! IsGene (fcp->featlist->data.ptrvalue)
        && IsCDS (last_fcp->featlist->data.ptrvalue)))
  {
    return FALSE;
  }
  if ((IsExon (fcp->featlist->data.ptrvalue) && !IsExon(last_fcp->featlist->data.ptrvalue))
      || (IsExon (last_fcp->featlist->data.ptrvalue) && !IsExon(fcp->featlist->data.ptrvalue))
      || (IsIntron (fcp->featlist->data.ptrvalue) && !IsIntron(last_fcp->featlist->data.ptrvalue))
      || (IsIntron (last_fcp->featlist->data.ptrvalue) && !IsIntron(fcp->featlist->data.ptrvalue)))
  {
    return FALSE;
  }
  return TRUE;
}

static void RemoveRedundantGeneFeatureFromConsolidatedClause (
  FeatureClausePtr fcp
)
{
  ValNodePtr featlist, prevfeat, tmpfeat;
  SeqFeatPtr feat1, feat2;

  prevfeat = NULL;
  featlist = fcp->featlist;
  while ( featlist != NULL
         && featlist->choice == DEFLINE_FEATLIST
         && featlist->next != NULL
         && featlist->next->choice == DEFLINE_FEATLIST)
  {
    feat1 = featlist->data.ptrvalue;
    feat2 = featlist->next->data.ptrvalue;
    if (feat1 == NULL || feat2 == NULL) return;
    if (IsGene (feat1) && ! IsGene (feat2))
    {
      if (prevfeat == NULL)
      {
        fcp->featlist = featlist->next;
        featlist->next = NULL;
        FreeListElement (featlist);
        featlist = fcp->featlist->next;
      }
      else
      {
        prevfeat->next = featlist->next;
        featlist->next = NULL;
        FreeListElement (featlist);
        featlist = prevfeat->next;
      }
    }
    else if ( !IsGene (feat1) && IsGene (feat2))
    {
      tmpfeat = featlist->next;
      featlist->next = tmpfeat->next;
      tmpfeat->next = NULL;
      FreeListElement (tmpfeat);
    }
    else
    {
      featlist = featlist->next;
    }
  }
}

static void PluralizeConsolidatedClauseDescription (
  FeatureClausePtr fcp
)
{
  CharPtr new_desc;

  /* prevent crash */
  if (fcp == NULL || fcp->feature_label_data.description == NULL) return;

  /* don't pluralize tRNA names */
  if (StringNCmp (fcp->feature_label_data.description, "tRNA-", 5) ==0) return;

  new_desc = MemNew (StringLen (fcp->feature_label_data.description) + 2);
  if (new_desc == NULL) return;

  StringCpy (new_desc, fcp->feature_label_data.description);
  StringCat (new_desc, "s");
  MemFree (fcp->feature_label_data.description);
  fcp->feature_label_data.description = new_desc;
}

static void ConsolidateClauses (
  ValNodePtr PNTR list,
  BioseqPtr  bsp,
  Uint1      biomol,
  Boolean    delete_now,
  DeflineFeatureRequestListPtr rp)
{
  ValNodePtr       vnp;
  FeatureClausePtr fcp;
  FeatureClausePtr last_cds_fcp;
  CharPtr          last_desc = NULL, new_desc;
  Boolean          last_partial, new_partial, partial3, partial5;
  SeqLocPtr        new_loc;

  if (list == NULL || *list == NULL) return;
  last_cds_fcp = NULL;
  for (vnp = *list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != DEFLINE_CLAUSEPLUS
      || (fcp = vnp->data.ptrvalue) == NULL
      || fcp->featlist == NULL
      || fcp->featlist->choice != DEFLINE_FEATLIST)
    {
      continue;
    }

    ConsolidateClauses (&(fcp->featlist), bsp, biomol, FALSE, rp);

    if (last_cds_fcp == NULL)
    {
      last_cds_fcp = fcp;
      if (fcp->feature_label_data.description == NULL)
      {
        last_desc = GetGeneProtDescription (fcp, bsp, rp);
      }
      else
      {
        last_desc = StringSave (fcp->feature_label_data.description);
      }
      CheckSeqLocForPartial (fcp->slp, &partial5, &partial3);
      if (partial5 || partial3) 
      {
        last_partial = TRUE;
      }
      else
      {
        last_partial = FALSE;
      }
    }
    else 
    {
      if (fcp->feature_label_data.description == NULL)
      {
        new_desc = GetGeneProtDescription (fcp, bsp, rp);
      }
      else
      {
        new_desc = StringSave (fcp->feature_label_data.description);
      }
      CheckSeqLocForPartial (fcp->slp, &partial5, &partial3);
      if (partial5 || partial3) 
      {
        new_partial = TRUE;
      }
      else
      {
        new_partial = FALSE;
      }
      if ( OkToConsolidate (last_desc, new_desc,
                            last_partial, new_partial,
                            last_cds_fcp, fcp))
      {
        /* two clauses have identical descriptions - combine them */
        MoveSubclauses (last_cds_fcp, fcp);
        RemoveRedundantGeneFeatureFromConsolidatedClause (last_cds_fcp);
        fcp->featlist = NULL;
        fcp->delete_me = TRUE;
        new_loc = SeqLocMerge (bsp, last_cds_fcp->slp, fcp->slp,
                                         FALSE, TRUE, FALSE);
        last_cds_fcp->slp = SeqLocFree (last_cds_fcp->slp);
        last_cds_fcp->slp = new_loc;
        /* if we have two clauses that are really identical instead of
         * just sharing a "prefix", make the description plural
         */
        if (StringCmp (last_cds_fcp->interval, fcp->interval) == 0)
        {
          last_cds_fcp->make_plural = TRUE;
/*          PluralizeConsolidatedClauseDescription (last_cds_fcp); */
        }

        /* Recalculate the interval */
        if (last_cds_fcp->interval != NULL)
        {
          MemFree (last_cds_fcp->interval);
        }
        last_cds_fcp->interval =
                  GetGenericInterval (last_cds_fcp, biomol, bsp, rp);
        MemFree (new_desc);
      }
      else
      {
        MemFree (last_desc);
        last_desc = new_desc;
        last_cds_fcp = fcp;
        last_partial = new_partial;
      }
    }  
  }   
  last_desc = MemFree (last_desc);
  if (delete_now) 
  {
    DeleteFeatureClauses (list);
  }
} 

static void CountUnknownGenes (
  ValNodePtr PNTR clause_list,
  BioseqPtr       bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp, new_fcp;
  ValNodePtr vnp, new_vnp;
  CharPtr gene_name;
  Int4 num_unknown_genes;
  
  num_unknown_genes = 0;
  vnp = *clause_list;
  new_vnp = NULL;
  new_fcp = NULL;
  while (vnp != NULL)
  {
    if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (fcp = vnp->data.ptrvalue) != NULL
      && ! fcp->is_unknown
      && fcp->grp == NULL) 
    {
      CountUnknownGenes (&(fcp->featlist), bsp, rp);
      gene_name = GetGeneProtDescription (fcp, bsp, rp);
      if (StringCmp (gene_name, "unknown") == 0
        && fcp->featlist != NULL
        && fcp->featlist->choice == DEFLINE_FEATLIST)
      {
        if (new_fcp == NULL)
        {
          new_vnp = ValNodeNew (*clause_list);
          if (new_vnp == NULL) return;
          new_fcp = NewFeatureClause (fcp->featlist->data.ptrvalue, 
                                      bsp, rp);
          new_fcp->is_unknown = TRUE;
          new_vnp->choice = DEFLINE_CLAUSEPLUS;
          new_vnp->data.ptrvalue = new_fcp;
        }
        else
        {
          new_vnp = ValNodeNew (new_fcp->featlist);
          if (new_vnp == NULL) return;
          new_vnp->choice = DEFLINE_FEATLIST;
          new_vnp->data.ptrvalue = fcp->featlist->data.ptrvalue;
        }  
        num_unknown_genes ++;
        fcp->delete_me = TRUE;
      }
      gene_name = MemFree (gene_name);
    }
    vnp = vnp->next;
  }
   
  if (num_unknown_genes > 0)
  {
    DeleteFeatureClauses (clause_list);
    if (num_unknown_genes > 1)
    {
      new_fcp->feature_label_data.typeword = StringSave ("genes");
    }
  }
}

NLM_EXTERN void ReplaceDefinitionLine (
  SeqEntryPtr sep,
  CharPtr defline
)
{
  ValNodePtr ttl;
  if (sep == NULL || defline == NULL) return;

  ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  if (ttl == NULL)
    ttl = CreateNewDescriptor (sep, Seq_descr_title);
  if (ttl != NULL) {
    MemFree (ttl->data.ptrvalue);
    ttl->data.ptrvalue = defline;
    defline = NULL;
  }
  MemFree (defline);
}

FeatureClausePtr NewFeatureClause 
( SeqFeatPtr sfp,
  BioseqPtr  bsp,
  DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  Boolean          partial5, partial3;

  fcp = (FeatureClausePtr) MemNew (sizeof (FeatureClauseData));
  if (fcp == NULL) return NULL;

  fcp->feature_label_data.typeword = NULL;
  fcp->feature_label_data.description = NULL;
  fcp->feature_label_data.productname = NULL;
  fcp->feature_label_data.pluralizable = FALSE;
  fcp->feature_label_data.is_typeword_first = FALSE;
  fcp->allelename = NULL;
  fcp->interval = NULL;
  fcp->featlist = NULL;
  fcp->delete_me = FALSE;
  fcp->clause_info_only = FALSE;
  fcp->make_plural = FALSE;
  fcp->is_unknown = FALSE;
  fcp->grp = NULL;
  if (sfp == NULL) return fcp;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  fcp->slp = SeqLocMerge (bsp, sfp->location, NULL,
                                 FALSE, TRUE, FALSE);
  SetSeqLocPartial (fcp->slp, partial5, partial3);
  
  if (sfp->data.choice == SEQFEAT_GENE)
  {
    fcp->grp = sfp->data.value.ptrvalue;
  }
  else
  { 
    fcp->grp = SeqMgrGetGeneXref (sfp);
  }
  if (( IsCDS (sfp) || IsExon (sfp) || IsNoncodingProductFeat (sfp))
    && StringStr (sfp->comment, "alternatively spliced") != NULL) 
  {
    fcp->is_alt_spliced = TRUE;
  }
  else
  {
    fcp->is_alt_spliced = FALSE;
  }
  if (IsCDS (sfp))
  {
    fcp->feature_label_data.productname = GetProductName (sfp, bsp, rp);
  }
  fcp->featlist = ValNodeNew (NULL);
  if (fcp->featlist == NULL)
  {
    MemFree (fcp);
    return NULL;
  }

  fcp->featlist->data.ptrvalue = sfp;
  fcp->featlist->choice = DEFLINE_FEATLIST;
  
  return fcp;
}

static ValNodePtr GetFeatureList (BioseqPtr bsp, DeflineFeatureRequestListPtr rp)
{
  ValNodePtr        head, vnp;
  SeqFeatPtr        sfp;
  FeatureClausePtr  fcp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL) return NULL;

  /* get list of all recognized features */
  head = NULL;
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL)
  {
    if (IsRecognizedFeature (sfp))
    {
      fcp = NewFeatureClause (sfp, bsp, rp);
      if (fcp == NULL) return NULL;
      fcp->numivals = fcontext.numivals;
      fcp->ivals = fcontext.ivals;
      vnp = ValNodeNew (head);
      if (head == NULL) head = vnp;
      if (vnp == NULL) return NULL;
      vnp->data.ptrvalue = fcp;
      vnp->choice = DEFLINE_CLAUSEPLUS;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
  return head;
}

static void ExtractSegmentClauses (
  ValNodePtr segment_features,
  ValNodePtr parent_features,
  ValNodePtr PNTR segment_clauses
);

static Boolean FeatureIsOnSegment (
  SeqFeatPtr sfp, 
  ValNodePtr segment_features
)
{
  ValNodePtr vnp, featclause;
  FeatureClausePtr fcp;

  for (vnp = segment_features; vnp != NULL; vnp = vnp->next)
  {
    fcp = vnp->data.ptrvalue;
    if (fcp != NULL)
    {
      for (featclause = fcp->featlist;
           featclause != NULL;
           featclause = featclause->next)
      {
        if (featclause->data.ptrvalue == sfp) return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean FeatureClauseIsOnSegment (
  FeatureClausePtr fcp,
  ValNodePtr segment_features
)
{
  ValNodePtr vnp;

  if (fcp == NULL || fcp->featlist == NULL) return FALSE;
  for (vnp = fcp->featlist; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST 
      && FeatureIsOnSegment (vnp->data.ptrvalue, segment_features))
    {
      return TRUE;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS)
    {
      if (FeatureClauseIsOnSegment (vnp->data.ptrvalue, segment_features))
      {
        return TRUE;
      }     
    }
  }
  return FALSE;
}

static FeatureClausePtr CopyMatchingClauses (
  FeatureClausePtr fcp,
  ValNodePtr segment_features
)
{
  FeatureClausePtr new_fcp, copy_clause;
  ValNodePtr       featlist, new_feat;
  Boolean          found_feat_on_segment;
  Boolean          partial5, partial3;

  new_fcp = (FeatureClausePtr) MemNew (sizeof (FeatureClauseData));
  if (new_fcp == NULL) return NULL;
  new_fcp->feature_label_data.pluralizable = 
    fcp->feature_label_data.pluralizable;
  new_fcp->feature_label_data.is_typeword_first = 
    fcp->feature_label_data.is_typeword_first;
  new_fcp->feature_label_data.typeword =
    StringSave (fcp->feature_label_data.typeword);
  new_fcp->feature_label_data.description = 
    StringSave (fcp->feature_label_data.description);
  new_fcp->feature_label_data.productname =
    StringSave (fcp->feature_label_data.productname);
  new_fcp->allelename = StringSave (fcp->allelename);
  new_fcp->interval = StringSave (fcp->interval);
  new_fcp->is_alt_spliced = fcp->is_alt_spliced;
  CheckSeqLocForPartial (fcp->slp, &partial5, &partial3);
  new_fcp->slp = (SeqLocPtr) AsnIoMemCopy (fcp->slp, (AsnReadFunc) SeqLocAsnRead,
                                           (AsnWriteFunc) SeqLocAsnWrite);
  SetSeqLocPartial (new_fcp->slp, partial5, partial3);
  new_fcp->grp = fcp->grp;
  new_fcp->has_mrna = fcp->has_mrna;
  new_fcp->delete_me = FALSE;
  new_fcp->clause_info_only = fcp->clause_info_only;
  new_fcp->featlist = NULL;
  found_feat_on_segment = FALSE;
  for (featlist = fcp->featlist; featlist != NULL; featlist = featlist->next)
  {
    new_feat = NULL;
    if (featlist->choice == DEFLINE_FEATLIST)
    {
      new_feat = ValNodeNew (new_fcp->featlist);
      if (new_feat == NULL) return NULL;
      new_feat->data.ptrvalue = featlist->data.ptrvalue;
      new_feat->choice = DEFLINE_FEATLIST;

      /* some portions of the clause are present for product and gene info */
      /* if they aren't actually on the segment */
      if ( segment_features == NULL
        || FeatureIsOnSegment (new_feat->data.ptrvalue, segment_features))
      {
        found_feat_on_segment = TRUE;
      }
    }
    else if (featlist->choice == DEFLINE_CLAUSEPLUS)
    {
      copy_clause = featlist->data.ptrvalue;
      if ( segment_features == NULL
        || FeatureClauseIsOnSegment ( copy_clause, segment_features))
      {
        new_feat = ValNodeNew (new_fcp->featlist);
        if (new_feat == NULL) return NULL;
        new_feat->data.ptrvalue = CopyMatchingClauses ( copy_clause,
                                                        segment_features);
        new_feat->choice = DEFLINE_CLAUSEPLUS;
      }
    }
    if (new_feat != NULL && new_fcp->featlist == NULL) 
    {
      new_fcp->featlist = new_feat;
    }
  }
  if (found_feat_on_segment)
  {
    new_fcp->clause_info_only = FALSE;
  }
  else
  {
    new_fcp->clause_info_only = TRUE;
  }
  return new_fcp; 
}

static void CopyFeatureList (
  ValNodePtr match_features,
  ValNodePtr parent_features,
  ValNodePtr PNTR new_list
)
{
  ValNodePtr vnp, addvnp;

  for (vnp = parent_features;
       vnp != NULL;
       vnp = vnp->next)
  {
    if (vnp->choice == DEFLINE_FEATLIST 
      && (match_features == NULL
        || FeatureIsOnSegment ( vnp->data.ptrvalue, match_features)))
    {
      addvnp = ValNodeNew (*new_list);
      if (addvnp == NULL) return;
      addvnp->data.ptrvalue = vnp->data.ptrvalue;
      addvnp->choice = DEFLINE_FEATLIST;
      if (*new_list == NULL) *new_list = addvnp;
    }
    else if (vnp->choice == DEFLINE_CLAUSEPLUS
      && (match_features == NULL
        || FeatureClauseIsOnSegment ( vnp->data.ptrvalue, match_features)))
    {
      addvnp = ValNodeNew (*new_list);
      if (addvnp == NULL) return;
      addvnp->data.ptrvalue = CopyMatchingClauses ( vnp->data.ptrvalue,
                                                    match_features);
      addvnp->choice = DEFLINE_CLAUSEPLUS;
      if (*new_list == NULL) *new_list = addvnp;
    }
  }
    
}

static void ExtractSegmentClauses (
  ValNodePtr segment_features,
  ValNodePtr parent_features,
  ValNodePtr PNTR segment_clauses
)
{
  CopyFeatureList (segment_features, parent_features, segment_clauses);
}

typedef struct segmentdeflinedata {
  BioseqPtr  parent_bsp;
  ValNodePtr parent_feature_list;
  Uint1      molecule_type;
  DeflineFeatureRequestList PNTR feature_requests;
  ModifierCombinationPtr m;
  ModifierItemLocalPtr modList;
  OrganismDescriptionModifiersPtr odmp;
  Int2 product_flag;
} SegmentDefLineData, PNTR SegmentDefLinePtr;

typedef struct segmentdeflinefeatureclausedata {
  BioseqPtr  parent_bsp;
  ValNodePtr parent_feature_list;
  Uint1      molecule_type;
  DeflineFeatureRequestList PNTR feature_requests;
  Int2            product_flag;
  Boolean         alternate_splice_flag;
  Boolean         gene_cluster_opp_strand;
  ValNodePtr PNTR list;
} SegmentDefLineFeatureClauseData, PNTR SegmentDefLineFeatureClausePtr;

typedef struct deflinefeatclause {
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  CharPtr     clauselist;
} DefLineFeatClauseData, PNTR DefLineFeatClausePtr;

NLM_EXTERN void DefLineFeatClauseListFree (ValNodePtr vnp)
{
  DefLineFeatClausePtr deflist;

  if (vnp == NULL) return;
  DefLineFeatClauseListFree (vnp->next);
  vnp->next = NULL;
  deflist = vnp->data.ptrvalue;
  if (deflist != NULL)
  {
    MemFree (deflist->clauselist);
    MemFree (deflist);
  }
  ValNodeFree (vnp);
}


static Boolean IntervalIntersectsIvals 
(Int2    numivals,
 Int4Ptr ivals,
 SeqMgrSegmentContextPtr context)
{
  Int2 idx;
  Int4 start, stop;

  if (numivals == 0 || ivals == NULL || context == NULL) return FALSE;

  for (idx = 0; idx < numivals; idx ++) {
    start = ivals [idx * 2];
    stop = ivals [idx * 2 + 1];
    if ( start <= context->cumOffset + context->to - context->from
         && stop >= context->cumOffset)
    {
      return TRUE;
    }
  }
  return FALSE;
}


/* if there are no features at all on this segment, select the genes that 
 * traverse the segment.
 */
static ValNodePtr GrabTraversingGenes 
(ValNodePtr              parent_feature_list,
 SeqMgrSegmentContextPtr context,
 BioseqPtr               parent_bsp,
 DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr  fcp, new_fcp;
  ValNodePtr        clause;
  SeqFeatPtr        sfp;
  ValNodePtr        segment_feature_list;
  ValNodePtr        vnp;

  segment_feature_list = NULL;
  for (clause = parent_feature_list;
       clause != NULL;
       clause = clause->next)
  {
    fcp = clause->data.ptrvalue;

    if (fcp != NULL  &&  fcp->featlist != NULL
        &&  (sfp = fcp->featlist->data.ptrvalue) != NULL
        &&  sfp->idx.subtype == FEATDEF_GENE
        &&  fcp->ivals != NULL && fcp->numivals > 0)
    {
      if (IntervalIntersectsIvals (fcp->numivals, fcp->ivals, context)) {
        new_fcp = NewFeatureClause (fcp->featlist->data.ptrvalue, parent_bsp,
                                    rp);
        if (new_fcp == NULL) return FALSE;
        vnp = ValNodeNew (segment_feature_list);
        if (vnp == NULL) return FALSE;
        vnp->data.ptrvalue = new_fcp;
        vnp->choice = DEFLINE_CLAUSEPLUS;
        if (segment_feature_list == NULL) segment_feature_list = vnp;
      }
    } 
  }
  return segment_feature_list;
}


static CharPtr BuildFeatureClauses (
  BioseqPtr   bsp,
  Uint1       molecule_type,
  SeqEntryPtr sep,
  ValNodePtr  PNTR feature_list,
  Boolean     isSegment,
  ValNodePtr  PNTR seg_feature_list,
  Int2        product_flag,
  Boolean     alternate_splice_flag,
  Boolean     gene_cluster_opp_strand,
  DeflineFeatureRequestList PNTR feature_requests
);

static Boolean LIBCALLBACK GetFeatureClauseForSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context)
{
  SegmentDefLineFeatureClausePtr sdlp;
  ValNodePtr        clause, tmp_parent_list;
  FeatureClausePtr  fcp, new_fcp;
  Int2              idx;
  Int4              start, stop;
  ValNodePtr        segment_feature_list, vnp;
  SeqIdPtr          sip;
  BioseqPtr         bsp;
  Uint2             entityID;
  SeqLocPtr         loc;
  DefLineFeatClausePtr deflist;

  if (slp == NULL || context == NULL) return FALSE;
  sdlp = (SegmentDefLineFeatureClausePtr) context->userdata;

  sip = SeqLocId (slp);
  
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqFind (sip);

  if (bsp == NULL) return TRUE;


  segment_feature_list = NULL;
  for (clause = sdlp->parent_feature_list;
       clause != NULL;
       clause = clause->next)
  {
    fcp = clause->data.ptrvalue;

    if (fcp != NULL && fcp->ivals != NULL && fcp->numivals > 0)
    {
      idx = (fcp->numivals - 1) * 2;
      start = fcp->ivals [idx];
      stop = fcp->ivals [idx + 1];
      if ( stop <= context->cumOffset + context->to - context->from
           && stop >= context->cumOffset)
      {
        new_fcp = NewFeatureClause (fcp->featlist->data.ptrvalue,
                                    sdlp->parent_bsp, 
                                    sdlp->feature_requests);
        if (new_fcp == NULL) return FALSE;
        vnp = ValNodeNew (segment_feature_list);
        if (vnp == NULL) return FALSE;
        vnp->data.ptrvalue = new_fcp;
        vnp->choice = DEFLINE_CLAUSEPLUS;
        if (segment_feature_list == NULL) segment_feature_list = vnp;
      }
    } 
  }

  if (segment_feature_list == NULL) {
    segment_feature_list = GrabTraversingGenes (sdlp->parent_feature_list,
                                                context, sdlp->parent_bsp,
                                                sdlp->feature_requests);
  }

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  tmp_parent_list = NULL;
  CopyFeatureList (NULL, sdlp->parent_feature_list, &tmp_parent_list);

  deflist = (DefLineFeatClausePtr) MemNew (sizeof (DefLineFeatClauseData));
  if (deflist == NULL) return TRUE;
  deflist->sep = SeqMgrGetSeqEntryForData (bsp);
  deflist->bsp = bsp;
  deflist->clauselist = BuildFeatureClauses (sdlp->parent_bsp,
                            sdlp->molecule_type,
                            SeqMgrGetSeqEntryForData (bsp),
                            &tmp_parent_list,
                            TRUE,
                            &segment_feature_list,
                            sdlp->product_flag,
                            sdlp->alternate_splice_flag,
                            sdlp->gene_cluster_opp_strand,
                            sdlp->feature_requests);
  vnp = ValNodeNew (*(sdlp->list));
  if (vnp == NULL) return TRUE;
  if (*(sdlp->list) == NULL) *(sdlp->list) = vnp;
  vnp->data.ptrvalue = deflist;

  FreeListElement (tmp_parent_list);
  FreeListElement (segment_feature_list);
  DeleteMarkedObjects (entityID, 0, NULL);
  return TRUE;
}

static Boolean HasAnyPromoters (BioseqPtr bsp)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_promoter, &fcontext);
  if (sfp == NULL) {
    return FALSE;
  } else {
    return TRUE;
  }
}

static void AddFakePromoterClause (ValNodePtr PNTR feature_list, BioseqPtr bsp, DeflineFeatureRequestListPtr rp)
{
  FeatureClausePtr fcp;
  SeqFeatPtr       sfp = NULL;
  
  /* create fake promoter */
  sfp = CreateNewFeature (SeqMgrGetSeqEntryForData (bsp), NULL,
                          SEQFEAT_IMP, NULL);

  sfp->location = SeqLocIntNew (0, bsp->length - 1, Seq_strand_plus, SeqIdDup (bsp->id));
  sfp->data.choice = SEQFEAT_IMP;
  sfp->idx.subtype = FEATDEF_promoter;
  /* mark promoter for deletion */
  sfp->idx.deleteme = TRUE;
  
  fcp = NewFeatureClause (sfp, bsp, rp);
  if (fcp == NULL) return;
  fcp->numivals = 1;
  fcp->ivals = (Int4Ptr) MemNew (sizeof (Int4) * 2);
  fcp->ivals[0] = 0;
  fcp->ivals[1] = bsp->length - 1;
  ValNodeAddPointer (feature_list, DEFLINE_CLAUSEPLUS, fcp);
  
}

static Boolean IsInGenProdSet (BioseqPtr bsp)
{
  BioseqSetPtr bssp;
  if (bsp == NULL || bsp->idx.parentptr == NULL || bsp->idx.parenttype != OBJ_BIOSEQSET) {
    return FALSE;
  }
  bssp = (BioseqSetPtr) bsp->idx.parentptr;
  if (bssp->_class != BioseqseqSet_class_nuc_prot || bssp->idx.parentptr == NULL || bsp->idx.parenttype != OBJ_BIOSEQSET) {
    return FALSE;
  }
  bssp = bssp->idx.parentptr;
  if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    return TRUE;
  } else {
    return FALSE;
  }
}


NLM_EXTERN CharPtr BuildNonFeatureListClause (BioseqPtr bsp, DefLineType feature_list_type)
{
  CharPtr      str = NULL;
  BioSourcePtr biop;
  SeqDescrPtr  sdp;
  SeqMgrDescContext context;
  Char         ending_str [200];
  CharPtr      mol_name;
  MolInfoPtr   molinfo;

  ending_str [0] = 0;

  if (feature_list_type == DEFLINE_COMPLETE_SEQUENCE
      || feature_list_type == DEFLINE_PARTIAL_SEQUENCE
      || feature_list_type == DEFLINE_COMPLETE_GENOME
      || feature_list_type == DEFLINE_PARTIAL_GENOME)
  {
    biop = GetBiopForBsp (bsp);
    if (biop != NULL)
    {
      switch (biop->genome) {
        case GENOME_macronuclear :
          sprintf (ending_str, "macronuclear");
          break;
        case GENOME_nucleomorph :
          sprintf (ending_str, "nucleomorph");
          break;
        case GENOME_mitochondrion :
          sprintf (ending_str, "mitochondrion");
          break;
        case GENOME_apicoplast :
        case GENOME_chloroplast :
        case GENOME_chromoplast :
        case GENOME_kinetoplast :
        case GENOME_plastid :
        case GENOME_cyanelle :
        case GENOME_leucoplast :
        case GENOME_proplastid :
        case GENOME_hydrogenosome :
        case GENOME_chromatophore :
          sprintf (ending_str, "%s", organelleByGenome [biop->genome]);
          break;
      }
    }
  }

  if (feature_list_type == DEFLINE_COMPLETE_SEQUENCE)
  {
    StringCat (ending_str, ", complete sequence.");
  }
  else if (feature_list_type == DEFLINE_PARTIAL_SEQUENCE)
  {
    StringCat (ending_str, ", partial sequence.");
  }
  else if (feature_list_type == DEFLINE_COMPLETE_GENOME)
  {
    StringCat (ending_str, ", complete genome.");
  }
  else if (feature_list_type == DEFLINE_PARTIAL_GENOME)
  {
    StringCat (ending_str, ", partial genome.");
  }
  else if (feature_list_type == DEFLINE_SEQUENCE) 
  {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
    if (sdp == NULL || (molinfo = sdp->data.ptrvalue) == NULL) 
    {
      mol_name = NULL;
    } 
    else 
    {
      mol_name = BiomolNameFromBiomol (molinfo->biomol);
    }
    if (mol_name == NULL) 
    {
      sprintf (ending_str + StringLen (ending_str), " sequence");
    } 
    else 
    {
      sprintf (ending_str + StringLen (ending_str), "%s sequence", mol_name);
    }
  }

  str = StringSave (ending_str);

  return str;
}


/* NOTE: under some circumstances this function will create features that
 * are marked for deletion, so DeleteMarkedObjects should always be called
 * at some later point.
 */
static CharPtr BuildFeatureClauses (
  BioseqPtr   bsp,
  Uint1       molecule_type,
  SeqEntryPtr sep,
  ValNodePtr  PNTR feature_list,
  Boolean     isSegment,
  ValNodePtr  PNTR seg_feature_list,
  Int2        product_flag,
  Boolean     alternate_splice_flag,
  Boolean     gene_cluster_opp_strand,
  DeflineFeatureRequestList PNTR feature_requests
)
{
  ValNodePtr   strings = NULL;
  ValNodePtr   clause;
  CharPtr      str = NULL;
  ValNodePtr   tmp_feat_list;

  if ((feature_requests->feature_list_type == DEFLINE_USE_FEATURES
       || (IsmRNASequence(bsp) && IsInGenProdSet(bsp)))
      && (! isSegment || (seg_feature_list != NULL && *seg_feature_list != NULL)))
  {
    /* remove features that indexer has chosen to suppress before they are grouped
     * with other features or used to determine loneliness etc.
     */
    RemoveSuppressedFeatures (feature_list, feature_requests->suppressed_feature_list);
  
    GroupmRNAs (feature_list, bsp, feature_requests);

    /* genes are added to other clauses */
    GroupGenes (feature_list, feature_requests->suppress_locus_tags);

    if (! feature_requests->suppress_alt_splice_phrase)
    {
      /* find alt-spliced CDSs */
      FindAltSplices (*feature_list, bsp, feature_requests);
    }

    GroupAltSplicedExons (feature_list, bsp, TRUE);
    
    if (!isSegment)
    {
       /* group CDSs that have the same name and are under the same gene together */
      GroupSegmentedCDSs (feature_list, bsp, TRUE, feature_requests);
    }

    /* per Susan's request, if promoters are requested and no promoters are found, add a promoter */
    if (feature_requests->keep_items[RemovablePromoter] 
        && feature_requests->add_fake_promoters
        && !HasAnyPromoters (bsp)) {
      AddFakePromoterClause (feature_list, bsp, feature_requests);
    }

    /* now group clauses */
    GroupAllClauses ( feature_list, gene_cluster_opp_strand, bsp );

    ExpandAltSplicedExons (*feature_list, bsp, feature_requests);

    FindGeneProducts (*feature_list, bsp, feature_requests);

    if (seg_feature_list != NULL && *seg_feature_list != NULL)
    {
      tmp_feat_list = NULL; 
      ExtractSegmentClauses ( *seg_feature_list, *feature_list, &tmp_feat_list);
      FreeListElement (*feature_list);
      DeleteMarkedObjects (bsp->idx.entityID, 0, NULL);
      *feature_list = tmp_feat_list;
    }
   
    /* remove exons and other unwanted features */
    RemoveUnwantedFeatures (feature_list, bsp, isSegment, feature_requests);

    RemoveGenesMentionedElsewhere (feature_list, *feature_list, TRUE,
                                   feature_requests->suppress_locus_tags);

    if (feature_requests->remove_subfeatures)
    {
      DeleteSubfeatures (feature_list, TRUE);
    }

    DeleteOperonAndGeneClusterSubfeatures (feature_list, TRUE);

    CountUnknownGenes (feature_list, bsp, feature_requests);

    if (feature_requests->misc_feat_parse_rule == 1)
    {
      RenameMiscFeats (*feature_list, molecule_type);
    }
    else
    {
      RemoveUnwantedMiscFeats (feature_list, TRUE);
    }

    ReplaceRNAClauses (feature_list, bsp, feature_requests);

    /* take any exons on the minus strand */
    /* and reverse their order within the clause */
    ReverseClauses (feature_list, IsExonOrIntron);

    RenameExonSequences ( feature_list, bsp, TRUE);

    LabelClauses (*feature_list, molecule_type, bsp, 
                  feature_requests);

    /* parse lists of tRNA and intergenic spacer clauses in misc_feat notes */
    /* need to do this after LabelClauses, since LabelClauses labels intergenic
     * spacers with more relaxed restrictions.  The labels from LabelClauses
     * for intergenic spacers are the default values.
     */
    ReplaceIntergenicSpacerClauses (feature_list, bsp, feature_requests);

    ConsolidateClauses (feature_list, bsp, molecule_type, TRUE,
                        feature_requests);

    /* this allows genes to be listed together even if they are from */
    /* separate sequences */
/*    SmashTallClauses (feature_list, TRUE); */

    clause = *feature_list;
    ListClauses (clause, &strings, TRUE, FALSE, feature_requests->suppress_allele);

    AutoDef_AddEnding (clause, &strings, bsp, 
                       product_flag, alternate_splice_flag);
    str = MergeValNodeStrings (strings, FALSE);
      ValNodeFreeData (strings);
  } else {
    str = BuildNonFeatureListClause(bsp, feature_requests->feature_list_type);
  }
  
  return str;
}

/* This function looks at the product names for the CDSs on the Bioseq,
 * and sets the flag for the "nuclear genes for X products" ending
 * based on the contents of the CDS products. */
static Int2 GetProductFlagFromCDSProductNames (BioseqPtr bsp)
{
  SeqMgrFeatContext context;
  SeqFeatPtr        cds = NULL;
  Int2              product_flag;
  Int2              i;
  CharPtr           found;
  Char              ch_before, ch_after;

  product_flag = 0;
  for (cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &context);
       cds != NULL && product_flag == 0;
       cds = cds->next)
  {
    for (i = 1; organelleByPopup[i] != NULL && product_flag == 0; i++)
    {
      found = StringStr (context.label, organelleByPopup[i]);
      if (found != NULL)
      {
        if (found == context.label) {
          ch_before = ' ';
        } else {
          ch_before = *(found - 1);
        }
        ch_after = *(found + StringLen (organelleByPopup[i]));
        if (ch_before == ' ' && (ch_after == 0 || ch_after == ' '))
        {
          product_flag = i;
        }
      }
    }
  }

  return product_flag;
}


static void BuildFeatClauseListForSegSet (
  BioseqPtr bsp,
  Uint2 entityID,
  DeflineFeatureRequestList PNTR feature_requests,
  Int2 product_flag,
  Boolean alternate_splice_flag,
  Boolean gene_cluster_opp_strand,
  ValNodePtr PNTR list
)
{
  SegmentDefLineFeatureClauseData sdld;
  DefLineFeatClausePtr deflist;
  ValNodePtr    vnp;

  /* get default product flag if necessary */
  if (product_flag == -1 || product_flag == DEFAULT_ORGANELLE_CLAUSE) {
    product_flag = GetProductFlagFromCDSProductNames (bsp);
  }

  sdld.parent_bsp = bsp;
  sdld.molecule_type = GetMoleculeType (bsp, entityID);
  sdld.parent_feature_list = GetFeatureList (bsp, 
                                              feature_requests);

  sdld.feature_requests =  feature_requests;
  sdld.product_flag = product_flag;
  sdld.alternate_splice_flag = alternate_splice_flag;
  sdld.gene_cluster_opp_strand = gene_cluster_opp_strand;
  sdld.list = list;
  SeqMgrExploreSegments (bsp, (Pointer) &sdld, GetFeatureClauseForSeg);
  deflist = (DefLineFeatClausePtr) MemNew (sizeof (DefLineFeatClauseData));
  if (deflist == NULL) return;
  deflist->sep = SeqMgrGetSeqEntryForData (bsp),
  deflist->bsp = bsp;

  deflist->clauselist = BuildFeatureClauses (bsp,
                        sdld.molecule_type,
                        SeqMgrGetSeqEntryForData (bsp),
                        &sdld.parent_feature_list,
                        FALSE,
                        NULL,
                        product_flag,
                        alternate_splice_flag,
                        gene_cluster_opp_strand,
                        feature_requests);
  vnp = ValNodeNew (*list);
  if (vnp == NULL) return;
  if (*list == NULL) *list = vnp;
  vnp->data.ptrvalue = deflist;
  FreeListElement (sdld.parent_feature_list);
}


static Boolean Is5SList (ValNodePtr feature_list)
{
  FeatureClausePtr  fcp;
  Boolean           is_5s_list = TRUE;
  SeqFeatPtr        sfp;
  RnaRefPtr         rrp;

  if (feature_list == NULL || feature_list->next == NULL) {
    return FALSE;
  }

  while (feature_list != NULL && is_5s_list) {
    if (feature_list->choice != DEFLINE_CLAUSEPLUS) {
      is_5s_list = FALSE;
    } else if ((fcp = (FeatureClausePtr) feature_list->data.ptrvalue) == NULL) {
      is_5s_list = FALSE;
    } else if (fcp->featlist->choice != DEFLINE_FEATLIST 
               || (sfp = (SeqFeatPtr) fcp->featlist->data.ptrvalue) == NULL) {
      is_5s_list = FALSE;
    } else if (sfp->idx.subtype == FEATDEF_rRNA) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) {
        is_5s_list = FALSE;
      } else if (rrp->ext.choice != 1) {
        is_5s_list = FALSE;
      } else if (StringCmp (rrp->ext.value.ptrvalue, "5S ribosomal RNA") != 0) {
        is_5s_list = FALSE;
      }
    } else if (sfp->idx.subtype == FEATDEF_misc_feature) {
      if (StringCmp (sfp->comment, "nontranscribed spacer") != 0
          && StringCmp (sfp->comment, "contains 5S ribosomal RNA and nontranscribed spacer") != 0) {
        is_5s_list = FALSE;
      }
    } else {
      is_5s_list = FALSE;
    }
    feature_list = feature_list->next;
  }
  return is_5s_list;
}


static void BuildOneFeatClauseList (
  SeqEntryPtr sep,
  Uint2 entityID,
  DeflineFeatureRequestList PNTR feature_requests,
  Int2 product_flag,
  Boolean alternate_splice_flag,
  Boolean gene_cluster_opp_strand,
  ValNodePtr PNTR list
)
{
  BioseqPtr    bsp;
  ValNodePtr    head;
  Uint1      molecule_type;
  SeqEntryPtr   nsep;
  DefLineFeatClausePtr deflist;
  ValNodePtr    vnp;

  nsep = FindNucSeqEntry (sep);
  if (nsep != NULL)
  {
    bsp = (BioseqPtr) nsep->data.ptrvalue;
    if ( SpecialHandlingForSpecialTechniques (bsp)) 
    {
      return;
    }
    if (bsp != NULL && bsp->repr == Seq_repr_seg &&
      bsp->seq_ext != NULL && bsp->seq_ext_type == 1) 
    {
      BuildFeatClauseListForSegSet (bsp, entityID, feature_requests, product_flag, alternate_splice_flag, gene_cluster_opp_strand, list);
    }
  }

  if (nsep != NULL && nsep != sep)
    sep = nsep;


  if (! IS_Bioseq (sep)) return;

  /* get list of all recognized features */
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if ( SpecialHandlingForSpecialTechniques (bsp)) 
  {
    return;
  }
  molecule_type = GetMoleculeType (bsp, entityID);
  head = GetFeatureList (bsp, feature_requests);

  /* get default product flag if necessary */
  if (product_flag == -1 || product_flag == DEFAULT_ORGANELLE_CLAUSE) {
    product_flag = GetProductFlagFromCDSProductNames (bsp);
  }

  deflist = (DefLineFeatClausePtr) MemNew (sizeof (DefLineFeatClauseData));
  if (deflist == NULL) return;
  deflist->sep = SeqMgrGetSeqEntryForData (bsp),
  deflist->bsp = bsp;
  if (Is5SList(head)) {
    deflist->clauselist = StringSave ("5S ribosomal RNA gene region");
  } else {
    deflist->clauselist = BuildFeatureClauses (bsp,
                                               molecule_type,
                                               SeqMgrGetSeqEntryForData (bsp),
                                               &head,
                                               FALSE,
                                               NULL,
                                               product_flag,
                                               alternate_splice_flag,
                                               gene_cluster_opp_strand,
                                               feature_requests);
  }
    
  vnp = ValNodeNew (*list);
  if (vnp == NULL) return;
  if (*list == NULL) *list = vnp;
  vnp->data.ptrvalue = deflist;
  FreeListElement (head);
}


static void RecurseForBuildingFeatClauseLists(
  SeqEntryPtr sep,
  Uint2 entityID,
  DeflineFeatureRequestList PNTR feature_requests,
  Int2 product_flag,
  Boolean alternate_splice_flag,
  Boolean gene_cluster_opp_strand,
  ValNodePtr PNTR list
)
{
  BioseqSetPtr    bssp;

  if ( IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if ( bssp->_class == 7 || IsPopPhyEtcSet (bssp->_class)
        || bssp->_class == BioseqseqSet_class_gen_prod_set
        || bssp->_class == BioseqseqSet_class_not_set)
    {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        RecurseForBuildingFeatClauseLists (sep, entityID, feature_requests,
                                    product_flag, alternate_splice_flag,
                                    gene_cluster_opp_strand, list);
      }
      return;
    }
  }
  BuildOneFeatClauseList (sep, entityID, feature_requests, product_flag, alternate_splice_flag, gene_cluster_opp_strand, list);
}


NLM_EXTERN void BuildDefLineFeatClauseList (
  SeqEntryPtr sep,
  Uint2 entityID,
  DeflineFeatureRequestList PNTR feature_requests,
  Int2 product_flag,
  Boolean alternate_splice_flag,
  Boolean gene_cluster_opp_strand,
  ValNodePtr PNTR list
)
{
  RecurseForBuildingFeatClauseLists (sep, entityID, feature_requests, product_flag, alternate_splice_flag, gene_cluster_opp_strand, list);
  DeleteMarkedObjects (entityID, 0, NULL);
}

static Boolean IdenticalExceptForPartialComplete (CharPtr str1, CharPtr str2)
{
    CharPtr cp, word_in_first, word_in_second;
    Int4    first_len, second_len, compare_len;
    
    if (StringHasNoText (str1) && StringHasNoText (str2)) {
        return TRUE;
    } else if (StringHasNoText (str1) || StringHasNoText (str2)) {
        return FALSE;
    }
    
    word_in_first = StringISearch (str1, "partial");
    cp = StringISearch (str1, "complete");
    if (word_in_first == NULL || (cp != NULL && word_in_first > cp)) {
        word_in_first = cp;
        first_len = 8;
    } else {
        first_len = 7;
    }
    
    word_in_second = StringISearch (str2, "partial");
    cp = StringISearch (str2, "complete");
    if (word_in_second == NULL || (cp != NULL && word_in_second > cp)) {
        word_in_second = cp;
        second_len = 8;
    } else {
        second_len = 7;
    }
    
    if (word_in_first == NULL && word_in_second == NULL) {
        if (StringCmp (str1, str2) == 0) {
            return TRUE;
        } else {
            return FALSE;
        }
    } else if (word_in_first == NULL || word_in_second == NULL) {
        return FALSE;
    } else if ((compare_len = word_in_first - str1) != word_in_second - str2) {
        return FALSE;
    } else if (StringNCmp (str1, str2, compare_len) != 0) {
        return FALSE;
    } else {
        return IdenticalExceptForPartialComplete (word_in_first + first_len, word_in_second + second_len);
    }    
}


static CharPtr GetTaxnameForBsp (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  BioSourcePtr      biop;
  CharPtr           taxname = NULL;

  if (bsp != NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
    if (sdp != NULL && (biop = sdp->data.ptrvalue) != NULL
      && biop->org != NULL) {
      taxname = biop->org->taxname;
      }
  }
  return taxname;
}


NLM_EXTERN Boolean AreFeatureClausesUnique (ValNodePtr list)
{
  ValNodePtr vnp1, vnp2;
  DefLineFeatClausePtr deflist1, deflist2;
  CharPtr              taxname1;

  for (vnp1 = list; vnp1 != NULL && vnp1->next != NULL; vnp1 = vnp1->next)
  {
    deflist1 = vnp1->data.ptrvalue;
    if (deflist1 == NULL || deflist1->clauselist == NULL) return FALSE;
    taxname1 = GetTaxnameForBsp (deflist1->bsp);
    
    for (vnp2 = vnp1->next; vnp2 != NULL; vnp2 = vnp2->next)
    {
      deflist2 = vnp2->data.ptrvalue;
      if (deflist2 == NULL || deflist2->clauselist == NULL
        || (StringCmp (taxname1, GetTaxnameForBsp (deflist2->bsp)) == 0 
            && IdenticalExceptForPartialComplete (deflist1->clauselist, deflist2->clauselist)))
      {
        return FALSE;
      }
    }
  }
  return TRUE;
}


NLM_EXTERN CharPtr GetKeywordPrefix (SeqEntryPtr sep)
{
  ValNodePtr vnp;
  GBBlockPtr gbp;
  MolInfoPtr mip;
  
  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_genbank, NULL);
  if (vnp == NULL) {
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL && mip->tech == MI_TECH_tsa) {
        return "TSA: ";
      }
    }
    } else {
      gbp = (GBBlockPtr) vnp->data.ptrvalue;
      if (gbp != NULL)
      {
        for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next)
        {
          if (StringICmp((CharPtr)vnp->data.ptrvalue, "TPA:inferential") == 0)
          {
            return "TPA_inf: ";
          }
          else if (StringICmp((CharPtr)vnp->data.ptrvalue, "TPA:experimental") == 0)
          {
            return "TPA_exp: ";
          }
        }
      }
    }
    return "";
}


NLM_EXTERN CharPtr 
BuildOneDefinitionLine 
(SeqEntryPtr sep, 
 BioseqPtr   bsp,
 CharPtr  featclause_list, 
  ModifierItemLocalPtr modList,
  ValNodePtr modifier_indices,
  OrganismDescriptionModifiersPtr odmp)
{
  CharPtr    org_desc, tmp_str, keyword_prefix;

  keyword_prefix = GetKeywordPrefix (sep);
  
  org_desc = GetOrganismDescription (bsp, modList, modifier_indices, odmp);
  tmp_str = (CharPtr) MemNew (StringLen (keyword_prefix) 
                              + StringLen (org_desc) 
                              + StringLen (featclause_list) + 2);
  if (tmp_str == NULL) return NULL;
  tmp_str [0] = 0;
  if (keyword_prefix != NULL)
  {
    StringCat (tmp_str, keyword_prefix);
  }
  StringCat (tmp_str, org_desc);
  MemFree (org_desc);
  if (featclause_list != NULL
    && featclause_list [0] != ','
    && featclause_list [0] != '.'
    && featclause_list [0] != 0)
  {
    StringCat (tmp_str, " ");
  }
  StringCat (tmp_str, featclause_list);
  tmp_str [0] = toupper (tmp_str [0]);

  return tmp_str;
}


typedef struct deflineclauseoptions {
  DeflineFeatureRequestList feature_requests;
  Int2 product_flag;
  Boolean alternate_splice_flag;
  Boolean gene_cluster_opp_strand;
} DefLineClauseOptions, PNTR DefLineClauseOptionsPtr;

static DefLineClauseOptionsPtr DefLineClauseOptionsNew (void)
{
  DefLineClauseOptionsPtr clause_options;

  clause_options = MemNew (sizeof (DefLineClauseOptions));
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->product_flag = DEFAULT_ORGANELLE_CLAUSE;
  clause_options->alternate_splice_flag = FALSE;
  clause_options->gene_cluster_opp_strand = FALSE;

  return clause_options;
}


static DefLineClauseOptionsPtr DefLineClauseOptionsFree (DefLineClauseOptionsPtr clause_options)
{
  if (clause_options != NULL) {
    clause_options = MemFree (clause_options);
  }
  return clause_options;
}


typedef void (*Nlm_SetFeatureRequestsProc) PROTO ((DefLineClauseOptionsPtr));

static void DefaultClauseOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->product_flag = DEFAULT_ORGANELLE_CLAUSE;
  clause_options->alternate_splice_flag = FALSE;
  clause_options->gene_cluster_opp_strand = FALSE;
}


static void CompleteSequenceClauseOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->feature_requests.feature_list_type = DEFLINE_COMPLETE_SEQUENCE;
  clause_options->product_flag = DEFAULT_ORGANELLE_CLAUSE;
  clause_options->alternate_splice_flag = FALSE;
  clause_options->gene_cluster_opp_strand = FALSE;
}


static void CompleteGenomeClauseOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->feature_requests.feature_list_type = DEFLINE_COMPLETE_GENOME;
  clause_options->product_flag = DEFAULT_ORGANELLE_CLAUSE;
  clause_options->alternate_splice_flag = FALSE;
  clause_options->gene_cluster_opp_strand = FALSE;
}


static void SequenceClauseOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->feature_requests.feature_list_type = DEFLINE_SEQUENCE;
  clause_options->product_flag = DEFAULT_ORGANELLE_CLAUSE;
  clause_options->alternate_splice_flag = FALSE;
  clause_options->gene_cluster_opp_strand = FALSE;
}


static void MiscFeatNonCodingOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  
  clause_options->feature_requests.keep_items[RemovableNoncodingProductFeat] = TRUE;
}


static void MiscFeatSemicolonOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  
  clause_options->feature_requests.keep_items[RemovableNoncodingProductFeat] = TRUE;
  clause_options->feature_requests.misc_feat_parse_rule = 1;
 
}


static void MitochondrialProductClauseOptions (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->product_flag = 2;
  clause_options->alternate_splice_flag = FALSE;
  clause_options->gene_cluster_opp_strand = FALSE;
}


static void RequestPromoterAndExon (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->feature_requests.keep_items[RemovableExon] = TRUE;
  clause_options->feature_requests.keep_items[RemovablePromoter] = TRUE;
}


static void RequestIntronAndExon (DefLineClauseOptionsPtr clause_options)
{
  InitFeatureRequests (&(clause_options->feature_requests));
  clause_options->feature_requests.keep_items[RemovableExon] = TRUE;
  clause_options->feature_requests.keep_items[RemovableIntron] = TRUE;
}


static Nlm_SetFeatureRequestsProc ClauseOptionSetList[] = {
  DefaultClauseOptions,
  RequestPromoterAndExon,
  RequestIntronAndExon,
  CompleteSequenceClauseOptions,
  CompleteGenomeClauseOptions,
  SequenceClauseOptions,
  MiscFeatNonCodingOptions,
  MiscFeatSemicolonOptions,
  MitochondrialProductClauseOptions,
  NULL };

typedef Boolean (*Nlm_SetOrgModifiersProc) PROTO ((OrganismDescriptionModifiersPtr, ValNodePtr PNTR, ModifierItemLocalPtr));


static void DefaultOrgOptions (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list)
{
  odmp->use_modifiers = TRUE;
  odmp->use_labels = TRUE;
}


static Boolean SubstituteMod (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available, Int4 mod1, Int4 mod2) 
{
  ValNodePtr vnp;
  Boolean has_one = FALSE;
  Boolean has_two = FALSE;
  Boolean rval = FALSE;

  if (mod_list == NULL || *mod_list == NULL || available == NULL || !available[mod2].any_present) {
    return FALSE;
  }
  DefaultOrgOptions(odmp, mod_list);

  for (vnp = *mod_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->data.intvalue == mod1) {
      has_one = TRUE;
    } else if (vnp->data.intvalue == mod2) {
      has_two = TRUE;
    }
  }

  if (has_one && !has_two) {
    for (vnp = *mod_list; vnp != NULL; vnp = vnp->next) {
      if (vnp->data.intvalue == mod1) {
        vnp->data.intvalue = mod2;
        rval = TRUE;
      }
    }
  }
  return rval;
}


static Boolean IsolateInsteadOfClone (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available) 
{
  return SubstituteMod (odmp, mod_list, available, DEFLINE_POS_Clone, DEFLINE_POS_Isolate);
} 


static Boolean  HaplotypeInsteadOfVoucher (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available) 
{
  return SubstituteMod (odmp, mod_list, available, DEFLINE_POS_Specimen_voucher, DEFLINE_POS_Haplotype);
} 


static Boolean  VoucherInsteadOfIsolate (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available) 
{
  return SubstituteMod (odmp, mod_list, available, DEFLINE_POS_Isolate, DEFLINE_POS_Specimen_voucher);
} 


static Boolean UseNone (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  if (mod_list == NULL || *mod_list == NULL) {
    return FALSE;
  }

  DefaultOrgOptions(odmp, mod_list);
  odmp->use_modifiers = FALSE;
  return TRUE;
}


static Boolean AddOneAdjustmentMod (ValNodePtr PNTR mod_list, ModifierItemLocalPtr available, Int4 specific)
{
  Boolean already_has_sv = FALSE;
  Boolean rval = FALSE;
  ValNodePtr vnp, vnp_new;

  if (mod_list == NULL || available == NULL || !available[specific].any_present) {
    return FALSE;
  }

  for (vnp = *mod_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->data.intvalue == specific) {
      already_has_sv = TRUE;
    }
  }
  if (!already_has_sv) {
    vnp_new = ValNodeNew (NULL);
    vnp_new->choice = 0;
    vnp_new->data.intvalue = specific;
    vnp_new->next = *mod_list;
    *mod_list = vnp_new;
    rval = TRUE;
  }
  return rval;
}


static Boolean UseOneSpecific (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available, Int4 specific)
{
  DefaultOrgOptions(odmp, mod_list);
  return AddOneAdjustmentMod (mod_list, available, specific);
}


static Boolean UseStrain (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  return UseOneSpecific (odmp, mod_list, available, DEFLINE_POS_Strain);
}


static Boolean UseSpecimenVoucher (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  return UseOneSpecific (odmp, mod_list, available, DEFLINE_POS_Specimen_voucher);
}


static Boolean UseHaplotype (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  return UseOneSpecific (odmp, mod_list, available, DEFLINE_POS_Haplotype);
}


static Boolean UseAutoDefId (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  size_t index;
  ValNodePtr vnp;
  Boolean found, changed = FALSE;
  ModifierItemLocalPtr cpy;

  if (mod_list == NULL) {
    return FALSE;
  }

  cpy = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (cpy == NULL) return FALSE;
  MemCpy (cpy, available, NumDefLineModifiers () * sizeof (ModifierItemLocalData));

  DefaultOrgOptions(odmp, mod_list);
  SetAutoDefIDModifiers (cpy);
  for (index = 0; index < NumDefLineModifiers(); index++) {
    if (cpy[index].required) {
      found = FALSE;
      if (*mod_list != NULL) {
        for (vnp = *mod_list; vnp != NULL && !found; vnp = vnp->next) {
          if (vnp->data.intvalue == index) {
            found = TRUE;
          }
        }
      }
      if (!found) {
        ValNodeAddInt (mod_list, 0, index);
        changed = TRUE;
      }
    }
  }
  cpy = MemFree (cpy);
  
  return changed;
}


static Boolean DontExcludeSp (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  if (odmp == NULL || !odmp->exclude_sp) {
    return FALSE;
  }
  DefaultOrgOptions(odmp, mod_list);
  odmp->exclude_sp = FALSE;
  return TRUE;
}


static Boolean UseCountryAndIsolate (OrganismDescriptionModifiersPtr odmp, ValNodePtr PNTR mod_list, ModifierItemLocalPtr available)
{
  Boolean add_country, add_isolate;
  Boolean rval = FALSE;

  DefaultOrgOptions(odmp, mod_list);
  add_country =  AddOneAdjustmentMod (mod_list, available, DEFLINE_POS_Country);
  add_isolate = AddOneAdjustmentMod (mod_list, available, DEFLINE_POS_Isolate);
  if (add_country || add_isolate) {
    rval = TRUE;
  }
  return rval;
}


static Nlm_SetOrgModifiersProc OrgModifiersSetList[] = {
  IsolateInsteadOfClone,
  UseNone,
  UseStrain,
  UseSpecimenVoucher,
  UseHaplotype,
  HaplotypeInsteadOfVoucher,
  VoucherInsteadOfIsolate,
  UseAutoDefId,
  DontExcludeSp,
  UseCountryAndIsolate,
  NULL };


typedef Boolean (*Nlm_CompareDeflinesProc) PROTO ((CharPtr, CharPtr));

static Boolean MatchExceptOldProduct (CharPtr old_str, CharPtr new_str)
{
  Int4 old_len, new_len, pattern_len, new_pattern_len, organelle_len, i;
  CharPtr pattern_fmt = "%s %s gene%s for %s product%s.";
  CharPtr other_pattern_fmt = "%s %s gene%s for %s RNA.";
  CharPtr new_pattern_fmt = "; %s.";
  CharPtr pattern, new_pattern;
  Boolean rval = FALSE;
  CharPtr delimiters[] = {";", ",", NULL};
  CharPtr plural[] = {"", "s", NULL};
  Int4 j, k;

  old_len = StringLen (old_str);
  new_len = StringLen (new_str);

  for (i = 2; organelleByGenome[i] != NULL && !rval; i++) {
    if (StringHasNoText (organelleByGenome[i])) {
      continue;
    }
    organelle_len = StringLen (organelleByGenome[i]);
    new_pattern = (CharPtr) MemNew (sizeof (Char) * (StringLen (new_pattern_fmt) + organelle_len));
    sprintf (new_pattern, new_pattern_fmt, organelleByGenome[i]);
    new_pattern_len = StringLen (new_pattern);

    if (new_len > new_pattern_len
        && StringCmp (new_str + new_len - new_pattern_len, new_pattern) == 0) {
      pattern = (CharPtr) MemNew (sizeof (Char) * (StringLen (pattern_fmt) + (2 * organelle_len)));
      for (j = 0; delimiters[j] != NULL; j++) {
        for (k = 0; plural[k] != NULL; k++) {
          /* first pattern */
          sprintf (pattern, pattern_fmt, delimiters[j], organelleByGenome[i], plural[k], organelleByGenome[i], plural[k]);
          pattern_len = StringLen (pattern);     
          if (old_len > pattern_len && new_len > organelle_len + 2
              && StringCmp (old_str + old_len - pattern_len, pattern) == 0) {
            rval = TRUE;
          }
          /* other pattern */
          sprintf (pattern, other_pattern_fmt, delimiters[j], organelleByGenome[i], plural[k], organelleByGenome[i]);
          pattern_len = StringLen (pattern);     
          if (old_len > pattern_len && new_len > organelle_len + 2
              && StringCmp (old_str + old_len - pattern_len, pattern) == 0) {
            rval = TRUE;
          }
        }
      }

      pattern = MemFree (pattern);
    }   
    new_pattern = MemFree (new_pattern);
  }
  return rval;
}


static Boolean NoSequenceForIntron (CharPtr old_str, CharPtr new_str)
{
  Int4 old_len, new_len;

  old_len = StringLen (old_str);
  new_len = StringLen (new_str);
  
  if (new_len < old_len && StringNCmp (old_str, new_str, new_len - 1) == 0
      && (StringCmp (old_str + new_len - 7, "intron, partial sequence.") == 0
          || StringCmp (old_str + new_len - 7, "intron, complete sequence.") == 0)) {
    return TRUE;
  } else {
    return FALSE;
  }

}


static Boolean SkipPair (CharPtr PNTR a, CharPtr PNTR b, CharPtr a_start, CharPtr b_start, CharPtr val1, CharPtr val2)
{
  Int4 len1, len2;
  Char end1, end2;
  Boolean rval = FALSE;

  if (a == NULL || *a == NULL
      || b == NULL || *b == NULL 
      || a_start == NULL || b_start == NULL
      || val1 == NULL || val2 == NULL) {
    return FALSE;
  }
  len1 = StringLen (val1);
  len2 = StringLen (val2);

  end1 = *(val1 + len1 - 1);
  end2 = *(val2 + len2 - 1);

  if (**a == end1 && **b == end2 && *a - a_start > len1 && *b - b_start > len2
            && StringNCmp (*a - len1 + 1, val1, len1) == 0
            && StringNCmp (*b - len2 + 1, val2, len2) == 0) {
    *a -= len1 - 1;
    *b -= len2 - 1;
    rval = TRUE;
  } else if (**a == end2 && **b == end1 && *a - a_start > len2 && *b - b_start > len2
              && StringNCmp (*a - len2 + 1, val2, len2) == 0
              && StringNCmp (*b - len1 + 1, val1, len1) == 0) {
    *a -= len2 - 1;
    *b -= len1 - 1;
    rval = TRUE;
  }
  return rval;
}


static Boolean AdjustForSpace (CharPtr PNTR a, CharPtr PNTR b, CharPtr a_start, CharPtr b_start) 
{
  Boolean rval = FALSE;

  if (a == NULL || *a == NULL 
      || b == NULL || *b == NULL) {
    rval = FALSE;
  } else if (**a == ' ' && *a > a_start && *((*a) - 1) == **b) {
    (*a)--;
    rval = TRUE;
  } else if (**b == ' ' && *b > b_start && *((*b) - 1) == **a) {
    (*b)--;
    rval = TRUE;
  }
  return rval;
}


static Boolean AdjustForCharBeforePhrase (CharPtr PNTR a, CharPtr PNTR b, CharPtr phrase, Char ch)
{
  Boolean rval = FALSE;
  Int4    phrase_len = StringLen (phrase);

  if (a == NULL || *a == NULL 
      || b == NULL || *b == NULL) {
    rval = FALSE;
  } else if (**a == ch && StringNCmp ((*a) + 1, phrase, phrase_len) == 0
             && **b != ch && StringNCmp ((*b) + 1, phrase, phrase_len) == 0) {
    (*a)--;
    rval = TRUE;
  } else if (**b == ch && StringNCmp ((*b) + 1, phrase, phrase_len) == 0
             && **a != ch && StringNCmp ((*a) + 1, phrase, phrase_len) == 0) {
    (*b)--;
    rval = TRUE;
  }
  return rval;
}


static Boolean AdjustForCommaBeforePhrase (CharPtr PNTR a, CharPtr PNTR b, CharPtr phrase)
{
  return AdjustForCharBeforePhrase (a, b, phrase, ',');
}


static Boolean AdjustForCommaBeforeAnd (CharPtr PNTR a, CharPtr PNTR b)
{
  return AdjustForCommaBeforePhrase (a, b, " and ");
}


static Boolean AdjustForSkippableWord (CharPtr PNTR a, CharPtr PNTR b, CharPtr str1, CharPtr str2, CharPtr word)
{
  Int4 len;
  Boolean rval = FALSE;

  if (a == NULL || *a == NULL || b == NULL || *b == NULL) {
    return FALSE;
  }
  len = StringLen (word);
  if (*a - str1 > len && StringNCmp (*a - len + 1, word, len) == 0) {
    *a -= len - 1;
    *b += 1;
    rval = TRUE;
  }
  if (*b - str2 > len && StringNCmp (*b - len + 1, word, len) == 0) {
    *b -= len - 1;
    *a += 1;
    rval = TRUE;
  }
  return rval;
}


static Boolean AdjustForKnownDiffs (CharPtr PNTR a, CharPtr PNTR b, CharPtr str1, CharPtr str2)
{
  Boolean rval = SkipPair (a, b, str1, str2, " pseudogene, partial sequence", " gene, partial cds")
                || SkipPair (a, b, str1, str2, " pseudogene, complete sequence", " gene, partial cds")
                || SkipPair (a, b, str1, str2, " pseudogene, partial sequence", " gene, complete cds")
                || SkipPair (a, b, str1, str2, " pseudogene, complete sequence", " gene, complete cds")
                || SkipPair (a, b, str1, str2, " partial", " complete")
                || SkipPair (a, b, str1, str2, "cds", "sequence")
                || SkipPair (a, b, str1, str2, "gene, partial sequence;", "gene")
                || SkipPair (a, b, str1, str2, "gene, complete sequence;", "gene")
                || SkipPair (a, b, str1, str2, "spacer, partial sequence;", "spacer")
                || SkipPair (a, b, str1, str2, "spacer, complete sequence;", "spacer")
                || SkipPair (a, b, str1, str2, "gene", "genes")
                || AdjustForSkippableWord (a, b, str1, str2, "-like")
                || AdjustForSkippableWord (a, b, str1, str2, " truncated")
                || AdjustForSkippableWord (a, b, str1, str2, " nonfunctional")
                || AdjustForSkippableWord (a, b, str1, str2, " pseudo")
                || AdjustForSkippableWord (a, b, str1, str2, " and")
                || AdjustForSkippableWord (a, b, str1, str2, " partial sequence")
                || AdjustForSkippableWord (a, b, str1, str2, " partial cds")
                || AdjustForSkippableWord (a, b, str1, str2, " complete sequence")
                || AdjustForSkippableWord (a, b, str1, str2, " complete cds")
                || AdjustForSkippableWord (a, b, str1, str2, " gene")
                || AdjustForSkippableWord (a, b, str1, str2, " genes")
                || AdjustForSkippableWord (a, b, str1, str2, " pseudogene")
                || AdjustForSkippableWord (a, b, str1, str2, " pseudogenes")
                || AdjustForSkippableWord (a, b, str1, str2, ";")
                || AdjustForSkippableWord (a, b, str1, str2, ",")
                || AdjustForSpace (a, b, str1, str2);
  return rval;
}


static Boolean AdjustForCommaBeforeGenomicSequence (CharPtr PNTR a, CharPtr PNTR b)
{
  return AdjustForCommaBeforePhrase (a, b, " genomic sequence");
}


static CharPtr defline_skippable_words[] = {
  " transposon",
  " mitochondrial",
  NULL};

static Boolean SkipWord (CharPtr PNTR a, CharPtr a_start, CharPtr PNTR b)
{
  Int4 index, len;
  Boolean rval = FALSE;

  if (a == NULL || *a == NULL || *a <= a_start || b == NULL || *b == NULL) {
     return FALSE;
  }

  for (index = 0; defline_skippable_words[index] != NULL && !rval; index++) {
    len = StringLen (defline_skippable_words[index]);
    if (*a - a_start >= len && StringNCmp (*a - len + 1, defline_skippable_words[index], len) == 0) {
      *a -= len - 1;
      (*b)++;
      rval = TRUE;
    }
  }
  return rval;
}
   

static Boolean MatchWithPhraseExceptions (CharPtr str1, CharPtr str2)
{
  Int4 len_curr, len_new;
  CharPtr a, b;
  CharPtr mitochondrial = "; mitochondrial";
  Int4 len_mito = StringLen (mitochondrial);

  len_curr = StringLen (str1);
  len_new = StringLen (str2);
  b = str2;
  if (len_curr > 0) {
    a = str1 + len_curr - 1;
    b = str2 + len_new - 1;
    if (*a == '.') {
      a--;
    }
    if (*b == '.') {
      b--;
    }
    /* skip trailing "mitochondrial" */
    if (a - str1 > len_mito && StringNCmp (a - len_mito + 1, mitochondrial, len_mito) == 0) {
      a -= len_mito;
    }
    if (b - str2 > len_mito && StringNCmp (b - len_mito + 1, mitochondrial, len_mito) == 0) {
      b -= len_mito;
    }

    while (a >= str1 && b >= str2 
            && (toupper (*a) == toupper(*b)
                || SkipPair (&a, &b, str1, str2, " gene", " genes")
                || SkipPair (&a, &b, str1, str2, "; and", ";")
                || AdjustForCommaBeforeAnd (&a, &b)
                || AdjustForCommaBeforeGenomicSequence (&a, &b)
                || AdjustForSpace (&a, &b, str1, str2)
                || SkipWord (&a, str1, &b)
                || SkipWord (&b, str2, &a))) {
      a--;
      b--;
    }
  }
  if (a + 1 == str1 && b + 1 == str2) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Nlm_CompareDeflinesProc CompareDeflinesList[] = {
  MatchExceptOldProduct,
  NoSequenceForIntron,
  MatchWithPhraseExceptions,
  NULL };


static Boolean DeflinesMatch (CharPtr old_str, CharPtr new_str)
{
  Int4 index;
  Boolean match = FALSE;

  if (StringCmp (old_str, new_str) == 0) {
    return TRUE;
  }

  for (index = 0; CompareDeflinesList[index] != NULL && !match; index++) {
    match = CompareDeflinesList[index](old_str, new_str);
  }
  return match;
}


static void ReplaceOldName (BioseqPtr bsp, CharPtr PNTR old_title)
{
  SeqDescPtr                    sdp;
  SeqMgrDescContext             dcontext;
  BioSourcePtr                  biop;
  OrgModPtr                     omp;

  if (old_title == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL && (biop = (BioSourcePtr)sdp->data.ptrvalue) != NULL && biop->org != NULL && biop->org->orgname != NULL) {
    omp = biop->org->orgname->mod;
    while (omp != NULL) {
      if (omp->subtype == ORGMOD_old_name) {
        FindReplaceString (old_title, omp->subname, biop->org->taxname, TRUE, TRUE);
      }
      omp = omp->next;
    }
  }
}


static void AddMissingPeriod (CharPtr PNTR old_title)
{
  Int4 len;
  CharPtr new_title;

  if (old_title == NULL || *old_title == NULL) {
    return;
  }

  len = StringLen (*old_title);
  if ((*old_title)[len - 1] != '.') {
    new_title = (CharPtr) MemNew (sizeof (Char) * (len + 2));
    sprintf (new_title, "%s.", *old_title);
    *old_title = MemFree (*old_title);
    *old_title = new_title;
  }
}


static void ChangeExonList (CharPtr PNTR old_title)
{
  CharPtr exon_start, first_and, second_and;
  CharPtr src, dst;

  if (old_title == NULL || StringHasNoText (*old_title) 
      || (exon_start = StringStr (*old_title, "exons ")) == NULL
      || (first_and = StringStr (exon_start, " and ")) == NULL
      || (second_and = StringStr (first_and + 1, " and ")) == NULL) {
    return;
  }

  if (StringCmp (second_and, " and partial cds.") != 0 && StringCmp (second_and, " and complete cds.") != 0) {
    return;
  }

  *first_and = ',';
  src = first_and + 4,
  dst = first_and + 1;
  while (*src != 0) {
    *dst = *src;
    dst++;
    src++;
  }
  *dst = 0;
}


NLM_EXTERN DefLineClauseOptionsPtr MakeFeatureRequestsMatchExpectedTitle (BioseqPtr bsp)
{
  SeqEntryPtr                   sep;
  ValNodePtr                    defline_clauses = NULL;
  size_t                        index;
  Int4                          mod_index;
  ValNodePtr                    best_modifier_indices, default_modifier_indices, modifier_indices = NULL, tmp_mod_list;
  ValNodePtr                    vnp;
  ModifierItemLocalPtr          modList;
  OrganismDescriptionModifiers  odmp;
  SeqEntryPtr                   oldscope;
  CharPtr                       clause_list = NULL, attempt = NULL, old_title;
  DefLineFeatClausePtr          df;
  DefLineClauseOptionsPtr       clause_options;
  SeqDescPtr                    sdp;
  SeqMgrDescContext             dcontext;
  Boolean                       found_match = FALSE;
  Boolean                       default_exclude_sp;

  if (bsp == NULL) {
    return NULL;
  }
  sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  if (sep == NULL) return NULL;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return NULL;

  sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_title, &dcontext);
  if (sdp == NULL || (old_title = (CharPtr) sdp->data.ptrvalue) == NULL) {
    return NULL;
  }

  old_title = StringSave (old_title);
  /* replacements for old RNA style */
  FindReplaceString (&old_title, "RNA gene, mitochondrial gene for mitochondrial RNA, partial sequence", "RNA gene, partial sequence; mitochondrial", TRUE, TRUE);
  FindReplaceString (&old_title, "16S ribosomal RNA, complete sequence", "16S ribosomal RNA gene, complete sequence", TRUE, TRUE);
  FindReplaceString (&old_title, "16S ribosomal RNA, partial sequence", "16S ribosomal RNA gene, partial sequence", TRUE, TRUE);
  FindReplaceString (&old_title, "subunit ribosomal RNA, complete sequence", "subunit ribosomal RNA gene, complete sequence", TRUE, TRUE);
  FindReplaceString (&old_title, "subunit ribosomal RNA, partial sequence", "subunit ribosomal RNA gene, partial sequence", TRUE, TRUE);
  FindReplaceString (&old_title, "specimen-voucher", "voucher", TRUE, TRUE);

  ReplaceOldName (bsp, &old_title);
  AddMissingPeriod (&old_title);
  ChangeExonList (&old_title);

  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  InitOrganismDescriptionModifiers (&odmp, sep);
  default_exclude_sp = odmp.exclude_sp;
  odmp.use_modifiers = TRUE;

  oldscope = SeqEntrySetScope (sep);

  clause_options = DefLineClauseOptionsNew();

  best_modifier_indices = FindBestModifiersEx (sep, modList, TRUE);
  default_modifier_indices = GetModifierIndicesFromModList (modList);

  /* start loop here */
  for (index = 0;
       ClauseOptionSetList[index] != NULL && !found_match;
       index++) {
    ClauseOptionSetList[index](clause_options);
    BuildDefLineFeatClauseList (sep, bsp->idx.entityID, 
                                &(clause_options->feature_requests),
                                clause_options->product_flag, 
                                clause_options->alternate_splice_flag, 
                                clause_options->gene_cluster_opp_strand,
                                &defline_clauses);
    if (AreFeatureClausesUnique (defline_clauses)) {
      modifier_indices = ValNodeDupIntList(default_modifier_indices);
    } else {
      modifier_indices = ValNodeDupIntList(best_modifier_indices);
    }
    for (vnp = defline_clauses; vnp != NULL; vnp = vnp->next) {
      df = (DefLineFeatClausePtr) vnp->data.ptrvalue;
      if (df->bsp == bsp) {
        clause_list = df->clauselist;
        break;
      }
    }

    attempt = BuildOneDefinitionLine (sep, bsp, clause_list, modList, modifier_indices, &odmp);
    if (DeflinesMatch(old_title, attempt)) {
      found_match = TRUE;
    }
    attempt = MemFree (attempt);

    if (!found_match) {
      odmp.use_labels = FALSE;
      attempt = BuildOneDefinitionLine (sep, bsp, clause_list, modList, modifier_indices, &odmp);
      if (DeflinesMatch(old_title, attempt)) {
        found_match = TRUE;
      }
      attempt = MemFree (attempt);
      odmp.use_labels = TRUE; 
    }

    if (!found_match) {
      for (mod_index = 0; OrgModifiersSetList[mod_index] != NULL && !found_match; mod_index++) {
        tmp_mod_list = ValNodeDupIntList (modifier_indices);
        if (OrgModifiersSetList[mod_index](&odmp, &tmp_mod_list, modList)) {
          attempt = BuildOneDefinitionLine (sep, bsp, clause_list, modList, tmp_mod_list, &odmp);
          if (DeflinesMatch(old_title, attempt)) {
            found_match = TRUE;
          }
          attempt = MemFree (attempt);
          if (!found_match) {
            odmp.use_labels = FALSE;
            attempt = BuildOneDefinitionLine (sep, bsp, clause_list, modList, tmp_mod_list, &odmp);
            if (DeflinesMatch(old_title, attempt)) {
              found_match = TRUE;
            }
            attempt = MemFree (attempt);
            odmp.use_labels = TRUE; 
          }
        }
        tmp_mod_list = ValNodeFree (tmp_mod_list);
        odmp.exclude_sp = default_exclude_sp;
      }
    }

    DefLineFeatClauseListFree (defline_clauses);
    defline_clauses = NULL;
    modifier_indices = ValNodeFree (modifier_indices);
  }
  /* end loop here */
  best_modifier_indices = ValNodeFree (best_modifier_indices);
  default_modifier_indices = ValNodeFree (default_modifier_indices);

  old_title = MemFree (old_title);

  if (modList != NULL) {
    for (index = 0; index < NumDefLineModifiers (); index++) {
      ValNodeFree (modList [index].values_seen);
    }
    MemFree (modList);
  }

  SeqEntrySetScope (oldscope);

  if (!found_match) {
    clause_options = DefLineClauseOptionsFree (clause_options);
  }

  return clause_options;
}


static Int4 MatchlenForAutodef (CharPtr str1, CharPtr str2)
{
  Int4 len_curr, len_new;
  CharPtr a, b;

  len_curr = StringLen (str1);
  len_new = StringLen (str2);
  b = str2;
  if (len_curr > 0) {
    a = str1 + len_curr - 1;
    b = str2 + len_new - 1;
    if (*a == '.') {
      a--;
    }
    if (*b == '.') {
      b--;
    }
    if (a - str1 > 13 && StringNCmp (a - 12, "; chloroplast", 13) == 0) {
      a-= 13;
    }
    if (b - str2 > 13 && StringNCmp (b - 12, "; chloroplast", 13) == 0) {
      b -= 13;
    }
    while (a >= str1 && b >= str2 
            && (*a == *b || AdjustForKnownDiffs(&a, &b, str1, str2))) {
      a--;
      b--;
    }
  }
  return len_new - (b - str2 + 1);
}


static void FindCommonTitleCallback (BioseqPtr bsp, Pointer data)
{
  Int4 len_curr, len_new;
  CharPtr new_title;
  CharPtr PNTR current_title;
  CharPtr a, b;
  CharPtr src, dst;
  SeqMgrDescContext context;
  SeqDescPtr sdp;

  if (bsp == NULL || ISA_aa (bsp->mol)) {
    return;
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);

  if (sdp == NULL || sdp->choice != Seq_descr_title 
      || (new_title = (CharPtr) sdp->data.ptrvalue) == NULL
      || (current_title = (CharPtr PNTR) data) == NULL) {
    return;
  }

  if (*current_title == NULL) {
    *current_title = StringSave (new_title);
  } else {
    len_curr = StringLen (*current_title);
    len_new = StringLen (new_title);
    if (len_curr > 0) {
      a = (*current_title) + len_curr - 1;
      b = new_title + len_new - 1;
      if (*a == '.') {
        a--;
      }
      if (*b == '.') {
        b--;
      }
      if (a - *current_title > 13 && StringNCmp (a - 12, "; chloroplast", 13) == 0) {
        a-= 13;
      }
      if (b - new_title > 13 && StringNCmp (b - 12, "; chloroplast", 13) == 0) {
        b -= 13;
      }
      while (a >= *current_title && b >= new_title 
            && (*a == *b || AdjustForKnownDiffs(&a, &b, *current_title, new_title))) {
        a--;
        b--;
      }
      /* don't copy part of a word */
      while (*(a + 1) != 0 && !isspace (*(a + 1))) {
        a++;
      }
      if (a + 1 > *current_title) {
        src = a + 1;
        dst = *current_title;
        while (*src != 0) {
          *dst = *src;
          dst++;
          src++;
        }
        *dst = 0;
      }
    }
  }
}


typedef struct verifycommonfeatureclause {
  CharPtr common_clause;
  Boolean is_ok;
} VerifyCommonFeatureClauseData, PNTR VerifyCommonFeatureClausePtr;


static void PrintBioSource (BioSourcePtr biop)
{
  OrgModPtr mod;
  SubSourcePtr ssp;

  printf ("Taxname: %s", biop->org->taxname);
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    printf ("\tOrgMod%d:%s", mod->subtype, mod->subname);
  }
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    printf ("\tSubSource%d:%s", ssp->subtype, ssp->name);
  }
  printf ("\n");
}


static void VerifyCommonFeatureClauseCallback (BioseqPtr bsp, Pointer data)
{
  VerifyCommonFeatureClausePtr v;
  SeqMgrDescContext context;
  SeqDescPtr sdp;
  CharPtr title;
  BioSourcePtr biop;
  Int4 common_clause_len, title_len, last_word_len = 0;
  SourceConstraint source_constraint;
  CharPtr last_word, cp;
  Char id_txt[100];

  if (bsp == NULL || ISA_aa (bsp->mol) 
      || (v = (VerifyCommonFeatureClausePtr) data) == NULL
      || !v->is_ok) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);

  if (sdp == NULL || sdp->choice != Seq_descr_title || (title = (CharPtr) sdp->data.ptrvalue) == NULL) {
    return;
  }
  title_len = StringLen (title);
  common_clause_len = MatchlenForAutodef (v->common_clause, title);

  cp = title + title_len - common_clause_len - 1;
  while (cp > title && !isspace (*cp)) {
    --cp;
    last_word_len++;
  }
  if (*(cp + last_word_len) == ',') {
    last_word_len--;
  }

  /* note - allocate extra space, in case we need to try adding a semicolon */
  last_word = (CharPtr) MemNew (sizeof (Char) * (last_word_len + 2));
  StringNCpy (last_word, cp + 1, last_word_len);
  last_word[last_word_len] = 0;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);

  if (sdp == NULL || sdp->choice != Seq_descr_source || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL) {
    return;
  }

  MemSet (&source_constraint, 0, sizeof (SourceConstraint));
  source_constraint.constraint = StringConstraintNew ();
  source_constraint.constraint->match_text = last_word;
  source_constraint.constraint->match_location = String_location_ends;
  source_constraint.constraint->case_sensitive = FALSE;
  source_constraint.constraint->whole_word = TRUE;
  source_constraint.constraint->not_present = FALSE;

  if (!DoesBiosourceMatchConstraint(biop, &source_constraint)) {
    StringCat (last_word, ";");
    source_constraint.constraint->match_location = String_location_contains;
    if (!DoesBiosourceMatchConstraint(biop, &source_constraint)) {
      /* try country with colon */
      last_word[last_word_len] = ':';
      source_constraint.constraint->match_location = String_location_contains;
      source_constraint.field1 = ValNodeNew (NULL);
      source_constraint.field1->choice = SourceQualChoice_textqual;
      source_constraint.field1->data.intvalue = Source_qual_country;
      if (!DoesBiosourceMatchConstraint(biop, &source_constraint)) {
        v->is_ok = FALSE;
        SeqIdWrite (bsp->id, id_txt, PRINTID_TEXTID_ACC_ONLY, sizeof (id_txt) - 1);
        printf ("For %s, common title fails, last word is %s\n", id_txt, last_word);
        PrintBioSource (biop);
      }
      source_constraint.field1 = ValNodeFree (source_constraint.field1);
    }
  }
  source_constraint.constraint = StringConstraintFree (source_constraint.constraint);
}


static CharPtr GetCommonFeatureClause (SeqEntryPtr sep)
{
  CharPtr common_clause = NULL;
  VerifyCommonFeatureClauseData v;
  CharPtr src, dst;

  VisitBioseqsInSep (sep, &common_clause, FindCommonTitleCallback);
  if (StringHasNoText (common_clause)) {
    common_clause = MemFree (common_clause);
  } else {
    MemSet (&v, 0, sizeof (VerifyCommonFeatureClauseData));
    v.common_clause = common_clause;
    v.is_ok = TRUE;
    VisitBioseqsInSep (sep, &v, VerifyCommonFeatureClauseCallback);
    if (!v.is_ok) {
      if (StringNCmp (common_clause, " chloroplast", 12) == 0) {
        dst = common_clause;
        for (src = common_clause + 12; *src != 0; src++) {
          *dst = *src;
        }
        *dst = 0;
        v.is_ok = TRUE;
        VisitBioseqsInSep (sep, &v, VerifyCommonFeatureClauseCallback);
      }
    }
    if (!v.is_ok) {
      common_clause = MemFree (common_clause);
    }
  }
  return common_clause;
}


NLM_EXTERN void BuildDefinitionLinesFromFeatureClauseLists (
  ValNodePtr list,
  ModifierItemLocalPtr modList,
  ValNodePtr modifier_indices,
  OrganismDescriptionModifiersPtr odmp
)
{
  ValNodePtr vnp;
  DefLineFeatClausePtr defline_featclause;
  CharPtr    tmp_str;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL)
    {
      defline_featclause = vnp->data.ptrvalue;
      tmp_str = BuildOneDefinitionLine (defline_featclause->sep, defline_featclause->bsp,
                                        defline_featclause->clauselist,
                                        modList, modifier_indices, odmp);
      ReplaceDefinitionLine (defline_featclause->sep, tmp_str);
    }
  }
}

NLM_EXTERN void BuildDefLinesFromFeatClauseListsForOneBsp (
  ValNodePtr list,
  ModifierItemLocalPtr modList,
  ValNodePtr modifier_indices,
  OrganismDescriptionModifiersPtr odmp,
  BioseqPtr bsp
)
{
  ValNodePtr vnp;
  DefLineFeatClausePtr defline_featclause;
  CharPtr    tmp_str;

  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue != NULL)
    {
      defline_featclause = vnp->data.ptrvalue;
      if (defline_featclause == NULL) continue;
      if (defline_featclause->bsp != bsp) continue;
      tmp_str = BuildOneDefinitionLine (defline_featclause->sep, defline_featclause->bsp,
                                        defline_featclause->clauselist,
                                        modList, modifier_indices, odmp);
      ReplaceDefinitionLine (defline_featclause->sep, tmp_str);
    }
  }
}


/* This removes redundant titles on nuc-prot sets, which will not be
 * visible in the flat file if all sequences in the nuc-prot set have
 * their own title.
 */
NLM_EXTERN void RemoveNucProtSetTitles (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  this_sep;
  SeqDescrPtr  sdp, prev = NULL, sdp_next;
  
  if (sep == NULL || ! IS_Bioseq_set (sep))
  {
    return;
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  for (this_sep = bssp->seq_set; this_sep != NULL; this_sep = this_sep->next)
  {
    RemoveNucProtSetTitles (this_sep);
  }
  
  if (bssp->_class != BioseqseqSet_class_nuc_prot) 
  {
    return;
  }
  for (sdp = bssp->descr; sdp != NULL; sdp = sdp_next)
  {
    sdp_next = sdp->next;
    if (sdp->choice == Seq_descr_title)
    {
      if (prev == NULL)
      {
        bssp->descr = sdp->next;
      }
      else
      {
        prev->next = sdp->next;
      }
      sdp->next = NULL;
      SeqDescrFree (sdp);
    }
    else
    {
      prev = sdp;
    }
  }
}


static void ProtTitleRemoveProc (BioseqPtr bsp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title && sdp->extended) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}

NLM_EXTERN void RemoveProteinTitles (SeqEntryPtr sep)

{
  Uint2  entityID;

  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, ProtTitleRemoveProc);
  entityID = ObjMgrGetEntityIDForChoice (sep);
  DeleteMarkedObjects (entityID, 0, NULL);
}

static void MRnaTitleRemoveProc (BioseqPtr bsp, Pointer userdata)

{
  MolInfoPtr     mip;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp == NULL) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL || mip->biomol != MOLECULE_TYPE_MRNA) return;

  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title && sdp->extended) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}

NLM_EXTERN void RemoveMRnaTitles (SeqEntryPtr sep)

{
  Uint2  entityID;

  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, MRnaTitleRemoveProc);
  entityID = ObjMgrGetEntityIDForChoice (sep);
  DeleteMarkedObjects (entityID, 0, NULL);
}


typedef struct popsetdefline {
 DeflineFeatureRequestListPtr feature_requests;
 Int2 product_flag;
 Boolean alternate_splice_flag;
 Boolean gene_cluster_opp_strand;
} PopsetDeflineData, PNTR PopsetDeflinePtr;



NLM_EXTERN Boolean GetsDocsumTitle(Uint1 set_class)
{
  if (set_class == BioseqseqSet_class_pop_set
      || set_class == BioseqseqSet_class_phy_set
      || set_class == BioseqseqSet_class_eco_set
      || set_class == BioseqseqSet_class_mut_set) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean HasTitle(SeqDescrPtr descr)
{
  while (descr != NULL) {
    if (descr->choice == Seq_descr_title) {
      return TRUE;
    }
    descr = descr->next;
  }
  return FALSE;
}


NLM_EXTERN void AddPopsetDeflineWithClause (BioseqSetPtr bssp, CharPtr clause)
{
  SeqEntryPtr set_sep;
  OrgRefPtr   org;
  CharPtr     org_desc = NULL, keyword_prefix = NULL, tmp_str = NULL;

  set_sep = SeqMgrGetSeqEntryForData (bssp);
  org = GetCommonOrgRefForSeqEntry (set_sep);
  if (org == NULL || StringHasNoText (org->taxname)) {
    return;
  }
  org_desc = StringSave (org->taxname);
  org = OrgRefFree (org);
    
  keyword_prefix = GetKeywordPrefix (set_sep);
  
  tmp_str = (CharPtr) MemNew (StringLen (keyword_prefix) 
                              + StringLen (org_desc) 
                              + StringLen (clause) + 2);
  if (tmp_str == NULL) return;
  tmp_str [0] = 0;
  if (keyword_prefix != NULL)
  {
    StringCat (tmp_str, keyword_prefix);
  }
  StringCat (tmp_str, org_desc);
  if (clause != NULL
    && clause[0] != ','
    && clause[0] != '.'
    && clause[0] != ' '
    && clause[0] != 0)
  {
    StringCat (tmp_str, " ");
  }
  StringCat (tmp_str, clause);
  tmp_str [0] = toupper (tmp_str [0]);
  ReplaceDefinitionLine (set_sep, tmp_str);
  MemFree (org_desc);
}


static void AddPopsetCallback (BioseqSetPtr bssp, Pointer data)
{
  SeqEntryPtr          set_sep, first_sep;
  OrgRefPtr            org;
  ValNodePtr           defline_clauses = NULL;
  DefLineFeatClausePtr defline_featclause;
  CharPtr              org_desc = NULL, keyword_prefix = NULL, tmp_str = NULL;
  PopsetDeflinePtr     pop;

  if (bssp == NULL 
      || !GetsDocsumTitle(bssp->_class) 
      || HasTitle(bssp->descr)
      || (pop = (PopsetDeflinePtr)data) == NULL) {
    return;
  }

  set_sep = SeqMgrGetSeqEntryForData (bssp);
  org = GetCommonOrgRefForSeqEntry (set_sep);
  if (org == NULL || StringHasNoText (org->taxname)) {
    return;
  }
  org_desc = StringSave (org->taxname);
  org = OrgRefFree (org);

  /* make title here */
  first_sep = bssp->seq_set;

  BuildDefLineFeatClauseList (first_sep, bssp->idx.entityID,
                              pop->feature_requests,
                              pop->product_flag, pop->alternate_splice_flag, pop->gene_cluster_opp_strand,
                              &defline_clauses);

  if (defline_clauses != NULL) {
    defline_featclause = defline_clauses->data.ptrvalue;
    
    keyword_prefix = GetKeywordPrefix (defline_featclause->sep);
    
    tmp_str = (CharPtr) MemNew (StringLen (keyword_prefix) 
                                + StringLen (org_desc) 
                                + StringLen (defline_featclause->clauselist) + 2);
    if (tmp_str == NULL) return;
    tmp_str [0] = 0;
    if (keyword_prefix != NULL)
    {
      StringCat (tmp_str, keyword_prefix);
    }
    StringCat (tmp_str, org_desc);
    if (defline_featclause->clauselist != NULL
      && defline_featclause->clauselist [0] != ','
      && defline_featclause->clauselist [0] != '.'
      && defline_featclause->clauselist [0] != 0)
    {
      StringCat (tmp_str, " ");
    }
    StringCat (tmp_str, defline_featclause->clauselist);
    tmp_str [0] = toupper (tmp_str [0]);
    ReplaceDefinitionLine (set_sep, tmp_str);
  }
  MemFree (org_desc);
}


NLM_EXTERN void AddPopsetTitles 
(SeqEntryPtr sep,
 DeflineFeatureRequestListPtr feature_requests,
 Int2 product_flag,
 Boolean alternate_splice_flag,
 Boolean gene_cluster_opp_strand)
{
  PopsetDeflineData pop;

  pop.feature_requests = feature_requests;
  /* forcibly suppress alleles in popset titles */
  pop.feature_requests->suppress_allele = TRUE;
  pop.product_flag = product_flag;
  pop.alternate_splice_flag = alternate_splice_flag;
  pop.gene_cluster_opp_strand = gene_cluster_opp_strand;

  VisitSetsInSep (sep, &pop, AddPopsetCallback);
}


static void RemovePopsetTitlesCallback(BioseqSetPtr bssp, Pointer data)
{
  SeqDescrPtr   sdp;
  ObjValNodePtr ovp;

  if (bssp == NULL || !GetsDocsumTitle(bssp->_class)) {
    return;
  }
  for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title
        && sdp->extended) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}


NLM_EXTERN void RemovePopsetTitles(SeqEntryPtr sep)
{
  Uint2 entityID;

  VisitSetsInSep (sep, NULL, RemovePopsetTitlesCallback);
  entityID = ObjMgrGetEntityIDForChoice (sep);
  DeleteMarkedObjects (entityID, 0, NULL);
}


NLM_EXTERN void 
AutoDefForSeqEntry 
(SeqEntryPtr sep,
 Uint2 entityID,
 OrganismDescriptionModifiersPtr odmp,
 ModifierItemLocalPtr modList,
 ValNodePtr modifier_indices,
 DeflineFeatureRequestListPtr feature_requests,
 Int2 product_flag,
 Boolean alternate_splice_flag,
 Boolean gene_cluster_opp_strand)
{
  ValNodePtr defline_clauses = NULL;

  if (sep == NULL) return;

  RemoveNucProtSetTitles (sep);
  
  SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID,
                              feature_requests,
                              product_flag, alternate_splice_flag,
                              gene_cluster_opp_strand,
                              &defline_clauses);
                              
/*  dlfp->feature_requests.suppressed_feature_list = ValNodeFree (dlfp->feature_requests.suppressed_feature_list);                               */

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, odmp);
  DefLineFeatClauseListFree (defline_clauses);
  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);

  AddPopsetTitles (sep, feature_requests, product_flag, 
                   alternate_splice_flag, gene_cluster_opp_strand);
}


/* Retro PopSet Title Functions */
static SeqDescPtr BioseqHasTitleOrNucProtSetHasTitle (BioseqPtr bsp)
{
  SeqDescPtr sdp = NULL;
  SeqMgrDescContext context;
  BioseqSetPtr bssp;
  ObjValNodePtr ovn;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &context);
  while (sdp != NULL) {
    if (sdp->extended) {
      ovn = (ObjValNodePtr) sdp;
      if (ovn->idx.parentptr == bsp) {
        return sdp;
      } else if (ovn->idx.parenttype == OBJ_BIOSEQSET
                 && (bssp = (BioseqSetPtr) ovn->idx.parentptr) != NULL
                 && bssp->_class == BioseqseqSet_class_nuc_prot) {
        return sdp;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_title, &context);
  }
  return NULL;
}


static BioseqPtr FindFirstNucBioseqWithTitle (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sep == NULL || sep->data.ptrvalue == NULL ||
      /* sep->choice < 0 || */ sep->choice > 2) return NULL;
  if (IS_Bioseq (sep) ) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (!ISA_aa (bsp->mol) && BioseqHasTitleOrNucProtSetHasTitle(bsp) != NULL) {
      return bsp;
    } else {
      return NULL;
    }
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    bsp = FindFirstNucBioseqWithTitle (sep);
    if (bsp != NULL) return bsp;
  }
  return NULL;
}


static void RetroPopSetAutoDefCallback (BioseqSetPtr bssp, Pointer data)
{
  SeqEntryPtr set_sep;
  BioseqPtr   first_bsp;
  DefLineClauseOptionsPtr clause_options = NULL;
  CharPtr                    common_title = NULL;
  SeqEntryPtr oldscope, sep;
  PopSetRetroStatPtr stat = (PopSetRetroStatPtr) data;

  if (bssp == NULL || !GetsDocsumTitle (bssp->_class) || HasTitle (bssp->descr)) {
    return;
  }

  set_sep = SeqMgrGetSeqEntryForData (bssp);
  first_bsp = FindFirstNucBioseqWithTitle (set_sep);
  if (first_bsp == NULL) {
    if (stat != NULL) {
      stat->uncalculatable++;
    }
    return;
  }
  sep = GetTopSeqEntryForEntityID (bssp->idx.entityID);

  /* infer feature list from first title */
  clause_options = MakeFeatureRequestsMatchExpectedTitle (first_bsp);

  if (clause_options == NULL) {
    common_title = GetCommonFeatureClause (set_sep);
    if (common_title == NULL) {
      if (stat != NULL) {
        stat->uncalculatable++;
      }
    } else {
      AddPopsetDeflineWithClause (bssp, common_title);
      common_title = MemFree (common_title);
      if (stat != NULL) {
        stat->common_title++;
        stat->title_added = TRUE;
      }
    }
  } else {
    oldscope = SeqEntrySetScope (sep);

    AddPopsetTitles (sep, &(clause_options->feature_requests),
                     clause_options->product_flag,
                     clause_options->alternate_splice_flag,
                     clause_options->gene_cluster_opp_strand);

    SeqEntrySetScope (oldscope);
    clause_options = DefLineClauseOptionsFree (clause_options);
    if (stat != NULL) {
      stat->feature_clause++;
      stat->title_added = TRUE;
    }
  }
}


NLM_EXTERN void PopSetAutoDefRetro (SeqEntryPtr sep, PopSetRetroStatPtr stat)

{
  SeqEntryPtr oldscope;

  if (sep == NULL) return;

  oldscope = SeqEntrySetScope (sep);

  VisitSetsInSep (sep, stat, RetroPopSetAutoDefCallback);

  SeqEntrySetScope (oldscope);
}




/* functions for editing seq-locs */
NLM_EXTERN Int4 ExtendSeqLocToEnd (SeqLocPtr slp, BioseqPtr bsp, Boolean end5)
{
  Uint1          strand;
  SeqLocPtr      slp_to_change, slp_index;
  Int4           extent_to_change;
  Int4           start, stop;
  SeqIdPtr       sip;
  Int4           start_diff = 0;
  
  if (slp == NULL || bsp == NULL) return 0;

  slp_to_change = NULL;
  strand = SeqLocStrand (slp);
  switch (slp->choice)
  {
    case SEQLOC_INT:
      slp_to_change = slp;
      break;
    case SEQLOC_MIX:
      case SEQLOC_PACKED_INT:
      sip = SeqLocId (slp);
      if (sip == NULL) return 0; /* can only process if all on one bioseq */
      slp_to_change = NULL;
      if ((strand == Seq_strand_minus && end5)
        || (strand != Seq_strand_minus && !end5))
      {
        extent_to_change = 0;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          stop = GetOffsetInBioseq (slp_index, bsp, SEQLOC_STOP);
          if (stop > extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = stop;
          }
        }
      }
      else
      {
        extent_to_change = bsp->length;
        for (slp_index = (SeqLocPtr)slp->data.ptrvalue; slp_index != NULL; slp_index = slp_index->next)
        {
          start = GetOffsetInBioseq (slp_index, bsp, SEQLOC_START);
          if (start < extent_to_change)
          {
            slp_to_change = slp_index;
            extent_to_change = start;
          }
        }
      }
      break;
  }

  if (slp_to_change != NULL)
  {
    if ((strand == Seq_strand_minus && end5)
      || (strand != Seq_strand_minus && !end5))
    {
      start = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_START);
      stop = bsp->length - 1;
    }
    else
    {
      start = 0;
      stop = GetOffsetInBioseq (slp_to_change, bsp, SEQLOC_STOP);
    }
    if (end5) {
        if (strand == Seq_strand_minus) {
            start_diff = bsp->length - 1 - GetOffsetInBioseq(slp_to_change, bsp, SEQLOC_START);
        } else {
            start_diff = GetOffsetInBioseq(slp_to_change, bsp, SEQLOC_START);
        }
    }
    
    expand_seq_loc (start, stop, strand, slp_to_change);
  }
  return start_diff;
}


/* functions for feature conversion.  shared by sequin5 and macroapi */

NLM_EXTERN Boolean IsBioseqSetInGPS (BioseqSetPtr bssp)
{
  if (bssp == NULL) return FALSE;
  if (bssp->_class == BioseqseqSet_class_gen_prod_set) return TRUE;
  if (bssp->idx.parentptr == NULL || bssp->idx.parenttype != OBJ_BIOSEQSET) return FALSE;
  return IsBioseqSetInGPS ((BioseqSetPtr) bssp->idx.parentptr);
}


NLM_EXTERN Boolean IsBioseqInGPS (BioseqPtr bsp)
{
  if (bsp == NULL || bsp->idx.parentptr == NULL || bsp->idx.parenttype != OBJ_BIOSEQSET)
  {
    return FALSE;
  }
  else
  {
    return IsBioseqSetInGPS ((BioseqSetPtr) bsp->idx.parentptr);
  }
}


NLM_EXTERN Boolean IsFeatInGPS (SeqFeatPtr sfp)
{
  if (sfp == NULL) return FALSE;
  return IsBioseqInGPS (BioseqFindFromSeqLoc (sfp->location));
}


NLM_EXTERN RnaRefPtr RnaRefFromLabel (Uint2 featdef_to, CharPtr label, BoolPtr add_label_to_comment)
{
  RnaRefPtr rrp;
  tRNAPtr   trp = NULL;
  Boolean   just_trna_text;
  Int4      j;

  rrp = RnaRefNew ();
  if (NULL == rrp)
    return NULL;

  *add_label_to_comment = FALSE;

  switch (featdef_to) 
  {
    case FEATDEF_preRNA :
      rrp->type = 1;
      break;
    case FEATDEF_mRNA :
      rrp->type = 2;
      break;
    case FEATDEF_tRNA :
      rrp->type = 3;
      break;
    case FEATDEF_rRNA :
      rrp->type = 4;
      break;
    case FEATDEF_snRNA :
      rrp->type = 8;
      break;
    case FEATDEF_scRNA :
      rrp->type = 8;
      break;
    case FEATDEF_snoRNA :
      rrp->type = 8;
      break;
    case FEATDEF_ncRNA:
      rrp->type = 8;
      break;
    case FEATDEF_tmRNA:
      rrp->type = 9;
      break;
    case FEATDEF_misc_RNA:
      rrp->type = 10;
      break;
    case FEATDEF_otherRNA :
      rrp->type = 255;
      break;
    default :
      break;
  }

  if (featdef_to == FEATDEF_tRNA) 
  {
    trp = (tRNAPtr) MemNew (sizeof (tRNA));
    rrp->ext.choice = 2;
    rrp->ext.value.ptrvalue = (Pointer) trp;
    trp->aa = ParseTRnaString (label, &just_trna_text, NULL, FALSE);
    trp->aatype = 2;
    for (j = 0; j < 6; j++) {
        trp->codon [j] = 255;
    }
    if (!just_trna_text)
    {
      *add_label_to_comment = TRUE;
    }
  } 
  else 
  {
    SetRNARefProductString (rrp, NULL, label, ExistingTextOption_replace_old);
  }
  return rrp;
}


NLM_EXTERN Boolean ConvertProtToProtFunc 
(SeqFeatPtr sfp,
 Uint2      featdef_to)
{
  ProtRefPtr prp;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (NULL == prp)
    return FALSE;

  switch (featdef_to) {
    case FEATDEF_PROT :
      prp->processed = 0;
      break;
    case FEATDEF_preprotein :
      prp->processed = 1;
      break;
    case FEATDEF_mat_peptide_aa :
      prp->processed = 2;
      break;
    case FEATDEF_sig_peptide_aa :
      prp->processed = 3;
      break;
    case FEATDEF_transit_peptide_aa :
      prp->processed = 4;
      break;
    default :
      break;
  }
  return TRUE;
}


NLM_EXTERN void 
ApplyCDSOptionsToFeature
(SeqFeatPtr sfp,
 Boolean remove_mRNA,
 Boolean remove_gene,
 Boolean remove_transcript_id,
 Boolean keep_original)
{
  BioseqPtr         product_bsp;
  SeqFeatPtr        gene, mrna;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL) return;

  if (sfp->product != NULL) {
    product_bsp = BioseqFindFromSeqLoc (sfp->product);
    if (product_bsp != NULL && !keep_original)
    {
      product_bsp->idx.deleteme = TRUE;
    }
    if (!IsFeatInGPS (sfp) || remove_transcript_id)
    {
      sfp->product = SeqLocFree (sfp->product);
    }
  }

  if (remove_gene)
  {
    gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
    if (gene != NULL)
    {
      gene->idx.deleteme = TRUE;
    }
  }
  
  if (remove_mRNA)
  {
    mrna = SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
    if (mrna != NULL)
    {
      mrna->idx.deleteme = TRUE;
    }
  }

}


NLM_EXTERN Boolean 
ConvertCDSToRNA 
(SeqFeatPtr  sfp,
 Uint2       rna_type)
{
  Char                   label [256];
  CharPtr                new_comment;
  Int4                   comment_len = 0;
  Boolean                add_label_to_comment = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;

  FeatDefLabel (sfp, label, sizeof (label), OM_LABEL_CONTENT);

  sfp->data.value.ptrvalue =
      CdRegionFree ((CdRegionPtr) sfp->data.value.ptrvalue);

  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = RnaRefFromLabel (rna_type, label, &add_label_to_comment);
  
  if (add_label_to_comment && StringCmp (label, sfp->comment) != 0)
  {
    if (StringHasNoText (sfp->comment)) 
    {
      new_comment = StringSave (label);
    }
    else
    {
      comment_len = StringLen (sfp->comment) + StringLen (label) + 3;
      new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
      sprintf (new_comment, "%s; %s", sfp->comment, label);
    }
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;
  }
  /* change subtype so that feature will be reindexed */
  sfp->idx.subtype = 0;  
  
  return TRUE;
}


NLM_EXTERN Boolean ConvertGeneToRNA (SeqFeatPtr sfp, Uint2 featdef_to)
{
  Char                   label [256];
  GeneRefPtr grp;
  Boolean    add_label_to_comment = FALSE;
  CharPtr    new_comment;
  Int4       comment_len = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return FALSE;

  FeatDefLabel (sfp, label, sizeof (label), OM_LABEL_CONTENT);

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;

  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = RnaRefFromLabel (featdef_to, label, &add_label_to_comment);
  
  if (add_label_to_comment)
  {
    comment_len += StringLen (label) + 2;
  }
  if (grp != NULL)
  {
    if (!StringHasNoText (grp->locus) && StringCmp (grp->locus, label) != 0)
    {
      comment_len += StringLen (grp->locus) + 2;
    }
    if (!StringHasNoText (grp->allele) && StringCmp (grp->allele, label) != 0)
    {
      comment_len += StringLen (grp->allele) + 2;
    }
    if (!StringHasNoText (grp->desc) && StringCmp (grp->desc, label) != 0)
    {
      comment_len += StringLen (grp->desc) + 2;
    }
    if (!StringHasNoText (grp->maploc) && StringCmp (grp->maploc, label) != 0)
    {
      comment_len += StringLen (grp->maploc) + 2;
    }
    if (!StringHasNoText (grp->locus_tag) && StringCmp (grp->locus_tag, label) != 0)
    {
      comment_len += StringLen (grp->locus_tag) + 2;
    }
  }
  if (comment_len > 0)
  {
    comment_len += StringLen (sfp->comment) + 3;
    new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
    if (!StringHasNoText (sfp->comment))
    {
      StringCat (new_comment, sfp->comment);
      StringCat (new_comment, "; ");
    }
    if (add_label_to_comment)
    {
      StringCat (new_comment, label);
      StringCat (new_comment, "; ");
    }
    /* append unused gene qualifiers */
    if (grp != NULL)
    {
      if (!StringHasNoText (grp->locus) && StringCmp (grp->locus, label) != 0)
      {
        StringCat (new_comment, grp->locus);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->allele) && StringCmp (grp->allele, label) != 0)
      {
        StringCat (new_comment, grp->allele);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->desc) && StringCmp (grp->desc, label) != 0)
      {
        StringCat (new_comment, grp->desc);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->maploc) && StringCmp (grp->maploc, label) != 0)
      {
        StringCat (new_comment, grp->maploc);
        StringCat (new_comment, "; ");
      }
      if (!StringHasNoText (grp->locus_tag) && StringCmp (grp->locus_tag, label) != 0)
      {
        StringCat (new_comment, grp->locus_tag);
        StringCat (new_comment, "; ");
      }
    }
    /* remove last semicolon */
    new_comment[StringLen (new_comment) - 2] = 0;
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;

  }
  
  /* free generef */
  grp = GeneRefFree (grp);      
  
  return TRUE;
}


/* These functions are used for converting features on nucleotide sequences to
 * features on protein sequences */

/* copied from seqport.c, for the benefit of load_fuzz_to_DNA */
static Boolean add_fuzziness_to_loc (SeqLocPtr slp, Boolean less)
{
    IntFuzzPtr ifp;
    SeqIntPtr sint;
    SeqPntPtr spnt;    

    sint = NULL;
    spnt = NULL;

    if(slp->choice == SEQLOC_INT)
        sint = (SeqIntPtr) slp->data.ptrvalue;
    else
    {
        if(slp->choice == SEQLOC_PNT)
            spnt = (SeqPntPtr) slp->data.ptrvalue;
        else
            return FALSE;
    }
    ifp = IntFuzzNew();
    ifp->choice = 4;
    ifp->a = less ? 2 : 1;

    if(spnt != NULL)
        spnt->fuzz = ifp;
    else if (sint != NULL)
    {
        if(less)
            sint->if_from = ifp;
        else
            sint->if_to = ifp;
    }

    return TRUE;
}


/* copied from seqport.c, for the benefit of MYdnaLoc_to_aaLoc */
static Boolean load_fuzz_to_DNA(SeqLocPtr dnaLoc, SeqLocPtr aaLoc, Boolean 
first)
{
    Uint1 strand;
    SeqPntPtr spnt;
    SeqIntPtr sint;
    IntFuzzPtr ifp;
    Boolean load, less;

    load = FALSE;
    strand = SeqLocStrand(aaLoc);
    if(aaLoc->choice == SEQLOC_INT)
    {
        sint = (SeqIntPtr) aaLoc->data.ptrvalue;
        if((first && strand != Seq_strand_minus ) || 
            (!first && strand == Seq_strand_minus))    /*the first 
Seq-loc*/
        {
            ifp = sint->if_from;
            if(ifp && ifp->choice == 4 )
                load = (ifp->a == 2);
        }
        else
        {
            ifp = sint->if_to;
            if(ifp && ifp->choice == 4)
                load = (ifp->a == 1);
        }
    }
    else if(aaLoc->choice == SEQLOC_PNT)
    {
        spnt = (SeqPntPtr) aaLoc->data.ptrvalue;
        ifp = spnt->fuzz;
        if(ifp && ifp->choice == 4)
        {
            if(first)
                load = (ifp->a == 2);
            else
                load = (ifp->a == 1);
        }
    }

    if(load)
    {
        if(SeqLocStrand(dnaLoc) == Seq_strand_minus)
            less = (first == FALSE);
        else
            less = first;
        add_fuzziness_to_loc (dnaLoc, less);
        return TRUE;
    }
    else
        return FALSE;
}    


static SeqLocPtr MYdnaLoc_to_aaLoc(SeqFeatPtr sfp, 
                                   SeqLocPtr location_loc, 
                                   Boolean merge, 
                                   Int4Ptr frame, 
                                   Boolean allowTerminator)
{
    SeqLocPtr aa_loc = NULL, loc;
    CdRegionPtr crp;
    Int4 product_len, end_pos, frame_offset;
    GatherRange gr;
    Int4 a_left = 0, a_right, last_aa = -20, aa_from, aa_to;
    SeqLocPtr slp, slp1, slp2;
    Int2 cmpval;
    SeqIdPtr aa_sip;
    BioseqPtr bsp;

    if ((sfp == NULL) || (location_loc == NULL)) return aa_loc;
    if (sfp->data.choice != 3) return aa_loc;
    if (sfp->product == NULL) return aa_loc;

    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if(crp == NULL) return aa_loc;

    /* each interval of location_loc must be equal to or contained in
     * an interval of sfp->location
     */
    slp1 = SeqLocFindNext (sfp->location, NULL);
    slp2 = SeqLocFindNext (location_loc, NULL);
    while (slp2 != NULL && slp1 != NULL) {
      cmpval = SeqLocCompare (slp2, slp1);
      if (cmpval == SLC_A_IN_B || cmpval == SLC_A_EQ_B) {
        slp2 = SeqLocFindNext (location_loc, slp2);
      } else {
        slp1 = SeqLocFindNext (sfp->location, slp1);
      }
    }
    if (slp1 == NULL) return aa_loc;
      
    aa_sip = SeqLocId(sfp->product);
    if (aa_sip == NULL) return aa_loc;
    bsp = BioseqLockById(aa_sip);
    if (bsp == NULL) return aa_loc;
    end_pos = bsp->length - 1;
    BioseqUnlock(bsp);

    if(crp->frame == 0)
        frame_offset = 0;
    else
        frame_offset = (Int4)crp->frame-1;

    slp = NULL;
    product_len = 0;
    loc = NULL;
    while ((slp = SeqLocFindNext(sfp->location, slp))!=NULL)
    {
       if (SeqLocOffset(location_loc, slp, &gr, 0))
       {
            SeqLocOffset(slp, location_loc, &gr, 0);
        
            a_left = gr.left + product_len - frame_offset;
            a_right = gr.right + product_len - frame_offset;

            aa_from = a_left / 3;
            aa_to = a_right / 3;

            if (aa_from < 0)
                aa_from = 0;
            if (aa_to > end_pos)
                aa_to = end_pos;

            if (merge)
            {
                if (aa_from <= last_aa)  /* overlap due to 
codons */
                    aa_from = last_aa+1;  /* set up to merge 
*/
            }

            if (aa_from <= aa_to || (allowTerminator && aa_from == aa_to + 1))
            {
                if(loc != NULL)
                {
                    if(aa_loc == NULL)
                        load_fuzz_to_DNA(loc, location_loc, TRUE);
                    SeqLocAdd(&aa_loc, loc, merge, FALSE);
                }
                loc = SeqLocIntNew(aa_from, aa_to, 0, aa_sip);
                last_aa = aa_to;
            }
         }

         product_len += SeqLocLen(slp);        
    }

    if(loc != NULL)
    {
        if(aa_loc == NULL)
            load_fuzz_to_DNA(loc, location_loc, TRUE);
        load_fuzz_to_DNA(loc, location_loc, FALSE);
        SeqLocAdd(&aa_loc, loc, merge, FALSE);
    }
    if (frame != NULL)
        *frame = a_left % 3;

    return SeqLocPackage(aa_loc);
}


NLM_EXTERN SeqLocPtr BuildProtLoc (SeqFeatPtr overlapping_cds, SeqLocPtr slp, Int4Ptr frame)
{
  SeqLocPtr tmp_loc, aa_loc = NULL, prot_loc = NULL, last_loc = NULL, next_loc;
  Boolean   partial5, partial3;
  BioseqPtr prot_bsp;
  Boolean   is_ordered = FALSE;
  Boolean   first = TRUE;
  
  prot_bsp = BioseqFindFromSeqLoc (overlapping_cds->product);
  if (prot_bsp == NULL) {
    return NULL;
  }
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  tmp_loc = SeqLocFindNext (slp, NULL);
  while (tmp_loc != NULL) {
    if (tmp_loc->choice == SEQLOC_NULL) {
      is_ordered = TRUE;
    } else {
      if (first) {
        aa_loc = MYdnaLoc_to_aaLoc (overlapping_cds, tmp_loc, FALSE, frame, FALSE);
        first = FALSE;
      } else {
        aa_loc = MYdnaLoc_to_aaLoc (overlapping_cds, tmp_loc, FALSE, NULL, FALSE);
      }
    }
    if (last_loc == NULL) {
      prot_loc = aa_loc;
    } else {
      last_loc->next = aa_loc;
    }
    last_loc = aa_loc;
    tmp_loc = SeqLocFindNext (slp, tmp_loc);
  }
  if (prot_loc != NULL && prot_loc->next != NULL) {
    tmp_loc = NULL;
    for (aa_loc = prot_loc; aa_loc != NULL; aa_loc = next_loc) {
      next_loc = aa_loc->next;
      aa_loc->next = NULL;
      
      last_loc = SeqLocMerge (prot_bsp, tmp_loc, aa_loc, FALSE, TRUE, is_ordered);
      tmp_loc = SeqLocFree (tmp_loc);
      aa_loc = SeqLocFree (aa_loc);
      tmp_loc = last_loc;
      last_loc = NULL;

      aa_loc = next_loc;
    }
    prot_loc = tmp_loc;
  }
  SetSeqLocPartial (prot_loc, partial5, partial3);
  return prot_loc;
}

NLM_EXTERN Boolean ConvertRegionToProtFunc (SeqFeatPtr sfp, Uint2 featdef_to)
{
  BioseqPtr  bsp;
  ProtRefPtr prp;
  SeqFeatPtr cds;
  SeqLocPtr  location;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION)
  {
    return FALSE; 
  }
  
  /* only convert features that are on protein sequences */
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (!ISA_aa (bsp->mol))
  {
    cds = SeqMgrGetOverlappingCDS (sfp->location, NULL);
    if (cds == NULL)
    {
      return FALSE;
    }
    else
    {
      location = BuildProtLoc (cds, sfp->location, NULL);
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = location;
    }
  }
  
  prp = ProtRefNew ();
  if (prp != NULL)
  {
    prp->name = ValNodeNew(NULL);
    if (prp->name != NULL)
    {
      /* use region name for protein name */
      prp->name->data.ptrvalue = sfp->data.value.ptrvalue;
      switch (featdef_to) 
      {
        case FEATDEF_PROT :
          prp->processed = 0;
          break;
        case FEATDEF_preprotein :
          prp->processed = 1;
          break;
        case FEATDEF_mat_peptide_aa :
          prp->processed = 2;
          break;
        case FEATDEF_sig_peptide_aa :
          prp->processed = 3;
          break;
        case FEATDEF_transit_peptide_aa :
          prp->processed = 4;
          break;
        default :
          break;
      }

      sfp->data.value.ptrvalue = prp;
      sfp->data.choice = SEQFEAT_PROT;
    }
  }  
  return TRUE;
}


NLM_EXTERN SeqLocPtr GetProteinLocationForNucleotideFeatureConversion (SeqLocPtr nuc_slp, BoolPtr no_cds)
{
  SeqFeatPtr cds;
  SeqMgrFeatContext cds_context;
  SeqLocPtr  prot_slp;

  cds = SeqMgrGetOverlappingCDS (nuc_slp, &cds_context);
  if (cds == NULL) {
    if (no_cds != NULL) {
      *no_cds = TRUE;
    }
    return NULL;
  } else if (no_cds != NULL) {
    *no_cds = FALSE;
  }

  prot_slp = BuildProtLoc (cds, nuc_slp, NULL);
  return prot_slp;
}



/*---------------------------------------------------------------------*/
/*                                                                     */
/* ConvertImpToProt () - Convert a given import feature to a    */
/*                           protein feature.                          */
/*                                                                     */
/*    Note : Any of the Import feature's gbquals that can be converted */
/*           to protein fields are caught in the automatic cleanup     */
/*           called during reindexing, so they don't need to be        */
/*           converted here.                                           */
/*                                                                     */
/*---------------------------------------------------------------------*/

NLM_EXTERN Boolean ConvertImpToProtFunc 
(SeqFeatPtr  sfp,
 Uint2       featdef_to)
{
  ImpFeatPtr ifp;
  SeqFeatPtr cds;
  SeqLocPtr  slp;
  SeqFeatPtr newSfp;
  Int4       frame;
  ProtRefPtr prp;
  SeqIdPtr   sip;
  BioseqPtr  bsp;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP)
  {
    return FALSE;
  }
  /* Get the Import Feature */

  ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
  if (NULL == ifp)
  {
    return FALSE;
  }

  /* Convert the location to a protein location */
  cds = SeqMgrGetOverlappingCDS (sfp->location, &fcontext);
  if (cds == NULL)
  {
    return FALSE;
  }

  slp = BuildProtLoc (cds, sfp->location, &frame);
  if (slp == NULL)
  {
    return FALSE;
  }

  /* Create a new generic feature */

  sip = SeqLocId (cds->product);
  if (sip == NULL)
  {
    slp = SeqLocFree (slp);
    return FALSE;
  }

  bsp = BioseqLockById (sip);
  if (bsp == NULL)
  {
    slp = SeqLocFree (slp);
    return FALSE;
  }

  newSfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_PROT, slp);
  BioseqUnlock (bsp);
  if (newSfp == NULL)
  {
    slp = SeqLocFree (slp);
    return FALSE;
  }

  /* Make it into a protein feature */

  prp = ProtRefNew ();
  newSfp->data.value.ptrvalue = (Pointer) prp;
  if (NULL == prp)
  {
    slp = SeqLocFree (slp);
    newSfp = SeqFeatFree (newSfp);
    return FALSE;
  }

  switch (featdef_to) {
    case FEATDEF_mat_peptide_aa :
      prp->processed = 2;
      break;
    case FEATDEF_sig_peptide_aa :
      prp->processed = 3;
      break;
    case FEATDEF_transit_peptide_aa :
      prp->processed = 4;
      break;
  }

  /* Transfer unchanged fields from old feature */

  newSfp->partial     = sfp->partial;
  newSfp->excpt       = sfp->excpt;
  newSfp->exp_ev      = sfp->exp_ev;
  newSfp->pseudo      = sfp->pseudo;
  newSfp->comment     = sfp->comment;
  newSfp->qual        = sfp->qual;
  newSfp->title       = sfp->title;
  newSfp->ext         = sfp->ext;
  newSfp->cit         = sfp->cit;
  newSfp->xref        = sfp->xref;
  newSfp->dbxref      = sfp->dbxref;
  newSfp->except_text = sfp->except_text;

  /* Null out pointers to transferred fields from old feature  */
  /* so that they don't get deleted when the feature does,     */

  sfp->comment     = NULL;
  sfp->qual        = NULL;
  sfp->title       = NULL;
  sfp->ext         = NULL;
  sfp->cit         = NULL;
  sfp->xref        = NULL;
  sfp->dbxref      = NULL;
  sfp->except_text = NULL;

  /* Mark the old feature for deletion */

  sfp->idx.deleteme = TRUE;
  return TRUE;
}


NLM_EXTERN SeqLocPtr FindNucleotideLocationForProteinFeatureConversion (SeqLocPtr slp)
{
  SeqMgrFeatContext context;
  SeqFeatPtr cds;
  SeqLocPtr  slp_nuc;
  Boolean    partial5, partial3;

  cds = SeqMgrGetCDSgivenProduct (BioseqFindFromSeqLoc (slp), &context);
  if (NULL == cds)
  {
    return NULL;
  }

  slp_nuc = aaLoc_to_dnaLoc (cds, slp);
  if (slp_nuc == NULL) 
  {
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    if (partial5 && partial3) {
      slp_nuc = AsnIoMemCopy (cds->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
    }
  }
  return slp_nuc;
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* ConvertProtToImp () -                                        */
/*                                                                     */
/*---------------------------------------------------------------------*/

NLM_EXTERN Boolean ConvertProtToImpFunc (SeqFeatPtr  sfp, Uint2 featdef_to)
{
  ProtRefPtr    prp;
  SeqLocPtr     slp;
  ImpFeatPtr    ifp;
  CharPtr       name;
  CharPtr       ec;
  CharPtr       activity;
  ValNodePtr    vnp;
  GBQualPtr     gbqual = NULL;
  GBQualPtr     prevGbq;
  GBQualPtr     topOfGbqList;
  DbtagPtr      dbt;
  Char          idStr[64];
  ObjectIdPtr   oip;
  Uint2         entityID;
  BioseqPtr     bsp;
  SeqAnnotPtr   old_sap, new_sap;
  SeqFeatPtr    tmp, tmp_prev;

  /* Make sure that we have a matching peptide feature */

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT)
  {
    return FALSE;
  }
  entityID = sfp->idx.entityID;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (NULL == prp)
  {
    return FALSE;
  }

  switch (sfp->idx.subtype) {
    case FEATDEF_mat_peptide_aa :
      if (2 != prp->processed)
        return FALSE;
      break;
    case FEATDEF_sig_peptide_aa :
      if (3 != prp->processed)
        return FALSE;
      break;
    case FEATDEF_transit_peptide_aa :
      if (4 != prp->processed)
        return FALSE;
      break;
  }

  /* Convert the location from the protein */
  /* to the nucleotide Bioseq.             */

  slp = FindNucleotideLocationForProteinFeatureConversion (sfp->location);
  if (NULL == slp)
    return FALSE;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  /* move feature to correct annot */
  if (sfp->idx.parenttype == OBJ_SEQANNOT
      && (old_sap = (SeqAnnotPtr) sfp->idx.parentptr) != NULL
      && old_sap->type == 1
      && (bsp = BioseqFindFromSeqLoc (sfp->location)) != NULL)
  {
    tmp_prev = NULL;
    tmp = old_sap->data;
    while (tmp != sfp) {
      tmp_prev = tmp;
      tmp = tmp->next;
    }
    if (tmp != NULL) {
      if (tmp_prev == NULL) {
        old_sap->data = tmp->next;
      } else {
        tmp_prev->next = tmp->next;
      }
      if (old_sap->data == NULL) {
        old_sap->idx.deleteme = TRUE;
      }
      new_sap = bsp->annot;
      while (new_sap != NULL && new_sap->type != 1) {
        new_sap = new_sap->next;
      }
      if (new_sap == NULL) {
        new_sap = SeqAnnotNew ();
        new_sap->type = 1;
        new_sap->next = bsp->annot;
        bsp->annot = new_sap;
      }
      sfp->next = new_sap->data;
      new_sap->data = sfp;
      sfp->idx.parentptr = new_sap;
    }
  }
        

  /* Create a new import feature and */
  /* attach it to the feature.       */

  ifp = ImpFeatNew ();
  if (NULL == ifp)
  {
    return FALSE;
  }

  /* set key */
  ifp->key = StringSave (GetFeatureNameFromFeatureType(featdef_to));

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;

  /* Store the protein fields as  */
  /* gbqual qualifier/value pairs */

  name = NULL;
  vnp = prp->name;
  if (vnp != NULL)
    name = vnp->data.ptrvalue;
  if (name == NULL) 
    name = prp->desc;

  if (name != NULL) {
    gbqual = GBQualNew ();
    if (NULL == gbqual)
      return FALSE;
    topOfGbqList = gbqual;
    gbqual->qual = StringSave ("product");
    gbqual->val = StringSave (name);
  }

  prevGbq = gbqual;

  ec = NULL;
  vnp = prp->ec;
  if (vnp != NULL)
    ec = (CharPtr) vnp->data.ptrvalue;
  
  if (ec != NULL) {
    gbqual = GBQualNew ();
    if (NULL == gbqual)
      return FALSE;
    prevGbq->next = gbqual;
    gbqual->qual = StringSave ("EC_number");
    gbqual->val = StringSave (ec);
  }

  prevGbq = gbqual;

  activity = NULL;
  vnp = prp->activity;
  if (vnp != NULL)
    activity = (CharPtr) vnp->data.ptrvalue;
  
  if (NULL != activity) {
    gbqual = GBQualNew ();
    if (NULL == gbqual)
      return FALSE;
    prevGbq->next = gbqual;
    gbqual->qual = StringSave ("function");
    gbqual->val = StringSave (activity);
  }

  prevGbq = gbqual;

  for (vnp = prp->db; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (NULL == dbt ) 
      continue;
    if (! StringHasNoText (dbt->db)) {
      gbqual = GBQualNew ();
      if (NULL == gbqual)
        continue;
      prevGbq->next = gbqual;
      oip = dbt->tag;
      if (oip->str != NULL && (! StringHasNoText (oip->str))) {
        sprintf (idStr, "%s:%s", (CharPtr)dbt->tag, oip->str);
        gbqual->qual = StringSave ("db_xref");
        gbqual->val = StringSave (idStr);
      } else {
        sprintf (idStr, "%s:%ld", (CharPtr)dbt->tag, (long) oip->id);
        gbqual->qual = StringSave ("db_xref");
        gbqual->val = StringSave (idStr);
      }
      prevGbq = gbqual;
    }
  }

  /* Insert the new qualifiers in front of any existing ones */

  gbqual->next = sfp->qual;
  sfp->qual = topOfGbqList;

  /* Free the obsolete Protein reference */

  ProtRefFree (prp);
  return TRUE;
}


/* functions for converting from biosource */
NLM_EXTERN CharPtr SubSourceText (BioSourcePtr biop, Uint1 subtype, BoolPtr found)
{
  Int4 subtype_len = 0;
  SubSourcePtr ssp;
  CharPtr subtype_txt = NULL;
  
  if (biop == NULL || biop->subtype == NULL) return NULL;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == subtype) {
      if (found != NULL) *found = TRUE;
      if (!StringHasNoText (ssp->name)) {
        subtype_len += StringLen (ssp->name) + 1;
      }
    }
  }
  if (subtype_len == 0) return NULL;
  subtype_txt = (CharPtr) MemNew (sizeof (Char) * subtype_len);
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == subtype && !StringHasNoText (ssp->name)) {
      if (!StringHasNoText (subtype_txt)) {
        StringCat (subtype_txt, ";");
      }
      StringCat (subtype_txt, ssp->name);
    }
  }
  return subtype_txt;
}

NLM_EXTERN CharPtr OrgModText (BioSourcePtr biop, Uint1 subtype, BoolPtr found)
{
  Int4 subtype_len = 0;
  OrgModPtr omp;
  CharPtr subtype_txt = NULL;
  
  if (biop == NULL
     || biop->org == NULL 
     || biop->org->orgname == NULL 
     || biop->org->orgname->mod == NULL) {
    return NULL;
  }
     
  for (omp = biop->org->orgname->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == subtype) {
      if (found != NULL) *found = TRUE;
      if (!StringHasNoText (omp->subname)) {
        subtype_len += StringLen (omp->subname) + 1;
      }
    }
  }
  if (subtype_len == 0) return NULL;
  subtype_txt = (CharPtr) MemNew (sizeof (Char) * subtype_len);
  for (omp = biop->org->orgname->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype == subtype && !StringHasNoText (omp->subname)) {
      if (!StringHasNoText (subtype_txt)) {
        StringCat (subtype_txt, ";");
      }
      StringCat (subtype_txt, omp->subname);
    }
  }
  return subtype_txt;
}

NLM_EXTERN CharPtr NoteText (BioSourcePtr biop, CharPtr comment)
{
  CharPtr orgmod_note, subsource_note;
  Int4    text_len = 0;
  CharPtr note_text = NULL;
  
  orgmod_note = OrgModText(biop, ORGMOD_other, NULL);
  if (!StringHasNoText (orgmod_note)) {
    text_len += StringLen (orgmod_note) + 1;
  }
  subsource_note = SubSourceText (biop, SUBSRC_other, NULL);
  if (!StringHasNoText (subsource_note)) {
    text_len += StringLen (subsource_note) + 1;
  }
  if (!StringHasNoText (comment)) {
    text_len += StringLen (comment) + 1;
  }
  
  if (text_len == 0) return NULL;  
  
  note_text = (CharPtr) MemNew (sizeof(Char) * text_len);
  if (!StringHasNoText (orgmod_note)) {
    StringCat (note_text, orgmod_note);
  }
  orgmod_note = MemFree (orgmod_note);
  if (!StringHasNoText (subsource_note)) {
    if (!StringHasNoText (note_text)) {
      StringCat (note_text, ";");
    }
    StringCat (note_text, subsource_note);
  }
  subsource_note = MemFree (subsource_note);
  
  if (!StringHasNoText (comment)) {
    if (!StringHasNoText (note_text)) {
      StringCat (note_text, ";");
    }
    StringCat (note_text, comment);
  }
  return note_text;
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* ConvertBioSrcToRepeatRegion ()                                  */
/*                                                                     */
/* 9/28/2004: Changed to convert all BioSource features with notes     */
/* instead of ones with transposon or insertion_seq qualifiers.        */
/*---------------------------------------------------------------------*/

NLM_EXTERN Boolean ConvertBioSrcToRepeatRegion (SeqFeatPtr sfp, Uint2 featdef_to)
{
  BioSourcePtr  biop;
  GBQualPtr     gbqual;
  ImpFeatPtr    ifp;
  CharPtr       transposon_txt, insertion_seq_txt, note_txt;
  Boolean       is_transposon = FALSE, is_insertion_seq = FALSE;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_BIOSRC) return FALSE;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  
  transposon_txt = SubSourceText (biop, SUBSRC_transposon_name, &is_transposon);
  insertion_seq_txt = SubSourceText (biop, SUBSRC_insertion_seq_name, &is_insertion_seq);
  note_txt = NoteText (biop, sfp->comment);
  
  
  /* Create a new Import Feature */

  ifp = ImpFeatNew ();
  if (NULL == ifp)
    return FALSE;
  ifp->key = StringSave ("repeat_region");

  /* Copy relevant info from the BioSource */
  /* feature to the Import feature.        */

    
  /* Delete the old BioSource feature */

  sfp->data.value.ptrvalue = BioSourceFree (biop);

  /* Attach the new Import feature in its place */

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = ifp;
  
  if (is_transposon) {
    gbqual = GBQualNew ();
    gbqual->qual = StringSave ("mobile_element");
    gbqual->val = (CharPtr) MemNew (sizeof(Char) * (StringLen (transposon_txt) + 12));
    StringCat (gbqual->val, "transposon:");
    StringCat (gbqual->val, transposon_txt);
    gbqual->next = sfp->qual;
    sfp->qual = gbqual;
  }
  transposon_txt = MemFree (transposon_txt); 

  if (is_insertion_seq) {
    gbqual = GBQualNew ();
    gbqual->qual = StringSave ("mobile_element");
    gbqual->val = (CharPtr) MemNew (sizeof(Char) * (StringLen (insertion_seq_txt) + 19));
    StringCat (gbqual->val, "insertion sequence:");
    StringCat (gbqual->val, insertion_seq_txt);
    gbqual->next = sfp->qual;
    sfp->qual = gbqual;
  }
  insertion_seq_txt = MemFree (insertion_seq_txt); 
  
  sfp->comment = MemFree (sfp->comment);
  sfp->comment = note_txt;
  return TRUE;
}


NLM_EXTERN Boolean 
ConvertNonPseudoCDSToMiscFeat 
(SeqFeatPtr sfp,
 Boolean viral)
{
  CdRegionPtr          cdrp;
  ImpFeatPtr           ifp;
  CharPtr              noteStr;
  BioseqPtr            protBsp;
  SeqMgrFeatContext    protContext;
  CharPtr              protName = NULL;
  SeqFeatPtr           protSfp;
  ProtRefPtr           prp;
  ValNodePtr           vnp;
  Int4                 note_len = 0;
  CharPtr              viral_fmt = "nonfunctional %s due to mutation";
  CharPtr              similar_fmt = "similar to %s";

  if (sfp == NULL 
      || sfp->data.choice != SEQFEAT_CDREGION
      || sfp->product == NULL)
  {
    return FALSE;
  }
      
  /* Get the CD region part of the feature, and */
  /* the associated protein bioseq.             */

  cdrp = (CdRegionPtr) sfp->data.value.ptrvalue;
  protBsp = BioseqFindFromSeqLoc (sfp->product);

  if (protBsp == NULL) return FALSE;

  /* Convert the CDS feature to a misc_feat */

  CdRegionFree (cdrp);
  sfp->data.value.ptrvalue = NULL;

  ifp = ImpFeatNew ();
  if (NULL == ifp) return FALSE;
  ifp->key = StringSave ("misc_feature");

  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;

  /* Add a name key to the misc_feature */

  protSfp = SeqMgrGetBestProteinFeature (protBsp, &protContext);
  if (protSfp != NULL)
  {
    prp = (ProtRefPtr) protSfp->data.value.ptrvalue;

    if (prp != NULL) 
    {
      note_len = StringLen (sfp->comment) + StringLen (prp->desc) + 5;
          
      vnp = prp->name;
      if (NULL != vnp)
      {
        protName = (CharPtr) vnp->data.ptrvalue;
        if (NULL != protName) 
        {
          if (viral) {
            note_len += StringLen (viral_fmt) + StringLen (protName);
          } else {   
            note_len += StringLen (similar_fmt) + StringLen (protName);
          }
        }
      }  
      noteStr = (CharPtr) MemNew (sizeof (Char) * note_len);
        
      if (NULL != protName) {
        if (viral) {
          sprintf (noteStr, viral_fmt, protName);
        } else {
          sprintf (noteStr, similar_fmt, protName);
        }
      }
      if (!StringHasNoText (prp->desc)) {
        if (!StringHasNoText (noteStr)) {
          StringCat (noteStr, "; ");
        }
        StringCat (noteStr, prp->desc);
      }
      if (!StringHasNoText (sfp->comment)) {
        if (!StringHasNoText (noteStr)) {
          StringCat (noteStr, "; ");
        }
        StringCat (noteStr, sfp->comment);
      }
      sfp->comment = MemFree (sfp->comment);
      sfp->comment = noteStr;
    }
  }

  /* set the subtype to zero so that it will be reindexed */
  sfp->idx.subtype = 0;
  return TRUE;
}


NLM_EXTERN Uint1 RnaTypeFromFeatdef (Uint2 featdef)
{
  switch (featdef) 
  {
    case FEATDEF_preRNA:
      return 1;
      break;
    case FEATDEF_mRNA:
      return 2;
      break;
    case FEATDEF_tRNA:
      return 3;
      break;
    case FEATDEF_rRNA:
      return 4;
      break;
    case FEATDEF_snRNA:
      return 8;
      break;
    case FEATDEF_scRNA:
      return 8;
      break;
    case FEATDEF_snoRNA:
      return 8;
      break;
    case FEATDEF_ncRNA:
      return 8;
      break;
    case FEATDEF_tmRNA:
      return 9; 
      break;
    case FEATDEF_otherRNA:
    default:
      return 255;
      break;
  }
}


NLM_EXTERN Boolean ConvertRegionToRNAFunc 
(SeqFeatPtr sfp,
 Uint2      featdef_to)
{
  RnaRefPtr  rrp;
  CharPtr    str, new_comment;
  Boolean    add_to_comment = FALSE;
  Int4       len;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION)
  {
    return FALSE;
  }

  str = (CharPtr) sfp->data.value.ptrvalue;
  rrp = RnaRefFromLabel (featdef_to, str, &add_to_comment);

  sfp->data.choice = SEQFEAT_RNA;
  sfp->data.value.ptrvalue = (Pointer) rrp;

  if (add_to_comment) {
    if (sfp->comment == NULL) {
      sfp->comment = str;
      str = NULL;
    } else {      
      len = StringLen (sfp->comment) + StringLen (str) + 3;
      new_comment = MemNew (sizeof (Char) * len);
      sprintf (new_comment, "%s; %s", sfp->comment, str);
      sfp->comment = MemFree (sfp->comment);
      str = MemFree (str);
      sfp->comment = new_comment;
    }
  } else {
    str = MemFree (str);
  }
  return TRUE;
  
  return TRUE;
}


NLM_EXTERN CharPtr GetImportFeatureName (Uint2 featdef_key)
{
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label = NULL;

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) 
  {
    if (featdef_key == key)
    {
      return curr->typelabel;
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  return NULL;
}


NLM_EXTERN Boolean ConvertRegionToImpFunc (SeqFeatPtr sfp, Uint2 featdef_to)
{
  GBQualPtr          gbqual;
  ImpFeatPtr         ifp;
  CharPtr            str;
  CharPtr            featname_to;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_REGION) return FALSE;
  featname_to = GetImportFeatureName (featdef_to);
  ifp = ImpFeatNew ();
  if (NULL == ifp)
    return FALSE;

  str = (CharPtr) sfp->data.value.ptrvalue;
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  if (featname_to == NULL)
  {
    ifp->key = StringSave ("misc_feature");
  }
  else
  {
    ifp->key = StringSave (featname_to);
  }

  if (! StringHasNoText (str)) {
    gbqual = GBQualNew ();
    if (gbqual != NULL) {
      gbqual->qual = StringSave ("note");
      gbqual->val = str;
      gbqual->next = sfp->qual;
      sfp->qual = gbqual;
    }
  }
  return TRUE;
}


NLM_EXTERN Boolean ConvertImpToImpFunc (SeqFeatPtr sfp, Uint2 featdef_to)
{
  ImpFeatPtr         ifp;
  CharPtr            featname;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_IMP) return FALSE;

  ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
  if (NULL == ifp)
    return FALSE;
  
  featname = GetImportFeatureName (featdef_to);
  ifp->key = MemFree (ifp->key);
  if (featname == NULL)
  {
    ifp->key = StringSave ("misc_feature");
  }
  else
  {
    ifp->key = StringSave (featname);
  }
  sfp->idx.subtype = 0;

  return TRUE;
}


NLM_EXTERN Boolean ConvertGeneToMiscFeatFunc 
(SeqFeatPtr sfp,
 Uint2      featdef_to)
{
  ImpFeatPtr  ifp;
  CharPtr     new_comment;
  GeneRefPtr  grp;
  Int4        comment_len = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE)
  {
    return FALSE;
  }
  ifp = ImpFeatNew ();
  if (NULL == ifp)
  {
    return FALSE;
  }

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp != NULL) 
  {
    if (!StringHasNoText (grp->locus))
    {
      comment_len += StringLen (grp->locus) + 2;
    }
    if (!StringHasNoText (grp->desc))
    {
      comment_len += StringLen (grp->desc) + 2;
    }
  }
  if (comment_len == 0) 
  {
    /* nothing to add to comment */
  }
  else
  {
    /* add one for terminating NULL */
    comment_len++;
    if (!StringHasNoText (sfp->comment))
    {
      comment_len += StringLen (sfp->comment) + 2;
    }

    new_comment = (CharPtr) MemNew (sizeof (Char) * comment_len);
    /* NOTE - I don't have to check for grp == NULL because
     * comment_len would only have been > 0 if grp had existed
     * and had nonempty fields.
     */
    if (!StringHasNoText (grp->desc))
    {
      StringCat (new_comment, grp->desc);
      StringCat (new_comment, "; ");    
    }
    if (!StringHasNoText (grp->locus))
    {
      StringCat (new_comment, grp->locus);
      StringCat (new_comment, "; ");    
    }
    if (!StringHasNoText (sfp->comment))
    {
      StringCat (new_comment, sfp->comment);
      StringCat (new_comment, "; ");    
    }
    /* remove last semicolon */
    new_comment[StringLen (new_comment) - 2] = 0;
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = new_comment;
  }

  sfp->data.value.ptrvalue =
    GeneRefFree ((GeneRefPtr) sfp->data.value.ptrvalue);
  sfp->data.choice = SEQFEAT_IMP;
  sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave ("misc_feature");
  return TRUE;
}




/* For mat-peptide instantiation */
static SeqIdPtr MakeMatPeptideProductId (SeqLocPtr mat_peptide_loc)
{
  ObjectIdPtr oip;
  SeqIdPtr    sip;
  Char        id_str[500];
  Int4        len;
  BioseqPtr   bsp;

  if (mat_peptide_loc == NULL) return NULL;

  bsp = BioseqFindFromSeqLoc (mat_peptide_loc);
  if (bsp == NULL) return NULL;
  sip = bsp->id;
  while (sip != NULL && sip->choice != SEQID_OTHER) {
    sip = sip->next;
  }
  if (sip == NULL) {
    sip = bsp->id;
    while (sip != NULL && sip->choice != SEQID_GENBANK) {
      sip = sip->next;
    }
  }
  if (sip == NULL) {
    sip = SeqIdFindBest (bsp->id, 0);
  }
  SeqIdWrite (sip, id_str, PRINTID_TEXTID_ACC_ONLY, sizeof (id_str) - 41);
  len = StringLen (id_str);
  sprintf (id_str + len, ":%d-%d", SeqLocStart (mat_peptide_loc) + 1, SeqLocStop (mat_peptide_loc) + 1);
  
  oip = ObjectIdNew ();
  oip->str = StringSave (id_str);
  sip = ValNodeNew (NULL);
  sip->choice = SEQID_LOCAL;
  sip->data.ptrvalue = oip;  
  return sip;
}


static void InstantiateMatPeptideProductForProteinFeature (SeqFeatPtr sfp, Pointer data)
{
  BioseqPtr mat_bsp, prot_bsp;
  SeqEntryPtr master, sep, old;
  Int4        i, len;
  Int2        residue;
  SeqLocPtr   slp;
  ByteStorePtr src, dst;
  ProtRefPtr   prp_orig, prp_mat;
  SeqFeatPtr   sfp_mat;
  SeqDescrPtr  sdp;
  MolInfoPtr   mip;
  Boolean      partial5, partial3;
  Char         defline_buf[1024];

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_mat_peptide_aa || sfp->product != NULL) {
    return;
  }

  prp_orig = sfp->data.value.ptrvalue;

  prot_bsp = BioseqFindFromSeqLoc (sfp->location);
  if (prot_bsp == NULL) {
    return;
  }
  master = GetBestTopParentForData (sfp->idx.entityID, prot_bsp);
  if (master == NULL) return;

  src = (ByteStorePtr) prot_bsp->seq_data;

  mat_bsp = BioseqNew ();
  if (mat_bsp == NULL) {
     return;
  }
  mat_bsp->mol = Seq_mol_aa;
  mat_bsp->repr = Seq_repr_raw;
  mat_bsp->seq_data_type = Seq_code_ncbieaa;
  mat_bsp->length = SeqLocLen (sfp->location);
  dst = BSNew (0);
  mat_bsp->seq_data = (SeqDataPtr) dst;
  BSSeek (dst, 0, SEEK_SET);

  for (slp = SeqLocFindNext (sfp->location, NULL);
       slp != NULL;
       slp = SeqLocFindNext (sfp->location, slp)) {
    BSSeek (src, SeqLocStart (slp), SEEK_SET);
    len = SeqLocLen (slp);
    for (i = 0; i < len; i++) {
      residue = BSGetByte (src);
      BSPutByte (dst, residue);
    }
  }

  old = SeqEntrySetScope (master);
  
  /*mat_bsp->id = MakeNewProteinSeqId (sfp->location, NULL); */
  mat_bsp->id = MakeMatPeptideProductId (sfp->location);
  SeqMgrAddToBioseqIndex (mat_bsp);
  SeqEntrySetScope (old);
  sep = SeqEntryNew ();
  if (sep != NULL) {
    sep->choice = 1;
    sep->data.ptrvalue = (Pointer) mat_bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) mat_bsp, sep);
  }
  SetSeqFeatProduct (sfp, mat_bsp);
  if (sep != NULL) {
    AddSeqEntryToSeqEntry (master, sep, TRUE);
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  /* set molinfo for new protein sequence */
  sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
  if (sdp != NULL) {
    mip = MolInfoNew ();
    sdp->data.ptrvalue = (Pointer) mip;
    if (mip != NULL) {
      mip->biomol = 8;
      mip->tech = 13;
      if (partial5 && partial3) {
        mip->completeness = 5;
      } else if (partial5) {
        mip->completeness = 3;
      } else if (partial3) {
        mip->completeness = 4;
      /*
      } else if (partial) {
        mip->completeness = 2;
      */
      } else {
        mip->completeness = 0;
      }
    }
  }


  /* create protein feature for mat_bsp */
  sfp_mat = CreateNewFeature (SeqMgrGetSeqEntryForData (mat_bsp), NULL,
                          SEQFEAT_PROT, NULL);

  sfp_mat->location = SeqLocIntNew (0, mat_bsp->length - 1, Seq_strand_plus, SeqIdDup (mat_bsp->id));
  SetSeqLocPartial (sfp_mat->location, partial5, partial3);
  if (partial5 || partial3) {
    sfp_mat->partial = TRUE;
  }

  prp_mat = AsnIoMemCopy (prp_orig, (AsnReadFunc) ProtRefAsnRead, (AsnWriteFunc) ProtRefAsnWrite);
  prp_mat->processed = 0;
  sfp_mat->data.value.ptrvalue = prp_mat;
  if (sfp->comment != NULL) {
    sfp_mat->comment = StringSave (sfp->comment);
  }

  /* add title */
  SeqMgrIndexFeatures (0, mat_bsp);
  if (NewCreateDefLineBuf (NULL, mat_bsp, defline_buf, sizeof (defline_buf), FALSE, FALSE)
      && !StringHasNoText (defline_buf)) {
    sdp = CreateNewDescriptor (sep, Seq_descr_title);
    sdp->data.ptrvalue = StringSave (defline_buf);
  }

}


NLM_EXTERN void ExtraCDSCreationActions (SeqFeatPtr cds, SeqEntryPtr parent_sep)
{
  ByteStorePtr       bs;
  CharPtr            prot, ptr;
  BioseqPtr          bsp;
  Char               ch;
  Int4               i;
  SeqEntryPtr        psep, nsep;
  MolInfoPtr         mip;
  ValNodePtr         vnp, descr;
  SeqFeatPtr         prot_sfp;
  ProtRefPtr         prp;
  Boolean            partial5, partial3;

  if (cds == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  /* Create corresponding protein sequence data for the CDS */

  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
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
  i = StringLen (prot);
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
  bsp->id = MakeNewProteinSeqId (cds->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  
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
  if (partial5 && partial3) {
    mip->completeness = 5;
  } else if (partial5) {
    mip->completeness = 3;
  } else if (partial3) {
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
  SetSeqFeatProduct (cds, bsp);
  
  prp = ProtRefNew ();
  
  if (prp != NULL) {
    prot_sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (prot_sfp != NULL) {
      prot_sfp->data.value.ptrvalue = (Pointer) prp;
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    }
  }
}


NLM_EXTERN SeqFeatPtr GetProtFeature (BioseqPtr protbsp)
{
  SeqMgrFeatContext fcontext;
  SeqAnnotPtr sap;
  SeqFeatPtr prot_sfp;
  ProtRefPtr prp;

  if (protbsp == NULL) return NULL;

  prot_sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &fcontext);
  if (prot_sfp == NULL) {
    sap = protbsp->annot;
    while (sap != NULL && prot_sfp == NULL) {
      if (sap->type == 1) {
        prot_sfp = sap->data;
        while (prot_sfp != NULL
               && (prot_sfp->data.choice != SEQFEAT_PROT
                   || (prp = prot_sfp->data.value.ptrvalue) == NULL
                   || prp->processed != 0)) {
          prot_sfp = prot_sfp->next;
        }
      }
      sap = sap->next;
    }
  }
  return prot_sfp;
}


NLM_EXTERN Boolean ConvertMiscFeatToGene (SeqFeatPtr sfp)
{
  GeneRefPtr grp;
  CharPtr    cp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_misc_feature) {
    return FALSE;
  }
  sfp->data.value.ptrvalue = ImpFeatFree (sfp->data.value.ptrvalue);
  grp = GeneRefNew ();
  sfp->data.value.ptrvalue = grp;
  sfp->data.choice = SEQFEAT_GENE;
  sfp->idx.subtype = 0;

  if (!StringHasNoText (sfp->comment)) {
    cp = StringChr (sfp->comment, ';');
    if (cp != NULL) {
      *cp = 0;
    }
    grp->locus = StringSave (sfp->comment);
    if (cp != NULL) {
      cp = StringSave (cp + 1);
    }
    sfp->comment = MemFree (sfp->comment);
    sfp->comment = cp;
  }
  return TRUE;
}


NLM_EXTERN Boolean ConvertMiscFeatToCodingRegion (SeqFeatPtr sfp)
{
  BioseqPtr bsp, prot_bsp;
  SeqFeatPtr prot;
  ProtRefPtr prp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_misc_feature) {
    return FALSE;
  }

  sfp->data.value.ptrvalue = ImpFeatFree (sfp->data.value.ptrvalue);
  sfp->data.value.ptrvalue = CdRegionNew ();
  sfp->data.choice = SEQFEAT_CDREGION;
  sfp->idx.subtype = 0;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp != NULL) {
    ExtraCDSCreationActions (sfp, GetBestTopParentForData (bsp->idx.entityID, bsp));
    if (!StringHasNoText (sfp->comment)) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      prot = GetProtFeature (prot_bsp);
      if (prot != NULL) {
        prp = prot->data.value.ptrvalue;
        ValNodeAddPointer (&prp->name, 0, sfp->comment);
        sfp->comment = NULL;
      }
    }
  }

  return TRUE;
}


NLM_EXTERN void InstantiateMatPeptideProducts (SeqEntryPtr sep)
{
  VisitFeaturesInSep (sep, NULL, InstantiateMatPeptideProductForProteinFeature);
}


typedef struct tsa_id_callback {
  SeqEntryPtr top_sep;
  CharPtr suffix;
} TSAIdCallbackData, PNTR TSAIdCallbackPtr;


static void MakeTSAIdsCallback (BioseqPtr bsp, Pointer data)
{
  Int4            gpid = 0;
  SeqIdPtr        sip, sip_tsa = NULL, sip_local = NULL;
  TSAIdCallbackPtr t;
  SeqIdPtr        sip_new = NULL;
  DbtagPtr        dbtag, dbtag_old = NULL;
  CharPtr         id_fmt = "gpid:%d";
  ObjectIdPtr     oip = NULL;
  Int4            id_num = 0;
  CharPtr         id_str = NULL;
  CharPtr         cp;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL || !IsTSA (bsp)) {
    return;
  }

  t = (TSAIdCallbackPtr) data;

  for (sip = bsp->id;
       sip != NULL;
       sip = sip->next) 
  {
    if (sip->choice == SEQID_LOCAL) 
    {
      sip_local = sip;
    } 
    else if (sip->choice == SEQID_GENERAL
               && (dbtag_old = (DbtagPtr) sip->data.ptrvalue) != NULL
               && StringNCmp (dbtag_old->db, "gpid:", 5) == 0) 
    {
      sip_tsa = sip;
    }
  }

  if (sip_tsa == NULL && sip_local == NULL) {
    return;
  }
  gpid = GetGenomeProjectID (bsp);
  if (gpid <= 0) {
    return;
  }

  dbtag = DbtagNew ();
  dbtag->db = MemNew (sizeof (Char) * (StringLen (id_fmt) + 15));
  sprintf (dbtag->db, id_fmt, gpid);
  dbtag->tag = ObjectIdNew ();
  
  if (sip_tsa != NULL && (dbtag_old = (DbtagPtr) sip_tsa->data.ptrvalue) != NULL) {
    oip = dbtag_old->tag;
    sip = sip_tsa;
  } else if (sip_local != NULL) {
    oip = sip_local->data.ptrvalue;
    sip = sip_local;
  }

  if (oip == NULL) {
    return;
  }

  if (oip->str == NULL) {
    id_num = oip->id;
  } else {
    id_str = StringSave (oip->str);
    if (sip == sip_tsa && (cp = StringRChr (id_str, '.')) != NULL) {
      *cp = 0;
    }
  }

  if (t->suffix == NULL) {
    if (id_str == NULL) {
      dbtag->tag->id = id_num;
    } else {
      dbtag->tag->str = StringSave (id_str);
    }
  } else {
    if (id_str == NULL) {
      dbtag->tag->str = (CharPtr) MemNew (sizeof (Char) * (16 + StringLen (t->suffix)));
      sprintf (dbtag->tag->str, "%d.%s", id_num, t->suffix);
    } else {
      dbtag->tag->str = (CharPtr) MemNew (sizeof (Char) * (StringLen (oip->str) + StringLen (t->suffix) + 2));
      sprintf (dbtag->tag->str, "%s.%s", id_str, t->suffix);
    }
  }
  id_str = MemFree (id_str);
  sip_new = ValNodeNew (NULL);
  sip_new->choice = SEQID_GENERAL;
  sip_new->data.ptrvalue = dbtag;
  ReplaceSeqIdWithSeqId (sip, sip_new, t->top_sep);
}     
          
  
NLM_EXTERN void ConvertLocalIdsToTSAIds (SeqEntryPtr sep, CharPtr suffix)
{
  TSAIdCallbackData t;

  t.top_sep = sep;
  t.suffix = suffix;
  VisitBioseqsInSep (sep, &t, MakeTSAIdsCallback);
}


static void ConvertLocalIdsToBarcodeIdsCallback (BioseqPtr bsp, Pointer data)
{
  SeqIdPtr        sip_local = NULL;
  SeqEntryPtr     top_sep;
  SeqIdPtr        sip_new;
  DbtagPtr        dbtag;
  ObjectIdPtr     oip = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  top_sep = (SeqEntryPtr) data;

  for (sip_local = bsp->id;
       sip_local != NULL && sip_local->choice != SEQID_LOCAL;
       sip_local = sip_local->next)
  {}
  if (sip_local == NULL) return;
  oip = sip_local->data.ptrvalue;
  if (oip == NULL) return;

  dbtag = DbtagNew ();
  dbtag->db = StringSave ("uoguelph");
  dbtag->tag = ObjectIdNew ();
  if (oip->str == NULL) {
    dbtag->tag->id = oip->id;
  } else {
    dbtag->tag->str = StringSave (oip->str);
  }
  sip_new = ValNodeNew (NULL);
  sip_new->choice = SEQID_GENERAL;
  sip_new->data.ptrvalue = dbtag;
  ReplaceSeqIdWithSeqId (sip_local, sip_new, top_sep);
}     
          
  
NLM_EXTERN void ConvertLocalIdsToBarcodeIds (SeqEntryPtr sep)
{
  VisitBioseqsInSep (sep, sep, ConvertLocalIdsToBarcodeIdsCallback);
}


NLM_EXTERN Int4 GetDeflinePosForFieldType (ValNodePtr field)
{
  Int4    i, rval = -1;
  CharPtr name;

  name = SummarizeFieldType (field);
  if (StringICmp (name, "specimen-voucher") == 0) {
    rval = DEFLINE_POS_Specimen_voucher;
  } else {
    for (i = 0; i < numDefLineModifiers; i++) {
      if (StringICmp (name, DefLineModifiers[i].name) == 0) {
        rval = i;
        break;
      }
    }
  }
  name = MemFree (name);
  return rval;
}


static void RemoveUnusedFieldTypes (FieldTypePtr PNTR orig_list)
{
  ValNodePtr vnp, prev = NULL, vnp_next;

  if (orig_list == NULL || *orig_list == NULL) {
    return;
  }
  for (vnp = *orig_list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (GetDeflinePosForFieldType (vnp) < 0) {
      if (prev == NULL) {
        *orig_list = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = FieldTypeFree (vnp);
    } else {
      prev = vnp;
    }
  }  
}


static Boolean RemoveMatchingFieldType (FieldTypePtr PNTR orig_list, FieldTypePtr match)
{
  ValNodePtr vnp, prev = NULL, vnp_next;
  Boolean    rval = FALSE;

  if (orig_list == NULL || *orig_list == NULL || match == NULL) {
    return rval;
  }

  for (vnp = *orig_list; vnp != NULL && !rval; vnp = vnp_next) {
    vnp_next = vnp->next;
    if (CompareFieldTypes (vnp, match) == 0) {
      if (prev == NULL) {
        *orig_list = vnp->next;
      } else {
        prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = FieldTypeFree (vnp);
      rval = TRUE;
    } else {
      prev = vnp;
    }
  }
  return rval;
}


static Boolean ListHasMatchingFieldType (FieldTypePtr list, FieldTypePtr match)
{
  Boolean rval = FALSE;

  if (list == NULL || match == NULL) return rval;

  while (list != NULL && !rval) {
    if (CompareFieldTypes (list, match) == 0) {
      rval = TRUE;
    } else {
      list = list->next;
    }
  }
  return rval;
}


static Int4 DefLineFieldTypeSortOrder [] = {
  Source_qual_strain,
  Source_qual_isolate,
  Source_qual_clone,
  Source_qual_haplotype,
  Source_qual_cultivar,
  Source_qual_specimen_voucher,
  Source_qual_ecotype,
  Source_qual_type,
  Source_qual_serotype,
  Source_qual_authority,
  Source_qual_breed
};


static int CompareFieldTypeByImportance (FieldTypePtr field1, FieldTypePtr field2)
{
  int rval = 0;
  Int4 index, num_defline_qual_sort_order;
  SourceQualChoicePtr scp1, scp2;

  if (field1 == NULL && field2 == NULL) {
    rval = 0;
  } else if (field1 == NULL) {
    rval = 1;
  } else if (field2 == NULL) {
    rval = -1;
  } else if (field1->choice == FieldType_source_qual && field2->choice == FieldType_source_qual) {
    scp1 = field1->data.ptrvalue;
    scp2 = field2->data.ptrvalue;
    num_defline_qual_sort_order = sizeof (DefLineFieldTypeSortOrder) / sizeof (Int4);
    for (index = 0; index < num_defline_qual_sort_order; index++)
    {
      if (scp1->data.intvalue == DefLineFieldTypeSortOrder [ index ]) return -1;
      if (scp2->data.intvalue == DefLineFieldTypeSortOrder [ index ]) return 1;
    }
    rval = CompareFieldTypes (field1, field2);
  } else {
    rval = CompareFieldTypes (field1, field2);
  }
  return rval;
}

static int LIBCALLBACK SortFieldTypeByImportance (
  VoidPtr ptr1,
  VoidPtr ptr2
)
{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL && ptr2 == NULL) return 0;
  
  if (ptr1 == NULL && ptr2 != NULL) return -1;
  if (ptr1 != NULL && ptr2 == NULL) return 1;
 
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  if (vnp1->data.ptrvalue == NULL || vnp2->data.ptrvalue == NULL) return 0;

  return CompareFieldTypeByImportance (vnp1, vnp2);
}


typedef struct uniqbiosource {
  BioSourcePtr biop;
  ValNodePtr   available_fields;
  ValNodePtr   strings;
} UniqBioSourceData, PNTR UniqBioSourcePtr;

static Boolean AddQualToUniqBioSource (
  UniqBioSourcePtr u,
  FieldTypePtr     field
)
{
  CharPtr             val = NULL, tmp;
  SourceQualChoicePtr q;
  Boolean             rval = FALSE;

  if (u == NULL || field == NULL) return FALSE;
  if (field->choice != FieldType_source_qual) return FALSE;
  q = (SourceQualChoicePtr) field->data.ptrvalue;
  if (q == NULL || q->choice != SourceQualChoice_textqual) return FALSE;

  val = GetSourceQualFromBioSource (u->biop, q, NULL);
  if (StringHasNoText (val)) {
    val = MemFree (val);
  } else if (q->data.intvalue == Source_qual_specimen_voucher && StringNICmp (val, "personal:", 9) == 0) {
    tmp = StringSave (val + 9);
    val = MemFree (val);
    val = tmp;
  } else if (IsNonTextSourceQual (q->data.intvalue)) {
    val = MemFree (val);
    val = StringSave (GetSourceQualName (q->data.intvalue));
  }
  if (val != NULL) {
    ValNodeAddPointer (&(u->strings), 0, val);
    rval = TRUE;
  }
  RemoveMatchingFieldType (&(u->available_fields), field);

  u->available_fields = ValNodeSort (u->available_fields, SortFieldTypeByImportance);
  return rval;
}
 

static UniqBioSourcePtr UniqBioSourceNew (BioSourcePtr biop)
{
  UniqBioSourcePtr u;

  if (biop == NULL) return NULL;
  u = (UniqBioSourcePtr) MemNew (sizeof (UniqBioSourceData));
  u->biop = biop;
  u->available_fields = GetSourceQualFieldListFromBioSource (u->biop);
  u->strings = NULL;

  /* add tax name as first string */
  AddQualToUniqBioSource (u, u->available_fields);
  RemoveUnusedFieldTypes (&(u->available_fields));

  return u;
}

static UniqBioSourcePtr UniqBioSourceFree (UniqBioSourcePtr u)
{
  if (u != NULL) {
    u->available_fields = FieldTypeListFree (u->available_fields);
    u->strings = ValNodeFreeData (u->strings);
    u = MemFree (u);
  }
  return u;
}


static UniqBioSourcePtr UniqBioSourceCopy (UniqBioSourcePtr u)
{
  UniqBioSourcePtr u2;
  ValNodePtr       vnp;

  if (u == NULL) return NULL;
  u2 = (UniqBioSourcePtr) MemNew (sizeof (UniqBioSourceData));
  u2->biop = u->biop;
  u2->available_fields = FieldTypeListCopy (u->available_fields);
  for (vnp = u->strings; vnp != NULL; vnp = vnp->next) {
    ValNodeAddPointer (&(u2->strings), 0, StringSave (vnp->data.ptrvalue));
  }
  return u2;
}


/* The CompareOrganismDescriptors function compares the contents of the
 * lists of strings for each BioSrcDesc item.
 * The function returns:
 *     -1 if org1 < org2
 *      0 if org1 = org2
 *      1 if org1 > org2
 */
static int CompareUniqBioSource (
  UniqBioSourcePtr org1,
  UniqBioSourcePtr org2
)
{
  ValNodePtr vnp1, vnp2;
  int cmpval;

  vnp1 = org1->strings;
  vnp2 = org2->strings;

  while (vnp1 != NULL && vnp2 != NULL)
  {
    cmpval = StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    if (cmpval != 0) return cmpval;

    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return 0;
  }
  else if (vnp1 != NULL && vnp2 == NULL)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}


static Boolean RemoveFieldFromUniqBioSource (UniqBioSourcePtr u, FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (u != NULL) {
    rval = RemoveMatchingFieldType (&(u->available_fields), field);
  }
  return rval;
}


static ValNodePtr UniqBioSourceListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = UniqBioSourceFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list= list_next;
  }
  return list;
}


static ValNodePtr UniqBioSourceListCopy (ValNodePtr orig)
{
  ValNodePtr list = NULL, prev = NULL, vnp;
  UniqBioSourcePtr u;

  while (orig != NULL) {
    u = (UniqBioSourcePtr) orig->data.ptrvalue;
    if (u != NULL && u->biop != NULL) {
      vnp = ValNodeNew (prev);
      vnp->choice = 0;
      vnp->data.ptrvalue = UniqBioSourceCopy (u);
      if (prev == NULL) {
        list = vnp;
      }
      prev = vnp;
    }
    orig = orig->next;
  }
  return list;
}


static int LIBCALLBACK SortUniqBioSource (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return CompareUniqBioSource (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static ValNodePtr UniqBioSourceListSort (ValNodePtr orig)
{
  orig = ValNodeSort (orig, SortUniqBioSource);
  return orig;
}


static Boolean AddQualToUniqBioSourceList (ValNodePtr list, FieldTypePtr field)
{
  Boolean rval = FALSE;
  while (list != NULL) {
    rval |= AddQualToUniqBioSource (list->data.ptrvalue, field);
    list = list->next;
  }
  return rval;
}


static Boolean RemoveFieldFromUniqBioSourceList (ValNodePtr list, FieldTypePtr field)
{
  Boolean rval = FALSE;

  while (list != NULL) {
    rval |= RemoveFieldFromUniqBioSource (list->data.ptrvalue, field);
    list = list->next;
  }
  return rval;
}


/* The UniqBioSrcGrp structure contains a list of UniqBioSrc items
 * for which the contents of the descriptive strings list are identical,
 * i.e., all the organisms in the group would have the same description
 * if you used the modifiers used to generate this list of strings.
 * The structure also contains the number of organisms in the list
 * so that it will be easy to tell that the UniqBioSrcGrp now contains a
 * single organism with a unique description.
 */
typedef struct uniqbiosrcgrp {
  Int4 num_biop;
  ValNodePtr biop_list;
} UniqBioSrcGrpData, PNTR UniqBioSrcGrpPtr;


static UniqBioSrcGrpPtr UniqBioSrcGrpNew (ValNodePtr biop_list)
{
  UniqBioSrcGrpPtr g;

  g = (UniqBioSrcGrpPtr) MemNew (sizeof (UniqBioSrcGrpData));
  g->num_biop = ValNodeLen (biop_list);
  g->biop_list = UniqBioSourceListCopy (biop_list);
  return g;
}


static UniqBioSrcGrpPtr UniqBioSrcGrpFree (UniqBioSrcGrpPtr g)
{
  if (g != NULL) {
    g->biop_list = UniqBioSourceListFree (g->biop_list);
    g = MemFree (g);
  }
  return g;
}


static UniqBioSrcGrpPtr UniqBioSrcGrpCopy (UniqBioSrcGrpPtr orig)
{
  UniqBioSrcGrpPtr g;

  if (orig == NULL) {
    return NULL;
  }

  g = UniqBioSrcGrpNew (orig->biop_list);
  return g;
}


static Boolean AddQualToUniqBioSrcGrp (UniqBioSrcGrpPtr g, FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (g != NULL) {
    rval = AddQualToUniqBioSourceList (g->biop_list, field);
    if (rval) {
      g->biop_list = UniqBioSourceListSort (g->biop_list);
    }
  }
  return rval;
}


static Boolean RemoveFieldFromUniqBioSrcGrp (UniqBioSrcGrpPtr g, FieldTypePtr field)
{
  Boolean rval = FALSE;

  if (g != NULL) {
    rval = RemoveFieldFromUniqBioSourceList (g->biop_list, field);
  }
  return rval;
}


static FieldTypePtr GetAllPresentQualsForGroup (UniqBioSrcGrpPtr g)
{
  ValNodePtr vnp;
  FieldTypePtr match_list = NULL, ft, ft_next, ft_prev;
  UniqBioSourcePtr u;

  if (g == NULL || g->num_biop < 2) {
    return NULL;
  }

  u = g->biop_list->data.ptrvalue;
  match_list = FieldTypeListCopy (u->available_fields);
  for (vnp = g->biop_list->next; vnp != NULL && match_list != NULL; vnp = vnp->next) {
    u = vnp->data.ptrvalue;
    ft = match_list;
    ft_prev = NULL;
    while (ft != NULL) {
      ft_next = ft->next;
      if (ListHasMatchingFieldType (u->available_fields, ft)) {
        ft_prev = ft;
      } else{
        if (ft_prev == NULL) {
          match_list = ft->next;
        } else {
          ft_prev->next = ft->next;
        }
        ft->next = NULL;
        ft = FieldTypeFree (ft);
      }
      ft = ft_next;
    }
  }
  return match_list;
}


static FieldTypePtr GetAllQualsForGroup (UniqBioSrcGrpPtr g)
{
  ValNodePtr vnp, tmp;
  FieldTypePtr field_list_head = NULL;
  FieldTypePtr field_list_tail = NULL;
  UniqBioSourcePtr u;

  if (g == NULL || g->num_biop < 2) {
    return NULL;
  }

  for (vnp = g->biop_list; vnp != NULL; vnp = vnp->next) {
    u = vnp->data.ptrvalue;
    if (u != NULL && u->available_fields != NULL) {
      tmp = ValNodeLink (&field_list_tail, FieldTypeListCopy (u->available_fields));
      if (field_list_head == NULL) {
        field_list_head = tmp;
      }
      if (tmp != NULL) {
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
      }
      field_list_tail = tmp;
    }
  }
  SortUniqueFieldTypeList (&field_list_head);
  return field_list_head;
}


static ValNodePtr UniqBioSrcGrpListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = UniqBioSrcGrpFree (list->data.ptrvalue);
    list = ValNodeFree (list); 
    list = list_next;
  }
  return list;
}


static ValNodePtr UniqBioSrcGrpListCopy (ValNodePtr orig)
{
  ValNodePtr list = NULL, prev = NULL, vnp;
  UniqBioSrcGrpPtr g;

  while (orig != NULL) {
    g = (UniqBioSrcGrpPtr) orig->data.ptrvalue;
    if (g != NULL) {
      vnp = ValNodeNew (prev);
      vnp->choice = 0;
      vnp->data.ptrvalue = UniqBioSrcGrpCopy (g);
      if (prev == NULL) {
        list = vnp;
      }
      prev = vnp;
    }
    orig = orig->next;
  }
  return list;
}


/* NOTE - we want to sort groups from most biops to least biops */
static int LIBCALLBACK SortUniqBioSrcGrp (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  UniqBioSrcGrpPtr g1, g2;
  int              rval = 0;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL && vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      g1 = (UniqBioSrcGrpPtr) vnp1->data.ptrvalue;
      g2 = (UniqBioSrcGrpPtr) vnp2->data.ptrvalue;
      if (g1->num_biop > g2->num_biop) {
        rval = -1;
      } else if (g1->num_biop < g2->num_biop) {
        rval = 1;
      }
    }
  }
  return rval;
}


static ValNodePtr BioSrcGrpListSort (ValNodePtr orig)
{
  orig = ValNodeSort (orig, SortUniqBioSrcGrp);
  return orig;
}


static Boolean RemoveFieldFromUniqBioSrcGrpList (ValNodePtr list, FieldTypePtr field)
{
  Boolean rval = FALSE;

  while (list != NULL) {
    rval |= RemoveFieldFromUniqBioSrcGrp (list->data.ptrvalue, field);
    list = list->next;
  }
  return rval;
}


static void ReGroupUniqBioSrcGrpList (ValNodePtr list)
{
  ValNodePtr list_next, vnp;
  UniqBioSrcGrpPtr g, g2;

  while (list != NULL) {
    g = (UniqBioSrcGrpPtr) list->data.ptrvalue;
    vnp = g->biop_list;
    while (vnp != NULL && vnp->next != NULL 
           && CompareUniqBioSource (vnp->data.ptrvalue, vnp->next->data.ptrvalue) == 0) {
      vnp = vnp->next;
    }
    if (vnp != NULL && vnp->next != NULL) {
      g2 = UniqBioSrcGrpNew (NULL);
      g2->biop_list = vnp->next;
      g2->num_biop = ValNodeLen (vnp->next);
      vnp->next = NULL;
      g->num_biop -= g2->num_biop;
      list_next = ValNodeNew (NULL);
      list_next->data.ptrvalue = g2;
      list_next->next = list->next;
      list->next = list_next;
    }
    list = list->next;
  }
}


static Int4 FindMaxOrgsInUniqBioSrcGrpList (ValNodePtr list)
{
  Int4 max = 0;
  UniqBioSrcGrpPtr g;

  while (list != NULL) {
    g = (UniqBioSrcGrpPtr) list->data.ptrvalue;
    if (g != NULL && g->num_biop > max) {
      max = g->num_biop;
    }
    list = list->next;
  }
  return max;
}


static Int4 CountUniqueOrgsInUniqBioSrcGrpList (ValNodePtr list)
{
  Int4 count = 0;
  UniqBioSrcGrpPtr g;

  while (list != NULL) {
    g = (UniqBioSrcGrpPtr) list->data.ptrvalue;
    if (g != NULL && g->num_biop == 1) {
      count++;
    }
    list = list->next;
  }
  return count;
}    


static Boolean AddQualToUniqBioSrcGrpList (ValNodePtr list, FieldTypePtr field)
{
  Boolean    rval = FALSE;
  ValNodePtr vnp;

  vnp = list;
  while (vnp != NULL) {
    rval |= AddQualToUniqBioSrcGrp (vnp->data.ptrvalue, field);
    vnp = vnp->next;
  }
  if (rval) {
    /* regroup */
    ReGroupUniqBioSrcGrpList (list);
  }
  return rval;
}


typedef struct qualcombo {
  Int4         num_groups;
  Int4         num_mods;
  Int4         max_orgs_in_group;
  Int4         num_unique_orgs;
  FieldTypePtr field_list;
  ValNodePtr   group_list;
} QualComboData, PNTR QualComboPtr;


/* This function creates a new ModifierCombination item using the supplied
 * OrgGroup list.  It calculates the number of groups, maximum number of 
 * organisms in any one group, and number of unique organisms.
 * Initially there are no modifiers.
 */
static QualComboPtr QualComboNew (ValNodePtr grp_list)
{
  QualComboPtr newm;

  newm = (QualComboPtr) MemNew (sizeof (QualComboData));
  if (newm == NULL) return NULL;

  newm->num_mods = 0;
  newm->field_list = NULL;
  
  /* copy groups */
  newm->group_list = UniqBioSrcGrpListCopy (grp_list);
 
  newm->max_orgs_in_group = FindMaxOrgsInUniqBioSrcGrpList (newm->group_list);
  newm->num_unique_orgs = CountUniqueOrgsInUniqBioSrcGrpList (newm->group_list);
  newm->num_groups = ValNodeLen (newm->group_list);

  return newm; 
}


/* The CopyQualCombo creates a copy of a QualCombo item.
 * This includes creating a copy of the number and list of modifiers
 * and a copy of the number and list of OrgGroups, as well as copying the
 * maximum number of organisms in any one group and the number of unique
 * organism descriptions produced by this combination of modifiers.
 */
static QualComboPtr QualComboCopy (
  QualComboPtr m
)
{
  QualComboPtr newm;

  newm = QualComboNew (m->group_list);
  if (newm == NULL) return NULL;

  newm->field_list = FieldTypeListCopy (m->field_list);
  newm->num_mods = m->num_mods;
  
  return newm; 
}

/* This function frees the memory associated with a list of
 * ModifierCombination items.
 */
static QualComboPtr QualComboFree (
  QualComboPtr m
)
{
  if (m != NULL) {
    m->group_list = UniqBioSrcGrpListFree (m->group_list);
    m->field_list = FieldTypeListFree (m->field_list);
    m = MemFree (m);
  }
  return m;
}


static void TESTDisplayQualCombo (QualComboPtr q)
{
  ValNodePtr vnp_t, vnp_b, vnp_q, vnp_f;
  UniqBioSrcGrpPtr g;
  UniqBioSourcePtr bio;
  Int4             val;

  for (vnp_t = q->group_list; vnp_t != NULL; vnp_t = vnp_t->next) {
    g = (UniqBioSrcGrpPtr) vnp_t->data.ptrvalue;
    for (vnp_b = g->biop_list; vnp_b != NULL; vnp_b = vnp_b->next) {
      bio = (UniqBioSourcePtr) vnp_b->data.ptrvalue;
      for (vnp_q = bio->available_fields; vnp_q != NULL; vnp_q = vnp_q->next) {
        vnp_f = (ValNodePtr) vnp_q->data.ptrvalue;
        val = vnp_f->data.intvalue;
      }
    }
  }
}


static Boolean AddQualToQualCombo (
  QualComboPtr m,
  FieldTypePtr field
)
{
  Boolean    rval = FALSE;

  if (m == NULL || field == NULL) return rval;
  
  
  if (AddQualToUniqBioSrcGrpList (m->group_list, field)) {
    ValNodeLink (&(m->field_list), AsnIoMemCopy (field, (AsnReadFunc) FieldTypeAsnRead, (AsnWriteFunc) FieldTypeAsnWrite));
    m->field_list = ValNodeSort (m->field_list, SortFieldTypeByImportance);
    m->num_mods ++;
    m->group_list = BioSrcGrpListSort (m->group_list);
    m->max_orgs_in_group = FindMaxOrgsInUniqBioSrcGrpList (m->group_list);
    m->num_unique_orgs = CountUniqueOrgsInUniqBioSrcGrpList (m->group_list);
    m->num_groups = ValNodeLen (m->group_list);
    rval = TRUE;
  }
  return rval;
}


static ValNodePtr LIBCALLBACK QualComboListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = QualComboFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


/* NOTE - we want to sort groups from most unique organisms to least unique organisms */
/* secondary sort - most groups to least groups */
/* tertiary sort - fewer max orgs in group to most max orgs in group */
/* fourth sort - least mods to most mods */
static int LIBCALLBACK SortQualCombo (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr    vnp1, vnp2;
  QualComboPtr  g1, g2;
  FieldTypePtr  field1, field2;
  int           rval = 0;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL && vnp1->data.ptrvalue != NULL && vnp2->data.ptrvalue != NULL) {
      g1 = (QualComboPtr) vnp1->data.ptrvalue;
      g2 = (QualComboPtr) vnp2->data.ptrvalue;
      if (g1->num_unique_orgs > g2->num_unique_orgs) {
        rval = -1;
      } else if (g1->num_unique_orgs < g2->num_unique_orgs) {
        rval = 1;
      } else if (g1->num_groups > g2->num_groups) {
        rval = -1;
      } else if (g1->num_groups < g2->num_groups) {
        rval = 1;
      } else if (g1->max_orgs_in_group < g2->max_orgs_in_group) {
        rval = -1;
      } else if (g1->max_orgs_in_group > g2->max_orgs_in_group) {
        rval = 1;
      } else if (g1->num_mods < g2->num_mods) {
        rval = -1;
      } else if (g1->num_mods > g2->num_mods) {
        rval = 1;
      } else {
        /* compare modifiers */
        field1 = g1->field_list;
        field2 = g2->field_list;
        while (field1 != NULL && field2 != NULL 
               && (rval = CompareFieldTypeByImportance (field1, field2)) == 0) {
          field1 = field1->next;
          field2 = field2->next;
        }
        if (rval == 0) {
          if (field1 == NULL && field2 != NULL) {
            rval = -1;
          } else if (field1 != NULL && field2 == NULL) {
            rval = 1;
          }
        }
      }
    }
  }
  return rval;
}


static ValNodePtr QualComboListSort (ValNodePtr orig)
{
  orig = ValNodeSort (orig, SortQualCombo);
  ValNodeUnique (&orig, SortQualCombo, QualComboListFree);
  return orig;
}



static ValNodePtr ExpandOneComboListUsingAllPresentQuals (QualComboPtr q)
{
  ValNodePtr new_list = NULL, vnp, vnp_m;
  FieldTypePtr match_list;
  UniqBioSrcGrpPtr g;
  QualComboPtr q_new;
  Boolean      found_group_improvement = FALSE;

  if (q == NULL) return NULL;
  for (vnp = q->group_list; vnp != NULL && !found_group_improvement; vnp = vnp->next) {
    g = (UniqBioSrcGrpPtr) vnp->data.ptrvalue;
    if (g->num_biop == 1) break;
    match_list = GetAllPresentQualsForGroup (g);
    for (vnp_m = match_list; vnp_m != NULL; vnp_m = vnp_m->next) {
      q_new = QualComboCopy (q);
      if (AddQualToQualCombo (q_new, vnp_m) && q_new->num_groups > q->num_groups) {
        ValNodeAddPointer (&new_list, 0, q_new);
        found_group_improvement = TRUE;
      } else {
        q_new = QualComboFree (q_new);
        RemoveFieldFromUniqBioSrcGrp (g, vnp_m);
      }
    }
    match_list = FieldTypeListFree (match_list);
  }
  return new_list;
}


static void TESTDisplayList (ValNodePtr new_list)
{
  QualComboPtr q;
  ValNodePtr vnp;

  for (vnp = new_list; vnp != NULL; vnp = vnp->next) {
    q = (QualComboPtr) vnp->data.ptrvalue;
    TESTDisplayQualCombo (q);
  }

}


static Boolean IsQualOkForDefline (ValNodePtr vnp)
{
  ValNodePtr scp;
  
  if (vnp == NULL || vnp->choice != FieldType_source_qual
      || (scp = (SourceQualChoicePtr) vnp->data.ptrvalue) == NULL
      || scp->choice != SourceQualChoice_textqual) {
    return FALSE;
  }
  if (scp->data.intvalue == Source_qual_map || scp->data.intvalue == Source_qual_nat_host) {
    return FALSE;
  } else {
    return TRUE;
  }
}
  
      
    
static ValNodePtr ExpandOneComboListUsingAnyPresentQuals (QualComboPtr q)
{
  ValNodePtr new_list = NULL, vnp, vnp_m;
  FieldTypePtr match_list;
  UniqBioSrcGrpPtr g;
  QualComboPtr q_new;
  Boolean      found_group_improvement = FALSE;

  if (q == NULL) return NULL;
  for (vnp = q->group_list; vnp != NULL && !found_group_improvement; vnp = vnp->next) {
    g = (UniqBioSrcGrpPtr) vnp->data.ptrvalue;
    if (g->num_biop == 1) break;
    match_list = GetAllQualsForGroup (g);
    for (vnp_m = match_list; vnp_m != NULL; vnp_m = vnp_m->next) {
      if (!IsQualOkForDefline(vnp_m)) {
        RemoveFieldFromUniqBioSrcGrp (g, vnp_m);
      } else {
        q_new = QualComboCopy (q);
        if (AddQualToQualCombo (q_new, vnp_m) && q_new->num_groups > q->num_groups) {
          ValNodeAddPointer (&new_list, 0, q_new);
          found_group_improvement = TRUE;
        } else {
          q_new = QualComboFree (q_new);
          RemoveFieldFromUniqBioSrcGrp (g, vnp_m);
        }
      }
    }
    match_list = FieldTypeListFree (match_list);
  }

  return new_list;
}


static Boolean ExpandComboList (ValNodePtr PNTR list)
{
  QualComboPtr  q;
  ValNodePtr    new_list, vnp, vnp_next, prev = NULL;
  Boolean       any_expansion = FALSE;

  if (*list == NULL) return FALSE;
  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    q = (QualComboPtr) vnp->data.ptrvalue;
    new_list = ExpandOneComboListUsingAnyPresentQuals (q);

    if (new_list == NULL) {
      prev = vnp;
    } else {
      if (prev == NULL) {
        *list = new_list; 
      } else {
        prev->next = new_list;
      }
      prev = new_list;
      while (prev->next != NULL) {
        prev = prev->next;
      }
      ValNodeLink (&new_list, vnp->next);
      vnp->next = NULL;
      vnp = QualComboListFree (vnp);
      any_expansion = TRUE;
    }
    vnp = vnp_next;
  }
  return any_expansion;
}


static void BuildUniqBioSrcList (
  BioSourcePtr biop,
  Pointer userdata
)
{
  UniqBioSourcePtr u;
  ValNodeBlockPtr  vnbp;
  ValNodePtr       vnp;

  u = UniqBioSourceNew (biop);
  vnbp = (ValNodeBlockPtr) userdata;
  vnp = ValNodeAddPointer (&vnbp->tail, 0, u);
  if (vnbp->head == NULL) {
    vnbp->head = vnp;
  }
  vnbp->tail = vnp;
}


/* The function FindBestQualCombo tries to find the best combination of modifiers
 * to create unique organism descriptions.  This is accomplished by
 * creating a list of required modifiers, and then creating a list of
 * combinations of modifiers by adding modifiers one at a time
 * to see if the additional modifiers provide any more differentiation in
 * the list.
 * In order to do this, I start with a list of required modifiers, and 
 * then create copies of this list.  For each copy I add one of the modifiers
 * that are present in the bio sources and not already on the list.
 * If adding the modifier increases the differentiation, I add that copy to
 * the list of possible combinations, otherwise I discard it.
 * The function then makes copies of all of the new items added to the list,
 * starting with the item pointed to by start_of_expand, and adds another
 * modifier to each combination, keeping the combinations that increase
 * the differentiation and discarding the rest.
 * This process continues until I have a combination that produces completely
 * differentiated bio sources, or I run out of possible combinations.
 * If the list of possible combinations is exhausted before each organism
 * has a unique description, the function selects the combination from the
 * list with the largest number of unique organism descriptions.  If more 
 * than one combination produces the largest number of unique organisms,
 * the combination with the largest number of unique organisms and the
 * largest number of groups will be selected.
 */
static QualComboPtr FindBestQualComboEx(ValNodePtr PNTR biop_list, ModifierItemLocalPtr ItemList)
{
  QualComboPtr initial_combo = NULL, best_combo = NULL;
  ValNodePtr   group_list = NULL, combo_list = NULL;
  UniqBioSrcGrpPtr g;
  SourceQualChoice scd;
  FieldType ft;
  Int4      i, qual;

  if (biop_list == NULL || *biop_list == NULL) {
    return NULL;
  }

  /* sort organisms */
  *biop_list = UniqBioSourceListSort (*biop_list);

  /* create group list */
  g = UniqBioSrcGrpNew (*biop_list);
  ValNodeAddPointer (&group_list, 0, g);  

  ReGroupUniqBioSrcGrpList (group_list);
  group_list = BioSrcGrpListSort (group_list);

  /* create combo with just the org groups */
  initial_combo = QualComboNew (group_list);
  group_list = UniqBioSrcGrpListFree (group_list);
  if (initial_combo == NULL) return NULL;

  /* add required quals */
  ft.choice = FieldType_source_qual;
  ft.data.ptrvalue = &scd;
  ft.next = NULL;
  scd.choice = SourceQualChoice_textqual;
  if (ItemList == NULL) {
    /* endogenous virus name */
    scd.data.intvalue = Source_qual_endogenous_virus_name;
    AddQualToQualCombo (initial_combo, &ft);
    /* plasmid name */
    scd.data.intvalue = Source_qual_plasmid_name;
    AddQualToQualCombo (initial_combo, &ft);
    /* transgenic */
    scd.data.intvalue = Source_qual_transgenic;
    AddQualToQualCombo (initial_combo, &ft);
  } else {
    for (i = 0; i < numDefLineModifiers; i++) {
      if (ItemList[i].required && ItemList[i].any_present) {
        qual = GetSrcQualFromSubSrcOrOrgMod (DefLineModifiers[i].subtype, DefLineModifiers[i].isOrgMod);
        if (qual > -1) {
          scd.data.intvalue = qual;
          AddQualToQualCombo (initial_combo, &ft);
        }
      }
    }
  }  

  if (initial_combo->max_orgs_in_group == 1)
  {
    /* we're done - they're all unique */
    return initial_combo;
  }

  /* they're not unique yet.  Need to find a combination of modifiers that will make this as unique as possible */
  ValNodeAddPointer (&combo_list, 0, initial_combo);
  best_combo = initial_combo;
  while (ExpandComboList (&combo_list)) {
    /* sort after expansion */  
    combo_list = QualComboListSort (combo_list);
    best_combo = combo_list->data.ptrvalue;
    if (best_combo->max_orgs_in_group == 1) {
      break;
    } 
  }
  best_combo = QualComboCopy (best_combo);
  combo_list = QualComboListFree (combo_list);
  return best_combo;
}


static QualComboPtr FindBestQualCombo(SeqEntryPtr sep, ModifierItemLocalPtr ItemList)
{
  QualComboPtr  best_combo;
  ValNodeBlock  vnb;

  /* first, get list of organisms */
  vnb.head = NULL;
  vnb.tail = NULL;
  VisitBioSourcesInSep (sep, &vnb, BuildUniqBioSrcList);

  best_combo = FindBestQualComboEx (&vnb.head, ItemList);

  UniqBioSourceListFree (vnb.head);
  return best_combo;
}


static ModifierCombinationPtr ModifierCombinationFromQualCombo (QualComboPtr q)
{
  ModifierCombinationPtr m;
  FieldTypePtr           field;
  Int4                   i;

  if (q == NULL) {
    return NULL;
  }

  m = (ModifierCombinationPtr) MemNew (sizeof (ModifierCombinationData));
  m->num_groups = q->num_groups;
  m->num_mods = q->num_mods;
  m->max_orgs_in_group = q->max_orgs_in_group;
  m->num_unique_orgs = q->num_unique_orgs;
  m->next = NULL;
  m->group_list = NULL;
  m->modifier_indices = NULL;
  for (field = q->field_list; field != NULL; field = field->next) {
    i = GetDeflinePosForFieldType (field);
    if (i >= 0) {
      ValNodeAddInt (&(m->modifier_indices), 0, i);
    }
  }
  return m;
}


NLM_EXTERN ValNodePtr FindBestModifiersForDeflineClauseList (
  ValNodePtr defline_clauses,
  ModifierItemLocalPtr ItemList
)

{
  QualComboPtr best_combo;
  ValNodePtr   biop_list = NULL, vnp;
  DefLineFeatClausePtr df;
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  UniqBioSourcePtr  u;
  ModifierCombinationPtr m;
  ValNodePtr modifier_indices = NULL;

  /* first, create list of organisms */
  for (vnp = defline_clauses; vnp != NULL; vnp = vnp->next) {
    df = (DefLineFeatClausePtr) vnp->data.ptrvalue;
    if (df != NULL) {
      sdp = SeqMgrGetNextDescriptor (df->bsp, NULL, Seq_descr_source, &context);
      if (sdp != NULL && sdp->data.ptrvalue != NULL) {
        u = UniqBioSourceNew (sdp->data.ptrvalue);
        ValNodeAddPointer (&(u->strings), 0, StringSave (df->clauselist));
        ValNodeAddPointer (&biop_list, 0, u);
      }
    }
  }

  best_combo = FindBestQualComboEx (&biop_list, ItemList);

  biop_list = UniqBioSourceListFree (biop_list);
  m = ModifierCombinationFromQualCombo (best_combo);
  if (m != NULL) {
    modifier_indices = CopyModifierIndices (m->modifier_indices);
    FreeModifierCombo (m);
  }
  return modifier_indices;
}


NLM_EXTERN ValNodePtr FindBestModifiersEx(
  SeqEntryPtr sep,
  ModifierItemLocalPtr ItemList,
  Boolean use_new
)

{
  ModifierCombinationPtr m;
  QualComboPtr q;
  ValNodePtr modifier_indices = NULL;

  if (use_new) {
    q = FindBestQualCombo (sep, ItemList);
    m = ModifierCombinationFromQualCombo (q);
    q = QualComboFree (q);
  } else {
    m = FindBestCombo (sep, ItemList);
  }
  modifier_indices = CopyModifierIndices (m->modifier_indices);
  FreeModifierCombo (m);
  return modifier_indices;
}


NLM_EXTERN ValNodePtr FindBestModifiers(
  SeqEntryPtr sep,
  ModifierItemLocalPtr ItemList
)

{
  return FindBestModifiersEx (sep, ItemList, FALSE);
}


/* In this test function, we create a list of biosources with various combinations of modifiers,
 * and then calculate the best combination to use for the organism description.
 */
static CharPtr strings1[] = {"a", "b", "c"};
static CharPtr strings2[] = {"foo", "bar", "baz"};
static CharPtr strings3[] = {"d", "e", "f"};

static void SetBiopQual (BioSourcePtr biop, Int4 qual, CharPtr val)
{
  OrgModPtr mod;
  SubSourcePtr ssp;

  if (DefLineModifiers[qual].isOrgMod)
  {
    mod = OrgModNew ();
    mod->subtype = DefLineModifiers[qual].subtype;
    mod->subname = StringSave (val);
    mod->next = biop->org->orgname->mod;
    biop->org->orgname->mod = mod;
  } else {
    ssp = SubSourceNew ();
    ssp->subtype = DefLineModifiers[qual].subtype;
    ssp->name = StringSave (val);
    ssp->next = biop->subtype;
    biop->subtype = ssp;
  }
}


static void ClearBiopQuals (BioSourcePtr biop)
{
  biop->org->orgname->mod = OrgModSetFree (biop->org->orgname->mod);
  biop->subtype = SubSourceSetFree (biop->subtype);
}


static void PrintBiopQuals (BioSourcePtr biop, FILE *fp)
{
  OrgModPtr mod;
  SubSourcePtr ssp;

  fprintf (fp, "Taxname: %s", biop->org->taxname);
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next) {
    fprintf (fp, "\tOrgMod%d:%s", mod->subtype, mod->subname);
  }
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    fprintf (fp, "\tSubSource%d:%s", ssp->subtype, ssp->name);
  }
  fprintf (fp, "\n");
}


static void PrintModifiers (ValNodePtr modifiers, FILE *fp)
{
  ValNodePtr vnp;

  if (modifiers == NULL) {
    fprintf (fp, "\tNo combo");
  } 
  for (vnp = modifiers; vnp != NULL; vnp = vnp->next) {
    fprintf (fp, "\t%s:%d", DefLineModifiers[vnp->data.intvalue].isOrgMod ? "OrgMod" : "SubSource", 
                            DefLineModifiers[vnp->data.intvalue].subtype);
  }
  fprintf (fp, "\n");
}


static Boolean IsNonTextDeflineQual (Int4 srcqual)
{
  if (srcqual == DEFLINE_POS_Transgenic
      || srcqual == DEFLINE_POS_Germline
      || srcqual == DEFLINE_POS_Metagenomic
      || srcqual == DEFLINE_POS_Environmental_sample
      || srcqual == DEFLINE_POS_Rearranged)
  {
    return TRUE;  
  }
  else
  {
    return FALSE;
  }
}


static void CreateOneTest (FILE *fp, Int4 i, Int4 j, Int4 k, BioSourcePtr PNTR biops, Int4 num_biops, Boolean vary1, Boolean vary2, Boolean vary3)
{
  Int4 n;
  ValNodePtr uniq_biop_list = NULL;
  QualComboPtr q;
  ModifierCombinationPtr m;

  for (n = 0; n < num_biops; n++) {
    if (i < numDefLineModifiers) {
      if (vary1) {
        SetBiopQual (biops[n], i, strings1[n]);
      } else {
        SetBiopQual (biops[n], i, strings1[0]);
      }
    }
    if (j < numDefLineModifiers) {
      if (vary2) {
        SetBiopQual (biops[n], j, strings2[n]);
      } else {
        SetBiopQual (biops[n], j, strings2[0]);
      }
    }
    if (k < numDefLineModifiers) {
      if (vary3) {
        SetBiopQual (biops[n], k, strings3[n]);
      } else {
        SetBiopQual (biops[n], k, strings3[0]);
      }
    }
    ValNodeAddPointer (&uniq_biop_list, 0, UniqBioSourceNew (biops[n]));
  }
  q = FindBestQualComboEx (&uniq_biop_list, NULL);
  m = ModifierCombinationFromQualCombo (q);
  /* print results */
  for (n = 0; n < num_biops; n++) {
    PrintBiopQuals (biops[n], fp);
  }
  PrintModifiers (m->modifier_indices, fp);
  q = QualComboFree (q);
  FreeModifierCombo (m);
  uniq_biop_list = UniqBioSourceListFree (uniq_biop_list);

  /* clear quals */
  for (n = 0; n < num_biops; n++) {
    ClearBiopQuals (biops[n]);
  }
  
}


extern void TestFindBestQualCombo (FILE *fp)
{
  BioSourcePtr biops[3];
  Int4         num_biops = 3;
  Int4         i, j, k;

  for (i = 0; i < num_biops; i++) {
    biops[i] = BioSourceNew ();
    biops[i]->org = OrgRefNew ();
    biops[i]->org->orgname = OrgNameNew ();
  }
  
  /* first try with all organism names the same */
  for (i = 0; i < num_biops; i++) {
    biops[i]->org->taxname = StringSave ("Homo sapiens");
  }

  for (i = 0; i <= numDefLineModifiers; i++) {
    if (IsNonTextDeflineQual (i)) {
      continue;
    }
    for (j = i + 1; j <= numDefLineModifiers; j++) {
      if (IsNonTextDeflineQual (j)) {
        continue;
      }
      for (k = j + 1; k <= numDefLineModifiers; k++) {
        if (IsNonTextDeflineQual (k)) {
          continue;
        }
        if (k != numDefLineModifiers) {
          /* try all the same but 1*/
          CreateOneTest (fp, i, j, k, biops, num_biops, FALSE, FALSE, TRUE);
        }
        if (j != numDefLineModifiers || k != numDefLineModifiers) {
          /* try 2 different 1 same*/
          CreateOneTest (fp, i, j, k, biops, num_biops, FALSE, TRUE, TRUE);
        }
        if (i != numDefLineModifiers || j != numDefLineModifiers || k != numDefLineModifiers) {
          /* try all different */
          CreateOneTest (fp, i, j, k, biops, num_biops, TRUE, TRUE, TRUE);
        }
      }
    }
  }
  for (i = 0; i < num_biops; i++) {
    biops[i] = BioSourceFree (biops[i]);
  }
 
}




/* collection_date has a controlled format.  
 * It is YYYY or Mmm-YYYY or DD-Mmm-YYYY where Mmm = Jan, Feb, Mar, Apr, May, 
 *                                                   Jun, Jul, Aug, Sep, Oct, 
 *                                                   Nov, Dec
 * This function will convert other formats  to this format.
 * For instance, September 12, 2004 should be converted to 12-Sep-2004
 * 12/15/2003 should be converted to 15-Dec-2003.  
 * 
 * If the date supplied is ambiguous (01/03/05), can you allow the indexer to choose which field goes in Mmm and which in DD.
 */

NLM_EXTERN Int4 ReadNumberFromToken (CharPtr token, Int4 token_len)
{
  Int4 val = 0;
  
  if (token == NULL || !isdigit (*token))
  {
    return val;
  }
  while (token_len > 0)
  {
    val *= 10;
    val += *token - '0';
    token++;
    token_len--;
  }
  
  return val;
}

static Int4 GetYearFromNumber(Int4 year)
{
    Nlm_DayTime dt;

  if (year < 1000)
  {
    GetDayTime (&dt);
    if (year + 2000 > dt.tm_year + 1901)
    {
      year += 1900;
    }
    else
    {
      year += 2000;
    }
  }
  return year;
}

NLM_EXTERN Int4 GetYearFromToken (CharPtr token, Int4 token_len)
{
  Int4        year = 0;
  
  if (token == NULL || token_len == 0 || token_len > 4)
  {
    return 0;
  }
  
  year = GetYearFromNumber(ReadNumberFromToken (token, token_len));
  
  return year;
}

static CharPtr month_abbrevs [12] =
{
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", 
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};


NLM_EXTERN CharPtr GetMonthAbbrev (Int4 n)
{
  if (n > 0 && n <= 12) {
    return month_abbrevs[n - 1];
  } else {
    return NULL;
  }
}
  

static Int4 days_in_month [12] =
{
  31, 29, 31, 30, 31, 30,
  31, 31, 30, 31, 30, 31
};


NLM_EXTERN Int4 GetDaysInMonth (Int4 n)
{
  if (n > 0 && n <= 12) {
    return days_in_month[n - 1];
  } else {
    return 0;
  }
}


NLM_EXTERN Int4 GetMonthNumFromAbbrev (CharPtr month_abbrev) 
{
  Int4 i;

  for (i = 0; i < 12; i++) {
    if (StringICmp (month_abbrev, month_abbrevs[i]) == 0) {
      return i;
    }
  }
  return -1;
}

static Int4 GetDaysInMonthByName (CharPtr month)
{
  Int4 month_num;
  
  for (month_num = 0; month_num < 12; month_num++)
  {
    if (StringCmp (month, month_abbrevs [month_num]) == 0)
    {
      return days_in_month [month_num];
    }
  }
  return 0;
}

NLM_EXTERN CharPtr GetMonthFromToken (CharPtr token, Int4 token_len)
{
  Int4    month_num;
  
  if (token == NULL || token_len == 0)
  {
    return NULL;
  }
  
  if (isdigit (*token)) 
  {
    if (token_len > 2)
    {
      return NULL;
    }
    else
    {
      month_num = ReadNumberFromToken (token, token_len);
      if (month_num == 0 || month_num > 12)
      {
        return NULL;
      }
      else
      {
        return month_abbrevs [month_num - 1];
      }
    }
  }
  else
  {
    for (month_num = 0; month_num < 12; month_num++)
    {
      if (StringNICmp (token, month_abbrevs[month_num], 3) == 0)
      {
        return month_abbrevs[month_num];
      }
    }
    return NULL;
  }
}

static Boolean
ChooseDayAndYear 
(Int4    num_1,
 Int4    num_2,
 CharPtr month,
 Boolean year_first,
 Int4Ptr day,
 Int4Ptr year)
{  
  if (day == NULL || year == NULL)
  {
    return FALSE;
  }
  
  if (num_1 == 0 && num_2 == 0)
  {
    return FALSE;
  }
  else if (num_1 == 0)
  {
    *year = 2000;
    *day = num_2;
  }
  else if (num_2 == 0)
  {
    *year = 2000;
    *day = num_1;
  }
  else if (num_1 > GetDaysInMonthByName (month))
  {
    *year = num_1;
    *day = num_2;
  }
  else if (num_2 > GetDaysInMonthByName (month))
  {
    *year = num_2;
    *day = num_1;
  }
  else if (year_first)
  {
    *year = num_1;
    *day = num_2;
  }
  else
  {
    *year = num_2;
    *day = num_1;
  }
  
  return TRUE;
}

static Boolean 
ChooseMonthAndYear
(Int4    num_1,
 Int4    num_2,
 Boolean month_first,
 CharPtr PNTR month,
 Int4Ptr year,
 BoolPtr month_ambiguous)
{
  if (year == NULL || month == NULL 
      || (num_1 == 0 && num_2 == 0)
      || (num_1 > 12 && num_2 > 12)
      || (num_1 == 0 && num_2 > 12)
      || (num_2 == 0 && num_1 > 12))
  {
    return FALSE;
  }
  
  if (num_1 == 0)
  {
    *year = 2000;
    *month = month_abbrevs[num_2 - 1];
  }
  else if (num_2 == 0)
  {
    *year = 2000;
    *month = month_abbrevs[num_1 - 1];
  }
  else if (num_1 > 12)
  {
    *year = GetYearFromNumber(num_1);
    *month = month_abbrevs [num_2 - 1];
  }
  else if (num_2 > 12)
  {
    *year = GetYearFromNumber(num_2);
    *month = month_abbrevs [num_1 - 1];
  }
  else if (month_first)
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *year = GetYearFromNumber(num_2);
    *month = month_abbrevs [num_1 - 1];
  }
  else
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *year = GetYearFromNumber(num_1);
    *month = month_abbrevs [num_2 - 1];
  }
  return TRUE;
}


static Boolean ChooseMonthAndDay 
(Int4    num_1,
 Int4    num_2,
 Boolean month_first,
 CharPtr PNTR month,
 Int4Ptr day,
 BoolPtr month_ambiguous)
{
  if (day == NULL || month == NULL || num_1 == 0 || num_2 == 0
      || (num_1 > 12 && num_2 > 12))
  {
    return FALSE;
  }
  
  if (num_1 > 12)
  {
    *day = num_1;
    *month = month_abbrevs [num_2 - 1];
  }
  else if (num_2 > 12)
  {
    *day = num_2;
    *month = month_abbrevs [num_1 - 1];
  }
  else if (month_first)
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *day = num_2;
    *month = month_abbrevs [num_1 - 1];
  }
  else
  {
    if (month_ambiguous != NULL) 
    {
      *month_ambiguous = TRUE;
    }
    *day = num_1;
    *month = month_abbrevs [num_2 - 1];
  }
  return TRUE;
}

NLM_EXTERN CharPtr ReformatDateStringEx (CharPtr orig_date, Boolean month_first, BoolPtr month_ambiguous)
{
  CharPtr reformatted_date = NULL, cp;
  Int4    year = 0, day = 0;
  CharPtr month = NULL;
  CharPtr token_list[3];
  Int4    token_lens[3];
  CharPtr numbers = "0123456789";
  CharPtr letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  Int4    num_tokens = 0;
  Int4    token_len;
  Int4    month_token = -1;
  Boolean is_num;
  Int4    num_1, num_2, num_3;
  
  if (StringHasNoText (orig_date))
  {
    return NULL;
  }
  
  /* divide our original date into tokens */
  /* skip over any leading spaces */
  cp = orig_date;
  while (*cp != 0 && num_tokens < 3)
  {
    is_num = FALSE;
    token_len = StringSpn (cp, numbers);  
    if (token_len == 0)
    {
      token_len = StringSpn (cp, letters);
    }
    else
    {
      is_num = TRUE;
    }
    if (token_len == 0)
    {
      cp++;
    }
    else
    {
      if (!is_num)
      {
        if (month_token == -1)
        {
          month_token = num_tokens;
        }
        else
        {
          /* already found a month string */
          return NULL;
        }
      }
      token_list [num_tokens] = cp;
      token_lens [num_tokens] = token_len;
      num_tokens ++;
      cp += token_len;
    }
  }
 
  if (num_tokens == 0 || *cp != 0)
  {
    return NULL;
  }

  if (num_tokens == 1)
  {
    if (month_token == 0)
    {
      return NULL;
    }
    year = GetYearFromToken (token_list [0], token_lens [0]);
  }
  else if (num_tokens == 2)
  {
    if (month_token == 0)
    {
      month = GetMonthFromToken (token_list [0], token_lens [0]);
      year = GetYearFromToken (token_list [1], token_lens [1]);
    }
    else if (month_token == 1)
    {
      month = GetMonthFromToken (token_list [1], token_lens [1]);
      year = GetYearFromToken (token_list [0], token_lens [0]);
    }
    else
    {
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      if (! ChooseMonthAndYear (num_1, num_2, month_first, &month, &year, month_ambiguous))
      {
        return NULL;
      }
    }
  }
  else if (num_tokens == 3)
  {
    if (month_token == 0)
    {
      month = GetMonthFromToken (token_list [0], token_lens [0]);
      num_1 = ReadNumberFromToken (token_list [1], token_lens [1]);
      num_2 = ReadNumberFromToken (token_list [2], token_lens [2]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else if (month_token == 1)
    {
      month = GetMonthFromToken (token_list [1], token_lens [1]);
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [2], token_lens [2]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else if (month_token == 2)
    {
      month = GetMonthFromToken (token_list [2], token_lens [2]);
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else
    {
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      num_3 = ReadNumberFromToken (token_list [2], token_lens [2]);
      
      if (num_1 > 31 || num_1 == 0)
      {
        year = num_1;
        if (! ChooseMonthAndDay (num_2, num_3, month_first, &month, &day, month_ambiguous))
        {
          return NULL;
        }
      }
      else if (num_2 > 31 || num_2 == 0)
      {
        year = num_2;
        if (! ChooseMonthAndDay (num_1, num_3, month_first, &month, &day, month_ambiguous))
        {
          return NULL;
        }
      }
      else if (num_3 > 31 || num_3 == 0)
      {
        year = num_3;
        if (! ChooseMonthAndDay (num_1, num_2, month_first, &month, &day, month_ambiguous))
        {
          return NULL;
        }
      }
      else if (num_1 > 0 && num_1 < 13 && num_2 > days_in_month [num_1] && num_3 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_1 - 1];
        year = num_2;
        day = num_3;
      }
      else if (num_1 > 0 && num_1 < 13 && num_3 > days_in_month [num_1] && num_2 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_1 - 1];
        year = num_3;
        day = num_2;
      }
      else if (num_2 > 0 && num_2 < 13 && num_1 > days_in_month [num_2] && num_3 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_2 - 1];
        year = num_1;
        day = num_3;
      }
      else if (num_2 > 0 && num_2 < 13 && num_3 > days_in_month [num_2] && num_1 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_2 - 1];
        year = num_3;
        day = num_1;
      }
      else if (num_3 > 0 && num_3 < 13 && num_1 > days_in_month [num_3] && num_2 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_3 - 1];
        year = num_1;
        day = num_2;
      }
      else if (num_3 > 0 && num_3 < 13 && num_2 > days_in_month [num_3] && num_1 <= days_in_month [num_1])
      {
        month = month_abbrevs [num_3 - 1];
        year = num_2;
        day = num_1;
      }
      else
      {
        year = num_3;
        if (! ChooseMonthAndDay (num_1, num_2, month_first, &month, &day, month_ambiguous))
        {
          year = num_1;
          if (!ChooseMonthAndDay (num_2, num_3, month_first, &month, &day, month_ambiguous))
          {
            return NULL;
          }
        }
      }
                
    }
    year = GetYearFromNumber(year);
  }
  
  if (month == NULL && day > 0)
  {
    return NULL;
  }
  
  reformatted_date = (CharPtr) MemNew (sizeof (Char) * 12);
  if (reformatted_date == NULL)
  {
    return NULL;
  }
   
  if (month == NULL)
  {
    sprintf (reformatted_date, "%d", year);
  }
  else if (day == 0)
  {
    sprintf (reformatted_date, "%s-%d", month, year);
  }
  else
  {
    sprintf (reformatted_date, "%02d-%s-%d", day, month, year);
  }
  return reformatted_date;
}


NLM_EXTERN CharPtr ReformatDateWithMonthNames (CharPtr orig_date)
{
  CharPtr reformatted_date = NULL, cp;
  Int4    year = 0, day = 0;
  CharPtr month = NULL;
  CharPtr token_list[3];
  Int4    token_lens[3];
  CharPtr numbers = "0123456789";
  CharPtr letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  Int4    num_tokens = 0;
  Int4    token_len;
  Int4    month_token = -1;
  Boolean is_num;
  Int4    num_1, num_2;
  
  if (StringHasNoText (orig_date))
  {
    return NULL;
  }
  
  /* divide our original date into tokens */
  /* skip over any leading spaces */
  cp = orig_date;
  while (*cp != 0 && num_tokens < 3)
  {
    is_num = FALSE;
    token_len = StringSpn (cp, numbers);  
    if (token_len == 0)
    {
      token_len = StringSpn (cp, letters);
    }
    else
    {
      is_num = TRUE;
    }
    if (token_len == 0)
    {
      cp++;
    }
    else
    {
      if (!is_num)
      {
        if (month_token == -1)
        {
          month_token = num_tokens;
        }
        else
        {
          /* already found a month string */
          return NULL;
        }
      }
      token_list [num_tokens] = cp;
      token_lens [num_tokens] = token_len;
      num_tokens ++;
      cp += token_len;
    }
  }
 
  if (num_tokens == 0 || *cp != 0 || month_token == -1 || num_tokens < 2)
  {
    return NULL;
  }

  if (num_tokens == 2)
  {
    if (month_token == 0)
    {
      month = GetMonthFromToken (token_list [0], token_lens [0]);
      year = GetYearFromToken (token_list [1], token_lens [1]);
    }
    else if (month_token == 1)
    {
      month = GetMonthFromToken (token_list [1], token_lens [1]);
      year = GetYearFromToken (token_list [0], token_lens [0]);
    }
    else
    {
      return NULL;
    }
  }
  else if (num_tokens == 3)
  {
    if (month_token == 0)
    {
      month = GetMonthFromToken (token_list [0], token_lens [0]);
      num_1 = ReadNumberFromToken (token_list [1], token_lens [1]);
      num_2 = ReadNumberFromToken (token_list [2], token_lens [2]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else if (month_token == 1)
    {
      month = GetMonthFromToken (token_list [1], token_lens [1]);
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [2], token_lens [2]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else if (month_token == 2)
    {
      month = GetMonthFromToken (token_list [2], token_lens [2]);
      num_1 = ReadNumberFromToken (token_list [0], token_lens [0]);
      num_2 = ReadNumberFromToken (token_list [1], token_lens [1]);
      if (!ChooseDayAndYear (num_1, num_2, month, FALSE, &day, &year))
      {
        return NULL;
      }
    }
    else
    {
      return NULL;                
    }
    year = GetYearFromNumber(year);
  }
  
  if (month == NULL && day > 0)
  {
    return NULL;
  }
  
  reformatted_date = (CharPtr) MemNew (sizeof (Char) * 12);
  if (reformatted_date == NULL)
  {
    return NULL;
  }
   
  if (month == NULL)
  {
    sprintf (reformatted_date, "%d", year);
  }
  else if (day == 0)
  {
    sprintf (reformatted_date, "%s-%d", month, year);
  }
  else
  {
    sprintf (reformatted_date, "%02d-%s-%d", day, month, year);
  }
  return reformatted_date;
}


NLM_EXTERN Boolean CreateMatPeptideFromCDS (SeqFeatPtr sfp)
{
  SeqFeatPtr        orig_prot, new_prot;
  SeqMgrFeatContext fcontext;
  SeqLocPtr         prot_loc;
  ProtRefPtr        prp;
  BioseqPtr         prot_bsp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }

  prot_bsp = BioseqFindFromSeqLoc (sfp->product);
  if (prot_bsp != NULL) {
    prot_loc = SeqLocIntNew (0, prot_bsp->length - 1, Seq_strand_plus, SeqLocId (sfp->product));
    new_prot = CreateNewFeatureOnBioseq (prot_bsp, SEQFEAT_PROT, prot_loc);
    orig_prot = SeqMgrGetNextFeature (prot_bsp, NULL, 0, FEATDEF_PROT, &fcontext);
    if (orig_prot != NULL) {
      prp = AsnIoMemCopy (orig_prot->data.value.ptrvalue, (AsnReadFunc) ProtRefAsnRead, (AsnWriteFunc) ProtRefAsnWrite);
    } else {
      prp = ProtRefNew ();
    }
    prp->processed = 2;
    new_prot->data.value.ptrvalue = prp;
    return TRUE;
  } else {
    return FALSE;
  }
}


NLM_EXTERN Boolean ConvertCDSToMatPeptideForOverlappingCDS (SeqFeatPtr sfp, SeqFeatPtr top_cds, Boolean remove_original)
{
  BioseqPtr prot_bsp;
  CdRegionPtr crp;
  SeqLocPtr   prot_loc;
  SeqFeatPtr  new_sfp, orig_prot;
  SeqMgrFeatContext prot_context;
  ProtRefPtr prp;
  Int4 frame;
  Boolean rval = FALSE;
  Boolean partial5, partial3;

  if (sfp == NULL || top_cds == NULL || sfp->data.choice != SEQFEAT_CDREGION || top_cds->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }

  prot_bsp = BioseqFindFromSeqLoc (top_cds->product);
  if (prot_bsp != NULL)
  {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    prot_loc = dnaLoc_to_aaLoc(top_cds, sfp->location, TRUE, &frame, !partial3);
    if (prot_loc != NULL)
    {
      /* Create new feature on prot_bsp */
      new_sfp = CreateNewFeatureOnBioseq (prot_bsp, SEQFEAT_PROT, prot_loc);
      if (new_sfp != NULL)
      {
        prot_bsp = BioseqFindFromSeqLoc (sfp->product);
        orig_prot = SeqMgrGetNextFeature (prot_bsp, NULL, 0, FEATDEF_PROT, &prot_context);
        if (orig_prot == NULL) {
          prp = ProtRefNew ();
        } else {          
          prp = AsnIoMemCopy (orig_prot->data.value.ptrvalue, (AsnReadFunc) ProtRefAsnRead, (AsnWriteFunc) ProtRefAsnWrite);
        }
        prp->processed = 2;
        new_sfp->data.value.ptrvalue = prp;
        prot_bsp->idx.deleteme = TRUE;
        
        rval = TRUE;
      }
      /* mark old feature for deletion */
      sfp->idx.deleteme = TRUE;
    }
  }
  return rval;
}


NLM_EXTERN Boolean AutoConvertCDSToMiscFeat (SeqFeatPtr cds, Boolean remove_original)
{
  BioseqPtr bsp;
  SeqFeatPtr top_cds = NULL;
  SeqMgrFeatContext fcontext;
  Boolean rval = FALSE;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) {
    return FALSE;
  }

  bsp = BioseqFindFromSeqLoc (cds->location);

  /* find overlapping coding region */
  for (top_cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       top_cds != NULL && (top_cds == cds || SeqLocCompare (top_cds->location, cds->location) != SLC_B_IN_A);
       top_cds = SeqMgrGetNextFeature (bsp, top_cds, SEQFEAT_CDREGION, 0, &fcontext)) {
  }

  if (top_cds == NULL) {
    if (remove_original) {
      rval = FALSE;
    } else {
      rval = CreateMatPeptideFromCDS (cds);
    }
  } else {
    rval = ConvertCDSToMatPeptideForOverlappingCDS (cds, top_cds, remove_original);
  }
  if (rval) {
    /* have to remove CDS because ConvertFeature has already created a duplicate if the feature is going to be kept */
    cds->idx.deleteme = TRUE;
  }
  return rval;
}


NLM_EXTERN SeqEntryPtr GetBestSeqEntryForItem (ValNodePtr vnp)
{
  SeqFeatPtr       sfp;
  SeqDescrPtr      sdp;
  ObjValNodePtr    ovp;
  SeqEntryPtr      sep = NULL;
  BioseqPtr        bsp;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  if (vnp->choice == OBJ_SEQFEAT) {
    sfp = vnp->data.ptrvalue;
    sep = GetBestTopParentForData (sfp->idx.entityID, BioseqFindFromSeqLoc (sfp->location));
  } else if (vnp->choice == OBJ_SEQDESC) {
    sdp = vnp->data.ptrvalue;
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
        sep = SeqMgrGetSeqEntryForData (ovp->idx.parentptr);
      } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
        sep = GetBestTopParentForData (ovp->idx.entityID, ovp->idx.parentptr);
      }
    }
  } else if (vnp->choice == OBJ_BIOSEQ) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
  } else if (vnp->choice == OBJ_SEQENTRY) {
    sep = vnp->data.ptrvalue;
  }
  return sep;
}

NLM_EXTERN void AddNewUniqueDescriptors (SeqDescrPtr PNTR new_set, SeqDescrPtr parent_set)
{
  SeqDescrPtr sdp, sdp_next, sdp_tmp, sdp_tmp_next;
  Boolean     found_match;

  if (new_set == NULL || parent_set == NULL) return;

  if (*new_set == NULL) {
    ValNodeLink (new_set,
                 AsnIoMemCopy ((Pointer) parent_set,
                               (AsnReadFunc) SeqDescrAsnRead,
                               (AsnWriteFunc) SeqDescrAsnWrite));
  } else {
    sdp = parent_set;
    while (sdp != NULL) {
      sdp_next = sdp->next;
      sdp->next = NULL;
      for (sdp_tmp = *new_set, found_match = FALSE;
           sdp_tmp != NULL && !found_match; 
           sdp_tmp = sdp_tmp->next) {
        sdp_tmp_next = sdp_tmp->next;
        sdp_tmp->next = NULL;
        if (AsnIoMemComp (sdp, sdp_tmp, (AsnWriteFunc) SeqDescrAsnWrite)) {
          found_match = TRUE;
        }
        sdp_tmp->next = sdp_tmp_next;
      }
      if (!found_match) {
        ValNodeLink (new_set,
                     AsnIoMemCopy ((Pointer) sdp,
                                   (AsnReadFunc) SeqDescrAsnRead,
                                   (AsnWriteFunc) SeqDescrAsnWrite));
      }
      sdp->next = sdp_next;
      sdp = sdp->next;
    }
  }
}

static void AddNewUniqueDescriptorsToSeqEntry (SeqEntryPtr sep, SeqDescrPtr parent_set)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;

  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL) {
      AddNewUniqueDescriptors (&(bsp->descr), parent_set);
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      AddNewUniqueDescriptors (&(bssp->descr), parent_set);
    }
  }
}

NLM_EXTERN void AddNewUniqueAnnotations (SeqAnnotPtr PNTR new_set, SeqAnnotPtr parent_set)
{
  SeqAnnotPtr sap, sap_next, sap_tmp, sap_tmp_next, sap_copy, last_sap;
  Boolean     found_match;

  if (new_set == NULL || parent_set == NULL) return;

  sap = parent_set;
  while (sap != NULL) {
    sap_next = sap->next;
    sap->next = NULL;
    last_sap = NULL;
    for (sap_tmp = *new_set, found_match = FALSE;
         sap_tmp != NULL && !found_match;
         sap_tmp = sap_tmp->next) {
      sap_tmp_next = sap_tmp->next;
      sap_tmp->next = NULL;
      if (AsnIoMemComp (sap, sap_tmp, (AsnWriteFunc) SeqAnnotAsnWrite)) {
        found_match = TRUE;
      }
      sap_tmp->next = sap_tmp->next;
      last_sap = sap_tmp;
    }
    if (!found_match) {        
      sap_copy = (SeqAnnotPtr) AsnIoMemCopy ((Pointer) sap, (AsnReadFunc) SeqAnnotAsnRead, (AsnWriteFunc) SeqAnnotAsnWrite); 
      if (last_sap == NULL) {
        *new_set = sap_copy;
      } else {
        last_sap->next = sap_copy;
      }
    }
    sap->next = sap_next;
    sap = sap->next;
  }
}


static void AddItemListToSet (ValNodePtr item_list, BioseqSetPtr newset, Boolean for_segregate)
{
  ValNodePtr vnp_item;
  SeqEntryPtr sep, last_sep, prev_sep, remove_sep;
  BioseqSetPtr bssp, orig_parent;
  BioseqPtr bsp;

  if (newset == NULL || item_list == NULL) return;

  last_sep = newset->seq_set;
  while (last_sep != NULL && last_sep->next != NULL) {
    last_sep = last_sep->next;
  }

  for (vnp_item = item_list; vnp_item != NULL; vnp_item = vnp_item->next) {
    sep = GetBestSeqEntryForItem (vnp_item);
    if (sep == NULL || sep->data.ptrvalue == NULL) continue;
    orig_parent = NULL;
    bsp = NULL;
    bssp = NULL;
    if (IS_Bioseq (sep)) {
      bsp = sep->data.ptrvalue;
      if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
        orig_parent = bsp->idx.parentptr;
        bsp->idx.parentptr = NULL;
      }
    } else if (IS_Bioseq_set (sep)) {
      bssp = sep->data.ptrvalue;
      if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
        orig_parent = bssp->idx.parentptr;
        bssp->idx.parentptr = NULL;
      }
    } else {
      continue;
    }
       
    if (orig_parent != NULL) {
      /* remove this seq-entry from the original parent */
      prev_sep = NULL;
      for (remove_sep = orig_parent->seq_set;
           remove_sep != NULL && remove_sep != sep;
           remove_sep = remove_sep->next) {
        prev_sep = remove_sep;
      }
      if (remove_sep == sep) {
        if (prev_sep == NULL) {
          orig_parent->seq_set = orig_parent->seq_set->next;
          if (orig_parent->seq_set == NULL) {
            orig_parent->idx.deleteme = TRUE; 
          }
        } else {
          prev_sep->next = sep->next;
        }
      }
      /* set class type if not already set */
      if (newset->_class == BioseqseqSet_class_genbank && for_segregate) {
        newset->_class = orig_parent->_class;
      }
    }
    if (orig_parent != NULL) {
      if (for_segregate) {
        /* add descriptors from the orig_parent to the new parent */
        AddNewUniqueDescriptors (&(newset->descr), orig_parent->descr);

        /* add annotations from the orig_parent to the new parent */
        AddNewUniqueAnnotations (&(newset->annot), orig_parent->annot);
      } else {
        /* add descriptors from the orig_parent to the bioseq itself (or nuc-prot-set if that's what moved) */
        if (bsp != NULL) {
          AddNewUniqueDescriptors (&(bsp->descr), orig_parent->descr);
        } else if (bssp != NULL) {
          AddNewUniqueDescriptors (&(bssp->descr), orig_parent->descr);
        }
      }
    }

    /* add to new parent */
    sep->next = NULL;
    if (last_sep == NULL) {
      newset->seq_set = sep;
    } else {
      last_sep->next = sep;
    }
    last_sep = sep;
    SeqMgrLinkSeqEntry (sep, OBJ_BIOSEQSET, newset);
  }

}


static void AddCategorySeqEntriesToSet (BioseqSetPtr newset, ClickableItemPtr category)
{
  ValNodePtr vnp_item;

  if (newset == NULL || category == NULL || category->item_list == NULL) return;

  if (category->chosen) {
    AddItemListToSet (category->item_list, newset, TRUE);
  } else {
    for (vnp_item = category->subcategories; vnp_item != NULL; vnp_item = vnp_item->next) {
      AddCategorySeqEntriesToSet (newset, vnp_item->data.ptrvalue);
    }
  }
}

static Boolean NeedsNewSet (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  while (sep != NULL) {
    if (IS_Bioseq (sep)) {
      return TRUE;
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL 
          && (bssp->_class == BioseqseqSet_class_nuc_prot
          || bssp->_class == BioseqseqSet_class_segset)) {
        return TRUE;
      }
    }
    sep = sep->next;
  }
  return FALSE;
}


static Boolean IsSingletonSet (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  SeqAnnotPtr sap;

  if (sep == NULL 
      || !IS_Bioseq_set(sep) 
      || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL
      || bssp->seq_set == NULL
      || bssp->seq_set->next != NULL) {
    return FALSE;
  }

  /* not a singleton set if it has an alignment annotation */
  for (sap = bssp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 2) {
      return FALSE;
    }
  }
  return TRUE;
}


static void AddAnnotsToSeqEntry (SeqEntryPtr sep, SeqAnnotPtr sap)
{
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  SeqAnnotPtr last_sap;

  if (sep == NULL) {
    return;
  }
  if (IS_Bioseq(sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp->annot == NULL) {
      bsp->annot = sap;
    } else {
      last_sap = bsp->annot;
      while (last_sap->next != NULL) {
        last_sap = last_sap->next;
      }
      last_sap->next = sap;
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->annot == NULL) {
      bssp->annot = sap;
    } else {
      last_sap = bssp->annot;
      while (last_sap->next != NULL) {
        last_sap = last_sap->next;
      }
      last_sap->next = sap;
    }
  }
}


static void PromoteSingletonSetsInSet (SeqEntryPtr sep)
{
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  BioseqSetPtr      bssp, child_bssp;
  SeqEntryPtr       sep_child, child_next, child_prev = NULL;
  ValNodePtr        titles;
  Uint2             top_parenttype;
  Pointer           top_parentptr;

  if (sep == NULL || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL) {
    return;
  }

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &top_parentptr, &top_parenttype);

  for (sep_child = bssp->seq_set; sep_child != NULL; sep_child = child_next) {
    child_next = sep_child->next;
    if (IsSingletonSet(sep_child)) {
      child_bssp = (BioseqSetPtr) sep_child->data.ptrvalue;
      /* remove set title if any */
      titles = ValNodeExtractList (&(child_bssp->descr), Seq_descr_title);
      titles = SeqDescrFree (titles);
      /* propagate remaining descriptors */
      SetDescriptorPropagate (child_bssp);
      /* push down annotation */
      AddAnnotsToSeqEntry (child_bssp->seq_set, child_bssp->annot);
      /* replace in list */
      if (child_prev == NULL) {
        bssp->seq_set = child_bssp->seq_set;
      } else {
        child_prev->next = child_bssp->seq_set;
      }
      child_bssp->seq_set->next = child_next;
      child_prev = child_bssp->seq_set;
      child_bssp->seq_set = NULL;
      sep_child = SeqEntryFree (sep_child);
    } else {
      child_prev = sep_child;
    }
  }
      
  SeqMgrLinkSeqEntry (sep, top_parenttype, top_parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
}


NLM_EXTERN BioseqSetPtr MakeGroupsForUniqueValues
(BioseqSetPtr bssp, 
 ValNodePtr   value_lists)
{
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  BioseqSetPtr parent_set;
  SeqEntryPtr  sep, first_new_sep = NULL;
  SeqEntryPtr  tmp;
  BioseqSetPtr newset;
  ValNodePtr   vnp;
  ClickableItemPtr cip;
  Uint1        child_class;
  Uint2        entityID;
  Boolean      child_became_parent = FALSE;

  if (bssp == NULL) return NULL;
  entityID = bssp->idx.entityID;
  
  sep = SeqMgrGetSeqEntryForData (bssp);

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  parent_set = (BioseqSetPtr)(bssp->idx.parentptr);

  if (parent_set == NULL || parent_set->seq_set == NULL) {
    /* this set has no parent, so make it the parent set, class GenBank,
     * and create two new sets using the original set class as members of this set
     */
    parent_set = bssp;
    child_class = parent_set->_class;
    child_became_parent = TRUE;
  } else {
    /* we already have a parent set. */
    child_class = bssp->_class;
  }

  for (vnp = value_lists; vnp != NULL; vnp = vnp->next) {      
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip == NULL || (!cip->chosen && ! AnyDiscrepanciesChosen (cip->subcategories))) {
      continue;
    }

    newset = BioseqSetNew ();
    newset->_class = child_class;
    /* add SeqEntries for this category here */
    AddCategorySeqEntriesToSet (newset, cip);

    /* check - is any member of the newset also a set? */
    tmp = newset->seq_set;
    while (tmp != NULL && IS_Bioseq (tmp)) {
      tmp = tmp->next;
    }

    /* add to members for parent_set */
    tmp = SeqEntryNew ();
    tmp->choice = 2;
    tmp->data.ptrvalue = (Pointer) newset;

    tmp->next = parent_set->seq_set;
    parent_set->seq_set = tmp;

    if (parent_set->idx.deleteme) {
      parent_set->idx.deleteme = FALSE;
    }
    if (first_new_sep == NULL) {
      first_new_sep = tmp;
    }

    SeqMgrLinkSeqEntry (tmp, OBJ_BIOSEQSET, parent_set);
  }

  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 
  DeleteMarkedObjects (entityID, 0, NULL);

  if (child_became_parent) {
    if (first_new_sep != NULL && NeedsNewSet (first_new_sep->next)) {
      /* make new set to hold what wasn't moved */
      newset = BioseqSetNew ();
      /* contents of set are unmoved sequences */
      newset->seq_set = first_new_sep->next;
      /* set class type */
      newset->_class = parent_set->_class;
      /* add descriptors from the parent */
      ValNodeLink (&(newset->descr),
                   AsnIoMemCopy ((Pointer) parent_set->descr,
                                 (AsnReadFunc) SeqDescrAsnRead,
                                 (AsnWriteFunc) SeqDescrAsnWrite));


      tmp = SeqEntryNew();
      tmp->choice = 2;
      tmp->data.ptrvalue = (Pointer) newset;
      tmp->next = NULL;
      first_new_sep->next = tmp;
      SeqMgrLinkSeqEntry (tmp, OBJ_BIOSEQSET, parent_set);
    } else {
      sep = first_new_sep->next;
      while (sep != NULL) {
        AddNewUniqueDescriptorsToSeqEntry (sep, parent_set->descr);
        sep = sep->next;
      }
    }
    /* set parent class to GenBank set */
    parent_set->_class = BioseqseqSet_class_genbank;
    /* remove descriptors on parent, they will all have been copied down */
    parent_set->descr = SeqDescrFree (parent_set->descr);
  }

  sep = SeqMgrGetSeqEntryForData (parent_set);
  PromoteSingletonSetsInSet (sep);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);

  return parent_set;
}


static void RemoveBioseqFromAlignmentsCallback (SeqAnnotPtr sap, Pointer data)
{
  BioseqPtr   bsp;
  SeqAlignPtr salphead, salp, salp_next, prev_salp, remove_salp, last_remove;
  SeqIdPtr    sip, tmpsip;
  Uint4       seqid_order;

  if (sap == NULL || sap->type != 2 
      || (bsp = (BioseqPtr) data) == NULL
      || (salphead = (SeqAlignPtr) sap->data) == NULL) {
    return;
  }
  
  salp = salphead;
  prev_salp = NULL;
  remove_salp = NULL;
  last_remove = NULL;
  while (salp != NULL)
  {
    salp_next = salp->next;
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = 0;
    for (sip = bsp->id; sip != NULL && seqid_order == 0; sip = sip->next) {
      seqid_order = SeqIdOrderInBioseqIdList(sip, tmpsip);
    }
    if (seqid_order == 0)
    {
      /* do nothing for this subalignment */
      prev_salp = salp;
    }
    else if (salp->dim == 2 || salphead->segtype ==1)
    {
      /* This is for a pairwise alignment or a DENDIAG alignment */
      if (prev_salp == NULL)
      {
          salphead = salp->next;
      }
      else
      {
          prev_salp->next = salp->next;
      }
      /* save the alignments that we want to free in a list and get rid of them
       * at the end - freeing them beforehand causes problems with listing the
       * IDs in the alignment.
       */
      salp->next = NULL;
      if (remove_salp == NULL)
      {
          remove_salp = salp;
      }
      else
      {
          last_remove->next = salp;
      }
      last_remove = salp;
    }
    else 
    {
      SeqAlignBioseqDeleteById (salphead, sip);  
      prev_salp = salp;
    }
    salp = salp_next;
  }
  /* Now we can free the alignment */
  SeqAlignFree (remove_salp);
  
  sap->data = salphead;
  if (sap->data == NULL) {
    sap->idx.deleteme = TRUE;
  }
}


/* expect that list is a valnode list with choice OBJ_BIOSEQ and data.ptrvalue a bioseq */
NLM_EXTERN void MoveSequencesFromSetToWrapper (ValNodePtr list, Uint2 entityID)
{
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  BioseqSetPtr      bssp;
  SeqEntryPtr       sep;
  ValNodePtr        vnp;

  if (list == NULL) return;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL 
      || !IS_Bioseq_set (sep) 
      || (bssp = sep->data.ptrvalue) == NULL
      || bssp->_class != BioseqseqSet_class_genbank) {
    return;
  }

  /* first, propagate descriptors */
  SetDescriptorPropagate (bssp);

  /* pull sequences out of current positions and add to top-level set */
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  AddItemListToSet (list, bssp, FALSE);

  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata); 

  /* remove sequences from alignments */
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    VisitAnnotsInSep (sep, vnp->data.ptrvalue, RemoveBioseqFromAlignmentsCallback);
  }
  DeleteMarkedObjects (entityID, 0, NULL);

  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);

}


static void GetBioseqListCallback (BioseqPtr bsp, Pointer userdata)
{
  if (bsp != NULL && userdata != NULL && ! ISA_aa (bsp->mol))
  {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_BIOSEQ, bsp);
  }
}

NLM_EXTERN ValNodePtr PrepareSequenceListForSegregateByNumberOfSets (Int4 num_sets, SeqEntryPtr sep)
{
  ValNodePtr cip_list = NULL;
  ValNodePtr         vnp;
  CharPtr            fmt = "set contains %d sequences";
  ValNodePtr         bsp_list = NULL, v_next;
  Int4               num_bioseqs, num_per_set, num_in_set;

  VisitBioseqsInSep (sep, &bsp_list, GetBioseqListCallback);
  num_bioseqs = ValNodeLen (bsp_list);
  num_per_set = num_bioseqs / num_sets;
  if (num_bioseqs % num_sets > 0) {
    num_per_set ++;
  }

  num_in_set = 0;
  
  for (vnp = bsp_list; vnp != NULL; vnp = v_next) {
    v_next = vnp->next;
    num_in_set ++;
    if (num_in_set == num_per_set) {
      vnp->next = NULL;
      ValNodeAddPointer (&cip_list, 0, NewClickableItem (0, fmt, bsp_list));
      bsp_list = v_next;
      num_in_set = 0;
    }
  }
  if (bsp_list != NULL) {
    ValNodeAddPointer (&cip_list, 0, NewClickableItem (0, fmt, bsp_list));
    bsp_list = NULL;
  }
  return cip_list;
}


NLM_EXTERN ValNodePtr PrepareSequenceListForSegregateByNumberPerSet (Int4 num_per_set, SeqEntryPtr sep)
{
  ValNodePtr cip_list = NULL;
  ValNodePtr         vnp;
  CharPtr            fmt = "set contains %d sequences";
  ValNodePtr         bsp_list = NULL, v_next;
  Int4               num_in_set;

  VisitBioseqsInSep (sep, &bsp_list, GetBioseqListCallback);

  num_in_set = 0;
  
  for (vnp = bsp_list; vnp != NULL; vnp = v_next) {
    v_next = vnp->next;
    num_in_set ++;
    if (num_in_set == num_per_set) {
      vnp->next = NULL;
      ValNodeAddPointer (&cip_list, 0, NewClickableItem (0, fmt, bsp_list));
      bsp_list = v_next;
      num_in_set = 0;
    }
  }
  if (bsp_list != NULL) {
    ValNodeAddPointer (&cip_list, 0, NewClickableItem (0, fmt, bsp_list));
    bsp_list = NULL;
  }
  return cip_list;
}


NLM_EXTERN void SegregateSetsByNumber (SeqEntryPtr sep, Int4 num_sets)
{
  ValNodePtr set_list;
  BioseqSetPtr bssp;

  if (sep == NULL || !IS_Bioseq_set (sep) || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL) {
    return;
  }

  set_list = PrepareSequenceListForSegregateByNumberOfSets (num_sets, sep);
  ChooseAllDiscrepancies(set_list);

  MakeGroupsForUniqueValues (bssp, set_list);

  FreeClickableList (set_list);

}


NLM_EXTERN void SegregateSetsByNumberPerSet (SeqEntryPtr sep, Int4 num_per_set)
{
  ValNodePtr set_list;
  BioseqSetPtr bssp;

  if (sep == NULL || !IS_Bioseq_set (sep) || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL) {
    return;
  }

  set_list = PrepareSequenceListForSegregateByNumberPerSet (num_per_set, sep);
  ChooseAllDiscrepancies(set_list);

  MakeGroupsForUniqueValues (bssp, set_list);

  FreeClickableList (set_list);

}


static CharPtr s_PlantGroupList[] = {
  "Anthocerotophyta",
  "Bryophyta",
  "Charophyceae",
  "Chlorokybophyceae",
  "Chlorophyceae",
  "Coleochaetophyceae",
  "Coniferophyta",
  "Cycadophyta",
  "Ginkgophyta",
  "Gnetophyta",
  "Klebsormidiophyceae",
  "Lycopodiophyta",
  "Magnoliophyta",
  "Mamiellophyceae",
  "Marchantiophyta",
  "Mesostigmatophyceae",
  "Moniliformopses",
  "Pedinophyceae",
  "Picocystophyceae",
  "Prasinophyceae",
  "Trebouxiophyceae",
  "Ulvophyceae",
  "Zygnemophyceae",
  NULL
};


static void AssignBioseqToPlantGroup (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;
  SeqMgrDescContext context;
  BioSourcePtr biop;
  CharPtr      lineage = NULL;
  ValNodePtr   list, vnp, vnp_last = NULL;
  ClickableItemPtr cip;

  if (bsp == NULL || ISA_aa (bsp->mol) || (list = (ValNodePtr) data) == NULL) {
    return;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL 
      && (biop = (BioSourcePtr) sdp->data.ptrvalue) != NULL
      && biop->org != NULL
      && biop->org->orgname != NULL) {
    lineage = biop->org->orgname->lineage;
  }

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    vnp_last = vnp;
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (vnp->next == NULL
        || StringISearch (lineage, cip->description) != NULL) {
      ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
      break;
    }
  }
}


static ValNodePtr MakePlantGroupList (SeqEntryPtr sep)
{
  ClickableItemPtr cip;
  Int4 i;
  ValNodePtr list = NULL, vnp, vnp_prev = NULL, vnp_next;

  for (i = 0; s_PlantGroupList[i] != NULL; i++) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    cip->description = StringSave(s_PlantGroupList[i]);
    ValNodeAddPointer (&list, 0, cip);
  }
  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->description = StringSave("Other");
  ValNodeAddPointer (&list, 0, cip);
  

  VisitBioseqsInSep (sep, list, AssignBioseqToPlantGroup);

  for (vnp = list; vnp != NULL; vnp = vnp_next) {
    vnp_next = vnp->next;
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip->item_list == NULL) {
      if (vnp_prev == NULL) {
        list = vnp_next;
      } else {
        vnp_prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp->data.ptrvalue = ClickableItemFree (vnp->data.ptrvalue);
      vnp = ValNodeFree (vnp);
    } else {
      vnp_prev = vnp;
    }
  }
  return list;
}


NLM_EXTERN void SegregateSetsByPlantGroup (SeqEntryPtr sep)
{
  ValNodePtr set_list;
  BioseqSetPtr bssp;

  if (sep == NULL || !IS_Bioseq_set (sep) || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL) {
    return;
  }

  set_list = MakePlantGroupList(sep);
  ChooseAllDiscrepancies(set_list);

  MakeGroupsForUniqueValues (bssp, set_list);

  FreeClickableList (set_list);
}


typedef  Boolean  (*Nlm_ParseProc) PROTO ((CharPtr, Pointer));


static Boolean SkipToken (CharPtr cp, Pointer data)
{
  return TRUE;
}


static Boolean ParseLineOfTokens (CharPtr line, Nlm_ParseProc PNTR token_funcs, Pointer data)
{
  CharPtr cp, cp_next;
  Char    ch_was;
  Int4    token_num = 0;
  Boolean rval = TRUE;

  if (StringHasNoText (line) || token_funcs == NULL) {
    return FALSE;
  }

  cp = line;
  cp_next = StringChr (cp, '\t');
  while (cp_next != NULL && rval && token_funcs[token_num] != NULL) {
    ch_was = *cp_next;
    *cp_next = 0;
    rval = token_funcs[token_num] (cp, data);
    *cp_next = ch_was;
    cp = cp_next + 1;
    cp_next = StringChr (cp, '\t');
    token_num++;
  }

  if (rval && token_funcs[token_num] != NULL) {
    /* last token_func for end of line */
    rval = token_funcs[token_num](cp, data);
    token_num++;
    while (token_funcs[token_num] != NULL && rval) {
      rval = token_funcs[token_num](NULL, data);
      token_num++;
    }
  }
  return rval;
}


/* Output from Fungal ITS sequence extractor:
 * first column is ID (after position in set) and length
 * next column is ITS1 or ---- (---- means it's not there)
 * next column is ITS2 or ---- (---- means it's not there)
 * next column is range for ITS1 (or -----)
 * next column is range for ITS2 (or -----)
 * next column (if present) indicates reverse complement.
 */

typedef struct extractorinfo {
  CharPtr id;
  Int4    length;
  Boolean has_its1;
  Boolean has_its2;
  CharPtr its1_range;
  CharPtr its2_range;
  Boolean is_complement;
} ExtractorInfoData, PNTR ExtractorInfoPtr;


static ExtractorInfoPtr ExtractorInfoNew ()
{
  ExtractorInfoPtr ep = (ExtractorInfoPtr) MemNew (sizeof (ExtractorInfoData));
  MemSet (ep, 0, sizeof (ExtractorInfoData));
  return ep;
}


static ExtractorInfoPtr ExtractorInfoFree (ExtractorInfoPtr ep)
{
  if (ep != NULL) {
    ep->id = MemFree (ep->id);
    ep->its1_range = MemFree (ep->its1_range);
    ep->its2_range = MemFree (ep->its2_range);
    ep = MemFree (ep);
  }
  return ep;
}


static Boolean ParseExtractorIdAndLength (CharPtr cp, Pointer data)
{
  ExtractorInfoPtr ep;
  CharPtr div, id_start, id_end;
  Int4    len;
  Char    ch_was;

  if (StringHasNoText (cp) || (ep = (ExtractorInfoPtr) data) == NULL) {
    return FALSE;
  }
  /* separate length and ID */
  len = StringLen (cp);
  if (len < 4) {
    return FALSE;
  }
  if (StringCmp (cp + len - 4, " bp.") != 0) {
    return FALSE;
  }

  div = cp + len - 5;
  while (div > cp && isdigit (*div)) {
    div--;
  }
  if (!isdigit (*(div + 1))) {
    return FALSE;
  }
  ep->length = atoi (div + 1);

  /* skip over the part that indicates the position of the sequence (1 of N, 2 of N, etc.) */
  id_start = StringChr (cp, ')');
  if (id_start == NULL) {
    return FALSE;
  }
  id_start++;
  while (isspace (*id_start)) {
    id_start++;
  }

  if (id_start >= div) {
    return FALSE;
  }

  /* if we have a list of IDs, truncate after just the first one */
  id_end = StringChr (id_start, '|');
  if (id_end != NULL && id_end < div) {
    id_end = StringChr (id_end + 1, '|');
    if (id_end != NULL && id_end < div) {
      div = id_end;
    }
  }
  
  ch_was = *div;
  *div = 0;
  ep->id = StringSave (id_start);  
  *div = ch_was;
  /* trim spaces from end of ID */
  cp = ep->id + StringLen (ep->id) - 1;
  while (cp > ep->id && isspace (*cp)) {
    cp--;
  }
  *(cp + 1) = 0;
  return TRUE;
}


static Boolean ParseHasITS1 (CharPtr cp, Pointer data)
{
  ExtractorInfoPtr ep;

  if (StringHasNoText (cp) || (ep = (ExtractorInfoPtr) data) == NULL) {
    return FALSE;
  }

  if (StringCmp (cp, "ITS1") == 0) {
    ep->has_its1 = TRUE;
  } else if (StringCmp (cp, "----") == 0) {
    ep->has_its1 = FALSE;
  } else {
    return FALSE;
  }

  return TRUE;
}


static Boolean ParseHasITS2 (CharPtr cp, Pointer data)
{
  ExtractorInfoPtr ep;

  if (StringHasNoText (cp) || (ep = (ExtractorInfoPtr) data) == NULL) {
    return FALSE;
  }

  if (StringCmp (cp, "ITS2") == 0) {
    ep->has_its2 = TRUE;
  } else if (StringCmp (cp, "----") == 0) {
    ep->has_its2 = FALSE;
  } else {
    return FALSE;
  }

  return TRUE;
}


static Boolean ParseITS1Range (CharPtr cp, Pointer data)
{
  ExtractorInfoPtr ep;
  Boolean          rval = TRUE;

  if (StringHasNoText (cp) || (ep = (ExtractorInfoPtr) data) == NULL) {
    return FALSE;
  }

  if (*cp == '-') {
    if (ep->has_its1) {
      rval = FALSE;
    } else {
      rval = TRUE;
    }
  } else if (StringNCmp (cp, "ITS1: ", 6) == 0) {
    if (ep->has_its1) {
      ep->its1_range = StringSave (cp + 6);
    } else {
      rval = FALSE;
    }
  } else {
    rval = FALSE;
  }
  return rval;
}


static Boolean ParseITS2Range (CharPtr cp, Pointer data)
{
  ExtractorInfoPtr ep;
  Boolean          rval = TRUE;

  if (StringHasNoText (cp) || (ep = (ExtractorInfoPtr) data) == NULL) {
    return FALSE;
  }

  if (*cp == '-') {
    if (ep->has_its2) {
      rval = FALSE;
    } else {
      rval = TRUE;
    }
  } else if (StringNCmp (cp, "ITS2: ", 6) == 0) {
    if (ep->has_its2) {
      ep->its2_range = StringSave (cp + 6);
    } else {
      rval = FALSE;
    }
  } else {
    rval = FALSE;
  }
  return rval;
}


static Boolean ParseIsComplement (CharPtr cp, Pointer data)
{
  ExtractorInfoPtr ep;
  Boolean rval = TRUE;

  if ((ep = (ExtractorInfoPtr) data) == NULL) {
    return FALSE;
  }
  if (StringHasNoText (cp)) {
    ep->is_complement = FALSE;
  } else if (StringNCmp (cp, "Reverse complementary", 21) == 0) {
    ep->is_complement = TRUE;
  } else {
    rval = FALSE;
  }
  return rval;
}


static Nlm_ParseProc token_parsers[] = {
  ParseExtractorIdAndLength,
  ParseHasITS1,
  ParseHasITS2,
  SkipToken,
  SkipToken,
  ParseITS1Range,
  ParseITS2Range,
  ParseIsComplement,
  NULL};


typedef enum {
  eExtractorFeat18S = 0,
  eExtractorFeatITS1,
  eExtractorFeat58S,
  eExtractorFeatITS2,
  eExtractorFeat28S
} EExtractorFeat;

CharPtr extractor_feature_labels[] = {
  "18S ribosomal RNA",
  "internal transcribed spacer 1",
  "5.8S ribosomal RNA",
  "internal transcribed spacer 2",
  "28S ribosomal RNA"
};


static CharPtr MakeLabelFromExtractorInfo (ExtractorInfoPtr ep)
{
  Boolean feat_present[5];
  CharPtr cp, label;
  Int4 len, i, num_feat = 0, feat_num = 0;

  if (ep == NULL) {
    return NULL;
  }

  MemSet (feat_present, 0, sizeof (feat_present));
  if (ep->has_its1) {
    feat_present[eExtractorFeatITS1] = TRUE;
    if (StringNCmp (ep->its1_range, "1-", 2) == 0) {
      feat_present[eExtractorFeat18S] = FALSE;
    } else {
      feat_present[eExtractorFeat18S] = TRUE;
    }
    if (ep->has_its2) {
      feat_present[eExtractorFeat58S] = TRUE;
      feat_present[eExtractorFeatITS2] = TRUE;
      cp = StringChr (ep->its2_range, '-');
      if (cp != NULL && StringCmp (cp + 1, "end") == 0) {
        feat_present[eExtractorFeat28S] = FALSE;
      } else {
        feat_present[eExtractorFeat28S] = TRUE;
      }
    } else {
      cp = StringChr (ep->its1_range, '-');
      if (cp != NULL && StringCmp (cp + 1, "end") == 0) {
        feat_present[eExtractorFeat58S] = FALSE;
      } else {
        feat_present[eExtractorFeat58S] = TRUE;
      }
      feat_present[eExtractorFeatITS2] = FALSE;
      feat_present[eExtractorFeat28S] = FALSE;
    }
  } else {
    feat_present[eExtractorFeat18S] = FALSE;
    feat_present[eExtractorFeatITS1] = FALSE;
    if (StringNCmp (ep->its2_range, "1-", 2) == 0) {
      feat_present[eExtractorFeat58S] = FALSE;
    } else {
      feat_present[eExtractorFeat58S] = TRUE;
    }
    feat_present[eExtractorFeatITS2] = TRUE;
    cp = StringChr (ep->its2_range, '-');
    if (cp != NULL && StringCmp (cp + 1, "end") == 0) {
      feat_present[eExtractorFeat28S] = FALSE;
    } else {
      feat_present[eExtractorFeat28S] = TRUE;
    }
  }

  len = 15;
  for (i = 0; i < 5; i++) {
    if (feat_present[i]) {
      len += StringLen (extractor_feature_labels[i]) + 2;
      num_feat++;
    } else if (num_feat > 0) {
      break;
    }
  }
  label = (CharPtr) MemNew (sizeof (Char) * len);
  sprintf (label, "contains ");
  for (i = 0; i < 5; i++) {
    if (feat_present[i]) {
      if (feat_num > 0) {
        if (feat_num == num_feat - 1) {
          if (num_feat == 2) {
            StringCat (label, " and ");
          } else {
            StringCat (label, ", and ");
          }
        } else {
          StringCat (label, ", ");
        }
      }
      StringCat (label, extractor_feature_labels[i]);
      feat_num++;
    } else if (feat_num > 0) {
      break;
    }
  }
  return label;
}


NLM_EXTERN void RevCompFeats (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  bsp = (BioseqPtr) mydata;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        RevCompOneFeatForBioseq (sfp, bsp);
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}


static SeqFeatPtr ParseExtractorResultRowToFeatures (CharPtr line, SeqEntryPtr sep)
{
  ExtractorInfoPtr ep;
  SeqFeatPtr sfp = NULL;
  CharPtr  label;
  Int4     len;
  SeqIdPtr sip;
  BioseqPtr bsp;
  RnaRefPtr rrp;
  RNAGenPtr rgp;

  if (StringHasNoText (line)) {
    return NULL;
  }

  ep = ExtractorInfoNew ();
  if (!ParseLineOfTokens(line, token_parsers, ep)) {
    ep = ExtractorInfoFree (ep);
    Message (MSG_POSTERR, "Unable to parse extractor line %s", line);
    return NULL;
  }
  if (!ep->has_its1 && !ep->has_its2) {
    ep = ExtractorInfoFree (ep);
    Message (MSG_POSTERR, "Unable to determine feature list for line %s", line);
    return NULL;
  }

  /* figure out ID */
  len = StringLen (ep->id);
  if (len > 3 && ep->id[len - 1] == '.' && ep->id[len - 2] == '.' && ep->id[len - 3] == '.') {
    ep = ExtractorInfoFree (ep);
    Message (MSG_POSTERR, "ID was truncated for line %s", line);
    return NULL;
  }
  sip = CreateSeqIdFromText (ep->id, sep);
  bsp = BioseqFind (sip);
  sip = SeqIdFree (sip);
  if (bsp == NULL) {
    ep = ExtractorInfoFree (ep);
    Message (MSG_POSTERR, "ID for sequence not present in record in line %s", line);
    return NULL;
  }


  /* calculate label */
  label = MakeLabelFromExtractorInfo(ep);

  if (ep->is_complement) {
    BioseqRevComp (bsp);
    SeqEntryExplore (sep, (Pointer) bsp, RevCompFeats);
  }

  /* make feature and attach to appropriate annots */
  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_RNA, NULL);
  rrp = RnaRefNew ();
  sfp->data.value.ptrvalue = rrp;
  rrp->type = 255;
  rgp = RNAGenNew ();
  rrp->ext.choice = 3;
  rrp->ext.value.ptrvalue = rgp;
  sfp->comment = label;
  SetSeqLocPartial (sfp->location, TRUE, TRUE);
  
  ep = ExtractorInfoFree (ep);
  return sfp;
}

NLM_EXTERN void ParseExtractorResultsTableToFeatures (FILE *fp, SeqEntryPtr sep)
{
  ReadBufferData  rbd;
  CharPtr         line;

  rbd.fp = fp;
  rbd.current_data = NULL;
  line = AbstractReadFunction (&rbd); 
  while (line != NULL && line[0] != EOF) {
    /* TODO: skip intro lines */
    ParseExtractorResultRowToFeatures(line, sep);
    line = MemFree (line);
    line = AbstractReadFunction (&rbd);
  }
}

extern CharPtr latlon_onedegree [];
CharPtr latlon_onedegree [] = {
  "1",
  "Afghanistan",
  "\t39\t69\t72",
  "\t38\t63\t75",
  "\t37\t62\t75",
  "\t36\t60\t75",
  "\t35\t59\t75",
  "\t34\t59\t72",
  "\t33\t59\t72",
  "\t32\t59\t71",
  "\t31\t59\t70",
  "\t30\t59\t70",
  "\t29\t59\t67",
  "\t28\t59\t67",
  "Albania",
  "\t43\t18\t21",
  "\t42\t18\t21",
  "\t41\t18\t22",
  "\t40\t18\t22",
  "\t39\t18\t22",
  "\t38\t18\t21",
  "Algeria",
  "\t38\t5\t8",
  "\t37\t-1\t9",
  "\t36\t-3\t9",
  "\t35\t-3\t9",
  "\t34\t-3\t9",
  "\t33\t-3\t10",
  "\t32\t-4\t10",
  "\t31\t-6\t10",
  "\t30\t-9\t10",
  "\t29\t-9\t10",
  "\t28\t-9\t10",
  "\t27\t-9\t10",
  "\t26\t-9\t11",
  "\t25\t-9\t12",
  "\t24\t-7\t12",
  "\t23\t-5\t12",
  "\t22\t-4\t12",
  "\t21\t-2\t12",
  "\t20\t-1\t10",
  "\t19\t0\t8",
  "\t18\t1\t7",
  "\t17\t2\t4",
  "American Samoa",
  "\t-10\t-172\t-170",
  "\t-11\t-172\t-170",
  "\t-12\t-172\t-170",
  "\t-13\t-171\t-167",
  "\t-14\t-171\t-167",
  "\t-15\t-171\t-167",
  "Andorra",
  "\t43\t0\t2",
  "\t42\t0\t2",
  "\t41\t0\t2",
  "Angola",
  "\t-3\t11\t14",
  "\t-4\t11\t17",
  "\t-5\t11\t17\t19\t21",
  "\t-6\t11\t22",
  "\t-7\t11\t22",
  "\t-8\t11\t23",
  "\t-9\t11\t25",
  "\t-10\t11\t25",
  "\t-11\t11\t25",
  "\t-12\t11\t25",
  "\t-13\t11\t25",
  "\t-14\t10\t25",
  "\t-15\t10\t23",
  "\t-16\t10\t24",
  "\t-17\t10\t24",
  "\t-18\t10\t24",
  "\t-19\t19\t22",
  "Anguilla",
  "\t19\t-64\t-61",
  "\t18\t-64\t-61",
  "\t17\t-64\t-61",
  "Antarctica",
  "\t-59\t-47\t-43",
  "\t-60\t-59\t-53\t-47\t-43",
  "\t-61\t-62\t-53\t-47\t-43",
  "\t-62\t-62\t-53",
  "\t-63\t-65\t-54",
  "\t-64\t-66\t-54\t51\t56\t99\t104\t110\t114",
  "\t-65\t-69\t-56\t47\t58\t86\t117\t119\t144",
  "\t-66\t-70\t-59\t42\t70\t79\t147",
  "\t-67\t-91\t-89\t-73\t-59\t31\t35\t38\t71\t76\t156",
  "\t-68\t-91\t-89\t-76\t-60\t31\t161",
  "\t-69\t-91\t-89\t-77\t-60\t-11\t168",
  "\t-70\t-103\t-95\t-77\t-59\t-13\t171",
  "\t-71\t-106\t-87\t-81\t-79\t-77\t-58\t-15\t171",
  "\t-72\t-128\t-117\t-115\t-112\t-106\t-58\t-22\t-19\t-17\t171",
  "\t-73\t-137\t-109\t-106\t-58\t-23\t171",
  "\t-74\t-147\t-58\t-27\t170",
  "\t-75\t-150\t-59\t-32\t166",
  "\t-76\t-159\t-62\t-48\t-44\t-36\t170",
  "\t-77\t-165\t-65\t-51\t-42\t-37\t170",
  "\t-78\t-165\t-64\t-62\t-58\t-52\t-41\t-37\t170",
  "\t-79\t-165\t-58\t-55\t168",
  "\t-80\t-165\t-58\t-55\t164",
  "\t-81\t-175\t-170\t-164\t169",
  "\t-82\t-175\t177",
  "\t-83\t-180\t180",
  "\t-84\t-180\t180",
  "\t-85\t-180\t180",
  "\t-86\t-180\t180",
  "\t-87\t-180\t180",
  "\t-88\t-180\t180",
  "\t-89\t-180\t180",
  "\t-90\t-180\t180",
  "\t-90\t-180\t180",
  "Antigua and Barbuda",
  "\t18\t-62\t-60",
  "\t17\t-62\t-60",
  "\t16\t-62\t-60",
  "\t15\t-62\t-60",
  "Argentina",
  "\t-20\t-67\t-61",
  "\t-21\t-68\t-60",
  "\t-22\t-68\t-59",
  "\t-23\t-69\t-57",
  "\t-24\t-69\t-52",
  "\t-25\t-69\t-52",
  "\t-26\t-70\t-52",
  "\t-27\t-70\t-52",
  "\t-28\t-71\t-52",
  "\t-29\t-71\t-54",
  "\t-30\t-71\t-55",
  "\t-31\t-71\t-56",
  "\t-32\t-71\t-56",
  "\t-33\t-71\t-56",
  "\t-34\t-71\t-56",
  "\t-35\t-72\t-55",
  "\t-36\t-72\t-55",
  "\t-37\t-72\t-55",
  "\t-38\t-72\t-55",
  "\t-39\t-72\t-56",
  "\t-40\t-72\t-60",
  "\t-41\t-73\t-61",
  "\t-42\t-73\t-61",
  "\t-43\t-73\t-62",
  "\t-44\t-73\t-63",
  "\t-45\t-73\t-64",
  "\t-46\t-73\t-64",
  "\t-47\t-74\t-64",
  "\t-48\t-74\t-64",
  "\t-49\t-74\t-64",
  "\t-50\t-74\t-66",
  "\t-51\t-74\t-66",
  "\t-52\t-73\t-66",
  "\t-53\t-71\t-62",
  "\t-54\t-69\t-62",
  "\t-55\t-69\t-62",
  "\t-56\t-67\t-65",
  "Armenia",
  "\t42\t42\t46",
  "\t41\t42\t46",
  "\t40\t42\t47",
  "\t39\t42\t47",
  "\t38\t43\t47",
  "\t37\t45\t47",
  "Aruba",
  "\t13\t-71\t-68",
  "\t12\t-71\t-68",
  "\t11\t-71\t-68",
  "Ashmore and Cartier Islands",
  "\t-11\t122\t124",
  "\t-12\t122\t124",
  "\t-13\t122\t124",
  "Australia",
  "\t-8\t141\t143",
  "\t-9\t131\t133\t141\t143",
  "\t-10\t129\t137\t140\t144",
  "\t-11\t129\t137\t140\t144",
  "\t-12\t124\t137\t140\t144",
  "\t-13\t123\t137\t140\t146",
  "\t-14\t123\t138\t140\t146",
  "\t-15\t121\t146",
  "\t-16\t121\t147",
  "\t-17\t120\t147",
  "\t-18\t118\t149",
  "\t-19\t114\t150",
  "\t-20\t112\t151",
  "\t-21\t112\t151",
  "\t-22\t112\t152",
  "\t-23\t112\t154",
  "\t-24\t111\t154",
  "\t-25\t111\t154",
  "\t-26\t111\t154",
  "\t-27\t112\t154",
  "\t-28\t112\t154",
  "\t-29\t113\t154",
  "\t-30\t113\t154\t158\t160",
  "\t-31\t113\t154\t158\t160",
  "\t-32\t113\t154\t158\t160",
  "\t-33\t113\t129\t131\t153",
  "\t-34\t113\t125\t133\t152",
  "\t-35\t113\t124\t134\t152",
  "\t-36\t115\t119\t134\t151",
  "\t-37\t135\t151",
  "\t-38\t138\t151",
  "\t-39\t139\t149",
  "\t-40\t142\t149",
  "\t-41\t142\t149",
  "\t-42\t143\t149",
  "\t-43\t144\t149",
  "\t-44\t144\t149",
  "\t-53\t157\t159",
  "\t-54\t157\t159",
  "\t-55\t157\t159",
  "Australia: Australian Capital Territory",
  "\t-34\t147\t150",
  "\t-35\t147\t150",
  "\t-36\t147\t150",
  "Australia: Jervis Bay Territory",
  "\t-34\t149\t151",
  "\t-35\t149\t151",
  "\t-36\t149\t151",
  "Australia: New South Wales",
  "\t-27\t147\t154",
  "\t-28\t140\t154",
  "\t-29\t140\t154",
  "\t-30\t140\t154",
  "\t-31\t140\t154",
  "\t-32\t140\t154",
  "\t-33\t140\t153",
  "\t-34\t140\t152",
  "\t-35\t140\t152",
  "\t-36\t142\t151",
  "\t-37\t143\t151",
  "\t-38\t147\t151",
  "Australia: Northern Territory",
  "\t-9\t131\t133",
  "\t-10\t129\t137",
  "\t-11\t129\t137",
  "\t-12\t128\t137",
  "\t-13\t128\t137",
  "\t-14\t128\t138",
  "\t-15\t128\t139",
  "\t-16\t128\t139",
  "\t-17\t128\t139",
  "\t-18\t128\t139",
  "\t-19\t128\t139",
  "\t-20\t128\t139",
  "\t-21\t128\t139",
  "\t-22\t128\t139",
  "\t-23\t128\t139",
  "\t-24\t128\t139",
  "\t-25\t128\t139",
  "\t-26\t128\t139",
  "\t-27\t128\t139",
  "Australia: Queensland",
  "\t-9\t141\t143",
  "\t-10\t140\t144",
  "\t-11\t140\t144",
  "\t-12\t140\t144",
  "\t-13\t140\t146",
  "\t-14\t140\t146",
  "\t-15\t137\t146",
  "\t-16\t137\t147",
  "\t-17\t137\t147",
  "\t-18\t137\t149",
  "\t-19\t137\t150",
  "\t-20\t137\t151",
  "\t-21\t137\t151",
  "\t-22\t137\t152",
  "\t-23\t137\t154",
  "\t-24\t137\t154",
  "\t-25\t137\t154",
  "\t-26\t137\t154",
  "\t-27\t137\t154",
  "\t-28\t140\t154",
  "\t-29\t140\t154",
  "\t-30\t140\t152",
  "Australia: South Australia",
  "\t-25\t128\t142",
  "\t-26\t128\t142",
  "\t-27\t128\t142",
  "\t-28\t128\t142",
  "\t-29\t128\t142",
  "\t-30\t128\t142",
  "\t-31\t128\t142",
  "\t-32\t128\t142",
  "\t-33\t131\t142",
  "\t-34\t133\t142",
  "\t-35\t134\t142",
  "\t-36\t134\t141",
  "\t-37\t135\t141",
  "\t-38\t138\t141",
  "\t-39\t139\t141",
  "Australia: Tasmania",
  "\t-38\t142\t149",
  "\t-39\t142\t149",
  "\t-40\t142\t149",
  "\t-41\t142\t149",
  "\t-42\t143\t149",
  "\t-43\t144\t149",
  "\t-44\t144\t149",
  "Australia: Victoria",
  "\t-32\t139\t141",
  "\t-33\t139\t144",
  "\t-34\t139\t148",
  "\t-35\t139\t149",
  "\t-36\t139\t150",
  "\t-37\t139\t150",
  "\t-38\t139\t150",
  "\t-39\t139\t148",
  "\t-40\t145\t147",
  "Australia: Western Australia",
  "\t-12\t124\t128",
  "\t-13\t123\t130",
  "\t-14\t123\t130",
  "\t-15\t121\t130",
  "\t-16\t121\t130",
  "\t-17\t120\t130",
  "\t-18\t118\t130",
  "\t-19\t114\t130",
  "\t-20\t112\t130",
  "\t-21\t112\t130",
  "\t-22\t112\t130",
  "\t-23\t112\t130",
  "\t-24\t111\t130",
  "\t-25\t111\t130",
  "\t-26\t111\t130",
  "\t-27\t112\t130",
  "\t-28\t112\t130",
  "\t-29\t113\t130",
  "\t-30\t113\t130",
  "\t-31\t113\t130",
  "\t-32\t113\t130",
  "\t-33\t113\t129",
  "\t-34\t113\t125",
  "\t-35\t113\t124",
  "\t-36\t115\t119",
  "Austria",
  "\t50\t13\t16",
  "\t49\t11\t18",
  "\t48\t8\t18",
  "\t47\t8\t18",
  "\t46\t8\t18",
  "\t45\t8\t17",
  "Azerbaijan",
  "\t42\t43\t50",
  "\t41\t43\t51",
  "\t40\t43\t51",
  "\t39\t43\t51",
  "\t38\t43\t50",
  "\t37\t44\t50",
  "Bahamas",
  "\t27\t-79\t-76",
  "\t26\t-80\t-75",
  "\t25\t-80\t-73",
  "\t24\t-80\t-72",
  "\t23\t-80\t-71",
  "\t22\t-80\t-71",
  "\t21\t-76\t-71",
  "\t20\t-74\t-71",
  "\t19\t-74\t-72",
  "Bahrain",
  "\t27\t49\t51",
  "\t26\t49\t51",
  "\t25\t49\t51",
  "\t24\t49\t51",
  "Baker Island",
  "\t1\t-177\t-175",
  "\t0\t-177\t-175",
  "\t-1\t-177\t-175",
  "Bangladesh",
  "\t27\t87\t90",
  "\t26\t87\t93",
  "\t25\t87\t93",
  "\t24\t87\t93",
  "\t23\t87\t93",
  "\t22\t87\t93",
  "\t21\t87\t93",
  "\t20\t88\t93",
  "\t19\t91\t93",
  "Barbados",
  "\t14\t-60\t-58",
  "\t13\t-60\t-58",
  "\t12\t-60\t-58",
  "Bassas da India",
  "\t-20\t38\t40",
  "\t-21\t38\t40",
  "\t-22\t38\t40",
  "Belarus",
  "\t57\t26\t30",
  "\t56\t25\t32",
  "\t55\t23\t32",
  "\t54\t22\t33",
  "\t53\t22\t33",
  "\t52\t22\t33",
  "\t51\t22\t32",
  "\t50\t22\t31",
  "Belgium",
  "\t52\t1\t6",
  "\t51\t1\t7",
  "\t50\t1\t7",
  "\t49\t1\t7",
  "\t48\t3\t6",
  "Belize",
  "\t19\t-90\t-86",
  "\t18\t-90\t-86",
  "\t17\t-90\t-86",
  "\t16\t-90\t-86",
  "\t15\t-90\t-87",
  "\t14\t-90\t-87",
  "Benin",
  "\t13\t1\t4",
  "\t12\t-1\t4",
  "\t11\t-1\t4",
  "\t10\t-1\t4",
  "\t9\t-1\t4",
  "\t8\t0\t4",
  "\t7\t0\t3",
  "\t6\t0\t3",
  "\t5\t0\t3",
  "Bermuda",
  "\t33\t-65\t-63",
  "\t32\t-65\t-63",
  "\t31\t-65\t-63",
  "Bhutan",
  "\t29\t88\t92",
  "\t28\t87\t93",
  "\t27\t87\t93",
  "\t26\t87\t93",
  "\t25\t87\t93",
  "Bolivia",
  "\t-8\t-67\t-64",
  "\t-9\t-70\t-64",
  "\t-10\t-70\t-63",
  "\t-11\t-70\t-61",
  "\t-12\t-70\t-59",
  "\t-13\t-70\t-59",
  "\t-14\t-70\t-59",
  "\t-15\t-70\t-57",
  "\t-16\t-70\t-56",
  "\t-17\t-70\t-56",
  "\t-18\t-70\t-56",
  "\t-19\t-70\t-56",
  "\t-20\t-69\t-56",
  "\t-21\t-69\t-56",
  "\t-22\t-69\t-61",
  "\t-23\t-69\t-61",
  "Borneo",
  "\t6\t113\t116",
  "\t5\t113\t116",
  "\t4\t113\t116",
  "\t3\t113\t116",
  "Borneo",
  "\t5\t114\t118",
  "\t4\t107\t109\t114\t118",
  "\t3\t107\t110\t113\t119",
  "\t2\t107\t120",
  "\t1\t107\t120",
  "\t0\t107\t120",
  "\t-1\t107\t120",
  "\t-2\t107\t118",
  "\t-3\t109\t117",
  "\t-4\t109\t117",
  "\t-5\t113\t117",
  "Borneo",
  "\t8\t115\t118",
  "\t7\t115\t119",
  "\t6\t114\t120",
  "\t5\t112\t120",
  "\t4\t111\t120",
  "\t3\t108\t119",
  "\t2\t108\t116",
  "\t1\t108\t116",
  "\t0\t108\t115",
  "\t-1\t109\t112",
  "Bosnia and Herzegovina",
  "\t46\t14\t19",
  "\t45\t14\t20",
  "\t44\t14\t20",
  "\t43\t14\t20",
  "\t42\t15\t20",
  "\t41\t16\t19",
  "Botswana",
  "\t-16\t22\t26",
  "\t-17\t19\t26",
  "\t-18\t19\t27",
  "\t-19\t19\t28",
  "\t-20\t19\t30",
  "\t-21\t18\t30",
  "\t-22\t18\t30",
  "\t-23\t18\t30",
  "\t-24\t18\t28",
  "\t-25\t18\t27",
  "\t-26\t19\t26",
  "\t-27\t19\t23",
  "Bouvet Island",
  "\t-53\t2\t4",
  "\t-54\t2\t4",
  "\t-55\t2\t4",
  "Brazil",
  "\t6\t-61\t-58",
  "\t5\t-65\t-58\t-52\t-50",
  "\t4\t-65\t-58\t-53\t-49",
  "\t3\t-69\t-49",
  "\t2\t-70\t-48",
  "\t1\t-71\t-45\t-30\t-28",
  "\t0\t-71\t-43\t-30\t-28",
  "\t-1\t-71\t-38\t-30\t-28",
  "\t-2\t-70\t-37\t-33\t-31",
  "\t-3\t-73\t-35\t-33\t-31",
  "\t-4\t-74\t-31",
  "\t-5\t-74\t-33",
  "\t-6\t-75\t-33",
  "\t-7\t-75\t-33",
  "\t-8\t-75\t-33",
  "\t-9\t-74\t-33",
  "\t-10\t-74\t-34",
  "\t-11\t-73\t-35",
  "\t-12\t-71\t-36",
  "\t-13\t-65\t-36",
  "\t-14\t-63\t-37",
  "\t-15\t-61\t-37",
  "\t-16\t-61\t-37",
  "\t-17\t-61\t-37",
  "\t-18\t-59\t-38",
  "\t-19\t-59\t-38\t-30\t-27",
  "\t-20\t-59\t-38\t-30\t-27",
  "\t-21\t-59\t-39\t-30\t-27",
  "\t-22\t-58\t-39",
  "\t-23\t-58\t-39",
  "\t-24\t-56\t-42",
  "\t-25\t-56\t-45",
  "\t-26\t-56\t-46",
  "\t-27\t-57\t-47",
  "\t-28\t-58\t-47",
  "\t-29\t-58\t-47",
  "\t-30\t-58\t-48",
  "\t-31\t-58\t-49",
  "\t-32\t-57\t-49",
  "\t-33\t-54\t-51",
  "\t-34\t-54\t-51",
  "British Indian Ocean Territory",
  "\t-4\t70\t73",
  "\t-5\t70\t73",
  "\t-6\t70\t73",
  "\t-7\t70\t73",
  "\t-8\t71\t73",
  "British Virgin Islands",
  "\t19\t-65\t-63",
  "\t18\t-65\t-63",
  "\t17\t-65\t-63",
  "Brunei",
  "\t6\t113\t116",
  "\t5\t113\t116",
  "\t4\t113\t116",
  "\t3\t113\t116",
  "Bulgaria",
  "\t45\t21\t28",
  "\t44\t21\t29",
  "\t43\t21\t29",
  "\t42\t21\t29",
  "\t41\t21\t29",
  "\t40\t21\t29",
  "Burkina Faso",
  "\t16\t-1\t1",
  "\t15\t-3\t1",
  "\t14\t-5\t2",
  "\t13\t-5\t3",
  "\t12\t-6\t3",
  "\t11\t-6\t3",
  "\t10\t-6\t3",
  "\t9\t-6\t1",
  "\t8\t-5\t-1",
  "Burundi",
  "\t-1\t27\t31",
  "\t-2\t27\t31",
  "\t-3\t27\t31",
  "\t-4\t28\t31",
  "\t-5\t28\t31",
  "Cambodia",
  "\t15\t101\t108",
  "\t14\t101\t108",
  "\t13\t101\t108",
  "\t12\t101\t108",
  "\t11\t101\t108",
  "\t10\t101\t107",
  "\t9\t102\t107",
  "Cameroon",
  "\t14\t13\t15",
  "\t13\t13\t16",
  "\t12\t12\t16",
  "\t11\t12\t16",
  "\t10\t11\t16",
  "\t9\t11\t16",
  "\t8\t9\t16",
  "\t7\t8\t16",
  "\t6\t7\t16",
  "\t5\t7\t16",
  "\t4\t7\t16",
  "\t3\t7\t17",
  "\t2\t8\t17",
  "\t1\t8\t17",
  "\t0\t14\t17",
  "Canada",
  "\t84\t-78\t-67",
  "\t83\t-90\t-60",
  "\t82\t-96\t-60",
  "\t81\t-101\t-60",
  "\t80\t-106\t-61",
  "\t79\t-115\t-108\t-106\t-66",
  "\t78\t-121\t-69",
  "\t77\t-124\t-73",
  "\t76\t-124\t-74",
  "\t75\t-125\t-76",
  "\t74\t-125\t-75",
  "\t73\t-126\t-73",
  "\t72\t-126\t-69",
  "\t71\t-132\t-66",
  "\t70\t-142\t-65",
  "\t69\t-142\t-63",
  "\t68\t-142\t-61",
  "\t67\t-142\t-60",
  "\t66\t-142\t-60",
  "\t65\t-142\t-60",
  "\t64\t-142\t-61",
  "\t63\t-142\t-62",
  "\t62\t-142\t-63",
  "\t61\t-142\t-89\t-84\t-63",
  "\t60\t-142\t-91\t-81\t-62",
  "\t59\t-142\t-91\t-81\t-61",
  "\t58\t-140\t-88\t-81\t-60",
  "\t57\t-138\t-86\t-81\t-59",
  "\t56\t-133\t-57",
  "\t55\t-134\t-56",
  "\t54\t-134\t-54",
  "\t53\t-134\t-54",
  "\t52\t-134\t-54",
  "\t51\t-133\t-54",
  "\t50\t-132\t-52",
  "\t49\t-129\t-51",
  "\t48\t-128\t-51",
  "\t47\t-126\t-51",
  "\t46\t-90\t-51",
  "\t45\t-86\t-58\t-56\t-51",
  "\t44\t-84\t-58",
  "\t43\t-84\t-73\t-67\t-58",
  "\t42\t-84\t-75\t-67\t-63\t-61\t-58",
  "\t41\t-84\t-77",
  "\t40\t-84\t-80",
  "Canada: Alberta",
  "\t61\t-121\t-109",
  "\t60\t-121\t-109",
  "\t59\t-121\t-109",
  "\t58\t-121\t-109",
  "\t57\t-121\t-109",
  "\t56\t-121\t-109",
  "\t55\t-121\t-109",
  "\t54\t-121\t-109",
  "\t53\t-121\t-109",
  "\t52\t-121\t-109",
  "\t51\t-119\t-109",
  "\t50\t-118\t-109",
  "\t49\t-116\t-109",
  "\t48\t-115\t-109",
  "\t47\t-115\t-109",
  "Canada: British Columbia",
  "\t61\t-140\t-119",
  "\t60\t-140\t-119",
  "\t59\t-140\t-119",
  "\t58\t-140\t-119",
  "\t57\t-138\t-119",
  "\t56\t-134\t-119",
  "\t55\t-134\t-119",
  "\t54\t-134\t-117",
  "\t53\t-134\t-116",
  "\t52\t-134\t-114",
  "\t51\t-133\t-113",
  "\t50\t-132\t-113",
  "\t49\t-129\t-113",
  "\t48\t-128\t-113",
  "\t47\t-126\t-113",
  "Canada: Manitoba",
  "\t61\t-103\t-93",
  "\t60\t-103\t-93",
  "\t59\t-103\t-91",
  "\t58\t-103\t-88",
  "\t57\t-103\t-87",
  "\t56\t-103\t-87",
  "\t55\t-103\t-87",
  "\t54\t-103\t-89",
  "\t53\t-102\t-90",
  "\t52\t-102\t-92",
  "\t51\t-102\t-93",
  "\t50\t-102\t-94",
  "\t49\t-102\t-94",
  "\t48\t-102\t-94",
  "\t47\t-102\t-94",
  "Canada: New Brunswick",
  "\t49\t-67\t-63",
  "\t48\t-70\t-63",
  "\t47\t-70\t-62",
  "\t46\t-70\t-62",
  "\t45\t-68\t-62",
  "\t44\t-68\t-63",
  "\t43\t-67\t-65",
  "Canada: Newfoundland and Labrador",
  "\t61\t-65\t-63",
  "\t60\t-65\t-62",
  "\t59\t-65\t-61",
  "\t58\t-65\t-60",
  "\t57\t-65\t-59",
  "\t56\t-68\t-57",
  "\t55\t-68\t-56",
  "\t54\t-68\t-54",
  "\t53\t-68\t-54",
  "\t52\t-68\t-54",
  "\t51\t-68\t-54",
  "\t50\t-65\t-63\t-59\t-52",
  "\t49\t-60\t-51",
  "\t48\t-60\t-51",
  "\t47\t-60\t-51",
  "\t46\t-60\t-51",
  "\t45\t-56\t-51",
  "Canada: Northwest Territories",
  "\t79\t-115\t-109",
  "\t78\t-121\t-109",
  "\t77\t-124\t-109",
  "\t76\t-124\t-109",
  "\t75\t-125\t-109",
  "\t74\t-125\t-109",
  "\t73\t-126\t-109",
  "\t72\t-126\t-109",
  "\t71\t-132\t-109",
  "\t70\t-136\t-109",
  "\t69\t-137\t-109",
  "\t68\t-137\t-115\t-113\t-111",
  "\t67\t-137\t-113",
  "\t66\t-137\t-108",
  "\t65\t-135\t-100",
  "\t64\t-134\t-100",
  "\t63\t-133\t-100",
  "\t62\t-131\t-100",
  "\t61\t-130\t-101",
  "\t60\t-129\t-101",
  "\t59\t-127\t-101",
  "Canada: Nova Scotia",
  "\t48\t-61\t-59",
  "\t47\t-65\t-58",
  "\t46\t-66\t-58",
  "\t45\t-67\t-58",
  "\t44\t-67\t-58",
  "\t43\t-67\t-58",
  "\t42\t-67\t-63\t-61\t-58",
  "Canada: Nunavut",
  "\t84\t-78\t-67",
  "\t83\t-90\t-60",
  "\t82\t-96\t-60",
  "\t81\t-101\t-60",
  "\t80\t-106\t-61",
  "\t79\t-111\t-108\t-106\t-66",
  "\t78\t-111\t-69",
  "\t77\t-111\t-73",
  "\t76\t-111\t-74",
  "\t75\t-111\t-76",
  "\t74\t-111\t-75",
  "\t73\t-111\t-73",
  "\t72\t-111\t-69",
  "\t71\t-118\t-66",
  "\t70\t-121\t-65",
  "\t69\t-121\t-63",
  "\t68\t-121\t-61",
  "\t67\t-121\t-60",
  "\t66\t-121\t-60",
  "\t65\t-118\t-60",
  "\t64\t-114\t-61",
  "\t63\t-110\t-62",
  "\t62\t-103\t-63",
  "\t61\t-103\t-89\t-84\t-63",
  "\t60\t-103\t-91\t-81\t-77\t-69\t-63",
  "\t59\t-103\t-93\t-81\t-76\t-69\t-63",
  "\t58\t-81\t-75",
  "\t57\t-81\t-75",
  "\t56\t-81\t-75",
  "\t55\t-82\t-75",
  "\t54\t-83\t-76",
  "\t53\t-83\t-77",
  "\t52\t-83\t-77",
  "\t51\t-82\t-77",
  "\t50\t-80\t-78",
  "Canada: Ontario",
  "\t57\t-90\t-86",
  "\t56\t-92\t-81",
  "\t55\t-94\t-81",
  "\t54\t-95\t-81",
  "\t53\t-96\t-79",
  "\t52\t-96\t-78",
  "\t51\t-96\t-78",
  "\t50\t-96\t-78",
  "\t49\t-96\t-78",
  "\t48\t-96\t-78",
  "\t47\t-95\t-76",
  "\t46\t-90\t-73",
  "\t45\t-86\t-73",
  "\t44\t-84\t-73",
  "\t43\t-84\t-73",
  "\t42\t-84\t-75",
  "\t41\t-84\t-77",
  "\t40\t-84\t-80",
  "Canada: Prince Edward Island",
  "\t48\t-65\t-62",
  "\t47\t-65\t-60",
  "\t46\t-65\t-60",
  "\t45\t-65\t-60",
  "\t44\t-63\t-61",
  "Canada: Quebec",
  "\t63\t-79\t-71",
  "\t62\t-79\t-68",
  "\t61\t-79\t-68\t-66\t-63",
  "\t60\t-79\t-68\t-66\t-63",
  "\t59\t-79\t-62",
  "\t58\t-79\t-62",
  "\t57\t-79\t-62",
  "\t56\t-79\t-62",
  "\t55\t-80\t-62",
  "\t54\t-80\t-62",
  "\t53\t-80\t-56",
  "\t52\t-80\t-56",
  "\t51\t-80\t-56",
  "\t50\t-80\t-56",
  "\t49\t-80\t-57",
  "\t48\t-80\t-60",
  "\t47\t-80\t-60",
  "\t46\t-80\t-65\t-63\t-60",
  "\t45\t-80\t-68",
  "\t44\t-78\t-69",
  "\t43\t-75\t-73",
  "Canada: Saskatchewan",
  "\t61\t-111\t-101",
  "\t60\t-111\t-101",
  "\t59\t-111\t-101",
  "\t58\t-111\t-101",
  "\t57\t-111\t-101",
  "\t56\t-111\t-100",
  "\t55\t-111\t-100",
  "\t54\t-111\t-100",
  "\t53\t-111\t-100",
  "\t52\t-111\t-100",
  "\t51\t-111\t-100",
  "\t50\t-111\t-100",
  "\t49\t-111\t-100",
  "\t48\t-111\t-100",
  "\t47\t-111\t-100",
  "Canada: Yukon",
  "\t70\t-142\t-136",
  "\t69\t-142\t-135",
  "\t68\t-142\t-132",
  "\t67\t-142\t-131",
  "\t66\t-142\t-131",
  "\t65\t-142\t-129",
  "\t64\t-142\t-128",
  "\t63\t-142\t-127",
  "\t62\t-142\t-125",
  "\t61\t-142\t-122",
  "\t60\t-142\t-122",
  "\t59\t-142\t-122",
  "Cape Verde",
  "\t18\t-26\t-23",
  "\t17\t-26\t-21",
  "\t16\t-26\t-21",
  "\t15\t-26\t-21",
  "\t14\t-25\t-21",
  "\t13\t-25\t-22",
  "Cayman Islands",
  "\t20\t-82\t-78",
  "\t19\t-82\t-78",
  "\t18\t-82\t-78",
  "Central African Republic",
  "\t12\t21\t23",
  "\t11\t20\t24",
  "\t10\t18\t24",
  "\t9\t17\t25",
  "\t8\t14\t26",
  "\t7\t13\t27",
  "\t6\t13\t28",
  "\t5\t13\t28",
  "\t4\t13\t28",
  "\t3\t13\t26",
  "\t2\t14\t19",
  "\t1\t14\t17",
  "Chad",
  "\t24\t13\t17",
  "\t23\t13\t19",
  "\t22\t13\t21",
  "\t21\t13\t23",
  "\t20\t14\t24",
  "\t19\t14\t24",
  "\t18\t14\t24",
  "\t17\t13\t24",
  "\t16\t12\t24",
  "\t15\t12\t24",
  "\t14\t12\t24",
  "\t13\t12\t23",
  "\t12\t12\t23",
  "\t11\t13\t23",
  "\t10\t12\t23",
  "\t9\t12\t23",
  "\t8\t12\t22",
  "\t7\t13\t20",
  "\t6\t14\t18",
  "Chile",
  "\t-16\t-70\t-68",
  "\t-17\t-71\t-67",
  "\t-18\t-71\t-67",
  "\t-19\t-71\t-67",
  "\t-20\t-71\t-67",
  "\t-21\t-71\t-66",
  "\t-22\t-71\t-66",
  "\t-23\t-71\t-66",
  "\t-24\t-71\t-66",
  "\t-25\t-106\t-104\t-81\t-78\t-71\t-66",
  "\t-26\t-110\t-108\t-106\t-104\t-81\t-78\t-72\t-67",
  "\t-27\t-110\t-108\t-106\t-104\t-81\t-78\t-72\t-67",
  "\t-28\t-110\t-108\t-72\t-67",
  "\t-29\t-72\t-68",
  "\t-30\t-72\t-68",
  "\t-31\t-72\t-68",
  "\t-32\t-81\t-77\t-72\t-68",
  "\t-33\t-81\t-77\t-73\t-68",
  "\t-34\t-81\t-77\t-73\t-68",
  "\t-35\t-74\t-68",
  "\t-36\t-74\t-69",
  "\t-37\t-74\t-69",
  "\t-38\t-74\t-69",
  "\t-39\t-74\t-69",
  "\t-40\t-75\t-70",
  "\t-41\t-75\t-70",
  "\t-42\t-75\t-70",
  "\t-43\t-76\t-70",
  "\t-44\t-76\t-70",
  "\t-45\t-76\t-70",
  "\t-46\t-76\t-70",
  "\t-47\t-76\t-70",
  "\t-48\t-76\t-70",
  "\t-49\t-76\t-71",
  "\t-50\t-76\t-69",
  "\t-51\t-76\t-67",
  "\t-52\t-76\t-67",
  "\t-53\t-76\t-66",
  "\t-54\t-75\t-65",
  "\t-55\t-74\t-65",
  "\t-56\t-72\t-65",
  "China",
  "\t54\t119\t126",
  "\t53\t119\t127",
  "\t52\t118\t127",
  "\t51\t118\t128",
  "\t50\t85\t88\t115\t130",
  "\t49\t84\t90\t114\t135",
  "\t48\t81\t91\t114\t135",
  "\t47\t81\t92\t114\t135",
  "\t46\t79\t94\t110\t135",
  "\t45\t78\t96\t110\t135",
  "\t44\t78\t96\t109\t134",
  "\t43\t78\t132",
  "\t42\t75\t132",
  "\t41\t72\t132",
  "\t40\t72\t129",
  "\t39\t72\t127",
  "\t38\t72\t125",
  "\t37\t72\t123",
  "\t36\t73\t123",
  "\t35\t73\t123",
  "\t34\t74\t121",
  "\t33\t77\t122",
  "\t32\t77\t122",
  "\t31\t77\t123",
  "\t30\t77\t123",
  "\t29\t78\t123",
  "\t28\t81\t123",
  "\t27\t83\t122",
  "\t26\t84\t93\t96\t122",
  "\t25\t96\t121",
  "\t24\t96\t120",
  "\t23\t96\t120",
  "\t22\t96\t118",
  "\t21\t98\t117",
  "\t20\t98\t102\t105\t114",
  "\t19\t107\t112",
  "\t18\t107\t112",
  "\t17\t107\t111",
  "China: Hainan",
  "\t21\t108\t111",
  "\t20\t107\t112",
  "\t19\t107\t112",
  "\t18\t107\t112",
  "\t17\t107\t111",
  "Christmas Island",
  "\t-9\t104\t106",
  "\t-10\t104\t106",
  "\t-11\t104\t106",
  "Clipperton Island",
  "\t11\t-110\t-108",
  "\t10\t-110\t-108",
  "\t9\t-110\t-108",
  "Cocos Islands",
  "\t-11\t95\t97",
  "\t-12\t95\t97",
  "\t-13\t95\t97",
  "Colombia",
  "\t14\t-82\t-80",
  "\t13\t-82\t-80\t-73\t-70",
  "\t12\t-82\t-80\t-75\t-70",
  "\t11\t-82\t-80\t-76\t-70",
  "\t10\t-77\t-70",
  "\t9\t-78\t-71",
  "\t8\t-78\t-69",
  "\t7\t-78\t-66",
  "\t6\t-78\t-66",
  "\t5\t-78\t-66",
  "\t4\t-78\t-66",
  "\t3\t-79\t-66",
  "\t2\t-80\t-65",
  "\t1\t-80\t-65",
  "\t0\t-80\t-65",
  "\t-1\t-79\t-68",
  "\t-2\t-75\t-68",
  "\t-3\t-74\t-68",
  "\t-4\t-71\t-68",
  "\t-5\t-71\t-68",
  "Comoros",
  "\t-10\t42\t44",
  "\t-11\t42\t45",
  "\t-12\t42\t45",
  "\t-13\t42\t45",
  "Cook Islands",
  "\t-7\t-159\t-156",
  "\t-8\t-159\t-156",
  "\t-9\t-166\t-164\t-162\t-156",
  "\t-10\t-166\t-164\t-162\t-159",
  "\t-11\t-166\t-164\t-162\t-159",
  "\t-17\t-160\t-158",
  "\t-18\t-160\t-156",
  "\t-19\t-160\t-156",
  "\t-20\t-160\t-156",
  "\t-21\t-160\t-156",
  "\t-22\t-160\t-156",
  "Coral Sea Islands",
  "\t-15\t146\t151",
  "\t-16\t146\t151",
  "\t-17\t146\t151",
  "\t-18\t147\t149",
  "\t-20\t152\t156",
  "\t-21\t152\t156",
  "\t-22\t152\t156",
  "\t-23\t154\t156",
  "Costa Rica",
  "\t12\t-86\t-83",
  "\t11\t-86\t-82",
  "\t10\t-86\t-81",
  "\t9\t-86\t-81",
  "\t8\t-86\t-81",
  "\t7\t-84\t-81",
  "\t6\t-88\t-86",
  "\t5\t-88\t-86",
  "\t4\t-88\t-86",
  "Cote d'Ivoire",
  "\t11\t-9\t-3",
  "\t10\t-9\t-1",
  "\t9\t-9\t-1",
  "\t8\t-9\t-1",
  "\t7\t-9\t-1",
  "\t6\t-9\t-1",
  "\t5\t-9\t-1",
  "\t4\t-8\t-1",
  "\t3\t-8\t-4",
  "Croatia",
  "\t47\t14\t18",
  "\t46\t12\t20",
  "\t45\t12\t20",
  "\t44\t12\t20",
  "\t43\t12\t20",
  "\t42\t14\t19",
  "\t41\t15\t19",
  "Cuba",
  "\t24\t-84\t-79",
  "\t23\t-85\t-76",
  "\t22\t-85\t-74",
  "\t21\t-85\t-73",
  "\t20\t-85\t-73",
  "\t19\t-80\t-73",
  "\t18\t-78\t-73",
  "Curacao",
  "\t13\t-70\t-67",
  "\t12\t-70\t-67",
  "\t11\t-70\t-67",
  "Cyprus",
  "\t36\t31\t35",
  "\t35\t31\t35",
  "\t34\t31\t35",
  "\t33\t31\t35",
  "Cyprus",
  "\t36\t31\t35",
  "\t35\t31\t35",
  "\t34\t31\t35",
  "Cyprus",
  "\t35\t31\t34",
  "\t34\t31\t34",
  "\t33\t31\t34",
  "Cyprus",
  "\t36\t32\t34",
  "\t35\t32\t34",
  "\t34\t32\t34",
  "\t33\t32\t34",
  "Czech Republic",
  "\t52\t13\t16",
  "\t51\t11\t19",
  "\t50\t11\t19",
  "\t49\t11\t19",
  "\t48\t11\t19",
  "\t47\t12\t18",
  "Democratic Republic of the Congo",
  "\t6\t18\t20\t23\t28",
  "\t5\t17\t31",
  "\t4\t17\t31",
  "\t3\t17\t32",
  "\t2\t16\t32",
  "\t1\t16\t32",
  "\t0\t15\t32",
  "\t-1\t15\t31",
  "\t-2\t14\t30",
  "\t-3\t11\t30",
  "\t-4\t11\t30",
  "\t-5\t11\t31",
  "\t-6\t11\t31",
  "\t-7\t11\t13\t15\t31",
  "\t-8\t15\t31",
  "\t-9\t16\t31",
  "\t-10\t20\t29",
  "\t-11\t21\t30",
  "\t-12\t21\t30",
  "\t-13\t25\t30",
  "\t-14\t27\t30",
  "Denmark",
  "\t58\t7\t12",
  "\t57\t7\t13",
  "\t56\t7\t16",
  "\t55\t7\t16",
  "\t54\t7\t16",
  "\t53\t7\t13",
  "Djibouti",
  "\t13\t41\t44",
  "\t12\t40\t44",
  "\t11\t40\t44",
  "\t10\t40\t44",
  "\t9\t40\t43",
  "Dominica",
  "\t16\t-62\t-60",
  "\t15\t-62\t-60",
  "\t14\t-62\t-60",
  "Dominican Republic",
  "\t20\t-72\t-67",
  "\t19\t-73\t-67",
  "\t18\t-73\t-67",
  "\t17\t-73\t-67",
  "\t16\t-72\t-70",
  "East Timor",
  "\t-7\t123\t128",
  "\t-8\t123\t128",
  "\t-9\t123\t128",
  "\t-10\t123\t127",
  "Ecuador",
  "\t2\t-80\t-77",
  "\t1\t-81\t-74",
  "\t0\t-81\t-74",
  "\t-1\t-82\t-74",
  "\t-2\t-82\t-74",
  "\t-3\t-82\t-74",
  "\t-4\t-81\t-76",
  "\t-5\t-81\t-77",
  "\t-6\t-80\t-78",
  "Ecuador: Galapagos",
  "\t2\t-93\t-90",
  "\t1\t-93\t-88",
  "\t0\t-93\t-88",
  "\t-1\t-92\t-88",
  "\t-2\t-92\t-88",
  "Egypt",
  "\t32\t23\t35",
  "\t31\t23\t35",
  "\t30\t23\t35",
  "\t29\t23\t35",
  "\t28\t23\t35",
  "\t27\t23\t35",
  "\t26\t23\t35",
  "\t25\t23\t36",
  "\t24\t23\t36",
  "\t23\t23\t37",
  "\t22\t23\t37",
  "\t21\t23\t37",
  "\t20\t23\t37",
  "El Salvador",
  "\t15\t-90\t-87",
  "\t14\t-91\t-86",
  "\t13\t-91\t-86",
  "\t12\t-91\t-86",
  "Equatorial Guinea",
  "\t4\t7\t9",
  "\t3\t7\t12",
  "\t2\t7\t12",
  "\t1\t8\t12",
  "\t0\t4\t6\t8\t12",
  "\t-1\t4\t6\t8\t10",
  "\t-2\t4\t6",
  "Eritrea",
  "\t19\t37\t39",
  "\t18\t35\t40",
  "\t17\t35\t41",
  "\t16\t35\t41",
  "\t15\t35\t42",
  "\t14\t35\t43",
  "\t13\t35\t44",
  "\t12\t39\t44",
  "\t11\t40\t44",
  "Estonia",
  "\t60\t21\t29",
  "\t59\t20\t29",
  "\t58\t20\t29",
  "\t57\t20\t28",
  "\t56\t20\t28",
  "Ethiopia",
  "\t15\t35\t41",
  "\t14\t35\t42",
  "\t13\t34\t43",
  "\t12\t33\t43",
  "\t11\t33\t44",
  "\t10\t33\t44",
  "\t9\t32\t47",
  "\t8\t31\t48",
  "\t7\t31\t48",
  "\t6\t31\t48",
  "\t5\t33\t47",
  "\t4\t33\t46",
  "\t3\t34\t46",
  "\t2\t36\t42",
  "Europa Island",
  "\t-21\t39\t41",
  "\t-22\t39\t41",
  "\t-23\t39\t41",
  "Falkland Islands (Islas Malvinas)",
  "\t-50\t-62\t-56",
  "\t-51\t-62\t-56",
  "\t-52\t-62\t-56",
  "\t-53\t-62\t-57",
  "Faroe Islands",
  "\t63\t-8\t-5",
  "\t62\t-8\t-5",
  "\t61\t-8\t-5",
  "\t60\t-7\t-5",
  "Fiji",
  "\t-11\t176\t178",
  "\t-12\t176\t178",
  "\t-13\t176\t178",
  "\t-15\t-180\t-178\t176\t180",
  "\t-16\t-180\t-177\t176\t180",
  "\t-17\t-180\t-177\t176\t180",
  "\t-18\t-180\t-177\t176\t180",
  "\t-19\t-180\t-177\t176\t180",
  "\t-20\t-180\t-177\t173\t180",
  "\t-21\t173\t175",
  "\t-22\t173\t175",
  "Finland",
  "\t71\t26\t28",
  "\t70\t19\t30",
  "\t69\t19\t30",
  "\t68\t19\t31",
  "\t67\t19\t31",
  "\t66\t22\t31",
  "\t65\t22\t31",
  "\t64\t20\t32",
  "\t63\t20\t32",
  "\t62\t20\t32",
  "\t61\t20\t32",
  "\t60\t20\t31",
  "\t59\t20\t29",
  "\t58\t21\t25",
  "Finland",
  "\t61\t18\t22",
  "\t60\t18\t22",
  "\t59\t18\t22",
  "\t58\t19\t21",
  "France",
  "\t52\t0\t3",
  "\t51\t0\t5",
  "\t50\t-2\t8",
  "\t49\t-6\t9",
  "\t48\t-6\t9",
  "\t47\t-6\t9",
  "\t46\t-5\t8",
  "\t45\t-3\t8",
  "\t44\t-2\t8",
  "\t43\t-2\t8",
  "\t42\t-2\t8",
  "\t41\t-2\t7",
  "France: Corsica",
  "\t44\t8\t10",
  "\t43\t7\t10",
  "\t42\t7\t10",
  "\t41\t7\t10",
  "\t40\t7\t10",
  "France: Saint Barthelemy",
  "\t18\t-63\t-61",
  "\t17\t-63\t-61",
  "\t16\t-63\t-61",
  "France: Saint Martin",
  "\t19\t-64\t-62",
  "\t18\t-64\t-62",
  "\t17\t-64\t-62",
  "French Guiana",
  "\t6\t-55\t-51",
  "\t5\t-55\t-50",
  "\t4\t-55\t-50",
  "\t3\t-55\t-50",
  "\t2\t-55\t-50",
  "\t1\t-55\t-51",
  "French Polynesia",
  "\t-6\t-141\t-139",
  "\t-7\t-141\t-138",
  "\t-8\t-141\t-137",
  "\t-9\t-141\t-137",
  "\t-10\t-141\t-137",
  "\t-11\t-140\t-137",
  "\t-13\t-149\t-140",
  "\t-14\t-149\t-139",
  "\t-15\t-152\t-139",
  "\t-16\t-152\t-137",
  "\t-17\t-152\t-135",
  "\t-18\t-150\t-148\t-146\t-135",
  "\t-19\t-142\t-135",
  "\t-20\t-142\t-134",
  "\t-21\t-152\t-150\t-141\t-134",
  "\t-22\t-152\t-146\t-141\t-133",
  "\t-23\t-152\t-146\t-136\t-133",
  "\t-24\t-150\t-146\t-136\t-133",
  "\t-26\t-145\t-143",
  "\t-27\t-145\t-143",
  "\t-28\t-145\t-143",
  "French Southern and Antarctic Lands",
  "\t-10\t46\t48",
  "\t-11\t46\t48",
  "\t-12\t46\t48",
  "\t-14\t53\t55",
  "\t-15\t53\t55",
  "\t-16\t41\t43\t53\t55",
  "\t-17\t41\t43",
  "\t-18\t41\t43",
  "\t-20\t38\t40",
  "\t-21\t38\t41",
  "\t-22\t38\t41",
  "\t-23\t39\t41",
  "\t-36\t76\t78",
  "\t-37\t76\t78",
  "\t-38\t76\t78",
  "\t-39\t76\t78",
  "\t-45\t49\t52",
  "\t-46\t49\t52",
  "\t-47\t49\t52\t67\t70",
  "\t-48\t67\t71",
  "\t-49\t67\t71",
  "\t-50\t67\t71",
  "Gabon",
  "\t3\t10\t14",
  "\t2\t8\t15",
  "\t1\t7\t15",
  "\t0\t7\t15",
  "\t-1\t7\t15",
  "\t-2\t7\t15",
  "\t-3\t8\t15",
  "\t-4\t9\t12",
  "Gambia",
  "\t14\t-17\t-12",
  "\t13\t-17\t-12",
  "\t12\t-17\t-12",
  "Gaza Strip",
  "\t32\t33\t35",
  "\t31\t33\t35",
  "\t30\t33\t35",
  "Georgia",
  "\t44\t38\t44",
  "\t43\t38\t47",
  "\t42\t38\t47",
  "\t41\t39\t47",
  "\t40\t40\t47",
  "Germany",
  "\t56\t7\t9",
  "\t55\t7\t15",
  "\t54\t5\t15",
  "\t53\t5\t15",
  "\t52\t4\t16",
  "\t51\t4\t16",
  "\t50\t4\t16",
  "\t49\t4\t15",
  "\t48\t5\t14",
  "\t47\t6\t14",
  "\t46\t6\t14",
  "Ghana",
  "\t12\t-2\t1",
  "\t11\t-3\t1",
  "\t10\t-3\t1",
  "\t9\t-3\t1",
  "\t8\t-4\t1",
  "\t7\t-4\t2",
  "\t6\t-4\t2",
  "\t5\t-4\t2",
  "\t4\t-4\t2",
  "\t3\t-3\t0",
  "Gibraltar",
  "\t37\t-6\t-4",
  "\t36\t-6\t-4",
  "\t35\t-6\t-4",
  "Glorioso Islands",
  "\t-10\t46\t48",
  "\t-11\t46\t48",
  "\t-12\t46\t48",
  "Greece",
  "\t42\t20\t27",
  "\t41\t19\t27",
  "\t40\t18\t27",
  "\t39\t18\t27",
  "\t38\t18\t28",
  "\t37\t19\t29",
  "\t36\t19\t29",
  "\t35\t20\t29",
  "\t34\t22\t28",
  "\t33\t23\t26",
  "Greenland",
  "\t84\t-47\t-23",
  "\t83\t-60\t-18",
  "\t82\t-65\t-10",
  "\t81\t-68\t-10",
  "\t80\t-69\t-10",
  "\t79\t-74\t-13",
  "\t78\t-74\t-16",
  "\t77\t-74\t-16",
  "\t76\t-73\t-16",
  "\t75\t-72\t-16",
  "\t74\t-68\t-65\t-61\t-16",
  "\t73\t-58\t-16",
  "\t72\t-57\t-19",
  "\t71\t-57\t-20",
  "\t70\t-56\t-20",
  "\t69\t-56\t-20",
  "\t68\t-55\t-21",
  "\t67\t-54\t-24",
  "\t66\t-54\t-31",
  "\t65\t-54\t-32",
  "\t64\t-54\t-34",
  "\t63\t-53\t-39",
  "\t62\t-52\t-39",
  "\t61\t-51\t-40",
  "\t60\t-50\t-41",
  "\t59\t-49\t-41",
  "\t58\t-45\t-42",
  "Grenada",
  "\t13\t-62\t-60",
  "\t12\t-62\t-60",
  "\t11\t-62\t-60",
  "Guadeloupe",
  "\t17\t-62\t-59",
  "\t16\t-62\t-59",
  "\t15\t-62\t-59",
  "\t14\t-62\t-60",
  "Guam",
  "\t14\t143\t145",
  "\t13\t143\t145",
  "\t12\t143\t145",
  "Guatemala",
  "\t18\t-92\t-88",
  "\t17\t-92\t-88",
  "\t16\t-93\t-87",
  "\t15\t-93\t-87",
  "\t14\t-93\t-87",
  "\t13\t-93\t-88",
  "\t12\t-92\t-89",
  "Guernsey",
  "\t50\t-3\t-1",
  "\t49\t-3\t-1",
  "\t48\t-3\t-1",
  "Guinea",
  "\t13\t-14\t-7",
  "\t12\t-15\t-7",
  "\t11\t-16\t-6",
  "\t10\t-16\t-6",
  "\t9\t-16\t-6",
  "\t8\t-14\t-6",
  "\t7\t-11\t-6",
  "\t6\t-10\t-7",
  "Guinea-Bissau",
  "\t13\t-17\t-12",
  "\t12\t-17\t-12",
  "\t11\t-17\t-12",
  "\t10\t-17\t-12",
  "\t9\t-16\t-13",
  "Guyana",
  "\t9\t-61\t-58",
  "\t8\t-61\t-57",
  "\t7\t-62\t-56",
  "\t6\t-62\t-56",
  "\t5\t-62\t-56",
  "\t4\t-62\t-56",
  "\t3\t-61\t-55",
  "\t2\t-61\t-55",
  "\t1\t-61\t-55",
  "\t0\t-60\t-55",
  "Haiti",
  "\t21\t-73\t-71",
  "\t20\t-74\t-70",
  "\t19\t-75\t-70",
  "\t18\t-75\t-70",
  "\t17\t-75\t-70",
  "Heard Island and McDonald Islands",
  "\t-51\t72\t74",
  "\t-52\t72\t74",
  "\t-53\t72\t74",
  "\t-54\t72\t74",
  "Honduras",
  "\t18\t-84\t-82",
  "\t17\t-87\t-82",
  "\t16\t-90\t-82",
  "\t15\t-90\t-82",
  "\t14\t-90\t-82",
  "\t13\t-90\t-82",
  "\t12\t-89\t-84",
  "\t11\t-88\t-86",
  "Hong Kong",
  "\t23\t112\t115",
  "\t22\t112\t115",
  "\t21\t112\t115",
  "Howland Island",
  "\t1\t-177\t-175",
  "\t0\t-177\t-175",
  "\t-1\t-177\t-175",
  "Hungary",
  "\t49\t16\t23",
  "\t48\t15\t23",
  "\t47\t15\t23",
  "\t46\t15\t23",
  "\t45\t15\t22",
  "\t44\t16\t20",
  "Iceland",
  "\t67\t-24\t-13",
  "\t66\t-25\t-12",
  "\t65\t-25\t-12",
  "\t64\t-25\t-12",
  "\t63\t-25\t-12",
  "\t62\t-23\t-15",
  "India",
  "\t36\t76\t79",
  "\t35\t72\t79",
  "\t34\t72\t80",
  "\t33\t72\t80",
  "\t32\t72\t80",
  "\t31\t72\t82",
  "\t30\t71\t82\t93\t97",
  "\t29\t69\t82\t87\t89\t91\t98",
  "\t28\t68\t98",
  "\t27\t68\t98",
  "\t26\t68\t98",
  "\t25\t67\t96",
  "\t24\t67\t96",
  "\t23\t67\t95",
  "\t22\t67\t95",
  "\t21\t67\t94",
  "\t20\t68\t93",
  "\t19\t69\t88",
  "\t18\t71\t87",
  "\t17\t71\t85",
  "\t16\t72\t84",
  "\t15\t72\t83",
  "\t14\t72\t82\t91\t95",
  "\t13\t73\t81\t91\t95",
  "\t12\t71\t81\t91\t95",
  "\t11\t71\t81\t91\t94",
  "\t10\t71\t80\t91\t94",
  "\t9\t71\t80\t91\t94",
  "\t8\t72\t80\t91\t94",
  "\t7\t72\t79\t92\t94",
  "\t6\t92\t94",
  "\t5\t92\t94",
  "Indonesia",
  "\t6\t94\t98\t125\t127",
  "\t5\t94\t99\t106\t109\t114\t118\t125\t128",
  "\t4\t94\t100\t104\t109\t114\t118\t124\t128",
  "\t3\t94\t102\t104\t110\t113\t119\t124\t129",
  "\t2\t94\t132",
  "\t1\t94\t137",
  "\t0\t96\t139",
  "\t-1\t96\t141",
  "\t-2\t97\t141",
  "\t-3\t98\t141",
  "\t-4\t99\t141",
  "\t-5\t101\t141",
  "\t-6\t101\t116\t118\t141",
  "\t-7\t104\t141",
  "\t-8\t105\t132\t136\t141",
  "\t-9\t109\t132\t136\t141",
  "\t-10\t115\t126\t139\t141",
  "\t-11\t119\t125",
  "Iran",
  "\t40\t43\t49",
  "\t39\t43\t49\t54\t58",
  "\t38\t43\t61",
  "\t37\t43\t62",
  "\t36\t43\t62",
  "\t35\t43\t62",
  "\t34\t44\t62",
  "\t33\t44\t62",
  "\t32\t44\t62",
  "\t31\t45\t62",
  "\t30\t46\t62",
  "\t29\t46\t63",
  "\t28\t47\t64",
  "\t27\t49\t64",
  "\t26\t50\t64",
  "\t25\t52\t64",
  "\t24\t53\t62",
  "Iraq",
  "\t38\t41\t45",
  "\t37\t40\t46",
  "\t36\t40\t47",
  "\t35\t39\t47",
  "\t34\t37\t47",
  "\t33\t37\t48",
  "\t32\t37\t48",
  "\t31\t37\t49",
  "\t30\t39\t49",
  "\t29\t41\t49",
  "\t28\t42\t49",
  "Ireland",
  "\t56\t-9\t-5",
  "\t55\t-11\t-5",
  "\t54\t-11\t-5",
  "\t53\t-11\t-4",
  "\t52\t-11\t-4",
  "\t51\t-11\t-4",
  "\t50\t-11\t-6",
  "Isle of Man",
  "\t55\t-5\t-3",
  "\t54\t-5\t-3",
  "\t53\t-5\t-3",
  "Israel",
  "\t34\t34\t36",
  "\t33\t33\t36",
  "\t32\t33\t36",
  "\t31\t33\t36",
  "\t30\t33\t36",
  "\t29\t33\t36",
  "\t28\t33\t36",
  "Italy",
  "\t48\t10\t13",
  "\t47\t6\t14",
  "\t46\t5\t14",
  "\t45\t5\t14",
  "\t44\t5\t14",
  "\t43\t5\t16",
  "\t42\t6\t18",
  "\t41\t7\t19",
  "\t40\t7\t19",
  "\t39\t7\t19",
  "\t38\t7\t19",
  "\t37\t7\t18",
  "\t36\t10\t17",
  "\t35\t10\t16",
  "\t34\t11\t13",
  "Jamaica",
  "\t19\t-79\t-75",
  "\t18\t-79\t-75",
  "\t17\t-79\t-75",
  "\t16\t-78\t-75",
  "Jan Mayen",
  "\t72\t-9\t-6",
  "\t71\t-10\t-6",
  "\t70\t-10\t-6",
  "\t69\t-10\t-7",
  "Japan",
  "\t46\t139\t143",
  "\t45\t139\t146",
  "\t44\t139\t146",
  "\t43\t138\t146",
  "\t42\t138\t146",
  "\t41\t138\t146",
  "\t40\t138\t144",
  "\t39\t137\t143",
  "\t38\t135\t143",
  "\t37\t131\t142",
  "\t36\t131\t142",
  "\t35\t128\t141",
  "\t34\t128\t141",
  "\t33\t127\t140",
  "\t32\t127\t141",
  "\t31\t127\t134\t138\t141",
  "\t30\t128\t132\t139\t141",
  "\t29\t128\t132\t139\t141",
  "\t28\t126\t131\t139\t143",
  "\t27\t125\t131\t139\t143",
  "\t26\t122\t132\t139\t143",
  "\t25\t121\t132\t140\t143\t152\t154",
  "\t24\t121\t126\t130\t132\t140\t142\t152\t154",
  "\t23\t121\t126\t140\t142\t152\t154",
  "Jarvis Island",
  "\t1\t-161\t-159",
  "\t0\t-161\t-159",
  "\t-1\t-161\t-159",
  "Jersey",
  "\t50\t-3\t-1",
  "\t49\t-3\t-1",
  "\t48\t-3\t-1",
  "Johnston Atoll",
  "\t17\t-170\t-168",
  "\t16\t-170\t-168",
  "\t15\t-170\t-168",
  "Jordan",
  "\t34\t37\t39",
  "\t33\t34\t40",
  "\t32\t34\t40",
  "\t31\t34\t40",
  "\t30\t33\t39",
  "\t29\t33\t38",
  "\t28\t33\t38",
  "Juan de Nova Island",
  "\t-16\t41\t43",
  "\t-17\t41\t43",
  "\t-18\t41\t43",
  "Kazakhstan",
  "\t56\t67\t71",
  "\t55\t60\t77",
  "\t54\t59\t79",
  "\t53\t59\t79",
  "\t52\t48\t84",
  "\t51\t46\t86",
  "\t50\t45\t88",
  "\t49\t45\t88",
  "\t48\t45\t88",
  "\t47\t45\t87",
  "\t46\t46\t86",
  "\t45\t47\t86",
  "\t44\t48\t83",
  "\t43\t48\t56\t58\t81",
  "\t42\t49\t56\t60\t81",
  "\t41\t50\t56\t64\t81",
  "\t40\t51\t56\t65\t71",
  "\t39\t66\t69",
  "Kenya",
  "\t6\t33\t36",
  "\t5\t32\t42",
  "\t4\t32\t42",
  "\t3\t32\t42",
  "\t2\t33\t42",
  "\t1\t32\t42",
  "\t0\t32\t42",
  "\t-1\t32\t42",
  "\t-2\t32\t42",
  "\t-3\t34\t42",
  "\t-4\t36\t41",
  "\t-5\t37\t40",
  "Kerguelen Archipelago",
  "\t-47\t67\t70",
  "\t-48\t67\t71",
  "\t-49\t67\t71",
  "\t-50\t67\t71",
  "Kingman Reef",
  "\t7\t-163\t-161",
  "\t6\t-163\t-161",
  "\t5\t-163\t-161",
  "Kiribati",
  "\t5\t-161\t-159",
  "\t4\t-161\t-158\t171\t173",
  "\t3\t-161\t-156\t171\t173",
  "\t2\t-160\t-156\t171\t174",
  "\t1\t-158\t-156\t171\t175",
  "\t0\t-158\t-156\t171\t177",
  "\t-1\t-172\t-170\t171\t177",
  "\t-2\t-172\t-170\t173\t177",
  "\t-3\t-173\t-170\t-156\t-153\t174\t177",
  "\t-4\t-173\t-171\t-156\t-153",
  "\t-5\t-173\t-171\t-156\t-153",
  "\t-10\t-152\t-150",
  "\t-11\t-152\t-150",
  "\t-12\t-152\t-150",
  "Kosovo",
  "\t44\t19\t22",
  "\t43\t19\t22",
  "\t42\t19\t22",
  "\t41\t19\t22",
  "\t40\t19\t21",
  "Kuwait",
  "\t31\t46\t49",
  "\t30\t45\t49",
  "\t29\t45\t49",
  "\t28\t45\t49",
  "\t27\t46\t49",
  "Kyrgyzstan",
  "\t44\t72\t75",
  "\t43\t69\t81",
  "\t42\t69\t81",
  "\t41\t68\t81",
  "\t40\t68\t80",
  "\t39\t68\t78",
  "\t38\t68\t74",
  "Laos",
  "\t23\t100\t103",
  "\t22\t99\t104",
  "\t21\t99\t105",
  "\t20\t99\t105",
  "\t19\t99\t106",
  "\t18\t99\t107",
  "\t17\t99\t108",
  "\t16\t99\t108",
  "\t15\t103\t108",
  "\t14\t104\t108",
  "\t13\t104\t108",
  "\t12\t104\t107",
  "Latvia",
  "\t59\t23\t26",
  "\t58\t20\t28",
  "\t57\t19\t29",
  "\t56\t19\t29",
  "\t55\t19\t29",
  "\t54\t24\t28",
  "Lebanon",
  "\t35\t34\t37",
  "\t34\t34\t37",
  "\t33\t34\t37",
  "\t32\t34\t37",
  "Lesotho",
  "\t-27\t26\t30",
  "\t-28\t26\t30",
  "\t-29\t26\t30",
  "\t-30\t26\t30",
  "\t-31\t26\t29",
  "Liberia",
  "\t9\t-11\t-8",
  "\t8\t-12\t-7",
  "\t7\t-12\t-6",
  "\t6\t-12\t-6",
  "\t5\t-12\t-6",
  "\t4\t-11\t-6",
  "\t3\t-10\t-6",
  "Libya",
  "\t34\t10\t12",
  "\t33\t9\t16\t19\t25",
  "\t32\t9\t26",
  "\t31\t8\t26",
  "\t30\t8\t26",
  "\t29\t8\t25",
  "\t28\t8\t25",
  "\t27\t8\t25",
  "\t26\t8\t25",
  "\t25\t8\t25",
  "\t24\t8\t25",
  "\t23\t9\t25",
  "\t22\t10\t25",
  "\t21\t12\t25",
  "\t20\t17\t25",
  "\t19\t20\t25",
  "\t18\t21\t25",
  "Liechtenstein",
  "\t48\t8\t10",
  "\t47\t8\t10",
  "\t46\t8\t10",
  "Lithuania",
  "\t57\t20\t26",
  "\t56\t19\t27",
  "\t55\t19\t27",
  "\t54\t19\t27",
  "\t53\t21\t27",
  "\t52\t22\t25",
  "Luxembourg",
  "\t51\t4\t7",
  "\t50\t4\t7",
  "\t49\t4\t7",
  "\t48\t4\t7",
  "Macau",
  "\t23\t112\t114",
  "\t22\t112\t114",
  "\t21\t112\t114",
  "Macedonia",
  "\t43\t19\t23",
  "\t42\t19\t24",
  "\t41\t19\t24",
  "\t40\t19\t24",
  "\t39\t19\t22",
  "Madagascar",
  "\t-10\t48\t50",
  "\t-11\t47\t50",
  "\t-12\t46\t51",
  "\t-13\t46\t51",
  "\t-14\t44\t51",
  "\t-15\t43\t51",
  "\t-16\t42\t51",
  "\t-17\t42\t51",
  "\t-18\t42\t50",
  "\t-19\t42\t50",
  "\t-20\t42\t50",
  "\t-21\t42\t49",
  "\t-22\t42\t49",
  "\t-23\t42\t48",
  "\t-24\t42\t48",
  "\t-25\t42\t48",
  "\t-26\t43\t48",
  "Malawi",
  "\t-8\t31\t35",
  "\t-9\t31\t35",
  "\t-10\t31\t35",
  "\t-11\t31\t35",
  "\t-12\t31\t36",
  "\t-13\t31\t36",
  "\t-14\t31\t36",
  "\t-15\t31\t36",
  "\t-16\t33\t36",
  "\t-17\t33\t36",
  "\t-18\t34\t36",
  "Malaysia",
  "\t8\t115\t118",
  "\t7\t98\t103\t115\t119",
  "\t6\t98\t104\t114\t120",
  "\t5\t98\t104\t112\t120",
  "\t4\t99\t104\t111\t120",
  "\t3\t99\t105\t108\t119",
  "\t2\t99\t105\t108\t116",
  "\t1\t100\t105\t108\t116",
  "\t0\t101\t105\t108\t115",
  "\t-1\t109\t112",
  "Maldives",
  "\t8\t71\t73",
  "\t7\t71\t74",
  "\t6\t71\t74",
  "\t5\t71\t74",
  "\t4\t71\t74",
  "\t3\t71\t74",
  "\t2\t71\t74",
  "\t1\t71\t74",
  "\t0\t71\t74",
  "\t-1\t71\t74",
  "Mali",
  "\t25\t-7\t-2",
  "\t24\t-7\t0",
  "\t23\t-7\t1",
  "\t22\t-7\t2",
  "\t21\t-7\t3",
  "\t20\t-7\t5",
  "\t19\t-7\t5",
  "\t18\t-7\t5",
  "\t17\t-6\t5",
  "\t16\t-12\t5",
  "\t15\t-13\t5",
  "\t14\t-13\t4",
  "\t13\t-13\t1",
  "\t12\t-13\t-1",
  "\t11\t-12\t-3",
  "\t10\t-12\t-3",
  "\t9\t-9\t-4",
  "Malta",
  "\t37\t13\t15",
  "\t36\t13\t15",
  "\t35\t13\t15",
  "\t34\t13\t15",
  "Marshall Islands",
  "\t15\t167\t170",
  "\t14\t167\t170",
  "\t13\t167\t170",
  "\t12\t164\t167",
  "\t11\t164\t167\t169\t171",
  "\t10\t164\t167\t169\t171",
  "\t9\t166\t171",
  "\t8\t166\t172",
  "\t7\t166\t173",
  "\t6\t167\t173",
  "\t5\t167\t173",
  "\t4\t167\t170",
  "\t3\t167\t169",
  "Martinique",
  "\t15\t-62\t-59",
  "\t14\t-62\t-59",
  "\t13\t-62\t-59",
  "Mauritania",
  "\t28\t-9\t-7",
  "\t27\t-9\t-5",
  "\t26\t-13\t-3",
  "\t25\t-13\t-3",
  "\t24\t-14\t-3",
  "\t23\t-14\t-3",
  "\t22\t-18\t-5",
  "\t21\t-18\t-5",
  "\t20\t-18\t-4",
  "\t19\t-18\t-4",
  "\t18\t-17\t-4",
  "\t17\t-17\t-4",
  "\t16\t-17\t-4",
  "\t15\t-17\t-4",
  "\t14\t-17\t-4",
  "\t13\t-13\t-10",
  "Mauritius",
  "\t-9\t55\t57",
  "\t-10\t55\t57",
  "\t-11\t55\t57",
  "\t-18\t56\t58\t62\t64",
  "\t-19\t56\t58\t62\t64",
  "\t-20\t56\t58\t62\t64",
  "\t-21\t56\t58",
  "Mayotte",
  "\t-11\t44\t46",
  "\t-12\t44\t46",
  "\t-13\t44\t46",
  "Mexico",
  "\t33\t-118\t-112",
  "\t32\t-118\t-104",
  "\t31\t-118\t-103",
  "\t30\t-119\t-99",
  "\t29\t-119\t-98",
  "\t28\t-119\t-98",
  "\t27\t-119\t-96",
  "\t26\t-116\t-96",
  "\t25\t-115\t-96",
  "\t24\t-113\t-96",
  "\t23\t-113\t-96",
  "\t22\t-111\t-96\t-91\t-85",
  "\t21\t-111\t-95\t-91\t-85",
  "\t20\t-111\t-109\t-107\t-94\t-92\t-85",
  "\t19\t-115\t-109\t-106\t-85",
  "\t18\t-115\t-109\t-106\t-86",
  "\t17\t-115\t-109\t-105\t-86",
  "\t16\t-103\t-87",
  "\t15\t-101\t-89",
  "\t14\t-98\t-90",
  "\t13\t-93\t-91",
  "Micronesia",
  "\t10\t137\t139",
  "\t9\t137\t139\t148\t151\t153\t155",
  "\t8\t137\t139\t148\t155",
  "\t7\t148\t159",
  "\t6\t148\t154\t156\t159\t161\t164",
  "\t5\t148\t150\t152\t154\t156\t159\t161\t164",
  "\t4\t152\t154\t156\t158\t161\t164",
  "Midway Islands",
  "\t29\t-178\t-176",
  "\t28\t-178\t-176",
  "\t27\t-178\t-176",
  "Moldova",
  "\t49\t25\t29",
  "\t48\t25\t30",
  "\t47\t25\t31",
  "\t46\t26\t31",
  "\t45\t27\t31",
  "\t44\t27\t29",
  "Monaco",
  "\t44\t6\t8",
  "\t43\t6\t8",
  "\t42\t6\t8",
  "Mongolia",
  "\t53\t97\t100",
  "\t52\t96\t103",
  "\t51\t88\t108\t112\t117",
  "\t50\t86\t117",
  "\t49\t86\t119",
  "\t48\t86\t120",
  "\t47\t86\t120",
  "\t46\t88\t120",
  "\t45\t89\t120",
  "\t44\t89\t117",
  "\t43\t91\t115",
  "\t42\t94\t112",
  "\t41\t95\t111",
  "\t40\t102\t106",
  "Montenegro",
  "\t44\t17\t21",
  "\t43\t17\t21",
  "\t42\t17\t21",
  "\t41\t17\t21",
  "\t40\t18\t20",
  "Montserrat",
  "\t17\t-63\t-61",
  "\t16\t-63\t-61",
  "\t15\t-63\t-61",
  "Morocco",
  "\t36\t-7\t-1",
  "\t35\t-7\t0",
  "\t34\t-9\t0",
  "\t33\t-10\t0",
  "\t32\t-10\t0",
  "\t31\t-10\t0",
  "\t30\t-11\t-1",
  "\t29\t-13\t-2",
  "\t28\t-14\t-4",
  "\t27\t-15\t-7",
  "\t26\t-15\t-7",
  "\t25\t-16\t-8",
  "\t24\t-17\t-11",
  "\t23\t-17\t-11",
  "\t22\t-18\t-12",
  "\t21\t-18\t-13",
  "\t20\t-18\t-13",
  "Mozambique",
  "\t-9\t38\t41",
  "\t-10\t33\t41",
  "\t-11\t33\t41",
  "\t-12\t33\t41",
  "\t-13\t29\t41",
  "\t-14\t29\t41",
  "\t-15\t29\t41",
  "\t-16\t29\t41",
  "\t-17\t29\t41",
  "\t-18\t31\t39",
  "\t-19\t31\t37",
  "\t-20\t30\t36",
  "\t-21\t30\t36",
  "\t-22\t30\t36",
  "\t-23\t30\t36",
  "\t-24\t30\t36",
  "\t-25\t30\t36",
  "\t-26\t30\t34",
  "\t-27\t31\t33",
  "Myanmar",
  "\t29\t96\t99",
  "\t28\t94\t99",
  "\t27\t94\t99",
  "\t26\t93\t99",
  "\t25\t92\t99",
  "\t24\t92\t100",
  "\t23\t91\t100",
  "\t22\t91\t102",
  "\t21\t91\t102",
  "\t20\t91\t102",
  "\t19\t91\t101",
  "\t18\t91\t100",
  "\t17\t92\t99",
  "\t16\t93\t99",
  "\t15\t92\t99",
  "\t14\t92\t100",
  "\t13\t92\t94\t96\t100",
  "\t12\t96\t100",
  "\t11\t96\t100",
  "\t10\t96\t100",
  "\t9\t96\t100",
  "\t8\t97\t99",
  "Namibia",
  "\t-15\t12\t14",
  "\t-16\t10\t26",
  "\t-17\t10\t26",
  "\t-18\t10\t26",
  "\t-19\t10\t25",
  "\t-20\t11\t21",
  "\t-21\t12\t21",
  "\t-22\t12\t21",
  "\t-23\t13\t21",
  "\t-24\t13\t20",
  "\t-25\t13\t20",
  "\t-26\t13\t20",
  "\t-27\t13\t20",
  "\t-28\t14\t20",
  "\t-29\t14\t20",
  "Nauru",
  "\t1\t165\t167",
  "\t0\t165\t167",
  "\t-1\t165\t167",
  "Navassa Island",
  "\t19\t-76\t-74",
  "\t18\t-76\t-74",
  "\t17\t-76\t-74",
  "Nepal",
  "\t31\t79\t83",
  "\t30\t79\t85",
  "\t29\t79\t87",
  "\t28\t79\t89",
  "\t27\t79\t89",
  "\t26\t80\t89",
  "\t25\t83\t89",
  "Netherlands",
  "\t54\t3\t8",
  "\t53\t3\t8",
  "\t52\t2\t8",
  "\t51\t2\t8",
  "\t50\t2\t7",
  "\t49\t4\t7",
  "Netherlands Antilles",
  "\t18\t-64\t-61",
  "\t17\t-64\t-61",
  "\t16\t-64\t-61",
  "\t13\t-69\t-67",
  "\t12\t-69\t-67",
  "\t11\t-69\t-67",
  "New Caledonia",
  "\t-18\t158\t160\t162\t164",
  "\t-19\t158\t160\t162\t168",
  "\t-20\t158\t160\t162\t169",
  "\t-21\t162\t169",
  "\t-22\t163\t169",
  "\t-23\t165\t168",
  "New Zealand",
  "\t-7\t-173\t-171",
  "\t-8\t-173\t-170",
  "\t-9\t-173\t-170",
  "\t-10\t-172\t-170",
  "\t-28\t-178\t-176",
  "\t-29\t-178\t-176",
  "\t-30\t-178\t-176",
  "\t-33\t171\t174",
  "\t-34\t171\t175",
  "\t-35\t171\t176",
  "\t-36\t172\t179",
  "\t-37\t172\t179",
  "\t-38\t172\t179",
  "\t-39\t171\t179",
  "\t-40\t170\t179",
  "\t-41\t169\t177",
  "\t-42\t-177\t-175\t167\t177",
  "\t-43\t-177\t-175\t166\t175",
  "\t-44\t-177\t-175\t165\t174",
  "\t-45\t-177\t-175\t165\t172",
  "\t-46\t165\t172",
  "\t-47\t165\t171",
  "\t-48\t165\t169\t177\t179",
  "\t-49\t164\t167\t177\t179",
  "\t-50\t164\t167\t177\t179",
  "\t-51\t164\t170",
  "\t-52\t168\t170",
  "\t-53\t168\t170",
  "Nicaragua",
  "\t16\t-84\t-82",
  "\t15\t-87\t-81",
  "\t14\t-88\t-81",
  "\t13\t-88\t-81",
  "\t12\t-88\t-82",
  "\t11\t-88\t-82",
  "\t10\t-87\t-82",
  "\t9\t-85\t-82",
  "Niger",
  "\t24\t10\t14",
  "\t23\t8\t16",
  "\t22\t6\t16",
  "\t21\t5\t16",
  "\t20\t3\t16",
  "\t19\t3\t16",
  "\t18\t3\t16",
  "\t17\t2\t16",
  "\t16\t0\t16",
  "\t15\t-1\t16",
  "\t14\t-1\t15",
  "\t13\t-1\t14",
  "\t12\t-1\t14",
  "\t11\t0\t10",
  "\t10\t1\t4",
  "Nigeria",
  "\t14\t3\t15",
  "\t13\t2\t15",
  "\t12\t2\t15",
  "\t11\t2\t15",
  "\t10\t1\t15",
  "\t9\t1\t14",
  "\t8\t1\t14",
  "\t7\t1\t13",
  "\t6\t1\t13",
  "\t5\t1\t12",
  "\t4\t4\t10",
  "\t3\t4\t9",
  "Niue",
  "\t-17\t-170\t-168",
  "\t-18\t-170\t-168",
  "\t-19\t-170\t-168",
  "\t-20\t-170\t-168",
  "Norfolk Island",
  "\t-27\t166\t168",
  "\t-28\t166\t168",
  "\t-29\t166\t168",
  "\t-30\t166\t168",
  "North Korea",
  "\t44\t128\t130",
  "\t43\t127\t131",
  "\t42\t125\t131",
  "\t41\t123\t131",
  "\t40\t123\t131",
  "\t39\t123\t130",
  "\t38\t123\t129",
  "\t37\t123\t129",
  "\t36\t123\t127",
  "Northern Mariana Islands",
  "\t21\t143\t146",
  "\t20\t143\t146",
  "\t19\t143\t146",
  "\t18\t144\t146",
  "\t17\t144\t146",
  "\t16\t144\t146",
  "\t15\t144\t146",
  "\t14\t144\t146",
  "\t13\t144\t146",
  "Norway",
  "\t72\t22\t29",
  "\t71\t17\t32",
  "\t70\t14\t32",
  "\t69\t11\t32",
  "\t68\t11\t31",
  "\t67\t11\t26",
  "\t66\t10\t18",
  "\t65\t8\t17",
  "\t64\t6\t15",
  "\t63\t3\t15",
  "\t62\t3\t13",
  "\t61\t3\t13",
  "\t60\t3\t13",
  "\t59\t3\t13",
  "\t58\t4\t13",
  "\t57\t4\t12",
  "\t56\t5\t8",
  "Oman",
  "\t27\t55\t57",
  "\t26\t55\t57",
  "\t25\t54\t58",
  "\t24\t54\t60",
  "\t23\t54\t60",
  "\t22\t54\t60",
  "\t21\t54\t60",
  "\t20\t51\t60",
  "\t19\t50\t59",
  "\t18\t50\t58",
  "\t17\t50\t58",
  "\t16\t51\t57",
  "\t15\t51\t55",
  "Pakistan",
  "\t38\t73\t75",
  "\t37\t70\t77",
  "\t36\t70\t78",
  "\t35\t68\t78",
  "\t34\t68\t78",
  "\t33\t68\t78",
  "\t32\t65\t76",
  "\t31\t65\t76",
  "\t30\t59\t75",
  "\t29\t59\t75",
  "\t28\t59\t74",
  "\t27\t60\t73",
  "\t26\t60\t72",
  "\t25\t60\t72",
  "\t24\t60\t72",
  "\t23\t65\t72",
  "\t22\t66\t69",
  "Palau",
  "\t8\t133\t135",
  "\t7\t133\t135",
  "\t6\t131\t135",
  "\t5\t131\t135",
  "\t4\t130\t133",
  "\t3\t130\t132",
  "\t2\t130\t132",
  "\t1\t130\t132",
  "Palmyra Atoll",
  "\t6\t-163\t-161",
  "\t5\t-163\t-161",
  "\t4\t-163\t-161",
  "Panama",
  "\t10\t-83\t-76",
  "\t9\t-84\t-76",
  "\t8\t-84\t-76",
  "\t7\t-84\t-76",
  "\t6\t-82\t-76",
  "Papua New Guinea",
  "\t0\t141\t143\t145\t151",
  "\t-1\t139\t143\t145\t153",
  "\t-2\t139\t155",
  "\t-3\t139\t155",
  "\t-4\t139\t156",
  "\t-5\t139\t156",
  "\t-6\t139\t156",
  "\t-7\t139\t156",
  "\t-8\t139\t154",
  "\t-9\t139\t154",
  "\t-10\t139\t155",
  "\t-11\t146\t155",
  "\t-12\t152\t155",
  "Paracel Islands",
  "\t18\t110\t112",
  "\t17\t110\t113",
  "\t16\t110\t113",
  "\t15\t110\t113",
  "Paraguay",
  "\t-18\t-62\t-57",
  "\t-19\t-63\t-56",
  "\t-20\t-63\t-56",
  "\t-21\t-63\t-54",
  "\t-22\t-63\t-53",
  "\t-23\t-63\t-53",
  "\t-24\t-62\t-53",
  "\t-25\t-61\t-53",
  "\t-26\t-59\t-53",
  "\t-27\t-59\t-53",
  "\t-28\t-59\t-54",
  "Peru",
  "\t1\t-76\t-73",
  "\t0\t-76\t-72",
  "\t-1\t-78\t-69",
  "\t-2\t-81\t-69",
  "\t-3\t-82\t-68",
  "\t-4\t-82\t-68",
  "\t-5\t-82\t-68",
  "\t-6\t-82\t-71",
  "\t-7\t-82\t-71",
  "\t-8\t-80\t-69",
  "\t-9\t-80\t-68",
  "\t-10\t-79\t-68",
  "\t-11\t-79\t-67",
  "\t-12\t-78\t-67",
  "\t-13\t-78\t-67",
  "\t-14\t-77\t-67",
  "\t-15\t-77\t-67",
  "\t-16\t-76\t-67",
  "\t-17\t-75\t-67",
  "\t-18\t-73\t-68",
  "\t-19\t-71\t-68",
  "Philippines",
  "\t22\t120\t122",
  "\t21\t120\t123",
  "\t20\t120\t123",
  "\t19\t119\t123",
  "\t18\t119\t123",
  "\t17\t118\t123",
  "\t16\t118\t123",
  "\t15\t118\t125",
  "\t14\t118\t125",
  "\t13\t118\t126",
  "\t12\t118\t126",
  "\t11\t117\t127",
  "\t10\t116\t127",
  "\t9\t115\t127",
  "\t8\t115\t127",
  "\t7\t115\t127",
  "\t6\t115\t127",
  "\t5\t117\t127",
  "\t4\t118\t126",
  "\t3\t118\t120",
  "Pitcairn Islands",
  "\t-22\t-131\t-129",
  "\t-23\t-131\t-127\t-125\t-123",
  "\t-24\t-131\t-127\t-125\t-123",
  "\t-25\t-131\t-127\t-125\t-123",
  "\t-26\t-131\t-129",
  "Poland",
  "\t55\t13\t24",
  "\t54\t13\t24",
  "\t53\t13\t24",
  "\t52\t13\t24",
  "\t51\t13\t25",
  "\t50\t13\t25",
  "\t49\t13\t25",
  "\t48\t16\t24",
  "\t47\t21\t23",
  "Portugal",
  "\t43\t-9\t-7",
  "\t42\t-9\t-5",
  "\t41\t-9\t-5",
  "\t40\t-10\t-5",
  "\t39\t-10\t-5",
  "\t38\t-10\t-5",
  "\t37\t-10\t-5",
  "\t36\t-9\t-6",
  "\t35\t-8\t-6",
  "Portugal: Azores",
  "\t40\t-32\t-26",
  "\t39\t-32\t-26",
  "\t38\t-32\t-24",
  "\t37\t-29\t-24",
  "\t36\t-26\t-24",
  "\t35\t-26\t-24",
  "Portugal: Madeira",
  "\t34\t-17\t-15",
  "\t33\t-18\t-15",
  "\t32\t-18\t-15",
  "\t31\t-18\t-14",
  "\t30\t-17\t-14",
  "\t29\t-17\t-14",
  "Puerto Rico",
  "\t19\t-68\t-64",
  "\t18\t-68\t-64",
  "\t17\t-68\t-64",
  "\t16\t-68\t-64",
  "Qatar",
  "\t27\t50\t52",
  "\t26\t49\t52",
  "\t25\t49\t52",
  "\t24\t49\t52",
  "\t23\t49\t52",
  "Republic of the Congo",
  "\t4\t15\t19",
  "\t3\t12\t19",
  "\t2\t12\t19",
  "\t1\t12\t19",
  "\t0\t11\t19",
  "\t-1\t10\t18",
  "\t-2\t10\t18",
  "\t-3\t10\t17",
  "\t-4\t10\t17",
  "\t-5\t10\t16",
  "\t-6\t10\t13",
  "Reunion",
  "\t-19\t54\t56",
  "\t-20\t54\t56",
  "\t-21\t54\t56",
  "\t-22\t54\t56",
  "Romania",
  "\t49\t21\t28",
  "\t48\t20\t29",
  "\t47\t19\t29",
  "\t46\t19\t30",
  "\t45\t19\t30",
  "\t44\t19\t30",
  "\t43\t20\t30",
  "\t42\t21\t29",
  "Russia",
  "\t82\t49\t51\t53\t66\t88\t97",
  "\t81\t35\t37\t43\t66\t77\t81\t88\t100",
  "\t80\t35\t37\t43\t66\t75\t81\t88\t105",
  "\t79\t35\t37\t43\t66\t75\t81\t89\t108",
  "\t78\t49\t52\t57\t60\t66\t68\t75\t78\t88\t108\t155\t157",
  "\t77\t59\t70\t88\t114\t136\t143\t147\t153\t155\t157",
  "\t76\t54\t70\t80\t114\t134\t153\t155\t157",
  "\t75\t53\t70\t78\t117\t134\t153",
  "\t74\t52\t130\t134\t151",
  "\t73\t50\t61\t67\t130\t134\t151",
  "\t72\t-180\t-174\t50\t59\t65\t159\t177\t180",
  "\t71\t-180\t-174\t50\t61\t65\t163\t167\t172\t177\t180",
  "\t70\t-180\t-174\t27\t37\t47\t180",
  "\t69\t-180\t-175\t27\t180",
  "\t68\t-180\t-171\t27\t180",
  "\t67\t-180\t-168\t27\t180",
  "\t66\t-180\t-167\t28\t180",
  "\t65\t-180\t-167\t28\t180",
  "\t64\t-180\t-167\t28\t180",
  "\t63\t-176\t-171\t28\t180",
  "\t62\t27\t180",
  "\t61\t25\t180",
  "\t60\t25\t175",
  "\t59\t25\t173",
  "\t58\t26\t167\t169\t171",
  "\t57\t26\t143\t150\t165",
  "\t56\t19\t23\t26\t141\t154\t167",
  "\t55\t18\t23\t26\t144\t154\t169",
  "\t54\t18\t23\t27\t144\t154\t169",
  "\t53\t18\t23\t29\t144\t154\t163\t165\t169",
  "\t52\t30\t63\t71\t144\t154\t161",
  "\t51\t30\t63\t77\t121\t124\t145\t154\t159",
  "\t50\t33\t62\t78\t121\t125\t145\t153\t159",
  "\t49\t34\t50\t53\t62\t78\t99\t101\t120\t126\t145\t152\t157",
  "\t48\t36\t49\t83\t90\t94\t98\t106\t120\t126\t145\t151\t156",
  "\t47\t36\t50\t129\t145\t149\t155",
  "\t46\t35\t50\t129\t144\t146\t154",
  "\t45\t35\t50\t129\t153",
  "\t44\t35\t49\t129\t138\t141\t151",
  "\t43\t36\t49\t129\t137\t144\t148",
  "\t42\t38\t49\t129\t136\t144\t147",
  "\t41\t42\t49\t129\t135",
  "\t40\t45\t49",
  "Rwanda",
  "\t0\t28\t31",
  "\t-1\t27\t31",
  "\t-2\t27\t31",
  "\t-3\t27\t31",
  "Saint Helena",
  "\t-6\t-15\t-13",
  "\t-7\t-15\t-13",
  "\t-8\t-15\t-13",
  "\t-14\t-6\t-4",
  "\t-15\t-6\t-4",
  "\t-16\t-6\t-4",
  "\t-17\t-6\t-4",
  "\t-36\t-13\t-11",
  "\t-37\t-13\t-11",
  "\t-38\t-13\t-11",
  "\t-39\t-11\t-8",
  "\t-40\t-11\t-8",
  "\t-41\t-11\t-8",
  "Saint Kitts and Nevis",
  "\t18\t-63\t-61",
  "\t17\t-63\t-61",
  "\t16\t-63\t-61",
  "Saint Lucia",
  "\t15\t-62\t-59",
  "\t14\t-62\t-59",
  "\t13\t-62\t-59",
  "\t12\t-62\t-59",
  "Saint Pierre and Miquelon",
  "\t48\t-57\t-55",
  "\t47\t-57\t-55",
  "\t46\t-57\t-55",
  "\t45\t-57\t-55",
  "Saint Vincent and the Grenadines",
  "\t14\t-62\t-60",
  "\t13\t-62\t-60",
  "\t12\t-62\t-60",
  "\t11\t-62\t-60",
  "Samoa",
  "\t-12\t-173\t-170",
  "\t-13\t-173\t-170",
  "\t-14\t-173\t-170",
  "\t-15\t-172\t-170",
  "San Marino",
  "\t44\t11\t13",
  "\t43\t11\t13",
  "\t42\t11\t13",
  "Sao Tome and Principe",
  "\t2\t6\t8",
  "\t1\t5\t8",
  "\t0\t5\t8",
  "\t-1\t5\t7",
  "Saudi Arabia",
  "\t33\t37\t40",
  "\t32\t35\t43",
  "\t31\t35\t44",
  "\t30\t33\t48",
  "\t29\t33\t49",
  "\t28\t33\t50",
  "\t27\t33\t51",
  "\t26\t33\t51",
  "\t25\t34\t52",
  "\t24\t35\t53",
  "\t23\t36\t56",
  "\t22\t37\t56",
  "\t21\t37\t56",
  "\t20\t37\t56",
  "\t19\t38\t56",
  "\t18\t39\t55",
  "\t17\t40\t52",
  "\t16\t40\t48",
  "\t15\t40\t48",
  "Senegal",
  "\t17\t-17\t-12",
  "\t16\t-17\t-11",
  "\t15\t-18\t-10",
  "\t14\t-18\t-10",
  "\t13\t-18\t-10",
  "\t12\t-17\t-10",
  "\t11\t-17\t-10",
  "Serbia",
  "\t47\t18\t21",
  "\t46\t17\t22",
  "\t45\t17\t23",
  "\t44\t17\t23",
  "\t43\t17\t23",
  "\t42\t18\t23",
  "\t41\t19\t23",
  "Seychelles",
  "\t-2\t54\t56",
  "\t-3\t54\t56",
  "\t-4\t52\t56",
  "\t-5\t51\t56",
  "\t-6\t51\t57",
  "\t-7\t51\t53\t55\t57",
  "\t-8\t45\t48\t51\t53\t55\t57",
  "\t-9\t45\t48",
  "\t-10\t45\t48",
  "Sierra Leone",
  "\t10\t-14\t-9",
  "\t9\t-14\t-9",
  "\t8\t-14\t-9",
  "\t7\t-14\t-9",
  "\t6\t-13\t-9",
  "\t5\t-12\t-10",
  "Singapore",
  "\t2\t102\t105",
  "\t1\t102\t105",
  "\t0\t102\t105",
  "Sint Maarten",
  "\t19\t-64\t-62",
  "\t18\t-64\t-62",
  "\t17\t-64\t-62",
  "Slovakia",
  "\t50\t16\t23",
  "\t49\t15\t23",
  "\t48\t15\t23",
  "\t47\t15\t23",
  "\t46\t16\t19",
  "Slovenia",
  "\t47\t12\t17",
  "\t46\t12\t17",
  "\t45\t12\t17",
  "\t44\t12\t16",
  "Solomon Islands",
  "\t-5\t154\t158",
  "\t-6\t154\t161",
  "\t-7\t154\t163",
  "\t-8\t154\t163\t166\t168",
  "\t-9\t155\t168",
  "\t-10\t158\t168",
  "\t-11\t158\t169",
  "\t-12\t158\t161\t165\t169",
  "\t-13\t167\t169",
  "Somalia",
  "\t12\t47\t52",
  "\t11\t47\t52",
  "\t10\t47\t52",
  "\t9\t47\t52",
  "\t8\t45\t51",
  "\t7\t44\t51",
  "\t6\t44\t50",
  "\t5\t40\t50",
  "\t4\t40\t49",
  "\t3\t39\t49",
  "\t2\t39\t48",
  "\t1\t39\t47",
  "\t0\t39\t46",
  "\t-1\t39\t44",
  "\t-2\t40\t42",
  "Somalia",
  "\t12\t42\t44\t46\t49",
  "\t11\t41\t49",
  "\t10\t41\t49",
  "\t9\t41\t49",
  "\t8\t42\t49",
  "\t7\t43\t49",
  "\t6\t45\t48",
  "South Africa",
  "\t-21\t26\t32",
  "\t-22\t25\t32",
  "\t-23\t18\t21\t24\t32",
  "\t-24\t18\t32",
  "\t-25\t18\t33",
  "\t-26\t18\t33",
  "\t-27\t15\t33",
  "\t-28\t15\t33",
  "\t-29\t15\t33",
  "\t-30\t15\t32",
  "\t-31\t16\t31",
  "\t-32\t16\t31",
  "\t-33\t16\t30",
  "\t-34\t16\t28",
  "\t-35\t17\t26",
  "\t-45\t36\t38",
  "\t-46\t36\t38",
  "\t-47\t36\t38",
  "South Georgia and the South Sandwich Islands",
  "\t-52\t-43\t-36",
  "\t-53\t-43\t-33",
  "\t-54\t-43\t-33",
  "\t-55\t-39\t-33\t-29\t-26",
  "\t-56\t-35\t-33\t-29\t-25",
  "\t-57\t-29\t-25",
  "\t-58\t-28\t-25",
  "\t-59\t-28\t-25",
  "\t-60\t-28\t-25",
  "South Korea",
  "\t39\t125\t129",
  "\t38\t123\t131",
  "\t37\t123\t131",
  "\t36\t123\t131",
  "\t35\t124\t130",
  "\t34\t124\t130",
  "\t33\t124\t129",
  "\t32\t125\t127",
  "South Sudan",
  "\t13\t31\t34",
  "\t12\t31\t34",
  "\t11\t23\t34",
  "\t10\t23\t35",
  "\t9\t23\t35",
  "\t8\t23\t35",
  "\t7\t23\t35",
  "\t6\t23\t36",
  "\t5\t25\t36",
  "\t4\t25\t36",
  "\t3\t26\t35",
  "\t2\t29\t34",
  "Spain",
  "\t44\t-10\t0",
  "\t43\t-10\t4",
  "\t42\t-10\t4",
  "\t41\t-10\t5",
  "\t40\t-9\t5",
  "\t39\t-8\t5",
  "\t38\t-8\t5",
  "\t37\t-8\t2",
  "\t36\t-8\t1",
  "\t35\t-7\t0",
  "\t34\t-6\t-1",
  "Spain: Canary Islands",
  "\t30\t-14\t-12",
  "\t29\t-19\t-12",
  "\t28\t-19\t-12",
  "\t27\t-19\t-12",
  "\t26\t-19\t-14",
  "Spratly Islands",
  "\t12\t113\t115",
  "\t11\t112\t115",
  "\t10\t112\t116",
  "\t9\t110\t116",
  "\t8\t110\t116",
  "\t7\t110\t116",
  "\t6\t112\t114",
  "\t5\t112\t114",
  "Sri Lanka",
  "\t10\t78\t81",
  "\t9\t78\t82",
  "\t8\t78\t82",
  "\t7\t78\t82",
  "\t6\t78\t82",
  "\t5\t78\t82",
  "\t4\t79\t81",
  "Sudan",
  "\t23\t30\t32",
  "\t22\t23\t38",
  "\t21\t23\t38",
  "\t20\t22\t38",
  "\t19\t22\t39",
  "\t18\t22\t39",
  "\t17\t22\t39",
  "\t16\t21\t39",
  "\t15\t21\t38",
  "\t14\t20\t37",
  "\t13\t20\t37",
  "\t12\t20\t37",
  "\t11\t20\t37",
  "\t10\t21\t36",
  "\t9\t21\t35",
  "\t8\t22\t35",
  "\t7\t22\t25",
  "Suriname",
  "\t7\t-56\t-54",
  "\t6\t-58\t-52",
  "\t5\t-59\t-52",
  "\t4\t-59\t-52",
  "\t3\t-59\t-52",
  "\t2\t-59\t-52",
  "\t1\t-58\t-53",
  "\t0\t-57\t-54",
  "Svalbard",
  "\t81\t15\t28\t30\t34",
  "\t80\t9\t34",
  "\t79\t9\t34",
  "\t78\t9\t31",
  "\t77\t9\t31",
  "\t76\t12\t26",
  "\t75\t14\t20\t23\t26",
  "\t74\t17\t20",
  "\t73\t17\t20",
  "Swaziland",
  "\t-24\t30\t33",
  "\t-25\t29\t33",
  "\t-26\t29\t33",
  "\t-27\t29\t33",
  "\t-28\t29\t32",
  "Sweden",
  "\t70\t19\t21",
  "\t69\t16\t24",
  "\t68\t15\t24",
  "\t67\t13\t25",
  "\t66\t13\t25",
  "\t65\t11\t25",
  "\t64\t10\t25",
  "\t63\t10\t22",
  "\t62\t10\t21",
  "\t61\t11\t19",
  "\t60\t10\t20",
  "\t59\t10\t20",
  "\t58\t10\t20",
  "\t57\t10\t20",
  "\t56\t10\t20",
  "\t55\t11\t19",
  "\t54\t11\t15",
  "Switzerland",
  "\t48\t5\t10",
  "\t47\t4\t11",
  "\t46\t4\t11",
  "\t45\t4\t11",
  "\t44\t5\t10",
  "Syria",
  "\t38\t39\t43",
  "\t37\t35\t43",
  "\t36\t34\t43",
  "\t35\t34\t43",
  "\t34\t34\t42",
  "\t33\t34\t42",
  "\t32\t34\t40",
  "\t31\t34\t39",
  "Taiwan",
  "\t26\t120\t123",
  "\t25\t117\t123",
  "\t24\t117\t123",
  "\t23\t117\t122",
  "\t22\t118\t122",
  "\t21\t119\t122",
  "\t20\t119\t121",
  "Tajikistan",
  "\t42\t69\t71",
  "\t41\t67\t71",
  "\t40\t66\t74",
  "\t39\t66\t75",
  "\t38\t66\t76",
  "\t37\t66\t76",
  "\t36\t66\t76",
  "\t35\t66\t73",
  "Tanzania",
  "\t1\t29\t31",
  "\t0\t29\t36",
  "\t-1\t29\t38",
  "\t-2\t29\t39",
  "\t-3\t28\t40",
  "\t-4\t28\t40",
  "\t-5\t28\t40",
  "\t-6\t28\t40",
  "\t-7\t28\t40",
  "\t-8\t29\t40",
  "\t-9\t29\t41",
  "\t-10\t30\t41",
  "\t-11\t33\t41",
  "\t-12\t33\t40",
  "Thailand",
  "\t21\t98\t101",
  "\t20\t96\t102",
  "\t19\t96\t105",
  "\t18\t96\t105",
  "\t17\t96\t106",
  "\t16\t96\t106",
  "\t15\t97\t106",
  "\t14\t97\t106",
  "\t13\t97\t106",
  "\t12\t97\t103",
  "\t11\t97\t103",
  "\t10\t96\t103",
  "\t9\t96\t101",
  "\t8\t96\t101",
  "\t7\t97\t103",
  "\t6\t97\t103",
  "\t5\t98\t103",
  "\t4\t99\t102",
  "Togo",
  "\t12\t-1\t1",
  "\t11\t-1\t2",
  "\t10\t-1\t2",
  "\t9\t-1\t2",
  "\t8\t-1\t2",
  "\t7\t-1\t2",
  "\t6\t-1\t2",
  "\t5\t-1\t2",
  "Tokelau",
  "\t-8\t-172\t-170",
  "\t-9\t-172\t-170",
  "\t-10\t-172\t-170",
  "Tonga",
  "\t-14\t-176\t-172",
  "\t-15\t-176\t-172",
  "\t-16\t-176\t-172",
  "\t-17\t-175\t-172",
  "\t-18\t-176\t-172",
  "\t-19\t-176\t-172",
  "\t-20\t-176\t-173",
  "\t-21\t-177\t-173",
  "\t-22\t-177\t-173",
  "\t-23\t-177\t-175",
  "Trinidad and Tobago",
  "\t12\t-61\t-59",
  "\t11\t-62\t-59",
  "\t10\t-62\t-59",
  "\t9\t-62\t-59",
  "Tromelin Island",
  "\t-14\t53\t55",
  "\t-15\t53\t55",
  "\t-16\t53\t55",
  "Tunisia",
  "\t38\t7\t12",
  "\t37\t7\t12",
  "\t36\t7\t12",
  "\t35\t6\t12",
  "\t34\t6\t12",
  "\t33\t6\t12",
  "\t32\t6\t12",
  "\t31\t7\t12",
  "\t30\t8\t11",
  "\t29\t8\t11",
  "Turkey",
  "\t43\t25\t28\t32\t36",
  "\t42\t25\t44",
  "\t41\t24\t45",
  "\t40\t24\t45",
  "\t39\t24\t45",
  "\t38\t24\t45",
  "\t37\t25\t45",
  "\t36\t26\t45",
  "\t35\t26\t41\t43\t45",
  "\t34\t34\t37",
  "Turkmenistan",
  "\t43\t51\t61",
  "\t42\t51\t62",
  "\t41\t51\t63",
  "\t40\t51\t64",
  "\t39\t51\t67",
  "\t38\t51\t67",
  "\t37\t52\t67",
  "\t36\t52\t67",
  "\t35\t59\t65",
  "\t34\t60\t65",
  "Turks and Caicos Islands",
  "\t22\t-73\t-70",
  "\t21\t-73\t-70",
  "\t20\t-73\t-70",
  "Tuvalu",
  "\t-4\t175\t177",
  "\t-5\t175\t178",
  "\t-6\t175\t179",
  "\t-7\t175\t180",
  "\t-8\t176\t180",
  "\t-9\t177\t180",
  "\t-10\t178\t180",
  "Uganda",
  "\t5\t32\t35",
  "\t4\t29\t35",
  "\t3\t29\t35",
  "\t2\t29\t36",
  "\t1\t28\t36",
  "\t0\t28\t36",
  "\t-1\t28\t35",
  "\t-2\t28\t34",
  "Ukraine",
  "\t53\t29\t35",
  "\t52\t22\t36",
  "\t51\t22\t39",
  "\t50\t21\t41",
  "\t49\t21\t41",
  "\t48\t21\t41",
  "\t47\t21\t41",
  "\t46\t21\t40",
  "\t45\t27\t38",
  "\t44\t27\t37",
  "\t43\t32\t36",
  "United Arab Emirates",
  "\t27\t55\t57",
  "\t26\t53\t57",
  "\t25\t50\t57",
  "\t24\t50\t57",
  "\t23\t50\t57",
  "\t22\t50\t56",
  "\t21\t51\t56",
  "United Kingdom",
  "\t61\t-3\t1",
  "\t60\t-4\t1",
  "\t59\t-8\t1",
  "\t58\t-14\t-12\t-9\t0",
  "\t57\t-14\t-12\t-9\t0",
  "\t56\t-14\t-12\t-9\t0",
  "\t55\t-9\t1",
  "\t54\t-9\t1",
  "\t53\t-9\t2",
  "\t52\t-6\t2",
  "\t51\t-6\t2",
  "\t50\t-7\t2",
  "\t49\t-7\t1",
  "\t48\t-7\t-4",
  "Uruguay",
  "\t-29\t-58\t-54",
  "\t-30\t-59\t-52",
  "\t-31\t-59\t-52",
  "\t-32\t-59\t-52",
  "\t-33\t-59\t-52",
  "\t-34\t-59\t-52",
  "\t-35\t-59\t-52",
  "USA",
  "\t72\t-158\t-153",
  "\t71\t-163\t-141",
  "\t70\t-164\t-140",
  "\t69\t-167\t-140",
  "\t68\t-167\t-140",
  "\t67\t-167\t-140",
  "\t66\t-169\t-140",
  "\t65\t-169\t-140",
  "\t64\t-172\t-140",
  "\t63\t-172\t-140",
  "\t62\t-172\t-140",
  "\t61\t-174\t-138",
  "\t60\t-174\t-171\t-168\t-133",
  "\t59\t-174\t-171\t-168\t-132",
  "\t58\t-171\t-131",
  "\t57\t-171\t-168\t-163\t-150\t-138\t-129",
  "\t56\t-171\t-168\t-164\t-151\t-137\t-128",
  "\t55\t-170\t-152\t-136\t-128",
  "\t54\t-170\t-154\t-135\t-128\t171\t173",
  "\t53\t-177\t-158\t-134\t-129\t171\t180",
  "\t52\t-180\t-165\t171\t180",
  "\t51\t-180\t-167\t171\t180",
  "\t50\t-180\t-174\t-96\t-93\t176\t180",
  "\t49\t-125\t-86",
  "\t48\t-125\t-84\t-70\t-66",
  "\t47\t-125\t-82\t-71\t-66",
  "\t46\t-125\t-81\t-75\t-66",
  "\t45\t-125\t-81\t-77\t-65",
  "\t44\t-125\t-65",
  "\t43\t-125\t-65",
  "\t42\t-125\t-68",
  "\t41\t-125\t-68",
  "\t40\t-125\t-68",
  "\t39\t-125\t-71",
  "\t38\t-124\t-73",
  "\t37\t-124\t-73",
  "\t36\t-124\t-74",
  "\t35\t-123\t-74",
  "\t34\t-122\t-74",
  "\t33\t-121\t-75",
  "\t32\t-121\t-76",
  "\t31\t-119\t-78",
  "\t30\t-114\t-79",
  "\t29\t-106\t-79",
  "\t28\t-105\t-79",
  "\t27\t-174\t-172\t-104\t-94\t-90\t-88\t-83\t-79",
  "\t26\t-174\t-166\t-100\t-95\t-83\t-79",
  "\t25\t-174\t-166\t-100\t-96\t-83\t-79",
  "\t24\t-172\t-160\t-98\t-96\t-83\t-79",
  "\t23\t-165\t-158\t-83\t-79",
  "\t22\t-165\t-155",
  "\t21\t-161\t-154",
  "\t20\t-161\t-153",
  "\t19\t-158\t-153",
  "\t18\t-157\t-153",
  "\t17\t-156\t-154",
  "USA: Alabama",
  "\t36\t-89\t-86",
  "\t35\t-89\t-84",
  "\t34\t-89\t-84",
  "\t33\t-89\t-83",
  "\t32\t-89\t-83",
  "\t31\t-89\t-83",
  "\t30\t-89\t-83",
  "\t29\t-89\t-86",
  "USA: Alaska",
  "\t72\t-158\t-153",
  "\t71\t-163\t-141",
  "\t70\t-164\t-140",
  "\t69\t-167\t-140",
  "\t68\t-167\t-140",
  "\t67\t-167\t-140",
  "\t66\t-169\t-140",
  "\t65\t-169\t-140",
  "\t64\t-172\t-140",
  "\t63\t-172\t-140",
  "\t62\t-172\t-140",
  "\t61\t-174\t-138",
  "\t60\t-174\t-171\t-168\t-133",
  "\t59\t-174\t-171\t-168\t-132",
  "\t58\t-171\t-131",
  "\t57\t-171\t-168\t-163\t-150\t-138\t-129",
  "\t56\t-171\t-168\t-164\t-151\t-137\t-128",
  "\t55\t-170\t-152\t-136\t-128",
  "\t54\t-170\t-154\t-135\t-128\t171\t173",
  "\t53\t-177\t-158\t-134\t-129\t171\t180",
  "\t52\t-180\t-165\t171\t180",
  "\t51\t-180\t-167\t171\t180",
  "\t50\t-180\t-174\t176\t180",
  "USA: Alaska, Aleutian Islands",
  "\t60\t-154\t-149\t-147\t-145",
  "\t59\t-162\t-159\t-154\t-149\t-147\t-145",
  "\t58\t-171\t-169\t-162\t-159\t-155\t-149\t-147\t-145",
  "\t57\t-171\t-168\t-162\t-150",
  "\t56\t-171\t-168\t-164\t-151",
  "\t55\t-170\t-152",
  "\t54\t-170\t-154",
  "\t53\t-177\t-158",
  "\t52\t-180\t-165",
  "\t51\t-180\t-167",
  "\t50\t-180\t-174",
  "USA: Arizona",
  "\t38\t-115\t-108",
  "\t37\t-115\t-108",
  "\t36\t-115\t-108",
  "\t35\t-115\t-108",
  "\t34\t-115\t-108",
  "\t33\t-115\t-108",
  "\t32\t-115\t-108",
  "\t31\t-115\t-108",
  "\t30\t-114\t-108",
  "USA: Arkansas",
  "\t37\t-95\t-88",
  "\t36\t-95\t-88",
  "\t35\t-95\t-88",
  "\t34\t-95\t-88",
  "\t33\t-95\t-89",
  "\t32\t-95\t-90",
  "USA: California",
  "\t43\t-125\t-119",
  "\t42\t-125\t-119",
  "\t41\t-125\t-119",
  "\t40\t-125\t-119",
  "\t39\t-125\t-117",
  "\t38\t-124\t-116",
  "\t37\t-124\t-115",
  "\t36\t-124\t-113",
  "\t35\t-123\t-113",
  "\t34\t-122\t-113",
  "\t33\t-121\t-113",
  "\t32\t-121\t-113",
  "\t31\t-119\t-113",
  "USA: Colorado",
  "\t42\t-110\t-101",
  "\t41\t-110\t-101",
  "\t40\t-110\t-101",
  "\t39\t-110\t-101",
  "\t38\t-110\t-101",
  "\t37\t-110\t-101",
  "\t36\t-110\t-101",
  "USA: Connecticut",
  "\t43\t-74\t-70",
  "\t42\t-74\t-70",
  "\t41\t-74\t-70",
  "\t40\t-74\t-70",
  "USA: Delaware",
  "\t40\t-76\t-74",
  "\t39\t-76\t-74",
  "\t38\t-76\t-74",
  "\t37\t-76\t-74",
  "USA: District of Columbia",
  "\t39\t-78\t-75",
  "\t38\t-78\t-75",
  "\t37\t-78\t-75",
  "USA: Florida",
  "\t32\t-88\t-84",
  "\t31\t-88\t-80",
  "\t30\t-88\t-79",
  "\t29\t-88\t-79",
  "\t28\t-86\t-79",
  "\t27\t-83\t-79",
  "\t26\t-83\t-79",
  "\t25\t-83\t-79",
  "\t24\t-83\t-79",
  "\t23\t-83\t-79",
  "USA: Georgia",
  "\t36\t-84\t-82",
  "\t35\t-86\t-81",
  "\t34\t-86\t-80",
  "\t33\t-86\t-79",
  "\t32\t-86\t-79",
  "\t31\t-86\t-79",
  "\t30\t-86\t-79",
  "\t29\t-86\t-80",
  "USA: Hawaii",
  "\t27\t-174\t-172",
  "\t26\t-174\t-166",
  "\t25\t-174\t-166",
  "\t24\t-172\t-160",
  "\t23\t-165\t-158",
  "\t22\t-165\t-155",
  "\t21\t-161\t-154",
  "\t20\t-161\t-153",
  "\t19\t-158\t-153",
  "\t18\t-157\t-153",
  "\t17\t-156\t-154",
  "USA: Idaho",
  "\t49\t-118\t-115",
  "\t48\t-118\t-114",
  "\t47\t-118\t-113",
  "\t46\t-118\t-112",
  "\t45\t-118\t-110",
  "\t44\t-118\t-110",
  "\t43\t-118\t-110",
  "\t42\t-118\t-110",
  "\t41\t-118\t-110",
  "USA: Illinois",
  "\t43\t-91\t-86",
  "\t42\t-92\t-86",
  "\t41\t-92\t-86",
  "\t40\t-92\t-86",
  "\t39\t-92\t-86",
  "\t38\t-92\t-86",
  "\t37\t-91\t-86",
  "\t36\t-91\t-87",
  "\t35\t-90\t-88",
  "USA: Indiana",
  "\t42\t-88\t-83",
  "\t41\t-88\t-83",
  "\t40\t-88\t-83",
  "\t39\t-89\t-83",
  "\t38\t-89\t-83",
  "\t37\t-89\t-83",
  "\t36\t-89\t-85",
  "USA: Iowa",
  "\t44\t-97\t-90",
  "\t43\t-97\t-89",
  "\t42\t-97\t-89",
  "\t41\t-97\t-89",
  "\t40\t-97\t-89",
  "\t39\t-96\t-89",
  "USA: Kansas",
  "\t41\t-103\t-94",
  "\t40\t-103\t-93",
  "\t39\t-103\t-93",
  "\t38\t-103\t-93",
  "\t37\t-103\t-93",
  "\t36\t-103\t-93",
  "USA: Kentucky",
  "\t40\t-85\t-83",
  "\t39\t-87\t-81",
  "\t38\t-90\t-80",
  "\t37\t-90\t-80",
  "\t36\t-90\t-80",
  "\t35\t-90\t-81",
  "USA: Louisiana",
  "\t34\t-95\t-90",
  "\t33\t-95\t-89",
  "\t32\t-95\t-88",
  "\t31\t-95\t-87",
  "\t30\t-95\t-87",
  "\t29\t-94\t-87",
  "\t28\t-94\t-87",
  "\t27\t-90\t-88",
  "USA: Maine",
  "\t48\t-70\t-66",
  "\t47\t-71\t-66",
  "\t46\t-72\t-66",
  "\t45\t-72\t-65",
  "\t44\t-72\t-65",
  "\t43\t-72\t-65",
  "\t42\t-71\t-68",
  "USA: Maryland",
  "\t40\t-80\t-74",
  "\t39\t-80\t-74",
  "\t38\t-80\t-74",
  "\t37\t-78\t-74",
  "\t36\t-77\t-74",
  "USA: Massachusetts",
  "\t43\t-74\t-69",
  "\t42\t-74\t-68",
  "\t41\t-74\t-68",
  "\t40\t-72\t-68",
  "USA: Michigan",
  "\t49\t-90\t-86",
  "\t48\t-91\t-84",
  "\t47\t-91\t-82",
  "\t46\t-91\t-81",
  "\t45\t-91\t-81",
  "\t44\t-89\t-81",
  "\t43\t-88\t-81",
  "\t42\t-88\t-81",
  "\t41\t-88\t-81",
  "\t40\t-88\t-82",
  "USA: Minnesota",
  "\t50\t-96\t-93",
  "\t49\t-98\t-88",
  "\t48\t-98\t-88",
  "\t47\t-98\t-88",
  "\t46\t-98\t-88",
  "\t45\t-97\t-90",
  "\t44\t-97\t-90",
  "\t43\t-97\t-90",
  "\t42\t-97\t-90",
  "USA: Mississippi",
  "\t36\t-91\t-88",
  "\t35\t-92\t-87",
  "\t34\t-92\t-87",
  "\t33\t-92\t-87",
  "\t32\t-92\t-87",
  "\t31\t-92\t-87",
  "\t30\t-92\t-87",
  "\t29\t-90\t-87",
  "USA: Missouri",
  "\t41\t-96\t-90",
  "\t40\t-96\t-89",
  "\t39\t-96\t-89",
  "\t38\t-96\t-88",
  "\t37\t-95\t-88",
  "\t36\t-95\t-88",
  "\t35\t-95\t-88",
  "\t34\t-91\t-88",
  "USA: Montana",
  "\t49\t-117\t-103",
  "\t48\t-117\t-103",
  "\t47\t-117\t-103",
  "\t46\t-117\t-103",
  "\t45\t-116\t-103",
  "\t44\t-115\t-103",
  "\t43\t-114\t-110",
  "USA: Nebraska",
  "\t44\t-105\t-97",
  "\t43\t-105\t-95",
  "\t42\t-105\t-94",
  "\t41\t-105\t-94",
  "\t40\t-105\t-94",
  "\t39\t-103\t-94",
  "USA: Nevada",
  "\t43\t-121\t-113",
  "\t42\t-121\t-113",
  "\t41\t-121\t-113",
  "\t40\t-121\t-113",
  "\t39\t-121\t-113",
  "\t38\t-121\t-113",
  "\t37\t-120\t-113",
  "\t36\t-119\t-113",
  "\t35\t-118\t-113",
  "\t34\t-116\t-113",
  "USA: New Hampshire",
  "\t46\t-72\t-70",
  "\t45\t-73\t-69",
  "\t44\t-73\t-69",
  "\t43\t-73\t-69",
  "\t42\t-73\t-69",
  "\t41\t-73\t-69",
  "USA: New Jersey",
  "\t42\t-76\t-73",
  "\t41\t-76\t-72",
  "\t40\t-76\t-72",
  "\t39\t-76\t-72",
  "\t38\t-76\t-73",
  "\t37\t-75\t-73",
  "USA: New Mexico",
  "\t38\t-110\t-102",
  "\t37\t-110\t-102",
  "\t36\t-110\t-102",
  "\t35\t-110\t-102",
  "\t34\t-110\t-102",
  "\t33\t-110\t-102",
  "\t32\t-110\t-102",
  "\t31\t-110\t-102",
  "\t30\t-110\t-105",
  "USA: New York",
  "\t46\t-75\t-72",
  "\t45\t-77\t-72",
  "\t44\t-80\t-72",
  "\t43\t-80\t-72",
  "\t42\t-80\t-70",
  "\t41\t-80\t-70",
  "\t40\t-79\t-70",
  "\t39\t-75\t-71",
  "USA: North Carolina",
  "\t37\t-83\t-74",
  "\t36\t-85\t-74",
  "\t35\t-85\t-74",
  "\t34\t-85\t-74",
  "\t33\t-85\t-75",
  "\t32\t-79\t-76",
  "USA: North Dakota",
  "\t49\t-105\t-96",
  "\t48\t-105\t-95",
  "\t47\t-105\t-95",
  "\t46\t-105\t-95",
  "\t45\t-105\t-95",
  "\t44\t-105\t-95",
  "USA: Ohio",
  "\t43\t-82\t-79",
  "\t42\t-85\t-79",
  "\t41\t-85\t-79",
  "\t40\t-85\t-79",
  "\t39\t-85\t-79",
  "\t38\t-85\t-79",
  "\t37\t-85\t-80",
  "USA: Oklahoma",
  "\t38\t-104\t-93",
  "\t37\t-104\t-93",
  "\t36\t-104\t-93",
  "\t35\t-104\t-93",
  "\t34\t-101\t-93",
  "\t33\t-101\t-93",
  "\t32\t-98\t-93",
  "USA: Oregon",
  "\t47\t-124\t-115",
  "\t46\t-125\t-115",
  "\t45\t-125\t-115",
  "\t44\t-125\t-115",
  "\t43\t-125\t-115",
  "\t42\t-125\t-115",
  "\t41\t-125\t-116",
  "USA: Pennsylvania",
  "\t43\t-81\t-77",
  "\t42\t-81\t-73",
  "\t41\t-81\t-73",
  "\t40\t-81\t-73",
  "\t39\t-81\t-73",
  "\t38\t-81\t-74",
  "USA: Rhode Island",
  "\t43\t-72\t-70",
  "\t42\t-72\t-70",
  "\t41\t-72\t-70",
  "\t40\t-72\t-70",
  "USA: South Carolina",
  "\t36\t-84\t-79",
  "\t35\t-84\t-77",
  "\t34\t-84\t-77",
  "\t33\t-84\t-77",
  "\t32\t-83\t-77",
  "\t31\t-82\t-78",
  "USA: South Dakota",
  "\t46\t-105\t-95",
  "\t45\t-105\t-95",
  "\t44\t-105\t-95",
  "\t43\t-105\t-95",
  "\t42\t-105\t-95",
  "\t41\t-99\t-95",
  "USA: Tennessee",
  "\t37\t-90\t-80",
  "\t36\t-91\t-80",
  "\t35\t-91\t-80",
  "\t34\t-91\t-81",
  "\t33\t-90\t-83",
  "USA: Texas",
  "\t37\t-104\t-99",
  "\t36\t-104\t-99",
  "\t35\t-104\t-97",
  "\t34\t-104\t-93",
  "\t33\t-107\t-93",
  "\t32\t-107\t-92",
  "\t31\t-107\t-92",
  "\t30\t-107\t-92",
  "\t29\t-106\t-92",
  "\t28\t-105\t-92",
  "\t27\t-104\t-94",
  "\t26\t-100\t-95",
  "\t25\t-100\t-96",
  "\t24\t-98\t-96",
  "USA: Utah",
  "\t43\t-115\t-110",
  "\t42\t-115\t-108",
  "\t41\t-115\t-108",
  "\t40\t-115\t-108",
  "\t39\t-115\t-108",
  "\t38\t-115\t-108",
  "\t37\t-115\t-108",
  "\t36\t-115\t-108",
  "USA: Vermont",
  "\t46\t-74\t-70",
  "\t45\t-74\t-70",
  "\t44\t-74\t-70",
  "\t43\t-74\t-70",
  "\t42\t-74\t-71",
  "\t41\t-74\t-71",
  "USA: Virginia",
  "\t40\t-79\t-76",
  "\t39\t-81\t-74",
  "\t38\t-83\t-74",
  "\t37\t-84\t-74",
  "\t36\t-84\t-74",
  "\t35\t-84\t-74",
  "USA: Washington",
  "\t49\t-125\t-116",
  "\t48\t-125\t-116",
  "\t47\t-125\t-115",
  "\t46\t-125\t-115",
  "\t45\t-125\t-115",
  "\t44\t-123\t-118",
  "USA: West Virginia",
  "\t41\t-81\t-79",
  "\t40\t-83\t-76",
  "\t39\t-83\t-76",
  "\t38\t-83\t-76",
  "\t37\t-83\t-77",
  "\t36\t-83\t-79",
  "USA: Wisconsin",
  "\t48\t-92\t-88",
  "\t47\t-93\t-87",
  "\t46\t-93\t-85",
  "\t45\t-93\t-85",
  "\t44\t-93\t-85",
  "\t43\t-93\t-85",
  "\t42\t-92\t-86",
  "\t41\t-92\t-86",
  "USA: Wyoming",
  "\t46\t-112\t-103",
  "\t45\t-112\t-103",
  "\t44\t-112\t-103",
  "\t43\t-112\t-103",
  "\t42\t-112\t-103",
  "\t41\t-112\t-103",
  "\t40\t-112\t-103",
  "Uzbekistan",
  "\t46\t55\t60",
  "\t45\t54\t62",
  "\t44\t54\t66",
  "\t43\t54\t67\t69\t72",
  "\t42\t54\t73",
  "\t41\t54\t74",
  "\t40\t54\t74",
  "\t39\t60\t74",
  "\t38\t61\t72",
  "\t37\t63\t69",
  "\t36\t65\t69",
  "Vanuatu",
  "\t-12\t165\t168",
  "\t-13\t165\t169",
  "\t-14\t165\t169",
  "\t-15\t165\t169",
  "\t-16\t165\t169",
  "\t-17\t166\t170",
  "\t-18\t167\t170",
  "\t-19\t167\t170",
  "\t-20\t168\t170",
  "\t-21\t168\t170",
  "Venezuela",
  "\t13\t-71\t-66",
  "\t12\t-73\t-62",
  "\t11\t-73\t-60",
  "\t10\t-74\t-59",
  "\t9\t-74\t-58",
  "\t8\t-74\t-58",
  "\t7\t-73\t-58",
  "\t6\t-73\t-59",
  "\t5\t-72\t-59",
  "\t4\t-68\t-59",
  "\t3\t-68\t-59",
  "\t2\t-68\t-61",
  "\t1\t-68\t-62",
  "\t0\t-68\t-62",
  "\t-1\t-67\t-64",
  "Viet Nam",
  "\t24\t103\t106",
  "\t23\t101\t107",
  "\t22\t101\t108",
  "\t21\t101\t108",
  "\t20\t101\t108",
  "\t19\t102\t108",
  "\t18\t102\t108",
  "\t17\t103\t109",
  "\t16\t104\t109",
  "\t15\t105\t110",
  "\t14\t106\t110",
  "\t13\t105\t110",
  "\t12\t104\t110",
  "\t11\t102\t110",
  "\t10\t102\t110",
  "\t9\t102\t109",
  "\t8\t103\t107",
  "\t7\t103\t107",
  "Virgin Islands",
  "\t19\t-66\t-63",
  "\t18\t-66\t-63",
  "\t17\t-66\t-63",
  "\t16\t-65\t-63",
  "Wake Island",
  "\t20\t165\t167",
  "\t19\t165\t167",
  "\t18\t165\t167",
  "Wallis and Futuna",
  "\t-12\t-177\t-175",
  "\t-13\t-179\t-175",
  "\t-14\t-179\t-175",
  "\t-15\t-179\t-176",
  "West Bank",
  "\t33\t33\t36",
  "\t32\t33\t36",
  "\t31\t33\t36",
  "\t30\t33\t36",
  "Western Sahara",
  "\t28\t-14\t-7",
  "\t27\t-15\t-7",
  "\t26\t-15\t-7",
  "\t25\t-16\t-7",
  "\t24\t-17\t-7",
  "\t23\t-17\t-11",
  "\t22\t-18\t-11",
  "\t21\t-18\t-12",
  "\t20\t-18\t-12",
  "\t19\t-18\t-16",
  "Yemen",
  "\t19\t47\t53",
  "\t18\t42\t53",
  "\t17\t41\t54",
  "\t16\t41\t54",
  "\t15\t41\t54",
  "\t14\t41\t53",
  "\t13\t41\t55",
  "\t12\t41\t49\t51\t55",
  "\t11\t42\t46\t51\t55",
  "Zambia",
  "\t-7\t27\t32",
  "\t-8\t27\t34",
  "\t-9\t22\t25\t27\t34",
  "\t-10\t22\t34",
  "\t-11\t22\t34",
  "\t-12\t20\t34",
  "\t-13\t20\t34",
  "\t-14\t20\t34",
  "\t-15\t20\t34",
  "\t-16\t20\t31",
  "\t-17\t20\t29",
  "\t-18\t21\t28",
  "\t-19\t24\t27",
  "Zimbabwe",
  "\t-14\t27\t32",
  "\t-15\t26\t33",
  "\t-16\t24\t34",
  "\t-17\t24\t34",
  "\t-18\t24\t34",
  "\t-19\t24\t34",
  "\t-20\t24\t34",
  "\t-21\t25\t34",
  "\t-22\t26\t33",
  "\t-23\t28\t32",
  NULL
};

extern CharPtr water_onedegree [];
CharPtr water_onedegree [] = {
  "1",
  "Adriatic Sea",
  "46\t11\t15",
  "\t45\t11\t16",
  "\t44\t11\t18",
  "\t43\t11\t20",
  "\t42\t11\t20",
  "\t41\t12\t20",
  "\t40\t14\t20",
  "\t39\t16\t20",
  "\t38\t17\t20",
  "Aegean Sea",
  "\t41\t21\t27",
  "\t40\t21\t27",
  "\t39\t21\t28",
  "\t38\t21\t29",
  "\t37\t21\t29",
  "\t36\t23\t29",
  "\t35\t23\t29",
  "Albemarle Sound",
  "\t37\t-77\t-74",
  "\t36\t-77\t-74",
  "\t35\t-77\t-74",
  "\t34\t-77\t-74",
  "Alboran Sea",
  "\t37\t-6\t-1",
  "\t36\t-6\t0",
  "\t35\t-6\t0",
  "\t34\t-6\t0",
  "Amundsen Gulf",
  "\t72\t-126\t-117",
  "\t71\t-128\t-116",
  "\t70\t-128\t-116",
  "\t69\t-128\t-116",
  "\t68\t-127\t-117",
  "Amundsen Sea",
  "\t-71\t-108\t-101",
  "\t-72\t-115\t-97",
  "\t-73\t-115\t-97",
  "\t-74\t-115\t-97",
  "\t-75\t-115\t-97",
  "\t-76\t-112\t-97",
  "Andaman Sea",
  "\t17\t93\t96",
  "\t16\t92\t97",
  "\t15\t91\t99",
  "\t14\t91\t99",
  "\t13\t91\t99",
  "\t12\t91\t99",
  "\t11\t91\t99",
  "\t10\t91\t99",
  "\t9\t91\t99",
  "\t8\t91\t99",
  "\t7\t91\t99",
  "\t6\t92\t99",
  "\t5\t92\t98",
  "\t4\t94\t96",
  "Antongila Bay",
  "\t-14\t48\t51",
  "\t-15\t48\t51",
  "\t-16\t48\t51",
  "\t-17\t48\t50",
  "Arabian Sea",
  "\t26\t60\t67",
  "\t25\t59\t68",
  "\t24\t59\t69",
  "\t23\t58\t70",
  "\t22\t57\t71",
  "\t21\t57\t72",
  "\t20\t56\t74",
  "\t19\t55\t74",
  "\t18\t53\t74",
  "\t17\t51\t74",
  "\t16\t50\t74",
  "\t15\t50\t75",
  "\t14\t50\t75",
  "\t13\t50\t75",
  "\t12\t50\t74",
  "\t11\t50\t72",
  "\t10\t50\t72",
  "\t9\t50\t72",
  "\t8\t52\t72",
  "\t7\t54\t73",
  "\t6\t56\t73",
  "\t5\t58\t73",
  "\t4\t60\t73",
  "\t3\t63\t73",
  "\t2\t65\t73",
  "\t1\t67\t74",
  "\t0\t69\t74",
  "\t-1\t71\t74",
  "Arafura Sea",
  "\t-2\t132\t135",
  "\t-3\t132\t138",
  "\t-4\t131\t139",
  "\t-5\t131\t139",
  "\t-6\t130\t141",
  "\t-7\t129\t141",
  "\t-8\t129\t142",
  "\t-9\t129\t143",
  "\t-10\t129\t143",
  "\t-11\t130\t143",
  "\t-12\t130\t143",
  "\t-13\t133\t142",
  "Aral Sea",
  "\t47\t58\t62",
  "\t46\t57\t62",
  "\t45\t57\t62",
  "\t44\t57\t61",
  "\t43\t57\t61",
  "Arctic Ocean",
  "\t90\t-180\t180",
  "\t89\t-180\t180",
  "\t88\t-180\t180",
  "\t87\t-180\t180",
  "\t86\t-180\t180",
  "\t85\t-180\t180",
  "\t84\t-180\t180",
  "\t83\t-180\t180",
  "\t82\t-180\t180",
  "\t81\t-180\t-69\t-18\t180",
  "\t80\t-180\t-75\t-4\t180",
  "\t79\t-180\t-75\t10\t50\t100\t180",
  "\t78\t-180\t-99\t-88\t-79\t108\t180",
  "\t77\t-180\t-107\t117\t180",
  "\t76\t-180\t-112\t125\t180",
  "\t75\t-180\t-119\t133\t148\t157\t180",
  "\t74\t-180\t-125\t161\t180",
  "\t73\t-180\t-132\t165\t180",
  "\t72\t-180\t-138\t169\t180",
  "\t71\t-180\t-145\t173\t180",
  "\t70\t-180\t-152\t177\t180",
  "Atlantic Ocean",
  "\t69\t-33\t-29",
  "\t68\t-34\t-27",
  "\t67\t-39\t-26",
  "\t66\t-42\t-24",
  "\t65\t-42\t-20\t-17\t-11",
  "\t64\t-43\t-8",
  "\t63\t-43\t-6",
  "\t62\t-43\t-3",
  "\t61\t-44\t0",
  "\t60\t-44\t0",
  "\t59\t-45\t0",
  "\t58\t-45\t0",
  "\t57\t-46\t-1",
  "\t56\t-47\t-5",
  "\t55\t-48\t-6",
  "\t54\t-48\t-7",
  "\t53\t-49\t-5",
  "\t52\t-50\t-5",
  "\t51\t-51\t-4",
  "\t50\t-51\t-4",
  "\t49\t-52\t-4",
  "\t48\t-56\t-4",
  "\t47\t-60\t-4",
  "\t46\t-62\t-5",
  "\t45\t-65\t-5",
  "\t44\t-66\t-6",
  "\t43\t-68\t-6",
  "\t42\t-74\t-6",
  "\t41\t-75\t-7",
  "\t40\t-75\t-7",
  "\t39\t-76\t-7",
  "\t38\t-76\t-5",
  "\t37\t-77\t-4",
  "\t36\t-77\t-4",
  "\t35\t-78\t-4",
  "\t34\t-80\t-4",
  "\t33\t-81\t-5",
  "\t32\t-82\t-5",
  "\t31\t-82\t-7",
  "\t30\t-82\t-8",
  "\t29\t-82\t-8",
  "\t28\t-82\t-8",
  "\t27\t-81\t-9",
  "\t26\t-81\t-11",
  "\t25\t-81\t-12",
  "\t24\t-81\t-13",
  "\t23\t-81\t-13",
  "\t22\t-81\t-14",
  "\t21\t-81\t-15",
  "\t20\t-78\t-15",
  "\t19\t-76\t-15",
  "\t18\t-74\t-15",
  "\t17\t-69\t-15",
  "\t16\t-62\t-15",
  "\t15\t-61\t-15",
  "\t14\t-61\t-14",
  "\t13\t-60\t-14",
  "\t12\t-61\t-14",
  "\t11\t-61\t-13",
  "\t10\t-62\t-12",
  "\t9\t-62\t-11",
  "\t8\t-62\t-10",
  "\t7\t-61\t-9",
  "\t6\t-59\t-8",
  "\t5\t-59\t-5",
  "\t4\t-58\t-2",
  "\t3\t-53\t2",
  "\t2\t-52\t5",
  "\t1\t-51\t7",
  "\t0\t-51\t7",
  "\t-1\t-51\t7",
  "Atlantic Ocean",
  "\t1\t-50\t9",
  "\t0\t-50\t10",
  "\t-1\t-50\t11",
  "\t-2\t-47\t12",
  "\t-3\t-44\t13",
  "\t-4\t-40\t14",
  "\t-5\t-39\t14",
  "\t-6\t-37\t14",
  "\t-7\t-36\t14",
  "\t-8\t-36\t14",
  "\t-9\t-37\t14",
  "\t-10\t-38\t14",
  "\t-11\t-39\t14",
  "\t-12\t-40\t14",
  "\t-13\t-40\t14",
  "\t-14\t-40\t13",
  "\t-15\t-40\t13",
  "\t-16\t-40\t13",
  "\t-17\t-40\t13",
  "\t-18\t-41\t13",
  "\t-19\t-41\t14",
  "\t-20\t-42\t14",
  "\t-21\t-45\t15",
  "\t-22\t-46\t15",
  "\t-23\t-48\t15",
  "\t-24\t-49\t15",
  "\t-25\t-49\t16",
  "\t-26\t-49\t16",
  "\t-27\t-50\t17",
  "\t-28\t-51\t17",
  "\t-29\t-51\t18",
  "\t-30\t-52\t19",
  "\t-31\t-53\t19",
  "\t-32\t-54\t19",
  "\t-33\t-55\t20",
  "\t-34\t-56\t20",
  "\t-35\t-57\t20",
  "\t-36\t-58\t20",
  "\t-37\t-62\t20",
  "\t-38\t-63\t20",
  "\t-39\t-63\t20",
  "\t-40\t-64\t20",
  "\t-41\t-66\t20",
  "\t-42\t-66\t20",
  "\t-43\t-66\t20",
  "\t-44\t-66\t20",
  "\t-45\t-66\t20",
  "\t-46\t-67\t20",
  "\t-47\t-68\t20",
  "\t-48\t-68\t20",
  "\t-49\t-68\t20",
  "\t-50\t-69\t20",
  "\t-51\t-69\t20",
  "\t-52\t-69\t20",
  "\t-53\t-70\t20",
  "\t-54\t-70\t20",
  "\t-55\t-70\t20",
  "\t-56\t-69\t20",
  "\t-57\t-69\t20",
  "\t-58\t-69\t20",
  "\t-59\t-69\t20",
  "\t-60\t-69\t20",
  "\t-61\t-69\t20",
  "Bab el Mandeb",
  "\t14\t42\t44",
  "\t13\t42\t44",
  "\t12\t42\t44",
  "\t11\t42\t44",
  "Baffin Bay",
  "\t79\t-77\t-71",
  "\t78\t-83\t-70",
  "\t77\t-83\t-66",
  "\t76\t-83\t-61",
  "\t75\t-81\t-55",
  "\t74\t-81\t-54",
  "\t73\t-81\t-53",
  "\t72\t-79\t-53",
  "\t71\t-78\t-53",
  "\t70\t-76\t-53",
  "\t69\t-73\t-53",
  "\t68\t-70\t-53",
  "Bahia Blanca",
  "\t-37\t-63\t-60",
  "\t-38\t-63\t-60",
  "\t-39\t-63\t-60",
  "\t-40\t-63\t-60",
  "Bahia de Campeche",
  "\t22\t-94\t-89",
  "\t21\t-98\t-89",
  "\t20\t-98\t-89",
  "\t19\t-98\t-89",
  "\t18\t-97\t-89",
  "\t17\t-96\t-90",
  "Bahia Grande",
  "\t-48\t-69\t-66",
  "\t-49\t-70\t-66",
  "\t-50\t-70\t-66",
  "\t-51\t-70\t-66",
  "\t-52\t-70\t-67",
  "\t-53\t-69\t-67",
  "Bahia Inutil",
  "\t-52\t-71\t-68",
  "\t-53\t-71\t-68",
  "\t-54\t-71\t-68",
  "\t-55\t-71\t-68",
  "Baia de Maputo",
  "\t-24\t31\t33",
  "\t-25\t31\t33",
  "\t-26\t31\t33",
  "\t-27\t31\t33",
  "Baia de Marajo",
  "\t1\t-49\t-47",
  "\t0\t-50\t-47",
  "\t-1\t-50\t-47",
  "\t-2\t-50\t-47",
  "\t-3\t-50\t-48",
  "Baia de Sao Marcos",
  "\t0\t-45\t-43",
  "\t-1\t-45\t-42",
  "\t-2\t-45\t-42",
  "\t-3\t-45\t-42",
  "\t-4\t-45\t-43",
  "Baird Inlet",
  "\t61\t-165\t-162",
  "\t60\t-165\t-162",
  "\t59\t-165\t-162",
  "Balearic Sea",
  "\t42\t0\t4",
  "\t41\t-1\t5",
  "\t40\t-1\t5",
  "\t39\t-1\t5",
  "\t38\t-1\t5",
  "\t37\t-1\t3",
  "Bali Sea",
  "\t-5\t114\t117",
  "\t-6\t113\t118",
  "\t-7\t113\t118",
  "\t-8\t113\t118",
  "\t-9\t113\t118",
  "\t-10\t115\t117",
  "Baltic Sea",
  "\t60\t16\t24",
  "\t59\t15\t24",
  "\t58\t15\t24",
  "\t57\t13\t23",
  "\t56\t11\t23",
  "\t55\t11\t22",
  "\t54\t11\t22",
  "\t53\t11\t21",
  "\t52\t13\t15",
  "Banda Sea",
  "\t1\t121\t124",
  "\t0\t120\t126",
  "\t-1\t119\t129",
  "\t-2\t119\t131",
  "\t-3\t119\t133",
  "\t-4\t119\t134",
  "\t-5\t119\t134",
  "\t-6\t119\t134",
  "\t-7\t119\t133",
  "\t-8\t119\t132",
  "\t-9\t121\t132",
  "Barents Sea",
  "\t82\t49\t66",
  "\t81\t16\t19\t26\t66",
  "\t80\t16\t67",
  "\t79\t16\t67",
  "\t78\t16\t68",
  "\t77\t16\t69",
  "\t76\t16\t69",
  "\t75\t16\t69",
  "\t74\t18\t61",
  "\t73\t20\t57",
  "\t72\t22\t55",
  "\t71\t24\t59",
  "\t70\t26\t61",
  "\t69\t26\t61",
  "\t68\t28\t61",
  "\t67\t36\t61",
  "\t66\t43\t50",
  "\t65\t44\t48",
  "Bass Strait",
  "\t-36\t143\t150",
  "\t-37\t142\t150",
  "\t-38\t142\t150",
  "\t-39\t142\t149",
  "\t-40\t142\t149",
  "\t-41\t142\t149",
  "\t-42\t144\t148",
  "Bathurst Inlet",
  "\t68\t-109\t-106",
  "\t67\t-109\t-106",
  "\t66\t-109\t-106",
  "\t65\t-109\t-106",
  "Bay of Bengal",
  "\t24\t89\t91",
  "\t23\t86\t92",
  "\t22\t85\t93",
  "\t21\t85\t94",
  "\t20\t83\t95",
  "\t19\t82\t95",
  "\t18\t81\t95",
  "\t17\t80\t95",
  "\t16\t79\t95",
  "\t15\t79\t95",
  "\t14\t79\t94",
  "\t13\t78\t93",
  "\t12\t78\t93",
  "\t11\t78\t93",
  "\t10\t78\t93",
  "\t9\t78\t93",
  "\t8\t79\t94",
  "\t7\t79\t95",
  "\t6\t80\t96",
  "\t5\t84\t96",
  "\t4\t91\t96",
  "Bay of Biscay",
  "\t49\t-6\t-3",
  "\t48\t-7\t0",
  "\t47\t-7\t0",
  "\t46\t-8\t1",
  "\t45\t-8\t1",
  "\t44\t-9\t1",
  "\t43\t-9\t0",
  "\t42\t-9\t0",
  "Bay of Fundy",
  "\t46\t-68\t-62",
  "\t45\t-68\t-62",
  "\t44\t-68\t-62",
  "\t43\t-68\t-64",
  "Bay of Plenty",
  "\t-35\t174\t177",
  "\t-36\t174\t179",
  "\t-37\t174\t179",
  "\t-38\t174\t179",
  "Beaufort Sea",
  "\t77\t-126\t-121",
  "\t76\t-133\t-121",
  "\t75\t-139\t-121",
  "\t74\t-146\t-122",
  "\t73\t-153\t-122",
  "\t72\t-157\t-122",
  "\t71\t-157\t-123",
  "\t70\t-157\t-124",
  "\t69\t-157\t-125",
  "\t68\t-145\t-127",
  "Bellingshausen Sea",
  "\t-67\t-74\t-70",
  "\t-68\t-80\t-70",
  "\t-69\t-86\t-68",
  "\t-70\t-92\t-68",
  "\t-71\t-96\t-68",
  "\t-72\t-96\t-68",
  "\t-73\t-96\t-73",
  "\t-74\t-96\t-73",
  "Bering Sea",
  "\t67\t-171\t-168",
  "\t66\t-173\t-165",
  "\t65\t-175\t-163",
  "\t64\t-177\t-163",
  "\t63\t-179\t-163\t174\t180",
  "\t62\t-180\t-163\t171\t180",
  "\t61\t-180\t-160\t165\t180",
  "\t60\t-180\t-160\t165\t180",
  "\t59\t-180\t-160\t163\t180",
  "\t58\t-180\t-160\t161\t180",
  "\t57\t-180\t-160\t161\t180",
  "\t56\t-180\t-160\t161\t180",
  "\t55\t-180\t-160\t161\t180",
  "\t54\t-180\t-160\t163\t180",
  "\t53\t-180\t-161\t165\t180",
  "\t52\t-180\t-163\t167\t180",
  "\t51\t-180\t-166\t169\t180",
  "\t50\t-180\t-171\t171\t180",
  "\t49\t178\t180",
  "Bering Strait",
  "\t67\t-171\t-168",
  "\t66\t-171\t-166",
  "\t65\t-171\t-166",
  "\t64\t-171\t-166",
  "Bight of Benin",
  "\t7\t0\t5",
  "\t6\t-1\t6",
  "\t5\t-1\t6",
  "\t4\t-1\t6",
  "\t3\t2\t6",
  "Bight of Biafra",
  "\t5\t5\t10",
  "\t4\t5\t10",
  "\t3\t5\t10",
  "\t2\t7\t10",
  "\t1\t8\t10",
  "Bismarck Sea",
  "\t0\t141\t148",
  "\t-1\t140\t152",
  "\t-2\t140\t153",
  "\t-3\t140\t153",
  "\t-4\t141\t153",
  "\t-5\t143\t153",
  "\t-6\t144\t152",
  "Black Sea",
  "\t48\t30\t32",
  "\t47\t29\t34",
  "\t46\t28\t37",
  "\t45\t27\t39",
  "\t44\t26\t41",
  "\t43\t26\t42",
  "\t42\t26\t42",
  "\t41\t26\t42",
  "\t40\t26\t42",
  "\t39\t37\t41",
  "Bo Hai",
  "\t41\t119\t123",
  "\t40\t116\t123",
  "\t39\t116\t123",
  "\t38\t116\t122",
  "\t37\t116\t122",
  "\t36\t117\t121",
  "Boca Grande",
  "\t10\t-62\t-59",
  "\t9\t-62\t-59",
  "\t8\t-62\t-59",
  "\t7\t-62\t-59",
  "Bohol Sea",
  "\t11\t122\t126",
  "\t10\t122\t126",
  "\t9\t122\t126",
  "\t8\t122\t126",
  "\t7\t122\t126",
  "Boknafjorden",
  "\t60\t4\t7",
  "\t59\t4\t7",
  "\t58\t4\t7",
  "\t57\t4\t7",
  "Bosporus",
  "\t42\t27\t30",
  "\t41\t27\t30",
  "\t40\t27\t30",
  "Bransfield Strait",
  "\t-60\t-58\t-53",
  "\t-61\t-63\t-53",
  "\t-62\t-63\t-53",
  "\t-63\t-64\t-53",
  "\t-64\t-64\t-54",
  "\t-65\t-64\t-59",
  "\t-66\t-64\t-62",
  "Bristol Bay",
  "\t60\t-161\t-155",
  "\t59\t-163\t-155",
  "\t58\t-163\t-155",
  "\t57\t-163\t-155",
  "\t56\t-163\t-156",
  "\t55\t-162\t-157",
  "\t54\t-162\t-159",
  "Bristol Channel",
  "\t52\t-7\t-1",
  "\t51\t-7\t-1",
  "\t50\t-7\t-1",
  "\t49\t-6\t-3",
  "Caribbean Sea",
  "\t23\t-84\t-79",
  "\t22\t-88\t-77",
  "\t21\t-88\t-73",
  "\t20\t-88\t-71",
  "\t19\t-89\t-60",
  "\t18\t-89\t-60",
  "\t17\t-89\t-59",
  "\t16\t-89\t-59",
  "\t15\t-88\t-58",
  "\t14\t-87\t-58",
  "\t13\t-84\t-58",
  "\t12\t-84\t-58",
  "\t11\t-84\t-58",
  "\t10\t-84\t-59",
  "\t9\t-84\t-59",
  "\t8\t-84\t-74\t-63\t-59",
  "\t7\t-83\t-75",
  "Caspian Sea",
  "\t48\t49\t52",
  "\t47\t47\t54",
  "\t46\t46\t54",
  "\t45\t45\t54",
  "\t44\t45\t54",
  "\t43\t45\t53",
  "\t42\t46\t53",
  "\t41\t46\t53",
  "\t40\t47\t54",
  "\t39\t47\t54",
  "\t38\t47\t54",
  "\t37\t47\t55",
  "\t36\t47\t55",
  "\t35\t49\t55",
  "Celebes Sea",
  "\t8\t121\t125",
  "\t7\t120\t126",
  "\t6\t117\t126",
  "\t5\t116\t126",
  "\t4\t116\t126",
  "\t3\t116\t126",
  "\t2\t116\t126",
  "\t1\t116\t126",
  "\t0\t116\t126",
  "\t-1\t117\t124",
  "Ceram Sea",
  "\t0\t124\t133",
  "\t-1\t124\t134",
  "\t-2\t124\t134",
  "\t-3\t124\t134",
  "\t-4\t124\t126\t129\t134",
  "\t-5\t130\t134",
  "\t-6\t132\t134",
  "Chaun Bay",
  "\t70\t167\t171",
  "\t69\t167\t171",
  "\t68\t167\t171",
  "\t67\t168\t171",
  "Chesapeake Bay",
  "\t40\t-77\t-74",
  "\t39\t-78\t-74",
  "\t38\t-78\t-74",
  "\t37\t-78\t-74",
  "\t36\t-78\t-74",
  "\t35\t-77\t-75",
  "Chukchi Sea",
  "\t72\t-179\t-155\t177\t179",
  "\t71\t-180\t-155\t175\t180",
  "\t70\t-180\t-155\t174\t180",
  "\t69\t-180\t-156\t174\t180",
  "\t68\t-180\t-161\t174\t180",
  "\t67\t-180\t-162\t179\t180",
  "\t66\t-176\t-162",
  "\t65\t-175\t-163",
  "\t64\t-169\t-165",
  "Cook Inlet",
  "\t62\t-152\t-148",
  "\t61\t-154\t-148",
  "\t60\t-155\t-148",
  "\t59\t-155\t-148",
  "\t58\t-155\t-150",
  "\t57\t-154\t-151",
  "Cook Strait",
  "\t-39\t173\t176",
  "\t-40\t173\t176",
  "\t-41\t173\t176",
  "\t-42\t173\t176",
  "Coral Sea",
  "\t-7\t142\t147",
  "\t-8\t141\t148\t164\t168",
  "\t-9\t141\t153\t161\t168",
  "\t-10\t141\t168",
  "\t-11\t141\t168",
  "\t-12\t141\t168",
  "\t-13\t142\t169",
  "\t-14\t142\t169",
  "\t-15\t142\t169",
  "\t-16\t144\t169",
  "\t-17\t144\t170",
  "\t-18\t144\t170",
  "\t-19\t145\t170",
  "\t-20\t145\t170",
  "\t-21\t147\t170",
  "\t-22\t148\t169",
  "\t-23\t148\t168",
  "\t-24\t149\t167",
  "\t-25\t150\t166",
  "\t-26\t151\t165",
  "\t-27\t152\t164",
  "\t-28\t152\t162",
  "\t-29\t152\t161",
  "\t-30\t152\t160",
  "Cordova Bay",
  "\t56\t-134\t-131",
  "\t55\t-134\t-131",
  "\t54\t-134\t-131",
  "\t53\t-133\t-131",
  "Cumberland Sound",
  "\t67\t-69\t-63",
  "\t66\t-69\t-62",
  "\t65\t-69\t-62",
  "\t64\t-69\t-62",
  "\t63\t-67\t-62",
  "\t62\t-65\t-63",
  "Dardanelles",
  "\t41\t25\t27",
  "\t40\t25\t27",
  "\t39\t25\t27",
  "\t38\t25\t27",
  "Darnley Bay",
  "\t70\t-125\t-122",
  "\t69\t-125\t-122",
  "\t68\t-125\t-122",
  "Davao Gulf",
  "\t8\t124\t126",
  "\t7\t124\t127",
  "\t6\t124\t127",
  "\t5\t124\t127",
  "\t4\t124\t126",
  "Davis Sea",
  "\t-62\t90\t104",
  "\t-63\t86\t111",
  "\t-64\t84\t113",
  "\t-65\t83\t113",
  "\t-66\t82\t113",
  "\t-67\t82\t111",
  "\t-68\t82\t87",
  "Davis Strait",
  "\t70\t-70\t-52",
  "\t69\t-70\t-50",
  "\t68\t-70\t-49",
  "\t67\t-70\t-49",
  "\t66\t-67\t-49",
  "\t65\t-64\t-49",
  "\t64\t-66\t-48",
  "\t63\t-66\t-47",
  "\t62\t-66\t-44",
  "\t61\t-66\t-43",
  "\t60\t-65\t-43",
  "\t59\t-65\t-43",
  "Delaware Bay",
  "\t40\t-76\t-73",
  "\t39\t-76\t-73",
  "\t38\t-76\t-73",
  "\t37\t-76\t-73",
  "Denmark Strait",
  "\t71\t-23\t-21",
  "\t70\t-26\t-19",
  "\t69\t-31\t-18",
  "\t68\t-31\t-16",
  "\t67\t-31\t-15",
  "\t66\t-30\t-15",
  "\t65\t-28\t-15",
  "\t64\t-27\t-16",
  "\t63\t-25\t-22",
  "Disko Bay",
  "\t71\t-55\t-49",
  "\t70\t-55\t-49",
  "\t69\t-55\t-49",
  "\t68\t-54\t-49",
  "\t67\t-54\t-49",
  "Dixon Entrance",
  "\t55\t-134\t-130",
  "\t54\t-134\t-130",
  "\t53\t-134\t-130",
  "Dmitriy Laptev Strait",
  "\t74\t138\t144",
  "\t73\t138\t144",
  "\t72\t138\t144",
  "\t71\t139\t144",
  "Drake Passage",
  "\t-53\t-67\t-62",
  "\t-54\t-69\t-61",
  "\t-55\t-69\t-60",
  "\t-56\t-69\t-58",
  "\t-57\t-69\t-57",
  "\t-58\t-69\t-56",
  "\t-59\t-69\t-55",
  "\t-60\t-69\t-54",
  "\t-61\t-69\t-54",
  "\t-62\t-69\t-54",
  "\t-63\t-69\t-57",
  "\t-64\t-69\t-61",
  "\t-65\t-69\t-62",
  "\t-66\t-69\t-64",
  "\t-67\t-69\t-65",
  "East China Sea",
  "\t34\t124\t127\t129\t131",
  "\t33\t122\t131",
  "\t32\t120\t131",
  "\t31\t120\t131",
  "\t30\t120\t131",
  "\t29\t119\t131",
  "\t28\t119\t131",
  "\t27\t118\t130",
  "\t26\t118\t129",
  "\t25\t118\t129",
  "\t24\t118\t128",
  "\t23\t120\t127",
  "\t22\t122\t125",
  "East Korea Bay",
  "\t41\t127\t129",
  "\t40\t126\t129",
  "\t39\t126\t129",
  "\t38\t126\t129",
  "\t37\t126\t129",
  "East Siberian Sea",
  "\t78\t147\t158",
  "\t77\t137\t162",
  "\t76\t137\t166",
  "\t75\t137\t170",
  "\t74\t138\t174",
  "\t73\t138\t178",
  "\t72\t138\t180",
  "\t71\t142\t180",
  "\t70\t147\t180",
  "\t69\t150\t155\t157\t178",
  "\t68\t158\t176",
  "\t67\t159\t162",
  "Eclipse Sound",
  "\t74\t-81\t-79",
  "\t73\t-82\t-76",
  "\t72\t-82\t-76",
  "\t71\t-82\t-76",
  "\t70\t-81\t-77",
  "English Channel",
  "\t52\t0\t2",
  "\t51\t-6\t2",
  "\t50\t-7\t2",
  "\t49\t-7\t2",
  "\t48\t-7\t2",
  "\t47\t-6\t0",
  "Eskimo Lakes",
  "\t70\t-134\t-130",
  "\t69\t-134\t-130",
  "\t68\t-134\t-130",
  "\t67\t-134\t-131",
  "Estrecho de Magellanes",
  "\t-51\t-75\t-67",
  "\t-52\t-75\t-67",
  "\t-53\t-75\t-67",
  "\t-54\t-74\t-69",
  "\t-55\t-72\t-69",
  "Finger Lakes",
  "\t43\t-78\t-75",
  "\t42\t-78\t-75",
  "\t41\t-78\t-75",
  "Flores Sea",
  "\t-4\t118\t121",
  "\t-5\t117\t121",
  "\t-6\t116\t122",
  "\t-7\t116\t123",
  "\t-8\t116\t123",
  "\t-9\t116\t123",
  "Foxe Basin",
  "\t71\t-80\t-76",
  "\t70\t-83\t-73",
  "\t69\t-83\t-72",
  "\t68\t-85\t-71",
  "\t67\t-87\t-71",
  "\t66\t-87\t-71",
  "\t65\t-87\t-71",
  "\t64\t-86\t-72",
  "\t63\t-84\t-74",
  "\t62\t-81\t-78",
  "Franklin Bay",
  "\t70\t-126\t-124",
  "\t69\t-126\t-124",
  "\t68\t-126\t-124",
  "Frobisher Bay",
  "\t64\t-69\t-64",
  "\t63\t-69\t-64",
  "\t62\t-69\t-64",
  "\t61\t-68\t-64",
  "Fury and Hecla Strait",
  "\t71\t-86\t-82",
  "\t70\t-86\t-81",
  "\t69\t-86\t-81",
  "\t68\t-86\t-81",
  "Garabogaz Bay",
  "\t43\t52\t54",
  "\t42\t51\t55",
  "\t41\t51\t55",
  "\t40\t51\t55",
  "\t39\t51\t55",
  "Geographe Bay",
  "\t-29\t114\t116",
  "\t-30\t114\t116",
  "\t-31\t114\t116",
  "\t-32\t114\t116",
  "\t-33\t114\t116",
  "\t-34\t114\t116",
  "George VI Sound",
  "\t-68\t-70\t-67",
  "\t-69\t-70\t-66",
  "\t-70\t-74\t-65",
  "\t-71\t-75\t-65",
  "\t-72\t-75\t-65",
  "\t-73\t-75\t-65",
  "\t-74\t-75\t-67",
  "Goldsmith Channel",
  "\t74\t-108\t-104",
  "\t73\t-108\t-104",
  "\t72\t-108\t-104",
  "\t71\t-106\t-104",
  "Golfe du Lion",
  "\t44\t2\t6",
  "\t43\t2\t6",
  "\t42\t2\t6",
  "\t41\t2\t5",
  "\t40\t2\t4",
  "Golfo Corcovado",
  "\t-40\t-74\t-71",
  "\t-41\t-74\t-71",
  "\t-42\t-74\t-71",
  "\t-43\t-75\t-71",
  "\t-44\t-75\t-71",
  "\t-45\t-75\t-71",
  "\t-46\t-74\t-71",
  "Golfo de California",
  "\t32\t-115\t-112",
  "\t31\t-115\t-111",
  "\t30\t-115\t-111",
  "\t29\t-115\t-110",
  "\t28\t-115\t-108",
  "\t27\t-114\t-108",
  "\t26\t-113\t-107",
  "\t25\t-113\t-106",
  "\t24\t-112\t-105",
  "\t23\t-111\t-105",
  "\t22\t-110\t-105",
  "Golfo de Guayaquil",
  "\t-1\t-81\t-78",
  "\t-2\t-81\t-78",
  "\t-3\t-81\t-78",
  "\t-4\t-81\t-78",
  "Golfo de Panama",
  "\t10\t-80\t-78",
  "\t9\t-81\t-76",
  "\t8\t-81\t-76",
  "\t7\t-81\t-76",
  "\t6\t-81\t-77",
  "Golfo de Penas",
  "\t-45\t-76\t-73",
  "\t-46\t-76\t-73",
  "\t-47\t-76\t-73",
  "\t-48\t-76\t-73",
  "Golfo de Tehuantepec",
  "\t17\t-96\t-92",
  "\t16\t-97\t-92",
  "\t15\t-97\t-92",
  "\t14\t-97\t-92",
  "Golfo de Uraba",
  "\t9\t-78\t-75",
  "\t8\t-78\t-75",
  "\t7\t-78\t-75",
  "\t6\t-77\t-75",
  "Golfo San Jorge",
  "\t-43\t-67\t-65",
  "\t-44\t-68\t-64",
  "\t-45\t-68\t-64",
  "\t-46\t-68\t-64",
  "\t-47\t-68\t-64",
  "\t-48\t-67\t-64",
  "Golfo San Matias",
  "\t-39\t-66\t-63",
  "\t-40\t-66\t-62",
  "\t-41\t-66\t-62",
  "\t-42\t-66\t-62",
  "\t-43\t-65\t-62",
  "Great Australian Bight",
  "\t-30\t127\t133",
  "\t-31\t123\t135",
  "\t-32\t118\t136",
  "\t-33\t117\t136",
  "\t-34\t117\t140",
  "\t-35\t117\t140",
  "\t-36\t117\t141",
  "\t-37\t119\t144",
  "\t-38\t123\t144",
  "\t-39\t126\t145",
  "\t-40\t129\t146",
  "\t-41\t133\t146",
  "\t-42\t136\t147",
  "\t-43\t139\t147",
  "\t-44\t143\t147",
  "Great Barrier Reef",
  "\t-8\t141\t146",
  "\t-9\t141\t146",
  "\t-10\t141\t146",
  "\t-11\t141\t146",
  "\t-12\t141\t147",
  "\t-13\t142\t148",
  "\t-14\t142\t148",
  "\t-15\t142\t148",
  "\t-16\t144\t149",
  "\t-17\t144\t150",
  "\t-18\t144\t151",
  "\t-19\t145\t151",
  "\t-20\t145\t152",
  "\t-21\t147\t154",
  "\t-22\t148\t154",
  "\t-23\t148\t154",
  "\t-24\t149\t154",
  "\t-25\t150\t154",
  "\t-26\t151\t154",
  "Great Bear Lake",
  "\t68\t-121\t-118",
  "\t67\t-126\t-116",
  "\t66\t-126\t-116",
  "\t65\t-126\t-116",
  "\t64\t-125\t-116",
  "\t63\t-123\t-119",
  "Great Salt Lake",
  "\t42\t-114\t-110",
  "\t41\t-114\t-110",
  "\t40\t-114\t-110",
  "\t39\t-113\t-110",
  "Great Slave Lake",
  "\t63\t-117\t-108",
  "\t62\t-118\t-108",
  "\t61\t-118\t-108",
  "\t60\t-118\t-110",
  "\t59\t-117\t-113",
  "Greenland Sea",
  "\t84\t-32\t-17",
  "\t83\t-33\t-3",
  "\t82\t-33\t11",
  "\t81\t-33\t18",
  "\t80\t-30\t-27\t-25\t18",
  "\t79\t-24\t18",
  "\t78\t-22\t18",
  "\t77\t-23\t18",
  "\t76\t-23\t18",
  "\t75\t-23\t17",
  "\t74\t-28\t14",
  "\t73\t-28\t10",
  "\t72\t-28\t5",
  "\t71\t-27\t0",
  "\t70\t-26\t-4",
  "\t69\t-29\t-7",
  "\t68\t-29\t-9",
  "\t67\t-29\t-10",
  "\t66\t-27\t-10",
  "\t65\t-25\t-11",
  "\t64\t-24\t-12",
  "Guba Gusinaya",
  "\t73\t144\t148",
  "\t72\t143\t148",
  "\t71\t143\t148",
  "\t70\t143\t147",
  "Gulf of Aden",
  "\t16\t49\t52",
  "\t15\t46\t52",
  "\t14\t44\t52",
  "\t13\t42\t52",
  "\t12\t41\t52",
  "\t11\t41\t52",
  "\t10\t41\t52",
  "\t9\t42\t47",
  "Gulf of Alaska",
  "\t61\t-150\t-138",
  "\t60\t-152\t-137",
  "\t59\t-156\t-135",
  "\t58\t-157\t-135",
  "\t57\t-159\t-135",
  "\t56\t-164\t-139",
  "\t55\t-164\t-145",
  "\t54\t-164\t-152",
  "\t53\t-164\t-158",
  "Gulf of Anadyr",
  "\t67\t-180\t-177",
  "\t66\t-180\t-174",
  "\t65\t-180\t-172",
  "\t64\t-180\t-172",
  "\t63\t-180\t-172",
  "\t62\t-180\t-174",
  "\t61\t-180\t-176",
  "\t60\t-180\t-178",
  "Gulf of Anadyr",
  "\t66\t175\t180",
  "\t65\t173\t180",
  "\t64\t173\t180",
  "\t63\t173\t180",
  "\t62\t177\t180",
  "\t61\t178\t180",
  "Gulf of Aqaba",
  "\t30\t33\t35",
  "\t29\t33\t35",
  "\t28\t33\t35",
  "\t27\t33\t35",
  "\t26\t33\t35",
  "Gulf of Boothia",
  "\t72\t-90\t-88",
  "\t71\t-93\t-84",
  "\t70\t-93\t-83",
  "\t69\t-93\t-83",
  "\t68\t-93\t-83",
  "\t67\t-91\t-83",
  "\t66\t-89\t-85",
  "Gulf of Bothnia",
  "\t66\t20\t26",
  "\t65\t20\t26",
  "\t64\t17\t26",
  "\t63\t16\t26",
  "\t62\t16\t24",
  "\t61\t16\t24",
  "\t60\t16\t24",
  "\t59\t16\t24",
  "\t58\t17\t24",
  "Gulf of Buli",
  "\t2\t127\t130",
  "\t1\t127\t130",
  "\t0\t127\t130",
  "\t-1\t127\t130",
  "Gulf of Carpentaria",
  "\t-11\t135\t142",
  "\t-12\t134\t142",
  "\t-13\t134\t142",
  "\t-14\t134\t142",
  "\t-15\t134\t142",
  "\t-16\t134\t142",
  "\t-17\t136\t142",
  "\t-18\t138\t141",
  "Gulf of Finland",
  "\t61\t23\t31",
  "\t60\t21\t31",
  "\t59\t21\t31",
  "\t58\t21\t31",
  "Gulf of Gabes",
  "\t36\t9\t12",
  "\t35\t9\t12",
  "\t34\t9\t12",
  "\t33\t9\t12",
  "\t32\t9\t11",
  "Gulf of Guinea",
  "\t6\t-6\t3",
  "\t5\t-8\t8",
  "\t4\t-8\t9",
  "\t3\t-8\t10",
  "\t2\t-6\t10",
  "\t1\t-3\t11",
  "\t0\t1\t11",
  "\t-1\t4\t11",
  "Gulf of Honduras",
  "\t18\t-89\t-87",
  "\t17\t-89\t-86",
  "\t16\t-89\t-85",
  "\t15\t-89\t-85",
  "\t14\t-89\t-85",
  "Gulf of Kamchatka",
  "\t57\t161\t164",
  "\t56\t160\t164",
  "\t55\t160\t164",
  "\t54\t160\t164",
  "\t53\t160\t163",
  "Gulf of Kau",
  "\t3\t127\t129",
  "\t2\t126\t130",
  "\t1\t126\t130",
  "\t0\t126\t130",
  "\t-1\t126\t128",
  "Gulf of Khambhat",
  "\t23\t71\t73",
  "\t22\t71\t74",
  "\t21\t69\t74",
  "\t20\t69\t74",
  "\t19\t69\t73",
  "\t18\t71\t73",
  "Gulf of Kutch",
  "\t24\t67\t71",
  "\t23\t67\t71",
  "\t22\t67\t71",
  "\t21\t67\t71",
  "Gulf of Maine",
  "\t45\t-70\t-65",
  "\t44\t-71\t-64",
  "\t43\t-71\t-64",
  "\t42\t-71\t-64",
  "\t41\t-71\t-65",
  "\t40\t-70\t-67",
  "Gulf of Mannar",
  "\t10\t77\t80",
  "\t9\t76\t80",
  "\t8\t76\t80",
  "\t7\t76\t80",
  "\t6\t78\t80",
  "Gulf of Martaban",
  "\t18\t95\t98",
  "\t17\t94\t98",
  "\t16\t94\t98",
  "\t15\t94\t98",
  "\t14\t94\t98",
  "\t13\t96\t98",
  "Gulf of Masira",
  "\t21\t56\t59",
  "\t20\t56\t59",
  "\t19\t56\t59",
  "\t18\t56\t58",
  "Gulf of Mexico",
  "\t31\t-90\t-83",
  "\t30\t-96\t-81",
  "\t29\t-98\t-81",
  "\t28\t-98\t-81",
  "\t27\t-98\t-80",
  "\t26\t-98\t-79",
  "\t25\t-98\t-78",
  "\t24\t-98\t-78",
  "\t23\t-98\t-78",
  "\t22\t-98\t-82",
  "\t21\t-98\t-82",
  "\t20\t-98\t-83",
  "\t19\t-98\t-93",
  "Gulf of Ob",
  "\t73\t71\t76",
  "\t72\t70\t76",
  "\t71\t70\t76",
  "\t70\t70\t77",
  "\t69\t71\t78",
  "\t68\t70\t79",
  "\t67\t68\t79",
  "\t66\t68\t79",
  "\t65\t68\t74",
  "Gulf of Olenek",
  "\t74\t117\t124",
  "\t73\t117\t124",
  "\t72\t117\t124",
  "\t71\t118\t124",
  "Gulf of Oman",
  "\t27\t55\t58",
  "\t26\t55\t62",
  "\t25\t55\t62",
  "\t24\t55\t62",
  "\t23\t55\t61",
  "\t22\t56\t61",
  "\t21\t58\t60",
  "Gulf of Papua",
  "\t-6\t142\t146",
  "\t-7\t141\t147",
  "\t-8\t141\t147",
  "\t-9\t141\t147",
  "Gulf of Riga",
  "\t60\t22\t24",
  "\t59\t21\t25",
  "\t58\t20\t25",
  "\t57\t20\t25",
  "\t56\t20\t25",
  "\t55\t22\t24",
  "Gulf of Sakhalin",
  "\t55\t138\t143",
  "\t54\t138\t143",
  "\t53\t138\t143",
  "\t52\t139\t143",
  "Gulf of Sidra",
  "\t33\t14\t20",
  "\t32\t14\t21",
  "\t31\t14\t21",
  "\t30\t14\t21",
  "\t29\t16\t21",
  "Gulf of St. Lawrence",
  "\t52\t-59\t-55",
  "\t51\t-65\t-55",
  "\t50\t-65\t-55",
  "\t49\t-67\t-56",
  "\t48\t-67\t-53",
  "\t47\t-67\t-53",
  "\t46\t-67\t-53",
  "\t45\t-65\t-54",
  "\t44\t-64\t-60",
  "Gulf of Suez",
  "\t30\t31\t34",
  "\t29\t31\t34",
  "\t28\t31\t35",
  "\t27\t31\t35",
  "\t26\t32\t35",
  "Gulf of Thailand",
  "\t14\t98\t101",
  "\t13\t98\t103",
  "\t12\t98\t104",
  "\t11\t98\t106",
  "\t10\t98\t106",
  "\t9\t98\t106",
  "\t8\t98\t106",
  "\t7\t98\t105",
  "\t6\t99\t104",
  "\t5\t99\t103",
  "Gulf of Tomini",
  "\t1\t119\t124",
  "\t0\t119\t124",
  "\t-1\t119\t124",
  "\t-2\t119\t122",
  "Gulf of Tonkin",
  "\t22\t105\t110",
  "\t21\t105\t111",
  "\t20\t104\t111",
  "\t19\t104\t111",
  "\t18\t104\t111",
  "\t17\t104\t109",
  "\t16\t105\t108",
  "Gulf of Yana",
  "\t76\t135\t138",
  "\t75\t135\t141",
  "\t74\t135\t141",
  "\t73\t133\t142",
  "\t72\t131\t142",
  "\t71\t131\t142",
  "\t70\t131\t140",
  "Gulf St. Vincent",
  "\t-31\t136\t138",
  "\t-32\t135\t138",
  "\t-33\t134\t139",
  "\t-34\t134\t139",
  "\t-35\t134\t139",
  "\t-36\t135\t139",
  "Hadley Bay",
  "\t74\t-109\t-107",
  "\t73\t-109\t-106",
  "\t72\t-109\t-106",
  "\t71\t-109\t-106",
  "\t70\t-109\t-106",
  "Hall Basin",
  "\t83\t-63\t-61",
  "\t82\t-69\t-60",
  "\t81\t-69\t-60",
  "\t80\t-69\t-60",
  "\t79\t-64\t-62",
  "Halmahera Sea",
  "\t1\t126\t131",
  "\t0\t126\t132",
  "\t-1\t126\t132",
  "\t-2\t126\t132",
  "Hamilton Inlet",
  "\t55\t-59\t-56",
  "\t54\t-61\t-56",
  "\t53\t-61\t-56",
  "\t52\t-61\t-57",
  "Hangzhou Bay",
  "\t31\t119\t123",
  "\t30\t119\t123",
  "\t29\t119\t123",
  "\t28\t120\t123",
  "Hecate Straight",
  "\t56\t-133\t-129",
  "\t55\t-133\t-128",
  "\t54\t-133\t-128",
  "\t53\t-133\t-128",
  "\t52\t-133\t-128",
  "\t51\t-132\t-129",
  "Hudson Bay",
  "\t67\t-87\t-84",
  "\t66\t-88\t-84",
  "\t65\t-94\t-81",
  "\t64\t-94\t-78",
  "\t63\t-94\t-77",
  "\t62\t-95\t-76",
  "\t61\t-95\t-76",
  "\t60\t-95\t-76",
  "\t59\t-95\t-76",
  "\t58\t-95\t-75",
  "\t57\t-95\t-75",
  "\t56\t-93\t-75",
  "\t55\t-93\t-75",
  "\t54\t-88\t-75",
  "\t53\t-83\t-77",
  "Hudson Strait",
  "\t65\t-79\t-71",
  "\t64\t-81\t-69",
  "\t63\t-81\t-64",
  "\t62\t-81\t-63",
  "\t61\t-79\t-63",
  "\t60\t-73\t-63",
  "\t59\t-71\t-63",
  "IJsselmeer",
  "\t54\t4\t6",
  "\t53\t3\t6",
  "\t52\t3\t6",
  "\t51\t3\t6",
  "Indian Ocean",
  "\t11\t49\t53",
  "\t10\t49\t55",
  "\t9\t49\t57",
  "\t8\t48\t59\t80\t85",
  "\t7\t48\t61\t79\t92",
  "\t6\t47\t64\t78\t96",
  "\t5\t46\t66\t77\t97",
  "\t4\t45\t68\t75\t98",
  "\t3\t44\t70\t74\t99",
  "\t2\t43\t99",
  "\t1\t40\t101",
  "\t0\t39\t101",
  "\t-1\t39\t102",
  "\t-2\t38\t103",
  "\t-3\t38\t104",
  "\t-4\t37\t105",
  "\t-5\t37\t107",
  "\t-6\t37\t111",
  "\t-7\t37\t119",
  "\t-8\t38\t120",
  "\t-9\t38\t123",
  "\t-10\t38\t125",
  "\t-11\t38\t126",
  "\t-12\t43\t127",
  "\t-13\t48\t127",
  "\t-14\t48\t127",
  "\t-15\t48\t127",
  "\t-16\t48\t126",
  "\t-17\t48\t125",
  "\t-18\t47\t124",
  "\t-19\t47\t123",
  "\t-20\t47\t122",
  "\t-21\t46\t120",
  "\t-22\t46\t117",
  "\t-23\t46\t115",
  "\t-24\t38\t114",
  "\t-25\t31\t114",
  "\t-26\t31\t115",
  "\t-27\t30\t115",
  "\t-28\t30\t115",
  "\t-29\t29\t116",
  "\t-30\t28\t116",
  "\t-31\t27\t116",
  "\t-32\t22\t116",
  "\t-33\t18\t117",
  "\t-34\t18\t120",
  "\t-35\t18\t124",
  "\t-36\t18\t127",
  "\t-37\t18\t130",
  "\t-38\t18\t134",
  "\t-39\t18\t137",
  "\t-40\t18\t140",
  "\t-41\t18\t144",
  "\t-42\t18\t148",
  "\t-43\t18\t151",
  "\t-44\t18\t153",
  "\t-45\t18\t156",
  "\t-46\t18\t159",
  "\t-47\t18\t161",
  "\t-48\t18\t164",
  "\t-49\t18\t167",
  "\t-50\t18\t167",
  "\t-51\t18\t167",
  "\t-52\t18\t167",
  "\t-53\t18\t167",
  "\t-54\t18\t167",
  "\t-55\t18\t167",
  "\t-56\t18\t167",
  "\t-57\t18\t167",
  "\t-58\t18\t167",
  "\t-59\t18\t167",
  "\t-60\t18\t167",
  "\t-61\t18\t167",
  "Inner Sea",
  "\t35\t129\t136",
  "\t34\t129\t136",
  "\t33\t129\t136",
  "\t32\t129\t136",
  "\t31\t130\t133",
  "Inner Seas",
  "\t59\t-7\t-4",
  "\t58\t-8\t-4",
  "\t57\t-8\t-3",
  "\t56\t-9\t-3",
  "\t55\t-9\t-3",
  "\t54\t-9\t-3",
  "\t53\t-8\t-4",
  "Internal Canada (B.C.) Waters",
  "\t54\t-130\t-126",
  "\t53\t-130\t-126",
  "\t52\t-130\t-126",
  "\t51\t-130\t-126",
  "Internal Canada (B.C.) Waters",
  "\t55\t-133\t-131",
  "\t54\t-133\t-131",
  "\t53\t-133\t-131",
  "\t52\t-133\t-131",
  "Internal Canada (B.C.) Waters",
  "\t56\t-131\t-128",
  "\t55\t-131\t-128",
  "\t54\t-131\t-128",
  "\t53\t-131\t-129",
  "Internal Canada Arctic Waters",
  "\t72\t-119\t-116",
  "\t71\t-119\t-116",
  "\t70\t-119\t-116",
  "Internal Canada Arctic Waters",
  "\t69\t-108\t-104",
  "\t68\t-108\t-104",
  "\t67\t-108\t-104",
  "Internal Canada Arctic Waters",
  "\t71\t-82\t-79",
  "\t70\t-82\t-78",
  "\t69\t-82\t-78",
  "\t68\t-82\t-78",
  "Internal Denmark Waters",
  "\t57\t9\t11",
  "\t56\t8\t12",
  "\t55\t8\t12",
  "\t54\t8\t12",
  "\t53\t8\t12",
  "Internal Philippines Waters",
  "\t11\t124\t127",
  "\t10\t124\t127",
  "\t9\t124\t127",
  "\t8\t124\t127",
  "Internal Philippines Waters",
  "\t11\t122\t124",
  "\t10\t122\t124",
  "\t9\t122\t124",
  "\t8\t122\t124",
  "Internal U.S. (Alaska) Waters",
  "\t60\t-138\t-134",
  "\t59\t-138\t-132",
  "\t58\t-138\t-131",
  "\t57\t-138\t-130",
  "\t56\t-136\t-129",
  "\t55\t-135\t-129",
  "\t54\t-134\t-129",
  "Ionian Sea",
  "\t41\t15\t18",
  "\t40\t15\t22",
  "\t39\t14\t24",
  "\t38\t14\t24",
  "\t37\t14\t24",
  "\t36\t14\t23",
  "\t35\t14\t23",
  "Irish Sea",
  "\t55\t-7\t-1",
  "\t54\t-7\t-1",
  "\t53\t-7\t-1",
  "\t52\t-7\t-1",
  "\t51\t-7\t-2",
  "\t50\t-7\t-4",
  "James Bay",
  "\t55\t-83\t-77",
  "\t54\t-83\t-77",
  "\t53\t-83\t-77",
  "\t52\t-83\t-77",
  "\t51\t-83\t-77",
  "\t50\t-81\t-77",
  "\t49\t-80\t-78",
  "Java Sea",
  "\t-1\t105\t114",
  "\t-2\t104\t117",
  "\t-3\t104\t119",
  "\t-4\t103\t120",
  "\t-5\t103\t120",
  "\t-6\t103\t120",
  "\t-7\t104\t119",
  "\t-8\t111\t118",
  "Jones Sound",
  "\t77\t-92\t-77",
  "\t76\t-92\t-77",
  "\t75\t-92\t-77",
  "\t74\t-91\t-78",
  "Joseph Bonaparte Gulf",
  "\t-12\t126\t130",
  "\t-13\t126\t130",
  "\t-14\t126\t130",
  "\t-15\t126\t130",
  "\t-16\t127\t130",
  "Kaliningrad",
  "\t56\t19\t22",
  "\t55\t19\t22",
  "\t54\t19\t22",
  "\t53\t19\t22",
  "Kane Basin",
  "\t81\t-73\t-63",
  "\t80\t-79\t-63",
  "\t79\t-79\t-63",
  "\t78\t-79\t-63",
  "\t77\t-79\t-67",
  "Kangertittivaq",
  "\t72\t-29\t-23",
  "\t71\t-30\t-20",
  "\t70\t-30\t-20",
  "\t69\t-30\t-20",
  "\t68\t-28\t-26",
  "Kara Sea",
  "\t82\t64\t96",
  "\t81\t64\t98",
  "\t80\t64\t103",
  "\t79\t64\t103",
  "\t78\t65\t103",
  "\t77\t65\t102",
  "\t76\t59\t102",
  "\t75\t56\t102",
  "\t74\t55\t100",
  "\t73\t54\t88",
  "\t72\t54\t88",
  "\t71\t54\t80",
  "\t70\t54\t69\t74\t80",
  "\t69\t55\t70\t77\t80",
  "\t68\t59\t70",
  "\t67\t65\t70",
  "Karaginskiy Gulf",
  "\t61\t162\t167",
  "\t60\t161\t167",
  "\t59\t160\t167",
  "\t58\t160\t167",
  "\t57\t160\t166",
  "\t56\t161\t164",
  "Karskiye Strait",
  "\t71\t56\t60",
  "\t70\t56\t60",
  "\t69\t56\t60",
  "Kattegat",
  "\t59\t10\t12",
  "\t58\t9\t13",
  "\t57\t9\t13",
  "\t56\t9\t13",
  "\t55\t9\t13",
  "\t54\t10\t12",
  "Kennedy Channel",
  "\t82\t-67\t-63",
  "\t81\t-68\t-63",
  "\t80\t-68\t-63",
  "\t79\t-68\t-63",
  "Khatanga Gulf",
  "\t76\t111\t114",
  "\t75\t108\t114",
  "\t74\t105\t114",
  "\t73\t104\t114",
  "\t72\t104\t113",
  "\t71\t104\t107",
  "Korea Strait",
  "\t37\t128\t131",
  "\t36\t126\t133",
  "\t35\t125\t133",
  "\t34\t125\t133",
  "\t33\t125\t133",
  "\t32\t125\t131",
  "\t31\t126\t130",
  "Kotzebue Sound",
  "\t68\t-164\t-160",
  "\t67\t-165\t-159",
  "\t66\t-165\t-159",
  "\t65\t-165\t-159",
  "Kronotskiy Gulf",
  "\t55\t158\t162",
  "\t54\t158\t162",
  "\t53\t158\t162",
  "\t52\t158\t161",
  "La Perouse Strait",
  "\t47\t140\t142",
  "\t46\t140\t143",
  "\t45\t140\t143",
  "\t44\t140\t143",
  "Labrador Sea",
  "\t61\t-65\t-43",
  "\t60\t-65\t-42",
  "\t59\t-65\t-42",
  "\t58\t-65\t-42",
  "\t57\t-64\t-43",
  "\t56\t-63\t-43",
  "\t55\t-63\t-44",
  "\t54\t-62\t-45",
  "\t53\t-60\t-46",
  "\t52\t-58\t-46",
  "\t51\t-57\t-47",
  "\t50\t-57\t-48",
  "\t49\t-57\t-49",
  "\t48\t-57\t-49",
  "\t47\t-55\t-50",
  "\t46\t-54\t-51",
  "Laccadive Sea",
  "\t15\t73\t75",
  "\t14\t70\t75",
  "\t13\t70\t76",
  "\t12\t70\t76",
  "\t11\t70\t77",
  "\t10\t70\t77",
  "\t9\t70\t79",
  "\t8\t70\t80",
  "\t7\t70\t81",
  "\t6\t70\t81",
  "\t5\t71\t81",
  "\t4\t71\t80",
  "\t3\t71\t79",
  "\t2\t71\t78",
  "\t1\t71\t76",
  "\t0\t71\t75",
  "\t-1\t71\t74",
  "Lago de Maracaibo",
  "\t11\t-72\t-70",
  "\t10\t-73\t-70",
  "\t9\t-73\t-70",
  "\t8\t-73\t-70",
  "Lake Baikal",
  "\t56\t107\t110",
  "\t55\t107\t110",
  "\t54\t105\t110",
  "\t53\t104\t110",
  "\t52\t102\t110",
  "\t51\t102\t109",
  "\t50\t102\t107",
  "Lake Chad",
  "\t14\t13\t15",
  "\t13\t13\t15",
  "\t12\t13\t15",
  "\t11\t13\t15",
  "Lake Champlain",
  "\t46\t-74\t-72",
  "\t45\t-74\t-72",
  "\t44\t-74\t-72",
  "\t43\t-74\t-72",
  "\t42\t-74\t-72",
  "Lake Erie",
  "\t44\t-80\t-77",
  "\t43\t-84\t-77",
  "\t42\t-84\t-77",
  "\t41\t-84\t-77",
  "\t40\t-84\t-79",
  "Lake Huron",
  "\t47\t-82\t-80",
  "\t46\t-82\t-78",
  "\t45\t-82\t-78",
  "\t44\t-82\t-78",
  "\t43\t-82\t-78",
  "Lake Huron",
  "\t47\t-85\t-80",
  "\t46\t-85\t-78",
  "\t45\t-85\t-78",
  "\t44\t-85\t-78",
  "\t43\t-84\t-78",
  "\t42\t-84\t-80",
  "Lake Huron",
  "\t47\t-84\t-80",
  "\t46\t-84\t-80",
  "\t45\t-84\t-80",
  "\t44\t-84\t-80",
  "Lake Huron",
  "\t45\t-84\t-82",
  "\t44\t-84\t-82",
  "\t43\t-84\t-82",
  "\t42\t-84\t-82",
  "Lake Malawi",
  "\t-8\t32\t35",
  "\t-9\t32\t35",
  "\t-10\t32\t35",
  "\t-11\t32\t35",
  "\t-12\t33\t36",
  "\t-13\t33\t36",
  "\t-14\t33\t36",
  "\t-15\t33\t36",
  "Lake Michigan",
  "\t47\t-86\t-84",
  "\t46\t-88\t-83",
  "\t45\t-89\t-83",
  "\t44\t-89\t-83",
  "\t43\t-89\t-84",
  "\t42\t-88\t-85",
  "\t41\t-88\t-85",
  "\t40\t-88\t-85",
  "Lake Okeechobee",
  "\t28\t-82\t-79",
  "\t27\t-82\t-79",
  "\t26\t-82\t-79",
  "\t25\t-82\t-79",
  "Lake Ontario",
  "\t45\t-78\t-74",
  "\t44\t-80\t-74",
  "\t43\t-80\t-74",
  "\t42\t-80\t-75",
  "Lake Pontchartrain",
  "\t31\t-91\t-88",
  "\t30\t-91\t-88",
  "\t29\t-91\t-88",
  "Lake Saint Clair",
  "\t43\t-84\t-81",
  "\t42\t-84\t-81",
  "\t41\t-84\t-81",
  "Lake Shasta",
  "\t41\t-123\t-121",
  "\t40\t-123\t-121",
  "\t39\t-123\t-121",
  "Lake Superior",
  "\t50\t-89\t-87",
  "\t49\t-90\t-84",
  "\t48\t-92\t-83",
  "\t47\t-93\t-83",
  "\t46\t-93\t-83",
  "\t45\t-93\t-83",
  "Lake Superior",
  "\t48\t-85\t-83",
  "\t47\t-86\t-83",
  "\t46\t-86\t-83",
  "\t45\t-86\t-83",
  "Lake Tahoe",
  "\t40\t-121\t-118",
  "\t39\t-121\t-118",
  "\t38\t-121\t-118",
  "\t37\t-121\t-118",
  "Lake Tanganyika",
  "\t-2\t28\t30",
  "\t-3\t28\t30",
  "\t-4\t28\t30",
  "\t-5\t28\t31",
  "\t-6\t28\t31",
  "\t-7\t28\t32",
  "\t-8\t28\t32",
  "\t-9\t29\t32",
  "Lake Victoria",
  "\t1\t30\t35",
  "\t0\t30\t35",
  "\t-1\t30\t35",
  "\t-2\t30\t35",
  "\t-3\t30\t34",
  "\t-4\t31\t33",
  "Lake Winnipeg",
  "\t55\t-99\t-96",
  "\t54\t-100\t-96",
  "\t53\t-100\t-95",
  "\t52\t-100\t-95",
  "\t51\t-99\t-95",
  "\t50\t-99\t-95",
  "\t49\t-97\t-95",
  "Laptev Sea",
  "\t82\t95\t101",
  "\t81\t95\t109",
  "\t80\t95\t118",
  "\t79\t95\t126",
  "\t78\t96\t134",
  "\t77\t101\t139",
  "\t76\t103\t139",
  "\t75\t104\t139",
  "\t74\t111\t138",
  "\t73\t111\t137",
  "\t72\t111\t137",
  "\t71\t112\t114\t126\t136",
  "\t70\t127\t134",
  "\t69\t129\t132",
  "Leyte Gulf",
  "\t12\t124\t126",
  "\t11\t124\t126",
  "\t10\t124\t126",
  "\t9\t124\t126",
  "Liddon Gulf",
  "\t76\t-115\t-110",
  "\t75\t-116\t-110",
  "\t74\t-116\t-110",
  "\t73\t-116\t-111",
  "Ligurian Sea",
  "\t45\t7\t10",
  "\t44\t6\t10",
  "\t43\t6\t10",
  "\t42\t6\t10",
  "Lincoln Sea",
  "\t84\t-70\t-36",
  "\t83\t-70\t-36",
  "\t82\t-70\t-36",
  "\t81\t-69\t-37",
  "\t80\t-54\t-48\t-46\t-43",
  "Long Island Sound",
  "\t42\t-74\t-71",
  "\t41\t-74\t-71",
  "\t40\t-74\t-71",
  "\t39\t-74\t-71",
  "Lutzow-Holm Bay",
  "\t-67\t32\t41",
  "\t-68\t32\t41",
  "\t-69\t32\t41",
  "\t-70\t32\t40",
  "\t-71\t37\t39",
  "Luzon Strait",
  "\t23\t119\t121",
  "\t22\t119\t122",
  "\t21\t119\t123",
  "\t20\t119\t123",
  "\t19\t119\t123",
  "\t18\t119\t123",
  "\t17\t119\t123",
  "M'Clure Strait",
  "\t77\t-123\t-119",
  "\t76\t-124\t-114",
  "\t75\t-125\t-113",
  "\t74\t-125\t-113",
  "\t73\t-125\t-113",
  "\t72\t-116\t-114",
  "Mackenzie Bay",
  "\t70\t-140\t-133",
  "\t69\t-140\t-133",
  "\t68\t-140\t-133",
  "\t67\t-138\t-134",
  "Makassar Strait",
  "\t2\t116\t122",
  "\t1\t116\t122",
  "\t0\t115\t122",
  "\t-1\t115\t121",
  "\t-2\t115\t120",
  "\t-3\t115\t120",
  "\t-4\t115\t120",
  "\t-5\t116\t120",
  "\t-6\t118\t120",
  "Marguerite Bay",
  "\t-66\t-70\t-65",
  "\t-67\t-71\t-65",
  "\t-68\t-71\t-65",
  "\t-69\t-71\t-65",
  "\t-70\t-71\t-65",
  "Massachusetts Bay",
  "\t43\t-72\t-69",
  "\t42\t-72\t-69",
  "\t41\t-72\t-69",
  "\t40\t-71\t-69",
  "Matochkin Shar Strait",
  "\t74\t53\t57",
  "\t73\t53\t57",
  "\t72\t53\t57",
  "McMurdo Sound",
  "\t-71\t165\t167",
  "\t-72\t163\t170",
  "\t-73\t161\t170",
  "\t-74\t159\t170",
  "\t-75\t159\t170",
  "\t-76\t159\t170",
  "\t-77\t161\t170",
  "\t-78\t161\t170",
  "\t-79\t162\t166",
  "Mecklenburger Bucht",
  "\t55\t9\t13",
  "\t54\t9\t13",
  "\t53\t9\t13",
  "\t52\t9\t12",
  "Mediterranean Sea",
  "\t38\t10\t15",
  "\t37\t9\t24\t26\t37",
  "\t36\t9\t37",
  "\t35\t9\t37",
  "\t34\t9\t36",
  "\t33\t9\t36",
  "\t32\t9\t36",
  "\t31\t11\t36",
  "\t30\t23\t35",
  "\t29\t27\t30",
  "Mediterranean Sea",
  "\t44\t4\t10",
  "\t43\t3\t10",
  "\t42\t2\t10",
  "\t41\t2\t10",
  "\t40\t2\t10",
  "\t39\t-1\t11",
  "\t38\t-2\t13",
  "\t37\t-3\t13",
  "\t36\t-3\t13",
  "\t35\t-3\t11",
  "\t34\t-2\t1",
  "Melville Bay",
  "\t77\t-68\t-59",
  "\t76\t-68\t-56",
  "\t75\t-68\t-55",
  "\t74\t-67\t-55",
  "\t73\t-62\t-55",
  "Minto Inlet",
  "\t72\t-119\t-114",
  "\t71\t-119\t-114",
  "\t70\t-119\t-114",
  "Molucca Sea",
  "\t5\t125\t127",
  "\t4\t124\t128",
  "\t3\t124\t129",
  "\t2\t123\t129",
  "\t1\t122\t129",
  "\t0\t122\t129",
  "\t-1\t122\t129",
  "\t-2\t122\t128",
  "Monterey Bay",
  "\t37\t-123\t-120",
  "\t36\t-123\t-120",
  "\t35\t-123\t-120",
  "Mozambique Channel",
  "\t-9\t39\t44",
  "\t-10\t39\t49",
  "\t-11\t39\t50",
  "\t-12\t39\t50",
  "\t-13\t39\t50",
  "\t-14\t39\t49",
  "\t-15\t38\t48",
  "\t-16\t35\t48",
  "\t-17\t34\t46",
  "\t-18\t33\t45",
  "\t-19\t33\t45",
  "\t-20\t33\t45",
  "\t-21\t33\t45",
  "\t-22\t34\t44",
  "\t-23\t33\t45",
  "\t-24\t31\t46",
  "\t-25\t31\t46",
  "\t-26\t31\t46",
  "\t-27\t31\t39",
  "Murchison Sound",
  "\t79\t-73\t-71",
  "\t78\t-73\t-65",
  "\t77\t-73\t-65",
  "\t76\t-73\t-65",
  "North Sea",
  "\t61\t-2\t7",
  "\t60\t-3\t7",
  "\t59\t-4\t8",
  "\t58\t-5\t10",
  "\t57\t-5\t10",
  "\t56\t-5\t10",
  "\t55\t-4\t10",
  "\t54\t-4\t10",
  "\t53\t-2\t10",
  "\t52\t-1\t10",
  "\t51\t-1\t5",
  "\t50\t-1\t5",
  "\t49\t0\t2",
  "Norton Sound",
  "\t65\t-165\t-159",
  "\t64\t-165\t-159",
  "\t63\t-165\t-159",
  "\t62\t-165\t-159",
  "Norwegian Sea",
  "\t77\t13\t19",
  "\t76\t9\t21",
  "\t75\t4\t23",
  "\t74\t-1\t25",
  "\t73\t-5\t27",
  "\t72\t-9\t28",
  "\t71\t-10\t28",
  "\t70\t-11\t28",
  "\t69\t-12\t28",
  "\t68\t-12\t24",
  "\t67\t-13\t18",
  "\t66\t-14\t15",
  "\t65\t-14\t15",
  "\t64\t-14\t14",
  "\t63\t-14\t12",
  "\t62\t-12\t10",
  "\t61\t-9\t9",
  "\t60\t-7\t7",
  "\t59\t-4\t6",
  "Oresund",
  "\t57\t11\t13",
  "\t56\t11\t13",
  "\t55\t11\t13",
  "\t54\t11\t13",
  "Ozero Mogotoyevo",
  "\t73\t143\t147",
  "\t72\t143\t147",
  "\t71\t143\t147",
  "Pacific Ocean",
  "\t59\t-140\t-135",
  "\t58\t-146\t-134",
  "\t57\t-153\t-132\t161\t164",
  "\t56\t-159\t-132\t161\t166",
  "\t55\t-164\t-131\t160\t168",
  "\t54\t-167\t-131\t157\t170",
  "\t53\t-172\t-130\t157\t172",
  "\t52\t-180\t-129\t156\t179",
  "\t51\t-180\t-126\t155\t180",
  "\t50\t-180\t-123\t154\t180",
  "\t49\t-180\t-123\t153\t180",
  "\t48\t-180\t-123\t152\t180",
  "\t47\t-180\t-122\t150\t180",
  "\t46\t-180\t-122\t148\t180",
  "\t45\t-180\t-122\t147\t180",
  "\t44\t-180\t-122\t143\t180",
  "\t43\t-180\t-123\t142\t180",
  "\t42\t-180\t-123\t141\t180",
  "\t41\t-180\t-123\t140\t180",
  "\t40\t-180\t-122\t140\t180",
  "\t39\t-180\t-121\t139\t180",
  "\t38\t-180\t-121\t139\t180",
  "\t37\t-180\t-120\t139\t180",
  "\t36\t-180\t-119\t138\t180",
  "\t35\t-180\t-117\t138\t180",
  "\t34\t-180\t-116\t138\t180",
  "\t33\t-180\t-115\t138\t180",
  "\t32\t-180\t-115\t138\t180",
  "\t31\t-180\t-114\t139\t180",
  "\t30\t-180\t-113\t139\t180",
  "\t29\t-180\t-113\t139\t180",
  "\t28\t-180\t-113\t140\t180",
  "\t27\t-180\t-111\t141\t180",
  "\t26\t-180\t-111\t141\t180",
  "\t25\t-180\t-110\t141\t180",
  "\t24\t-180\t-105\t140\t180",
  "\t23\t-180\t-104\t140\t180",
  "\t22\t-180\t-104\t140\t180",
  "\t21\t-180\t-104\t140\t180",
  "\t20\t-180\t-103\t140\t180",
  "\t19\t-180\t-101\t143\t180",
  "\t18\t-180\t-99\t144\t180",
  "\t17\t-180\t-97\t145\t180",
  "\t16\t-180\t-91\t145\t180",
  "\t15\t-180\t-90\t144\t180",
  "\t14\t-180\t-86\t144\t180",
  "\t13\t-180\t-85\t143\t180",
  "\t12\t-180\t-84\t141\t180",
  "\t11\t-180\t-84\t140\t180",
  "\t10\t-180\t-82\t138\t180",
  "\t9\t-180\t-80\t136\t180",
  "\t8\t-180\t-76\t135\t180",
  "\t7\t-180\t-76\t133\t180",
  "\t6\t-180\t-76\t132\t180",
  "\t5\t-180\t-76\t130\t180",
  "\t4\t-180\t-76\t128\t180",
  "\t3\t-180\t-76\t127\t180",
  "\t2\t-180\t-76\t127\t180",
  "\t1\t-180\t-76\t127\t180",
  "\t0\t-180\t-77\t128\t180",
  "\t-1\t-180\t-78\t128\t180",
  "Pacific Ocean",
  "\t4\t171\t173",
  "\t3\t-93\t-90\t170\t174",
  "\t2\t-93\t-89\t170\t174",
  "\t1\t-180\t-79\t130\t180",
  "\t0\t-180\t-79\t130\t180",
  "\t-1\t-180\t-79\t130\t180",
  "\t-2\t-180\t-79\t133\t180",
  "\t-3\t-180\t-79\t133\t142\t145\t180",
  "\t-4\t-180\t-79\t133\t136\t150\t180",
  "\t-5\t-180\t-78\t152\t180",
  "\t-6\t-180\t-78\t153\t180",
  "\t-7\t-180\t-77\t154\t180",
  "\t-8\t-180\t-77\t156\t180",
  "\t-9\t-180\t-76\t158\t180",
  "\t-10\t-180\t-76\t160\t180",
  "\t-11\t-180\t-75\t160\t180",
  "\t-12\t-180\t-75\t166\t180",
  "\t-13\t-180\t-74\t166\t180",
  "\t-14\t-180\t-73\t166\t180",
  "\t-15\t-180\t-71\t166\t180",
  "\t-16\t-180\t-69\t167\t180",
  "\t-17\t-180\t-69\t167\t180",
  "\t-18\t-180\t-69\t167\t180",
  "\t-19\t-180\t-69\t168\t180",
  "\t-20\t-180\t-69\t167\t180",
  "\t-21\t-180\t-69\t166\t180",
  "\t-22\t-180\t-69\t165\t180",
  "\t-23\t-180\t-69\t164\t180",
  "\t-24\t-180\t-69\t163\t180",
  "\t-25\t-180\t-69\t161\t180",
  "\t-26\t-180\t-69\t160\t180",
  "\t-27\t-180\t-69\t159\t180",
  "\t-28\t-180\t-69\t158\t180",
  "\t-29\t-180\t-70\t158\t180",
  "\t-30\t-180\t-70\t158\t180",
  "\t-31\t-180\t-70\t158\t180",
  "\t-32\t-180\t-70\t158\t180",
  "\t-33\t-180\t-70\t161\t180",
  "\t-34\t-180\t-70\t165\t180",
  "\t-35\t-180\t-70\t169\t180",
  "\t-36\t-180\t-71\t172\t180",
  "\t-37\t-180\t-71\t173\t180",
  "\t-38\t-180\t-72\t174\t180",
  "\t-39\t-180\t-72\t175\t180",
  "\t-40\t-180\t-72\t173\t180",
  "\t-41\t-180\t-72\t172\t180",
  "\t-42\t-180\t-72\t170\t180",
  "\t-43\t-180\t-72\t169\t180",
  "\t-44\t-180\t-72\t169\t180",
  "\t-45\t-180\t-72\t166\t180",
  "\t-46\t-180\t-72\t166\t180",
  "\t-47\t-180\t-72\t165\t180",
  "\t-48\t-180\t-72\t165\t180",
  "\t-49\t-180\t-72\t165\t180",
  "\t-50\t-180\t-71\t165\t180",
  "\t-51\t-180\t-71\t165\t180",
  "\t-52\t-180\t-71\t165\t180",
  "\t-53\t-180\t-68\t165\t180",
  "\t-54\t-180\t-67\t165\t180",
  "\t-55\t-180\t-67\t165\t180",
  "\t-56\t-180\t-67\t165\t180",
  "\t-57\t-180\t-67\t165\t180",
  "\t-58\t-180\t-67\t165\t180",
  "\t-59\t-180\t-67\t165\t180",
  "\t-60\t-180\t-67\t165\t180",
  "\t-61\t-180\t-67\t165\t180",
  "Palk Strait",
  "\t11\t78\t80",
  "\t10\t77\t81",
  "\t9\t77\t81",
  "\t8\t77\t81",
  "\t7\t78\t80",
  "Pamlico Sound",
  "\t36\t-78\t-74",
  "\t35\t-78\t-74",
  "\t34\t-78\t-74",
  "\t33\t-77\t-75",
  "Peacock Sound",
  "\t-71\t-103\t-94",
  "\t-72\t-103\t-94",
  "\t-73\t-103\t-94",
  "\t-74\t-102\t-94",
  "Persian Gulf",
  "\t31\t46\t51",
  "\t30\t46\t51",
  "\t29\t46\t52",
  "\t28\t46\t57",
  "\t27\t47\t58",
  "\t26\t47\t58",
  "\t25\t48\t58",
  "\t24\t49\t57",
  "\t23\t49\t55",
  "\t22\t50\t53",
  "Philippine Sea",
  "\t36\t135\t139",
  "\t35\t135\t140",
  "\t34\t132\t141",
  "\t33\t130\t141",
  "\t32\t129\t141",
  "\t31\t129\t141",
  "\t30\t129\t142",
  "\t29\t128\t143",
  "\t28\t127\t143",
  "\t27\t126\t143",
  "\t26\t126\t143",
  "\t25\t120\t143",
  "\t24\t120\t143",
  "\t23\t119\t143",
  "\t22\t119\t144",
  "\t21\t119\t146",
  "\t20\t119\t146",
  "\t19\t120\t147",
  "\t18\t121\t147",
  "\t17\t120\t147",
  "\t16\t120\t147",
  "\t15\t120\t147",
  "\t14\t120\t147",
  "\t13\t120\t146",
  "\t12\t121\t146",
  "\t11\t123\t145",
  "\t10\t124\t144",
  "\t9\t124\t142",
  "\t8\t125\t141",
  "\t7\t125\t139",
  "\t6\t124\t137",
  "\t5\t124\t136",
  "\t4\t124\t134",
  "\t3\t124\t133",
  "\t2\t124\t131",
  "\t1\t127\t129",
  "Porpoise Bay",
  "\t-65\t125\t131",
  "\t-66\t125\t131",
  "\t-67\t125\t131",
  "\t-68\t126\t130",
  "Prince ALbert Sound",
  "\t71\t-118\t-110",
  "\t70\t-118\t-110",
  "\t69\t-118\t-110",
  "Prince of Wales Strait",
  "\t74\t-118\t-113",
  "\t73\t-121\t-113",
  "\t72\t-121\t-113",
  "\t71\t-121\t-115",
  "\t70\t-121\t-117",
  "Prince William Sound",
  "\t62\t-149\t-145",
  "\t61\t-149\t-144",
  "\t60\t-149\t-144",
  "\t59\t-149\t-144",
  "\t58\t-148\t-146",
  "Prydz Bay",
  "\t-66\t68\t75",
  "\t-67\t68\t80",
  "\t-68\t66\t80",
  "\t-69\t66\t80",
  "\t-70\t66\t78",
  "\t-71\t65\t74",
  "\t-72\t65\t72",
  "\t-73\t65\t71",
  "\t-74\t65\t68",
  "Puget Sound",
  "\t49\t-123\t-121",
  "\t48\t-124\t-121",
  "\t47\t-124\t-121",
  "\t46\t-124\t-121",
  "Qiongzhou Strait",
  "\t21\t108\t111",
  "\t20\t108\t111",
  "\t19\t108\t111",
  "\t18\t108\t110",
  "Queen Charlotte Sound",
  "\t54\t-130\t-128",
  "\t53\t-132\t-127",
  "\t52\t-132\t-126",
  "\t51\t-132\t-126",
  "\t50\t-132\t-126",
  "\t49\t-130\t-126",
  "Queen Charlotte Straight",
  "\t52\t-128\t-125",
  "\t51\t-128\t-123",
  "\t50\t-128\t-123",
  "\t49\t-128\t-123",
  "Ragay Gulf",
  "\t14\t121\t124",
  "\t13\t121\t124",
  "\t12\t121\t124",
  "\t11\t122\t124",
  "Red Sea",
  "\t29\t33\t36",
  "\t28\t32\t36",
  "\t27\t32\t37",
  "\t26\t32\t38",
  "\t25\t32\t39",
  "\t24\t33\t39",
  "\t23\t34\t40",
  "\t22\t34\t40",
  "\t21\t34\t41",
  "\t20\t35\t42",
  "\t19\t36\t42",
  "\t18\t36\t43",
  "\t17\t36\t43",
  "\t16\t37\t43",
  "\t15\t38\t44",
  "\t14\t38\t44",
  "\t13\t39\t44",
  "\t12\t40\t44",
  "\t11\t41\t44",
  "Richard Collinson Inlet",
  "\t74\t-115\t-113",
  "\t73\t-115\t-112",
  "\t72\t-115\t-112",
  "\t71\t-115\t-112",
  "Rio de la Plata",
  "\t-31\t-59\t-57",
  "\t-32\t-59\t-57",
  "\t-33\t-59\t-53",
  "\t-34\t-59\t-53",
  "\t-35\t-59\t-53",
  "\t-36\t-58\t-54",
  "\t-37\t-58\t-55",
  "Robeson Channel",
  "\t83\t-63\t-56",
  "\t82\t-63\t-55",
  "\t81\t-63\t-55",
  "\t80\t-62\t-55",
  "Ronne Entrance",
  "\t-70\t-76\t-74",
  "\t-71\t-76\t-72",
  "\t-72\t-76\t-72",
  "\t-73\t-76\t-72",
  "Ross Sea",
  "\t-70\t169\t180",
  "\t-71\t167\t180",
  "\t-72\t167\t180",
  "\t-73\t167\t180",
  "\t-74\t168\t180",
  "\t-75\t168\t180",
  "\t-76\t165\t180",
  "\t-77\t160\t180",
  "\t-78\t158\t180",
  "\t-79\t157\t180",
  "\t-80\t157\t180",
  "\t-81\t157\t180",
  "\t-82\t159\t180",
  "\t-83\t160\t180",
  "\t-84\t166\t180",
  "\t-85\t176\t180",
  "Ross Sea",
  "\t-70\t-180\t-101",
  "\t-71\t-180\t-101",
  "\t-72\t-180\t-101",
  "\t-73\t-180\t-101",
  "\t-74\t-180\t-107",
  "\t-75\t-180\t-130\t-125\t-113",
  "\t-76\t-180\t-135",
  "\t-77\t-180\t-144",
  "\t-78\t-180\t-147",
  "\t-79\t-180\t-147",
  "\t-80\t-180\t-147",
  "\t-81\t-180\t-147",
  "\t-82\t-180\t-149",
  "\t-83\t-180\t-152",
  "\t-84\t-180\t-155",
  "\t-85\t-180\t-155",
  "\t-86\t-159\t-155",
  "Salton Sea",
  "\t34\t-117\t-114",
  "\t33\t-117\t-114",
  "\t32\t-117\t-114",
  "Samar Sea",
  "\t14\t122\t124",
  "\t13\t122\t125",
  "\t12\t122\t126",
  "\t11\t122\t126",
  "\t10\t123\t126",
  "San Francisco Bay",
  "\t39\t-123\t-120",
  "\t38\t-123\t-120",
  "\t37\t-123\t-120",
  "\t36\t-123\t-121",
  "Sargasso Sea",
  "\t36\t-68\t-51",
  "\t35\t-69\t-50",
  "\t34\t-70\t-49",
  "\t33\t-70\t-49",
  "\t32\t-71\t-49",
  "\t31\t-71\t-49",
  "\t30\t-71\t-49",
  "\t29\t-71\t-49",
  "\t28\t-71\t-49",
  "\t27\t-71\t-49",
  "\t26\t-71\t-49",
  "\t25\t-71\t-49",
  "\t24\t-71\t-49",
  "\t23\t-70\t-49",
  "\t22\t-69\t-49",
  "\t21\t-68\t-49",
  "\t20\t-67\t-49",
  "\t19\t-64\t-50",
  "Savu Sea",
  "\t-7\t117\t126",
  "\t-8\t117\t126",
  "\t-9\t117\t126",
  "\t-10\t117\t125",
  "\t-11\t119\t124",
  "Scotia Sea",
  "\t-50\t-59\t-53",
  "\t-51\t-60\t-46",
  "\t-52\t-60\t-40",
  "\t-53\t-60\t-35",
  "\t-54\t-60\t-35",
  "\t-55\t-59\t-35",
  "\t-56\t-59\t-36",
  "\t-57\t-58\t-37",
  "\t-58\t-58\t-39",
  "\t-59\t-57\t-40",
  "\t-60\t-57\t-41",
  "\t-61\t-56\t-43",
  "\t-62\t-56\t-49",
  "Sea of Azov",
  "\t48\t36\t40",
  "\t47\t33\t40",
  "\t46\t33\t40",
  "\t45\t33\t39",
  "\t44\t33\t39",
  "Sea of Crete",
  "\t39\t22\t24",
  "\t38\t21\t25",
  "\t37\t21\t29",
  "\t36\t21\t29",
  "\t35\t22\t29",
  "\t34\t22\t28",
  "Sea of Japan",
  "\t52\t139\t143",
  "\t51\t139\t143",
  "\t50\t139\t143",
  "\t49\t138\t143",
  "\t48\t137\t143",
  "\t47\t137\t143",
  "\t46\t135\t143",
  "\t45\t134\t143",
  "\t44\t130\t143",
  "\t43\t129\t142",
  "\t42\t128\t142",
  "\t41\t127\t141",
  "\t40\t127\t141",
  "\t39\t127\t141",
  "\t38\t127\t141",
  "\t37\t127\t140",
  "\t36\t126\t139",
  "\t35\t125\t138",
  "\t34\t125\t137",
  "\t33\t125\t133",
  "\t32\t125\t131",
  "\t31\t126\t130",
  "Sea of Marmara",
  "\t42\t26\t30",
  "\t41\t25\t30",
  "\t40\t25\t30",
  "\t39\t25\t30",
  "Sea of Okhotsk",
  "\t60\t141\t156",
  "\t59\t139\t156",
  "\t58\t137\t157",
  "\t57\t137\t157",
  "\t56\t137\t157",
  "\t55\t136\t157",
  "\t54\t136\t157",
  "\t53\t136\t157",
  "\t52\t136\t139\t141\t158",
  "\t51\t142\t158",
  "\t50\t142\t158",
  "\t49\t141\t157",
  "\t48\t141\t156",
  "\t47\t141\t155",
  "\t46\t140\t154",
  "\t45\t140\t153",
  "\t44\t140\t151",
  "\t43\t141\t149",
  "\t42\t143\t148",
  "Selat Bali",
  "\t-7\t113\t116",
  "\t-8\t113\t116",
  "\t-9\t113\t116",
  "Selat Dampier",
  "\t1\t128\t132",
  "\t0\t128\t132",
  "\t-1\t128\t132",
  "\t-2\t130\t132",
  "Seno de Skyring",
  "\t-51\t-74\t-70",
  "\t-52\t-74\t-70",
  "\t-53\t-74\t-70",
  "\t-54\t-74\t-71",
  "Seno Otway",
  "\t-51\t-72\t-70",
  "\t-52\t-73\t-70",
  "\t-53\t-73\t-70",
  "\t-54\t-73\t-70",
  "Shark Bay",
  "\t-23\t112\t114",
  "\t-24\t112\t115",
  "\t-25\t112\t115",
  "\t-26\t112\t115",
  "\t-27\t112\t115",
  "Shelikhova Gulf",
  "\t63\t162\t166",
  "\t62\t155\t166",
  "\t61\t153\t166",
  "\t60\t153\t165",
  "\t59\t153\t164",
  "\t58\t153\t162",
  "\t57\t154\t160",
  "\t56\t155\t158",
  "Sherman Basin",
  "\t69\t-99\t-97",
  "\t68\t-99\t-96",
  "\t67\t-99\t-96",
  "\t66\t-99\t-96",
  "Sibuyan Sea",
  "\t14\t120\t123",
  "\t13\t120\t124",
  "\t12\t120\t124",
  "\t11\t120\t124",
  "\t10\t120\t124",
  "Skagerrak",
  "\t60\t8\t12",
  "\t59\t6\t12",
  "\t58\t6\t12",
  "\t57\t6\t12",
  "\t56\t6\t11",
  "\t55\t7\t9",
  "Smith Sound",
  "\t53\t-128\t-125",
  "\t52\t-129\t-125",
  "\t51\t-129\t-125",
  "\t50\t-129\t-125",
  "Sognefjorden",
  "\t62\t3\t8",
  "\t61\t3\t8",
  "\t60\t3\t8",
  "\t59\t4\t8",
  "Solomon Sea",
  "\t-3\t151\t155",
  "\t-4\t146\t155",
  "\t-5\t145\t157",
  "\t-6\t145\t160",
  "\t-7\t145\t161",
  "\t-8\t146\t162",
  "\t-9\t147\t163",
  "\t-10\t147\t163",
  "\t-11\t148\t163",
  "\t-12\t152\t162",
  "South China Sea",
  "\t24\t112\t121",
  "\t23\t112\t121",
  "\t22\t109\t122",
  "\t21\t108\t123",
  "\t20\t108\t123",
  "\t19\t107\t123",
  "\t18\t105\t123",
  "\t17\t105\t123",
  "\t16\t105\t121",
  "\t15\t106\t121",
  "\t14\t107\t121",
  "\t13\t108\t121",
  "\t12\t107\t121",
  "\t11\t104\t121",
  "\t10\t104\t120",
  "\t9\t103\t120",
  "\t8\t102\t119",
  "\t7\t101\t118",
  "\t6\t101\t117",
  "\t5\t101\t117",
  "\t4\t101\t117",
  "\t3\t102\t116",
  "\t2\t102\t114",
  "\t1\t101\t113",
  "\t0\t101\t112",
  "\t-1\t101\t111",
  "\t-2\t103\t111",
  "\t-3\t103\t111",
  "\t-4\t105\t107",
  "Southern Ocean",
  "\t-59\t-180\t180",
  "\t-60\t-180\t180",
  "\t-61\t-180\t180",
  "\t-62\t-180\t180",
  "\t-63\t-180\t180",
  "\t-64\t-180\t180",
  "\t-65\t-180\t91\t103\t180",
  "\t-66\t-180\t87\t110\t180",
  "\t-67\t-180\t52\t54\t85\t112\t180",
  "\t-68\t-180\t51\t54\t84\t112\t122\t141\t180",
  "\t-69\t-180\t-69\t-66\t44\t74\t80\t145\t180",
  "\t-70\t-180\t-73\t-64\t33\t154\t180",
  "\t-71\t-180\t-79\t-63\t33\t159\t180",
  "\t-72\t-180\t-85\t-62\t2\t24\t27\t161\t163\t166\t180",
  "\t-73\t-100\t-91",
  "St. Helena Bay",
  "\t-30\t16\t19",
  "\t-31\t16\t19",
  "\t-32\t16\t19",
  "\t-33\t16\t19",
  "St. Lawrence River",
  "\t51\t-67\t-63",
  "\t50\t-69\t-63",
  "\t49\t-72\t-63",
  "\t48\t-72\t-63",
  "\t47\t-74\t-66",
  "\t46\t-75\t-68",
  "\t45\t-75\t-69",
  "\t44\t-75\t-72",
  "Stettiner Haff",
  "\t55\t12\t14",
  "\t54\t12\t15",
  "\t53\t12\t15",
  "\t52\t12\t15",
  "Storfjorden",
  "\t79\t17\t22",
  "\t78\t16\t22",
  "\t77\t15\t22",
  "\t76\t15\t22",
  "\t75\t15\t19",
  "Strait of Belle Isle",
  "\t53\t-56\t-54",
  "\t52\t-58\t-54",
  "\t51\t-58\t-54",
  "\t50\t-58\t-54",
  "Strait of Georgia",
  "\t51\t-126\t-122",
  "\t50\t-126\t-121",
  "\t49\t-126\t-121",
  "\t48\t-126\t-121",
  "\t47\t-124\t-121",
  "Strait of Gibraltar",
  "\t37\t-7\t-4",
  "\t36\t-7\t-4",
  "\t35\t-7\t-4",
  "\t34\t-6\t-4",
  "Strait of Juan de Fuca",
  "\t49\t-125\t-121",
  "\t48\t-125\t-121",
  "\t47\t-125\t-121",
  "Strait of Malacca",
  "\t9\t97\t99",
  "\t8\t97\t100",
  "\t7\t95\t101",
  "\t6\t94\t101",
  "\t5\t94\t101",
  "\t4\t94\t102",
  "\t3\t96\t103",
  "\t2\t97\t104",
  "\t1\t98\t104",
  "\t0\t99\t104",
  "\t-1\t101\t104",
  "Strait of Singapore",
  "\t2\t102\t105",
  "\t1\t102\t105",
  "\t0\t102\t105",
  "Straits of Florida",
  "\t27\t-81\t-77",
  "\t26\t-82\t-77",
  "\t25\t-84\t-77",
  "\t24\t-84\t-77",
  "\t23\t-84\t-77",
  "\t22\t-84\t-78",
  "Sulu Sea",
  "\t13\t118\t122",
  "\t12\t118\t123",
  "\t11\t118\t123",
  "\t10\t117\t124",
  "\t9\t116\t124",
  "\t8\t115\t124",
  "\t7\t115\t124",
  "\t6\t115\t123",
  "\t5\t115\t123",
  "\t4\t116\t121",
  "Sulzberger Bay",
  "\t-75\t-153\t-144",
  "\t-76\t-159\t-144",
  "\t-77\t-159\t-144",
  "\t-78\t-159\t-144",
  "Surigao Strait",
  "\t11\t124\t126",
  "\t10\t124\t126",
  "\t9\t124\t126",
  "\t8\t124\t126",
  "Taiwan Strait",
  "\t26\t117\t122",
  "\t25\t116\t122",
  "\t24\t116\t122",
  "\t23\t116\t121",
  "\t22\t116\t121",
  "Tasman Sea",
  "\t-28\t152\t160",
  "\t-29\t152\t160",
  "\t-30\t151\t162",
  "\t-31\t150\t166",
  "\t-32\t150\t170",
  "\t-33\t149\t174",
  "\t-34\t149\t174",
  "\t-35\t148\t175",
  "\t-36\t148\t175",
  "\t-37\t147\t175",
  "\t-38\t146\t176",
  "\t-39\t146\t176",
  "\t-40\t146\t176",
  "\t-41\t146\t176",
  "\t-42\t145\t175",
  "\t-43\t145\t172",
  "\t-44\t145\t171",
  "\t-45\t147\t169",
  "\t-46\t150\t168",
  "\t-47\t152\t168",
  "\t-48\t155\t168",
  "\t-49\t158\t167",
  "\t-50\t160\t167",
  "\t-51\t163\t167",
  "Tatar Strait",
  "\t54\t139\t142",
  "\t53\t139\t142",
  "\t52\t139\t142",
  "\t51\t140\t142",
  "\t50\t140\t142",
  "Tayabas Bay",
  "\t14\t119\t123",
  "\t13\t119\t123",
  "\t12\t119\t123",
  "The North Western Passages",
  "\t81\t-101\t-95",
  "\t80\t-108\t-90\t-88\t-82",
  "\t79\t-114\t-80",
  "\t78\t-117\t-80",
  "\t77\t-120\t-80",
  "\t76\t-120\t-81",
  "\t75\t-120\t-78",
  "\t74\t-120\t-76",
  "\t73\t-106\t-76",
  "\t72\t-106\t-76",
  "\t71\t-118\t-116\t-106\t-83",
  "\t70\t-119\t-112\t-108\t-83",
  "\t69\t-119\t-85",
  "\t68\t-119\t-92",
  "\t67\t-118\t-92",
  "\t66\t-116\t-106\t-104\t-94",
  "\t65\t-97\t-94",
  "Timor Sea",
  "\t-7\t125\t131",
  "\t-8\t123\t131",
  "\t-9\t121\t132",
  "\t-10\t121\t133",
  "\t-11\t121\t133",
  "\t-12\t122\t133",
  "\t-13\t124\t133",
  "\t-14\t125\t131",
  "Torres Strait",
  "\t-8\t140\t144",
  "\t-9\t140\t144",
  "\t-10\t140\t144",
  "\t-11\t141\t143",
  "\t-12\t141\t143",
  "Trondheimsfjorden",
  "\t65\t10\t12",
  "\t64\t7\t12",
  "\t63\t7\t12",
  "\t62\t7\t12",
  "Tsugaru Strait",
  "\t42\t139\t142",
  "\t41\t139\t142",
  "\t40\t139\t142",
  "\t39\t139\t142",
  "Tyrrhenian Sea",
  "\t45\t8\t11",
  "\t44\t8\t11",
  "\t43\t8\t12",
  "\t42\t8\t14",
  "\t41\t8\t16",
  "\t40\t8\t17",
  "\t39\t7\t17",
  "\t38\t7\t17",
  "\t37\t7\t17",
  "\t36\t10\t14",
  "Uchiura Bay",
  "\t43\t139\t144",
  "\t42\t139\t144",
  "\t41\t139\t144",
  "\t40\t139\t143",
  "\t39\t140\t142",
  "Uda Bay",
  "\t57\t136\t139",
  "\t56\t134\t139",
  "\t55\t134\t139",
  "\t54\t134\t139",
  "\t53\t134\t139",
  "\t52\t135\t138",
  "Ungava Bay",
  "\t61\t-71\t-63",
  "\t60\t-71\t-63",
  "\t59\t-71\t-63",
  "\t58\t-71\t-64",
  "\t57\t-71\t-64",
  "\t56\t-70\t-66",
  "Uummannaq Fjord",
  "\t73\t-54\t-52",
  "\t72\t-55\t-50",
  "\t71\t-55\t-49",
  "\t70\t-55\t-49",
  "\t69\t-55\t-49",
  "Vestfjorden",
  "\t69\t12\t18",
  "\t68\t11\t18",
  "\t67\t11\t18",
  "\t66\t11\t17",
  "\t65\t12\t14",
  "Vil'kitskogo Strait",
  "\t79\t99\t106",
  "\t78\t99\t106",
  "\t77\t99\t106",
  "\t76\t99\t106",
  "\t75\t99\t101",
  "Vincennes Bay",
  "\t-65\t103\t111",
  "\t-66\t103\t111",
  "\t-67\t103\t111",
  "Visayan Sea",
  "\t13\t122\t124",
  "\t12\t121\t125",
  "\t11\t121\t125",
  "\t10\t121\t125",
  "\t9\t121\t125",
  "Viscount Melville Sound",
  "\t76\t-110\t-103",
  "\t75\t-115\t-103",
  "\t74\t-116\t-103",
  "\t73\t-116\t-103",
  "\t72\t-116\t-104",
  "\t71\t-114\t-107",
  "Waddenzee",
  "\t54\t3\t7",
  "\t53\t3\t7",
  "\t52\t3\t7",
  "\t51\t3\t6",
  "Wager Bay",
  "\t66\t-92\t-86",
  "\t65\t-92\t-86",
  "\t64\t-92\t-86",
  "Weddell Sea",
  "\t-70\t-62\t-9",
  "\t-71\t-63\t-9",
  "\t-72\t-63\t-9",
  "\t-73\t-64\t-10",
  "\t-74\t-66\t-13",
  "\t-75\t-78\t-14",
  "\t-76\t-84\t-17",
  "\t-77\t-84\t-25",
  "\t-78\t-84\t-22",
  "\t-79\t-84\t-22",
  "\t-80\t-82\t-22",
  "\t-81\t-79\t-23",
  "\t-82\t-70\t-36",
  "\t-83\t-66\t-50\t-48\t-42",
  "\t-84\t-62\t-57",
  "White Sea",
  "\t69\t37\t45",
  "\t68\t30\t33\t37\t45",
  "\t67\t30\t45",
  "\t66\t30\t45",
  "\t65\t31\t45",
  "\t64\t33\t41",
  "\t63\t33\t41",
  "\t62\t35\t38",
  "Wrigley Gulf",
  "\t-72\t-131\t-124",
  "\t-73\t-135\t-123",
  "\t-74\t-135\t-123",
  "\t-75\t-135\t-123",
  "Wynniat Bay",
  "\t73\t-112\t-109",
  "\t72\t-112\t-109",
  "\t71\t-112\t-109",
  "Yellow Sea",
  "\t41\t123\t125",
  "\t40\t120\t126",
  "\t39\t120\t126",
  "\t38\t119\t127",
  "\t37\t119\t127",
  "\t36\t118\t127",
  "\t35\t118\t127",
  "\t34\t118\t127",
  "\t33\t118\t127",
  "\t32\t119\t127",
  "\t31\t119\t125",
  "\t30\t120\t123",
  "Yellowstone Lake",
  "\t45\t-111\t-109",
  "\t44\t-111\t-109",
  "\t43\t-111\t-109",
  "Yenisey Gulf",
  "\t74\t77\t81",
  "\t73\t77\t83",
  "\t72\t77\t84",
  "\t71\t77\t84",
  "\t70\t79\t84",
  "\t69\t81\t84",
  "Yucatan Channel",
  "\t23\t-86\t-84",
  "\t22\t-88\t-83",
  "\t21\t-88\t-83",
  "\t20\t-88\t-83",
  NULL
};


static Uint4 sqn_binary_search_on_uint4_list(Uint4Ptr list, Uint4 pos, Uint4 listlen)
{
   Uint4  L;
   Uint4  mid;
   Uint4  R;

   if (list == NULL || listlen == 0)
      return 0;
   L = 0;
   R = listlen - 1;
   while (L < R)
   {
      mid = (L+R)/2;
      if (list[mid + 1] <= pos)
      {
         L = mid + 1;
      } else
      {
         R = mid;
      }
   }
   return R;
}


static Int4 MapRowCoordsSpecial(SeqAlignPtr sap, Uint4 pos, Int4 row, Boolean is_left_end)
{
   DenseSegPtr  dsp;
   Int4         idx;
   Int4         offset;
   SAIndexPtr   saip;
   Int4         start;

   if (sap == NULL || row < 0)
      return -1;
   if (sap->saip == NULL)
      return -1;
   saip = (SAIndexPtr)sap->saip;
   dsp = (DenseSegPtr)sap->segs;
   start = sqn_binary_search_on_uint4_list(saip->aligncoords, pos, dsp->numseg);
   offset = pos - saip->aligncoords[start];
   idx = (dsp->dim*start) + row - 1;
   if (dsp->starts[idx] == -1)
   {
      if (is_left_end)
      {
         /* round up */
         while (start < dsp->numseg) {
            idx = (dsp->dim*start) + row - 1;
            if (dsp->starts[idx] != -1)
               return (dsp->starts[idx]);
            start++;
         }
         return -2;
      }
      else
      {
         /* round down */
         while (start >= 0) {
            idx = (dsp->dim*start) + row - 1;
            if (dsp->starts[idx] != -1)
               return (dsp->starts[idx] + dsp->lens[start] - 1);
            start--;
         }
         return -2;
      }
   } else
   {
      idx = (dsp->dim*start) + row - 1;
      if (dsp->strands[idx] != Seq_strand_minus)
         return (dsp->starts[idx] + offset);
      else
         return (dsp->starts[idx] + dsp->lens[start] - 1 - offset);
   }
   return -1;
}


static Int4 MapBioseqToBioseqSpecial(SeqAlignPtr sap, Int4 begin, Int4 fin, Int4 pos, Boolean is_left_end)
{
   Int4  bspos;
   Int4  sapos;
   Int4  start1;
   Int4  start2;
   Int4  stop1;
   Int4  stop2;

   if (sap == NULL || sap->saip == NULL)
      return -2;
   AlnMgr2GetNthSeqRangeInSA(sap, begin, &start1, &stop1);
   AlnMgr2GetNthSeqRangeInSA(sap, fin, &start2, &stop2);
   /* check to see whether the position is outside the alignment */
   if (pos < start1)
      return (start2 - (start1 - pos));
   else if (pos > stop1)
      return (stop2 + (pos-stop1));
   sapos = AlnMgr2MapBioseqToSeqAlign(sap, pos, begin);
   bspos = MapRowCoordsSpecial(sap, sapos, fin, is_left_end);
   if (bspos >= 0)
      return bspos;
   else if (is_left_end)
      return (start2-1);
   else
      return (stop2+1);
}


/* This function adjusts the endpoints of a location, as long as the
 * endpoints are in the area represented by the alignment.
 * When we are adjusting locations for an alignment of a part, we will 
 * be looking at all features indexed on the main segment, but we only
 * want to adjust feature endpoints located on the segment that we are
 * updating.
 */
static Int4 AdjustEndpoint 
(SeqAlignPtr salp,
 SeqLocPtr   slp, 
 Int4        max_length, 
 Int4        begin, 
 Int4        fin,
 Int4        endpoint,
 Boolean     is_left_end)
{
  BioseqPtr            slp_bsp, parent_bsp, old_bsp;
  SeqMgrSegmentContext segcontext; 
  SeqIdPtr             old_sip; 
  Int4                 pt;
  
  if (slp == NULL || salp == NULL)
  {
    return endpoint;
  }
  
  old_sip = AlnMgr2GetNthSeqIdPtr (salp, begin);
  old_bsp = BioseqFind (old_sip);

  parent_bsp = SeqMgrGetParentOfPart (old_bsp, &segcontext);
  
  slp_bsp = BioseqFind (SeqLocId (slp));
  if (slp_bsp == old_bsp
      || (slp_bsp == parent_bsp 
          && endpoint >= segcontext.cumOffset + segcontext.from
          && endpoint < segcontext.cumOffset + segcontext.to))
  {
    if (slp_bsp == parent_bsp)
    {
      endpoint -= segcontext.cumOffset + segcontext.from;
    }
    pt = MapBioseqToBioseqSpecial (salp, begin, fin, endpoint, is_left_end);
    if (pt < 0) {
      pt = 0;
    } else if (pt >= max_length) {
      pt = max_length - 1;
    }
    if (slp_bsp == parent_bsp)
    {
      pt += segcontext.cumOffset + segcontext.from;
    }
  }
  else
  {
    pt = endpoint;
  }
  
  return pt;  
}


static void ReplaceLocation (SeqAlignPtr salp, SeqLocPtr slp, Int4 length, Int4 begin, Int4 fin)

{
  PackSeqPntPtr  psp;
  SeqIntPtr      sinp;
  SeqPntPtr      spp;
  Uint1          used;

  if (slp == NULL) return;
  switch (slp->choice) {
    case SEQLOC_INT :
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sinp->from = AdjustEndpoint (salp, slp, length, begin, fin,
                                     sinp->from, TRUE);
        sinp->to = AdjustEndpoint (salp, slp, length, begin, fin,
                                   sinp->to, FALSE);
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        spp->point = AdjustEndpoint (salp, slp, length, begin, fin,
                                     spp->point, TRUE);
      }
      break;
    case SEQLOC_PACKED_PNT :
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        for (used = 0; used < psp->used; used++) {
          psp->pnts [used] = AdjustEndpoint (salp, slp, length, begin, fin,
                                             psp->pnts [used], TRUE);
        }
      }
      break;
    default :
      break;
  }
}


/* this function iterates through the pieces of a complex location
 * and calls ReplaceLocation for each one.  ReplaceLocation will only
 * act on SEQLOC_INT, SEQLOC_PNT, and SEQLOC_PACKED_PNT and will ignore
 * other types.
 */
NLM_EXTERN void 
ReplaceComplexLocation 
(SeqLocPtr   slp,
 SeqAlignPtr salp,
 Int4        new_len,
 Int4        begin,
 Int4        fin)
{
  SeqLocPtr subslp;
  
  if (slp == NULL || salp == NULL)
  {
    return;
  }
  
  subslp = SeqLocFindNext (slp, NULL);
  while (subslp != NULL) {
    ReplaceLocation (salp, subslp, new_len, begin, fin);
    subslp = SeqLocFindNext (slp, subslp);
  }
}


static void UpdateOneFeatureForSequenceReplace 
(SeqFeatPtr  sfp, 
 SeqAlignPtr salp,
 BioseqPtr   oldbsp,
 Int4        new_len)
{
  CodeBreakPtr cbp;
  CdRegionPtr  crp;
  RnaRefPtr    rrp;
  tRNAPtr      trp;
  
  if (sfp == NULL || salp == NULL)
  {
    return;
  }
  
  ReplaceComplexLocation (sfp->location, salp, new_len, 1, 2);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) 
      {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) 
        {
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


static void UpdateLocationsForSequenceReplace 
(SeqAlignPtr salp,
 BioseqPtr oldbsp,
 BioseqPtr newbsp)
{
  BioseqPtr    parentbsp;
  SeqMgrFeatContext  context;
  SeqMgrSegmentContext segcontext;
  SeqFeatPtr           sfp;
  
  if (salp == NULL || oldbsp == NULL || newbsp == NULL)
  {
    return;
  }
    
  /* if this sequence is a part, the features will be indexed on
   * the parent.
   */
  parentbsp = SeqMgrGetParentOfPart (oldbsp, &segcontext);
  if (parentbsp == NULL)
  {
    sfp = SeqMgrGetNextFeature (oldbsp, NULL, 0, 0, &context);
    while (sfp != NULL)
    {
      UpdateOneFeatureForSequenceReplace (sfp, salp, oldbsp, newbsp->length);
      sfp = SeqMgrGetNextFeature (oldbsp, sfp, 0, 0, &context);
                                  
    }
  }
  else
  {
    sfp = SeqMgrGetNextFeature (parentbsp, NULL, 0, 0, &context);
    while (sfp != NULL)
    {
      UpdateOneFeatureForSequenceReplace (sfp, salp, oldbsp, newbsp->length);
      sfp = SeqMgrGetNextFeature (parentbsp, sfp, 0, 0, &context);
    }
  }  
}


NLM_EXTERN void 
ReplaceOneSequence 
(SeqAlignPtr salp,
 BioseqPtr oldbsp,
 BioseqPtr newbsp)
{
  SeqDataPtr         bs;
  Int4               len, len_change;
  Uint1              seq_data_type, seq_ext_type;
  Pointer            seq_ext;
  Uint1              repr;
  BioseqPtr          parent_bsp;
  SeqMgrSegmentContext context;
  
  if (oldbsp == NULL || newbsp == NULL)
  {
    return;
  }

  UpdateLocationsForSequenceReplace (salp, oldbsp, newbsp);
  len_change = newbsp->length - oldbsp->length;

  /* switch bioseqs to finish update */

  bs = oldbsp->seq_data;
  oldbsp->seq_data = newbsp->seq_data;
  newbsp->seq_data = bs;
  len = oldbsp->length;
  oldbsp->length = newbsp->length;
  newbsp->length = len;
  seq_data_type = oldbsp->seq_data_type;
  oldbsp->seq_data_type = newbsp->seq_data_type;
  newbsp->seq_data_type = seq_data_type;  
  /* also move seq_ext, for delta sequences */
  seq_ext_type = oldbsp->seq_ext_type;
  seq_ext = oldbsp->seq_ext;
  oldbsp->seq_ext_type = newbsp->seq_ext_type;
  oldbsp->seq_ext = newbsp->seq_ext;
  newbsp->seq_ext_type = seq_ext_type;
  newbsp->seq_ext = seq_ext;
  
  /* swap repr */
  repr = oldbsp->repr;
  oldbsp->repr = newbsp->repr;
  newbsp->repr = repr;
  
  /* if this was part of a segmented set, update the parent length */
  parent_bsp = SeqMgrGetParentOfPart (oldbsp, &context);
  if (parent_bsp != NULL)
  {
    parent_bsp->length += len_change;
  }
}


NLM_EXTERN Boolean AreSequenceResiduesIdentical (BioseqPtr bsp1, BioseqPtr bsp2)
{
  SeqPortPtr    spp1, spp2;
  Uint1         seqcode;
  Int4          buf_len = 255;
  Char          buf1[255], buf2[255];
  Int4          ctr1, ctr2, offset;
  Boolean       rval;
  
  if (bsp1 == NULL && bsp2 == NULL)
  {
    return TRUE;
  }
  else if (bsp1 == NULL || bsp2 == NULL)
  {
    return FALSE;
  }
  else if (bsp1->length != bsp2->length)
  {
    return FALSE;
  }
  else if (ISA_na (bsp1->mol) && ! ISA_na (bsp2->mol))
  {
    return FALSE;
  }
  else if (!ISA_na (bsp1->mol) && ISA_na (bsp2->mol))
  {
    return FALSE;
  }

  if (ISA_na (bsp1->mol))
  {
    seqcode = Seq_code_iupacna;
  }
  else
  {
    seqcode = Seq_code_iupacaa;
  }

  
  spp1 = SeqPortNew (bsp1, 0, bsp1->length - 1, Seq_strand_plus, seqcode);
  spp2 = SeqPortNew (bsp2, 0, bsp2->length - 1, Seq_strand_plus, seqcode);
  
  ctr1 = SeqPortRead (spp1, (Uint1Ptr)buf1, buf_len - 1);
  ctr2 = SeqPortRead (spp2, (Uint1Ptr)buf2, buf_len - 1);
  buf1 [ctr1] = 0;
  buf2 [ctr2] = 0;
  offset = ctr1;

  while (ctr1 == ctr2 && StringCmp (buf1, buf2) == 0 && offset < bsp1->length)
  {
    ctr1 = SeqPortRead (spp1, (Uint1Ptr)buf1, buf_len - 1);
    ctr2 = SeqPortRead (spp2, (Uint1Ptr)buf2, buf_len - 1);
    buf1 [ctr1] = 0;
    buf2 [ctr2] = 0;
    offset += ctr1;
  }
  
  if (ctr1 != ctr2 || StringCmp (buf1, buf2) != 0 || offset < bsp1->length)
  {
    rval = FALSE;
  }
  else
  {
    rval = TRUE;
  }
  
  spp1 = SeqPortFree (spp1);
  spp2 = SeqPortFree (spp2);
  
  return rval;
}


static Boolean FindBestCitSubCallback (GatherContextPtr gcp)

{
  CitSubPtr   best;
  CitSubPtr   PNTR bestp;
  CitSubPtr   csp;
  PubdescPtr  pdp;
  ValNodePtr  sdp;
  ValNodePtr  vnp;

  if (gcp == NULL) return TRUE;
  bestp = (CitSubPtr PNTR) gcp->userdata;
  if (bestp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQDESC) return TRUE;
  sdp = (ValNodePtr) gcp->thisitem;
  if (sdp == NULL || sdp->choice != Seq_descr_pub) return TRUE;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return TRUE;
  vnp = pdp->pub;
  if (vnp == NULL || vnp->choice != PUB_Sub) return TRUE;
  csp = (CitSubPtr) vnp->data.ptrvalue;
  if (csp == NULL) return TRUE;
  if (*bestp == NULL) {
    *bestp = csp;
    return TRUE;
  }
  best = *bestp;
  if (DateMatch (best->date, csp->date, FALSE) == -1) {
    *bestp = csp;
    return TRUE;
  }
  return TRUE;
}


NLM_EXTERN ValNodePtr CreateUpdateCitSubFromBestTemplate (
  SeqEntryPtr top_sep,
  SeqEntryPtr upd_sep,
  CharPtr     update_txt
)
{
  CitSubPtr    best;
  CitSubPtr    csp;
  DatePtr      dp;
  GatherScope  gs;
  PubdescPtr   pdp;
  ValNodePtr   sdp;
  ValNodePtr   vnp;
  ValNodePtr   err_list = NULL;

  best = NULL;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.scope = top_sep;
  GatherSeqEntry (top_sep, (Pointer) (&best), FindBestCitSubCallback, &gs);
  if (best == NULL) {
    ValNodeAddPointer (&err_list, 0, "There is no earlier cit-sub template");
    return err_list;
  }
  dp = DateCurr ();
  if (dp != NULL) {
    if (DateMatch (best->date, dp, FALSE) == 0) {
      DateFree (dp);
	    dp = NULL;
      if (upd_sep == top_sep
          && StringICmp (best->descr, update_txt) == 0)
      {
        ValNodeAddPointer (&err_list, 0, "There already exists an update on today's date");
        return err_list;
      } else if (best->descr == NULL) {
        best->descr = StringSave (update_txt);
        ValNodeAddPointer (&err_list, 1, "Adding update indication to existing cit-sub");
        return err_list;
      }
    }
    DateFree (dp);
  }
  sdp = CreateNewDescriptor (upd_sep, Seq_descr_pub);
  if (sdp == NULL) return FALSE;
  pdp = PubdescNew ();
  if (pdp == NULL) return FALSE;
  sdp->data.ptrvalue = (Pointer) pdp;
  vnp = ValNodeNew (NULL);
  csp = AsnIoMemCopy ((Pointer) best,
                      (AsnReadFunc) CitSubAsnRead,
                      (AsnWriteFunc) CitSubAsnWrite);
  pdp->pub = vnp;
  vnp->choice = PUB_Sub;
  vnp->data.ptrvalue = csp;
  csp->date = DateFree (csp->date);
  csp->date = DateCurr ();
  if (!StringHasNoText (update_txt)) {
    csp->descr = MemFree (csp->descr);
    csp->descr = StringSave (update_txt);
  }

  if (top_sep == upd_sep)
  {
    ValNodeAddPointer (&err_list, 1, "The update Cit-sub has been placed on the top Seq-entry");
  }
  return err_list;
}                               


CharPtr kSubmitterUpdateText = "Sequence update by submitter";

NLM_EXTERN void AddCitSubToUpdatedSequence (BioseqPtr upd_bsp, Uint2 input_entityID, CharPtr update_txt)
{
  SeqEntryPtr top_sep, upd_sep;
  ValNodePtr  msg_list;

  upd_sep = SeqMgrGetSeqEntryForData (upd_bsp);
  if (upd_sep == NULL) return;
  top_sep = GetTopSeqEntryForEntityID ( input_entityID);
  if (top_sep == NULL) return;
  msg_list = CreateUpdateCitSubFromBestTemplate (top_sep, upd_sep, update_txt);
  if (msg_list != NULL) {
    if (msg_list->choice == 0) {
      Message (MSG_ERROR, (CharPtr) msg_list->data.ptrvalue);
    } else {
      Message (MSG_OK, (CharPtr) msg_list->data.ptrvalue);
    }
    msg_list = ValNodeFree (msg_list);
  }
}


static void ListPhrapGraphsCallback (SeqGraphPtr sgp, Pointer userdata)
{
  ValNodePtr PNTR vnpp;
  
  if (sgp == NULL || userdata == NULL) return;
  if (StringICmp (sgp->title, "Phrap Quality") == 0)
  {
    vnpp = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (vnpp, 0, sgp);
  }
}


NLM_EXTERN void RemoveQualityScores 
(BioseqPtr bsp, 
 FILE      *log_fp,
 BoolPtr   data_in_log)

{
  ValNodePtr score_list = NULL, vnp;
  Char          acc_str [256];
  SeqGraphPtr   sgp;

  if (bsp == NULL) return;

  VisitGraphsOnBsp (bsp, &score_list, ListPhrapGraphsCallback);
  if (score_list == NULL) return;

  if (log_fp != NULL && data_in_log != NULL)
  {
    SeqIdWrite (bsp->id, acc_str, PRINTID_REPORT, sizeof (acc_str));
    fprintf (log_fp, "Quality scores cleared for %s\n", acc_str);
    *data_in_log = TRUE;
  }

  for (vnp = score_list; vnp != NULL; vnp = vnp->next) 
  {
    sgp = vnp->data.ptrvalue;
    if (sgp != NULL)
    {
      sgp->idx.deleteme = TRUE;
    }
  }

  score_list = ValNodeFree (score_list);
  DeleteMarkedObjects (0, OBJ_BIOSEQ, (Pointer) bsp);
}


static Char GetNextCharacterFromFile (FILE *fp, BoolPtr pIsASN)
{
  FileCache    fc;
  CharPtr      str;
  Char         special_symbol;
  Char         line [128];
  Int4         pos;
  
  /* look ahead to see what character caused inability to interpret line */
  FileCacheSetup (&fc, fp);
  /* pos = FileCacheTell (&fc); */
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  if (str != NULL && StringDoesHaveText (str)) {
    TrimSpacesAroundString (str);
  }
  special_symbol = line [0];
  if (pIsASN != NULL)
  {
    if (StringStr (line, "::=") != NULL)
    {
      *pIsASN = TRUE;
    }
    else
    {
      *pIsASN = FALSE;
    }
  }
  /* seek to start of next line after one that could not be interpreted */
  pos = FileCacheTell (&fc);
  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, pos);
  fseek (fp, pos, SEEK_SET);
  return special_symbol;
}


NLM_EXTERN void ReplaceFakeIDWithIDFromTitle (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  CharPtr     title_txt, new_id_str;
  Int4        id_len;
  Boolean     remove_punct = FALSE;
  
  if (bsp == NULL)
  {
    return;
  }
  
  bsp->id = SeqIdFree (bsp->id);
  
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp != NULL && !StringHasNoText (sdp->data.ptrvalue))
  {
    title_txt = sdp->data.ptrvalue;
    /* skip any leading spaces */
    title_txt += StringSpn (title_txt, " \t");
    /* look for local IDs surrounded by quotes - the real way to have an ID with a space */
    if (*title_txt == '\'')
    {
      title_txt ++;
      id_len = StringCSpn (title_txt, "\'");
      if (title_txt [id_len] == '\'')
      {
        remove_punct = TRUE;
      }
      else
      {
        id_len = 0;
      }
    }
    else if (*title_txt == '\"')
    {
      title_txt ++;
      id_len = StringCSpn (title_txt, "\"");
      if (title_txt [id_len] == '\"')
      {
        remove_punct = TRUE;
      }
      else
      {
        id_len = 0;
      }
    }
    else
    {
      id_len = StringCSpn (title_txt, " \t");
    }
    
    if (id_len > 0)
    {
      new_id_str = (CharPtr) MemNew ((id_len + 1) * sizeof (Char));
      if (new_id_str != NULL)
      {
        StringNCpy (new_id_str, title_txt, id_len);
        new_id_str [id_len] = 0;
        bsp->id = MakeSeqID (new_id_str);
        new_id_str = MemFree (new_id_str);
        /* remove id from title */
        title_txt += id_len;
        if (remove_punct)
        {
          title_txt ++;
        }
        title_txt += StringSpn (title_txt, " \t");
        title_txt = StringSave (title_txt);
        sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
        sdp->data.ptrvalue = title_txt;
      }
    }
  }  
}


static void PutDeflineIDBackInTitle (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  CharPtr     id_txt;
  CharPtr     title_txt;
  
  if (bsp == NULL || bsp->id == NULL)
  {
    return;
  }

  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_title)
  {
    sdp = sdp->next;
  }
  if (sdp == NULL)
  {
    sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_title);
  }
  if (sdp == NULL)
  {
    return;
  }
  
  id_txt = SeqIdWholeLabel (bsp->id, PRINTID_REPORT);
  
  if (StringHasNoText (sdp->data.ptrvalue))
  {
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = id_txt;
  }
  else
  {
    title_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (id_txt) + StringLen (sdp->data.ptrvalue) + 2));
    StringCpy (title_txt, id_txt);
    StringCat (title_txt, " ");
    StringCat (title_txt, sdp->data.ptrvalue);
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = title_txt;
    id_txt = MemFree (id_txt);
  }
  
  bsp->id = SeqIdFree (bsp->id);  
}



static SeqEntryPtr 
ImportOneNucBioseq 
(FILE   *fp, 
 BoolPtr this_chars_stripped, 
 Boolean parse_id, 
 CharPtr supplied_id_txt, 
 CharPtr lastchar, 
 BoolPtr read_from_delta,
 BoolPtr isAsn,
 ValNodePtr PNTR err_msg_list)
{
  SeqEntryPtr sep = NULL;
  BioseqPtr   bsp = NULL;
  ErrSev      oldsev;

  if (feof (fp)) {
    lastchar = 0;
    return NULL;
  }
  oldsev = ErrSetMessageLevel (SEV_MAX);
  /* first, try to read as delta sequence */
  bsp = ReadDeltaFastaEx (fp, NULL, this_chars_stripped);
  ErrSetMessageLevel (oldsev);
  if (bsp == NULL) 
  {
    *read_from_delta = FALSE;
    /* if failed, try regular FASTA */
    /* note - we pass in FALSE for parse_id in SequinFastaToSeqEntryEx
     * because we do not want to use Sequin's auto-generated sequence IDs.
     * We then parse the sequence ID from the title ourselves using
     * ReplaceFakeIDWithIDFromTitle if parse_id is TRUE, or leave the ID
     * as blank to force the user to select a real ID later.
     */

    bsp = ReadFastaOnly (fp,
                         TRUE, FALSE,
                         this_chars_stripped,
                         lastchar);
    if (bsp == NULL)
    {
/*
      if (*lastchar == '[')
      {
        sep = ReadOneSegSet (fp, parse_id, err_msg_list, this_chars_stripped);
      } 
      else  */
      if (lastchar == 0) 
      {
        *lastchar = GetNextCharacterFromFile(fp, isAsn);
      }
    }
    else 
    {
      if (parse_id)
      {
        ReplaceFakeIDWithIDFromTitle (bsp);
      }
      else
      {
        bsp->id = SeqIdFree (bsp->id);
        if (!StringHasNoText (supplied_id_txt))
        {
        	bsp->id = MakeSeqID (supplied_id_txt);
        }
      }
      sep = SeqMgrGetSeqEntryForData (bsp);
      SeqEntryPack (sep);
    }
  } 
  else 
  {
    /* successfully read delta sequence */
    *read_from_delta = TRUE;

    if (!parse_id) 
    {
      PutDeflineIDBackInTitle (bsp);
      if (!StringHasNoText (supplied_id_txt))
      {
        bsp->id = MakeSeqID (supplied_id_txt);
      }
    }
    sep = SeqMgrGetSeqEntryForData (bsp);
  }

  return sep;
}


static Boolean HasGapID (SeqEntryPtr sep)
{
  BioseqPtr bsp;
  Char      id_str [128];
  Int4      j;
  
  if (sep == NULL || ! IS_Bioseq (sep) || (bsp = sep->data.ptrvalue) == NULL)
  {
    return FALSE;
  }
  
  SeqIdWrite (bsp->id, id_str, PRINTID_REPORT, sizeof (id_str));
  
  if (id_str [0] != '?')
  {
    return FALSE;
  }
  if (StringICmp (id_str + 1, "unk100") == 0)
  {
    return TRUE;
  }
  else 
  {
    /* make sure there are only numbers after the question mark */
    j = 1;
    while (isdigit (id_str [j]))
    {
      j++;
    }
    if (id_str [j] == 0)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
}


static Boolean HasNoSeqID (SeqEntryPtr sep)
{
  BioseqPtr bsp;

  if (sep == NULL || ! IS_Bioseq (sep) || (bsp = sep->data.ptrvalue) == NULL)
  {
    return FALSE;
  }
  if (bsp->id == NULL)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static Int4 FindLineForStartOfBadRead (FILE *fp, Int4 pos)
{
  FileCache    fc;
  Int4         line_num = 0;
  Char         line [4096];
  CharPtr      str;
  
  if (fp == NULL || pos == 0) {
    return 0;
  }

  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, 0);
  fseek (fp, 0, SEEK_SET);
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL && FileCacheTell (&fc) < pos) {
      line_num++;
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }
  return line_num;
}


static Int4 FindLineForBadReadChar (FILE *fp, Char badchar)
{
  FileCache    fc;
  Int4         line_num = 0;
  Char         line [4096];
  CharPtr      str;
  
  if (fp == NULL) {
    return 0;
  }

  FileCacheSetup (&fc, fp);
  FileCacheSeek (&fc, 0);
  fseek (fp, 0, SEEK_SET);
  str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  while (str != NULL && StringChr (str, badchar) == NULL) {
      line_num++;
      str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
  }
  return line_num;
}


NLM_EXTERN SeqEntryPtr 
ImportNucleotideFASTASequencesFromFile
(FILE           *fp, 
 Boolean         parse_id,
 CharPtr         supplied_id_txt,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped,
 Boolean         allow_char_stripping)
{
  Int4          count;
  SeqEntryPtr   last;
  Char          lastchar;
  SeqEntryPtr   nextsep;
  BioseqPtr     bsp = NULL;
  SeqEntryPtr   new_sep_list = NULL;
  Boolean       read_from_delta;
  Int4          pos, last_no_id_start = -1;
  Boolean       this_chars_stripped = FALSE;
  Boolean       isASN = FALSE, isOnlyFASTA = FALSE;
  Int4          bad_start = 0, bad_line = 0;
  ValNodePtr    special_list = NULL;
  
  if (chars_stripped != NULL)
  {
    *chars_stripped = FALSE;
  }
    
  count = 0;
  
  new_sep_list = NULL;
  last = NULL;
  
  pos = ftell (fp);

  while ((nextsep = ImportOneNucBioseq (fp, &this_chars_stripped, parse_id, supplied_id_txt, &lastchar, &read_from_delta, &isASN, err_msg_list)) != NULL)
  {
    if (new_sep_list == NULL) {
      new_sep_list = nextsep;
    }
    if (chars_stripped != NULL)
    {
      *chars_stripped |= this_chars_stripped;
    }
    if (last_no_id_start > -1)
    {
      if (HasGapID (nextsep))
      {
        nextsep = SeqEntryFree (nextsep);
        bsp = last->data.ptrvalue;
        SeqMgrDeleteFromBioseqIndex (bsp);
        bsp = BioseqFree (bsp);
        fseek (fp, last_no_id_start, SEEK_SET);
        bsp = ReadDeltaFastaWithEmptyDefline (fp, NULL, chars_stripped);
        last->data.ptrvalue = bsp;
        bsp->id = SeqIdFree (bsp->id);
        last_no_id_start = -1;
      }
      else if (HasNoSeqID (nextsep))
      {
        last_no_id_start = pos;
      }
      else
      {
        last_no_id_start = -1;
      }
    }
    else if (HasNoSeqID (nextsep))
    {
      last_no_id_start = pos;
    }

    /* add to list */
    if (nextsep != NULL)
    {
      if (last == NULL) 
      {
        new_sep_list = nextsep;
        last = nextsep;
      }
      else 
      {
        last->next = nextsep;
        last = nextsep;
      }
    } 
    pos = ftell (fp);
    if (!allow_char_stripping && this_chars_stripped)
    {
      break;
    }
  }
  
  if ((!allow_char_stripping && this_chars_stripped) || (lastchar != (Char) EOF && lastchar != NULLB && lastchar != (Char) 255))
  {
    if (!this_chars_stripped && !isASN) {
      bad_start = FindLineForStartOfBadRead (fp, pos);
      bad_line = FindLineForBadReadChar (fp, lastchar);
      Message (MSG_ERROR, "Unable to read file, starting at line %d (found bad character '%c' at line %d)!", bad_start, lastchar, bad_line);
    }
    new_sep_list = SeqEntryFree (new_sep_list);
  }
  return new_sep_list;
}


static void StripStopCodons (SeqEntryPtr sep_list)
{
  BioseqPtr pbsp;
  Char      prot_str[2];

  prot_str[1] = 0;
  while (sep_list != NULL) {
    if (IS_Bioseq (sep_list) && (pbsp = (BioseqPtr) sep_list->data.ptrvalue) != NULL
        && ISA_aa(pbsp->mol) && pbsp->repr == Seq_repr_raw) {     
      SeqPortStreamInt (pbsp, pbsp->length - 1, pbsp->length - 1, Seq_strand_plus, EXPAND_GAPS_TO_DASHES, (Pointer) (prot_str), NULL);
      if (prot_str[0] == '*') {
        /* trim trailing stop codon */
        BSSeek ((ByteStorePtr) pbsp->seq_data, -1, SEEK_END);
        BSDelete ((ByteStorePtr) pbsp->seq_data, 1);
        pbsp->length -= 1;
      }
    }
    sep_list = sep_list->next;
  }
}


NLM_EXTERN SeqEntryPtr ImportProteinFASTASequences 
(FILE            *fp,
 Boolean         parse_id,
 CharPtr         supplied_id_txt,
 ValNodePtr PNTR err_msg_list,
 BoolPtr         chars_stripped)
{
  Pointer       dataptr;
  Uint2         datatype;
  SeqEntryPtr   new_list = NULL, nextsep, lastsep = NULL;
  Boolean       error_reading = FALSE;
  Int4          pos, bad_start;
  BioseqPtr     bsp;

  pos = ftell (fp);  
  dataptr = ReadAsnFastaOrFlatFileEx (fp, &datatype, NULL, FALSE, TRUE,
                                    FALSE, FALSE, chars_stripped);
  while (dataptr != NULL) 
  {
    bsp = NULL;
    if (datatype == OBJ_SEQENTRY)
    {
      nextsep = (SeqEntryPtr) dataptr;
      if (IS_Bioseq(nextsep)) {
        bsp = nextsep->data.ptrvalue;
      } else {
        nextsep = NULL;
      }      
    }
    else if (datatype == OBJ_BIOSEQ)
    {
      bsp = dataptr;
      nextsep = SeqMgrGetSeqEntryForData(bsp);
    }
    if (bsp != NULL)
    {
      if (parse_id) 
      {
        ReplaceFakeIDWithIDFromTitle ((BioseqPtr) nextsep->data.ptrvalue);
      } 
      else 
      {
        if (!StringHasNoText (supplied_id_txt))
        {
          bsp->id = MakeSeqID (supplied_id_txt);
        }
      } 
      SeqEntryPack (nextsep);
      if (lastsep == NULL) {
        new_list = nextsep;
      } else {
        lastsep->next = nextsep;
      }
      lastsep = nextsep;
      pos = ftell (fp);  
      dataptr = ReadAsnFastaOrFlatFileEx (fp, &datatype, NULL, FALSE, TRUE,
                                          FALSE, FALSE, chars_stripped);
    }
    else
    {
      if (dataptr != NULL) {
        error_reading = TRUE;
      }
      dataptr = NULL;
    }
  }
  if (*chars_stripped || error_reading)
  {
    bad_start = FindLineForStartOfBadRead (fp, pos);
    Message (MSG_ERROR, "Unable to read file, starting at line %d!", bad_start);
    new_list = SeqEntryFree (new_list);
  }
  
  StripStopCodons (new_list);

  return new_list;
  
}


NLM_EXTERN void AddUniqueUpdateSequenceIDs (SeqEntryPtr sep)
{
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL)
  {
    return;
  }
  else if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && bsp->id == NULL)
    {
      bsp->id = MakeUniqueSeqID ("UpdateSequence");
    }
  }
  else if (IS_Bioseq_set (sep))
  {
    /* we could add IDs to segmented sets, but maybe we should just remove them? */
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL)
    {
      AddUniqueUpdateSequenceIDs (bssp->seq_set);
    }
  }
  
  AddUniqueUpdateSequenceIDs (sep->next);
}


NLM_EXTERN void 
ListBioseqsInSeqEntry 
(SeqEntryPtr     sep, 
 Boolean         is_na,
 Int4Ptr         seq_num, 
 ValNodePtr PNTR bioseq_list)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  
  if (sep == NULL || bioseq_list == NULL || seq_num == NULL)
  {
    return;
  }
  if (IS_Bioseq (sep) && sep->data.ptrvalue != NULL)
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_na (bsp->mol))
    {
      if (is_na)
      {
        ValNodeAddPointer (bioseq_list, *seq_num, bsp);
        (*seq_num)++;
      }
    }
    else if (!is_na)
    {
      ValNodeAddPointer (bioseq_list, *seq_num, bsp);
      (*seq_num)++;
    } 
  }
  else if (IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL)
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    ListBioseqsInSeqEntry (bssp->seq_set, is_na, seq_num, bioseq_list);
  }
  
  ListBioseqsInSeqEntry (sep->next, is_na, seq_num, bioseq_list);
}


static Boolean SeqIdListsOverlap (SeqIdPtr sip1, SeqIdPtr sip2)
{ 
  SeqIdPtr sip_next;
  Char     tmp_id_str [255];
  Boolean  rval = FALSE;
  StringConstraint scd;

  MemSet (&scd, 0, sizeof (StringConstraint));

  while (sip1 != NULL && !rval)
  {
    if (SeqIdIn (sip1, sip2))
    {
      rval = TRUE;
    }
    else if (sip1->choice == SEQID_LOCAL)
    {
      sip_next = sip1->next;
      sip1->next = NULL;
      SeqIdWrite (sip1, tmp_id_str, PRINTID_REPORT, sizeof (tmp_id_str) - 1);     
      sip1->next = sip_next;
      scd.match_text = tmp_id_str;
      scd.match_location = String_location_equals;
      rval = DoesSeqIDListMeetStringConstraint (sip2, &scd);
    }
    sip1 = sip1->next;
  }
  return rval;
}


NLM_EXTERN ValNodePtr ShuffleUpdateBioseqList (ValNodePtr PNTR update_bioseq_list, ValNodePtr orig_bioseq_list)
{
  ValNodePtr unmatched_list = NULL;
  ValNodePtr orig_vnp, update_vnp;
  ValNodePtr new_update_list = NULL;
  Int4       bsp_pos = 0, update_pos, pos;
  BioseqPtr  orig_bsp, update_bsp;
  
  if (update_bioseq_list == NULL || *update_bioseq_list == NULL)
  {
    return NULL;
  }
  else if (orig_bioseq_list == NULL)
  {
    unmatched_list = *update_bioseq_list;
    *update_bioseq_list = NULL;
  }

  for (orig_vnp = orig_bioseq_list; orig_vnp != NULL; orig_vnp = orig_vnp->next)
  {
    if (orig_vnp->data.ptrvalue == NULL)
    {
      ValNodeAddPointer (&new_update_list, bsp_pos, NULL);
      bsp_pos ++;
      continue;
    }
    orig_bsp = (BioseqPtr) orig_vnp->data.ptrvalue;
    update_pos = -1;
    for (update_vnp = *update_bioseq_list, pos = 0;
         update_vnp != NULL && update_pos < 0; 
         update_vnp = update_vnp->next, pos++)
    {
      if (update_vnp->data.ptrvalue != NULL)
      {
        update_bsp = (BioseqPtr) update_vnp->data.ptrvalue;
        if (SeqIdListsOverlap (update_bsp->id, orig_bsp->id))
        {
          update_pos = pos;
        }
      }
    }
    if (update_pos >= 0)
    {
      update_vnp = ExtractNthValNode (update_bioseq_list, update_pos);
      update_vnp->choice = bsp_pos;
      ValNodeLink (&new_update_list, update_vnp);
    }
    else
    {
      ValNodeAddPointer (&new_update_list, bsp_pos, NULL);
    }
    bsp_pos++;
  }
  
  unmatched_list = *update_bioseq_list;
  *update_bioseq_list = new_update_list;
  
  /* renumber unmatched_list */
  for (update_vnp = unmatched_list, update_pos = 0;
       update_vnp != NULL;
       update_vnp = update_vnp->next, update_pos++)
  {
    update_vnp->choice = update_pos;
  }
  
  return unmatched_list;
}


/* This function compares the text from a local ID against the 
 * report string from non-local IDs in sip_list, useful when
 * comparing values from a file in which the user did not specify
 * the version or the gb| in the sequence ID.
 */
NLM_EXTERN Boolean RelaxedSeqIdIn (SeqIdPtr sip, SeqIdPtr sip_list)
{
  SeqIdPtr sip_next;
  Char     id_txt1 [128], id_txt2 [128];
  CharPtr  ptr;
  Int4     len;
  DbtagPtr dp;
  Int4     bankit_num;
  
  if (sip == NULL || sip_list == NULL || sip->choice != SEQID_LOCAL)
  {
    return FALSE;
  }
  
  SeqIdWrite (sip, id_txt1, PRINTID_REPORT, sizeof (id_txt1) - 1);
  
  while (sip_list != NULL)
  {
    if (sip_list->choice != SEQID_LOCAL)
    {
      sip_next = sip_list->next;
      sip_list->next = NULL;
      SeqIdWrite (sip_list, id_txt2, PRINTID_REPORT, sizeof (id_txt1) - 1);
      sip_list->next = sip_next;
      if (StringCmp (id_txt1, id_txt2) == 0)
      {
        return TRUE;
      }
      ptr = StringChr (id_txt2, '.');
      if (ptr != NULL)  /* ID in list has version */
      {
        len = StringLen (id_txt1);
        if (len == ptr - id_txt2 && StringNCmp (id_txt1, id_txt2, len) == 0)
        {
          return TRUE;
        }
      }
      if (StringNICmp (id_txt1, "bankit", 6) == 0
          && sip_list->choice == SEQID_GENERAL
          && sip_list->data.ptrvalue != NULL) {
        bankit_num = atoi (id_txt1 + 6);
        dp = (DbtagPtr) sip_list->data.ptrvalue;
        if(StringICmp(dp->db, "BankIt") == 0 && dp->tag != NULL && dp->tag->id == bankit_num) {
          return TRUE;
        }
      }
    }
    sip_list = sip_list->next;
  }
  return FALSE;
}


NLM_EXTERN BioseqPtr FindBioseqInList (ValNodePtr bioseq_list, SeqIdPtr sip, Int4Ptr position)
{
  ValNodePtr vnp;
  BioseqPtr  bsp = NULL;
  Int4       vnp_pos;
  
  if (position != NULL)
  {
    *position = -1;
  }
  if (bioseq_list == NULL)
  {
    return NULL;
  }
  
  for (vnp = bioseq_list, vnp_pos = 0;
       vnp != NULL && bsp == NULL;
       vnp = vnp->next, vnp_pos++)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (SeqIdIn (sip, bsp->id) || RelaxedSeqIdIn (sip, bsp->id))
    {
      if (position != NULL)
      {
        *position = vnp_pos;
      }
    }
    else
    {
      bsp = NULL;
    }
  }
  return bsp;
}


/* This function should find all update Bioseqs that have colliding sequence IDs and
 * replace the colliding IDs with new sequence IDs.
 */
NLM_EXTERN void ReplaceCollidingUpdateIDs (ValNodePtr update_bioseq_list, ValNodePtr orig_bioseq_list)
{
  ValNodePtr vnp, orig_vnp;
  SeqIdPtr   replace_sip, sip;
  BioseqPtr  bsp;
  Char       id_txt [128];
  Int4       orig_pos;
  
  if (update_bioseq_list == NULL || orig_bioseq_list == NULL)
  {
    return;
  }
  
  for (vnp = update_bioseq_list; vnp != NULL; vnp = vnp->next)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      if (FindBioseqInList (orig_bioseq_list, bsp->id, &orig_pos))
      {
        orig_vnp = GetNthValNode (orig_bioseq_list, orig_pos);
        if (orig_vnp != NULL && orig_vnp->data.ptrvalue != NULL)
        {
          replace_sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
          SeqIdWrite (replace_sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
          StringCat (id_txt, "_update");
          sip = MakeUniqueSeqID (id_txt);
          SeqMgrDeleteFromBioseqIndex (bsp);
          BioseqReplaceID (bsp, sip);
          sip = SeqIdFree (sip);   
        
          SeqMgrReplaceInBioseqIndex (orig_vnp->data.ptrvalue);
          SeqMgrReplaceInBioseqIndex (bsp);
        }
      }
    }
  }
  
}


NLM_EXTERN void RemoveSequencesWithoutUpdates (ValNodePtr PNTR orig_bioseq_list, ValNodePtr PNTR update_bioseq_list)
{
  ValNodePtr orig_prev = NULL, update_prev = NULL;
  ValNodePtr orig_next = NULL, update_next = NULL;
  ValNodePtr orig_vnp, update_vnp;
  Int4       seq_num;
  
  if (orig_bioseq_list == NULL || update_bioseq_list == NULL
      || *orig_bioseq_list == NULL || *update_bioseq_list == NULL)
  {
    return;
  }
  
  orig_vnp = *orig_bioseq_list;
  update_vnp = *update_bioseq_list;
  
  while (orig_vnp != NULL && update_vnp != NULL)
  {
    orig_next = orig_vnp->next;
    update_next = update_vnp->next;
    if (orig_vnp->data.ptrvalue == NULL || update_vnp->data.ptrvalue == NULL)
    {
      if (orig_prev == NULL || update_prev == NULL)
      {
        *orig_bioseq_list = orig_vnp->next;
        *update_bioseq_list = update_vnp->next;
      }
      else
      {
        orig_prev->next = orig_vnp->next;
        update_prev->next = update_vnp->next;
      }
      orig_vnp->next = NULL;
      update_vnp->next = NULL;
      ValNodeFree (orig_vnp);
      ValNodeFree (update_vnp);
    }
    else
    {
      orig_prev = orig_vnp;
      update_prev = update_vnp;
    }
    
    orig_vnp = orig_next;
    update_vnp = update_next;
  }
  
  for (orig_vnp = *orig_bioseq_list, update_vnp = *update_bioseq_list, seq_num = 0;
       orig_vnp != NULL && update_vnp != NULL;
       orig_vnp = orig_vnp->next, update_vnp = update_vnp->next, seq_num++)
  {
    orig_vnp->choice = seq_num;
    update_vnp->choice = seq_num;
  }
}


NLM_EXTERN SeqAlignPtr AlignForSequenceUpdate (BioseqPtr bsp1, BioseqPtr bsp2, BoolPtr revcomp, GlobalAlignFunc align_func)
{
  SeqIdPtr    old_id;
  SeqAlignPtr salp = NULL;
  Boolean     put_id_back = FALSE;

  if (bsp1 == NULL || bsp2 == NULL || align_func == NULL) {
    return NULL;
  }

  if (SeqIdIn (bsp2->id, bsp1->id)) {
    old_id = bsp2->id;
    bsp2->id = MakeUniqueSeqID ("UpdateSequence");
    SeqMgrReplaceInBioseqIndex (bsp2);
    put_id_back = TRUE;
  }
  salp = align_func (bsp1, bsp2, revcomp);
  if (put_id_back) {
    bsp2->id = SeqIdFree (bsp2->id);
    bsp2->id = old_id;
    SeqMgrReplaceInBioseqIndex (bsp2);
  }
  return salp;
}




