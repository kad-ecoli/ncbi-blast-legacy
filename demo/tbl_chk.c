/*   tbl_chk.c
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
* File Name:  tbl_chk.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   8/05/10
*
* $Revision: 1.3 $
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
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <explore.h>
#include <pmfapi.h>
#define NLM_GENERATED_CODE_PROTO
#include <asnmacro.h>
#include <objmacro.h>
#include <macroapi.h>
#ifdef INTERNAL_NCBI_TBL_CHK
#include <accpubseq.h>
#endif

#define TBL_CHK_APP_VER "1.0"

CharPtr TBL_CHK_APPLICATION = TBL_CHK_APP_VER;

static Boolean debug_mode = FALSE;

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


static SeqIdPtr SmartGuessMakeId (CharPtr str)
{
  CharPtr id_txt;
  SeqIdPtr sip = NULL;

  if (StringHasNoText (str)) {
    return NULL;
  } else if (StringChr (str, '|') != NULL) {
    sip = MakeSeqID (str);
  } else if (IsAllDigits (str)) {
    id_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 4));
    sprintf (id_txt, "gi|%s", str);
    sip = MakeSeqID (id_txt);
    id_txt = MemFree (id_txt);
  } else if (StringChr (str, '_') != NULL) {
    id_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 5));
    sprintf (id_txt, "oth|%s", str);
    sip = MakeSeqID (id_txt);
    id_txt = MemFree (id_txt);
  } else {
    id_txt = (CharPtr) MemNew (sizeof (Char) * (StringLen (str) + 4));
    sprintf (id_txt, "gb|%s", str);
    sip = MakeSeqID (id_txt);
    id_txt = MemFree (id_txt);
  }
  return sip;
}


typedef struct fetchitem {
  SeqIdPtr sip;
  CharPtr  id_txt;
  Int4     index_pos;
  ValNodePtr field_values;
} FetchItemData, PNTR FetchItemPtr;


static FetchItemPtr FetchItemNew (CharPtr accn) 
{
  FetchItemPtr item;
  
  item = (FetchItemPtr) MemNew (sizeof (FetchItemData));
  item->id_txt = StringSave (accn);
  item->sip = SmartGuessMakeId(accn);
  item->field_values = NULL;
  item->index_pos = -1;
  return item;
}


static ValNodePtr FetchItemFieldValuesFree (ValNodePtr field_values)
{
  ValNodePtr vnp;

  for (vnp = field_values; vnp != NULL; vnp = vnp->next) {
    vnp->data.ptrvalue = ValNodeFreeData (vnp->data.ptrvalue);
  }
  field_values = ValNodeFree (field_values);
  return field_values;
}


static FetchItemPtr FetchItemFree (FetchItemPtr item)
{
  if (item != NULL) {
    item->id_txt = MemFree (item->id_txt);
    item->sip = SeqIdFree (item->sip);
    item->field_values = FetchItemFieldValuesFree (item->field_values);
    item = MemFree (item);
  }
  return item;
}


static ValNodePtr FetchItemListFree (ValNodePtr list)
{
  ValNodePtr list_next;

  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = FetchItemFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}


static int FetchItemCompare (FetchItemPtr f1, FetchItemPtr f2)
{
  if (f1 == NULL || f2 == NULL) {
    return 0;
  } else {
    return StringCmp (f1->id_txt, f2->id_txt);
  }
}


static ValNodePtr FetchItemListFromTable (ValNodePtr table)
{
  ValNodePtr       list = NULL, prev = NULL, this_item, row, col;

  if (table == NULL || table->next == NULL) {
    return NULL;
  }
  for (row = table->next; row != NULL; row = row->next) {
    col = (ValNodePtr) row->data.ptrvalue;
    if (col != NULL) {
      this_item = ValNodeNew (NULL);
      this_item->data.ptrvalue = FetchItemNew (col->data.ptrvalue);
      if (prev == NULL) {
        list = this_item;
      } else {
        prev->next = this_item;
      }
      prev = this_item;
    }
  }
  return list;
}


static void ResetFetchItemListIndex (ValNodePtr list)
{
  FetchItemPtr fetch_item;

  while (list != NULL) {
    fetch_item = (FetchItemPtr) list->data.ptrvalue;
    if (fetch_item != NULL) {
      fetch_item->index_pos = -1;
    } 
    list = list->next;
  }
}


static void FetchItemFieldValuesPrint (FILE *fp, ValNodePtr field_values)
{
  ValNodePtr vnp, vnp_val;

  for (vnp = field_values; vnp != NULL; vnp = vnp->next) {
    for (vnp_val = vnp->data.ptrvalue; vnp_val != NULL; vnp_val = vnp_val->next) {
      fprintf (fp, "%s%s", vnp_val->data.ptrvalue == NULL ? "" : (char *) vnp_val->data.ptrvalue, vnp_val->next == NULL ? "\n" : "\t");
    }
  }
}


static int LIBCALLBACK SortVnpByFetchItem (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return FetchItemCompare (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static FetchItemPtr PNTR sFetchIndex = NULL;
static Int4 sFetchIndexSize = 0;

static void MakeFetchItemIndex(ValNodePtr fetch_list)
{
  Int4 num;
  ValNodePtr tmp_list = NULL, prev = NULL, vnp, vnp_new;

  if (sFetchIndex != NULL) {
    sFetchIndex = MemFree (sFetchIndex);
    sFetchIndexSize = 0;
  }

  if (fetch_list == NULL) {
    return;
  }

  tmp_list = ValNodeNew (NULL);
  tmp_list->data.ptrvalue = fetch_list->data.ptrvalue;
  prev = tmp_list;
  for (vnp = fetch_list->next; vnp != NULL; vnp = vnp->next) {
    vnp_new = ValNodeNew (prev);
    vnp_new->data.ptrvalue = vnp->data.ptrvalue;
    prev = vnp_new;
  }

  tmp_list = ValNodeSort (tmp_list, SortVnpByFetchItem);
  sFetchIndexSize = ValNodeLen (tmp_list);
  sFetchIndex = (FetchItemPtr PNTR) MemNew (sizeof (FetchItemPtr) * sFetchIndexSize);
  for (vnp = tmp_list, num = 0; vnp != NULL; vnp = vnp->next, num++) {
    sFetchIndex[num] = vnp->data.ptrvalue;
  }
  tmp_list = ValNodeFree (tmp_list);
}


static FetchItemPtr FindInFetchIndex (CharPtr str)
{
  Int4 imin = 0;
  Int4 imax = sFetchIndexSize - 1;
  Int4 num = -1, i, j;
  CharPtr tmp;

  while (imax >= imin)
  {
    i = (imax + imin)/2;
    tmp = sFetchIndex[i]->id_txt;
    if ((j = StringCmp(tmp, str)) > 0)
      imax = i - 1;
    else if (j < 0)
      imin = i + 1;
    else
    {
      num = i;
      break;
    }
  }
  if (num == -1) {
    return NULL;
  } else {
    return sFetchIndex[num];
  }
}


#ifdef INTERNAL_NCBI_TBL_CHK

static CharPtr smartfetchproc = "SmartBioseqFetch";

static CharPtr srcchkfetchcmd = NULL;
static CharPtr smartfetchcmd = NULL;


static Int4 NumInSet (BioseqSetPtr bssp)
{
  SeqEntryPtr  sep;
  Int4         num = 0;

  if (bssp != NULL) {
    sep = bssp->seq_set;
    while (sep != NULL) {
      num++;
      sep = sep->next;
    }
  }
  return num;
}
   

static BioseqPtr FetchBioseqFromSmartNotId (CharPtr accn, Uint2Ptr pEntityID)
{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  Char              err_path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;
  Int4              gi = 0;
  ValNodePtr        vnp;
  time_t            t1, t2;

  if (srcchkfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TBL_CHK", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	srcchkfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (srcchkfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	srcchkfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (srcchkfetchcmd == NULL) return NULL;

  TmpNam (path);

  t1 = time(NULL);
#ifdef OS_UNIX
  sprintf (err_path, "%s.err", path);
  sprintf (cmmd, "csh %s %s > %s 2>%s", srcchkfetchcmd, accn, path, err_path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", srcchkfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
#ifdef OS_UNIX
    FileRemove (err_path);
#endif
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
#ifdef OS_UNIX
  FileRemove (err_path);
#endif

  if (dataptr == NULL) return NULL;

  sep = GetTopSeqEntryForEntityID (entityID);

  if (sep == NULL) return NULL;
  sip = SmartGuessMakeId (accn);
  bsp = BioseqFindInSeqEntry (sip, sep);
  sip = SeqIdFree (sip);
  if (debug_mode) {
    t2 = time(NULL);
    if (t2 - t1 > 1) {
      printf("Time to download %s from SMART:%d\n", accn, t2 - t1);
    }
  }
  if (pEntityID != NULL) {
    *pEntityID = entityID;
  }
  return bsp;
}


static Int2 LIBCALLBACK SmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  Char              err_path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;
  Int4              gi = 0;
  ValNodePtr        vnp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (err_path, "%s.err", path);
  sprintf (cmmd, "csh %s %s > %s 2>%s", smartfetchcmd, tsip->accession, path, err_path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
#ifdef OS_UNIX
    FileRemove (err_path);
#endif
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
#ifdef OS_UNIX
  FileRemove (err_path);
#endif

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);

  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean SmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, smartfetchproc, smartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  SmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr tpasmartfetchproc = "TPASmartBioseqFetch";

static CharPtr tpasmartfetchcmd = NULL;

extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];
  Char     err_path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (err_path, "%s.err", path);
  sprintf (cmmd, "csh %s %s > %s 2>%s", tpasmartfetchcmd, accn, path, err_path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
#ifdef OS_UNIX
  FileRemove (err_path);
#endif
  return dataptr;
}


static Int2 LIBCALLBACK TPASmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  Char              err_path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_TPG) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (err_path, "%s.err", path);
  sprintf (cmmd, "csh %s %s > %s 2>%s", tpasmartfetchcmd, tsip->accession, path, err_path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
#ifdef OS_UNIX
    FileRemove (err_path);
#endif
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
#ifdef OS_UNIX
  FileRemove (err_path);
#endif

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean TPASmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, tpasmartfetchproc, tpasmartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  TPASmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}
#endif


static ValNodePtr CollectFieldList(BioseqPtr bsp)
{
  BioSourcePtr biop;
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  ValNodePtr list = NULL, vnp;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    vnp = GetSourceQualFieldListFromBioSource (biop);
    ValNodeLink (&list, vnp);
  }
  return list;
}


static void PrintHeader (FILE *fp, ValNodePtr field_list)
{
  CharPtr txt;

  if (fp == NULL || field_list == NULL) {
    return;
  }
  /* first field accession, second field GI, third field tax ID */
  fprintf (fp, "\t\tTaxID");
  while (field_list != NULL) {
    txt = SummarizeFieldType (field_list);
    fprintf (fp, "\t%s", txt);
    txt = MemFree (txt);
    field_list = field_list->next;
  }
  fprintf (fp, "\n");
}


static Int4 GetTaxIdFromOrgRef (OrgRefPtr orp)
{
  Int4       tax_id = -1;
  ValNodePtr vnp;
  DbtagPtr   d;

  if (orp != NULL)
  {
    for (vnp = orp->db; vnp != NULL; vnp = vnp->next) 
    {
      d = (DbtagPtr) vnp->data.ptrvalue;
      if (StringCmp(d->db, "taxon") == 0) 
      {
        tax_id = d->tag->id;
        break;
      }
    }
  }
  return tax_id;
}


static ValNodePtr CollectBioSourceValues (BioSourcePtr biop, ValNodePtr field_list)
{
  Char       taxid_buf[30];
  ValNodePtr field_values = NULL;
  CharPtr    txt;

  sprintf (taxid_buf, "%d", GetTaxIdFromOrgRef(biop->org));
  ValNodeAddPointer (&field_values, 0, StringSave (taxid_buf));
 
  while (field_list != NULL) {
    txt = GetSourceQualFromBioSource (biop, field_list->data.ptrvalue, NULL);
    ValNodeAddPointer (&field_values, 0, txt);
    field_list = field_list->next;
  }
  return field_values;
}


static ValNodePtr CollectBioseqLineValues (BioseqPtr bsp, ValNodePtr field_list, Boolean want_gi)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Char              id_txt[255], id_txt2[255];
  SeqIdPtr          sip, sip_gi = NULL, sip_gb = NULL;
  ValNodePtr        line_list = NULL, line_values;

  if (bsp == NULL) {
    return NULL;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK
        || (sip->choice == SEQID_EMBL && sip_gb == NULL)
        || (sip->choice == SEQID_SWISSPROT && sip_gb == NULL)
        || (sip->choice == SEQID_DDBJ && sip_gb == NULL)
        || (sip->choice == SEQID_PIR && sip_gb == NULL)) {
      sip_gb = sip;
    } else if (sip->choice == SEQID_GI) {
      sip_gi = sip;
    }
  }

  if (sip_gb == NULL && sip_gi == NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    id_txt2[0] = 0;
  } else {
    if (sip_gb == NULL) {
      id_txt[0] = 0;
    } else {
      SeqIdWrite (sip_gb, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    }
    if (sip_gi == NULL) {
      id_txt2[0] = 0;
    } else {
      SeqIdWrite (sip_gi, id_txt2, PRINTID_REPORT, sizeof (id_txt2) - 1);
    }
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
    line_values = NULL;
    ValNodeAddPointer (&line_values, 0, StringSave (id_txt));
    if (want_gi) {
      ValNodeAddPointer (&line_values, 0, StringSave (id_txt2));
    }
    ValNodeLink (&line_values, CollectBioSourceValues (sdp->data.ptrvalue, field_list));
    ValNodeAddPointer (&line_list, 0, line_values);
  }
  return line_list;
}


#ifdef INTERNAL_NCBI_TBL_CHK

typedef struct populate {
  ValNodePtr field_list;
  Boolean    want_gi;
} PopulateData, PNTR PopulatePtr;

static void PopulateFetchItemCallback (BioseqPtr bsp, Pointer data)
{
  PopulatePtr pp;
  SeqIdPtr    sip;
  TextSeqIdPtr tsip;
  Char         buffer[15];
  FetchItemPtr fetch_item = NULL;

  if (bsp == NULL || ISA_aa(bsp->mol) || (pp = (PopulatePtr)data) == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL && fetch_item == NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI:
        printf (buffer, "%d", sip->data.intvalue);
        fetch_item = FindInFetchIndex(buffer);
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
      case SEQID_OTHER :
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL) {
          fetch_item = FindInFetchIndex(tsip->accession);
        }
        break;
      default :
        break;
    }
  }

  if (fetch_item != NULL && fetch_item->index_pos < 0) {
    /* collect field values */
    fetch_item->field_values = CollectBioseqLineValues (bsp, pp->field_list, pp->want_gi);
    fetch_item->index_pos = 0;
  }
}


static void PopulateFetchItemCache (Uint2 entityID, ValNodePtr fetch_list, ValNodePtr field_list, Boolean want_gi)
{
  SeqEntryPtr sep;
  BioseqSetPtr bssp;
  ValNodePtr   vnp;
  FetchItemPtr fetch_item;
  BioseqPtr    bsp;
  time_t       t1, t2;
  Int4         num_in_list, num_in_set;
  PopulateData pd;

  t1 = time(NULL);
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL && IS_Bioseq_set(sep) && (bssp = (BioseqSetPtr)sep->data.ptrvalue) != NULL
      && bssp->_class != BioseqseqSet_class_nuc_prot) {
    num_in_set = NumInSet (bssp);
    if (num_in_set < sFetchIndexSize) {
      for (vnp = fetch_list; vnp != NULL; vnp = vnp->next) {
        fetch_item = (FetchItemPtr) vnp->data.ptrvalue;
        if (fetch_item->index_pos < 0) {
          bsp = BioseqFindInSeqEntry (fetch_item->sip, sep);
          if (bsp != NULL) {
            /* collect field values */
            fetch_item->field_values = CollectBioseqLineValues (bsp, field_list, want_gi);
            fetch_item->index_pos = 0;
          }
        }
      }
    } else {
      pd.field_list = field_list;
      pd.want_gi = want_gi;
      VisitBioseqsInSep (sep, &pd, PopulateFetchItemCallback);
    }     
  }
  if (debug_mode) {
    t2 = time(NULL);
    if (t2 - t1 > 2) {
      printf ("Time to populate cache: %d\n", t2 - t1);
    }
  }
}

#endif


static void PrintBioSourceLine (FILE *fp, BioSourcePtr biop, ValNodePtr field_list)
{
  CharPtr txt;

  if (fp == NULL || biop == NULL || field_list == NULL) {
    return;
  }

  fprintf (fp, "\t%d", GetTaxIdFromOrgRef(biop->org));
 
  while (field_list != NULL) {
    txt = GetSourceQualFromBioSource (biop, field_list->data.ptrvalue, NULL);
    fprintf (fp, "\t%s", txt == NULL ? "" : txt);
    txt = MemFree (txt);
    field_list = field_list->next;
  }
}


static void PrintBioseqLines (FILE *fp, BioseqPtr bsp, ValNodePtr field_list, Boolean want_gi)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Char              id_txt[255], id_txt2[255];
  SeqIdPtr          sip, sip_gi = NULL, sip_gb = NULL;

  if (fp == NULL || bsp == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK
        || (sip->choice == SEQID_EMBL && sip_gb == NULL)
        || (sip->choice == SEQID_SWISSPROT && sip_gb == NULL)
        || (sip->choice == SEQID_DDBJ && sip_gb == NULL)
        || (sip->choice == SEQID_PIR && sip_gb == NULL)) {
      sip_gb = sip;
    } else if (sip->choice == SEQID_GI) {
      sip_gi = sip;
    }
  }

  if (sip_gb == NULL && sip_gi == NULL) {
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    id_txt2[0] = 0;
  } else {
    if (sip_gb == NULL) {
      id_txt[0] = 0;
    } else {
      SeqIdWrite (sip_gb, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);
    }
    if (sip_gi == NULL) {
      id_txt2[0] = 0;
    } else {
      SeqIdWrite (sip_gi, id_txt2, PRINTID_REPORT, sizeof (id_txt2) - 1);
    }
  }

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
    if (want_gi) {
      fprintf (fp, "%s\t%s", id_txt, id_txt2);
    } else {
      fprintf (fp, "%s", id_txt);
    }
    PrintBioSourceLine (fp, sdp->data.ptrvalue, field_list);
    fprintf (fp, "\n");
  }
}


static void PrintBioseqErrorLine (FILE *fp, SeqIdPtr sip)
{
  Char              id_txt[255];

  if (fp == NULL || sip == NULL) {
    return;
  }

  SeqIdWrite (sip, id_txt, PRINTID_REPORT, sizeof (id_txt) - 1);

  if (sip->choice == SEQID_GI) {
    fprintf (fp, "\t%s\n", id_txt);
  } else {
    fprintf (fp, "%s\t\n", id_txt);
  }
}


/* Args structure contains command-line arguments */

#define i_argInputFile         0
#define o_argOutputFile        1
#define D_argDebugMode         2
#define m_argMetaMode          3

Args myargs [] = {
  {"Input File", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Debug Mode", "F", NULL, NULL,
    TRUE, 'D', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Metadata Mode", "F", NULL, NULL,
    TRUE, 'm', ARG_BOOLEAN, 0.0, 0, NULL}
};


static void SortFieldListForSrcChk (ValNodePtr PNTR field_list)
{
  ValNodePtr vnp, vnp_s, vnp_prev = NULL;

  if (field_list == NULL || *field_list == NULL) return;

  SortUniqueFieldTypeList (field_list);

  /* move taxname to front of list */
  for (vnp = *field_list; vnp != NULL; vnp_prev = vnp, vnp = vnp->next) {
    if (vnp->choice == FieldType_source_qual) {
      vnp_s = vnp->data.ptrvalue;
      if (vnp_s != NULL
          && vnp_s->choice == SourceQualChoice_textqual
          && vnp_s->data.intvalue == Source_qual_taxname) {
        /* only need to move if not already at front of list */
        if (vnp_prev != NULL) {
          vnp_prev->next = vnp->next;
          vnp->next = *field_list;
          *field_list = vnp;
        }
        break;
      }
    }
  }       


}


static ValNodePtr FieldsFromFieldListString (CharPtr str)
{
  CharPtr cpy, val, comma;
  Int4    qual;
  ValNodePtr field_list = NULL, qc;

  if (StringHasNoText (str)) {
    return NULL;
  }
  cpy = StringSave (str);
  val = cpy;
  comma = StringChr(val, ',');
  while (comma != NULL) {
    *comma = 0;
    qual = GetSourceQualTypeByName(val);
    if (qual < 0) {
      Message (MSG_ERROR, "%s is not a recognized source field name", val);
    } else {
      qc = ValNodeNew (NULL);
      qc->choice = SourceQualChoice_textqual;
      qc->data.intvalue = qual;
      ValNodeAddPointer (&field_list, FieldType_source_qual, qc);
    }
    *comma = ',';
    val = comma + 1;
    comma = StringChr (val, ',');
  }

  qual = GetSourceQualTypeByName(val);
  if (qual < 0) {
    Message (MSG_ERROR, "%s is not a recognized source field name", val);
  } else {
    qc = ValNodeNew (NULL);
    qc->choice = SourceQualChoice_textqual;
    qc->data.intvalue = qual;
    ValNodeAddPointer (&field_list, FieldType_source_qual, qc);
  }

  cpy = MemFree (cpy);
  return field_list;
}


static Boolean PvtStrToLong (CharPtr str, Int4Ptr longval)

{
   Char     ch;
   Int2     i;
   Int2     len;
   Char     local [64];
   Boolean  nodigits;
   Boolean  rsult;
   long int     val;

   rsult = FALSE;
   if (longval != NULL) {
     *longval = (Int4) 0;
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
       if (sscanf (local, "%ld", &val) == 1) {
         *longval = val;
       }
     }
   }
   return rsult;
}

static Int4 AccessionToGi (CharPtr string)

{
    Int4      gi;
    SeqIdPtr  sip;

    sip = SeqIdFromAccessionDotVersion (string);
    if (sip == NULL) return 0;
    gi = GetGIForSeqId (sip);
    SeqIdFree (sip);
    return gi;
}


static BioseqPtr FetchOnlyBioseqFromID (CharPtr str)
{
   BioseqPtr    bsp = NULL;
   Uint2        entityID;
   SeqEntryPtr  sep;
   SeqIdPtr     sip = NULL;
   Int4         uid;
   time_t       t1, t2;

   t1 = time (NULL);
   uid = 0;
   TrimSpacesAroundString (str);
   if (IsAllDigits (str)) {
     if (PvtStrToLong (str, &uid)) {
       sip = ValNodeNew(NULL);
       sip->choice = SEQID_GI;
       sip->data.intvalue = uid;
     } else {
      uid = 0;
     }
   } else {
     sip = SeqIdFromAccessionDotVersion (str);
     uid = AccessionToGi (str);
   }
   sep = NULL;
   if (uid > 0) {
     sep = PubSeqSynchronousQuery (uid, 3, 0); /* retcode was 0 */
     if (sep != NULL) {
       bsp = BioseqFindInSeqEntry (sip, sep);
       entityID = ObjMgrGetEntityIDForChoice (sep); 
     }
   }
   sip = SeqIdFree (sip);
   if (debug_mode) {
    t2 = time (NULL);
    if (t2 - t1 > 1) {
      printf ("Time to download %s from ID:%d\n", str, (int) (t2 - t1));
    }
   }
   return bsp;
}


static void CollectFieldsInFetchItemCache (Uint2 entityID, ValNodePtr PNTR field_list, ValNodePtr fetch_list)
{
  SeqEntryPtr sep;
  BioseqSetPtr bssp;
  ValNodePtr   vnp;
  FetchItemPtr fetch_item;
  BioseqPtr    bsp;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL && IS_Bioseq_set(sep) && (bssp = (BioseqSetPtr)sep->data.ptrvalue) != NULL
      && bssp->_class != BioseqseqSet_class_nuc_prot) {
    for (vnp = fetch_list; vnp != NULL; vnp = vnp->next) {
      fetch_item = (FetchItemPtr) vnp->data.ptrvalue;
      if (fetch_item->index_pos < 0) {
        bsp = BioseqFindInSeqEntry (fetch_item->sip, sep);
        if (bsp != NULL) {
          ValNodeLink (field_list, CollectFieldList(bsp));
          fetch_item->index_pos = 0;
        }
      }
    }
  }
}


static ValNodePtr CollectBioseqFieldsWithCaching (ValNodePtr fetch_list) 
{
  ValNodePtr   vnp;
  ValNodePtr   field_list = NULL;
  BioseqPtr    bsp;
  FetchItemPtr fetch_item;
  Uint2        entityID;

  for (vnp = fetch_list; vnp != NULL; vnp = vnp->next) {
    fetch_item = (FetchItemPtr) vnp->data.ptrvalue;
    if (fetch_item->index_pos < 0) {
      entityID = 0;
#ifdef INTERNAL_NCBI_TBL_CHK
      bsp = FetchBioseqFromSmartNotId (fetch_item->id_txt, &entityID);
      if (bsp != NULL) {
        CollectFieldsInFetchItemCache (entityID, &field_list, vnp->next);
      } else {
        bsp = FetchOnlyBioseqFromID (fetch_item->id_txt);
      }
#else
      bsp = FetchOnlyBioseqFromID (fetch_item->id_txt);
#endif
      ValNodeLink (&field_list, CollectFieldList(bsp));
    }
  }
  ResetFetchItemListIndex (fetch_list);
  return field_list;
}


static Boolean HasMismatch
(ValNodePtr values,
 ValNodePtr col)
{
  Boolean rval = FALSE;

  while (!rval && values != NULL && col != NULL) {
    if (StringHasNoText (values->data.ptrvalue)
        && StringHasNoText (col->data.ptrvalue)) {
      /* both empty, ignore */
    } else if (StringCmp (values->data.ptrvalue, col->data.ptrvalue)) {
      rval = TRUE;
    }
    values = values->next;
    col = col->next;
  }
  if (values != NULL || col != NULL) {
    rval = TRUE;
  }
  return rval;
}
 
static void PrintMisMatches(FILE *out, CharPtr id_txt, ValNodePtr val, ValNodePtr col, ValNodePtr field)
{
  CharPtr label = NULL;

  while (field != NULL && (val != NULL || col != NULL)) {
    if (val == NULL) {
      if (!StringHasNoText (col->data.ptrvalue)) {
        label = SummarizeFieldType (field);
        fprintf (out, "%s\t%s\t\t%s\n", id_txt, label, (char *) col->data.ptrvalue);
      }
    } else if (col == NULL) {
      if (!StringHasNoText (val->data.ptrvalue)) {
        label = SummarizeFieldType (field);
        fprintf (out, "%s\t%s\t%s\t\n", id_txt, label, (char *) val->data.ptrvalue);
      }
    } else if (!StringHasNoText (val->data.ptrvalue) && !StringHasNoText (col->data.ptrvalue)
               && StringCmp (val->data.ptrvalue, col->data.ptrvalue) != 0) {
      label = SummarizeFieldType (field);
      fprintf (out, "%s\t%s\t%s\t%s\n", id_txt, label, (char *) val->data.ptrvalue, (char *) col->data.ptrvalue);
    }
    label = MemFree (label);

    if (val != NULL) {
      val = val->next;
    }
    if (col != NULL) {
      col = col->next;
    }
    field = field->next;
  }
}

static Int2 
ProcessBioseqsWithCaching 
(ValNodePtr fetch_list,
 ValNodePtr field_list,
 ValNodePtr table,
 Boolean    meta_mode,
 FILE *out)
{
  ValNodePtr vnp, row, col, line, val;
  BioseqPtr  bsp;
  FetchItemPtr fetch_item;
  Uint2        entityID;
  Int4         num_processed = 0;
  time_t       t, last_t;

  last_t = time(NULL);

  for (vnp = fetch_list, row = table->next; 
       vnp != NULL && row != NULL; 
       vnp = vnp->next, row = row->next) {
    fetch_item = (FetchItemPtr) vnp->data.ptrvalue;
    if (fetch_item->index_pos < 0) {
      entityID = 0;
#ifdef INTERNAL_NCBI_TBL_CHK
      bsp = FetchBioseqFromSmartNotId (fetch_item->id_txt, &entityID);
      if (bsp != NULL) {
        PopulateFetchItemCache (entityID, vnp->next, field_list, TRUE);
      } else {
        bsp = FetchOnlyBioseqFromID (fetch_item->id_txt);
        if (bsp == NULL) {
          printf ("Unable to download Bioseq for %s\n", fetch_item->id_txt);
        }
      }
#else
      bsp = FetchOnlyBioseqFromID (fetch_item->id_txt);
#endif
      fetch_item->field_values = CollectBioseqLineValues (bsp, field_list, FALSE);
    }
    if (fetch_item->field_values == NULL) {
      fprintf (out, "%s: Bioseq not found\n", fetch_item->id_txt);
    } else {
      for (line = fetch_item->field_values; line != NULL; line = line->next) {
        val = line->data.ptrvalue;
        col = row->data.ptrvalue;
        if (meta_mode) {
          PrintMisMatches(out, fetch_item->id_txt, val->next->next, col->next, field_list);
        } else {
          if (HasMismatch(val->next->next, col->next)) {
            fprintf (out, "%s\n", fetch_item->id_txt);
          }
        }
      }
      fetch_item->field_values = FetchItemFieldValuesFree (fetch_item->field_values);
    }
    fflush(out);
    if (debug_mode) {
      num_processed++;
      if (num_processed == 10) {
        t = time (NULL);
        printf ("Time to process %d: %d\n", num_processed, (int) (t - last_t)); 
        num_processed = 0;
        last_t = t;
      }    
    }
  }

  return 0;
}


static ValNodePtr GetFieldListFromHeader(ValNodePtr col)
{
  ValNodePtr val, list = NULL, field, prev = NULL;

  /* note - first item is accession */
  if (col == NULL || col->next == NULL) {
    return NULL;
  }
  for (val = col->next; val != NULL; val = val->next) {
    field = FieldTypeFromString (val->data.ptrvalue);
    if (field == NULL) {
      Message (MSG_ERROR, "Unable to match field name for %s", val->data.ptrvalue);
      list = FieldTypeListFree (list);
      return NULL;
    } else if (prev == NULL) {
      list = field;
    } else {
      prev->next = field;
    }
    prev = field;
  }
  return list;
}


Int2 Main(void)
{
  Char             app [64];
  Int4             rval = 0;
  CharPtr          id_file;
  ValNodePtr       fetch_list = NULL;
  ValNodePtr       field_list = NULL;
  ValNodePtr       table;
  Boolean          meta_mode = FALSE;
  FILE *fp;

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

#ifdef INTERNAL_NCBI_TBL_CHK
    SmartFetchEnable ();
    TPASmartFetchEnable ();

    if (! PUBSEQBioseqFetchEnable ("tbl_chk", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    } 
#else
  PubSeqFetchEnable ();
#endif

  /* process command line arguments */

  sprintf (app, "tbl_chk %s", TBL_CHK_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  id_file = (CharPtr) myargs [i_argInputFile].strvalue;
  debug_mode = (Boolean) myargs [D_argDebugMode].intvalue;
  meta_mode = (Boolean) myargs [m_argMetaMode].intvalue;

  fp = FileOpen (id_file, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", id_file);
    return 1;
  }

  table = ReadTabTableFromFile (fp);
  FileClose (fp);
  if (table == NULL || table->next == NULL) {
    Message (MSG_ERROR, "Table must have at least two rows, one header and one data");
    return 1;
  }

  field_list = GetFieldListFromHeader(table->data.ptrvalue);
  if (field_list == NULL) {
    Message (MSG_ERROR, "Unable to read table header");
    table = FreeTabTable (table);
    return 1;
  }
  fetch_list = FetchItemListFromTable (table);
  MakeFetchItemIndex(fetch_list);

  fp = FileOpen ((CharPtr) myargs [o_argOutputFile].strvalue, "w");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", (CharPtr) myargs [o_argOutputFile].strvalue);
    rval = 1;
  } else {
    if (meta_mode) {
      fprintf (fp, "Accession\tField\tOld Value\tNew Value\n");
    }
    ProcessBioseqsWithCaching (fetch_list, field_list, table, meta_mode, fp);
  }
  FileClose (fp);
  field_list = FieldTypeListFree (field_list);
  fetch_list = FetchItemListFree (fetch_list);
  MakeFetchItemIndex(NULL);

  return rval;
}

