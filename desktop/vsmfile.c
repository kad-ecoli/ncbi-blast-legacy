/*   vsmfile.c
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
* File Name:  vsm.c
*
* Author:  Jim Ostell
*
* Version Creation Date:   11-29-94
*
* $Revision: 6.27 $
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

#include <vsmpriv.h>
#include <tofasta.h>
#include <objmgr.h>
#include <dlogutil.h>
#include <asn2gnbk.h>
#include <explore.h>
#include <asn2gnbp.h>

/*****************************************************************************
*
*   VSMFileInit()
*        Initialize VSM file I/O routines
*
*****************************************************************************/
Boolean LIBCALL VSMFileInit(void)
{
                         /** register default functions */
                         /*** OPEN ***/
    ObjMgrProcLoad(OMPROC_OPEN, "Read FASTA Protein File","FASTA protein file", 0,0,0,0,NULL,
                           VSMFastaProtOpen, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_OPEN, "Read FASTA Nucleotide File","FASTA nucleotide file", 0,0,0,0,NULL,
                           VSMFastaNucOpen, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_OPEN, "Read Binary ASN1 File","ASN.1 binary file", 0,0,0,0,NULL,
                           VSMGenericBinAsnOpen, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_OPEN, "Read Text ASN1 File","ASN.1 text file", 0,0,0,0,NULL,
                           VSMGenericTextAsnOpen, PROC_PRIORITY_DEFAULT);
   
                         /*** SAVE ***/   
  ObjMgrProcLoad(OMPROC_SAVE, "Save Descriptor File", "Descriptor file", 0,0,0,0,NULL,
                           VSMDescriptorAsnSave, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_SAVE, "Export Protein Feature Table","Protein Feature Table", 0,0,0,0,NULL,
                           VSMExportProteinFeatureTable, PROC_PRIORITY_DEFAULT);

  ObjMgrProcLoadEx (OMPROC_SAVE,  
                    "Export Nucleotide Feature Table, Selected Features Only, Suppress Protein IDs",
                    "Selected Features Only, Suppress Protein IDs",
                    0,0,0,0,NULL,
                    VSMExportNucFeatureTableSelectedFeaturesSuppressProteinIDs, 
                    PROC_PRIORITY_DEFAULT, "Nucleotide Feature Table");

  ObjMgrProcLoadEx (OMPROC_SAVE,  
                    "Export Nucleotide Feature Table, Selected Features Only",
                    "Selected Features Only",
                    0,0,0,0,NULL,
                    VSMExportNucFeatureTableSelectedFeatures, 
                    PROC_PRIORITY_DEFAULT, "Nucleotide Feature Table");                   

    ObjMgrProcLoadEx(OMPROC_SAVE, 
                     "Export Nucleotide Feature Table Without Sources, Suppress Protein IDs", 
                     "Without Sources, Suppress Protein IDs",
                     0,0,0,0,NULL,
                   VSMExportNucFeatureTableWithoutSourcesSuppressProteinIDs, 
                   PROC_PRIORITY_DEFAULT, "Nucleotide Feature Table");

    ObjMgrProcLoadEx(OMPROC_SAVE, 
                     "Export Nucleotide Feature Table Without Sources", 
                     "Without Sources",
                     0,0,0,0,NULL,
                   VSMExportNucFeatureTableWithoutSources, 
                   PROC_PRIORITY_DEFAULT, "Nucleotide Feature Table");

  ObjMgrProcLoadEx(OMPROC_SAVE, 
                   "Export Nucleotide Feature Table, Suppress Protein IDs",
                   "Suppress Protein IDs",
                     0,0,0,0,NULL,
                   VSMExportNucFeatureTableSuppressProteinIDs,
                   PROC_PRIORITY_DEFAULT, "Nucleotide Feature Table");

  ObjMgrProcLoadEx(OMPROC_SAVE, 
                   "Export Nucleotide Feature Table",
                   "Do Not Suppress Protein IDs",
                     0,0,0,0,NULL,
                   VSMExportNucFeatureTable,
                   PROC_PRIORITY_DEFAULT, "Nucleotide Feature Table");

    ObjMgrProcLoad(OMPROC_SAVE, "Export Nucleotide Feature Table, Suppress Protein IDs","Nucleotide Feature Table, Suppress Protein IDs", 0,0,0,0,NULL,
                           VSMExportNucFeatureTableSuppressProteinIDs, PROC_PRIORITY_DEFAULT);

    ObjMgrProcLoad(OMPROC_SAVE, "Export Nucleotide Feature Table","Nucleotide Feature Table", 0,0,0,0,NULL,
                           VSMExportNucFeatureTable, PROC_PRIORITY_DEFAULT);

  ObjMgrProcLoad(OMPROC_SAVE, "Save Entries Within Selected Sets", "ASN.1 text files for each entry in selected sets", 0,0,0,0,NULL,
                           VSMSaveSetsAsFiles, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_SAVE, "Save As FASTA Protein File With Products","FASTA protein file with products", 0,0,0,0,NULL,
                           VSMFastaProtSaveWithProduct, PROC_PRIORITY_DEFAULT);

    ObjMgrProcLoad(OMPROC_SAVE, "Save As Sorted FASTA Protein File","Sorted FASTA protein file", 0,0,0,0,NULL,
                           VSMFastaSortedProtSave, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_SAVE, "Save As FASTA Protein File","FASTA protein file", 0,0,0,0,NULL,
                           VSMFastaProtSave, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_SAVE, "Save As FASTA Nucleotide File","FASTA nucleotide file", 0,0,0,0,NULL,
                           VSMFastaNucSave, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_SAVE, "Save As Binary ASN1 File","ASN.1 binary file", 0,0,0,0,NULL,
                           VSMGenericBinAsnSave, PROC_PRIORITY_DEFAULT);
                           
    ObjMgrProcLoad(OMPROC_SAVE, "Save As Text ASN1 File","ASN.1 text file", 0,0,0,0,NULL,
                           VSMGenericTextAsnSave, PROC_PRIORITY_DEFAULT);

    return TRUE;
}

static void PromoteToSeqEntry (Uint2 entityID, Uint2 datatype, Pointer dataptr)

{
    BioseqPtr     bsp;
    BioseqSetPtr  bssp;
    SeqEntryPtr   sep;

    sep = GetTopSeqEntryForEntityID (entityID);
    if (sep != NULL) return;
    sep = SeqEntryNew ();
    if (sep == NULL) return;
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

Int2 LIBCALLBACK VSMGenericTextAsnOpen ( Pointer data )
{
    Char filename[255];
    Pointer ptr = NULL;
    Uint2 entityID, datatype;
    Int2 retval = OM_MSG_RET_ERROR;
    OMProcControlPtr ompcp;

    ompcp = (OMProcControlPtr)data;

    filename[0] = '\0';
    if (GetInputFileName(filename, (size_t)254, NULL, NULL))
    {
        WatchCursor();
        ptr = ObjMgrGenericAsnTextFileRead (filename, &datatype, &entityID);
        ArrowCursor();
        if (ptr == NULL) goto erret;

        ompcp->output_data = ptr;
        ompcp->output_entityID = entityID;
        PromoteToSeqEntry (entityID, datatype, ptr);

        retval = OM_MSG_RET_DONE;
    }
    else
        retval = OM_MSG_RET_OK;

ret:
    return retval;
erret:
    goto ret;
}

typedef struct vsmreadbinstr {
    Boolean do_it,
            window_done;
    Int2 the_type;
    PopuP p;
} VSMReadBinStr, PNTR VSMReadBinStrPtr;

static void AsnBinAcceptProc (ButtoN b)
{
    WindoW w;
    VSMReadBinStrPtr vrp;
    
    w = ParentWindow(b);
    vrp = (VSMReadBinStrPtr) GetWindowExtra(w);
    vrp->do_it = TRUE;
    vrp->the_type = GetValue(vrp->p);
    Remove(w);
    vrp->window_done = TRUE;
    return;    
}

static void AsnBinCancelProc (ButtoN b)
{
    WindoW w;
    VSMReadBinStrPtr vrp;
    
    w = ParentWindow(b);
    vrp = (VSMReadBinStrPtr) GetWindowExtra(w);
    vrp->do_it = FALSE;
    Remove(w);
    vrp->window_done = TRUE;
    return;    
}

Int2 LIBCALLBACK VSMGenericBinAsnOpen ( Pointer data )
{
    Char filename[255];
    AsnIoPtr aip;
    Pointer ptr;
    Uint2 entityID;
    Int2 ct, i, retval = OM_MSG_RET_OK;
    ObjMgrPtr omp;
    ObjMgrTypePtr omtp = NULL;
    OMProcControlPtr ompcp;
    WindoW w;
    GrouP g;
    ButtoN b;
    VSMReadBinStr vrb;
    VSeqMgrPtr vsmp;

    ompcp = (OMProcControlPtr)data;
    vsmp = VSeqMgrGet();
    omp = vsmp->omp;

    filename[0] = '\0';
    if (GetInputFileName(filename, (size_t)254, NULL, NULL))
    {
        vrb.do_it = FALSE;
        vrb.window_done = FALSE;
        vrb.the_type = 0;
        w = ModalWindow (-50, -33, -10, -10, NULL);
        SetWindowExtra(w, &vrb, NULL);
        g = HiddenGroup(w, 0, 2, NULL);
        StaticPrompt(g, "Select ASN.1 type:", 0,0,systemFont,'l');
        vrb.p = PopupList(g, TRUE, NULL);
        i = 0;
        ct = 0;
        omtp = NULL;
        while ((omtp = ObjMgrTypeFindNext(omp, omtp)) != NULL)
        {
            if (omtp->asnname != NULL)
            {
                i++;
                PopupItem(vrb.p, omtp->asnname);
                if (! StringCmp(vsmp->lastASNtype, omtp->asnname))
                    ct = i;
            }
        }
        
        if (! i)
        {
            ErrPostEx(SEV_ERROR,0,0, "No ASN.1 types are registered");
            Remove(w);
            return OM_MSG_RET_ERROR;
        }
        
        if (! ct)
            ct = 1;
        SetValue(vrb.p, ct);
        
        g = HiddenGroup(w, 2, 0, NULL);
        DefaultButton(g, "Accept", AsnBinAcceptProc);
        b = PushButton(g, "Cancel", AsnBinCancelProc);
        
        Show(w);
        Nlm_WaitForCondition (! vrb.window_done);
        ProcessAnEvent();
        
        if (! vrb.do_it)
            return retval;
        
        i = 0;
        omtp = NULL;
        while ((omtp = ObjMgrTypeFindNext(omp, omtp)) != NULL)
        {
            i++;
            if (i == vrb.the_type)
                break;
        }
        
        if (omtp == NULL)
        {
            ErrPostEx(SEV_ERROR,0,0,"Couldn't find vrb.the_type");
            return OM_MSG_RET_ERROR;
        }
        
        StringMove(vsmp->lastASNtype, omtp->asnname);    
        WatchCursor();
        aip = AsnIoOpen(filename, "rb");
        ptr = (*(omtp->asnread))(aip, NULL);
        AsnIoClose(aip);
        if (ptr == NULL)
        {
            ErrPostEx(SEV_ERROR,0,0,"Couldn't read [%s], type [%s]", filename, omtp->asnname);
            retval = OM_MSG_RET_ERROR;
        }
        else
        {
            entityID = ObjMgrRegister(omtp->datatype, ptr);
            ompcp->output_data = ptr;
            ompcp->output_entityID = entityID;
            PromoteToSeqEntry (entityID, omtp->datatype, ptr);

            retval = OM_MSG_RET_DONE;
        }
        ArrowCursor();
    }
    else
        retval = OM_MSG_RET_OK;

    return retval;
}

static Int2 LIBCALLBACK VSMGenericFastaOpen ( Boolean is_na )
{
    Char filename[255];
    FILE * fp;
    SeqEntryPtr sep;
    
    filename[0] = '\0';
    if (GetInputFileName(filename, (size_t)254, NULL, NULL))
    {
        WatchCursor();
        fp = FileOpen(filename, "r");
        while ((sep = FastaToSeqEntry(fp, is_na)) != NULL)
            ObjMgrRegister(OBJ_SEQENTRY, (Pointer)sep);
        FileClose(fp);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;
}

Int2 LIBCALLBACK VSMFastaProtOpen ( Pointer data )
{
    return VSMGenericFastaOpen(FALSE);
}

Int2 LIBCALLBACK VSMFastaNucOpen ( Pointer data )
{
    return VSMGenericFastaOpen(TRUE);
}


static Int2 LIBCALLBACK VSMGenericFastaSave (OMProcControlPtr ompcp, Boolean is_na)
{
    Char filename[255];
    FILE * fp;
    ValNode vn;
    SeqEntryPtr sep = NULL;
    SeqFeatPtr sfp;
    SeqLocPtr slp;
    Uint1 code;
    BioseqPtr bsp;
    
    sfp = NULL;
    switch(ompcp->input_itemtype)
    {
        case OBJ_SEQENTRY:
        case OBJ_BIOSEQ:
        case OBJ_BIOSEQSET:
            break;
        case OBJ_SEQFEAT:
            sfp = (SeqFeatPtr) ompcp->input_data;
            if (sfp == NULL) return OM_MSG_RET_ERROR;
            break;
        default:
            ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can only write Seq-entry, Bioseq, or Bioseq-set");
            return OM_MSG_RET_ERROR;
    }
    if (sfp != NULL) {
    } else if (ompcp->input_choicetype == OBJ_SEQENTRY)
        sep = (SeqEntryPtr)(ompcp->input_choice);
    else
    {
        vn.next = NULL;
        vn.data.ptrvalue = ompcp->input_data;
        if (ompcp->input_itemtype == OBJ_BIOSEQ)
            vn.choice = 1;
        else
            vn.choice = 2;
        sep = &vn;
    }
        
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
#ifdef WIN_MAC
        fp = FileOpen (filename, "r");
        if (fp != NULL) {
            FileClose (fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif
        fp = FileOpen(filename, "w");
        if (sfp != NULL) {
            if (is_na)
                code = Seq_code_iupacna;
            else
                code = Seq_code_ncbieaa;
            slp = sfp->location;
            if (sfp->data.choice == SEQFEAT_CDREGION && (! is_na)) {
              slp = sfp->product;
            }
            /*
            spp = SeqPortNewByLoc (slp, code);
            */
            bsp = GetBioseqGivenSeqLoc (slp, ompcp->input_entityID);
            /*
            if (spp != NULL && bsp != NULL) {
                while (FastaSeqLine(spp, buf, 70, is_na))
                    FastaFileFunc(bsp, FASTA_SEQLINE, buf, sizeof (buf), (Pointer)fp);
                SeqPortFree(spp);
                FastaFileFunc(bsp, FASTA_EOS, buf, sizeof (buf), (Pointer)fp);
            }
            */
            if (slp != NULL && bsp != NULL) {
                SeqLocFastaStream (slp, fp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, 70, 0, 0);
            }
        } else {
            /*
            SeqEntryToFasta(sep, fp, is_na);
            */
            SeqEntryFastaStreamEx (sep, fp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, 70, 0, 0, is_na, !is_na, FALSE, !is_na, FALSE);
        }
        FileClose(fp);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;
}

Int2 LIBCALLBACK VSMFastaProtSave ( Pointer data )
{
    return VSMGenericFastaSave((OMProcControlPtr)data, FALSE);
}

Int2 LIBCALLBACK VSMFastaNucSave ( Pointer data )
{
    return VSMGenericFastaSave((OMProcControlPtr)data, TRUE);
}

typedef struct fastaexportoptions {
  FILE *fp;
  StreamFlgType  flags;
  Int2           linelen;
  Int2           blocklen;
  Int2           grouplen;
} FastaExportOptionsData, PNTR FastaExportOptionsPtr;

static void WriteOneProteinWithProduct (BioseqPtr bsp, Pointer data)
{
    FastaExportOptionsPtr fe;
    SeqFeatPtr            prot;
    SeqMgrFeatContext     fcontext;
    Char                  id [128];

    if (bsp != NULL && ISA_aa (bsp->mol) && (fe = (FastaExportOptionsPtr) data) != NULL) {
        prot = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_PROT, &fcontext);
        if (prot == NULL) {
            BioseqFastaStreamEx (bsp, fe->fp, fe->flags, fe->linelen, fe->blocklen, fe->grouplen,
                                 TRUE, FALSE, FALSE);
        } else {
            
            SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
            fprintf (fe->fp, ">%s [prot=%s]\n", id, fcontext.label);
            BioseqFastaStreamEx (bsp, fe->fp, fe->flags, fe->linelen, fe->blocklen, fe->grouplen,
                                 FALSE, FALSE, FALSE);
        }
    }
}


Int2 LIBCALLBACK VSMFastaProtSaveWithProduct ( Pointer data )
{

    Char filename[255];
    FastaExportOptionsData fe;
    ValNode vn;
    SeqEntryPtr sep = NULL;
    SeqFeatPtr sfp;
    SeqLocPtr slp;
    Uint1 code;
    BioseqPtr bsp;    
    OMProcControlPtr ompcp;

    if ((ompcp = (OMProcControlPtr) data) == NULL) {
      return OM_MSG_RET_ERROR;
    }

    sfp = NULL;
    switch(ompcp->input_itemtype)
    {
        case OBJ_SEQENTRY:
        case OBJ_BIOSEQ:
        case OBJ_BIOSEQSET:
            break;
        case OBJ_SEQFEAT:
            sfp = (SeqFeatPtr) ompcp->input_data;
            if (sfp == NULL || (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_PROT)) return OM_MSG_RET_ERROR;
            break;
        default:
            ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can only write Seq-entry, Bioseq, or Bioseq-set");
            return OM_MSG_RET_ERROR;
    }
    if (sfp != NULL) {
    } else if (ompcp->input_choicetype == OBJ_SEQENTRY)
        sep = (SeqEntryPtr)(ompcp->input_choice);
    else
    {
        vn.next = NULL;
        vn.data.ptrvalue = ompcp->input_data;
        if (ompcp->input_itemtype == OBJ_BIOSEQ)
            vn.choice = 1;
        else
            vn.choice = 2;
        sep = &vn;
    }
        
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
#ifdef WIN_MAC
        fe.fp = FileOpen (filename, "r");
        if (fe.fp != NULL) {
            FileClose (fe.fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif
        fe.fp = FileOpen(filename, "w");
        if (sfp != NULL) {
            code = Seq_code_ncbieaa;
            slp = sfp->location;
            if (sfp->data.choice == SEQFEAT_CDREGION) {
              slp = sfp->product;
            }
            bsp = GetBioseqGivenSeqLoc (slp, ompcp->input_entityID);
            if (slp != NULL && bsp != NULL) {
                SeqLocFastaStream (slp, fe.fp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, 70, 0, 0);
            }
        } else {
            fe.flags = STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL;
            fe.linelen = 70;
            fe.blocklen = 0;
            fe.grouplen = 0;
            VisitBioseqsInSep (sep, &fe, WriteOneProteinWithProduct);
        }
        FileClose(fe.fp);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;
}

typedef struct alphaprot
{
  BioseqPtr bsp;
  CharPtr   prot_name;  
} AlphaProtData, PNTR AlphaProtPtr;

static void GetProtListCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR   pList;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  ProtRefPtr        prp;
  AlphaProtPtr      app;
  
  if (bsp == NULL || userdata == NULL || ! ISA_aa (bsp->mol)) return;
  pList = (ValNodePtr PNTR) userdata;
  app = (AlphaProtPtr) MemNew (sizeof (AlphaProtData));
  if (app == NULL) return;
  app->bsp = bsp;
  app->prot_name = NULL;
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &fcontext);
  if (sfp != NULL && sfp->data.value.ptrvalue != NULL) 
  {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp->name != NULL)
    {
      app->prot_name = StringSave (prp->name->data.ptrvalue);
    }
    else
    {
      app->prot_name = StringSave (fcontext.label);
    }
  }
  ValNodeAddPointer (pList, 0, app);
}

static int LIBCALLBACK SortProtByName (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1, vnp2;
  AlphaProtPtr app1, app2;
  
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      app1 = (AlphaProtPtr) vnp1->data.ptrvalue;
      app2 = (AlphaProtPtr) vnp2->data.ptrvalue;
      if (app1 != NULL && app2 != NULL) {
        return StringCmp (app1->prot_name, app2->prot_name);
      }
    }
  }
  return 0;
}

static void WriteSortedProteinsToFile (FILE *fp, SeqEntryPtr sep)
{
    AlphaProtPtr app;
  ValNodePtr   prot_list = NULL, vnp;
  CharPtr      prev_name = NULL;

  if (fp == NULL || sep == NULL) {
    return;
  }

  VisitBioseqsInSep (sep, &prot_list, GetProtListCallback);
  prot_list = ValNodeSort (prot_list, SortProtByName);
    for (vnp = prot_list; vnp != NULL; vnp = vnp->next)
  {
      app = (AlphaProtPtr) vnp->data.ptrvalue;
      if (app == NULL) continue;
      if (prev_name != NULL && StringCmp (prev_name, app->prot_name) != 0) {
        fprintf (fp, "\n");
      }

    sep = SeqMgrGetSeqEntryForData (app->bsp);
    SeqEntryFastaStreamEx (sep, fp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, 70, 0, 0, FALSE, TRUE, FALSE, TRUE, TRUE);

    prev_name = app->prot_name;
    }
    for (vnp = prot_list; vnp != NULL; vnp = vnp->next)
    {
      app = (AlphaProtPtr) vnp->data.ptrvalue;
      if (app == NULL) continue;
      app->prot_name = MemFree (app->prot_name);
      vnp->data.ptrvalue = MemFree (app);
    }
  prot_list = ValNodeFree (prot_list);        
}

Int2 LIBCALLBACK VSMFastaSortedProtSave (Pointer data)
{
  OMProcControlPtr ompcp;
    Char filename[255];
    FILE * fp;
    ValNode vn;
    SeqEntryPtr sep = NULL;

  ompcp = (OMProcControlPtr)data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  
    switch(ompcp->input_itemtype)
    {
        case OBJ_SEQENTRY:
        case OBJ_BIOSEQ:
        case OBJ_BIOSEQSET:
            break;
        default:
            ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can only write Seq-entry, Bioseq, or Bioseq-set");
            return OM_MSG_RET_ERROR;
    }
    if (ompcp->input_choicetype == OBJ_SEQENTRY)
        sep = (SeqEntryPtr)(ompcp->input_choice);
    else
    {
        vn.next = NULL;
        vn.data.ptrvalue = ompcp->input_data;
        if (ompcp->input_itemtype == OBJ_BIOSEQ)
            vn.choice = 1;
        else
            vn.choice = 2;
        sep = &vn;
    }
        
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
#ifdef WIN_MAC
        fp = FileOpen (filename, "r");
        if (fp != NULL) {
            FileClose (fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif
        fp = FileOpen(filename, "w");
        
    WriteSortedProteinsToFile (fp, sep);

    FileClose(fp);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;  
}


NLM_EXTERN void ViewSortedProteins (SeqEntryPtr sep)
{
  Char         path [PATH_MAX];
  FILE        *fp;

  if (sep == NULL) return;
  TmpNam (path);
  fp = FileOpen(path, "w");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", path);
    return;
  }
        
  WriteSortedProteinsToFile (fp, sep);

  FileClose(fp);

  LaunchGeneralTextViewer (path, "Sorted Proteins");
  FileRemove (path);
}


typedef struct featuretableexport
{
  FILE *fp;
  Boolean show_nucs;
  Boolean show_prots;
  Boolean hide_sources;
  Boolean export_only_selected;
  Boolean suppress_protein_ids;
} FeatureTableData, PNTR FeatureTablePtr;

static void ExciseProteinIDLine (CharPtr line)
{
  CharPtr protein_id_line_start = NULL, protein_id_line_end;
  
  if (StringHasNoText (line))
  {
    return;
  }
  
  protein_id_line_start = StringStr (line, "\n\t\t\tprotein_id\t");
  if (protein_id_line_start == NULL)
  {
    return;
  }
  protein_id_line_end = StringChr (protein_id_line_start + 1, '\n');
  if (protein_id_line_end == NULL)
  {
    return;
  }
  
  while (*protein_id_line_end != 0)
  {
    *protein_id_line_start = *protein_id_line_end;
    protein_id_line_start ++;
    protein_id_line_end ++;
  }
  *protein_id_line_start = 0;
}

static Boolean IsBaseBlockFeatureSelected (BaseBlockPtr bbp)
{
  SelStructPtr     ssp = NULL;
  Boolean          rval = FALSE;
  
  if (bbp == NULL 
      || (bbp->blocktype != FEATURE_BLOCK && bbp->blocktype != SOURCEFEAT_BLOCK))
  {
    return FALSE;
  }
  ssp = ObjMgrGetSelected ();
  while (ssp != NULL)
  {
      if (ssp->entityID == bbp->entityID
          && ssp->itemtype == bbp->itemtype
            && ssp->itemID == bbp->itemID)
      {
        rval = TRUE;
      }
    
    ssp = ssp->next;    
  }
  return rval;
}

static Boolean BioseqHasSelectedFeatures (Asn2gbJobPtr ajp, Boolean hide_sources)
{
  Int4            index;
  BaseBlockPtr    bbp;

  if (ajp == NULL)
  {
    return FALSE;
  }

  for (index = 0; index < ajp->numParagraphs; index++) 
  {
    bbp = ajp->paragraphArray [index];
    if ((bbp->blocktype == FEATURE_BLOCK 
          || (bbp->blocktype == SOURCEFEAT_BLOCK && ! hide_sources))
        && IsBaseBlockFeatureSelected (bbp))
    {
      return TRUE;
    }
  }
  return FALSE;
}

static void VSMExportFeatureTableBioseqCallback (BioseqPtr bsp, Pointer userdata)

{
  FeatureTablePtr ftp;
  CstType         custom_flags = 0;
  Asn2gbJobPtr    ajp;
  BaseBlockPtr    bbp;
  XtraBlock       extra;
  Int4            index;
  CharPtr         string;
  
  if (bsp == NULL || userdata == NULL) return;
  
  ftp = (FeatureTablePtr) userdata;
  if (ftp->fp == NULL) return;
  if (!ftp->show_nucs && ISA_na (bsp->mol))
  {
    return;
  }
  if (!ftp->show_prots && ISA_aa (bsp->mol))
  {
    return;
  }
  if (ftp->hide_sources)
  {
    custom_flags |= HIDE_SOURCE_FEATS;
  }
  MemSet ((Pointer) &extra, 0, sizeof (XtraBlock));
  ajp = asn2gnbk_setup (bsp, NULL, NULL, FTABLE_FMT, DUMP_MODE, NORMAL_STYLE,
                           0, 0, custom_flags, &extra);
                           
  if (ftp->export_only_selected 
      && ! BioseqHasSelectedFeatures (ajp, ftp->hide_sources))
  {
    /* nothing to export */
  }                    
  else if (ajp != NULL) 
  {
    for (index = 0; index < ajp->numParagraphs; index++) 
    {
      bbp = ajp->paragraphArray [index];
      if (bbp->blocktype == FEATURE_BLOCK)
      {
        if (!ftp->export_only_selected || IsBaseBlockFeatureSelected (bbp))
        {
          string = asn2gnbk_format (ajp, (Int4) index);
          if (ftp->suppress_protein_ids)
          {
            ExciseProteinIDLine (string);
          }
          fprintf (ftp->fp, "%s", string);
          MemFree (string);
        }
      }
      else if (bbp->blocktype == SOURCEFEAT_BLOCK)
      {
        if (!ftp->hide_sources 
            && (!ftp->export_only_selected || IsBaseBlockFeatureSelected (bbp)))
        {
          string = asn2gnbk_format (ajp, (Int4) index);
          fprintf (ftp->fp, "%s", string);
          MemFree (string);
        }
      }
      else
      {
        string = asn2gnbk_format (ajp, (Int4) index);
        fprintf (ftp->fp, "%s", string);
        MemFree (string);
      }
    }
  }
  asn2gnbk_cleanup (ajp);
}

static Int4 FindNecessaryBioseqLength (SeqFeatPtr sfp)
{
  Int4 len = 0;
  
  while (sfp != NULL) {
    len = MAX (SeqLocStop (sfp->location), len);
    len = MAX (SeqLocStart (sfp->location), len);
    sfp = sfp->next;
  }
  return len + 1;
}

extern void ExportSeqAnnotFeatureTable (FILE *fp, SeqAnnotPtr sap) 
{
  FeatureTableData ftd;
  BioseqPtr fake_bsp;
  SeqFeatPtr first_feat;
  SeqEntryPtr sep;
  
  if (fp == NULL || sap == NULL) return;

  /* create fake bioseq to hold annotation */
  fake_bsp = BioseqNew();
  fake_bsp->annot = sap;
  
  /* create SeqEntry for temporary bioseq to live in */
  sep = SeqEntryNew ();
  sep->choice = 1;
  sep->data.ptrvalue = fake_bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) fake_bsp, sep);  

  ftd.export_only_selected = FALSE;
  ftd.fp = fp;
  ftd.show_nucs = TRUE;
  ftd.show_prots = TRUE;
  ftd.suppress_protein_ids = FALSE;

  /* show sources if there are any */
  first_feat = (SeqFeatPtr) sap->data;
  while (first_feat != NULL && first_feat->data.choice != SEQFEAT_BIOSRC) {
    first_feat = first_feat->next;
  }
  if (first_feat != NULL) {
    ftd.hide_sources = FALSE;
  } else {
    ftd.hide_sources = TRUE;
  }
      
  first_feat = (SeqFeatPtr) sap->data;
  fake_bsp->id = SeqIdDup (SeqLocId (first_feat->location));
  if (first_feat->data.choice == SEQFEAT_PROT) {
    fake_bsp->mol = Seq_mol_aa;
  } else {
    fake_bsp->mol = Seq_mol_dna;
  }
  fake_bsp->repr = Seq_repr_raw;
  fake_bsp->length = FindNecessaryBioseqLength (first_feat);
  
  VisitBioseqsInSep (sep, &ftd, VSMExportFeatureTableBioseqCallback);

  fake_bsp->annot = NULL;
  fake_bsp->idx.deleteme = TRUE;
  DeleteMarkedObjects (fake_bsp->idx.entityID, 0, NULL);
 
}

static Int2 
VSMExportFeatureTable 
(OMProcControlPtr ompcp, 
 Boolean show_nucs, 
 Boolean show_prots,
 Boolean hide_sources,
 Boolean suppress_protein_ids,
 Boolean export_only_selected)

{
    Char filename[255];
    ValNode vn;
    SeqEntryPtr sep = NULL;
    FeatureTableData ftd;
    SeqAnnotPtr sap = NULL;
    BioseqPtr fake_bsp = NULL;
    SeqFeatPtr first_feat = NULL;
    
    switch(ompcp->input_itemtype)
    {
        case OBJ_SEQENTRY:
        case OBJ_BIOSEQ:
        case OBJ_BIOSEQSET:
      case OBJ_SEQFEAT:
            break;
        case OBJ_SEQANNOT:
            sap = (SeqAnnotPtr) ompcp->input_data;
            if (sap== NULL) {
              ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can't write NULL Seq-annot");
              return OM_MSG_RET_ERROR;
            } else if (sap->type != 1) {
                ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can only write Feature Table Seq-annot");
              return OM_MSG_RET_ERROR;
            } else if (sap->data == NULL) {
                ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can't write empty Feature Table Seq-annot");
              return OM_MSG_RET_ERROR;
            }
            break;
        default:
            ErrPostEx(SEV_ERROR, 0,0,"ToFasta: Can only write Seq-entry, Feature Table Seq-annot, Bioseq, or Bioseq-set");
            return OM_MSG_RET_ERROR;
    }
    
    ftd.show_nucs = show_nucs;
    ftd.show_prots = show_prots;
    ftd.hide_sources = hide_sources;
    ftd.suppress_protein_ids = suppress_protein_ids;
    ftd.export_only_selected = export_only_selected;
    
  if (ompcp->input_itemtype == OBJ_SEQFEAT)
  {
    sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  }
  else if (ompcp->input_choicetype == OBJ_SEQENTRY)
  {
        sep = (SeqEntryPtr)(ompcp->input_choice);
  }
  else if (ompcp->input_itemtype == OBJ_SEQANNOT) {
    sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
    if (sep == NULL) {
      if (sap != NULL) {
        fake_bsp = BioseqNew();
        first_feat = (SeqFeatPtr) sap->data;
        fake_bsp->id = SeqIdDup (SeqLocId (first_feat->location));
        fake_bsp->annot = sap;
        if (first_feat->data.choice == SEQFEAT_PROT) {
          fake_bsp->mol = Seq_mol_aa;
        } else {
          fake_bsp->mol = Seq_mol_dna;
        }
        fake_bsp->repr = Seq_repr_raw;
        fake_bsp->length = FindNecessaryBioseqLength (first_feat);

        /* create SeqEntry for temporary bioseq to live in */
        sep = SeqEntryNew ();
        sep->choice = 1;
        sep->data.ptrvalue = fake_bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) fake_bsp, sep);
      }
    }
  }
    else
    {
        vn.next = NULL;
        vn.data.ptrvalue = ompcp->input_data;
        if (ompcp->input_itemtype == OBJ_BIOSEQ)
            vn.choice = 1;
        else
            vn.choice = 2;
        sep = &vn;
    }
        
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
#ifdef WIN_MAC
        ftd.fp = FileOpen (filename, "r");
        if (ftd.fp != NULL) {
            FileClose (ftd.fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif
        ftd.fp = FileOpen(filename, "w");
        VisitBioseqsInSep (sep, &ftd, VSMExportFeatureTableBioseqCallback);
        FileClose(ftd.fp);
        ArrowCursor();
    }
    
    if (fake_bsp != NULL) {
      fake_bsp->annot = NULL;
      fake_bsp->idx.deleteme = TRUE;
      DeleteMarkedObjects (fake_bsp->idx.entityID, 0, NULL);
    }
    
    return OM_MSG_RET_DONE;  
}

Int2 LIBCALLBACK VSMExportNucFeatureTable ( Pointer data )
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, TRUE, FALSE, FALSE, FALSE, FALSE);
}

Int2 LIBCALLBACK VSMExportNucFeatureTableSuppressProteinIDs ( Pointer data )
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, TRUE, FALSE, FALSE, TRUE, FALSE);
}

Int2 LIBCALLBACK VSMExportNucFeatureTableWithoutSources ( Pointer data )
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, TRUE, FALSE, TRUE, FALSE, FALSE);
}

Int2 LIBCALLBACK VSMExportNucFeatureTableWithoutSourcesSuppressProteinIDs ( Pointer data )
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, TRUE, FALSE, TRUE, TRUE, FALSE);
}

Int2 LIBCALLBACK VSMExportNucFeatureTableSelectedFeatures ( Pointer data )
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, TRUE, FALSE, FALSE, FALSE, TRUE);
}

Int2 LIBCALLBACK VSMExportNucFeatureTableSelectedFeaturesSuppressProteinIDs ( Pointer data )
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, TRUE, FALSE, FALSE, TRUE, TRUE);
}

Int2 LIBCALLBACK VSMExportProteinFeatureTable (Pointer data)
{
  return VSMExportFeatureTable ((OMProcControlPtr)data, FALSE, TRUE, FALSE, FALSE, FALSE);
}

typedef struct selectedsave
{
    AsnIoPtr     aip;
  SelStructPtr ssp;    
    ObjMgrPtr    omp;
} SelectedSaveData, PNTR SelectedSavePtr;

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));
static Boolean SaveOneSelectedItem (GatherObjectPtr gop)

{
  SelectedSavePtr ssp;
  SelStructPtr    sel;
    ObjMgrTypePtr   omtp;

  if (gop == NULL || gop->dataptr == NULL) return TRUE;
  ssp = (SelectedSavePtr) gop->userdata;
  if (ssp == NULL || ssp->aip == NULL || ssp->ssp == NULL) return TRUE;

  sel = ssp->ssp;
  while (sel != NULL 
         && (sel->entityID != gop->entityID 
             || sel->itemtype != gop->itemtype 
             || sel->itemID != gop->itemID))
  {
    sel = sel->next;
  }

  if (sel == NULL) return TRUE;
   
     omtp = ObjMgrTypeFind(ssp->omp, sel->itemtype, NULL, NULL);
    if (omtp == NULL)
    {
        ErrPostEx(SEV_ERROR,0,0,"Can't locate type record for [%d]", (int)sel->itemtype);
        return TRUE;
    }    
        
  (*(omtp->asnwrite))(gop->dataptr, ssp->aip, NULL);
  AsnPrintNewLine (ssp->aip);
  AsnIoFlush (ssp->aip);
 
  return TRUE;
}

static void SaveSeqLocEntity (Uint2 entityID, SelectedSavePtr ssp)
{
  ObjMgrDataPtr  omdp;
    ObjMgrTypePtr  omtp;

  if (entityID == 0 || ssp == NULL)
  {
    return;
  }
  
  omdp = ObjMgrGetDataStruct (ssp->omp, entityID);
  if (omdp == NULL) return;
  if (omdp->choicetype == OBJ_SEQENTRY || omdp->datatype != OBJ_SEQLOC) 
  {
    return; /* not seqloc */  
  }
  
     omtp = ObjMgrTypeFind(ssp->omp, OBJ_SEQLOC, NULL, NULL);
    if (omtp == NULL)
    {
        ErrPostEx(SEV_ERROR,0,0,"Can't locate type record for [%d]", (int)OBJ_SEQLOC);
        return;
    }    
        
  (*(omtp->asnwrite))(omdp->dataptr, ssp->aip, NULL);
  AsnPrintNewLine (ssp->aip);
  AsnIoFlush (ssp->aip);
}

static Int2 LIBCALLBACK VSMGenericAsnSave (OMProcControlPtr ompcp, CharPtr mode )
{
    Char filename[255];
    SelStructPtr  ssp, sel;
#ifdef WIN_MAC
    FILE * fp;
#endif
  ValNodePtr entity_list = NULL, vnp;
  SelectedSaveData ssd;

  ssp = ObjMgrGetSelected();
  if (ssp == NULL)
    {
        return OM_MSG_RET_ERROR;
    }
    
    for (sel = ssp; sel != NULL; sel = sel->next)
    {
      for (vnp = entity_list;
           vnp != NULL && vnp->data.intvalue != sel->entityID;
           vnp = vnp->next)
      {}
      if (vnp == NULL)
      {
        ValNodeAddInt (&entity_list, 0, sel->entityID);
      }
    }

    ssd.omp = ObjMgrGet();

  /* get file name to use */    
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
#ifdef WIN_MAC
        fp = FileOpen (filename, "r");
        if (fp != NULL) {
            FileClose (fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif

        ssd.aip = AsnIoOpen(filename, mode);
        ssd.ssp = ssp;
    
      for (vnp = entity_list; vnp != NULL; vnp = vnp->next)
      {
        SaveSeqLocEntity (vnp->data.intvalue, &ssd);
      GatherObjectsInEntity (vnp->data.intvalue, 0, NULL, SaveOneSelectedItem, (Pointer) &ssd, NULL);
      }

    ValNodeFree (entity_list);
        AsnIoClose(ssd.aip);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;
}


Int2 LIBCALLBACK VSMGenericTextAsnSave ( Pointer data )
{
    return VSMGenericAsnSave((OMProcControlPtr)data, "w");
}

Int2 LIBCALLBACK VSMGenericBinAsnSave ( Pointer data )
{
    return VSMGenericAsnSave((OMProcControlPtr)data, "wb");
}


typedef struct setsave {
  SelStructPtr sel;
  CharPtr file_base;
  Int4 file_num;
} SetSaveData, PNTR SetSavePtr;

static Boolean SaveSetsInOneSelectedSet (GatherObjectPtr gop)

{
  SetSavePtr      ssp;
  SelStructPtr    sel;
  BioseqSetPtr    bssp;
  CharPtr         filename, file_fmt = "%s_%d.sqn";
  SeqEntryPtr     sep;
  AsnIoPtr        aip;
#ifdef WIN_MAC
  FILE *fp;
#endif

  if (gop == NULL || gop->dataptr == NULL || gop->itemtype != OBJ_BIOSEQSET) return TRUE;
  ssp = (SetSavePtr) gop->userdata;
  if (ssp == NULL || StringHasNoText (ssp->file_base)) return TRUE;

  sel = ssp->sel;
  while (sel != NULL 
         && (sel->entityID != gop->entityID 
             || sel->itemtype != gop->itemtype 
             || sel->itemID != gop->itemID))
  {
    sel = sel->next;
  }

  if (sel == NULL) return TRUE;

  bssp = gop->dataptr;
   
  filename = (CharPtr) MemNew (sizeof (Char) + (StringLen (ssp->file_base) + StringLen (file_fmt) + 15));
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    sprintf (filename, file_fmt, ssp->file_base, ssp->file_num);
    ssp->file_num++;
#ifdef WIN_MAC
        fp = FileOpen (filename, "r");
        if (fp != NULL) {
            FileClose (fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif

        aip = AsnIoOpen(filename, "w");
        SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);
  }
  filename = MemFree (filename);
  return TRUE;
}


Int2 LIBCALLBACK VSMSaveSetsAsFiles (Pointer data)
{
    Char filename[255];
    SelStructPtr  ssp, sel;
  ValNodePtr entity_list = NULL, vnp;
  SetSaveData ssd;
  OMProcControlPtr ompcp;

  ompcp = (OMProcControlPtr) data;
  ssp = ObjMgrGetSelected();
  if (ssp == NULL)
    {
    Message (MSG_ERROR, "Nothing selected!");
        return OM_MSG_RET_ERROR;
    }

  /* get file name to use */    
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
    ssd.file_base = filename;
    ssd.file_num = 1;

    ssd.sel = ssp;

    /* get list of entity IDs */
      for (sel = ssp; sel != NULL; sel = sel->next)
      {
      if (sel->itemtype != OBJ_BIOSEQSET) {
        continue;
      }
        for (vnp = entity_list;
            vnp != NULL && vnp->data.intvalue != sel->entityID;
            vnp = vnp->next)
        {}
        if (vnp == NULL)
        {
          ValNodeAddInt (&entity_list, 0, sel->entityID);
        GatherObjectsInEntity (sel->entityID, 0, NULL, SaveSetsInOneSelectedSet, (Pointer) &ssd, NULL);
        }
      }
    

    ValNodeFree (entity_list);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;
}


static Boolean AddToSaveList (GatherContextPtr gcp)

{
  ValNodePtr PNTR list;

  list = (ValNodePtr PNTR) gcp->userdata;
  if (list == NULL) return TRUE;
  ValNodeAddPointer (list, gcp->thistype, gcp->thisitem);
  return TRUE;
}


typedef struct savesetsdp {
  AsnIoPtr aip;
  Boolean  already_have_molinfo;  
} SaveSetSdpData, PNTR SaveSetSdpPtr;


static void SaveSetDescriptors (SeqDescrPtr sdp, Pointer userdata)
{
  SaveSetSdpPtr sp;
  MolInfoPtr    mip;

  if (sdp == NULL || userdata == NULL) return;
  if (sdp->choice != Seq_descr_source
      && sdp->choice != Seq_descr_pub
      && sdp->choice != Seq_descr_molinfo
      && sdp->choice != Seq_descr_comment) return;

  sp = (SaveSetSdpPtr) userdata;
  if (sdp->choice == Seq_descr_molinfo
      && (mip = (MolInfoPtr) sdp->data.ptrvalue) != NULL)
  {
    if (mip->biomol == MOLECULE_TYPE_PEPTIDE || sp->already_have_molinfo) 
    {
      return;
    }
    else
    {
      sp->already_have_molinfo = TRUE;
    }
  }
  SeqDescAsnWrite (sdp, sp->aip, NULL);
  AsnPrintNewLine (sp->aip);
  AsnIoFlush (sp->aip);
}

Int2 LIBCALLBACK VSMDescriptorAsnSave (Pointer data)
{
  OMProcControlPtr ompcp;
  Char filename[255];
  SelStructPtr  ssp, sel;
  SeqDescrPtr   sdp;
  BioseqPtr     bsp;
  SeqEntryPtr   sep;
  SeqMgrDescContext dcontext;
#ifdef WIN_MAC
  FILE * fp;
#endif
  ValNodePtr vnp;
  SaveSetSdpData sd;
  ValNodePtr obj_list = NULL;

  ompcp = (OMProcControlPtr)data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;

  sd.already_have_molinfo = FALSE;
  ssp = ObjMgrGetSelected();
  if (ssp == NULL)
    {
    Message (MSG_ERROR, "You must select a sequence or set from which descriptors should be saved");
    return OM_MSG_RET_DONE;
  } else {
      for (sel = ssp; sel != NULL; sel = sel->next)
      {
      GatherItem (sel->entityID, sel->itemID, sel->itemtype,
                (Pointer) &obj_list, AddToSaveList);
    }
  }

  /* get file name to use */    
    filename[0] = '\0';
    if (GetOutputFileName(filename, (size_t)254, NULL))
    {
        WatchCursor();
#ifdef WIN_MAC
        fp = FileOpen (filename, "r");
        if (fp != NULL) {
            FileClose (fp);
        } else {
            FileCreate (filename, "TEXT", "ttxt");
        }
#endif

        sd.aip = AsnIoOpen(filename, "w");
      for (vnp = obj_list; vnp != NULL; vnp = vnp->next)
      {
      switch (vnp->choice) {
        case OBJ_SEQDESC:
          sdp = (SeqDescrPtr) vnp->data.ptrvalue;
          SeqDescAsnWrite (sdp, sd.aip, NULL);
          AsnPrintNewLine (sd.aip);
          AsnIoFlush (sd.aip);
          break;
        case OBJ_BIOSEQ:
          bsp = (BioseqPtr) vnp->data.ptrvalue;
          for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
               sdp != NULL;
               sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext)) {
            if (sdp->choice != Seq_descr_source
                && sdp->choice != Seq_descr_pub
                && sdp->choice != Seq_descr_molinfo
                && sdp->choice != Seq_descr_comment) continue;
            SeqDescAsnWrite (sdp, sd.aip, NULL);
            AsnPrintNewLine (sd.aip);
            AsnIoFlush (sd.aip);
          }
          break;
        case OBJ_BIOSEQSET:
          sep = SeqMgrGetSeqEntryForData (vnp->data.ptrvalue);
          VisitDescriptorsInSep (sep, &sd, SaveSetDescriptors);
          break;
      }              
    }   
        AsnIoClose(sd.aip);
        ArrowCursor();
    }
    
    return OM_MSG_RET_DONE;  
}


