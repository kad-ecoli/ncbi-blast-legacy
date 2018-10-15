/*  objfeat.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  objfeat.c
*
* Author:  James Ostell
*   
* Version Creation Date: 4/1/91
*
* $Revision: 6.71 $
*
* File Description:  Object manager for module NCBI-SeqFeat
*
* Modifications:  
* --------------------------------------------------------------------------
* Date       Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

/** for ErrPostEx() ****/

static char *this_module = "ncbiobj";
#define THIS_MODULE this_module
static char *this_file = __FILE__;
#define THIS_FILE this_file

/**********************/

#include <objfeat.h>           /* the features interface */
#include <asnfeat.h>        /* the AsnTool header */
#include <objmgr.h>
#include <objfdef.h>
#include <objseq.h>

static Boolean loaded = FALSE;

/*****************************************************************************
*
*   SeqFeat ObjMgr Routines
*
*****************************************************************************/
static Pointer LIBCALLBACK SeqFeatNewFunc (void)
{
    return (Pointer) SeqFeatNew();
}

static Pointer LIBCALLBACK SeqFeatFreeFunc (Pointer data)
{
    return (Pointer) SeqFeatFree ((SeqFeatPtr) data);
}

static Boolean LIBCALLBACK SeqFeatAsnWriteFunc (Pointer data, AsnIoPtr aip, AsnTypePtr atp)
{
    return SeqFeatAsnWrite((SeqFeatPtr)data, aip, atp);
}

static Pointer LIBCALLBACK SeqFeatAsnReadFunc (AsnIoPtr aip, AsnTypePtr atp)
{
    return (Pointer) SeqFeatAsnRead (aip, atp);
}

static Int2 LIBCALLBACK SeqFeatLabelFunc ( Pointer data, CharPtr buffer, Int2 buflen, Uint1 content)
{
    Int2 rsult = 0;

    switch (content)
    {
        case OM_LABEL_TYPE:
        case OM_LABEL_CONTENT:
        case OM_LABEL_BOTH:
            rsult = FeatDefLabel((SeqFeatPtr)data, buffer, buflen, content);
            break;
        case OM_LABEL_SUMMARY:
        default:
            rsult = FeatDefLabel((SeqFeatPtr)data, buffer, buflen, OM_LABEL_BOTH);
            break;
    }
    return rsult;
}

static Uint2 LIBCALLBACK SeqFeatSubTypeFunc (Pointer ptr)
{
    if (ptr == NULL)
        return 0;
    return (Uint2)FindFeatDefType((SeqFeatPtr)ptr);
}

/*****************************************************************************
*
*   SeqFeatAsnLoad()
*      requires SeqAsnLoad() to be called first
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatAsnLoad (void)
{
    if (loaded)
        return TRUE;
    loaded = TRUE;

    if (! GeneralAsnLoad())
    {
        loaded = FALSE;
        return FALSE;
    }
    if (! SeqLocAsnLoad())      /* SeqLoc refers to Seqfeat and vice versa */
    {
        loaded = FALSE;
        return FALSE;
    }
    if (! PubAsnLoad())
    {
        loaded = FALSE;
        return FALSE;
    }
    if ( ! AsnLoad())
    {
        loaded = FALSE;
        return FALSE;
    }

    ObjMgrTypeLoad(OBJ_SEQFEAT, "Seq-feat", "SeqFeat", "Sequence Feature",
        SEQ_FEAT, SeqFeatNewFunc, SeqFeatAsnReadFunc, SeqFeatAsnWriteFunc,
        SeqFeatFreeFunc, SeqFeatLabelFunc, SeqFeatSubTypeFunc);

    return TRUE;
}

/*****************************************************************************
*
*   SeqFeat Routines
*
*****************************************************************************/
/*****************************************************************************
*
*   SeqFeatNew()
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL SeqFeatNew (void)
{
    return (SeqFeatPtr)MemNew(sizeof(SeqFeat));
}

/*****************************************************************************
*
*   SeqFeatFree(sfp)
*       Frees one SeqFeat and associated data
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL SeqFeatFree (SeqFeatPtr sfp)
{
    ValNodePtr vnp, next=NULL;
  UserObjectPtr uop, unp = NULL;

    if (sfp == NULL)
        return (SeqFeatPtr)NULL;

    SeqFeatIdFree(&sfp->id);
    SeqFeatDataFree(&sfp->data);
    MemFree(sfp->title);
    MemFree(sfp->comment);
    UserObjectFree(sfp->ext);
    SeqLocFree(sfp->product);
    SeqLocFree(sfp->location);
    GBQualFree(sfp->qual);
    PubSetFree(sfp->cit);
    SeqFeatXrefFree(sfp->xref);
    for (vnp = sfp->dbxref; vnp != NULL; vnp = next)
    {
        next = vnp->next;
        DbtagFree((DbtagPtr)(vnp->data.ptrvalue));
        MemFree(vnp);
    }
    MemFree(sfp->except_text);
    for (vnp = sfp->ids; vnp != NULL; vnp = next)
    {
        next = vnp->next;
        SeqFeatIdFree((ChoicePtr)(vnp->data.ptrvalue));
        MemFree(vnp);
    }
    /*
    AsnGenericChoiceSeqOfFree (sfp->ids, (AsnOptFreeFunc) SeqFeatIdFree);
    */
    for (uop = sfp->exts; uop != NULL; uop = unp)
    {
        unp = uop->next;
        UserObjectFree(uop);
    }
    /*
    AsnGenericUserSeqOfFree (sfp->exts, (AsnOptFreeFunc) UserObjectFree);
    */
    SeqFeatSupportFree (sfp->support);

    ObjMgrDelete(OBJ_SEQFEAT, (Pointer)sfp);

    return (SeqFeatPtr)MemFree(sfp);
}

/*****************************************************************************
*
*   SeqFeatAsnWrite(sfp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqFeat ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatAsnWrite (SeqFeatPtr sfp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    GBQualPtr gbp;
    Boolean retval = FALSE;
    SeqFeatXrefPtr sfxp;
    ValNodePtr vnp;
  UserObjectPtr uop;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, SEQ_FEAT);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (sfp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)sfp))
        goto erret;

    if (sfp->id.choice)              /* feat-id present */
    {
        if (! SeqFeatIdAsnWrite(&sfp->id, aip, SEQ_FEAT_id))
            goto erret;
    }

    if (! SeqFeatDataAsnWrite(&sfp->data, aip, SEQ_FEAT_data)) goto erret;


    av.boolvalue = TRUE;
    if (sfp->partial)
    {
        if (! AsnWrite(aip, SEQ_FEAT_partial, &av)) goto erret;
    }
    if (sfp->excpt)
    {
        if (! AsnWrite(aip, SEQ_FEAT_except, &av)) goto erret;
    }
    if (sfp->comment != NULL)
    {
        av.ptrvalue = sfp->comment;
        if (! AsnWrite(aip, SEQ_FEAT_comment, &av)) goto erret;
    }
    if (sfp->product != NULL)
    {
        if (! SeqLocAsnWrite(sfp->product, aip, SEQ_FEAT_product))
            goto erret;
    }
    if (sfp->location != NULL)
    {
        if (! SeqLocAsnWrite(sfp->location, aip, SEQ_FEAT_location))
            goto erret;
    }
    if (sfp->qual != NULL)
    {                                  /* SEQUENCE OF */
        if (! AsnOpenStruct(aip, SEQ_FEAT_qual, (Pointer)sfp->qual))
            goto erret;
        gbp = sfp->qual;
        while (gbp != NULL)
        {
            if (! GBQualAsnWrite(gbp, aip, SEQ_FEAT_qual_E))
                goto erret;
            gbp = gbp->next;
        }
        if (! AsnCloseStruct(aip, SEQ_FEAT_qual, (Pointer)sfp->qual))
            goto erret;
    }
    if (sfp->title != NULL)
    {
        av.ptrvalue = (Pointer) sfp->title;
        if (! AsnWrite(aip, SEQ_FEAT_title, &av)) goto erret;
    }
    if (sfp->ext != NULL)
    {
        if (! UserObjectAsnWrite(sfp->ext, aip, SEQ_FEAT_ext))
            goto erret;
    }

    if (sfp->cit != NULL)
    {
        if (! PubSetAsnWrite(sfp->cit, aip, SEQ_FEAT_cit))
            goto erret;
    }

    if (sfp->exp_ev)
    {
        av.intvalue = (Int4)sfp->exp_ev;
        if (! AsnWrite(aip, SEQ_FEAT_exp_ev, &av)) goto erret;
    }

    if (sfp->xref != NULL)
    {
        if (! AsnOpenStruct(aip, SEQ_FEAT_xref, (Pointer)sfp->xref))
            goto erret;
        sfxp = sfp->xref;
        while (sfxp != NULL)
        {
            if (! SeqFeatXrefAsnWrite(sfxp, aip, SEQ_FEAT_xref_E))
                goto erret;
            sfxp = sfxp->next;
        }
        if (! AsnCloseStruct(aip, SEQ_FEAT_xref, (Pointer)sfp->xref))
            goto erret;
    }

    if (sfp->dbxref != NULL)
    {
        if (! AsnOpenStruct(aip, SEQ_FEAT_dbxref, (Pointer)sfp->dbxref))
            goto erret;
        vnp = sfp->dbxref;
        while (vnp != NULL)
        {
            if (! DbtagAsnWrite((DbtagPtr)(vnp->data.ptrvalue), aip, SEQ_FEAT_dbxref_E))
                goto erret;
            vnp = vnp->next;
        }
        if (! AsnCloseStruct(aip, SEQ_FEAT_dbxref, (Pointer)sfp->dbxref))
            goto erret;
    }
    if (sfp->pseudo)
    {
        av.boolvalue = TRUE;
        if ((aip->spec_version == 0) || (aip->spec_version >= 6))
        {
        if (! AsnWrite(aip, SEQ_FEAT_pseudo, &av)) goto erret;
        }
    }
    if (sfp->except_text != NULL)
    {
        av.ptrvalue = sfp->except_text;
        if ((aip->spec_version == 0) || (aip->spec_version >= 6))
        {
        if (! AsnWrite(aip, SEQ_FEAT_except_text, &av)) goto erret;
        }
    }
    if (sfp->ids != NULL)
    {
        if (! AsnOpenStruct(aip, SEQ_FEAT_ids, (Pointer)sfp->ids))
            goto erret;
        vnp = sfp->ids;
        while (vnp != NULL)
        {
            if (! SeqFeatIdAsnWrite((ChoicePtr)(vnp->data.ptrvalue), aip, SEQ_FEAT_ids_E))
                goto erret;
            vnp = vnp->next;
        }
        if (! AsnCloseStruct(aip, SEQ_FEAT_ids, (Pointer)sfp->ids))
            goto erret;
        /*
        AsnGenericChoiceSeqOfAsnWrite (sfp->ids, (AsnWriteFunc) SeqFeatIdAsnWrite, aip, SEQ_FEAT_ids, SEQ_FEAT_ids_E);
        */
    }
    if (sfp->exts != NULL)
    {
        if (! AsnOpenStruct(aip, SEQ_FEAT_exts, (Pointer)sfp->exts))
            goto erret;
        uop = sfp->exts;
        while (uop != NULL)
        {
            if (! UserObjectAsnWrite(uop, aip, SEQ_FEAT_exts_E))
                goto erret;
            uop = uop->next;
        }
        if (! AsnCloseStruct(aip, SEQ_FEAT_exts, (Pointer)sfp->exts))
            goto erret;
        /*
        AsnGenericUserSeqOfAsnWrite (sfp->exts, (AsnWriteFunc) UserObjectAsnWrite, aip, SEQ_FEAT_exts, SEQ_FEAT_exts_E);
        */
    }
    if (sfp->support != NULL)
    {
        if (! SeqFeatSupportAsnWrite(sfp->support, aip, SEQ_FEAT_support))
            goto erret;
    }
    if (! AsnCloseStruct(aip, atp, (Pointer)sfp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SeqFeatAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL SeqFeatAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldtype;
    SeqFeatPtr sfp;
    GBQualPtr gbq, gbqlast = NULL;
    SeqFeatXrefPtr sfxp, sfxplast = NULL;
    ValNodePtr vnp, vnplast = NULL;
    DbtagPtr dbtp;
  Choice id;
  ChoicePtr cp;
    UserObjectPtr  uop, last = NULL;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (SeqFeatPtr)NULL;
    }

    if (aip == NULL)
        return (SeqFeatPtr)NULL;

    if (orig == NULL)           /* SeqFeat ::= (self contained) */
        atp = AsnReadId(aip, amp, SEQ_FEAT);
    else
        atp = AsnLinkType(orig, SEQ_FEAT);    /* link in local tree */
    oldtype = atp;
    if (atp == NULL)
        return (SeqFeatPtr)NULL;

    sfp = SeqFeatNew();
    if (sfp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) != oldtype)
    {
        if (atp == NULL)
            goto erret;
        if (atp == SEQ_FEAT_id)
        {
            if (! SeqFeatIdAsnRead(aip, atp, &sfp->id))
                goto erret;
        }
        else if (atp == SEQ_FEAT_data)
        {
            if (! SeqFeatDataAsnRead(aip, atp, &sfp->data))
                goto erret;
        }
        else if (atp == SEQ_FEAT_product)
        {
            sfp->product = SeqLocAsnRead(aip, atp);
            if (sfp->product == NULL) goto erret;
        }
        else if (atp == SEQ_FEAT_location)
        {
            sfp->location = SeqLocAsnRead(aip, atp);
            if (sfp->location == NULL) goto erret;
        }
        else if (atp == SEQ_FEAT_qual_E)
        {
            gbq = GBQualAsnRead(aip, atp);
            if (gbq == NULL) goto erret;
            if (gbqlast == NULL)
                sfp->qual = gbq;
            else
                gbqlast->next = gbq;
            gbqlast = gbq;
        }
        else if (atp == SEQ_FEAT_ext)
        {
            sfp->ext = UserObjectAsnRead(aip, atp);
            if (sfp->ext == NULL)
                goto erret;
        }
        else if (atp == SEQ_FEAT_cit)
        {
            sfp->cit = PubSetAsnRead(aip, atp);
            if (sfp->cit == NULL)
                goto erret;
        }
        else if (atp == SEQ_FEAT_xref_E)
        {
            sfxp = SeqFeatXrefAsnRead(aip, atp);
            if (sfxp == NULL) goto erret;
            if (sfxplast == NULL)
                sfp->xref = sfxp;
            else
                sfxplast->next = sfxp;
            sfxplast = sfxp;
        }
        else if (atp == SEQ_FEAT_dbxref_E)
        {
            dbtp = DbtagAsnRead(aip, atp);
            if (dbtp == NULL) goto erret;
            vnp = ValNodeNew(vnplast);
            if (vnplast == NULL)
                sfp->dbxref = vnp;
            vnplast = vnp;
            vnp->data.ptrvalue = (Pointer)dbtp;
        }
        else if (atp == SEQ_FEAT_ids)
        {
        atp = AsnLinkType (SEQ_FEAT_ids, FEAT_ID);
        if (atp == NULL) goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the SEQUENCE */
            while ((atp = AsnReadId(aip, amp, atp)) == SEQ_FEAT_ids_E) {
                if (! SeqFeatIdAsnRead(aip, atp, &id))
                    goto erret;
                cp = (ChoicePtr) MemNew (sizeof (Choice));
                if (cp == NULL) goto erret;
                cp->choice = id.choice;
                if (id.choice) {
                    cp->value.intvalue = id.value.intvalue;
                } else {
                    cp->value.ptrvalue = id.value.ptrvalue;
                }
                ValNodeAddPointer (&(sfp->ids), 0, (Pointer) cp);
            }
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the END STRUCT */
        AsnUnlinkType (SEQ_FEAT_ids);
            /*
            sfp->ids = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqFeatIdAsnRead, (AsnOptFreeFunc) SeqFeatIdFree);
            if (isError && sfp->ids == NULL) {
                goto erret;
            }
            */
        }
        else if (atp == SEQ_FEAT_exts)
        {
            atp = AsnLinkType (SEQ_FEAT_exts, FEAT_ID);
            if (atp == NULL) goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the SEQUENCE */
                while ((atp = AsnReadId(aip, amp, atp)) == SEQ_FEAT_exts_E) {
                uop = UserObjectAsnRead (aip, atp);
                if (uop == NULL) goto erret;
                if (last != NULL) {
                  last->next = uop;
                }
                if (sfp->exts == NULL) {
                  sfp->exts = uop;
                }
                last = uop;
            }
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the END STRUCT */
            AsnUnlinkType (SEQ_FEAT_exts);
            /*
            sfp->exts = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) UserObjectAsnRead, (AsnOptFreeFunc) UserObjectFree);
            if (isError && sfp->exts == NULL) {
                goto erret;
            }
            */
        }
        else if (atp == SEQ_FEAT_support)
        {
            sfp->support = SeqFeatSupportAsnRead(aip, atp);
            if (sfp->support == NULL)
                goto erret;
        }
        else
        {
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            if (atp == SEQ_FEAT_partial)
                sfp->partial = av.boolvalue;
            else if (atp == SEQ_FEAT_except)
                sfp->excpt = av.boolvalue;
            else if (atp == SEQ_FEAT_comment)
                sfp->comment = (CharPtr) av.ptrvalue;
            else if (atp == SEQ_FEAT_title)
                sfp->title = (CharPtr) av.ptrvalue;
            else if (atp == SEQ_FEAT_exp_ev)
                sfp->exp_ev = (Uint1) av.intvalue;
            else if (atp == SEQ_FEAT_pseudo)
                sfp->pseudo = av.boolvalue;
            else if (atp == SEQ_FEAT_except_text)
                sfp->except_text = (CharPtr)(av.ptrvalue);
        }
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return sfp;
erret:
    aip->io_failure = TRUE;
    sfp = SeqFeatFree(sfp);
    goto ret;
}

/*****************************************************************************
*
*   SeqFeatSetAsnWrite(sfp, aip, set, element)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatSetAsnWrite (SeqFeatPtr sfp, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    return SeqFeatSetAsnWriteExtra(sfp, aip, set, element, (ValNodePtr) NULL);
}

/*****************************************************************************
*
*   SeqFeatSetAsnWriteExtra(sfp, aip, set, element, extras)
*      this is to accomodate extra features from SeqEntryAsnOut()
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatSetAsnWriteExtra (SeqFeatPtr sfp, AsnIoPtr aip,
                    AsnTypePtr set, AsnTypePtr element, ValNodePtr extras)
{
    AsnTypePtr atp;
    SeqFeatPtr oldsfp, tsfp;
    Boolean retval = FALSE;
    Int2 ctr = 0;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(element, SEQ_FEAT);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    oldsfp = sfp;
    if (extras != NULL)   /* this is an extra set of feats */
        oldsfp = (SeqFeatPtr)(extras->data.ptrvalue);

    if (! AsnOpenStruct(aip, set, (Pointer)oldsfp))
        goto erret;

    while (sfp != NULL)
    {
        if (! SeqFeatAsnWrite(sfp, aip, atp))
            goto erret;
        sfp = sfp->next;
        ctr++;
        if (ctr == 20)
        {
            if (! ProgMon("Write SeqFeat"))
                goto erret;
            ctr = 0;
        }

    }

    while (extras != NULL)
    {
    tsfp = (SeqFeatPtr)(extras->data.ptrvalue);
    if (! SeqFeatAsnWrite(tsfp, aip, atp))
        goto erret;
    extras = extras->next;
    }

    if (! AsnCloseStruct(aip, set, (Pointer)oldsfp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(element);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SeqFeatSetAsnRead(aip, set, element)
*
*****************************************************************************/
NLM_EXTERN SeqFeatPtr LIBCALL SeqFeatSetAsnRead (AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    DataVal av;
    AsnTypePtr atp;
    SeqFeatPtr sfp, curr = NULL, first = NULL;
    Int2 ctr = 0;

    if (aip == NULL)
        return (SeqFeatPtr)NULL;

    AsnLinkType(element, SEQ_FEAT);    /* link in local tree */
    atp = set;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) == element)
    {
        sfp = SeqFeatAsnRead(aip, atp);
        if (sfp == NULL)
            goto erret;
        if (first == NULL)
            first = sfp;
        else
            curr->next = sfp;
        curr = sfp;
        ctr++;
        if (ctr == 20)
        {
            if (! ProgMon("Read SeqFeat"))
                goto erret;
            ctr = 0;
        }

    }
    if (atp == NULL)
        goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    if (first == NULL)
        ErrPostEx(SEV_WARNING, 0, 0, "Empty SET OF Seq-feat. line %ld", (long) aip->linenumber);
ret:
    AsnUnlinkType(element);       /* unlink local tree */
    return first;
erret:
    while (first != NULL)
    {
        curr = first;
        first = curr->next;
        SeqFeatFree(curr);
    }
    goto ret;
}

/*****************************************************************************
*
*   SeqFeatDataFree(cp)
*
*****************************************************************************/
NLM_EXTERN void LIBCALL SeqFeatDataFree(ChoicePtr cp)
{
    Pointer pnt;

    pnt = cp->value.ptrvalue;
    switch (cp->choice)
    {
        case SEQFEAT_GENE :
            GeneRefFree((GeneRefPtr)pnt);
            break;
        case SEQFEAT_ORG :
            OrgRefFree((OrgRefPtr)pnt);
            break;
        case SEQFEAT_CDREGION :
            CdRegionFree((CdRegionPtr)pnt);
            break;
        case SEQFEAT_PROT :
            ProtRefFree((ProtRefPtr)pnt);
            break;
        case SEQFEAT_RNA :
            RnaRefFree((RnaRefPtr)pnt);
            break;
        case SEQFEAT_PUB :
            PubdescFree((PubdescPtr)pnt);
            break;
        case SEQFEAT_SEQ :         /* seq */
            SeqLocFree((ValNodePtr)pnt);
            break;
        case SEQFEAT_IMP :
            ImpFeatFree((ImpFeatPtr)pnt);
            break;
        case SEQFEAT_REGION:        /* region */
        case SEQFEAT_NON_STD_RESIDUE :       /* non-std residue */
        case SEQFEAT_HET:       /* heterogen */
            MemFree(pnt);
            break;
        case SEQFEAT_COMMENT:       /* comment - nothing */
        case SEQFEAT_BOND:       /* bond - int */
        case SEQFEAT_SITE:       /* site - int */
        case SEQFEAT_PSEC_STR:       /* psec_struct - int */
            break;
        case SEQFEAT_RSITE:
            RsiteRefFree((RsiteRefPtr)pnt);
            break;
        case SEQFEAT_USER:
            UserObjectFree((UserObjectPtr)pnt);
            break;
        case SEQFEAT_TXINIT:
            TxinitFree((TxinitPtr)pnt);
            break;
        case SEQFEAT_NUM:
            NumberingFree((NumberingPtr)pnt);
            break;
        case SEQFEAT_BIOSRC:
            BioSourceFree((BioSourcePtr)pnt);
            break;
        case SEQFEAT_CLONEREF:
            CloneRefFree((CloneRefPtr)pnt);
            break;
        case SEQFEAT_VARIATIONREF:
            VariationRefFree((VariationRefPtr)pnt);
            break;
    }
    return;
}

/*****************************************************************************
*
*   SeqFeatDataAsnWrite(cp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqFeatData ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatDataAsnWrite (ChoicePtr cp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    Pointer pnt;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, SEQFEATDATA);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (cp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    pnt = cp->value.ptrvalue;
    av.ptrvalue = pnt;
    if (! AsnWriteChoice(aip, atp, (Int2)cp->choice, &av)) goto erret;

    switch (cp->choice)
    {
        case SEQFEAT_GENE :
            if (! GeneRefAsnWrite((GeneRefPtr)pnt, aip, SEQFEATDATA_gene))
                goto erret;
            break;
        case SEQFEAT_ORG :
            if (! OrgRefAsnWrite((OrgRefPtr)pnt, aip, SEQFEATDATA_org))
                goto erret;
            break;
        case SEQFEAT_CDREGION :
            if (! CdRegionAsnWrite((CdRegionPtr)pnt, aip, SEQFEATDATA_cdregion))
                goto erret;
            break;
        case SEQFEAT_PROT :
            if (! ProtRefAsnWrite((ProtRefPtr)pnt, aip, SEQFEATDATA_prot))
                goto erret;
            break;
        case SEQFEAT_RNA :
            if (! RnaRefAsnWrite((RnaRefPtr)pnt, aip, SEQFEATDATA_rna))
                goto erret;
            break;
        case SEQFEAT_PUB:
            if (! PubdescAsnWrite((PubdescPtr)pnt, aip, SEQFEATDATA_pub))
                goto erret;
            break;
        case SEQFEAT_SEQ:         /* seq reference */
            if (! SeqLocAsnWrite((ValNodePtr)pnt, aip, SEQFEATDATA_seq))
                goto erret;
            break;
        case SEQFEAT_IMP:
            if (! ImpFeatAsnWrite((ImpFeatPtr)pnt, aip, SEQFEATDATA_imp))
                goto erret;
            break;
        case SEQFEAT_REGION:        /* region */
            av.ptrvalue = pnt;
            if (! AsnWrite(aip, SEQFEATDATA_region, &av)) goto erret;
            break;
        case SEQFEAT_COMMENT:       /* comment - nothing */
            if (! AsnWrite(aip, SEQFEATDATA_comment, &av)) goto erret;
            break;
        case SEQFEAT_BOND:       /* bond - int */
            av.intvalue = cp->value.intvalue;
            if (! AsnWrite(aip, SEQFEATDATA_bond, &av)) goto erret;
            break;
        case SEQFEAT_SITE:       /* site - int */
            av.intvalue = cp->value.intvalue;
            if (! AsnWrite(aip, SEQFEATDATA_site, &av)) goto erret;
            break;
        case SEQFEAT_RSITE:       /* rsite-ref */
            if (! RsiteRefAsnWrite((RsiteRefPtr)pnt, aip, SEQFEATDATA_rsite))
                goto erret;
            break;
        case SEQFEAT_USER:       /* user object */
            if (! UserObjectAsnWrite((UserObjectPtr)pnt, aip, SEQFEATDATA_user))
                goto erret;
            break;
        case SEQFEAT_TXINIT:
            if (! TxinitAsnWrite((TxinitPtr)pnt, aip, SEQFEATDATA_txinit))
                goto erret;
            break;
        case SEQFEAT_NUM:
            if (! NumberingAsnWrite((NumberingPtr)pnt, aip, SEQFEATDATA_num))
                goto erret;
            break;
        case SEQFEAT_PSEC_STR:       /* psec-struct - int */
            av.intvalue = cp->value.intvalue;
            if (! AsnWrite(aip, SEQFEATDATA_psec_str, &av)) goto erret;
            break;
        case SEQFEAT_NON_STD_RESIDUE:        /* non-std residue */
            av.ptrvalue = pnt;
            if (! AsnWrite(aip, SEQFEATDATA_non_std_residue, &av)) goto erret;
            break;
        case SEQFEAT_HET:        /* heterogen */
            av.ptrvalue = pnt;
            if (! AsnWrite(aip, SEQFEATDATA_het, &av)) goto erret;
            break;
        case SEQFEAT_BIOSRC:
            if (! BioSourceAsnWrite((BioSourcePtr)pnt, aip, SEQFEATDATA_biosrc))
                goto erret;
            break;
        case SEQFEAT_CLONEREF:
            if (! CloneRefAsnWrite((CloneRefPtr)pnt, aip, SEQFEATDATA_clone))
                goto erret;
            break;
        case SEQFEAT_VARIATIONREF:
            if (! VariationRefAsnWrite((VariationRefPtr)pnt, aip, SEQFEATDATA_variation))
                goto erret;
            break;
    }

    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SeqFeatDataAsnRead(aip, atp, cp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatDataAsnRead (AsnIoPtr aip, AsnTypePtr orig, ChoicePtr cp)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return retval;
    }

    if (aip == NULL)
        return retval;

    if (orig == NULL)           /* SeqFeatData ::= (self contained) */
        atp = AsnReadId(aip, amp, SEQFEATDATA);
    else
        atp = AsnLinkType(orig, SEQFEATDATA);    /* link in local tree */
    if (atp == NULL)
        return retval;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read CHOICE */
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    if (atp == SEQFEATDATA_gene)
    {
        cp->choice = SEQFEAT_GENE ;
        cp->value.ptrvalue = (Pointer)GeneRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_org)
    {
        cp->choice = SEQFEAT_ORG ;
        cp->value.ptrvalue = (Pointer)OrgRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_cdregion)
    {
        cp->choice = SEQFEAT_CDREGION ;
        cp->value.ptrvalue = (Pointer)CdRegionAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_prot)
    {
        cp->choice = SEQFEAT_PROT ;
        cp->value.ptrvalue = (Pointer)ProtRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_rna)
    {
        cp->choice = SEQFEAT_RNA ;
        cp->value.ptrvalue = (Pointer)RnaRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_pub)
    {
        cp->choice = SEQFEAT_PUB ;
        cp->value.ptrvalue = (Pointer)PubdescAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_seq)
    {
        cp->choice = SEQFEAT_SEQ ;
        cp->value.ptrvalue = (Pointer)SeqLocAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_imp)
    {
        cp->choice = SEQFEAT_IMP ;
        cp->value.ptrvalue = (Pointer)ImpFeatAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_region)
    {
        cp->choice = SEQFEAT_REGION ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        cp->value.ptrvalue = av.ptrvalue;
    }
    else if (atp == SEQFEATDATA_comment)
    {
        cp->choice = SEQFEAT_COMMENT ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
    }
    else if (atp == SEQFEATDATA_bond)
    {
        cp->choice = SEQFEAT_BOND ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        cp->value.intvalue = av.intvalue;
    }
    else if (atp == SEQFEATDATA_site)
    {
        cp->choice = SEQFEAT_SITE ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        cp->value.intvalue = av.intvalue;
    }
    else if (atp == SEQFEATDATA_rsite)
    {
        cp->choice = SEQFEAT_RSITE ;
        cp->value.ptrvalue = (Pointer)RsiteRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_user)
    {
        cp->choice = SEQFEAT_USER ;
        cp->value.ptrvalue = (Pointer)UserObjectAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_txinit)
    {
        cp->choice = SEQFEAT_TXINIT ;
        cp->value.ptrvalue = (Pointer)TxinitAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_num)
    {
        cp->choice = SEQFEAT_NUM ;
        cp->value.ptrvalue = (Pointer)NumberingAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_psec_str)
    {
        cp->choice = SEQFEAT_PSEC_STR ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        cp->value.intvalue = av.intvalue;
    }
    else if (atp == SEQFEATDATA_non_std_residue)
    {
        cp->choice = SEQFEAT_NON_STD_RESIDUE ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        cp->value.ptrvalue = av.ptrvalue;
    }
    else if (atp == SEQFEATDATA_het)
    {
        cp->choice = SEQFEAT_HET ;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        cp->value.ptrvalue = av.ptrvalue;
    }
    else if (atp == SEQFEATDATA_biosrc)
    {
        cp->choice = SEQFEAT_BIOSRC ;
        cp->value.ptrvalue = (Pointer)BioSourceAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_clone)
    {
        cp->choice = SEQFEAT_CLONEREF ;
        cp->value.ptrvalue = (Pointer)CloneRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == SEQFEATDATA_variation)
    {
        cp->choice = SEQFEAT_VARIATIONREF ;
        cp->value.ptrvalue = (Pointer)VariationRefAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }

    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   GBQualNew()
*
*****************************************************************************/
NLM_EXTERN GBQualPtr LIBCALL GBQualNew (void)
{
    return (GBQualPtr)MemNew(sizeof(GBQual));
}

/*****************************************************************************
*
*   GBQualFree(gbp)
*       Frees a chain of GBQual and associated data
*
*****************************************************************************/
NLM_EXTERN GBQualPtr LIBCALL GBQualFree (GBQualPtr gbp)
{
    GBQualPtr next;

    while (gbp != NULL)
    {
        next = gbp->next;
        MemFree(gbp->qual);
        MemFree(gbp->val);
        MemFree(gbp);
        gbp = next;
    }
    return (GBQualPtr)NULL;
}

/*****************************************************************************
*
*   GBQualAsnWrite(gbp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (GBQual ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL GBQualAsnWrite (GBQualPtr gbp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, GB_QUAL);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (gbp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)gbp))
        goto erret;

    av.ptrvalue = (Pointer)gbp->qual;
    if (! AsnWrite(aip, GB_QUAL_qual, &av)) goto erret;
    av.ptrvalue = (Pointer)gbp->val;
    if (! AsnWrite(aip, GB_QUAL_val, &av)) goto erret;

    if (! AsnCloseStruct(aip, atp, (Pointer)gbp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   GBQualAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN GBQualPtr LIBCALL GBQualAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    GBQualPtr gbp=NULL;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return gbp;
    }

    if (aip == NULL)
        return gbp;

    if (orig == NULL)           /* GBQual ::= (self contained) */
        atp = AsnReadId(aip, amp, GB_QUAL);
    else
        atp = AsnLinkType(orig, GB_QUAL);    /* link in local tree */
    if (atp == NULL)
        return gbp;

    gbp = GBQualNew();
    if (gbp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
    gbp->qual = (CharPtr)av.ptrvalue;
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
    gbp->val = (CharPtr)av.ptrvalue;
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return gbp;
erret:
    gbp = GBQualFree(gbp);
    goto ret;
}

/*****************************************************************************
*
*   SeqFeatXrefNew()
*
*****************************************************************************/
NLM_EXTERN SeqFeatXrefPtr LIBCALL SeqFeatXrefNew (void)
{
    return (SeqFeatXrefPtr)MemNew(sizeof(SeqFeatXref));
}

/*****************************************************************************
*
*   SeqFeatToXref(sfp)
*       Copy a SeqFeat to an Xref
*       Currently only works on Gene-ref and Prot-ref
*
*****************************************************************************/
NLM_EXTERN SeqFeatXrefPtr LIBCALL SeqFeatToXref (SeqFeatPtr sfp)
{
    SeqFeatXrefPtr sfxp;

    if (sfp == NULL)
        return (SeqFeatXrefPtr)NULL;

    sfxp = SeqFeatXrefNew();
    if (sfxp == NULL) return NULL;
    SeqFeatIdDup(&sfxp->id, &sfp->id);
    switch (sfp->data.choice)
    {
        case SEQFEAT_GENE :
            sfxp->data.value.ptrvalue = (Pointer) GeneRefDup((GeneRefPtr)sfp->data.value.ptrvalue);
            sfxp->data.choice = sfp->data.choice;
            break;
        case SEQFEAT_PROT :
            sfxp->data.value.ptrvalue = (Pointer) ProtRefDup((ProtRefPtr)sfp->data.value.ptrvalue);
            sfxp->data.choice = sfp->data.choice;
            break;
        default:
            ErrPostEx(SEV_WARNING, 0, 0, "Unsupported SeqFeatToXref");
            break;
    }

    return sfxp;
}

/*****************************************************************************
*
*   SeqFeatXrefFree(sfxp)
*       Frees a chain of SeqFeatXref and associated data
*
*****************************************************************************/
NLM_EXTERN SeqFeatXrefPtr LIBCALL SeqFeatXrefFree (SeqFeatXrefPtr sfxp)
{
    SeqFeatXrefPtr next;

    while (sfxp != NULL)
    {
        next = sfxp->next;
        SeqFeatIdFree(&sfxp->id);
        SeqFeatDataFree(&sfxp->data);
        MemFree(sfxp);
        sfxp = next;
    }
    return (SeqFeatXrefPtr)NULL;
}

/*****************************************************************************
*
*   SeqFeatXrefAsnWrite(sfxp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqFeatXref ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatXrefAsnWrite (SeqFeatXrefPtr sfxp, AsnIoPtr aip, AsnTypePtr orig)
{
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, SEQFEATXREF);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (sfxp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)sfxp))
        goto erret;

    if (sfxp->id.choice)              /* feat-id present */
    {
        if (! SeqFeatIdAsnWrite(&sfxp->id, aip, SEQFEATXREF_id))
            goto erret;
    }

    if (sfxp->data.choice)            /* data present */
    {
        if (! SeqFeatDataAsnWrite(&sfxp->data, aip, SEQFEATXREF_data))
            goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)sfxp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SeqFeatXrefAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN SeqFeatXrefPtr LIBCALL SeqFeatXrefAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    SeqFeatXrefPtr sfxp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (SeqFeatXrefPtr)NULL;
    }

    if (aip == NULL)
        return (SeqFeatXrefPtr)NULL;

    if (orig == NULL)           /* SeqFeatXref ::= (self contained) */
        atp = AsnReadId(aip, amp, SEQFEATXREF);
    else
        atp = AsnLinkType(orig, SEQFEATXREF);    /* link in local tree */
    if (atp == NULL)
        return (SeqFeatXrefPtr)NULL;

    sfxp = SeqFeatXrefNew();
    if (sfxp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;

    if (atp == SEQFEATXREF_id)
    {
        if (! SeqFeatIdAsnRead(aip, atp, &sfxp->id))
            goto erret;
        atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    }

    if (atp == SEQFEATXREF_data)
    {
        if (! SeqFeatDataAsnRead(aip, atp, &sfxp->data))
            goto erret;
        atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    }

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return sfxp;
erret:
    sfxp = SeqFeatXrefFree(sfxp);
    goto ret;
}



/**************************************************
*
*    ModelEvidenceSupportNew()
*
**************************************************/
NLM_EXTERN 
ModelEvidenceSupportPtr LIBCALL
ModelEvidenceSupportNew(void)
{
   ModelEvidenceSupportPtr ptr = MemNew((size_t) sizeof(ModelEvidenceSupport));

   ptr -> full_length = 0;
   ptr -> supports_all_exon_combo = 0;
   return ptr;

}


/**************************************************
*
*    ModelEvidenceSupportFree()
*
**************************************************/
NLM_EXTERN 
ModelEvidenceSupportPtr LIBCALL
ModelEvidenceSupportFree(ModelEvidenceSupportPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> method);
   AsnGenericValNodeSetFree(ptr -> mrna, (AsnOptFreeFunc) ModelEvidenceItemFree);
   AsnGenericValNodeSetFree(ptr -> est, (AsnOptFreeFunc) ModelEvidenceItemFree);
   AsnGenericValNodeSetFree(ptr -> protein, (AsnOptFreeFunc) ModelEvidenceItemFree);
   SeqIdFree(ptr -> identification);
   AsnGenericValNodeSetFree(ptr -> dbxref, (AsnOptFreeFunc) DbtagFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ModelEvidenceSupportAsnRead()
*
**************************************************/
NLM_EXTERN 
ModelEvidenceSupportPtr LIBCALL
ModelEvidenceSupportAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ModelEvidenceSupportPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ModelEvidenceSupport ::= (self contained) */
      atp = AsnReadId(aip, amp, MODELEVIDENCESUPPORT);
   } else {
      atp = AsnLinkType(orig, MODELEVIDENCESUPPORT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ModelEvidenceSupportNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MODELEVIDENCESUPPORT_method) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> method = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_mrna) {
      ptr -> mrna = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ModelEvidenceItemAsnRead, (AsnOptFreeFunc) ModelEvidenceItemFree);
      if (isError && ptr -> mrna == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_est) {
      ptr -> est = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ModelEvidenceItemAsnRead, (AsnOptFreeFunc) ModelEvidenceItemFree);
      if (isError && ptr -> est == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_protein) {
      ptr -> protein = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ModelEvidenceItemAsnRead, (AsnOptFreeFunc) ModelEvidenceItemFree);
      if (isError && ptr -> protein == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_identification) {
      ptr -> identification = SeqIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_dbxref) {
      ptr -> dbxref = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      if (isError && ptr -> dbxref == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_exon_count) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> exon_count = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_exon_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> exon_length = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_full_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> full_length = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCESUPPORT_supports_all_exon_combo) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> supports_all_exon_combo = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ModelEvidenceSupportFree(ptr);
   goto ret;
}



/**************************************************
*
*    ModelEvidenceSupportAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ModelEvidenceSupportAsnWrite(ModelEvidenceSupportPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MODELEVIDENCESUPPORT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> method != NULL) {
      av.ptrvalue = ptr -> method;
      retval = AsnWrite(aip, MODELEVIDENCESUPPORT_method,  &av);
   }
   AsnGenericValNodeSetAsnWrite(ptr -> mrna, (AsnWriteFunc) ModelEvidenceItemAsnWrite, aip, MODELEVIDENCESUPPORT_mrna, MODELEVIDENCESUPPORT_mrna_E);
   AsnGenericValNodeSetAsnWrite(ptr -> est, (AsnWriteFunc) ModelEvidenceItemAsnWrite, aip, MODELEVIDENCESUPPORT_est, MODELEVIDENCESUPPORT_est_E);
   AsnGenericValNodeSetAsnWrite(ptr -> protein, (AsnWriteFunc) ModelEvidenceItemAsnWrite, aip, MODELEVIDENCESUPPORT_protein, MODELEVIDENCESUPPORT_protein_E);
   if (ptr -> identification != NULL) {
      if ( ! SeqIdAsnWrite(ptr -> identification, aip, MODELEVIDENCESUPPORT_identification)) {
         goto erret;
      }
   }
   AsnGenericValNodeSetAsnWrite(ptr -> dbxref, (AsnWriteFunc) DbtagAsnWrite, aip, MODELEVIDENCESUPPORT_dbxref, MODELEVIDENCESUPPORT_dbxref_E);
   if (ptr -> exon_count || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> exon_count;
      retval = AsnWrite(aip, MODELEVIDENCESUPPORT_exon_count,  &av);
   }
   if (ptr -> exon_length || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> exon_length;
      retval = AsnWrite(aip, MODELEVIDENCESUPPORT_exon_length,  &av);
   }
   if (ptr -> full_length) {
      av.boolvalue = ptr -> full_length;
      retval = AsnWrite(aip, MODELEVIDENCESUPPORT_full_length,  &av);
   }
   if (ptr -> supports_all_exon_combo) {
      av.boolvalue = ptr -> supports_all_exon_combo;
      retval = AsnWrite(aip, MODELEVIDENCESUPPORT_supports_all_exon_combo,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ExperimentSupportNew()
*
**************************************************/
NLM_EXTERN 
ExperimentSupportPtr LIBCALL
ExperimentSupportNew(void)
{
   ExperimentSupportPtr ptr = MemNew((size_t) sizeof(ExperimentSupport));

   return ptr;

}


/**************************************************
*
*    ExperimentSupportFree()
*
**************************************************/
NLM_EXTERN 
ExperimentSupportPtr LIBCALL
ExperimentSupportFree(ExperimentSupportPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> explanation);
   AsnGenericBaseSeqOfFree(ptr -> pmids ,ASNCODE_INTVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> dois ,ASNCODE_PTRVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    ExperimentSupportAsnRead()
*
**************************************************/
NLM_EXTERN 
ExperimentSupportPtr LIBCALL
ExperimentSupportAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   ExperimentSupportPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ExperimentSupport ::= (self contained) */
      atp = AsnReadId(aip, amp, EXPERIMENTSUPPORT);
   } else {
      atp = AsnLinkType(orig, EXPERIMENTSUPPORT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ExperimentSupportNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == EXPERIMENTSUPPORT_category) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> category = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EXPERIMENTSUPPORT_explanation) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> explanation = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EXPERIMENTSUPPORT_pmids) {
      ptr -> pmids = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> pmids == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EXPERIMENTSUPPORT_dois) {
      ptr -> dois = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> dois == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ExperimentSupportFree(ptr);
   goto ret;
}



/**************************************************
*
*    ExperimentSupportAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ExperimentSupportAsnWrite(ExperimentSupportPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, EXPERIMENTSUPPORT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> category || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> category;
      retval = AsnWrite(aip, EXPERIMENTSUPPORT_category,  &av);
   }
   if (ptr -> explanation != NULL) {
      av.ptrvalue = ptr -> explanation;
      retval = AsnWrite(aip, EXPERIMENTSUPPORT_explanation,  &av);
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> pmids ,ASNCODE_INTVAL_SLOT, aip, EXPERIMENTSUPPORT_pmids, EXPERIMENTSUPPORT_pmids_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> dois ,ASNCODE_PTRVAL_SLOT, aip, EXPERIMENTSUPPORT_dois, EXPERIMENTSUPPORT_dois_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ProgramIdNew()
*
**************************************************/
NLM_EXTERN 
ProgramIdPtr LIBCALL
ProgramIdNew(void)
{
   ProgramIdPtr ptr = MemNew((size_t) sizeof(ProgramId));

   return ptr;

}


/**************************************************
*
*    ProgramIdFree()
*
**************************************************/
NLM_EXTERN 
ProgramIdPtr LIBCALL
ProgramIdFree(ProgramIdPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   MemFree(ptr -> version);
   return MemFree(ptr);
}


/**************************************************
*
*    ProgramIdAsnRead()
*
**************************************************/
NLM_EXTERN 
ProgramIdPtr LIBCALL
ProgramIdAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   ProgramIdPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ProgramId ::= (self contained) */
      atp = AsnReadId(aip, amp, PROGRAM_ID);
   } else {
      atp = AsnLinkType(orig, PROGRAM_ID);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ProgramIdNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PROGRAM_ID_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PROGRAM_ID_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> version = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ProgramIdFree(ptr);
   goto ret;
}



/**************************************************
*
*    ProgramIdAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ProgramIdAsnWrite(ProgramIdPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PROGRAM_ID);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, PROGRAM_ID_name,  &av);
   }
   if (ptr -> version != NULL) {
      av.ptrvalue = ptr -> version;
      retval = AsnWrite(aip, PROGRAM_ID_version,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    EvidenceBasisNew()
*
**************************************************/
NLM_EXTERN 
EvidenceBasisPtr LIBCALL
EvidenceBasisNew(void)
{
   EvidenceBasisPtr ptr = MemNew((size_t) sizeof(EvidenceBasis));

   return ptr;

}


/**************************************************
*
*    EvidenceBasisFree()
*
**************************************************/
NLM_EXTERN 
EvidenceBasisPtr LIBCALL
EvidenceBasisFree(EvidenceBasisPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericValNodeSetFree(ptr -> programs, (AsnOptFreeFunc) ProgramIdFree);
   AsnGenericChoiceSeqOfFree(ptr -> accessions, (AsnOptFreeFunc) SeqIdFree);
   return MemFree(ptr);
}


/**************************************************
*
*    EvidenceBasisAsnRead()
*
**************************************************/
NLM_EXTERN 
EvidenceBasisPtr LIBCALL
EvidenceBasisAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   EvidenceBasisPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* EvidenceBasis ::= (self contained) */
      atp = AsnReadId(aip, amp, EVIDENCEBASIS);
   } else {
      atp = AsnLinkType(orig, EVIDENCEBASIS);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = EvidenceBasisNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == EVIDENCEBASIS_programs) {
      ptr -> programs = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ProgramIdAsnRead, (AsnOptFreeFunc) ProgramIdFree);
      if (isError && ptr -> programs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EVIDENCEBASIS_accessions) {
      ptr -> accessions = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) SeqIdAsnRead, (AsnOptFreeFunc) SeqIdFree);
      if (isError && ptr -> accessions == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = EvidenceBasisFree(ptr);
   goto ret;
}



/**************************************************
*
*    EvidenceBasisAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
EvidenceBasisAsnWrite(EvidenceBasisPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, EVIDENCEBASIS);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericValNodeSetAsnWrite(ptr -> programs, (AsnWriteFunc) ProgramIdAsnWrite, aip, EVIDENCEBASIS_programs, EVIDENCEBASIS_programs_E);
   AsnGenericChoiceSeqOfAsnWrite(ptr -> accessions, (AsnWriteFunc) SeqIdAsnWrite, aip, EVIDENCEBASIS_accessions, EVIDENCEBASIS_accessions_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    InferenceSupportNew()
*
**************************************************/
NLM_EXTERN 
InferenceSupportPtr LIBCALL
InferenceSupportNew(void)
{
   InferenceSupportPtr ptr = MemNew((size_t) sizeof(InferenceSupport));

   ptr -> type = 0;
   ptr -> same_species = 0;
   return ptr;

}


/**************************************************
*
*    InferenceSupportFree()
*
**************************************************/
NLM_EXTERN 
InferenceSupportPtr LIBCALL
InferenceSupportFree(InferenceSupportPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> other_type);
   EvidenceBasisFree(ptr -> basis);
   AsnGenericBaseSeqOfFree(ptr -> pmids ,ASNCODE_INTVAL_SLOT);
   AsnGenericBaseSeqOfFree(ptr -> dois ,ASNCODE_PTRVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    InferenceSupportAsnRead()
*
**************************************************/
NLM_EXTERN 
InferenceSupportPtr LIBCALL
InferenceSupportAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   InferenceSupportPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* InferenceSupport ::= (self contained) */
      atp = AsnReadId(aip, amp, INFERENCESUPPORT);
   } else {
      atp = AsnLinkType(orig, INFERENCESUPPORT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = InferenceSupportNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == INFERENCESUPPORT_category) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> category = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INFERENCESUPPORT_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INFERENCESUPPORT_other_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> other_type = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INFERENCESUPPORT_same_species) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> same_species = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INFERENCESUPPORT_basis) {
      ptr -> basis = EvidenceBasisAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INFERENCESUPPORT_pmids) {
      ptr -> pmids = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> pmids == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == INFERENCESUPPORT_dois) {
      ptr -> dois = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> dois == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = InferenceSupportFree(ptr);
   goto ret;
}



/**************************************************
*
*    InferenceSupportAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
InferenceSupportAsnWrite(InferenceSupportPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, INFERENCESUPPORT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> category || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> category;
      retval = AsnWrite(aip, INFERENCESUPPORT_category,  &av);
   }
   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, INFERENCESUPPORT_type,  &av);
   if (ptr -> other_type != NULL) {
      av.ptrvalue = ptr -> other_type;
      retval = AsnWrite(aip, INFERENCESUPPORT_other_type,  &av);
   }
   av.boolvalue = ptr -> same_species;
   retval = AsnWrite(aip, INFERENCESUPPORT_same_species,  &av);
   if (ptr -> basis != NULL) {
      if ( ! EvidenceBasisAsnWrite(ptr -> basis, aip, INFERENCESUPPORT_basis)) {
         goto erret;
      }
   }
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> pmids ,ASNCODE_INTVAL_SLOT, aip, INFERENCESUPPORT_pmids, INFERENCESUPPORT_pmids_E);
   retval = AsnGenericBaseSeqOfAsnWrite(ptr -> dois ,ASNCODE_PTRVAL_SLOT, aip, INFERENCESUPPORT_dois, INFERENCESUPPORT_dois_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ModelEvidenceItemNew()
*
**************************************************/
NLM_EXTERN 
ModelEvidenceItemPtr LIBCALL
ModelEvidenceItemNew(void)
{
   ModelEvidenceItemPtr ptr = MemNew((size_t) sizeof(ModelEvidenceItem));

   ptr -> full_length = 0;
   ptr -> supports_all_exon_combo = 0;
   return ptr;

}


/**************************************************
*
*    ModelEvidenceItemFree()
*
**************************************************/
NLM_EXTERN 
ModelEvidenceItemPtr LIBCALL
ModelEvidenceItemFree(ModelEvidenceItemPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SeqIdFree(ptr -> id);
   return MemFree(ptr);
}


/**************************************************
*
*    ModelEvidenceItemAsnRead()
*
**************************************************/
NLM_EXTERN 
ModelEvidenceItemPtr LIBCALL
ModelEvidenceItemAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   ModelEvidenceItemPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ModelEvidenceItem ::= (self contained) */
      atp = AsnReadId(aip, amp, MODELEVIDENCEITEM);
   } else {
      atp = AsnLinkType(orig, MODELEVIDENCEITEM);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ModelEvidenceItemNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == MODELEVIDENCEITEM_id) {
      ptr -> id = SeqIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCEITEM_exon_count) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> exon_count = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCEITEM_exon_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> exon_length = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCEITEM_full_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> full_length = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == MODELEVIDENCEITEM_supports_all_exon_combo) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> supports_all_exon_combo = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ModelEvidenceItemFree(ptr);
   goto ret;
}



/**************************************************
*
*    ModelEvidenceItemAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ModelEvidenceItemAsnWrite(ModelEvidenceItemPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, MODELEVIDENCEITEM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> id != NULL) {
      if ( ! SeqIdAsnWrite(ptr -> id, aip, MODELEVIDENCEITEM_id)) {
         goto erret;
      }
   }
   if (ptr -> exon_count || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> exon_count;
      retval = AsnWrite(aip, MODELEVIDENCEITEM_exon_count,  &av);
   }
   if (ptr -> exon_length || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> exon_length;
      retval = AsnWrite(aip, MODELEVIDENCEITEM_exon_length,  &av);
   }
   if (ptr -> full_length) {
      av.boolvalue = ptr -> full_length;
      retval = AsnWrite(aip, MODELEVIDENCEITEM_full_length,  &av);
   }
   if (ptr -> supports_all_exon_combo) {
      av.boolvalue = ptr -> supports_all_exon_combo;
      retval = AsnWrite(aip, MODELEVIDENCEITEM_supports_all_exon_combo,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SeqFeatSupportNew()
*
**************************************************/
NLM_EXTERN 
SeqFeatSupportPtr LIBCALL
SeqFeatSupportNew(void)
{
   SeqFeatSupportPtr ptr = MemNew((size_t) sizeof(SeqFeatSupport));

   return ptr;

}


/**************************************************
*
*    SeqFeatSupportFree()
*
**************************************************/
NLM_EXTERN 
SeqFeatSupportPtr LIBCALL
SeqFeatSupportFree(SeqFeatSupportPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericValNodeSetFree(ptr -> experiment, (AsnOptFreeFunc) ExperimentSupportFree);
   AsnGenericValNodeSetFree(ptr -> inference, (AsnOptFreeFunc) InferenceSupportFree);
   AsnGenericValNodeSetFree(ptr -> model_evidence, (AsnOptFreeFunc) ModelEvidenceSupportFree);
   return MemFree(ptr);
}


/**************************************************
*
*    SeqFeatSupportAsnRead()
*
**************************************************/
NLM_EXTERN 
SeqFeatSupportPtr LIBCALL
SeqFeatSupportAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SeqFeatSupportPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SeqFeatSupport ::= (self contained) */
      atp = AsnReadId(aip, amp, SEQFEATSUPPORT);
   } else {
      atp = AsnLinkType(orig, SEQFEATSUPPORT);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SeqFeatSupportNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == SEQFEATSUPPORT_experiment) {
      ptr -> experiment = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ExperimentSupportAsnRead, (AsnOptFreeFunc) ExperimentSupportFree);
      if (isError && ptr -> experiment == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQFEATSUPPORT_inference) {
      ptr -> inference = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) InferenceSupportAsnRead, (AsnOptFreeFunc) InferenceSupportFree);
      if (isError && ptr -> inference == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == SEQFEATSUPPORT_model_evidence) {
      ptr -> model_evidence = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ModelEvidenceSupportAsnRead, (AsnOptFreeFunc) ModelEvidenceSupportFree);
      if (isError && ptr -> model_evidence == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SeqFeatSupportFree(ptr);
   goto ret;
}



/**************************************************
*
*    SeqFeatSupportAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SeqFeatSupportAsnWrite(SeqFeatSupportPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, SEQFEATSUPPORT);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   AsnGenericValNodeSetAsnWrite(ptr -> experiment, (AsnWriteFunc) ExperimentSupportAsnWrite, aip, SEQFEATSUPPORT_experiment, SEQFEATSUPPORT_experiment_E);
   AsnGenericValNodeSetAsnWrite(ptr -> inference, (AsnWriteFunc) InferenceSupportAsnWrite, aip, SEQFEATSUPPORT_inference, SEQFEATSUPPORT_inference_E);
   AsnGenericValNodeSetAsnWrite(ptr -> model_evidence, (AsnWriteFunc) ModelEvidenceSupportAsnWrite, aip, SEQFEATSUPPORT_model_evidence, SEQFEATSUPPORT_model_evidence_E);
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}




/*****************************************************************************
*
*   CdRegionNew()
*
*****************************************************************************/
NLM_EXTERN CdRegionPtr LIBCALL CdRegionNew (void)
{
    return (CdRegionPtr)MemNew(sizeof(CdRegion));
}

/*****************************************************************************
*
*   CdRegionFree(cdp)
*       Frees a CdRegion and associated data
*
*****************************************************************************/
NLM_EXTERN CdRegionPtr LIBCALL CdRegionFree (CdRegionPtr cdp)
{
    if (cdp == NULL)
        return cdp;

    GeneticCodeFree(cdp->genetic_code);
    CodeBreakFree(cdp->code_break);
    return (CdRegionPtr)MemFree(cdp);
}

/*****************************************************************************
*
*   CdRegionAsnWrite(cdp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (CdRegion ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL CdRegionAsnWrite (CdRegionPtr cdp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    CodeBreakPtr cbp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, CDREGION);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (cdp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)cdp))
        goto erret;
    
    if (cdp->orf)
    {
        av.boolvalue = cdp->orf;
        if (! AsnWrite(aip, CDREGION_orf, &av)) goto erret;
    }
    if (cdp->frame)
    {
        av.intvalue = cdp->frame;
        if (! AsnWrite(aip, CDREGION_frame, &av)) goto erret;
    }
    if (cdp->conflict)
    {
        av.boolvalue = cdp->conflict;
        if (! AsnWrite(aip, CDREGION_conflict, &av)) goto erret;
    }
    if (cdp->gaps)
    {
        av.intvalue = cdp->gaps;
        if (! AsnWrite(aip, CDREGION_gaps, &av)) goto erret;
    }
    if (cdp->mismatch)
    {
        av.intvalue = cdp->mismatch;
        if (! AsnWrite(aip, CDREGION_mismatch, &av)) goto erret;
    }
    if (cdp->genetic_code != NULL)
    {
        if (! GeneticCodeAsnWrite(cdp->genetic_code, aip, CDREGION_code))
            goto erret;
    }
    if (cdp->code_break != NULL)
    {
        if (! AsnOpenStruct(aip, CDREGION_code_break, (Pointer)cdp->code_break))
            goto erret;
        cbp = cdp->code_break;
        while (cbp != NULL)
        {
            if (! CodeBreakAsnWrite(cbp, aip, CDREGION_code_break_E))
                goto erret;
            cbp = cbp->next;
        }
        if (! AsnCloseStruct(aip, CDREGION_code_break, (Pointer)cdp->code_break))
            goto erret;
    }
    if (cdp->stops)
    {
        av.intvalue = cdp->stops;
        if (! AsnWrite(aip, CDREGION_stops, &av)) goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)cdp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   CdRegionAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN CdRegionPtr LIBCALL CdRegionAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    CdRegionPtr cdp;
    CodeBreakPtr curr, next;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (CdRegionPtr)NULL;
    }

    if (aip == NULL)
        return (CdRegionPtr)NULL;

    if (orig == NULL)           /* CdRegion ::= (self contained) */
        atp = AsnReadId(aip, amp, CDREGION);
    else
        atp = AsnLinkType(orig, CDREGION);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (CdRegionPtr)NULL;

    cdp = CdRegionNew();
    if (cdp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (atp == CDREGION_code)
        {
            cdp->genetic_code = GeneticCodeAsnRead(aip, atp);
            if (cdp->genetic_code == NULL)
                goto erret;
        }
        else
        {
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            if (atp == CDREGION_orf)
                cdp->orf = av.boolvalue;
            else if (atp == CDREGION_conflict)
                cdp->conflict = av.boolvalue;
            else if (atp == CDREGION_gaps)
            {
                if (av.intvalue > 255)
                    av.intvalue = 255;
                cdp->gaps = (Uint1)av.intvalue;
            }
            else if (atp == CDREGION_frame)
            {
                if (av.intvalue > 255)
                    av.intvalue = 255;
                cdp->frame = (Uint1)av.intvalue;
            }
            else if (atp == CDREGION_mismatch)
            {
                if (av.intvalue > 255)
                    av.intvalue = 255;
                cdp->mismatch = (Uint1)av.intvalue;
            }
            else if (atp == CDREGION_stops)
            {
                if (av.intvalue > 255)
                    av.intvalue = 255;
                cdp->stops = (Uint1)av.intvalue;
            }
            else if (atp == CDREGION_code_break)
            {
                curr = NULL;
                while ((atp = AsnReadId(aip, amp, atp)) == CDREGION_code_break_E)
                {
                    next = CodeBreakAsnRead(aip, atp);
                    if (next == NULL)
                        goto erret;
                    if (cdp->code_break == NULL)
                        cdp->code_break = next;
                    else
                        curr->next = next;
                    curr = next;
                }
                if (atp == NULL)
                    goto erret;
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            }

        }
    }

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return cdp;
erret:
    cdp = CdRegionFree(cdp);
    goto ret;
}

static GeneticCodePtr genetic_codes = NULL;    /* array of ValNodes */
static Boolean GeneticCodeAdd PROTO((GeneticCodePtr gc));

/*****************************************************************************
*
*   GeneticCodeNew()
*       returns a ValNodePtr for the head of the chain
*       first ValNode has the choice 254, 
*       The chain of GeneticCode elements is attached to data.ptrvalue of
*         this first node.
*
*****************************************************************************/
NLM_EXTERN GeneticCodePtr LIBCALL GeneticCodeNew (void)
{
    GeneticCodePtr gc;

    gc = ValNodeNew(NULL);
    if (gc == NULL) return gc;
    gc->choice = 254;

    return gc;
}

/*****************************************************************************
*
*   GeneticCodeAdd()
*       Adds to genetic_codes global list
*       This should only be used by GeneticCodeTableRead()
*
*****************************************************************************/
static Boolean GeneticCodeAdd (GeneticCodePtr gc)
{
    GeneticCodePtr tmp;

    if (gc == NULL)
        return FALSE;

    if (genetic_codes == NULL)
        genetic_codes = gc;
    else
    {
        tmp = genetic_codes;
        while (tmp->next != NULL)
            tmp = tmp->next;
        tmp->next = gc;
    }

    return TRUE;
}

/*****************************************************************************
*
*   GeneticCodeFree(gcp)
*       Frees a GeneticCode and associated data
*
*****************************************************************************/
NLM_EXTERN GeneticCodePtr LIBCALL GeneticCodeFree (GeneticCodePtr gcp)
{
    ValNodePtr prev, next;

    if (gcp == NULL)
        return (GeneticCodePtr)NULL;

    if (gcp->choice == 254)   /* head of chain */
    {
        next = genetic_codes;
        prev = NULL;
        while ((next != NULL) && (gcp != next))
        {
            prev = next;
            next = prev->next;
        }

        if (next != NULL)      /* found it */
        {
            if (prev != NULL)   /* not the first one */
                prev->next = next->next;
            else
                genetic_codes = next->next;
            gcp->next = NULL;
        }

        prev = gcp;
        gcp = (GeneticCodePtr) prev->data.ptrvalue;
        MemFree(prev);
    }

    while (gcp != NULL)
    {
        next = gcp->next;
        switch (gcp->choice)
        {
            case 1:   /* name */
            case 3:   /* ncbieaa */
            case 6:   /* sncbieaa */
                MemFree(gcp->data.ptrvalue);
                break;
            case 4:   /* ncbi8aa */
            case 5:      /* ncbistdaa */
            case 7:   /* sncbi8aa */
            case 8:   /* sncbistdaa */
                BSFree((ByteStorePtr)gcp->data.ptrvalue);
                break;
            default:
                break;
        }
        MemFree(gcp);
        gcp = next;
    }
    return (GeneticCodePtr)NULL;
}

/*****************************************************************************
*
*   GeneticCodeAsnWrite(gcp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (GeneticCode ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL GeneticCodeAsnWrite (GeneticCodePtr gcp, AsnIoPtr aip, AsnTypePtr orig)
{
    GeneticCodePtr tmp;
    DataVal av;
    AsnTypePtr atp, gtp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, GENETIC_CODE);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (gcp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)gcp))
        goto erret;

    tmp = (ValNodePtr) gcp->data.ptrvalue;
    while (tmp != NULL)
    {   
        switch (tmp->choice)
        {
            case 1:   /* name */
                gtp = GENETIC_CODE_E_name;
                break;
            case 2:
                gtp = GENETIC_CODE_E_id;
                break;
            case 3:   /* ncbieaa */
                gtp = GENETIC_CODE_E_ncbieaa;
                break;
            case 4:      /* ncbi8aa */
                gtp = GENETIC_CODE_E_ncbi8aa;
                break;
            case 5:   /* ncbistdaa */
                gtp = GENETIC_CODE_E_ncbistdaa;
                break;
            case 6:   /* sncbieaa */
                gtp = GENETIC_CODE_E_sncbieaa;
                break;
            case 7:      /* sncbi8aa */
                gtp = GENETIC_CODE_E_sncbi8aa;
                break;
            case 8:   /* sncbistdaa */
                gtp = GENETIC_CODE_E_sncbistdaa;
                break;
            default:
                ErrPostEx(SEV_WARNING, 0, 0, "Unsupported GeneticCode Type = %d", (int)gcp->choice);
                goto erret;
        }
        
        av.ptrvalue = tmp;
        if (! AsnWriteChoice(aip, GENETIC_CODE_E, (Int2)tmp->choice, &av)) goto erret;
        if (! AsnWrite(aip, gtp, &tmp->data)) goto erret;
        tmp = tmp->next;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)gcp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}
/*****************************************************************************
*
*   GeneticCodePtr GeneticCodeAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN GeneticCodePtr LIBCALL GeneticCodeAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    GeneticCodePtr gcp = NULL, tmp = NULL;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return 0;
    }

    if (aip == NULL)
        return 0;

    if (orig == NULL)           /* GeneticCode ::= (self contained) */
        atp = AsnReadId(aip, amp, GENETIC_CODE);
    else
        atp = AsnLinkType(orig, GENETIC_CODE);    /* link in local tree */
    if (atp == NULL)
        return 0;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    gcp = GeneticCodeNew();
    if (gcp == NULL) goto erret;

    while((atp = AsnReadId(aip, amp, atp)) == GENETIC_CODE_E)
    {
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;    /* what kind? */

        tmp = ValNodeNew(tmp);
        if (tmp == NULL) goto erret;
        if (gcp->data.ptrvalue == NULL)
            gcp->data.ptrvalue = tmp;

        if (AsnReadVal(aip, atp, &tmp->data) <= 0) goto erret;

        if (atp == GENETIC_CODE_E_name)
            tmp->choice = 1;
        else if (atp == GENETIC_CODE_E_id)
            tmp->choice = 2;
        else if (atp == GENETIC_CODE_E_ncbieaa)
            tmp->choice = 3;
        else if (atp == GENETIC_CODE_E_ncbi8aa)
            tmp->choice = 4;
        else if (atp == GENETIC_CODE_E_ncbistdaa)
            tmp->choice = 5;
        else if (atp == GENETIC_CODE_E_sncbieaa)
            tmp->choice = 6;
        else if (atp == GENETIC_CODE_E_sncbi8aa)
            tmp->choice = 7;
        else if (atp == GENETIC_CODE_E_sncbistdaa)
            tmp->choice = 8;
        else
            tmp->choice = 255;   /* unrecognized */
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    goto ret;

erret:
    gcp = GeneticCodeFree(gcp);
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return gcp;
}

/*****************************************************************************
*
*   genCodeTblMemStr as last resort embedded version of gc.prt
*
*****************************************************************************/

#ifndef WIN16
static CharPtr genCodeTblMemStr = "Genetic-code-table ::= {\n" \
"{ name \"Standard\" , name \"SGC0\" , id 1 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"---M---------------M---------------M----------------------------\" } ,\n" \
"{ name \"Vertebrate Mitochondrial\" , name \"SGC1\" , id 2 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"--------------------------------MMMM---------------M------------\" } ,\n" \
"{ name \"Yeast Mitochondrial\" , name \"SGC2\" , id 3 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"----------------------------------MM----------------------------\" } ,\n" \
"{ name \"Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate\n" \
"Mitochondrial; Mycoplasma; Spiroplasma\" ,\n" \
"name \"SGC3\" , id 4 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"--MM---------------M------------MMMM---------------M------------\" } ,\n" \
"{ name \"Invertebrate Mitochondrial\" , name \"SGC4\" , id 5 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"---M----------------------------MMMM---------------M------------\" } ,\n" \
"{ name \"Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear\" ,\n" \
"name \"SGC5\" , id 6 ,\n" \
"ncbieaa  \"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M----------------------------\" } ,\n" \
"{ name \"Echinoderm Mitochondrial; Flatworm Mitochondrial\" , name \"SGC8\" , id 9 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M---------------M------------\" } ,\n" \
"{ name \"Euplotid Nuclear\" , name \"SGC9\" , id 10 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M----------------------------\" } ,\n" \
"{ name \"Bacterial, Archaeal and Plant Plastid\" , id 11 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"---M---------------M------------MMMM---------------M------------\" } ,\n" \
"{ name \"Alternative Yeast Nuclear\" , id 12 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-------------------M---------------M----------------------------\" } ,\n" \
"{ name \"Ascidian Mitochondrial\" , id 13 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"---M------------------------------MM---------------M------------\" } ,\n" \
"{ name \"Alternative Flatworm Mitochondrial\" , id 14 ,\n" \
"ncbieaa  \"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M----------------------------\" } ,\n" \
"{ name \"Blepharisma Macronuclear\" , id 15 ,\n" \
"ncbieaa  \"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M----------------------------\" } ,\n" \
"{ name \"Chlorophycean Mitochondrial\" , id 16 ,\n" \
"ncbieaa  \"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M----------------------------\" } ,\n" \
"{ name \"Trematode Mitochondrial\" , id 21 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M---------------M------------\" } ,\n" \
"{ name \"Scenedesmus obliquus Mitochondrial\" , id 22 ,\n" \
"ncbieaa  \"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"-----------------------------------M----------------------------\" } ,\n" \
"{ name \"Thraustochytrium Mitochondrial\" , id 23 ,\n" \
"ncbieaa  \"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"--------------------------------M--M---------------M------------\" } ,\n" \
"{ name \"Pterobranchia Mitochondrial\" , id 24 ,\n" \
"ncbieaa  \"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG\",\n" \
"sncbieaa \"---M---------------M---------------M---------------M------------\" } };\n";
#endif

/*****************************************************************************
*
*   GeneticCodePtr GeneticCodeTableLoad()
*       loads all current genetic codes
*       looks for "gc.val" in the "data" directory
*       now uses genCodeTblMemStr if unable to find "data" directory
*
*****************************************************************************/
static Boolean LoadGenCodeFromLocalString (void)

{
#ifndef WIN16
  AsnIoMemPtr aimp;

  aimp = AsnIoMemOpen ("r", (BytePtr) genCodeTblMemStr, (Int4) StringLen (genCodeTblMemStr));
  if (aimp == NULL || aimp->aip == NULL) return FALSE;
  GeneticCodeTableAsnRead (aimp->aip, NULL);
  AsnIoMemClose (aimp);
#endif
  return (Boolean) (genetic_codes != NULL);
}

NLM_EXTERN GeneticCodePtr LIBCALL GeneticCodeTableLoad (void)
{
    Char buf[PATH_MAX];
    AsnIoPtr aip;

    if (genetic_codes != NULL)
        return genetic_codes;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (GeneticCodePtr)NULL;
    }

#ifdef OS_UNIX
    if (getenv ("USE_GENCODE_FILE") == NULL) {
          if (LoadGenCodeFromLocalString ()) {
              return genetic_codes;
          }
    }
#endif

    if (! FindPath("ncbi", "ncbi", "data", buf, sizeof (buf)))
    {

        if (LoadGenCodeFromLocalString ()) {
            return genetic_codes;
        }

        ErrPostEx(SEV_WARNING, 0, 0, "FindPath failed in GeneticCodeTableLoad - ncbi configuration file missing or incorrect");
        return genetic_codes;
    }

    StringCat(buf, "gc.val");
    if ((aip = AsnIoOpen(buf, "rb")) == NULL)
    {

        if (LoadGenCodeFromLocalString ()) {
            return genetic_codes;
        }

        ErrPostEx(SEV_WARNING, 0, 0, "Couldn't open [%s]", buf);
        return genetic_codes;
    }

    GeneticCodeTableAsnRead(aip, NULL);

    AsnIoClose(aip);
    return genetic_codes;
}

/*****************************************************************************
*
*   GeneticCodeFind (id, name)
*       if (name == NULL)  finds by id
*       returns NULL, if can't find anything
*
*****************************************************************************/
NLM_EXTERN GeneticCodePtr LIBCALL GeneticCodeFind (Int4 id, CharPtr name)
{
    GeneticCodePtr gcp;
    ValNodePtr tmp;

    if (GeneticCodeTableLoad() == NULL)
        return (GeneticCodePtr)NULL;

    for (gcp = genetic_codes; gcp != NULL; gcp = gcp->next)
    {
        for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next)
        {
            switch (tmp->choice)
            {
                case 1:   /* name */
                    if (name != NULL)
                    {
                        if (! StringICmp(name, (CharPtr)tmp->data.ptrvalue))
                            return gcp;
                    }
                    break;
                case 2:   /* id */
                    if (name == NULL)
                    {
                        if (id == tmp->data.intvalue)
                            return gcp;
                    }
                    break;
                default:
                    break;
            }
        }
    }

    return (GeneticCodePtr)NULL;
}


/*****************************************************************************
*
*   GeneticCodeTableAsnWrite(gcp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (Genetic-code-table ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL GeneticCodeTableAsnWrite (GeneticCodePtr gcp, AsnIoPtr aip, AsnTypePtr orig)
{
    GeneticCodePtr tmp;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, GENETIC_CODE_TABLE);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (gcp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)gcp))
        goto erret;
    
    tmp = gcp;
    while (tmp != NULL)
    {
        if (! GeneticCodeAsnWrite(tmp, aip, GENETIC_CODE_TABLE_E))
            goto erret;
        tmp = tmp->next;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)gcp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   GeneticCodePtr GeneticCodeTableAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN GeneticCodePtr LIBCALL GeneticCodeTableAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    GeneticCodePtr gcp = NULL, tmp = NULL;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return 0;
    }

    if (aip == NULL)
        return 0;

    if (orig == NULL)           /* Genetic-code-table ::= (self contained) */
        atp = AsnReadId(aip, amp, GENETIC_CODE_TABLE);
    else
        atp = AsnLinkType(orig, GENETIC_CODE_TABLE);    /* link in local tree */
    if (atp == NULL)
        return 0;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while((atp = AsnReadId(aip, amp, atp)) == GENETIC_CODE_TABLE_E)
    {
        if ((tmp = GeneticCodeAsnRead(aip, atp)) == NULL) goto erret;
        if (gcp == NULL)     /* keep first one */
            gcp = tmp;
        GeneticCodeAdd(tmp);   /* add to global genetic_codes */
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    goto ret;

erret:
    while (gcp != NULL)
    {
        tmp = gcp->next;
        gcp = GeneticCodeFree(gcp);
        gcp = tmp;
    }
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return gcp;
}

/*****************************************************************************
*
*   CodeBreakNew()
*
*****************************************************************************/
NLM_EXTERN CodeBreakPtr LIBCALL CodeBreakNew (void)
{
    return (CodeBreakPtr)MemNew(sizeof(CodeBreak));
}

/*****************************************************************************
*
*   CodeBreakFree(cbp)
*       Frees a chain of CodeBreak and associated data
*
*****************************************************************************/
NLM_EXTERN CodeBreakPtr LIBCALL CodeBreakFree (CodeBreakPtr cbp)
{
    CodeBreakPtr next;

    while (cbp != NULL)
    {
        next = cbp->next;
        SeqLocFree(cbp->loc);
        MemFree(cbp);
        cbp = next;
    }
    return (CodeBreakPtr)NULL;
}

/*****************************************************************************
*
*   CodeBreakAsnWrite(cbp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (CodeBreak ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL CodeBreakAsnWrite (CodeBreakPtr cbp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, tmp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, CODE_BREAK);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (cbp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)cbp))
        goto erret;

    if (! SeqLocAsnWrite(cbp->loc, aip, CODE_BREAK_loc))
        goto erret;
    av.intvalue = cbp->aa.value.intvalue;
    if (! AsnWriteChoice(aip, CODE_BREAK_aa, (Int2)cbp->aa.choice, &av)) goto erret;
    switch(cbp->aa.choice)
    {
        case 1:
            tmp = CODE_BREAK_aa_ncbieaa;
            break;
        case 2:
            tmp = CODE_BREAK_aa_ncbi8aa;
            break;
        case 3:
            tmp = CODE_BREAK_aa_ncbistdaa;
            break;
    }
    if (! AsnWrite(aip, tmp, &av)) goto erret;
    if (! AsnCloseStruct(aip, atp, (Pointer)cbp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   CodeBreakAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN CodeBreakPtr LIBCALL CodeBreakAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    CodeBreakPtr cbp;
    Uint1 choice = 0;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (CodeBreakPtr)NULL;
    }

    if (aip == NULL)
        return (CodeBreakPtr)NULL;

    if (orig == NULL)           /* CodeBreak ::= (self contained) */
        atp = AsnReadId(aip, amp, CODE_BREAK);
    else
        atp = AsnLinkType(orig, CODE_BREAK);    /* link in local tree */
    if (atp == NULL)
        return (CodeBreakPtr)NULL;

    cbp = CodeBreakNew();
    if (cbp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;   /* location */
    cbp->loc = SeqLocAsnRead(aip, atp);
    if (cbp->loc == NULL)
        goto erret;
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;    /* read the CHOICE */
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;    /* the aa */
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
    cbp->aa.value.intvalue = av.intvalue;
    if (atp == CODE_BREAK_aa_ncbieaa)
        choice = 1;
    else if (atp == CODE_BREAK_aa_ncbi8aa)
        choice = 2;
    else if (atp == CODE_BREAK_aa_ncbistdaa)
        choice = 3;
    cbp->aa.choice = choice;

    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return cbp;
erret:
    cbp = CodeBreakFree(cbp);
    goto ret;
}

/*****************************************************************************
*
*   ImpFeatNew()
*
*****************************************************************************/
NLM_EXTERN ImpFeatPtr LIBCALL ImpFeatNew (void)
{
    return (ImpFeatPtr)MemNew(sizeof(ImpFeat));
}

/*****************************************************************************
*
*   ImpFeatFree(ifp)
*       Frees an ImpFeat and associated data
*
*****************************************************************************/
NLM_EXTERN ImpFeatPtr LIBCALL ImpFeatFree (ImpFeatPtr ifp)
{
    if (ifp == NULL)
        return (ImpFeatPtr)NULL;

    MemFree(ifp->key);
    MemFree(ifp->loc);
    MemFree(ifp->descr);
    return (ImpFeatPtr)MemFree(ifp);
}

/*****************************************************************************
*
*   ImpFeatAsnWrite(ifp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (ImpFeat ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL ImpFeatAsnWrite (ImpFeatPtr ifp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, IMP_FEAT);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (ifp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)ifp))
        goto erret;

    av.ptrvalue = ifp->key;
    if (! AsnWrite(aip, IMP_FEAT_key, &av)) goto erret;
    if (ifp->loc != NULL)
    {
        av.ptrvalue = ifp->loc;
        if (! AsnWrite(aip, IMP_FEAT_loc, &av)) goto erret;
    }
    if (ifp->descr != NULL)
    {
        av.ptrvalue = ifp->descr;
        if (! AsnWrite(aip, IMP_FEAT_descr, &av)) goto erret;
    }
    if (! AsnCloseStruct(aip, atp, (Pointer)ifp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   ImpFeatAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN ImpFeatPtr LIBCALL ImpFeatAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    ImpFeatPtr ifp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (ImpFeatPtr)NULL;
    }

    if (aip == NULL)
        return (ImpFeatPtr)NULL;

    if (orig == NULL)           /* ImpFeat ::= (self contained) */
        atp = AsnReadId(aip, amp, IMP_FEAT);
    else
        atp = AsnLinkType(orig, IMP_FEAT);    /* link in local tree */
    oldatp = atp;

    if (atp == NULL)
        return (ImpFeatPtr)NULL;

    ifp = ImpFeatNew();
    if (ifp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == IMP_FEAT_key)
            ifp->key = (CharPtr)av.ptrvalue;
        else if (atp == IMP_FEAT_loc)
            ifp->loc = (CharPtr)av.ptrvalue;
        else if (atp == IMP_FEAT_descr)
            ifp->descr = (CharPtr)av.ptrvalue;
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return ifp;
erret:
    ifp = ImpFeatFree(ifp);
    goto ret;
}

/*****************************************************************************
*
*   RnaRefNew()
*
*****************************************************************************/
NLM_EXTERN RnaRefPtr LIBCALL RnaRefNew (void)
{
    return (RnaRefPtr)MemNew(sizeof(RnaRef));
}

/*****************************************************************************
*
*   RnaRefFree(rrp)
*       Frees an RnaRef and associated data
*
*****************************************************************************/
NLM_EXTERN RnaRefPtr LIBCALL RnaRefFree (RnaRefPtr rrp)
{
    tRNAPtr trna;
    RNAGenPtr rgp;

    if (rrp == NULL)
        return (RnaRefPtr)NULL;
    
    if (rrp->ext.choice == 1)  /* name */
        MemFree(rrp->ext.value.ptrvalue);
    else if (rrp->ext.choice == 2)   /* tRNA */
    {
        trna = (tRNAPtr)(rrp->ext.value.ptrvalue);
        SeqLocFree(trna->anticodon);
        MemFree(trna);
    }
    else if (rrp->ext.choice == 3)   /* gen */
    {
        rgp = (RNAGenPtr)(rrp->ext.value.ptrvalue);
        RNAGenFree(rgp);
    }
    return (RnaRefPtr)MemFree(rrp);
}

/*****************************************************************************
*
*   RnaRefAsnWrite(rrp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (RnaRef ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL RnaRefAsnWrite (RnaRefPtr rrp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, tmp = NULL;
    tRNAPtr trna;
    RNAGenPtr rgp;
    Int2 i;
    Boolean first = TRUE, retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, RNA_REF);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (rrp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)rrp))
        goto erret;

    av.intvalue = rrp->type;
    if (! AsnWrite(aip, RNA_REF_type, &av)) goto erret;
    if (rrp->pseudo)
    {
        av.boolvalue = rrp->pseudo;
        if (! AsnWrite(aip, RNA_REF_pseudo, &av)) goto erret;
    }
    if (rrp->ext.choice)
    {
        av.ptrvalue = rrp->ext.value.ptrvalue;
        if (! AsnWriteChoice(aip, RNA_REF_ext, (Int2)rrp->ext.choice, &av)) goto erret;
        if (rrp->ext.choice == 1)   /* name */
        {
            if (! AsnWrite(aip, RNA_REF_ext_name, &av)) goto erret;
        }
        else if (rrp->ext.choice == 2)   /* trna */
        {
            trna = (tRNAPtr)rrp->ext.value.ptrvalue;
            if (! AsnOpenStruct(aip, RNA_REF_ext_tRNA, (Pointer)trna))
                goto erret;
            if (trna->aatype)
            {
                av.intvalue = trna->aa;
                if (! AsnWriteChoice(aip, TRNA_EXT_aa, (Int2)trna->aatype, &av)) goto erret;
                switch (trna->aatype)
                {
                    case 1:
                        tmp = TRNA_EXT_aa_iupacaa;
                        break;
                    case 2:
                        tmp = TRNA_EXT_aa_ncbieaa;
                        break;
                    case 3:
                        tmp = TRNA_EXT_aa_ncbi8aa;
                        break;
                    case 4:
                        tmp = TRNA_EXT_aa_ncbistdaa;
                        break;
                }
                if (! AsnWrite(aip, tmp, &av)) goto erret;
            }
            for (i = 0; i < 6; i++)
            {
                if (trna->codon[i] != 255)
                {
                    if (first)
                    {
                        if (! AsnOpenStruct(aip, TRNA_EXT_codon, (Pointer)trna->codon))
                            goto erret;
                        first = FALSE;
                    }
                    av.intvalue = trna->codon[i];
                    if (! AsnWrite(aip, TRNA_EXT_codon_E, &av)) goto erret;
                }
            }
            if (! first)    /* wrote one */
                if (! AsnCloseStruct(aip, TRNA_EXT_codon, (Pointer)trna->codon))
                    goto erret;
            if (trna->anticodon != NULL)
                if (! SeqLocAsnWrite(trna->anticodon, aip, TRNA_EXT_anticodon))
                    goto erret;

            if (! AsnCloseStruct(aip, RNA_REF_ext_tRNA, (Pointer)trna))
                goto erret;
        }
        else if (rrp->ext.choice == 3)   /* gen */
        {
            rgp = (RNAGenPtr)(rrp->ext.value.ptrvalue);
            if (! RNAGenAsnWrite (rgp, aip, RNA_REF_ext_gen))
                goto erret;
        }
    }
    if (! AsnCloseStruct(aip, atp, (Pointer)rrp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   RnaRefAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN RnaRefPtr LIBCALL RnaRefAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    RnaRefPtr rrp;
    tRNAPtr trna = NULL;
    Int2 i = 0;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (RnaRefPtr)NULL;
    }

    if (aip == NULL)
        return (RnaRefPtr)NULL;

    if (orig == NULL)           /* RnaRef ::= (self contained) */
        atp = AsnReadId(aip, amp, RNA_REF);
    else
        atp = AsnLinkType(orig, RNA_REF);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (RnaRefPtr)NULL;

    rrp = RnaRefNew();
    if (rrp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (atp == TRNA_EXT_anticodon)
        {
            trna->anticodon = SeqLocAsnRead(aip, atp);
            if (trna->anticodon == NULL) goto erret;
        }
        else
        {
            if (atp == RNA_REF_type)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                rrp->type = (Uint1)av.intvalue;
            }
            else if (atp == RNA_REF_pseudo)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                rrp->pseudo = av.boolvalue;
            }
            else if (atp == RNA_REF_ext_name)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                rrp->ext.choice = 1;
                rrp->ext.value.ptrvalue = av.ptrvalue;
            }
            else if (atp == RNA_REF_ext_tRNA)   /* its a trna */
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                if (trna == NULL)   /* first time */
                {
                    rrp->ext.choice = 2;
                    trna = (tRNAPtr)MemNew(sizeof(tRNA));
                    if (trna == NULL)
                        goto erret;
                    rrp->ext.value.ptrvalue = trna;
                    for (i = 0; i < 6; i++)
                        trna->codon[i] = 255;   /* not used */
                    i = 0;
                }
            }
            else if (atp == TRNA_EXT_aa_iupacaa)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                trna->aatype = 1;
                trna->aa = (Uint1)av.intvalue;
            }
            else if (atp == TRNA_EXT_aa_ncbieaa)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                trna->aatype = 2;
                trna->aa = (Uint1)av.intvalue;
            }
            else if (atp == TRNA_EXT_aa_ncbi8aa)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                trna->aatype = 3;
                trna->aa = (Uint1)av.intvalue;
            }
            else if (atp == TRNA_EXT_aa_ncbistdaa)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                trna->aatype = 4;
                trna->aa = (Uint1)av.intvalue;
            }
            else if (atp == TRNA_EXT_codon_E)
            {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                if (i < 6) {
                    trna->codon[i] = (Uint1)av.intvalue;
                }
                i++;
            }
            else if (atp == RNA_REF_ext_gen)
            {
                rrp->ext.choice = 3;
                rrp->ext.value.ptrvalue = RNAGenAsnRead (aip, atp);
            } else {
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            }
        }
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return rrp;
erret:
    rrp = RnaRefFree(rrp);
    goto ret;
}



/**************************************************
*
*    RNAQualNew()
*
**************************************************/
NLM_EXTERN 
RNAQualPtr LIBCALL
RNAQualNew(void)
{
   RNAQualPtr rqp = MemNew((size_t) sizeof(RNAQual));

   return rqp;

}


/**************************************************
*
*    RNAQualFree()
*
**************************************************/
NLM_EXTERN 
RNAQualPtr LIBCALL
RNAQualFree(RNAQualPtr rqp)
{

   if(rqp == NULL) {
      return NULL;
   }
   MemFree(rqp -> qual);
   MemFree(rqp -> val);
   return MemFree(rqp);
}


/**************************************************
*
*    RNAQualAsnRead()
*
**************************************************/
NLM_EXTERN 
RNAQualPtr LIBCALL
RNAQualAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   RNAQualPtr rqp;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RNAQual ::= (self contained) */
      atp = AsnReadId(aip, amp, RNA_QUAL);
   } else {
      atp = AsnLinkType(orig, RNA_QUAL);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   rqp = RNAQualNew();
   if (rqp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == RNA_QUAL_qual) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      rqp -> qual = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == RNA_QUAL_val) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      rqp -> val = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return rqp;

erret:
   aip -> io_failure = TRUE;
   rqp = RNAQualFree(rqp);
   goto ret;
}



/**************************************************
*
*    RNAQualAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RNAQualAsnWrite(RNAQualPtr rqp, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, RNA_QUAL);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (rqp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) rqp)) {
      goto erret;
   }

   if (rqp -> qual != NULL) {
      av.ptrvalue = rqp -> qual;
      retval = AsnWrite(aip, RNA_QUAL_qual,  &av);
   }
   if (rqp -> val != NULL) {
      av.ptrvalue = rqp -> val;
      retval = AsnWrite(aip, RNA_QUAL_val,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)rqp)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    RNAQualSetFree()
*
**************************************************/
NLM_EXTERN 
RNAQualSetPtr LIBCALL
RNAQualSetFree(RNAQualSetPtr rqp)
{

   if(rqp == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(rqp,  (AsnOptFreeFunc) RNAQualFree);
   return NULL;
}


/**************************************************
*
*    RNAQualSetAsnRead()
*
**************************************************/
NLM_EXTERN 
RNAQualSetPtr LIBCALL
RNAQualSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   RNAQualSetPtr rqp;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RNAQualSet ::= (self contained) */
      atp = AsnReadId(aip, amp, RNA_QUAL_SET);
   } else {
      atp = AsnLinkType(orig, RNA_QUAL_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   rqp  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) RNAQualAsnRead, (AsnOptFreeFunc) RNAQualFree);
   if (isError && rqp  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return rqp;

erret:
   aip -> io_failure = TRUE;
   rqp = RNAQualSetFree(rqp);
   goto ret;
}



/**************************************************
*
*    RNAQualSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RNAQualSetAsnWrite(RNAQualSetPtr rqp, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, RNA_QUAL_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (rqp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(rqp , (AsnWriteFunc) RNAQualAsnWrite, aip, atp, RNA_QUAL_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    RNAGenNew()
*
**************************************************/
NLM_EXTERN 
RNAGenPtr LIBCALL
RNAGenNew(void)
{
   RNAGenPtr rgp = MemNew((size_t) sizeof(RNAGen));
   return rgp;

}


/**************************************************
*
*    RNAGenFree()
*
**************************************************/
NLM_EXTERN 
RNAGenPtr LIBCALL
RNAGenFree(RNAGenPtr rgp)
{

   if(rgp == NULL) {
      return NULL;
   }
   MemFree(rgp -> _class);
   MemFree(rgp -> product);
   RNAQualSetFree(rgp -> quals);
   return MemFree(rgp);
}


/**************************************************
*
*    RNAGenAsnRead()
*
**************************************************/
NLM_EXTERN 
RNAGenPtr LIBCALL
RNAGenAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   RNAGenPtr rgp;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* RNAGen ::= (self contained) */
      atp = AsnReadId(aip, amp, RNA_GEN);
   } else {
      atp = AsnLinkType(orig, RNA_GEN);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   rgp = RNAGenNew();
   if (rgp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == RNA_GEN_class) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      rgp -> _class = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == RNA_GEN_product) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      rgp -> product = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == RNA_GEN_quals) {
      rgp -> quals = RNAQualSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return rgp;

erret:
   aip -> io_failure = TRUE;
   rgp = RNAGenFree(rgp);
   goto ret;
}



/**************************************************
*
*    RNAGenAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
RNAGenAsnWrite(RNAGenPtr rgp, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, RNA_GEN);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (rgp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) rgp)) {
      goto erret;
   }

   if (rgp -> _class != NULL) {
      av.ptrvalue = rgp -> _class;
      retval = AsnWrite(aip, RNA_GEN_class,  &av);
   }
   if (rgp -> product != NULL) {
      av.ptrvalue = rgp -> product;
      retval = AsnWrite(aip, RNA_GEN_product,  &av);
   }
   if (rgp -> quals != NULL) {
      if ( ! RNAQualSetAsnWrite(rgp -> quals, aip, RNA_GEN_quals)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)rgp)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}




/**************************************************
*
*    GeneNomenclatureNew()
*
**************************************************/
NLM_EXTERN 
GeneNomenclaturePtr LIBCALL
GeneNomenclatureNew (void)
{
   GeneNomenclaturePtr gnp = MemNew((size_t) sizeof(GeneNomenclature));

   return gnp;

}


/**************************************************
*
*    GeneNomenclatureFree()
*
**************************************************/
NLM_EXTERN 
GeneNomenclaturePtr LIBCALL
GeneNomenclatureFree (GeneNomenclaturePtr gnp)
{

   if (gnp == NULL) {
      return NULL;
   }
   MemFree(gnp -> symbol);
   MemFree(gnp -> name);
   DbtagFree(gnp -> source);
   return MemFree(gnp);
}


/**************************************************
*
*    GeneNomenclatureAsnRead()
*
**************************************************/
NLM_EXTERN 
GeneNomenclaturePtr LIBCALL
GeneNomenclatureAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   GeneNomenclaturePtr gnp;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* GeneNomenclature ::= (self contained) */
      atp = AsnReadId(aip, amp, GENE_NOMENCLATURE);
   } else {
      atp = AsnLinkType(orig, GENE_NOMENCLATURE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   gnp = GeneNomenclatureNew();
   if (gnp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == GENE_NOMENCLATURE_status) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      gnp -> status = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GENE_NOMENCLATURE_symbol) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      gnp -> symbol = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GENE_NOMENCLATURE_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      gnp -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == GENE_NOMENCLATURE_source) {
      gnp -> source = DbtagAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return gnp;

erret:
   aip -> io_failure = TRUE;
   gnp = GeneNomenclatureFree(gnp);
   goto ret;
}


/**************************************************
*
*    GeneNomenclatureAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
GeneNomenclatureAsnWrite (GeneNomenclaturePtr gnp, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, GENE_NOMENCLATURE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (gnp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) gnp)) {
      goto erret;
   }

   av.intvalue = gnp -> status;
   retval = AsnWrite(aip, GENE_NOMENCLATURE_status,  &av);
   if (gnp -> symbol != NULL) {
      av.ptrvalue = gnp -> symbol;
      retval = AsnWrite(aip, GENE_NOMENCLATURE_symbol,  &av);
   }
   if (gnp -> name != NULL) {
      av.ptrvalue = gnp -> name;
      retval = AsnWrite(aip, GENE_NOMENCLATURE_name,  &av);
   }
   if (gnp -> source != NULL) {
      if ( ! DbtagAsnWrite(gnp -> source, aip, GENE_NOMENCLATURE_source)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)gnp)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/*****************************************************************************
*
*   GeneRefNew()
*
*****************************************************************************/
NLM_EXTERN GeneRefPtr LIBCALL GeneRefNew (void)
{
    return (GeneRefPtr)MemNew(sizeof(GeneRef));
}

/*****************************************************************************
*
*   GeneRefFree(grp)
*       Frees an GeneRef and associated data
*
*****************************************************************************/
NLM_EXTERN GeneRefPtr LIBCALL GeneRefFree (GeneRefPtr grp)
{
    ValNodePtr anp, next;

    if (grp == NULL)
        return (GeneRefPtr)NULL;

    MemFree(grp->locus);
    MemFree(grp->allele);
    MemFree(grp->desc);
    MemFree(grp->maploc);
    anp = grp->db;
    while (anp != NULL)
    {
        next = anp->next;
        DbtagFree((DbtagPtr)anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    anp = grp->syn;
    while (anp != NULL)
    {
        next = anp->next;
        MemFree(anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    MemFree (grp->locus_tag);
    GeneNomenclatureFree (grp->formal_name);
    return (GeneRefPtr)MemFree(grp);
}

/*****************************************************************************
*
*   GeneRefDup(grp)
*       Duplicates an GeneRef and associated data
*
*****************************************************************************/
NLM_EXTERN GeneRefPtr LIBCALL GeneRefDup (GeneRefPtr grp)
{
    ValNodePtr anp, vnp;
    GeneRefPtr ngp;

    if (grp == NULL)
        return (GeneRefPtr)NULL;

    ngp = GeneRefNew();
    if (ngp == NULL) return ngp;

    if (grp->locus != NULL)
        ngp->locus = StringSave(grp->locus);
    if (grp->allele != NULL)
        ngp->allele = StringSave(grp->allele);
    if (grp->desc != NULL)
        ngp->desc = StringSave(grp->desc);
    if (grp->maploc != NULL)
        ngp->maploc = StringSave(grp->maploc);
    anp = grp->db;
    while (anp != NULL)
    {
        vnp = ValNodeNew(ngp->db);
        if (vnp == NULL) return GeneRefFree(ngp);
        if (ngp->db == NULL)
            ngp->db = vnp;
        vnp->data.ptrvalue = (Pointer)DbtagDup((DbtagPtr)anp->data.ptrvalue);
        anp = anp->next;
    }
    anp = grp->syn;
    while (anp != NULL)
    {
        vnp = ValNodeNew(ngp->syn);
        if (vnp == NULL) return GeneRefFree(ngp);
        if (ngp->syn == NULL)
            ngp->syn = vnp;
        vnp->data.ptrvalue = (Pointer)StringSave((CharPtr)anp->data.ptrvalue);
        anp = anp->next;
    }
    if (grp->locus_tag != NULL)
        ngp->locus_tag = StringSave(grp->locus_tag);
    if (grp->formal_name != NULL) {
      ngp = AsnIoMemCopy ((Pointer) grp->formal_name,
                          (AsnReadFunc) GeneNomenclatureAsnRead,
                          (AsnWriteFunc) GeneNomenclatureAsnWrite);
    }
    return ngp;
}

/*****************************************************************************
*
*   GeneRefAsnWrite(grp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (GeneRef ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL GeneRefAsnWrite (GeneRefPtr grp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    ValNodePtr anp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, GENE_REF);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (grp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)grp))
        goto erret;

    if (grp->locus != NULL)
    {
        av.ptrvalue = grp->locus;
        if (! AsnWrite(aip, GENE_REF_locus, &av)) goto erret;
    }
    if (grp->allele != NULL)
    {
        av.ptrvalue = grp->allele;
        if (! AsnWrite(aip, GENE_REF_allele, &av)) goto erret;
    }
    if (grp->desc != NULL)
    {
        av.ptrvalue = grp->desc;
        if (! AsnWrite(aip, GENE_REF_desc, &av)) goto erret;
    }
    if (grp->maploc != NULL)
    {
        av.ptrvalue = grp->maploc;
        if (! AsnWrite(aip, GENE_REF_maploc, &av)) goto erret;
    }
    if (grp->pseudo)
    {
        av.boolvalue = grp->pseudo;
        if (! AsnWrite(aip, GENE_REF_pseudo, &av)) goto erret;
    }
    if (grp->db != NULL)
    {
        anp = grp->db;
        if (! AsnOpenStruct(aip, GENE_REF_db, (Pointer)grp->db))
            goto erret;
        while (anp != NULL)
        {
            if (! DbtagAsnWrite((DbtagPtr)anp->data.ptrvalue, aip, GENE_REF_db_E))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, GENE_REF_db, (Pointer)grp->db))
            goto erret;
    }
    if (grp->syn != NULL)
    {
        anp = grp->syn;
        if (! AsnOpenStruct(aip, GENE_REF_syn, (Pointer)grp->syn))
            goto erret;
        while (anp != NULL)
        {
            if (! AsnWrite(aip, GENE_REF_syn_E, &anp->data)) goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, GENE_REF_syn, (Pointer)grp->syn))
            goto erret;
    }
    if (grp->locus_tag != NULL)
    {
        av.ptrvalue = grp->locus_tag;
        if (! AsnWrite(aip, GENE_REF_locus_tag, &av)) goto erret;
    }
   if (grp->formal_name != NULL) {
      if ( ! GeneNomenclatureAsnWrite(grp->formal_name, aip, GENE_REF_formal_name)) {
         goto erret;
      }
   }
    if (! AsnCloseStruct(aip, atp, (Pointer)grp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   GeneRefAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN GeneRefPtr LIBCALL GeneRefAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    GeneRefPtr grp;
    ValNodePtr anp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (GeneRefPtr)NULL;
    }

    if (aip == NULL)
        return (GeneRefPtr)NULL;

    if (orig == NULL)           /* GeneRef ::= (self contained) */
        atp = AsnReadId(aip, amp, GENE_REF);
    else
        atp = AsnLinkType(orig, GENE_REF);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (GeneRefPtr)NULL;

    grp = GeneRefNew();
    if (grp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (atp == GENE_REF_locus) {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            grp->locus = (CharPtr)av.ptrvalue;
        } else if (atp == GENE_REF_allele) {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            grp->allele = (CharPtr)av.ptrvalue;
        } else if (atp == GENE_REF_desc) {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            grp->desc = (CharPtr)av.ptrvalue;
        } else if (atp == GENE_REF_maploc) {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            grp->maploc = (CharPtr)av.ptrvalue;
        } else if (atp == GENE_REF_pseudo) {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            grp->pseudo = av.boolvalue;
        } else if (atp == GENE_REF_db)
        {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            while ((atp = AsnReadId(aip, amp, atp)) == GENE_REF_db_E)
            {
                anp = ValNodeNew(grp->db);
                if (anp == NULL)
                    goto erret;
                if (grp->db == NULL)
                    grp->db = anp;
                anp->data.ptrvalue = (Pointer)DbtagAsnRead(aip, atp);
                if (anp->data.ptrvalue == NULL)
                    goto erret;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
        }
        else if (atp == GENE_REF_syn)
        {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            while ((atp = AsnReadId(aip, amp, atp)) == GENE_REF_syn_E)
            {
                anp = ValNodeNew(grp->syn);
                if (anp == NULL)
                    goto erret;
                if (grp->syn == NULL)
                    grp->syn = anp;
                if (AsnReadVal(aip, atp, &anp->data) <= 0) goto erret;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
        }
        else if (atp == GENE_REF_locus_tag) {
           if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
           grp->locus_tag = (CharPtr)av.ptrvalue;
        } else if (atp == GENE_REF_formal_name) {
            grp->formal_name = GeneNomenclatureAsnRead(aip, atp);
            if (aip -> io_failure) goto erret;
        }
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return grp;
erret:
    grp = GeneRefFree(grp);
    goto ret;
}

/*****************************************************************************
*
*   OrgRefNew()
*
*****************************************************************************/
NLM_EXTERN OrgRefPtr LIBCALL OrgRefNew (void)
{
    return (OrgRefPtr)MemNew(sizeof(OrgRef));
}

/*****************************************************************************
*
*   OrgRefFree(orp)
*       Frees an OrgRef and associated data
*
*****************************************************************************/
NLM_EXTERN OrgRefPtr LIBCALL OrgRefFree (OrgRefPtr orp)
{
    ValNodePtr anp, next;

    if (orp == NULL)
        return (OrgRefPtr)NULL;

    orp->taxname = MemFree(orp->taxname);
    orp->common = MemFree(orp->common);
    ValNodeFreeData(orp->mod);
    orp->mod = NULL;
    anp = orp->db;
    while (anp != NULL)
    {
        next = anp->next;
        DbtagFree((DbtagPtr)anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    orp->db = NULL;
    ValNodeFreeData(orp->syn);
    orp->syn = NULL;
    OrgNameFree(orp->orgname);
    orp->orgname = NULL;
    return (OrgRefPtr)MemFree(orp);
}

/*****************************************************************************
*
*   OrgRefAsnWrite(orp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (OrgRef ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgRefAsnWrite (OrgRefPtr orp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    ValNodePtr anp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, ORG_REF);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (orp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)orp))
        goto erret;

    if (orp->taxname != NULL)
    {
        av.ptrvalue = orp->taxname;
        if (! AsnWrite(aip, ORG_REF_taxname, &av)) goto erret;
    }
    if (orp->common != NULL)
    {
        av.ptrvalue = orp->common;
        if (! AsnWrite(aip, ORG_REF_common, &av)) goto erret;
    }
    if (orp->mod != NULL)
    {
        anp = orp->mod;
        if (! AsnOpenStruct(aip, ORG_REF_mod, (Pointer)orp->mod))
            goto erret;
        while (anp != NULL)
        {
            av.ptrvalue = anp->data.ptrvalue;
            if (! AsnWrite(aip, ORG_REF_mod_E, &av)) goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, ORG_REF_mod, (Pointer)orp->mod))
            goto erret;
    }
    if (orp->db != NULL)
    {
        anp = orp->db;
        if (! AsnOpenStruct(aip, ORG_REF_db, (Pointer)orp->db))
            goto erret;
        while (anp != NULL)
        {
            if (! DbtagAsnWrite((DbtagPtr)anp->data.ptrvalue, aip, ORG_REF_db_E))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, ORG_REF_db, (Pointer)orp->db))
            goto erret;
    }
    if (orp->syn != NULL)
    {
        anp = orp->syn;
        if (! AsnOpenStruct(aip, ORG_REF_syn, (Pointer)orp->syn))
            goto erret;
        while (anp != NULL)
        {
            av.ptrvalue = anp->data.ptrvalue;
            if (! AsnWrite(aip, ORG_REF_syn_E, &av)) goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, ORG_REF_syn, (Pointer)orp->syn))
            goto erret;
    }
    if (orp->orgname != NULL)
    {
        if (! OrgNameAsnWrite(orp->orgname, aip, ORG_REF_orgname))
            goto erret;
    }
    if (! AsnCloseStruct(aip, atp, (Pointer)orp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   OrgRefAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN OrgRefPtr LIBCALL OrgRefAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    OrgRefPtr orp;
    ValNodePtr anp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (OrgRefPtr)NULL;
    }

    if (aip == NULL)
        return (OrgRefPtr)NULL;

    if (orig == NULL)           /* OrgRef ::= (self contained) */
        atp = AsnReadId(aip, amp, ORG_REF);
    else
        atp = AsnLinkType(orig, ORG_REF);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (OrgRefPtr)NULL;

    orp = OrgRefNew();
    if (orp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (atp == ORG_REF_orgname)
        {
            orp->orgname = OrgNameAsnRead(aip, atp);
            if (orp->orgname == NULL) goto erret;
        } else {
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            if (atp == ORG_REF_taxname)
                orp->taxname = (CharPtr)av.ptrvalue;
            else if (atp == ORG_REF_common)
                orp->common = (CharPtr)av.ptrvalue;
            else if (atp == ORG_REF_mod_E)
            {
                anp = ValNodeNew(orp->mod);
                if (anp == NULL)
                    goto erret;
                anp->data.ptrvalue = av.ptrvalue;
                if (orp->mod == NULL)
                    orp->mod = anp;
            }
            else if (atp == ORG_REF_db)
            {
                while ((atp = AsnReadId(aip, amp, atp)) == ORG_REF_db_E)
                {
                    anp = ValNodeNew(orp->db);
                    if (anp == NULL)
                        goto erret;
                    if (orp->db == NULL)
                        orp->db = anp;
                    anp->data.ptrvalue = (Pointer)DbtagAsnRead(aip, atp);
                    if (anp->data.ptrvalue == NULL)
                        goto erret;
                }
                if (atp == NULL)
                    goto erret;
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
            }
            else if (atp == ORG_REF_syn_E)
            {
                anp = ValNodeNew(orp->syn);
                if (anp == NULL)
                    goto erret;
                anp->data.ptrvalue = av.ptrvalue;
                if (orp->syn == NULL)
                    orp->syn = anp;
            }
        }
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return orp;
erret:
    orp = OrgRefFree(orp);
    goto ret;
}

/*****************************************************************************
*
*   ValNodeStringMatch(vnp1, vnp2)
*
*   This function is used by OrgRefMatch.
*****************************************************************************/
static Boolean ValNodeStringMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return TRUE;
  }
  else if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  else if (vnp1->choice != vnp2->choice)
  {
    return FALSE;
  }
  else if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0)
  {
    return FALSE;
  }
  else
  {
    return ValNodeStringMatch (vnp1->next, vnp2->next);
  }
}


/*****************************************************************************
*
*   ValNodeDbtagMatch(vnp1, vnp2)
*
*   This function is used by OrgRefMatch.
*****************************************************************************/
static Boolean ValNodeDbtagMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return TRUE;
  }
  else if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  else if (vnp1->choice != vnp2->choice)
  {
    return FALSE;
  }
  else if (!DbtagMatch (vnp1->data.ptrvalue, vnp2->data.ptrvalue))
  {
    return FALSE;
  }
  else
  {
    return ValNodeDbtagMatch (vnp1->next, vnp2->next);
  }
}

/*****************************************************************************
*
*   OrgRefMatch(orp1, orp2)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgRefMatch (OrgRefPtr orp1, OrgRefPtr orp2)
{
  if (orp1 == NULL && orp2 == NULL)
  {
    return TRUE;
  }
  else if (orp1 == NULL || orp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (orp1->taxname, orp2->taxname) != 0
           || StringCmp (orp1->common, orp2->common) != 0)
  {
    return FALSE;
  }
  else if (!ValNodeStringMatch (orp1->syn, orp2->syn)
          || ! ValNodeDbtagMatch (orp1->db, orp2->db))
  {
    return FALSE;
  }
  else if (! OrgNameMatch (orp1->orgname, orp2->orgname))
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }  
}


/*****************************************************************************
*
*   OrgNameNew()
*
*****************************************************************************/
NLM_EXTERN OrgNamePtr LIBCALL OrgNameNew (void)
{
    return (OrgNamePtr)MemNew(sizeof(OrgName));
}

/*****************************************************************************
*
*   OrgNameFree(onp)
*       Frees an OrgName and associated data
*
*****************************************************************************/
NLM_EXTERN OrgNamePtr LIBCALL OrgNameFree (OrgNamePtr onp)
{
    if (onp == NULL)
        return (OrgNamePtr)NULL;

    switch (onp->choice)
    {
        case 1:      /* binomial name */
        case 4:      /* named hydrid */
            BinomialOrgNameFree((BinomialOrgNamePtr)(onp->data));
            break;
        case 2:      /* virus name */
            MemFree(onp->data);
            break;
        case 3:   /* hybrid */
            OrgNameSetFree((OrgNamePtr)(onp->data));
            break;
        case 5:   /* partially characterized name */
            TaxElementSetFree((TaxElementPtr)(onp->data));
            break;

    }
    MemFree(onp->attrib);
    OrgModSetFree(onp->mod);
     MemFree(onp->lineage);
     MemFree(onp->div);
    return (OrgNamePtr)MemFree(onp);
}

/*****************************************************************************
*
*   OrgNameAsnWrite(onp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (OrgName ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgNameAsnWrite (OrgNamePtr onp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, ORGNAME);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (onp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)onp))
        goto erret;

    if (onp->choice != 0)
    {
    if (! AsnWriteChoice(aip, ORGNAME_name, (Int2)onp->choice, &av)) goto erret;


    switch (onp->choice)
    {
        case 1:      /* binomial name */
            if (! BinomialOrgNameAsnWrite((BinomialOrgNamePtr)(onp->data), aip, ORGNAME_name_binomial))
                goto erret;
            break;
        case 2:      /* virus name */
            av.ptrvalue = onp->data;
            if (AsnWrite(aip, ORGNAME_name_virus, &av) <= 0)
                goto erret;
            break;
        case 3:   /* hybrid */
            if (! OrgNameSetAsnWrite((OrgNamePtr)(onp->data), aip, ORGNAME_name_hybrid, MULTIORGNAME_E))
                goto erret;
            break;
        case 4:      /* named hybrid */
            if (! BinomialOrgNameAsnWrite((BinomialOrgNamePtr)(onp->data), aip, ORGNAME_name_namedhybrid))
                goto erret;
            break;
        case 5:   /* partially characterized name */
            if (! TaxElementSetAsnWrite((TaxElementPtr)(onp->data), aip, ORGNAME_name_partial, PARTIALORGNAME_E))
                goto erret;
            break;

    }
    }

    if (onp->attrib != NULL)
    {
        av.ptrvalue = onp->attrib;
        if (! AsnWrite(aip, ORGNAME_attrib, &av)) goto erret;
    }
    if (onp->mod != NULL)
    {
        if (! OrgModSetAsnWrite(onp->mod, aip, ORGNAME_mod, ORGNAME_mod_E))
            goto erret;
    }

    if (onp->lineage != NULL)
    {
          if (aip->spec_version == 3)    /* ASN3 strip new value */
          {
              ErrPostEx(SEV_ERROR,0,0,"ASN3: OrgName.lineage stripped");
          }
          else
          {
        av.ptrvalue = onp->lineage;
        if (! AsnWrite(aip, ORGNAME_lineage, &av)) goto erret;
         }      
    }

    if (onp->gcode)
    {
          if (aip->spec_version == 3)    /* ASN3 strip new value */
          {
              ErrPostEx(SEV_ERROR,0,0,"ASN3: OrgName.gcode stripped");
          }
          else
          {
        av.intvalue = (Int4)(onp->gcode);
        if (! AsnWrite(aip, ORGNAME_gcode, &av)) goto erret;
         }      
    }

    if (onp->mgcode)
    {
          if (aip->spec_version == 3)    /* ASN3 strip new value */
          {
              ErrPostEx(SEV_ERROR,0,0,"ASN3:  OrgName.mgcode stripped");
          }
          else
          {
        av.intvalue = (Int4)(onp->mgcode);
        if (! AsnWrite(aip, ORGNAME_mgcode, &av)) goto erret;
         }      
    }

    if (onp->div != NULL)
    {
          if (aip->spec_version > 0 && aip->spec_version < 5)    /* ASN4 strip new value */
          {
              ErrPostEx(SEV_ERROR,0,0,"ASN3: OrgName.lineage stripped");
          }
          else
          {
        av.ptrvalue = onp->div;
        if (! AsnWrite(aip, ORGNAME_div, &av)) goto erret;
         }      
    }

    if (onp->pgcode)
    {
          if (aip->spec_version == 3)    /* ASN3 strip new value */
          {
              ErrPostEx(SEV_ERROR,0,0,"ASN3:  OrgName.pgcode stripped");
          }
          else
          {
        av.intvalue = (Int4)(onp->pgcode);
        if (! AsnWrite(aip, ORGNAME_pgcode, &av)) goto erret;
         }      
    }


    if (! AsnCloseStruct(aip, atp, (Pointer)onp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   OrgNameAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN OrgNamePtr LIBCALL OrgNameAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    OrgNamePtr onp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (OrgNamePtr)NULL;
    }

    if (aip == NULL)
        return (OrgNamePtr)NULL;

    if (orig == NULL)           /* OrgName ::= (self contained) */
        atp = AsnReadId(aip, amp, ORGNAME);
    else
        atp = AsnLinkType(orig, ORGNAME);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (OrgNamePtr)NULL;

    onp = OrgNameNew();
    if (onp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (atp == ORGNAME_name_binomial)
        {
            onp->choice = 1;
            onp->data = (Pointer)BinomialOrgNameAsnRead(aip, atp);
            if (onp->data == NULL) goto erret;
        }
        else if (atp == ORGNAME_name_hybrid)
        {
            onp->choice = 3;
            onp->data = (Pointer)OrgNameSetAsnRead(aip, atp, MULTIORGNAME_E);
            if (onp->data == NULL) goto erret;
        }
        else if (atp == ORGNAME_name_namedhybrid)
        {
            onp->choice = 4;
            onp->data = (Pointer)BinomialOrgNameAsnRead(aip, atp);
            if (onp->data == NULL) goto erret;
        }
        else if (atp == ORGNAME_name_partial)
        {
            onp->choice = 5;
            onp->data = (Pointer)TaxElementSetAsnRead(aip, atp, PARTIALORGNAME_E);
            if (onp->data == NULL) goto erret;
        }
        else if (atp == ORGNAME_mod)
        {
            onp->mod = OrgModSetAsnRead(aip, atp, ORGNAME_mod_E);
            if (onp->mod == NULL) goto erret;
        }
        else
        {
         if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
         if (atp == ORGNAME_name_virus)
            {
                onp->choice = 2;
                onp->data = av.ptrvalue;
            }
            else if (atp == ORGNAME_attrib)
                onp->attrib = (CharPtr)av.ptrvalue;
            else if (atp == ORGNAME_lineage)
            {
               if (aip->spec_version == 3)    /* ASN3 strip new value */
              {
                  ErrPostEx(SEV_ERROR,0,0,"ASN3:  OrgName.lineage stripped");
                MemFree((CharPtr)(av.ptrvalue));
              }
              else
                onp->lineage = (CharPtr)av.ptrvalue;
           }
            else if (atp == ORGNAME_gcode)
            {
               if (aip->spec_version == 3)    /* ASN3 strip new value */
              {
                  ErrPostEx(SEV_ERROR,0,0,"ASN3:  OrgName.gcode stripped");
              }
              else
                onp->gcode = (Uint1)av.intvalue;
           }
            else if (atp == ORGNAME_mgcode)
            {
               if (aip->spec_version == 3)    /* ASN3 strip new value */
              {
                  ErrPostEx(SEV_ERROR,0,0,"ASN3:  OrgName.mgcode stripped");
              }
              else
                onp->mgcode = (Uint1)av.intvalue;
           }
            else if (atp == ORGNAME_div)
            {
               if (aip->spec_version > 0 && aip->spec_version < 5)    /* ASN4 strip new value */
              {
                  ErrPostEx(SEV_ERROR,0,0,"ASN4:  OrgName.div stripped");
              }
              else
                onp->div = (CharPtr)(av.ptrvalue);
           }
            else if (atp == ORGNAME_pgcode)
            {
               if (aip->spec_version == 3)    /* ASN3 strip new value */
              {
                  ErrPostEx(SEV_ERROR,0,0,"ASN3:  OrgName.pgcode stripped");
              }
              else
                onp->pgcode = (Uint1)av.intvalue;
           }
        }
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return onp;
erret:
    onp = OrgNameFree(onp);
    goto ret;
}

/*****************************************************************************
*
*   OrgNameSetFree (onp)
*
*****************************************************************************/
NLM_EXTERN OrgNamePtr LIBCALL OrgNameSetFree (OrgNamePtr onp)
{
    OrgNamePtr next;

    while (onp != NULL)
    {
        next = onp->next;
        OrgNameFree(onp);
        onp = next;
    }
    return onp;
}

/*****************************************************************************
*
*   OrgNameSetAsnWrite(onp, aip, set, element)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgNameSetAsnWrite (OrgNamePtr onp, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    AsnTypePtr atp;
    OrgNamePtr oldonp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(element, ORGNAME);   /* link local tree */
    if (atp == NULL)
        return FALSE;
    oldonp = onp;

    if (! AsnOpenStruct(aip, set, (Pointer)oldonp))
        goto erret;

    while (onp != NULL)
    {
        if (! OrgNameAsnWrite(onp, aip, atp))
            goto erret;
        onp = onp->next;
    }

    if (! AsnCloseStruct(aip, set, (Pointer)oldonp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(element);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   OrgNameSetAsnRead(aip, set, element)
*
*****************************************************************************/
NLM_EXTERN OrgNamePtr LIBCALL OrgNameSetAsnRead (AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    DataVal av;
    AsnTypePtr atp;
    OrgNamePtr onp, curr = NULL, first = NULL;

    if (aip == NULL)
        return (OrgNamePtr)NULL;

    AsnLinkType(element, ORGNAME);    /* link in local tree */
    atp = set;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) == element)
    {
        onp = OrgNameAsnRead(aip, atp);
        if (onp == NULL)
            goto erret;
        if (first == NULL)
            first = onp;
        else
            curr->next = onp;
        curr = onp;
    }
    if (atp == NULL)
        goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    if (first == NULL)
        ErrPostEx(SEV_WARNING, 0, 0, "Empty SET OF OrgName. line %ld", (long) aip->linenumber);
ret:
    AsnUnlinkType(element);       /* unlink local tree */
    return first;
erret:
    first = OrgNameSetFree(first);
    goto ret;
}

/*****************************************************************************
*
*   OrgNameMatch (onp1, onp2)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgNameMatch (OrgNamePtr onp1, OrgNamePtr onp2)
{
    if (onp1 == NULL && onp2 == NULL)
  {
    return TRUE;
  }
  else if (onp1 == NULL || onp2 == NULL)
  {
    return FALSE;
  }
  else if (onp1->choice != onp2->choice
           || onp1->gcode != onp2->gcode
           || onp2->mgcode != onp2->mgcode
           || StringCmp (onp1->attrib, onp2->attrib) != 0
           || StringCmp (onp1->lineage, onp2->lineage) != 0
           || StringCmp (onp1->div, onp2->div) != 0
           || ! OrgModSetMatch (onp1->mod, onp2->mod)
           || ! OrgNameMatch (onp1->next, onp2->next))
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }  
}

/*****************************************************************************
*
*   ProtRefNew()
*
*****************************************************************************/
NLM_EXTERN ProtRefPtr LIBCALL ProtRefNew (void)
{
    return (ProtRefPtr)MemNew(sizeof(ProtRef));
}

/*****************************************************************************
*
*   ProtRefFree(prp)
*       Frees an ProtRef and associated data
*
*****************************************************************************/
NLM_EXTERN ProtRefPtr LIBCALL ProtRefFree (ProtRefPtr prp)
{
    ValNodePtr anp, next;

    if (prp == NULL)
        return (ProtRefPtr)NULL;

    MemFree(prp->desc);
    anp = prp->name;
    while (anp != NULL)
    {
        next = anp->next;
        MemFree(anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    anp = prp->ec;
    while (anp != NULL)
    {
        next = anp->next;
        MemFree(anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    anp = prp->activity;
    while (anp != NULL)
    {
        next = anp->next;
        MemFree(anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    anp = prp->db;
    while (anp != NULL)
    {
        next = anp->next;
        DbtagFree((DbtagPtr)anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    return (ProtRefPtr)MemFree(prp);
}

/*****************************************************************************
*
*   ProtRefDup(prp)
*       Duplicates an ProtRef and associated data
*
*****************************************************************************/
NLM_EXTERN ProtRefPtr LIBCALL ProtRefDup (ProtRefPtr prp)
{
    ValNodePtr anp, vnp;
    ProtRefPtr npp;

    if (prp == NULL)
        return (ProtRefPtr)NULL;

    npp = ProtRefNew();
    if (npp == NULL) return npp;

    if (prp->desc != NULL)
        npp->desc = StringSave(prp->desc);
    anp = prp->name;
    while (anp != NULL)
    {
        vnp = ValNodeNew(npp->name);
        if (vnp == NULL) return ProtRefFree(prp);
        if (npp->name == NULL)
            npp->name = vnp;
        vnp->data.ptrvalue = StringSave((CharPtr)anp->data.ptrvalue);
        anp = anp->next;
    }
    anp = prp->ec;
    while (anp != NULL)
    {
        vnp = ValNodeNew(npp->ec);
        if (vnp == NULL) return ProtRefFree(prp);
        if (npp->ec == NULL)
            npp->ec = vnp;
        vnp->data.ptrvalue = StringSave((CharPtr)anp->data.ptrvalue);
        anp = anp->next;
    }
    anp = prp->activity;
    while (anp != NULL)
    {
        vnp = ValNodeNew(npp->activity);
        if (vnp == NULL) return ProtRefFree(prp);
        if (npp->activity == NULL)
            npp->activity = vnp;
        vnp->data.ptrvalue = StringSave((CharPtr)anp->data.ptrvalue);
        anp = anp->next;
    }
    anp = prp->db;
    while (anp != NULL)
    {
        vnp = ValNodeNew(npp->db);
        if (vnp == NULL) return ProtRefFree(prp);
        if (npp->db == NULL)
            npp->db = vnp;
        vnp->data.ptrvalue = DbtagDup((DbtagPtr)anp->data.ptrvalue);
        anp = anp->next;
    }
    npp->processed = prp->processed;
    return npp;
}

/*****************************************************************************
*
*   ProtRefAsnWrite(prp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (ProtRef ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL ProtRefAsnWrite (ProtRefPtr prp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    ValNodePtr anp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, PROT_REF);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (prp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)prp))
        goto erret;

    if (prp->name != NULL)
    {
        if (! AsnOpenStruct(aip, PROT_REF_name, (Pointer)prp->name))
            goto erret;
        anp = prp->name;
        while (anp != NULL)
        {
            if (! AsnWrite(aip, PROT_REF_name_E, &anp->data))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, PROT_REF_name, (Pointer)prp->name))
            goto erret;
    }
    if (prp->desc != NULL)
    {
        av.ptrvalue = prp->desc;
        if (! AsnWrite(aip, PROT_REF_desc, &av)) goto erret;
    }
    if (prp->ec != NULL)
    {
        if (! AsnOpenStruct(aip, PROT_REF_ec, (Pointer)prp->ec))
            goto erret;
        anp = prp->ec;
        while (anp != NULL)
        {
            if (! AsnWrite(aip, PROT_REF_ec_E, &anp->data))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, PROT_REF_ec, (Pointer)prp->ec))
            goto erret;
    }
    if (prp->activity != NULL)
    {
        if (! AsnOpenStruct(aip, PROT_REF_activity, (Pointer)prp->activity))
            goto erret;
        anp = prp->activity;
        while (anp != NULL)
        {
            if (! AsnWrite(aip, PROT_REF_activity_E, &anp->data))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, PROT_REF_activity, (Pointer)prp->activity))
            goto erret;
    }
    if (prp->db != NULL)
    {
        anp = prp->db;
        if (! AsnOpenStruct(aip, PROT_REF_db, (Pointer)prp->db))
            goto erret;
        while (anp != NULL)
        {
            if (! DbtagAsnWrite((DbtagPtr)anp->data.ptrvalue, aip, PROT_REF_db_E))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, PROT_REF_db, (Pointer)prp->db))
            goto erret;
    }
    if (prp->processed)
    {
        if ((prp->processed < 3) || (aip->spec_version == 0 || aip->spec_version > 4))
        {
            av.intvalue = (Int4)(prp->processed);
            if (AsnWrite(aip, PROT_REF_processed, &av) <= 0) goto erret;
        }
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)prp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   ProtRefAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN ProtRefPtr LIBCALL ProtRefAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    ProtRefPtr prp;
    ValNodePtr anp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (ProtRefPtr)NULL;
    }

    if (aip == NULL)
        return (ProtRefPtr)NULL;

    if (orig == NULL)           /* ProtRef ::= (self contained) */
        atp = AsnReadId(aip, amp, PROT_REF);
    else
        atp = AsnLinkType(orig, PROT_REF);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (ProtRefPtr)NULL;

    prp = ProtRefNew();
    if (prp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == PROT_REF_desc)
            prp->desc = (CharPtr)av.ptrvalue;
        else if (atp == PROT_REF_name)
        {
            while ((atp = AsnReadId(aip, amp, atp)) == PROT_REF_name_E)
            {
                anp = ValNodeNew(prp->name);
                if (anp == NULL)
                    goto erret;
                if (prp->name == NULL)
                    prp->name = anp;
                if (AsnReadVal(aip, atp, &anp->data) <= 0) goto erret;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
        }
        else if (atp == PROT_REF_ec)
        {
            while ((atp = AsnReadId(aip, amp, atp)) == PROT_REF_ec_E)
            {
                anp = ValNodeNew(prp->ec);
                if (anp == NULL)
                    goto erret;
                if (prp->ec == NULL)
                    prp->ec = anp;
                if (AsnReadVal(aip, atp, &anp->data) <= 0) goto erret;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
        }
        else if (atp == PROT_REF_activity)
        {
            while ((atp = AsnReadId(aip, amp, atp)) == PROT_REF_activity_E)
            {
                anp = ValNodeNew(prp->activity);
                if (anp == NULL)
                    goto erret;
                if (prp->activity == NULL)
                    prp->activity = anp;
                if (AsnReadVal(aip, atp, &anp->data) <= 0) goto erret;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
        }
        else if (atp == PROT_REF_db)
        {
            while ((atp = AsnReadId(aip, amp, atp)) == PROT_REF_db_E)
            {
                anp = ValNodeNew(prp->db);
                if (anp == NULL)
                    goto erret;
                if (prp->db == NULL)
                    prp->db = anp;
                anp->data.ptrvalue = (Pointer)DbtagAsnRead(aip, atp);
                if (anp->data.ptrvalue == NULL)
                    goto erret;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;  /* end SET OF */
        }
        else if (atp == PROT_REF_processed)
            prp->processed = (Uint1)(av.intvalue);
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return prp;
erret:
    prp = ProtRefFree(prp);
    goto ret;
}

/*****************************************************************************
*
*   RsiteRefFree(rrp)
*       Frees an RsiteRef and associated data
*
*****************************************************************************/
NLM_EXTERN RsiteRefPtr LIBCALL RsiteRefFree (RsiteRefPtr rrp)
{
    if (rrp == NULL)
        return (RsiteRefPtr)NULL;

    if (rrp->choice == 1)   /* str */
        MemFree(rrp->data.ptrvalue);
    else
        DbtagFree((DbtagPtr)rrp->data.ptrvalue);
    return (RsiteRefPtr)MemFree(rrp);
}

/*****************************************************************************
*
*   RsiteRefAsnWrite(rrp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (RsiteRef ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL RsiteRefAsnWrite (RsiteRefPtr rrp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, RSITE_REF);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (rrp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnWrite(aip, atp, &av)) goto erret;     /* CHOICE */
    switch (rrp->choice)
    {
        case 1:
            if (! AsnWrite(aip, RSITE_REF_str, &rrp->data))
                goto erret;
            break;
        case 2:
            if (! DbtagAsnWrite((DbtagPtr)rrp->data.ptrvalue, aip, RSITE_REF_db))
                goto erret;
            break;
    }
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   RsiteRefAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN RsiteRefPtr LIBCALL RsiteRefAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    RsiteRefPtr rrp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (RsiteRefPtr)NULL;
    }

    if (aip == NULL)
        return (RsiteRefPtr)NULL;

    if (orig == NULL)           /* RsiteRef ::= (self contained) */
        atp = AsnReadId(aip, amp, RSITE_REF);
    else
        atp = AsnLinkType(orig, RSITE_REF);    /* link in local tree */
    if (atp == NULL)
        return (RsiteRefPtr)NULL;

    rrp = ValNodeNew(NULL);
    if (rrp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the CHOICE */
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;   /* which CHOICE */
    if (atp == NULL)
        goto erret;
    if (atp == RSITE_REF_str)
    {
        rrp->choice = 1;
        if (AsnReadVal(aip, atp, &rrp->data) <= 0)
            goto erret;
    }
    else if (atp == RSITE_REF_db)
    {
        rrp->choice = 2;
        rrp->data.ptrvalue = (Pointer) DbtagAsnRead(aip, atp);
        if (rrp->data.ptrvalue == NULL)
            goto erret;
    }
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return rrp;
erret:
    rrp = RsiteRefFree(rrp);
    goto ret;
}

/*****************************************************************************
*
*   SeqFeatIdFree(cp)
*       Frees one SeqFeatId's associated data
*       Does NOT free the Choice itself
*
*****************************************************************************/
NLM_EXTERN void LIBCALL SeqFeatIdFree (ChoicePtr cp)
{
    Pointer pnt;

    if (cp == NULL)
        return;

    pnt = cp->value.ptrvalue;
    switch (cp->choice)
    {
        case 1:     /* gibb id */
            cp->value.intvalue = 0;
            break;
        case 2:     /* giim id */
            cp->value.ptrvalue = GiimFree((GiimPtr)pnt);
            break;
        case 3:
            cp->value.ptrvalue = ObjectIdFree((ObjectIdPtr)pnt);
            break;
        case 4:
            cp->value.ptrvalue = DbtagFree((DbtagPtr)pnt);
            break;
    }
    cp->choice = 0;
    return;
}

/*****************************************************************************
*
*   SeqFeatIdDup(dest, src)
*       Duplicates one SeqFeatId's associated data
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatIdDup (ChoicePtr dest, ChoicePtr src)
{
    GiimPtr ga, gb;

    if (src == NULL)
        return FALSE;

    dest->choice = src->choice;
    switch (src->choice)
    {
        case 1:     /* gibb */
            dest->value.intvalue = src->value.intvalue;
            break;
        case 2:     /* giim id */
            ga = (GiimPtr) src->value.ptrvalue;
            gb = GiimNew();
            if (gb == NULL) return FALSE;
            gb->id = ga->id;
            gb->db = StringSave(ga->db);
            gb->release = StringSave(ga->release);
            dest->value.ptrvalue = gb;
            break;
        case 3:
            dest->value.ptrvalue = ObjectIdDup((ObjectIdPtr)src->value.ptrvalue);
            break;
        case 4:
            dest->value.ptrvalue = DbtagDup((DbtagPtr)src->value.ptrvalue);
            break;
    }
    return TRUE;
}

/*****************************************************************************
*
*   SeqFeatIdAsnWrite(cp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqFeatId ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatIdAsnWrite (ChoicePtr cp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, FEAT_ID);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (cp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (cp->choice)              /* feat-id present */
    {
        if (! AsnWrite(aip, atp, &av)) goto erret;   /* CHOICE */
        if (cp->choice == 1)           /* gibb */
        {
            av.intvalue = cp->value.intvalue;
            if (! AsnWrite(aip, FEAT_ID_gibb, &av)) goto erret;
        }
        else if (cp->choice == 2)      /* giim */
        {
            if (! GiimAsnWrite((GiimPtr)cp->value.ptrvalue, aip, FEAT_ID_giim))
                goto erret;
        }
        else if (cp->choice == 3)      /* local */
        {
            if (! ObjectIdAsnWrite((ObjectIdPtr)cp->value.ptrvalue, aip, FEAT_ID_local))
                goto erret;
        }
        else if (cp->choice == 4)      /* general */
        {
            if (! DbtagAsnWrite((DbtagPtr)cp->value.ptrvalue, aip, FEAT_ID_general))
                goto erret;
        }
    }
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SeqFeatIdAsnRead(aip, atp, cp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SeqFeatIdAsnRead (AsnIoPtr aip, AsnTypePtr orig, ChoicePtr cp)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return retval;
    }

    if (aip == NULL)
        return retval;

    if (orig == NULL)           /* SeqFeatId ::= (self contained) */
        atp = AsnReadId(aip, amp, FEAT_ID);
    else
        atp = AsnLinkType(orig, FEAT_ID);    /* link in local tree */
    if (atp == NULL)
        return retval;

    if (AsnReadVal(aip, atp, &av) <= 0)    /* read CHOICE */
        goto erret;
    atp = AsnReadId(aip, amp, atp); if (atp == NULL) goto erret;
    if (atp == FEAT_ID_gibb)
    {
        cp->choice = 1;
        if (AsnReadVal(aip, atp, &av) <= 0)
            goto erret;
        cp->value.intvalue = av.intvalue;
    }
    else if (atp == FEAT_ID_giim)
    {
        cp->choice = 2;
        cp->value.ptrvalue = (Pointer)GiimAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == FEAT_ID_local)
    {
        cp->choice = 3;
        cp->value.ptrvalue = (Pointer)ObjectIdAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    else if (atp == FEAT_ID_general)
    {
        cp->choice = 4;
        cp->value.ptrvalue = (Pointer)DbtagAsnRead(aip, atp);
        if (cp->value.ptrvalue == NULL)
            goto erret;
    }
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   TxinitNew()
*
*****************************************************************************/
NLM_EXTERN TxinitPtr LIBCALL TxinitNew (void)
{
    TxinitPtr txp;

    txp = (TxinitPtr)MemNew(sizeof(Txinit));
    if (txp == NULL) return txp;
    txp->inittype = 255;
    return txp;
}

/*****************************************************************************
*
*   TxinitFree(txp)
*       Frees a Txinit and associated data
*
*****************************************************************************/
NLM_EXTERN TxinitPtr LIBCALL TxinitFree (TxinitPtr txp)
{
    ValNodePtr anp, next;
    TxEvidencePtr txep, txnext;

    if (txp == NULL)
        return (TxinitPtr)NULL;

    MemFree(txp->name);
    ValNodeFreeData(txp->syn);
    anp = txp->gene;
    while (anp != NULL)
    {
        next = anp->next;
        GeneRefFree((GeneRefPtr)anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    anp = txp->protein;
    while (anp != NULL)
    {
        next = anp->next;
        ProtRefFree((ProtRefPtr)anp->data.ptrvalue);
        MemFree(anp);
        anp = next;
    }
    ValNodeFreeData(txp->rna);
    MemFree(txp->expression);
    MemFree(txp->txdescr);
    OrgRefFree(txp->txorg);
    txep = txp->evidence;
    while (txep != NULL)
    {
        txnext = txep->next;
        MemFree(txep);
        txep = txnext;
    }

    return (TxinitPtr)MemFree(txp);
}

/*****************************************************************************
*
*   TxinitAsnWrite(txp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (Txinit ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL TxinitAsnWrite (TxinitPtr txp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    ValNodePtr anp;
    TxEvidencePtr txep;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, TXINIT);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (txp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)txp))
        goto erret;

    av.ptrvalue = (Pointer)txp->name;
    if (! AsnWrite(aip, TXINIT_name, &av)) goto erret;
    if (txp->syn != NULL)
    {
        anp = txp->syn;
        if (! AsnOpenStruct(aip, TXINIT_syn, (Pointer)txp->syn))
            goto erret;
        while (anp != NULL)
        {
            av.ptrvalue = anp->data.ptrvalue;
            if (! AsnWrite(aip, TXINIT_syn_E, &av)) goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, TXINIT_syn, (Pointer)txp->syn))
            goto erret;
    }
    if (txp->gene != NULL)
    {
        anp = txp->gene;
        if (! AsnOpenStruct(aip, TXINIT_gene, (Pointer)txp->gene))
            goto erret;
        while (anp != NULL)
        {
            if (! GeneRefAsnWrite((GeneRefPtr)anp->data.ptrvalue, aip, TXINIT_gene_E))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, TXINIT_gene, (Pointer)txp->gene))
            goto erret;
    }
    if (txp->protein != NULL)
    {
        anp = txp->protein;
        if (! AsnOpenStruct(aip, TXINIT_protein, (Pointer)txp->protein))
            goto erret;
        while (anp != NULL)
        {
            if (! ProtRefAsnWrite((ProtRefPtr)anp->data.ptrvalue, aip, TXINIT_protein_E))
                goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, TXINIT_protein, (Pointer)txp->protein))
            goto erret;
    }
    if (txp->rna != NULL)
    {
        anp = txp->rna;
        if (! AsnOpenStruct(aip, TXINIT_rna, (Pointer)txp->rna))
            goto erret;
        while (anp != NULL)
        {
            av.ptrvalue = anp->data.ptrvalue;
            if (! AsnWrite(aip, TXINIT_rna_E, &av)) goto erret;
            anp = anp->next;
        }
        if (! AsnCloseStruct(aip, TXINIT_rna, (Pointer)txp->rna))
            goto erret;
    }
    if (txp->expression != NULL)
    {
        av.ptrvalue = (Pointer)txp->expression;
        if (! AsnWrite(aip, TXINIT_expression, &av)) goto erret;
    }

    av.intvalue = (Int4)txp->txsystem;
    if (! AsnWrite(aip, TXINIT_txsystem, &av)) goto erret;

    if (txp->txdescr != NULL)
    {
        av.ptrvalue = (Pointer)txp->txdescr;
        if (! AsnWrite(aip, TXINIT_txdescr, &av)) goto erret;
    }

    if (txp->txorg != NULL)
    {
        if (! OrgRefAsnWrite(txp->txorg, aip, TXINIT_txorg))
            goto erret;
    }

    if (txp->mapping_precise)
    {
        av.boolvalue = TRUE;
        if (! AsnWrite(aip, TXINIT_mapping_precise, &av)) goto erret;
    }

    if (txp->location_accurate)
    {
        av.boolvalue = TRUE;
        if (! AsnWrite(aip, TXINIT_location_accurate, &av)) goto erret;
    }

    if (txp->inittype != 255)
    {
        av.intvalue = (Int4) txp->inittype;
        if (! AsnWrite(aip, TXINIT_inittype, &av)) goto erret;
    }

    if (txp->evidence != NULL)
    {
        txep = txp->evidence;
        if (! AsnOpenStruct(aip, TXINIT_evidence, (Pointer)txp->evidence))
            goto erret;
        while (txep != NULL)
        {
            if (! AsnOpenStruct(aip, TXINIT_evidence_E, (Pointer)txep))
                goto erret;
            av.intvalue = (Int4) txep->exp_code;
            if (! AsnWrite(aip, TX_EVIDENCE_exp_code, &av)) goto erret;
            if (txep->exp_sys != 1)
            {
                av.intvalue = (Int4) txep->exp_sys;
                if (! AsnWrite(aip, TX_EVIDENCE_expression_system, &av)) goto erret;
            }
            if (txep->low_prec_data)
            {
                av.boolvalue = TRUE;
                if (! AsnWrite(aip, TX_EVIDENCE_low_prec_data, &av)) goto erret;
            }
            if (txep->from_homolog)
            {
                av.boolvalue = TRUE;
                if (! AsnWrite(aip, TX_EVIDENCE_from_homolog, &av)) goto erret;
            }
            if (! AsnCloseStruct(aip, TXINIT_evidence_E, (Pointer)txep))
                goto erret;
            txep = txep->next;
        }
        if (! AsnCloseStruct(aip, TXINIT_evidence, (Pointer)txp->evidence))
            goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)txp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   TxinitAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN TxinitPtr LIBCALL TxinitAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    TxinitPtr txp;
    ValNodePtr curr, last;
    TxEvidencePtr txep, txlast;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (TxinitPtr)NULL;
    }

    if (aip == NULL)
        return (TxinitPtr)NULL;

    if (orig == NULL)           /* Txinit ::= (self contained) */
        atp = AsnReadId(aip, amp, TXINIT);
    else
        atp = AsnLinkType(orig, TXINIT);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (TxinitPtr)NULL;

    txp = TxinitNew();
    if (txp == NULL)
        return txp;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == TXINIT_name)
            txp->name = (CharPtr)av.ptrvalue;
        else if (atp == TXINIT_syn)
        {
            last = NULL;
            while ((atp = AsnReadId(aip, amp, atp)) == TXINIT_syn_E)
            {
                curr = ValNodeNew(last);
                if (curr == NULL)
                    goto erret;
                if (AsnReadVal(aip, atp, &curr->data) <= 0) goto erret;
                if (last == NULL)
                    txp->syn = curr;
                last = curr;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        }
        else if (atp == TXINIT_gene)
        {
            last = NULL;
            while ((atp = AsnReadId(aip, amp, atp)) == TXINIT_gene_E)
            {
                curr = ValNodeNew(last);
                if (curr == NULL)
                    goto erret;
                curr->data.ptrvalue = (Pointer)GeneRefAsnRead(aip, atp);
                if (curr->data.ptrvalue == NULL)
                    goto erret;
                if (last == NULL)
                    txp->gene = curr;
                last = curr;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        }
        else if (atp == TXINIT_protein)
        {
            last = NULL;
            while ((atp = AsnReadId(aip, amp, atp)) == TXINIT_protein_E)
            {
                curr = ValNodeNew(last);
                if (curr == NULL)
                    goto erret;
                curr->data.ptrvalue = (Pointer)ProtRefAsnRead(aip, atp);
                if (curr->data.ptrvalue == NULL)
                    goto erret;
                if (last == NULL)
                    txp->protein = curr;
                last = curr;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        }
        else if (atp == TXINIT_rna)
        {
            last = NULL;
            while ((atp = AsnReadId(aip, amp, atp)) == TXINIT_rna_E)
            {
                curr = ValNodeNew(last);
                if (curr == NULL)
                    goto erret;
                if (AsnReadVal(aip, atp, &curr->data) <= 0) goto erret;
                if (last == NULL)
                    txp->rna = curr;
                last = curr;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        }
        else if (atp == TXINIT_expression)
            txp->expression = (CharPtr)av.ptrvalue;
        else if (atp == TXINIT_txsystem)
            txp->txsystem = (Uint1)av.intvalue;
        else if (atp == TXINIT_txdescr)
            txp->txdescr = (CharPtr)av.ptrvalue;
        else if (atp == TXINIT_txorg)
        {
            txp->txorg = OrgRefAsnRead(aip, atp);
            if (txp->txorg == NULL)
                goto erret;
        }
        else if (atp == TXINIT_mapping_precise)
            txp->mapping_precise = av.boolvalue;
        else if (atp == TXINIT_location_accurate)
            txp->location_accurate = av.boolvalue;
        else if (atp == TXINIT_inittype)
            txp->inittype = (Uint1)av.intvalue;
        else if (atp == TXINIT_evidence)
        {
            txlast = NULL;
            txep = NULL;
            while ((atp = AsnReadId(aip, amp, atp)) != TXINIT_evidence)
            {
                if (atp == NULL)
                    goto erret;
                if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
                if ((atp == TXINIT_evidence_E) && (av.intvalue == START_STRUCT))
                {
                    txep = (TxEvidencePtr)MemNew(sizeof(TxEvidence));
                    if (txep == NULL)
                        goto erret;
                    txep->exp_sys = 1;    /* default */
                    if (txlast != NULL)
                        txlast->next = txep;
                    else
                        txp->evidence = txep;
                    txlast = txep;
                }
                else if (atp == TX_EVIDENCE_exp_code)
                    txep->exp_code = (Uint1) av.intvalue;
                else if (atp == TX_EVIDENCE_expression_system)
                    txep->exp_sys = (Uint1) av.intvalue;
                else if (atp == TX_EVIDENCE_low_prec_data)
                    txep->low_prec_data = av.boolvalue;
                else if (atp == TX_EVIDENCE_from_homolog)
                    txep->from_homolog = av.boolvalue;
            }
            if (atp == NULL)
                goto erret;
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
        }
    }
    if (atp == NULL)
        goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return txp;
erret:
    txp = TxinitFree(txp);
    goto ret;
}


/*****************************************************************************
*
*   OrgModNew()
*
*****************************************************************************/
NLM_EXTERN OrgModPtr LIBCALL OrgModNew (void)
{
    return (OrgModPtr)MemNew(sizeof(OrgMod));
}

/*****************************************************************************
*
*   OrgModFree(omp)
*       Frees an OrgMod and associated data
*
*****************************************************************************/
NLM_EXTERN OrgModPtr LIBCALL OrgModFree (OrgModPtr omp)
{
    if (omp == NULL)
        return (OrgModPtr)NULL;

    MemFree(omp->subname);
    MemFree(omp->attrib);
    return (OrgModPtr)MemFree(omp);
}

/*****************************************************************************
*
*   OrgModAsnWrite(omp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (OrgMod ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgModAsnWrite (OrgModPtr omp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, ORGMOD);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (omp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)omp))
        goto erret;

    av.intvalue = (Int4)(omp->subtype);
    if (AsnWrite(aip, ORGMOD_subtype, &av) <= 0) goto erret;

    if (omp->subname != NULL)
    {
        av.ptrvalue = omp->subname;
        if (! AsnWrite(aip, ORGMOD_subname, &av)) goto erret;
    }

    if (omp->attrib != NULL)
    {
        av.ptrvalue = omp->attrib;
        if (! AsnWrite(aip, ORGMOD_attrib, &av)) goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)omp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   OrgModAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN OrgModPtr LIBCALL OrgModAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    OrgModPtr omp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (OrgModPtr)NULL;
    }

    if (aip == NULL)
        return (OrgModPtr)NULL;

    if (orig == NULL)           /* OrgMod ::= (self contained) */
        atp = AsnReadId(aip, amp, ORGMOD);
    else
        atp = AsnLinkType(orig, ORGMOD);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (OrgModPtr)NULL;

    omp = OrgModNew();
    if (omp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == ORGMOD_subtype)
            omp->subtype = (Uint1)(av.intvalue);
        else if (atp == ORGMOD_subname)
            omp->subname = (CharPtr)(av.ptrvalue);
        else if (atp == ORGMOD_attrib)
            omp->attrib = (CharPtr)(av.ptrvalue);
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return omp;
erret:
    omp = OrgModFree(omp);
    goto ret;
}

/*****************************************************************************
*
*   OrgModSetFree (omp)
*
*****************************************************************************/
NLM_EXTERN OrgModPtr LIBCALL OrgModSetFree (OrgModPtr omp)
{
    OrgModPtr next;

    while (omp != NULL)
    {
        next = omp->next;
        OrgModFree(omp);
        omp = next;
    }
    return omp;
}

/*****************************************************************************
*
*   OrgModSetAsnWrite(omp, aip, set, element)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgModSetAsnWrite (OrgModPtr omp, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    AsnTypePtr atp;
    OrgModPtr oldomp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(element, ORGMOD);   /* link local tree */
    if (atp == NULL)
        return FALSE;
    oldomp = omp;

    if (! AsnOpenStruct(aip, set, (Pointer)oldomp))
        goto erret;

    while (omp != NULL)
    {
        if (! OrgModAsnWrite(omp, aip, atp))
            goto erret;
        omp = omp->next;
    }

    if (! AsnCloseStruct(aip, set, (Pointer)oldomp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(element);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   OrgModSetAsnRead(aip, set, element)
*
*****************************************************************************/
NLM_EXTERN OrgModPtr LIBCALL OrgModSetAsnRead (AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    DataVal av;
    AsnTypePtr atp;
    OrgModPtr omp, curr = NULL, first = NULL;

    if (aip == NULL)
        return (OrgModPtr)NULL;

    AsnLinkType(element, ORGMOD);    /* link in local tree */
    atp = set;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) == element)
    {
        omp = OrgModAsnRead(aip, atp);
        if (omp == NULL)
            goto erret;
        if (first == NULL)
            first = omp;
        else
            curr->next = omp;
        curr = omp;
    }
    if (atp == NULL)
        goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    if (first == NULL)
        ErrPostEx(SEV_WARNING, 0, 0, "Empty SET OF OrgMod. line %ld", (long) aip->linenumber);
ret:
    AsnUnlinkType(element);       /* unlink local tree */
    return first;
erret:
    first = OrgModSetFree(first);
    goto ret;
}

/*****************************************************************************
*
*   OrgModSetMatch (mod1, mod2)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL OrgModSetMatch (OrgModPtr mod1, OrgModPtr mod2)
{
  if (mod1 == NULL && mod2 == NULL)
  {
    return TRUE;
  }
  else if (mod1 == NULL || mod2 == NULL)
  {
    return FALSE;
  }
  else if (mod1->subtype != mod2->subtype
           || StringCmp (mod1->attrib, mod2->attrib) != 0
           || StringCmp (mod1->subname, mod2->subname) != 0)
  {
    return FALSE;
  }
  else
  {
    return OrgModSetMatch (mod1->next, mod2->next);
  }
}

/*****************************************************************************
*
*   TaxElementNew()
*
*****************************************************************************/
NLM_EXTERN TaxElementPtr LIBCALL TaxElementNew (void)
{
    return (TaxElementPtr)MemNew(sizeof(TaxElement));
}

/*****************************************************************************
*
*   TaxElementFree(tep)
*       Frees an TaxElement and associated data
*
*****************************************************************************/
NLM_EXTERN TaxElementPtr LIBCALL TaxElementFree (TaxElementPtr tep)
{
    if (tep == NULL)
        return (TaxElementPtr)NULL;

    MemFree(tep->level);
    MemFree(tep->name);
    return (TaxElementPtr)MemFree(tep);
}

/*****************************************************************************
*
*   TaxElementAsnWrite(tep, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (TaxElement ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL TaxElementAsnWrite (TaxElementPtr tep, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, TAXELEMENT);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (tep == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)tep))
        goto erret;

        av.intvalue = (Int4)(tep->fixed_level);
        if (AsnWrite(aip, TAXELEMENT_fixed_level, &av) <= 0) goto erret;

    if (tep->level != NULL)
    {
        av.ptrvalue = tep->level;
        if (! AsnWrite(aip, TAXELEMENT_level, &av)) goto erret;
    }

    if (tep->name != NULL)
    {
        av.ptrvalue = tep->name;
        if (! AsnWrite(aip, TAXELEMENT_name, &av)) goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)tep))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   TaxElementAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN TaxElementPtr LIBCALL TaxElementAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    TaxElementPtr tep;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (TaxElementPtr)NULL;
    }

    if (aip == NULL)
        return (TaxElementPtr)NULL;

    if (orig == NULL)           /* TaxElement ::= (self contained) */
        atp = AsnReadId(aip, amp, TAXELEMENT);
    else
        atp = AsnLinkType(orig, TAXELEMENT);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (TaxElementPtr)NULL;

    tep = TaxElementNew();
    if (tep == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == TAXELEMENT_fixed_level)
            tep->fixed_level = (Uint1)(av.intvalue);
        else if (atp == TAXELEMENT_level)
            tep->level = (CharPtr)(av.ptrvalue);
        else if (atp == TAXELEMENT_name)
            tep->name = (CharPtr)(av.ptrvalue);
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return tep;
erret:
    tep = TaxElementFree(tep);
    goto ret;
}

/*****************************************************************************
*
*   TaxElementSetFree (tep)
*
*****************************************************************************/
NLM_EXTERN TaxElementPtr LIBCALL TaxElementSetFree (TaxElementPtr tep)
{
    TaxElementPtr next;

    while (tep != NULL)
    {
        next = tep->next;
        TaxElementFree(tep);
        tep = next;
    }
    return tep;
}

/*****************************************************************************
*
*   TaxElementSetAsnWrite(tep, aip, set, element)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL TaxElementSetAsnWrite (TaxElementPtr tep, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    AsnTypePtr atp;
    TaxElementPtr oldtep;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(element, TAXELEMENT);   /* link local tree */
    if (atp == NULL)
        return FALSE;
    oldtep = tep;

    if (! AsnOpenStruct(aip, set, (Pointer)oldtep))
        goto erret;

    while (tep != NULL)
    {
        if (! TaxElementAsnWrite(tep, aip, atp))
            goto erret;
        tep = tep->next;
    }

    if (! AsnCloseStruct(aip, set, (Pointer)oldtep))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(element);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   TaxElementSetAsnRead(aip, set, element)
*
*****************************************************************************/
NLM_EXTERN TaxElementPtr LIBCALL TaxElementSetAsnRead (AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    DataVal av;
    AsnTypePtr atp;
    TaxElementPtr tep, curr = NULL, first = NULL;

    if (aip == NULL)
        return (TaxElementPtr)NULL;

    AsnLinkType(element, TAXELEMENT);    /* link in local tree */
    atp = set;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) == element)
    {
        tep = TaxElementAsnRead(aip, atp);
        if (tep == NULL)
            goto erret;
        if (first == NULL)
            first = tep;
        else
            curr->next = tep;
        curr = tep;
    }
    if (atp == NULL)
        goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    if (first == NULL)
        ErrPostEx(SEV_WARNING, 0, 0, "Empty SET OF TaxElement. line %ld", (long) aip->linenumber);
ret:
    AsnUnlinkType(element);       /* unlink local tree */
    return first;
erret:
    first = TaxElementSetFree(first);
    goto ret;
}

/*****************************************************************************
*
*   BinomialOrgNameNew()
*
*****************************************************************************/
NLM_EXTERN BinomialOrgNamePtr LIBCALL BinomialOrgNameNew (void)
{
    return (BinomialOrgNamePtr)MemNew(sizeof(BinomialOrgName));
}

/*****************************************************************************
*
*   BinomialOrgNameFree(bop)
*       Frees an BinomialOrgName and associated data
*
*****************************************************************************/
NLM_EXTERN BinomialOrgNamePtr LIBCALL BinomialOrgNameFree (BinomialOrgNamePtr bop)
{
    if (bop == NULL)
        return (BinomialOrgNamePtr)NULL;

    MemFree(bop->genus);
    MemFree(bop->species);
    MemFree(bop->subspecies);
    return (BinomialOrgNamePtr)MemFree(bop);
}

/*****************************************************************************
*
*   BinomialOrgNameAsnWrite(bop, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (BinomialOrgName ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BinomialOrgNameAsnWrite (BinomialOrgNamePtr bop, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, BINOMIALORGNAME);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (bop == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)bop))
        goto erret;

    if (bop->genus != NULL)
    {
        av.ptrvalue = bop->genus;
        if (! AsnWrite(aip, BINOMIALORGNAME_genus, &av)) goto erret;
    }

    if (bop->species != NULL)
    {
        av.ptrvalue = bop->species;
        if (! AsnWrite(aip, BINOMIALORGNAME_species, &av)) goto erret;
    }

    if (bop->subspecies != NULL)
    {
        av.ptrvalue = bop->subspecies;
        if (! AsnWrite(aip, BINOMIALORGNAME_subspecies, &av)) goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)bop))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   BinomialOrgNameAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN BinomialOrgNamePtr LIBCALL BinomialOrgNameAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    BinomialOrgNamePtr bop;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (BinomialOrgNamePtr)NULL;
    }

    if (aip == NULL)
        return (BinomialOrgNamePtr)NULL;

    if (orig == NULL)           /* BinomialOrgName ::= (self contained) */
        atp = AsnReadId(aip, amp, BINOMIALORGNAME);
    else
        atp = AsnLinkType(orig, BINOMIALORGNAME);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (BinomialOrgNamePtr)NULL;

    bop = BinomialOrgNameNew();
    if (bop == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == BINOMIALORGNAME_genus)
            bop->genus = (CharPtr)(av.ptrvalue);
        else if (atp == BINOMIALORGNAME_species)
            bop->species = (CharPtr)(av.ptrvalue);
        else if (atp == BINOMIALORGNAME_subspecies)
            bop->subspecies = (CharPtr)(av.ptrvalue);
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return bop;
erret:
    bop = BinomialOrgNameFree(bop);
    goto ret;
}


/*****************************************************************************
*
*   BioSourceNew()
*
*****************************************************************************/
NLM_EXTERN BioSourcePtr LIBCALL BioSourceNew (void)
{
    return (BioSourcePtr)MemNew(sizeof(BioSource));
}

/*****************************************************************************
*
*   BioSourceFree(bsp)
*       Frees an BioSource and associated data
*
*****************************************************************************/
NLM_EXTERN BioSourcePtr LIBCALL BioSourceFree (BioSourcePtr bsp)
{
    if (bsp == NULL)
        return (BioSourcePtr)NULL;

    OrgRefFree(bsp->org);
    SubSourceSetFree(bsp->subtype);
    PCRReactionSetFree(bsp->pcr_primers);
    return (BioSourcePtr)MemFree(bsp);
}

/*****************************************************************************
*
*   BioSourceAsnWrite(bsp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (BioSource ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioSourceAsnWrite (BioSourcePtr bsp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, BIOSOURCE);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (bsp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)bsp))
        goto erret;

    if (bsp->genome)
    {
        if ((bsp->genome < 14) || (aip->spec_version == 0 || aip->spec_version > 4))
        {
            av.intvalue = (Int4)bsp->genome;
            if (! AsnWrite(aip, BIOSOURCE_genome, &av)) goto erret;
        }
    }

    if (bsp->origin)
    {
        av.intvalue = (Int4)bsp->origin;
        if (! AsnWrite(aip, BIOSOURCE_origin, &av)) goto erret;
    }

    if (bsp->org != NULL)
    {
        if (! OrgRefAsnWrite(bsp->org, aip, BIOSOURCE_org))
            goto erret;
    }

    if (bsp->subtype != NULL)
    {
        if (! SubSourceSetAsnWrite(bsp->subtype, aip, BIOSOURCE_subtype, BIOSOURCE_subtype_E))
            goto erret;
    }

    if ((bsp->is_focus) && ((aip->spec_version == 0) || (aip->spec_version >= 6)))
    {
        av.ptrvalue = NULL;
        if (! AsnWrite(aip, BIOSOURCE_is_focus, &av))
            goto erret;
    }

    if (bsp->pcr_primers != NULL) {
        if (! PCRReactionSetAsnWrite(bsp->pcr_primers, aip, BIOSOURCE_pcr_primers)) {
            goto erret;
        }
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)bsp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   BioSourceAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN BioSourcePtr LIBCALL BioSourceAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    BioSourcePtr bsp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (BioSourcePtr)NULL;
    }

    if (aip == NULL)
        return (BioSourcePtr)NULL;

    if (orig == NULL)           /* BioSource ::= (self contained) */
        atp = AsnReadId(aip, amp, BIOSOURCE);
    else
        atp = AsnLinkType(orig, BIOSOURCE);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (BioSourcePtr)NULL;

    bsp = BioSourceNew();
    if (bsp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (atp == BIOSOURCE_org)
        {
            bsp->org = OrgRefAsnRead(aip, atp);
            if (bsp->org == NULL) goto erret;
        }
        else if (atp == BIOSOURCE_subtype)
        {
            bsp->subtype = SubSourceSetAsnRead(aip, atp, BIOSOURCE_subtype_E);
            if (bsp->subtype == NULL) goto erret;
        }
        else if (atp == BIOSOURCE_pcr_primers)
        {
            bsp->pcr_primers = PCRReactionSetAsnRead(aip, atp);
            if (bsp->pcr_primers == NULL) goto erret;
        }
        else
        {
            if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
            if (atp == BIOSOURCE_genome)
                bsp->genome = (Uint1)(av.intvalue);
            else if (atp == BIOSOURCE_origin)
                bsp->origin = (Uint1)(av.intvalue);
            else if (atp == BIOSOURCE_is_focus)
                bsp->is_focus = TRUE;
        }
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return bsp;
erret:
    bsp = BioSourceFree(bsp);
    goto ret;
}

/*****************************************************************************
*
*   BioSourceMatch(biop1, biop2)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL BioSourceMatch (BioSourcePtr biop1, BioSourcePtr biop2)
{
  if (biop1 == NULL && biop2 == NULL)
  {
    return TRUE;
  }
  else if (biop1 == NULL || biop2 == NULL)
  {
    return FALSE;
  }
  else if (biop1->origin != biop2->origin
           || biop1->is_focus != biop2->is_focus
           || ! OrgRefMatch (biop1->org, biop2->org)
           || ! SubSourceSetMatch (biop1->subtype, biop2->subtype))
  {
    return FALSE;
  }
  else if (biop1->genome == biop2->genome
           || (biop1->genome == 0 && biop2->genome == 1)
           || (biop1->genome == 1 && biop2->genome == 0))
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/*****************************************************************************
*
*   SubSourceNew()
*
*****************************************************************************/
NLM_EXTERN SubSourcePtr LIBCALL SubSourceNew (void)
{
    return (SubSourcePtr)MemNew(sizeof(SubSource));
}

/*****************************************************************************
*
*   SubSourceFree(ssp)
*       Frees an SubSource and associated data
*
*****************************************************************************/
NLM_EXTERN SubSourcePtr LIBCALL SubSourceFree (SubSourcePtr ssp)
{
    if (ssp == NULL)
        return (SubSourcePtr)NULL;

    MemFree(ssp->name);
    MemFree(ssp->attrib);
    return (SubSourcePtr)MemFree(ssp);
}

/*****************************************************************************
*
*   SubSourceAsnWrite(ssp, aip, atp)
*       atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SubSource ::=)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SubSourceAsnWrite (SubSourcePtr ssp, AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(orig, SUBSOURCE);   /* link local tree */
    if (atp == NULL)
        return FALSE;

    if (ssp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

    if (! AsnOpenStruct(aip, atp, (Pointer)ssp))
        goto erret;

    av.intvalue = (Int4)(ssp->subtype);
    if (AsnWrite(aip, SUBSOURCE_subtype, &av) <= 0) goto erret;

    if (ssp->name != NULL)
    {
        av.ptrvalue = ssp->name;
        if (! AsnWrite(aip, SUBSOURCE_name, &av)) goto erret;
    }

    if (ssp->attrib != NULL)
    {
        av.ptrvalue = ssp->attrib;
        if (! AsnWrite(aip, SUBSOURCE_attrib, &av)) goto erret;
    }

    if (! AsnCloseStruct(aip, atp, (Pointer)ssp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SubSourceAsnRead(aip, atp)
*       atp is the current type (if identifier of a parent struct)
*            assumption is readIdent has occurred
*       if atp == NULL, then assumes it stands alone and read ident
*            has not occurred.
*
*****************************************************************************/
NLM_EXTERN SubSourcePtr LIBCALL SubSourceAsnRead (AsnIoPtr aip, AsnTypePtr orig)
{
    DataVal av;
    AsnTypePtr atp, oldatp;
    SubSourcePtr ssp;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (SubSourcePtr)NULL;
    }

    if (aip == NULL)
        return (SubSourcePtr)NULL;

    if (orig == NULL)           /* SubSource ::= (self contained) */
        atp = AsnReadId(aip, amp, SUBSOURCE);
    else
        atp = AsnLinkType(orig, SUBSOURCE);    /* link in local tree */
    oldatp = atp;
    if (atp == NULL)
        return (SubSourcePtr)NULL;

    ssp = SubSourceNew();
    if (ssp == NULL)
        goto erret;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */

    while ((atp = AsnReadId(aip, amp, atp)) != oldatp)
    {
        if (atp == NULL)
            goto erret;
        if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
        if (atp == SUBSOURCE_subtype)
            ssp->subtype = (Uint1)(av.intvalue);
        else if (atp == SUBSOURCE_name)
            ssp->name = (CharPtr)(av.ptrvalue);
        else if (atp == SUBSOURCE_attrib)
            ssp->attrib = (CharPtr)(av.ptrvalue);
    }
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
ret:
    AsnUnlinkType(orig);       /* unlink local tree */
    return ssp;
erret:
    ssp = SubSourceFree(ssp);
    goto ret;
}

/*****************************************************************************
*
*   SubSourceSetFree (ssp)
*
*****************************************************************************/
NLM_EXTERN SubSourcePtr LIBCALL SubSourceSetFree (SubSourcePtr ssp)
{
    SubSourcePtr next;

    while (ssp != NULL)
    {
        next = ssp->next;
        SubSourceFree(ssp);
        ssp = next;
    }
    return ssp;
}

/*****************************************************************************
*
*   SubSourceSetAsnWrite(ssp, aip, set, element)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SubSourceSetAsnWrite (SubSourcePtr ssp, AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    AsnTypePtr atp;
    SubSourcePtr oldssp;
    Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

    if (aip == NULL)
        return FALSE;

    atp = AsnLinkType(element, SUBSOURCE);   /* link local tree */
    if (atp == NULL)
        return FALSE;
    oldssp = ssp;

    if (! AsnOpenStruct(aip, set, (Pointer)oldssp))
        goto erret;

    while (ssp != NULL)
    {
        if (! SubSourceAsnWrite(ssp, aip, atp))
            goto erret;
        ssp = ssp->next;
    }

    if (! AsnCloseStruct(aip, set, (Pointer)oldssp))
        goto erret;
    retval = TRUE;
erret:
    AsnUnlinkType(element);       /* unlink local tree */
    return retval;
}

/*****************************************************************************
*
*   SubSourceSetAsnRead(aip, set, element)
*
*****************************************************************************/
NLM_EXTERN SubSourcePtr LIBCALL SubSourceSetAsnRead (AsnIoPtr aip, AsnTypePtr set, AsnTypePtr element)
{
    DataVal av;
    AsnTypePtr atp;
    SubSourcePtr ssp, curr = NULL, first = NULL;

    if (aip == NULL)
        return (SubSourcePtr)NULL;

    AsnLinkType(element, SUBSOURCE);    /* link in local tree */
    atp = set;

    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;    /* read the start struct */
    while ((atp = AsnReadId(aip, amp, atp)) == element)
    {
        ssp = SubSourceAsnRead(aip, atp);
        if (ssp == NULL)
            goto erret;
        if (first == NULL)
            first = ssp;
        else
            curr->next = ssp;
        curr = ssp;
    }
    if (atp == NULL)
        goto erret;
    if (AsnReadVal(aip, atp, &av) <= 0) goto erret;   /* end struct */
    if (first == NULL)
        ErrPostEx(SEV_WARNING, 0, 0, "Empty SET OF SubSource. line %ld", (long) aip->linenumber);
ret:
    AsnUnlinkType(element);       /* unlink local tree */
    return first;
erret:
    first = SubSourceSetFree(first);
    goto ret;
}

/*****************************************************************************
*
*   SubSourceSetMatch(ssp1, ssp2)
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL SubSourceSetMatch (SubSourcePtr ssp1, SubSourcePtr ssp2)
{
  while (ssp1 != NULL && ssp2 != NULL)
  {
    if (ssp1->subtype != ssp2->subtype)
    {
      return FALSE;
    }
    else if (StringCmp (ssp1->name, ssp2->name) != 0)
    {
      return FALSE;
    }
    else if (StringCmp (ssp1->attrib, ssp2->attrib) != 0)
    {
      return FALSE;
    }
    
    ssp1 = ssp1->next;
    ssp2 = ssp2->next;
  }
 
  if (ssp1 == NULL && ssp2 == NULL)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/**************************************************
*
*    PCRPrimerNew()
*
**************************************************/
NLM_EXTERN 
PCRPrimerPtr LIBCALL
PCRPrimerNew(void)
{
   PCRPrimerPtr ptr = MemNew((size_t) sizeof(PCRPrimer));

   return ptr;

}

/**************************************************
*
*    PCRPrimerFree()
*
**************************************************/
NLM_EXTERN 
PCRPrimerPtr LIBCALL
PCRPrimerFree(PCRPrimerPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> seq);
   MemFree(ptr -> name);
   return MemFree(ptr);
}

/**************************************************
*
*    PCRPrimerAsnRead()
*
**************************************************/
NLM_EXTERN 
PCRPrimerPtr LIBCALL
PCRPrimerAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   PCRPrimerPtr ptr;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (PCRPrimerPtr)NULL;
    }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PCRPrimer ::= (self contained) */
      atp = AsnReadId(aip, amp, PCRPRIMER);
   } else {
      atp = AsnLinkType(orig, PCRPRIMER);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PCRPrimerNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PCRPRIMER_seq) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> seq = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PCRPRIMER_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PCRPrimerFree(ptr);
   goto ret;
}

/**************************************************
*
*    PCRPrimerAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PCRPrimerAsnWrite(PCRPrimerPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PCRPRIMER);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> seq != NULL) {
      av.ptrvalue = ptr -> seq;
      retval = AsnWrite(aip, PCRPRIMER_seq,  &av);
   }
   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, PCRPRIMER_name,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

/**************************************************
*
*    PCRPrimerSetFree()
*
**************************************************/
NLM_EXTERN 
PCRPrimerSetPtr LIBCALL
PCRPrimerSetFree(PCRPrimerSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) PCRPrimerFree);
   return NULL;
}

/**************************************************
*
*    PCRPrimerSetAsnRead()
*
**************************************************/
NLM_EXTERN 
PCRPrimerSetPtr LIBCALL
PCRPrimerSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PCRPrimerSetPtr ptr;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (PCRPrimerSetPtr)NULL;
    }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PCRPrimerSet ::= (self contained) */
      atp = AsnReadId(aip, amp, PCRPRIMERSET);
   } else {
      atp = AsnLinkType(orig, PCRPRIMERSET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) PCRPrimerAsnRead, (AsnOptFreeFunc) PCRPrimerFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PCRPrimerSetFree(ptr);
   goto ret;
}

/**************************************************
*
*    PCRPrimerSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PCRPrimerSetAsnWrite(PCRPrimerSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PCRPRIMERSET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) PCRPrimerAsnWrite, aip, atp, PCRPRIMERSET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

/**************************************************
*
*    PCRReactionNew()
*
**************************************************/
NLM_EXTERN 
PCRReactionPtr LIBCALL
PCRReactionNew(void)
{
   PCRReactionPtr ptr = MemNew((size_t) sizeof(PCRReaction));

   return ptr;

}

/**************************************************
*
*    PCRReactionFree()
*
**************************************************/
NLM_EXTERN 
PCRReactionPtr LIBCALL
PCRReactionFree(PCRReactionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   PCRPrimerSetFree(ptr -> forward);
   PCRPrimerSetFree(ptr -> reverse);
   return MemFree(ptr);
}

/**************************************************
*
*    PCRReactionAsnRead()
*
**************************************************/
NLM_EXTERN 
PCRReactionPtr LIBCALL
PCRReactionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   PCRReactionPtr ptr;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (PCRReactionPtr)NULL;
    }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PCRReaction ::= (self contained) */
      atp = AsnReadId(aip, amp, PCRREACTION);
   } else {
      atp = AsnLinkType(orig, PCRREACTION);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PCRReactionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PCRREACTION_forward) {
      ptr -> forward = PCRPrimerSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PCRREACTION_reverse) {
      ptr -> reverse = PCRPrimerSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PCRReactionFree(ptr);
   goto ret;
}

/**************************************************
*
*    PCRReactionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PCRReactionAsnWrite(PCRReactionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PCRREACTION);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> forward != NULL) {
      if ( ! PCRPrimerSetAsnWrite(ptr -> forward, aip, PCRREACTION_forward)) {
         goto erret;
      }
   }
   if (ptr -> reverse != NULL) {
      if ( ! PCRPrimerSetAsnWrite(ptr -> reverse, aip, PCRREACTION_reverse)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

/**************************************************
*
*    PCRReactionSetFree()
*
**************************************************/
NLM_EXTERN 
PCRReactionSetPtr LIBCALL
PCRReactionSetFree(PCRReactionSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) PCRReactionFree);
   return NULL;
}

/**************************************************
*
*    PCRReactionSetAsnRead()
*
**************************************************/
NLM_EXTERN 
PCRReactionSetPtr LIBCALL
PCRReactionSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PCRReactionSetPtr ptr;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return (PCRReactionSetPtr)NULL;
    }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PCRReactionSet ::= (self contained) */
      atp = AsnReadId(aip, amp, PCRREACTIONSET);
   } else {
      atp = AsnLinkType(orig, PCRREACTIONSET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) PCRReactionAsnRead, (AsnOptFreeFunc) PCRReactionFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PCRReactionSetFree(ptr);
   goto ret;
}

/**************************************************
*
*    PCRReactionSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PCRReactionSetAsnWrite(PCRReactionSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

    if (! loaded)
    {
        if (! SeqFeatAsnLoad())
            return FALSE;
    }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PCRREACTIONSET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) PCRReactionAsnWrite, aip, atp, PCRREACTIONSET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    CloneRefNew()
*
**************************************************/
NLM_EXTERN 
CloneRefPtr LIBCALL
CloneRefNew(void)
{
   CloneRefPtr ptr = MemNew((size_t) sizeof(CloneRef));

   ptr -> concordant = 0;
   ptr -> unique = 0;
   return ptr;

}


/**************************************************
*
*    CloneRefFree()
*
**************************************************/
NLM_EXTERN 
CloneRefPtr LIBCALL
CloneRefFree(CloneRefPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> name);
   MemFree(ptr -> library);
   CloneSeqSetFree(ptr -> clone_seq);
   return MemFree(ptr);
}


/**************************************************
*
*    CloneRefAsnRead()
*
**************************************************/
NLM_EXTERN 
CloneRefPtr LIBCALL
CloneRefAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   CloneRefPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CloneRef ::= (self contained) */
      atp = AsnReadId(aip, amp, CLONE_REF);
   } else {
      atp = AsnLinkType(orig, CLONE_REF);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CloneRefNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CLONE_REF_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_REF_library) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> library = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_REF_concordant) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> concordant = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_REF_unique) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> unique = av.boolvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_REF_placement_method) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> placement_method = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_REF_clone_seq) {
      ptr -> clone_seq = CloneSeqSetAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CloneRefFree(ptr);
   goto ret;
}



/**************************************************
*
*    CloneRefAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CloneRefAsnWrite(CloneRefPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CLONE_REF);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, CLONE_REF_name,  &av);
   }
   if (ptr -> library != NULL) {
      av.ptrvalue = ptr -> library;
      retval = AsnWrite(aip, CLONE_REF_library,  &av);
   }
   av.boolvalue = ptr -> concordant;
   retval = AsnWrite(aip, CLONE_REF_concordant,  &av);
   av.boolvalue = ptr -> unique;
   retval = AsnWrite(aip, CLONE_REF_unique,  &av);
   if (ptr -> placement_method || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> placement_method;
      retval = AsnWrite(aip, CLONE_REF_placement_method,  &av);
   }
   if (ptr -> clone_seq != NULL) {
      if ( ! CloneSeqSetAsnWrite(ptr -> clone_seq, aip, CLONE_REF_clone_seq)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CloneSeqSetFree()
*
**************************************************/
NLM_EXTERN 
CloneSeqSetPtr LIBCALL
CloneSeqSetFree(CloneSeqSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericUserSeqOfFree(ptr,  (AsnOptFreeFunc) CloneSeqFree);
   return NULL;
}


/**************************************************
*
*    CloneSeqSetAsnRead()
*
**************************************************/
NLM_EXTERN 
CloneSeqSetPtr LIBCALL
CloneSeqSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   CloneSeqSetPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CloneSeqSet ::= (self contained) */
      atp = AsnReadId(aip, amp, CLONE_SEQ_SET);
   } else {
      atp = AsnLinkType(orig, CLONE_SEQ_SET);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   func = NULL;

   ptr  = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) CloneSeqAsnRead, (AsnOptFreeFunc) CloneSeqFree);
   if (isError && ptr  == NULL) {
      goto erret;
   }



ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CloneSeqSetFree(ptr);
   goto ret;
}



/**************************************************
*
*    CloneSeqSetAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CloneSeqSetAsnWrite(CloneSeqSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CLONE_SEQ_SET);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   retval = AsnGenericUserSeqOfAsnWrite(ptr , (AsnWriteFunc) CloneSeqAsnWrite, aip, atp, CLONE_SEQ_SET_E);
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    CloneSeqNew()
*
**************************************************/
NLM_EXTERN 
CloneSeqPtr LIBCALL
CloneSeqNew(void)
{
   CloneSeqPtr ptr = MemNew((size_t) sizeof(CloneSeq));

   return ptr;

}


/**************************************************
*
*    CloneSeqFree()
*
**************************************************/
NLM_EXTERN 
CloneSeqPtr LIBCALL
CloneSeqFree(CloneSeqPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SeqLocFree(ptr -> location);
   SeqLocFree(ptr -> seq);
   DbtagFree(ptr -> align_id);
   return MemFree(ptr);
}


/**************************************************
*
*    CloneSeqAsnRead()
*
**************************************************/
NLM_EXTERN 
CloneSeqPtr LIBCALL
CloneSeqAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   CloneSeqPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* CloneSeq ::= (self contained) */
      atp = AsnReadId(aip, amp, CLONE_SEQ);
   } else {
      atp = AsnLinkType(orig, CLONE_SEQ);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = CloneSeqNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == CLONE_SEQ_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_SEQ_confidence) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> confidence = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_SEQ_location) {
      ptr -> location = SeqLocAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_SEQ_seq) {
      ptr -> seq = SeqLocAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == CLONE_SEQ_align_id) {
      ptr -> align_id = DbtagAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = CloneSeqFree(ptr);
   goto ret;
}



/**************************************************
*
*    CloneSeqAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
CloneSeqAsnWrite(CloneSeqPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, CLONE_SEQ);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, CLONE_SEQ_type,  &av);
   if (ptr -> confidence || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> confidence;
      retval = AsnWrite(aip, CLONE_SEQ_confidence,  &av);
   }
   if (ptr -> location != NULL) {
      if ( ! SeqLocAsnWrite(ptr -> location, aip, CLONE_SEQ_location)) {
         goto erret;
      }
   }
   if (ptr -> seq != NULL) {
      if ( ! SeqLocAsnWrite(ptr -> seq, aip, CLONE_SEQ_seq)) {
         goto erret;
      }
   }
   if (ptr -> align_id != NULL) {
      if ( ! DbtagAsnWrite(ptr -> align_id, aip, CLONE_SEQ_align_id)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    SomaticOriginNew()
*
**************************************************/
NLM_EXTERN 
SomaticOriginPtr LIBCALL
SomaticOriginNew(void)
{
   SomaticOriginPtr ptr = MemNew((size_t) sizeof(SomaticOrigin));

   return ptr;

}


/**************************************************
*
*    SomaticOriginConditionNew()
*
**************************************************/
NLM_EXTERN 
SomaticOriginConditionPtr LIBCALL
SomaticOriginConditionNew(void)
{
   SomaticOriginConditionPtr ptr = MemNew((size_t) sizeof(SomaticOriginCondition));

   return ptr;

}


/**************************************************
*
*    VariationRefNew()
*
**************************************************/
NLM_EXTERN 
VariationRefPtr LIBCALL
VariationRefNew(void)
{
   VariationRefPtr ptr = MemNew((size_t) sizeof(VariationRef));

   return ptr;

}


/**************************************************
*
*    Consequence_loss_of_heterozygosityNew()
*
**************************************************/
NLM_EXTERN 
Consequence_loss_of_heterozygosityPtr LIBCALL
Consequence_loss_of_heterozygosityNew(void)
{
   Consequence_loss_of_heterozygosityPtr ptr = MemNew((size_t) sizeof(Consequence_loss_of_heterozygosity));

   return ptr;

}


/**************************************************
*
*    SomaticOriginFree()
*
**************************************************/
NLM_EXTERN 
SomaticOriginPtr LIBCALL
SomaticOriginFree(SomaticOriginPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   SubSourceFree(ptr -> source);
   SomaticOriginConditionFree(ptr -> condition);
   return MemFree(ptr);
}


/**************************************************
*
*    SomaticOriginConditionFree()
*
**************************************************/
NLM_EXTERN 
SomaticOriginConditionPtr LIBCALL
SomaticOriginConditionFree(SomaticOriginConditionPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> description);
   AsnGenericValNodeSetFree(ptr -> object_id, (AsnOptFreeFunc) DbtagFree);
   return MemFree(ptr);
}


/**************************************************
*
*    VariationRefFree()
*
**************************************************/
NLM_EXTERN 
VariationRefPtr LIBCALL
VariationRefFree(VariationRefPtr ptr)
{
   if(ptr == NULL) {
      return NULL;
   }
   DbtagFree(ptr -> id);
   DbtagFree(ptr -> parent_id);
   ObjectIdFree(ptr -> sample_id);
   AsnGenericValNodeSetFree(ptr -> other_ids, (AsnOptFreeFunc) DbtagFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> other_ids, (AsnOptFreeFunc) DbtagFree);
   */
   MemFree(ptr -> name);
   /*
   AsnGenericValNodeSetFree(ptr -> synonyms, (AsnOptFreeFunc) MemFree);
   */
   AsnGenericBaseSeqOfFree(ptr -> synonyms ,ASNCODE_PTRVAL_SLOT);
   MemFree(ptr -> description);
   AsnGenericValNodeSetFree(ptr -> phenotype, (AsnOptFreeFunc) PhenotypeFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> phenotype, (AsnOptFreeFunc) PhenotypeFree);
   */
   AsnGenericBaseSeqOfFree(ptr -> method ,ASNCODE_INTVAL_SLOT);
   AsnGenericValNodeSetFree(ptr -> population_data, (AsnOptFreeFunc) PopulationDataFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> population_data, (AsnOptFreeFunc) PopulationDataFree);
   */
   VariantPropertiesFree(ptr -> variant_prop);
   AsnGenericValNodeSetFree(ptr -> clinical_test, (AsnOptFreeFunc) DbtagFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> clinical_test, (AsnOptFreeFunc) DbtagFree);
   */
   PubFree(ptr -> pub);
   VarRefDataFree(ptr -> data);
   ConsequenceFree(ptr -> consequence);
   SeqLocFree(ptr -> location);
   AsnGenericValNodeSetFree(ptr -> ext_locs, (AsnOptFreeFunc) ExtLocFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> ext_locs, (AsnOptFreeFunc) ExtLocFree);
   */
   UserObjectFree(ptr -> ext);
   AsnGenericValNodeSetFree(ptr -> somatic_origin, (AsnOptFreeFunc) SomaticOriginFree);
   return MemFree(ptr);
}


/**************************************************
*
*    ConsequenceFree()
*
**************************************************/
NLM_EXTERN 
ConsequencePtr LIBCALL
ConsequenceFree(ValNodePtr anp)
{

   if (anp == NULL) {
      return NULL;
   }

   AsnGenericChoiceSeqOfFree(anp, (AsnOptFreeFunc) Consequence_elementFree);    
   return NULL;
}


/**************************************************
*
*    SomaticOriginAsnRead()
*
**************************************************/
NLM_EXTERN 
SomaticOriginPtr LIBCALL
SomaticOriginAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   SomaticOriginPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SomaticOrigin ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_somatic_origin);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_somatic_origin);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SomaticOriginNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_REF_somatic_origin_E_source) {
      ptr -> source = SubSourceAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_somatic_origin_E_condition) {
      ptr -> condition = SomaticOriginConditionAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SomaticOriginFree(ptr);
   goto ret;
}



/**************************************************
*
*    SomaticOriginConditionAsnRead()
*
**************************************************/
NLM_EXTERN 
SomaticOriginConditionPtr LIBCALL
SomaticOriginConditionAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   SomaticOriginConditionPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* SomaticOriginCondition ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_somatic_origin_E_condition);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_somatic_origin_E_condition);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = SomaticOriginConditionNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_REF_somatic_origin_E_condition_description) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> description = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_somatic_origin_E_condition_object_id) {
      ptr -> object_id = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      if (isError && ptr -> object_id == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = SomaticOriginConditionFree(ptr);
   goto ret;
}



/**************************************************
*
*    VariationRefAsnRead()
*
**************************************************/
NLM_EXTERN 
VariationRefPtr LIBCALL
VariationRefAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   VariationRefPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* VariationRef ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = VariationRefNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_REF_id) {
      ptr -> id = DbtagAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_parent_id) {
      ptr -> parent_id = DbtagAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_sample_id) {
      ptr -> sample_id = ObjectIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_other_ids) {
      ptr -> other_ids = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      /*
      ptr -> other_ids = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      */
      if (isError && ptr -> other_ids == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_synonyms) {
      /*
      ptr -> synonyms = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) NULL, (AsnOptFreeFunc) MemFree);
      */
      ptr -> synonyms = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_PTRVAL_SLOT, &isError);
      if (isError && ptr -> synonyms == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_description) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> description = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_phenotype) {
      ptr -> phenotype = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) PhenotypeAsnRead, (AsnOptFreeFunc) PhenotypeFree);
      /*
      ptr -> phenotype = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) PhenotypeAsnRead, (AsnOptFreeFunc) PhenotypeFree);
      */
      if (isError && ptr -> phenotype == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_method) {
      ptr -> method = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> method == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_population_data) {
      ptr -> population_data = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) PopulationDataAsnRead, (AsnOptFreeFunc) PopulationDataFree);
      /*
      ptr -> population_data = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) PopulationDataAsnRead, (AsnOptFreeFunc) PopulationDataFree);
      */
      if (isError && ptr -> population_data == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_variant_prop) {
      ptr -> variant_prop = VariantPropertiesAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_validated) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> validated = av.boolvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_clinical_test) {
      ptr -> clinical_test = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      /*
      ptr -> clinical_test = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      */
      if (isError && ptr -> clinical_test == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_allele_origin) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_origin = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_allele_state) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_state = av.intvalue;
      ptr -> OBbits__ |= 1<<2;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_allele_frequency) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_frequency = av.realvalue;
      ptr -> OBbits__ |= 1<<3;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_is_ancestral_allele) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_ancestral_allele = av.boolvalue;
      ptr -> OBbits__ |= 1<<4;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_pub) {
      ptr -> pub = PubAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_data) {
      ptr -> data = VarRefDataAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_consequence) {
      ptr -> consequence =  ConsequenceAsnRead(aip, atp);
      if (isError && ptr -> consequence == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_location) {
      ptr -> location = SeqLocAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_ext_locs) {
      ptr -> ext_locs = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) ExtLocAsnRead, (AsnOptFreeFunc) ExtLocFree);
      /*
      ptr -> ext_locs = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ExtLocAsnRead, (AsnOptFreeFunc) ExtLocFree);
      */
      if (isError && ptr -> ext_locs == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_ext) {
      ptr -> ext = UserObjectAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_somatic_origin) {
      ptr -> somatic_origin = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) SomaticOriginAsnRead, (AsnOptFreeFunc) SomaticOriginFree);
      if (isError && ptr -> somatic_origin == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = VariationRefFree(ptr);
   goto ret;
}



/**************************************************
*
*    ConsequenceAsnRead()
*
**************************************************/
NLM_EXTERN 
ConsequencePtr LIBCALL
ConsequenceAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{


   AsnTypePtr atp;
   ValNodePtr anp;
   Boolean isError = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Consequence_element ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_consequence);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_consequence);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp =
   AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError,
   (AsnReadFunc) Consequence_elementAsnRead, (AsnOptFreeFunc) Consequence_elementFree);
   if (isError) 
   goto erret;


ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    SomaticOriginAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SomaticOriginAsnWrite(SomaticOriginPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_REF_somatic_origin);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> source != NULL) {
      if ( ! SubSourceAsnWrite(ptr -> source, aip, VARIATION_REF_somatic_origin_E_source)) {
         goto erret;
      }
   }
   if (ptr -> condition != NULL) {
      if ( ! SomaticOriginConditionAsnWrite(ptr -> condition, aip, VARIATION_REF_somatic_origin_E_condition)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    SomaticOriginConditionAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
SomaticOriginConditionAsnWrite(SomaticOriginConditionPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_REF_somatic_origin_E_condition);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> description != NULL) {
      av.ptrvalue = ptr -> description;
      retval = AsnWrite(aip, VARIATION_REF_somatic_origin_E_condition_description,  &av);
   }
   if (ptr -> object_id != NULL) {
     AsnGenericValNodeSetAsnWrite(ptr -> object_id, (AsnWriteFunc) DbtagAsnWrite, aip, VARIATION_REF_somatic_origin_E_condition_object_id, VARIATION_REF_somatic_origin_E_condition_object_id_E);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    VariationRefAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
VariationRefAsnWrite(VariationRefPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_REF);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> id != NULL) {
      if ( ! DbtagAsnWrite(ptr -> id, aip, VARIATION_REF_id)) {
         goto erret;
      }
   }
   if (ptr -> parent_id != NULL) {
      if ( ! DbtagAsnWrite(ptr -> parent_id, aip, VARIATION_REF_parent_id)) {
         goto erret;
      }
   }
   if (ptr -> sample_id != NULL) {
      if ( ! ObjectIdAsnWrite(ptr -> sample_id, aip, VARIATION_REF_sample_id)) {
         goto erret;
      }
   }
   if (ptr -> other_ids != NULL) {
     AsnGenericValNodeSetAsnWrite (ptr -> other_ids, (AsnWriteFunc) DbtagAsnWrite, aip, VARIATION_REF_other_ids, VARIATION_REF_other_ids_E);
     /*
     AsnGenericUserSeqOfAsnWrite(ptr -> other_ids, (AsnWriteFunc) DbtagAsnWrite, aip, VARIATION_REF_other_ids, VARIATION_REF_other_ids_E);
     */
   }
   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, VARIATION_REF_name,  &av);
   }
   if (ptr -> synonyms != NULL) {
     retval = AsnGenericBaseSeqOfAsnWrite(ptr -> synonyms ,ASNCODE_PTRVAL_SLOT, aip, VARIATION_REF_synonyms, VARIATION_REF_synonyms_E);
   }
   if (ptr -> description != NULL) {
      av.ptrvalue = ptr -> description;
      retval = AsnWrite(aip, VARIATION_REF_description,  &av);
   }
   if (ptr -> phenotype != NULL) {
     AsnGenericValNodeSetAsnWrite (ptr -> phenotype, (AsnWriteFunc) PhenotypeAsnWrite, aip, VARIATION_REF_phenotype, VARIATION_REF_phenotype_E);
     /*
     AsnGenericUserSeqOfAsnWrite(ptr -> phenotype, (AsnWriteFunc) PhenotypeAsnWrite, aip, VARIATION_REF_phenotype, VARIATION_REF_phenotype_E);
     */
   }
   if (ptr -> method != NULL) {
     retval = AsnGenericBaseSeqOfAsnWrite(ptr -> method ,ASNCODE_INTVAL_SLOT, aip, VARIATION_REF_method, VARIATION_REF_method_E);
   }
   if (ptr -> population_data != NULL) {
     AsnGenericValNodeSetAsnWrite (ptr -> population_data, (AsnWriteFunc) PopulationDataAsnWrite, aip, VARIATION_REF_population_data, VARIATION_REF_population_data_E);
     /*
     AsnGenericUserSeqOfAsnWrite(ptr -> population_data, (AsnWriteFunc) PopulationDataAsnWrite, aip, VARIATION_REF_population_data, VARIATION_REF_population_data_E);
     */
   }
   if (ptr -> variant_prop != NULL) {
      if ( ! VariantPropertiesAsnWrite(ptr -> variant_prop, aip, VARIATION_REF_variant_prop)) {
         goto erret;
      }
   }
   if (ptr -> validated || (ptr -> OBbits__ & (1<<0) )){   av.boolvalue = ptr -> validated;
      retval = AsnWrite(aip, VARIATION_REF_validated,  &av);
   }
   if (ptr -> clinical_test != NULL) {
     AsnGenericValNodeSetAsnWrite (ptr -> clinical_test, (AsnWriteFunc) DbtagAsnWrite, aip, VARIATION_REF_clinical_test, VARIATION_REF_clinical_test_E);
     /*
     AsnGenericUserSeqOfAsnWrite(ptr -> clinical_test, (AsnWriteFunc) DbtagAsnWrite, aip, VARIATION_REF_clinical_test, VARIATION_REF_clinical_test_E);
     */
   }
   if (ptr -> allele_origin || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> allele_origin;
      retval = AsnWrite(aip, VARIATION_REF_allele_origin,  &av);
   }
   if (ptr -> allele_state || (ptr -> OBbits__ & (1<<2) )){   av.intvalue = ptr -> allele_state;
      retval = AsnWrite(aip, VARIATION_REF_allele_state,  &av);
   }
   if (ptr -> allele_frequency || (ptr -> OBbits__ & (1<<3) )){   av.realvalue = ptr -> allele_frequency;
      retval = AsnWrite(aip, VARIATION_REF_allele_frequency,  &av);
   }
   if (ptr -> is_ancestral_allele || (ptr -> OBbits__ & (1<<4) )){   av.boolvalue = ptr -> is_ancestral_allele;
      retval = AsnWrite(aip, VARIATION_REF_is_ancestral_allele,  &av);
   }
   if (ptr -> pub != NULL) {
      if ( ! PubAsnWrite(ptr -> pub, aip, VARIATION_REF_pub)) {
         goto erret;
      }
   }
   if (ptr -> data != NULL) {
      if ( ! VarRefDataAsnWrite(ptr -> data, aip, VARIATION_REF_data)) {
         goto erret;
      }
   }
   if (ptr -> consequence != NULL) {
     retval = ConsequenceAsnWrite(ptr -> consequence, aip,VARIATION_REF_consequence);
   }
   if (ptr -> location != NULL) {
      if ( ! SeqLocAsnWrite(ptr -> location, aip, VARIATION_REF_location)) {
         goto erret;
      }
   }
   if (ptr -> ext_locs != NULL) {
     AsnGenericValNodeSetAsnWrite (ptr -> ext_locs, (AsnWriteFunc) ExtLocAsnWrite, aip, VARIATION_REF_ext_locs, VARIATION_REF_ext_locs_E);
     /*
     AsnGenericUserSeqOfAsnWrite(ptr -> ext_locs, (AsnWriteFunc) ExtLocAsnWrite, aip, VARIATION_REF_ext_locs, VARIATION_REF_ext_locs_E);
     */
   }
   if (ptr -> ext != NULL) {
      if ( ! UserObjectAsnWrite(ptr -> ext, aip, VARIATION_REF_ext)) {
         goto erret;
      }
   }
   if (ptr -> somatic_origin != NULL) {
     AsnGenericValNodeSetAsnWrite (ptr -> somatic_origin, (AsnWriteFunc) SomaticOriginAsnWrite, aip, VARIATION_REF_somatic_origin, VARIATION_REF_somatic_origin_E);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ConsequenceAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ConsequenceAsnWrite(ValNodePtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, VARIATION_REF_consequence);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   retval = AsnGenericChoiceSeqOfAsnWrite(anp, 
   (AsnWriteFunc) Consequence_elementAsnWrite, aip, atp, VARIATION_REF_consequence_E);
erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    Consequence_elementAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
Consequence_elementAsnWrite(Consequence_elementPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, VARIATION_REF_consequence_E);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case Consequence_unknown:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_unknown, &av);
      break;
   case Consequence_splicing:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_splicing, &av);
      break;
   case Consequence_note:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_note, &av);
      break;
   case Consequence_variation:
      writetype = VARIATION_REF_consequence_E_variation;
      func = (AsnWriteFunc) VariationRefAsnWrite;
      break;
   case Consequence_Consequence_Frameshift:
      writetype = VARIATION_REF_consequence_E_frameshift;
      func = (AsnWriteFunc) Consequence_frameshiftAsnWrite;
      break;
   case Consequence_Consequence_LossOfHeterozygosity:
      writetype = VARIATION_REF_consequence_E_loss_of_heterozygosity;
      func = (AsnWriteFunc) Consequence_loss_of_heterozygosityAsnWrite;
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    Consequence_loss_of_heterozygosityAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
Consequence_loss_of_heterozygosityAsnWrite(Consequence_loss_of_heterozygosityPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_REF_consequence_E_loss_of_heterozygosity);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> reference != NULL) {
      av.ptrvalue = ptr -> reference;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_loss_of_heterozygosity_reference,  &av);
   }
   if (ptr -> test != NULL) {
      av.ptrvalue = ptr -> test;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_loss_of_heterozygosity_test,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Consequence_frameshiftAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
Consequence_frameshiftAsnWrite(Consequence_frameshiftPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_REF_consequence_E_frameshift);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> phase || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> phase;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_frameshift_phase,  &av);
   }
   if (ptr -> x_length || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> x_length;
      retval = AsnWrite(aip, VARIATION_REF_consequence_E_frameshift_x_length,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Data_setAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
VarRefDataSetAsnWrite(VarRefDataSetPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_REF_data_set);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, VARIATION_REF_data_set_type,  &av);
   AsnGenericValNodeSetAsnWrite(ptr -> variations, (AsnWriteFunc) VariationRefAsnWrite, aip, VARIATION_REF_data_set_variations, VARIATION_REF_data_set_variations_E);
   /*
   AsnGenericUserSeqOfAsnWrite(ptr -> variations, (AsnWriteFunc) VariationRefAsnWrite, aip, VARIATION_REF_data_set_variations, VARIATION_REF_data_set_variations_E);
   */
   if (ptr -> name != NULL) {
      av.ptrvalue = ptr -> name;
      retval = AsnWrite(aip, VARIATION_REF_data_set_name,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}


/**************************************************
*
*    VarRefDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
VarRefDataAsnWrite( VarRefDataPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, VARIATION_REF_data);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case VarRefData_unknown:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, VARIATION_REF_data_unknown, &av);
      break;
   case VarRefData_note:
      av.ptrvalue = anp->data.ptrvalue;
      retval = AsnWrite(aip, VARIATION_REF_data_note, &av);
      break;
   case VarRefData_uniparental_disomy:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, VARIATION_REF_data_uniparental_disomy, &av);
      break;
   case VarRefData_instance:
      writetype = VARIATION_REF_data_instance;
      func = (AsnWriteFunc) VariationInstAsnWrite;
      break;
   case VarRefData_set:
      writetype = VARIATION_REF_data_set;
      func = (AsnWriteFunc) VarRefDataSetAsnWrite;
      break;
   case VarRefData_complex:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, VARIATION_REF_data_complex, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Consequence_elementAsnRead()
*
**************************************************/
NLM_EXTERN 
Consequence_elementPtr LIBCALL
Consequence_elementAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Consequence_element ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_consequence_E);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_consequence_E);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == VARIATION_REF_consequence_E_unknown) {
      choice = Consequence_unknown;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == VARIATION_REF_consequence_E_splicing) {
      choice = Consequence_splicing;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == VARIATION_REF_consequence_E_note) {
      choice = Consequence_note;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == VARIATION_REF_consequence_E_variation) {
      choice = Consequence_variation;
      func = (AsnReadFunc) VariationRefAsnRead;
   }
   else if (atp == VARIATION_REF_consequence_E_frameshift) {
      choice = Consequence_Consequence_Frameshift;
      func = (AsnReadFunc) Consequence_frameshiftAsnRead;
   }
   else if (atp == VARIATION_REF_consequence_E_loss_of_heterozygosity) {
      choice = Consequence_Consequence_LossOfHeterozygosity;
      func = (AsnReadFunc) Consequence_loss_of_heterozygosityAsnRead;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    Consequence_loss_of_heterozygosityAsnRead()
*
**************************************************/
NLM_EXTERN 
Consequence_loss_of_heterozygosityPtr LIBCALL
Consequence_loss_of_heterozygosityAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   Consequence_loss_of_heterozygosityPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Consequence_loss_of_heterozygosity ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_consequence_E_loss_of_heterozygosity);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_consequence_E_loss_of_heterozygosity);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = Consequence_loss_of_heterozygosityNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_REF_consequence_E_loss_of_heterozygosity_reference) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> reference = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_consequence_E_loss_of_heterozygosity_test) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> test = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = Consequence_loss_of_heterozygosityFree(ptr);
   goto ret;
}



/**************************************************
*
*    Consequence_frameshiftAsnRead()
*
**************************************************/
NLM_EXTERN 
Consequence_frameshiftPtr LIBCALL
Consequence_frameshiftAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   Consequence_frameshiftPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Consequence_frameshift ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_consequence_E_frameshift);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_consequence_E_frameshift);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = Consequence_frameshiftNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_REF_consequence_E_frameshift_phase) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> phase = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_consequence_E_frameshift_x_length) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> x_length = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = Consequence_frameshiftFree(ptr);
   goto ret;
}



/**************************************************
*
*    Data_setFree()
*
**************************************************/
NLM_EXTERN 
VarRefDataSetPtr LIBCALL
VarRefDataSetFree(VarRefDataSetPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericValNodeSetFree (ptr -> variations, (AsnOptFreeFunc) VariationRefFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> variations, (AsnOptFreeFunc) VariationRefFree);
   */
   MemFree(ptr -> name);
   return MemFree(ptr);
}


/**************************************************
*
*    VarRefDataFree()
*
**************************************************/
NLM_EXTERN 
VarRefDataPtr LIBCALL
VarRefDataFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case VarRefData_note:
      MemFree(anp -> data.ptrvalue);
      break;
   case VarRefData_instance:
      VariationInstFree(anp -> data.ptrvalue);
      break;
   case VarRefData_set:
      VarRefDataSetFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    Data_setNew()
*
**************************************************/
NLM_EXTERN 
VarRefDataSetPtr LIBCALL
VarRefDataSetNew(void)
{
   VarRefDataSetPtr ptr = MemNew((size_t) sizeof(VarRefData));

   return ptr;

}


/**************************************************
*
*    Data_setAsnRead()
*
**************************************************/
NLM_EXTERN 
VarRefDataSetPtr LIBCALL
VarRefDataSetAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   VarRefDataSetPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Data_set ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_data_set);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_data_set);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = VarRefDataSetNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_REF_data_set_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_data_set_variations) {
      ptr -> variations = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) VariationRefAsnRead, (AsnOptFreeFunc) VariationRefFree);
      /*
      ptr -> variations = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) VariationRefAsnRead, (AsnOptFreeFunc) VariationRefFree);
      */
      if (isError && ptr -> variations == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_REF_data_set_name) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> name = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = VarRefDataSetFree(ptr);
   goto ret;
}


/**************************************************
*
*    VarRefDataAsnRead()
*
**************************************************/
NLM_EXTERN 
VarRefDataPtr LIBCALL
VarRefDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* data ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_REF_data);
   } else {
      atp = AsnLinkType(orig, VARIATION_REF_data);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == VARIATION_REF_data_unknown) {
      choice = VarRefData_unknown;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == VARIATION_REF_data_note) {
      choice = VarRefData_note;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.ptrvalue = av.ptrvalue;
   }
   else if (atp == VARIATION_REF_data_uniparental_disomy) {
      choice = VarRefData_uniparental_disomy;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   else if (atp == VARIATION_REF_data_instance) {
      choice = VarRefData_instance;
      func = (AsnReadFunc) VariationInstAsnRead;
   }
   else if (atp == VARIATION_REF_data_set) {
      choice = VarRefData_set;
      func = (AsnReadFunc) VarRefDataSetAsnRead;
   }
   else if (atp == VARIATION_REF_data_complex) {
      choice = VarRefData_complex;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}



/**************************************************
*
*    Consequence_elementFree()
*
**************************************************/
NLM_EXTERN 
Consequence_elementPtr LIBCALL
Consequence_elementFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case Consequence_note:
      MemFree(anp -> data.ptrvalue);
      break;
   case Consequence_variation:
      VariationRefFree(anp -> data.ptrvalue);
      break;
   case Consequence_Consequence_Frameshift:
      Consequence_frameshiftFree(anp -> data.ptrvalue);
      break;
   case Consequence_Consequence_LossOfHeterozygosity:
      Consequence_loss_of_heterozygosityFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    Consequence_loss_of_heterozygosityFree()
*
**************************************************/
NLM_EXTERN 
Consequence_loss_of_heterozygosityPtr LIBCALL
Consequence_loss_of_heterozygosityFree(Consequence_loss_of_heterozygosityPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> reference);
   MemFree(ptr -> test);
   return MemFree(ptr);
}


/**************************************************
*
*    Consequence_frameshiftFree()
*
**************************************************/
NLM_EXTERN 
Consequence_frameshiftPtr LIBCALL
Consequence_frameshiftFree(Consequence_frameshiftPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   return MemFree(ptr);
}


/**************************************************
*
*    Consequence_frameshiftNew()
*
**************************************************/
NLM_EXTERN 
Consequence_frameshiftPtr LIBCALL
Consequence_frameshiftNew(void)
{
   Consequence_frameshiftPtr ptr = MemNew((size_t) sizeof(Consequence_frameshift));

   return ptr;

}


/**************************************************
*
*    VariationInstNew()
*
**************************************************/
NLM_EXTERN 
VariationInstPtr LIBCALL
VariationInstNew(void)
{
   VariationInstPtr ptr = MemNew((size_t) sizeof(VariationInst));

   return ptr;

}


/**************************************************
*
*    VariationInstFree()
*
**************************************************/
NLM_EXTERN 
VariationInstPtr LIBCALL
VariationInstFree(VariationInstPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericValNodeSetFree (ptr -> delta, (AsnOptFreeFunc) DeltaItemFree);
   /*
   AsnGenericUserSeqOfFree(ptr -> delta, (AsnOptFreeFunc) DeltaItemFree);
   */
   return MemFree(ptr);
}


/**************************************************
*
*    VariationInstAsnRead()
*
**************************************************/
NLM_EXTERN 
VariationInstPtr LIBCALL
VariationInstAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   VariationInstPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* VariationInst ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIATION_INST);
   } else {
      atp = AsnLinkType(orig, VARIATION_INST);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = VariationInstNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIATION_INST_type) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> type = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_INST_delta) {
      ptr -> delta = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) DeltaItemAsnRead, (AsnOptFreeFunc) DeltaItemFree);
      /*
      ptr -> delta = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DeltaItemAsnRead, (AsnOptFreeFunc) DeltaItemFree);
      */
      if (isError && ptr -> delta == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIATION_INST_observation) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> observation = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = VariationInstFree(ptr);
   goto ret;
}



/**************************************************
*
*    VariationInstAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
VariationInstAsnWrite(VariationInstPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIATION_INST);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> type;
   retval = AsnWrite(aip, VARIATION_INST_type,  &av);
   AsnGenericValNodeSetAsnWrite (ptr -> delta, (AsnWriteFunc) DeltaItemAsnWrite, aip, VARIATION_INST_delta, VARIATION_INST_delta_E);
   /*
   AsnGenericUserSeqOfAsnWrite(ptr -> delta, (AsnWriteFunc) DeltaItemAsnWrite, aip, VARIATION_INST_delta, VARIATION_INST_delta_E);
   */
   if (ptr -> observation || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> observation;
      retval = AsnWrite(aip, VARIATION_INST_observation,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    VariantPropertiesNew()
*
**************************************************/
NLM_EXTERN 
VariantPropertiesPtr LIBCALL
VariantPropertiesNew(void)
{
   VariantPropertiesPtr ptr = MemNew((size_t) sizeof(VariantProperties));

   return ptr;

}


/**************************************************
/**************************************************
*
*    VariantPropertiesFree()
*
**************************************************/
NLM_EXTERN 
VariantPropertiesPtr LIBCALL
VariantPropertiesFree(VariantPropertiesPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   AsnGenericBaseSeqOfFree(ptr -> project_data, ASNCODE_INTVAL_SLOT);
   return MemFree(ptr);
}


/**************************************************
*
*    VariantPropertiesAsnRead()
*
**************************************************/
NLM_EXTERN 
VariantPropertiesPtr LIBCALL
VariantPropertiesAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   VariantPropertiesPtr ptr;
   Boolean isError = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* VariantProperties ::= (self contained) */
      atp = AsnReadId(aip, amp, VARIANTPROPERTIES);
   } else {
      atp = AsnLinkType(orig, VARIANTPROPERTIES);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = VariantPropertiesNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == VARIANTPROPERTIES_version) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> version = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_resource_link) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> resource_link = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_gene_location) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> gene_location = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_effect) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> effect = av.intvalue;
      ptr -> OBbits__ |= 1<<2;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_mapping) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> mapping = av.intvalue;
      ptr -> OBbits__ |= 1<<3;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_map_weight) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> map_weight = av.intvalue;
      ptr -> OBbits__ |= 1<<4;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_frequency_based_validation) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> freq_based_validation = av.intvalue;
      ptr -> OBbits__ |= 1<<5;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_genotype) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> genotype = av.intvalue;
      ptr -> OBbits__ |= 1<<6;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_project_data) {
      ptr -> project_data = AsnGenericBaseSeqOfAsnRead(aip, amp, atp, ASNCODE_INTVAL_SLOT, &isError);
      if (isError && ptr -> project_data == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_quality_check) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> quality_check = av.intvalue;
      ptr -> OBbits__ |= 1<<7;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_confidence) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> confidence = av.intvalue;
      ptr -> OBbits__ |= 1<<8;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_other_validation) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> other_validation = av.boolvalue;
      ptr -> OBbits__ |= 1<<9;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_allele_origin) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_origin = av.intvalue;
      ptr -> OBbits__ |= 1<<10;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_allele_state) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_state = av.intvalue;
      ptr -> OBbits__ |= 1<<11;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_allele_frequency) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_frequency = av.realvalue;
      ptr -> OBbits__ |= 1<<12;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == VARIANTPROPERTIES_is_ancestral_allele) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> is_ancestral_allele = av.boolvalue;
      ptr -> OBbits__ |= 1<<13;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = VariantPropertiesFree(ptr);
   goto ret;
}



/**************************************************
*
*    VariantPropertiesAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
VariantPropertiesAsnWrite(VariantPropertiesPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, VARIANTPROPERTIES);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   av.intvalue = ptr -> version;
   retval = AsnWrite(aip, VARIANTPROPERTIES_version,  &av);
   if (ptr -> resource_link || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> resource_link;
      retval = AsnWrite(aip, VARIANTPROPERTIES_resource_link,  &av);
   }
   if (ptr -> gene_location || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> gene_location;
      retval = AsnWrite(aip, VARIANTPROPERTIES_gene_location,  &av);
   }
   if (ptr -> effect || (ptr -> OBbits__ & (1<<2) )){   av.intvalue = ptr -> effect;
      retval = AsnWrite(aip, VARIANTPROPERTIES_effect,  &av);
   }
   if (ptr -> mapping || (ptr -> OBbits__ & (1<<3) )){   av.intvalue = ptr -> mapping;
      retval = AsnWrite(aip, VARIANTPROPERTIES_mapping,  &av);
   }
   if (ptr -> map_weight || (ptr -> OBbits__ & (1<<4) )){   av.intvalue = ptr -> map_weight;
      retval = AsnWrite(aip, VARIANTPROPERTIES_map_weight,  &av);
   }
   if (ptr -> freq_based_validation || (ptr -> OBbits__ & (1<<5) )){   av.intvalue = ptr -> freq_based_validation;
      retval = AsnWrite(aip, VARIANTPROPERTIES_frequency_based_validation,  &av);
   }
   if (ptr -> genotype || (ptr -> OBbits__ & (1<<6) )){   av.intvalue = ptr -> genotype;
      retval = AsnWrite(aip, VARIANTPROPERTIES_genotype,  &av);
   }
   if (ptr -> project_data != NULL) {
      retval = AsnGenericBaseSeqOfAsnWrite(ptr -> project_data ,ASNCODE_INTVAL_SLOT, aip, VARIANTPROPERTIES_project_data, VARIANTPROPERTIES_project_data_E);
   }
   if (ptr -> quality_check || (ptr -> OBbits__ & (1<<7) )){   av.intvalue = ptr -> quality_check;
      retval = AsnWrite(aip, VARIANTPROPERTIES_quality_check,  &av);
   }
   if (ptr -> confidence || (ptr -> OBbits__ & (1<<8) )){   av.intvalue = ptr -> confidence;
      retval = AsnWrite(aip, VARIANTPROPERTIES_confidence,  &av);
   }
   if (ptr -> other_validation || (ptr -> OBbits__ & (1<<9) )){   av.boolvalue = ptr -> other_validation;
      retval = AsnWrite(aip, VARIANTPROPERTIES_other_validation,  &av);
   }
   if (ptr -> allele_origin || (ptr -> OBbits__ & (1<<10) )){   av.intvalue = ptr -> allele_origin;
      retval = AsnWrite(aip, VARIANTPROPERTIES_allele_origin,  &av);
   }
   if (ptr -> allele_state || (ptr -> OBbits__ & (1<<11) )){   av.intvalue = ptr -> allele_state;
      retval = AsnWrite(aip, VARIANTPROPERTIES_allele_state,  &av);
   }
   if (ptr -> allele_frequency || (ptr -> OBbits__ & (1<<12) )){   av.realvalue = ptr -> allele_frequency;
      retval = AsnWrite(aip, VARIANTPROPERTIES_allele_frequency,  &av);
   }
   if (ptr -> is_ancestral_allele || (ptr -> OBbits__ & (1<<13) )){   av.boolvalue = ptr -> is_ancestral_allele;
      retval = AsnWrite(aip, VARIANTPROPERTIES_is_ancestral_allele,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    PhenotypeNew()
*
**************************************************/
NLM_EXTERN 
PhenotypePtr LIBCALL
PhenotypeNew(void)
{
   PhenotypePtr ptr = MemNew((size_t) sizeof(Phenotype));

   return ptr;

}


/**************************************************
*
*    PhenotypeFree()
*
**************************************************/
NLM_EXTERN 
PhenotypePtr LIBCALL
PhenotypeFree(PhenotypePtr ptr)
{
   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> source);
   MemFree(ptr -> term);
   /*
   AsnGenericUserSeqOfFree(ptr -> xref, (AsnOptFreeFunc) DbtagFree);
   */
   AsnGenericValNodeSetFree (ptr->xref, (AsnOptFreeFunc) DbtagFree);
   return MemFree(ptr);
}


/**************************************************
*
*    PhenotypeAsnRead()
*
**************************************************/
NLM_EXTERN 
PhenotypePtr LIBCALL
PhenotypeAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   PhenotypePtr ptr;
   Boolean isError = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Phenotype ::= (self contained) */
      atp = AsnReadId(aip, amp, PHENOTYPE);
   } else {
      atp = AsnLinkType(orig, PHENOTYPE);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PhenotypeNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == PHENOTYPE_source) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> source = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PHENOTYPE_term) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> term = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PHENOTYPE_xref) {
      ptr->xref = AsnGenericValNodeSetAsnRead (aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      /*
      ptr -> xref = AsnGenericUserSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) DbtagAsnRead, (AsnOptFreeFunc) DbtagFree);
      */
      if (isError && ptr -> xref == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == PHENOTYPE_clinical_significance) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> clinical_significance = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PhenotypeFree(ptr);
   goto ret;
}



/**************************************************
*
*    PhenotypeAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PhenotypeAsnWrite(PhenotypePtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, PHENOTYPE);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> source != NULL) {
      av.ptrvalue = ptr -> source;
      retval = AsnWrite(aip, PHENOTYPE_source,  &av);
   }
   if (ptr -> term != NULL) {
      av.ptrvalue = ptr -> term;
      retval = AsnWrite(aip, PHENOTYPE_term,  &av);
   }
   AsnGenericValNodeSetAsnWrite (ptr->xref, (AsnWriteFunc) DbtagAsnWrite, aip, PHENOTYPE_xref, PHENOTYPE_xref_E);
   /*
   AsnGenericUserSeqOfAsnWrite(ptr -> xref, (AsnWriteFunc) DbtagAsnWrite, aip, PHENOTYPE_xref, PHENOTYPE_xref_E);
   */
   if (ptr -> clinical_significance || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> clinical_significance;
      retval = AsnWrite(aip, PHENOTYPE_clinical_significance,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    PopulationDataNew()
*
**************************************************/
NLM_EXTERN 
PopulationDataPtr LIBCALL
PopulationDataNew(void)
{
   PopulationDataPtr ptr = MemNew((size_t) sizeof(PopulationData));

   return ptr;

}


/**************************************************
*
*    PopulationDataFree()
*
**************************************************/
NLM_EXTERN 
PopulationDataPtr LIBCALL
PopulationDataFree(PopulationDataPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   MemFree(ptr -> population);
   AsnGenericValNodeSetFree (ptr -> sample_ids, (AsnOptFreeFunc) ObjectIdFree);
   /*
   AsnGenericChoiceSeqOfFree(ptr -> sample_ids, (AsnOptFreeFunc) ObjectIdFree);
   */
   return MemFree(ptr);
}


/**************************************************
*
*    PopulationDataAsnRead()
*
**************************************************/
NLM_EXTERN 
PopulationDataPtr LIBCALL
PopulationDataAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean isError = FALSE;
   AsnReadFunc func;
   PopulationDataPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* PopulationData ::= (self contained) */
      atp = AsnReadId(aip, amp, POPULATION_DATA);
   } else {
      atp = AsnLinkType(orig, POPULATION_DATA);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = PopulationDataNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == POPULATION_DATA_population) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> population = av.ptrvalue;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == POPULATION_DATA_genotype_frequency) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> genotype_frequency = av.realvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == POPULATION_DATA_chromosomes_tested) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> chromosomes_tested = av.intvalue;
      ptr -> OBbits__ |= 1<<1;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == POPULATION_DATA_sample_ids) {
      ptr -> sample_ids = AsnGenericValNodeSetAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ObjectIdAsnRead, (AsnOptFreeFunc) ObjectIdFree);
      /*
      ptr -> sample_ids = AsnGenericChoiceSeqOfAsnRead(aip, amp, atp, &isError, (AsnReadFunc) ObjectIdAsnRead, (AsnOptFreeFunc) ObjectIdFree);
      */
      if (isError && ptr -> sample_ids == NULL) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == POPULATION_DATA_allele_frequency) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> allele_frequency = av.realvalue;
      ptr -> OBbits__ |= 1<<2;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == POPULATION_DATA_flags) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> flags = av.intvalue;
      ptr -> OBbits__ |= 1<<3;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = PopulationDataFree(ptr);
   goto ret;
}



/**************************************************
*
*    PopulationDataAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
PopulationDataAsnWrite(PopulationDataPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, POPULATION_DATA);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> population != NULL) {
      av.ptrvalue = ptr -> population;
      retval = AsnWrite(aip, POPULATION_DATA_population,  &av);
   }
   if (ptr -> genotype_frequency || (ptr -> OBbits__ & (1<<0) )){   av.realvalue = ptr -> genotype_frequency;
      retval = AsnWrite(aip, POPULATION_DATA_genotype_frequency,  &av);
   }
   if (ptr -> chromosomes_tested || (ptr -> OBbits__ & (1<<1) )){   av.intvalue = ptr -> chromosomes_tested;
      retval = AsnWrite(aip, POPULATION_DATA_chromosomes_tested,  &av);
   }
   AsnGenericValNodeSetAsnWrite(ptr -> sample_ids, (AsnWriteFunc) ObjectIdAsnWrite, aip, POPULATION_DATA_sample_ids, POPULATION_DATA_sample_ids_E);
   /*
   AsnGenericChoiceSeqOfAsnWrite(ptr -> sample_ids, (AsnWriteFunc) ObjectIdAsnWrite, aip, POPULATION_DATA_sample_ids, POPULATION_DATA_sample_ids_E);
   */
   if (ptr -> allele_frequency || (ptr -> OBbits__ & (1<<2) )){   av.realvalue = ptr -> allele_frequency;
      retval = AsnWrite(aip, POPULATION_DATA_allele_frequency,  &av);
   }
   if (ptr -> flags || (ptr -> OBbits__ & (1<<3) )){   av.intvalue = ptr -> flags;
      retval = AsnWrite(aip, POPULATION_DATA_flags,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    ExtLocNew()
*
**************************************************/
NLM_EXTERN 
ExtLocPtr LIBCALL
ExtLocNew(void)
{
   ExtLocPtr ptr = MemNew((size_t) sizeof(ExtLoc));

   return ptr;

}


/**************************************************
*
*    ExtLocFree()
*
**************************************************/
NLM_EXTERN 
ExtLocPtr LIBCALL
ExtLocFree(ExtLocPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   ObjectIdFree(ptr -> id);
   SeqLocFree(ptr -> location);
   return MemFree(ptr);
}


/**************************************************
*
*    ExtLocAsnRead()
*
**************************************************/
NLM_EXTERN 
ExtLocPtr LIBCALL
ExtLocAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   ExtLocPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* ExtLoc ::= (self contained) */
      atp = AsnReadId(aip, amp, EXT_LOC);
   } else {
      atp = AsnLinkType(orig, EXT_LOC);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = ExtLocNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == EXT_LOC_id) {
      ptr -> id = ObjectIdAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == EXT_LOC_location) {
      ptr -> location = SeqLocAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = ExtLocFree(ptr);
   goto ret;
}



/**************************************************
*
*    ExtLocAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
ExtLocAsnWrite(ExtLocPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, EXT_LOC);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> id != NULL) {
      if ( ! ObjectIdAsnWrite(ptr -> id, aip, EXT_LOC_id)) {
         goto erret;
      }
   }
   if (ptr -> location != NULL) {
      if ( ! SeqLocAsnWrite(ptr -> location, aip, EXT_LOC_location)) {
         goto erret;
      }
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    DeltaItemNew()
*
**************************************************/
NLM_EXTERN 
DeltaItemPtr LIBCALL
DeltaItemNew(void)
{
   DeltaItemPtr ptr = MemNew((size_t) sizeof(DeltaItem));

   return ptr;

}


/**************************************************
*
*    DeltaItemFree()
*
**************************************************/
NLM_EXTERN 
DeltaItemPtr LIBCALL
DeltaItemFree(DeltaItemPtr ptr)
{

   if(ptr == NULL) {
      return NULL;
   }
   Seq_seqFree(ptr -> Seq_seq);
   IntFuzzFree(ptr -> multiplier_fuzz);
   return MemFree(ptr);
}


/**************************************************
*
*    Seq_seqFree()
*
**************************************************/
NLM_EXTERN 
Seq_seqPtr LIBCALL
Seq_seqFree(ValNodePtr anp)
{
   Pointer pnt;

   if (anp == NULL) {
      return NULL;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   default:
      break;
   case Seq_seq_literal:
      SeqLitFree(anp -> data.ptrvalue);
      break;
   case Seq_seq_loc:
      SeqLocFree(anp -> data.ptrvalue);
      break;
   }
   return MemFree(anp);
}


/**************************************************
*
*    DeltaItemAsnRead()
*
**************************************************/
NLM_EXTERN 
DeltaItemPtr LIBCALL
DeltaItemAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   AsnReadFunc func;
   DeltaItemPtr ptr;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* DeltaItem ::= (self contained) */
      atp = AsnReadId(aip, amp, DELTA_ITEM);
   } else {
      atp = AsnLinkType(orig, DELTA_ITEM);
   }
   /* link in local tree */
   if (atp == NULL) {
      return NULL;
   }

   ptr = DeltaItemNew();
   if (ptr == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the start struct */
      goto erret;
   }

   atp = AsnReadId(aip,amp, atp);
   func = NULL;

   if (atp == DELTA_ITEM_seq) {
      ptr -> Seq_seq = Seq_seqAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == DELTA_ITEM_multiplier) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> multiplier = av.intvalue;
      ptr -> OBbits__ |= 1<<0;
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == DELTA_ITEM_multiplier_fuzz) {
      ptr -> multiplier_fuzz = IntFuzzAsnRead(aip, atp);
      if (aip -> io_failure) {
         goto erret;
      }
      atp = AsnReadId(aip,amp, atp);
   }
   if (atp == DELTA_ITEM_action) {
      if ( AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      ptr -> action = av.intvalue;
      atp = AsnReadId(aip,amp, atp);
   }

   if (AsnReadVal(aip, atp, &av) <= 0) {
      goto erret;
   }
   /* end struct */

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return ptr;

erret:
   aip -> io_failure = TRUE;
   ptr = DeltaItemFree(ptr);
   goto ret;
}



/**************************************************
*
*    Seq_seqAsnRead()
*
**************************************************/
NLM_EXTERN 
Seq_seqPtr LIBCALL
Seq_seqAsnRead(AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   ValNodePtr anp;
   Uint1 choice;
   Boolean nullIsError = FALSE;
   AsnReadFunc func;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return NULL;
      }
   }

   if (aip == NULL) {
      return NULL;
   }

   if (orig == NULL) {         /* Seq_seq ::= (self contained) */
      atp = AsnReadId(aip, amp, DELTA_ITEM_seq);
   } else {
      atp = AsnLinkType(orig, DELTA_ITEM_seq);    /* link in local tree */
   }
   if (atp == NULL) {
      return NULL;
   }

   anp = ValNodeNew(NULL);
   if (anp == NULL) {
      goto erret;
   }
   if (AsnReadVal(aip, atp, &av) <= 0) { /* read the CHOICE or OpenStruct value (nothing) */
      goto erret;
   }

   func = NULL;

   atp = AsnReadId(aip, amp, atp);  /* find the choice */
   if (atp == NULL) {
      goto erret;
   }
   if (atp == DELTA_ITEM_seq_literal) {
      choice = Seq_seq_literal;
      func = (AsnReadFunc) SeqLitAsnRead;
   }
   else if (atp == DELTA_ITEM_seq_loc) {
      choice = Seq_seq_loc;
      func = (AsnReadFunc) SeqLocAsnRead;
   }
   else if (atp == DELTA_ITEM_seq_this) {
      choice = Seq_seq_this__;
      if (AsnReadVal(aip, atp, &av) <= 0) {
         goto erret;
      }
      anp->data.boolvalue = av.boolvalue;
   }
   anp->choice = choice;
   if (func != NULL)
   {
      anp->data.ptrvalue = (* func)(aip, atp);
      if (aip -> io_failure) goto erret;

      if (nullIsError && anp->data.ptrvalue == NULL) {
         goto erret;
      }
   }

ret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return anp;

erret:
   anp = MemFree(anp);
   aip -> io_failure = TRUE;
   goto ret;
}


/**************************************************
*
*    DeltaItemAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
DeltaItemAsnWrite(DeltaItemPtr ptr, AsnIoPtr aip, AsnTypePtr orig)
{
   DataVal av;
   AsnTypePtr atp;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad()) {
         return FALSE;
      }
   }

   if (aip == NULL) {
      return FALSE;
   }

   atp = AsnLinkType(orig, DELTA_ITEM);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (ptr == NULL) { AsnNullValueMsg(aip, atp); goto erret; }
   if (! AsnOpenStruct(aip, atp, (Pointer) ptr)) {
      goto erret;
   }

   if (ptr -> Seq_seq != NULL) {
      if ( ! Seq_seqAsnWrite(ptr -> Seq_seq, aip, DELTA_ITEM_seq)) {
         goto erret;
      }
   }
   if (ptr -> multiplier || (ptr -> OBbits__ & (1<<0) )){   av.intvalue = ptr -> multiplier;
      retval = AsnWrite(aip, DELTA_ITEM_multiplier,  &av);
   }
   if (ptr -> multiplier_fuzz != NULL) {
      if ( ! IntFuzzAsnWrite(ptr -> multiplier_fuzz, aip, DELTA_ITEM_multiplier_fuzz)) {
         goto erret;
      }
   }
   if (ptr -> action != 0) {
      av.intvalue = ptr -> action;
      retval = AsnWrite(aip, DELTA_ITEM_action,  &av);
   }
   if (! AsnCloseStruct(aip, atp, (Pointer)ptr)) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}



/**************************************************
*
*    Seq_seqAsnWrite()
*
**************************************************/
NLM_EXTERN Boolean LIBCALL 
Seq_seqAsnWrite(Seq_seqPtr anp, AsnIoPtr aip, AsnTypePtr orig)

{
   DataVal av;
   AsnTypePtr atp, writetype = NULL;
   Pointer pnt;
   AsnWriteFunc func = NULL;
   Boolean retval = FALSE;

   if (! loaded)
   {
      if (! SeqFeatAsnLoad())
      return FALSE;
   }

   if (aip == NULL)
   return FALSE;

   atp = AsnLinkType(orig, DELTA_ITEM_seq);   /* link local tree */
   if (atp == NULL) {
      return FALSE;
   }

   if (anp == NULL) { AsnNullValueMsg(aip, atp); goto erret; }

   av.ptrvalue = (Pointer)anp;
   if (! AsnWriteChoice(aip, atp, (Int2)anp->choice, &av)) {
      goto erret;
   }

   pnt = anp->data.ptrvalue;
   switch (anp->choice)
   {
   case Seq_seq_literal:
      writetype = DELTA_ITEM_seq_literal;
      func = (AsnWriteFunc) SeqLitAsnWrite;
      break;
   case Seq_seq_loc:
      writetype = DELTA_ITEM_seq_loc;
      func = (AsnWriteFunc) SeqLocAsnWrite;
      break;
   case Seq_seq_this__:
      av.boolvalue = anp->data.boolvalue;
      retval = AsnWrite(aip, DELTA_ITEM_seq_this, &av);
      break;
   }
   if (writetype != NULL) {
      retval = (* func)(pnt, aip, writetype);   /* write it out */
   }
   if (!retval) {
      goto erret;
   }
   retval = TRUE;

erret:
   AsnUnlinkType(orig);       /* unlink local tree */
   return retval;
}

