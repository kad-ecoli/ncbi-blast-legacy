/*  objgen.h
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
* File Name:  objgen.h
*
* Author:  James Ostell
*   
* Version Creation Date: 1/1/91
*
* $Revision: 6.15 $
*
* File Description:  Object manager interface for module NCBI-General
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBI_General_
#define _NCBI_General_

#ifndef _ASNTOOL_
#include <asn.h>
#endif

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* gather/objmgr fields to be added to bioseq, features, descriptors, etc. */

#define EXTRA_OBJMGR_FIELDS \
  Uint2             entityID;   \
  Uint2             itemtype;   \
  Uint1             subtype;    \
  Uint1             deleteme;   \
  Uint2             parenttype; \
  Uint4             itemID;     \
  Pointer           parentptr;  \
  Pointer PNTR      prevlink;   \
  Pointer           scratch;

/* structure containing gather/objmgr fields to add as a block to other structures */

typedef struct gatherindex {
  EXTRA_OBJMGR_FIELDS
} GatherIndex, PNTR GatherIndexPtr;

/* extended valnode for linking seqdesc, perhaps seqid and seqloc */

typedef struct objvalnode {
  ValNode      vn;
  GatherIndex  idx;
} ObjValNode, PNTR ObjValNodePtr;

/* ValNode equivalent functions that allocate ObjValNode size */

NLM_EXTERN ValNodePtr LIBCALL SeqDescrNew (ValNodePtr vnp);
NLM_EXTERN ValNodePtr LIBCALL SeqDescrAdd (ValNodePtr PNTR head);
NLM_EXTERN ValNodePtr LIBCALL SeqDescrAddPointer (ValNodePtr PNTR head, Int2 choice, VoidPtr value);

/*****************************************************************************
*
*   loader
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL GeneralAsnLoad PROTO((void));

/*****************************************************************************
*
*   internal structures for NCBI-General objects
*
*****************************************************************************/

/*****************************************************************************
*
*   Date, Date-std share the same structure
*      any data[2] or data[3] values = 0 means not set or not present
*   data [0] - CHOICE of date ,0=str, 1=std
*        [1] - year (- 1900)
*        [2] - month (1-12)  optional
*   	 [3] - day (1-31)	 optional
*        [4] - hour (0-23) optional 255=not set
*        [5] - minute (0-59) optional 255=not set
*        [6] - second (0-59) optional 255=not set
*        [7] - not currently used
*
*****************************************************************************/

#define NOT_SET 255

typedef struct date {
	Uint1 data[8];      /* see box above */
	CharPtr str;		/* str or season or NULL */
} NCBI_Date, PNTR NCBI_DatePtr;
#define DatePtr NCBI_DatePtr


NLM_EXTERN NCBI_DatePtr LIBCALL DateNew PROTO((void));
NLM_EXTERN NCBI_DatePtr LIBCALL DateClean PROTO((NCBI_DatePtr dp));
NLM_EXTERN NCBI_DatePtr LIBCALL DateFree PROTO((NCBI_DatePtr dp));
NLM_EXTERN Boolean      LIBCALL DateWrite PROTO((NCBI_DatePtr dp, Int2 year, Int2 month, Int2 day, CharPtr season));
NLM_EXTERN Boolean      LIBCALL DateRead PROTO((NCBI_DatePtr dp, Int2Ptr year, Int2Ptr month, Int2Ptr day, CharPtr season));
NLM_EXTERN Boolean      LIBCALL DatePrint PROTO((NCBI_DatePtr dp, CharPtr buf));
NLM_EXTERN NCBI_DatePtr LIBCALL DateCurr PROTO((void)); /* time fields not set */
NLM_EXTERN NCBI_DatePtr LIBCALL DateTimeCurr PROTO((void)); /* fills time fields too */
NLM_EXTERN NCBI_DatePtr LIBCALL DateDup PROTO((NCBI_DatePtr dp));
NLM_EXTERN Boolean      LIBCALL DateAsnWrite PROTO((NCBI_DatePtr dp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN NCBI_DatePtr LIBCALL DateAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Int2         LIBCALL DateMatch PROTO((DatePtr a, DatePtr b, Boolean all));
NLM_EXTERN DatePtr      LIBCALL DateParse (CharPtr str);

/*****************************************************************************
*
*   DateCheck (dp)
*   	Checks the date and month values in a date structure
*       returns:
*      -4 = NULL pointer passed in
*      -3 = string date, can't be checked
*      -2 = month not set (but otherwise ok)
*      -1 = day not set	 (but otherwise ok)
*       0 = date ok, month,day,year all set
*       1 = day invalid
*       2 = month invalid
*   	3 = year not set (required for date)
*
*****************************************************************************/
NLM_EXTERN Int2         LIBCALL DateCheck PROTO((DatePtr dp));

/*****************************************************************************
*
*   Object-id stuff
*
*****************************************************************************/
typedef struct objid {
	Int4 id;
	CharPtr str;
} ObjectId, PNTR ObjectIdPtr;

NLM_EXTERN ObjectIdPtr LIBCALL ObjectIdNew PROTO((void));
NLM_EXTERN ObjectIdPtr LIBCALL ObjectIdFree PROTO(( ObjectIdPtr oid));
NLM_EXTERN ObjectIdPtr LIBCALL ObjectIdAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean     LIBCALL ObjectIdAsnWrite PROTO((ObjectIdPtr oid, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean     LIBCALL ObjectIdMatch PROTO((ObjectIdPtr a, ObjectIdPtr b));
NLM_EXTERN Boolean     LIBCALL ObjectIdMatchEx PROTO((ObjectIdPtr a, ObjectIdPtr b, Boolean case_sensitive));
NLM_EXTERN ObjectIdPtr LIBCALL ObjectIdDup PROTO((ObjectIdPtr oldid));

/*****************************************************************************
*
*   DBtag stuff
*
*****************************************************************************/
typedef struct dbtag {
	CharPtr db;
	ObjectIdPtr tag;
} Dbtag, PNTR DbtagPtr;

NLM_EXTERN DbtagPtr LIBCALL DbtagNew PROTO((void));
NLM_EXTERN DbtagPtr LIBCALL DbtagFree PROTO(( DbtagPtr dbt));
NLM_EXTERN DbtagPtr LIBCALL DbtagAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean  LIBCALL DbtagAsnWrite PROTO((DbtagPtr dbt, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean LIBCALL DbtagMatchEx PROTO((DbtagPtr a, DbtagPtr b, Boolean case_sensitive));
NLM_EXTERN Boolean  LIBCALL DbtagMatch PROTO((DbtagPtr a, DbtagPtr b));
NLM_EXTERN DbtagPtr LIBCALL DbtagDup PROTO((DbtagPtr oldtag));
NLM_EXTERN Int2     LIBCALL DbtagLabel PROTO((DbtagPtr dbt, CharPtr buf, Int2 buflen));

/*****************************************************************************
*
*   Name-std
*   names[0] = last
*        [1] = first
*        [2] = middle
*        [3] = full
*        [4] = initials
*        [5] = suffix
*        [6] = title
*
*****************************************************************************/
typedef struct namestd {
	CharPtr names[7];
} NameStd, PNTR NameStdPtr;

NLM_EXTERN NameStdPtr LIBCALL NameStdNew PROTO((void));
NLM_EXTERN NameStdPtr LIBCALL NameStdFree PROTO(( NameStdPtr nsp));
NLM_EXTERN NameStdPtr LIBCALL NameStdAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean    LIBCALL NameStdAsnWrite PROTO((NameStdPtr nsp, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean LIBCALL NameStdMatch (NameStdPtr nsp1, NameStdPtr nsp2);

/*****************************************************************************
*
*   Person-id
*     choice = 0 = not set
*              1 = dbtag
*              2 = name
*              3 = ml
*              4 = str
*              5 = consortium
*
*****************************************************************************/
typedef struct personid {
	Uint1 choice;         /* which CHOICE, see above */
	Pointer data;         /* points to appropriate data structure */
} PersonId, PNTR PersonIdPtr;

NLM_EXTERN PersonIdPtr	LIBCALL PersonIdNew PROTO((void));
NLM_EXTERN PersonIdPtr	LIBCALL PersonIdFree PROTO(( PersonIdPtr pid));
NLM_EXTERN PersonIdPtr	LIBCALL PersonIdAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean     	LIBCALL PersonIdAsnWrite PROTO((PersonIdPtr pid, AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean LIBCALL PersonIdMatch (PersonIdPtr pip1, PersonIdPtr pip2);

/*****************************************************************************
*
*   PersonIdLabel(pid, buf, buflen, format)
*   	Makes a short label, lastname then initials if it can
*   	format = PIDLABEL_GENBANK   last,initials
*   		   PIDLABEL_EMBL      last initials
*
*       Modeled from GBGetAuthNames in asn2ff
*   	returns number of bytes in buf
*   	buf MUST be at least (buflen + 1) long
*
*****************************************************************************/
NLM_EXTERN Int2 LIBCALL PersonIdLabel PROTO((PersonIdPtr pid, CharPtr buf, Int2 buflen, Int2 format));
#define PIDLABEL_GENBANK 1
#define PIDLABEL_EMBL 2

/*****************************************************************************
*
*   Int-fuzz
*
*****************************************************************************/
typedef struct intfuzz {
	Uint1 choice;       /* 1=p-m, 2=range, 3=pct, 4=lim 5=alt */
	Int4 a, b;          /* a=p-m,max,pct,orlim, b=min */
	Int4Ptr alt;		/* alternate positions, a=num alts, b=array size */
} IntFuzz, PNTR IntFuzzPtr;

NLM_EXTERN IntFuzzPtr	LIBCALL IntFuzzNew PROTO((void));
NLM_EXTERN IntFuzzPtr	LIBCALL IntFuzzFree PROTO(( IntFuzzPtr ifp));
NLM_EXTERN IntFuzzPtr	LIBCALL IntFuzzAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean     	LIBCALL IntFuzzAsnWrite PROTO((IntFuzzPtr ifp, AsnIoPtr aip, AsnTypePtr atp));

/*****************************************************************************
*
*   User-field
*      data is an DataVal where:
*    choice    asn1              data. =
        1 = str VisibleString ,  ptrvalue = CharPtr
        2 = int INTEGER ,        intvalue
        3 = real REAL ,          realvalue
        4 = bool BOOLEAN ,       boolvalue
        5 = os OCTET STRING ,    ptrvalue = ByteStorePtr
        6 = object User-object ,   ptrvalue = UserObjectPtr
        7 = strs SEQUENCE OF VisibleString ,  ptrvalue = CharPtr PNTR
        8 = ints SEQUENCE OF INTEGER ,        ptrvalue = Int4Ptr
        9 = reals SEQUENCE OF REAL ,          ptrvalue = FloatHiPtr
        10 = oss SEQUENCE OF OCTET STRING ,   ptrvalue = ByteStorePtr PNTR
        11 = fields SEQUENCE OF User-field ,  ptrvalue = UserFieldPtr
        12 = objects SEQUENCE OF User-object } }  ptrvalue = UserObjectPtr

*   User-object
*
*****************************************************************************/
typedef struct userfield {
    ObjectIdPtr label;
    Int4 num;
    Uint1 choice;
    DataVal data;
    struct userfield PNTR next;
} UserField, PNTR UserFieldPtr;

NLM_EXTERN UserFieldPtr LIBCALL UserFieldNew PROTO((void));
NLM_EXTERN UserFieldPtr LIBCALL UserFieldFree PROTO(( UserFieldPtr ufp));
NLM_EXTERN UserFieldPtr LIBCALL UserFieldAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean		 LIBCALL UserFieldAsnWrite PROTO((UserFieldPtr ufp, AsnIoPtr aip, AsnTypePtr atp));

typedef struct userobj {
    CharPtr _class;
    ObjectIdPtr type;
    UserFieldPtr data;
    struct userobj PNTR next;   /* for SEQUENCE OF User-object */
} UserObject, PNTR UserObjectPtr;

NLM_EXTERN UserObjectPtr LIBCALL UserObjectNew PROTO((void));
NLM_EXTERN UserObjectPtr LIBCALL UserObjectFree PROTO(( UserObjectPtr uop));
NLM_EXTERN UserObjectPtr LIBCALL UserObjectAsnRead PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean       LIBCALL UserObjectAsnWrite PROTO((UserObjectPtr uop, AsnIoPtr aip, AsnTypePtr atp));

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif

