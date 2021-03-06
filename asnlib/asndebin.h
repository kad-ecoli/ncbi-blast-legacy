/* asndebin.h
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
* File Name:  asndebin.h
*
* Author:  James Ostell
*
* Version Creation Date: 3/4/91
*
* $Revision: 6.1 $
*
* File Description:
*   Typedefs and prototypes for internal routines of asndebin.c
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* 3/4/91   Kans        Stricter typecasting for GNU C and C++
* 3/4/91   Kans        AsnDeBinReadBoolean returns Boolean
*
* $Log: asndebin.h,v $
* Revision 6.1  2000/12/12 15:56:11  ostell
* added support BigInt
*
* Revision 6.0  1997/08/25 18:09:48  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/12/03 21:43:48  vakatov
* Adopted for 32-bit MS-Windows DLLs
*
 * Revision 5.0  1996/05/28  14:00:29  ostell
 * Set to revision 5.0
 *
 * Revision 4.0  1995/07/26  13:47:38  ostell
 * force revision to 4.0
 *
 * Revision 2.1  1995/05/15  18:38:28  ostell
 * added Log line
 *
*
* ==========================================================================
*/

#ifndef _ASNDEBIN_
#define _ASNDEBIN_

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


	/* states in AsnDeBinTWord and AsnDeBinWord */
#define INDEFINATE -1              /* reading indefinate length item */

/*****************************************************************************
*
*   prototypes
*   	AsnBinReadVal
*       AsnBinReadId  prototyped in asngen.h
*
*****************************************************************************/
NLM_EXTERN Boolean AsnDeBinDecr PROTO((Int4 used, AsnIoPtr aip));
NLM_EXTERN Int4 AsnDeBinScanTag PROTO((AsnIoPtr aip));
NLM_EXTERN AsnTypePtr AsnDeBinFindType PROTO((AsnIoPtr aip, AsnModulePtr amp));
NLM_EXTERN AsnTypePtr AsnDeBinFindElement PROTO((AsnIoPtr aip, AsnTypePtr atp));

NLM_EXTERN Int4 AsnDeBinReadInteger PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Int8 AsnDeBinReadBigInt PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN Boolean AsnDeBinReadBoolean PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN CharPtr AsnDeBinReadString PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN void AsnDeBinSkipString PROTO((AsnIoPtr aip));
NLM_EXTERN ByteStorePtr AsnDeBinReadOctets PROTO((AsnIoPtr aip));
NLM_EXTERN void AsnDeBinSkipOctets PROTO((AsnIoPtr aip));
NLM_EXTERN FloatHi AsnDeBinReadReal PROTO((AsnIoPtr aip));
NLM_EXTERN void AsnDeBinReadNull PROTO((AsnIoPtr aip, AsnTypePtr atp));
NLM_EXTERN void AsnDeBinSkipStruct PROTO((AsnIoPtr aip));


#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
