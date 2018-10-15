/*  ncbibs.h
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
* File Name:  ncbibs.h
*
* Author:  Jim Ostell
*
* Version Creation Date:  1/1/91
*
* $Revision: 6.5 $
*
* File Description:
*   ByteStore typedefs, prototypes, and defines
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#ifndef _NCBIBS_
#define _NCBIBS_

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif

typedef struct bsunit {             /* for building multiline strings */
	Nlm_Handle str;            /* the string piece */
	Nlm_Int2 len_avail,
		 len;
	struct bsunit PNTR next; }       /* the next one */
Nlm_BSUnit, PNTR Nlm_BSUnitPtr;

typedef struct bytestore {
	Nlm_Int4 seekptr,       /* current position */
		totlen,             /* total stored data length in bytes */
		chain_offset;       /* offset in ByteStore of first byte in curchain */
	Nlm_BSUnitPtr chain,       /* chain of elements */
		curchain;           /* the BSUnit containing seekptr */
} Nlm_ByteStore, PNTR Nlm_ByteStorePtr;

NLM_EXTERN Nlm_VoidPtr LIBCALL Nlm_BSMerge PROTO((Nlm_ByteStorePtr ssp, Nlm_VoidPtr dest));
NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSNew PROTO((Nlm_Int4 len));
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSSeek PROTO((Nlm_ByteStorePtr bsp, Nlm_Int4 offset, Nlm_Int2 origin));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSTell PROTO((Nlm_ByteStorePtr bsp));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSDelete PROTO((Nlm_ByteStorePtr bsp, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSWrite PROTO((Nlm_ByteStorePtr bsp, Nlm_VoidPtr ptr, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSInsert PROTO((Nlm_ByteStorePtr bsp, Nlm_VoidPtr ptr, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSInsertFromBS PROTO((Nlm_ByteStorePtr bsp, Nlm_ByteStorePtr bsp2, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSRead PROTO((Nlm_ByteStorePtr bsp, Nlm_VoidPtr ptr, Nlm_Int4 len));
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSGetByte PROTO((Nlm_ByteStorePtr bsp));
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSPutByte PROTO((Nlm_ByteStorePtr bsp, Nlm_Int2 value));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSLen PROTO((Nlm_ByteStorePtr ssp));
NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSFree PROTO((Nlm_ByteStorePtr ssp));
NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSDup PROTO((Nlm_ByteStorePtr source));
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_BSEqual PROTO((Nlm_ByteStorePtr bs1, Nlm_ByteStorePtr bs2));
/*****************************************************************************
*
*   Int4 BSAdd(bsp, len, use_min_size)
*   	adds len bytes BEFORE current bsp->seekptr
*       bsp->seekptr returned pointing at first added byte
*   	returns bytes added
*       if (use_min_size) then does not add anything smaller than MIN_BSALLOC
*
*       BSAdd does not change the length of the byte store.  When inserting in
*       the middle of the chain, it is up to the calling function to adjust the
*       length after writing the bytes.
*
*****************************************************************************/
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSAdd PROTO((Nlm_ByteStorePtr bsp, Nlm_Int4 len, Nlm_Boolean use_min_size));

/****************************************************************************
*
*   Integer storage utilities
*      These assume integers are store in BIG_ENDIAN order in the ByteStore
*      They read and write Uint2 or Uint4 in the NATIVE endian order
*      All work with UNSIGNED 2 or 4 byte integers
*         (you should cast to make them signed)
*      These are just helper functions. They do no internal consistency
*         checking.
*      They are primarily to facilitate encoding SEQUENCE OF INTEGER as
*         OCTET STRING for ASN.1
*
****************************************************************************/

NLM_EXTERN Nlm_Uint2 LIBCALL Nlm_BSGetUint2 (Nlm_ByteStorePtr bsp);
NLM_EXTERN Nlm_Uint4 LIBCALL Nlm_BSGetUint4 (Nlm_ByteStorePtr bsp);
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSPutUint2 (Nlm_ByteStorePtr bsp, Nlm_Uint2 value);
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSPutUint4 (Nlm_ByteStorePtr bsp, Nlm_Uint4 value);

       /* In functions below, size is 2 or 4
          Integers are converted from ByteStore to native endian (BSUintXRead)
	        or from native endian to ByteStore (BSUintXWrite)
	      len is number of INTEGERS, not number of bytes
	      returns count of integers put in ptr
	      WARNING: On LITTLE_ENDIAN machines the data in ptr is changed to BIG_ENDIAN in the
	        XXWrite functions and LEFT THAT WAY
	   */
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSUint4Read (Nlm_ByteStorePtr bsp, Nlm_Uint4Ptr ptr, Nlm_Int4 len);
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSUint4Write (Nlm_ByteStorePtr bsp, Nlm_Uint4Ptr ptr, Nlm_Int4 len);
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSUint2Read (Nlm_ByteStorePtr bsp, Nlm_Uint2Ptr ptr, Nlm_Int4 len);
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSUint2Write (Nlm_ByteStorePtr bsp, Nlm_Uint2Ptr ptr, Nlm_Int4 len);

/* for copying and swapping of UID lists passed over network to BIG_ENDIAN server */

NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSDupAndSwapUint4 (Nlm_ByteStorePtr source);


#define ByteStore Nlm_ByteStore
#define ByteStorePtr Nlm_ByteStorePtr
#define BSMerge Nlm_BSMerge
#define BSNew Nlm_BSNew
#define BSLen Nlm_BSLen
#define BSFree Nlm_BSFree
#define BSSeek Nlm_BSSeek
#define BSTell Nlm_BSTell
#define BSDelete Nlm_BSDelete
#define BSWrite Nlm_BSWrite
#define BSInsert Nlm_BSInsert
#define BSInsertFromBS Nlm_BSInsertFromBS
#define BSDup Nlm_BSDup
#define BSEqual Nlm_BSEqual
#define BSRead Nlm_BSRead
#define BSGetByte Nlm_BSGetByte
#define BSPutByte Nlm_BSPutByte
#define BSDupAndSwapUint4 Nlm_BSDupAndSwapUint4

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

