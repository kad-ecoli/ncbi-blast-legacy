/*   urlquery.h
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
 * File Name:  urlquery.h
 *
 * Author:  Jonathan Kans
 *
 * Version Creation Date:   4/16/98
 *
 * $Revision: 6.22 $
 *
 * File Description: 
 *
 * Modifications:  
 * --------------------------------------------------------------------------
 *
 * ==========================================================================
 */

#ifndef _URLQUERY_
#define _URLQUERY_

#include <ncbi.h>
#include <connect/ncbi_http_connector.h>
#include <connect/ncbi_service_connector.h>
#include <asn.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif


/* URL QUERY CONVENIENCE FUNCTIONS */

/*
  Initializes a POST connector based on URL query parameters.  For example:

  QUERY_OpenUrlQuery (
    "www.ncbi.nlm.nih.gov", 80,
    "/cgi-bin/Sequin/testcgi.cgi",
    "request=seg&window=12&lowcut=2.3&hicut=2.6",
    "Sequin", 30, eMIME_T_NcbiData, eMIME_Fasta, eENCOD_Url,
    fHTTP_SureFlush | fHTTP_UrlDecodeInput | fHTTP_UrlEncodeOutput
  );

  The returned CONN value is then passed data before being sent to the cgi.
*/
NLM_EXTERN CONN QUERY_OpenUrlQuery (
  const char* host_machine,
  Nlm_Uint2 host_port,
  const char* host_path,
  const char* arguments,
  const char* appName,
  Nlm_Uint4 timeoutsec,
  EMIME_Type type,
  EMIME_SubType subtype,
  EMIME_Encoding encoding,
  THTTP_Flags flags
);

/*
  Returns connection to NCBI named service.  Pass parameters (if any)
  via the connection to the service (if successful).  Optionally,
  arguments may be passed to services that support connections via
  URLs (note that unlike the service parameters there is no guarantee
  that the arguments will be passed since services can be implemented
  in non-URL-compatible manners, like e.g. via standalone servers).
  Return NULL on error.
*/

NLM_EXTERN CONN QUERY_OpenServiceQueryEx (
  const char* service,
  const char* parameters,
  Nlm_Uint4 timeoutsec,
  const char* arguments
);

/*
  Same as QUERY_OpenServiceQueryEx(service, parameters, timeoutsec, 0);
*/
NLM_EXTERN CONN QUERY_OpenServiceQuery (
  const char* service,
  const char* parameters,
  Nlm_Uint4 timeoutsec
);

/*
  Calculates length of data written to connection, writes HTML header, calls
  CONN_Wait with a timeout of 0 to send the query without waiting for response.
  This should be called after all data are copied to the connection, typically
  either with CONN_Read, QUERY_CopyFileToQuery, or with QUERY_AsnIoConnOpen
  followed by a specific AsnWrite function.
*/
NLM_EXTERN EIO_Status QUERY_SendQuery (
  CONN conn
);

/*
  Copies file to CONN by repeated calls to CONN_Write.  Writing of data must be
  done before the query is sent.  For the above seg example, a FASTA file of
  protein sequence data is appropriate.
*/
NLM_EXTERN void QUERY_CopyFileToQuery (
  CONN conn,
  FILE *fp
);

/*
  Copies results from CONN by repeated calls to CONN_Read.  Reading of data is
  only done when the connection is ready (CONN_Wait returns eCONN_Success).  In
  the query system below, the user's completion routine is called when results
  are ready to be read.  For the seg example, the results are in FASTA format
  with lower case x characters replacing amino acids in low-complexity regions.
*/
NLM_EXTERN void QUERY_CopyResultsToFile (
  CONN conn,
  FILE *fp
);

NLM_EXTERN CharPtr QUERY_CopyResultsToString (
  CONN conn
);

/* Structure for direct AsnIo-CONN link */
typedef struct asnioconn {    /* for AsnIo to and from a connection */
  AsnIoPtr aip;
  CONN conn;
} AsnIoConn, PNTR AsnIoConnPtr;

/*
  Creates structure containing AsnIoPtr attached to CONN.  Pass aicp->aip
  to AsnRead and AsnWrite functions.
*/
NLM_EXTERN AsnIoConnPtr QUERY_AsnIoConnOpen (
  const char* mode,
  CONN conn
);

/*
  Frees AsnIoPtr and structure, CONN will be closed later by QUERY_CheckQueue.
*/
NLM_EXTERN AsnIoConnPtr QUERY_AsnIoConnClose (
  AsnIoConnPtr aicp
);

/* FUNCTIONS FOR MAINTAINING A QUEUE OF PENDING URL QUERIES */

/* Callback type for queued queries */
typedef Nlm_Boolean (LIBCALLBACK *QueryResultProc) (
  CONN conn, Nlm_VoidPtr userdata, EIO_Status status
);

/* Opaque handle type.
   Variable must be kept by application and initialized to NULL.
*/
struct SQueueTag;
typedef struct SQueueTag* QUEUE;  /* queue handle */

/*
  Return number of currently queued connections.
*/
NLM_EXTERN Nlm_Int4 QUERY_QueueSize (
  QUEUE q
);

/*
  Records connection, completion routine, and user data in queue.
*/
NLM_EXTERN void QUERY_AddToQueue (
  QUEUE* q,
  CONN conn,
  QueryResultProc resultproc,
  Nlm_VoidPtr userdata,
  Nlm_Boolean closeConn
);

/*
  Checks queued connections (with CONN_Wait), calls completion routine, then
  removes query from queue and closes connection (if closeConn was TRUE).
  Application is responsible for calling QUERY_CheckQueue() every once in a
  while, typically with a timer.  Returns the number of pending connections.
*/
NLM_EXTERN Nlm_Int4 QUERY_CheckQueue (
  QUEUE* q
);

/*
  Forces all connections to be closed, removes all queries from queue.
*/
NLM_EXTERN void QUERY_CloseQueue (
  QUEUE* q
);

#ifdef OS_MAC
NLM_EXTERN void QUERY_WaitForNextMacEvent (
  void
);
#endif

/* Simple query functions for access to EUtils services */

NLM_EXTERN CharPtr QUERY_UrlSynchronousQuery (
  CharPtr machine,
  Uint2 port,
  CharPtr path,
  CharPtr prefix,
  CharPtr arguments,
  CharPtr suffix,
  CharPtr post
);

NLM_EXTERN Boolean QUERY_UrlAsynchronousQuery (
  CharPtr machine,
  Uint2 port,
  CharPtr path,
  CharPtr prefix,
  CharPtr arguments,
  CharPtr suffix,
  CharPtr post,
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _URLQUERY_ */

