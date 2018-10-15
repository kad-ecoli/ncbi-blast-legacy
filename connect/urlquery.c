/* $Id: urlquery.c,v 6.61 2011/10/26 15:24:55 lavr Exp $
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
 * File Name:  urlquery.c
 *
 * Author:  Jonathan Kans
 *
 * Version Creation Date:   4/16/98
 *
 * $Revision: 6.61 $
 *
 * File Description: 
 *
 * Modifications:  
 * --------------------------------------------------------------------------
 *
 * ==========================================================================
 */

#include "asnbuild.h"
#include <urlquery.h>


#ifdef OS_MAC
#include <Events.h>

NLM_EXTERN void QUERY_WaitForNextMacEvent (void)
{
  EventRecord  currEvent;

  WaitNextEvent (0, &currEvent, 0, NULL);
}
#endif


/* Set HTTP user header */
static void x_SetupUserHeader (
  SConnNetInfo*  net_info,
  const char*    appName,
  EMIME_Type     type,
  EMIME_SubType  subtype,
  EMIME_Encoding encoding
)
{
  const char* userAgentName = NULL;
  char        user_header [MAX_CONTENT_TYPE_LEN + 80];

  /* content-type if specified */
  if (type < eMIME_T_Unknown) {
    VERIFY( MIME_ComposeContentTypeEx (type, subtype, encoding,
                                       user_header, MAX_CONTENT_TYPE_LEN) );
    ConnNetInfo_OverrideUserHeader (net_info, user_header);
  }

  /* allow the user to specify a prog. name, otherwise get it from elsewhere */
  if (StringHasNoText (appName)) {
    const char* progName = GetProgramName ();
    if (StringHasNoText (progName)) {
      char path [PATH_MAX];
      Nlm_ProgramPath (path, sizeof (path));
      userAgentName = StringRChr (path, DIRDELIMCHR);
      if (userAgentName)
        ++userAgentName;
    } else
      userAgentName = progName;
  } else
    userAgentName = appName;
  if (StringDoesHaveText (userAgentName)) {
    sprintf (user_header, "User-Agent: %.80s\r\n", userAgentName);
    ConnNetInfo_ExtendUserHeader (net_info, user_header);
  }
}


NLM_EXTERN CONN QUERY_OpenUrlQuery (
  const char*     host_machine,
  Nlm_Uint2       host_port,
  const char*     host_path,
  const char*     arguments,
  const char*     appName,
  Nlm_Uint4       timeoutsec,
  EMIME_Type      type,
  EMIME_SubType   subtype,
  EMIME_Encoding  encoding,
  THTTP_Flags     flags
)
{
  CONN           conn;
  CONNECTOR      connector;
  SConnNetInfo*  net_info;
  EIO_Status     status;

  if (StringHasNoText (host_path))
    return NULL;

  /* fill in connection info fields and create the connection */
  net_info = ConnNetInfo_Create (0);
  ASSERT ( net_info );

  x_SetupUserHeader (net_info, appName, type, subtype, encoding);

  if (StringDoesHaveText (host_machine)) {
    StringNCpy_0 (net_info->host, host_machine, sizeof (net_info->host));
  }
  if ( host_port ) {
    net_info->port = host_port;
  }
  StringNCpy_0 (net_info->path, host_path, sizeof (net_info->path));
  if (StringDoesHaveText (arguments)) {
    StringNCpy_0 (net_info->args, arguments, sizeof (net_info->args));
  }

  if (timeoutsec == (Nlm_Uint4)(-1L)) {
    net_info->timeout  = kInfiniteTimeout;
  } else if ( timeoutsec ) {
    net_info->tmo.sec  = timeoutsec;
    net_info->tmo.usec = 0;
    net_info->timeout  = &net_info->tmo;
  }

  connector = HTTP_CreateConnector (net_info, NULL, flags);

  ConnNetInfo_Destroy (net_info);

  if (connector == NULL) {
    ErrPostEx (SEV_ERROR, 0, 0, "QUERY_OpenUrlQuery failed in HTTP_CreateConnector");
    conn = NULL;
  } else if ((status = CONN_Create (connector, &conn)) != eIO_Success) {
    ErrPostEx (SEV_ERROR, 0, 0, "QUERY_OpenUrlQuery failed in CONN_Create: %s",
              IO_StatusStr (status));
    ASSERT (conn == NULL);
  }

  return conn;
}


NLM_EXTERN CONN QUERY_OpenServiceQueryEx (
  const char* service,
  const char* parameters,
  Nlm_Uint4   timeoutsec,
  const char* arguments
)
{
  CONN           conn;
  CONNECTOR      connector;
  SConnNetInfo*  net_info;
  size_t         n_written;
  EIO_Status     status;

  /* fill in connection info fields and create the connection */
  net_info = ConnNetInfo_Create (service);
  ASSERT ( net_info );

  /* let the user agent be set with a program name */
  x_SetupUserHeader (net_info,
                     NULL, eMIME_T_Unknown, eMIME_Unknown, eENCOD_None);

  if (timeoutsec == (Nlm_Uint4)(-1L)) {
    net_info->timeout  = kInfiniteTimeout;
  } else if ( timeoutsec ) {
    net_info->tmo.sec  = timeoutsec;
    net_info->tmo.usec = 0;
    net_info->timeout  = &net_info->tmo;
  }

  ConnNetInfo_PostOverrideArg (net_info, arguments, 0);

  connector = SERVICE_CreateConnectorEx (service, fSERV_Any, net_info, 0);

  ConnNetInfo_Destroy (net_info);

  if (connector == NULL) {
    ErrPostEx (SEV_ERROR, 0, 0, "QUERY_OpenServiceQuery failed in SERVICE_CreateConnectorEx");
    conn = NULL;
  } else if ((status = CONN_Create (connector, &conn)) != eIO_Success) {
    ErrPostEx (SEV_ERROR, 0, 0, "QUERY_OpenServiceQuery failed in CONN_Create:"
               " %s", IO_StatusStr (status));
    ASSERT (conn == NULL);
  } else if (StringDoesHaveText (parameters)) {
    status = CONN_Write (conn, parameters, StringLen (parameters),
                         &n_written, eIO_WritePersist);
    if (status != eIO_Success) {
      ErrPostEx (SEV_ERROR, 0, 0, "QUERY_OpenServiceQuery failed to write service parameters in CONN_Write: %s", IO_StatusStr (status));
      CONN_Close (conn);
      conn = NULL;
    }
  }

  return conn;
}


NLM_EXTERN CONN QUERY_OpenServiceQuery (
  const char* service,
  const char* parameters,
  Nlm_Uint4   timeoutsec
)
{
  return QUERY_OpenServiceQueryEx (service, parameters, timeoutsec, 0);
}


NLM_EXTERN EIO_Status QUERY_SendQuery (
  CONN conn
)

{
  static const STimeout kPollTimeout = { 0 };
  EIO_Status            status;

  if (conn == NULL) return eIO_Closed;

  /* flush buffer, sending query, without waiting for response */
  status = CONN_Wait (conn, eIO_Read, &kPollTimeout);
  return status == eIO_Timeout ? eIO_Success : status;
}


#define URL_QUERY_BUFLEN  4096

NLM_EXTERN void QUERY_CopyFileToQuery (
  CONN conn,
  FILE *fp
)
{
  Char         buffer [URL_QUERY_BUFLEN + 4];
  size_t       ct;
  size_t       n_written;
  EIO_Status   status;

  if (conn == NULL || fp == NULL) return;

  while ((ct = FileRead (buffer, 1, URL_QUERY_BUFLEN, fp)) > 0) {
    status = CONN_Write (conn, (const void *) buffer, ct,
                         &n_written, eIO_WritePersist);
    if (status != eIO_Success) break;
  }
}


NLM_EXTERN void QUERY_CopyResultsToFile (
  CONN conn,
  FILE *fp
)
{
  Char         buffer [URL_QUERY_BUFLEN + 4];
  size_t       n_read;
  EIO_Status   status;

  if (conn == NULL || fp == NULL) return;

  do {
    status = CONN_Read (conn, buffer, URL_QUERY_BUFLEN, &n_read, eIO_ReadPlain);
    if ( n_read ) {
      FileWrite (buffer, 1, n_read, fp);
    }
  } while (status == eIO_Success);
}


NLM_EXTERN CharPtr QUERY_CopyResultsToString (
  CONN conn
)

{
  Char        buffer [URL_QUERY_BUFLEN + 4];
  ValNodePtr  head = NULL, last = NULL, vnp;
  size_t      n_read;
  EIO_Status  status;
  CharPtr     str;

  if (conn == NULL) return NULL;

  do {
    status = CONN_Read (conn, buffer, URL_QUERY_BUFLEN, &n_read, eIO_ReadPlain);
    if ( n_read ) {
      buffer [n_read] = '\0';
      vnp = ValNodeCopyStr (&last, 0, buffer);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
    }
  } while (status == eIO_Success);

  if (head == NULL) return NULL;

  str = ValNodeMergeStrs (head);
  ValNodeFreeData (head);

  return str;
}


static Nlm_Int2 LIBCALL AsnIoConnWrite (
  Pointer ptr,
  CharPtr buf,
  Nlm_Uint2 count
)
{
  size_t        bytes;
  AsnIoConnPtr  aicp;

  aicp = (AsnIoConnPtr) ptr;
  if (aicp == NULL || aicp->conn == NULL) return 0;
  CONN_Write (aicp->conn, buf, (size_t) count, &bytes, eIO_WritePersist);
  return (Nlm_Int2) bytes;
}


static Nlm_Int2 LIBCALL AsnIoConnRead (
  Pointer ptr,
  CharPtr buf,
  Nlm_Uint2 count
)
{
  size_t        bytes;
  AsnIoConnPtr  aicp;

  aicp = (AsnIoConnPtr) ptr;
  if (aicp == NULL || aicp->conn == NULL) return 0;
  CONN_Read (aicp->conn, buf, (size_t) count, &bytes, eIO_ReadPlain);
  return (Nlm_Int2) bytes;
}


NLM_EXTERN AsnIoConnPtr QUERY_AsnIoConnOpen (
  const char* mode,
  CONN conn
)
{
  Int1          type;
  AsnIoConnPtr  aicp;

  if (strcmp(mode, "r") == 0)
    type = (ASNIO_IN | ASNIO_TEXT);
  else if (strcmp(mode, "rb") == 0)
    type = (ASNIO_IN | ASNIO_BIN);
  else if (strcmp(mode, "w") == 0)
    type = (ASNIO_OUT | ASNIO_TEXT);
  else if (strcmp(mode, "wb") == 0)
    type = (ASNIO_OUT | ASNIO_BIN);
  else {
    AsnIoErrorMsg (NULL, 81, mode);
    return NULL;
  }

  aicp = (AsnIoConnPtr) MemNew (sizeof (AsnIoConn));
  aicp->aip = AsnIoNew (type, NULL, (Pointer) aicp, AsnIoConnRead, AsnIoConnWrite);
  aicp->conn = conn;
  return aicp;
}


NLM_EXTERN AsnIoConnPtr QUERY_AsnIoConnClose (
  AsnIoConnPtr aicp
)
{
  if (aicp == NULL) return NULL;
  AsnIoClose (aicp->aip);
  aicp->conn = NULL;
  return (AsnIoConnPtr) MemFree (aicp);
}


typedef struct SQueueTag {
  CONN                conn;
  QueryResultProc     resultproc;
  Nlm_VoidPtr         userdata;
  Nlm_Boolean         closeConn;
  Nlm_Boolean         protect;
  struct SQueueTag*   next;
} SConnQueue, PNTR QueuePtr;


NLM_EXTERN void QUERY_AddToQueue (
  QUEUE* queue,
  CONN conn,
  QueryResultProc resultproc,
  Nlm_VoidPtr userdata,
  Nlm_Boolean closeConn
)
{
  QueuePtr       cqp;
  QueuePtr PNTR  qptr;

  ASSERT (queue != NULL);

  if (conn == NULL || resultproc == NULL) return;

  /* allocate queue element */

  cqp = (QueuePtr) MemNew (sizeof (SConnQueue));
  if (cqp == NULL) return;

  cqp->conn = conn;
  cqp->resultproc = resultproc;
  cqp->userdata = userdata;
  cqp->closeConn = closeConn;
  cqp->protect = FALSE;

  /* add to polling queue */

  for (qptr = (QueuePtr PNTR) queue;  *qptr != NULL;  qptr = &(*qptr)->next);
  *qptr = cqp;
}


static void QUERY_RemoveFromQueue (
  QUEUE* queue,
  CONN conn
)
{
  QueuePtr       curr;
  QueuePtr       next;
  QueuePtr PNTR  prev;
  QueuePtr PNTR  qptr;

  qptr = (QueuePtr PNTR) queue;
  if (qptr == NULL || *qptr == NULL || conn == NULL) return;

  prev = qptr;
  curr = *qptr;

  while (curr != NULL) {
    next = curr->next;
    if (curr->conn == conn) {
      *prev = next;
      curr->next = NULL;
      MemFree (curr);
    } else {
      prev = &(curr->next);
    }
    curr = next;
  }
}


NLM_EXTERN Nlm_Int4 QUERY_CheckQueue (
  QUEUE* queue
)
{
  static const STimeout kPollTimeout = { 0 };
  Nlm_Int4       count = 0;
  QueuePtr       curr;
  QueuePtr       next;
  QueuePtr PNTR  qptr;
  EIO_Status     status;

  qptr = (QueuePtr PNTR) queue;
  if (qptr == NULL || *qptr == NULL) return 0;

  curr = *qptr;

  while (curr != NULL) {
    next = curr->next;

    if (curr->conn != NULL && (! curr->protect)) {
      status = CONN_Wait (curr->conn, eIO_Read, &kPollTimeout);

      if (status == eIO_Success || status == eIO_Closed) {
        /* protect against reentrant calls if resultproc is GUI and processes timer */
        curr->protect = TRUE;
        if (curr->resultproc != NULL) {
          /* result could eventually be used to reconnect on timeout */
          curr->resultproc (curr->conn, curr->userdata, status);
        }
        if (curr->closeConn) {
          CONN_Close (curr->conn);
        }
        QUERY_RemoveFromQueue (queue, curr->conn);
      } else {
        count++;
      }
    }

    curr = next;
  }

  return count;
}


NLM_EXTERN Nlm_Int4 QUERY_QueueSize (
  QUEUE queue
)
{
  Nlm_Int4  count = 0;
  QueuePtr  qptr;

  for (qptr = (QueuePtr) queue;  qptr;  qptr = qptr->next) {
    count++;
  }

  return count;
}


NLM_EXTERN void QUERY_CloseQueue (
  QUEUE* queue
)
{
  QueuePtr       curr;
  QueuePtr       next;
  QueuePtr PNTR  qptr;

  qptr = (QueuePtr PNTR) queue;
  if (qptr == NULL || *qptr == NULL) return;

  curr = *qptr;

  while (curr != NULL) {
    next = curr->next;

    if (curr->conn != NULL) {
      CONN_Close (curr->conn);
      QUERY_RemoveFromQueue (queue, curr->conn);
    }

    curr = next;
  }
}

static CONN QUERY_SendSimpleUrlQuery (
  CharPtr machine,
  Uint2 port,
  CharPtr path,
  CharPtr prefix,
  CharPtr arguments,
  CharPtr suffix,
  CharPtr post
)

{
  Char     ch;
  CONN     conn;
  size_t   len;
  size_t   n_written;
  CharPtr  ptr;
  CharPtr  query = NULL;

  len = StringLen (prefix) + StringLen (arguments) + StringLen (suffix);
  if (len > 0) {
    query = (CharPtr) MemNew (len + 2);
    if (query != NULL) {
      StringCpy (query, prefix);
      StringCat (query, arguments);
      StringCat (query, suffix);
      ptr = query;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n') {
          *ptr = '+';
        }
        ptr++;
        ch = *ptr;
      }
    }
  }

  conn = QUERY_OpenUrlQuery (machine, port, path, query, NULL,
                             30, eMIME_T_Text, eMIME_Xml, eENCOD_Url,
                             fHTTP_SureFlush | /* fHTTP_UrlDecodeInput | */ fHTTP_UrlEncodeOutput);

  MemFree (query);
  if (conn == NULL) return NULL;

  if (StringDoesHaveText (post)) {
    CONN_Write (conn, (const void *) post, StringLen (post), &n_written, eIO_WritePersist);
  }

  QUERY_SendQuery (conn);

  return conn;
}

NLM_EXTERN CharPtr QUERY_UrlSynchronousQuery (
  CharPtr machine,
  Uint2 port,
  CharPtr path,
  CharPtr prefix,
  CharPtr arguments,
  CharPtr suffix,
  CharPtr post
)

{
  CONN        conn;
  time_t      currtime, starttime, max = 0;
  EIO_Status  status;
  CharPtr     str = NULL;
  STimeout    timeout;

  conn = QUERY_SendSimpleUrlQuery (machine, port, path, prefix, arguments, suffix, post);

  if (conn == NULL) return NULL;

#ifdef OS_MAC 
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 100;
  timeout.usec = 0;
#endif


  starttime = GetSecs ();
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) != eIO_Success && max < 300) {
    currtime = GetSecs ();
    max = currtime - starttime;
  }

  if (status == eIO_Success) {
    str = QUERY_CopyResultsToString (conn);
  }

  CONN_Close (conn);

  return str;
}

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
)

{
  CONN  conn;

  conn = QUERY_SendSimpleUrlQuery (machine, port, path, prefix, arguments, suffix, post);

  if (conn == NULL) return FALSE;

  QUERY_AddToQueue (queue, conn, resultproc, userdata, TRUE);

  return TRUE;
}

