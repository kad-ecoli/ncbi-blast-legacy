#ifndef CONNECT___NCBI_SENDMAIL__H
#define CONNECT___NCBI_SENDMAIL__H

/* $Id: ncbi_sendmail.h,v 6.21 2010/04/05 13:59:32 kazimird Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   Send mail (in accordance with RFC821 [protocol] and RFC822 [headers])
 *
 */

#include <connect/ncbi_types.h>


/** @addtogroup Sendmail
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/* Options apply to various fields of SSendMailInfo structure, below */
typedef enum {
    fSendMail_NoMxHeader       = (1 << 0), /* Don't add standard mail header,
                                            * just use what user provided    */
    fSendMail_StripNonFQDNHost = (1 << 8)  /* Strip host part in "from" field
                                            * if it does not look like an FQDN
                                            * (i.e. doesn't have at least two
                                            * domain name labels separated by a
                                            * dot); leave only username part */
} ESendMailOptions;
typedef unsigned int TSendMailOptions;     /* Bitwise OR of ESendMailOption  */


/* Define optional parameters for communication with sendmail
 */
typedef struct {
    unsigned int     magic_number;  /* Filled in by SendMailInfo_Init        */
    const char*      cc;            /* Carbon copy recipient(s)              */
    const char*      bcc;           /* Blind carbon copy recipient(s)        */
    char             from[1024];    /* Originator address                    */
    const char*      header;        /* Custom header fields ('\n'-separated) */
    size_t           body_size;     /* Message body size (if specified)      */
    const char*      mx_host;       /* Host to contact an MTA at             */
    short            mx_port;       /* Port to contact an MTA at             */
    STimeout         mx_timeout;    /* Timeout for all network transactions  */
    TSendMailOptions mx_options;    /* From the above                        */
} SSendMailInfo;


/* NOTE about recipient lists:
 * They are not parsed; valid recipient (according to the standard)
 * can be specified in the form "Name" <address>; recipients should
 * be separated by commas.  In case of address-only recipients (with no
 * "Name" part above), angle brackets around the address may be omitted.
 *
 * NOTE about message body size:
 * If not specified (0), and by default, the message body size is calculated
 * as strlen() of passed body argument, which thus must be '\0'-terminated.
 * Otherwise, exactly "body_size" bytes are read from the location pointed
 * to by "body" parameter, and are sent in the message.
 *
 * NOTE about MX header:
 * A message header is automatically prepended to a message that has the
 * fSendMail_NoMxHeader flag cleared in "mx_options" (default).  Otherwise,
 * only "header" part (if any) is sent to the MTA (mail tranfser agent), and
 * the rest of the header and the message body is passed "as is" following it.
 * Message header separator and proper message formatting is then
 * the caller's responsibility.
 */

/* Initialize SSendMailInfo structure, setting:
 *   'magic_number' to a proper value (verified by CORE_SendMailEx()!);
 *   'cc', 'bcc', 'header' to NULL (means no recipients/additional headers);
 *   'from' filled out using either the provided (non-empty) user name
 *          or the name of the current user if discovered, 'anonymous'
 *          otherwise, and host in the form: username@hostname; may be later
 *          reset by the application to "" for sending no-return messages
 *          (aka MAILER-DAEMON messages);
 *   'mx_*' filled out in accordance with some hard-coded defaults, which are
 *          very NCBI-specific; an outside application is likely to choose and
 *          use different values (at least for mx_host).
 * Return value equals the argument passed in.
 * Note: This call is the only valid way to properly init SSendMailInfo.
 */
extern NCBI_XCONNECT_EXPORT SSendMailInfo* SendMailInfo_InitEx
(SSendMailInfo*       info,
 const char*          user
 );

#define SendMailInfo_Init(info)  SendMailInfo_InitEx(info, 0)

/* Send a simple message to recipient(s) defined in 'to',
 * and having subject 'subject', which may be empty (both NULL and "" treated
 * equally as empty subjects), and message body 'body' (may be NULL/empty,
 * if not empty, lines are separated by '\n', must be '\0'-terminated).
 * Return value 0 means success; otherwise descriptive error message
 * gets returned.  Communicaiton parameters for connection with sendmail
 * are set using default values as described in SendMailInfo_Init().
 * NOTE:  Use of this call in out-of-house applications is discouraged as
 *        it is likely to fail since MTA communication parameters set
 *        to their defaults (which are NCBI-specific) may not be suitable.
 */
extern NCBI_XCONNECT_EXPORT const char* CORE_SendMail
(const char*          to,
 const char*          subject,
 const char*          body
 );

/* Send a message as in CORE_SendMail() but by explicitly specifying
 * all additional parameters of the message and the communication via
 * argument 'info'. In case of 'info' == NULL, the call is completely
 * equivalent to CORE_SendMail().
 * NB: Body is not neccessarily '\0'-terminated if 'info->body_size' specifies
 * non-zero message body size (see SSendMailInfo::body_size above).
 *
 * NOTE: Body can have additional header part if fSendMail_NoMxHeader is
 *       set in 'info->mx_options'.  Even if no additional headers are
 *       supplied, the body must provide the proper header / message text
 *       delimiter (an empty line), which will not be automatically inserted
 *       in the no-header (aka "as-is") mode.
 */
extern NCBI_XCONNECT_EXPORT const char* CORE_SendMailEx
(const char*          to,
 const char*          subject,
 const char*          body,
 const SSendMailInfo* info
 );


#ifdef __cplusplus
}  /* extern "C" */
#endif


/* @} */

#endif /* CONNECT___NCBI_SENDMAIL__H */
