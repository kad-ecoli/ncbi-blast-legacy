/*  $Id: hspfilter_queue.c,v 1.1 2009/06/01 13:54:56 maning Exp $
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
 * Author:  Ning Ma
 *
 */

/** @file hspfilter_queue.c
 * Default implementation of the BlastHSPWriter interface to save hits from
 * a BLAST search, and subsequently return them in sorted order.
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: hspfilter_queue.c,v 1.1 2009/06/01 13:54:56 maning Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */


#include <ncbithr.h>
#include <algo/blast/core/blast_hspstream.h>
#include <algo/blast/api/hspfilter_queue.h>
//#include <algo/blast/core/blast_util.h>

/** Data structure used by the writer */
typedef struct BlastHSPQueueData {
   ListNode * start;      /**< First element of the queue */
   ListNode * end;        /**< First element of the queue */
   Boolean    writeDone;  /**< Has writing to this stream been finished? */
   TNlmMutex  lock;       /**< reading/writing lock */
   TNlmSemaphore sema;    /**< Semaphore for reading */
} BlastHSPQueueData;

/*************************************************************/
/** The following are implementations for BlastHSPWriter ADT */

/** Perform pre-run stage-specific initialization 
 * @param data The internal data structure [in][out]
 * @param results The HSP results to operate on  [in]
 */ 
static int 
s_BlastHSPQueueInit(void* data, BlastHSPResults* results)
{
   BlastHSPQueueData * q_data = data;
   return 0;
}

/** Perform post-run clean-ups
 * @param data The buffered data structure [in]
 * @param results The HSP results to propagate [in][out]
 */ 
static int 
s_BlastHSPQueueFinal(void* data, BlastHSPResults* results)
{
   BlastHSPQueueData * q_data = data;

   NlmMutexLockEx(&q_data->lock);
   q_data->writeDone = TRUE;
   NlmSemaPost(q_data->sema);
   NlmMutexUnlock(q_data->lock);

   return 0;
}

/** Perform writing task
 * ownership of the HSP list and sets the dereferenced pointer to NULL.
 * @param data To store results to [in][out]
 * @param hsp_list Pointer to the HSP list to save in the queue. [in]
 */
static int 
s_BlastHSPQueueRun(void* data, BlastHSPList* hsp_list)
{
   BlastHSPQueueData * q_data = data;

   if (!hsp_list)
      return 0;

   if (hsp_list->hspcnt == 0) {
      Blast_HSPListFree(hsp_list);
      return 0;
   }

   if (q_data->writeDone)
      return -1;

   NlmMutexLockEx(&q_data->lock);
   q_data->end = ListNodeAddPointer(&q_data->end, 0, (void *)hsp_list);
   if (!q_data->start)
      q_data->start = q_data->end;
   hsp_list = NULL;
   NlmSemaPost(q_data->sema);
   NlmMutexUnlock(q_data->lock);

   return 0; 
}

/** Free the writer 
 * @param writer The writer to free [in]
 * @return NULL.
 */
static
BlastHSPWriter*
s_BlastHSPQueueFree(BlastHSPWriter* writer) 
{
   ListNode * p;
   BlastHSPQueueData *q_data = writer->data;

   NlmSemaDestroy(q_data->sema);
   NlmMutexDestroy(q_data->lock);

   for(p = q_data->start; p; p = p->next) {
      p->ptr = (void *) Blast_HSPListFree((BlastHSPList*) p->ptr);
   }
   q_data->start = ListNodeFree(q_data->start);
   sfree(writer->data);
   sfree(writer);
   return NULL;
}

/** create the writer
 * @param params Pointer to the hit paramters [in]
 * @param query_info BlastQueryInfo (not used) [in]
 * @return writer
 */
static
BlastHSPWriter* 
s_BlastHSPQueueNew(void* params, BlastQueryInfo* query_info)
{
   BlastHSPWriter * writer = NULL;
   BlastHSPQueueData * data = NULL;

   /* allocate space for writer */
   writer = malloc(sizeof(BlastHSPWriter));

   /* fill up the function pointers */
   writer->InitFnPtr   = &s_BlastHSPQueueInit;
   writer->FinalFnPtr  = &s_BlastHSPQueueFinal;
   writer->FreeFnPtr   = &s_BlastHSPQueueFree;
   writer->RunFnPtr    = &s_BlastHSPQueueRun;

   /* allocate for data structure */
   writer->data = calloc(1,sizeof(BlastHSPQueueData));
   data = writer->data;
   data->sema = NlmSemaInit(0);
    
   return writer;
}

/*************************************************************/
/** The following are exported functions to be used by APP   */

BlastHSPQueueParams*
BlastHSPQueueParamsNew()
{
    return NULL;
}

BlastHSPQueueParams*
BlastHSPQueueParamsFree(BlastHSPQueueParams* opts)
{
    return NULL;
}

BlastHSPWriterInfo*
BlastHSPQueueInfoNew(BlastHSPQueueParams* params) {
   BlastHSPWriterInfo * writer_info =
                        malloc(sizeof(BlastHSPWriterInfo)); 
   writer_info->NewFnPtr = &s_BlastHSPQueueNew;
   writer_info->params = params;
   return writer_info;
}

/************************************************************/
/** The follwoing is added to support queue implementation  */

int BlastHSPQueueRead(void* data, BlastHSPList** hsp_list_out) 
{
   BlastHSPQueueData* q_data = (BlastHSPQueueData*) data;
   int status = kBlastHSPStream_Error;

   /* Lock the mutex */
   NlmMutexLockEx(&q_data->lock);

   if (!q_data->writeDone) {
      while (!q_data->writeDone && !q_data->start) {
         /* Decrement the semaphore count to 0, then wait for it to be 
          * incremented. Note that mutex must be locked whenever the 
          * contents of the stream are checked, but it must be unlocked
          * for the semaphore wait. */
         NlmMutexUnlock(q_data->lock);
         NlmSemaWait(q_data->sema);
         NlmMutexLockEx(&q_data->lock);
      }
   }

   if (!q_data->start) {
      /* Nothing in the queue, but no more writing to the queue is expected. */
      *hsp_list_out = NULL;
      status =  kBlastHSPStream_Eof;
   } else {
      ListNode* start_node = q_data->start;

      *hsp_list_out = (BlastHSPList*) start_node->ptr;

      q_data->start = start_node->next;
      start_node->next = NULL;
      ListNodeFree(start_node);
      if (!q_data->start)
         q_data->end = NULL;
      status = kBlastHSPStream_Success;
   }

   NlmMutexUnlock(q_data->lock);

   return status;
}

