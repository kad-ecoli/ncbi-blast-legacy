/*  $Id: hspfilter_queue.h,v 1.1 2009/06/01 13:54:56 maning Exp $
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

/** @file hspfilter_queue.h
 * Implementation of HSPWriter in queue mode.
 */

#ifndef ALGO_BLAST_CORE__HSPFILTER_QUEUE__H
#define ALGO_BLAST_CORE__HSPFILTER_QUEUE__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hspfilter.h>
#include <algo/blast/core/blast_hits.h>
#include <connect/ncbi_core.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Probably need no options */
typedef struct BlastHSPQueueParams {
   EBlastProgramType program;/**< program type */
} BlastHSPQueueParams;

/** Sets up parameter set for use by queue.
 * @param program Blast program type.[in]
 * @param hit_options field hitlist_size and hsp_num_max needed, a pointer to 
 *      this structure will be stored on resulting structure.[in]
 * @param ext_options field compositionBasedStats needed here. [in]
 * @param scoring_options gapped_calculation needed here. [in]
 * @return the pointer to the allocated parameter
 */
NCBI_XBLAST_EXPORT
BlastHSPQueueParams*
BlastHSPQueueParamsNew();

/** Deallocates the BlastHSPQueueParams structure passed in
 * @param params structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
BlastHSPQueueParams*
BlastHSPQueueParamsFree(BlastHSPQueueParams* params);

/** WriterInfo to create a default writer: the collecter
 * @param params The queue parameters.
 * @return pointer to WriterInfo
 */
NCBI_XBLAST_EXPORT
BlastHSPWriterInfo* 
BlastHSPQueueInfoNew(BlastHSPQueueParams* params);


/** Read one HSP list from a queue of HSP lists. If the queue is empty, this 
 * function waits for more results to be written, unless results queue is 
 * already closed for writing.
 * @param hsp_stream HSP list stream to read from [in]
 * @param hsp_list_out The read HSP list. NULL, if there is nothing left 
 *                     in the queue to read.
 * @return Status: success, error or end of reading.
 */
NCBI_XBLAST_EXPORT
int 
BlastHSPQueueRead(void* data, BlastHSPList** hsp_list_out);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__HSPFILTER_QUEUE__H */
