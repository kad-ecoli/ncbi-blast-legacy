/* $Id: seqpanel.h,v 6.24 2010/12/06 17:16:46 bollin Exp $
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
* Author:  Vlad Lebedev
**
* ==========================================================================
*/

#ifndef seqpanel_h
#define seqpanel_h
extern PaneL CreateSeqViewPanel (GrouP g, Int2 w, Int2 h);
extern void UpdateSeqViewPanel (PaneL pnl);

extern void WriteAlignmentContiguousToFile
(SeqAlignPtr salp,
 FILE        *fp,
 Int4        seq_chars_per_row,
 Boolean     show_substitutions);
extern void 
WriteAlignmentInterleaveToFile 
(SeqAlignPtr salp,
 FILE        *fp,
 Int4        seq_chars_per_row,
 Boolean     show_substitutions); 
extern void 
WriteAlignmentInterleaveToFileEx 
(SeqAlignPtr salp,
 FILE        *fp,
 Int4        seq_chars_per_row,
 Boolean     show_substitutions,
 Boolean     show_coordinates);

extern ForM CreateSeqEditorWindow (Int2 left, Int2 top, CharPtr windowname, BioseqPtr bsp);

NLM_EXTERN SeqAlignPtr Sequin_GlobalAlign2Seq (BioseqPtr bsp1, BioseqPtr bsp2, BoolPtr revcomp);
extern ForM CreateAlnEditorWindow (Int2 left, Int2 top, CharPtr windowname, SeqAlignPtr salp, Uint2 entityID);
extern ForM 
CreateAlnEditorWindowEx 
(Int2 left, Int2 top, CharPtr windowname, SeqAlignPtr salp, Uint2 entityID,
 Nlm_ChangeNotifyProc on_close_func, Pointer on_close_data);
NLM_EXTERN void RemoveSeqEdCloseFunc (WindoW w);
NLM_EXTERN void SeqAlnWindowScrollToAlnPos (WindoW w, Int4 aln_pos, Int4 seq_num);

NLM_EXTERN CharPtr FeatureLocationAlignment (SeqFeatPtr sfp, SeqAlignPtr salp, Int4 begin, Int4 fin);
extern void FlipAlignment (SeqAlignPtr salp);

NLM_EXTERN Boolean CloseAlignmentEditor(Uint2 entityID, Uint4 itemID);

#endif
