/* $Id: ncbi_base64.c,v 1.1 2010/06/08 15:24:52 lavr Exp $
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
 *   BASE-64 Encoding/Decoding
 *
 */

#include "ncbi_base64.h"


extern void BASE64_Encode
(const void* src_buf,
 size_t      src_size,
 size_t*     src_read,
 void*       dst_buf,
 size_t      dst_size,
 size_t*     dst_written,
 size_t*     line_len)
{
    static const char syms[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ" /*26*/
        "abcdefghijklmnopqrstuvwxyz" /*52*/
        "0123456789+/";              /*64*/
    const size_t max_len = line_len ? *line_len : 76;
    const size_t max_src =
        ((dst_size - (max_len ? dst_size/(max_len + 1) : 0)) >> 2) * 3;
    unsigned char* src = (unsigned char*) src_buf;
    unsigned char* dst = (unsigned char*) dst_buf;
    size_t len = 0, i = 0, j = 0;
    unsigned char temp = 0, c;
    unsigned char shift = 2;
    if (!max_src  ||  !src_size) {
        *src_read    = 0;
        *dst_written = 0;
        if (dst_size > 0) {
            *dst = '\0';
        }
        return;
    }
    if (src_size > max_src) {
        src_size = max_src;
    }
    c = src[0];
    for (;;) {
        unsigned char bits = (c >> shift) & 0x3F;
        if (max_len  &&  len >= max_len) {
            dst[j++] = '\n';
            len = 0;
        }
        _ASSERT((size_t)(temp | bits) < sizeof(syms) - 1);
        dst[j++] = syms[temp | bits];
        len++;
        if (i >= src_size) {
            break;
        }
        shift += 2;
        shift &= 7;
        temp = (c << (8 - shift)) & 0x3F;
        if (shift) {
            c = ++i < src_size ? src[i] : 0;
        } else if (i + 1 == src_size) {
            i++;
        }
    }
    _ASSERT(j <= dst_size);
    *src_read = i;
    for (i = 0; i < (3 - src_size % 3) % 3; i++) {
        if (max_len  &&  len >= max_len) {
            dst[j++] = '\n';
            len = 0;
        }
        dst[j++] = '=';
        len++;
    }
    _ASSERT(j <= dst_size);
    *dst_written = j;
    if (j < dst_size) {
        dst[j] = '\0';
    }
}


extern int/*bool*/ BASE64_Decode
(const void* src_buf,
 size_t      src_size,
 size_t*     src_read,
 void*       dst_buf,
 size_t      dst_size,
 size_t*     dst_written)
{
    unsigned char* src = (unsigned char*) src_buf;
    unsigned char* dst = (unsigned char*) dst_buf;
    size_t i = 0, j = 0, k = 0, l;
    unsigned int temp = 0;
    if (src_size < 4  ||  dst_size < 3) {
        *src_read    = 0;
        *dst_written = 0;
        return 0/*false*/;
    }
    for (;;) {
        int/*bool*/  ok = i < src_size ? 1/*true*/ : 0/*false*/;
        unsigned char c = ok ? src[i++] : '=';
        if (c == '=') {
            c  = 64; /*end*/
        } else if (c >= 'A'  &&  c <= 'Z') {
            c -= 'A';
        } else if (c >= 'a'  &&  c <= 'z') {
            c -= 'a' - 26;
        } else if (c >= '0'  &&  c <= '9') {
            c -= '0' - 52;
        } else if (c == '+') {
            c  = 62;
        } else if (c == '/') {
            c  = 63;
        } else {
            continue;
        }
        temp <<= 6;
        temp  |= c & 0x3F;
        if (!(++k & 3)  ||  c == 64) {
            if (c == 64) {
                if (k < 2) {
                    if (ok) {
                        /* pushback leading '=' */
                        --i;
                    }
                    break;
                }
                switch (k) {
                case 2:
                    temp >>= 4;
                    break;
                case 3:
                    temp >>= 10;
                    break;
                case 4:
                    temp >>= 8;
                    break;
                default:
                    _ASSERT(0);
                    break;
                }
                l = 4 - k;
                while (l > 0) {
                    /* eat up '='-padding */
                    if (i >= src_size)
                        break;
                    if (src[i] == '=')
                        l--;
                    else if (src[i] != '\r'  &&  src[i] != '\n')
                        break;
                    i++;
                }
            } else {
                k = 0;
            }
            switch (k) {
            case 0:
                dst[j++] = (temp & 0xFF0000) >> 16;
                /*FALLTHRU*/;
            case 4:
                dst[j++] = (temp & 0xFF00) >> 8;
                /*FALLTHRU*/
            case 3:
                dst[j++] = (temp & 0xFF);
                break;
            default:
                break;
            }
            if (j + 3 >= dst_size  ||  c == 64) {
                break;
            }
            temp = 0;
        }
    }
    *src_read    = i;
    *dst_written = j;
    return i  &&  j ? 1/*true*/ : 0/*false*/;
}
