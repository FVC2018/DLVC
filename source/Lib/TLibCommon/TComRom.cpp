/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2015, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TComRom.cpp
    \brief    global variables & functions
*/

#include "TComRom.h"
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#if COM16_C806_EMT || COM16_C806_T64
#include <math.h>
#endif
#include <iomanip>
#include <assert.h>
#include "TComDataCU.h"
#include "Debug.h"
// ====================================================================================================================
// Initialize / destroy functions
// ====================================================================================================================

//! \ingroup TLibCommon
//! \{
#if GuoLei
Bool g_new = false;
#endif
#if JVET_C0024_AMAX_BT
UInt g_uiBlkSize[ 10 ];
UInt g_uiNumBlk[ 10 ];
#if JVET_C0024_AMAX_BT_FIX
UInt g_uiPrevISlicePOC = 0;
Bool g_bInitAMaxBT = false;
#endif
#endif

#if VCEG_AZ08_KLT_COMMON
short **g_ppsEigenVector[USE_MORE_BLOCKSIZE_DEPTH_MAX];
#define MAX_KLTAREA (1<<(((USE_MORE_BLOCKSIZE_DEPTH_MAX)<<1) + 2))
#if VCEG_AZ08_INTER_KLT
Bool g_bEnableCheck = true;
#endif
Void reOrderCoeff(TCoeff *pcCoef, const UInt *scan, UInt uiWidth, UInt uiHeight)
{
    TCoeff coeff[MAX_KLTAREA];
    UInt uiMaxNumCoeff = uiWidth * uiHeight;
    memcpy(coeff, pcCoef, uiMaxNumCoeff*sizeof(TCoeff));

    for (UInt i = 0; i < uiMaxNumCoeff; i++)
    {
        pcCoef[scan[i]] = coeff[i];
    }
}
Void recoverOrderCoeff(TCoeff *pcCoef, const UInt *scan, UInt uiWidth, UInt uiHeight)
{
    TCoeff coeff[MAX_KLTAREA];
    UInt uiMaxNumCoeff = uiWidth * uiHeight;
    memcpy(coeff, pcCoef, uiMaxNumCoeff*sizeof(TCoeff));
    for (UInt i = 0; i < uiMaxNumCoeff; i++)
    {
        pcCoef[i] = coeff[scan[i]];
    }
}
#endif
#if VCEG_AZ08_INTRA_KLT
Int getZorder(Int iLCUX, Int iLCUY, Int NumInRow)
{
    //get raster id
    Int rasterId = (iLCUY >> 2)*NumInRow + (iLCUX >> 2);
    Int zOrder = g_auiRasterToZscan[rasterId];
    return zOrder;
}
#endif

const Char* nalUnitTypeToString(NalUnitType type)
{
  switch (type)
  {
  case NAL_UNIT_CODED_SLICE_TRAIL_R:    return "TRAIL_R";
  case NAL_UNIT_CODED_SLICE_TRAIL_N:    return "TRAIL_N";
  case NAL_UNIT_CODED_SLICE_TSA_R:      return "TSA_R";
  case NAL_UNIT_CODED_SLICE_TSA_N:      return "TSA_N";
  case NAL_UNIT_CODED_SLICE_STSA_R:     return "STSA_R";
  case NAL_UNIT_CODED_SLICE_STSA_N:     return "STSA_N";
  case NAL_UNIT_CODED_SLICE_BLA_W_LP:   return "BLA_W_LP";
  case NAL_UNIT_CODED_SLICE_BLA_W_RADL: return "BLA_W_RADL";
  case NAL_UNIT_CODED_SLICE_BLA_N_LP:   return "BLA_N_LP";
  case NAL_UNIT_CODED_SLICE_IDR_W_RADL: return "IDR_W_RADL";
  case NAL_UNIT_CODED_SLICE_IDR_N_LP:   return "IDR_N_LP";
  case NAL_UNIT_CODED_SLICE_CRA:        return "CRA";
  case NAL_UNIT_CODED_SLICE_RADL_R:     return "RADL_R";
  case NAL_UNIT_CODED_SLICE_RADL_N:     return "RADL_N";
  case NAL_UNIT_CODED_SLICE_RASL_R:     return "RASL_R";
  case NAL_UNIT_CODED_SLICE_RASL_N:     return "RASL_N";
  case NAL_UNIT_VPS:                    return "VPS";
  case NAL_UNIT_SPS:                    return "SPS";
  case NAL_UNIT_PPS:                    return "PPS";
  case NAL_UNIT_ACCESS_UNIT_DELIMITER:  return "AUD";
  case NAL_UNIT_EOS:                    return "EOS";
  case NAL_UNIT_EOB:                    return "EOB";
  case NAL_UNIT_FILLER_DATA:            return "FILLER";
  case NAL_UNIT_PREFIX_SEI:             return "Prefix SEI";
  case NAL_UNIT_SUFFIX_SEI:             return "Suffix SEI";
  default:                              return "UNK";
  }
}

class ScanGenerator
{
private:
  UInt m_line, m_column;
  const UInt m_blockWidth, m_blockHeight;
  const UInt m_stride;
  const COEFF_SCAN_TYPE m_scanType;

public:
  ScanGenerator(UInt blockWidth, UInt blockHeight, UInt stride, COEFF_SCAN_TYPE scanType)
    : m_line(0), m_column(0), m_blockWidth(blockWidth), m_blockHeight(blockHeight), m_stride(stride), m_scanType(scanType)
  { }

  UInt GetCurrentX() const { return m_column; }
  UInt GetCurrentY() const { return m_line; }

  UInt GetNextIndex(UInt blockOffsetX, UInt blockOffsetY)
  {
    Int rtn=((m_line + blockOffsetY) * m_stride) + m_column + blockOffsetX;

    //advance line and column to the next position
    switch (m_scanType)
    {
      //------------------------------------------------

      case SCAN_DIAG:
        {
          if ((m_column == (m_blockWidth - 1)) || (m_line == 0)) //if we reach the end of a rank, go diagonally down to the next one
          {
            m_line   += m_column + 1;
            m_column  = 0;

            if (m_line >= m_blockHeight) //if that takes us outside the block, adjust so that we are back on the bottom row
            {
              m_column += m_line - (m_blockHeight - 1);
              m_line    = m_blockHeight - 1;
            }
          }
          else
          {
            m_column++;
            m_line--;
          }
        }
        break;

      //------------------------------------------------

      case SCAN_HOR:
        {
          if (m_column == (m_blockWidth - 1))
          {
            m_line++;
            m_column = 0;
          }
          else
          {
            m_column++;
          }
        }
        break;

      //------------------------------------------------

      case SCAN_VER:
        {
          if (m_line == (m_blockHeight - 1))
          {
            m_column++;
            m_line = 0;
          }
          else
          {
            m_line++;
          }
        }
        break;

      //------------------------------------------------

      default:
        {
          std::cerr << "ERROR: Unknown scan type \"" << m_scanType << "\"in ScanGenerator::GetNextIndex" << std::endl;
          exit(1);
        }
        break;
    }

    return rtn;
  }
};

#if JVET_E0077_MMLM
Int g_aiLMDivTable[] = {
    65536, 32768, 21845, 16384, 13107, 10923, 9362, 8192, 7282, 6554, 5958, 5461, 5041, 4681, 4369, 4096, 3855, 3641, 3449, 3277, 3121, 2979, 2849, 2731, 2621, 2521, 2427, 2341, 2260, 2185, 2114, 2048, 1986, 1928, 1872, 1820, 1771, 1725, 1680, 1638, 1598, 1560, 1524, 1489, 1456, 1425, 1394, 1365, 1337, 1311, 1285, 1260, 1237, 1214, 1192, 1170, 1150, 1130, 1111, 1092, 1074, 1057, 1040, 1024, 1008, 993, 978, 964, 950, 936, 923, 910, 898, 886, 874, 862, 851, 840, 830, 819, 809, 799, 790, 780, 771, 762, 753, 745, 736, 728, 720, 712, 705, 697, 690, 683, 676, 669, 662, 655, 649, 643, 636, 630, 624, 618, 612, 607, 601, 596, 590, 585, 580, 575, 570, 565, 560, 555, 551, 546, 542, 537, 533, 529, 524, 520, 516, 512, 508, 504, 500, 496, 493, 489, 485, 482, 478, 475, 471, 468, 465, 462, 458, 455, 452, 449, 446, 443, 440, 437, 434, 431, 428, 426, 423, 420, 417, 415, 412, 410, 407, 405, 402, 400, 397, 395, 392, 390, 388, 386, 383, 381, 379, 377, 374, 372, 370, 368, 366, 364, 362, 360, 358, 356, 354, 352, 350, 349, 347, 345, 343, 341, 340, 338, 336, 334, 333, 331, 329, 328, 326, 324, 323, 321, 320, 318, 317, 315, 314, 312, 311, 309, 308, 306, 305, 303, 302, 301, 299, 298, 297, 295, 294, 293, 291, 290, 289, 287, 286, 285, 284, 282, 281, 280, 279, 278, 277, 275, 274, 273, 272, 271, 270, 269, 267, 266, 265, 264, 263, 262, 261, 260, 259, 258, 257, 256, 255, 254, 253, 252, 251, 250, 249, 248, 247, 246, 245, 245, 244, 243, 242, 241, 240, 239, 238, 237, 237, 236, 235, 234, 233, 232, 232, 231, 230, 229, 228, 228, 227, 226, 225, 224, 224, 223, 222, 221, 221, 220, 219, 218, 218, 217, 216, 216, 215, 214, 213, 213, 212, 211, 211, 210, 209, 209, 208, 207, 207, 206, 205, 205, 204, 204, 203, 202, 202, 201, 200, 200, 199, 199, 198, 197, 197, 196, 196, 195, 194, 194, 193, 193, 192, 192, 191, 191, 190, 189, 189, 188, 188, 187, 187, 186, 186, 185, 185, 184, 184, 183, 183, 182, 182, 181, 181, 180, 180, 179, 179, 178, 178, 177, 177, 176, 176, 175, 175, 174, 174, 173, 173, 172, 172, 172, 171, 171, 170, 170, 169, 169, 168, 168, 168, 167, 167, 166, 166, 165, 165, 165, 164, 164, 163, 163, 163, 162, 162, 161, 161, 161, 160, 160, 159, 159, 159, 158, 158, 158, 157, 157, 156, 156, 156, 155, 155, 155, 154, 154, 153, 153, 153, 152, 152, 152, 151, 151, 151, 150, 150, 150, 149, 149, 149, 148, 148, 148, 147, 147, 147, 146, 146, 146, 145, 145, 145, 144, 144, 144, 143, 143, 143, 142, 142, 142, 142, 141, 141, 141, 140, 140, 140, 139, 139, 139, 139, 138, 138, 138, 137, 137, 137, 137, 136, 136, 136, 135, 135, 135, 135, 134, 134, 134, 133, 133, 133, 133, 132, 132, 132, 132, 131, 131, 131, 131, 130, 130, 130, 130, 129, 129, 129, 129, 128, 128,
};

Int g_aiLMDivTableHigh[] = {
    65536, 32768, 21845, 16384, 13107, 10922, 9362, 8192, 7281, 6553, 5957, 5461, 5041, 4681, 4369, 4096, 3855, 3640, 3449, 3276, 3120, 2978, 2849, 2730, 2621, 2520, 2427, 2340, 2259, 2184, 2114, 2048, 1985, 1927, 1872, 1820, 1771, 1724, 1680, 1638, 1598, 1560, 1524, 1489, 1456, 1424, 1394, 1365, 1337, 1310, 1285, 1260, 1236, 1213, 1191, 1170, 1149, 1129, 1110, 1092, 1074, 1057, 1040, 1024, 1008, 992, 978, 963, 949, 936, 923, 910, 897, 885, 873, 862, 851, 840, 829, 819, 809, 799, 789, 780, 771, 762, 753, 744, 736, 728, 720, 712, 704, 697, 689, 682, 675, 668, 661, 655, 648, 642, 636, 630, 624, 618, 612, 606, 601, 595, 590, 585, 579, 574, 569, 564, 560, 555, 550, 546, 541, 537, 532, 528, 524, 520, 516, 512, 508, 504, 500, 496, 492, 489, 485, 481, 478, 474, 471, 468, 464, 461, 458, 455, 451, 448, 445, 442, 439, 436, 434, 431, 428, 425, 422, 420, 417, 414, 412, 409, 407, 404, 402, 399, 397, 394, 392, 390, 387, 385, 383, 381, 378, 376, 374, 372, 370, 368, 366, 364, 362, 360, 358, 356, 354, 352, 350, 348, 346, 344, 343, 341, 339, 337, 336, 334, 332, 330, 329, 327, 326, 324, 322, 321, 319, 318, 316, 315, 313, 312, 310, 309, 307, 306, 304, 303, 302, 300, 299, 297, 296, 295, 293, 292, 291, 289, 288, 287, 286, 284, 283, 282, 281, 280, 278, 277, 276, 275, 274, 273, 271, 270, 269, 268, 267, 266, 265, 264, 263, 262, 261, 260, 259, 258, 257, 256, 255, 254, 253, 252, 251, 250, 249, 248, 247, 246, 245, 244, 243, 242, 241, 240, 240, 239, 238, 237, 236, 235, 234, 234, 233, 232, 231, 230, 229, 229, 228, 227, 226, 225, 225, 224, 223, 222, 222, 221, 220, 219, 219, 218, 217, 217, 216, 215, 214, 214, 213, 212, 212, 211, 210, 210, 209, 208, 208, 207, 206, 206, 205, 204, 204, 203, 202, 202, 201, 201, 200, 199, 199, 198, 197, 197, 196, 196, 195, 195, 194, 193, 193, 192, 192, 191, 191, 190, 189, 189, 188, 188, 187, 187, 186, 186, 185, 185, 184, 184, 183, 183, 182, 182, 181, 181, 180, 180, 179, 179, 178, 178, 177, 177, 176, 176, 175, 175, 174, 174, 173, 173, 172, 172, 172, 171, 171, 170, 170, 169, 169, 168, 168, 168, 167, 167, 166, 166, 165, 165, 165, 164, 164, 163, 163, 163, 162, 162, 161, 161, 161, 160, 160, 159, 159, 159, 158, 158, 157, 157, 157, 156, 156, 156, 155, 155, 154, 154, 154, 153, 153, 153, 152, 152, 152, 151, 151, 151, 150, 150, 149, 149, 149, 148, 148, 148, 147, 147, 147, 146, 146, 146, 145, 145, 145, 144, 144, 144, 144, 143, 143, 143, 142, 142, 142, 141, 141, 141, 140, 140, 140, 140, 139, 139, 139, 138, 138, 138, 137, 137, 137, 137, 136, 136, 136, 135, 135, 135, 135, 134, 134, 134, 134, 133, 133, 133, 132, 132, 132, 132, 131, 131, 131, 131, 130, 130, 130, 130, 129, 129, 129, 129, 128, 128, 128, 128,
};
Int g_aiLMDivTableLow[] = {
    0, 0, 21845, 0, 13107, 43690, 18724, 0, 50972, 39321, 53620, 21845, 15123, 9362, 4369, 0, 3855, 58254, 17246, 52428, 49932, 59578, 25644, 43690, 28835, 40329, 16990, 37449, 56496, 34952, 4228, 0, 61564, 34695, 29959, 29127, 15941, 41391, 26886, 26214, 28771, 24966, 6096, 29789, 23301, 45590, 25098, 21845, 30761, 47185, 1285, 20164, 34622, 41263, 36938, 18724, 49439, 61016, 51095, 17476, 23635, 2114, 16644, 0, 16131, 63550, 9781, 50115, 52238, 14979, 2769, 14563, 49376, 40738, 53302, 20695, 7660, 13443, 37330, 13107, 5663, 14385, 38689, 12483, 771, 3048, 18832, 47662, 23563, 11650, 11522, 22795, 45100, 12549, 55878, 43690, 41213, 48148, 64212, 23592, 57100, 33410, 17815, 10082, 9986, 17311, 31849, 53399, 16233, 51237, 27159, 9362, 63216, 57487, 57557, 63276, 8962, 25547, 47362, 8738, 40621, 11817, 53281, 33825, 18874, 8322, 2064, 0, 2032, 8065, 18009, 31775, 49275, 4890, 29612, 57825, 23918, 58887, 31589, 7489, 52056, 34152, 19248, 7281, 63728, 57456, 53944, 53137, 54979, 59419, 868, 10347, 22273, 36598, 53274, 6721, 27967, 51433, 11540, 39321, 3663, 35599, 4020, 39960, 12312, 52112, 28255, 6241, 51575, 33153, 16479, 1524, 53792, 42184, 32206, 23831, 17031, 11781, 8054, 5825, 5069, 5761, 7878, 11397, 16295, 22550, 30139, 39042, 49238, 60707, 7891, 21845, 37012, 53374, 5377, 24074, 43912, 64874, 21406, 44564, 3260, 28550, 54882, 16705, 45075, 8907, 39258, 5041, 37314, 4993, 39135, 8655, 44613, 15924, 53648, 26699, 604, 40884, 16458, 58386, 35585, 13579, 57895, 37449, 17767, 64376, 46192, 28743, 12019, 61546, 46244, 31638, 17720, 4481, 57448, 45541, 34288, 23681, 13710, 4369, 61185, 53078, 45578, 38676, 32366, 26640, 21491, 16912, 12896, 9437, 6527, 4161, 2331, 1032, 257, 0, 255, 1016, 2277, 4032, 6277, 9004, 12210, 15887, 20031, 24637, 29699, 35213, 41173, 47574, 54411, 61680, 3840, 11959, 20494, 29443, 38801, 48562, 58724, 3744, 14693, 26028, 37746, 49844, 62316, 9624, 22834, 36408, 50342, 64632, 13737, 28728, 44063, 59740, 10219, 26568, 43249, 60257, 12055, 29709, 47682, 434, 19033, 37941, 57155, 11136, 30953, 51067, 5938, 26637, 47624, 3360, 24916, 46751, 3328, 25716, 48376, 5770, 28967, 52428, 10616, 34599, 58840, 17799, 42547, 2010, 27256, 52748, 12947, 38924, 65140, 26056, 52743, 14127, 41277, 3120, 30726, 58555, 21072, 49344, 12300, 41007, 4394, 33530, 62876, 26896, 56659, 21092, 51264, 16103, 46678, 11915, 42886, 8515, 39875, 5890, 37632, 4027, 36145, 2912, 35400, 2534, 35385, 2880, 36089, 3939, 37500, 5698, 39605, 8147, 42395, 11275, 45857, 15069, 49982, 19521, 54758, 24619, 60175, 30353, 688, 36713, 7357, 43690, 14639, 51274, 22522, 59455, 30999, 2688, 40059, 12037, 49693, 21956, 59894, 32437, 5117, 43471, 16425, 55050, 28273, 1630, 40655, 14275, 53561, 27441, 1449, 41120, 15382, 55305, 29818, 4453, 44748, 19629, 60166, 35288, 10529, 51425, 26902, 2496, 43742, 19567, 61042, 37095, 13261, 55074, 31463, 7962, 50106, 26824, 3649, 46117, 23157, 302, 43088, 20442, 63436, 40997, 18660, 61961, 39826, 17792, 61393, 39557, 17819, 61715, 40171, 18724, 62908, 41651, 20489, 64956, 43980, 23096, 2304, 47139, 26529, 6009, 51115, 30773, 10519, 55890, 35811, 15819, 61448, 41628, 21892, 2240, 48208, 28724, 9322, 55538, 36301, 17144, 63604, 44608, 25692, 6855, 53632, 34952, 16349, 63360, 44911, 26539, 8242, 55557, 37410, 19338, 1340, 48951, 31099, 13320, 61149, 43513, 25949, 8456, 56569, 39216, 21932, 4718, 53109, 36031, 19022, 2080, 50741, 33933, 17191, 516, 49441, 32896, 16416, 0,
};

#endif


#if RRC_DWDIP

DECL_DWDIP_TABLE_SET(weights)
DECL_SELECTOR_FUN(weights)

DECL_DWDIP_TABLE_SET(SideDecision)
DECL_SELECTOR_FUN(SideDecision)

#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
DECL_DWDIP_TABLE_SET(SideDeltaPos)
DECL_SELECTOR_FUN(SideDeltaPos)
#endif

#endif


#if JVET_E0077_ENHANCED_LM

Bool   IsLMMode(UInt uiIntraMode)
{
    if (uiIntraMode == LM_CHROMA_IDX)
    {
        return true;
    }
#if JVET_E0077_MMLM
    else if (uiIntraMode == MMLM_CHROMA_IDX)
    {
        return true;
    }
#endif

#if JVET_E0077_LM_MF
    else if (uiIntraMode >= LM_CHROMA_F1_IDX && uiIntraMode < (LM_CHROMA_F1_IDX + LM_FILTER_NUM))
    {
        return true;
    }
#endif
    return false;
}


const Int g_aiNonLMPosThrs[] = { 3, 1, 0 };

const Int g_aiLAP_MinSize[] = { 11, 31 };
const Int g_aiMFLM_MinSize[] = { 0, 0 };
const Int g_aiMMLM_MinSize[] = { 0, 0 };
Int g_aiLMCodeWord[LM_SYMBOL_NUM][16];

#endif


#if RRC_DWDIP

Int getIntPart(Int pos)
{
  return(pos >> 5);//32 = 2^5  
}

Int getFractPart(Int pos)
{
  return(pos & 0x1F);//0x1F = 31 = 32 - 1  
}


Void xGetLUTsRow(
    Int y
  , Int deltaPos
  , Int width
  , Int height
  , Int * const piWeight
  , Int * const pbDecision
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
  , Int * const piSideDeltaPos
#endif 
)
{
  Int distGen, intGen, fractGen;

  const Int deltaInt = deltaPos >> 5;
  const Int deltaFract = deltaPos & (32 - 1);

  {
    Int  x, refMainIndex;
    for (x = 0; x < width; x++)
    {
      refMainIndex = x + deltaInt + 1;
      { //params derivation

        //redefine actual main pixel in certain case:
        Bool bIsRightSideReferenced = false;
        Int intRightGen = -1;

        Double dSideDeltaPos;
        Int sideDeltaInt;

        //Bool bUseSideSamples = xGetRefSideParams(x, y, height, deltaPos, dSideDeltaPos, sideDeltaInt, bIsRightSideReferenced);
        Bool bUseSideSamples = false;
        {
          if (deltaPos)
          {
            // case '/': works for positive deltaPos only (derives left boundary offset from top boundary offset (deltaPos) )
            dSideDeltaPos = static_cast<Double>(((y + 1) << 5) * ((x+1) << 5)) / deltaPos;
            sideDeltaInt = getIntPart(static_cast<Int>(round(dSideDeltaPos)));
          }
          else
          {
            dSideDeltaPos = 0;
            sideDeltaInt = 0;
          }

          if (deltaPos >= 0)
          {
            //switch adjacent/opposite sides processing:
            bUseSideSamples = !((y << 5) + dSideDeltaPos > (height << 5) || dSideDeltaPos <= 0);
          }
        }

        Int refSideIndex = y + sideDeltaInt + 1; /// cASE: /

        Double dHeightWeight;
        Double dHeightRatio;

        //switch adjacent/opposite sides processing:
        //if ( refSideIndex >= (2*height+1) || dSideDeltaPos<=0 )
        if (!bUseSideSamples)
        {
          dHeightWeight = static_cast<Double> (y + 1) / height;
          dHeightRatio = static_cast<Double> (y + 1) / (height - y);

          Int iProjectedDeltaPos = deltaPos >= 0 ? (x + 1) : (width - x);

          Double dGenDeltaPos = deltaPos ? -static_cast<Double>(((y + 1) << 5) * (iProjectedDeltaPos << 5)) / deltaPos : 0;
          if ((y << 5) + dGenDeltaPos > (height << 5) || dGenDeltaPos <= 0) // bottom side referenced
          {
            distGen = (x << 5) + static_cast<Int>(-round(deltaPos / dHeightRatio));

            Double iUpperPart = ((x + 1) << 5) + deltaPos >= 0 ? (y + 1) << 5 : (((x + 1) << 5)*((y + 1) << 5) / abs(deltaPos));
            dHeightWeight = static_cast<Double> (iUpperPart) / (iUpperPart + ((height - y) << 5));
            dHeightRatio = static_cast<Double> (iUpperPart) / ((height - y) << 5);
          }
          else
          {
            distGen = ((width + height) << 5) - ((y << 5) + static_cast<Int>(dGenDeltaPos)); //right side referenced

            Int iUpperPart = (x << 5) + deltaPos >= 0 ? y : ((x << 5)*(y << 5) / abs(deltaPos)) >> 5;
            dHeightWeight = static_cast<Double> (iUpperPart + 1) / (iUpperPart + 1 + (dGenDeltaPos / 32));
            dHeightRatio = static_cast<Double> (iUpperPart + 1) / (1 + (dGenDeltaPos / 32));

          }
          intGen = getIntPart(distGen);
          fractGen = getFractPart(distGen);

          // adjust ratios to calculate distances correctly
          if (bIsRightSideReferenced)
          {
            dHeightWeight = static_cast<Double> (y - intRightGen) / (height - intRightGen);
            dHeightRatio = static_cast<Double> (y - intRightGen) / (height - y);
          }

          refMainIndex = intGen; //"+1" - no need to add 1, because first refgen pixel is not at -1
          assert(refMainIndex >= -1);
          assert(refMainIndex <= width + height);
        }
        else
        {
          if (bIsRightSideReferenced)
          {
            //right and left are referenced
            dHeightWeight = static_cast<Double> ((width - x)) / (width + 1);
          }
          else
            //top and left are referenced
            //no need to adjust dHeightWeight in this case
            dHeightWeight = static_cast<Double> ((y + 1) << 5) / std::min<Double>(height << 5, ((y + 1) << 5) + dSideDeltaPos);

          intGen = refSideIndex;
          fractGen = getFractPart(static_cast<Int>(round(dSideDeltaPos)));
        }
        const Int dWeightPrecisionBits = 6;

        piWeight[x] = static_cast<Int>(round(dHeightWeight * (1 << dWeightPrecisionBits)));
        pbDecision[x] = bUseSideSamples;

#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
        piSideDeltaPos[x] = (bUseSideSamples ? static_cast<Int>(round(dSideDeltaPos)) : distGen);
#endif
      }// param derivation

    }
  }
}

#endif // RRC_DWDIP


// initialize ROM variables
Void initROM()
{
  Int i, c;

  // g_aucConvertToBit[ x ]: log2(x/4), if x=4 -> 0, x=8 -> 1, x=16 -> 2, ...
  ::memset( g_aucConvertToBit,   -1, sizeof( g_aucConvertToBit ) );
  c=0;
#if JVET_C0024_QTBT
  for ( i=1<<MIN_CU_LOG2; i<=MAX_CU_SIZE; i*=2 )
#else
  for ( i=4; i<=MAX_CU_SIZE; i*=2 )
#endif
  {
    g_aucConvertToBit[ i ] = c;
    c++;
  }
#if JVET_C0024_QTBT
  if (MIN_CU_LOG2>=1)
  {
    c = -1;
    for ( i=(1<<(MIN_CU_LOG2-1)); i>=1; i>>=1 )
    {
      g_aucConvertToBit[ i ] = c;
      c--;
    }
  }
#endif

#if COM16_C806_EMT || COM16_C806_T64
#if JVET_C0024_QTBT
  c = 2;
#else
  c = 4;
#endif

  const Double PI = 3.14159265358979323846;

#if JVET_C0024_QTBT
  for ( i=0; i<7; i++ ) 
#else
  for ( i=0; i<5; i++ )
#endif
  {
    TMatrixCoeff *iT = NULL;
    const Double s = sqrt((Double)c) * ( 64 << COM16_C806_TRANS_PREC );

    switch(i)
    {
#if JVET_C0024_QTBT
      case 0: iT = g_aiTr2 [0][0]; break; 
      case 1: iT = g_aiTr4 [0][0]; break;
      case 2: iT = g_aiTr8[0][0]; break;
      case 3: iT = g_aiTr16[0][0]; break;
      case 4: iT = g_aiTr32[0][0]; break;
      case 5: iT = g_aiTr64[0][0]; break;
      case 6: iT = g_aiTr128[0][0]; break;
      case 7: exit(0); break;
#else
    case 0: iT = g_aiTr4 [0][0]; break;
    case 1: iT = g_aiTr8 [0][0]; break;
    case 2: iT = g_aiTr16[0][0]; break;
    case 3: iT = g_aiTr32[0][0]; break;
    case 4: iT = g_aiTr64[0][0]; break;
    case 5: exit(0); break;
#endif
    }

    for( Int k=0; k<c; k++ )
    {
      for( Int n=0; n<c; n++ )
      {
        Double w0, w1, v;

        // DCT-II
        w0 = k==0 ? sqrt(0.5) : 1;
        v = cos(PI*(n+0.5)*k/c ) * w0 * sqrt(2.0/c);
        iT[DCT2*c*c + k*c + n] = (Short) ( s * v + ( v > 0 ? 0.5 : -0.5) );

        // DCT-V
        w0 = ( k==0 ) ? sqrt(0.5) : 1.0;
        w1 = ( n==0 ) ? sqrt(0.5) : 1.0;
        v = cos(PI*n*k/(c-0.5)) * w0 * w1 * sqrt(2.0/(c-0.5));
        iT[DCT5*c*c + k*c + n] = (Short) ( s * v + ( v > 0 ? 0.5 : -0.5) );

        // DCT-VIII
        v = cos(PI*(k+0.5)*(n+0.5)/(c+0.5) ) * sqrt(2.0/(c+0.5));
        iT[DCT8*c*c + k*c + n] = (Short) ( s * v + ( v > 0 ? 0.5 : -0.5) );

        // DST-I
        v = sin(PI*(n+1)*(k+1)/(c+1)) * sqrt(2.0/(c+1));
        iT[DST1*c*c + k*c + n] = (Short) ( s * v + ( v > 0 ? 0.5 : -0.5) );

        // DST-VII
        v = sin(PI*(k+0.5)*(n+1)/(c+0.5)) * sqrt(2.0/(c+0.5));
        iT[DST7*c*c + k*c + n] = (Short) ( s * v + ( v > 0 ? 0.5 : -0.5) );
      }
    }
    c <<= 1;
  }
#endif

#if VCEG_AZ07_INTRA_4TAP_FILTER
  for( i=17; i<32; i++ )
  {
    for( c=0; c<4; c++ )
    {
      g_aiIntraCubicFilter[i][c] = g_aiIntraCubicFilter[32-i][3-c];
      g_aiIntraGaussFilter[i][c] = g_aiIntraGaussFilter[32-i][3-c];
    }
  }
#if DIMD_NUM_INTRA_DIR_INC
  for (i = 33; i < 64; i++)
  {
    for (c = 0; c < 4; c++)
    {
      g_aiExtIntraCubicFilter[i][c] = g_aiExtIntraCubicFilter[64 - i][3 - c];
      g_aiExtIntraGaussFilter[i][c] = g_aiExtIntraGaussFilter[64 - i][3 - c];
    }
  }
#endif
#endif

#if JVET_D0120_NSST_IMPROV
  const Double scPi = 2 * PI / NSST_HYGT_PTS;
  for ( i=0; i<NSST_HYGT_PTS; i++ )
  {
    g_tabSinCos[i].c = Int(floor(cos(i * scPi) * 1024.0 + 0.5));
    g_tabSinCos[i].s = Int(floor(sin(i * scPi) * 1024.0 + 0.5));
  }
#endif

  // initialise scan orders
#if COM16_C806_T64
  for(UInt log2BlockHeight = 0; log2BlockHeight < MAX_LOG2_TU_SIZE_PLUS_ONE; log2BlockHeight++)
  {
    for(UInt log2BlockWidth = 0; log2BlockWidth < MAX_LOG2_TU_SIZE_PLUS_ONE; log2BlockWidth++)
#else
  for(UInt log2BlockHeight = 0; log2BlockHeight < MAX_CU_DEPTH; log2BlockHeight++)
  {
    for(UInt log2BlockWidth = 0; log2BlockWidth < MAX_CU_DEPTH; log2BlockWidth++)
#endif
    {
      const UInt blockWidth  = 1 << log2BlockWidth;
      const UInt blockHeight = 1 << log2BlockHeight;
      const UInt totalValues = blockWidth * blockHeight;

      //--------------------------------------------------------------------------------------------------

      //non-grouped scan orders

      for (UInt scanTypeIndex = 0; scanTypeIndex < SCAN_NUMBER_OF_TYPES; scanTypeIndex++)
      {
        const COEFF_SCAN_TYPE scanType = COEFF_SCAN_TYPE(scanTypeIndex);

        g_scanOrder[SCAN_UNGROUPED][scanType][log2BlockWidth][log2BlockHeight] = new UInt[totalValues];

#if VCEG_AZ07_CTX_RESIDUALCODING
        if (scanType == SCAN_VER && log2BlockWidth == 1 && log2BlockHeight == 1) 
        {
          for (UInt scanPosition = 0; scanPosition < totalValues; scanPosition++)
          {
            g_scanOrder[SCAN_UNGROUPED][scanTypeIndex][log2BlockWidth][log2BlockHeight][scanPosition] = g_scanOrder[SCAN_UNGROUPED][scanTypeIndex - 1][log2BlockWidth][log2BlockHeight][scanPosition];
          }
          continue;
        }
#endif
        ScanGenerator fullBlockScan(blockWidth, blockHeight, blockWidth, scanType);

        for (UInt scanPosition = 0; scanPosition < totalValues; scanPosition++)
        {
          g_scanOrder[SCAN_UNGROUPED][scanType][log2BlockWidth][log2BlockHeight][scanPosition] = fullBlockScan.GetNextIndex(0, 0);
        }
      }

      //--------------------------------------------------------------------------------------------------

      //grouped scan orders
#if VCEG_AZ07_CTX_RESIDUALCODING || JVET_C0024_QTBT
      UInt  groupWidth                 = 1           << MLS_CG_LOG2_WIDTH;
      UInt  groupHeight                = 1           << MLS_CG_LOG2_HEIGHT;
      UInt  widthInGroups              = blockWidth  >> MLS_CG_LOG2_WIDTH;
      UInt  heightInGroups             = blockHeight >> MLS_CG_LOG2_HEIGHT;
#else
      const UInt  groupWidth           = 1           << MLS_CG_LOG2_WIDTH;
      const UInt  groupHeight          = 1           << MLS_CG_LOG2_HEIGHT;
      const UInt  widthInGroups        = blockWidth  >> MLS_CG_LOG2_WIDTH;
      const UInt  heightInGroups       = blockHeight >> MLS_CG_LOG2_HEIGHT;
#endif

#if JVET_C0024_QTBT
      if (log2BlockWidth==1 || log2BlockHeight==1)
      {
        groupWidth                 = 1           << (MLS_CG_LOG2_WIDTH-1);
        groupHeight                = 1           << (MLS_CG_LOG2_HEIGHT-1);
        widthInGroups              = blockWidth  >> (MLS_CG_LOG2_WIDTH-1);
        heightInGroups             = blockHeight >> (MLS_CG_LOG2_HEIGHT-1);
      }
#endif

      const UInt  groupSize            = groupWidth    * groupHeight;
      const UInt  totalGroups          = widthInGroups * heightInGroups;
      for (UInt scanTypeIndex = 0; scanTypeIndex < SCAN_NUMBER_OF_TYPES; scanTypeIndex++)
      {
        const COEFF_SCAN_TYPE scanType = COEFF_SCAN_TYPE(scanTypeIndex);

        g_scanOrder[SCAN_GROUPED_4x4][scanType][log2BlockWidth][log2BlockHeight] = new UInt[totalValues];
#if VCEG_AZ07_CTX_RESIDUALCODING
        Bool bHorVerCGScan = (scanType && log2BlockWidth == 3 && log2BlockHeight == 3) ;
        if ( bHorVerCGScan ) 
        {
          for (UInt scanPosition = 0; scanPosition < totalValues; scanPosition++)
          {
            g_scanOrder[SCAN_GROUPED_4x4][scanType][log2BlockWidth][log2BlockHeight][scanPosition] = g_scanOrder[SCAN_UNGROUPED][scanType][log2BlockWidth][log2BlockHeight][scanPosition];
          }
        }
        else
        {
#endif
        ScanGenerator fullBlockScan(widthInGroups, heightInGroups, groupWidth, scanType);

        for (UInt groupIndex = 0; groupIndex < totalGroups; groupIndex++)
        {
          const UInt groupPositionY  = fullBlockScan.GetCurrentY();
          const UInt groupPositionX  = fullBlockScan.GetCurrentX();
          const UInt groupOffsetX    = groupPositionX * groupWidth;
          const UInt groupOffsetY    = groupPositionY * groupHeight;
          const UInt groupOffsetScan = groupIndex     * groupSize;

          ScanGenerator groupScan(groupWidth, groupHeight, blockWidth, scanType);

          for (UInt scanPosition = 0; scanPosition < groupSize; scanPosition++)
          {
            g_scanOrder[SCAN_GROUPED_4x4][scanType][log2BlockWidth][log2BlockHeight][groupOffsetScan + scanPosition] = groupScan.GetNextIndex(groupOffsetX, groupOffsetY);
          }

          fullBlockScan.GetNextIndex(0,0);
        }
#if VCEG_AZ07_CTX_RESIDUALCODING
        }
#endif
      }

      //--------------------------------------------------------------------------------------------------
    }
  }

#if RRC_DWDIP
  //Initialize weights & shifts
#if COM16_C806_T64
  for(UInt log2BlockHeight = 0; log2BlockHeight < MAX_LOG2_TU_SIZE_PLUS_ONE; log2BlockHeight++)
  {
    for(UInt log2BlockWidth = 0; log2BlockWidth < MAX_LOG2_TU_SIZE_PLUS_ONE; log2BlockWidth++)
#else
  for (UInt log2BlockHeight = 0; log2BlockHeight < MAX_CU_DEPTH; log2BlockHeight++)
  {
    for (UInt log2BlockWidth = 0; log2BlockWidth < MAX_CU_DEPTH; log2BlockWidth++)
#endif
    {
      const UInt width = 1 << log2BlockWidth;
      const UInt height = 1 << log2BlockHeight;
      //const UInt totalValues = blockWidth * blockHeight;

      if (isEnabledForDWDIP(CHANNEL_TYPE_LUMA, width, height))
      {
        for (Int iModeIdx = 0, dirMode = DIA_IDX; dirMode <= VDIA_IDX; ++iModeIdx, ++dirMode)
        {
          const Int        intraPredAngleMode = (Int)dirMode - VER_IDX;
          const Int        absAngMode = abs(intraPredAngleMode);
          const Int        signAng = intraPredAngleMode < 0 ? -1 : 1;

          // Set bitshifts and scale the angle parameter to block size
#if VCEG_AZ07_INTRA_65ANG_MODES
          static const Int angTable[17] = { 0,    1,    2,    3,    5,    7,    9,   11,   13,   15,   17,   19,   21,   23,   26,   29,   32 };
          static const Int invAngTable[17] = { 0, 8192, 4096, 2731, 1638, 1170,  910,  745,  630,  546,  482,  431,  390,  356,  315,  282,  256 }; // (256 * 32) / Angle
#else
          static const Int angTable[9] = { 0,    2,    5,   9,  13,  17,  21,  26,  32 };
          static const Int invAngTable[9] = { 0, 4096, 1638, 910, 630, 482, 390, 315, 256 }; // (256 * 32) / Angle
#endif
          Int invAngle = invAngTable[absAngMode];
          Int absAng = angTable[absAngMode];
          Int intraPredAngle = signAng * absAng;

          Int * pWeights = getDWDIPweightsGen(width, height, dirMode);
          Int * pbDecision = getDWDIPSideDecisionGen(width, height, dirMode);
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
          Int * piSideDeltaPos = getDWDIPSideDeltaPosGen(width, height, dirMode);
#endif // RRC_DWDIP_DIVISION_FREE_ADDSHIFT

          for (Int y = 0, deltaPos = intraPredAngle; y < height; y++, deltaPos += intraPredAngle
            , pWeights += width
            , pbDecision += width
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
            , piSideDeltaPos += width
#endif
            )
          {
            xGetLUTsRow(y, deltaPos, width, height, pWeights, pbDecision
#if RRC_DWDIP_DIVISION_FREE_ADDSHIFT
              , piSideDeltaPos
#endif // RRC_DWDIP_DIVISION_FREE_ADDSHIFT
            );
          }
        } // for (...)
      } // if (isEnabledForDWDIP(... ))
    } // width loop
  } // height loop
#endif
}

Void destroyROM()
{
  for(UInt groupTypeIndex = 0; groupTypeIndex < SCAN_NUMBER_OF_GROUP_TYPES; groupTypeIndex++)
  {
    for (UInt scanOrderIndex = 0; scanOrderIndex < SCAN_NUMBER_OF_TYPES; scanOrderIndex++)
    {
#if COM16_C806_T64
      for (UInt log2BlockWidth = 0; log2BlockWidth < MAX_LOG2_TU_SIZE_PLUS_ONE; log2BlockWidth++)
      {
        for (UInt log2BlockHeight = 0; log2BlockHeight < MAX_LOG2_TU_SIZE_PLUS_ONE; log2BlockHeight++)
#else
      for (UInt log2BlockWidth = 0; log2BlockWidth < MAX_CU_DEPTH; log2BlockWidth++)
      {
        for (UInt log2BlockHeight = 0; log2BlockHeight < MAX_CU_DEPTH; log2BlockHeight++)
#endif
        {
          delete [] g_scanOrder[groupTypeIndex][scanOrderIndex][log2BlockWidth][log2BlockHeight];
        }
      }
    }
  }
}

#if BLOCK_BASED_BACKGROUND_REFERENCE
UInt g_uiMaxCUWidth = MAX_CU_SIZE;
UInt g_uiMaxCUHeight = MAX_CU_SIZE;
#endif

// ====================================================================================================================
// Data structure related table & variable
// ====================================================================================================================
UInt g_auiZscanToRaster [ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ] = { 0, };
UInt g_auiRasterToZscan [ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ] = { 0, };
UInt g_auiRasterToPelX  [ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ] = { 0, };
UInt g_auiRasterToPelY  [ MAX_NUM_PART_IDXS_IN_CTU_WIDTH*MAX_NUM_PART_IDXS_IN_CTU_WIDTH ] = { 0, };

#if !JVET_C0024_QTBT
const UInt g_auiPUOffset[NUMBER_OF_PART_SIZES] = { 0, 8, 4, 4, 2, 10, 1, 5};
#endif
#if VCEG_AZ07_CTX_RESIDUALCODING
const UInt g_auiGoRiceRange[MAX_GR_ORDER_RESIDUAL] =
{
  6, 5, 6, COEF_REMAIN_BIN_REDUCTION, COEF_REMAIN_BIN_REDUCTION, COEF_REMAIN_BIN_REDUCTION, COEF_REMAIN_BIN_REDUCTION, COEF_REMAIN_BIN_REDUCTION, COEF_REMAIN_BIN_REDUCTION, COEF_REMAIN_BIN_REDUCTION
};
#if !COM16_C806_T64
const UInt g_uiLastCtx[ 28 ]    =         //!!!!to be modified for when COM16_C806_T64 = 1  
{
  0,   1,  2,  2,                         // 4x4    4
  3,   4,  5,  5, 2,  2,                  // 8x8    6  
  6,   7,  8,  8, 9,  9, 2, 2,            // 16x16  8
  10, 11, 12, 12, 13, 13, 14, 14, 2, 2    // 32x32  10
                                          // 64x64  12    
};
#endif
#endif

#if PVC_AQS // prepare error sensitivity map
//initial the look-up table for error sensitivity functions
int  g_aiFactorContrast[256];  // look-up table from std to contrast masking factor
int  g_aiFactorLuminance[256];  // look-up table from avg to luminance masking factor
Int  g_ESMFactor = 256;
Int  g_ESMFactorIdx = 0;

Void initFactorValue()
{
  Int i = 0;
  Double dTemp;
  for (i = 0; i<60; i++)
  {
    dTemp = 1.0 / (1.0 + (60.0 - i) / 150.0);
    g_aiFactorLuminance[i] = (Int) (dTemp * (1 << ESM_SHIFT) + 0.5);
  }

  for (i = 60; i<170; i++)
  {
    dTemp = 1.0;
    g_aiFactorLuminance[i] = (Int) (dTemp * (1 << ESM_SHIFT) + 0.5);
  }

  for (i = 170; i<256; i++)
  {
    dTemp = 1.0 / (1.0 + (i - 170.0) / 425.0);
    g_aiFactorLuminance[i] = (Int) (dTemp * (1 << ESM_SHIFT) + 0.5);
  }

  for (i = 0; i<256; i++)
  {
    Int k = 5; // threshold of the zero-slope line
    // with pedestal effect
    if (i<k)
      dTemp = 1.0;
    else
      dTemp = 10.0 / (Double) (i - k + 10.0);

    g_aiFactorContrast[i] = (Int) (dTemp * (1 << ESM_SHIFT) + 0.5);
  }
}
#endif

Void initZscanToRaster ( Int iMaxDepth, Int iDepth, UInt uiStartVal, UInt*& rpuiCurrIdx )
{
  Int iStride = 1 << ( iMaxDepth - 1 );

  if ( iDepth == iMaxDepth )
  {
    rpuiCurrIdx[0] = uiStartVal;
    rpuiCurrIdx++;
  }
  else
  {
    Int iStep = iStride >> iDepth;
    initZscanToRaster( iMaxDepth, iDepth+1, uiStartVal,                     rpuiCurrIdx );
    initZscanToRaster( iMaxDepth, iDepth+1, uiStartVal+iStep,               rpuiCurrIdx );
    initZscanToRaster( iMaxDepth, iDepth+1, uiStartVal+iStep*iStride,       rpuiCurrIdx );
    initZscanToRaster( iMaxDepth, iDepth+1, uiStartVal+iStep*iStride+iStep, rpuiCurrIdx );
  }
}

Void initRasterToZscan ( UInt uiMaxCUWidth, UInt uiMaxCUHeight, UInt uiMaxDepth )
{
  UInt  uiMinCUWidth  = uiMaxCUWidth  >> ( uiMaxDepth - 1 );
  UInt  uiMinCUHeight = uiMaxCUHeight >> ( uiMaxDepth - 1 );

  UInt  uiNumPartInWidth  = (UInt)uiMaxCUWidth  / uiMinCUWidth;
  UInt  uiNumPartInHeight = (UInt)uiMaxCUHeight / uiMinCUHeight;

  for ( UInt i = 0; i < uiNumPartInWidth*uiNumPartInHeight; i++ )
  {
    g_auiRasterToZscan[ g_auiZscanToRaster[i] ] = i;
  }
}

Void initRasterToPelXY ( UInt uiMaxCUWidth, UInt uiMaxCUHeight, UInt uiMaxDepth )
{
  UInt    i;

  UInt* uiTempX = &g_auiRasterToPelX[0];
  UInt* uiTempY = &g_auiRasterToPelY[0];

  UInt  uiMinCUWidth  = uiMaxCUWidth  >> ( uiMaxDepth - 1 );
  UInt  uiMinCUHeight = uiMaxCUHeight >> ( uiMaxDepth - 1 );

  UInt  uiNumPartInWidth  = uiMaxCUWidth  / uiMinCUWidth;
  UInt  uiNumPartInHeight = uiMaxCUHeight / uiMinCUHeight;

  uiTempX[0] = 0; uiTempX++;
  for ( i = 1; i < uiNumPartInWidth; i++ )
  {
    uiTempX[0] = uiTempX[-1] + uiMinCUWidth; uiTempX++;
  }
  for ( i = 1; i < uiNumPartInHeight; i++ )
  {
    memcpy(uiTempX, uiTempX-uiNumPartInWidth, sizeof(UInt)*uiNumPartInWidth);
    uiTempX += uiNumPartInWidth;
  }

  for ( i = 1; i < uiNumPartInWidth*uiNumPartInHeight; i++ )
  {
    uiTempY[i] = ( i / uiNumPartInWidth ) * uiMinCUWidth;
  }
}

const Int g_quantScales[SCALING_LIST_REM_NUM] =
{
  26214,23302,20560,18396,16384,14564
};

const Int g_invQuantScales[SCALING_LIST_REM_NUM] =
{
  40,45,51,57,64,72
};

//--------------------------------------------------------------------------------------------------

//structures

#define DEFINE_DST4x4_MATRIX(a,b,c,d) \
{ \
  {  a,  b,  c,  d }, \
  {  c,  c,  0, -c }, \
  {  d, -a, -c,  b }, \
  {  b, -d,  c, -a }, \
}

#define DEFINE_DCT4x4_MATRIX(a,b,c) \
{ \
  { a,  a,  a,  a}, \
  { b,  c, -c, -b}, \
  { a, -a, -a,  a}, \
  { c, -b,  b, -c}  \
}

#define DEFINE_DCT8x8_MATRIX(a,b,c,d,e,f,g) \
{ \
  { a,  a,  a,  a,  a,  a,  a,  a}, \
  { d,  e,  f,  g, -g, -f, -e, -d}, \
  { b,  c, -c, -b, -b, -c,  c,  b}, \
  { e, -g, -d, -f,  f,  d,  g, -e}, \
  { a, -a, -a,  a,  a, -a, -a,  a}, \
  { f, -d,  g,  e, -e, -g,  d, -f}, \
  { c, -b,  b, -c, -c,  b, -b,  c}, \
  { g, -f,  e, -d,  d, -e,  f, -g}  \
}

#define DEFINE_DCT16x16_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
{ \
  { a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a}, \
  { h,  i,  j,  k,  l,  m,  n,  o, -o, -n, -m, -l, -k, -j, -i, -h}, \
  { d,  e,  f,  g, -g, -f, -e, -d, -d, -e, -f, -g,  g,  f,  e,  d}, \
  { i,  l,  o, -m, -j, -h, -k, -n,  n,  k,  h,  j,  m, -o, -l, -i}, \
  { b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b}, \
  { j,  o, -k, -i, -n,  l,  h,  m, -m, -h, -l,  n,  i,  k, -o, -j}, \
  { e, -g, -d, -f,  f,  d,  g, -e, -e,  g,  d,  f, -f, -d, -g,  e}, \
  { k, -m, -i,  o,  h,  n, -j, -l,  l,  j, -n, -h, -o,  i,  m, -k}, \
  { a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a}, \
  { l, -j, -n,  h, -o, -i,  m,  k, -k, -m,  i,  o, -h,  n,  j, -l}, \
  { f, -d,  g,  e, -e, -g,  d, -f, -f,  d, -g, -e,  e,  g, -d,  f}, \
  { m, -h,  l,  n, -i,  k,  o, -j,  j, -o, -k,  i, -n, -l,  h, -m}, \
  { c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c}, \
  { n, -k,  h, -j,  m,  o, -l,  i, -i,  l, -o, -m,  j, -h,  k, -n}, \
  { g, -f,  e, -d,  d, -e,  f, -g, -g,  f, -e,  d, -d,  e, -f,  g}, \
  { o, -n,  m, -l,  k, -j,  i, -h,  h, -i,  j, -k,  l, -m,  n, -o}  \
}

#define DEFINE_DCT32x32_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E) \
{ \
  { a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a}, \
  { p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,  A,  B,  C,  D,  E, -E, -D, -C, -B, -A, -z, -y, -x, -w, -v, -u, -t, -s, -r, -q, -p}, \
  { h,  i,  j,  k,  l,  m,  n,  o, -o, -n, -m, -l, -k, -j, -i, -h, -h, -i, -j, -k, -l, -m, -n, -o,  o,  n,  m,  l,  k,  j,  i,  h}, \
  { q,  t,  w,  z,  C, -E, -B, -y, -v, -s, -p, -r, -u, -x, -A, -D,  D,  A,  x,  u,  r,  p,  s,  v,  y,  B,  E, -C, -z, -w, -t, -q}, \
  { d,  e,  f,  g, -g, -f, -e, -d, -d, -e, -f, -g,  g,  f,  e,  d,  d,  e,  f,  g, -g, -f, -e, -d, -d, -e, -f, -g,  g,  f,  e,  d}, \
  { r,  w,  B, -D, -y, -t, -p, -u, -z, -E,  A,  v,  q,  s,  x,  C, -C, -x, -s, -q, -v, -A,  E,  z,  u,  p,  t,  y,  D, -B, -w, -r}, \
  { i,  l,  o, -m, -j, -h, -k, -n,  n,  k,  h,  j,  m, -o, -l, -i, -i, -l, -o,  m,  j,  h,  k,  n, -n, -k, -h, -j, -m,  o,  l,  i}, \
  { s,  z, -D, -w, -p, -v, -C,  A,  t,  r,  y, -E, -x, -q, -u, -B,  B,  u,  q,  x,  E, -y, -r, -t, -A,  C,  v,  p,  w,  D, -z, -s}, \
  { b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b}, \
  { t,  C, -y, -p, -x,  D,  u,  s,  B, -z, -q, -w,  E,  v,  r,  A, -A, -r, -v, -E,  w,  q,  z, -B, -s, -u, -D,  x,  p,  y, -C, -t}, \
  { j,  o, -k, -i, -n,  l,  h,  m, -m, -h, -l,  n,  i,  k, -o, -j, -j, -o,  k,  i,  n, -l, -h, -m,  m,  h,  l, -n, -i, -k,  o,  j}, \
  { u, -E, -t, -v,  D,  s,  w, -C, -r, -x,  B,  q,  y, -A, -p, -z,  z,  p,  A, -y, -q, -B,  x,  r,  C, -w, -s, -D,  v,  t,  E, -u}, \
  { e, -g, -d, -f,  f,  d,  g, -e, -e,  g,  d,  f, -f, -d, -g,  e,  e, -g, -d, -f,  f,  d,  g, -e, -e,  g,  d,  f, -f, -d, -g,  e}, \
  { v, -B, -p, -C,  u,  w, -A, -q, -D,  t,  x, -z, -r, -E,  s,  y, -y, -s,  E,  r,  z, -x, -t,  D,  q,  A, -w, -u,  C,  p,  B, -v}, \
  { k, -m, -i,  o,  h,  n, -j, -l,  l,  j, -n, -h, -o,  i,  m, -k, -k,  m,  i, -o, -h, -n,  j,  l, -l, -j,  n,  h,  o, -i, -m,  k}, \
  { w, -y, -u,  A,  s, -C, -q,  E,  p,  D, -r, -B,  t,  z, -v, -x,  x,  v, -z, -t,  B,  r, -D, -p, -E,  q,  C, -s, -A,  u,  y, -w}, \
  { a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a}, \
  { x, -v, -z,  t,  B, -r, -D,  p, -E, -q,  C,  s, -A, -u,  y,  w, -w, -y,  u,  A, -s, -C,  q,  E, -p,  D,  r, -B, -t,  z,  v, -x}, \
  { l, -j, -n,  h, -o, -i,  m,  k, -k, -m,  i,  o, -h,  n,  j, -l, -l,  j,  n, -h,  o,  i, -m, -k,  k,  m, -i, -o,  h, -n, -j,  l}, \
  { y, -s, -E,  r, -z, -x,  t,  D, -q,  A,  w, -u, -C,  p, -B, -v,  v,  B, -p,  C,  u, -w, -A,  q, -D, -t,  x,  z, -r,  E,  s, -y}, \
  { f, -d,  g,  e, -e, -g,  d, -f, -f,  d, -g, -e,  e,  g, -d,  f,  f, -d,  g,  e, -e, -g,  d, -f, -f,  d, -g, -e,  e,  g, -d,  f}, \
  { z, -p,  A,  y, -q,  B,  x, -r,  C,  w, -s,  D,  v, -t,  E,  u, -u, -E,  t, -v, -D,  s, -w, -C,  r, -x, -B,  q, -y, -A,  p, -z}, \
  { m, -h,  l,  n, -i,  k,  o, -j,  j, -o, -k,  i, -n, -l,  h, -m, -m,  h, -l, -n,  i, -k, -o,  j, -j,  o,  k, -i,  n,  l, -h,  m}, \
  { A, -r,  v, -E, -w,  q, -z, -B,  s, -u,  D,  x, -p,  y,  C, -t,  t, -C, -y,  p, -x, -D,  u, -s,  B,  z, -q,  w,  E, -v,  r, -A}, \
  { c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c}, \
  { B, -u,  q, -x,  E,  y, -r,  t, -A, -C,  v, -p,  w, -D, -z,  s, -s,  z,  D, -w,  p, -v,  C,  A, -t,  r, -y, -E,  x, -q,  u, -B}, \
  { n, -k,  h, -j,  m,  o, -l,  i, -i,  l, -o, -m,  j, -h,  k, -n, -n,  k, -h,  j, -m, -o,  l, -i,  i, -l,  o,  m, -j,  h, -k,  n}, \
  { C, -x,  s, -q,  v, -A, -E,  z, -u,  p, -t,  y, -D, -B,  w, -r,  r, -w,  B,  D, -y,  t, -p,  u, -z,  E,  A, -v,  q, -s,  x, -C}, \
  { g, -f,  e, -d,  d, -e,  f, -g, -g,  f, -e,  d, -d,  e, -f,  g,  g, -f,  e, -d,  d, -e,  f, -g, -g,  f, -e,  d, -d,  e, -f,  g}, \
  { D, -A,  x, -u,  r, -p,  s, -v,  y, -B,  E,  C, -z,  w, -t,  q, -q,  t, -w,  z, -C, -E,  B, -y,  v, -s,  p, -r,  u, -x,  A, -D}, \
  { o, -n,  m, -l,  k, -j,  i, -h,  h, -i,  j, -k,  l, -m,  n, -o, -o,  n, -m,  l, -k,  j, -i,  h, -h,  i, -j,  k, -l,  m, -n,  o}, \
  { E, -D,  C, -B,  A, -z,  y, -x,  w, -v,  u, -t,  s, -r,  q, -p,  p, -q,  r, -s,  t, -u,  v, -w,  x, -y,  z, -A,  B, -C,  D, -E}  \
}

#if COM16_C806_EMT
#if FAST_EMT_v1
Int g_aiTrSubsetIntra[3] = { DCT2,  DCT2,  DCT2 };
#else
Int g_aiTrSubsetIntra[3][2] = { {DST7, DCT8}, {DST7, DST1}, {DST7, DCT5} };
#endif
#if FAST_EMT_v1 //EMT_SIGNAL_INTER
Int g_aiTrSubsetInter[4] = { DST7,DCT2 };
#else
Int g_aiTrSubsetInter   [4] =   {DCT8, DST7};
#endif
#if VCEG_AZ07_INTRA_65ANG_MODES
const UChar g_aucTrSetVert[NUM_INTRA_MODE-1] =
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66
   2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0
};
const UChar g_aucTrSetHorz[NUM_INTRA_MODE-1] =
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66
   2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0
};
#else
const UChar g_aucTrSetVert[35] =
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
   2, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2, 2, 1, 0, 1, 0, 1, 0
};
const UChar g_aucTrSetHorz[35] =
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
   2, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2, 2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0
};
#endif
#if FAST_EMT_v1
const UInt g_iEmtSigNumThr = 2;
#else
const UInt g_iEmtSigNumThr = 2;
#endif
#endif

#if COM16_C806_EMT || COM16_C806_T64
#if JVET_C0024_QTBT
TMatrixCoeff g_aiTr2 [NUM_TRANS_TYPE][ 2][ 2];
TMatrixCoeff g_aiTr128 [NUM_TRANS_TYPE][ 128][ 128];
#endif
TMatrixCoeff g_aiTr4 [NUM_TRANS_TYPE][ 4][ 4];
TMatrixCoeff g_aiTr8 [NUM_TRANS_TYPE][ 8][ 8];
TMatrixCoeff g_aiTr16[NUM_TRANS_TYPE][16][16];
TMatrixCoeff g_aiTr32[NUM_TRANS_TYPE][32][32];
TMatrixCoeff g_aiTr64[NUM_TRANS_TYPE][64][64];
#endif

#if COM16_C1044_NSST
#if JVET_D0120_NSST_IMPROV
#if VCEG_AZ07_INTRA_65ANG_MODES
const UChar g_NsstLut[NUM_INTRA_MODE-1] = 
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2
};
#else
const UChar g_NsstLut[NUM_INTRA_MODE-1] = 
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
   0, 1, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28,30,32,34,32,30,28,26,24,22,20,18,16,14,12,10, 8, 6, 4, 2
};
#endif
#else
#if VCEG_AZ07_INTRA_65ANG_MODES
const UChar g_NsstLut[NUM_INTRA_MODE-1] = 
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66
   0, 0, 1, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 8, 9, 8, 9,10,11,10,11,10,11,10,11,10,11,10, 9, 8, 9, 8, 9, 8, 7, 7, 7, 6, 6, 6, 5, 5, 5, 4, 3, 4, 3, 4, 3, 2, 1, 2, 1, 2, 1
};
#else
const UChar g_NsstLut[NUM_INTRA_MODE-1] = 
{//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
   0, 0, 1, 1, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 9,10,10,11,11,11,10,10, 9, 8, 7, 7, 6, 5, 5, 4, 3, 2, 2, 1, 1
};
#endif
#endif

#if VCEG_AZ07_CTX_RESIDUALCODING
const UInt g_auiCoefScanFirstCG8x8[ 3 ][ 16 ] = {
  { 0, 8,  1, 16, 9, 2, 24, 17, 10,  3, 25, 18, 11, 26, 19, 27 },
  { 0, 1,  2,  3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27 },
  { 0, 8, 16, 24, 1, 9, 17, 25,  2, 10, 18, 26,  3, 11, 19, 27 },
};
#endif

#if JVET_D0120_NSST_IMPROV
#if JVET_C0024_QTBT
const UInt g_auiCoefTopLeftDiagScan8x8[5][64] = 
{
  {0,  8,  1, 16,  9,  2, 24, 17, 10,  3, 25, 18, 11, 26, 19, 27,    32, 40, 33, 48, 41, 34, 56, 49, 42, 35, 57, 50, 43, 58, 51, 59,     4, 12,  5, 20, 13,  6, 28, 21, 14,  7, 29, 22, 15, 30, 23, 31,    36, 44, 37, 52, 45, 38, 60, 53, 46, 39, 61, 54, 47, 62, 55, 63}, //   8
  {0, 16,  1, 32, 17,  2, 48, 33, 18,  3, 49, 34, 19, 50, 35, 51,    64, 80, 65, 96, 81, 66,112, 97, 82, 67,113, 98, 83,114, 99,115,     4, 20,  5, 36, 21,  6, 52, 37, 22,  7, 53, 38, 23, 54, 39, 55,    68, 84, 69,100, 85, 70,116,101, 86, 71,117,102, 87,118,103,119}, //  16
  {0, 32,  1, 64, 33,  2, 96, 65, 34,  3, 97, 66, 35, 98, 67, 99,   128,160,129,192,161,130,224,193,162,131,225,194,163,226,195,227,     4, 36,  5, 68, 37,  6,100, 69, 38,  7,101, 70, 39,102, 71,103,   132,164,133,196,165,134,228,197,166,135,229,198,167,230,199,231}, //  32
  {0, 64,  1,128, 65,  2,192,129, 66,  3,193,130, 67,194,131,195,   256,320,257,384,321,258,448,385,322,259,449,386,323,450,387,451,     4, 68,  5,132, 69,  6,196,133, 70,  7,197,134, 71,198,135,199,   260,324,261,388,325,262,452,389,326,263,453,390,327,454,391,455}, //  64
  {0,128,  1,256,129,  2,384,257,130,  3,385,258,131,386,259,387,   512,640,513,768,641,514,896,769,642,515,897,770,643,898,771,899,     4,132,  5,260,133,  6,388,261,134,  7,389,262,135,390,263,391,   516,644,517,772,645,518,900,773,646,519,901,774,647,902,775,903}, // 128
};
#else
const UInt g_auiCoefTopLeftDiagScan8x8[3][64] = 
{
  {0, 16,  1, 32, 17,  2, 48, 33, 18,  3, 49, 34, 19, 50, 35, 51,    64, 80, 65, 96, 81, 66,112, 97, 82, 67,113, 98, 83,114, 99,115,     4, 20,  5, 36, 21,  6, 52, 37, 22,  7, 53, 38, 23, 54, 39, 55,    68, 84, 69,100, 85, 70,116,101, 86, 71,117,102, 87,118,103,119}, // 16x16
  {0, 32,  1, 64, 33,  2, 96, 65, 34,  3, 97, 66, 35, 98, 67, 99,   128,160,129,192,161,130,224,193,162,131,225,194,163,226,195,227,     4, 36,  5, 68, 37,  6,100, 69, 38,  7,101, 70, 39,102, 71,103,   132,164,133,196,165,134,228,197,166,135,229,198,167,230,199,231}, // 32x32
  {0, 64,  1,128, 65,  2,192,129, 66,  3,193,130, 67,194,131,195,   256,320,257,384,321,258,448,385,322,259,449,386,323,450,387,451,     4, 68,  5,132, 69,  6,196,133, 70,  7,197,134, 71,198,135,199,   260,324,261,388,325,262,452,389,326,263,453,390,327,454,391,455}, // 64x64
};
#endif
#endif

#if JVET_D0120_NSST_IMPROV
tabSinCos g_tabSinCos[NSST_HYGT_PTS]; 

const int g_nsstHyGTPar4x4[35][3][64] = 
{
  {  //   0
    {
      233,124, 86, 87, 81, 89, 85, 87,247,254,  9, 10,141, 11,139, 13, 49, 72,255, 51, 64,197, 74,  0,251,251,  0,251,245, 25,254,  0,
       38,129, 88, 58, 33, 59,166,  7,  2,  2,252,128,  3,128,  2,  1,241,177,  8,107,125,197,152, 88,  0, 14,139,  9,168,195,116,246,
    },
    {
      212,205, 70,200,220,199,177, 45,  7,  7,245,247,  4, 26,116,118, 17,  9,220,188, 77,  2,122,135, 65,196,199,216, 24,  3,251,  3,
       59, 56, 54, 58,150, 32,247,205, 10,126,  7,  4,251,254,129,254,187, 75, 49,222,  7,255,130,123, 39, 69,237, 44, 98, 62,152,240,
    },
    {
      254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,
      254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,254,252,
    },
  },
  {  //   1
    {
       12, 45, 16, 42, 11, 39, 14, 50,  5, 12,  8, 11, 10,142, 10, 20, 52, 38, 46, 46, 81,251,201, 43,  8, 13,246, 11,  5,  2,245,250,
      247, 30,  0, 33,254,229,  4,165,252,  3,253,  8,  6,132,251,232,210,222, 87,151,249,  7,134,199,251,206,  6,197,  8,  6,  3,  6,
    },
    {
        4,255,198,  0, 12,  1, 74, 63, 59,  6,253,  0,  8,194,136,202,238, 42, 63,241,197,163,133,246,  3,  1,194,128,  9,246,  0, 77,
      199, 52,195,193, 71,224, 59,132,253,254,253,183,195,193,  0,252,115,198,112, 76, 14, 53,129, 92,241,254, 60,248, 10,183,250,254,
    },
    {
      254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,
      254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,254,  3,
    },
  },
  {  //   2
    {
        7, 73,212, 27, 46,121,221,239,255, 25, 31,  3,255, 30, 40, 43,210,221, 53, 71,216, 37,  3, 47,196, 57, 30,  0,253,248, 31,211,
       12,198, 11, 80,248,172, 15, 57,215, 49,  2,  0,253,  8, 30, 10, 66,240,212, 15, 30,254,242, 14, 13,  3,202,224, 49,  4,246,225,
    },
    {
      147, 17, 25, 12,219,231,205, 42,255, 63,  0,181,210,246,189,126,113,215,217,204, 92, 53, 10, 11,185, 32,102,237, 68,132, 66, 62,
       66,118,117,152,193, 70,242,215, 88,161,219, 68, 20, 57,220,226,244,192, 81,251, 66, 23, 11,254,111,132,211,  9, 58,177, 61, 58,
    },
    {
      197,197,195,193,222,225, 40,183, 67, 69, 67, 76, 79, 80, 59,167,  1,212,171,246, 75, 84,122, 70,186,202,215,191, 21, 34,191, 47,
       22, 26, 52,174, 59,  7, 38,238,255,210, 12,186,108, 63,203, 61,245,234,173,114, 70,  1, 72,220, 64,193,126, 42,175,  3, 72, 21,
    },
  },
  {  //   3
    {
       26,202,196, 21,233, 15,191, 66,192,179,213, 66,158,  9,230,208,199, 14,203,222,244, 23,237,  6, 40, 41,  2,  4,127, 18,255,  2,
      123, 46,  7, 46,135, 57, 18,214,212,151,196, 60,  9, 58,207,222,189,  7,126,246,225, 33, 25,127, 12,101,210, 31,253, 65,186,230,
    },
    {
      127, 17,195, 36,255, 12, 67,212,195,186,200, 64,188, 79, 99, 84,243,239,245,163,244, 17, 72,252,219, 38, 73,192,239,224,240,219,
      118,  6,228,255, 17,  6,104, 25,  3, 12,187,181,  6,  7,251,204,242, 79,129,230,202,195,254,  3, 49,155,200,110,149, 13, 73,209,
    },
    {
      199, 23, 13, 91,100,168, 22, 95,149,128,  2,157, 28,  1,  4,160, 25,118,  0,200,146,102,155,225,255,  1, 42,218, 40, 72,134,130,
       56,245,193,252, 49,  5,198,119,190, 68, 23,245,228, 21, 66,128,178,122, 66, 70,136,250,158,184, 16,216,211,191,122,148, 59,208,
    },
  },
  {  //   4
    {
      253,245,  4,252, 17,254,254,243,  3,  3, 16,  4, 12,  4,144,  7, 58, 34, 36,103,254,200, 20,175, 18,210,238, 15,255,  8, 58,111,
      250,249,196,  9, 37,  8, 26,  4,248,  3,  0,  0, 67,254,193,  0,213,192,173, 80,115, 39, 43,211, 10, 59,113, 97,  2,232,181,213,
    },
    {
        0, 55,223,198,  8,244,  3,166,  2,  0,  5, 19, 28,  4,244,207,  8,225,149,227,240, 51, 85, 70,  4,255,  3,252, 72,174, 58,244,
        2,245,243, 29,248,205, 23, 17,  3,239,195,206, 17,206,  2,255,255, 43,240, 12,226, 63, 57,214,  2, 64,110,174,223,  3,151,243,
    },
    {
       20, 84, 11, 64, 17, 63, 13, 75,  0, 15,  0, 18,  1, 19,  0, 22,193, 78, 45, 74,212, 35,255,  5,226,201, 17, 11, 39, 35, 47,240,
      254,205,  6,189,  5,205,  5,197,251,  3,  8,  4,252,  3,253,  9, 38,249,209,214, 81, 44,222, 12,216,  6,245,218, 10,199,  5,219,
    },
  },
  {  //   5
    {
       77, 12,212,  5,190, 45,210,236,228,  0, 19,  1,254,244,227,164,233, 55,218,230,199, 18,219,147, 64, 48,223,255,240,  0, 17, 66,
       70,  3,247, 12, 28,182, 31, 17, 14,  9,221,  0,254,254,224, 27,221, 64,  3,  5,244,249, 26,220, 68,252, 44, 38,224, 62, 39, 36,
    },
    {
       72,248,186, 10,205, 21,158,  1, 25,  1, 21,  0, 39, 10,129, 30,248,245,237, 12, 42,198,172,219, 60, 18,191,240, 10,242, 61,254,
      191,154, 23,250,191, 62,182, 27,255,207, 31,  0,  1,223,  8, 26,  8,  4,107,244,  6, 41, 85, 29,237,253, 14,108, 34, 32,197, 15,
    },
    {
       27,144,156,143,215, 91,160,  3,202,  1,182,251,111,188,206,175,185, 67,201,179,242,  5,  0,  1,248,122, 64, 64,193, 65, 78,193,
       88,132,167,191,105,195,162, 66, 70,227, 60,191, 15, 62, 18,192,216, 59,  2,  4,254,184,  0,130,174, 51,108, 83,  1,191,192,130,
    },
  },
  {  //   6
    {
       68, 27,235,248,201,110, 92,253,  6,  5,  7,  4,240,253, 67,252,247, 68,174,199,241, 20,223,151,229,193,  7,212,  5, 36, 33, 64,
      221, 26,252,  4,156,247,  0, 95,123,  0,252,  0,158,  1,  3,252, 14,  5,223, 68,174,  1,204, 97,195,169, 27,205, 12,220, 50, 35,
    },
    {
      193,133, 72,193,152,243,195, 55,125,  2,127, 69,231, 34, 55, 73, 67,  0, 41,255,  2, 88,173,182,149,232,192,192,129,235,139,251,
      198,192,141,190, 51,211,196,180,  2,236,173, 49,  3,192, 33,203,117,195,187,111,248,125, 58,  3, 24,181,130, 16,149,219,127,249,
    },
    {
      255,  2,255,  2,  0,  3,  2,  4, 65,  2, 66,  4,198,  8, 77,  8,213,149,235,223, 33,209, 94,234,141,219,252,251, 82,181,215, 55,
       64, 65, 65, 65, 10,193,  1,192, 63,189,196,192,196,235,194,151,  7, 65, 41, 11,245,245, 47,202,242,  1,232,250, 28,190,194,  7,
    },
  },
  {  //   7
    {
       15,205,198,  3, 11,227,169,182,129,122,  6,  0,253,246,202,209,  6, 19, 65, 44,196, 77,159, 33,192, 61,210,199,228, 49,188,178,
        6, 93,255,221, 12,196,249, 61, 59,104,121,255,  0, 72,  1,255, 59, 92,183, 32,246,  8, 36,252,248,238,127,201,237,185, 17, 52,
    },
    {
       78,129,174, 89, 36, 15, 23,133, 10,253,253,249,198,  3, 84,251,209,244,212,117,156,206,108,171,250,  0, 64, 65,241,161,174, 10,
      201, 63, 42, 90,145, 33,209,250,  0,175,  0, 32, 66, 94, 38,  5,247,146,127, 14,255, 58,229,223,253, 87,253,  0,206,251,105,183,
    },
    {
       77,201,255,135, 68, 68,  7,255,246, 41, 91,245, 67,130,  7,181, 19,150,129, 24,  2, 59, 45,163,  1,  1,131,250, 14,225, 25,123,
       79,194,200, 60, 60, 65, 55, 66, 29,120,125, 97,250,119, 59,186,119, 23, 85,  3,149, 57,241,137,254,203,194,229,172, 50,197, 79,
    },
  },
  {  //   8
    {
      144, 72,141, 68,155, 98,242,141, 63,106, 63, 96,187,252, 34,253, 87, 76,202, 19,121,108,  0, 40, 96, 48,116,125,251,255,130, 26,
      181,252, 63,192,  6, 64, 99, 59,141,201, 68,128, 84, 38,125,104,129,112, 18,  0,176,139,107,180, 10,203,  1,251,141, 65, 52,248,
    },
    {
      131, 77, 22, 41,205,211,  6,204, 64,239, 64,  8, 70,237, 61, 33,211, 60,248,  3,236,232, 51,211, 67, 34,  3,146, 42,228,233,204,
      193,  2, 64,238,182, 21, 68,  3, 57,196, 93,189, 10,231,218,215,182, 94,219,209, 18,228,194, 43, 20, 67,226, 31, 76, 34,229,255,
    },
    {
       60,130,137,  5,201,251, 10,  7,  5,254,249,252, 24,255, 29,136,210,252,239,238, 42, 91,145,137,188,197,171,178,216,224,143,202,
      195, 10, 85,  2,212,203,138,201,159, 62,182,  1,  1,  5,127,200, 73,193, 20, 68, 71,  1,137, 43,176,142,117,255, 68, 75,254,236,
    },
  },
  {  //   9
    {
       52,193,198,196,218,192,235,192,124,  7,  4, 13,  6, 27, 10, 23, 27, 28, 90, 89,169,168, 23,235,191,181,176,213, 10, 22, 19,253,
       13,251, 31,130,250,255, 10,216,  8,112,253, 61,251,  2,  7,  3, 58,228,208,115, 53, 10,  5, 69,199,244,122,147, 98, 53, 42,  9,
    },
    {
       95,223, 33,203, 64,  6, 18,205,125,  5,255,242, 44,  2, 60,  5, 21, 13, 44,220,202, 17,214,248,  0,251,199,248, 59,249,187,159,
       88, 46,210,246,  1, 34,196,232, 59,246, 95,221,244,254,  4,255, 11, 13,250,200,  4, 20,251,187,130, 17,  1,218,111, 58,247,122,
    },
    {
        4, 14,  8,189,202,195,231,201, 38,  0,255,229,  0, 19, 66, 14,  1, 48,  0, 64, 63, 27,194,110,  0, 33,  5,174, 15,238,109,  3,
      251,254,246,229,192,215, 54, 24,217,217,255,  1,232, 63,255, 12, 15,238,244,141,229, 25,  5, 41,252, 41, 52, 39,250, 26,127,142,
    },
  },
  {  //  10
    {
        6, 69, 21,207, 93, 94, 61,  6,207, 59, 47,194,221, 59,238,199, 23,201,236, 46, 24,182,153,  7,242, 64,  9,218,205, 67, 71,196,
        0,118, 33,  9,239, 21,243,140,207, 61, 40,223,214,  1,209, 52,255, 42,250,  7,236, 54, 79, 62, 17,  4,178,254, 51,245,218, 21,
    },
    {
       20,204,250,202,225, 23,211,  9,128, 43,  2, 49, 11, 67,  6,194,170, 22,136,140,153,242, 30, 32, 19,255, 68, 35,185,  9, 46,174,
      234, 67, 22,198,  5, 66,  1,195, 76, 81, 62,223, 82, 40,202,195,239,200,185,246, 88, 98, 55,  9, 96,247,202,138, 94, 62, 44,249,
    },
    {
       24, 81, 55,211, 41,205, 88,202,  0, 22,249,251,122,242,123,120,254,143,245,135,180,207, 70,115,  0,250,253,254,226,124, 14,181,
      203,126,148,185,228, 55, 27, 56,251, 70,247,196, 32,193,212, 59,  2,132, 67,185,182,195,247,189,250,254,232, 87,247,251,192,189,
    },
  },
  {  //  11
    {
      191, 52,182,201,246,203,249, 46,253,254,130,251,255, 15,126,  9, 35,254,225,244,253,172,234, 87,215, 61, 73,245, 32,235,203, 76,
      249,  5,244, 49,250, 41, 14,247,254, 15,189,  1,  4,  1,  3,243,244,  6,133,246, 77,131, 64,188,242, 62, 76,235, 60, 95,199, 58,
    },
    {
       42,  1, 18, 31, 66, 13, 12, 20,241,  1,252,  8,252,  0,  7,  3, 36, 16, 14,144,253, 18,197,227,  4, 20,253,  3, 85,141,156,231,
      224,254,109,186,187,176,  0,205, 12,  1,  1, 80,  1, 49,  0, 63,249,193,244, 60,235,254, 51,241,239,222,242, 28, 16,  9,113,249,
    },
    {
       72,139,197,  7,219, 13, 13, 11,  3,  0,249,255, 11,251, 70,251,211, 18,230,101,233,239,190,236, 17,  1, 59,125,  0,196,140,185,
      124,128,196, 62, 20,127,  1,244,168,  0,181, 23,250,  7,236,255,  3, 73, 40,190, 11, 49,228, 89, 13,179,234, 79, 10, 57, 22, 78,
    },
  },
  {  //  12
    {
        9,195,211,129,  5, 66,205, 68,123,232,242,255,254, 15,122,  6,250, 41,  0,  5, 43, 30,183,221,149,159,171,175, 27, 36,248,  1,
      109,239, 25, 80,238,242,200, 11,189,251,194,248,201,  7, 64, 17,230,121,250, 89,202, 31, 15, 13, 86, 85, 23, 29, 33, 23,226, 71,
    },
    {
       64, 70, 73, 74, 60,202, 59,195,191,  0,194,  0, 61,254,187,255,241,220, 33, 11,204, 84, 86,119,182, 53, 36,175,155, 12,  0,233,
       62,245, 64,  7,191,249, 66,243, 17,195, 71, 62,224, 59, 19,197, 29, 70, 59,229,225, 37,204,236,184, 31,  0,171, 14,  5, 76,240,
    },
    {
      241, 46,214, 15, 71, 15,177,  3,255,  9, 21,  1,171,  1,156,251, 74, 60,182, 69, 68, 10,249, 14,213,193,114, 39, 63, 66, 59, 64,
       53,  9, 48,249, 62,252, 32,218,125,  2,147,  6,117,254,  4, 23, 82, 85,227,  1,203, 60,  0,195,  0,190,156, 67,182,231,253, 87,
    },
  },
  {  //  13
    {
       71,  4, 13,201, 23,186,233,209,241,251,  0,237,255,247,  4,143, 53,  6,134,  6, 11, 34, 41,180,250, 64, 65,191,  3,183,  8,248,
        7,  4, 52,252,189,199, 11,107,  1,243,255, 10,250,193,  2,204, 16,247, 13,145,  2,179, 12,240,205,175,184,250,207,188,212,204,
    },
    {
      194,164,194,  6,244,128,196, 71,  0,254,143,  1,253,255,254,249, 37, 86, 13, 52,142, 48,199,120, 30,167,  1,163, 59,  1,218,131,
      158,  6, 65,202,235,  4,188, 43,255,194,  3,114,  1, 63,242,186, 64,  0,222, 72, 14,108,  1,207, 56,194, 67,239,238,187,202,102,
    },
    {
       41,  0, 61,  1,252,  3,147,  6,255,  0,252,  0, 10,  1,236,254, 95, 28, 26, 25, 18,232, 95,166,227, 51, 66,235, 79, 58,133, 45,
      202,  0,254,251, 10,255, 69,252,218,  0, 17,  4,255,  0,  0,  4,234, 15,107,250,110,174,  7,245, 49, 69,249,174,249, 54,211,229,
    },
  },
  {  //  14
    {
      201,200,211,217,215,214,167,226,186,  6,202,252,200,250,188,219,  7, 11,198,192,  2,253, 20,235,199,200, 29, 15,187,193,189,195,
      253,221,183,230,253,214,254,222,195,224,175, 16,115,  2,184, 54, 14,191,244,  9, 24,171, 56,248,199,  0,238,181,235,250,241,  3,
    },
    {
      120,252,121, 57,232,254,218, 62,254,255,  0,  7,  3,255,254,246, 15,140,188,127, 70, 42,187,146,255,147,138,138,237, 16, 56, 16,
        8,183, 54, 29, 88,126, 17,  1,255,254,244,  0,  1,255,  1,248,195, 71, 60, 60,176,102,192,172, 49,194, 41, 72,254,  8, 64,120,
    },
    {
      204,139,204,139,205,  9,202,  2, 22,194, 24,195,101,189,221,189,225,226,223,217,157, 99, 27,220,210,208,216,206, 13,244,  4,249,
      120,  4,248,  3,187,190,201,184,227,  4,157,254,253,146,253,147,248,201,  8,252,236, 81,227,151, 11, 14,  1,103,253,131, 27,251,
    },
  },
  {  //  15
    {
      216, 73,194, 71, 74, 72, 76,199,252, 68,125, 61,132,194,124, 63,122,135,125,240,197, 47,234,112,205,206, 78, 39,151,215,208,184,
      179,129,180,129, 48,126, 84,  0,  4,137,  7,154,119, 15,138,253,225, 82,181, 54, 57,203,125, 91, 67, 93,135, 20, 13,206,200,176,
    },
    {
      124,  4,136,141,  2,  9,  1,  4,254,193,254,194,  1, 65,  2, 66,242,166,215, 53, 98, 29,147, 41, 66,220,193, 41,145,177,140, 14,
      129, 51,128, 53,  0, 75,127, 56, 17,255,243,129,251,  8,252,130,246,107,156,211,222,179,173,247, 51,120, 69,140, 32,211, 22,187,
    },
    {
      117,254,113,  2,110,  2,217,126,254,  0,251,  0,122,255,249,254,173,158, 27, 92,255,156, 56,125,156, 49,225,115,  4,247, 17,235,
       49,195,208,  2,111,  5,229,122, 65,254,202,  0, 68,255, 67,  2,139,122, 53, 64, 20,229, 56, 91,190,124,  3,134,  5,130, 51,152,
    },
  },
  {  //  16
    {
      207,254,207,254, 82,253, 81,251, 13,  0, 14,  0,244,255,249,254,232,227, 18,144, 23, 20, 46, 37, 19, 27, 49,208,251,253, 20,241,
      242,240,246,245,181, 11,238,139,197,  0,190,  0,196,  3,194,  0,171, 51,241,176,  4,  1,191,  4, 48,  3, 39,144, 73,244, 23, 58,
    },
    {
      243,196,253,197,  6,120,244,242, 11,  6,  1,  9,255,253,  2,134, 62, 23, 64, 27, 33, 38,146,207,225, 56,  5, 66,  1,  2,254,  0,
        5,245,  7,132,236,185,238, 70,  0,  2, 10, 10,  2,  2,  2,251, 19,212,103,  1, 72,191,192, 95,209,221, 88,117,210,196,206,212,
    },
    {
      209,200,205,218,206,219,209,218,187,254,204,  2,204,  1,202,253,186,203,215,207,213,203,215,204,199,196,206,202,  5,  1,254,  6,
      254,217,255,220,255,223,255,214,182,  7, 21,249,  4, 51,249, 63,205,  1,  0,190, 16,217,225, 16,226,150,  1,176, 13,251, 56,  6,
    },
  },
  {  //  17
    {
      213, 21,221, 29,204, 17,212,183,165,  4,194,  6,157,  2,174,249,251,  2, 14, 17,  0, 10, 19,237,244,246,  2,  4,252,151, 69,187,
      243, 51, 12, 44,236,113,254, 10, 65,254, 66,254, 58,  0,254,255,192,253,  9,  1,196,  9,241, 59,217, 14, 16,253,  5,110,224,196,
    },
    {
       14, 57, 16, 52, 22, 53,  9, 59, 64,204, 65,202, 66,203,192,196, 58,189, 22, 24,248, 16,  9,249,199,190,252,250,163, 38,218,162,
       43,203,176,243, 84, 22, 50, 13,194, 64,189,196,194, 63, 60,189,  0,  0,206,106, 35, 33,227,219,235,254,216,185,244,248,114,217,
    },
    {
       22, 19,244,241,208, 71,235,111,  5,  8,  6,132,124,250,120,123, 66, 67,252,  6,224,244,157, 50,  8,236,  1,127,255,255,123,246,
      242, 85, 44, 13, 16,143,149,140,  1, 20,  3,  0,254,255,133,132,235,  1,  7,122,226, 84, 96,211,  1,192, 55, 52,247,248,  5,132,
    },
  },
  {  //  18
    {
       14, 66, 16, 65, 19,190, 14, 66,  0, 19,  1, 20,254,236,  1, 16,244,239,231,232,214, 73, 63, 65,  2,254,250,250, 14,236, 27, 30,
      250,209,248, 18,  7, 79, 14,209,  1,  1, 64,254,  1,254,253,  2,  0,251,254,179, 34,168,109, 46,  6,  6,  2,  3,  0,  1,  0,  2,
    },
    {
      201, 15,204, 11,207, 14, 72,127, 28,  0, 28,  0, 33,  1, 11,  1,  9,  9, 11, 13, 49, 47,170, 45,  1,  1,255,  0, 17, 14, 44,206,
      201, 13,  6, 19,133, 49, 67, 48,251,253,  5,  0,126,  1,  0,236,253, 80,196,  4, 65,227,119, 77,  3, 10,  4,116,  8,250, 65, 37,
    },
    {
        0,  0,255,  1,254,  3,253,128, 14,211, 23,217,241, 49,  9,164,229,225,225,103,168,217,206, 89,  6,248,249,  6, 82, 32, 39,224,
        0,249,  0,255,  1,129,  0,  2, 53,210, 82,212, 13, 85,238, 42,248,250, 65,  4,245,185, 43, 30,248,252,253, 48, 64,179,251,235,
    },
  },
  {  //  19
    {
       17, 90, 23, 40, 22, 46, 18, 44,125, 10,  3, 10,  3, 12,  2,  7, 54,206, 60, 65,250,251,179,180,237, 25,110,236,  8,  9,252,  2,
       82, 31,243, 78,237,167,239,215,122,  0,  1,244,254,  2,254,254,  2,248, 47, 16, 25, 28,155,164,251, 64,251,100,249,185, 31,  5,
    },
    {
      251,182,250, 54,250,183,250, 58,  0, 20,  0,242,  0,164,  0,243,185, 28,201,219, 72,237, 69,254, 56, 54,119,240, 56,200, 57, 69,
        2, 69,  0, 61,253,254,  1,178,  0,190,  0,191,  0,251,255,255,184, 11, 57,165,  7, 30,129,191,188,191,  5, 60,188, 63,100, 44,
    },
    {
       16,164,142,161,  8, 29,252,148,  0,175,  1, 78,  4, 74,248, 59, 35,164, 89,217,115, 18, 22,103,211,123, 23, 28, 61, 64, 62, 65,
        2, 89,128, 44,  0,208,254,220, 14,178, 17, 65,  2, 57, 14,203,  4,188,252,132,158, 34, 28, 32, 21,139, 30, 78,202, 91,243,  7,
    },
  },
  {  //  20
    {
      214, 10, 87,  7, 11,198, 15,192, 28,  2,229,253,252,231,252,239, 35,234, 17, 82, 49,250, 76,195, 38,  9,  1, 32, 42,  8,  9, 56,
      128,190, 15,228, 25, 43,254,248,253,  1,  6,  4,252,  1,  3,252,234,244,  6,192, 14, 10, 91,253,105,  8,221, 68, 19,246, 20, 19,
    },
    {
       46, 51, 43, 50,166,172, 28, 32, 12, 16, 12, 18, 17, 21, 20, 23, 23, 19,151, 22, 40,178, 38, 49,  5,254,121,  3,210, 38, 51,214,
      211,209,216,171, 95,100,233,145,244,248,241, 13,209, 18, 76,240,199, 13, 78, 51,235, 83, 37, 32,  0,  1,154,117,  0,  9,255,191,
    },
    {
      209, 18,204, 19,209,144, 68, 13, 36,  2, 31,  2, 99,254,107,255, 22, 16, 84, 78,190, 64,208,109,  7,  7,238,231, 25, 17,  6,  3,
       79, 14, 16, 12,211, 74, 78,  5,  2,  2,239,  0,  0,  7,  4,  0,  4,241,196,202, 15, 26,253,  5,199,  8,129,  7,239,250, 39, 53,
    },
  },
  {  //  21
    {
       15, 77, 17,204, 17, 69,  4, 66,  0, 18,  0,239,  1, 15,  1, 14, 24,152,102,223,229,115,131,127,252,251,247, 78, 34, 37,225, 99,
      193,126,196,125, 10, 62,  3,187,  5,  1,255,251,  1,254,  3,  2,  1, 60,253, 28, 45,227,126,107, 21,  7,202,203,133,  3,185,123,
    },
    {
      127,249,253,250,247,251,  1,  0,251,124,  5,  3,  8,  3,  7,  4,233, 95, 27,230,242, 18,245,123,  5,212,248,248, 30,  6, 25,162,
       62,253,255,255,  4,254,252,254, 35,128,253, 64,251, 65,253, 66,218,  0,254,190,229, 65, 35,218, 54,  0, 64,179,196,226, 68,251,
    },
    {
      227,208, 14, 58, 32,184, 39, 50,  2,  3,  0,  8,255,250,  3,  5, 58,196,  6,248,255,238, 32,186,191, 54, 37, 30,  2,  1, 91, 29,
       33,238,227,250,243,235,215,140,  0,248,  0,  3,  0,252,  0,236, 56,236,243,166, 12,208,250,126,  1,207,105, 97,  0, 69,  9,232,
    },
  },
  {  //  22
    {
       23, 41, 29, 22, 23, 23, 30, 20,255, 86,247, 69,  3,198,253, 67,234,235,177,175, 22,231,156,123,249,  5,  8,235,  6,  8,  2,  5,
        4, 34,  4, 38,253, 20,  4,100,  5,198,243,247,  5,254,251,132, 26,254, 74,179,231,158, 36, 28, 10, 62,251, 41,254, 36,  2, 62,
    },
    {
       10, 49,133,252,  5,  4,252,252,191,  8, 62,255, 66,  3,196,  5,163, 63, 79,219, 44, 30,253,238,158, 88,152,128,127,131,255,  6,
       71,255,198,  0, 62,255,208,130, 58, 15, 73,169, 24,252,126,247,250, 20, 17,229,107,201, 81, 29, 55, 27, 50, 42,  0, 13, 83,125,
    },
    {
       10, 62,  7, 58,155, 69, 11, 69,  0,  6,  0,  2,  0,245,  0, 17,102, 43, 85, 85, 81,224,  1,100,195,229,202, 27,  6,113, 25,242,
       63, 73, 75,205,  2, 50, 41,180, 23,  1,  1,  0,  1,254, 10,  0, 23,185,111,230,189, 72,119, 67,203,172,149, 55, 56,239,170,250,
    },
  },
  {  //  23
    {
      254,198,135, 71, 11,199, 32,196,  0,253,  0,  0,255,227,253,241, 49, 53,228,174, 17,254,226, 98,217, 29,249,103,247,244,240,  5,
       66,200, 14,129,253,130, 17,245, 11,255,  0,237,  0,233,255,245, 66,251,245,237, 56,170, 88, 32,  8, 15, 73,193,192, 54,181, 74,
    },
    {
      253,234,252,236,252,242,  2,251,191, 54,192, 55,193, 57,196,186,222,221,227,221, 76,184, 43,190,222,228,231,226, 16,244,  8,245,
      238,255,223,248, 43,194, 85,192,  0,  7,  3, 12, 55,  0, 54,253,  1, 65,246, 63,229,148,170,237,170,128,255, 69,236,  6,  0, 30,
    },
    {
      250,242,249,103,240,227,236,232,  2,  3,253,254,  9,  1,  7,  7,  8,180, 22, 19,  5,184,238,  1, 62, 28, 67,183, 16, 57,214,225,
      235, 31,242,219,253,184,234, 80,251, 12,  1, 18,  0,250,  2, 11, 13, 49, 78, 59,  1, 19, 26,227,200,234,222, 37,247,221,152,226,
    },
  },
  {  //  24
    {
      146, 32, 95,251,221,  3,238,137,  0,252,120,255, 16,  1,242,251,200,198, 21, 17,  0,151,198, 27, 77,229,106,102,252,116, 80,228,
      237,238,200, 14,233,113, 99, 11,217,  0,207,128,184,129,112,255,211, 11,150, 14,190, 87,192,137, 98,137,241, 60,239,114,220, 82,
    },
    {
        1,207,255, 10,  0,  5,253,  0,  5,249,  4,255, 10,  0, 14,  4,253, 65,254,  8,124,187,228, 38,  2, 31,  1,  2, 59, 23, 24, 76,
       68, 54, 97, 10, 87,  6,207,247, 74,252, 53,  1,172,255,  3,255,  0,  0,213,217,253,  7, 55, 13, 43,246, 55, 17, 47,  5, 95,101,
    },
    {
        2, 76, 15, 69,  9, 69,  8,191,255, 22,  1, 19,  0, 16,255,238,203,190,243,251, 20,218,  4,254, 47, 32, 37, 41,237,  7,  6,  7,
       13, 16, 27,206,253,205,251,143,255, 14,  1,  0,  5,  3,  0, 14, 47, 38,  7, 60,245, 22, 69,247,220,223,171,233, 46, 14, 93, 25,
    },
  },
  {  //  25
    {
       11,251,254,219,  4,203,  2,239,  1, 16,  2,  7,  1, 34,  3,  8, 18,209,210,220, 94,232, 41,  6,  2,  7, 33, 50, 90,  7,  2,  2,
      249,244, 60, 33,251,227,  6,225,  0,219,233,254,254,220,  0,  3,  3,251, 54, 21, 80,167, 15,210,255,247,249,245,254,251, 16,234,
    },
    {
        1,  3,248,  0,  3,254,  9,  5,  2,  3,  3,  4,  5,  4,  7,  5,225,219, 42,105, 94,154, 46,106,246,208,249, 27,207,249, 98, 53,
      241,  3, 11,  5,233, 61,  0, 63,  2,  1,  0,  1,  0,233,255,227, 33,  5, 16,233,184, 84,177, 90,201,194,252,206,190,252,123,  0,
    },
    {
       20, 25, 77, 17, 18, 78,200,141, 57,  8,173,255,193, 20, 45,  0,239,253, 59, 66,208,  5, 71, 11, 17, 48, 52,201,246,248,252,  2,
      223,154,193,254, 73,254,194,  0,253,178,251, 65,252,195,184,108, 52,  6,  0,228, 70,212,125,250,246,193,152, 20, 32, 39,255,131,
    },
  },
  {  //  26
    {
      137, 68, 68, 79, 11,186,218, 71, 74, 62,100,192, 65, 49,251, 28,254,214,245,246, 64,209,  3, 39,128, 32, 58, 58,253,118,  0, 16,
      253,252,239,  2,116,251,205,  1,130,243,225, 61, 63,205, 72,191,100,185,103, 76,255, 23,224,197, 52,131,238, 11,253,  9,242, 38,
    },
    {
      200, 61, 99,191, 51,179,175,184, 59,189, 62,182,190, 54,194, 79, 39, 16,108,242,207,181,233, 46,156,228, 56,192,198,148,199,170,
      183,165,240,176,136, 78, 59,220, 91,193,191,254,192,193, 78,195,249,205,251,190, 66,157,193,  5,236,208, 36,228,253,205, 15,192,
    },
    {
      131, 75,139, 72, 86,  8, 32,204,  1,246,  2,236, 34,  1, 52,227, 55, 31, 21,120,247,241,240, 17,197,246,254,  7,255,252,  1,236,
      248, 59, 66,178, 17,  2, 66, 88,236,  1,  9,  1,242,254,100,204, 51, 54,217, 54,253,254, 46,167,237, 10,165,244, 29,235,246,177,
    },
  },
  {  //  27
    {
      212,209,184, 45,140, 35,  7, 64,  9,  5,123,250,253,244,  1, 38,103,195, 47, 79,239,247,187, 69, 24,235,241,253, 85,195,153,155,
      250, 55,205, 32, 34, 53,242, 45,197,254, 64,  4,209,255,216,  0, 28,237,255,144,231,212, 22,211,254,234, 58,231, 70, 55, 95,169,
    },
    {
       81, 20, 78,  6,205, 10,213,  5, 95,133,234,134, 30,136,162,249, 42,165,231,170,246,180,233,210, 35,135,  4,118,127,227,181,232,
      133, 11, 22,125,121,143,  3, 14,211, 70, 66, 72,  6, 73,230,185,188,240,137,118,172,230,148, 79,177,105, 53,178,134,175,184,  0,
    },
    {
      240,190,241,194,241,209,  1, 24,195,213, 62,209, 66, 78, 67, 75,159, 82, 94, 22, 93,156, 79, 72, 55,224,198,189, 36,123, 30,210,
        1, 15,  1,234,247,248, 82,  4, 62,242,193,142,195,192,193,255, 78,206,245,117,250,248,248, 58, 51, 79, 96,224, 64,121,220,128,
    },
  },
  {  //  28
    {
      217,217, 87,252, 88,129,226,124, 76, 63,178,192,209,196,175, 60,249,126,234,255,151, 78,192,229,  3,127,  0,130,138, 30, 34,  5,
       62, 12,147,117,114,115,147,113,210, 61, 86, 65,238,193,156,195, 66, 56,167, 92, 29, 35,128,  8,  5,186, 24,109,200,199, 70, 40,
    },
    {
       28, 66, 12,  2,153,198, 70,253,  5, 32,  1, 58,  4,164,233,192, 13,184, 14,236,187,196,223, 65,  0,247,251,  4, 27,  8, 12,109,
      255,216,254, 93,  1, 58, 12, 41,244,175,  8,  3, 45,250,243, 66,122,245, 83,225,114,215, 16,121,255,246,  6, 22, 48,255,  4,110,
    },
    {
       38, 80,255, 76,222, 25,213,141,  0,211,255,209,  8,195, 12, 62,220, 43,152,254, 37, 37,220, 26,254,  3,  3,242, 26,  6, 57,219,
      183,  1,219,247,203,  7,201,243,128,172,102,187,193, 65, 18,208, 47, 31,194, 74,202,  9, 56, 57, 45,  8,  5,213, 38, 18,235,  0,
    },
  },
  {  //  29
    {
      203, 15,213,134, 24,195,244,133, 16,  0,243,  0,  0,239,239,247, 89, 56, 51,153,185, 34,225, 52,  2, 60,233,200, 36,  1,189,  9,
      254,187,200, 82,252,226,235, 52,255, 48,  0,255,217,255,214,255,200,175,111,171,231,237,188,114,255, 69,255,162,237, 27, 58, 90,
    },
    {
       27, 77,243, 46, 92, 79,138,  6, 62, 70, 64, 63, 15, 38,188,139, 68,178,183, 54,218, 97, 24,138,247,  8, 71, 83,126, 14,254,  1,
      164,205,115, 92, 14,117,247,  7,241,121,  1, 68,135,100, 49,157,186,  0,223, 42,214, 78, 88, 21,194,173,100,207,106,199,153, 25,
    },
    {
      175,160,182,164,193, 74,207,201,  4,  2,  4,  3,252,248,  4, 12, 16, 10,213, 83, 80, 96, 93, 45,  0,255, 54,204, 46, 91,249,233,
       17,103,222,152,253,211,255, 55,  1,  0, 22,  1,219,252, 43,255, 70,251, 75,  0, 21, 78,164,168,121,244, 15,220, 12, 65, 31,  6,
    },
  },
  {  //  30
    {
      105,194,103,198,246,239, 28,175,254,243,255,246,  8,  7,251,101, 29, 25,229, 25,195, 45,191, 12,255,  3,239,247,223, 89,101,  4,
      211, 60,207,242,166,241, 60,164,184,255,  6,  4,152,251,170,  0, 87, 55,  9, 51, 49,249,178,226,180,217,  5, 25,191,189, 60, 27,
    },
    {
       15, 68, 15, 67, 17, 66, 23, 63,  0, 22,  1, 24,  1, 25,  3, 20, 24, 84,211, 76,228, 14, 96, 88,253,218,226, 13, 36,254, 76,178,
       59,176,194, 43,252,139,253,212,  3,  0,  3,  1,  1,243,  3,  0, 57,195,233,113, 46, 16,192, 64, 35,  9, 46,142,255,  0,128,190,
    },
    {
      136,253, 28, 51,202,113, 12, 35,244,138,252,244, 25,  3,121,215,183,199,204,164,243,  6, 90, 18,251,113,189,213,128,129,237, 24,
       39, 56,201,182,149,203, 30, 77,188, 10, 78,118,224,199,242, 57,139, 70,158, 68, 62,137, 97,  7,113,133,195,238, 86, 29,169,241,
    },
  },
  {  //  31
    {
      208, 14, 11,196, 19,214,200, 21,  9,  0,  1,243, 68, 44, 75, 54,237, 16,197,200,235,238, 28,218,  6,  4,237, 23,245,213, 15,247,
       56,  5,192,252, 52,244,198,133, 65, 19,192,254,230,199, 28,233, 88,253,252,250,207, 28, 38,  8,206,  1,  6, 19,229, 28, 43,243,
    },
    {
      244, 53, 11,186, 72, 24,248, 34,127,230,254,205, 10,223,240, 61,119, 42, 27, 44, 53,231,225, 85,255,  3, 79,201, 42,189,  0, 19,
      253,215,211,  4,205,252,243,225,128,203,183,226,112,  2, 11, 33, 65, 51, 90, 22,241,  3,222, 24, 58, 66,161, 38, 73, 59, 81,249,
    },
    {
       24, 87,  9,136, 33, 95,155,210,  2, 27,252,246,  4, 93,251, 62,235, 22,149, 71,  9, 88, 29, 65,  5,152, 40,143, 51, 71, 62,234,
        0, 27, 22,190, 13, 80, 27,220,193,184,248,185,252, 70,240,180, 12, 90,  7, 26,204,229,147,200, 58, 38,  0,243,196,182,142,239,
    },
  },
  {  //  32
    {
      220, 29,221, 32,224, 43, 22,134,212,  3,160,205,222,  8, 73, 37,142,254,184, 94, 40,  3,230, 24,188,196, 37, 54,110,253,  6,134,
      185,226,241,242,136,140,131,247, 76,254,108,243, 19,  2,123,245,115, 45,250, 56,248,  0,216,122,216,193,  1, 28,  4,172,251,118,
    },
    {
      228, 57,220,196, 76, 61,227,241,255,250,  1, 10,  6, 16, 30, 10,  4, 19,221, 27, 42, 33, 78,242,  0,  5, 11,176,246,212, 44,  1,
      221, 39,220, 59,247,208,218,182,247,  1, 31,255, 36,251, 60,255, 12, 11,119,206,196,191,184,194,  1, 58, 35, 88, 14,204,249,188,
    },
    {
       20,211,252, 73, 32, 51,160, 18, 44, 63,252, 61,187,235, 79,140,194, 46,192, 53, 29,  6,  7, 78,241,100, 80,121,  2,254,  0,246,
      191,177,  2, 14, 38,181,137, 48,204,145,212,174,143,157,211,221, 58,199,189, 46, 22,184,221, 50,194,191,133,191, 23,178, 46,232,
    },
  },
  {  //  33
    {
      254,139, 36, 33,229,206,117,116,255,146,253,102,214,225,244,216,  6, 22,158, 19,133,188,173,231,  0,192,157,  1,193,255,155,250,
        1,212,194, 55,192,181,236,  5,  1,210, 70,196,231, 36,230,206,  3,180,214,212, 62,242, 39,232,200,239,234, 31,226,190,248,207,
    },
    {
      198,  8, 36, 72,108,136,182, 72, 76,  0,189,128, 68, 80, 89, 16, 22,118,242,  2, 44,136,199,208,136,130,128,242, 75,252, 20,212,
      129,128,221, 20,147,170,144, 90,190,127,172, 94, 79,108,196, 56, 89,183,143,  0, 25, 18, 27,222,  1, 69,210, 76,144, 65,224,183,
    },
    {
      140, 16, 38, 78, 55, 69,193, 87,231, 61, 72, 67,255,197, 16,135,182,189, 14, 37,  4,245, 77,250,  3, 80,254, 12,  1,  5,226, 30,
      139, 11,  0,222, 27,193, 24, 34,114,120, 63, 65,253, 40, 65,186, 69, 30,228,  2, 33,136, 18, 12, 94,114,238, 19,253, 34,255, 79,
    },
  },
  {  //  34
    {
      156, 20, 27, 25,158, 19,153, 17,200,246,185, 14,200,244,186,244, 23,232, 27,239, 15, 15, 17, 24,189,183,225,224,254,  4,254,  5,
       95,252, 42,  3, 87,252,237,254,203,  6, 66,249,234,  2,210, 99, 50, 25, 27,106,255,  2,253,127, 13,186, 58, 73, 20, 25,166,192,
    },
    {
       20, 41, 23, 43, 24, 49, 21, 42,  0,  8,255, 11,  5, 11,  7, 21, 18, 16, 33, 44,212,208, 77, 26,254,255,255,245,216,221,220,229,
      246, 32,237, 26,241,208,226,232,  1,  4,  3,251,  7, 15,  8,  1, 15,210,180,178, 81,214,222,120,  1, 55, 14,228,191,  0,  9,  4,
    },
    {
      242,238,  0,231,151,151,203, 97, 67,101,222,247, 57, 29,  5,226,140,100, 53,221,  5, 81, 14, 48,239, 40,214, 28,123, 36,  1, 11,
        5,235,232,  3,112, 26,159, 82, 30, 30,251,235,  1,214, 49,226, 41,186,189, 24, 90, 74, 57,  2,130,251,  0,223, 96,218,237, 72,
    },
  },
};

const int g_nsstHyGTPermut4x4[35][3][16] = 
{
  {  //   0
    {  1,12, 0,13, 2, 4, 3, 9, 6, 8,10, 5, 7,11,14,15,},
    {  1, 8, 0, 9, 2, 4,12,10, 3, 5, 6,13,11,14, 7,15,},
    {  4,13, 5,10, 8, 6, 2,14, 9, 1, 0,12,15,11, 7, 3,},
  },
  {  //   1
    {  0, 9, 4, 7, 5,12,11, 1, 8, 2,13, 3, 6,15,14,10,},
    {  1, 0, 4, 3,13, 9, 8, 5,10, 7, 6,11,15, 2,12,14,},
    {  4, 6, 8, 5,10, 1,14, 2,13,12,11, 9, 7, 3, 0,15,},
  },
  {  //   2
    {  0, 5, 9, 4, 1,14, 8,10, 3, 2,13, 7,12, 6,15,11,},
    { 11, 6, 5, 8, 2, 7, 3, 9, 0,13,14, 4,10,15,12, 1,},
    {  1,13, 0, 6, 4, 2,10,14,15, 3,12, 9,11, 5, 8, 7,},
  },
  {  //   3
    { 13, 4, 0, 3, 1,10,11, 8, 9,12,15,14, 2, 7, 5, 6,},
    {  2,11, 4, 9, 3, 8, 6,12, 7, 5,10, 1,14,13, 0,15,},
    { 14,15,13, 2, 8, 7, 6, 4, 1, 5,11,10, 9, 0,12, 3,},
  },
  {  //   4
    {  9, 4, 0, 1,12, 6, 5, 8,15,10,14,11,13, 2, 3, 7,},
    {  0, 1, 8, 5, 7, 4, 9,12,11,14, 2,15, 6,13, 3,10,},
    {  0, 5, 4, 6,13,12, 9,10, 7, 8, 2, 1,11,14,15, 3,},
  },
  {  //   5
    {  1, 0, 9,14, 4, 6, 8,12,15, 7, 5,10, 3,13, 2,11,},
    {  0, 5,12, 9,13, 4,14,15, 3,10,11, 1, 8, 2, 7, 6,},
    {  1, 5,15,13,11,12, 3, 9, 2, 0, 7, 6, 8,10, 4,14,},
  },
  {  //   6
    {  5, 4, 1, 9, 2, 6, 8,10, 7,12,14,13,15, 0, 3,11,},
    {  0, 4,11, 2,12,13, 9, 5, 6, 7, 8, 3, 1,10,14,15,},
    { 13,11, 5, 3,14, 1, 0, 6, 9,15, 8,10, 2, 7, 4,12,},
  },
  {  //   7
    {  8, 5, 3,12, 0, 1,10,11, 9, 4, 2,13, 7,14,15, 6,},
    {  0, 3,12, 9,14, 5,15, 4,10,11, 6,13, 1, 7, 2, 8,},
    { 14,15, 9,12, 6, 1, 0, 4,13, 7, 2, 3,11, 5,10, 8,},
  },
  {  //   8
    {  7,14, 2,15, 9, 1, 6, 0, 5,13,11, 4, 8, 3,12,10,},
    { 12, 0, 2, 4, 5, 6, 9, 8,14, 3,15,10,13,11, 1, 7,},
    {  4, 9, 6, 1, 0,11,12, 8,15, 5, 2,14, 7,13,10, 3,},
  },
  {  //   9
    {  9, 8, 7, 0,15, 5, 3, 2,12, 1,11, 6, 4,13,14,10,},
    {  2, 5, 6, 9, 4, 1, 3, 7,13,12,14,11, 0,10, 8,15,},
    {  0, 4,12,13, 5, 8, 9, 2, 1,11,10, 3,14, 7,15, 6,},
  },
  {  //  10
    {  0, 4, 5, 9,15, 8,11, 1, 7,13,10,12, 6, 2,14, 3,},
    {  2, 8, 6,12, 3,14, 4,13, 0,10, 9, 7, 1,11,15, 5,},
    { 15,11, 0, 3, 4, 6, 1, 7,12, 8,10, 2,14, 5, 9,13,},
  },
  {  //  11
    {  1, 5, 0, 4,13,14,12,11, 9, 8, 6, 7,10, 3,15, 2,},
    {  0, 5,12, 2,11, 8, 4, 3,13,15, 1, 7,10, 9, 6,14,},
    { 13, 1, 9, 8, 5, 6,11, 2,14, 0,12, 3, 7,10, 4,15,},
  },
  {  //  12
    {  2, 1, 6, 3, 9, 5,13, 7,12,10,14,11, 0,15, 4, 8,},
    {  2, 0, 6,12,14, 4,10, 1,13, 8, 5,15, 9,11, 3, 7,},
    {  8, 4, 1,12,10,13, 5, 0, 3, 6, 2, 9,14, 7,11,15,},
  },
  {  //  13
    {  0, 8,11,13, 7,12, 5, 4, 9,15, 1, 2, 3,14, 6,10,},
    {  8,15, 4,11, 0, 7,14,13, 3,12, 6,10, 5, 1, 9, 2,},
    { 12, 9, 0, 5,13,14, 8, 4, 7, 2, 1,15, 6,10,11, 3,},
  },
  {  //  14
    { 13,15, 9,11, 5, 3,12, 4, 7, 8,14, 6, 0,10, 2, 1,},
    {  4, 9, 0,12, 5, 8, 2, 1,13,11,10,14, 3, 7,15, 6,},
    {  7,14, 6, 4,15,10, 3, 9,11, 1, 2,13,12, 0, 8, 5,},
  },
  {  //  15
    {  7,11, 5,15, 3, 9,12, 1, 0,13, 8,10, 4, 6,14, 2,},
    { 12, 0,14, 4, 6, 8,10,11,15, 2, 3, 9, 5, 7, 1,13,},
    {  5, 1, 2, 6,10,14, 4, 9,12, 8,13, 7,15, 0, 3,11,},
  },
  {  //  16
    {  1, 6,14, 0, 9, 5,11,12,13,10,15, 2, 8, 4, 3, 7,},
    { 12, 8, 4, 5, 9, 0, 6, 1, 7,13,14, 2,15,11,10, 3,},
    {  5, 3,11,15,13,14, 1, 9, 6, 7, 8,10, 4, 0,12, 2,},
  },
  {  //  17
    {  1, 5, 2, 6, 9, 3,10,13,14,15,12, 4,11, 0, 7, 8,},
    {  1, 0, 8, 9, 2, 4,12, 7, 6, 5,11,14,13,10, 3,15,},
    {  5, 0, 4,13,11, 9, 8,10, 6,12, 3, 7, 1, 2,14,15,},
  },
  {  //  18
    {  0, 6, 1, 8, 5, 2,12, 9, 3, 7,14,13, 4,10,11,15,},
    {  0, 1, 5, 8, 4,14, 9,13,11, 3, 2,12, 6, 7,15,10,},
    {  6, 4, 1, 2, 0,13, 8,11, 9,12,10, 5, 3,14,15, 7,},
  },
  {  //  19
    {  4, 9,13, 8, 2, 0,12, 7,14, 1, 6, 5,11,10, 3,15,},
    {  4, 0, 8, 1,12, 5,14,15, 6, 3,11, 7, 9, 2,13,10,},
    {  8, 4, 2, 0,14,13,12,10, 9, 6,11, 5,15, 1, 3, 7,},
  },
  {  //  20
    {  1, 5,12, 4, 8, 7, 9, 0, 6,13,15,11,10,14, 2, 3,},
    {  4, 0, 1, 9, 6, 5,13,15, 2,10,14, 3,11, 7,12, 8,},
    {  8, 5, 1,12, 4, 2, 0,13,14, 9, 3,11, 6,15,10, 7,},
  },
  {  //  21
    {  5, 8, 0, 1,13,12, 2, 4,11,10, 9, 6, 3, 7,14,15,},
    {  1, 4,11,12, 8,10, 0, 7, 2,15, 3,13,14, 6, 9, 5,},
    {  0, 8, 4,13, 5,12, 3, 9,15, 1, 7, 2,14,11,10, 6,},
  },
  {  //  22
    {  0, 4, 2, 8, 6,11, 3,12, 9, 7,10,15,14,13, 5, 1,},
    {  4, 8,14, 0, 6, 2, 7,10,12, 9, 5, 1, 3,13,11,15,},
    {  4, 1,12, 9, 0, 5,10, 8, 2,13,15, 6,11,14, 3, 7,},
  },
  {  //  23
    {  0, 4,13, 1,12, 9, 5,15, 8, 3,11,10, 7, 2,14, 6,},
    {  6,11, 3, 4,12, 2, 1,13, 7,15, 0,10, 8, 9, 5,14,},
    {  0, 4, 8, 5, 1,12, 3,13, 7, 9, 2,11, 6,15,10,14,},
  },
  {  //  24
    {  5, 1,14,13, 2, 4, 9, 8, 6, 0,15,12,10, 3,11, 7,},
    {  1,13, 9, 5, 4, 6,14,10, 8,12, 0, 2,15, 3, 7,11,},
    {  0, 1, 8, 4, 6, 9, 5,12, 7,13,11,15,14,10, 2, 3,},
  },
  {  //  25
    {  0, 4, 8, 2, 1, 5, 9, 3,12, 7,13, 6,15,11,10,14,},
    {  0, 5, 8, 2, 9,12,15,10,13, 4, 1, 7, 3, 6,11,14,},
    {  6, 4, 0,10, 3, 8,12,14, 7, 2, 9,11,13, 5, 1,15,},
  },
  {  //  26
    {  0, 1, 2, 9, 6,13,10, 3, 7,14,15, 5, 8, 4,11,12,},
    { 12, 5, 1, 4,13, 8, 9, 3,10, 0,11, 7, 2, 6,15,14,},
    {  1, 8, 0, 9,14, 5, 6,12,13, 2, 4,10, 3, 7,15,11,},
  },
  {  //  27
    {  9,10, 5, 1, 0, 2, 4, 3, 6,13, 8,15,14,12, 7,11,},
    { 11,14, 7, 4,10, 3,12, 2, 0,15, 6, 1, 8, 9, 5,13,},
    {  0, 1, 4, 7, 9, 8,13,14, 3, 5, 2,12,11,15,10, 6,},
  },
  {  //  28
    { 13, 9,12, 5, 0,14, 4, 1, 6, 2,15,10, 7, 8, 3,11,},
    {  0,12, 5, 2, 4, 6,15, 7, 8,14,10, 1, 3, 9,11,13,},
    {  7, 5, 3, 9, 8,11, 1, 2,10, 4,13,15,12, 0, 6,14,},
  },
  {  //  29
    {  0, 1, 4, 9,13,10, 7,12, 5, 6,15, 8,11, 3, 2,14,},
    { 13,14, 5,10, 6, 3, 4, 2, 0,11, 7, 1,12,15, 9, 8,},
    {  1, 4, 9, 2,13, 0,12,10,14, 7, 5, 3, 6, 8,11,15,},
  },
  {  //  30
    { 13, 2, 5,12, 8, 9,14, 0,10,11,15, 1, 3, 4, 6, 7,},
    {  5, 8, 0, 4, 2, 1, 6, 9,11,12,10,13,15, 7, 3,14,},
    {  1,14,10,13,11, 6, 3, 5, 4, 8,12,15, 2, 9, 0, 7,},
  },
  {  //  31
    {  2, 5, 0, 1,13,10,14, 3, 9, 8,15, 6, 4,11,12, 7,},
    { 12, 3, 4, 7, 5, 9, 0,15, 6,10,11, 1,13, 2,14, 8,},
    {  2, 0,11, 4, 5, 7,12, 6, 3, 8,14, 1,15,10, 9,13,},
  },
  {  //  32
    {  1,11, 6, 4, 2, 0, 8,10,12,14, 7,13, 3, 9, 5,15,},
    {  9,13, 0, 1, 6, 2, 4, 7,10,15,12, 5,14, 3, 8,11,},
    {  8,13, 9, 4, 0, 7,11,14,12, 1, 2,15, 5, 3, 6,10,},
  },
  {  //  33
    {  8,12, 0, 2,14, 7, 1, 3,10, 4,13,15,11, 6, 5, 9,},
    { 13, 0, 9,15, 4, 6, 5, 7,12,14,10, 3, 1,11, 8, 2,},
    {  6, 8, 3,13,10,12, 2,11, 1, 7, 0,14, 4,15, 9, 5,},
  },
  {  //  34
    {  0, 8, 2, 1,10, 4, 9,14,15, 5,12, 3, 6,11,13, 7,},
    {  0, 5,12, 7, 9, 8, 6,10, 1,14, 3, 4,13,15,11, 2,},
    {  2, 3, 1, 6, 8, 0, 7, 4, 5,10,14,11, 9,12,13,15,},
  },
};


const int g_nsstHyGTPar8x8[35][3][768] = 
{
  {  //   0
    {
      228,231,217,226,226,232,218,225,230,232,218,227,227,228,213,223,219,224,215,224,230,232,212,224,226,230,217,220,215,225,219,207,215,205,  5,255,213,206,  5,255,210,207,  7,  2,212,210,  5,  1,209,206, 10,  3,205,208,  9,  0,214,211,  9,255,209,210,  9,254,254,254,255,255,  0,  3,255,254,  0,  1,  0,255,255,  2,  0,255,  3,  5,255,255,  2, 10,253,  0,  3,  4,253,  1,  5,  3,251,  0,225,227,225,241,220,219,225,237,231,232,218,222,230,230,233,229,241,250,250,252,245,249,243,249,212,213,215,214,216,209,218,227, 59, 68, 53, 60, 44, 48, 55, 46,237,244,250,240,227,230,241,  0,255,  0,  2,  0,  1,254,  1,255,248,245,248,252,247,246,243, 37,245,236,246,235,248,245,247,247,250,254,254,255,251,250,252,253,  1,  1,  1,  1,  2,  0,  2,  3,255,  1,  0,  1,  0,255,254,  7,
        2, 30, 21, 16,252,  0, 24,255,  5, 25, 22, 10,254,252, 22,  0,  3, 15, 25, 16,243,245, 25,255,  5, 12, 17,251,255,  1, 20,  2,  1,255,251,252,  0,  2,250,251,250,250,247,251,  1,  3,253,252,255,254,250,253,  0,  1,251,  1,249,247,  3,252,241,246,245,  5,214,204,192,192,210,196,193,196,203,197,193,194,217,206,192,194,206,199,193,193,208,198,195,192,203,199,192,193,203,198,198,192,  3,  4,  6,  0,  4,  3,252,  2,  2,  3,  2, 10,  6,  8,  9,  6,251,250,250,252,254,249,249,248, 17, 22,233,248, 15, 15, 15, 19,186,187,175,176,184,180,198,200,246,246,240,245,242,228,241,244,  0,  0,  1,253,  1,  0,253,254,255,253,246,198,  1,  1,252,249,208,207,203,200,205,203,205,207,208,209,204,199,211,215,204,209,191,191,194,191,191,191,193,192,192,192,192,189,192,192,192,192,
      242,225,242,230,236,251,255,248,240,  8,248,247,239,255,254,254,245,238,251,242,240,254,  8,253,244,241,253,218,246,  0,  6,255,190,196,214,220,197,198,218,213,179,201,214,218,205,203,222,217,186,190,213,213,195,195,218,208,201,203,213,211,196,197,216,205,255,255,252,255,255,255,253,253,  4,253,254,  2,252,254,  0,  2,255,254,253,253,255,253,253,250,  0,  1,251,254,  0,254,251,249,251,240,251,247,249,239,250,239,  7,242,231,237,236,239,238,239,252,242,249,254,250,249,251,251,  6,  8, 11, 17,  9, 11,  1,  7,250,  1,252,  0,252,254,253,255, 17, 13,  5,  3,  1,255,255,255,252,251,  0,253,251,241,239,234,253,226,  7,255,255,251,250,249,252,255,247,248,249,248,248,247,  1,254,252,249,249,240,247,231,255,  2,  1,  0,  0,  1,  1,  1,  2,  5,255,  0,  0,255,  0,255,
      234,220,228,232,229,231,238,236,224,225,221,225,229,229,223,225,224,223,222,223,223,227,229,220,226,228,233,225,227,227,234,221,  1,  6, 41, 40,252,  3, 40, 44,  7,  8, 44, 43,251,  4, 45, 46,  6, 12, 47, 48,  0,  5, 40, 38,249,255, 37, 38,254,  3, 33, 34,  3,  2,  3,  2,  2,  2,  2,  2,  3,  2,  1,  1,  4,  3,  0,  0,  2,  1,  3,  3,  1,  1,  2,  1,  0,255,  3,  1,  0,  0,  2,  1,212,223,218,218,221,228,220,226,231,244,244,244,243,244,240,242,221,225,231,231,228,230,228,226,217,217,218,218,227,230,236,238,202,200,202,200,199,198,200,199,201,197,199,199,199,198,198,197,237,237,231,229,233,233,215,216,244,247,241,239,245,245,227,225,254,255,  0,  1,  1,  1,  0,255,254,254,  0,  0,  1,  1,  0,  0,  0,  0,  1,  2,255,255,255,255,  3,  2,  3,  3,  1,  1,255,  1,
    },
    {
      236,237,236,240,226,232,232,233,204,218,211,211,222,229,220,221,229,226,217,217,227,226,221,227,226,224,222,217,222,218,214,206,  6,185, 13, 18, 15,193,141,250,  8,  6,187,189,  7,  2,186,190,  9, 17,189,193,  9, 17,190,165, 13, 20,195, 64,  9, 14,196,151,  6,252,  0,254,251,254,  0,  0,  0,  0,254,253,  0,  1,252,250,  0,  1,250,238,  2,  0,239,252,  1,255,218,206,  2,254,191, 41, 42,  9,  4, 23,249,  4,245,  2,242,245,238,243,232,239,236,241,179,180,186,195,185,196,189,195,249,235,195,196,122,190,122,193,252,  3,  0,  0,  0,  0,  0,  0,253,  3,  2,  0,  1,  0,  4,254,185,184,193,191,200,193,200,191,194,201,190,191,193,191,192,191,  0,254,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  1,254,255,  0,255,  0,255,  1,255,  0,243,239,244,247,238,240,240,242,
      252,  8,  1,254,  0,  2,255,  3, 29, 17, 18, 19,255,248,  3,  3,251,254,  1,  0,  4,  1,254,  3,252,254,254,  1,250,  4, 18, 14,226, 41,195,171,217, 29, 62,198,232,219, 19,  6,236,  8, 23, 12,191,194,  2, 37,197,190,  1,125,246,255, 73, 30,214,254, 58,  5,199,190,242,244,199,194,242,247,199,195,244,249,198,193,244,242,194,228,255,252,192,205,253,  2,196,193,  1,235,208,192,254,250,195,205,214,208,204,197,200,208,240,249,247,254,241,249,247,251, 52,  5, 66, 39, 60,246,127,  6, 52, 53, 61, 56, 56,  8, 66, 58,249,228,248,234,247,228,246,230,215,190,212,193,218,191,214,191, 64, 64, 69, 67, 68, 64, 62, 63, 64, 74, 68, 69, 64, 64, 62, 60,254,  0,  0,255,  0,  0,  0,254,255,254,  1,255,  1,244,  0,255,  1,253,  0,  2,  0,  0,  0,  0,255,255,  0,254,  0,  0,255,255,
        0,245,252,246,  0,  1,  0,254,254,255,254,255,255,255,  1,255,  0,  1,  0,255,  0,  2,254,254,  2,252,  1,255, 17, 11,  1,254,252,  0,207,214,250,252,202,195,251,255,207,207,250,  1,203,200,252,254,195,193,252,  1,192,194,  1,255,176,250,255,255,167, 13, 55, 66,215,205, 50, 61,231,254, 54, 64,224,216, 49, 65,225,231, 13, 64,207,248,254, 70,242,  2, 64, 72,203,193, 62, 61,190,192,255,252,253,252,253,251,241,242,  0,  2,  6, 10,  4,  5,235,247,212,213,214,219,198,206,237, 40,251, 66,252, 54, 11,248,227, 13,246,  3,244,  4,239,  4,237,  0, 41, 57, 35, 50, 40, 58, 11, 18,249,247,240,245,  3,  4,247,244,  2,  3,  1,  2,  8,  6,253,  2,204,206,206,207,204,206,202,203,205,205,204,204,206,207,199,199,194,194,194,194,191,194,192,191,193,192,192,192,193,192,191,193,
      225,225,218,222,211,217,207,220,223,223,237,221,227,231,209,210,224,225,216,224,212,219,206,223,219,223,216,206,230,232,225,241,202,200,  0,  0,206,207,  0,  0,205,203,  0,  1,217,217,255,255,204,202,  0,  0,206,207,  0,  0,203,201,  0,  0,203,199,  0,  0,190,191,  2,  3,191,191,255,255,192,192,  1,  0,192,192,  1,  1,190,191,  1,  2,191,191,  0,255,190,190,  0,  0,190,191,  1,  0,232,234,229,232,229,231,218,216,230,239,183,183,186,186,149,152,231,232,228,231,230,232,223,222,196,194,196,195,197,199,214,224,195,195,197,195,195,194,192,191,191,193,203,201,203,202,197,196,223,221,220,218,219,218,246,244,221,223,224,220,227,225,226,223,  1,  2,  1,  1,  1,  1,  2,  1,  2,  3,255,  0,254,254,  0,  0,247,247,247,247,246,246,255,254,247,248,249,248,250,249,251,250,
    },
    {
       13, 26, 34, 23, 20, 30, 28, 21, 29, 27, 28, 28, 30, 25, 27, 28, 31, 21, 27, 30, 32, 19, 26, 28, 33, 26, 22, 24, 38, 21, 20, 24,  5,  2,239,249,  5,253,231,250,  5,  4,237,244,  5,  3,238,245,238,  3,239,245,  7,  4,239,247,  8,  3,239,240,  9,  2,239,239,  1,  1,  0,255,  1,  1,255,255,  1,  1,  0,255,  2,  1,  0,255,  2,  1,  0,  0,  2,  1,  0,  0,  2,  1,253,252,  2,  0,253,251,145,  7,  9,  4,  1,  4,  1, 27,  9, 47, 47, 40, 41, 39, 41, 39,255, 60,  0, 55, 25, 11, 26,250, 18,252, 36, 33, 40, 40, 40, 42,  0,  4,  0,  0,  0,  0,  0,  1,253, 36,  3,195,  1,  5,  1,255,  1,240,  0,254,  1,130,  1,254,253,255,253,  3,253,254,249,  1,  0,  0,  1,  0,  0,255,  0,  2,  0,253,  0,255,254,255,255,254,255,  0,252,  4,186,  6,186,255,254,  0,254,  0,251,253,252,252,
       51, 66,205,201,  8,255,254,249,252,  5,  0,  1,  8,  1,248,250,  1,  0,  0,254, 17,  1,253,255, 54,245,  0,255,211,  0, 16, 12,  2,  4, 66,  3,  1,  0, 83, 63,  0,  0,204, 53,  2,251,209,126, 25,254, 73,193,  1,  0, 72, 60,  0,  0, 74, 58,  0,255,205, 72,  0,255,245,  1,  0,  0,250,  0,  0,  0,211,  3,  0,  0,245,  0,  0,  0,238,249,  0,  0,  1,250,  0,  0,229,252,  0,  0,  2,254, 99,253,  5,  0,  1,255,  3,  0,  5,204,  0,254,  0,255,255,  0,249,197,  6,240,244,238,251,249,  1, 15,  1,  7,255,253,  0,  1,211, 43,  0,  0,  0,  0,  0,  0,  5,  8,246, 62,252,  1,248,  1, 35,226,  4,  3,204,124,  3,  1,  1,253,130,254,224, 26,253,  2,255,  0,  0,  0,  0,  0,  0,255,255,  3,  3,  0,255,  1,255,  0,  1,255,  7, 68,127, 67,  9,193,194, 11,  0,248,255,  0,  0,  1,
       42,190, 63, 65, 19,252,  9,252,247,  2, 33,  0,233,  0,254,  0,225,  4,255,254,251,255,237,248, 32,  1,  1,  1, 49,  0, 20, 16,  1,  2, 70,254,  0,  5,254,253,  0,  0,  0,  0,255,  1,125, 70,  0,  1,  1,  0,  1,  0,  1,  1,  0,  1,  0,244,  0,  0,129,255,  0,  0,  1,249,  0,  0,  8,  2,  0,  0, 46,251,  0,  0,249,254,  0,  0, 22,  6,  0,  0,252,  7,  0,  0, 26,  0,  0,  0,249,  1,  7,241,254,255,  6,238,  3,  1,  1,255,250,  0,  0,  1,255,  0, 26,244,255,224,  1,  0,253, 10,  0,  2,252,255,  3,  0,255,  0,244,  0,  0,  0,  0,  0,  0,  0, 10,247,  5,253,246,  5,243, 10,175, 51,253,254,197,249,  3,  7,  0,255,128,254,228,225,253,252,  5,254, 22, 23,  0,254,  0,255,254,251, 26, 23,255,254,  0,254,  3,  0,252,184,190,192, 55, 61,  8,213,  1,  7,255,  2,  0,  2,
        0,231,225,212,235,232,231,234,232,224,225,227, 12,228, 25, 20,254,231,223,226,224,232,244,239,133,239,193,232,209,228,235,231,  0,  0,204,196,  0,  0,201,204,  0,  0,203,202,  1,255,199, 72,  0,  0,203,203,  1,  1,203,204,  0,  0, 73, 73,  0,  1,201,199,  0,  1,  1,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  4,  0,  0,  4,  4,  0,244,251,246,247,249,247,243,255,253,215,215,212,212,211,210,  5,  3,251,246,239,236,238,237,  3,254,215,213,215,213,211,209, 54,207,  6,  7,  6,  6,  6,  5,245,247,  4,  5,  7,  7,  7,  7,205, 19,  1,  1,  3,  4,  2,  2,  2,  4,  1,  2,  2,  2,  2,  0,255,  1,234,234,  2,  1,  2,  1,  4,  3,234,234,  4,  4,  4,  4,  6,255,  2,  2,  5,  5,  6,  6, 43, 10,  1,  1,  5,  5,  6,  4,
    },
  },
  {  //   1
    {
      219,211,  1,205,222,207,255,208,217,208,255,205,219,208,  2,211,213,210,252,218,213,212,254,234,210,211,252,213, 60,215,  1,221, 12,  6, 63, 51,  6,  6, 64, 54,  7, 11, 63, 55,  5, 12, 64, 57,  7,  9, 64, 60,  9,  6, 67, 62, 10,  6, 64, 55,249,250, 68, 57,  0,  0, 20,254,  0,  0, 15,255,  0,  0, 16,254,  0,  0, 14,254,255,  0, 19,  0,  1,  1, 21,  1,255,255, 25,  1,  1,  2, 22,255, 38, 27, 25, 33, 19, 18, 18, 19, 24, 17, 26, 25, 29, 33, 25, 24, 20, 17, 24, 22, 29, 25, 21, 21, 72, 61, 51, 57, 28, 37, 44, 37, 49, 50, 47, 52, 33, 36, 32, 36, 56, 56, 51, 60, 43, 40, 33, 38,255,  0,  2,  2,  7,  6,  4,  8, 10, 12, 13, 13,  6, 16, 12, 15, 30, 27, 24, 22, 30, 22, 21, 23, 25, 27, 19, 17, 21, 17, 20, 17,  0,  0,  0,  0, 10,  5,  3,  4,  1,  0,  2,  1,255,  3,  1,  1,
        8, 34, 23,243,  4, 32, 13,240,  6, 44, 22,245,  0, 49,  8,242,  4, 32,  1,241,  3, 31,245,237,162, 20, 13,232,  8, 36, 14,233, 23,  7,  3,  1, 21,  7,  1,  1, 24,  6,  4,  7, 19,  5,  2,  6, 18,  6,252,252, 18,  6,249,249, 19,  8,  4, 12, 14,  5,  0,  7,  8,  9, 20, 28,  8,  9, 22, 30,  9,  9,  6, 12,  9,  9,  6, 10,  7,  6, 19, 33,  6,  7, 21, 33, 11, 14, 20, 29, 11, 12, 12, 23, 20, 19, 21, 24, 20, 18, 21, 25, 25, 27,  8,  9, 25, 27, 12, 13, 31, 30, 30, 34, 30, 27, 26, 32,248,245, 64, 71,247,242, 67, 70,246,242,252,  1,246,242,255,  8,254,255,  0,  4,254,254,  7, 18,  0,  0,254,254,253,  2,255,255,  0,  1,  1,  0,  1,  0,  1,  1, 21, 22, 44, 42, 22, 25, 39, 30, 12,  9, 33, 31, 14, 16, 33, 21, 62, 62, 64, 64, 61, 62, 64, 64, 62, 62, 64, 64, 62, 61, 63, 63,
      241,213,241,219,237,216,237,222,231,218,232,225,236,229,237,226,243,215,244,217,240,218,241,219,230,209,230,216,233,205,233,215, 67, 73, 62, 64, 66, 82, 56, 68,  1,252,254,254,  2,  2,  7,  3, 21, 20,  7, 10, 21, 25,  6,  9,  5,253,  5,252,254,252,252,251, 40, 45, 55, 54, 39, 44, 55, 54, 53, 52, 46, 51, 53, 52, 53, 48, 55, 55, 49, 55, 55, 56, 46, 52, 55, 55, 55, 62, 55, 55, 49, 55,250,  1,250,254,254,  3,249,  0,166,181,200,206,240,240,210,207,244,  3,254,250,254,  7,  0,  4,248,239,  0,250,235,235,254,251,  2,  1,255,255, 10, 12,255,  0,  1,  1,253,252,  3,  7,254,254,  2,247,  3, 14,  1,  7, 28, 35,253,251, 11, 26,  2,  9, 38, 48,245,248,254,  1,  1,254, 20, 17,246,245,  0,  3,  2,  2, 12,  5,  1,  0,254,255,  4,  4,  0,  0,  1,  2,255,255,  3,  4,255,255,
       62, 43, 60, 64, 64, 57, 61, 66, 74, 64, 67,117,122, 63, 66, 49, 61, 57, 63, 59, 73, 60, 64, 59, 57, 67, 70, 58, 67, 70, 70, 58,196,200,246,250,190,199,233,250,252,246, 64, 65,  2,  2, 61, 57,255,  1, 34, 40,  4,  4, 37, 43,  4,249, 48, 47,  3,  2, 46, 50,  2,  1,  4,  3,  2,  4,  7,  4,252,255,  1,254,255,  5,  0,  1,  1,  1,254,  6,  1,  1,  2,  8,255,  0,254,  1,  0,  0,  0,  2, 35, 35, 26, 32, 22, 22, 35, 34, 49, 62, 11, 24,  1,  1, 39, 45, 32, 42, 43, 36, 26, 29, 35, 34, 31, 23, 35, 31, 35, 34, 35, 28, 63, 65, 62, 61, 65, 62, 72, 67, 60, 59, 64, 62, 62, 64, 73, 70, 51, 42, 18, 26,246,254, 10, 21, 51, 43, 22, 28,248,  0, 10, 18,  1,251,  0,  1,  0,  0,  0,  0,  2,  5,  0,  1,  0,  1,  0,  0, 26, 25,  7,  7,  2,  2, 14, 14, 23, 83,  7,  5,  1,  2, 12, 10,
    },
    {
       13, 10,168,221,249,253,  6,  1,  4,  0,255,255,254,  0,255,  0,  4,  5,  1,  0,  2,  1,255,254,  4,  1,  0,  1,253,252,251,246, 12,239, 29,226, 11, 14,234,250,  7,216, 98,  3, 16,  7,198,253,  9,  9,192,246, 13, 14,186,  1,250, 19,197, 50, 11, 17, 69,191,  0,253,  1,254,255,  1,255,  0,254,  2,  1,  0,  1,  2,250,  0,  1,  2, 33,  1,  1,  2,248,  0,  1,255,251,  4,255,  0, 12,248,  8,186, 12,194,250,249,  8,253,  9,191, 10,193,254,192,255,194,  0,192,  0,192,251,190,  1,194,  5,187,254,185,  4,200,  1,192, 32,217,237,252,  6,  1,249,  2,211,186,249, 27,  2,  0,173, 13,251,196,198,197,212,189,194,242,249,  3,255, 77, 51,  0, 71,  2, 23,  0,  2,  1,  1,255,  0,255,  6, 15,  3,  0,  0,  0,255,  1,254,  5,  0,  8,  1, 28,255,  3,  0,255,  0,  5,248,  1,247,  0,
       63, 49, 93, 35,199, 50,248,  1,253, 60,251,255, 27, 23,  3,  3, 17,251, 12,  8,254,212,  0,  2,205, 29,  7,  0,179,236,  3,247,239,  7, 34,166, 46,  7, 24,197,255,  4,220,192,214,244,128,194,  1,  0, 46, 13,255,  3,121,194,  0,250,128,196,  8,252,208,144,  3,  0,  1,255,253,  0, 12,255,  1,  0, 13,  0,255,  0,255,  1,  0,  0, 10,255,  1,  2,253,  4,255,  0, 43,241,253,  0,  1,  4,  4,246,251,236,249,194,  3,184,  5,  2,255,246,255,  4,248,255,252,  3,  7,255,253,  2,255,255,  0,253,255,  4,  0,  2,255,254,173, 34, 17, 27,245, 35, 12, 10, 32, 65,192,211,  4,255,222,227, 56,190,  3, 60,118,220,  2,  1,242,193,  3,129, 63,192, 79,188,252,  3,  4,  0,254,  3,  1,  0,255,  0,  1,  8,  4,254,  1,  4, 26,  4,232,  0,  3,245,245,208,  1,252,255,  2,  0,  4,  1,254,
      221,199,  2, 14,243,252,  2,  4, 50,201,  0,249, 48,208,  1,  0,196,247, 15,252,191,142,252,143,202,175,250,  0,  6,232,253,  2,234,246, 29,  2,219,  1,239,254,  8,  2,  9,  3,239,252,  7,  2,247,  9, 27, 47,  2,191,193,253,  1,  0, 58,255,242,252, 30, 15,255,  0,  6,252,  4,  1,252,253,  1,  0,  2,254,  0,  0,252,250,  3,  1,244,254,  1,  0, 14,  0,  2,  0,251,252,  2,254,251,248,244,251,255,245, 10,  5,255,235, 76,255,248,251,253,248,  4,  1, 16,  0, 11,255,  3,124,  2,249,251,  2,253,253,  0,  1,255,  2,  5,208,242,230, 16,193,245, 12, 20,189, 69,230,246,212,208,237,197, 68, 70,209,193,187,244,252,172,223,205,255, 53,248,247,246,  0,217,  0,251,  1, 13,  1,  0,254,248,  1,251,  1,251,253,  0,254,241,  4,242,247,241,  0, 15,  0,255,247,233, 15,233,  1,216,
        4, 11,  0,250,  4,  0,253,252,252, 10,255,  4,  2, 27,248,249,182,177,  0,244,252, 69,254, 10,249,232,254,251,  5, 10,255,  0,  7,199,226,  0, 58,208, 40,250,243,200,191,  5,  0,197,189,249,  1,189,255, 28, 71,188, 73,254, 10,192, 51,  2,255,207,254,252,  2,192,251,193,250,193,255,194,  0,193,  2,192,  0,194,  2,192,  1,192,  2,192,  1,192,  2,192,255,193,254,191,  1,194,  2,195,252,244,  0, 20,245,  9,246, 13,184,248,253,249,238,253,250,  0,  7, 86,254, 11,254, 10,  5,  4,238,250,  5,  2,  3,  0,255,  3, 31,230, 18,186, 21,213,255,214, 14,240,  0,245, 20,251,  1, 17, 10,  0,254,  7,254,  8, 23,253,  3,254,  8,255,  8,  0, 86, 35,242,194,255,198,  2,189,  1,198,  0,190,  2,194,254,191,  0,194,  1,206,252,192,254,197,252,198,  1,213, 13,233,  1,210,  1,199,
    },
    {
        1,  7, 29,248,255,255,254,  6,248,  9, 15, 16,254,  3, 14, 14,253,  1,247,  3, 30, 42, 13, 13,  0, 20,  4,253, 30, 37, 13, 15,  3, 37,239,  7,255,  4,245,  9,  6,249,235,247,  5,203,251,244, 11, 41,253, 62, 13,  5, 57,251, 17,  0,242,176, 20,  8, 59,  0,255,246,251,255,  2,255,  7,248,  3,254,250,255, 71,  1,  3,  0, 27,  0,226,  6,  0,  0,  5,  0, 12,  0,241,255,  0,  0,  5,  0,244,255,  6,  0,  9,255,244,254,  0,  0,255,253,  0,  0,  0,116,  0,  0,  0,  0,  0,  1,  0,  1,  0,254,  0,252,  0,  0,255,  1,  2,245,252,255,253,253,251,  2,  0,255,254,  1,  0,  1,252,  0,  1,  2,  1, 10,  1,250,  1,254,254,254,252,253,  6,218, 17, 60,  0,254,  0,  0,  1,  0,  0,255,  0,  0,249,  1,  0,  2,  0,  0,  0,  2,  0,  1,255,  0,  0,  0,  1,  0,  2,  0,  0,  0,  0,255,
      255,252,247,252,  0,  8,254,250,  9,250,246,  1,254,  1,  0,  4,  1,  1,  0,  2,222,245,  0,  0,252,254,  0,253,221,244,255,  0,  2,210,  8,250,  0,  1,255,  3,247,245, 19,  1,  2, 22,  0,250,249,199,  9,194,249,254,207,  3,242,252, 19, 80,245,250,197,251,  2,193,  5,  2,255,  0,250,  6,253,  4,  5,255,185,255,252,  0,229,252, 31,  0,  0,  0,  0,  0,243,  1, 17,  0,255,  1,  0,  0, 11,  0,  0,251,253,  0, 24, 12,249,250,  1, 14,239,  9, 18,157,255,  1,  0,255,254, 11,  5,  2,250,  1,254,  0,  4, 10,  3,198,209,  3,245,254,  3,  0,  4,  2,  1,  6,255,248,253,  4,  9,254, 59,  0,255,213, 30,230,171,208, 58,  0,  0,  0, 70,126,181,133,  2,  0,  2,  0,255,  0,  0,  0,  4,  0,  8,  0,  0,  2,  0,  0,239, 11,  0,  0,  3,252,251,254, 60,  0,253,  0,  1,  0,  1,  0,
        0,250,  2, 11,  0,253,  5,  0,  0,  0, 18,245,  0,255,243,245,  7,  5,  4,  3,  0,218,249,238,254,  2,250,252,  0,226,242,  1,249,  7,  7,  0,  1,251,  1,244,  0,255,  1,  4,251, 12,254,  5,251,245,  2,253,  2,  2,255,254,  0,255,  6,  3,  0,229,  6,  0,254, 68,  0,255,255,  5,250,  2,  0,255,  0,  0,  1,  1,  0,  0,  2,255,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,255,  1,254,  6,252,  0,250,247,255,  2,  0,240, 17,243,239,235,254,254,247,  1, 59,218,253,189,254,  0,250,  0,193,255,251,209, 52,254, 32,229, 54, 55, 26,228,250,251,  5, 17,  9,  0,255,209,  3,251,  6, 43, 68, 70, 64,254,252,254,  6,  5,234,187,  3, 44,  1,255,253,  3,  0,  1,  6,  0,253,  1,  1,254,  0,  0,  0,  4,  8,244,  0,255,  0,  0,252,  0,195,  0,255,  0,  4,253,255,  1,
        1,  2,233,  3,  0,253,  0,  1,253,252,237,253,255,251,255,  0,249,249,250,  0,  0,  4,254,  1,  2,231, 10,  0,  0,  4,252,241,255,  6,  4,246,  2,  2, 12,  0,  8, 26,  1,  3,254, 22,  8,  3,  3, 24,255,255,254,248,  0,  0,  0,  6,255,  0,  0, 19,255,254,  0,  0,  0,  0,  0,252,  4,  0,  2,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,254,  0,  0,  0,  0,  0,  2,254,250,252,252,251,246,  0,  9,  3,255,  5,254,255,252,252, 11,  2,  9,255,190,244,255, 68,  3,  2,  8,  2, 60,207,252,251,252, 13,241, 25,206,207,233, 31,  1,  1,  0,249,253,255,255, 52,194,  2,249,  1,203,214, 20, 10,200,  4,253,255,237,234, 66, 68,240,  2,  1,  0,  1,  1,252,  0,250,255,253,  1,  1,  1,  1,  0,  2,255,  0,  0,255,  2,  3,  0,  4,  0,  2,  0,  0,  1,  0,  0,
    },
  },
  {  //   2
    {
       18, 33, 21, 21,  1, 36, 15,  5,224,211,  5,229, 69, 57, 40, 14,200,  5,215,  8,216,232, 42, 57, 85, 36,246,184,108,  6, 86, 30,243,216,249,  1,215,213,253,250,214, 19,254,235,219,201,250,253,250,196, 16,  1, 10,194, 59,253, 26,197, 63,195,139,170, 67,253,253,  8,255,  0,  1,  3,255,  0, 30,  2,255,  2, 28,  8,255,  0,234,234,  5,255,255,244, 30,255,  4,247,246,  0,  0,238,254,252, 64,250,240,236,249,251,  5,125,240, 12,192,  5,200,162,242,180, 17,248,  3,180,  2,240,239,  1, 22,201,215,191,175,237,217,216,250, 57, 64,249,195,  4,251, 32,  0,  2,245,255,210, 12,240,247,  5,  3,254,  5,254,253,  0, 43,253,255, 27, 18, 18, 54,244,249,  0,254,231,255,  1, 64, 68,247,248,255,  1,254,255,  1,224,  4,255,  0,255,  0,  0,  2,  3,233,  0,  0,  1,254,  0,  5,251,  2,
       56, 15, 24,234,  2,213,207,229,  2,  1,  3,153, 29,245, 18, 52, 39,169,200,255,198,244,195,237,184,152, 10,255, 15,  6, 52,211,250, 55,  4,255, 54,111,253,  1,  5,  7,253,254, 12, 53,  8,240, 11,248, 21,252, 32,246,255,218,  5,136,172, 13,217,210,235,231,  3,254,  0,  0,  2,242,  0,  0, 12,  1,  0,  1,193,  1,255,  1,254,245,  0,  0,  0,  2,  3,  0, 17, 12,  2,  1,248,251,  5, 10,244, 62,253,  8, 16,166,233,118,230, 27,221, 39,  6, 39, 70, 82,243, 77,177, 66,  1, 54,167, 58,199,114, 75, 51,239,  8,208,  9,  7,229,200,230, 69, 28, 16, 40,  4, 87,252,  0,226,208,225, 95,  6,255, 87,  8,  1,  3,  1,138,  1,249,255,  2, 51,  7,252,103,  7,  2,255,253,253,255,192,248,  2,249,254,250,245, 43,227,251,  1,  1,255,255,  1,  2,  1, 26,  1,  0,  0,  0,255,  0,254,  2,
      254,246, 10, 38,247, 20, 38,201, 42,244,235, 61,236,  6,238, 81,184,  9,239, 15,239,151, 25,249, 13,238,208,202,121, 35,225,201,221,252,  2,  6,222,  4,253,252, 31,229,  3, 25, 21,  0,250, 63,  3, 54,253, 24,  9,237, 28, 47,193, 33, 34, 96,241, 91, 38, 67,244,  7,  0,  1,254,255,  0,  0,253,  1,  1,255,  6,249,  1,  9,252,  2,252,  2,247,255,  3,  1,248,254,254,  3,245,  2,253,242,253,255, 50,194,  1,189,255,  2, 25,255, 55,226, 34,  0,105, 32, 54,245,255, 80,248,223,  9,254,255,250, 36,101, 36, 37,  0, 15,214,237,  1,255,252,239,241,255,  0,187,241, 87,232, 65,124, 11, 58,188, 38,  0, 63, 68,  8, 62,  0,  1,  3,254,  8,250,254,248,254,  0,255,  4,  0,254, 30,251, 14,254,252,235,  8,255, 46, 34,254,  0,  1,255,  3,  1,  9,254,  0,255,  1,  1,  7,  0,254,  2,
        0, 44, 11,222,  3,244,243,240,193, 55, 22,195, 33, 30, 52,192, 94,  3, 38,203, 17,227,152, 56, 60,229, 72, 18, 54,195,244, 70, 63, 59,253,248,  4,253,251,  2,242,227, 54,255,  6, 33, 13,240,182,247, 69,196,246,  9, 41,  8,249,216,  7,  7,167,183, 67,164,  2,  0,  0,  2,  2,249,  1,255, 44,  1,254,  0, 50,  0,  2,  1,248,224,  3,  2,  2,239,253,  1,  3,255,246,  0,246,127,  0,231,234,248, 94,231, 46, 81,159,150,190,251,241, 31,244,180,224, 13,  0,212, 67,186, 63, 21, 67,186, 17,  8, 21,254,239,  1, 34,  2,255, 20,190,239, 52, 85,227,  8,  1,223,222, 88, 89,253,238, 16, 65, 11,255, 19,247,189, 20, 54,  0,  6, 13,254,255,  1, 50,  4,255,  0,248,237,221,239,217,  0,246,  0,  0,  4,  0,242,211,254,254,  0,  0,  3,  1, 45, 10, 34,255,  0,  0,  1,  2,255,255,  0,
    },
    {
      248, 39,254,128, 55, 92,255,255,203, 27,  6, 19,176,254,  0, 52,203, 29,228, 27, 42, 28, 13, 64, 38, 51,  3,240,221, 14,113,  8, 25, 17, 11,  2,224, 42,  4, 24,220,250,249,251,208,251,  2,  9, 28, 20, 75,245, 44,245,255,  9,244, 63, 53,  7,185,196,119,255,  6,  0,  1,255,249,  1,  0,  1,253,  0,252,  4, 20,255,  0,  1, 14,  6, 25,254,  9, 27,197,  6,249,251, 24,  4,  5,123,129,221, 75,230, 16,  8,  1,226, 16,250, 14, 12, 59, 33, 16,232,193,183, 18, 21,  1, 20,235, 58,  1, 31,  9, 18, 42,218,193,  6,251, 36, 18,222,253,  6,203,210, 52,119,  3, 16,  2,  0,211,209,234, 57,239, 26, 62,  4,  2, 47,  4, 46,196,245, 17, 50,250,254,234,201,  3,  3,253,252,  7,  6,  6,109,254,  2,  0,  0,229,252,252,247,  1,255, 13,  2,255,  0,241,244,240,252,  1, 17,  0,255,  6,  2,
        1, 33,  6,  7,249, 42, 36,198,  8,245, 44,193, 19,189, 34, 94, 25, 48, 45,  4,  0, 18,  0,252, 11,244, 87,  5, 19, 56, 98,141, 25, 42, 67,248,226,  0,  2,  5, 43,255,195,253, 50,  0, 16,124,211,221,250,244,  8,242, 66,  6,  3,  0,180,251, 61, 69, 21, 58,  0,253,204,255,254,  1,254,  1,  0,255,  4,255,  9,  1,  1,253,  1,  8,218,  2,253,247, 64,253,  9,248,246,233,254,132,  2,  0,252,  5,250,179,  0, 95,225,179,  5,242,  0,251,254,  9, 48,220,  3, 44, 26,239, 14,206,255,254,  8,238, 10,253,233, 51, 44,188,215, 44,227, 10, 49, 40,  4,175,  0,254, 39, 15,  8,238, 45,183,202,232,  1, 10, 66,238,245, 10, 13,  3, 29, 51,  2, 10,240, 47,254,253, 35,255, 41,  7,  6,  9,  0,253,  7,253, 28,  3,  1, 17,247,  2,251,255, 12,254,254,254, 27,249,255,  1,254,  0,252,  3,
       17, 17,170,245,198,246,220, 61,241,240,218, 16, 32,199, 16,  2, 23,248, 65,230,202,245,245, 23,219, 11, 38,244,250, 44,248,139,  3, 18,176,233,242, 13,242,205,249,  8, 32,168,184,  3,  0, 52, 10, 10, 56, 35, 24, 28,252,246,  3,189, 73,238,203,150,212,127,  3,  1, 10,  1,255,  2,  0,251,  0,  0,253,255, 80,255,255,  3,  7,224, 30,240, 13,244,250,  0,244,222,244,  0,239,249,250,145,210, 89,225, 57,255,197,255,247, 13,  1, 28,235,  0, 25,184,240, 44,241, 19,222, 43,233,252,170,239,  0,239,245,211,120,152, 10,241, 39, 15,  4,249,199, 69,237,251,253, 36,228, 68, 71,243,198, 39, 86, 48,208,222,  9,196,250,  0,245, 36, 31,230,241,246,220,252,  1,203,245,  4, 10,  1,  2,251,  3,254,  2,254,  6,  0,250,  3,255, 33,  1, 13,253,  4,  7,253,243,252,255, 11, 11,240,  6,
      205,  5,  7, 13,250,254,  1,  4, 57, 51, 56,132,208,240, 26,188, 62, 10,210,240, 11, 32,  0, 34, 20,253,248,254, 40, 41,191, 14,  0,248,223, 77, 27,255,  1, 50,253,  6, 23,225,220, 55, 34,105,254, 19,  0,  5,242,244,  9,  6,251, 67, 87, 65,  4,244,221, 90, 18,  0,255,  1,230,255,255,254,  1,255,  9,  0,243,255,245,255, 39,  3,  1,251,245,  4, 10, 15,  5,252,249,189,244, 63,  1,  6, 19, 17, 30,213, 53, 11, 15, 48,252, 35,249,  5, 53, 13,243, 88, 13, 13,  8, 45, 38,222, 33, 56,219,247,252,179,  2,198,  0,233, 30,  9, 18, 36, 49, 11,222,161,250, 58,245,246,188, 52,131,219,187,255,229, 22,202,209,193,131, 11,195,230,251, 31,254,226,241,  0,244, 13,251,250,186, 56,254,  4,255,  2,249,  0,  7,247,253,255,254,222,  3,  3, 34,  7,173,255,  5,255,  0,255,  1,  1,  9,
    },
    {
      198,  5,255, 47,243,  5,  6,  6,  6,215, 40,254,247, 29,252, 45,  8, 23,247,236,  5, 43, 25,246,246, 71, 69, 54,225, 89, 88,114,  7,251,247,240,221,  2,  9,248, 15,247,  1,  1,  4,209,255,243,254,212,  4,  4,132,237, 22,  2,239,250, 12, 34,239,  2, 41, 47,  5,255,254,  2,  0,  4,249,254,252, 53,  3,255,  3,  8,  0,  0,  2,  3,  4,  0,254,249,249,  0,251,230,  1, 13,255,253, 20,  6,  3,244,220,184, 12, 40,207, 34,240,238,221,193, 46,  0,171,255, 84,220, 25,238,241,198, 95, 22,231,236,246,240,242,248, 73,253, 11,  0,  1,  5,249,  2,  1, 11, 23,  1,254,170,229, 51, 67,255, 69, 53, 57, 55, 56, 36,216, 51,  1,  3,  4,  1,  1,  6,  1, 48,  0,  0,  0,  0,  3,254,  0,255,  6,  0,  2,  6,254,249,  8,  0,249, 11,201,222, 44,254,251,254,253,  0,  0,  0,  0,254,254,  8,
       74,240,247,246,241,255,234, 27,252,194,  7,232, 72,203, 10,193,  8,215,144, 38,  4,255, 53,  7,251, 51, 13, 52, 98, 79,202,172,  8,  0, 19, 18, 26,  4,  4, 15,246, 12,  0,  1,  3,229,254,253, 71,191, 14, 40, 66, 45,243, 70,  1,  8,255,254,252,194,235,252,  3,250,  6,253, 12,  4,255,255, 22,206,  0,  0,252,  9,  0,255,247,254,253, 25,240,253,  1,255,243, 19,  0,  2,247,  1, 20,255,253,234, 49,227,250, 10,254,246, 25,206,219,  6,253,250, 71,188,219,211,186,213,209,238,230, 65,  2,  0,  0,  0,  1,  0,228, 16,  0,  4,218,217,  1,254,253, 21,228, 24,252,114, 35,161,  7,  0,125,197, 38,246,204,  3, 45, 45,248,  9,  1,241, 60,243, 18,  8,  0,  1,251,244,250,  0,255,255,253,  0,  0,  2,255,  1,254,  0,241,  3,211,221,  0, 11, 79,248,  0,  0,  0,255,  3,  0,238,255,
      209, 55,224,221,  3,239, 21,236,244,  2,254,193, 15,202, 10,131, 73,252, 68,134,  0, 48, 23, 27, 43,147,234, 35,176,184,184, 80,255,  1,240,244, 13,  9,  1, 11,244, 56,255,  0,  0, 10,  0,  1,149, 30,189,254,239, 87, 72, 15,254,251, 52,  6,209,253,202,131,247,247,254,  1,  5,249,  2,252,222, 33,255,  0, 23,  1,  0,  0,  0, 32,255,248,240,253,  1,  7,  6, 66,  2,  0,  1,  4,250, 40,  3,252,252,  2,253,252,246,158,214,249, 36, 16,170, 12,193,255,  6,  1,248,240,253,  0,246, 26,250,254,  0,  0,198,  4,  1,193,  4,  2, 50,233,250,250,254,255, 21,245,253,  9,234,243,186,  8,204,169,203,149, 37,197, 20,250,  7,243, 13,254, 25,151,  1,201,  0,  0, 47,  8,  0,253,  1,  6,  1,254,255,  0,254,  0,  0,  4,  2,250, 23,247,  0, 11, 62,188,  2,254,255,  0,245,255,255,  1,
       36,250, 33,226, 35, 31,215,233,244, 26,221,251, 31,  0, 13,  2,217,  6,192, 64, 54,184,109,  2, 42,251,  0, 58, 81, 65,233, 15,  5,  2,  0,253, 23,247,245,255, 22,  5,  3,255, 44,195,  1,253,  0, 21,254,193, 17,232,241,190,  3, 39, 50,  1,  5,119,  3, 74,  7, 16,254,  1,248, 11,  4,253,239,237,  0,255, 14,  5,255,  0,191,249,  1,  1,249,  0,252,255,  3,  4,  0,  1,  2,  2,252,  2, 13,211,214,251, 26,250,  4, 57,252, 11, 26,215,241,251,241,251,155,222,  5,248,248,210, 63,210,  1,255,  2,255, 60,255,151,  6,255, 16, 35, 71,241,223,254,  0,  0,  1,243, 26, 19, 62, 50,221, 42, 62, 13, 29,253,187,195, 59,  0,255,255,241, 39, 10,  0,  0,  0,  3,244,  5,245,  7,  0,  0,254,  0,  1,252,  3,  1,244,  0,  1,  2,254,246,  9,251,254,  1,  4,  0,  1,255,255,  3,  0,  7,
    },
  },
  {  //   3
    {
       33, 20, 60,  7, 41, 17, 19,  7,  6,  6, 59, 65,213,236,251, 67, 33, 45, 67, 72,253,253,  1,  5, 49,180, 66,219, 19,  2,  7, 48,  7,205,  4,  1, 34,223,248,254, 12, 34,  1,245, 87,  8,  2, 13,  3,221,251,248, 14,197, 57,249, 79, 62,202,250, 40,199, 64,  9,255, 23,255,  0,  0,  8,  0,255,  0,251,255,  4,244,255,  0,  0,  5, 37,253,255,125, 15, 28,255,  2, 28,  3,  0,254, 19, 16,  4,243,183, 38,235, 21,  0,230, 32,229,251, 20,216, 47,223, 67,248,233, 26,226,249, 22,  3, 61,  9,193, 79, 38,164, 81, 13,176,  5, 11, 27, 20, 11,199,  3, 69,192, 45, 58, 43,  1, 57,198, 17, 60,  9,254, 19,  4,  1,182, 77,254,  8,217,252,243,  6,248, 42, 48,255,  2,  3,  2, 21,254,251, 16,250,254,  3,253,250,  1,255,250,255,  0,  0,255,255,228,  3,  8,255,  0,  2,  0,  0,255,  8,247,
       18,  8,252,193,219,236,228,101, 25, 24,  3,  3,219,237,255,187,  1,  0,  2,243, 16,254,250,246,  0, 15, 63,255,251,253,127,133,  7, 47,251,255,  8, 93,255,  2,229,230, 18,254,224,244,  0,253,253, 41,250, 46,  9, 28,  6, 56,  9,180,219,205, 52, 78,222,224,255,252,  0,  0,252,  1,255,  1,255,  3,  1,  0,248,252,  0,254,252,  1,  0,250,132,  2,250,252, 15,  0, 35,254,  0,253,  1,  0, 43,241, 17, 42,254,244, 63,253, 60, 49,235,  4,162, 80,199,120, 30,240, 50,253,212,213, 23,230, 95,170,218,187,246,149, 42,154,  4, 19, 35,209, 62,254,198,202,247,  1,190,  2,199,220, 65,127, 10,240,247,255, 73, 82, 29,  3,249,254, 88,  0,  0,  4, 18, 62,255,255,220, 20,254,255,  1,  4,  8,  4,253,  3,254,  5,253,190,  0,  0,255,  2,240,253,  0,250,  2,  1,  6,  0,  0,  0,  0,  1,
       47,246,249, 19, 36, 10, 41,198, 36,244, 11,193,212,244,224,251,240, 13, 66,253, 12,250,161, 74,211,236, 46,209, 51,  0,152,247, 33,  0,  0,  2,249,222,  2,255, 18, 48,  3, 22,214, 16, 16,  3,  9,193,194,243,252, 80,221,225, 44, 19, 49,235,  5, 69, 52,126,253,  0,  1,  0,  0,254,  0,  0,255,254,  0,254, 70,254,  0,  2,  0,255,  0,  9,  1,255,251,254,255, 15,  3,252,  0,244,  1,213,  8, 59,237, 32, 11, 15,250,255,  3, 40, 29, 50, 53, 56, 16,144,214,211,249,223,  8,237,185,  3,178,247,  7, 50,  9,219,222, 52,184,253,204, 40,  2, 16,134, 20,254,242, 43,240,187,207,252, 52,  0,255,  2,255, 77,181,218,124,253,  0,212, 24,254,240,253, 53,253,249, 11,230,  1,  1,253,  1,  0,254, 33, 26,  3,254,149,  9,255,  1,  3,  0,  0,252, 65,  3,  0,  0,253,  1,252,255,252,254,
       18,  0, 45, 67,228, 15,243,  5, 20, 41,  1, 73, 37, 67,224, 48,  1,  8,172, 83,  6,  8,128,  8,173, 62,242,193,251,214, 17, 47,214,253,254,  3,253,  1,  0,241, 26,239,  0, 13,  6, 41,250, 73,227, 86,179, 58,  1, 71,  1,229,147,246,  6,  4, 18, 70,254,169,  4,  3,255,  1,  3,255,  0,  2,210,  2,  1,  1,225,254,255,  7,  2,255,  0,  0,255, 20,  0,  2,  3,  5,214,  1, 19,  5,252,243,217,215,247, 12,140,221,  4,141,243,252,158,224,214, 15,157,186,255,  6, 17, 33, 10, 55, 67,130,232,195, 27,252,210,245,134,218,249,  4, 21, 31, 31,209,215,  3,  3, 10,251,  6,225,244,132,  0,  1,  8,  2,  5, 29,232, 73,252,  4,252,  0,254,247,  8, 57, 63,254,234, 18,220,  4, 13,  3, 64,253,  0,228, 33,  9,  1,180,138,255,  0,  2, 16,  6,247,249,254,254,  1,254,  0,255,254,  1,  1,
    },
    {
      251,  2, 28,146, 20,  1, 62,148, 26,  0, 21, 44,234,  7,  3,147,240,  7, 10,148, 18,216,  5, 65, 63, 34,251, 70,245,240, 25, 27, 26,  7, 34,237,252, 30,194,254,  8, 16,235,  7,240,239,  7, 31,  8,235,236,  2,  6,237, 94, 33,255,248,242, 18, 98, 29, 60,154,  1,  3,253,255,  0,  9,  2,254,  1,  4,254,  3, 32,  0,254,  0,242,206,250,  9, 23,242, 26,234, 86,  7,  2,251, 12, 39,250,239, 37, 10,240,245, 12,  3,  4,253,  3,  8,219,  6,224,  8, 13,248,  7,243, 14,204,252, 57, 16, 45, 62, 45, 60,  7,204,150, 49,163,  3, 11,  7, 18, 12, 20,224,245,188, 22, 60,228,254,  0,239, 10, 50,251,198, 29, 69,  8,  0,196,  7,247,232,252,  9,254,246,136,  0,254,  2,  0,  4,  1,  7,  0, 17,255,  1,246,251,252,246,  8,249,255, 15, 27, 26,224, 10,251,  2,  2,252,  2, 13,220, 36,246,
        9, 44,251,235,  6, 20, 60,124,206, 12,226, 10,253,232,249, 85,185,251,178,199,  2,236,192,253,240,222,243, 68,205, 19,159,237,192, 64,243, 66,248,232,225,  9,249, 14, 96,190,  6, 35,220,  9,224,  8, 32,237,  0,  8,245, 39,247,  1, 43,198, 77,183, 33,214,249,  0,255,  0,  3,254,  0,255,  0,  1,253, 12,255,  0,248,  0,  7,247, 16,223,217,  0,224,238,  9,218,236, 11, 27, 93, 37, 16,  2,240,252,  7,254,244,  3,  6, 11, 50, 62,237, 73,202,220, 15,  1,  4,  6,  7, 62, 34,228,219,  3, 31,  2,226, 41, 42,200, 11,  6,  8,  0,250, 11,253,249,  3, 63, 15, 25,221,204,245, 59, 58,249, 59,250,248,  4,230, 41, 52,201, 31,  4,  9, 33, 12, 18, 22,251,255,  0,  0,248,  5,253,  7,  0,  4,  4, 72,225,  4,233,242,255,253,234,244, 24,243,255, 67, 43,251,  0,  2,253,191,  3,231,
        3,218,218,236, 22, 27,  6,197, 22,255, 32,247,232,244,  6,  3, 46, 13, 18,  0,  6,230, 16, 11,196,233,247,173, 55,  7,115,253, 59,193,253, 19, 40, 11,  4, 62,230, 60,206, 48,237,223, 91,244,191,235, 32,  2,190, 44,252,201,255, 21, 66, 54,  7,  2, 56, 71,245, 10,  1,249,252,  3,  1,254,  2,  9,  0,  3,241,254,  1,  1, 22,222, 42, 20,254, 14,193,242,221,  5,  7,  2,227,254,206,255,254,254,251,254,  1,  6,  7,  6, 32, 53, 20,242, 52,221, 70, 16,255,249,219, 28,199,228, 15, 41,252,238,250, 31,245,217,226, 28,252,253, 14,252,254,255,  4,  9,255, 26, 13, 53,246,212,  8,239, 52,244,252,  9,215, 15,190,198,209,250,206,228, 92,203,228, 42,  0,  0,254,255,  3,254,  5,  1,255,252,  4,255,212,169, 73,244,253,252,  7,244,  0,190,246,231,254, 12,  1, 12, 56, 57,254,  0,
      236,219, 22,197, 24,232,234,100,220, 29,249, 67,254,247, 19, 35,250,  3,  3,126,  1, 29, 26, 38,200,232,185, 25,191, 51, 72,222,237,254,244,135, 12, 17,217,254, 35,194,252, 61, 18, 30,253, 23, 24, 30,250,165, 38,213, 21,190,251,248, 30,252,224, 57,229,209,204,  9,  4,  1,196,253,254,  0,  9,  2,248,  0,197,  2,253,254,219,  0,217, 12,226,  8,251,  4, 15,254,254,  4,223, 37, 25, 29,  3,  7, 19,  4,  1, 62,  9,237,253, 24,216, 37,248, 85,181,205, 25,241,  6, 39, 51,211, 42,107,253,235,230,229,237,  0,  9,254,  1,244,  4, 58,  0, 47,  2, 41,240,237,244,  8,236, 61, 87, 49, 55,223,251, 31, 15,234,  4,217, 31, 47,247,200, 25, 54,242,192,  4,192,  1,251,255, 60, 59,231,143,  1,252, 30,254, 14,204,138,253,  1, 12,255,195, 39,188,167,216,250,236,  4,216,233,224,244,
    },
    {
       18,215,231,199,219,205,252,227,231,210,247,  3,220,218,  2,192,253,  8,244,170,250,195, 39,201,242,248,237,127,138,141,193,142,195,212,255,240,204, 41,253,238,197,216,251,255,218,251,255,245,207,209,254,  2,188,248,251,252,191, 50,  9,  1,218,193,219, 42,248,254,  0,  1,  8,  1,255,  0, 10,  2,  0,  0, 19,  0,  1,  2, 13,  6,  0,  1,  4,248,  0,  1, 16, 10,  4,  0,132, 23,255,253,225,236, 31,  8, 35, 24, 32,  5,226,230,  3,  0, 10,242,246,246, 70,172, 42,210, 31, 30, 53,107, 15,189,  3,197,254,216,165,133,208,222,188,206,224,229,193,194,209,199,193,194,206,203, 65, 64, 10,254,252,246,252,252,225,  1,254,255,  7,  4,254,246,215,190,208,203,200,204,195,196,201,200,216,199,230,  4,210,180,180, 25,192,194,194,192,196,193,191,192,195,195,191,195,192,192,193,194,
      242,245,222,  9,253, 18,238,172,248,194, 27,164,245,  6, 11,191,  0,238, 29, 10,251,204, 10,214,245,194,  4,198,247,255, 13,224, 11,235,249,244,  3,248, 12,255,247,218,  2, 12,  8,181, 72, 36,  5,241,  1,252, 25,238,254, 10,  7,226,  0,225,  2,179,  2, 45,  3,254,239,  0,253,227,238,244,254,254,  1,  1,  0,249,  0,  2,231,240,238,255, 12,241,239,244,  1,241,241,241,248,241,242,239,209,236,216,  3,206,204,215,  6, 40,240,215, 12, 11,231,241,210,243,230,  5,  0,250,235,  4,  1, 70,246,228,  8,244,247,228, 10,254,254,254,255,161,251,  0,255,253, 14,  5,  0,  0, 18,  1,  0,255,254,  2,  1,  3,252,  2,  0,248, 12,254,248,238, 11,253,249, 27,245,235,245,  2,246,236,251,242,251,221,252,244,251,223,252,253,253,  1,254,  4,255,  0,  0,  0,253,  1,253,  0,255,  0,  0,
      238,203, 21,196,249,183,203,182,247,253,111, 17, 35, 51, 37, 34,192,255,187,  0,247,  1,  1,  1,232, 86,238, 82,229,253,229,253,  3,255,247,  0,208,254,197,  0, 27, 42,254,  9,248,252,246,254,240,254,238,  0,201,194,180,190,231,230,234,230,250,235,250,237,192,152,197,211,204,204,191,208,184,217,192,194,162,197,194,192,206,223,210,193,210,205,174,206,207,206,207,207,207,206,206,208,247,247,  0,181,  0,226, 46,234,154, 19, 95,111,243, 72,234,187,255,255,  6,  0,246,180,190,216,  6,244,249, 15,230,207,225,247, 36,251,137,207,191,  7,192,  6,198,191,255, 64,254,252,226,197,166,252,  9,250,218,  1,  1,254,  0,  3,255, 15,255,255,226,236,  4, 19,254,255,255,252,  2, 16,  6,188,  7,254,198,249,253,  0,  8,255,255,  3,253,255,255,  0,  0,252,  1,  0,  0,  0,  0,255,
      249,212,255,253, 71,255,  5,205,188,135,191,  3,  1,253,  1,194,244,247,233,  8,191,195,214,178, 14,211,230,153,220,204,193,251,212,196,217,  1,234,  1, 70, 26,201,187, 13,  1,188,159,  2, 11,253,191, 71, 25,254, 11,254, 19,199,228,231,190,149,203,224,207,157,  1,  4,248,  2, 12,254,  0,254,255,  1,  3,  1,  4,255, 16,  6,255,255,  5,  2,  0,  0,254,  0,  1,  0,  4,255,  0,  2,  1,249,  6,225,192,255,242,  3,252,143,233,231,236,202,193,198,151,  0,195,235,244, 15, 64,  2,217, 58, 99,235,237,218,187,254,201,251, 87,215,  1,  6,197,255,197,  1,232,239,187,197,248,204, 70,252,126,248, 26,  8,255,241,247,126, 62,195,208,  2,  2,227,  6,  1, 34,251,255,  1,  0,255,254, 77, 75,255,255,245,252,255,  0,  0,  0,222,248,  1,  0,  1,238,  0,  2,  7,220,  1,255,251,192,
    },
  },
  {  //   4
    {
      208, 25,231,  2,224,237,216,210,251, 10,246,189,255, 20,250, 50,248, 35,250,  7,196,218,243,191,  9,  9,  0,246, 15, 50, 52,  2,  6,254,196,  1, 43,  0,  3,250,245,241,  0,  8,250,227,253,247,  2,234, 64,236, 23, 15, 61,244, 12,225, 71,  4,  6,238, 82,248,  0,  0,251,  1,  0,  1,  0,254,  5,  4,  0,255,  6,  4,255,  1,  0,  4, 15,  0,  1,255, 23,  3,252,231, 30,255,  0,  2, 21,253,242,251, 31,246,255,232, 36,220,244, 22, 36, 64, 41, 25,222,  2,  4,208, 10, 15, 87,246, 49,254,242,242,234,242, 44,172,187,  8, 11,255,224,  0,  8,247,  0,229, 30,246,251, 41,213,218,248,246,211,195, 46, 55, 69, 64,117, 70,  1,253,  3,231, 82, 54,245,245,  4,  0,  3,  8,  0,  1,  3,252,255,  0,  2,  0,254,238,251,  4, 16,255,243, 10,  0, 56, 27,191,254,253,253,254,250,  0,  6,251,
       58,233,225,152,245,  6,212,151,250,227, 11,  4, 46,201,  2,192,199,180,  2,252, 45, 17, 65,240,236,219,196,200,247,227,192,  1,  1,  9, 78,247,226,252, 66,  1,252,245, 62,238,255,240,190, 11,248, 10, 50,243,247,222,216, 82,249,255, 19, 85, 44, 94,212,247,  0,254,  1,  2,  2,253,  1,  1,  0,  1,  7,  0,  1,  0,252,  0,235,  2,  0,  1,  4,251,247,  6,251,255,  1,  0,241,244,  6,  5,249,237,240,185, 32,  6,241, 20, 48,237,  6, 24,175, 62,237,  0,199,251,236,229,173,  4,142,128,253,  1,  0,254,231,  1, 11, 13,254,253, 75,197,249,  3,246,250, 36,  5,223,230,254,  7,196, 14,223, 43, 34,255,  4,177,217,244,  0,251, 18,255, 65,219, 30,  4,255, 61,248,  5,253,255,  6,  2,254,  7,254,255,252,255, 43,  1,217,253,249, 29, 60,  3,  2,227,  0,255,  0,255,215, 40, 16, 23,
        2,227, 46, 10, 21, 77, 71,  5,228, 16, 50,  2,241, 32,196,  0,  0,  7,  4,248, 34, 90,143,253,  2,  0,  1, 70,117,217, 50,211,  0, 11, 35, 96,246, 16, 58,152,238,245,134,  1, 16,  2,  7, 17,235,  2,  3, 64,  8,  0,178,237, 14,  1, 83,200, 70, 73,197,230,253,  0,  1,  0,255,255, 25,  0,  0,  6,245,  0,255,  1,252,  0,235,  0,255,  1, 14,  0,248,251,  0,255,  0,  1,  9,  7,  1,  3,  5,  3, 55,248,  0,192,252,226,205,216,221, 51, 61,135,201,  7,  9,255,245,252,  0, 60, 60,222,  1, 10,250,253,226,180, 37,255, 10, 20, 18, 53,254,  1,251,253,236, 10, 67, 16, 54,  4,238,  2,  1,240,238,237,253, 51, 13,180,  2,254,230,  2,  0,125, 19,217,  1,201, 29,237,255,254,  3,  2,248,251,  2,  0, 24,  3,  0,  0,246,  0, 47,254,  9,228, 65, 46,  3,  4,  3,254,250,233,  1,225,
      255, 49,239, 45, 37, 15,192,  2,  3,209, 16,  2,245,213, 35,  2, 69,213,249,241, 41, 32,  1,208,248,243, 44,253,147, 32,  9,207,  1,251,254,255,  4,252,193,253, 30,255, 29,250, 24,238, 57, 18,  6,253,  8,154,230,244,241,  4, 14,255, 29,181,221,226,206, 55,  3,  2,  2,  2,  3,  0,  1,  0,255,  1,255,  1,243,232,255,  2, 58,  4,254,  6, 26,  0,110,253,  8,  0,  3,  0,252,238,250,  0, 14,  1, 27, 30,229,244, 71,240,241,246, 37,226,241, 10, 84,  3,249,245,246,  5,252,246, 88,255,255,  0,  2,249,255,188, 19,234,  4, 58,247, 19,  1,254,191,190,  4,  2, 59,218, 25,  6,  7, 80, 37,  8,223,217,  1,253, 34,131,  4,  5, 58,222,  0,197,202, 26,  1,239,231, 73, 61,  1,  1, 53,  0,  8,252,  7,253,  3,252, 31,255, 20, 17,  5, 52,131, 32,255,  3,  6,  4,  6,  1,  1,144, 59,
    },
    {
        7, 39, 10, 26, 41, 70, 67, 65,  1,  0,  6, 43,196, 66, 64, 65,251,  3,  1,  4,193, 62, 64, 57,237,250,254,  2,195, 63, 57, 64, 58, 85,246,175, 53, 80,189, 23, 51, 45,246, 35, 88,251,185,232,214, 10, 52, 12, 47,222, 35, 38,212, 18,254,226, 75, 31,227, 23, 12,  0,  0,  2,  2,  8,  3,  4,255, 13,  0,252,253,  3,250,247,255, 44,245,  2,255,224,235, 19, 18,  7,246,  4, 56,210, 44,194,  3,251,250, 62, 60, 60, 15,  6, 46,233,255, 43,  4,250,189,245, 69,173,  0, 30,209, 39, 31, 41,119, 87,190, 81,  7, 29, 40, 76,  3,244,228,250,199,246, 61, 14,254, 34,  7, 20, 31, 14,146,244,202, 16, 50, 59, 14, 82,192, 13,250,225,246,235,232, 80,243,  9,254,252,  1,254, 29, 37,196,237,237,235,251,  3,  1, 39,  1,245,242,  5,  4,244,253,255,234,217,246,255,  2,  9,251,196,241, 37,
       88,230,  4,  0,252,237,  2,250,255,  8,  3,254,  3,  0,  1,  5,251,244,  2,  0,  0,  3,249,  1,253,  0,  2,  2,  3,  0,253,249,  4, 24,252,246, 35, 42,250, 41, 44,  2, 48, 16, 38, 32,137, 77, 16, 29, 59,223, 50, 60,251,241, 17, 88, 29, 15, 27, 36, 32,238,255, 11,  4, 16,  0, 21,  5, 15,  1, 10,  1,  0,251,  9,  2,  6, 35, 23, 68, 86, 10,253,103, 29, 22,219, 13,213, 19, 28, 94, 46, 25, 45,248,  6, 26, 41,250,  4, 10, 19, 77,241, 50, 16,232,  3,198,251, 85, 43,203,194,  6,105, 69,  5,251, 13, 25, 56, 82, 19, 48,231, 33, 44, 91,237, 62, 66, 21, 40, 37, 20, 28, 37,251, 50,241,  2,  6, 10, 15, 37, 44,235,226, 26, 36, 62, 46, 73,233, 60, 59,253, 23,253, 34,252,  0,250,249, 17,  1,255, 26,  6,  2, 38, 34, 38, 57, 42, 13,240, 47, 10, 19, 23, 54, 31, 56, 21, 14, 69,
        0,253,  1,  0, 29,255,  2,254,251,  1,  3,  3,  2,255,  6,254,  0,  0,  4,255,255,250,253,251,246,  2,254,255,  0,255,  2,  0,  0, 13,236,  2, 15,  3, 28,242, 54, 12, 47, 22, 22, 40, 57, 17,249,  1,  2,  1, 14,  6, 15, 13, 55,255, 17,234,  8, 31, 59,178, 54, 47, 65, 62, 63,  5, 61, 58, 53, 39, 60, 36, 63, 53, 46, 71, 62, 27, 70,  7, 45, 51,  0,  7, 61, 79, 62, 63, 58, 80, 63, 68,207, 19, 15, 56,  8, 49, 20, 40, 52,104, 22,240, 11,229,  0,  6, 56, 48,255, 71,214,249,  8,  2, 13,232,  0, 41,255,231,  1,  2,239, 11, 53,124,216,247, 51,231, 49,254, 14,254, 11,254,253, 18,224,228, 50, 32,  8,251,248,  4,  8,221, 57, 33, 29, 18,243,  2, 44,  5,104,241, 28, 33,  8, 10, 28,  4, 70, 25, 62, 41, 40, 31,196,206, 12,208,  4, 34, 45, 31, 33, 90, 65, 87, 31, 25,  0, 40,
        5,  2,254,  8,252, 10,  1,  3,  4,  4,254,253,  0,  3,  0,  3,  2,  3,  1,252,254,255, 19,  2,254,230,  9,  0,254,254,  3,  1,233, 78, 46, 43, 47, 65, 16, 23,236, 55, 56,210,208, 55,249, 23,227, 63, 38, 35, 12,  8, 27,255, 64,196, 13, 15, 25,235, 32, 33,  5,  3,  7,  6,252,253,  9,238,  3,  1, 65,227,  0,210, 35,225,  6, 14,248, 14,  0,253,  6, 10,251,  0,  2,253,  0,250, 36,248,141,131,  7,237,231, 27, 72, 26,171,117,213,236,234, 40, 65,226, 99, 46,210, 17,  3,235, 43, 64,139,  1,245, 85, 36, 39,212, 21, 19, 23,250,198, 33, 31, 18, 14,104,212, 16, 14, 38, 55, 65, 66,201,  8,236, 60, 18, 65, 42,248, 98,201,  0,217, 17,184,107, 59,234, 93,230,183,246,  4, 20,  6,105, 21,231, 15,254,251,  1, 14, 37,242,  4,237, 14, 23,250,  7,224, 22,214,  1, 33, 66, 42, 63,
    },
    {
       52, 63, 37, 28,206, 45,  5,  0, 11,  9,  9,  4,  2,253, 11,  1,221,224,229,241,220,220,254,203,228,224, 20, 73,234,223, 43,207,  0, 24, 55,250, 13,  8,  1,255,  6, 29,  1,181, 18,254,  0,247,  9, 19,242,239, 12, 51,252, 83, 11, 23,  2,191, 15, 14,  0, 54,  0,  3,248,  1,  7,  1,255,  1,  2,  0,  2,253,  4,  1,  2,255,  2,  5,  2,253,  4,255,  0,249,  5,  2,  1, 24,  4,  5,  3,228,  7,218,231,216,193,248,  8,238,229,237,209, 54,231,249,189,195, 52, 18,240, 20, 42, 52, 59,231,235,218,230,240,237, 41, 29,117,249, 20,  4,243, 32,198,  1,177,218,251, 13,  6,201, 15,150,202,231, 62,186, 22,224,197, 30,234,  4, 12,  9,  1,  0,205,251, 22,254,  2,254,  1,193,250,  1,  6, 39,  3,  3,  5,247,  1,248,248, 37,239, 16,215,  3,  2,252,  3,255,250,  2,  1,  3, 14,255,  5,
      252,216,249,  0, 61,184, 10, 32,255,249, 19,209,203,246,193, 10,155, 28, 65,249, 47,207,249, 50,  8, 31,199, 50,169,220, 51, 65,246,226,254,186, 30,  8,  5,153, 20,243,253,  0,234,215,122,237,247, 26, 46,252, 11,  8,204,205,239, 21,215,126,195,101,252,  6,  0,  0,251,255,255,253,254,  5,  0,  3,  0,  2,251,255,255,  1,  0,  4,253, 13,  5,254, 30, 10,  0,  3,  0,255,  1,247,  4,253,222,241, 14, 83,195,202, 23,161, 57,228, 81,252,  3,  3,216,243,198,203,209,222,  2,112,190, 24,255,254,255,  0, 68,  5, 29,  6,249,199,228,  2,183, 76,251, 61,206, 51, 16,220,199, 81,255,211, 13,233, 17,255,  0, 15,  1, 25,  1,254,  0,  6,235,245,239,193,251, 26,252,  6,  1,  8,250,  7,255,  0,255,255,  1,  1,  2,255,221,252,235, 24, 24,  8, 13,252,  0,  1,254,  1,  0,247,249,  5,
      250, 26, 12,192,243, 38, 53,203,230, 51,251,  4, 21, 15,251, 12,244,  5,  2,197,200, 50,226,194, 12,224,250,219,  3, 61,118,210,255,193,196, 33,211,  6, 20, 26, 12,210,194, 35,249, 10,181, 27, 10,255,  2,148, 20, 10,188,  2, 11,213,213,209,231,165,134,143,  0, 10,241,253,  0, 10,  0,255,  0,  1, 35,  1,255,255,  0,254,  1,  6,255,219,  0,  0, 19,251,255,254,  2,228,227,  0, 31,  2,246,226, 32,253,241,175,192,236,242,193,  8,226,252,115,255,239,239,  1,212,231, 30,178,206,255,  1,255,  5,  0, 20,  8,100,253,  7,206,200,244, 59,245,  3,213, 98,243,236,244,  7,192,  1, 41,253, 21, 25,  7,187,167, 11,110,  2, 17,237,183,196,252,247, 56,  4,203, 47,  3, 19,250,  1,253,  6,253,249,  0, 11,  1,  1, 32,242,247, 15,253,255,239,  6,  3,  0,  1,254,254,  1,255, 20,255,
      255,245,  2,187, 25,236,238,255,229,214,220,244,253,  5,193,207,208,  3,255,  1,140,254, 67,222,246, 17, 59,248,188,253,  4,  8,  1, 27, 66,254, 27, 73,190,254, 28, 21, 56, 64,237,  3,  7,215,  1,  1,198, 74,246, 71,  5,231,  2, 37, 39,191,147, 66,  3,187,  1,252,  0,  1,  7,  0,  0,  1,  5,  4,255,254,  4,255,  1,254,  1,248,  0,  4,251, 24,250,  1,  2, 21,  0,254,  0,  1,237, 26,  3,252, 17, 36,  1,253, 20,198, 23, 27,222,221,241,  7,239,188,254,  5,242,248, 31,252,204,244,  4, 16,  3,  3, 69,  1,170, 55,  3, 31,  5,246,251,156, 20, 65,202, 32,231, 15,241, 75, 18, 24,227,236, 51,227,192,160,212,240, 12, 74,241,253,252, 67,212, 61,246, 13,205, 10,  1,246, 60, 32,  2,  2,253,255,  7,  0,111,  8,247,254,254,  6, 41,252,  4,  5,255,  0,  1,  0,255,187,246,250,
    },
  },
  {  //   5
    {
       14,  3, 16,188, 37,  9, 22, 46, 58, 65,  9,221,  2,205, 21,239,  8, 22, 26,185,104,215,  7, 42, 35, 51, 20,112,219, 43,248,171, 13,186,254,237,206,226,255,228,  9, 29,  1,234, 33, 22, 11, 10,247,212,252,244, 36,  2,  6,  5,202,225,252,255,  4,190,  5,144,  7,  0,  1,  0,  4,  2,255,  0,  4,  0,  3,254,  0,  2,255,  0,244, 11,255,  0,244,253,  5,254,249,252,  1,255,142, 16,138,128,240, 26, 29,233,227, 18,253, 14, 59, 64, 22,218,212, 17, 78,  5, 10,  4, 21, 12, 17,206, 68, 43, 28, 20,241, 68,193, 36,160,149,255,250,  6,  5,  7,217, 12, 13, 18,  7,  3, 15,249,  2,244,253,191, 62, 72, 56, 57, 22,242,122,  9,  5,  9,202,201, 43,171,228,  1,  0,  0,254,  0,250,  9,252,  3,249,  1,255,  4,247,255,  6,  4,233, 33,225,  8,  2,246,250,253,  1,255, 13,  0,241,  4,  6,
      213, 76, 11,110,  9, 25,252, 49, 43,204, 12,238,  2, 77,241, 47,227,235, 39, 19,248, 26,147,239,212,  9,  0,233, 43,194, 55, 85,253,225,249,  5, 36, 35,  0,166,207, 27, 11,239, 52,247,  8,207,201, 50,178, 77,  4, 48,126,242, 36, 17,  4,204, 66,189,225,  3, 68,  2,  0,  1, 13,255,  0,255,246,252,  2,  0, 22,  0,  1,254,248,  1,227,  7, 12, 10,  1,253,225,  1,  0,255, 12,210, 26,  4,  1, 11, 21,  4,  5,183,254,  4,  4,254,227,237, 20, 18, 25, 41,241, 53, 70, 45, 15, 62,115,196, 27,237,  0,255,226,240,  8,217,250,253,248,  6,255,  2, 14, 12,248, 19, 17,151,249,245,129,224, 24,247,  5, 28, 32, 69,250, 54, 37,244, 11, 36, 24, 66,207,178,250,  1,251,  0,  0,  9,244,255,124,  1,  0,255,253,  2,  0,246,245,  5, 23, 54,254, 24,248, 50, 12,255,  1,255,  8,  1,  4, 69,
      200,228,  0,246,254,  5,  0, 65,239, 26,205, 53,232, 24,254,166, 17,243, 68,208,251,214, 72,135,174, 44, 46,111, 23,242,254,189,251,198,  2,195,  9,251,  0,143, 67,248,  5,254,220,230,253, 70, 54,255,116,110, 11,217, 37,127, 27, 27,  2,233,195, 33, 20,215, 33,254,  1,  0, 16,  0,  1,254,  3, 15,  1,255, 14,  0,255,  0,  0, 34,  1,  0,  4, 22,232,  1,252, 89,  1,  3,  1,254,  1,  0, 20, 30,204,252, 57,202,218,163,191,212, 29,215,229,254,103,180, 59, 11,239, 18,112, 15,  5,  3,  2,247,238,253, 12,  3, 72, 66,198, 49,228, 27,245, 65,187,252,255,253,  3, 35, 88, 54,127,131,232,  8,  1,210,  0,177,131, 62,  1,155, 15,231,  2,177, 93,237,  8,  0,251,  7,  0,242,196,252,  6,254,  7,252,  1,  7,  3,117,255,  1,218, 14, 35, 11,132,232,255,250, 27,254,  1,  4, 73,252,
       42,253,  0,  4, 22,251,  2,141,231,245,226,232,217,  9,254,210,247,199, 63,145,250,  9,247,227,224, 13, 18,148,197, 18, 12, 70, 10,225,  0,  6,240,  1,  4,196,  3,  4,246,  6,252, 21,  0, 67, 48,226,255,230,229,235,229, 68,247,243, 31,  2,  9,201, 20,208,  1,239,  0,  2,255,255,255,  1,  0,252,  0,255, 12,251,  3,255, 15,129,  4,255,  0,125,  0,126,186, 23,253,255,184, 62,  2,135, 17,  3, 10, 34, 57,187,130,  0,  6,255, 34, 27,242,  0, 24, 58,252,234, 10, 16, 41,163,161, 81, 20,193,254,  4,144,241,180,179, 43,248,249,  6, 30, 55, 76, 94,250,  1, 27, 29,169,253,211,159,  9,193, 18, 89,212, 59,223,137,250,  9, 44,116,230,  4,201,243,  2,  2,222, 59,  4,247,246,131, 58,254,247,  2,253,  0, 18,253,255,  3, 36,193,  1, 55, 55,  3,255,251,246,  0,253,211,  9, 89,
    },
    {
       20,241,192,107,  9,  5,  1,247,255,231,  8, 38,254, 17,  7,232,254, 11, 20, 94,192,214,246,  1,240, 43, 22,  4,  3,213, 25,243,247, 11,  2,245, 12,233,252,  7,  1,  0,  9,  9,233,203,253, 15,253,212,  8,  4,253,  0, 57,253,  5,238,237,255,177,196, 23,217,  0,254,255,  1,  0,  2,  0,255,254,  5,  0,255,  6,  0,  1,255, 19,247,255,  2,254,  4, 29,  0, 18,240,249,254,205,221,  7,251, 26,225, 35,  3,253, 69,216,135, 13,  1, 39,227,  6,  6,  3, 18, 11,191,225, 27,  3,235,250,217,  5,240,239,250, 17,194, 55,  5,  1,229,229,222,  1,211,252,107,198, 64,210,192, 23,241, 24,138,228, 71,239, 12,  4, 98,  5, 17,252,  0,  3,252, 53, 13,201,111,255,  3,  5, 37,  3,242, 35, 66, 14, 17,245,  4,  2,  3,248,230,  3,  4,  0,  2,  5,  5,252, 45,  0,  0,  1,  0,198,232,  2,252,
        1, 57, 65,237, 63,211, 53, 20, 31,  0, 20,  7,208,231,240,215,248,  9,  0, 29, 42,217, 68, 12, 30,200,201,245,  4,210,242, 12,254,201,  7,  3, 21,240,  1,  6,  9,228,242,252,230,248, 30,  4,192,190,249,225,248,253,240,248,240,121, 31,241,229,203, 28, 27,  2,  4,  0,255,252, 16,255,  1,  4,  2,255,  1,  2,  0,  0,  1, 47,  1,252, 12, 41,252,255,  1,  1,255,  7,254, 60,248,205, 53,  0,251, 71, 32,222,  0,198, 13, 20, 20,  2,227, 69, 69, 74,190,  8, 36, 18,  4,243,249, 50, 38, 24, 23,252,199, 35, 39,222, 28,251, 28,  0, 23, 52,179, 82,153, 66,253,228, 19,  2, 20,153,188,241,  9, 12,194, 54, 78,  9, 88,237,  0,  8,246,241,113, 49,142,  1, 17, 42,255,238,232, 35,  5,  1,224,254,192, 18,255,239,255,253,  0,252,  0, 15,246,245, 93, 10,  3,  0,  0,  0,  0, 20,216,
        2, 36,237,241, 41,  5,  8,246,233, 41,110, 57,195, 22, 44, 40,  3,254, 36,214,239,219,219, 42, 18, 35, 58,242,117, 21,186,153,252,100,241,250,250,253,255,  2, 35, 21,  4,252,217, 41,  6,255, 79, 63, 67,230,216, 48,233,253, 77, 65,201,254,227, 61,213,  3,  0,  0,  0,255,253,236,  0,  1,  1,255,  0,254,  5, 20,  0,  0,  1,  4,  9,  0,  0,204,246,255,112,  1,  9,  0, 47, 18, 56,  1,  7,  2,163, 47, 53,253, 30, 63, 11, 56,207,252,245,212, 61,126,  0, 23,251,208, 43,235,251,179, 26, 36,205,224, 11, 62, 15,207, 12,235,  4,250, 36, 79,221, 19, 35, 22, 18, 32,  3,  2,219,164,255,232, 12,254, 82, 49, 20,113,  4,252, 17,241,202,226, 35,142,  2,  0,  1,177,  3,  2,  5,250,250,  0,223, 39,242,251, 11,242,254,254,255,254,247,223, 12,207,  0,  3,  1,248,255, 37,216,  0,
      252,196,239, 15,211, 26,196,  7,253,  6, 55, 38,231, 10, 60,  4,  4,207,238,234,  6, 32, 76,  7,250, 50,206,247,194, 27, 66, 38,  2,223,  2,254,  4,250,  0,246,  5,189,  3,250,  8,  4,250,252, 16, 10,198,  5, 32,255,242,248, 26,217,  5,252,203, 45,  0,214,254,251,  3,255, 19,255,  0,255,  1,  4,  0,  3,  0,  8,254,255,  8,223,  1,  1,228,224,  3,255,225,185,  1,243, 22, 59,  1,240,  6, 26,243,239, 80,126,232,196,220,209,  2, 24,216, 90, 13,203, 30,  7, 41, 59,224,238,105,248, 73, 25,119,220, 69, 49,187, 42,248, 34,  5,161,204,248, 81, 94,224,253,215,255, 16,119, 28,195, 20,240,  2,  3,212,  6,120, 52,247, 29, 50,110, 14, 79,255,145,  0,250,  7, 21, 29,  7,  5, 77,  1,  1,  6,216,  6,224,121,128,  1,  0,253,251,  0,  2, 93, 18,252,255,  9,  1, 16,239,221,234,
    },
    {
       22, 32, 19,221, 26, 24, 28,210, 24, 52, 19, 83,  4, 39,228, 90,231,251,222,  9,146,201,183,225, 13,  8, 71,132,  4,  3, 26,  3,  0,  2,  5,217,  4,245,  0,228, 44, 17,  2, 13,  2,250,247,  2,242,201, 60,  1,244,235,  0,  3,  1,234, 11,253,252, 55, 27,  3,  0,250,  3,251,255,255,  1, 59, 15,253,  0,251,  1,250,255,  3,255,231,  7,255,  8, 15,  4,  0,255,190,246,244,  3,113, 15,254,  8, 14,251,255,244,233,  7,  2, 35, 12, 18, 76,206,181,229, 70,207,228,201, 14,  2,209, 19,126, 63,231,251, 69,205,  7,130,190,  1,206,  9,243,  1, 29,252,198,  2,254,238,  6,197,249,253,  3,241,249,232,  1, 21,239, 16,189, 23,249,  2,252,228,255, 73, 54,255,254,255, 19,255, 13,254,  1,254,  3,179, 30,  2,235, 11,  0,  0,255,255,255,254,  1,  1,253,  0,255,255,  0,255,  1, 21,  0,
      254,243,215, 73,  2,  8,195,  2,253,  6,250,  7,197, 31,205, 10, 71,241, 59,  3, 49,  2, 57,118, 55,113,166, 60,242, 28,  3, 69,191, 45,254, 14,  8,  4,  2,237,182, 23,  9,250, 34,  3,  4,253, 95,139,235,232,  6,217,  2,253, 15, 18,  2,255,241,177,231,242,254, 25,  0,  0,254,  7,  0,191,250,249,  0,252,244,  5,254,  6,233,255,254,255, 40,250,  0,249,237,255,  1,251,254,  2, 14,245,254,241,  1, 16, 30,239, 34,179,  0, 59, 11,182, 25, 20,213,248,  9, 96,255, 61,  9,194,251, 68,255, 70, 30, 14, 12,171, 60, 98,230,246,253, 29, 92,242,255,199,  0,211, 42,215,241, 52,227, 11, 17,236,  8,252, 11,  1, 10,169,248,252,  1,  1,  1, 22,148,  2,  0,  0,  0,255,245,  2,249,246,252,251, 61,233,200,  2,248,  2,  0,  0,  0,  1,  9,254,209,  4,252,  0,254,255,  2,  2,  6,  2,
        5,  0,  5,188,  4,217,243, 51, 19,  8, 42,180, 52,105,241,185, 31, 44, 25, 68,229,  2,  3, 16,168,177, 52,192, 41,244,241,176, 61, 29,253,  4,252,216,  1,  5, 35,225,248,255, 21, 28,226,253,194, 29,249,106, 41,230,  0, 58, 12,  7, 67, 14, 33,232,254,249,  0,218,249,254,  0, 14,254,  0,  0,  6, 10,  1,  3,235,  5,  2,222, 44,247,255,207,255,  0,255,249,222, 10,240,  0,  0,236,  6, 58,253, 18,197, 35, 23,157, 90,169, 14,235, 16,148,243,234,207, 70, 55,252, 49, 61,253,187,236,212,214, 26,116, 39,252,171, 59, 11, 15, 16,247,184, 12, 27,215,255, 60,202,  6,237,  8,121, 66,  0, 10,  3,  4,  1, 73, 32,190,210, 14,  0,  6, 10, 28, 45,241,252,  8,  1,248, 25,253,254,  0,  0,253,235,  2, 24,255,  0, 16,254,253,255,  0,253,252,250,255,  2,255,255,  3,  0,255,254,248,
        4, 24,253,253,  0, 57, 76,248,205,245, 18,253,240,202, 38, 18,  8,194,255, 70,227,100,255,154,207,158,198,142,249, 44,172, 27, 11,207,240,253,  0, 19,  5,250,  8, 17,248,  0,231, 34, 14,251,234,226,  0,248,252,  5,253,205,248, 20,  0, 48,  1,245,244, 54,  3,228,250,255,255, 20,  3,255,  8,249,251,  0,249,  1,  4,  3, 24,  4,  0,  3,255, 14,  0,255,214,  4,242,155,254,174,252, 13,196, 14,254,181,250,  5,208,142, 39,234,221, 39,  7,195, 15, 55,187, 12,239,248,206, 25,215,138,188,121, 18, 42,  1,126,236, 50,  8, 11,  7,  7,225,255, 32, 62,  4,205,  8,209,242, 24,186,177,237,  7,  6,251, 15,208, 50,218,255, 11,252,254,  0, 35,127,212,  0,249,  3, 12,  8,  0,251,  1,  6,  1,250, 34,  0, 73,248,205,  0,255,255,  0,244,  6,255,252,251,254,  2,254,255,251,  7,  5,
    },
  },
  {  //   6
    {
      227,230,227,209,202,244,183, 65,198,239, 69,243,237,235,238,241, 13, 17, 46,249,244, 59,196,  2,208,226,213,235,194,251,241, 15,229,253,254,245, 15, 41,254,  5,223,  0,  8,  3,207,221,  2,255,221,210, 16,  3,202,  3,253,  0,211, 33,245,  5,226, 25,  2,255,  9,254,  0,254,254, 10,  0,  0,253,252,255,  0,255,  3,254,  0,  2, 15,  0,  0,  1,239,  0,  0, 18,255,  0,  3, 37,  5,  4,  0,247,203, 38,241, 29,255,167,  0,190,189,187,197,171,207,127,178,180,226,253,240, 22,  0, 52, 64,236,238,232,243,223,237,245,245, 12,206, 42, 64, 12,249, 28,  4,  9, 17,251,195, 57,223,203,111,225,193, 10,  0,200, 40,  8, 38,  0,240,  9,167, 22, 65,  8, 41,  2,252, 23, 19,  4,  5,250,  2,  3,255,255,255, 24,  1,244,129,233,255,253,  4,213,  0,  0,244,254,  1,  1,255,  0,  0,  0,250,
        5, 11, 30,229, 45,247, 68, 15,225, 53,186,239, 66, 36,204,191,251,239,181,245,194, 63,128, 68,244,219,224,183, 64, 62,193,188, 56, 11,251,  7,205,178,254,  2, 21,250,  5,254,242, 29,  2,  2,  7, 16,250,  1,  6,130,255, 47,239,  1,254,  4,  5,181,  2,  1,254,250,255,255, 19,254,  0,  0,255,253,  0,  0,255,  3,  0,  0,221, 14,255,  0,  5,  3,255,  0,250, 26,  0,  1,181,210,255,  2, 42,235,253,249, 74,166,189,105,247,250,229, 21, 70, 44,106,197, 19,255,255, 15,  7,  4,122,228,  6,255,247,255,246, 10, 87, 52, 37, 42,227, 59,210,188,151,  2,  3,241, 31,  3,255,240, 34,135, 63,255,  3,196,101,247, 51, 49, 30,253, 87,253,243,242,155, 51,243,  5,  6,  3,246,255,255,255,248, 24,  1,252,  1,  1,  1,  4,255, 38,  9,244,253,247,243,  1,255,  1,255,255, 14,250, 21,  5,
      251,241, 33, 25, 46,230, 62,175, 31,251, 19, 68,251,210, 68, 14,230,235, 48,158, 18, 43, 82,175,236, 20,215, 93,225, 16, 89,243,212,  5,  0,  4, 52, 73,  7,252, 22, 11,255,  0,  9,246,254,253,194,100,253,  1, 16,189, 31,  1,109,203,  2,  7,222, 91,254,  0,  0,242,254,255,255,254,254,255, 17,  7,255,  0,  0,255,  1,  0,  5,  1,  1,  0,255,  1,  2,  0,189,  8,253,255, 77,158,  5,  6, 17, 43,253, 38,205,240,147, 93, 23,  6, 26,252, 64,198,210,  6,193,242,202,212,242,194,255,132,247,254,  2,250,  1,235,205,167,215, 36, 30,188,  1,  7, 43,  8, 32, 28, 78, 67,202,184,249,245,250,228, 64,240, 29, 20,245,119,252,  7,254,252,165,250,101,  4,255,255,  4,  9,  8,128,  5,252,  3,251,241,  0,254,229,253,253,  6,254,242,250,235,117,237,232,  0,253,255,  1,  3,255,231,231,
       27,239, 41, 44, 66, 60, 19, 35,233, 55,228,178, 46,219,212,133, 60,201,246,130,176,235,108,255, 24, 34,243,111, 16,174, 19,157,  1, 58,  3,  9,206,  5,247,254,221, 18,255, 38, 36,254,254,  1,101, 43,  0,  0,254, 19, 55,  0,189, 30,  1, 10, 67, 84,  3,  1,250,255,254,255,253,  1,  0,  0,  0,  1,  0,  2, 10,  5,  0,  0,247,  8,  0,  0,252,  0,254,  0,229,254,254,254,  0,222,248,252, 27,215,211, 40,253,  4, 78,155,231,208,223, 44,251,204,134,  0, 45,206, 35, 53, 47,  4,  1,239,  6, 32,255,254,255,239,166,204,  3,216,222, 29, 11, 63,254,190,  6,221, 13,240, 26, 52, 43,255, 63,  7, 13,206, 73,225,194,  7,  1,  4,  9,252,  9,188,248,133,  1,  6,242,223, 10,126,  4,  0,255,  5,  2,  0,247,  6,191,194,  0,  3,229,  4, 15,145,  5,246,255,255,  2,248,  0,217, 15, 51,
    },
    {
       13,237,210, 68,  4,  5,252,190,241,245,250,252,191, 63,207,205,  2,254,  4,  3,187,191,195,194,251,  2,  2,  1,180, 60,195, 64,205,248, 11,  4,230, 59,239,240,198,236,255, 15,213, 32,251,248,217,215,255,  1,215,228,129,126,  0,212,253,  0,164, 79,254,246,  1,252,255,  1,  0,  5,  0,255,250,  2,255,254,  7, 10,  0,  1, 46, 11,  1,253,231,229,  3,  3,228, 31,  6,  3,242,212,252,250,252,205,  7,246,178, 56, 74, 11, 18, 50,213, 38, 28, 15,184,222,235,  1,216, 15,248, 19, 69,113, 29,169, 31, 32,249, 39,117,149, 73,193, 34, 13,191, 10,194, 10,  3,251, 18, 21, 33,  0, 33,181,185,  0,233, 31,231, 13, 18, 43,  5,239,255,241, 31,212,221,186,251, 19,254,  2,199, 89,248,  3,252,254,  0,  0,  7,252,  2,  9,  0,  1,254,  6,  1,  8,  1,  1,  0,  7,  1,255,254,250,251,  1,
      244,246,  0,  9, 15,241, 77,  1, 58,253, 95,227,  9,  2, 18,237,  2,  6,  1,250,  2,  2,250,  8,255,  2,  0,251,  1,246,252,247,  7, 14,254,  7,218, 17,250,253,  6,226, 19, 20,211,241,  5,  0,250,246,252,  5,  3,174,139,127, 22,  0,252,  7,137, 41,248,247,231,253,  0,  0,255,  4,  0,  1,252,  0,  1,  0,  2,250,  1,  0,217,227,254,255,220, 27,  0,  5,207,197,  6,253,  1, 59,  8,  0,251,247, 26,254,194, 58, 25, 39,250,239, 29,  8,237,242, 26, 20,218,  4, 46,250,177, 52,255,244,253,253,233,255,  2,187, 11,119,253, 40,212,252,226,197,205, 92,182, 10,230,213,126,236, 29,232,197, 62,194, 37, 43,234,162, 81,193, 45,208, 38, 14,126,218,110, 22,255,242,254,  0,  4,248,  8,  1,251,  0,254,214,243,255,  6,254,251,216,219,255,252,245,243,250,252,248,  0,  2,  5,121, 14,
        5,  9,203,247,  1,234,252,179,217,  7, 29,242,221, 11,230,  0,255,255,255,254,254,  2,250,  4,  0, 31,  0,  1,  0,251,255,  1,223,184,250,  1, 14, 26,254,242,247, 18,  0,245,  1, 20,  2,  1,  4,244,245, 10, 47, 74,  3,  6,255,255,  4,253,223, 59, 57,253,255,247,  0,255,  2,  1,255,  1,254,252,255,255,248,  1,254,  0,244,229,251,250,189,216,255,248,220,209,255,  0,254,254,253,252, 13,224,252,248,  3,251,  3, 12,248,  4, 42,218, 26,230, 52, 39,225, 59,252, 37,238,245,249,250,  0, 31, 18,254,135,  0, 65,233,202, 49,206, 87,208, 56,183, 35, 59,229,  0, 10, 61, 84,  0, 63,251,174,  4,  0, 35, 65, 24,206,240, 40,  4,238,194, 75, 14,155,222,224,227,226,208,  8,222, 14,192,  3,192,  1,  1,254,192, 31,228,  1,239,232,252,  1,252, 15,207,255,228,  7,250, 66,220, 12,
      249,  0,  1,  1,  1,  0, 78,  5,  2,247,211,242, 28,255, 13, 38, 36,  2,250,253,254,  0,253,252,  0,238,  3,254,  0,  0,  6,246,  1, 10,252,  7,211, 32,239,  2, 27, 34,  1,  4,  6, 52,  3,  0,  1,237,251,  0,  1,209,190,231, 26,252, 40, 10,  4,191, 13, 17,196,  4,193,  2,254,  2,196,254,187,  2,191,  1,180,255,194,  0,198,251,195,255,193,  9,192,  2,195,158,190,  1,194,193,195,  1,  8, 10, 81,  8,  1,151,198, 22,234, 27,252,193,208,218,195,254,249,254, 57, 14,241,  3,193,138,  7,231,  4,  4, 29, 52,187,132, 14,237,143,199,231, 60,254,244,  2,234,233, 29,245, 65, 13,  1, 55,252, 61,226,  4,248,192,133, 55, 89,242,  7, 16, 22,242,247,  6, 40,  2,210,222, 43,247,239,227,  0, 38,  3,  1,248,222,216, 22,255,201, 13,184,252,  1,249,  1,  4,223,  0,195,  6,222, 11,
    },
    {
       10,  7,191,212,254,  5,194,233,246,  8,118,240,230,  6,234,246,225,200,164,179,214,251,208,112,214, 17, 56,221, 94, 30, 49, 95,246,255,241,241,  4,247,241,250,243,235,  9,  2,236,249,245,245, 17, 17,254,246, 14,217,  6,  1,243,195, 15,  3, 10, 62,226,255,  0,252,255,253,255,254,255,255,252,254,255,  0,  6,  2,255,  0,  2,254,  0,  1,255,  4,254,  0,252,  6,  1,  0,  1,234,254,  2,230,228,252,  3,224,210,226,234,224,196,229,  4, 30,199,231,189,243,250, 70,253,219, 39, 75, 18,206,207,176,178,181,229,224,209,242,250,234,222, 29,254,240,254,247,247, 10,251,  7,  8, 13,247,192,191,190,192,192,173,182,191,225,214,191,213,204,191, 19,208,255,  0,  0,  0,254,  1,  1,  0,  0,255,254,255,  1,255,  6,251,227, 10,255,  6,255,246,  4,245,  2,  2, 10, 18,251,  1,244,254,
      215,219, 45,245,218,204, 32,158,212,227,250,242,202,176,172,180,164,241,194, 54,250, 24,177,251,191, 15,225,219,249,219,233,199,224,219,249,  5,  5,211,  3,247,220,242,  3,  6,  0,  2,  1,  7,230,217,  1, 39,178,244,243, 13,191,207,  3,  3,244,217,  2,255,  6,254,255,254,  1,254,253,254,252,255,254,254,  0,254,255,  0,255,255,  1,254,  6,  7,  4,253,239,  1,249,251,  4, 18,250,254, 12,248,249,239, 25,  1,245,244, 72,234,251,212,120,247,252,220,254,  5,252,246,  8,252,  2, 64,166,231,225,248,191,189,229,251,222,252,250,250,251,  4,239,252,234,248,235,247,  7, 43,246,247,254,254,  1,  0,181,255,255,255, 25,  5,208,  2,249, 30, 41,241,192,191,191,193,192,191,192,193,190,193,188,193,192,190,196,193,  0,226,250,249,251,251,248,252,195,205,191,205,197,192,206,214,
      156, 27,254, 77,207,196,208,190,212,216,193,205,234,  4, 10,255,236,228,185,234,240,  0,198,  5, 12, 13,248, 20,231,253,225,243,193,227,  7,239,233,251,217,  0,219,220,245,234,228,229,229,228,245,245,247, 23,201,239,220,237,223,241,221,240,196,210,237,234,190,190,196,191,191,183,189,186,193,197,193,194,195,194,196,195,201,195,194,194,196,197,193,194,198,194,199,197,192,194,191,193,185,  0,204, 12,126,195,251,  4,201,253,199,  0, 38, 31,  0,251,251,244,252,  7,164,193,  3,  2,254,245,250,243,193,223, 14, 11,  1,  0, 20,  0,254,255,  1,  0, 19,252,244,  3,255,254,208, 13,255,  0,  0,  3, 10,  4,251,228,220,173,211,172,237,246,200,247,255,  1,252,  0,254,  2,  0,255,  1,254,  0,254,  2,254,253,254,211, 23,216, 17,228,226,217,222,  1,252,250,250, 11, 35,255,243,
      250, 63,  3,228, 31,124,248, 65,  3,  5,248,250,112,103,180,248, 14,254,248,255, 32, 23,  3,254,239,233,234,250, 39, 53,240,246,189,130,243,201,209,183,253,238,224,224, 20,242,208,211,251,238,209,186, 11,253,223,212,224,247,227,220,  1,254,213,229,  8,  4,254,254,237,  2,  0,255,  0,252,  0,254,  0,255,255,  0, 24,253,  0,  2,255,  1,  1,  0,  2,  1,  0,  0,  0,255,  0,  0,254,  1,183, 52, 44, 42,119, 14, 67,244, 18,  0,255,  0,197,188,188,194, 16,  3, 10, 13, 60, 41,  5,224,249, 59,224,229, 45, 30,236,204,213,  3, 94,  3,  2,  1,123,251,250,251,  2,  4,  5,254,  1,  1,200,175,207,202,220,210,203,228, 45,  5, 16, 55,196,233,245,231,251,  2,250,255,255,  1,246,  0,187,  2,  0,  2, 22,255,  0,  0,  0,  1,  0,  5, 12,255,  0,  5, 25,  5, 16,199,207,229,254, 13,
    },
  },
  {  //   7
    {
       46, 13,181,230, 16,  1, 72,  0, 16, 10, 21, 19,197,  9,248,235,250, 10,  2, 35, 30,133, 59,130,  4,252,252, 33,170,  3, 66, 46,  8,  2, 52,252,  2,205,212,  3,253,  5,249,250,202,  1,252,253, 19,  5,  0,  4, 16,  6,  2,254,  7,224,  2,248,  0,163,253,252,255,254,  1,231,  1,  9,  2,  0,  0,  6,255,254,255,  0,255,255,  0,  1,255,  0,  2,255, 21,255,  5,  8,255,  3,  1,253,  1,  0,250,  0,237,  0,  5,  7, 12,255,254,  7, 76,  5,200,196, 17,154, 50, 60,  0,255,247, 32,255,237, 13,253,  1,255,245, 36,247, 17,  5,  3,252,  1,255,  2,  6,  0,  0,193, 73,200,249,232,191, 36,  8,  2, 64, 69, 70,199, 67, 65,  0,  1,242,120, 52,251,  0,243,  0,255,255,253,254,255,  0,  0,  0,241,  0, 23,  3,241,  2,248,249,  0, 23, 16, 11,250,  8,238,  0,255,  1,255,255,254,255,  0,
      219,  3, 12,249, 14,192,  4,182,248,253,201,  8, 63, 49,244,  3, 42, 13,252,226, 19,193,  0,190,  0,  0,  2,224, 54, 63,246, 12,  0,  5,202,234,244,  0,  5,  6, 13,192,  1,252, 58, 45,185,  6,252,235,187,251,212,235,  2, 24,  0,198, 63,  6,161,254,  2, 12,  0,252,255, 28,255,  0,  0,  0,255,235,255,  0,254,  2, 10,  0,  1, 25,233,253,  1,233,  0,234,  0, 12, 34,  3,251, 70,  1,  0,250,243,  0, 60, 19,252,  1,  2,213, 60,184,246, 11, 44, 58, 80, 10,235,255, 28,204,166,  3,  1, 13, 39,183,213,158, 66,  0,108, 62,  0,  9,243,248,  4,  1,254,213,  2, 17, 60, 53, 70,  9,218,237,226,103,230,226,253,  0,126,245,193,  7,247,241,201, 50, 60,253,  0,255,  0,255,  0,  2,  0,250,  0,210, 28,  0,240,  3,254,  0,  2,  3,255,  7,247,235,  3,  0,253,  6,  3, 27,251,247,253,
        1, 63, 74,226,  5,  8,245, 17,  2,  0, 11,252,  4,  1,255, 15,  2,  3,  0,  0,255,255,199,248,253,  0,247,255,254,255,  1,253,234, 27,247,235,  8,  2,  1,187,  0, 54,  9,252,250, 34,237,  5,226,247, 50,255, 48,253,  1, 11,  0,254,193,  0,167, 79,  6,241,  1,255,  0,253,  1,  0,  1,  5,  0,251,242,  0,  0,243,255,  0,243, 34,255,255,  1,255,  0,  1,238, 17,247,255,  2, 27,  1, 19, 66,219,220,250, 14,235,200,  0,186,255,241,211, 93, 14,  0,  5,194,251,125,180,253, 67, 70, 50, 20,250, 14,  9, 10,138,252,136,200, 26,242,  8,255, 22,  1,  3, 33,  5, 88,128, 15,255,244,255,246,201,  1,181,193,193,192,217, 62, 27,218,196,  2, 21,  7, 70,243,  1, 34,  2,  1,  3,231,  0,  0,  0,  1,244,  2,252, 19,  2,249,250,254, 15, 13,203,255,253,  1,254,  5,  0, 34,249, 11,  2,
       24,188,  0, 65,249,  0,234, 62,251,  0,  4,  2,  4,243, 64,240,  3,  4,251,  1,  3,  2,179,143,  0,  0,231,  8,  1,254,254, 12,  0,255, 20,  1,  0,  1,205,  3,  3,  4, 54,  7,  1, 41,  4,  6,230,128,  1,  1,221, 82,  0,237,252, 62,  0,  2, 16,233,251,250,255,  1,255,  0,  0,  0,  1,  2,  1,  0,  6,  1,  1,255,255,255,242,239,  1,  0,254,255,  2,244,253,  0,255,  1,  0,254,255,  1,218,247,207,  0,255, 26, 26, 60, 85, 25,225, 51,  3, 65,198,254,220,219,255,229,  2,223,245,248,251,240, 52,247,237, 68,201,236,233, 26, 14,197,212,252, 15,  8,242,  1,  5,167, 31,190, 24,222,245, 38,237, 64,  4,  9,  5, 98,216,220,  7, 12,  1,212,  0,255,  1,243,  9,  5, 29,  1,253,  3,  1,  2,  0,246, 19, 50,249,202,  2,  5,  1,  2,253,  1,  0,255,  4,255,  0,255,254, 59,  0,254,
    },
    {
       10,242, 74,190,  2,249, 65, 64,  4,246, 50, 47, 26, 16, 49, 46,239,254,  5,  7,217, 11,253,241,215,  1,233, 76, 37,  0,253, 76,  7,240,  9, 10, 26,224,247,246, 46, 53,246,249,  7, 19,247,249,224,  8,250,  1,224, 54,  0, 12, 51,204,249,251,230,214,253,  8,  1,  3,  0,  0,255,  6,  0,  0,  0,  5,  0,  0,253,  3,255,  0, 26, 13,255,  0, 15,244,  1,  2,246,  7,254,  4,  7,237,255,254, 45,241,253,  7,232, 27,181,229, 53, 46,  4,149, 30, 19,237,  0,204,242,  9, 17, 35, 13,247,  7, 13, 28,215, 62, 22,239,167,183, 18, 13, 25,230,239,220, 82,161,  3, 55, 11,  3,247,253, 64,239, 55, 13, 67,252,255, 20,180,  5,245,246, 25,217,252,  1, 92,224,  3,251,251,250,  0, 93,  3,  4,255,  2,246, 12,251,  7,254,  1,  1,  5, 20,  1, 24,  1, 70,  4,255,  2,252, 13,  1,255,  0,255,
        3, 71,254,239, 15,227,212,140,237, 29, 50,191, 11, 40,222,145,249, 66,236,  3, 33,224, 28, 55,220,252, 37, 57,248, 79, 16,  1,218,222,  1,255,248,225,  2,252,218, 52,254,252,241,254,  3,244, 41,250,255,253, 11,249,251,254,  4, 71,  2,  4,242,211,  5,  8,  0,  2,  0,  0,  1,  0,  0,  0,  4,253,  0,  0,  8,  1,  0,  0,250, 12,  0,  2, 20,  2,  1,  0,252, 32,  0,255, 27,133,  1,254,240, 23,234,227,239,235, 78,122, 20,209, 44,200, 21,  5,217,121,218,  2,253, 64,226,216,  1, 63,214, 39,249, 17,248,212, 59, 62,  8,245,253, 36, 15,  8, 46,225,  1,217, 19,  7,254, 12,184,111, 23, 69,  9,  2,212,  8,152,250,252,253,240,  7,243,213,218, 34,249,254,  1,  9,250,202,239,247, 10,246, 12,254,250,195,255,  3,  0,249,218, 11,248,194,250,  7,255,253,  2,  6,247,200,  4,250,
      246,253,  0, 37, 42, 17,234,108,218, 29,230,208,227,216, 10,178,254,205,242, 35,203, 66,235,228, 20,238,238, 40,208,127, 22,251, 18, 59,  2,248,232, 12,253,255,238,230,  6,  1, 12,251,  0,255,197, 42,  3,  1,  5,250,250,252,237,218,  8,  0,244, 55,247,  1,254,  0,  0,  0,  2,  1,255,  0,254,255,  0,  0,250,253,255,  0,250,240,  1,255,244,254,  5,  0,188,251,  0,255,218,221,255,  0,233, 13,205, 32, 97,190,249, 63, 38, 76, 84,  7, 58,203,237, 65,237, 28,254, 58,213, 17,219, 63,  4,  2, 14, 14, 49,191,234, 25,250, 35,210, 15, 31,234,150,221, 36,  7,  8, 33,242, 50,191,  9,  1,250,245,246,185, 67, 42, 10,  4,  1, 43,219, 13, 33,198,218,254,  0, 68,239,229, 80, 63,247,247,  9,  3,  1,181,244, 66, 79, 12,244,253,221, 87,212,252,  0,  1,  4,250,  4,255, 61,138, 53,
       21,230,204,  5, 11, 21, 11,  0,242,253, 66,241, 12,242, 91, 61,244,204, 11, 38,205,138,204,194,206,173,235,227, 94,  9,247, 12,  9, 56,  2,  0, 31,  7,254,  2,232, 12,  5,  4, 10, 25,  2, 17,230,  5,254,  0, 25, 66,244,  4,223, 39,  3,  2,  0,249,252,253,  1,  6,  0,  0,  2,  0,  0,  0,246,255,255,  1,252,  1,  0,  3,110,250,255,  0,255,109,  1,  1, 71,  6,  0,  0,227, 84,  0,  0, 38,201,231,224,190,216,242, 41,171,235,198, 13,187,253,255,199, 40,253,235, 23, 48,200,255, 26,211,250,  9, 10, 15,199,200, 62,  1,207,218, 24,184, 55,154, 59,254,221,250, 21,182,212,117, 15,124,229, 52, 56,  3, 86,244,236, 21, 75,221,  6,  5, 70, 65,217,255,232, 47,242,222,220, 70,215, 60,  1,253, 29, 54,220, 27,174,  1,251,236, 15,224,250, 13, 43,  5,  0,255,253,253, 54, 60,  8,
    },
    {
       41,254,213,235,196,217,203,226, 28, 11, 65,  2, 66, 66, 34,249,246,  1,178,195, 31, 45,247, 67,214,243,217,182,207, 36,249, 13,242,228,245,243, 52,223,  1,244,250,221,  0,254,253, 92,250,252,234,  1,253,253,202,235,  2, 18,219,187,  0,255,  9,243,  3,249,248,  1,253,255,  4,  4,  1,255,255,  7,  0,253,  2,251,255,254,252,  0,  0,  1, 20,234,255,  1, 12,  3,  0,  3, 59,  7,  0,  5,  0,254,199, 23, 16,233,186,111,239,186,254,226,205, 55,161,216,228,252,234, 13, 19,239, 28,121,222,255,245,243,155, 44, 45,231,  7,255, 23,  9,  5,  4,251,130, 23, 40,251,  7,251,218, 23,243,224, 49, 68,117, 57,255, 79,199,  0,246,246,242, 31, 69,187,  2,255,  1,253,254,254,  3,  7,254,  2,  0,  2,255,  3,249,254,  1,  3,231, 31,  1, 23,250, 18,  3,255,  0,251,  2,  3,  3,  2,251,
       14, 13, 28, 65,247,232,254, 26,  0, 16,175,230,251,255,197,250,  5,249, 57,101, 36,  0,231,  7,  2, 25,176, 42,196,239,148,243,  2, 19, 12,255,251, 58,124,  6,229,182,  1,  1, 26, 25,247,251,236, 41,254,244,  3, 23,  1,219,226,249, 36,249, 22,230,255,  2,244,  7,  0,255,255,254,254,  1, 11,253,  0,  3,250,253,255,  3,  9,253,  1,255, 73,252,  0,255,  3,  9,255,  0,235, 41,  0,254,254,222, 58,223,  4, 12,226,221, 85, 33, 58, 28,182, 27,240, 21,228, 16,232,224,186, 23,134, 35, 20,213,203,194,  3,231,193, 44,  6,  6,  3,  9,  3,249, 22,238, 18,252,222, 21,  5,  7, 11,246,  6, 62,207,213,197, 55,192, 15,252,245, 59, 48, 27, 52,209,241,  1,  0,  2,  1,  1,  2,  3,254,254,  2,251,  4, 16,  1,  7,255,254,  2, 41, 28, 12,126, 14,253,250,  0,  5,254,252,  1,254,  4,
      216,211, 31, 32, 81, 66,227, 34, 13,  4,246,  3,252,254,241,  6,240,243,  4,135, 50,233, 39,  5, 18,  7,240,247,  2,  1,116,  4,252,253,  2,  0,239, 19,  2,250,230, 26,198,255, 57, 29, 78,  1,205,251,  5,  4, 37,  8,252,  4, 50, 33,  2,254, 84,239,  0, 10,210,242,255,  1,252,246,  0,  0,255,255, 48,253,255,  5, 10,254,  4,247,  0, 17,  2, 13,  0,  0, 17,  1,  1,  0,  0,  1,  0,253,230,  2, 40,246,  0,248,170, 37,197, 31,  3,222,  1,  5,  6,196,210, 14, 25,151,229, 19,116,146,240, 44,245, 22,100,199,239, 68,247, 30,238,252,255,254,  1,  5,  6,245,  0,195, 15, 12,254,127,251,244,229, 35, 65,248,211,190, 11,188,241,239, 25,244,155, 57,252,  0,250,  1,255,253,  6,255,253,  7,  1,253,  1,  7,  6,  2,252,  5,  4,  3,252, 10,  1,253,249,  4,  4,255, 15, 16,255,  3,
      253, 58,249, 30,235,243, 82, 10,240, 51,  5,  1, 25, 26, 75,253, 14,248,  2, 35, 10, 26,212,  3,  8,236,103,192,253,  5,227,189,  8,254,  2,  2, 41,217,  2, 10,250,247, 63,  0,241,184,221,252, 52, 81,  4,  2,195, 17,250,  0,237,  2, 33,  0,  4, 20,249,  2,244,  0,255,255,253,  0,253,255,  9,254,  1,254,  1,254,253,  0, 20, 14,  2,  3,136,  4,255,  0,  0,  3,  6,  1,249,248,244,  2,  9,238,253,  4,  0, 18, 92, 41,225, 20,198, 30,229,190,  0,159,222, 43,155,185, 45,234,240,120, 24,217,  9,191, 57,172,230, 33,  9, 38,192, 40,255, 14, 42,  7,  0, 59,  6,228, 19, 68,188, 12, 30,232,252,250,249, 54, 59,140,234,  1,253,178,249, 24,171,114,  3,242, 12,  2,  1,  8, 13,  0,  0,  0,  3,  3,251, 19, 22,192,255,253,  0, 33,234,  5,253,  2,252,238,  2,  1,  0,128, 11,  2,
    },
  },
  {  //   8
    {
       11, 25, 23, 40, 24,249,  6, 10, 22, 18, 29, 11,252, 66, 16,198,234, 30, 27,207, 30,  3, 72,  5, 44,230,199, 21,221, 90,242,207, 15,253,251,226,255, 27,251,243,  1,  3,247,251, 59, 15,  0,221, 35,251,255,249,251, 35,237,  1, 62,239,244,  6,171, 14,255,254,  3,  1,  0,252,246,  3,  1,  0,  1,  2,  0,255, 48,  1,  1,253,  5,  2,  1,255,  2,244,  1,  2,  1,249,  2,  0, 71,252,255,  0,247, 18,234,247,252,  6,254,  0,  0, 20,242,211, 20,195, 46,252,242,244, 19,246,243, 62,236,250, 12,  2,235,  6, 51, 71,190, 28,  4,  1, 19,253, 14,249,  2,  0,  0,213,244,253,253,  0,245,  0,255, 13,  5, 54, 59,255,136, 68, 70,  3, 12,249,246, 68, 62, 43,  2,  0,  6,  1,255,  2,  1,255,  0,254,  0,  0,  0,253,254,  5,  1,  0,  8, 12, 19,  8,249,  9,243,254,253,  2,  3,236,  7,254,
      231,202,245, 17,  5, 80,254, 56,203,233, 27,134,  5,237, 25, 20, 52,  5,205, 42,125,221, 25, 59,  4,226,254,192,  3,174,  2,255,214,234,254, 40, 43, 57,  0, 24, 18,  2,  3,255,  6,225,  4, 18,236,  2,249,246,154,234,252,248,220, 24,  5,  6,237,  1,  1,  4, 32,  2,255,  0, 21,255,  0,  0,  1,  0,  0,  0,162,250,  0,  3,  4,  1,  0,255, 80,  7,254,255, 19,  6,  0,  0,  2,  8,252,  0,239, 27, 27,205, 10,240, 12,249,184, 53,233,165, 54,237, 48,198,  0, 29,251, 42,233,122,248,241,  5,  7,229, 70,  3, 26, 51, 45,249,253,255,246,252, 11,  9,255,248, 15,255,254,221,  1,246,  7, 11, 30, 55,253,  8,253,144, 92,  0, 66,246,255, 13, 55,  0, 56,  7,255,253,  4,  1,  1,  1,  0,  0,255,  0,  0,253,  6,253,  1,249,  7, 33,244, 14,254,255,  5,  3,  4,254,  3,  0,  0,246,254,
       61, 69,219,169,220,240,226, 28, 43,255,252, 67,181, 79,  5, 32,238, 71,217, 54, 78, 85,236, 20,214,254, 12,144,  6,127, 29,187,197, 96,  8,  1, 28,206,  0,252, 11,210,  9, 27, 23,  5,  1,254, 13, 10,250,  5,182,  0,254,255, 62,238,253,243,232,244, 14,196,255, 13,  0,  1,254,  2,  0,255,216,253,  1,  0,255,255,  0,  0,254, 17,  0,  1,236, 44,  1,  0,254,  5,251,254,  0,  1,129,  2, 10, 17,240,242,250,255,245, 23,  8,254,102, 49, 10, 64,243,186,253,255,  1,250, 38, 38, 83, 21,255,207,255, 60, 61,254,  9,  5,  3,252,  3,251,255,  1, 17,  0, 28,255, 10,  4, 25, 11,  2,238, 50,252,228,244,178, 23, 71, 13,255, 39, 12,252,235,  8, 15,193,237,  0,  1,254,  0,  1,254,  0,254,  1,  1,  0,255,  4,255,  6, 12, 31,252,254,  0,252,  0, 18,251,  0,  0,  3,242,215,205,254,
      252,254, 44, 39, 58,  0, 40,175, 33,189,  7,  3,240,  7,254, 61,207,207, 77,246, 88, 70, 27,199,235,  1, 36,190,  0,  5, 16,249, 52,211,  0,247,208, 15,  4,  6,209, 14,  1,  5, 27,245,246,254,232, 50,239,  8, 90,220, 12,  6,187,254,252,  0,248,194,  2,  3,244,253,  0,255,254,  0,  0,255, 32,  4,  0,255,  3,  0,255,255,240,253,255,  0,223, 17,  2,253,  5,  0,  0,255,  4,  7,127,  3, 24,  1, 11,241, 23, 41, 45,  1, 16,207,211,197,166,202, 35,  8, 10,252, 30,255,193, 81, 30, 68,240, 16,  1,242,183,248,250, 31,251, 13,245,246,237,253,237,253, 12,  0,252,253, 62, 42,  8,  5,233,191,201,243,233, 54,247,198,184,  7,253, 46,  1,196, 79,254,255,  0,231, 54,252,255,252,  2,  0,  2,255,  1,248,212, 57,253,  1,251,249, 24,246, 41,  0,  3,254,254,  0,  0,  0,  1,254,239,
    },
    {
      208,  2,201, 11,210,253,207,  9,219, 21, 13, 18,160, 17,231, 37, 37, 63, 34, 11, 46, 35, 17, 57,  6,  3, 15,206,248,  6, 38,144,238,195,  5,255,191,191, 21,  4,245,255,253,245,  4,211,  4,  7,202, 82,244,  0,206, 43,  1,  4,217, 28,  0,232, 21, 35,  1,  0,255,235,255,  1,  0,246,254,  1, 15,  0,  0,  1,  3,251,  0,  0,253,254,  0,  0,  2,255,  0,  5, 29,  4,  1,  1,249,245,253,  1,223,218,  8,221,235,216,233,227, 55,251, 53,248,  6,  9,254,101,  1,  7,  1, 19,177,205,  6,194, 29,209, 75, 25,210,245, 48,  3,  6,155,187,242,251,205, 35,173,  8,246,215,  4,212,223,181,  2,255,  4,  2,251,201,233,250,136,  2,222, 26, 13,  8,  5, 19, 62,233, 18, 51,  4,243,243,240,121,  8,  3,  4,  1,242, 17,  1,  3,  0,  1,255,254,246,  3,250,  1,255,254,  1,  4,237,  3,249,136,
       24,  9,224,209,243,255, 38,139,199, 16, 11, 24,249,212,  6, 44,225,249,237,144,221,238,231,201,232, 43,248,236,223, 76,240,234,223,189,250, 21, 26,  3,  4,  5, 36,185,  0,230, 23,  6,244,254,245,237,247,  2,  6,  5,  0,  1,194,237,254,239,174, 35,253,  1,255,  0,  0,  0, 14,253,  0,  1,  3,255,  0,  2,  1,  3,254,255,  1,  0,  1,255,240,255,  0,255,216,250,251,  3,190,  1,  1,254, 73,178,183,  6,253, 52,232, 53,  6, 29,249, 31,  9, 11, 12, 32,222, 45,244,253,251,244, 61,247, 10, 98, 36, 33,230, 52,142,156, 36,224,  1,254,  7, 88, 15, 53,231, 43,210,  1,249, 37,234,237,183,  7,255,252,254,190,205, 74, 12,  4, 13, 13,163, 20,255, 28,250, 23,241,219,250,251,250,  3,242, 23,238,254,  3,245,224,  3,253,  0,  0,250, 11,251,201,255,  0,  5,251,254, 43, 13,184,194,
       75,238,248,244,  5,  3,214, 25,219,215,255,253,254,  1,102,171,186, 24, 45,168,  1, 35,191, 65,  0,116,196,218,  1,  3,239,225,252,190,  9, 31,195,190,  4,252, 19,  8,  1,247, 11, 37,254,  3, 23, 16,  5,248,235,238,254,217, 64,133,253,252, 70,210, 42,248,254,  0,  0,  0,235,253,  1,254,  6,  1,  0,  0,254,  0,  0,  0,255,  8,  0,  0,252,255,255,252,  2,246,255,  2,189,254,242,255,187,250,235,  9,133, 10,187,215,232, 32,249,236,202,218,  5, 36, 18,208,222,  9, 52, 68,  8,191,  3,226,194,129,254,  3,255,222,244,230,210, 17,160,232, 60,251, 37,  4, 31,253,215, 19,133, 79, 62,  2,253,  1,168, 19,208,  4,  0,253,  5,  0,248,223, 65,186,254,  9,188, 26,196,  4, 11, 59,254,  0,237,248,  3,251, 64,222, 13,  2,253,254, 10,245,252,203,255,254,245,250, 25,  4,114,  0,
       60,  3, 62, 62,  9,238,220,253,253,  5,221,217,227, 44,250,255,255,229,  8,180,  7,250, 63,179,191, 45, 46,  6,  3,233,254, 43,219,197, 39,250, 32,  3,199,  2,237,  3,  7,248,235,248,  0,  9,255, 34, 43,  0,  4,  2,  1,  2,  2,255,251,251,255, 19,254,253,192,  0,195,  0, 64,  1,205,254,203,  0,193,  1,161,251,192,254, 63,254,189,  0,188,254,193,  2,189,  1,192,  0,  4,163,190,  1,180,248,231,254,235,  1,252,  5, 11,  9,  1, 40,205,  8,248, 28, 52,231,189, 19,253, 24, 37, 11,133,  2,100,209,255,198,233, 68,234,232,246,  9,245, 53, 22,190, 57,249,255, 58,192,225,  9,223,211, 33,  9,  8, 26,194,249,238,155,250,254,208,223, 34,  2,197, 31,189,253,167,193, 30,234,225,200,255,249,  2,190, 10,229,230,245,255, 61,231,192, 58,200, 44,204,255,  3, 36,  0,  5,189, 17,
    },
    {
       14,  4,249,255, 40, 35, 77, 64,255, 17,  7,251,208,  2,217,123,226, 24,252,100,247,224, 35,217, 10, 40,  3,247, 50, 96,224,  5, 11,252,251,236,  5,  4,  0,249, 17,  6,252,251,  3,252,  5,  3, 24, 11,  1,  0,249,239,  5,  3,253, 33,  7,  2,  3, 58,  4,253,  5,  8,  0,255,  4, 10,  0,  1,254,  8,255,  0,  4,251,  0,  1,  0,250,  0,  0, 51,250,  1,  0, 75,242,  3,  2,190, 36,251,  1, 17,213, 25,175, 45,213, 17,237, 18, 24, 31, 20,235,253,255,126,250,241, 18,251, 11, 12, 82,250, 13, 10,220, 14,183, 33,183,243,234,216,235,230, 79,222,243, 45,246,255,  3,247,226,249, 20, 59, 18,246,  0, 61, 94,228, 14,185, 59, 31, 90,236, 12,251,138,251,  3, 34, 48,255, 13,242,253,  4,  0,255,  0,  1,  4,255,  5,252,  1,255,  3, 19,252,  5,254,  4,  1,251,251, 10,  4,254,  3,  3,
        5, 11, 91,132,  7, 13,  4,205,208,  5,227, 51,249,249,244, 29, 37,  2, 55,237,  3, 46,  9,245,233,219, 59, 63,232,184,  3,244,255,231,  4,  4,229, 10,  1,  2,242,253,  4,252, 49, 52,255,251,228, 56,  2,255,202,209,252,250, 46,  0,254,253,229,  3,249,  3,254,240,  0,  1,252,244,  0,  0,254,251,  1,  0,  6, 12,255,255,  0, 45,  0,  0, 69, 35,255,254, 28,220,  3,255, 35,252,  3,255,251, 13,  1,251,244,  0,253,122, 55, 67,231,208,204, 66, 36,182,252,253,  6,231,255,160,253, 73,194,  9,  2,216,253,194,220,248,  1, 11, 37, 25,  9,224,116,  7, 28,251,233, 36, 94, 11,182,101,  0, 41,252, 12, 20, 41, 82,208,220,216,225, 49,238, 91,142,235,249,  5,229, 12,227,  7,253,237,255,253,255,  4,  1,  2,  8,  7,  0,  6,249,250,  1,248,  1,246,255,  1,  1,  6,  1,  4,253,  3,
      255,246, 18,  6,250,  6,219,228,210,226,189, 96, 22, 16, 56,131,255, 73, 14,  1,255,236, 28, 14,220,224, 14,210, 30, 32,209,243,255, 20,253,  2, 25,231,  4,  0, 44, 35,  1,253, 18,246,  0,254,232,  5,  5,254, 23,196,252,  0,232,206,  3,254, 26, 58,  0,  0,  2, 14,255,  2,  2,247,  1,  0, 15,  4,255,252, 14,254,  0,252,  7,  9,255,255,191,252,  9,  0,196,244,  4,  1,248,244,254,255, 47,244, 40,248, 30, 14,228,193,253,245, 12,237,  1, 37, 97,228, 17, 87,244,  7,242,242,216,155,236,194,154, 21,235,223,129,198,252, 29,252,227,  1, 54, 91,255,  9, 57, 31,203,248, 71,132, 68, 11,238,239,247,225, 37,105,  2,240,  6,251,192, 11, 53,235, 45,253, 55, 34,  4,239,  1, 13,244,  0,  2,255,253,250,248,  2,253,255,  2,231,251,  5,236,249,  4,254,255,247,  4,255,  2,246, 10,
        0,  6,192, 30, 19,  2, 24,136,234, 60,  8,182, 24, 21,230,168, 17,242, 67, 20, 12, 97, 15,239, 15,223,232,254, 21,237,198, 51,254, 60,255,255,250,245,251,243, 14,252,  1,251,  9,235,  0,255, 34, 13,245,  3, 80, 72,  2,  2,234,  3, 11,  0,201,195,  7,250,251,  9,  0,  0,252,  4,  0,  4,240,  2,253,  0,246, 11,  0,  0,250,251,  2,255, 43,239,254,253,250,  1,  2,255, 71, 31,255,253, 26,  2,248, 20,237, 68,228, 39,227, 65, 39, 97, 48,213,195,183,  6, 22,206,251,224,162,206,241,250,218, 48,192, 12,131,  2, 70,  0,250, 63, 37, 64,195,105,241,  0, 13,229, 55,220, 83,246, 24, 52, 37, 10,228, 56,170, 66, 48,120,  6,  0,  6,  9,  8,245,213,  1,255,251,248,255,243, 67,247,  0,  2,  6,254, 19, 21, 11,  1, 26,  0, 31, 22, 52,212,  2,250,254,  6,255, 50,253,238,  1,  6,
    },
  },
  {  //   9
    {
      216,196,246, 19,212,243,255,  0,214,223,241,238,239, 17,  6,255, 36, 23, 54,236, 22, 21,  3, 75,234,253,249, 55,244, 16,252,219,  5,  0,251,  5, 21,  8,251,240, 16,  9,255,253,248,  2,252,246,  9,244,  6,  0, 71,238,  1,252,203,226,  1,248,210,233,  1,255,  2,  0,  2,255,  2,  0,255,  0,255,  1,  0,  0,  3,  3,  1,  0,  1,  6,  1,255,250,  3,255,  0, 13,249,  1,  0, 22, 12,253,  0, 13,255,194,238,228, 12,220,255,211, 45,221,  4,248,251,225,240, 11, 59, 33, 31,238,230,200, 35,237,244, 39,222,196, 12,239,226, 24,252,250,243,246,217, 27,  9,  5,234,232,252,  9, 51,232, 73,224, 47,  0, 53, 40, 28,227, 14,245, 27,250, 13,232,249,232, 51,252,  0,249,252,  4,  3,253,238,  5,254, 14,  3,  7,  7,  8, 22,197,210,  9, 20, 49, 10,132,244,250,248,  9,248,204, 60,204,213,
       39, 35, 18,235, 41, 49, 10,229,227, 48, 17, 86,  9,224,  8, 19,248,226,  0, 11,235,214,188, 52,  0, 35,251, 23,  1, 59,  4,  4,241,  7, 51,247,248,  1,253,  6,240, 35,  0,249, 38,253,255,  3, 24,245,248,  2,238, 44,246,  6,  6, 16,  1,  0, 26,  2,  6,  3,253,  0,  6,  0,255,  1,  0,  1,  1,  0,  1,  0,  3,  0,  0,  0,  1, 12,255,254,255,  1,  0,254,255,  1,  0,255,  5,  0,255,  0,212,236, 11, 18,237,243,  0,231,210,213,243, 18, 79,241,214, 51,244,250,230, 40, 36,228, 23,  7,  2,227, 38,198, 26,  3,129, 65,  0,240,254,250,251,  7,  0, 90,  0,244, 15,  8,255,191, 88,143,228,194,254, 60,203, 14, 69,205, 25,232,235,240,  3, 63,251,120, 12,  2,  9, 46, 10,255,  0, 85, 25, 15,254, 13, 27, 59,243,231, 58,  2,253,184, 52, 36,229,227, 11, 17,  4,  4,229, 36, 24,241,
        7, 59,244,  2,  6,237,  0, 46, 43,  2, 32,143,  5,208, 30, 17,249,252,187, 44,  5,  5,  7,142,  3, 49,  6, 42,240,246, 48,131,  3,247,215,237,  3,  5,246,  4, 25,237,242,  1,  1, 19,196,  2, 65, 61, 13,  0,251,242,252,  5, 35,244,254,  0,236, 24,244,  0,  5,255,238,255,  0,  1,  0,  1,255,  0,  0,  1,  1,254,  1,  0,253, 18,  0,  0,  0,  0,  1,  1, 46,  9,  0,  1, 27,253,  1,255,220,251, 19,202, 19,225,252, 49, 31,241,244, 13,194, 13,243, 44,240,181,241, 42,241,  1,223,239,248,250,228,232,220,238,136, 58,254,240, 78, 23,  0,251,252, 15,252,236,253,  9,  7,  5,249, 49,251, 15, 10, 13,253,255,123,133,  0,247, 35, 19,228,253, 80,247,  1,247, 48, 12,  2,212,  2,255,245,239,255,  8,240,253, 16,236,248,  2,214,  4,  1,184,248, 10,  3,248,  0,246, 23, 39, 12,100,
      249,  1, 10,  0,  1, 48, 14,212,246,232,  7, 65, 40,249, 18,251,193, 23, 63, 21,249,239,250, 34,251,250,241,242, 60,240,244, 26, 11,254,251,254,242,236,  0,  3, 19, 23,196,245,  2, 36, 87,  2, 33, 19,255,  4,  9, 35,  2, 20,237,252,255,  2, 42, 71,253,253,252,  4,  3,  1,  1,255,255,  0,  2,255,255,254,  0,  3,236,254,  1,255,  0,255,  2,  0,  1,255,246, 17,  0,  2,252,253,254,  2, 55,243, 12, 49,242,  0,255,199,208,  7,253,246,  2, 25,191,243,255,242,  9,213,  0,203, 44,188,251, 13,245, 93, 23, 16, 49,204,255,227,210,254,  4,  2,174,237,253,223, 13,255,  8,  6,229, 14,  2,  0, 37,  1,197,193, 61,188,246,255, 11,188, 57,  7,224,255,  4, 23,236,237,  3,254,  2, 20,235, 51,  3,236,218, 70, 18,  3,  3,255,252,246,  8, 33,255, 37,240, 10,  7,228,209,238,  3, 11,
    },
    {
       25, 20, 66, 58, 21,248,151,175,243, 21, 51, 24,237,216,  9,  4, 25, 21, 31, 39,233,177,237,218,186, 71,201,191,255, 42,  4, 62,241,  3,230,248,  5, 10,250,253,236,  3,249,253,  5, 13,252,255,216,  3,254,253,221, 34,  0,127, 48, 21,254,255,  4,238,254,  2,253,  3,  0,255,  2,253,  0,255,254,  1,255,255,  0,  2,  1,  0,  5,  7,255,  1,  1,  0,255,255,  4,  2,  2,254, 54,  9,  3,  2, 17,235,238,  2, 24,108,229,192, 54,225, 34,238, 47, 15, 25, 42,  7, 16, 13,255, 15, 11, 92,106, 44, 32,241,240,209,  4,246, 79,251,192, 19,  3, 61,167,147, 20,234,  0,245, 62,196, 39,179,231,250,  6, 61, 35,249,  3,201,128,243,185,  1,225, 43,242,118,246,243, 32,  1, 12, 17,  0,234, 38,254,255,252,253, 21,252,  8,248,254,  0, 35, 25,  2,255, 18,254,  4,  5,252,215,245,  3,  8,117,
        5,255,211,216,  8, 50,154,166, 31,253,207,127, 45,221,243,235,249,159,222,141,241,220,247,167, 74,246,254, 28, 23,236,188,231, 48,213,  1,  2,230, 20,  4,  4, 75,222,  1,  0, 11,225,  3,  1,253,222,  0,254, 73, 19,  0,127,  3,  1,251,  2, 16,  6,251,  2,  2,237,  0,  0,  0,  0,  1,  1,  5,  2,  0,  0,  5,  0,255,  0,  4,251,  1,  1,  7,244,  0,  0,140,230,  0,  0,  3,  4,  3,255,  5,248,192,246, 45,254,169, 29,252,184,  7, 35, 23,128, 66, 66,231, 24,  1, 40,187,228,251,194,231,  1, 27, 78,236,170, 37, 44,243, 68,224,190,225,202, 15, 17, 65,195,250,194,190,155,142, 87,  2,255,198, 39,206, 13,250, 17,247, 30,194,224,187,234,118,243, 19,  8, 17,  5,245,255,242, 20,248, 15,  6,252,  0, 58, 28,  8,  2,  3,253, 22,229, 27, 43,253,  5,  1, 13,247,251,251,249,  3,
       24, 14,220, 67,  1,238, 27, 16,248, 15, 20,242,233, 73,  0,163,235, 29, 29, 11, 10,253,224,201,207,  0,198,184, 37,228,242,139,165, 41,  2,254,243,  3,254,  1,225,230,  0,  2,253,254,255,  0,253,  8,255,249,228, 72,  3, 26, 13, 29,  0,  4, 71, 35,255,255,254,  2,  0,  0,  5,254,  0,  0,252,253,  1,  1,  0,  1,255,  1,254,  9,  0,  0,252,247,  0,255,194,  7,  1,  1,240,244,  1,  0,251,205,218,217, 24, 24,253, 96,245,246,250,194, 60, 63,204, 43, 25,247,204, 42, 11, 40, 35,194,242, 28,236,243, 67, 60, 19,188,  6,237,  0,121,166, 12, 84,233,247, 45,229,247,135,250, 32, 93,  9,243,248,  5,238,  4,198,  6,214, 29,  3, 17,250,216,  2,187, 38,251,243,247,  1,253,  4,248,245,  0,255, 31,250,  1,210,231,  2,  2,  5,231,236,  1, 24, 66,  7,253,245, 17,  0,239,223, 24,
      242,244, 31,171,234,  1,185,220, 16,  0,231,249,207,240, 96, 85,  3,237,251,191,246,211,178,218, 67, 18, 18, 70,241,220,241, 23, 71, 41,  3,  2,  6,237,  3,254, 46,236,255,  1, 27,247,  1,  0,  0,  9,245,  4,255,  7,  4,  0, 11,250,  0,  0,196,204,  3,250,254,255,  1,  1,  5,  2,  0,  0,253,254,  0,255,255,  2,  1,  0,  2,248,  0,  3,  4,  3,  2,  0,186, 35,255,  0,235,  2,  1,  1, 11,201,205,192,199,143,128, 19, 52,250,  5,226,164, 70, 35,248,  6,241, 30,190,223,227,  7,254, 52,194, 47,242, 53,123,145,191,  6,  4,  5,253,109,208, 21,  2, 20, 33, 23,254, 34, 54,140,214,  4,232,249,244, 48, 84,252,187,255,  2,135,157, 71,181, 26, 55,  9,230,244,231,191,119,250, 18, 13,  2,237,150,247,247,249, 28,255, 20,249, 30, 48,198,248,200,242,255,  2,251,238,  7,241, 25,
    },
    {
       14,187, 50,200,230,200,243,244,242,185,240,  0,209,203,  0,251,224,219,252,248,239,244, 20, 22,234,226,237, 17, 33,  3,251,251,249,254,177, 10,  8,246,246,254,254,235,  2,254,252,228,254,253,219,238,250,254,204,230,255,  2, 34,213,255,  0, 31,221,  0,251,  1,  3,251,255,  2,254,  3,  0,  3,  0,254,  0,255,255,  0,  0,  8,  0,  0,255,  8,  3,255,  0,243,  8,  2,255,220,  7,  5,255,  4,246,247, 66,246,  4,  5,249,190,177, 31,250, 28,211,185, 54,241,219,244,253,205, 17, 60, 79,229,226,247,248,244, 41, 75,248,255,253,  3,  1, 14,250,250,  2,240, 15,  3,252,243,255,253,  0,230, 22, 29, 64, 86,  2, 83,227,  9, 45,  8,241,  3, 44,  7,187,  0,  4,255,251,  5,  0,  1,247,254,  0,242,253,216,254, 24, 25,244,237,  8, 17, 29,255,255,178,  2,  2,252,  1,  4,247,  0,203,
       45,113,139,229,250,  8,255,252, 11, 14,  7,  8, 32,249,  4,255,201,242, 58, 57,  6, 49,230,251, 25,237,  1, 87,187,  3, 42,  9,  5,  6, 20,153,253, 23, 10,  0,253, 42,255,  1,232, 37,  7,254,  9,247,254,  0, 53,  4,246,255, 34, 45,  0, 22,250,  1,255,  3,  3,  1,  0,  0,254,  2,  0,  0,254,  1,  0,  0,  5,  0,254,255,251,  5,  1,  1,244, 10,  0,  0,  2,  0,  1,  1,247, 10,  2,  0,  2,251,192,211, 13,218, 11,  1, 25,193,238,227,191,212, 41, 69, 10,203,234,250,  7,217,253, 13, 24,252,  3,235, 92,  7,139,211,  2,  1, 22,  1,250,  2,221,  2,240,252,255,240,223,243, 36,  3, 37,252, 36,243,191,182, 23,  2, 30,251, 30, 21, 88,191,252,245,  1,  0,251,  0,  0,  3,  1,254,  1,253,  0,253,241,209,242,224, 35,243, 60,225,200,200,  4,212,  0,  1,255,  8,  4,239, 23,  4,
      158, 56,215,179, 40,203, 29,  7, 43, 65,232,  0, 16, 18,  9, 11, 43,186,250,214,236,244,242,  4,231,244,253,238, 37,213,214, 38, 10,  3,251, 27,219,252,201,  4,250,204,253,  3, 24,225, 64,  3, 69,  0,  0, 59,207, 15, 58,  2,194,221,254,  6,205,209,255,  4,  2,  0,  0,255,255,255,254,  0,  0,  0,  0,  0,  0,  4, 21,  0,249,  0,  1,  1,254,  2, 30,  0,237,  1,  0,254,  3,  8,  2,  1,  0,254,247, 59,253, 19,254,253, 20, 22,236, 18,  4, 43,  0,145, 23,246, 37,235,  3,191,  1, 70,253,253,  1, 13,200,184,159,254, 10,  1,194,  1,249,  3, 30,255, 11,217,251, 17,206, 28, 26,177,254, 14, 14,  5,  8, 15,165, 51, 14,188, 19,162,254,250,252, 11,255,  0,255,255, 14, 19,  1,253,  0,254,255,254,203,  2, 16,236,251, 25,252, 64, 73, 91,  0,143,251, 26,247,231,255,202,255,228,
       49,189,  9,161,  9,255,246,  0, 19,  0,  3,  4,  6, 18,150,  1, 31, 21,209,236, 37, 25, 79,254, 59, 39,205, 10,255, 69,253, 16,  9,  0, 11, 26, 56,  3, 65,  0,253, 21,254,  2,  2,255,190,250,171, 58,250,255,  0, 30,192,208, 32, 43,255,253, 59,  0,  2,  0,255,254,252,  0,  3,  2,  0,  1,255,  1,255,  0,  1,  0,  2,  1, 39,252,  0,255, 23,  0,  0,255,  1,254,  0,  0, 14,  2,  1,  3,  1,249, 22,  0,219,244,247,254, 33, 23,253,222,145, 17,142, 86,232,  3,248,255,  1,  2,234,222,255, 30,244,  6,251,112, 57, 12,255,250, 66,  4,  0,  3,199,  2,  3, 53,219, 14, 54, 43, 35,221,235,206,209, 44, 19,247,236,244,222, 65,215,235,  7, 72,  3,116,  0,  3, 43,238,201,245,  4,253,  2,224,  0,  0, 18, 27,253,218,254,245,  0, 35,238,242,  4, 16,  3,  0,  0, 37, 60, 12,  6,  5,
    },
  },
  {  //  10
    {
      240,210, 14,  9,222, 33,  6,  5,  1, 18,248,253,253,  9,  9,  7,248, 48,  8,250,253, 10, 10,  4,  9, 15, 35,255,246,236,  1,  8,240,229,243,204, 27,234,248,251,244,  6,  4,254, 15,  9,253,251,237,232,246,254,253,241,247,253, 21,212,  5,254,208,190,255,  0,  3,  1,  0,  3,  2,  2,255,  0,  1,  0,255,  0,  4,  5,  1,  0,  4,  3,  1,  0,  6,  5,  2,255,252, 17,  1,  0, 21,254,  0,255,  1, 27,247, 33,  6, 31, 33,  2, 19, 43, 39,196, 13,179, 77,230,221,  4,  2,232, 34,203,160,246,229, 83, 69, 36, 21,217,210,151,  2, 12, 14,250,  4,243, 12,255, 21,249,239,254,231,245, 35,228, 17,218, 63, 37, 74,117, 71, 43,200,195,250,197, 13, 61,252,234,255,  2,  0,  0,255,  1,  1,  0,253,  0,253,  0,247,253,  3, 12,  2,  6, 14, 24,248,  3, 13, 18, 11,252,253,  0,  0,191,254,254,
       13, 43,240,249, 26,141,213,250,245, 61, 12,  9,247, 77,246,  3, 20, 16,198,  8,182,202, 13,254,200, 40,  1,254, 57,235,130,255, 98, 85, 33,218, 37,  8,248,  0,252,233,  0,  4, 58,238,253,  4,254,  7,  3,  3,253,249,250,254,212,232,243,252, 34,238,249,255,255,  1,255,  2,  2,254,  0,  0,  2,  0,255,  0,227,  0,  1,  0,252,253,  2,  0,198,  0,254,  0,  7,  0,  0,  0,  1,  0,  0,  0, 27,226,252,241,211,178,191,255,  0,235,  3,195, 11, 60,221,179, 33,248, 50, 41,  2, 19,241, 26,218,253,251,242,  1,197,253, 38,220, 15,253,244,250,213,236,  1,253, 15,245,241,251,  5,253, 97,238, 65, 17, 13, 11,135,250,223,206,217, 67,199,  4, 47, 75,  4,  2,249,  1,  1,  0,  1,255,  0,255,249,  2,  0,  1,  1,255,  8,223, 32,244, 66, 47,217,229,  1,243, 15,  2,  3,189, 10,222, 24,
       20, 45,222,  1, 61, 85, 60,193, 31,245,  8,252,  7, 15, 11,  3, 29, 79,194,  3,  3,235,251, 45,231, 10,  8,194,203, 14, 12,169,187,211, 16, 60, 44,217,  4, 37,  9, 53, 56,  2,  7, 47, 10,255, 48,  3,246,255, 41, 66,  1,  8, 18, 41,255,  1, 16,200,  1,255,  0,  3,  1,  0,  0,  2,  0,  2,  0,  1, 30,  0,253,253,  0,  0,  5,253,  1,  0,178,255,255,  0,  3,  1,255,  0,  4, 15,  0,255,250, 28, 27, 46,  0,191, 55, 80,247,  0,255, 26, 11,246,255, 60,250,197,239,198,251,219,115,205,249,  1, 16, 38, 12, 67, 17,100,  2,  1,  0,252,252,  4,255,254, 15, 61,252,241,  3, 69,249,  8, 11,218, 10,231,203, 52,174,200,243,238, 48,  2,235, 66,252, 53,255,  3,254,245,  3,246,  1,  3,  1,251,255,  1,  3, 39,  0,  0, 30,219, 40, 17,246,237,  0,234,  2,255,  4, 29, 63,  1,127,255,
      227,233, 59, 69,  1,223, 26, 12, 13,237,251,  0,184,243,251, 14,  6,197, 31, 75,218,246,245,193,235,224, 51,169,  6, 63,198,112,  7,  4, 14,254,182,245,248, 49,253,238,200,254, 21,221,254,253,193,  1,249,  1,  9,  2,  0,  0,  2,238,  5,  1, 12, 95,  4,  0,255,  0,  1,  0,254,  0,  1,250,253,  0,255,  0,254,  2,  0,  0,  8,  1,255,  1,129,  1,254,  0,244,  1,  0,  0,255,  1,255,  0,250,243, 16,199, 14,229,230,193,227,202,215,218, 10, 71, 25, 63,255,250, 17,240, 32, 68,136, 37,187,250, 31,248,180, 25,  1,150,250,  2, 79, 10, 11,211,  9,  7,236,222,250,228, 59,237,232,133,240,236,235,183,  5,230,  1,224,236,230, 17,210,250,241,212,124,253,  2, 31,243, 13,  1,252,  0,  1,245,  1,  5, 21,216,252,  2,255,  6, 13,194,253, 10, 10,  7,  0,  3, 18,248, 36, 96,128, 20,
    },
    {
      238, 36,249, 55, 11, 14, 22, 29,239,  0, 48, 22,201, 45,  1, 94,249, 18,255,208,  3,214,226, 32,249,250,252,213,249,234,  5,  0,  4, 24,220, 11,229,209,250,255,184,245,254,  4,196, 11,252,  2, 19,  6,  0,239,232,243,  4,  3,226,239,255,252,220,233,253,251,254,  1,  6,  0,251,  5,  0,  0,254,  1,  0,  1,250,  0,  0,255,  4,  1,  0,254,  4,  4,255,  0, 11,  6,  1,  0, 23,  3,  0,  0, 39,243,  6, 72, 43,240,  3,237, 71, 12,232, 46,232,217,230, 46, 22, 10,223, 52,239, 53, 58,252, 45, 12, 51, 34, 35, 81,198, 97,  4,223, 13,254,247,  0,  1,  0,235,246, 23, 20,230,229,  0,  9,  0,250, 38,204,187,  6,231,  4, 32, 21, 23,246,252, 60, 71, 41,252, 15,  7,248,250,254,  0,  0, 14,  0,248,  2,  8,255,253,  3,  4, 15, 12, 22,223,255,  1,  2,  2,242,251,251, 10,  5,246,  2,
       43,  4, 27,199,228,  1,187,220,  5,198,222,229,251,  3,234,218,227, 68,  3, 53,203,200,245,227, 45,214, 61, 18,221,244, 30, 51,229, 32,  3,214,208, 37,  4,  4, 19, 16,254,  0, 41,  0,  2,  6,176, 18,254,  4,246,194,  0,  1, 32, 30,  0,  1, 14, 85,  0,254,248,  1,  0,255,  2,  9,255,  1,255,254,  0,  0,  2,  0,  0,255,250,  4,255,  0,  4,  2,  0,  0,  0,  2,255,  0,254,254,  0,  0,239,203,  5,253,  5,235,233,204, 58,248, 32,245, 57, 16, 73,219, 21,241,245, 62,232,178,187,  2,225,229,239, 22,254,  0,242,130,233, 60,241, 12,239,236,244,194,250,244,244, 47, 84,190,207,250,  3,  9,193,232,  1, 23, 49, 26,227, 18,225,100, 51, 76,  7,133,252, 35,254,251,  8,  1,  2,240,  4,254, 39,253, 50, 43, 16,  0,  4,255, 49,243,253, 28,  0,  4,247,254,244,231, 12, 38,  0,233,
      254, 22,235,  3, 14,239,227,  4,189, 53,232, 60,  9,220, 37,250, 40,242, 50, 12,211,253,174,254,255,228,129, 40,215,222,223,215, 16, 16,  2,255, 66,210,  0,197,223,221,251,  0,  7,232,254,  0,234,210,  1,  0,213,238,254,255,194,244,250,255, 36,  3,254,255,199,253,255,  0,250,253,255,  0,255,255,  0,  0,  5,  2,255,  1,  0,252,255,  1,255,253,  1,  0,  0,  1,255,  0,249,254,  0,  0,  4,  7,254,234,255, 13, 44, 69,249, 33,  0, 13,213,215,195,250,251, 74,225,250,168,146, 58,174,210,  9, 16, 15, 60, 42,180,244, 25,196, 17,  2,195,205,186, 12,164, 14, 22, 24,214, 56, 29, 19,252, 65, 19, 11,242,154,  4,  7,215,254,  7,180, 64, 11,173,213,237,237,237, 24,228,231, 64,  3, 20,254, 28,242,252,198,  4,212, 10,240,252,  9,254,251,187,255, 20,246,247,244,210, 17,  2, 10,
        1,210,  7,245,223, 30,  7,253,  4,  7,  4,248,  6,244,235,188,218,226,217,254,245,  6,190,195, 36,232, 80,255,235, 14,253,245, 13,241,251,137,214, 56,253,  1, 58,231,  0,253,230, 26,  1, 11,255, 50, 18,  0,225,  3,  0, 17, 26, 17,  0,  0,240, 24,  1,  0,242,253,  0,  0,255,254,  0,  0, 20,  1,  0,  0,  3,  0,  1,  2,  3,  1,  1,  0,254,  1,  0,  0, 22,250,255,  0, 17,250,  1,  0,214, 26,210, 23,170, 20,254,  7, 16,222,251,  3,  0,249,249,145,  4,249,245,249,209,  0,  2,245,128,218,236, 13,223,189,  0, 84,247,235,249,  2, 63, 64, 12,245,253,218,  8,212,  3,  0,112, 65,248,244, 16, 60,217, 34,254,110,182, 83,240, 62, 66,194,238,248,249, 30,230, 49,  3,  9,192,252,  6,252,  0,255,208,  4,252,254, 19,250,212, 40, 19,  8, 66,254, 29,  4,232, 16, 12,255,255,247,
    },
    {
      248, 37,  2,  1,226,191,253,  1, 32,196, 39,  2,238,  0,241,  4, 16,247, 35,255, 47,  3,237,132,254,235,219,130,134,  0, 49,254,  0,246,250, 19,254,  5,246,  0,  2,238,236,255, 18,249,247,253, 15,255,244,253, 34,  9,254,  0, 19,  7,  3,  1, 17,  7,255,255,  1,255,  0,248,  3,  0,255,255,  1,  1,  2,255,  4,  2,  0,  0,  4,  1,  1,  0,  2,  1,254,  0,255,  6,  0,  0,  1,247,  1,  0,  6, 10,197,209, 54,196,225, 10, 17,221,246, 18, 50,246,114,221,201, 35, 26,179,  0, 51,239,193, 18,  6, 23,221, 59,156, 10,189,224,  6, 98,229,203,215, 59,  0,254,  4,192,  2, 16,  8,246,231,252,  8,248,207, 96,253,217,173, 21,239,  9,193, 13,235,202, 16,255,  4, 26,248,  5,  1,251,  1,  4,  0,212, 10,249,  6,  4,251,  0,255,254, 16,  0,  4,254,  3,253,  1,255,230,  3,235,  9,  3,
        9, 12,227,251, 33,250,155, 13,249,  9,242,  3, 51,252, 18,  5, 43, 64,  0,  1,222, 54, 23,122,251,229,  1,124, 43,250,195,  0,  0,247, 61, 41,254,242,  1,  0,  1, 22, 38,  3,  7, 12,254,  4,  2,  8,249,  0, 14, 29,254,  2,  3, 20, 13,  1,194, 10,  2,  2,  0,  4,  2,243,255,  2,  2,255,  1,  0,  3,  0,  2,254,249,  0,  4,  5,  0,  0,  1,  4,  0,  0,254,255,  3,255,  0,  1,254,  0,  1, 22, 18, 26,245,255,130,255, 41, 17,187,225,  3,247, 60,227,244, 29,233,232,161, 21, 33,205, 24,  2,  2, 94,170, 94, 27,233, 14, 11,193,250, 65,255,199,  1,  3,193, 55,243,171,251, 70, 31,223, 65,  4,231,183,196,236, 80, 14, 17,232, 24,156,  2, 70, 42,  0,252,  2,253,241,  6,  2,  0,253,221,235, 29, 12,  5,255,127,252,234,  4,254,254, 13,255, 20,  0,  1,242,235,252,254,251, 66,
      255, 13,252,247, 16, 46, 10,  1, 18, 47, 41,  2,  4,253, 55,  0,  4,246,253, 78, 41,  1, 68,  2,236,  8,209,  8, 16, 74,199,130,  0, 16, 10,209,  7,  0,206,  5,241,255,195,  0,255,254,249,  0, 52, 28,255,  4,193,  7,  2,  1, 26,232,252,253,  5,168,  1,  1,255,250,248,  8,  1,  0,248,  1,255,255,  2,  0,254,255,  2,  0,244,  1,  0,255,247,  1,  1,  0, 13,  1,255,  1, 17,  1,  1,  0, 37,249,249,253,225,247,254, 20,  3,242,221,219, 10,167,142,  2,228,254, 33, 41,205,171, 69,  4,180,248,189,193,193, 12, 37, 17,249,  4, 21,161, 50, 76,235,  6, 16, 24,237, 11,228, 39,189,  0, 18,253,217,255, 10, 30,249,121,251, 20, 12,247, 22, 93,  4, 62,  1,  3,221,201,249,  8, 29,255,255, 30, 28,201,  3,245,  1,118,255, 21,245,  0, 96,  1,  1,214,  1,  2,  7, 36,253,255,  0, 75,
       25,215, 19,  2,255, 42,253, 15,242,241,204,  3,206,  3,218,  0, 34,238,179,  0,254,241,  5, 63,207,  3,  7, 66, 19,252,228,  2,  3,250,  0,233,254,247,251, 18,  9,214,  1,  0,  6,  1, 53,  0,253,218,248,  1,  2,253,  0,  1, 39,206,254,  3, 20, 65,246,  0,  0,255,  0,252,255,255,  0,  1,  1,  1,255,255,  3,  2,254,  0,248,244,  0,  1, 20,  0,  0,  1,234,  5,  0,255,210,253,  1,  0,236,247,242,  9,249,127,209,  7,250,  6, 38,  1,106,  2,209,193, 42,205,254,  5,156,  0,  2,244,222,247,  1,  2,241,217, 70,128,254,195, 38, 29,253, 35,210,  0,252,234,234, 74,199, 60,195,254, 56,  2,252, 23, 10,  0,254, 39,225,236,254,252,233,163, 31, 14,  0, 14,239, 48, 78,243,250,  0,  0,  6,235, 65,249, 55,  2,  1,  2,  0, 29,253, 73,  6,  2,192,  1,253,250,252, 68, 11,  6,253,
    },
  },
  {  //  11
    {
        6,207,212,  2, 12,207,240,252,202,  0,186,249,255,250,251,255,254,227,251,247,253,234,252,247,251,223,253,253,248,  4,  9, 12,  2,246,249,  1,236,239,253,250,232,226,241,  0,  6,249,255,250,237, 19,  0,253,237,  9,247,249,241,237,254,255,187,249,255,249,  1,252,  0,  4,  1,255,255,  0,  2,  2,255,  0,  4,  0,255,254,  6,  2,255,  1,  9,  4,255,255, 12,  7,255,255,  4,  7,254,  0,255,242,248,  7, 19,196, 16, 61,232,232,196,240,227, 57, 18,104,231,  8,205,230,175,222,238,209,148,207, 11,203,148,211,216,181,  3,250,  0,231, 54, 36,  8,205,  0,221,235, 60, 58,255,  4, 63,233, 12,239,223, 12,212, 74,255,  0, 42, 56,219,  5,221,201, 23,  0,  5,  4,  7, 13,  6,254,248,  0,  5,255,252, 17, 16,254,165,  5,246,  5,236,252, 27, 38,254,247,  0,  6,241,  0,249,249,  1,
      248, 17, 36,191,222, 18,  8,119, 60, 59, 56,  2, 12,  9,220,233,251, 44,  0,224,  2, 37,  0, 76,  2,244,218, 67,252, 60, 16,192,  5, 15,  2,253, 26,242,255,  4, 38, 14,247,  8,247,206,  6, 10, 44,245,  1,  2,223,231,  1,253,239, 44,255,  3,  7, 35,  2,  4,  0,  0,  0,  1,  2,255,255,  2,  0,  0,  1,  0,  0, 10,255,  1,  1,  1,  0,  0,  7,  1,  0,255,  0,253,  0,  2, 11,  1,  0,254,222, 61, 62,252, 20, 87,246,  6,240,  4,209,225,  1, 18, 39, 65,243,  3, 12,255,211,234, 60, 33, 22,245,229, 35,130,229,226,250,  0, 18,236, 30,223,182,247, 54,255, 13, 29, 14,171,246, 11,185, 21, 29, 69,246, 49,100, 32,232,  8,209,198, 12,192, 57, 33,134,255,254,  4,251, 30,  6,247,  0,255,246,252,249,  3,242,  2,  8,  2,229,248,249, 62,247,  0,  0, 12,  1,  9,250,  1,249,254, 13,
      234, 20,170,253, 53,187, 64,189,  6,205, 12,243, 38, 17,  5, 22, 44, 15,  2,170, 65,173,241,255,195, 10, 13,191,231,255,253, 69,253,255,254,235,241,252,255,  1,228, 39,255,250, 30,  5, 16,  0, 45, 24,252,  4,241,245,252,254,  8,239,  9,  2, 38,238,  6,  0,  0,  0,  1,251,  4,  1,  0,253,  7,  2,  1,  0,237,252,  1,  0,  1,  1,  0,  1,254, 27,255,  1,  1,254,  0,  0,249,255,  0,  0,  0, 20,  9, 13, 13, 51, 34,  1,246,  1,232,240, 15,221, 19, 78,197,  0, 15, 16,140,  0,205,206,  2, 14, 85,199, 64,  4, 36, 85,255,238,  3,197,  0,242,  8,200,251,254,255,200,  2,219,  6, 24, 45, 76, 49, 30, 51,175,248, 68,237, 14, 20, 81,208,196,104, 65,  3, 25,220, 41,  8,252,251,255,  3,  0, 11,255, 11, 17,248,  0,  5, 32, 19, 15,127,255,254,  0,253,  2,246,236,  2,  0,  5, 60,
      254, 60, 77, 62,192, 40,168,191,239,237,  6, 25,252,237,247,  4, 64, 35,  7,241,  5, 93,247, 31, 31,234,  3,251, 67,248, 18,174,  2,237,255,  2,  1,248,250,  1, 13,  2,  1,213,197, 12,248,254,242,249,  0,  2,209,219,  3,  1,203,  7,  2,  8,238,  3,254,  0,  1,255,254,  0,  0,253,  0,255,  8,254,252,237,248,  2,  1,  0,  0,  0,  0,  0,  0,253,  1,  0,249,  1,255,  0,  9,238,  1,255, 23, 10,235, 17, 49,  2,175,204,  1,222,  3, 29, 22,250,227,127,239,253, 48,  2, 47, 27,158,233,233,253,211,217,236,185,137,253,  2,239,242,204, 44, 38,  1,188,  4,250,244,181, 62, 76,228, 51,240,180,  4, 29,126,215,225,185,239,227, 68, 30,252, 32,  0,194,  2,237,249, 24, 59,  1,254,255,255,  1, 10,234,241,245, 13,182,  1, 40,226,217,  6,  1,247,253,  3,237,252,235,  3,221, 28,133,
    },
    {
        7, 11,255,250, 26,241,  3,199,  8,254, 25, 71, 29,248,  5,204, 19,  9,  6,  9,250,  4,  5,  2, 13, 25,  4, 23,250,  4,  0, 11,229, 12,  5, 26, 14, 13,  2,251, 15, 18,255,  7, 31, 12,  1,  3,245,  0,254,253,246,227,254,251,225,247,254,251,  2,245,254,253,255,  0,  0,237,  0,  0,  1,254,  3,  0,  0,  1,  4,  2,  0,254,  0,  6,255,  0, 17,  9,255,  0, 17,  6,255,  1, 41,  6,  0,  0, 35,185,241, 10, 23,238,  4,251, 19, 44, 19,229,250, 96, 34,229, 67, 19,236,216, 23, 42, 14,  0, 51, 69,  6, 20,  8,  3, 78,220,236,203,  2, 61, 63,130, 34,  2,192,228,219,235,232,134, 56,252,  4, 11, 31,  0,243, 17,  0,244,  4,  0,249,  3,166,  5,139,237,254,235,249,240,205,250,248,255, 36, 10, 22, 25,248, 20,  2,255,242,  1, 12,  3, 32,134,255,255,247,254,255,251,253,  7,  8,249,
      229,190,  4,197, 25,  4, 23, 56,  1,239,243,100,233,240,238, 68,176, 29,239,250,238,252,  6, 98,213,  8,239, 64,119,233,  0,173,241,223,252,  2, 12,236,  1,255,242, 35,  1,253,  5, 10,  2,249, 10,249,  5,255,249, 29,  1,255, 14, 15,  0,  0,239,223,  1,255,  9,  1,253,  0,253,255,  0,  0,  2,  1,  0,  1,  2,  0,254,  0,240,251,  0,  0,252,255,  1,  1,  0,  2,255,  0,254,255,  1,  0,248,234, 19, 10,245,134,250, 11, 26, 18,  2,215, 14,247,  3,255,  6,218, 35,255,  9,250, 40,112, 20, 38, 43, 56,239,249, 13,131,244, 46,225,223,242,131,  3, 44,246, 30, 22, 24,  2, 26,195,187,242, 27,239,  3,  1,144, 66, 64,  1, 11,213,246, 33,249, 14,255,253,246,  5,252,  2,129,  0,254, 23,254,249,245, 23, 22, 57,249,  2,243,243,250,220,232, 33,207,  6,  2,  0,249,  9, 41,255,  2,
       19, 40, 10, 64,245,  0,  1,247, 17, 55,238,109,209, 19,  5, 19,232,  3,  7, 34,  2,  2, 10,254, 48, 32,240,146,211, 33,  0,189, 38, 21,255,209,207,  9,255,193,236,  7,255,  3,229, 18,  3,254, 15,235,  0,  0, 42, 57,  1,  3,248, 63,254,  2,249,255,  1,  0, 33,  4,  0,  2,  5,  0,  0,  3,251,255,  0,  0,250,  0,  0,  0, 17,  2,  0,  0,255,  4,  3,  2,250,  1,  0,  0,  5,254,255,254, 14,220,250,237,  0,  1,190,  0, 19, 24,254, 49, 66,255,255,  2, 39, 43,  6, 62,142,226,194, 62, 29, 27,251, 55,219, 13,  0,198,247,246,219,240,255,253,229,255,192,218,229,227,223,116, 94,255,  2,248,  3,255, 55, 38,218,255,244,  1,  9, 16,  8, 15, 18,246,250,243,  1,254,251,255,214,  1,218,231,250, 81, 18,253,  0,  0, 16,  7,  6, 31, 25, 21,222,  1,  4,249,  3,245,252,175,  1,  7,
      245,242,252,  4,239,  4,  0,247,245,222,244,  3, 37,212,245, 14, 57, 16,  3,229,218, 42, 31, 58,206,241,  4,243,227,  5,255, 62,216, 24,  1, 38,224,250,  0, 68, 25,209,254,254,  7,208,255,  1,246, 12,  0,  2,247,228,255,  2,251,237,254,255, 33,  4,  0,  2, 37,252,  4,224,252,  0,  1,  0,  2,  1,253,255,254,  0,  1,  0,240,255,  1,254,253,255,  0,  0,  6,  0,  0,255,233,  1,  0,  0, 36,228,246,  1,  3, 48, 48,  0,229,216,160,224,219, 57,235,237,236, 43, 19,250,  1, 73, 12,  4,220, 10,  7,214, 85,  3,113,113,208, 20,255,251,  0,245, 27,  1,  7, 26,241,253, 21, 72, 43,  3,213, 25, 44,222, 11,218,207,213,  7,231,  0,209,211,255, 55,208,252,235,  0,  0,  5,255,249,  0, 12, 37, 10,227,  8,  7,  6,252, 16,247,  5,  0,252,228,254,201,  2,  0,  5,  8,  6,  0,  1,254,
    },
    {
       18,189,  8,175,  2,  9,249,250, 20,233,  1,  0,239,218,  5,204,251,233, 28,  0,236,253,255,201,246, 20,255,  0,255, 43,253, 11,  8,241,  5,233,  3,  3,253,254,252,  5,252,255, 15,  4,  2,238, 37, 10,254,253, 36, 12,  0,253,254,237,253,248,  1, 22,253,242,  0,  0,  1,  4,  3,  1,  2,  1,  0,  2,  1,  0,  1,  3,  0,255,  4,  2,  0,  0,  9,  3,  1,  1, 19,  3,  0,255, 28,  4,  0,  0,229,249, 20,253, 25,227,  6,  1, 20, 12,251,225, 46,253, 83,  4,214, 15,218, 23, 18,240,  5,120, 30,  5,  2,223, 11, 51, 26, 41,240,  7,242,249,249,238,  2,  1,233,219,247,250, 28,245, 13,253,217, 13,184, 74,206,191,  0, 75,239, 26,230,  8, 68,217,182, 59,  0,  5,  4,246, 12,  0,251,  0, 33,254, 15,  2,251,  8,  2,  8,216, 29,  0,234, 38,  7,  5,  5,  0,247,  4,244, 14,253,  0,  1,
        3, 67,  0, 37, 22, 31, 17,254,218, 22, 32,238, 41, 24,246,236,201,246,249, 64,226, 42, 23, 17, 60,250,239,242,247,220,253,244,241,253,186, 44,251,  2, 15,243,252,249, 10,255,210,248,  4, 27,247,253,  2,197,254,  4,254,249, 37, 23,  0,  4, 28,228,254,  6,255,  2,  3,  0,  1,  0,255,255,  2,  0,  0,  1,  9,  1,255,  1,  1,  1,  0,254,  0,255,  0,  0,250,  2,  0,  0, 17,  0,  1,  0,255,201,  8,222,  6,236,  1,  2,198,233,224,234,187,245,190,192,246,201,251, 70,212,  1,205, 23, 17,208,243, 11, 43,245,  1,240,  1,  5,253,244, 46,249,  1,255,243,207, 20, 67, 13, 36, 23, 29, 34, 12, 82,239,186, 63, 91, 12,  2,219, 40,  9,244, 40, 51, 16,254, 26,248, 19,242,251,255,  0,226,238,254,  7,246,250,254,253, 35, 43, 19, 17,  0,196,254,  4,  3,255,  3, 20,  7,251,  0,250,
       12,  9,241,  7, 21,225,242,  8,246, 43,  3, 42,  5,241,238, 13, 18,  8, 10, 21,  0,221,193, 80,230,243,240, 23, 26, 11, 11,201, 12,  6, 58,162,  7,  1,249,  3, 38,246,255,165,251,  8,244,  6, 55,  8,  0,255,223,249,252,  1, 22, 10,  1,254, 24,  0,  1,255,  0,  0,247,249,254,  0,254,254,254,255,  0,  0,  0,  0,  0,  0, 14,  0,  0,  1,  1,253,  0,  1,255,  2,  0,  0,253,  2,255,  0, 10, 59, 31, 43,207,  4,254,249,244, 50, 20, 10,  1,  4,239, 65,243,185, 27,214,250, 44,173,209,188,  5,230,252,202,187,179,216,  6, 54,247,  3,206,254,  1,  3, 24,250, 53,252, 21, 56,241,249, 36,236, 31,238, 35,236, 16,201, 29,  3, 18, 57,246,208, 77,  1,  0,255,255, 12,  4,  0,  0,  0,  4,231,245,182,  3,249,255,253, 22, 28,253,243,232,249,  0,250,242, 55,253,254, 27,  8,  2,254,
      239,217, 51,250,233,  1,238, 21,  3, 48,255,252,  1,232, 11, 60,186,222,254,160, 21,  1, 50,200, 11,  0,252,254,223,  5, 20,  5,  0,  5, 50, 71,  2, 15,  1,253,221, 20,  3,  2,  0, 10,253,  4,223,209,  1,  0,242,252,  1,255,238,237,255,  2, 18,  9,  2,  1,255,  1,  2,  7,255,  2,  1,  0,  2,  1,  4,255,254,  2,255,  1, 14,  4,  0,  1,255,  0,  0,  0,214,  0,  0,  0,  1,255,  0,  0, 27,  5, 52,248,245,254,234,250, 28,  0, 11,238,255, 13,248,235,216,223,  0, 64, 60,237, 42, 32,253, 28, 66, 22,249, 59, 62,203,232, 37, 46, 14,  1,252,220,  3,248, 20, 25,246, 45,210,193,253,  3,232,221,194, 63,216, 21, 85,  5,241,240,189,201,220,212, 49,253, 20, 26, 32,246,  3,255,  1,252, 42,250,254, 88,237,229,  0,  0,229, 22,251,254, 64,255,247, 40,209, 21,249,  5,248,  1,  0,
    },
  },
  {  //  12
    {
      253, 27, 15,241, 20,  9,  6,242,241, 32,238, 41, 33,237,  6, 10, 30, 14,250, 14,194, 16,  9, 71,  3, 40,  1,  9, 11, 20,233,232, 48,  5,  0,248,240,  8,  1,254,242, 14,252,  0,253,  6,253,  0,190,  2,254,254,  5,250,254,253,255, 30,249,255,204,250,  2,249,  0,  0,  0,  0,  0,  1,  0,255,  1,  2,  0,255,255,  2,  0,255,  5,  3,  0,  0,247,255,  0,  1,  1,  1,  3,  1,  1,  3,  0,255, 42, 19, 20, 21, 22, 23, 31,  3, 24,229, 11,246,193,246, 32, 33,241,236, 15,214,249,163, 58, 14, 22,241, 26, 89, 63, 23, 60,191,243,252,  3,248,254,248, 14, 10, 29,239,254, 27, 33, 25, 14,255,254,  0,253, 47, 65, 84,  1,198,  9, 28,  4, 20,238, 28,211,192,249,254,248,255,  0,255,253,  0,  0, 10,  1, 27,249,  1,  4,  1,250,  1,253,247, 25,253,253,255,255,252,251,  8,  3,251,219,  7,
       18, 14, 42,207,253,  2,181, 19,247,242, 14,228,242, 32, 26,235, 24,  4, 26, 54, 10, 57,  8, 12, 32,249, 50,  8,251, 17, 14,251,253, 12,  8, 63, 43,  5, 74,  1,244,254,252,  4,  1,  9,255,  1, 49,  3,  7,252,  6,248,  3,255,  8,  9,255,255, 29, 12,255, 10,252,  3,  0, 19,255,255,  0,254,253,254,  0,255,  2,255,  0,254,  0,  0,  0,254,  9,254,  0,  0,  2,  1,  0,  2, 17,253,  1,  3, 24,  0, 33,255,  6,249,236,  0,  1,213,  7,  5, 17, 45,217,  3,  4,  0, 48, 65,234,223,193,188,235, 30, 45,194,198, 48,207, 35, 21,  6,248,  4,254,  8,215, 63,249, 25,251, 18,192,212, 53, 16, 74,215, 41, 78, 50,214,247,254,  9,234, 67,116,244,207, 93, 30,  0,  1,  0,255,  2,254,  4,  1,  2,223,  0,  2, 11,236,  4,255, 29,245,  5,  0,253,251,253,  0,  0,  4, 46,244,206,247,  0, 16,
        3,238,250, 17, 21,  2, 20,254, 27,246, 39,  1, 30,248,180,196,240, 87,239,251,249,247, 10,194,  0,  9, 52,215,238, 12,  5,255,255,238,  2,191,196,251,  4, 67, 34,  0,221,201, 23,  2,  4, 38,  5,242, 56,  3, 49,  5,  0,  1,  2,251,243,254,  8, 25,253,250,  3,  0,  1,251,255,  0,255,221,  0,253,  3,231,  1,255,  0,  4,  0,255,245,  0, 25,253,255,255,  0,  0,  0,  0,  2,  7,  0,254,  7,240,228,254,  3,  2,248,  1,250, 25, 32, 21,  9,254,  2,245,  8, 11, 18,254,255,253,253,201,249, 19, 21, 12, 66,229, 26,205,239, 16,218,250, 16,193, 41,251, 39, 50,  9,248, 39,255,238,  2,201, 20,247,167,248, 41,254, 64,  9,  1,213,216,194, 17,  3,197,  0,  7, 29,252,  1,  0, 10,253,  0,  9, 10,254, 36,  7,249,  0,  5,  1,  3,254,213,  5,253,  0,  0,  2,  4,  7,251,231, 10,  0,
        8,251,214, 74, 17, 12, 61, 71, 25, 14,230,  8,  9, 12, 32, 64,  9,  9,248,253, 66,  0,193, 56, 11,199,  1,116, 34,238,249,254,255, 21,  7,252, 11,  7,190, 92,240,247, 44, 62,245,  0,248,  0,252,231,191,  4,214,237,  1,255,227, 15,252,255,208, 62,255,  0,  0,  2,  0,254,254,  0,255,  2,  0,253,  0,  3,  0,  1,255,  2,  1,  1,  0,  1,252,  5,255,  1,  2,  0,  0,  0,252,216,  2,  0,245,218, 31, 18,  7,195,  9,247, 43,230,234,  2,  1,  1, 52,251, 20, 16,240, 21,  6, 19, 32, 28, 36,  1,245,  2,253,104,153,202, 27, 19, 42,  3,  1,247,  1,  8,242,240, 33,  9,203,251,233,  5,241, 13, 22, 40, 11, 59, 49,113,  6,  6, 10,201, 76, 49, 31,182,245,226,248,  7,254,219,  1,255,  1,  1,240,  1, 28,  0,  8,  0,  1,232,  0,  5, 49,255,254,  0,254,255, 14,230,  2,253,  5,210,
    },
    {
       27, 16, 21,253, 16, 17, 12, 94,  7, 13,  7, 53, 12,  9, 16,  2,242, 17,  6, 14,255, 36,  5,  9,248,209,  9, 10,255,222,  5, 13,247, 14,254,187, 10,  1,  0,  2, 13, 13,  2,255,238,  2,252,  1, 28,252,251,  2, 22, 11,249,254,240, 23,248,251,235, 26,252,253,255,254,  0,250,  1,  3,  0,  1,  2,  2,  0,  1,  4,  4,  0,  0,  1,  2,  0,  0,  3,  3,  0,  0, 14,  4,  0,  0, 15,  3,254,  0,228,253, 46, 10,250,231,246,249,  9,  6,232,194,248,245,205,236, 85,253,  2,212,245,182, 34, 76, 15,233,250, 14,250,250,189,  7,240,235, 58,251,251,214, 41,  2,243,196, 10,218, 70, 64, 40,225,224,245, 21, 21,  7, 25,250,228,252,243,  0,  4,  7,247,232,  9, 12,255,223, 15,254,241,248,  0,  2,232,  0,246, 22,247,247,  5,  3,243,  0, 13,244,253,  0,  0,  7,  0,252, 20,254,  0,  4,253,
      239, 40,252, 71,231, 15,231,192, 49,  1,248,230,217,255,250, 31,222,253,251,  1, 27,208,248,252,235, 61,249,247,191,251, 37,  1,  8,236,242,254,243,  4,  5,247,239,  1,  3,249, 50,  3,252,251,227, 14, 61,250,246,245, 66,  1, 20, 15, 26,  4,238, 56,255,  0,  1,255,  1,  1,255,  0,  1,  0,255,  2,255,  0,  3,  2,  0,  0,  0,  1,  6,  0,  0,253, 14,  0,255,255,  1,  0,  0,250,255,  0,251,213,218,238, 34,  8,252, 41,233,223,255,254,252,226, 86,100,224, 32, 71, 48,184,232,198, 18,  1, 24,  7,237, 28,  3, 55,247,195,221, 22,240, 28, 36,255, 26,240, 30,  2,251,223,201,227,  3,239,230, 33, 14,252,  2,  2,255, 12, 16,254,253,255,252,  2,242,225,  5,  5,231,255,  7,  1,  1,  2,  0,214, 18,  0,245,  1,251,  1, 19,236,250,255, 10,253,  1,  5,  8,255,228,255,  0,254,255,
      232,230,248,196, 27,237,  0, 16,203,  5, 12,189,247,  5,  9, 46, 38, 26,255,  0,199, 27,252,  1,  1,247, 29,250, 42, 11,230,  1, 40,255,  9,241, 22,253,253,  2, 11,251,255,  8,216,250,  7,247,203,252,193,  1,  9,  2,198,253, 47,243,255,251,230, 15,239,253,  1,  0,  1,250,253,  1,  0,  0,  1,  1,  0,  1,253,254,  0,  0,253,  0,  0,  0,251,  2,255,  0,253,  0,  0,  0,  0,  2,  2,  0,244,237,223, 73,255, 15,146,  5,212,241, 14,131, 52,232,197, 20,248, 33, 30, 15,255, 80, 60,180,254, 51,250, 49, 12, 32, 62, 26,  2,234,  7,233, 56,185, 48, 66, 10,228,251,119, 15, 17,253, 57,  5,243,  1,191,  0,  4,  1,  2,221, 77,238,  9,  0,247,252,246,228,244, 22,  1,  5,236,  1,  7,  8,243,  4,239, 25,225,  0,255,  1,211,213,  6,254,  0,253,  1, 32, 58, 15,255,196,  2,  0,  0,
       38,187,  3, 16,247,  2,  7, 76, 58, 33,242,  8,230,246,255, 39,  1,  0,243,247, 18,238,  0,249, 14, 26,221,249,245, 22,255,  3,186, 40, 10,  1,235,  9,  3,  0, 59,  4, 12,213, 21,254,255,250, 32,  5,252,253,250,  3,255,254,220,  9, 23,  2,  7,251,255,  4, 53,254,  0,  0,  2,254,  1,  0,  2,254,  0,255,  0,254,  1,  0,  0,  0,  0,  0,  7,255,  0,  0,  1,255,255,  0,255,  0,255,  0,255, 28, 74, 54, 10, 11, 32,  9,  6,221,250,132,253,  5,104,251,255,  5, 18, 30, 57, 36,254,255,  7, 18,223,  8, 59, 10,219,  3,  1,245,251, 40,187, 77,220,192,253, 31, 26, 79,250,211, 53,237,  8, 23,253, 50,222,  5,242,248,245,210, 30, 16, 44,  3,252,251,  3,  9,247, 15,  2,251, 36,255,249, 16, 18,177,221,  6,  4,255,217, 12,  4, 24,  4,227,252,  0,216,238,242, 25, 56,  5,254,  0,
    },
    {
      217, 15,  4, 27,238, 17,  1,254,209, 11,  1,251,244, 60,  9, 10,207, 17,226,173, 57,250, 14,  0,238,216,  1, 34, 49, 59, 52, 42,  4, 28,242, 32,251, 34,249,255, 13,254,246,250,  4, 25,246,253,  8,254,254,240, 45,  8,251,253,  2, 46,254,  2,209,  5,247,254,  5,  0,  2,  1,  1,  1,254,  0,  4,254,255,255,  1,  2,  1,  0,  2,254,  1,  1,255,254,  0,  0, 13,254,  2,  0,  7, 17,  0,  0,  1,244,236,242,239,231,  6,249,250, 83, 58,  9,185, 37,216,  2,243,121,235,208,217,  5,147,  3, 46, 66,230, 51, 16,239, 16, 44,  7,  0,  2,  2,  1,  2,  3,  0, 24,247,233,250, 24,227, 10,254,187, 49, 49, 21,218, 42,186,246,231,233, 13,235,239,133,111, 10,  1,  0,  2,255,  0,  0,255,  0,  2,254,255,255,254,  2,  2,  0,232,250, 19,  3,247,  3,254,  2,  1,  0,251,  6,255,255,  1,  1,
       60, 53,  7, 53,  7,240,202,151,245, 75, 17, 13,233, 51,223,  1,236,247, 59,180, 19, 55,  9,253,253,212,229, 27,244, 27, 36,165,226, 14, 49,189, 10,228,247, 69, 25,  7,244,  0,254,225, 12,253, 73, 23,121,255, 48,250,252,251,  3, 22,250,251,206, 10,  3,  3,  1,  1,232,  0,255,255,255,237,  2,254,  0,  0,255,  3,  0,  0,  0,  1,  0,  1,  0,  2,  0,255,  0,254,255,  0,242,  2,254,  0,255,238,253, 24,248,255,254,  6,  2,216,220,201, 46,202,255,229,225,171,255, 54, 21,194,217,129, 62,195, 48, 12,251,196,105,131,  2,254,255,249,  7,  1,  0,  5,216,  8,  3, 10,245,255,237,  0, 58, 89,197,230, 11,168, 18,140, 57,  9,  3, 27,210,133,195,214,  0,  0,  1,255,  0,  0,  0,  0,  8,  0,  4,255,255,  2,  3,  0,  6, 13,  3,230,255,  2,255,  1,  4,254, 80,  8, 11,253,  4,236,
      245,142, 80,252, 44, 61,147,255, 11,254,233,254,229,246, 16,254, 18,197,243,214,  9,174, 91,110, 19, 60,202,  4,  2,255,195,245,  2,251,197,245,248, 16,117,194,227, 49,204,  6,  1,  0,184,  0, 23,191,  0,128,198, 22,  1,  2, 64,  4,254,  5,247, 14,246,  3,  8,  3,254,  2,  2,254,  0,  2,  0,  0,  8,  0,  0,255,  2,  1,  0,  3,255,  0,244,255,  0,  0,254,245,255,  0,  0,  5,  0,  0,  4,244, 83,228,  1,247,236,234, 55,  0,235, 45,234, 74,197,194,181,225,250,215,145,210,164,  3,252,206,223,  7,108,186,169, 66,  2,240,255,252,  0,  6,255,  0, 13,247,254,237,252,252,  4,246,255,251, 32, 30,127, 25,181,  0,240,221,244, 29,  3, 10,149, 58,  0,254,  1,  2,  0,  0,  1,  1,255,  1,  3,  2,  1,255,254,  0, 13, 37, 17, 17,  6,  6,255,231,255,  2,196,  1,  1,248,  0,248,
      255, 40,255,253,239,205, 68,240, 12, 17,252,  1,214,190, 11,255,207, 85, 62,144, 49,  1,  5,229,248, 14,212,207,  4,150,201,  4,  4,  0,255, 14,255,253,  5,233,  1,252, 78,252, 20,  2, 57,  1, 49,212,  5,  0, 14,235,  2,255,218,250,  1,253, 11, 54,246,  2,246,254,255,  1,  2,  1,  1,255,255,  2,  0,  1,  3,255,  0,  0, 14,  7,  0,255,249,  0,  0,  0,  8,255,  0,  2,  7,254,  0,  1, 16, 32,244, 77,  4,  4, 63, 18,  9,255,193, 37, 57,160, 27, 82,195,200, 28, 74,171, 61,242,239,241,254, 47,254, 92, 87,  8,217,  2,  7,  1,  0,  4,  8,212,255,250,247, 14, 22,  1, 26,242,  4,  1,206,255, 40, 25,207,255,140,251,248, 30, 26,207, 51, 65, 59,  0,  7,  2,255,  0,  1,  8,255,  1,255,235,  8,255, 12,254,253,  3, 17,228, 43,  8,  6,  1,  7,  7,255,254,  7,216, 89,  4,  0,
    },
  },
  {  //  13
    {
      235, 50, 58,  9, 16, 67, 23, 20,  9,  1,242, 96,  5,  9, 71, 53,201,197,245,201,  1, 21,226,236, 21, 33,  8, 32, 51,165,255,  8,226,253, 72, 67, 56,  0,216,240,  4,  8,231,254, 37,254,232,249,240,247,  0,253, 40, 43,254,253, 64, 10,255,252,210,253,254,247, 17,  0,253,249, 19,  1,250,  0,254,255,255, 15,  2,  3,  1,254,  3,  1,255,  1,  3,  7,  1,  0,249,  3,  0,  0,253,  1,255,  0,  8,255,253,237,  6,233,248,  4, 27, 44, 29,252,246, 30,  4,  3,168,248,246, 24, 41, 70,  1,246,200, 11,188, 18, 46, 27,218,255,  0,  1,250,  0,254,  4,  5,  1,  0, 13,249,213,252,  6,  1,252, 58,243, 53,223, 15, 13,200, 46, 62,245,246,244, 21,186, 61,  2,255,  0,  1,  0,  0,254,  0,  0,254,  0,253,229,  1,  0,  0,  3,  0,250, 10,  1,  8,  0,  0,253,250,255,254,252,251,248,  0,  0,
      244,192,254,195,  4,  1,  1,  0,252,  3,  0,146, 21, 19,209,198,  2, 41,241,253, 38,237, 60,108, 21,245,  9,240,237,243, 44,  6,249,222, 17,234,201,  1, 38,  1,238,  0,240,239,213,255,252, 18, 46,  9,254,255, 63,202,255,  1, 10, 31,  0,  1,227, 21,  3,255,241,255, 10,  1,250,  0, 13,  0,  3,255,  0, 15,254,  2,  0,  1,  4,  2,  1,255,  2,251,  1,  0,251,  1,  0,  1,250,  5,  1,  0,  6,248,220,  0,249,251,  4,  0, 58,226, 35,  7,240, 52,140,221,229, 36, 24,209,189, 16,  3,185, 48,249, 46, 52, 34,  2, 65,191, 12,  3,248,  0,  5,  7,255,253, 18,202,242,241, 75,  6,  6,222,181, 86, 11,254, 69,240, 63,253,209,249,194,187, 53,245,200,212,  1,  1,254,  0,  0,  1,  0,  0,254,222,  0, 22, 19,255,255,  7,  0,  3, 13,  5,255,  1,  0,  0,  4,  1,254,  4,  0,250,  0,  1,
       32, 60,  9, 11,  0,255, 52,143, 37, 14,  6,186, 18,  8,194, 66, 55,231,223,184,249,221,239,176,205,227, 57,237, 13, 18,105, 62, 47, 37,230,195,249,254,251,252,253,239,236, 14, 69, 25,  7,253, 25, 12,  2,  2,248, 11,245,255,  3,  3,253,255,216, 68,255,  0,251,  2,247,252,255,  1,  0,255,252,  2,  0,252,255,  3,  0,  1,255,  6,  1,  2,248,253,255,  0,255,254,  0,  0,  4, 10,  0,  0, 12,253, 45,  0,244,194,  0, 11, 25, 18,  7,216,228,252,253, 36, 19,206,218, 45, 68,  5, 60,195, 11,241, 33, 22,201, 13,118,103,190,247,  3,  5,253, 40,255,  4,  3, 39,251,236, 13,239,195, 37, 49,242,  0,206,140,203,  2, 55,237,250, 17, 20,201,223,252,253,246,  0,  0,  2,  2,  1,  0,  0,230, 39,254,  0,254,  2,  3,  3,  1,  4, 22,253,253,253,  2,255,253,  1,253,255,  3,  5,255,  0,
        6,224,244,251,  3,253,221, 64,  2, 11, 61, 60,238,  0, 63,143, 33, 30,255,129,255, 57, 65,202, 41, 34, 25,183,227,216, 12, 36,  8,229,197, 16,  5,  1,  2,  6,  1, 22,243,  7,206,249,  1,248,  3,  8, 11,254,255,  0,254,  2, 35, 33,255,  0,200,183,255,254,255,  0,  2,  4,  1,  0,255,255,  0,  0,254,255,  1,  3,  2,  1,254,255,  1,  0,  0,  0,254,  0,254,253,  0,255,244,  5,  1,  1, 10,252,240, 10, 21, 14,248,245,  5,177,  8,230, 64, 45, 16, 22, 13, 57,230, 74, 38,168,254,176,248,212,188, 11, 42,255,  6, 40, 60, 11,249, 45, 19,207,246,  2,253, 18,191,205,251,217, 36,237, 15,215, 30,150,134, 59,249, 56,  0, 36, 85, 75, 36, 72,  0, 28,253,243,  0,253,  0,251,  2,  0, 24,255, 47,247,  0,247,  2,255,254,223,252,235,251,251,251,  0,  1,246,251, 23, 16,254,  0,  0,
    },
    {
       18,  4,  2,250,222,236,  4,  8, 13,  5, 22,248,219,  5,243, 14, 55, 14, 19, 13,244,227,  5, 11, 23, 18, 32, 21,229,251,238,245,249,  2,249,185, 24,  8,253,245,  5,  6,246,  0, 20,  2,245,  1,  2,  0,245,252, 21, 16,251,251,  7, 12,252,252,234, 14,250,251,253,254,255,  1,  7,  1,  1,  0,  1,  1,255,  0,  3,  1,255,255,  2,  1,255,255,  1,  1,  0,  1,  3,  4,  1,  0,  4,  3,  0,  0,255, 19,  9,  4,245,  5,246,253,  9, 50,  1, 21,  5, 14, 39,243, 19, 20,  8,214, 21,155,245, 19,230, 38, 23, 59,218, 29, 14,253,  7,243,148,247,253,243,  4,255,246,244,236,  0,  2,248,240, 12, 49, 23, 17,247, 41,211,  4, 29,213, 63, 25, 47,255, 17,  0,208,248, 20,247,250,  7,255,  0,  0,254,250,  4, 13,249,247,  1,  2,  4,  8,231, 19, 13,242,254,253,235,250,  4,105, 22,237,111,  0,
      248,245,  3,  8, 25, 17, 11, 15, 49,249,253,228, 37, 10, 35, 11,219,247, 18,244,190, 36, 12,249,186,222,250,233,244,249, 80,255,245,  1,  2, 83,239,240,254,  0,  2,  4,  4,255,233,250,255,254, 21,253,  9,253,250,248,247,  1,  2,  9,254,  1, 11,252,  2,  0,  6,  0,  2,118,254,255,255,  0,253,  0,255,254,  0,  3,  0,  1,  1,  1,  0,  1,  1,255,  0,  0,  2,  1,  0,255,  1,  1,  0,255,165,  6,  8, 34, 49,  7,  3,253,  5,245,  6, 13,  2, 13, 11,240,  9,255, 56, 41,210, 29, 34,167,206, 41,241,254,250,232, 29, 19,  2,249,237,  8,  5, 26,253,255,253,232, 24, 26,  1,252,255,243,  6,206, 28,200,243, 45,179, 16,246,225, 63, 43,250,248,252, 25,  4,234, 11,247,249, 10,255,  0,  2,195,236,154,253,192,  6,  2,237,  4, 73,246,245, 36,  1,  3, 49, 11, 31, 64, 30,240,152, 49,
      218, 19,  2,243,248, 28,255,228,252,249,246, 34, 17,247,236,255,211,251,  4,  0,  6,241,249,251, 20,254,246,  6, 43, 15,202,  5,252,251, 12,233,  2,250,255,  2,  1,251,252,  0,235,254,  4,255,  4,  2,  8,  2,  1,  1, 13,  0,229,  1,241,  2, 16,255,255,  8,  1,  0,253,138,  4,  0,  0,  0,  2,  0,  0,  1,255,254,  2,  1,  3,  0,  1,  0,  2,255,  0,  0,  1,  0,  0,  0,253,255,  0,255, 72,252,227, 15,  5,210,232,211,  0, 20,  7,  5, 29, 98,255,  5,  5,233, 33,240, 16, 67, 39,126,  8, 17,252, 15, 22,224, 12,251,219,224, 55, 48,253,156,  4,254,  5,  5,239, 14,  7,251,  1,  1,242, 19, 58, 79, 62,  7,203,211, 37, 15, 10,239, 39,113, 92,254,241,253,193,  7,231,228,238,254,246, 66, 18,  9,247,  8,  0,  0, 12,244,233,  5,249,217,246,247,254, 43,224,180,208,218, 69,204,
      211,  4, 19,255, 12,253,  4,252,226, 26,246,  1,232,234,251,246,190, 43,243,  1,211, 55, 12,242, 17,231,255,250, 66, 12, 14,254,  3,  2,249, 10,  2,237,  0,251,253,  0, 11,253,  5,  3,  2,  3,254,250,250,255,245,  5,  2,  2, 48,253, 11,255, 20,  9,250,254,  1,  0,  0,  0,  1,255,  0,255,  0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,254,255,  0,  1,254,  0,  2,  0,223,247, 17,249,225, 97, 44, 39,252,242,  1, 69,248,214, 65, 31,230,240, 15,214,238, 34, 10,223, 43, 15,220, 13, 28,225, 58,238,  4, 53, 21,219, 20,  1,  4,254,  5, 50,216,  1, 15,242, 63,  6, 21,  2,187, 87,  5,237,105, 20,218,241,167,  7, 51,176,176, 64,230,237, 13,253, 15,246,246,  0,253,229,219,255, 38, 17,  0,253,  1,246,  1, 19,253,  5, 53,  9, 16, 31, 10, 16, 43,219,252,  0,
    },
    {
      235, 14,201,255,239, 63,  8, 13,221,215,234,255,213,232,222,255,221,250,228, 63,244,  7,245, 52,228,219,254,232, 44, 31, 47, 64,  3,255,  2,  3,255,  5,246,246,  6,  3,245,255,  6,  1,  7,  0, 15,  6,  0,  0, 22,  4,  0,  0, 13,244,  1,254, 21, 40,250,252,  1,  0,  0,255,  4,255,254,255,  3,  1,  0,  0,  0,  1,  5,  0,  9,  1,  1,  3,  5,  4,  1,  1,  4,  2,  0,  0,  2,254,  0,  2, 24, 14,  2,  3,255,247,  3,194, 41,222,255,234,211,  0,  9,170,  1,252,201, 51, 24,170,180,232,216,233,203, 24,251,250,  1, 33,211,217,206,  1,241,231, 17,  8,239,227,  3,254, 20,204,  3,196,255,  8,193,248,244,235,  1,  4, 19, 13, 33, 30,  0,246,  2,244,  3,248,248,  0,255,  4,  3,243, 14, 20,249, 17, 29,229, 11,252,255,  5, 11,  4,  1,254,  2,223,  6,  1, 20,233,255,  1,  0,255,
       61,  0, 75, 75, 42,  2, 41,172,230, 31, 11,  2,239, 21, 10,253, 18,231, 17,197, 56,  3, 13, 92, 20, 12,253,200,  1, 20,209,238,253,  1,  8,240,253,  8,254,  1,232,  3,  9,255,246,255,  5,242,  3,  5,  7,254,253, 10,  1,  0, 51,  4,  2,  8,249,223,  3,  8,252,  0,255,  2,253,255,255,  0,  1,  0,  0,  0,  3,  0,250,255,249,254,254,  0,  2,254,  1,255,  3,  1,  2,  2,250,  0,  0,  2, 49, 12,246, 44, 55, 50, 31,207,216,248, 29, 31, 31,250,  7,254, 15,243,203,231,247, 64, 49,203,  1,225,  3,249, 10,165,254,243,244, 16,249,218,242,223,238, 57,242,  6,  3, 18,  0, 35,230,207,253,234,232,254,239, 63,227,  8,253,248, 41,203,  8,205,255, 52,240, 41, 13,  8,234,  0,243, 35,247,239,253,255,  9,255,229,  0,  5,254,  1,253,255,170,  2,  7,  1,255,  7, 39,245,  4,  0,251,
      231, 31, 18,150,233,235,240,109, 62,235, 59,198, 50,  7,243, 48, 34,  9,  4,253,235,  2,  0,123, 51,222, 14, 94,233,  1,244,201,255,250,250,  0,252,  0,  2,211,248,  8,251,  0,252,228,248,  6,246,238,  4,245,  1, 12,253,  3,240, 10,241,255, 19,254,254,  1,  2,  2,  0,255,  4,  1,  1,241,  0,  6,  0,  0,  1,253,255,  3,  2,  0,255,255,250,255,  2,254,  3,  1,  4,  2,  3,252,  1,  0,  0, 17, 31, 15,248,253,  0,  0,232,  8,209, 22,255,  7, 11, 42,215,  7,253,197,204,197, 55,231,221,230,236,193,234,254,252,252, 21,223, 16, 13,  9,236, 27,248, 22,251,254,247,207, 11,227,  1,232,233, 44,188,197,255,198, 56,193, 21,191,229,174,249,146, 54,244,  6,244,237, 10,252, 11,243,  4,  9,251,  0,245,  1,250,  0,253, 68, 15, 50, 41, 29,  7,255,250,224,  0,253,252,234,  1,253,
       23, 43,201, 76,246,236,251,  6,  2,  4,241, 83,  4, 13, 15,247,255, 11,  8,252,245,218,  0,252, 12,  7,212,  2, 14,234, 16, 15,  1,254,  0,249,  6,  4,  4, 62,  8,252, 11,253, 12, 32,  9,251,  4,  1, 14,  4,248,247,  1,  3,244,249,  0,  2,253,  7,  3,  0,  0,  1,  0,255,  0,  6,  0,  0,  0,  8,  3,255,  1,  3,  3,254,  0,255,  2,  0,  2,255,255,  0,255,253,  0,  0,  0,  4,  0,255,249,237,220, 10,241,  0, 70, 65, 60,224, 48,206,220,250,206,242,255,212,241,209,248,248, 16,191, 19,  4,245,  1,254,211,253,124, 52,232,253, 37, 40,  3,  6,241, 25,  7,  5,254, 43,251, 60, 15, 13, 10, 28, 36,107,205, 39, 59, 57,212, 64, 13, 25, 68,126,  2,237,233,  2,233, 14,221, 47,  1,  2,250,  0, 13,  5,  1,  2,241,241,204, 12,  6, 37,251,  1,  0,251, 10,246, 52,230,254,  5,250,
    },
  },
  {  //  14
    {
      211, 73, 26, 15,220, 74, 25, 10,228,206, 18,251,251,205, 13, 25,230,203, 14,  7,183,213,195, 14,  3, 11,229,240,112,154,254, 69,247,  2, 55,250,247, 13, 54,252,  8,237, 58,255,  6,230,186,254,  8,246, 62,254,  0,  0, 21,  2, 12, 29,244,254,247,235,249, 29,255,  0,236,  0,255,  0,239,  0,255,  1, 27,  0,  0,  2,246,  1,255,  0, 26,  0,  3,  2,255,  0,  1,  2,  1,  2,254,254,  0,  1, 26, 17, 24, 36, 27, 19,  5,  7,226, 54, 30,202,197,201, 18,  6,  8,240,  9,242,193,255,193, 17,254, 28, 70,219,250,247, 83, 43,254,  1,247,  0,  1,  1,  6, 21,230,  6,235, 19,242,240,223,249,236,195, 53,197, 80, 78,251,223,253,250, 24, 66,198,  1,254,170,  6,  0,  0,  0,  1,  0,  0,  0,  2,  1,  1,  0,  3,  4,250,  5,213,228,  9,230,254, 49,  0,253,255,  0,253,  2,236,  2,  0,251,
       31, 10,203,  3, 20,  8,206,  2,249, 13,199,251,221, 98,  9,251,  3,  1, 38, 32, 28, 31, 71,  0,241, 69, 63,168,192, 26,199,182,255,253, 23, 10,252,240, 26,255, 11,228, 59,  1, 14,  6,200,244,255,227, 63, 52,  0,252, 32, 56,  8,254,253,244,245, 41, 18,251,  1,  0,245,  1,  2,  0,245,  0,255,255,  0,  0,253,  0,255,  1,  0,  0,  8,254,255,  0, 35,  3,  1,  0,  0,  0,  0,255,  1,  0,254,237,234,228,253,232,241,  8, 59,217, 43,200, 47, 62,207, 20, 69,228, 59,245, 64,184,103, 46,  0,255,210,182,  1,209,252,119,255,255,  4,  2,  2,251,  4,246,  1,  4, 99,  4,243,244,252,255,252, 27, 54,  9, 37,247,132, 13, 64, 54,247,142, 63,252,  4,196,  0,  0,254,  0,  1,  1,255,  0,  1,  0,251,  0,  1,252,255,255,  8,  2,  5,  4,  3,  2,255,  0, 41,  3, 15,255,255,  3,  0,  3,
      223,251, 63,241, 38,255,  3,  1, 36, 37, 13,  0, 36, 13, 63,240,250,  6,  0, 10, 34, 92,243,252,  5,197,175, 20,246, 14, 61,  4,235, 66,210,239,  0, 28,237,  5,239,242,176,251,252, 24, 82, 10, 25,  1,247,244,  4, 23,193,220,  1,  2,200,  4,  3,  3,198,  7,  2,  1,255,  0,255,  0,  8,  0,254,  3,  4,254,  1,  0,255,  1,255,  0,  0,251, 12,  0,254,  1,  0,255,  2,255,  0,  0,  2,255,250,253,255,  6, 10,  1,  5,  1,244,  2,255, 41, 10,243,252,189,  4, 55,199,228,210,219,  1,185,211,239,135,246,  6,191, 24, 63, 14, 12, 13,  0,250,  4,255,250,251,217,234, 67,255,  0,253,  0, 49, 11,137, 14,193,224, 41,223,  1,170, 35,199, 63,  0,  6,  3,  5,  0,  0,  0,  0,  0,  0,  3,255,254,  0,243,  0,  2,255,253,  2,  1,  2,  2,247,  1,250,255,217,241,250,  1, 41, 23,252,  0,
       24,191, 30,252,  0,236, 69,254,222, 10, 19,  8,243,205, 63,198, 53,207,235,240, 88, 36,230, 49,212,254,  0,131,225,254, 34,  8, 65,237,239, 24,255,  3,195,  1, 24, 11,215,239,248,  2,249,250,254,254,205,197,240,  9, 61,207,  5,252, 33,223,245,255,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,254,254,  0,  0,255,  0,  1,255,255,255,  0,255,  0,  0,  0,255,  0,  6,255,  0,  1,  3,  1,  5,  5,  0,254, 22,241,  3,247,  7,248, 54,248,249, 16,248, 40, 12,250, 60,235, 75,244,249,247,224,251,252, 79, 39, 54, 18, 66,  6,  1,  1, 11,  3, 10, 15,253, 12, 44,228,184, 14, 18,  4, 22,203, 20, 62,198, 61,197,195,215,195, 68,188,244,222, 58,196,212,  0,255,  0,  0,  1,255,255,  0,254, 10,246, 10,254,  0,  0,  8,  0,  1,255,254,  0,  1,255,  1,  2,  1,  2, 24,248,  3,  2,246,
    },
    {
      238, 18,246,196,  1, 47,250,254,208,250,249, 59,225,  3,240,226,235,226,245,242,222,242,  7,246, 36,242, 14, 12,246,249,237,176,252,  1,246, 58, 13,  6,243,231, 11,250,  4, 37, 12,253,250,254, 10,  0,242,  1,  7,  1,243,  2, 13,250,252,252, 16,  2,  0,254,254,  1,255,255,  4,  0,254,  0,  4,255,  1,255,  3,  0,  0,255,  3,  0,  1,255,  5,254,  0,  1,  2,  1,255,  0,  3,  1,  0,254,252,  2,249,  5,  6,  2,246, 17,229,243, 26, 30, 16,236,  2,253, 12, 37, 21,250,155,239, 14, 34, 37, 26, 49,247, 15, 70,193,  4,254,  1, 12,  0,255,248,242,  2,248,245,225,243,  9,252,  3,  2,215, 41,  0,202, 50,198, 61,  0,255,250, 31, 20, 37,244,207,216,  2,254,  2,255,255,  0,  2,255,254,253,  1,  1,  1,  9,  1,  2,240,255,  0,242,225,247,239,  0,  0,253,  0,  0,253,  0,244,251,
       48, 22,190,246, 24, 24,255,  5, 35,214,219,217, 34,224, 44,232,  4,  4,  7, 34, 20, 20,222,254, 51,187, 54,  3, 23, 13,252,234,  0,251,  1,  3,238,  5, 11, 28,253, 13, 15,235,254,  7,  5,252,  1,254,  4,243,227,253,  4,248,  1,  0,241,  0, 28,  1, 14,252,  0,253,  2,  1,251,255,  1,  2,255,  0,  3,253,  1,255,  0,  1,  0,  0,  1,255,  1,  0,  2,  0,  2,254,254,  0,252,  0,254,  3,  2, 18,  5,241,250,253, 10,241, 33,206,253,249,249,240,  4,  2, 35,218,240,218, 96,249,247, 28, 10,  7, 20,234,  4,  2,  1,187, 15,250,  0,253,251,255, 35, 11,233, 24, 29, 12,231, 23,252, 31, 63, 31,226, 69, 83,255,  7, 52,255,  9, 22, 40, 47, 15, 93,  4,253,253,  0,  0,  4,  0,255,255,255,  8,  2,  0,  2,203,  2,254,  3,240, 12,  1, 15,242,246,  2,  0,  1,255,  0,  1,  1,  8,253,
      201,228, 16, 90,228,200,  4,168,252, 46, 11,218,  2, 25,243, 63,221,254,  6, 55, 51,253,248,255, 12,250, 49,  2,203,247,  9,255, 16,249,251, 57,251,  6,  1,  8,250,254,206,  6,243,248,  1,255,209,  0,251, 60,221,  0,246,255,  2,  1, 34,255,  4,236, 10,255,  1,  0,253,  5,  0,  0,255,  1,  0,  1,253,  3,254,  1,255,  0,  0,  0,  0,250,  6,  2,  1,252,  0,  0,  5,  0,  1,  1,254,  0,  2,  6,229,  6,  7,248,  0,254, 13, 52, 15,200,  5, 20,  0,  6, 14,242, 23,184,  5,  0, 41, 79,251,  2,  2, 52,  1,230,254,  4,  2,  9, 24,253,  8,  2, 12,  2, 37, 30,  7, 44, 43,  2, 38, 52,245,246,251,213, 14, 64,255,110,228,193,248,245, 48,  3,221,  0,  2,  3,239,  0,251,  3,  1,  0,  1,252,254,255,  1, 48,  3,  6,255, 16,228,252,251,  0, 12,253,222, 34, 28, 13, 13,254,  3,  1,
       33,250,  8,193,  4,  3,244,244, 32, 17, 36,246, 27,  1,241, 99, 31,254,253, 49, 23,  0,227,  6, 46,  2,255,  1,238,  2,  2,238,250, 11, 46,192,  3,  1,255,  0,  1,  0, 21, 15, 12,  3,  1,247, 45,254,  0,253, 66,  8,  7,232,  2,254,228, 30,247,  2,252,  0,  1,255,  0,255,  3,  2,255,  0,  0,  1,  1,255,  0,  0,255,255,254,  2,  0,  1,254,255,  2,252,  1,  1,255,  0,  0,  0,  0,  0, 13,  5, 49, 19,  6, 27,  0,253,  0, 12, 26, 17,  2,  9, 56,249,  8, 13, 60, 99,228, 47,  9,  6, 16, 60,239,212,246,190,  0,227,249, 10,  1,  6, 25,248,240, 13,  1, 14, 34,  6,252,  0,197,172,188,218,  1,244,200, 13,206, 11, 19, 63, 56,166,215,254,  0,252,253,  3,  1,255,  4,253,  0,  0,  3,248,237,251,249,  0,  4,  8,  1,249, 29,  8,254,  1,  2,255, 29,  6, 20,  6,248,  0,  2,  6,
    },
    {
       17, 29, 12, 17, 14, 16,142, 36, 23, 20, 13,  9,  0, 21,  8, 11, 12, 14,  6,  5,251, 12,249, 10,  6, 11,227, 19,229,  0,216, 24,  1,  9, 57,217,  3,  4,198, 36,254,  8,204,247, 10,  8,239,251,  5, 10,246,249, 12, 10,252,252,247,  5,  3,  3,233,  5,255,  4,  0,  0,247,  0,  0,255,249,  1,  1,  0,  3,  0,  2,  0,  0,  0,  3,  1,  0,255,  1,  1,255,  0,  1,  2,  3,255,  1,  2,  3,  0,234, 16,231, 39,230,168, 23, 36, 79, 19, 17, 13, 70, 47, 24,217,251, 19, 47, 17, 17, 66,209, 44,215, 26, 17, 83, 83, 93, 16,167,248,248,245,231,247,254,206,189,  9,  8,  4, 16,245,238,161,223,255,196,226, 67,193,206,239, 25,251,198, 56,234,207, 33,  2,241, 13, 36, 21,  0,  0,  3,  1, 23,246,251,254,255, 29,  2, 35,  4,  0, 15,252,223, 37,  2,  4,  1,255,247,226, 11,249,248,  0,255,
        1, 16,216,255,229,  5,133,235, 40,251,  1,  6,  5,  1,233,231,245,255, 48,  4,253, 37, 54,254,  2,  7, 22,  2,230,255, 19,243,254,  3,201,  0,  9,  0,250, 15,243,  2,196,  3, 12,  0,230, 14,  5,  0,130,  1,  1,  0,197,  2, 53,  4, 29,250, 34,147,254,255,  0,  1,  0,  1,  0,  1,  0,  0,253,  1,251,252,  0,  2,  1,  1,  0,  0,  1,  0,  0,  0,247,  0,250,  1, 14,  0,251,  2,  0,  1,235,  3, 39,246,209, 18, 19,145,251, 31,  3, 21, 13,225, 51, 36, 20,  0,197, 60,196,  2,  2, 24, 37,242, 53,207, 29,236, 67,214,223,202,250,  6,  2,239,164,  2,  6,251,254,245,224, 10, 72,177,240, 22,248,212,  7,191, 74,202,  2, 29,253,194, 38,135,246, 12,191,218,255,251, 34,252,  0,  2, 13,  1, 13, 29,  5,235,  5,255,254, 28, 17,254,  7,  0,255,  0,  4,246,  5,255,255,210,  0,255,
       22,  0, 57,254, 26, 57,214,242,248, 59,  9,226, 22, 12, 30, 11,253,  4, 26,229, 20,  5,206, 12,174,243,  1, 31, 13,251,243,255,226,  1, 66,251,  0, 16,254,252,  1, 24, 74,  0,238,  1,232,  4,  0,  1,255,255,214,247,222,  1, 33,  5,  6, 64,246,128,251,254,  0,  0,  1,  1,  0,255,235,255,255,255, 47,255,  0,  0,255,  0,  0,  0,  1,  0,255,254,255,  0,252,  0,  1,  0,  1,  0,240,  2,244, 30,  3, 17,  3,148,238,193,254, 26,  1,  4, 52,233,  0, 28,250,  1, 54, 77,198,227,  1,252,254,  4,  4,254,220, 62,243, 51, 20, 31,  9,216,239, 60,252,114,255,  7,  4,  6, 36,  5,253,241,249,219,211,  9,199, 10,194,249,253, 95,210,  1,242,248,249,  0, 36,248,253,250,240,248,255,  1,253, 53,  2,  0,246,  3,248, 19,  2,227,  9,  2,184, 33,  0,  2,253,  1,255,255,  0,255,252,254,
      229,250, 18,  0,  0,193,218,  3, 16,215, 81,  1,218, 36, 15, 17, 30,204, 23,  4,254, 37, 51,249, 66,  5, 33, 16, 32,  7,255,255, 14, 12,194,218,252,253,107,  4,  0,254,191, 23,  7,254, 12,249,255,  6,237,  1,  7,  2, 66,255,  8,242,  0,  0,254,  0, 23,241,  1,255,  0,  0,  0,  1,255,255,255,  1,  0,  0,  0,  0,  2,  2,  0,  1,  0,  1,  3,  1,254,  0,  1,  1,  0,  1,255,254,  5,255, 12,252, 10,250,239, 52, 62,  1,223,224,  0,  1, 61, 81, 23,216, 16, 33,253, 67, 57,254,156, 13,  4,  4,  6,251,252,189,196,209,248,  3,  5,242,247,244,  0,245,  4, 13,  6, 69,168,  1, 33, 45,  6,240, 60, 33, 87,189, 49, 77,  0,254, 57,255, 38,221,251, 13, 13,  9,247,  0,239,126,  0,  0,  0,203,254, 13,246,253,249, 12, 10,248,  2,  2, 20,  5,252,  1,255,  3,240,  0,255, 18,  0,251,
    },
  },
  {  //  15
    {
      239, 47,240,196,237, 48,238, 86,  4,215,246,202,254,220,240,243,235,214,242,195,216,205,208, 71,247,232,226, 90, 47,223,199, 56,255, 31,  4,246,  1, 36,  3,255,248,250,254,246,251,252,254,254,246,252,247,228,  1, 25,248,  2,250,243,179,  0,248,  0, 80,251,  2,252,  0,254,  2,253,255,  1,  1,  0,  0,  1,  2,255,  1,  0,  3,  2,  1,  1,  1,  3,  1,255,  0,  2,  0,  0,255,254,255,  0,252,221,203,217,215,219,  0,179,247,240,213,223,215,222,236,236,  4,243,  4,234,234,196,241,184, 20,231,  7,218, 84,194, 14,170,252,250,247,248,  0,253,212,  2,252,  2,248,255,  6,223,227,  0,192,198,198,196,194,202,216,171,197,193,201,195,202,208,169,195,  4,  1,  1,  0,  0,  1,  3,  1,254,255,  1,  0,  0,  1,  0,  0,249,252,250,  4,  1,246,  4,  1,  0,251,  2,254,251,252,  4,246,
        9,251, 29,150, 10,196,  2, 67,218,225,209, 15,225,245,  8,244,255, 15,124,  6, 20,  9, 37,125,  6,233,219,184,233,252,188,198,  1,236,190,191,  3,175,196,140,  4,247,190,191,  0,234,197,192,237,255, 64,220,  4,  3, 23,207,  7,235,190,201,233,201,214,207,199,199,203,195,199,190,185,193,193,192,197,202,193,192,198,216,193,194,192,191,190,193,194,193,192,191,193,201,191,192,192,194,229,253,249,248,234,  3,253,250,153,215,  1,255,216, 15, 15,254,227,238,198,128,168,226,196,238,  0,189, 49,  4,  6,239,249,247,  1,  0,  4,  0,253,  0, 10,  6,254,251,250,  1,251,255,247,  0,252,253,  2,255,  5,242,255,255,209,251,  0,248,240, 10,246,  1,  0,  0,  0,  0,  2,254,  1,  0,  0,255,  0,  0,254,255,255,  0,189,243,246,241,235,252,243,242, 12,255,  4,  0,  2,250,  1,253,
      241,237,243,199,253,254, 12,246,251, 62,  4, 45,184,252,251,254,192, 51,157,251,247, 52,233,231,  0, 69,242, 42, 35, 93, 18, 44,255,246,  0,244,  5, 48,250, 60,  0,250,  2,243,253,  1,  4,255,224, 94,253,247,246, 50, 14,  3,253,252,251,248,  9,  4, 32,  9, 57, 57, 50, 61, 57, 63, 62, 67, 64, 64, 45, 60, 63, 64, 53, 41, 64, 64, 58, 62, 63, 64, 65, 64, 64, 65, 46, 61, 64, 64, 53, 61,249,252,255,  5, 21,226,  2,  9, 38,220, 11,245,129, 18,  8,  3, 44, 16, 28, 16,127,139,250, 47,230,230,238,229,202,  8,251,234,212,212,210,207,188,200,215,212,196,205,199,197,141,182,206,240, 12,  2,  8, 10, 17,  8,  4,  3,  0, 12,  2,  8, 64,197,254,250,194,196,194,194,162,190,193,193,191,196,195,193,  6, 56,193,190,233,233,228,231, 24,235,238,244,229,229,223,221, 13,255,233,238,
      231,220,123,198,235,237,194,227,226,213,238,212,216,216,190,206,225,217, 13,209,234,249,234,202,223,209, 14,205,218,217,  2,204,201,209, 22,255,193,213,  2,  0,201,208,  7,  3,200,216,  1,252,202,212,  1,  2,201,219,  7,250,201,211,246,  4,203,224,  1,  0,192,192,192,193,194,191,192,193,191,192,194,197,193,193,194,197,192,192,194,195,195,192,193,193,192,192,193,196,193,192,193,196,250,  0, 38,249,221,242,245,232, 59,252,229,238,234,243,226,233,  2,171,228,230,220,248,219,223,206,  1,233,238,238,245,226,236,253,255,  0,  0,  2,  1,  1,  0,250, 75,255,  0,  3,  3,255,  0, 57,251, 17,  8, 26, 20, 22,  7, 64, 62, 21, 22, 24, 29, 20,  9,250,251,254,  1,  1,255,252,254,103,192,  0,254,  0,255,  2,254,180,215,204,215,220,221,216,219,252,195,209,213,225,224,218,222,
    },
    {
       19,  0, 37,224,239,  7, 14,245,243,248,249,252,249,244,246,  2,  8, 19, 10, 12,238,  7,240,254,243,  8,  4,228,  2,247,  5, 21,  7,254, 15, 17, 20,  3,252,252, 10,253,243,251,  8,  2,245,251,  0,  6,247,  6,  7,  1,237,255,  8,  6,249,240,  8,  5,254,253,  0,  1,255,253,  2,  0,  0,255,  2,  0,  0,  0,  1,  1,  0,  0,  2,  0,255,254,  2,  1,255,255,  2,  1,  0,  1,  1,  1,  1,255, 10, 11,  4, 16, 40, 10, 11,254, 57,  4, 16, 24,231,232,220, 74, 12,251, 82,218,190,212,196, 16,221,180, 20, 46,  0, 24,228, 56, 12,  9,248, 27,244,  0,253,  0,  4, 30, 29, 32, 39, 16,222, 29, 42,  4,210, 50, 23,234, 40, 62,227,  3,  0,  6,  2,  9,247,254,  6,253,  2,  0,254,  1,  2,  0, 30,230, 30, 41, 36,249,  0,  5, 20,  9,  1,253,252,  3,  1,224,  4,  6,253,252,  4,  0,  2,  4,
      241, 31, 25, 31, 40,  6,251,138, 11, 42, 18,204, 24, 14, 51,250,  4, 56, 15,  1,245, 16, 23,  6,242,252, 42,249, 21,191,251,219,252,  0,200, 45,226,  1, 91,194, 14, 10, 67,220,255,254, 58,194,254,250,201,  8,244,254, 35,254,254,255, 20, 30,249,  4, 67,241,  0,  0,255, 12,254,  0,254,230,  0,  0,  0,255,  0,  0,254,247,  0,  0, 12,  0,  0,  0,255,  0,  0,  0,255,  3,255,  0,  1,  1,  9, 29,  5,  1,247, 53, 29,246, 41,242,223,  1, 60,  0, 30,  1,238,249,  2, 45,188, 60,255,206,244, 46,  2, 30, 16,229,  2,199,217, 14,  7, 17,254,242,219,254,229,  0,214,  8,252,250, 17, 17,213,234,249, 59,235,202,246, 59,253,228,  2,192,226,220, 40, 58,  2,  0,  0,254,253,  0, 12,  0,217, 32,202,253,254,  1,248,255,250,249,255,193,  4, 34,  0,  4,250,244,254,250,250,  2,241,255,
       18, 39,  6,246,247, 62,250,213,  4,  4, 50,251,231,253, 61,  3, 53,187, 11,253, 15,254,251,  1, 13,185,205, 49,202, 65,215,  0,  3, 20,  2,193,  1,  6,237, 64,  2,  2,201,159, 12,  0,195, 60,254,255, 65,253,  3,  0, 65, 36,249, 15,251,  6,252,254,195,164,  0,  5,  6,  0,  1,  1,224,  0,  0,255,  1,255,  0,  3,255,255,255,  0,251,252,  0,  0,  0,  2,  0,  1,  0,255,  0,  0,255,254,249,232,235, 28, 31,250, 11, 44,250, 16, 29, 36, 14, 70,199,186, 51, 23, 60,215, 73,168,129, 15, 32, 34, 31,253,206,  0,  2, 31,230,191,199,250,207,  6,255,  0, 46,218,  2,253,  2,  8,249,193,226, 12,244,  9,216, 47,254, 20, 35,  1,189,229, 55, 51,  5,174,201, 26, 46,  2,246,  0,  5,255,  6, 22,238,254,172,  0,  1,254,  2,207,249,247,251,249,  6,  0, 50,216,  3,194, 33,254,253,  6,
       35,218,212,178, 11, 16,209,  0,  2,212,211,156, 11,  2,  2,234, 34,253,162,  0,223,254,255,248,252,131, 26,193,248, 21,  0, 28,245, 15,  8,  7,254,  4,225,239,  0,  2,253,  3,229,254,  5,156,  2,  0,220,247,  0,  0,  1, 38,  0,  0,  3,  0,255,253,255,209,254,  1,255,255,255,  0,  1,  0,  0,  0,  0,  1,255,255,  0,238,  0,  0,  0,  2,255,255,  0,  3,  0,  1,  0,255,  0,  0,  0,  0, 40, 23, 11,246,  6, 10, 11,  7,251, 26, 10,  0,254,250, 67,202,235,250,211, 30,251, 48,241, 13,243, 17,253,224, 64,255, 48,226, 31, 18,254,242,255,  5,243,  3,239,248,243,254,245,  3,  0,248,  1,224,215,235,223, 63,255,217,234, 34,  0,249,  6,  1,241, 58, 11,218,  4,253,  9,255,255,  0,229,252, 39,  0,248, 10,  0,255,246, 47,204,  0,  1,  0,  0,  5,212, 43,253, 63, 12,251,254,254,
    },
    {
      216, 18,240,248,211, 18,210,247,219, 13, 26,241,196,  3, 59,  6,200,253,  4,  2,194, 10,235,255,248,  6,255,183,  3,  2,254,144, 12,253,240,249, 10,255,227,  2,218,253,241,  3,  3,255,243,254,  1,254,255,250,247,255,254,  2,  5,  4,  3,236,  9,  4, 10,253,  3,  0,  1,  0,  1,  0,255,255,  1,255,  0,  0,  1,254,255,254,  0,  0,  0,255,  1,  1,  0,254,  1,  1,  0,  1,  2,  1,255,  3, 13, 43,141,186, 34,213,231,183, 23,236,253, 10,  7,249, 54,230,245, 56,181,225,170, 48, 56, 34,243,216,239, 29, 39,243,  0,244, 16, 27,252, 15,221,230,  5,227, 21, 64,245, 35,232,246,232, 17, 44, 29,232,245,  6, 17,254,  4, 23, 21,193,225,196,255,248,  4,254,253,  0, 10,  1,253,255,254,252, 37,  0,245,  1,255,250,253,253,252,255,241,255,  0,  0,250,  4,  4,223,246,229,255,  0,253,
       64,251,249, 11,233,249, 27, 15,236,  1,232,240,251,252,175,238,221,241, 61,  5,243,  2,  6,253,218,  0,194, 21, 49,  3,195,248,  0,254,194, 59,  1,  0, 11,242, 24,255,  0,252,  5,  1,141,  4, 15,252,199,  0,  2,255,221,  1,  2,  0, 26, 28,251,  3,247,255,  0,  0,254, 11,  2,  0,254,254,  2,  0,  0,253,  0,  0,  0,  1,255,  1,  1,255,  0,  0,255,  0,255,255,254,  1,  0,255,  0,255, 12,255, 87, 72, 44, 51, 13, 14,228,251,  2,242,227,226, 69, 54,227,177,  3,  4,186,205,144,207,  6, 21, 71,231,  7,254, 19, 14,  5, 36,100, 25, 63,238, 61, 65,245,198,247, 37, 34,254,137,116,252, 38, 39, 25,  2,  2,254, 17,246, 18, 74,232,252,247, 56,208,  0,  6,245, 13,241, 14,247,248, 10,  0,254,248,  0, 13,255,255,  6,  3, 13, 15,255,255,255,  0, 26,  8,  2,246,254, 14, 25,  4,
      246,251, 22,245,  6, 64,  4,239, 34,  7,235, 16, 44,252,230, 11, 50,250,237,255, 21,  3, 66,217, 50,254, 50, 48, 36,233, 62,  2,  3,  2,237,251,  2, 10, 46, 20,227,  5, 62,218,  6,  2,156,243,243,  2,229,243,252,  3,  3, 25, 22,  0,215,234,  0,252,  6,244,  0,  0,255,  1,  0,  0,250,  0,  0,  1,251, 22,  1,  0,  0,  2,  1,  1,255,  0,255,  0,  0,  1,  1,  1,255,  1,  3,  0,255,  0,250,241,255,251,  2,214,252,202,  2,  0,227,251,  9,248,255,251, 37,223,253,242,221,194,255,192,234,216,  1,103,234, 32, 13, 42,  4,211,197,252,216, 11,194,180, 15,252,251,250, 50,230,175,238, 24,229,240, 32,236, 18,244,165,  2,235,  0, 39, 15,125,191,  6,  1,223,217,248,211,  0,  1,  1,  2,  1,  0,  1,117, 36, 54,227, 19,  0,252,255,  4,  1,  6,  0,233,254,  0,252,255,  2,  1,253,
      255,252,251,199, 56,191,251,253, 37,  4, 12, 14, 19,  9, 10,  3,  3,247,  0,  2, 59,251, 10,  3,246,255,255,218,246,182,  2, 50,255,  4, 63,196, 11, 17,  0, 16, 11,  3, 31,249,  3,  1, 53, 11,253,253, 17, 13,  4,255, 51,240,232,  2,  3, 48, 70,  8,  0,  9,  0,255,  3,  2,  0,  0,  0,  1,  1,  1,  0,  6,255,  1,  2,254,255,  0,  1,  0,  0,255,254,255,  0,  0,  1,  1,252,  0,255,  0,196,209, 27, 10,251,252,  7,239,241,249,197,246, 18,248,225,241,238, 14,  6, 11,  9,190,109, 65,255,  7,212,146,216, 34,  1, 66,  5,251,232,251,254,231,251, 32, 11, 20, 58,  3, 22, 30, 19,223,198,  8,  2, 10,229,255,  1,251,252,255, 31, 71,255, 31,254, 24,  1, 22, 90, 41,252,  2,230,243,  5,231,  1,  2,148,255,209,  1,254,  0,  1,255, 39,  2,254,255,255,251,254,  9,  1,130,  1,251,
    },
  },
  {  //  16
    {
      210,216,  4,229,212,216,  5,217,211,217,  6,201,212,216,  4,196,212,218,255,204,217,218,253,225,207,210,239,221,200,196,212,200,252,  1,194, 50,248,  4,195, 48,253,  5,191, 46,252,  6,192, 55,245,  5,192, 66,245,  8,193, 60,250,  3,194, 68,253,  5,196, 73,  0,  1,227,255,  1,255,227,255,  1,255,233,  0,  1,255,234,  1,  1,255,235,  3,  0,255,231,  5,  0,255,240,  5,  0,255,246,  6,228,204,228,189,191,234,228,163,232,221,226,237,223,214,226,232,234,206,226,218,228,193,226,239, 27,207, 23,210, 22, 35, 15, 42,215,206,215,211,194,195,208,200,245,  8,  5, 35, 13,230,251, 43,217,172,214,201,216,217,209,212,  4,247,  3,248,  3,  0,  1, 12,219,192,225,193,207,185,214,236,239,251,247, 36,248,251,242,  9,242,251,240,250,254,253,244,248,192,200,191,199,191,192,192,185,
        1,  4, 19,  5,  1,239,201,254,  0,  0,248, 13,252,  1,155, 13,254,  2,153,  1,  0,250,156,  2,253,241,156,  0,  1,234,185,  0,204,  0,  1,  0,205,  0,252,  3,207,243,255,239,199,252,254,  6,205,249,254, 13,210,251,  0, 21,204,247,  1,  8,197,250,253, 16,192,191,199,200,205,192,209,192,195,192,202,191,193,192,194,194,194,191,201,192,193,191,206,192,192,192,202,192,203,192,200,192,  0,255,251,  9,250,255,246,  1,193,192,191,191,190,194,192,196,227,232,236,236,244,229,235,235,255,250,247,224,244,252,245,248,227, 31,  1,249, 65, 20,214,  3,  2,  8, 11, 19,  9, 47, 23,  4, 93, 47, 18, 70, 52,208,223,230,188, 31, 10, 69,201, 21,  8, 10,250, 13,232,243,  2, 16,224,  2,213,239,241,209,213,243,248,226, 52,233, 59, 66,244, 25, 31, 11,179,196,252,203, 24,206,235,222,
      246,255,  1,193,255,254,255,216, 96,203,  1,194,223,224,  0,206, 65,252,255,207,  2,219,254,194, 14,  5,  1,215,  2,227,  0,230,255,191,186,245,247,195,255,230,248,193,  3,205,202,185, 11,214,254,188,251,224,238,193,  1,209, 14,197,237,185,252,207,245,199,254,  0,224,  1,244,  0,251,  0,  1,  1,  1,255,  4,  0,  0,  0,253,  0,231,  0,244,255,  1,  0,  0,  0,255,255,  0,  0,  3,  1,  0, 20,  0,254,  0,  2,  1,  0,  1,253,  0,255,254,249,  1,245,  2,  0,  0,249,252,  2,255,  1,  0,  1,255,  2,  6,255,255,254, 99,195,185,174,194,250,  2,253, 10,204,  2,182,215,216,  9,208,254, 13,193,145,216,222,255, 53,247,234,239,251, 37,231, 49,188,  2,166,  4,  0,217,245,254,  6,232,242,215,189,243,219,216,  0,207,242,  3,227, 57, 11, 47,239,  4,254,  3, 10,221,152,218,180,
        0,185,  0,255, 19,183,  2,254,  1,192,255,  6,191,215,253,252, 12,231,  0,255,228,194,  0,254,  4,192,254,246,252,205,  4,  8,126, 18,200,221,215,233,253,246,193, 45,  8, 14,188, 55,218,253,183, 61,221,245,195, 28, 41,  1,244,254,251, 23,247,231, 33,229,  0,251,254,  0,  2,255,  0,  0,  1,  0,  0,  0,  1,  0,254,255,  0,  0,  2,  0,255,  1,  1,  1,  0,  0,  0,  2,  0,253,255,  1,254,  8,  3,  3,254,  0,253,  4,  3,  5,  1,  5,253,254,  5,  1,  1,  2,  2,  1,253,  1,255,  2,254,  0,254,  1,254,  0,255,252,150, 62,191,243,252,244,  8,242, 62, 10, 90,207,241, 19,247, 33,215,196,189, 11,198, 38, 10,212,158,194,  8,237,238, 18,176,  1,223,201, 11,239,  6, 36, 15,208,237,243,255,248, 12,  6,212,255,  1,190, 13,254,227,207,168,232,181,  7, 49,211,224, 28,219,249,
    },
    {
      221, 12,247, 20,222,  0,247, 11,214,252,241,251,217,251,242,252,237,248,241,245,251,248,241,245,238,244,243,245,241,249,248,255,  5,253,242,225,  7,252,235,234,  1,254,238,209,  8,  0,241,234,  5,  0,239,247,  5,  2,246,244,  7,  5,240,220,  3,  1,251,244,  2,255,255,251,  1,255,  0,254,  3,  0,254,251,  2,  0,255,252,  1,  0,  1,254,  1,  1,255,254,  1,  1,  0,251,  1,  1,255,253, 11, 11, 41, 22, 13, 19, 27, 11,  6,  3,254, 28, 65, 49,217, 24,209,192, 18,195,213,201,  0,183, 30, 21, 83,233,  2, 41, 12, 18, 19, 26,238,  3, 12,  5,246, 12,248, 11,249, 25, 20, 14,244, 40, 32, 15,253,203,211,205,195,217,251,250,159,  9,251,255,195,255,251,249,  7,  0,  1,  2,  1,  3, 26, 35,240, 19, 19, 18,  1,  8,249,253,250,  5, 20, 10,225,  4, 13,  0,252,  2,  3,250,251,255,
       43, 18, 14,  6, 34,242, 27, 13, 46,213,  7,  3, 43,255,  1,  2,243,  1,248,247, 35, 10, 17,  3, 42,236,255,253, 23,  4,253,255,252,228,230, 10,252,  6, 19,  7,  0,254,185,107,254,  3,241, 40,  0,  1, 21, 39,  0,  0,  1,249,  3,  1,104, 27,  4,  2, 39, 34,  0, 11,239,  2,  1,  3,252,  5,  0,  0,  8,255,  0,  1,254,  0,  0,  0,  0,  3,  1,  0,  0,  0,255,  0,243,  1,255,  0,  0,  0,  1, 30,  1,255, 34, 41,  2,251,251,  2,233,245,237,  2, 48,  2,  3,254, 15,  3,  8,235,250, 33,248,  7,  7,248,225,242,  7,209,234,221,251,254,  6,  7,254,249, 32, 69, 48, 60, 90, 64, 47, 58, 16,253, 58, 83, 58, 30, 67,  1,233,251, 74,  1, 33,  6, 22,251,249,252,249,  0,253,  0,253,254,209,254,203,235,229, 25,192,228,  3, 15, 13,  0,  9,250,  2,  1,  6,  4,  3,  0,  2,  0,253,  0,
      253,253,237,  4,245,  3, 12,254,255, 55,255, 14,253, 23,235,230, 28,  0,253,255,240,249,246,255,227, 20,  9,  7,251,239,255,  1,  1, 28, 43,255,  2,254,255,245,  0,  9, 52,179,  2,  1, 52,  4,  1,  2,129,239,  4,  6,246, 17,253,253,177,  9,  0,254,254,118,255,255, 14,  4,  0,255,  0,  1,  0,  1,  0,252,  0,255, 10,248,  0,  1,245,  2,  0,  0,  8,  2,  0,  1,  0,  1,  1,  0,  0,  1,  5,  3,  0,  0,  0,  3,255,254, 43, 33, 29, 30, 50,186,192, 25, 36, 71,  2,120, 43, 27,252,119,235,249,217, 44, 77,  2,224,245,251, 50, 11,  7,  2,215,242,243,222, 21, 25,205, 50,  0,255,190,248, 61,249,242,237, 39,217,170, 38, 24,250,136, 34,  7,  4,190, 14,248,  8,  2,255,255,  0,  0, 35,218, 51, 10,254,220, 65,  0, 18,242,  1,  1,251,  7,  0,255,220,253, 67,205,248,  1,  4,207,
        0,226,  6,225,  1, 22,246,243,  7,246,245,  5,250, 26, 44,248,  1,234,  2,  1,  3,  1,  1,247,  8,241,  6,  0,251,  4,  1,  3,  4, 21,  0, 26,255,  6, 36,250,  4,  2, 50,230,  0,  0,247,  3,  0,  0,  5, 14,254,255,254,  1,  0,  1,204,  0,252,  5, 15,192,  3,251,253,253,  1,255,255,  5,  1,255,254,  5,  0,  0,  2,  0,  0,255,  2,  3,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,255,255,253, 14,  3, 11,  0,251, 11,  1,234,  7, 13,  2,  6,118,104,253, 33, 23, 40, 46, 14, 60, 10, 82, 25,248,225, 12,162,251, 19,216,  0,221, 16,254, 12,248,  1,  6, 42,205, 51,  1,112,194, 29,  2,249,203,231,  4, 11, 17,233,109,  8,242, 50, 92, 58,244, 41,197,245,  4,255,253,  1,252,  2,  0,224,  6,175,253,254,251,252,  6,234,  4,  1,  0,255,  4,  3,  1, 16, 10,252,214,  1,252,  0,251,
    },
    {
       22, 57, 24, 18, 28, 59, 21, 15,  9, 64, 14,252, 19, 64, 11, 15,  4, 75, 21, 17, 23, 78, 20, 12, 19,207, 17,  2,255, 62, 13, 32,251, 12, 51,248,251, 11, 54,252,  1, 13,189,235,253, 10,190,255,255, 22,182,  1,254, 29,187,  2,  2,230,189,251,  4,  7,209, 37,  0,  1,239,  0,  0,  1,244,  0,  0,  1, 19,  1,  1,  2, 14,255,  0,  2, 25,  0,  1,  2, 29,255,  1,  1,230,  1,  2,  0, 12,255, 13, 32, 13, 19, 15, 16, 11, 11,243, 46, 17, 26,239,243,249,220, 60,192,197,212,209,214,  6, 30,225, 23,239,240,244, 19,201,166,232, 10, 12, 13, 13,245,238,247, 39, 50, 42, 41, 49, 73,223, 47,252, 95,235,248, 29, 28,  2, 12,234,249,255,  2,255,  3,240,239,  1,255,  0,255,  0,255,  1,  4,  6,  8,  8,  8,233,245,  4,254,  0,254,  0,254,  2,  1,  0,  7,  3,  1,255,255,  2,  2,  0,  0,
       27, 27,216,  6,250, 15, 40,254,  1, 22, 39,255,233,  2,223,254,  6, 16, 28,  3,242, 91, 32,  1,167, 72, 41, 12,221,  9,  6, 19,213, 34,112,241,  2, 10,145,209,255, 16,218,  8,  8,  3,247,255,250, 18,250,251,  7,  2, 91,  0,  8,  5, 52,238,  1,255, 68,192,  7,  0,255,  0,255,  0, 18,  1,  0,  0,  5,  0,255,  0,  4,255,  0,  1,225,255,  1,255,254,255,  0,  1,  0,  2,  1,  0,117,  0,  1,252,  0,255,252,253,  7,254,244,216,  2,247,  8,237, 32, 47,  4,247,254, 38,  7,249,219,176, 12,  5,247, 12,247,225,173,189,  2,243,  0,  2,  3,  6, 10,  6,221, 38, 34, 39, 91,189, 50, 53, 70,226, 45, 65,196,  0,254, 42, 17, 17,  1,  3,  9,251,253,250,  0,  0,  0,  0,  1,255,  5,  0,249, 33,  4,239,249, 18, 20, 28,  1,  1,250,  1,255,255,255,  1,  0,254,255,  2,  0,255,  0,  3,
      246, 22,  1,  5,  7, 11, 10,  1,249, 17,  2,  3,  2, 32,252,255,186,255, 24,  1, 21,196,246,  2,  7,183,247,193, 57,249, 15,239, 47,245,  9, 20,  3, 14,119,255,254,  7, 31,244,252, 13,  4,192,253,255,255, 19,  1,  8,233,253,  1,252,247, 61,  3,  1, 13, 62,254,  1, 27,  0,  1,  0,245,  0,  0,  0, 35,  0,255,  0,  9,  1,255,255, 69,  2,255,  0,  0,255,  0,255,255,  3,  0,  0,207,  1,251,  0,  3,  3,252,  0,  9,252, 33, 67, 19, 35, 40,  9,  6,240,  4, 36,  0,  4,187, 61,206, 69,  4,251, 61,197,215, 45,190,  5, 13,  4,231,242, 14,  8,235,246,202,225, 79, 16, 70, 69,  7,204,199,  6,193, 97,250,192,  1,221, 12,244,249,252,172,216,  1, 31,252,  1,  1,255,  6,  1,  3,  2, 24,248,252, 34,  5,  3,242,251,255,255,253,248,  1, 40,  0,251,  0,  8,227,252,255,  3,  0,254,
      252,200,245,  4,  2,245,  3,  0, 26,  1,  0,  1,243,191, 69,  0, 12,191,  1,  1,251,193, 10,  0,  6,250, 58,  3,234,219, 16,201,  1,  2, 56,  5,254,249,181,255,  3,  2, 38, 27,  2,  3, 66,  0,  0,255, 70,207,  0,  4, 65,224,  0,254,192, 10,255,248,  1,  3,  0,191,  2,  2,255,255,  2,  0,  0,  0,244,  0,  1,  0,  0,255,255,  0,228,252,255,  0,  2,  0,255,254,  3,  0,  1,255,  0,  0, 14,  0,255,  3,  7,  0,  8,  4, 26,  2, 17,255,254, 62,255,191,206, 49, 63, 66, 45, 12, 19,252, 30,  3,227,  5,252,  9,  0,253,  6,  0,  8,250,252,  0,236,  7, 53,240,195,137,213,224,122,202,  0,  2,  2,246,247, 45,  9, 42,248, 10,203,193,  2, 46,254,213,  8,  0,  1,  0,252,  0,254,  3,254,  7,  4, 11,251,  1,225,  1,  2,  0,254,  8,  0,254,  1,255,  2,245,  0,254,  2,  8,  0,  3,
    },
  },
  {  //  17
    {
       17, 17,221,210, 16, 14,222,208, 25, 14,220,207, 20, 12,221,208,220, 16,224,210,231, 15,225,212,218,  9,219,219,215, 10,222,216,233,199,191,240,235,194,191,239,227,195,188,242,232,191,187,240, 10, 67,180,241,  2, 66,190,244, 31, 66,186,254, 37, 63,182,248,  0,211,  0,  0,  0,209,  1,255,  0,216,  0,254,  0,210,  0,254,255, 46,  6,  1,254, 44,  6,  1,  0, 32,  4,  1,  0, 31,  5,  1,212,205,212,212,207,204,208,208,216,206,214,217,204,207,211,209,202,206,223,219, 69,195,208,205,227,216,231,226,254,196,224,194,213,215,221,227,203,204,213,214,205,205,206,206,210,211,208,210,  2,245,242,239,  5,252,249,255,251,253,250,239,  0,247,253,253,  2,250,  8,  3,  0,254,  2,  2,253,247,253,  1,254,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,255,  0,  0,  0,  0,  0,  0,  0,  0,
      206,231,253,227,189,209,  2,228,198,221,  1,234,188,244, 15,232,  4,242,  2, 78,  5,  1,250,215, 27,236,  8,  8, 27,243, 45,199,235,213,191,  2,  5, 21,249, 28, 36,  2, 36, 53, 41, 34,  8, 21,238,  0,246, 20,213,  1,244,230,244, 10,221,143,232, 15,253,254,255,251,  0,252,  0,250,  1,  8,252,250,255,252,252,250,252,254,254,248,  1,  0,255,249,254,249,252,  2,  1,254,254,255,  1,252,  2,  4,248,  1,252,  2,240,255,253,254,255,  4,253,253,255,  9, 13, 24,248, 13,208, 36,240,  8,245,240,230,240,  1,  1,245,246,254,  2,251,252,250,  6,253,251,  1,252,  3,255,254,251,  1,253,254,  0,  4,  0,253,  0,  2,254,254,  0,  5,254,254,  1,249,  5,  1,255,248,245,  0,  4,  3,251,249,253,251,250, 23,  2,237,254,  0,  0,  0,  1,  0,  0,  0,  0,253,253,251,253,255,  0,  1,  0,
      236,238,  0,  1,247, 15,  0,254,242,233,245,  0,224,254,228,252,251,235,227,238,235,  3,250,255,229,223,177,  5,223,251,213,255, 18,246, 64,231,240,194,199,204,201,198,202,177,196,172,200,201,222,190,233,222,250,199, 54,186,182,176, 18,208,198,176,186,225,192,193,194,196,192,202,192,198,193,207,197,197,196,217,197,196,194,194,192,190,193,197,192,186,192,195,194,192,192,195,193,192,240,255,  0, 15,243,241,  2,249,  2,248,  6,255,252,249,  0,  0,  3,249,250,244, 23,232, 12,  7,  6,  0, 13,243,  5,  1,  4,255,  2, 12,  6,254,253,245,241,239,  1,  0,  1,  1,249,247,251,247,  0,  0,  3,  1,253,  1,250,  1,253,250,  1,254,246, 11,  0,250,193,196,192,194,185,192,212,192,173,195,216,191,198,195,204,194,193,192,193,192,192,192,192,192,190,192,192,192,197,194,195,195,
      146, 23, 17,193,  1,247,241,194, 52, 10, 41,199, 15,  0,  8,201,  3,248, 40,251,  2,251, 36,253,  0,  0, 60, 15,  0, 10, 40, 31, 26, 36,212,242,135, 22,197,207,200, 31,190,208,  8,222,200,221,198,  5,194,  1, 11,  2,193,240,221,  4,207,  5,  0, 10,199,243,255,  0,  0,  1,253,  2,254,  2,253,254,253,255,251,  2,  0,  0,  0, 12,  0,255,  1,249,  0,  0,  0,  6,  0,  0,  1,243,  0,  1,226,245,166,224,214,243,202,233,252,217,176,  4,244,230,240,227,243,250,246,228,239,252,254,238,246,251,244,254,246,253,251,245,198,196,193,196,202,200,202,200,202,195,196,201,202,200,204,203,234,240,243,240,249,239,246,  0,242,245,244,245,249,254,255,255,  0,  0,  0,  0,  0,255,  1,  0,  1,  0,  0,  1,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,255,  0,  0,  0,  1,  0,
    },
    {
       12,246,247,246,  7,243,249,247, 13,236,250,247, 11,  0,254,247,  4,252,  1,248,  3, 15,255,247,  7,252,254,250,  1,  7,247,  8,  6,243,238,  7,  3,240,240, 12,  6,240,235,250,  6,248,232,252,  6,244,243,253,  3, 15,253, 55,  5,239,234, 63,  3,  9,189,198,  1,  0,  0,255,  1,  0,  0,255,  1,  1,  0,  0,  1,  1,  0,255,  1,  4,  1,254,  1,  0,  1,  2,  1,  0,  0,  2,  0,255,  4,252, 74, 33, 29, 33, 28, 29, 27, 39,226,234,225,233,214,229,218,203,167,213,216, 20,  7,255, 10, 11,231,183,229, 51,  6,225,179, 36,253,247,249,246,249,245,248,  2, 17, 21, 17, 22, 20,  8, 14,  4,  0, 26,  0,  2,  1,  3, 11,  2, 36, 11,  0,  4, 64, 53,201,228,236,  1,  0,  1,  3,  4,  3,  1, 23, 29, 38, 89, 64, 64, 66, 72,  0,  4,252,255,241,  0,252,  3,182,255,252,  0,  5,  2,255, 10,
        5, 26,  8,242,255,241,247,  5,249, 21,  1,  1,  7,  5,249,247,255,  8,252,  0,254, 21, 13,251,254,  9,253,254,  0,251,  2,  1,  4, 34,198,180,248,240,243,197,251, 28, 17,  5,252, 24,207, 77,  1, 10, 13, 44,253,  9, 21,238,253,251,187,213,  0,243,243, 52,  0,  0,244,205,  0,  0,  0,  0,  0,  0,  0,  0,  0,255,244,  4,  0,  0,  1,  4,  0,  1,  0,  1,  0,  0,251,  0,  0,  1,  0,  1,255,  3,255,  3,  0,255,  0,  0,  4, 41, 49, 36, 65, 65, 63, 62, 25, 44, 24, 38, 37,  1, 18,  1,197, 16, 42,138,145,118, 15, 15,243,  5,211,  9,212, 30,  5,255, 18, 21,  9, 20, 60,  2, 51, 64, 18,231,226,  3,255,  1, 19,  5,198, 24, 79,  0,245,198, 44, 28, 15,  1,  2,  0,254,  1,  0,  0,255,247,  0,180,212, 75,  2,249,202,253,  2,255,  4,  0,  2,253,  2, 34,  8,  5,253,224,215,194,
        4,252, 19,245,255,  5,  1,  4,187, 63,247,253,255,246,  0,  5,255,  1,  0,  6,  7,254,  0,  6,  0,  2,254,  4,  0,  7,251,  0,  1,236, 69,  0,  7, 31,198,  3, 36,195, 62,190,  4, 32,197,255,  0, 10,191,214,  1,252,198,236,255,255, 52,189,255, 13, 47,  0,255,  0,252, 40,  0,  0,249,  0,  8,  0,255,  0,  0,253,255,  1,  0,255,250,  0,  1,254,252,255,  0,255,252,255,  0,  1,252,  0,232,253, 12, 47,255,228, 11,  0,251,253,254,  0,  0,234,  0, 15,236,159,248, 82,215, 68,254, 56, 54, 28,222,190, 63,  0,230,199,253,  0,239,255, 52,195,112,  0, 51, 41, 16, 41,225,  3,  2,238, 16,242,191,230,222,199,  5, 67, 13,  6,174, 64, 51, 46,216,199, 11,  1,  4,252,254,  2,  8,  0,202,254,241,  0, 24,150,246,253,  5,  0,  2,  0,  0,249,  2,235,250,246,215,  5, 16, 16, 60, 16,
      255,  0,244,  0,253,249,254,255, 51,189,  0,255,255,  0,  2,255,  3,246,  0,  0,  2,254,251,255,  0,  1,  8,255,  1,  1,  5,255,  0,  4, 12, 64,  2,217, 59, 10,194, 32, 15,221,  0,221, 65,254,  0,250, 68,  0,  1,  1, 40,  2,  0, 27,251,251,  0,253,252,254,  0,237,253,  0,  0,255,  1,  1,  0,255,  3,  1,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,255,  5,  1,  0,  0,255,  1,254,  0,249,170,  2,  5,254,  1, 27, 12,  0, 11,  2, 43,  2, 10, 74,113,241,233, 10, 66, 39,205, 31, 49,254,187,203,208,121, 25, 30,  2, 48,  3,197,251,179,253,204,218,  0, 31,  8, 41,  4,243, 46,223,108,195,  3,171,227,195, 11, 15, 28, 43,  5,197,  4, 65, 13,  4,  4, 64,  2,  1,  0,  2, 18,255,  5,255,  8,  8,213,235, 18,  0,255,  1,255, 10,  3,254, 78, 10, 40,  1,  0,  3,  0,255,
    },
    {
      222,253,237,198,228,252,235,200,229,252,234,200,224,250,234,200,223,253,235,199,219,251,233,193,218,  4,240,184,222,242,229,206,218,242,250,226,218,244,250,230,215,241,251,228,211,241,253,231,214,241,254,225,218,238,253,223,214,238,  3,221,216,245,250,216,  3,  0,250,249,  0,  0,250,249,254,  1,250,249,  1,  1,250,249,255,  1,250,249,255,  1,250,251,  1,  1,249,250,252,  1,250,249,217,225,223,227,218,224,223,228,210,215,219,220,209,221,219,220,200,206,208,211,204,208,209,208,186,191,196,192,192,209,195,195, 35, 35, 29, 31, 32, 32, 30, 36, 52, 49, 50, 52, 52, 52, 50, 49,247,  0,246,  4,253,252,245,  3,  3,  3,252,  5,  8,  7,251,  7,236,237,241,238,236,238,241,236,220,221,230,222,219,227,230,221,247,247,248,248,252,249,248,248,250,251,255,251,251,252,255,249,
      249,250,  7,254,250,254,  4,  1,249,252,  6,255,249,254,  5,  0,249,250,  6,255,250,254,  3,  1,251,253,  8,255,248,254,  4,254, 10,239,  8, 21,  6,236,  5, 19,  7,234,  8, 17,  1,234,  4, 15, 11,237,  8, 19,  5,236,  4, 16,  8,233, 11, 17,  4,236,  3, 23,202,193,198,196,202,192,198,196,201,193,198,196,199,192,198,195,203,193,198,196,203,192,198,196,198,193,198,196,197,193,199,196,251,250,247,242,246,243,239,239, 14,  4,  4,254,  5,  2,  0,  1,  3,254,  0,251,252,253,247,250, 27,  7,  2,  4, 10,  9,  1,  9,230,233,233,222,223,224,228,223,230,222,231,226,221,222,223,220,249,249,248,243,242,241,244,241,253,249,251,251,246,248,251,247,210,208,208,210,209,210,207,211,215,208,205,214,214,213,207,211,195,198,194,199,197,202,198,202,196,199,189,197,199,199,193,198,
       10,249,240,230, 16,250,237,233, 11,243,207,210, 13,  2,230,209,242,  5,244,229,253,  2,241,222,222,251,239,223,244,253,239,221,191, 10,224,234,196, 21,228,227,191, 30,228,221,192, 16,228,236,  5, 19,227,227,  4, 19,237,232,253, 16,233,234,  3, 21,232,232,254,  3,251,255,  1,  0,253,255,255,  2,254,  0,  0,  1,255,255,  9,  4,  0,255,  6,  2,  0,255,  1,  6,255,255,  0,  3,  0,255,  3, 12,  2,  4,  5,  4, 33,248,250, 24,232,244,194,194,183,200,232,  6,  0,240, 17, 20,225,224,225,205,192,197,234,232,218,230, 15,238, 14, 19,  4,  9,  1,  3, 15,248, 15, 15,  3,  6, 15,  2,  2,254,  0,254,  0,  1,  4,  4,253,252,  6, 11,  7,  0,  0,201,  1,  2,  0,  0,  1,  1,  1,  0,  0,  3,255,252,255,  0,255,248,  0,  1,254,  0,249,249,249,251,  0,  0,  5,  2,  1,255,255,  0,
      187,201,215,  0,189,204,219,253,190,202,222, 19,190,202,246, 21,237,226,210,194,199,245,217,219,220,243,220,197,237,244,224,190,200,251,  3,252,191,250,  5,  7,183,250,  0,251,195,248,  1, 11,189,190,  5, 21,191,197,253,  0,193,185,  0,  6,199,188,  1,  4,  0,  1,  0,254,  2,255,255,  0,255,255,  0,  0,255,  0,  0,255,  0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,255,  0,  1,  0,  0,226,223,221,229,229,239,206,254,239,217,201,239, 40, 37, 40, 40,219,230,211,228,197,201,246,250, 21, 15,251,252,240,240,  0,  2,193,201,217,194,211,221,230,215,188,199,219,187,212,211,209,208,190,196,193,193,195,197,189,190,188,190,193,192,202,195,245,198,255,  1,255,255,  0,  1,254,  2,254,  1,  1,  2,252,254,  0,  1,  0,  0,255,  0,  1,  2,  1,  2,  0,  0,  0,  0,255,  1,255,  1,
    },
  },
  {  //  18
    {
      206,213, 52,  7,207,215, 50,  5,208,216, 44,  8,205,215, 43,  5,206,217, 38, 10,204,214, 38,  5,203,209, 34,  3,198,207, 36,  1,250,209,236,  0,250,212,239,  0,254,214,243,255,253,214,246,  0,  1,217,247,  0,254,218,247,  0,253,209,252,255,252,205,255,  0,254,252,  0,253,254,252,  0,253,254,253,254,253,254,252,255,253,253,253,253,253,254,253,254,253,253,254,254,253,253,254,  1,253,239,241,239,243,241,240,239,243,223,228,226,230,228,228,224,230,204,215,206,214,210,205,205,213,212,221,208,226,209,203,217,225,  5,  9,  6,  5, 13,  8,  6,  5,  9,  7,  8,  8, 16,  9, 12,  8,187,193,189,187,195,192,193,189,207,214,211,224,213,208,205,223,  2,  2,  1,  3,  0,  1,  1,  4,251,251,252,251,251,253,255,251,244,238,245,243,240,245,243,242,233,237,235,240,230,235,221,240,
      243,252,209,  4,245,252,214,  3,250,254,213,  2,246,253,214,  4,246,  0,214,  6,247,254,223,  1,244,253,225,  3,246,254,221,  2,213,248,  5,247,211,245,  5,246,210,249,  5,249,207,248,  5,249,203,249,  8,250,204,241,  6,244,206,241,  0,247,202,237,253,239,217,191,191,195,216,193,191,195,211,194,192,195,210,195,193,195,226,195,193,195,194,195,193,195,233,195,193,195,221,195,192,195,248,246,  0,247,214,247,230,245,242,238,254,243, 42,253,232,248,  7,  9, 17,  0,250,  5,  9, 11, 10, 21,254, 11, 24,  6,240,  0,217,212,236,229,237,234,220,210,214,215,197,211,255,219,220,205, 20,125, 74, 38, 18, 16, 58, 32, 46,  3,207,  7,217,216, 17,245,192,191,192,194,192,192,198,193,196,199,204,206,192,193,200,204,194,191,192,196,191,198,194,198,194,191,193,191,192,194,193,194,
        7,255,254,  7,254,  0,  4,255,251,  0,250,  2,  3,253,251,  6,219,212,  5,251,  9,253,  0,251,221,222,  0,253,248,231,  0,252,  4, 56,  1,252,251,212,  1,  9,255,197,  3,249,  4,210,  0,254,253,  0,253,242,247,217,247,247,255,  0,255,253,245,244,  0,251,234,255,  0,  0, 26,  0,255,  0,236,  0,  0,  0,244,  0,  0,  0,  0,  1,  0,  0,255,  0,  0,  0,  2,  4,  0,  0,255,  1,  0,  0,255,254,236,217,254,234,  6,234,248, 17,206,254,229,233,187,195,  4,  0,  1,253,246,254, 15,  5,  0,  1,254,253,194,243,  7,216,240,153,176,247, 13, 23,  4,  4,211,173, 47,  5, 63, 63,222,229,243,243,228,220,193,216,231,194, 14,  0, 17,  2,202,247,243,199,  0,  0,  4,  0,  0,  0,253,254,  6,255,  1,  1, 11, 16,255,254,  2,  0,  6,255,  1,  0,253,253,255,  0,  1,  2,254,  0,  1,254,
       65,243,250,248,192,253,252,251,  1,253,254,247,  1,  1,242,255,235,231,253,193,191,190,250,192,240,233,253,216,211,217,  7,213,142, 51,254,255, 15,207,255,248,200,145,  0,  0,198,  2,  4,  1,203,197,  0,201,239,204,  0,197,203,200,  2,202,214,209,  0,201,  0,  0,  0,  0,  0,  1,  0,  0,252,  1,  0,  0,253,255,  0,  0,  0,  0,  0,253,  0,  0,  0,254,  0,  0,  0,  1,  0,  0,  0,  2, 50, 64,101,229,237,233,247,237,197,143,213,134,214,219, 45, 40,210,139,216,216,221,222,252,207,204,199,205,206,254,219,216,213,192,133,  2,  1,248,244,199,233,  3,183,255, 87,215,208,246,  0,  1, 10,250,242,  6,  7, 54,  8,  3,253,253,250,  9,  5, 50,  5, 13,  0,  0,  1,  1,  1,  0, 63,  6,  0,  0,  0,  0,  0,  1, 66,  1,  0,  1,255,  1,  0,  0,  1,251,  0,255,  0,  5,  0,255,  1,
    },
    {
       11, 21,  4,199, 12, 21,  1,204,  7, 19, 18, 54,  7, 21, 17, 86,  7, 17, 15, 95,  8, 23, 21,207,  7, 15, 79,154,  0, 15, 66,187,  3,205,  0,238,254,203,254,241,  3,207,253,  0,  5,205,  3,  3,  6,203,  3,  2,254,199,255,251,  6,197, 12,  0,  3, 67, 17,249,  0,251,  0,  1,  1,253,  0,  0,  1,  7,  0,  0,  1,  7,  0,  0,  1,  7,  0,  1,  1,  4,  0,255,  0, 21,254,  0,  1,242,  0,  1, 24, 26,235,212, 22,207,133, 14, 33, 28,225,222, 24, 27,227, 25, 54, 66,253, 13, 60, 51, 62,194,242,244,238, 46,  7, 10,223,224, 19, 22, 17,  3, 16, 12,233,236,  5, 14,242,253, 13,255,253,248,250, 26, 64, 62, 70,217, 20, 10, 29,193,  5, 56,192,192, 60,191,  0,  0,255,255,  0,  3,254,255,  8,255,  5,  3,  0,  0,255,  1,250,  4, 19, 23,  6,252,255,  1, 33,229,255,  0,230,241, 31, 14,
        2,  0, 64, 76,  2, 17, 56, 81,  4,246, 62,209,  2,249, 69,180,  1, 24,125,198,  0, 16, 19, 40,254,209, 76,236,255, 10,204, 57,246,249, 24,255,254,  5, 60,254,250,  1, 14,254,253,  4, 23,  1,192,  5, 41,  0,  2,254, 44,  5,246,  0,  0,  0,  0,250,  0,253,255, 31,255,  0,  0, 22,254,  0,255, 17,  0,  0,255, 17,  0,  0,  0,  2,254,  0,  0,  1,  0,  0,255, 17,  0,  0,  0,247,  1,  0, 85,218,219, 26,253,232,  1,117,  0,217, 13, 41,  2,255,  2, 31,  0,234,  5,  3,  2,131,  1, 34,254, 92,  0,245,232, 12,190,188,226, 49, 21, 51, 38, 86, 22, 21,254, 16,  9,247,  9,250,  3, 15,  0,  1,251,  0,  0, 44,249,  5,239,198,136, 24,  0,  1,131, 56,  0,254,254,255,  1,253,254,  0,  7,  2,  1,252,  0,  0,  0,  1,  4,  2,  2,  3,  3,  2,255,225,230,254,255,  5, 26,  2,  1,253,
        1, 28,  2,  7,  0, 15,  4,  1,  1, 11,241,  2,  0, 30,128,252, 38,  1,251,  4,  0, 29, 15,247,  0, 59,218,114,  0,  7,244,238,249, 65,217,  1,247, 60,195,  0,  7, 65,245,  0,  2, 64, 10,254, 64, 56,201,  1,  1, 61, 40,247,254, 62,  5,255,  0,193,  1,  0,  0,253,  1,  1,  0,253,  0,255,  1,  3,  0,  0,  0,  3,  1,  0,254,  1,  1,  0,  0,  0,  0,  0,  1,  0,  2,255,255,  4,  1,255,203, 23,  8, 64, 71,210,226, 69, 10,235, 12, 20,  1,233,246, 17, 16, 11,178, 62,251,187,  3,188, 87,183,199, 29, 16,  6,251,145, 24,228,243,225, 47,213, 10,243, 26,  8,  3,172, 10,254,248,209,217,  7, 36,209,114,254, 10,196,207, 50,  0,245,  1,197,157,112,  6,  0,  1,254,  0,  0,255,  0,251,  6,  0,255,  0,254,  1,  6,252,250,  1,254,  3,  0,255,  5,247,253,248,255, 13,  4,  1,  0,
       19,235,164, 50,  7,255,200,  4,  1,242,214,  1,  3,233, 61, 70,255,242,199, 38, 19, 30,192, 13,  8,231,219,  0,  1,255, 21,  8,  1,  2,254,  8,254,255,252,  3,254,252,  8,  0,254,  1,  0, 22,  0,248,199,246,254,  1,  0,  9,254,  6,231,  0,254,  5,  1,255,  0,  1,  0,  1,  0,  0,  1,255,  0,  1,  0,  0,  0,  0,  1,  1,  0,  1,  1,253,  2,  0,  1,  0,  0,  0,254,  0,  0,  1,  0,  1,223, 22, 70,  4, 60,238, 58,  1,244, 41, 25,240,238,  1, 19, 10,180,194, 37,221,191,252,197,  5,204, 64, 40,245, 12, 23,253, 79,203,253, 11,247,176,218,214,253,247,  5,252, 52,250, 25, 50, 65, 59,245,199,212,255,205,189,115, 31, 76, 73,171, 63, 67,247, 21,248,  3,252,  0,  1,255,  3,  3,254,  2,255,254,  2,  3,  1,  0,253,  3,  1,255,  0,  1,  3,255,  4,254,  0,  1,  0,  1,  0,  0,
    },
    {
        0,252,228,203,248,245,226,201,253,249,233,200,253,250,233,197,246,243,235,200,253,245,233,196,250,245,238,195,248,244,238,196,208,221,  1,240,213,221,255,240,209,222,  1,245,211,223,  6,252,208,225,  3,254,208,218,  0,252,206,224,  3,  3,206,226,  3,  4,254,250,254,  2,  1,250,253,  2,  1,254,254,  1,  1,  0,254,  0,  2,255,254,  2,  0,255,254,  4,  3,  0,255,  2,  4,  0,255,  2,241,244,251,246, 22,244,  0,250,243,240,245,243,255,241,  1,244,253,252, 10,254,  3,253, 15,  3,231,233,231,234,233,232,233,234,228,223,230,224,229,231,234,231,  0,251,  1,253, 11,253,  3,254,220,215,223,217,224,218,214,221,  5,  3,  6,  3,  5,  5,  4,  2,246,247,252,248,241,246,252,247,244,246,247,246,245,244,248,244,251,250,250,252,245,247,252,246,253,254,254,255,  0,  1,  0,  1,
      236,250,  4,  1,240,254,  1,  2,239,252,255,254,239,247,  2,  1,238,251,254,  5,235,251,255,  6,242,253,254,255,239,246,111,  7,253,234,252, 17,252,234,  0, 21,253,236,  1, 16,247,232,253, 14,254,233,255, 12,252,237,  1, 15,253,230,  0, 14,237,223,229,239, 12,  6,  0,252,  2,  5,  1,251,254,  4,  0,253,253,  1,  0,254,  1,  4,255,253,  5,  6,  0,251,254,  5,  0,252,  1,254,  3,255,239,237,233,230,184,227,227,231,215,198,223,218,215,228,213,215,236,234,217,227,219,229,215,228,255,254,253,251,  0,255,  2,  1,243,248,243,244,221,254,233,251,240,243,242,246,250,237,243,249,231,231,227,231,227,246,225,241,254,  1,252,253,  0,  0,  3,  0,207,207,210,210,209,207,203,208,226,223,223,227,226,222,215,217,219,223,212,215,215,214,212,218,193,192,193,194,192,193,191,191,
      211,201,196,197,212,207,192,202,205,207,196,200,207,207,207,197,189,214,199,201,216,207,189,204,188,216,198,201,187,209,201,199,  9,251,246,252,  1,250,244,255,  5,253,249,252,  5,  4, 19, 27,254,255,241,252,249,247,242,254,254,  0,244,251,255,  0,249,252,197,197,193,193,192,190,194,193,193,198,194,194,188,192,191,191,199,209,194,194,190,190,194,194,197,209,193,194,192,201,193,194, 44, 61, 43, 48, 23, 22, 16,255,192,191,127,188,195,196,195,197,  5, 15,  1,  1,253,250,  3,250, 13,  9, 10,  8,248,237,  8,247,190,192,187,190,192,197,193,198,232,219,253,229,227,229,224,251,210,214,207,215,214,217,229,224,254,  6,  6,  5,  2,251,228,255,  0,  0,  1,  0,  0,  0,255,255,217,209,218,218,216,219,228,231,253,241,254,  0,  0,  2,254,  4,249,253,255,246,  4,252,246,244,
      150,191,  9,  3, 25, 42,  9,252, 20, 46,  8,  1, 19,  0, 11,  2, 21, 49, 28,244, 24, 44, 10,  3, 28, 45, 38,246, 20, 46, 34,  1,199,193,253,  5,198,199,251,  0,197,197,246,  3,191,192,243,253,199,200,  2,  0,198,198,252,255,203,208,  1,255,194,194,250,248,  0,  2,  1,  1,  1,  0,255,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,255,254,  2,  1,  0,  0,  1,  0,254,246,255,255,252,254,251,  3,189,196,189,188,214,217,221,180,189,189,191,192,204,207,205,228,  1,  1,  0,  0, 20, 17,  3,208,232,  9,  1,  1,  2,255,  3,255, 46, 46, 44, 45, 38, 38, 33, 14,253,255,  2,  3,250,252,251,232,251,248,  0,  0,215,220, 45,254,195, 14,  2,  2,  3,  3,  0,255,255,252,250,250,254,254,245,255,  0,  0,  0,  1,  0,255,255,254,226,226,236,232,241,239,  1,246,
    },
  },
  {  //  19
    {
      213,236,  2,221,214,220,197,206,212,233,242,202,216,246,244,237,213,238,245,197,207,207,255,170,218,240,183, 89,253,134,129,209,210, 68,  0,254,202, 67,  1,255,208, 66,  2,  3,206, 65,251,  1,207, 66,  2,  2,198, 65,  0,195,200, 68,250,192,188, 61,127,196, 16,254,255,  0,  3,253,  0,  0, 11,  1,  0,  0,  4,  0,255,  0, 10,254,  0,  0, 73,  0,  0,  0,254,  1,  0,255,249,255,255,  0,  2, 26,229, 39,  2, 47,212,182,  9, 14, 29, 29, 16, 63,244,221,  9, 40, 29, 69, 11, 64,187,  1,  5,246, 70,215,103,192, 73,127, 72,249,247,227, 12,196,159, 63, 13,  4, 43, 15, 11,  7,254, 10, 10, 73,201, 38,251, 73, 54,  1,  4,251,231, 46, 63, 62, 66,  2,240,249,  2,  8,  0,255,  0,238,  2,  4,254,250,  0,  0,255,  0,  0, 13,243,  1,255,  0,  0,  1,255,  0,246, 21,  0,249, 13,  2,
        8, 21,252,247,216,  0, 53,237,231,249, 13, 19, 11, 33,244,232,251, 13, 70,106,202, 36,  6,194,  5,238,  0,193, 69,246,  5,  3,  3,  2,  3,254,  5,  1,  4,  5,254,  6,  6,  0,  0,  4,  2,  5,254,  1,250,  1,  3,255, 56,222,  0,  4,188,173,  0,  2,  0,  0,204, 62,199, 57,193, 61,199, 57,195, 61,199, 57,196, 61,199, 57,192, 67,199, 57,191, 64,199, 57,188, 63,198, 57,  7, 63,198, 57,254, 20,  1,  1,231, 29,  5,  1,204,253,226, 49,200, 47,226, 48,186, 64,175, 52,216, 72,175, 52,  2,  0,162,250,225,248,166,250,191,114, 22,  3,134, 33, 22,  3,164, 65,207, 41,204, 56,207, 40, 15,201, 18,  9,205,192, 20,  9,191,182,254, 14,251,235,254, 16,202, 64,187, 64,200, 74,188, 64,191, 43,188, 60,191, 66,187, 60,193, 64,198, 49,218, 65,199, 49,192, 72,214, 51,191, 57,217, 51,
      191,234,225, 22,255,235,218, 22,255,  1,212,  0, 69,252,246,  6,237, 31,  5,225,251,243,  4, 13,236,244,234, 11, 12,  3,245,252,255,  0,254,  0,254,  1,  0,255,251,  0,  1,  0,  1,  0,  2,255,  1,  0,  4,  0,  8,252, 10,254, 14,  2,  0,253, 10,253,  5,251,  0,  1,249,  7,253,254,250,  7,251,  0,250,  7,251,  1,250,  7,254,254,249,  7,246,253,249,  7,251,  7,249,  7,255,  3,249,  7, 65,217,208, 67, 39, 23, 12, 52,197,  1,167, 83,  0,  4,213, 17,  1, 47,202, 98,242, 33,222, 46,156,254,156, 80,231, 16,208, 14,194,255,  5,180, 40,  5,241,  3,253, 60,253, 49,189, 65,189,  1,  1, 34,234,182, 51, 54,217, 42, 11, 72,253, 71,191, 43,219, 65, 40, 10,  5, 11, 22, 18, 16, 16,251,  1,  4,  9,255,  7,253, 39,  0,  0,234, 20,248,  7,245,  4,255, 11,239,255,  1,  3,  0,  8,
       28,130, 46,  3,253,241, 29,  4,192,  7, 44,246,208,  1,234,252, 38,240, 51,252,  3, 78,240, 16,  0,  0, 35,  0,  6,248, 19, 29,  0,195,251,255, 61,176,254,  0,  1,254,  7,255, 17,189,193,254,  3, 68,193,  0,193, 63,250,  1,185, 63,195,  0,178, 58,187,251,246,  0,  0,  0,  2,  0,  0,  0,  0,  3,  1,  0,  0,  0,  1,  0,255,  3,  1,  0,255,255,  1,  0,  0,  0,  2,  1,  0,  0,251,  0,  0, 48,136, 61,242, 55,222,  5,194,192, 97,198,  2, 38,  0,252,255, 64,248, 76,  9,199,189, 16, 37, 61,223, 56,191,247,199,234, 78,239,197, 67, 61, 56,208,227,199, 64, 34,  1,  1, 34,  2,233,  3, 66,186,219,201,242,215, 17,212, 74,166,242,  5,  1,175,247, 10,  0,  5,251,  1,253,247,234,  0,255,252,255,  2, 14,193,219, 75,  3,  0,  0,  1,  0,182,248,  3,  2,250,  1,251,254, 20,248,
    },
    {
       43,224, 41,176,245,  2,208,252, 48,251, 25,177, 25,255, 57,247, 25,249, 52,249, 12, 14, 51, 20, 51, 32, 49, 70,232, 70,162,185,254,  9,238,  8,  0, 21,237,254,  1,  2,251,  9,249,252, 20,  1, 11,  2,234,255, 18,  1,  0,  1,  3,  4,128,246,  5,  2,102,229,253, 13,  0,254,  6,247,  1,  1,254,  3,  1,253,254,255,  0,  0,  0,  4,  4,  0,255,  2,  1,  0,255,  1,  2,  4,  1,  0,255,  0,243, 11,240,  6,  7,251,  9,254,  1,  1,255,242,252,252, 16,251,254,  4,  0,239,234,248,  0, 13,229,248, 20,235,  4,250,  8,  7,216,  6,229, 96, 62,208, 72, 74, 12,201, 12,  4,245,255, 71,134,242,236,  2, 77,160,210,  2,230,237,249,220, 85,  4,192,194,187,224,227, 25,  4,  1,255,241,252,248, 19,253,255,  0,255,  1,255,254,220,  4,  1,255,251,  0,  1,249,  6,  4, 12,  0,245,  0, 20,
       10,254,234, 42, 26, 22, 32, 41,  2,234,235,255,184, 43,187,252,228,  4,172, 62, 65,242,196,183, 60,192,201, 21,237,140,  0,225,  6, 18,193,  0,255, 23, 57,255, 13,246,205, 10,  3,239,191,  3,252,  1,184,  3,254,  0, 68, 84,255,255, 62,  7, 12,251,250, 62,255,251,253,254,255,116,252,  0,255,253,  5,  0,  1,131,238,  0,  0,  0,255,254,  0,  0,  1,254,  0,  0, 10,  3,  0,  0,255,248,252,  5, 10,255,  9, 13,  2,  0,253,255,  0,255,252,255,  2,  1,255,255,  1,255,255,  3,252,  5,252,223,246,  3, 16,254,  1,254, 56, 16, 37,219,  8,  9,204,195,222,  0,200, 11,237,238, 50, 63,226,255,193,238,253, 91,212,  8,241,  7,234, 52,226,111,123,134, 10, 10, 17,  4, 12,  4,249,  2,239,255,  9,  2,  1,  0,  7,251,207,  4,  9,  8,  1,255,250,249,247,  5,  4,  3,255,  0,  0,190,
      223,227,198, 17, 22,229,  2,  1,246,246,200,255,195,231,246,  1,239, 74,  2,126,229, 83,127,127,186,251,202,133,244,223, 11, 53,252,180, 21,201,253,  9,186,  4,250,  6,206,  2,  1,253,117,  3,  2,  2,189, 40,255,244, 63,171,254,  0, 79,  3,254,  2,191, 44,  2,252, 10,  0,255,135,  0,255,254,  1,226,  1,  0,130, 15,  1,  1,  0,255,255,  1,255,  0,  0,  0,  0,255,  2,255,  0,  3,  0,  6,  4, 11,252,  3,152,  1,  0,  5,  0,226,252, 12,254,  3,249, 10,226, 18,  3,  3,  9,  1,251, 30,237,255,253, 59,  4,  4, 18, 22,  8, 10,247,  1,  4, 12,210,250, 54, 36,  5, 63,156,194,209, 12, 36, 73,244,253,253, 53,252, 13, 38,254,237, 67, 17, 20, 22,  8,251, 49,251, 16,227,  2, 18, 13,254,  0,  2,  2, 15,  0,131,  4, 32,  9,246,  3,254,  0,  6,  0,  4,  0,246,249,  1,  1,252,
      248,  0,  9,255,236,236,127,  2,245,244, 61,  0,240,212,249, 14, 31,242,249, 69,210,140,252,199,211,253,  1, 30, 63, 67,191,195,  1, 70,199,  4,  5,245, 62,  7,  2,  0, 63,247,249, 24,251,253,244,  4,  6,  2,253,  1,254,  1,255,  1,  0, 69,  0,  2,204, 85,253,  0,190,  1,  1,255,255,254, 12,  1,  0,  0,  2,  5,  0,  0,255,  1,  0,  0,  0,  0,  0,  0,  0,  0,254,  3,  0,  0,255,  4, 11,238,  5,245, 26,242,  0,247,  6,245, 17,255,  0,  5,  3,  7,250, 14,252, 17,243,248,  2,  1,245,  0,  8,232,253,  1,  9,120,224,248,  0,228, 71,233,  0, 23,  4,  0,239, 47,  3, 53,  6,145,248, 15,  5, 14,  2,253,  0,247, 56,  8,217, 16,254,255,202, 47,  8, 14,185,255,  2, 30,254,  2,  4,  4,254,  2,255,  6,  0,123,245,  8,249,  2,  0,  5,254,253, 10,255,  1,254,  0,254,255,  2,
    },
    {
       15,204,202, 79, 11,200,196,208,189,203, 21,  1, 26,194, 12, 10, 14,203,232,  7, 15,196, 28, 13,202, 14,119,204, 57,226, 15,192,  2,232, 13, 12,  2,242,242,248,  4,  4,194,246,  1,247,195,236,  3,235,197,251,253,247,218,255,  7,254,  0, 15,248,254,  1,113,  1,250,  1,  0,  0,252,  1,  0,  4,  0,239,  1,  0,  0,246,  1,  0,  1,252,  1,  1,255,252,  1,255,  0,255,  1,  5,  0,  0,254,174,233,221,219,208,223, 41, 14,255, 17,  0, 15,119,236,244, 24,238,251, 18, 50,209, 16,204,235,252, 43,184,225, 81,112,161, 35,  0,249,  4,247,  4,255,  6,  5, 12,253, 41, 15,  4,  3,252, 87,  4, 38, 29,199,208, 78,255,113,254,  3, 36,251, 42,131,192,  9,  6,  2,  0,  2,255,255,  5,  1,  1, 19,117, 19,  2,254,  4,252, 14,249,253,251,  3,227,  0,253,  1, 12, 17,  1,240,  2,242,  2,
      206,134, 72,213,  8,242, 65,201, 68,254, 33,255,236, 26, 25,134, 39, 61, 13,253, 11, 20,239, 74, 35, 18,155,189,190, 60, 46,252,255,247, 68, 61,255,246,254,247,  3,247,219,  2,255,  1, 48,180,  0,254, 47,  2,251,  0,251,253,254,252,254,  2,  0,250,197,  7,  1,  0,230,  1,  0,254,  2,  0,253,  6,  4,  0,255,  0,  1,  0,  2,  0,  1,254,255,254,254,  1,  1,  3,  0,  0,  0,  1,254,  0, 24,  2,237, 12,196, 47,255,191,246,214,176,213,157, 11,204,249, 33,252,235,209,240,  8,  5,180,250,226, 23, 16, 28, 10, 79,187,243,254,247,253, 53, 32,  2,106,254, 20, 10,  7,  6,  5, 36,254, 19,197, 95,  3, 16,127,143,121, 34,222,198,244, 12,180, 50, 73,252,  2,  0,  2,254,  2,  0,  3,253,  1,132,255,254,255,255,  0,  9, 15,252,  7,255, 12,255,  3, 59,220,  4,209,  7,233,254,  0,
       25,255,246,255,245,111, 16, 33, 14,108,242,  3,237, 81,180, 88,213,219,152,  4,229, 97, 24,203, 58,109,168, 38,  1,126,241,123,  0, 11,188,249,  1,251,234, 46,  0, 18,191,  2,254, 39,152,  3,  8,  1, 70,255,253,254,255,249,255,254, 65,252,  0,234,252,251,  0,  0,252,  0,  1,  0,  1,  2,  0,254, 42,255,  0,253,  2,  0,255,255,  0,  2,  1,  0,  0,  1,  0,  2,247,255,  0,  1,254,  0,222,250, 14,249,255,  9,249, 12,255, 65, 58, 47,231, 40,255,245, 21,248,212,217,248,188,  1, 26, 17,215, 29,192, 67,226,  0, 34,226,255,250,254,254,214,178,166,  0, 31,234,249,  3, 20,194,133,246,  3,148,254,252, 10, 28,181,191, 44, 47, 89,132,192, 19, 51,  0,253,255,255,255,254,  1,  3,245,216,247,  2,  2, 30,254,  3,199,250,224,  2, 14,191,211,252,254, 16,  0,253,254,254,  0,255,
       25,213,255,255,253,205,246, 57,  2,253, 34,255,248,203, 60, 13, 25,249, 27,246,254,215,130,243,254,  7,243,  0,237,236, 59, 70,  2,  1,  1,251,253,  6,253,232,255,211, 68,253, 25,  0,225, 42,  3,255,216, 31,  0,  2,  3,  2,  0, 24,182,  0,244,255,246, 15,  0,255,253,  0,255,  3,255,255,  1,255,  3,255,  0,  1,254,  2,  1,254, 12,254,  0,  1,  0,  0,  0,  1,  1,  0,  1,  0,  0,255, 83,  0, 17, 13,  3,219,  2,223,205, 10, 12,236,194,255,189,138,209, 15,201,196, 37,193,  1,230,  1, 12,215,254,127,  1,255,174,  0,253,  2,250,248,  0,  3,127,  3, 10,  0,253, 10, 18,253,233,240, 50,187, 71,193,245,183, 11, 21,  3,190,237,242,243, 71, 60,  4,  3,  9,  3,200, 12,  2,  0,  1, 28,248,253,  4,  2,  0,253, 61,240,  0,254,  1, 52,254,244,233,243,  2,  4,  2,  0,  0,  4,
    },
  },
  {  //  20
    {
      206,191,212, 12,202,195,210, 12,196,196,209, 20,193,194,207, 16,191,195,206,  5,190,192,209, 26,188,196,203,  9,192,198,216,248,199,254, 26,192,201,  0, 18,193,201,  2, 17,190,199,  2, 15,191,202,  4, 16,189,199,  2, 14,191,206,  5, 14,192,203,  5, 21,191,  1,  0,255,255,  3,  0,255,254,  2,  0,  0,252,  4,  0,  0,251,  5,  0,  0,248,  5,  0,  1,248,  6,  0,  0,244, 11,  0,  1,244,200,206,208,206,212,213,195,197,222,228,229,224,218,220,224,223,246,248,250,248,252,240,255,245,255,221,  7,  2,252, 15, 11,  2, 16, 15, 21, 19, 26, 31, 18, 17,232,250,240,240, 12, 14,231,237,201,197,201,202,205,207,206,205,196,208,196,198,199,190,197,198,255,255,255,255,253,254,255,254,  0,  0,  0,  1,  0,  2,  0,  0,222,236,224,222,205,224,223,222,247,249,247,247,248,244,245,244,
       52, 32,208,199, 48, 18,201,210, 43, 30,212,203, 48, 30,215,206, 49, 10,188,199, 79, 23,125,193, 48, 19,167,206, 45, 20,215,206,234, 31,192, 67,245,  8,190, 60,238, 18,193, 72,242, 30,193, 74,232,  4,194, 71, 11, 11,202, 25,236,  9,195, 64,245, 15,193, 63,194,194,194,193,192,196,194,192,193,195,194,193,193,196,194,193,194,192,192,191,194,195,196,194,192,195,192,192,194,195,191,192,  7,  7, 23,255,  4,  5, 17,255,250,252, 27, 31,251,251, 30, 14,198,199,197,192,183,202,198,179,232,233, 14,236,236,238,223,235,248,249,230,226,255,253,232,228,214,217,197,193,212,216,193,190, 67, 71, 76, 64, 63, 67, 73, 74,252,254,255, 26,244,247,250,244,192,190,195,192,192,192,189,192,194,193,192,198,192,192,191,194,191,193,191,193,193,192,191,191,202,199,191,191,200,201,198,199,
       13,166,212,188,  8,191,209,190,  4,142,178,181, 17,219,211,191,  6,209,208,209, 14,206,210,209, 12,204,216,189, 17,191,209,201,254,247,  3,240,254,213,  5,  2, 21,  8,255,  1,  1,252,  5, 10,253,254,250,253,  0,254,253,255,230,225,  5,229,229,  4, 21,  5,207,254,254,255,201,  0,253,  0,200,  1,  1,255,211,255,  0,  1,220,255,255,254,220,254,255,254,216,254,252,254,220,254,249,254, 16,252,240,  4,220,229,211,254,  3,  9,  2,  9,  8,  0,  3,  3,233,231,222,200,227,254,237,230,253,253,253,245,  1,255,248, 20,197,191,125,186,186,193,167,186,208,255,159, 21, 47,254, 45,244,240,208,  3,  1,227,218,  5,251, 49,221,  7,231,232,253,254, 23,  1,255, 11,252,  1,  0,  0,254,248,  0,252,  4,  1,  0,253,  1,  1,253,  1,243,255,253,255,253,254,  9,  0,255,255,249,  1,253,
      194,193,254,251,192,195,255,  4,194,195,  0,  2,194,218,  0,  1,193,197,255,230,196,194,  1,238,190,196,  0,  1,190,199,  0,  6,134,  3,245,179,236,195,248,199,184,193,248,218,197,192,248,192,203,197,251,211,125,198,188,198,206,208,251,222,137,219,190,191,  6,248,255,  0,255,  1,  0,255,  1,  0,255,  2,254,255,  1,254,255,  0,255,  0,255,  0,  1,  0,255,  0,255,  0,  1,255,  0,  0,229,  0,255,  2, 10,251,  8, 38,148, 19,240,227,228,245,224,229,254,254, 70,254,254,232,  5,240,250,252,248,  2,  3,220,  2,254,  1, 70,  0,  1,  6,  3,  6,229, 12, 71,  3, 62,210, 21,212,252,241,244,  1,224,242,247,235,229, 53,249,248,229, 35,  5, 37,253,252,  2,  0,  0,255,  0,  0,  3,  2,  1,251,  0,  0,  1,255,  2,254,  0,  3,255,  0,  2,255,  1,  0,255,  3,251,  0,  2,255,  1,
    },
    {
      214,220,226,231,223,220,224,230,226,224,225,223,218,215,213,231,204,220,214,170,198,209,204,164,198,201,190,153,195,198,185,182, 10,  2,191,205,  6,  1,194,201,  6,255,193,197,  4,255,191,192,  3,  1,195,191,  9,  3,198,191,  8,  6,202,192,  8,  3,215,197,  1,  0,246,254,  1,  0,246,254,  1,  0,247,252,  1,  0,249,252,  0,  0,249,250,  1,  0,250,249,  1,  0,250,248,  0,  0,250,250,226,244,219,216,211,215,218,220,239,236,231,250,230,220,235,248,199,242,209,215,218,216,210,215, 16, 58,151, 17, 68, 58, 67, 78,241,  9,236,232,208,219,230,233,234,215,195,232,243,195,237,203,227,247,163,214,213,217,195,194,  8, 12,255,255, 15,  7,  4,  9,  0,255,  1,  7, 12,250,253,255, 16,  9, 23, 12,  6, 32, 17, 12,204,201,196,196,191,192,195,200,255,  1,255,  1,  0,255,255,255,
       11, 22,  7, 10,239,251,  0,254,  5,215, 15,  0,237,  1,247,  3,254, 11, 17, 12,244,238, 65,244,231,248,253,245,  9, 11, 33, 11,195,193,174,171,196,199,  3,253,189,192,251,  3,194,193,  1,243,202,192,222,238,202,198,253,  2,195,193,254,251,195,190,252,254,195,197,192,192,196,191,199,194,192,192,191,193,192,189,196,196,192,192,191,192,190,192,200,197,194,192,198,192,193,194,190,193,  0,247,  1,255,  5,  0,  5,246,  0,  3,  1,  2,  4,254,248,255,254,252,  0,  2,  0,243,  6,233,196,192,191,184,191,191,205,185,  0,  2,255,  1,  6, 13,250,254, 79,105, 36, 58,114, 69, 45, 61,  1, 95,133,  4, 53,201,240,215, 18, 25,240, 27,255,253,238,225,  2,  3,  6,254,255,  2,251,249,221,220,208,241,209,207,204,229,232,211, 40,221, 11,205,242,220,230,193,222,245,203,221,220,  3,
      213,214,222,235,219,230,236,232,  2,231, 10,239,217,243,243,252, 11, 12, 12,248,210,250,239,231,254,246,195,239,248,244,238,  9, 82, 77,250,253,  1,  1,252,255,  0,255,255,  0,255,  2,245,253,243, 38,254,  0,254,249,  5,  4,253,  4,  1,255,  0,254,  4,  4,  1,  1,  0,  0,  1,  1,249,252,  0,  0,  0,  1,251,  0,249,255,254,  9,255,  0,  2,247,252,254,  0,252,255,255,  1,  2,253,  2,249, 15,254,254,  7,252,231, 17,  0,255,216,254,255,238,246,201,  0,  0,255,255,  1,  1,255,  3,253,  3,  0,  2,  1,255,249,  4,162, 14, 43, 97,221,132,192,247,194, 13,254,245,234,226,248,  6,226,205,127,114,  6,209,199, 14,189,184,226,193,184,197,218, 23,201, 43,  3,130,189,198,209,201,223,  5, 13,209, 56,243,238,244,190,  5,  7,251,232,124,195,191,175,128,  3,231, 21, 24,  1,228,
      253,234,233,242,158,249,253,214,  7,135,  0,248,220,243,  1,  0,246,252, 13,198,227,234,252,206,245, 25,194,188,253,236,  2,165,254,  0,253,254,  1,  0,  1,  2,255,254,255,  3,255, 16,  1,  0,255,112,  1,255,  2,  6,254,199,  0,247,  0,  2,  0,244,  5,  0,253,253,251,  1,255,  0,  1,  1,  0,  0,  0,  0,  1,  3,  0,  1,  1,  0,  0,  0,253,  0,  0,  6,254,255,  0,  0,  0,246,255,  2,253,194,251,  0,255,  1,  0,  2,253,  1,238,224,230,217,  4,245,  5,253,255,  0, 11,255,216,  0,206,255,  0,255,  0,252, 37,243,196,189, 63,125,  0,  1, 21,194,245, 66,254,  4,  7,  6,  7,252, 22,  5,144,141,  1, 89,218, 38,  0, 14,245,  1,245,  1,204,255,255,114,195,192,  5,  2,252,234,228,243,218,  0,195,  1,250,  1,  0,191,199, 57,193,231,255, 43, 58,122,250,  2,253,245,249,239,
    },
    {
      206,219,227,228,217,222,225,223,208,216,225,221,216,215,223,209,219,206,217,207,251,209,211, 19,179,197,216,154,202, 11, 17,202,202,250,202,234,201,247,207,238,200,250,202,237,199,246,183,238,196,247,188,243,196,253,222,  2,200, 10,170,  0,207,  1,195,240,221,  0,253,200,216,249,253,201,211,  2,254,199,191,252,254,199,192,  1,255,195,197,254,255,192,189,  0,  0,191,190,  1,  0,194,249,245, 10, 10,251,253,  4,  4,197,235,230,232,225,225,227,233,  6,232,231,  1,227,229, 18,  9,249,213,211,221,247,252,169,190,193,210,195,199,199,202,195,195,210,212,236,214,221,223,224,225,250,  4,255,  0,  7,  3,  4,  1,  4,229,  3, 19, 10,  1,251,  4,216,242,237,251,245,244,235,235,255,252,252,  2,250,254,255,  0,  2,  0,255,  0,251,254,  1,  0,  0,  4,  4,  0,  2,  0,  1,  1,
        0,194,238,248,254,195,236,248,  1,195,244,249,254,194,241,250,  0,193,238,251,  0,192,228,252,255,191,216,237,  0,191,196,251,253,250,243,  0,  1,250,207,  0,  9,  9,231,  4,227,  1,237,234, 10,  1,209,253,253, 12,203,226,174,198,201,249, 11,  1,174,209,249,194,218,225,188,204,231,217, 42,204, 31,219, 24,198, 31,234,  0,197,236,222,235,201,  1,214,  6,192, 10,190, 69,194,251,213,  1,237,250,255, 14,237,242,247, 10,  0, 13, 12,232,254,  9, 13, 29,218,129,235,164,216,221,203,132, 23,241, 50, 67,  7,196,239,255,254,  2,  3,248,247,  6,  3,226,247,227,222,233,216,196,229,  1,202,254,  3,  1,254,254,246, 21, 12,218,255,  3, 12,248,253,243,195,171,223,252,206,171,211,  3,  6,  0,243, 17,254,182,220,255,255,255,255,  0,255,  2,  2,255,236,  1,255,  2,250,  2,255,
        0,  0,  0,255,  7,  1,  3,255,  3,  1,  0,254,  5,255,  1,255,253,  4,  0,  4,  3,  0,  2,254,255,  0,  0,255,255,  1,255,250,205,247,191,243,226,250,  4,  3,235,229, 76,239,227,  2,243,208,252,  5, 51,213,250,254, 22,239,242,252, 19, 29,  2,  0,197,190,  1,255,  5,251,236,249,209,238,213,245,206, 15,  1,  1,251,254,252,252,167,255,248,250,254,249,182,249,216,  7,  2, 14,243,254,  3,  2,  2,249,250,  1, 34,  1, 12,  4,  7,241,253,248,227,223,164,227,  1,236, 13,232, 69, 16,  4,  2,130,212, 63,242, 17, 23,192,253, 64,255,254,  2, 34,254,192,203,201,204,228,236, 47,252,  3, 26, 10,254,233, 20,  1,249,168,166,137,192,191,191,  2,249, 62,168,186,230,  2,  9,235,248,190,210,192,191,192,193,190,238,128,198,125,192,192,194,191,191,195,187,194,216,238,207,192,194,
        2,  2,  2,254,  4,  0,  1,254,  0,  7,253,  0,  5,  0,  0,  1,  4,  0,  4,  2,254,253,253,  3,  0,  0,254,255,248,253,254,  0,245,237,250,  2,186, 49,  0,253,  4,  0,  3,226,188,233,194,  2,184,  8,246,  9, 55,235,242,  3,  7,244,219,217,176,227,255, 50,250,  7,  0,241, 77, 12,  4,255,239,253,201,255, 24,255,  1,249,236,  0, 15,250,215,  0,251,242,  0,253,  0,251,207,250,  1, 15,254,249,  7,245,192,241,255,240,196,255,  0,  2,  0,  0,  2,218, 22,255,234,  7, 67,  0,242,218,  0,230,255,222,250,201,250,196,  8,252,254,193,254,195,191,199,  1,  8,105,243,191,  2,201,204,193,  3,191,  3,221,254,245,247,245, 10,219, 26,183, 11,221,189,192, 91, 10,254,  0, 15,  0,247,  1,255,  2,  1,  3,254,246,254,255,250,  0,  1,  1,  0,230,  3,243,  2,  2,255,249,  2,254,203,
    },
  },
  {  //  21
    {
       70, 11, 24,  7, 70,  9, 28,  3, 69,  8, 26,  7, 69, 10, 25,  3, 65, 16, 18, 16, 66, 14, 27, 15, 68, 29, 10,192, 29,234, 11,234, 38, 56, 49, 56, 38, 55, 49, 58, 39, 55, 50, 60, 41, 55, 52, 59, 40, 54, 55, 63, 38, 54, 53, 62, 50, 44, 64, 51, 59, 39, 56, 41,  1,255,254,  3,  1,  0,255,  2,  0,255,  0,  2,  0,  0,  1,  3,  1,  0,  0,  2,  0,  0,  2,  3,  0,  1,  2,  1,  0,  0,  2,  0, 45, 39, 32, 40, 35, 32, 31, 35, 17, 19, 13, 18, 22, 22, 22, 27, 18, 21, 22, 20, 22, 83, 20, 27, 33, 19, 34, 32, 48, 63, 64, 73,226,226,230,224,229,232,229,220,  9, 19, 21, 18, 21,  8, 27, 18,251,  1,  2,  0,255,181,  3,  0,  8, 13, 17, 13,248, 78, 32, 12,  4,  4,  9,  5,  2,  2,  7,  1,253,254,252,254,248,  0,247,255,  2,  0,254,  0,  0,255,  0,255,255,  7,  1,  6,253,237,  8,  5,
       13,248, 80,231,235,220,  1,201,  8,  0, 80,234,  8,248, 89,234,  6,241, 63,234,  7,249, 88,238,  2,251,120,248,227,230, 66,203, 11, 13,251,  9,  8, 10,  3,248, 11, 12,251,  7,  7, 11,249,  4,  9, 14,  0, 10,  8, 15,252, 12,251, 26,252, 48, 18, 12,  4,251, 64, 63, 63, 62, 64, 64, 62, 62, 64, 65, 61, 62, 64, 65, 62, 62, 64, 64, 65, 62, 64, 61, 62, 62, 64, 63, 65, 64, 64, 64, 63, 62,254,  1,  3,  7,255,252,119, 25, 14,  9, 11, 12, 14, 13, 13,  6,253, 68, 67, 51, 63, 62, 54, 49,  2,246, 25,253,223, 11, 41,237, 26, 32, 35, 28, 22, 33,243, 22, 59, 67, 64, 61, 75, 63, 67, 52, 28, 72, 57, 60, 53, 57, 25, 62,182,254,255,  1,253,  2,250,  2, 60, 62, 63, 64, 62, 63, 64, 64, 10, 13, 16, 18, 18, 20, 20, 37, 63, 70, 73, 71, 70, 68, 66, 66, 66, 63, 59, 61, 63, 64, 64, 61,
      206,191, 58, 29,241,207,243,  0,  1, 12, 56,199,232,188, 59,238,188, 10, 66, 34,192,  8, 63, 11,192, 10, 64,211, 34, 44, 87,254,156,188,  2,  5,205, 56,  8,255,192,189,  2,247,194,193,  1,  0,200,193,254, 12,201,193,  1,  3,202,195,  2, 10,198,201,  3,  2,  2,  1,  0,  0,  1,115,  0,  0, 15,  1,  0,  0,252,250,  0,  0, 28,253,  0,  0, 20,  0,  0,  0, 29,255,  0,  0, 22,254,  0,  0,  0,  1,252,  5,  0, 13, 18,223,  2, 21, 55,  7,244,  4,204,  8, 98, 71, 58, 48,  9, 25, 65, 51,233, 35,  5,239, 75, 64,234, 17,  1,243, 37,253,221,230, 23,  8, 54, 66, 62, 65, 68, 72, 68, 55, 21,245, 19, 24,248,237,255,219, 92, 62, 62, 54, 63, 54,215, 69,  0,  1,  3, 10,  0,  2,  5,  2,  0,255,  2,  3,  1,  0,  3,  1,  0,246,  0,  2,254,253,  1,  2,  1,  1,  2,253,  4,  0,253,  0,
      136,203,  4,251,167,138, 49, 55, 66,252,254, 50,121,139, 56, 23, 67, 54, 18,  3, 60, 23, 16, 32, 68, 55,  6, 51, 85, 59,  2, 61, 87, 64,254,  2, 57, 69, 16,  1, 61, 72,  8, 11, 50, 62,  0,255, 64, 65,  5,  3, 62, 65,255,  1, 64, 62,  0,  2, 65, 67,255,255,  0,255,  1,  0,253,  2,255,  0,255,  1,  0,  2,130,255,  1,  0,  1,  0,  2,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  1,  1,160, 36, 12, 10, 57, 69, 38, 43, 67, 26, 31, 31, 29, 29, 18, 10, 37, 27, 43, 51, 78, 59, 24, 15, 24, 15,  4,247,221,224,236,237, 10,250,255,253,  5,252, 15,210,244, 16, 12,254, 20, 12,252, 36,219,231, 24,221, 57, 46,  0,  5, 10,  0,249,220, 47, 29,246,226,251,  3,  2,255,  1,255,254,253,254,  3,  5,  0,254,255,  0,  0, 10,255,  7,  0,  0,254,254,255,  1,253,  1,255,254,253,255,  0,
    },
    {
      235,251,232, 32,  1,  7,244,191,233,237,235,222, 26, 10,243,  9,234,  3, 14, 59,228,227, 17,232,  8,239,235,239,252,239, 15, 26, 12, 10,250, 20,  0,  7,240,219,  7, 11,  8, 19,  5,  4,245,249,254,  6,254,251, 13,  2, 15,243,  6,255,235, 69,  6,  7,251,143,  3,252,255,  1,  2,255,  0,254,  3,255,  8,  0,  1,255,255,  0,  1,  1,  0,255,  0,  0,  1,  0,  2,  0,252,  4,  2,254,  2,  0,226,  2, 20,254, 24,  9, 29,227, 28,242,253, 25, 14,231,  0,245,237, 74,177, 24, 26, 22, 57,194,222, 34,218,  2, 39,196, 46, 37, 16,249,  0,239,251,  0, 43,250, 48, 31, 34, 31, 38, 63,248, 52,221,207, 50,199, 54,215, 25, 32, 28, 53, 12,248, 13,  1, 42,  2,245,250,  0,  5,255,  0,255,  0, 20,  7, 27,250, 16,250, 15,241,251,249,252,253,  5,  8,  0,255,249,228,  2,  3,255,  0,  8,  3,
        5, 35, 33,204, 21, 49,245,238,253,  3,  4, 56,218,247,223, 79,242,203,219,224, 51,252,248,216, 40,151, 33, 85, 50,207,223,236, 19,240, 11,  0, 20,  2, 24,206,  6,233, 31, 11, 14,253, 26,246,  5,  4, 44,254,253,253,  4, 36,  1,246,197,252,  1,250, 85, 20,252,  6,  1,  1, 23,  1,  2,  0,  1,  1,  6,  2,  3,  0,  1,254,  0,  1,  1,  1,254,  0,  1,255,  1,  1,254,  1,  0,  0,254,254,228, 34,239,  4,251,  6, 20,248, 13, 11, 33,  5, 11,209,250,250, 16,228, 99,233,  6, 38, 52, 75,248, 36, 71, 44,236,194, 66,227,251, 27,246, 13,  4,  4,249, 40,195, 20,  7,  7,253,214, 10,246,251, 29,240,226,  0, 15,246,201,  7, 28, 42,  4,  5,210, 45,249, 17,254,252,  1,255,  2,  2,253,228,198,  6,  1,250,211,  2,  1,  6,  8,  1,  0,  1,  2,  0,  0, 18, 24,  1,229,  2,255,250,  8,
      249, 23,  7,248,245,240,  1,239, 59, 28, 59,225,229,214,253,193, 13,234,249,246,  6, 58,251,  5, 28,159,237,202, 19,165,182,205,240, 17, 66, 45,  1,  4, 61, 47,  0,249,  2, 20,  5,255,234, 88, 16,254,188,241,  0,  0,  0,247,255,253, 48, 30,252,  4, 13,198,254,254, 19,  8,  2,  3, 40,  0,  1,252, 11,254,  0,  1,  2,254,  1,  0,255,255,  1,  0,241,  0,255,  0,254,  0,  0,  0,  0,  0,252,  5, 35, 10,230,  3,  3,252, 15,  6,251, 19,247,  6,254,236, 18,249, 44,237,255,  0, 60,184,  1,222,250,255, 44,249, 14, 72, 10,228, 19,249,250,251,  9,255, 54, 19,254,238, 15,250,255, 92, 35,225,201,212, 69,187,  4, 12,209,  4,220,255,  6,197, 92, 47,254, 14,253,255,255,253,  3,  0,237, 44,  4,  9, 38,  1,192,255,248, 11,251,252,248,253,  0,250,227,233, 12,225,  0,217,254,248,
       74,220,208,233,211,205,229,255,212,238,182,  2, 29, 39,  6,253, 41,197, 47,  5,247,177,255,  0,251, 14,227, 65,223, 87, 57, 51,250,  2,194,214,  9,  6,  0, 13,  2, 19,231,250,251,255, 22, 36,  0,250, 64,  3,  3,  0, 65,  0,  3,254,183,219,  8,  2,  2,209,  0,254,  0,246,254,  0,  0,241,  2,255,242,  0,253,  0,254,  1,  1,  1,  1,  1, 27,255,  3,255,  1,  0,  1,  0,  0,  0,  0,  2,249,248,253, 66,250,  4,190,248,233, 21,  6,252,253,244,212, 22,238, 14,248, 76,254,138,  0,179,  5,  8, 43, 33,202,191,134, 12,253,240,254, 19,255,  3, 57,  9,239, 65,253, 17,252,255,  6,195,253, 69, 56, 33,  7,108,195,183,234,234,169,217,  0, 61,249,218,244,  1,  1, 25,  7,  1,223,253, 19,229,  4,  1,197,249,248,  2,  2, 10,  3,  4, 10,  1,183,253,245,245,216,248,255,  4,  0,254,
    },
    {
       13,  6, 54, 89, 20,245,253,197, 63, 46, 13,247,200,240, 96,199,  7,  8,  5,234, 34,248,185,216,223, 54,212,254,  9, 50, 53, 64,  9,230,  4,245,  1,251, 30, 22, 12,  1,245,255, 11,  0,251,251,  6,  8,126,242,243,  9, 83,252,234,246,249,254,  2, 42,253,  0,  3,  1,  0,255,  1,  1,254,255,255,  1,  0,255,253,  3,  0,255,254,  3,255,  2,  0,253,  0,  0,  3,  5,  0,  1,  2,  0,254,  0, 29,220, 42,208, 13,234, 11, 40,234,243,220,  1, 12,237, 38,205, 32, 16,218, 56,  6,  4,  3,242,  7,204,250, 21, 10,184,190,219,255,253, 29, 60,  3,240,255,212,254,  4,  6,247, 33,  0, 69,255,254,  2, 48,  5,221,177,  4,  6, 34,241, 54,238, 34,251,235,101,  0,  0,  7,  3,  1,248,  7,  8,  4,251,254,252,  1, 17,  3,  1,251,  0,245,  3, 20,  2,252,251,238,250,  5,  0,253,247,253,  0,
       12,  7,  4,166,246,247,  9,234,198,168, 17, 11, 47,255,200, 53,  9, 18, 35, 77,238, 15,  1, 67, 11,  0, 29, 63,251, 16,142, 54,  2, 31, 86,250,  2, 13, 14,236,  1, 24,  6, 16,255,239,251, 18, 10,247,150,213,  5, 15,249,185, 26, 50,  0,  0,254,243,  4, 18,  1,  3,254,  0,255,  3,255,255,  3,  9,  0,  1,  1,253,  0,  1,  0,  0,  0,  2,255,  1,  1,247,255,252,  0,255,  1,  2,  0,252,226, 23,216,249, 47,  0,255, 84,183,210,234,173, 47,  3,251,137, 18,200, 30,178, 17,220,246,210, 23,254, 24, 19,  0,191, 27,250, 17,230,246, 29,205,200,  1, 76,220,242,  9,173,187, 76,186, 99,220, 10,215,  8,  7, 52,253,255,182,  2,198, 44, 49,  0,  8,253,  2,  6, 11,  5,  6,  3,  0,247,254,  2,251, 13,  0,  5,254,251,241,  8,  0,  0, 12,250,  1,255,  4,244,  5,252,  1, 28,250,  0,
      251,246, 18,  0,  3, 31,  0, 19, 22, 39,246,248,  2, 95, 18,196, 23, 13,223,198,222,237,217,  7,249,220, 69,220,230,162,  0, 69,242,227, 15,252,  6,  5,  6,247, 16,249,  5,251,253, 17, 28,244,206,236,221,241, 11,242, 54, 55,250,230,  1,  0, 10,249, 13,  0,255,252,254,254,  1,  3,  1,  1,  1,253,255,255,  1, 11,  0,254,  1,255,  4,  2,  0,  2,  2,251,  1,  1,  0,255,  0,  0,252,  1,  2, 12,246, 45,  5,  5,  0,107,  1, 67, 84, 60,223, 66,212,  6,251,108,  3, 75,193,220,233,114, 19,  5,110, 31, 45,216, 54,  7,222, 42,226,202, 56,221,251, 92,238,238,235,217, 50,232,139,151, 10,213,  4,245, 67, 20,  5,251, 59,220,253, 47, 90,247,  0,249,253, 10,237, 10,  2,187,  1,  2,247,229,  0,  2,  8, 11, 19,246,254,252,255,255,  0,  5,  1,252,241,244,250, 13,243,  0,  4,  7,
      254,254,253,  4,250, 24, 20, 61, 16,  4,254,213,  5,206,  8, 86, 42, 21, 64, 79, 41, 39, 59, 24,254, 10,189,187,246,231, 20, 83,  2, 23,172,255,250,  2, 17,  6,  9, 14,253,  0,254, 22,247,  1, 52, 21,254,195,  4,  1,  7, 22,249,  3,  1,185,253,254,  0,  0,253,255,254,255,  2,  1,255,  1,251,253,255,  0,  1,  3,253,  0,245,  1,  1,  0,  0,255,255,  2,  1,255,  0,255,  1,  2,  0,  2,249, 13, 76,235,  9,251, 17,  1,228,233, 34,  8,197,173,251,248,193,194,229,251,250,230,169, 57,222,219,188,149,  9, 68, 41,209, 18, 13,242,229,208,189,230,237, 17,255,206, 16, 29,235, 62,102,251,252,  6,170,  1, 26,254,  6, 24, 22,  3, 44, 17,231,253,  8,255, 17,251,  1, 15, 10,253,  3,  9, 26,249,  1,200,252,239,  9, 18,251,246,224,255,191,  4,  0,  3,  8,253,253,243,  0,  0, 16,
    },
  },
  {  //  22
    {
        1, 91, 26,233,232, 63,247, 63,178,249, 87, 18,  5, 69, 61, 16,212, 10,209,143,243,105, 71,230,206, 14,211, 35, 62,174, 88,113,255, 10, 50,183,251, 12,  5, 19,  8,  0, 48,251,255, 15,217,254, 16,  1, 26,190,253, 15,187,191, 11,254, 20,  7,253,250,209,247,  0,  0,139,  2,  2,  0,  0,  0,  0,  4,  1,255,  0,  2,  0,  0,  1,  0,255, 14,  1,  2,  5,  6,  1,  2,  1,254,  1,  2,  5,252,244,  4,  0,253,254,254,  1,  2,247, 16,  3,247,253,  4,  3,252,220,127,  7,245, 56,  6,  8,249,212,192, 60,206,192,252,160,186,236,  4,205,  5,215,115, 63,246, 23,251,  0,  0,  3,240,241, 11,252,  2,255,  0,  4,254,  1,  0, 50,255, 14,  9,244,255,245,  0, 28,255, 11,  1,  0,  0, 12,  0,255,  0,  0,  0,  0,  0,  2,  3,255,255,  0,  0,  0,  0,  1,255,  3, 24,253,  0,  0,  2,252,  1,
        5, 79, 91,196, 20, 16, 13,176, 43,246,203,  1,252,  7, 83,246,  4,  4, 15, 52,220,186,207,203, 85,219,  0,202,240,226,209,206,252,  6,133,239,255,239, 73, 53,  1,253,241,  3,  2,  8, 54,195,255,254, 63, 75, 24,243, 13,147,252,  7,248, 76,253,254,196,171,  0,  0,255,  1,  0,  0,226,  0,252,  0,176,  1,  0,  0, 25,255,  0,  0,  0,  2,  1,255,254,  7,  2,  1,  5,245,  0,  0,254,254, 14,252,255,253,  1,250,  4,  0,248,  5,  0,  8,255,  0,  0,  0,  0,182,255,242,252,240,  3,195,188,251,187, 85,144,122,  2,248, 16, 61, 63, 88,206,201,  4,222,247, 10, 13,  6,241, 30,173,  6,  0,201,  0,249,  1,240,  9,  3,192,210,192,189, 69,244,196,244,232,255,  1,241,248,254,  0,252,  1, 10,  0,  0,  1,  6,255,  0,254,  2,255,255,  0,250,  0,255,243,223,238,254,  6,254,  3,  2,
      248,  6, 57,242, 26,  1,165,202, 55,226, 13,255, 14,  7,137,183,251,200,  3,187, 62,249, 69,254,239,170,186, 65,  1,181,243,  1,254,242, 41, 59,255, 37,122,  7,  9, 29,182,142,252,  3,124,124,255, 68,188,  2,255, 37,253,197,252, 18,113,126,  1,  3,  7, 55,255,  0,249, 27,  0,  0,243,255,252,  0,  1,254,  0,  0,  3,  3,  0,  0,255,  1,255,  1,243,  0,  0,  0,249,  0,  0,  0,229,255,  1,  7,254,251,  2,  2,253,  1,  5, 14,  0,243,  0,  8,255,  4, 50,195, 63,242, 25,146,  0,195, 26,213,224,  9,109,215,243, 47, 40,204, 63, 46,190,253,105, 63,253,254,  4, 24, 31,232, 15,163,245, 64, 46, 37, 83, 59,  1, 13,  0, 22,108, 61,198,  9, 69,151,  0,255,  0,  0,  0,255,  2,255,249,  1,  1,  0,250,  0,  1,  0,  2,  0,  0,255,  1, 11,  1,  0, 11,  6,252,255,153,  1,  2,  0,
        5,  3,  1,253,  0,214, 34, 34,  1,247,251,  6,252,201, 60, 58, 87, 55,201,  0,241,241,198,  1, 65, 61, 28,255, 87, 59, 61,255, 15,  9,252,221,  1,253,195, 70,  2,242,240, 65, 11,  2,187, 44,248,  4, 68,250,251,242, 64, 30, 73,  1, 63,241,242,  0, 58,  3,  7,  0,254,  1,  1,  0,  1,  0,255,  7, 17,255,  3,  0,  0,255,  2,  1,  2,  0,255,  0,  1,  0,254,255,  0,  1,254,  1,  1,  0, 20,  6,  1,  0,254,  7,  0,  9,  3,249,  2,  1,247,251, 13,  6,244,235,244,186,251,232,207,147, 76,  6,221,120,255,250,  3,  1,195, 37,239,164, 13, 54,127, 51,  9,  9,240,249, 15,  1,238, 70,250, 60,167,  2, 11,216,184,171, 11,222,126,210,254, 68, 18,215,255,253,  0,252,255,  1,254,255,  0,  1,255,  0,  2,  1,  1,  1,255,  2,  0,  0,  1,  0,  0,  2,  0,251,  0,  0,  2,  5,  1, 10,
    },
    {
      223,254,219,229,  8,205,251,222,223,234,199,250,242,233,242,248,211,233,173,192,249, 26,105,  1,214,243,219,214, 23, 11, 17, 19, 18,254,227,  4, 10,248, 10,208, 14,  3,253,255, 12,  6,237,215,  0,  0,  7,245,  5,  7, 12,  3, 28,  5,  0,248, 10,  7,251,244,  2,  0,  0,  0,  1,  0,  2,  2,  2,  0,  2,255,  2,  0,  0,  0,  1,  1,  0,255,254,255,255,255,  2,  3,  0,253,  3,  0,  0,255,253, 11, 14, 40, 21,  0,251, 29, 35,  7,245, 25, 35,  0,216,  1, 11,247,213, 82,216, 46, 91, 61,  3,244, 66,  9, 11,  4, 43,  9, 21, 17,236,249, 13, 29,250,  1,236, 44,221, 80,  4, 13,141, 41,239,241, 53, 40,248, 37,252,185,223,249, 32,  4,247,  6,254,  4,  3,  3,253,  2,  4,253,  1,  1,246,253,  7,227,250,233,  0,251,247,255, 15,254,  4,  3,  1, 50,232,243, 20,  0,254,254,254,  1,
       49, 38,222, 45,  6, 28,212,225, 51, 28,237,242,  6,231,226,  9, 48,224,115,  5,251,  0,157, 40,214, 14,  3,249, 20, 65, 19,210,  2,  1,224, 44,246,  1,255, 15,254,255,254,254,250,  3,208, 41,  2,  3,245, 28,  0,248,249,250, 14,  5,  7, 10,  0, 13, 68, 15,  0,  1,254,225,  0,  0,  2,254,255,  0,  1,  0,  1,  1,254,254,  1,  0,  0,251,  1,  0,  2,  0,  0,  0,  0,254,  0,  1,251,  0,  4, 72,215, 37,213,184,210,  1, 47,246, 18, 10, 29, 10,  3,253,  0, 24,237,100, 29,130,139, 40, 56, 23,255,242,  5, 38, 58,253,248, 47, 19, 78,170,237,201,191, 49, 16, 62, 44, 53,206,183, 21, 29,  5,181,172,  4, 25,211, 55, 43,223,253,254,  1,244,255, 42,248,213, 19,255,230,254,  4,255,  3,253,  3,250, 26, 60,246,  0,230,236,  2,  2,  6,254,255,  0, 22, 20,254,  6,  0,  0,  0,255,
        0, 49, 28,136, 11, 37,186, 39,217,  3, 59,  0,245,255, 75,201,247, 85,  9, 28, 33,192, 56,213,235,179, 96, 35,218,189,253, 57,  0,255, 70,255,254,248, 32, 16,  0, 41,195,  6,  2,  1, 53, 35,253,  4,253,239,  2, 20,255,251,  0,253, 17,202, 15,  0,192,249,  0,  2,  7,  0,  0,  0,  4,  0,  0,255,243,  4,  0,  2,  1,  0,  2,  0,254,  0,  2,  0,  1,  3,  0,  0,  1,253,  0,  1,  0,  0, 45,  3,251,222,255,  0,247,  3, 37,202, 20,252,222, 13,  1,254, 13, 52, 54,216,116, 60, 85,186, 22,243,  7,245,209, 29,190,208,246,225,254,219, 41, 49,  8,191,  6,  5,217, 14,184, 57,187,251, 35,206, 31,  8,254, 64,  4,  5,171,227,183,253, 73, 14,199, 77, 13, 58,  3,251, 22,254,  3,  0,249,240,  9,240,  5,  1,  6,  8, 31, 38,  1,251,  7, 45,  0,  0,251,211,  5,254, 58,255,  5,255,
      251,198, 29, 62,249,217, 62,225,230, 31,253,  1, 16, 16,240,  9, 36,239,  4,215,241,225,138,248,239, 31,191,254,  5,253,251,239,246,  2, 67,198,  5,  3,254, 34,  2,223, 58,  1,  7,  7,230,212,  0,255,  3, 10, 10,  1,251, 61,  3,  4,  0, 17,  8,  0,252,212,252,253,  2,  2,  1,  1,254,  0,  1,  1,  1,  0,  1,  1,  0,  2,  1,  0,251,255,255,  0,  0,  1,  1,  1,  0,  0,255,  1,  0,255,235,197,213, 20,252,197,255, 15,201,252, 24, 12,  3,232, 53, 65,  1,203, 60, 29,  1, 31,252, 77,253,233, 94,249,  3, 44, 65, 10,  2,  2,226,196,  3,231, 26, 26,247, 36,252,200, 20,  3, 27,245,207,  3,246,221,190, 22,246,  2, 18, 39,115,  0,235,195,  4,206,  2,  8,236,  3,  0,  3,254,  0,  8, 38,  0,240,253,255,243,235, 28,241,  1,  5,252,255,255,254,246, 44,245,  0,  5,  0,  2,  6,
    },
    {
       19, 13, 52,  1, 22, 19, 52, 42,231, 18,  4, 11, 34, 48,251,229, 12,248,  8,251, 23, 23, 17, 10, 19,  1, 33,247,237,251,  5,239,  1,  7,253,249,  2,  7,  3,247,251,255,253,244,  6, 13,243,233,  4,  9,242,246,241, 20,255,249,251,  8,239,252,252, 25,247,238,  0,  4,254,254,  0,  3,255,  0,  3,  0,253,255,  2,  2,  0,  1,  1,  1,255,255,  2,  1,255,  0,  4,  2,255,254,  4,  3,  0,  0, 13,255,226, 63,246,211,  3, 14,239,232, 69,239,241,  3,219,245,252, 54,242,194,211,189, 35,244,229,242,234, 15,234, 10,175,  4,  2,238,255,249, 24,252,246, 18,  4, 19, 18,249,253, 26,  9,245, 68, 41,224,  1,224, 20,248,108, 27,253,  4, 58, 76,205, 11, 35,  0,  2,250,253,  2,  0,  3,251,  1,  3,250,252,  0,  0, 18,  0,231,255,  6,  2,252,247,  3,  0,  7,255,  4,255,241,251, 11,  2,
      240,250,168, 22,242,176,116,227,246,  9,  0, 15, 40,242, 15,  5, 44, 65,200,  6, 29,253, 16,  4, 13, 23,248,172,252,221,  7,195,  1, 13,186, 15,  1,247,221,  9,  6,  3,  9, 17,248,242,251,246,249,233,  6,  2, 26,  1,196,244, 30, 13,251,  1,243,240,  1,247,  0,  0,253,255,  0,  1,  2,  0,255,  1,255,  0,254,254,253,255,  1,253,  1,  1,  0,255,  3,  1,  0,  1,  0,  1,  3,  1,  0,254,227,252,  3, 22,226, 81, 50,192, 50, 46,180, 59,190,  1,204,222, 20, 75,220,  9, 11,255, 75,181,255,244,  2,  7, 96,216, 68, 71,  3,  9,250,249,  9, 51,249,245, 14, 45,248, 76, 60,247, 10,242,222, 11, 48,215, 11,  2, 63,  5,243,  1,138, 17,110, 44, 14,  4,  0,  0,251,253,254,249,255,  0,  3,213,254,  0,  2,254,  0,  2,228,  4,252,206,  8, 60,247,  4,  2,255, 43,254,200,255,  1,  0,
        3,230, 57, 49,255, 94,200, 29, 39, 10, 17,253, 47,254,250,232,  3,  4, 65, 23, 37,238, 32,  6,247,255,252, 49,202,237,  3, 79,  1,236,250, 87,  0,  3, 42,244,246, 13, 11, 13,  1,253,194, 93, 36,  3,255, 17,255,  2,248,249,249,208,183,  5, 45,238,192,236,  0,255,  1,255,  0,  2, 11,254,255,  0,  1,  1,  0,  0,254,  1,  1,255,  0,  2,255,  1,255,  1,  0,255, 11,255,255,255,255,  0, 61,255, 69,201,254,236,  0,  6,215,194,249, 63,254,  6, 58,234,  8,172,250,208, 17,206,200,  2,  5, 51,  6,  3,126,  1,  6,183, 74, 24,189,244,223, 24,  2,  2,  1,245,255,251, 10,248, 66,237, 42,251,226, 13,249,253,  3, 44,255,  8, 82,  4,  0,237, 64,252, 56, 17,  3,255, 31,234,254,  2,  2,230,247,246,  0,  1, 54,  1, 44,253,237, 40,255,205,  8,  5,253, 21,  9, 20,255,  1,  0,  5,
      242,245, 32,183,255,220, 46,242, 33,217,227, 21,173, 37,  4,254,250,235,245,  4, 13,215,  3,  7,255, 32,  3,  0,  0,242,251,190,239,  4, 73,  7,  1,  5,  7, 39,  0,255, 22,242,  2,  3, 70,179,254,254,255,  7,  0,  7,250,242,  0, 52,194,252, 38, 30, 68, 28,  0,  1,  2,255,  0,  0,  0,  1,  0,  2,  1,  1,255,  1,  0,  0,255,255,  0,254,  1,  1,  0,255,  0,254,  0,255,  4,  0,254,  2,254,236,254, 38,190,195, 56, 87,214,231,223,184, 59,231, 57,255,238, 47, 22, 43,254, 45, 11,243,  0,225,249,231,  8,253,138, 68,227,226, 76, 32, 39,254,246,229,  5, 19,  9, 29,251, 51,  6, 17,195,241, 16, 21, 41,209,226, 26,254, 10, 12,217,237,  0,  0, 97,206,253,215, 23, 10,  6,236,  4,255, 46,  8,  4,  2,232,  0,252,224,  2, 13,230,  9,255,  2,  4,  0, 10, 16,252,252,  4,  4,246,
    },
  },
  {  //  23
    {
       16, 11, 43,237,  0, 13, 38, 18, 11,251,179, 82, 66,245,165,247, 35, 65, 54,239,231,255,  8, 12,223, 24,249,207,  8,  6, 10,127,231,  3,252,243,224,  9,246,252,253,  1,252,  2,  5,  9,211,  5,  1, 13,187,255, 24,  3,246, 60, 48,  5,255,232,220,210,  9,192,254,  0,255,250,253,  2,255,253,  0,254,  0,  3,255,  0,  2,255,  0,  1,  4,  1,  1,  0,255,250,253,  2,  1,  4, 13,  5,  0,  0, 10,  1, 12,  1, 16, 10, 20, 71,  7,253, 14, 42, 62, 46,194,197, 45, 51,239,249,194,  5,247, 48, 24,194, 10, 29,254,224,193, 54, 11, 12,  3,  0,255,252,250, 11,241,206,249,254,244,101,  5,  3,223, 66, 64, 65, 72, 69, 67,225,  2, 13, 56,250, 13, 15,  7, 35,  0,  0,  0,  1,255,  1,  0,254, 47,245,  1,255,  0,  0,  0,  1,255,  5,253,249,255,241,253,241,254,  0,254,252,254,  1,  1,  4,
        6,252,163,190, 21,243, 25,121, 25, 18,126,198,  2, 14, 30, 14,252, 21, 75,253,250,  7,  5,194,253,227, 11,201, 61,  8, 54,245, 22,  6, 80, 63,248,  3,174,247,237,  6,  2,112,225,  2,253,242, 52,236,233,205,  0,  1,254, 75,217,185,  0,233, 70,173, 18,  0,  4,  2,  1,239,  0,255,  1,  1,  1,  0,  0,251,  2,255,  0,  1,  1,  6,  0, 23,  0,  1,  1,  5,  7,252,253,  2,  2,  1,  0,  1,252,  0,  2,  5,243,  5,253,190, 36,247,  0, 84, 48,  9,193,131,  0,254, 26,182, 68, 27,  7,253, 73, 31,223,183,  0,  1, 53,  4,  2,  7, 12,  1,255,  3,  1, 10, 32, 49, 15,247,253,162,  0,255, 91,254, 56,  1,254, 59,212,251,254,  6, 34,252, 34,199, 33,  4,  0,  0,  0,255,  0,  1,  0,  2,224, 10,253,  0,  0,  4,  1,  0,  0,253,255,254,254,  1,  3,252,248,254,  2,225,254,255,240,  4,
       65, 10, 62, 62,240,  9,203,184,228,254,226, 54,192,210, 57,124,  1,252,126,  1, 41,  0,239, 64,211,  1, 69, 67,250, 67, 55, 41,228,255, 35,200,183,252, 69, 15, 32,255,113, 63, 11, 11,255, 56, 58,196,242, 21, 27,  0, 55,230, 15,238,188,255,217, 65, 17, 21,254,249,  0,  1,  0,  0,  1, 14,  1,255,253, 12,  1,  0,  0,  2,  3,160,  2, 24,254,  0,  1,  4,254,250,255,255,  0,  0,  0,254,  6,  1,242,180,180,  8,244,  8, 32, 20, 72,184,204, 39,191,251, 14, 63,  7,108,247,202, 20,244, 77,  2, 45, 52,129, 56,191,138,255,252,  6, 13,  2,  0,  0,  0,241, 19,  1, 65,254, 17, 22,254,219,192,254,121,185, 50, 37, 65,244,220,245, 10,  5, 18,  4,  8,  0,  0,  0,252,255,  0,  1,  0,  2,255,  9,  1,254,254,255,  1,  4, 18,222, 17,  1, 13,253,254,  5,  1,251,255,253,  1,255,248,
      250,252,252,  2,238,247, 34,198,244,215, 55,191, 17,255,255,194, 69, 27,244,208,200,252, 70,178, 72,  3,219, 56,250,250,198,245,  3,255,118,234,115,255,194, 64,  5,254, 56,253,243,  1, 23,186, 30,  8,  6,238,251,  0,122,255,238, 26,208, 16, 27,  1,  8,219,  0,  7,255,  1,  2,255,  0,  1,  1,  0,255,  0,  2,  1,255,  0,254,253,  0,  3,136,254,  2,254,  0,  2,  0,  0,252,250,254,  2, 14,250,232, 49,248,253,246,196,231, 31,  2,202,246, 25,  2,196, 29,197,190, 62, 65,243,195,100,252, 13,228,233,234,182,  9,169, 12,  9,250, 84, 11, 19, 16, 19,251,253,255, 10,  3, 41, 10, 31, 62, 13,220,  5,  9, 73,159,220,237,243, 57,228,  5, 60, 55, 64,  0,  0,  3,  1,255,  1,  7,255,  0,  2,  5,  4,  0,254,  4,  3,  0,238,239,252,  0,  6,  3,249,250,255,  2,  3, 12,254,  5,  3,
    },
    {
       18, 36, 41, 44,248, 66, 51, 46,  1, 11,251,251,226,246,217,219,244,  2,233,245,248, 10,243,249, 13,  5, 10, 18, 19, 77, 34,255,249,  4,244,246,  0, 12,  5,246, 14,  7,241,246, 18,  9,244,247, 13,252,241,245, 13, 13,249,250,241, 10,248,253,212, 27,249,185,255,  4,255,  0,  0,  2,255,  0,  2,  1,255,  0,  2,  2,  0,  0,  2,  1,  0,  0,  3,  2,  0,255,  5,  0,  0,255,254,  0,  0,252, 15,177,238,246,213,  1, 45, 16, 26, 17, 32, 30, 16, 53, 20, 56, 69, 15, 14, 21,246,199,255,200,214,247,  1, 53, 48,  3,223,252, 43,249, 13,212,207, 12,204, 29,236,249,245, 29,224, 10, 25, 14,194, 43,228,197,  8, 35,255,218,253, 21,215, 31,227,  0, 11,  3,  2,  5, 35,  2,246,  1,240,  1,244, 11,  8,  5,  5,  1,248,  4,255, 13,  8, 33,  4,251,  0,  3,  1,  6,  4,  6, 18,  0,234,254,
       28,231,156, 35, 10,207,255, 13, 30,  6, 38, 44, 57, 49, 20, 15,237,213, 20,254,186,239, 20, 27,255,250,171, 57,248,  7,  1,  8, 12,251,249,247, 21,  5,245,  2,  7,253,  6, 10, 14,252,  7,  4, 11,249,  9,  1,  9,  9,  9,  9,234,250, 10,  2, 24,254,248, 69,  1,  0,  0,255,  0,  3,  0,  0,254,255,  0,  0,255,  0,  0,  0,254,  1,  0,  0,252,255,255,  0,  0,255,  1,  0,  2,  0,  0,  0, 25, 32, 25, 35,203,229, 24, 14,245,214,254, 19,220,  4, 18, 24, 76,193, 19, 47, 60,  6,  5,253, 41,  0, 25,234,196,175, 66,216,221,221,225,240, 37,  5,213,253, 26,201, 27,174, 36, 43, 13,246,237, 19, 29,241,221,175, 12, 56, 28, 36, 46,253,186,183, 89,192,254,  2,  8,254, 16,  1,242,  0,  0,251, 37, 33,  2,  2,252,  7,  3, 33,253,231,251, 36, 11,  0,254, 27,248,214,187,238, 44,255,
       17, 20, 27, 54, 35, 23,249,183,239,  1,183,220,254,253,226,175, 59, 44,253,253, 28, 12,211,208,236, 17, 96,190,243,234, 40,  8,231,  1, 55, 63,253,251,  0,  3,227,249,255,251,232,254,  4,  2, 14,240,  5,253,244,  2,251,255,  5,  2,255,  9,252,255,  1,  0,  0,  0,  4,237,  1,  0,255,253,  0,  0,255,  0,255,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  0,255,255,  0,  0,255,255,  0,  0,  4,  7, 10,254,246,  3,254,  0,222, 24, 23, 20,246,233,158,138, 52,240,252,184,176,  8,250,196, 27, 29,231, 33,159,244,252, 17,254,  6,255,236,  0,  5,  1,  2, 15, 54, 23, 38,  8, 31,  0,191,240,175,  2,230, 27,182,199, 83,210, 21,251,217, 94,  1,  4, 68, 28,254,  6,  2,255,  1,  1,255,  1,248,206, 49,236, 49,183,243,228,235, 10, 13,164,254,  8,  0,  0,188, 59,249, 11, 51,248,  5,
      203,  4,243,193,250,  4,  2,  0,211,240, 54,179,252,  5,215,252,  6, 33,226,199,216,245,183,  6,216,229,  0,  0, 25, 45,247,  8,  9,  8, 71, 68,255,255,  0,255, 38,255,  6,  2,253,  3,  5,253, 27,  1,  1,255,  4, 10,254,255,223,254,  5,250,  1,255,  3,  0,  0,254,  1,254,  0,  0,255,  0,  2,  0,  0,  0,  2,  0,  1,252,  0,  1,  0,  1,253,  0,  1,255,  0,  0,  0,  0,  0,  0,  0,  0,239,248, 24,  1,206,254,  0, 12, 48, 19,  8,247, 64, 92, 75,244,249, 14,251, 25,232,178, 22,245,253, 50, 24, 17,202, 66,227, 36,243,  4,  7,229,  4,253, 10,  5,247, 25,  6,215,  6,186,193,112, 24, 33,242, 70, 57, 55, 38,106, 53,238,  0, 56,201,102, 16,200,247, 13,240,236, 69,  4,254,  1,242,242,  4,247,  2,250, 44,209,  2, 24,245,  4,  2,198,  4,  4,  0, 45,234, 81,  8, 22, 66, 18,
    },
    {
      251,243,250, 12,  2,  4, 21, 50,  0,245, 16,  0,  5,248,  2,255,230,249,  6,252,248,244, 10,119, 25,253, 24,  1,  2,  5,  6,134,  7,249,242, 45,255,  0,251,236,255,  1,245,  6,  5,255,248, 55,249,255,248,248, 13, 16,  4, 13, 35, 10,242,252, 38,  7, 68,  5,  1,255,  0,252,  0,  2,  0,  1,  0,  2,255,253,  1,  0,255,251,255,  0,255,255,  3,  2,  0,255,  2,  3,  0,254,  0,  1,  8,  2, 10, 17, 30,233, 16,212,217,226, 49,  6,204,237,181,  1, 33,167,205,255,169,  0,201,214,214, 56, 42, 19,248, 38,  3,218, 80,255,  3, 20,  1,  5,  2,  0, 31, 11, 13, 27, 90, 21,232,254,  0,  9, 20,188, 54, 44, 11, 21, 30, 30,253,196,  7, 20,  0,  7,197,155,  4,255,253,252,255,255,254,255,  1,  0,250,  0,  5,254,255,  0,254,232,255, 12,  1,  1,  2,255,  7,  9,235,115,  3,  2, 10,  1,
      238, 65, 61, 23,212, 16,198,142,225, 13, 67,245, 36,  5,255, 25, 10, 10, 51, 16, 46, 62,247,229,249, 33, 44,181, 27, 45,249,145,253, 30,157,177,234,252,227, 66,248,255,193, 57,255,  1,243, 27,251,254, 49,195,  9,242, 63, 68, 12,  0, 47, 50,  9,  3,134,  4,  0,254,  1,252,  0,  0,  0,245,  0,255,  2,245,  0,  0,255,254,  0,255,  4, 10,  1,  0, 16,  0,  0,  0,  0,  2,  0,  1,  1,  0, 17, 26,190,  6, 34,212, 23, 36,208,  9,  2,240,  0,209,248,253,  1, 75,195, 46,235, 94, 77,218, 11,251,246,  4, 13, 32, 81,  2,  2,249,231,  1, 17,252, 10,  6, 20,247,251,  8,252,235, 48,  2,237, 47,228, 75, 34, 65,252, 46, 26, 66,  4,  4,129, 21,247,  0,252,250,254,255,  0,  3,255,  0,255,  2,254,254,255,255,  0,  0,210, 14,255,  2,128,126,  0,  5,254,  3, 11, 79,253,255,255,  0,
       30,238, 27,  0, 43,238,222,  2, 15,248, 17, 23, 40,244,245,240,229,164,253,253, 18, 50,228,211,226, 51,141,251, 22, 68,197, 13,  7,237, 87, 31,254,244,  3,235,  4,255,193,197,238,  0,  6, 10, 57,  1,192, 25,231,245,194,242,244,253,242,192,249,245, 40, 71,  0,  1,  1,  2,  0,  0,  0,  6,  2,254,  0,  0,  1,  0,  0,  0,  5,255,  0, 12,254,  0,  1,247,  0,255,  0,  0,  0,  0,  0,  0,240,233, 62,243,215, 36, 41,204, 52, 49,201, 62, 37, 45,214,218, 45, 11,199, 63, 11,171, 16, 45,  4,183,243,254, 36,  3,219, 67,  5, 60,  7,254,230,  1,173,254,  1, 13,241,  1, 17,199,236,191,  9,  8, 70,237, 66, 70, 15,255, 32,167,205,232,190,251, 49,248,  0,243,  1,255,  0,  3,254,  0,255,  5,253,255,  1,255,  4,  4, 29,  1,  5,225, 23,  1,  0,191,237, 11, 18,252,155,  0,  4,255,
      251,196,193,199,  3, 13, 14,145,252,  5, 43,  1, 50,241, 56,198,191,193,195, 61, 12,216,206,  8,193,222,  5, 98,234,  6,  1, 29, 17,  0, 57, 12, 11, 15,248,243,249,  0,  8,255,241,  3, 27,  1,  7, 32,226, 17, 86, 10, 14,196, 21,255,254,254,  4,  0,  2,233,  1,  0,254,255,  1,255,  1,  1,  1,  0,255,  1,  2,  1,  0,255,  0,  0,250,  0,253,  1,255,255,253,  0,255,  0,  0,  0,  0,254,254, 20,  8, 57,245,253, 28, 26,236,  1, 35,253, 15,  8, 63, 67,  0,223, 26, 48, 75, 59,241, 93, 48, 41,  5, 34,250,160, 10,216,  6,  0,242,  4, 22, 15,218,252,  2,253, 10, 24, 48, 93,222,  1, 12,246, 27,134,201,253,  3, 56,244, 16,220,179,209,216, 28, 49,  2,252,  5,246,  2,254,  1,  8,  0,250,237,  1,  1,  0,  4,245,253,  0, 11,246,253,  2,252,  0,  7,212,254, 11,  9,255,  0,188,
    },
  },
  {  //  24
    {
      252, 84, 20, 12,233,240,254,  5,207,236, 12,  7,204,229,251,196, 16,239,  6,237,189, 78, 12, 12,250,  0,253, 21,151, 19,  3,119,  5,249,243,250,  5,246,246, 65, 11, 24,248,  0,235,244,  2,240,254,  9,248,241, 14,  1,249,249, 11,  7,254,255, 50,232,  4, 27,  8,  1,  1,  0,  6,  0,  0,252,  7,  3,  0,  0,  7,  3,255,255,  2,  2,  0,  2, 10,254,255,  0,  1,  4,  0,254,  0,254,  0,  3,242,205, 13,  9,248,  3, 29,255,243, 43,255, 60,243, 40,197,  9,251, 54,206,  8,198,245, 25,201,199, 24,  4,225,252, 54,183,120, 32,236,  1,251, 16, 99,206, 58,211,253,213, 19,207, 63, 69,187,253, 53, 91,126, 35, 59, 42,253,248, 16,244,  9, 13,157,  6,244, 35,250,  0,  8,255,193,248,  0,255,  3,214, 11,  1,  3, 11,227,  6,  6,  9,  2,249,255,  5,  0,  5,  3,  1,  0,251,  2,254,  4,
       20,233,  0, 24,242, 44, 33, 20,255, 50, 22,244,  6, 32,  3, 16,210,252, 12,243, 44,214,187,236,247, 52, 10, 53,200,245,186,212,217,  0,255, 59, 16,  5, 11,  1, 24, 59,  5,  3, 29,243,252, 68,212,  4,253,255, 24,  0, 40,250,193, 51,253,  8,250, 55,249,  1,237,  0,  0,  0,253,  1,252,  1,250,254,  0,  0,251,  0,  0,  0,  4,255,  0,  0,255, 10,255,  0, 36,  7,  0,254,  1, 12,255,  1,223, 47,  1,242,246,255,254, 64,  7,252, 25, 21,  2,222, 12, 86, 43, 62, 41, 80,241, 13,184,248, 12, 14,  1, 21, 21, 13, 39,222,227, 23, 36,249,241,190, 49,  1, 47, 15, 59, 64,  8,173, 80, 13,255,  1,227, 28,226, 41,215,  0,  0,  2,253,220, 63, 56, 99,231,230,249,247,252,  3,250,  6, 55,  0,  5,  0,254,249,  0,  0,  8,255,253,255,243,250,191, 20,252,250,  1,  0,232,  5,250,  2,239,
        1, 13,  6, 60, 42,254, 35,  6, 33, 29,225,217, 42,  7,207,240, 67, 18, 12, 40,242, 50, 33, 43, 15, 19,252,121,  0,197, 65,249, 26, 27,  1,  7,240,  0,254,  1,227,238,254,  0,  6, 27,  8,  0, 13,251,  4,251, 10,242,  0,  5, 51,249,  2, 69, 26,233,  3,251,  8,  0,  0,  0,254,  1,  0,  0,  2,254,  0,  0,  2,  2,  4,255,252,  0,  1,255,245,  0,255,255,251,  2,  1,  0,  3,249,  0,  1,  3,213,223,247,239,226,  0,  4,172, 66,251, 51,186,136, 11,197, 18,229,215,  8,255,  3,237,243, 43,232,243,241, 61,250,  5,123, 38,243,182,234,  0,226,235,226,  8, 32, 54,221, 27,207,237,130, 13,254,199,205,197,197,232,213, 55,  0, 96, 21,213, 29,163,128, 16, 75,  8,217, 23,  3,200,247,253, 50,250, 14,  7,  6,130,202,  2,  3,  5, 10, 32,236, 19,250,  2,254, 19,  0,  1,252,254, 16,
      227,226,251,197, 12, 34,204,253,246,214, 63, 75, 25, 17, 48,  4, 12, 55,  0, 29,254, 58,140,200,253, 24,240,172,233,255,  8, 60, 14,246,  3,251,249,249,248,  0,240, 20,244,255,239,213,255,254, 24,254,254, 13, 23,254, 10,188, 45,215,  7,253,250,223,  2,  9,252,  2,  0,255,  0,254,  0,  0,  1,253,  0,  0,251,  4,  0,  0,251,  0,  0,  0,  7,  0,254,  2,254,  0,  0,255,239,  0,255,254, 14,238,201,227, 24, 44,216,  2,  8,  4,248,  7, 34,255,139, 57,  2, 28,236,  1,  2,107,173,194,230,  9,  1, 20,114,207,101, 25,236,219,  4, 48,  2, 12,  7,145,247, 14,215, 12, 32, 30,232,188,  6,  5,  4, 35,244,254, 49,  6, 42,255,220,231, 67,201,186,243,240,230, 54,  1,234,  0,255, 53,  3,246,200,188, 13, 44,215,229,  0, 21, 30, 13, 13,  9, 18,  9,254,254,247,  2,253,  3,252,217,
    },
    {
       21,211,243,196,  0,245,149,102, 58, 27, 53,144, 26,  4,240,182, 10,  7, 30, 50,253,249,208, 47,198,  9,189,178,221, 29,250,214,  2, 54,254,  4, 13,  1,251,  5, 28,249, 14,  0,243, 11,  1,248, 22, 23,255,246, 12, 10, 10,243, 18, 54,232,  1,226, 11,  1,243,  2,  1,255,255,253,252,  0,  4,  3,  0,  1,  1,  3,  2,  1,  1,255,  3,  0,  1,  2,  3,  1,255,  1,  4,  0,  1,  3,  4,  0,255,254, 13,189,  1, 94, 33, 19,243,221,250,  2,105,  6, 71, 18, 57,230,228, 10, 78, 39,238, 54, 71,225,  4,219,230,117, 65, 44,203, 13,251,230,247,251,253, 11,241,255,  7,225, 40,231, 21,227,249, 77, 59, 12,232,  9,  5,144,  1, 13, 43, 42, 57,184,221, 30,171,247,255, 20,  0,  0,247, 10, 12,246,243,  5,119, 27,232,  1,157,  3,243,248,  5,245,  0, 60, 63,246, 17,  4, 38,242, 26, 26,247,
        4,248, 27, 33, 15, 10, 57,239,249,251,237,252,242,250, 43, 90, 66,237,227,182,232, 28, 27, 22, 20,239,255, 50, 64,243,233,221, 52,220,  3,246,  2,254,  6,  6,230,244,255,255,228,  2,252, 56,229,  2,242,  5,  4,  3,255, 66,255,  5, 51,144,227,246,254,  1,  3,  4,  0,  0,  2,254,255,254,255,  1,  0,  1,255,255,  0,201,  6,  1,254,  0,  1,  0,255,249,  1,  0,255,  0,  0,252,  1,  0,238,231, 64,250,207,135,222, 87, 35, 36,169, 16,  1,250,156,  0,223,211,220,177, 34,248,209,  1, 43,242,  9,241,135, 63,198, 48, 12,  3,236, 21, 40,243, 23, 39, 15,246,  8, 23, 80,254, 88,  0,206, 28, 56,241, 35,245,205, 11, 14, 20,252,  0, 34, 64,229, 57,251, 23,249,  0,238,234,231,  6,244, 51,  4,  8,  3,  0,237,163, 31, 42, 53,223,240, 47,215,184, 22,  2,248, 14,  2,192,250,200,
      223, 51, 51, 25,235,254,237,253,255,251,223, 27, 50,247,254,255, 29,  8,181, 16, 14, 37,202,  1, 10, 19, 52,226,  1, 56,251,  3, 61,244,  2,  5,230,  1,253,  5, 11, 10,254,  5,253, 20,253,194, 17,  1,254,  4, 32,255, 55,  1, 34,249,  2,  2, 47,245,  6,194,252,  1,  0,  0,  0,253,  0,  1,  3,254,  1,  0,254,255,  0,  0,  0,  2,  0,  1,253,  0,  0,  0,252,  0,  0,  0,  7,232,253,  1, 15,  0,223, 33,  5, 24,242,117,252, 38, 57, 52,238, 31,144,158,240, 10, 13,179,251,147,242,  4,252,  7,229,  1, 94,251, 37,132,228,248, 27,229,224, 26,172, 26,238,247,254, 96,175,  9, 66,241,240,175,235, 10, 40,227, 88,212,  4,254,244,245,140,254, 20,252, 12, 29,  5,229,230, 32,247,162,250,209,236, 40,  3, 14, 18,237,199, 19,246,227,230,250, 39,212,246,253, 40,215, 31, 77,182,203,
      215,250, 42,224,  7,254,145,189, 41,248, 41,245, 18, 12,194,201,244,255,198,178,166,248, 23,131,234,224, 60,241,198,203, 50, 42, 11,  7,  1,  1, 31,  1,251,  1, 17,  5,  2,  4,  9, 11,252,246,253,  2,252,  9,  6,251,121,229, 21,250,249,254,254,204,250,  0,253,244,  0,  2,254,  1,  0,254,  0,  1,  0,  0,  7,  1,  0,  0,255,255,  0,  2,253,  0,  0,  2,  0,  4,  0,255,  1,  3,254,  0,235,238,234, 23,230,242, 36, 30, 59,  7,238,  4,  2,225,109,245,225,254, 31,222, 44, 76, 75, 62,  7,  4, 10, 26,252,207,193,216,254,254,241, 20, 15,246,  9,  4,  8, 44,185,231,  1, 69,  2, 42, 23, 29,254, 13,183,201,149,254,254,203,237,250, 23, 45,228,137, 28,235,252, 21,254,250,246,209,230,237, 13, 28,254, 51,155, 70, 16,210,  6, 84,  1,  4, 68,188, 10, 12,238, 68, 24, 17,244,103,
    },
    {
        9, 20,253,239, 19, 13, 71,188,240,  4,231,  1, 63, 18, 57, 48, 34, 32, 54, 94, 25,  4, 14, 22, 48,197,215, 64,  0,173, 87,156, 11, 10, 47,195,  4,  6, 41,238, 13,  3,250,202,246,251, 50,  4,203,  4, 14,  4,186,253,181,164,  7, 40, 15, 12,  0, 95,209,136,  5,  1,247,  5,253,  3,254,249,  2,254,  1, 12,255,  3,254,  5,255,  3,255,254,255,  2,  2,  0,253,248,  0,  1,255,255,  3,  1, 17, 74,249,243,238,  2,255,235,222, 41,252,  5,227,236,161, 44, 32, 27,226, 46,200,202,  4, 72,229, 36, 37,246,246,229, 69,121, 10, 45,236,250,  6, 15,243, 19,244, 28, 55, 12,  2, 15,246,228,248, 11,  0,188, 31,186,133, 17,  0,250,252,253,255,  3,  0,  2,  2,247,  7, 17,255, 15,254,249,  2,  7,233, 31,  1,247,  7,248,240,  2,255,  2,  0,  3,  1,  2,  0,255,  1, 15,  0,255,255,  0,
       73,  4, 52, 45,  3,241,  0,243, 11, 16, 64,236,196,250,218, 13,  0, 13,160,170, 45,251,226, 46, 60, 44,110, 90, 26, 21, 74, 63,191,  2, 20, 13,251,  1,250,236,229,  1, 36, 51,  7,  5,243,205, 48,119,226,221, 47, 18,235,199,  0,255,215,243, 29,198, 71,215,252,252,254, 14,  2,251,  3,247,255,  0,  0, 12,255,  1,254,248,255,  0,  0,  1,  8,251,255,  2,  3,255,  1,254,  0,  1,253,253,  5, 59,241,  6, 32,186, 15,243, 21,252,226,223,219,223, 26, 33,  0,232,207,230, 55, 47, 18,250,212,199, 16, 24,128,  1, 11, 71, 17, 14,227, 20,202,  0, 20,213,208,243,  3,194,253, 56,  3, 82,185,246,253,125,187, 53,238,248,  2,180, 18,  8,124,254,247,  2,  1,  5, 21,245,254,  0,251,  2, 27,  2,233,221,254,  0,  3,222,237, 52,  1,254,254,  5,  0,254,  1, 26,  0,252,  0,255,254,  0,
       10,253,  7,209, 50, 17, 86,255,  4,254,  4, 21,  2,245, 49,243,254, 24,145, 53, 19,231, 93, 24, 25,219, 22,222,  1, 26, 75, 62, 54,251,  4,208,249,  2,108, 63, 28,248,220, 14,  3,  2,214, 23,250, 69,181,231,  3, 10,220,255,  2,  1,189,229,201,  0,122,127,  0,253,  3,232,  2,255,  4,  1,  2,253,  1,  6,  0,254,  2, 13,250,  0,  0,  3,255,  1,  6,249,  0,  0,214,253,  0,  0,254,  4,254, 12,250, 42,250,229, 56,252,246,  4, 59,  6,228, 24,232,212,238, 41, 74, 19,199,100, 60, 67,  2, 14,229,187, 45, 65,  1, 94, 61,  0, 25, 41, 50,198,250, 20, 44, 36,253, 12, 99, 39, 16, 20, 68, 12,  9, 62,  5,128,242,255,253, 58,251, 64,250, 52,245,191,239,  0,251,234,255,  5, 20,254,255,231,222, 19,253,  2,255,253,  0,198,  1, 15,255,  0,245,255,255,243,250, 19,  1,254,255,  3,
       37,254,227,  5, 36,  5, 33,221, 35,250, 19,190,205,247,182,166,  8,241,173,176,246,249, 34,242,209, 13, 53,  3, 63,242, 33,249, 34,  1,224, 51,  3,  2,227,196,  7,  5,194,177,225,  1,247,230,  2,  4, 74,253, 17,  0, 27,105,243,  2, 64,179, 21,195, 90, 62,254,  0,  2,  2,  1,  0,  0,253,253,  3,  1,247,  2,253,  2,  6,  0,  0,253,  3,255,  1,255,  1,252,  2,  0,  1,  1,  0,  2,255,  2, 14,255, 37,195, 59,167,234,226, 10,201,  8,  0, 30, 54, 16, 21, 29,168,217,  1, 79, 60, 19,  3, 43, 34, 71,  1, 80,135,166,220,229, 35, 22, 17,  7, 16, 74,251,207,227,187, 72,247, 54, 18, 12,229, 10,  8, 70,245,188,  6,  0,243,226, 40,253, 11,  2,198, 20,250, 21, 16,  2,248,255, 64,  3, 10,253,243,  9,  4,251,255,254,248,255,  4,240,254,  1,  2,247,255,  4, 15,253,  0,  4,  8,
    },
  },
  {  //  25
    {
       17,253,213,233,214,244,236,218,230,232,189,  4, 17, 25,201, 28,242,  1, 61,225,222, 95, 71,213,218,241,247,249,249,235,241, 24,253, 10,225,  1, 11,254,243, 14,253, 12, 54,  9, 36, 16, 43,234, 29,255, 19,  0, 16,254, 34,  0,219,245,241,254,242,100,  2,242,  1,  2,255,255,  2,255,  2,251,  0,255,  1,  1,  0,  2,  1,  2,  2, 10,  2,254,254,  3,  5,254,  9,  4,254,  1, 16,250,255,  1, 18,  8, 56, 36,181, 19, 32, 46,  8, 13,239, 28,231,239, 36, 32,186,  5, 13, 39,  8, 23, 37,211,196,249,172, 20,186, 45,227,254,  5, 57,253, 11, 91,252,221,212,  8, 13, 37,248,241,207,  0, 42,222, 76, 37, 31,233,  4,251, 66,253,  8,252, 27,245,242, 24,156,  3,251, 21,241,253, 40, 86,251,  6,  0, 68, 39,249, 31,242, 98,  3,247, 17, 42,239,  0,244, 10,246,255,  9,  3,238,  2,248,  0,
      216, 33, 74,199, 45, 29,106, 13,  8, 30,  8, 24,205,230, 45,157,203, 32,  6,252,  2, 14,  1,  0, 11,214,215,228,255,221, 30,221, 26,  5,250, 43,252,252, 27,251,239, 23, 31,248, 32,221,243,243, 73,215, 52,252,235,  1, 42,227,  3,  4,  1, 11, 45,208, 12,243,  1,  0,254,  6,255,  1,  1,  0,254,250,  0,253,  0,  3,  0,  2, 10, 15,254,254,255,  6,252,  2,  4,  1,  0,254,  4,  4,  1,  5,246, 27,  6,251, 15, 37, 16, 10,  7, 24,222, 10,  3, 19, 47,118,227,253,241, 45,  3,236,250,212,  8,187, 13,241, 85,231,  8,238, 13,253,228,255,  9, 18,250,  6, 11,240,252,  6,218,253,190, 85, 45, 12, 10, 30, 65, 30,253,182, 11, 33,242,197, 16, 11,  3, 59,253,253,  3, 24,203, 20,241,164,245,  6,189,189,  1, 18, 61, 83,245,  1,197,  7, 40, 17, 16,245,  0,231, 16,245,255,  6, 15, 29,
       14,  4,217, 19,  5,255,218,127,253,231,228, 15, 75, 12, 31,183, 65, 47, 27,200, 11,  9,233,100, 51,245,245,248, 30, 55,194,182,238,251, 69,249,252, 11,249,  4, 18,233,197, 11, 36,233, 65,250, 39,247,185,248, 15, 14,190, 63,  5, 10,244,238,214, 19,199,  1,255,  0,  6,  2,  0,255,  4,  0,  0,  7, 46,250,255,  1,253,255,  3, 11,255,255,  0,254,253,249,252,  6,  0,253,255,255,253,  3,240,253, 55, 18, 74,242,252,236, 21, 78,247,246, 53,  0,172,250, 28, 68,  6,  6,  8,  2,  0,245,234, 18,255,253,  4,251,182, 63, 12,249,252,193,255,192,  4,  0,  1, 29, 30, 42,220, 22,  0, 53, 11,212,  8, 26,216,234, 49,  5,246,212,245,254,186, 58, 10,238,  8,  6,  8,203,240, 42,236,  1,  8,237, 25,233,140,197,255,248,246,232, 53, 10, 22, 32,249,198, 14,253,244,  6,186,236,255, 33,
       41, 19,191, 14,253,235, 64, 55,  0,255,  1,254, 27,  2, 35,254,207, 18, 39,224,220, 19,  2, 40,  7, 20, 11, 65,  8,210,225,180,251,  6,230, 14,  2,237,211, 44, 18, 25, 63,  2,214,246,194, 13,255,  1,  4, 54,204,238, 56, 65, 23,  8,244,230,255, 53,248,239,  0,  2,250,254,  0,  0,255,252,254,255,  0,  2,244,  1,  2,  1,  5,  3,253,  9, 12,254,255,  0,252,247,254,254,253,254,  1,  2, 21,  2,191,231,  7,185,  0,138, 16,246, 59, 13,180,217, 63, 12,237,217,  0, 24,208,  1,127,161, 17,  7,235, 29,217,253, 13,136,247, 13, 57, 42,241,153,  5,215, 11, 10,215,  5,248, 53, 63, 20,236,178, 16,228,244,205, 62,204, 61,236, 27,205,  1,243,103,  7,  6,212, 33, 13, 54,254,244, 61,  4,245, 15, 30, 62,203, 37, 20, 12,216, 19,230,  5, 66,237,137,  7,252,218, 33,253,252, 10, 92,
    },
    {
       19, 28, 50, 53, 23, 28, 44, 17, 25, 34, 51, 95, 12,248,245,251,235,  1,203, 98, 29,237,  9, 21, 50, 28, 67, 68, 11, 65, 56,203,  2,  9,177,240,  2,  9, 40,250,  3, 21,212,232,  2, 10,223,244, 10,  6, 26, 58,228,244,233, 38, 44,248,253,249,204,  2,  7, 70,  0,  3,  2,  0,255,  3,253,255,  0,  1,  2,255,  1,  3,  0,  1,255,  0,  0,254,  4,  0,253,  0,  2,  8,  2,255,  1,255,255,  7, 33, 16, 33, 17,203,248, 15,216, 20,248, 24,  1, 81,  6, 46,  8,227,240,196,222, 19,  3, 82, 10,234,226,158,218,227,253,229, 53,242,255,249,250,247,246,251,241, 12, 25, 10,220,254, 16,233,199,  7,  0,220,195,211, 45, 18,232,  6, 38,237,199,245,211,251,223,  0,255,251,254,  3,  2,253, 19, 53,213,222,  6,204,255,  0,166,246,  1,254, 14, 11,251,255,235,254, 13,248, 49,253,  1,  2, 24,
        1, 32, 28,224,238,216,240, 12,244,245, 49,240,252,222,200,216, 50,253,246,213, 16, 18,242,147,253,215,208, 31,  0, 13,117, 64,  0,  2,  3,  8,253,  4,255,252,  5,247,245, 49,  7,251,  2, 79,242,255,243,191,255,  9, 85,133,  2,228, 29,169,198,250,232,118,255,  0,  1,  1,  2,255,253,252,255,  0,255,  7,253,  0,  1,254,  0,  0,254,  7,  0,  0,255,  1,  0,  2,  1,  0,  1,255,  1,  0,  0,196, 58,204, 24,230,213,225, 67,227,238,235, 29,170,255,238,  6,242, 10,235,240,244,252,223, 47, 32, 64,218, 23, 10,248,  0,  3,226,211,230,234, 19,249,240, 21, 36,251,248,231,247,221, 44,207, 25,  0,205,237,215, 48,167, 16,250,252, 41,  2, 20,216,255,255, 25,  3,230,  0,255,  2,208,235, 19,255,254,126,255, 51,163, 31,210,253,251, 18,235,216,137, 14,253,  7,228,226,  6, 17,  0,
      255,230,225,  2,208, 41,233, 57,180,239,134, 55, 62, 16,228,235,212,173,182, 45,240,182,250,200, 37,251, 25,230,244, 51,190,196,254,250, 16, 45,243,  2,243,190,  4,  3,237,148,  1,  3, 67,131, 35,  1,255,204, 52,247,189,143,186,239, 26,158,240,172,  3,206,  0,255,255,  5,  0,  0,  0,245,  6,  2,  0, 11,  0,  0,  1,  2,  0,254,242,  2,  1,  0,252,255,254,  0,254,253,  0,  1,  0,254,  0,  2, 48,255,226,230,188,253,  2,241, 14,216,  9,248,230,235,  4,239,  3,  3, 48,249,  7,136,253,172,234,210,247,245,255, 13,255,215,253,  8, 28,227,219,165, 44,216,  9,  0, 37,248,255, 22, 73,242,254, 28,224,123,247, 41,235,237,251,215,203, 58,255,234,  0, 25,  8,183,  1,245, 10,191,255,  1,252,240,198,247,254,253,221,  0,  0, 28,234,237, 63,255,246,246,245,255,  3,255,191,  6,
        3,  2,240,247,242, 13,249,224, 55,234, 30,182,213,202,229,192, 62,206, 95, 47,250,254,199,253, 67,  3, 15,255,164, 11, 52,229,  1, 58,  2,204,250,  1,  0,254, 19,253, 43,252,  4,251, 19,254,255,  2,184,247,209, 11,196, 43,198,  8,242, 15,238,238,196,249,  0,251,255,  1,  0,  0,  4,  3,  0,250,253,253,  0,255,250,  0,  3,  2,255,242,  4,  0,  0,  0,  3,255,  1,255,  0,255,  0,  2, 12,255, 25,224,137,131, 66,  4, 14,250, 65,250,236, 55, 44,126, 37,  3,  6, 34,250,214, 15, 44, 13,237, 40,101,  3, 22,241,251,  1, 27,246, 15, 38,214,234, 83, 65,  3,240,249,  3, 85, 69,  0,242,251, 10,239,236,234,234,187, 28,247, 25,226,  8,211, 15, 75,  2,187,252,240,232,255,  2,185, 27,249,  2,225, 55,171,243,200, 18, 28,209, 16,226,238,246,137,  0, 16,243, 16,223,230,  0,  1,
    },
    {
      209, 29, 20,251,  5,  4, 75, 19,225,253,251, 96,180, 74, 14, 61,  5, 27,226,  1,234,233,  9, 42, 11,208,247, 97, 37,244,253,197,250,251, 49,164, 11,  3, 42,  1, 16,252,  5, 18, 60, 13,250,  4, 12,  6, 62,251, 25,198,252, 60, 12, 13, 68,234,235,233, 60,  2,255,254,245,  4,254,  2,255,  1,  4,255,  0,255, 13,254,  0,  0,  0,  4,  4,  1,240,251,255,255,  0,  2,  0,255,  2,254, 14,  1, 19,  9, 18, 54, 90,192, 16,167,255,215,177,238,202, 14,145, 52,248,  4,247, 43, 64, 15, 56, 60,254, 68,218,  2, 98,220, 29, 59, 14,  3,  1, 15,  3,249,255,253,247,225, 39,251,  0,255,  0,132,201, 65, 64, 33,  1,121, 61,156,247,212, 40, 11,216, 28,228,200,  0,  0,  1,  3,  0,251,  2,  1,  0,  2,251,  2,249,  1,251,255,215, 21,253,  2,127, 12,255,246,255,  0,254,245,  2,202,  1,252,
       57, 48, 78,246, 14,  3,  6,227,  4, 10,185, 22, 41, 21,  0,140,  0, 60,187,202, 29, 31, 76,247,239,226,191, 97, 39, 34,248,243,251,171, 59,241, 13, 51,  2, 62,224,254,124,  9, 16,215,138,247,  3, 14,240,248, 12, 20,200,190,  0, 38,251,  2,239,233,181, 59,254,  1,245,253,  0,253,253,  0,  1,253,  1,  0,  3,  1,  0,255,  3,200,  0,245, 24,255,255,  1,  1,  0,  0,  1,  0,  2,  0,  7,252, 64,254, 46,  7, 63,245,196,235, 52,214, 10, 79,254, 13,232, 22, 53, 12, 64,174, 63, 47, 69,  0,239, 43,  8,141,233,205, 11,255,254,248,  2,  0,250,  1,  1, 14,247,248,246,  6,  1,255,133,255,220,203,207,118,180,122,210, 50, 20,249,240, 29,252,175,135,  1,247,252,250,254,  0,  0,  0,255,  0,  1,  1,  0,  4,255,  3, 33,222,228, 13,251,  1,  4,241,253,  6,  2,  3,  5,224,251,  2,
      252, 26, 58,178,238,219, 93,215,251,249, 59, 25, 87, 44,  1, 48,250, 61,254,232, 21,  2,203, 57, 16,253,240,149,252, 14,  5, 25,  7, 78, 74, 49,228,249, 87,162,  8, 64,254,254,235,252, 20,190,  2,243,253,190, 23, 16,192,245,254,227,  0,209, 12, 45,221, 55,253,  2,  3,244,253,  2,254,255,  1, 28,255,255,  0,255,255,  1,245,  0,  0,  1,  8, 64,244,255,  0,  3,  0,  1,  1,  7,  8,254,  6, 43, 36,250, 30,175, 77,  2,250,206, 82,225, 74,  4, 94, 12, 21,  6, 48,247,210,255,254, 93,  6, 12,253,197,146, 19,123, 32,  2,  6,243,254,  6,254,249,  0, 18,255, 14, 36,  7,  2,  2,191, 68, 17, 23,210,221, 71,216, 10,239,238,  1,  8,200,130,224, 90,  3,252,252,  4,252,254,  3,  0,  0,  0,  0,255,255,  1,  1,225,  3, 11,211,253, 11,  0, 36,255,  0,253,253,  3,109,255,254,255,
        5,251,228,200,  2, 55, 43, 63, 36, 43, 62,152,254,  0,250,112,209,247,241,196,  0,  4, 60,  9,248,216,253,199, 83, 23, 14,252,  0,  8, 70, 14,  4,226,128, 38,243,239,255,  1,250, 62,  3, 61, 48,  0,249,  4, 69, 72, 37,229,221,237,250,128, 37, 20, 44,255,254,  2,  0,  1,  0,  1,  0,253,249,  0,254,254,  1,  0,  0,  5,  0,255,  2,255,250,254,233,242,252,255,255,  3, 25,250,  2,  0,252,190, 26,246,231,227,192,151,229, 12, 45,211,229,220, 58, 29,231,255, 10, 16,  9,179,190,195,241,253, 38,237,191,182, 76, 63,  2,229, 41, 18, 18,237, 69, 17,251,  5,  7,201, 15,222,253,162,237,243, 36,239,  0, 73, 31,191,  3,230, 15,104,196, 62,223,134,  1, 16,245,252,255,244,201,  2,253,  1,253,239,255,254,254,122,  2,254,  5, 27, 12,248,193, 26,250,  5,  3,  2,236,254,252,218,
    },
  },
  {  //  26
    {
       31, 37, 41, 25, 18, 10, 42, 34, 20, 49, 26,126, 10, 21,183, 48, 25, 56, 12, 55,  3, 16, 74,133, 27, 83,208,214, 13, 29, 18,249,  2,233, 46,  1,  9,213, 60,243,  0,  9,205,  4,  7, 58,196, 14,  5, 27, 15,  1, 12,206,220,  2,  8,  9,129,248, 31,157, 51,186,  0,  9,249,254,255, 12,253,  0,255,  0,  8,255,  0, 10,249,  7,  1,  1,  0,  2,  0, 27,254,  0,254,  9,  0,  1,255,255, 10, 12,235, 20, 36, 16, 40, 17, 56,105, 12,  2,192,  1,234,245,212,  2,244,  2, 49, 63,197,127, 47, 18,  4,255,255,  4, 51,192, 54, 50, 12,255,  1,  8,250,  6,  1,241,  5,  0,245,  1,252,  0,  1,255,240,245, 48, 43,123, 66,189,  5,108, 17,244,  7,184,228, 67,192,  4,  3,  1,255,  2,  0,  1,254,254,  1,  4,  3,  1,  0,  0,  1,  4,  3, 11,  5,255,252,253,  1,255, 30, 10,  1,  0,  5,  1,  1,
       16, 20,255, 37, 16, 15, 83,149, 50, 14, 61,129,255, 22, 61, 40, 29,248,  7, 42,243, 20,  3,  8, 41, 27, 61,213,  1, 14, 45, 65,  7, 24,223,238,249, 57,233, 11, 13,246,234, 89,  1,193, 30,255,252,249,193, 76, 37, 61,  5,  1,226,246, 16,194, 39, 57,  1,239,255,252,  0,252,254,  0,  2,254,  2,  0,246,253,  1,255,239,  9,  0,254,252,254,253,255,  0,  1,253,254,255,255,  7, 44,  0,255, 25,242,236, 25, 21,250,  6,246, 19, 47, 19, 64,101, 34, 17, 49,198,207,  5, 20, 33, 75,204,235,233,206,238, 44,255,238,193,142, 14,255,253,  8,  0,254,255,  1,  6,231,206,  2,247,238,248,254, 11,239,246,247, 24, 13,250,114,196, 18,  8,246, 74,178,114, 77,255,254,254,  1,254,254,  2,  1,  0,237,244,  1,  2,255,  0,254,254,251, 23, 20,  8,  1,255,255, 11,  8,  0,253,  0,  7,250,254,
      227,254,208, 18,253,191,191,244,254,  2,195,  4,254,  5, 62,219, 13,  1,222,162,  4, 25, 12,215,  5,252, 65,127,249,  6, 41, 69,255, 27, 55,250,240,255,244,217,234, 37,192,240,  4, 13,198, 11,176,  6,  1,247, 48,206,  7,163,195, 88,250,  7, 31,245,243, 49,  1,  1,248, 10,  0,  0,251,255,  0,252,  0,254,  2,  1,  9,  8, 24,  0,  0,255, 40,  1,  0,  0,253,  3,  1,  1,  1,  1,  0,253,  6,251,248,  5,226,210, 11, 18, 13,128,  0,222, 46, 82,254,219,219,  7,235,249, 13,176,132,205,159,250,195,  3, 18, 11,159,186,234,  7,  5,  7,  0,  0,255,  6,  9,229,  4,231,255,254,248,252,243, 69,230,238, 50,244, 55, 59,  0,244, 10,241,190,236,176, 58,  1,  2,  0,  1,  4,255,  0,  0,  0,  5,  6,236,254,  3,  0,  0,252, 24,  1, 38,  0, 23,255,  1,240,  6,  1,  6,255,  7,  0,254,
        7,245,245,250,  0, 19,175,197,  2,  5,242, 70,247,  8, 11,253,  2,  6, 27,244,238, 72, 38,245,  3,248,229,196,  3,114,255,248,243,248,190, 50,  2,250, 48, 54,226, 24, 25,254,  1,248,  5, 32, 28, 49,250, 51,  1,253,  3,255,179, 42, 51,  1,253, 49,  0,233,  3,  2,  0,253,255,  0,  0,  0,  1,255,  1,  0,254,  0,246,247, 10,  0,255,  1, 53,251,255,  0,  2,253,  3,  0,244,255,254,  2,255, 57,227,  8, 17,188,254, 75,250, 52,173,200,172,190, 31,206, 12, 56, 63,246,214,201,203,  1, 49,227,252, 15,  5, 66, 59,  8,  9, 19, 20, 25, 33, 12, 13,  1,  1,251, 64,  7,202, 14,210, 20,194,194, 12, 12,232,156, 47, 65,  1,201,241,241, 67,176,193,143,  2,245,  8,253,254,  1,254,  0,  0,  3, 11,242,247,250,192,  0,255,255,246, 72,254,240,251, 11,  5,  4,  1, 33,  1,248,252,  4,
    },
    {
      245,221,  0,227, 22,211, 23,235,244,  9,224,242, 20,233,250,249,227,228,  3,  2,248, 22,251,  2,192, 14,218,219,245,233,254, 10,242,  8,254,248,251,247,250,254,229, 12,252,  0, 48,255,251,  3,236, 21,253,250,217, 16,253,254,241, 11,253,  3,239,245,  0,251,254,  0,254,255,  0,255,255,255,  3,  2,254,  0, 19,252,  2,255,  4,  0,  1,254, 14,  1,255,254,  8,  8,  3,  0,152,  7,  2,254,237, 33,241,242,234,239,  4,  3,239,  2,252,206,239, 67, 28,106, 76, 23,245,253,200, 17,218,224,225, 99,  3,254,219, 39, 32,202, 41,227, 40,232, 20,242, 68, 91, 40, 62,  5,  3,  7, 18, 46,226, 26,238, 58, 20,  4,188, 62,  7,  1,217,241, 28,204,239,143, 64,  0, 32,255, 23,244,223,252, 16,251,  7,234,  2, 33, 26,238,249,252,  4, 17,241,  2,216,249,  7,253,  1,255,233,209, 36, 12,  5,
       28,233, 18,  5,  2,254,195, 68,238,250,223,203,208, 38,252,  3, 10,246, 21, 33,  8,244, 57,141,231,  9, 13,193,254, 10,  8, 53, 30,  5,250,252,  7,241,  0, 19, 12,  8,255,254,252, 16,  6,  2, 66,252,  3,  1, 52,  6,  3,255, 65,246,  2,253,227,236,  0,  2,251,  1,  2,  0,  5,  1,  0,255,253,253,  0,255,245,  0,  2,  1,252,  1,  2,  3,  0,  5,255,253,  8,  3,255,253,203,  1,  1,  2, 28, 17,236,227,  1,238,199, 24, 11,240,247,206, 15,229,  2, 17,  2,209,194,203,188,223,213, 94, 16,177,206,210,228, 32, 42,246,238, 41,253,233,233,  7,221,222, 15,221,144, 12,255, 52,  0,216,185, 31, 66,215,189,245,244,230, 14,  4, 10, 17, 79, 62, 11, 12,  1,247, 12,  1,  6, 18,232,244, 22,236,224,244,  3,238,253,214,255,253,255,  0, 10, 31,118,212,  0,248,  2, 11,  0, 20,  0,255,
        4,226,201,225,252, 15, 15, 10,211, 28, 52,209, 28, 76,  8, 85,  9, 44,  2,253,233,251,231,240,194,245, 25,254,215,210,233, 22,247,253,  9,  0, 20, 11,245,  0, 10, 25,252,  0, 10,240,  4,242, 35,254,255,252, 19, 20,  2,  2,218, 25,  5, 61, 43, 36,253,  9,255,255,  1,  0,248,  2,  1,  0,  8,252,  0,  0, 14,  4,255,  0,244,248,  2,255,248,  3,  0,255,  0,254,  0,  4,  7,  1,  0,  2,243,178,241, 23,231,255,  1,223, 11, 28, 31,  2, 21,222,128,248, 62, 34,  1,236,231, 27, 42,212,251, 54,254, 17, 51,241,150,196,245, 26,232, 22, 13,250,  7, 14,225,226, 15, 45,246,227,184,154,255, 18,222, 18,250, 61,217, 13,253,235,239, 40,242,196,  0,192,  9,  3,  3, 13,239,125,249, 21,255,227, 15,243,  4, 22,  6,251, 11,  0, 14,  9,247,  0,189,115,254,248,  1,  2,  0,  1,  3, 18,
      247, 78,223,225,246,202, 34,239, 34,195,235, 80, 21,255,252,241,234,210,255,  0, 16,216,254, 72,248, 13, 91,102, 19, 65, 27, 48,  4,  7,252,255,  7,  0,240,  3, 39,  7,  1,  0,221, 16,253,251,237,206,252, 17, 15,  4,250,236, 43,  4,255,255,175,130,255,  0,  7,254,  0,  1,255,  2,254,  0,  0,  6,  1,  0,248,  2,  0,  0,252,254,255,  1, 53,254,  0,  3, 17,  9,  1,  1, 10,  7,  4,  0, 15,220,  2,184, 56, 76,204, 46, 13,189,219,248, 19,152,252, 67,  6,243,245, 22, 60,197,193,  2,214, 72,232,239,  7,  7,232,245,247,218,244,244, 30, 52,186,  2,224, 24, 33,237,254, 49,248,129,234, 60, 16,  5,  3, 12,  6,219,  7,  1,219,255, 25,250,221, 71,248,252, 16,246, 27,197,191,190,254,250, 12, 19,253,255,  6,254,253,253,249,  9, 34, 40,251,238,  4, 14,241,  4,237, 18,  6, 76,
    },
    {
        9,  4, 12, 51,  6,  9,222,232,  0, 27, 31,219,239,230,  5,253,243, 34, 50, 66,231,252,219, 88,254,253,  5, 27, 11, 10, 22,211,  4,  1,253, 16,  4, 11, 13,248, 22,  5,  1,  5, 31,  1,242, 19,  5, 11,250, 20,218,  2,  2,  3,243,240,  0,243,244,230,254,255,253,  4,  1,255,253,  4,255,253, 11,  0,255,  1,  3,  0,255,253,250,  3,  3,  1,  8, 11,  1,  2, 22,  3,  4,  0, 16, 23,254,  1,  4, 25, 38,  7, 10,236, 56, 37,252, 40, 38, 26,241,197,171,233,249,230, 30,247,  1, 59,222,  0, 69,  2, 18,214,189,215, 12,  7, 42, 46,222,202,227,214,197,199,253,250,244,159,  9,233, 11,248, 20, 24,255,243,232,191, 11,226,  9,229,  5,218, 58,  3,  4, 35,  3, 48, 60,  3,238,  0, 72,215, 56, 17, 28,253,241, 57,240, 59,  4, 30, 11,245, 22,229,224, 19,248,  5,232,248,237, 21,251,236,
       14, 74,  8, 25, 16, 51, 34, 84, 19, 27,221, 25,226,233,  3,  5, 44, 10, 49, 51, 10,  2, 45,242, 16, 49,221,238, 19,181,234,111,255,  9,  0,  0,231, 33,255,  1,  7, 15,250,  4, 35,240,252, 16,  0, 11, 32,251,239,253,  5, 11,208,234,244,  7, 12,  6,253, 10,  9,  1,  2,251,  1,254,  1,  1,255,  6,  2,  2,247,  1,255,255,  0, 17,  1,255,254,255,  1,252, 52,255,254,254,  0,  7,  1,  0,251,242, 44,  0,253, 66,237,208,234,  9,237, 25, 35, 25,239,145,215, 12,  1,243, 55,206,141, 26, 10,226, 31, 29, 15,249,  4, 51,245, 10, 18, 16, 24,196,255, 35,  6,  2,  9, 52,252, 38, 58,165,  0,228, 25,234,217, 11, 53, 19,253, 24,  6,245, 29, 30, 12,224, 45,  4,248, 30, 33, 41,250,252,210,220, 61,252, 42,194,255,239, 13,255,246,247, 23,102,193,199,  5,254,252,  9, 13,186, 77, 83,
      250, 31, 13,244,241, 19,249, 17, 31,241, 33, 18, 17, 16,240,194,  1, 43,208,225,  1,212,231,  7,190,210,227, 49, 11, 62,221,  8,250, 13,  1,250,250,243, 12,  1,222,234, 15,  4, 16, 48,  9,  0,248,241,  2,254,204,182,  1,246,251,247,  3,  1,249, 54,  2, 26,251,253,  1,  1,251,  0,  0,254,  0,252,  2,  0,  1,  0,  1,  0,251,253,253,253,246,231,  2,255, 27,252,  4,253,253,251,  0,  1,218, 38, 36,  9,234, 62,  0,237, 73,  3,225,  9,187,176,131,  5,  3, 11,248,  7,253, 20,187, 37,  4, 30, 46,244,246, 27,219, 32, 29, 22,216, 12,202, 27,167,  6,246, 23, 10,255,248,196,250, 69,252,252,250, 42,238, 73,188, 20,237,  7, 30,245,230,  3, 11, 35,177,  2,254,245, 39,183, 58, 41,  4,  3,238, 38,201,155, 41,140,  0, 18,241, 23,211,185, 33, 34,  4,239,238,243, 16, 65,232,253,
        9, 11,202, 11, 20,216, 58,251,217,200, 32, 71,248,  1, 28, 15,243,245,253,231, 58, 47,218, 41, 36, 39,254,238,214,244,214,  0, 10,249,  0,255,241,  9,243,  3, 15,233,252,  2, 14,214, 31,  6,229, 14, 16,  5, 74, 52,  8,253, 24, 23,249,251,251,233,  0,249,255,  1,255,253,  4,  4,  2,  1,252,  8,255,  1,  3,254,  4,  3,250, 13,253,  4,244,  9,  3,  0,246,  2,255,  2,  9,  8,  4,251,234, 17,245, 13,  0,223,  8,134,239,249,240,  9,  4,  1,205, 41,  5,224,225,  0,229,223,214, 51,221,  5,239,248, 39,212,210, 23,228,232,  2,  7, 74,240,  6,231,255, 24,255,  7, 62,233,255,171,247,  1,  2, 62,213,  1, 72,232, 18, 72, 55, 72,207, 40,  6,198, 11,224,224, 18,183,235, 58, 45,255, 17,221,238,249, 30,221,191, 13,234, 10, 24, 47,  2,  9,189,  0,249,253,254,  6,183,229, 19,
    },
  },
  {  //  27
    {
       38, 19,241,  8, 10, 17, 67,  1, 25,184, 43,224,242,179,246,  7,225,199, 14,  3, 14,  8, 62,188,255,172, 29,200,237,249,252,  9,229,218, 24,  0, 63,226,  8,251,249, 10, 49,  5,254,  9, 43,252,231,  1,248,254,  8,227,  4,  4,253, 26,  2,  5,251,238,  3,  0,253,  0,255,  0,  1,255,255,255,253,  2,  6,  0,254,  3,  4,  0, 13,  5,  0,  0,  4, 47,  1,255,187, 13,  0,  4,181, 20,  5,  1,249,230,  1,243,  8, 33,  6, 10, 62,219,196,  3, 61,255, 45,207, 68,  2, 11, 20, 61,230,245,198, 88, 23, 25,196,115, 10,182, 60,249,252,255,  8,  1,  4,  6,251,237, 40,  0,221,248,253,253,230,  0, 57, 40, 71,249,240,  0,253, 40, 21,242,249, 70,248, 52,148,254,  1,255,  0,252,  9,  0, 85,  5,255,253,  4,239,247,  2,253,  4, 33, 47,222, 14,234,  2,135,  2,255,254,  1,  0,214,255,251,
      253,250,214,254,251,  2,  4,  7,245,  9,251,254,254, 21,  0, 10,179, 46,241,246, 37,237,156,132,253,  8, 70, 22, 57,243,  4,237,  0, 36, 39,253,188, 24,  3,  8, 22, 64, 16,  0,245,254,235,251,237,251, 16,252,236, 11,129,255,220,222,  3,254,201,245,253,  7,  3,  0,  3,255,255,  3,255,  1,  1,255,  8,  0,  3,  0,  3,  0,245, 14,  1,253,255, 11,255,  3, 10,  8,255,  0, 66,252,254,  1,  8, 16,199,234,155, 13,243,251,252, 34,221, 24,232,  6,  0,  8,219,  3, 33,  1, 28, 60,  3,251,209,251, 44,  4, 68,244, 12,205, 13,248,210,  5,  6,  1,255,  7, 27,223,252, 31,  6,250,  5, 16, 18,217,255,  1,186, 20, 63,  4, 11,215,252, 53,228,237,245,233,254,253,238, 18, 12,253,  0,155,255,  0,  0,  0,248,248,  0, 10,255,199,  3, 32,250,195,196,253,  2,255,  2,  1,235,249,  0,250,
      230,  9, 50,254,250,246,197,251,197,  3,180,254, 35, 70,213, 33,174,239,116, 61,212, 17,244, 10,255,192,241,  1, 27, 42,205,  0,254,240,203, 10,227, 43,249,  8,254,232,251, 56, 52,244,241, 26,  7,251,237,  3, 15,254,127,  2,254,  1,225,  2,170,  1,252, 15,  3,  2,  2,254,  0,255,254,  0,251,  0,  0,  0,  2,255,  1,  0,249, 17,  0,255,250, 88,  0,255,224,  2,  0,251, 19, 14,  2,  2,220, 43,251,219,226, 47,  6,183, 36, 95, 20,254,  0,  0,249,219, 35, 26,  5, 21, 24, 14,132, 35,  8,  0,254, 71,  5,206, 65, 65,254, 22, 71,246,242,  0,255,255,254,  5,  9, 15,  9,252,250,  4,253,233,242, 26,  4, 26,  3,244,193, 72, 20, 20, 39,245,  0, 25,253,  1, 37,  0,255,  2,249,  1,  1,  0,255,  3,  4,  5,255,254,  0,  7,222,250, 11,255,255,107, 15,254,252,255, 32, 48,  6,  5,
       11,237,212,253,  8,250,  3,248,202,  1,255, 20, 77,211, 39,199, 61,252,230,253,163,228, 22, 69,249, 19, 89, 67,172,104, 60,250,254, 14, 46, 11, 37,232,255, 13,226, 70,  6, 59, 46, 39,243,  6,238,249,  2,252, 51,254,  1,  0,246, 58,243,254, 17,245,231, 13,  2,254,  0,  0,  1,254,  2,255,255,  0,  0,  0,  1,254,  1,254,  0,246,255,  1,  7,232,  1,  0, 31,  0,  1,  0,238,  6,246,251, 32,247, 23,105, 69,147, 33, 79, 22, 11,211,227,205,241, 35,  5,  0, 19,  5, 63,153,232,254, 55,255, 35,239,252,102,209,190, 68,250,241, 29, 20,242,255,206,251,  1, 30,  6, 57, 75, 13,209,225,255,255,  3,221, 35,237,190,196, 20, 40,252,205, 60,211, 54, 35,  0,  3,226,  9,236,247, 53, 11,255,253,  1, 16, 30,246,200,  2,255,  7,236,249, 49,  0, 10,  2,249,  0,251, 14,249, 37,253,249,
    },
    {
       31,  4, 48, 14,249, 10,  8,183, 19, 24, 14, 71,  3, 19, 68, 13, 17, 49, 35,196, 36, 73,  1, 11, 27, 58,108, 59,221,  4, 16,230,250,  1, 52, 68, 42,  9, 63, 55,249,  6, 62, 84, 67,  6, 55, 64, 33, 44,194, 58, 30, 44,248, 73,243,  9, 81, 53,  7,192, 60,245,  1,  1,251, 28,  2,  2,  4, 94,  0,255,251,  5,  1,255,  4, 27,  0,255, 19, 30, 11,255,250,  8,254,  0,  9,  8,  2, 11,  9,253,239,249, 19,244,  7,240, 19, 14,238, 33, 43, 32, 54,196, 17, 22, 21, 27,252,  0, 63,  2,254,254, 54,222, 17, 32,249,110,236,  0, 12,  7,  2,242,247,  2,237,251, 15,  3, 21, 13, 10,  1,252,244, 64, 61, 71, 65, 68, 55, 62, 62, 28, 67,235,  7, 78, 16,227,  2,  0,  1,  0,  0,  0,  0,  0,255,  1,  0,  3,  1,  0,  0,  0,255,  7,  6,250,  0,  0,254,  5,  4, 11,  8,252,255,229,255,  0,  1,
        1, 71,201, 21,  9, 59, 18, 46, 22, 62, 65, 22,  2, 50,219, 29, 21, 13, 21,169, 68, 32,212,214,190, 31,251, 30,210, 38,195, 59,  0, 28,  3,  2,254, 11,  2,  1,  4, 46,254,  5,245, 34,  3,  5, 41, 22, 19, 12, 17,  0,253,  6, 16,251, 58,252, 12,243, 79,126, 66, 67, 47, 71, 63, 62,225, 57, 64, 65, 42, 68, 63, 64, 38, 63, 65, 66, 40, 74, 62, 66, 70, 57, 55, 61, 54, 39, 65, 64, 62, 56,246, 10, 82, 28, 18,251, 42, 58,238, 13,174, 71, 56,  2, 47, 24,  2, 42, 19, 63, 70, 62, 60,  8,  0,248,254,255,250, 19,255,252,  0,  7,250,254,  0,255, 19,251, 11,  7,249,  5,  8,  7,253, 46,254,255,230,  5,  3,248, 34, 10,190, 67,185, 59,  7,  6,252, 63,  0,254,  3,253,  0,  4,255,255,254,251,255,  0,  2,  0,255,  0,253, 13, 28, 24, 11, 21, 16, 43,  1,  2,253,  0,  5,251,  1,  1,
       29, 18, 13,  1, 70, 14, 65,255, 24, 15, 64, 13,216,228,  6,249,225,205, 17, 12,  0,114, 37,202, 10, 26, 81,  4,198,233, 47,235,251,250,  0, 10,  4,255,249, 32,  1,253,251,254,  5,  1,255,  1,242,  2,  6,252, 20,127,225,243,250, 48,251,234, 11, 53,255,229,254,253, 12,  6,  1,255,250,  2,  0,  0,235,  5,  0,  0,246,  0,  1,  0,  5,249,  1,  1,  0,253,  7,  0,246,255,  0,  3,  2,255, 70, 74,223,213, 81, 48, 39, 10,223, 60,117,208,245, 24, 10,  2, 81,113, 72,203,235, 12, 73,255, 64, 83, 65, 63, 53, 66, 55, 69, 82, 64,254, 54, 25, 50, 69, 55,  2,255,231, 22,252, 10,  7,223, 35,  4,194,  6,  6,242, 41,210,  1, 65,247,  4,251,254,242,  2, 60, 74, 62, 64, 62, 21, 44, 57, 63, 68, 64, 66, 64, 62, 64, 64, 65, 61, 63, 64, 67, 61, 66, 63, 60, 82, 84, 64, 33, 45, 55, 37,
      226, 46, 27, 22,  3,200,122, 54,  2, 77,241, 16, 64, 12, 38, 24, 54,243,  3,250, 83, 93, 55, 16, 66,251,245,  1, 21, 18,252, 12,250,  4, 48, 48,  1,214,179, 12, 24, 54, 27, 47,  5, 49, 64, 52,  0,  4, 64, 26,254,254, 60, 55,255,  0, 64, 32,  3,  5, 66, 48,255,  0,  3,  0,  4,254,224, 15,  0,254,254,  3,  1,  1,254,  4,254,  1,  1,254,  4,  0,  1,  0,  1,  0,  1,251,255,255,  2,  2, 50,139, 38, 80, 10, 51,  3,  1,190,245,  3,  6,  4, 20,  3,  9,  4,250,200, 15,250,  1, 17,  0,203,212, 30, 77,  5, 49, 26, 29,  8, 62, 67, 54,109,  1,251, 44, 49, 45,247,  4,  3, 13,  1,253, 49, 25,118, 59, 12,254, 36, 12, 81, 89, 84, 28, 59, 50, 39, 29,227,  0,  1,  1,  0,  3, 22, 21,  0,243,  1,  4,  0,  0,  1,  0, 21,  6,255,  3,  4,  3,255,  7,255,  1,253,  0,  2,  6,  5,  2,
    },
    {
        2,245, 14,221, 11, 87, 23,  3,254,211, 17,254,234,232,247,  6,249,210,  0, 17,243,233,  8, 76, 11,241,  4, 72, 16,225,254,239,251,235,232,  5,243, 41,236,255,247, 12,209,  0, 25, 13,250,240,  3, 26,253,252, 58,249,  0,234,204,250,  3,  7,182,160,252, 19,  1,  2,  0,255,  1,  2,  0,  0,  2,  2,  5,  0,  2,  1,  0,  0,253,  3,  0,  0,252,255,255,  1,  8,  5,  0,  2, 20,254,  0,252, 20,242,  6,  3,244, 13, 89,249, 38,183,161,207,196, 46, 68, 11,248, 19,238, 15,221,190, 18,192,238,  6,226,198,  1, 70,178,250,  9,237,  1,  3,  0,  4, 22,249, 12, 20,239, 21,  6, 12, 15, 41,248,  2, 24, 59, 54,240,249, 29, 11,189, 27, 46, 63,  4, 39, 20,  1,  0,  0,  1,  2,  0, 12,250,  1,  0,249,  1,255,254,  2,255,255,  0,249, 36,  7,  3,  5,251,  1,  3,251,  1,  1,255,  0,  3,
        4, 18,195,254,235,165,225,243,220, 12,  2, 52, 45, 40, 29,  3,227,240,184,201,  6, 10,  3,239, 51,195, 12, 12, 12, 80,247,218, 33, 17, 60,238,254, 16,196, 42,252,238, 65, 15,227,238,156,245, 13,  0, 38,250, 32,186,254,  1, 49,  2,  2,224,239, 85,254,255,  0,254,  3,  0,255,  0,252,  0,  0,  0, 11,  1,  1,  1,248,  0,253,  1,255,  2,  2,  3,  0,  4,  7,254,254,  1, 71,  2,255,  0,242,222,  4, 27,  7, 11, 60,244,225, 40, 29, 41, 19,213,176,160, 37,252, 46, 69, 79,  1,140,221,  6, 14,210,248,  5,204, 28,228,  2,230,255,  3,234, 18, 39,183,251,222, 39, 19,251, 36, 10,196, 81,253, 25,178,218,192,248,239, 40, 52,  7,245,205,246,237, 20,  0, 59,  0,255,253,  2,  1, 15,  0,  1,  3,  1,255, 41,255,  3,  5,  0,220,  0,254,  1,  2,  0,  6,  0,255,  1,  8,  0,254,  0,
      255,  2, 26,  0,248, 28, 10,  2, 63, 19,  1,  8, 27, 20, 21,246,234,191,181, 47, 35,  7, 63,250,  8,255,209, 39, 57, 17,255, 41,237,236,189,  7, 11,230, 71, 20,  1,241,194,250,252,247, 67, 45, 20, 20, 78,243, 19,227,255,255,221, 28,251,  6, 37,233,  5,229,  0,  1,  0,  0,255,  0,  3,255,254,  1,  0,255,252,  1,  0,  0,  2,255,254,254,254, 23,  0,  6,  4,242,  0,  3,238,255,255,  0, 20, 53,  7, 35, 47,  0,209,105, 32, 33, 15,243,233, 31,183,230, 39, 12, 33,246,225,232,  2,  7,253, 38,248,189,215,  6, 89,225,192,  6,  8,245, 10,202,252, 53,  7, 54,  7,245, 62, 46, 29, 17,182, 13,219,250, 13, 31, 11,251, 29,  2,240, 12,  2,249, 46, 19,251,  0,255,  6,254,255,253,  2,  0,253,  2,  1,  5,252,  2,  2,  0,255, 25, 28,225,  4,  0,247,  0,  1,  4,255,250,253,  2,253,
      244,250, 45, 24,232,236, 11,243,  4,239, 66,228, 27,240,244, 21, 20, 76, 22,224, 17,246,186,246, 32,204,209,181,240, 70,254,  2,225, 25, 38,  1,243,  1,255, 11,255,  6, 20,244,220,225,  1, 12,237,253,162, 20,209, 19,  9,  8, 42, 57,  0,255,  2,239,  4, 35,  1,  0,  0,  0,255,255,255,  1,  0,  0,  0,  2,  7,  0,255,  0,255,  2,  3,255,244,253,253,254,244,  4,  0,254, 14,251,255,255,244, 22, 27,229, 15,252,252, 60,  1,232, 36, 46, 59, 24, 33,188,  2,239,249,  0, 69,194, 64, 18, 11,254, 10,249, 51,253,152,252, 56, 34,234,245,232,  1, 15,239,250,250, 25, 10, 25,213,241,228,239,250,222,  0, 48,229, 55,241,221, 61,248, 23,116,  1,  0,215,  0,251,254,  4,  0,  0,254,  1,  0,  1,255,249,255,  1,242, 23,  2,255,  3,247,  1,224,  1, 37,  1,  6,  2,  5,244,  2,  1,255,
    },
  },
  {  //  28
    {
       12, 16, 53, 84, 16, 12, 62,209, 23,  4, 31, 70,  6,225,  2,210,  3, 13, 25, 32,224,217,255,237, 12,255, 30, 35,202,220,253,246,253,209,  1,  2,212,208,  4,  6,248,249,254,252,205, 23,  2,241,228,251,255,  9,219,216,255,  2, 75,251,254, 38,188, 90,  0, 48,  0, 10,254,255,255, 12,  0,  0,255,  2,  0,  1, 11,  4,255,  1,  5,  1,  3,253, 28,  1,255,255,239,  2,  1,  3,168,255,254,254,  6, 23,253, 46, 16, 14, 81,209, 23, 25,235,  4,205, 40,244,229, 18, 64,246,244,197, 49,246,111,240,243,  6,  1,209,249,231,114, 13,250, 24,  7,  7,  6,251, 53, 22,245, 10, 28,255, 19,  5,235,246,  8, 44, 26, 36,  7, 76,243,198,255,  4,255, 53, 66, 37,  8,  2,254,  3,253,  3,  0,251,  1,  2, 46,  3,  4,250,  0,  1,254,254,  2, 34,  7,255,  1, 50,251, 19,  1,  2,255, 26, 23,250,254,
       23, 27,248,197,237,  1,217, 34,218, 51,  3,240,251,240, 31,171,221,255, 18,234,  0,226,187, 34,228, 26, 45,  3,255, 44,212,226,218, 46,254, 61, 35, 28, 12,218, 38,231,247, 11, 34,225,  0,255,  5, 23,248,246, 50,  8,211,  8,218,249, 21, 24,246, 60,250,252,254,  1,  0,  3,  0,255,255,  1,251,  0,254,  2,223,  1,255,  1,253,254,  1,255,  1,  4,  0,  1,  2, 34,255,  1,255,246,255,  1,252, 12,251, 14,212,191,201,200,240, 31,214,216,  4,164,173,251,  5,197,250,213,253,245, 52,175,188,180,246,  8,197,245,181, 48, 33,250,234,251,240, 10, 22,231,254, 27,238,242,  0,  0,239,  0, 41,  2,242, 16,221,251,204,235, 20,242,188,245,245,233,254,218,245,  3,  3,  2,243,  3,  8,  2,253,250,254,252,255,  1,252,254,  1,  0,252, 10,  1,120,  2,  0,250,  4,  0,  3,  0,248,  8,253,
      225, 18,173, 87, 13,241,254, 15, 89,242,238,236, 27,239,  3, 67,247,238,137,251,247,215, 66,242,252, 61,231, 74,  5,249,252, 29,  3, 33,238, 21, 66,  0, 27,194,231, 17,  3,246,252,245,  1,197, 52,  9,240,228,235,  1,253, 58,254, 66,253,254,210, 50,230,  3,  0,  1,255,255,  1,254,  1,  0,  7,  0,  0,255, 36,  3,  0,254,251,  1,  1,250,217,253,255,  2,248,  3,255,  2,253,244,  2,  4, 22, 13,249,238,117,243, 77,  4, 17, 11, 19, 53,243,225,241,207,  0,251,241,246,  4, 37,137,192, 65,  2,253, 70, 61,238, 30,132,201, 53,208,226,252,252, 69, 24,  9, 11,247,235,  6,255,  5,246, 56,204,  1,  2, 84,189, 80, 24,101,196,251,136,  0,202,231,135,  0,  5, 29,237,  2,252,  2,  2,  0,  1,  0,254, 15, 10,  4,  1,  4,224,  1, 45,  5, 87,253,253,255, 39,  6,  1,253,239,239,255,
       39,240,174, 80,  4,  0, 31, 72,  3,  1,  2,213, 37, 70,  9,  1,252, 19,117, 98, 41, 36, 88,252,239,241, 64, 64,245,  0, 14,133,252,246, 16,251,201, 31,  1,244, 48, 12, 18,195, 50,  6,255, 64, 57, 31,253,253, 45, 41,249,  2,  0,214,  0,236,252,255,  0,  9,254,254,  0,255,255,  1,  0,  1,  1,  0,  1,  4,250,  2,  1,253,  0,  9,  1,  0, 15,187,246, 11,  1,255,  0,  0,249,  2,  2,  0,254,  2,  4, 13,143, 10,  4, 83,254,  4,238,232, 31,112,253,179,247, 39, 20,250,204, 54,  1,196,250, 25, 37,  1, 59,198, 80,201, 53,235, 50, 54,219, 28, 77, 60, 16,252,221,253, 62,216,162,  2, 15, 59,208, 67,206, 61,121,254,190, 64, 76,146, 10, 40, 61,209,  0,234,241,228,227,  4,255, 11,255,  7,239,  5,254,  1,249, 65,  8, 10,  4,  3,254,193,252,  3,251,240,  2,  1,249, 11,197,  4,
    },
    {
       17, 21, 15, 15, 21, 82, 14, 28, 32, 29, 51,184, 55, 13, 45,168,232,  3, 12, 93,235, 46, 49, 27, 25, 45, 64, 27,105,233,217,197,206,  2,178,  4,255,218, 54,248,204, 10,135, 13, 13,  3,142, 15,251,  2, 64,247, 15, 22, 51,247, 10,232, 42,253, 23,  7,235, 57,  0,  1,  6,255,255,  3,234,  2,  1,  4,  0,255,251,  2,  1,  2,  1,255,  2,  3,  0,  0,251,  2,255, 27, 10,  0, 13,  0,253, 24,255,248,  0,  1,245,  7,251,251,248, 22,254, 51, 27, 24, 64, 45, 98,254, 20,  0,144, 57,  5,229, 53,199,216,239,223, 68,213, 55, 10, 14, 10,  0,  0,  3,  0,  2,  8,250,238,244,  1,252,255,  5,  3,  0, 61, 27, 62,  2,  3, 80,  7, 18,199,226, 15,100,230,197,  4, 67,  0,  2,  1,255,253,  3,254,  4,254,253,  0,  1,  0,  1,  1,  2,248, 24,  1,  0,251,248,  0,255,251,250,  4,251,255,253,
        4, 15, 37,191, 26,206,218,240,226, 64,162,157, 81,221, 30,108,200,223,  1,240, 72,208,194,254, 14,242,151, 58, 21,196, 91,156, 10,  8,  0,246,195, 18,198,244,247,236,185,245,191,203,255, 23,  0,  0,  6,197, 61,  1,239, 75,191, 75,104,255, 51,219,253, 16,253,  0,  2, 33,248,  0,  1,  1,  1,  1, 36,255,  8,  1,  0,  2,  0,255, 21,235,  1,255,  0,  0,255,  0,251,253,255,  2,  2,  4, 11,228,253,249,  7, 16,  2,  0, 23,213, 12, 51,205,245,  0,205,  1, 71,211,116,203,187, 67, 57,  0,  7, 58,167,252,  4,186,227,  1, 33,  2,255,  4,249,  0,  0, 24, 16,242,178,250,132,203, 17,251, 14,209,229, 63,109, 27,220, 22,255,249, 57,121, 22,  7,191,225,235,  1,  0,255,  0,254,  0,254,  0,  0,  0,  1,  7,  0,  0,226,251, 36,248,247,  0,  1,253,  0,  4,  0,  1,  9, 12,  1,  1,
      253,241,255, 66,249,209, 19, 21, 75,232, 74,196,250,255,209,227, 63,197,217,235, 13,227, 62, 13,245,249,124,118, 73, 96,168,230, 51,255, 57,192, 37,  0,247, 18,255, 35,  5,197,187,  6,252,  0,  9,  1,  5, 61,229,253,181, 21, 25,217,  4,139, 78,158, 13, 23, 22,  0,  2,253,255,  0,  0,  0,255,252,  2, 28,  0,  1,255,  0,255,  0,245,  0,248,255,254,251,  0,  0,  1,  0, 29, 14,250,240,  3,  7,  7,  5,255,255,  0,250,220, 57, 53,  2,  0,195,  2,  1,159,247,218,130,193,  6,197,158,  2,194, 20,  5,135, 59,105,199,  6,  4,252,  6,  0,  1,  2,  2,205, 20, 55, 15, 13,255, 67,135, 73,  0, 30,233,240,121,233,187,  2, 30,233,219,251,200,222,  7,  0,191,  0,  2,255,255,255,252,  0,  6,  1,  0,  0,  0,  1,254,  0,  0, 16, 27,222,217,255,  0,  0,  9,254,255,248,  1,255,  0,
       18, 20,  4, 59,217, 54,229, 15,254,242,233, 61, 82, 19, 70, 55, 13,203,  7, 60, 88,200, 37,  0,  0, 68,117,235, 63, 38,151,233,251,253,253, 59, 33, 28,248, 16,242,  4, 52,116,221,247,125,252,198, 22,233, 41,  4,255,197,201,246,255, 15,253, 11, 72, 19,121,243,  3,  1,  2,  1,  0,255,255,251,  0,  3,  0,253,255,  2,  0, 10,  2,254,249,254,255,254,  0,  0, 16,  6,255,132,  0,  2,  1, 14,251,  5,  3,255,  1,255,255, 72,249, 23,  8, 23,187,  3,123,102,251,254, 15,251, 17,116,253, 15, 54,249,252, 49, 25, 19, 91,  4,255, 59,254, 69,249,254,  3,  1,211,209, 62, 76,116, 36,214, 67, 15,  0, 18,228, 24,251,253,230,  0,211,167, 44,186,193, 24,  0,251,231,245,244,  2,193,252,254,  1,  1,  0,232,  0,  0,  7,  2,255,  0,255,  1,246,247,255,  4,255,  2,  0,255,249,  2,254,
    },
    {
       21,  5, 71, 17, 12,248, 73,  1, 78,  7,252,230,243, 28, 10, 38,  4, 57,189, 12, 43,  6, 28,  7, 28, 45, 42,206,  8, 18, 22, 23, 19,215,243,253, 32, 17,243,  0,233,241,  4,247,  9,  1,249,  3,253, 25,159,  0,233, 12,244,254,240,232,  1,254,255,216,  0,  1,  0, 13,255,  0,255,251,  1,  1,255,255,  1,  1,  2,255,  0,254,  1,253,  0,  0,  4,  4,  1,  0, 26, 39,  1,254, 43, 23,  6,  1, 37, 37,253,253,230,220, 85,190,222,212,233,199, 31,232,234, 32, 14,220,243,250,212, 21,210,245,229,207,233,229, 31,246,247,244, 26, 24, 11,247, 44, 55, 66,  4, 30,  5, 15,246,  0, 58, 57,187,249,  3, 41,194, 41,186, 55,  6, 11,255,182,125,  6,153,  4,181,219,  2,  1, 12,  2,  0,  8,231,  1,  4,  1,251, 15,251,246,231,  3,252, 30, 20,  5, 25,  8,  1,254,254, 12,  0,243,254,  4, 13,
        3, 44, 22,245,232,  8,199,231, 27,  5,223, 34,219,243,249,167,  2, 52,  4,102,230, 52,241,  4,222,  5,180,219,214,201,  2,139,230, 67,254,  1,239,236, 13,  4,242,  2,  2,  7, 23, 14,255,254,239, 16,255,253, 31,227, 46,250, 35,  2,  0,  6,239,158,253,243,251,254,  1,  0,  1,  2,  0,  0,250,  3,  0,253,240,  7,  0,  1,  0,246,  0,255,255,  9,  2,  2,241,229,253,  0,245,252,255,255,229, 38,208, 48,255,254,227,  4,  4,241,207, 45,  8,243,243,244,  1,  1,233, 10,212,248,  1,194,200, 28, 64,240,206,130, 61,  8,210, 61,  9,235,234,230,166,  2,254,231, 35,  5,252,255,209, 25,  0,  1,  0, 36, 18, 38,253,237,232,128, 97, 22,  1, 23,  8, 94, 30, 22,  7, 25,  6, 21,251,255,  1,  1,  6,252,206,252,252,180,253,  1,  3,255,  2,  5,230,  2, 10,  0,254,  4,231,  1, 43, 10,
        2,242,226,244, 52, 13, 30,237,245,242, 41, 55,252, 37,251,219, 29, 15,171, 44, 22,231,253,233,234,254,200, 52,251,247,  1,253, 22,202,  7,255,  9,239,246,255, 50, 48,  0,254, 26,222,  3,  1,234,  9,  3,252, 38, 26,110,255,184,217,  7,254, 19,223,  0,  7,  3,254,  1,  0,255,  1,  0,  0,  5,237,  0,  1, 17,247,  1,  1,  1,249,255,  1,  1,232,  1,255,133,250,  0,  1,238,  2,  0,  0,227,198, 34,236,241,  9, 80,166, 23,236, 15, 17, 17, 15,252,187,183,  5,246, 25,243,253, 14, 68, 63,223, 57, 23,121,254,  2, 15, 86,  2,242,  8,222, 16,253,195, 37,  3, 32,251, 19, 55,245,249, 12,  9,156,226,195,224,  7,248,255, 58, 55,240,247,  0, 21, 35,244, 20,228,248,184, 12, 10, 63,238,  5, 74, 25,230,129,161, 17, 11,244,244,250, 10, 13, 26,131,  1,  1,255,  2,  2,  2,218,  3,
      238,217, 56, 45,240,253, 15, 11,  5, 16,185, 12, 13,227,  0,  0,253, 13,240, 26,254, 12, 15,247,239, 33, 43,242,227,132,  1, 16, 26, 46,250,253,  0, 40,252,252, 17,233,  5,  1,237, 35,  5,252, 27,254,  1,252,238, 15, 19,246,223,230,248,255,241,225,  2,  6,253,  2,  1,  0,254,253,252,  1,254,252,  0,  0,  8,  1,  3,  1,250,  4,255,254,255, 10,255,  0, 10,238,253,254,251,  4,  3,  1,250, 21,255,247, 16,  2, 62,244,249, 13,251, 95,  8, 40,204,150, 52,250,246,191, 13,160, 71, 10,251,250,253,251, 59, 62, 72,196,204,  7,238, 18,215,195, 90, 18,236, 53,255, 26,235, 71, 76,201,228,246, 45,  4,  1, 76,189,204,  1, 11,  0, 16,200,235, 59, 28,253, 16,225, 26, 12,251,199,  6,  5,  1,180,214, 69,246, 17,129,254, 11,242,  8,254, 11, 14, 43,  1,  1,  3,250,249, 77,242,239,
    },
  },
  {  //  29
    {
       15, 29, 19, 34,222,255,205, 28,242,243,  0,232, 59,  6,131, 97,239, 13,221, 21,242,218, 31, 67,254,238, 25,  1,155, 15,215, 48,  1,202,255,  6,  1,205, 11,255,234, 38,255,250,  0,200,  3,250,249, 12,247,  3,249,245,  1,254,208,183,253,  0, 29, 30,  0,252,255,  3,255,255,  2,243,  0,  0, 14,  5,255,255,254,  2,254,  0,  4,  3,255,  0, 46,253,255,  0, 27,240,  1,254, 79,157,255,253,  0, 33,225, 14,245,204,239,149, 50, 12,243, 43, 44, 25, 50,222,  8,251,247,228, 58, 49, 56, 70,252, 72, 28,206,238, 44,  2, 69,  9, 25,240,234,219, 15,130, 35, 15,222, 25, 12,253,253,251,189,  6,255,202, 68,  8,247,251,193,  6,  5,222, 24,239, 50,204,232, 49, 43,  3,255,224,  5, 61,243,  0,250,  2,  1, 14,  2,  3, 86,  2,  0,246,255,  3,252,  1,240,  1,250,  0,254,250,252,  8,200,
       18,251, 21, 20, 62, 48, 25, 24,218, 49,192, 14,236, 54,204,184, 13, 13, 92, 27,202,219, 11, 86,202,  1,241,142,246, 22, 85,  7,249,231,  3,  0, 20,  2,242,  5,234, 18,  8,  0,226, 21,  1,253,255,  1,254,253, 41,226,253,254,  4, 14,238,245,  8, 41,119,251,239,  1,  0,255,252,255,255,255,252,  1,  2,  0,  0,222,  1,255,  3,  0,  3,255,232,251,  2,  1,  2,  9,255,  1, 16,  4,250,252, 18,236, 29,225,231,197, 32, 94,238, 15,236, 66,177,238,177,194,  4,  4,228, 53,228,225, 75,  3, 14, 60,251,248,225,192,129,130,  3,228, 41, 40,  6,  5, 21,220, 42,  3,224,230,218,239,138,156,196,195,  5,  7, 18, 92, 25, 57,  6, 23,218, 37,237,222,142, 46,  2,  3,239,  1,244,  1,114,147,  2,  1,  3,255,  4,  8,253,236,253,203,246,  3,252,  8,251,  0,  1,254,249,  1, 44,255,230,234,
        0,  1, 19,204,229,  1,248, 16,  0, 23, 79, 43,223,243,203,237,  1,  5, 17, 40,252,  6,  3, 43,  8, 42,247, 19, 71,236,  8, 65, 46,231,  6,  1,244,223,255,  0, 37, 36,  3,247,  9, 14,  8,255,  0,225,  1,219,229, 52,  0,254, 49, 32,136,252, 61, 65,253, 40,248,  1,  2,  1,  1,255,  0,  1,250,  0,  1,  0, 59, 11,  0,  1,186,251,255,249, 13, 70,252,251,  3,  1,  3,  0,  6,254,  2, 13,216, 55,  3, 31, 15,  3, 57,101,253,  3, 27, 16,251, 42,133,  4,215,237, 58,165, 68,205,191, 69,  0,  6,  1, 50, 57, 29,196,250,237, 27,254,197, 17,  6, 29,215,  8,251,  7, 55,248,175,219,  1, 49, 60, 61, 48, 62, 46,204,156,229,217,254,234, 47,226,126, 96,219,  1,205,  2,245, 10,247, 39,251,  2,248,  5,  3,  7, 16,128,  5,  4,254,241,  5,254,  2, 13,251,  4,252,  1,255,  1,254,137,
      252, 13,249,168,239,255,235,  0,188,250,231,154, 24,  1,240,231,230,255,242,152, 22,251, 50,201, 34,  0, 57, 23,161,252,  1,  5,208,116,  1,254, 32, 47,240,254,  2,235,  3,  1, 23,209,254,  0,203,  6,252, 46,230, 38, 31,233,244, 30, 42,  3, 37,  1,236, 64,198,255,255,  1,  0,  0,253,  0,  4,  9,  1,255,212,  1,  0,  0,250,  1,  1,  3,255,250,  2,251,189, 10,241,  0,250,232,253,252,  2, 99, 23,245,  1, 48, 62,117, 16,229,213,  7,235, 57, 33,210, 35,223,243,240, 20,  1,  9,207, 55,237,200,160, 55,216,189,  8,  5,251, 88, 14,253, 75,181,  2,249,255, 23,230,  8, 44,248,  6,240,  3,253,233, 61, 90, 68, 67,255,254,192, 46,249,208,247, 90, 40,  3,207,  2,250, 63,181,146,  4,  2,252, 30,  5,226,249,  0,  0,254, 19, 19,240,253,253,142,252,  3,251,  2,  2,  0,243,167,
    },
    {
      194,207,  0,109,250,255,253,124,187,193,251,116,250,  2,250,249,188,196,193, 89,250,  0,253,  6,183,192,193, 66,246,254,  2,254,222,195,219, 23,228,194,235,  3,240,186,224,248,  2, 13, 11, 12,227, 24,246, 14,  2,  8,  4,  0,233,229,  4,234,  1, 15, 77,253,255,251,  1,  3,255,  2,253,255,  0,  0,  0,  0,  4,  1,  2,  2, 19,241,  2,252, 58, 10,255,  2, 29, 21,  0,249,225, 42, 37,  2, 36,237,240,203,240,156,241,103,  5,  8, 25,250,251,241,  5, 11,240,252, 22,244,235,216,193,197,  5,  7, 33, 17, 34,107,  2,194,247, 67,235,228,208,245,214, 10,217,198,207,192,222,  4,168,232,250,222,241,  4, 13,194,  2,254,173,232,255,215, 24,241,255,250,255,  2,242,  5,185,252, 60,131,  2,247,172,208,211,254, 10,  3,  7,250,249,248,  8,255,255,254,246,244,  7,  0,253,  3,253,  0,
        3,252,250,  1,  1,  0,255,252, 10,251,202,201,  6,  0,  7,251,251,253, 22,183,  1,254,242,155,255,253,254,131,255,253,  2,  1,208,232, 31,205,247,214,244,253,208,  1,148, 42,  0,  3, 14,243,248,227,192,  5,255,220,241, 11,222,235,125,255,247,244,196,237,240,  1,  0,248,246,241,  1,  7,246,  9,235,255,242,237,255,243,  1,  4,  1,  3,199, 24,  2,252,  8, 13,255,  0, 61,225,  8,249,221,254,220,213,226,254,214,175,222,  7,201,  5,242,245,196,  5,203, 33, 12, 12, 62,216,  2, 65,255,238,242,  0,  4,162,  1,  7,234,172, 23,  2,245,195,  3,  2,234,249,252,253,245,237,254,  1,212,217,188,226, 16,202,252,253, 29,177,194,224,192,243,153,209,239,255, 21,247,249,253, 64,245, 11,218,  5, 10,241, 17, 11,253,  0,223, 14, 12,252,246,  5,  3,  3,249,253,  2,  0,250,255,  2,
      221,192,185,210,198,197,198,190,206,199,191,200,185,191,184,192,194,199,192,199,196,196,191,  4,211,196,203, 56,193,197,235,197,  7,241,244,  2, 27, 22, 12,243,232,252, 33,251,246,245,246,245,  5,249,255,252,227, 21,  3,254,244,199,  4,246,222,  5, 10,247,189,205,184,188,203,208,202,189,193,203,193,197,210,211,201,222,221,200,191,192,172,  9,186,190,190,191,191,192,194,192,200,186,239, 26, 66, 75,  0,241, 14, 36, 69, 29,132,  4,212,  1,240, 25,242,226, 96, 11,  9,  2,207,185,187, 79,  0,  2, 16,250,224,235,251, 48,  0, 72,242,  3,247,199,255,251,  5,255,250,249,  3,  9,255,196,191,  4,211,251,  2,204,255,192,181,  6,246,176, 44, 35,231,188,192,224,241,225,253,  9,  7,255,141,168,212,189,211,239, 62,190,203,186,209,204,195,200,185,193,190,188,186,199,191,195,
      251,254,246,250,246,180,  3,254,  0,197,  4,252,  2,149,  8,245,  0,203,250,255,  0, 32,247,252,248,253,157,252,255,  0,  7,  2,192,  7,238,237,197,190,237,231,193, 57,  2,  2, 34,180, 61,196,196,174,204,200,206,255,218,230,193, 55,242,254,233, 10,  9,254,  0,254,  9,  4,247,  0,254,235,  0,255,  0,  5,249,254,253,237, 67,  6,253,252,236,255,  4,255,  2,251,  2,249,  1,249,250,247,127,244, 36,  9,244,211,227,212,  2,239,  6,220,208,208,229,235,207,254,231,182,  0,254,255, 37, 41,252,203,192,241,252,  9,196, 96, 50,  0,225,251,227,199,  7,139,  2,  8,231, 27, 20,203,244,  2,194,  1,253, 59,216,188, 16,195, 73,246, 58,  6,  4,  9, 28,  2,247,220, 77,203,  3,247,225, 14,254,176, 84, 10,235,254,230,253,200,144,254, 60,151,227,  7, 23,230, 27,245, 12,  8,  2,226,
    },
    {
      222,244,238,250,221,246,250,253,237,240,241,253,224,244,231,247,  1,241,250, 15, 15, 39,  1,246, 11,212,248,244,238,225, 12,240,227,  4,240,253,246,  4, 12,251,235, 10,235,249,238, 26,246,  0,242, 39,253,251,229, 19,254,253,223, 57,252,  4,212, 51,  6,  0,  3,  0,  0,255,  1,  0,255,  1,  5,  2,  0,  0,  4,  0,  0,  0, 13,  5,255,  0, 17,253,  2,  1, 16,244,255,  2, 55, 12,255,  1, 20, 26, 16, 36,  2, 77, 83, 53,  7, 34, 34, 29,215, 92,232, 28, 19,  7,203, 48,207, 13,231, 52, 11, 67,  1, 15, 23,191, 11,192,217,255,203,234,204, 27, 23,226,  3, 30,251, 56,  0, 44, 41,181,243,  6,249, 15, 30,104, 21, 52,205,231,  3,247,251, 32,  8,100,  1,  0, 19,  1,  0,241, 14,255,  0,  1,  4,244,  5,  2,253,  8,  1,  0,  1, 11, 10,221,  8,  5,254,251,  1,252,  0,247,  2,  1,
      253,247,  7,  7, 27,240,253,254,  1,254,187,  5, 47,204,195,254,207, 65,254,255,197, 54, 68,230, 59, 20,222,235, 89,255, 17, 45, 35,253, 46,249, 15,245,231,253, 40,246,252,  2, 54,214, 15,  2, 81, 19,  0,  0,255,255,244,254, 10,209, 56,  3,206,247,  1,255,  0,  0,252,255,  0,254,  0,  0,  0,253,255,255,  7,  0,  4,255,237,  1,  1,  0,250,199,  3,  0,  2,244,  3,255,242,251,255,  0, 33, 13,  0,  2, 92,242, 38, 56,  1, 58,  7, 14, 92, 88,254,193, 27,  1, 15, 55, 31, 13, 79, 18, 26, 67, 45,  8, 42,  7, 68, 85, 60,  2,245,202,239, 58,229,190,255,234, 27,190,235,219, 23, 51,245,  0, 29, 35, 64,  8,248, 10,  4, 83,227,255,211,  2,  7,226,  1,255, 17, 63,  4,248,237,  5,  5,  2,  6, 43,  6,254,245,252,  0,  0,  0,251,250,  7,  5,151,  4,255,249,254,252,  2,  5,  4,
       44,  3,229,254,227,232,  2,  1, 22, 51,206,255,  3,238, 67,248, 12,  7,215,  7, 14,191, 28,  3, 38,255,243,144,220, 29,180,255,  0,246, 61,254,239, 10, 19,250,254,195,225,254,219, 12, 44,250,188,226, 94,253, 22, 11,236,255,245,  0,  0,  0, 27, 53,221,  8,  0,  1, 17,  0,  0,  0,240,  1,  2,254,255,  0,  5,  1,  8,255,  1,255,  4,  0,255,  3,255,  1,241, 56,254,  0, 10,254,236,  1,  8,249, 24,  3,213, 63,229,146, 64,252,246,  5,166, 20, 43, 67,237, 19, 72,  0, 70,234, 90, 11,201, 14,255,248,125,194,201,103,  1,  2,247,240,202,244,214, 59,250, 14, 71, 33,  1, 79, 57,250,249,250, 58, 17,240,166, 56, 79,  7,190,182,219,  6,210, 10,203,255,192,  1,208,246,  1,248,252,  2,255,204, 40,  2,188,  5,  0,  2,255,  1,234,  0,247,255,197,  0,253,  6,245,255,254,254, 23,
      255, 23, 65,252, 52,222,  5,  0, 61,  2,  3,  0, 18,235, 37,  0,  0,248, 22,248,  4, 52, 39,  7,  5,177,241,  6, 22, 14, 65,247,254, 23,214,237,  0,  1,203, 11,235, 47, 67,255,245,  4,212,  0, 42,  9, 23,  1,248,  5,159,  2, 12,  0,253,255,254, 37,252,  1,  0,255,252,  0,  0,255,255,  2, 15,  3,  0,  0,  5,  2,  8,255,  6,  0,  0,  1,254,  4,  2,255,234, 11,253,  1,  8,  4,254,  1,247,215,  0, 14,190,  3,  0,233,191,222,248,253, 30,250,219,232,105,  0,  2, 36,  6,  5,165, 78, 15,213,207,255, 15,178, 18,  3,255,  9, 69, 57,  1,187,  9,  8,  1, 40,187, 45, 29,255,239,254,227,  0,222,187,196,190,251, 75, 11,228, 43, 85, 57,198,242,202,  0, 58,  1,226,242,247,  8,239, 11,  1,232, 36,  2, 46, 85,  4,  0,  0, 13, 22,  0,  7, 68,212,  1,  2,209, 13,  4, 33,242,223,
    },
  },
  {  //  30
    {
      235, 29,235,236, 32, 12, 80,227, 16, 67, 36,118, 62, 53, 51, 46,  3,236,236,253, 21, 92,  2,235, 36,246,  5,255, 66,100,200,214,229,  3,237,251, 39,250,239,  0,  7,235,173, 10,247, 28,214, 56, 19,250,244,  2, 10,245,  0,245,232,  1, 68,254, 46, 99, 46, 15,  6,  0,  1,253,  1,  2,  0,  0,255,  7,254,  0,  0,  0,255, 14,  1,  2,  4,  0,253,  4,  0,254,  0,255, 29,  0,  7,  9, 38,  1,220,249, 48,242,  6,216,236, 48,  9,193,255, 61, 34,200, 83, 67, 21, 23,  9,248,227, 78,  7, 47,135,247, 53, 64,166, 43,236, 34, 16,253,248,249,251,255,239,  5,216, 49,186,221,215, 10,177, 66,  5,  4,192,189,  3,234,  8,238,255,  1,223,217,254,215, 15,202,  0,  2,  0,255,  6,255,  5,  5,254,  2,  3,  0,253,251,  3,254,  0,254,  4,243,255,245,254,  1,  0,255,  0,  1,  1,  5,255,254,
      231,  4, 47,  0,221,  1,226, 13,  4, 40, 13,252, 15, 32,245,163,220,243,255,  0, 35,205,245, 62,195, 19, 38,  6,  3,165, 79, 56, 21,  2, 82,  2, 13, 14, 14, 14,  9, 57,233, 57,223,234,251,220,255,252,200, 24,213,213,  1,151,253,  3,189, 33, 11,115,  1,218,  0,  2,  6,  2,  0,  4,  3,  0,  0,246,  1,  1,255,  0,255,  0,250,255,187,255, 37,  1,  0,  0, 27,  4,  0,  0,248,249,255, 14,  7, 10,246, 32, 10,179,  1,217, 34, 65,254,199,191, 75,185,184,230,183,254,  3, 15,174, 10,136, 69, 75,217, 25,108, 12,223,151,255,242, 21,229, 58, 40,250,249, 48,241,  5,220, 13, 19,234,250,  5,  0, 54, 58,249,  4, 70, 28,255,249, 46, 20, 21, 42, 82, 52,  0,  1,250,  6,  2,  0,  0,  0,  2,196,227,252,  9,255,  0,253,  0,255,  0,  0,255,  1,  1,  1,253,255,241,244,  3,251,255,243,
       84,250,223,  2,  2, 34,230,243, 41,235,243, 94,226,241,194,254, 20,229,169,  1,243,  9,  5,252, 50, 33,251,248, 42,255, 23,194,246,  5,213, 49,248, 35,229, 57, 11, 24, 76,234,234,  7,125,  3,  0,  9,252,132, 20,  3, 63,162, 72,200,249, 11,219,254,213,  0,255,  0,  0,244,255,255,  0,253,  2,  1,254,250,  8,  1,252,255,247,  6,  0,  0,  0,  0,  4,  2,247,  4,235, 18,  3,247,252,255,  9,  5, 18,253, 55, 36, 72,  1,  5,  7,240,248, 10,  5,175,246, 16,100,102, 65,205, 11, 60,208,244,235,252,  0,130, 57, 65,114,234, 18,198,235, 65,  7, 31,237,242,  6, 72, 37, 67, 30,171,189,246,  3,254,240,  5, 87, 73,223, 27, 72,  3,229,150,211,187,147,255,  1, 57,  0,  3,  6, 20,255,254,248, 53,248,  7,254,  1,  5,  1,  0,  1,255,254,254, 22,253,  9,  8,246, 20,  1,  4,  0,254,
      238,  3,224, 44,255,249,254,248,224,246,  1,  1,240,243,  3, 13, 46,  2,196,255,  3,247,240,247, 31,  6, 15,255,254,132, 82, 29,  0,  5, 56,189,254,151,245,190, 12,  7,237, 30,215, 27,148, 54, 32, 17, 53,198,245,238,132, 31,186, 16,230,  3,239,206, 47, 13,  0,  0,  0,  1,  1,  0,  0,  0,  0,254,  0,  4,249,253,255,  0,  1,253,255,  1,  6,  1,  1,253,  6,250, 35, 10,  5,253,  4,  0, 31,  4,203,  0,247,  7,230,182,  2, 24,  7, 45, 42,236, 27, 89,  8, 21,213,208, 84,204, 78, 96,253,201,255, 60,250,164, 83, 53, 33, 10, 37,233,186,  9, 31,  3,  5, 64,244,234,237,186,149, 22,  0,249,  1, 69,190, 42,167, 93,224, 33,  2,233,  7,241,195,154,  3,255,227,  0,  2,248,255, 35,  1,  0,203,199,208,  1,242,255,  0,  2,  3,  3,  4,254,  2,  1,250,253, 28,  2,  1, 28,222,  2,
    },
    {
       13, 62, 25, 30, 16, 14,  9, 33, 14, 13,118, 11, 13,  2,134, 41, 19, 12,233, 17, 19, 32,  5, 40, 20, 31, 64, 17, 30,  0,  7,131, 10,229,252,255,  0,216,252, 11,  4,  1,220, 57,254,169,151,219,246, 28,240,238,251, 20, 38, 27,248,220,  9,  3,245, 39, 23, 66,  0,255,  1,255,255,  2,  1,255,  0,  0,  4,  6,  1,  0,  0,252,  2,  1,255,255,255,  1,  3,  3,  3, 56,  7,  0,  8,255,  1,  8, 26, 40,  9, 27, 38,221,219,229, 51,248, 36,  6, 19, 36,249,  8,  6,  7,233, 57, 13, 17,246,240, 11, 35, 88, 53, 27,255, 39,  1,240,  3,  5,  1,  7,  9,253,251, 80,  2, 14,  6,  9,  2,230,  7, 59,247,249,210,197, 50, 30, 69,255,247,211, 51,246, 63,232,  2,  0,  5,  0,  0,  2,253,202, 10,250,  2,252,  2,255,  0, 10,227,245,  1,  2,252, 11,  3,  1,  0,  1,  4,231,255,254,218,  4, 15,
       13,213,182,247,240, 13, 84,242, 12,253,232, 49,  0, 27, 44,201,  3,  1,241,181,242,226,193,210,255, 16, 22,234,  2,  4,125,  4,241,232,250,  0,252, 41,244,253,  3,201, 13,202,  2, 28, 79, 27,251,196,  9,251,  1, 39,190,217,  1, 32,233, 18,255,178, 40,200,  0,226,  1,255,  3,  1,  0,  1,  1, 23,255,  0,  0, 40,  4,  1,  4,242,255,  1,255,255,  0,  1,  3,214,  3,253,  1,254,255,  2,  3,253, 31, 23, 28, 11,226,233, 33, 23, 25, 57,208, 66, 56,231,230,  1, 13, 43,247,230,  1,202,  5, 15,186, 27, 17,254,248,194,226,251,249,  6,255,  3,216,  6,228,203, 49,226,239,237,  9, 40,178,190,218,  8, 14,  6, 12,171, 15,247,244,251,201,255, 62, 72,251,235,217,  5,  3,255, 61, 11,  2, 21,255,251,  7,244, 16,253,251,  0,  9,255,  7, 24,  8,252,  0,  0,255,  1,246,254,  1, 38,
      253,243, 50,176,249, 61,246, 78,255, 77,  5, 57, 31,  3, 48,240,  1,  2, 29,208, 39,251,231, 32,  3, 14,  8,212,253,214,209,  2,  7, 70,253,  2,250, 62,252,220,  4, 58, 47,  1,  2, 35,241,  0,  4, 67, 21,241,254,189,  3,250,  3,251,  6,246, 18,  1,  0,250,252,  1,  1,  1,253,  6,  0,254,  0,243,254,134,  4,  2,255,255,255,  6,  1,  1,  1,  7,255,255,  0, 64,255,243,244,255,254,242,204,187, 13, 26,208,197, 25,255, 59,  2, 19,238,  1, 54,154, 12, 17, 44,  8,202,244,  5,121,229,  6, 34,254,204,246, 39,252,  5, 37,216,228,251,251,206,244, 69,253, 14,181,252, 43,190, 24,  0,249,254,221,  1, 60,252,246, 26,255,  2,242,253, 16,253,252,186,244, 34,231,  7,252,184, 12,254,  0,  6,  4,226,  2,233,200,245, 75, 11,224,  1,234,242,250,199,  1,  1,  1,255,  2,  1,  2, 18,
       54, 52,194, 80,254, 34,202,189,252, 10,240, 32,242,213,251,227, 14,255, 82, 16,199,  6, 52,207,239,255,176, 42,205,253, 98,  3,244,  6,227,252,251,243, 17,  1,252,  6, 32,  4,  0,234,  3,  9,  0,194,  4,245,  8,  9,252,  1,  0,  3,  0, 17, 25, 56,255,251,  0,254,255,254,  0,  4,  9,  0,  1,252,254,133,  0,255,  1,  2,254,  2,  1,  0,254,  5,  5,  0,  0,  1,  0,250, 55,221,  0,  0,255,208,217,254,244,182,250, 55,  4,220, 26,236,  0,  7, 49,124,243, 31,249,  9,196, 19,  5, 89,  1,  1,255,172, 52, 51,247,232,  4,245,  5,216,108,  2,157, 67,233, 14, 23,  9,253, 32,249, 30,230, 38, 38,  2, 35,204,255,211,  5,  0,239,248, 13,  0,  1,215,208,  6, 54,  3,  2, 60,241,250,  0,  1, 62,  4,206,  6,  5, 24,196,251, 19,254,216,230,252,221,  0,255,255,254,  1,  3,255,178,
    },
    {
       15, 16, 63,196,  8,  3,  0, 58,  8,  4, 11,247,  6,  4,249,249,251,  2,  1,254,251,248,  0,251,  9,  0,  3,254,240,  0, 10,129,243,  2,  8,  4,  3,  8,  1, 11,  2,236, 67,190,251,242,253,184,249,238,185,244, 10,  7,253,178, 10,238,196,  4, 13,246, 18, 62,255,  5,  0,253,  0,  0,  0,  2,254,  1,242,255,254,  0,255,  4,  9, 23,250, 15, 71,255,  9,  5,248, 42,231,252,  1,  4,234,253,  9,210,  0,197, 11,255,207,170, 16,203,  6,200,  1,126,  2,181,  0,197,252,191,  1,192,255,191,252,195,  2,195,  0, 64,  4, 67,  1,  6,  6,250,253,  0,  0,254,199,215, 10,251, 21,  1,254,  1, 27,  6,232,191,212, 54, 49, 27,250,238,222,216,  3,231, 67,231,  1,  1,  2,250,  0,  0,  0,254,  0,  5,  0,253,246,  2,  0,  4,253,240,234, 14,245,245,249,  5,255,  1,  8,253,224,255, 10,  0,
      247, 12,205,  4, 16, 21,251,250,  0,218, 17,  9,254,214,192,203,167,  0,191,226, 24,244,184, 31,223,  4,  4,250,  7,207,196,196, 16,204,  1,245,230,237,255, 15,255,215,185, 64,  8,199, 66, 13,  0,251, 65,212, 35,200,253, 14,251, 13, 50,203,  3,234,233,255,249,  3,  1,  1, 20,  3,255,  2,225,  6,  4,  1,231,250,  0,  0, 45,252,  0,  1,244, 34,239,  0, 28,177,  2,255, 55,239,  1, 35,  6,251,254,252,240, 50,222,240, 14,253,253,  1,  7,  5, 68,250,  3,252,  0,252,  1,255,  0,  3,  3,  7,  3,  1,253,  2,  4,249,  7,224, 23,252,242,244,201,  4, 37,233,230,238,246,241,253, 13,  2,253,185, 71,237,  9,191,254, 45,255, 43,200,239,173,238,216,  0,249,231,245,255,197,242,255,  0,244,251,  0,251,250,  2,  2, 12,  1,245,252,255,226, 66,245,  1,  4,255,  1,254,205,  0,210,
        2,220,  1,190, 12, 42,  9, 64,217,250,  1,192,  4,209,255,253, 43,236,184,  5,252, 56,  7,193,251, 19,226,  0,246,255,255,253,  2,237,  4,209, 76,249,250,  7,241,234,  0,192, 58,238,  3,239,251,  3, 50,187,248,250, 62,247, 47, 21,  0, 55,243, 38, 74,238,198,255,  1,255, 42,248,  5,254, 21,  0,255,242,247, 17,  1,  1,253,210,  8, 15,228,241,191,239,  8,  5,  0,234,  3,247, 42, 25,245,  1,  5,242, 10,  0,  5,245,  8,  2,  5,254, 10,  0,189,  4,  3,  1,  1,246,  1,  3,  3,  1,252,252,254,253,254,253,  0,  4,246,243,220,203,  9,255,236,244,  2,183,221,254,255,192,197,251,248,194, 12,196,  6,195,250,181,240,234,122,185, 28, 22, 74,169,253,186, 14,208,  0,168, 12,209,  2,235,252,  2,  0,202,  2,254,246,  0,255,246,246,255,  4,210,255,197,255,220,  0,194,253,  6,
        1,  2,  0,  0, 27,241,  0,254, 20,252,  0,  8,255,  0,  1,  2, 48,  3, 92,  6, 63,229, 62, 27,253,199,223,  0,254,254,255,  3,  0,244,255, 14,208,230,255,217,236,233,  8,  0,210,252,  3,245,248,224,184,188,201, 40,  2,201,244,  1,255,  5,201,233,188,188,128, 62,253,208,243,194,255,196, 19,192,253,198, 28,214,255,189,252,197,  3,194,235,222,  5,210, 10,244,  0,203,225,191,  8,196,249,251,242,195,  9,254,191,  0,254,  0,255, 64, 15,255, 64,252,253,255,252,186,  0,  7,  3,  0,255,  8,  0,194,252,240,  0,252,  4, 50,247,190,  2,254,202,253, 65,  0,244,126,  9,212, 38,206, 26,226, 75,  9, 52,240,120,228,238,192,244,226,252, 19,  5,  4, 10,182, 12,124,  0,254,  2,246,  0,249,  2,240, 11,190, 26,188,  9,207,252, 16,218,194,  1,189,255,191,  0, 29,  0, 44,  1,243,
    },
  },
  {  //  31
    {
       10,  3, 45, 12,217,217,251, 30,220,236,237,226,244,244,  3, 18,  3,  2,249,  8,254, 23, 15,221,251, 22, 10,255,224, 16,242,255, 46, 14,241,  1,245,252,249,248,225,246,  0,252,231,222,248,255,  4,238,244,  1,243,233,255,249,231,221, 53,  2,110, 19, 19,  1,  6,  4,  0,255,254,  5,  0,  0,  9,255,  0,  0,253,  4,  1,255,194, 11,255,  0, 26, 36,255,255,  3, 35,  4,  0,210,  6,125,  6,253, 33,  3, 31,181,  6, 37,210, 30, 72,240, 43,  6, 28, 49, 52,242,253,252,234,249,229,187, 99,236,224,250,197,227,  1,195, 55,255,  8,196,247,249,254, 16,241,185,  6,230,167,  7,  4,  5, 42,205, 46,  5, 57, 57, 58, 67,250,213, 14,248,  0,  1, 47,254, 76,  6,  1,251,  4,245,251,  1,237,246,  2,252,243,254,  1,254,247,  3,240,  1, 19,  0, 29,242, 28,  4,255,254,250,253,241,254,252,
      248, 25,194,  4, 31,228,  5, 17, 43, 56, 15, 49, 30,229,251,  1,251,253, 25,218,183,181, 22,255,250,241, 59,252,172, 15,255,  4,  5,  7,  2,  6,254,238,251,254, 19,253,255, 15,215, 16,  3,251,255,251,252,249,236,251,242,  5, 16,245,252,255,  2,245,251,253,253,  0,  0,  0,  0,249,  1,  0,  4,  2,  0,  0, 13,  3,  0,255,244,240,  1,  8, 28,249,251,255, 22,221,  4,  0,201,  0, 52,252, 57, 39,232, 62, 41, 48,202, 88,242,244,251,  5,195,244, 27,245, 44,217, 38, 20, 39,216,238, 10,  1,  2, 12, 65,219, 64,200, 75,  6,250, 35,190, 18,  0,201,255,  8, 21, 28, 16,  1,207,179,238, 16,236, 24, 15,  6,209, 45,230,  4, 11,  4, 16, 18,255, 24, 96,244,  8,  2,  0,  3,  2,  2,254,250,  0,242,  0,251, 17,136, 13,254, 26, 25,  0, 12,246,252,236,254,  1,253, 14,242,  5,251,  3,
       39, 53,  2,174, 21,250, 10, 36, 23, 15,197, 11,217, 34,219,246,222,241, 35, 12, 14, 11,153,212,157,  9, 40,240,  8,212,248,253,202,  2,254,  1,253,235,  1, 52,224,251,  2,  1, 38,232,252,  0, 10,245,247, 12, 63,  6,  7,219,237,242,249,254,255,213, 33,253,253, 11,  0,  0,  7,247,  0,  1,252, 21,255,  1,248,  5,254,255,251, 40,  5,255,252,250,254,  3,250, 10,  2,  3,181,253,  2,  1,237,208,212, 31,116, 32,250,  0,  8,  8,251, 73,222,255,237, 70,231, 35,247,215, 62,254,123,247,  2, 16, 20,250, 18,195,147,196,251,  3,220, 54, 27,  0,216,185, 31,  3, 16,  2,  5,254,252,  9,  3,  2,241, 40,239,230,144, 13,  2, 44,254, 37,  5, 17,179, 24,  2,254,  2, 29,249,246,  4, 18, 31,246,255,  4,  1,253,252,254,241,  2,247, 11,251,236,  7,239,  1,254,226,248,  8,200,  5, 20,
      255,  0,252,249,  1,252,254, 65, 38, 14, 36,167, 19, 28,104,220,  3,  3,201,224,243,218, 13, 28, 44, 12, 41,202, 74, 52,194, 21,  6,255,  3,  2, 46, 20,115,  4, 69, 13,252,254,  0,251,252,251,183, 13,255,  8,251, 10, 30, 23,  5,253,250,248, 17,  0,  3,253,254,243,  1,254,255,  6,  3,  0,253,247,  1,255,  0,  1,  2,  0,252,254,  2,250,197, 21,  5,254, 84, 84,  0,  8,  1,188,129,255,250, 63, 33,188,126,177,  3, 80,252, 57,251, 43, 39, 15, 17,  2, 15,216,243,247,  4, 72,  4,169,  8,251,246,196, 12, 30,220,204, 11,230,255, 44, 65,  4,190,158,  6,249,  3,247, 53,  3,  0,196, 62,248,  9,209, 40,197, 52, 87,  5, 76,246,  7,  1,178,126,205,  0, 50,236, 25,228,  1, 45,218,  0,  0, 13,251,254,  6, 74,  3,  1,253, 13, 33,249,247, 72, 61,  2,253,242,  6,194, 13,209,  1,
    },
    {
        3,226, 18, 93, 19, 53, 22,208,177, 22,244,  0, 19, 42, 35,  5,226,214,  4,239,246, 18,251,250, 24,233, 56,229,226,238, 32,  0,  3, 23,  2,249,  4,184,249,245, 23, 11,  1,251,  2,  5,248,246, 45,233, 62,250,237,230,  4,  9,216,234, 15,  1, 92,243,  0, 47,  2,  8,  1,  4,  2,254,255,255,  2,  3,  1,  0,255,253,  1,250,  3,246, 43,253, 53,  5,  0,254,  9, 20,  0,  0,233,246,199,235,218,252,  5, 24, 21,145, 23,232, 22, 33, 33, 16, 18,222,171, 19,243, 16,  6, 23, 30,254, 30, 15,244,247,247,190, 85,  2, 14,255,  5, 19,235, 13,217,235,  0,226,252,251,212,252,204, 15,239,242,  1,219,  1,  7, 23,  2,254,239, 51, 10, 82,224,  2, 32, 14,  8, 14,249,238,  5,  3,243,  8,243,  0,255,251,  0, 16,249,  6, 11,  3,  3,  8,  1,255,  4, 22,243,253,  2,  0,253,254,252,252,  0,
       10, 18,253, 69,  3, 14, 19, 37,248,216,241, 17,226,  7,215,  9, 35, 45,251,218, 46,238, 53,212, 14,244,248, 14, 19, 14, 65,210,  7,248,253,249,254, 54, 25,  6, 24, 11,247,  7, 26,  3,243,  4, 15,247,  1,  3, 15, 10,251,243, 40, 53, 41,253,255, 49, 14,238,249,  7,255,  0,  1,255,  3,  1,  1,251,  0,  2,  1,  4,  4,  5,254,  5,253,254,  7,197,245,  0, 12,223, 24,  1, 11,240,250, 19, 60,209, 18, 41, 22,  0, 58,211,254, 21, 35,196,206,242, 41,154,  3,  2,  5,252, 11,199,210,234,243, 28, 14,249,253,236, 59,209,  4, 79,  5,252,251,  4,224,  0,228, 22,246,  4, 57,  7, 21,232,  3,  2,254, 25,229, 57, 15,249,236,  5,171,220,215,  3, 23,232,254, 28, 46,254, 11,  4,254, 11,  1,  2,205,  0,249,246, 15,244,251,  0,253, 52,  1,248,  0,253,  0,  1,  1,  3,  2,253,  7, 13,
        2,  4, 21,244, 29,177,244,231, 20, 17,214, 22,219,168, 12, 10,233, 11,201, 39, 37, 65,177, 21,203,205, 52,  7,215,252, 24,255,  0, 28,244,251, 55,241,  1,249,198,219,238,  7,251,221,  6,  1, 24,232,172, 10, 43,251,  4,237, 54,180, 36,248,  1,155, 24, 18,192,  0,  0,  0,  4,252,  2,  3,235,252,  2,  1,255,251,  0,  0,  2,  1,243,  3,183,249,251,  2,255,217,213,241,250,240,221, 28,228,231,250,239, 67, 44,  5,243, 19,224,252,242,  5,211, 52,129, 12,  7,  5,  9, 18,185,  7,247, 42, 21, 16,224, 29,235,213, 39,243, 13,238,225,245, 16, 51, 20, 12,  0, 24,238,255,144,221,  9, 43,  2,  5,  4, 52,242,  5, 72, 60,  0, 39,  4, 57, 64, 13,  0,228,236,239,  9,  2,234,253, 42,232,  5,255,227,  0,  1,242, 14,242,253,  0,253,240,255,247,251,  1,246,253,251,226,241, 56,253,
      243,212,  0, 47, 12,255,255, 40,249,248,249,157,125,  2, 13, 57, 13,  9,250,  2,  3, 32,193, 15,228, 10,201, 11,  7,  1, 29, 21,229,244,  2, 25,192, 33,  0,  8,248, 37,243,235, 34, 32,  1,243,254,251, 39,247, 61, 39,  6,211,246, 68,  3,248, 16,100,  8, 13,248,  3,254,  2,249,254,254,  1,  6,  2,  1,  2,241,  6,255,  0,220, 10,  8,255, 58, 26,254,250, 46,253,  5,  0, 37,250,195,  2,242,245,202, 21,188, 20,  1,242, 63,  9, 28,  4,247, 34, 74,185,251, 29, 23, 29,  9, 54,250,231, 38,242,  4, 16, 39, 10,  1, 85,168,221, 17, 26,255, 27,120, 80,254,246,198, 59,201,250,166,168,234,  1,249,231,  7,231,  8,226,214, 12,254, 45,226,245,230, 48, 25,  1,  1, 19,255,237,214, 48,255,253,  9,243,  3,248,253, 65,  4,  3,255,227, 20,  9,193,204,  1, 13,254,252,  5, 13, 13,  4,
    },
    {
       26, 21,210,237, 41, 21,191,165,248,  5,197, 15,228,254,  5, 54,222,  2,215,201,214,229,215, 66,219,249,249,187, 75,216,228, 76,243,212,240,248,253,202,246,254,220,195, 18,  0,207,209,254, 15,225,197,251,250,221,189,176, 21,232,204,  5,239, 58,164,247,  2,  0, 10,  0,  1,255,251,  1,  0,  1,231,255,  0,236,  1,  1,  0,253,222,  0,255,253,252,  0,  0,  2,235,  4,  0, 13,  1,  0,255,243,240,237,219,253,251, 29,239,  1,194,223,218,191,251,183,195,238,236,194,212,197,246,245,  7, 20,238,108, 37,175,230, 63, 32,  4,248,219,238,  1,255,  0,240,255, 20, 13, 15,252, 10,232,172,193,192,176,195,207,214,191,197,198,198,195,192,230,209,182, 14,  0,254,254,254,  1,250,  1,  3,254,255,253,246,  1,  1,  2,254,247,242,242,245,  2,235,253, 29,241,232,  4,233,245,240, 44,251,
      253,243,251,154,248,253, 16,127,218,255,  6,215,222,231,  3,196, 28,250,255, 63,254,  3,176, 56,231,244,216, 33,234, 32,  0, 57,200,250,220,253,202,  3,243,255,  0,235,255,  4,247,248,  5,  4,253,255, 54,  0,201,  0,205,250,248,  0,252,251,  4,245,  5,  5,252,222,  0,  0,253,232,  0,  0,243,189,  0,255,246,248,255,255,248,191,255,255,255,226,  3,255,253,  3,255,  0,  1,250,  0,253,223,235,  2,252,222,213,223,113,216,255,  6,255,247,  0,222,  0,232,241, 10,233, 85,236,248, 13,226,255, 13,  5, 34,229, 31,  9,  1,  3,195,225,  0, 11, 23, 23,250,189,  3,249,246,251, 32,  1,  3, 57,190,255, 43, 42, 16, 59,255,195,  0,252,248,231,  2,  2,196,204,210,202,198,192,194,172,196,211,194,201,196,192,197,200,228,210,193,203,191,198,182,173,202,193,192,239,201,204,190,235,
      254,250, 52,246, 13,254,  0, 66,250,  6, 12,255,  0,252, 20,218,236,216,243, 62,  9,240,  2,255,222, 23, 19,254,195, 22,214, 33, 10,  2, 66,  5,250, 54,247,239,252,255,249,254,  0,255,239,241, 15,254,  3,255,252,255,  6,191,  0,  0,255,255,220,239,241,231,193,208,189,192,197,189,156,185,195,233,190,192,226,205,192,192,197,193,192,192,196,192,192,209,232,232,192,192,194,196,195,195, 15,  7,  0, 54, 30,  0,  0,  1,  8,  7, 14,164,225,238,201,241,255,255, 81,  2, 19,251, 13,246,228,251, 39,  4,240, 15,253,251,149,199, 48,255,212,212, 31,237,  1, 72,251, 46,  0, 98,196,253,  0,234,236,  1,219,243,250,236, 54,254,246,  1,230,  3,193,  4,  0,  3,  3,251,252,  2,242,  4,  0,  3,126,  0,255,  2,  0,255,  4,  0, 36,224,252,222,255,248,  1,  3,253,244,216,253,  0,229,
      246,240,219, 12, 16, 17,190, 30,  3,  1,  4, 19,198,167,191,255, 24, 47,224, 19,195,255,186,  0,238,240,196,227,  6, 20,239, 10, 54,239,  1,250, 32,188,201, 66, 60,192,233,  4,167,180,249, 21,200,199,252,252,151,252,231,227,253,189,255,236,198,233,255,230,  2, 14,255,  0,254,  0,255,255,  2,  0,  2,  0,228,  1,  7,255,  0,  0,  0,253,  1,  0,  1,246, 11,  0,  0,251,  0,  0,255,255,251,  3,134, 25,248, 21,193,191,225,199, 40, 19,241,242,254,227,255, 71,255, 16,199,  3,206,232,187, 13, 44,226,183,238,196,194, 11,212,145, 19,127,249,  6,  3,  3, 34, 19,  1,194,251, 64, 12, 78,188,183,206,255,226,255,222,  5,214,248,220,220, 45,  0,236,231,253,255,231,  3,255,  1,  1,248, 12,  1,  0,255,252,  0,248, 24,223,  0,222,248, 15,  1,250,241,254,246,254,  3,222,  2,254,
    },
  },
  {  //  32
    {
       17, 68, 19, 12, 16, 58,  8,209, 17, 71, 12, 65, 50, 47, 29, 19,217, 58, 11, 22,205,237, 18,180,161,172, 42, 32, 70, 17,254, 36,  1, 34, 52, 56,251, 21,191, 17,255, 24, 62, 55,200,245, 48,123, 21,  3, 59,205, 40,  1, 87, 36,255,  0, 69, 72, 25,249, 56,125,  0,  0,228,  5,  0,  1, 15,232,  0,255,248,  3,  2,252,248,  0,255,  0,241,  0,  3,254,221,  4,  2,  3,248,231,248,  0, 34,  1, 10, 60, 16, 64, 66,198, 55, 67,248, 52, 45,102,126,  4, 21,228,253, 23,  6,240,171, 32,196, 74,  4,153, 14, 15,211, 51,211,255,  9,254,253,255, 35, 17,253,  1,236,230, 30,  5,255,  1,  1,251,196,  2, 52,216,110,246, 59, 11,254,  7,  0,247,207, 67, 60,188,  0, 21,255,  5,  0,  1,  3,  0,255,254,  2,255,  1,  1,255,  1,  0,251,  8,  3,249,231,  8,230,  1,  1,  0,255, 21,  4,235,214,
        5,  0,247,235, 23, 54,197, 13,194,234,223,254,210, 81,213,244,254,181,212,  7,240,254,216, 54,  9, 63, 60,246, 50, 23,198,  3,254,  0,  2,251,247,  3,  2, 89, 34, 35,251, 20,248, 16,106,218, 41,253,  8, 14,242, 12,159, 19,228, 14,244,207,250,235,255,201,  0,  0, 10, 28,  0,  0,246,  1,  0,  0, 48,  0,  0,254,242,  0,  2,  0, 12,247,  0,  0, 78,251,255,  2, 85,  7,248,251,254,251,252, 16, 34,227,242, 76,  5,157,254, 17,232, 64,254,  2,178,120,223, 42,  4, 17, 20, 29,231,234,251,254,  0,  1,202,214,221,171,254, 31,253,  3,  9,252,254,253, 60,253, 20, 38,  2,  5,243,  0,221,183,232, 12,227,241,  7,225,  2,203, 42,  6,243,118,250,229,  0,  8,  0,235,250,  0,251,  4,  2,  0,  0,250,  3,  1,  1,  1,242, 11,248,249,  9, 73,230, 10,254,  4,  1,  0,213,241, 10,231,
        2,247, 61, 16, 28,192,253,244, 44, 68, 86,  5,236,183,201,  2,241,249, 71, 34,  4, 62,124,253,  9,150,244,  4,  2,248, 65, 37,  1,  0,158,121,222,246,253,222,249,250,  0, 34, 49,204,200,121,249, 27, 20,253,206,255,  3,236,253,254,  1,190,  8,  0, 21, 45,  1,  0,  1,253,  0,  1,  2,  1,  0,  1,232,  3,  7,254,254,  4,255,255, 42,  5,  1,  2,  2,205,  1,254, 63,  1,252,253, 87, 14,  2,254,199,246, 66,252, 68,218, 32,255, 14, 49,  3, 64,240,199,238,238,212,238, 13, 13,254,  2,  2,228,  3,  2, 24,223,122,101,  1,  4,237, 60,253, 11, 14,255, 22,  0, 56, 76,255, 58,  6,252, 36, 53,251,210, 59,254, 32,221,248, 50,240,185, 18,254,  2,247,  0,215,  5, 75,255,  3,253,  0,  6,230,249,  0,  1, 10,253,255,247,  2,243,  5,243,  4,  3, 77,  2,255,  0,  0,  3,  4,  8,  7,
        4, 15,247,252,  0,  0, 22, 31,176,249, 30,  2,252, 10,  7,250,250,196, 57,254,255,255, 41,252,255,255, 64, 63, 11,  4, 25, 21,255,253, 67,204, 97,  2, 36, 16, 37,252, 60,207,244,218,  9,  6,241,254,192,168, 22,  0,203, 10,232, 26,  5, 56, 14, 24, 13,248,  1,  0,  0,  0,  0,  1,  0,  0,241,  0,  1,  2,252,  6,  0,  0,215,  2,  3,  1,254,  0,254, 47,208,  2,254,  2,  2,  3, 37,124,  2,220,187,250,234,202,200,  7, 50,250,216, 42,219,  9, 56,  1,251,221, 20,  0,185, 53,121,170,255,211,  4,239,157,222,111,  9,  2, 19,  7,228,200, 13,  0,236,245,  3,238,215, 36,249,201,  3, 48,  0, 25, 22, 44,  0, 45,181,245,  6,  7,204, 27,213,  1, 49,  0, 19, 21,207, 34,  5,252,222,  8, 30,255,  2,228,  0,202, 79,  2,  1,  2,213,  2,253,250,152,  1,  0,255,  1, 14,253,244, 48,
    },
    {
        0,251, 37, 50,251, 12, 17,240,236,241,241, 13,201,  4,251,254, 15,225,197,248, 73,245,  0, 11,193, 54,197, 67,210,240,251,246,  1,  2,104,234,238,  9,  1,  4,233,204,248,250,  8, 67, 64, 10, 32,232, 65,249,186,202,241, 17,199,239, 20,210,214,  0,244, 15,  0,  2,  0,  4,  2,  5,254,  0, 10,  3,255,255,  2,  4, 16,  2,  1,219,  1,255, 28,245,  2,  3,255,  4,  6, 32, 35,  3,  0,  8, 12, 14,207,225, 60, 20, 35, 21,233, 23,239,166, 20,  3,  8, 91,215, 42,  0,246, 29,  0,229,208,232,  3,  8,239,223,238, 30,192,243, 51,226, 12,241,251,184,127,  6,  5, 10,226,247, 15,253, 12,  5,  8,250,  6, 74,253,  7, 38,  3,  0,239,  1, 45, 12,244, 69,  2,249,  8, 26,255,254,  0, 14,  1,  2,255,  3,253,254,249,254,255,  0,  1,  5,  2,255,  0,247,  2,  1,  1,  0,242,194,  0,233,
      252, 15, 44,122,236,237,186,210,240,  7,243, 23, 33, 14,  5,253,221, 29, 65,241,  0,  4,  7,220,197,  8,  7,  4,234,  5,238,132,  0,244,254,  6, 28,236,254,162, 34,250,254,  5,  2, 32,194,  1,  1, 28, 65,252, 14,236, 50, 18, 46, 70,  2, 52,251, 62,227,  7,253, 11,255,255,  1,254,  1,  0,255,  4,  0,  0, 17,  1,  2,254,  7,255, 38,255,  1, 31,238,  1,250,  3,  0,  1,254,  4,239,  6, 17,223,237,183,122, 65,229,223, 68, 43, 24, 86,233,195,  2,190, 26,253,  2, 15, 17, 77, 21,242, 10,  4,231,  5,196,205,196, 13,  0,194, 28,224, 23, 55,211, 35, 44, 21, 21, 51, 97,166, 29, 10, 63,253, 66, 49, 55,246,248,216,245,  0, 24,255,229,221, 11,219,  0,254,  2,  2,252, 11,255,  7,254,  0,  2,  1, 14,254,  3,247,225, 11,  1,  1,229, 10,252,240,  1,  0,254,  1,  2,  0,  6,243,
        9, 19,178, 30,232, 28,196, 15,227,220,214,210,  4,235,178,117,252,195,241,251,  0,218,253,215,232,194,206,  0,234,199, 86,224,241,233,  0,238,241, 34,  2,  6, 52, 56, 19,  0, 31, 10,  7,  2,249,247,228, 21, 86, 70,222, 17,236,218, 57,251,133, 22, 14,142,  0,  2,  1,  3,253,  1,  1,  0,254,254,  1,  0,247,  3,  0,253, 43,252,253,200,  6,  1,242, 13,  8,241,  1,  0, 13,217,245,  2, 34,241, 32,254,253, 20,191,144,246, 25, 54, 17, 53,202, 44, 54,208,249, 29, 51,206,237, 55, 56,  8,253,  4,251, 12,247, 49, 45,  7, 53, 13,229,242,230,135,122,249,247,220,223,203, 50, 21,222,165,243,212,226, 28,220,  7,206,  1,  1,251,246,  3, 50,252,  2,254,  9,227, 29, 52,  9,242,254,  0,  3,254,255,  3,  4,254,249,250,249,251, 12,  9,  0,217,  4,255,  0,  3,  0,  1,  0,  1, 27,
      233,  2,244,109,  5,253, 25, 74,236, 13, 14,212,254,222,182,179, 45,253, 38,120, 52,  5,217,  3,235,  0,242,201,219, 38,188, 84,253, 46,  2,242,254,245,252,  6, 17, 10,  3,  5, 36,245,249,  7, 41, 11, 17,255,244,255, 65,143, 24, 74, 29,249,138,253,129, 58,  0,251,255,  0,  0,254,  2,254,  2,  3,  0,255,  8,  5,  0,  1,230, 31,241,  0,232,250,246,204, 11,233,  2,253,  5,236,252,  2, 15,239,238, 13,243,218, 20, 37, 17, 14,218, 16, 13,  5, 50, 98, 10,  6, 11,  3,  5, 26, 58,254,  5,245,254,253,  8,243,  6, 59,  6, 23,218, 34,  4,208,244, 77,254, 41,236,250,215, 61, 47,222, 10, 40, 32,245,229, 22,241,244,  1,  0, 25,232,248, 29, 52,164,  1, 19, 52,247,251,248,199,  1,  0,255,  4, 12, 12,241,132,193,  3,  3,  8,223, 19, 14,  0,251,255,  3,  2,  2,  3,  3, 64, 71,
    },
    {
        4, 21, 11,205, 43,251, 35,  3, 32, 40, 25, 31, 16,227,248,169,220,241, 11,200,  2,218,255,224,238, 31, 36, 17,135,242,188,249,  7, 13,  4, 44, 15, 28, 12,254, 22,240,  9,251,217,245,  7,246,206,  4,251,249, 48,  3, 63,227,206, 27,  0, 23, 45,  9, 59,239,241,  1,254,  4,  4,  6,  0,  0,  4,  5,255,254, 26,  2,255,252, 38,  6,  3,255,  0,  4, 50,  2, 15,251,  5,  6,238, 50, 58,242, 42,252,  0,226,  4,120, 36,204,  4,245, 12,247, 63, 18,239, 30, 19,  3, 17,228,255,254,213, 58,204, 77,  5,247,226,241, 27,234,  2, 10, 38, 14,249,233, 27,247,251, 11,  4,217,232, 30, 13, 13,  0, 19, 20,250, 37, 39, 30, 16,  2,255,255,249,249,  3,248, 30,  3,  2,254,254,255,254, 20, 42,255,  8,  5,240, 16, 17, 28,253,251,  8,  4,  6,240, 35, 16,249,238,  0,  4, 25,  9, 19,254,216,
        7, 49, 10, 38,239,  5, 41, 71,  1, 31,248, 38, 61,  7,218,247,247,  9,209, 12, 38, 32, 12, 24, 21, 11, 15, 27, 31,243,236,228,238, 51,  4,  7, 24,244, 21, 19,249,  7,254,241, 45, 54, 12,206, 55,230, 29, 27, 24,244,246, 21,239,155, 16, 35,189,210,238,  1,253,  7,  1,  0,  4,  7, 10,  0,242,253,  3,  2, 25,  7,250,  2,251,  5,248,  8,  0,247,  1,252, 19,246,252,236, 11,  0,252,250, 30,214,227,255, 64,220,  3,241, 20,  2,253,234,188, 36,221, 43, 21, 33, 61, 38,200, 18, 33,156,237, 14, 18, 30, 28, 34,111,252,237, 32, 48, 23,249,251,246,233,  4,  0,  7,228,244,235,255, 49, 10, 23, 24, 72,250,  9, 18,  1,237,253, 45,251, 15,208,228,239,  7,243,253,  4, 12,240,249,227,  3,  6,248,225,  0,248, 13,141,  4,  7,  0, 14, 24, 21,225,247,251, 25,  1,248, 13,252,  9, 43,
        2,210, 13,146,255,  2,  7,120,221,231,242,186,  1,214,217, 17,235,  8,234,  3,244, 77, 21, 15, 74,  3,179,227,208,211, 85, 22,235,254,251,243,240, 18, 20,254,  3,227,  2,232,237,232,  6,255,191,254,223, 15, 24,224, 33,222, 31,  6,221,  3,121,241, 42,  1,246,250,  0,253,217,250, 15,249,215,  1,255,  0,244,235,  5,254,253,  1,253,254,223, 19,  1,  6,213,242,251,  5, 12, 49,118,  3,238,197,241,248,192,234,186,  9,  3,237, 19,252,254, 11,156,239,  6,214, 19, 49,  5, 16,251,250,147,  0, 61,238,254, 34, 56,199,248,253,169,215,  9,  2,255,  4, 19, 33,230, 32,251,255,253,238,  5, 24,248, 11,245, 11,235, 12,254,246,  1,112,  3, 17,  6,  4,  5,252, 13,  5,  0, 19,  1,234,255, 21,  9,236, 16,226,225,121,254, 69,  1, 15,  6,207, 52, 55,  1,251,  4, 11, 17,  0,234, 26,
       25,203,252, 13,253, 38,  8, 53,  4,  5,  8,235,  1,210, 53,207,233,228,  8, 21, 14, 40,253,112,255,255,240,246, 49,210,235,220,244,226,246, 10, 18, 13,248,245,  0,  3,  3,239, 49,233,242,237, 55,231,247, 22,218, 24, 35, 15,224,230, 36,  1,  0, 78,185, 58,236,253,  3,  0, 10,  5,249,  1, 57,233,  0,  1,188,  2,  5,253, 71, 21,  8,254,  6,249,244,  3,223,  2,133,  3, 60, 20,238, 53,215, 30,252, 33,  3, 10, 88, 39,255, 62,250, 31,219,250,113, 10,  0, 45, 53,239, 29,175, 56, 13, 64,  7,212,  3, 31,231,106, 10,  5, 12,  5,251,247,241, 73,230,255,  0,223,233,245,  7,211,126, 34,229,224,198,  7,213, 18,208,205,215,201,200,  5,249,222,253,254,243,251,232,  1,230, 12, 24,254,246,252,252,239,247,237,223,  0,226,  2, 36, 32, 62,  5,204,248,192,  2,  2,  1, 60,249, 40,
    },
  },
  {  //  33
    {
       17,  7,200, 12, 39, 68, 60, 33,  2,244,251, 11, 40, 61, 64,246,243,253,246,240,226, 92, 44,  7, 51,194,235, 43, 48, 57,216, 11,212,252,244,253, 56,254,255,250,208,  4,255,244,  4,222,  2,255,221,224, 11,255, 11, 11,  8,  7, 45, 11, 15,237,237,244,232,145,  3, 16,  0,  0,  7,  0,  0,  0, 10,  1,255,255,254, 13,  0,255, 28,  8,254,254, 86, 10,255,  2, 30,  0,  4,  0,234,162,254, 16, 10, 50,238, 44,255,252, 17,  4,244, 52,  3, 31, 82,236, 12,230, 14, 29,245,244, 21,238,  9,160, 10,  9, 11,234, 13,195, 18,  7, 12,  0, 13, 32,253,  1,243,250,230, 13,232,197,  2,215,  2,210,203,216, 16, 10, 46,252,190, 14,  4,254, 20, 22,  7, 18, 32,179,255,253,  0,  1,  0,  2,  5,  8,  0,  2,  0,  1,  9,239,  1,253,253,250,255,254,254,  2,236,248,255,255,  0,254,  1,255,216,  1,
      247,  1,249,  2,  1,235,  9,231,  1,253, 26,219,  5, 10, 24, 46,221,236,219, 49, 58,222,241,239,248,250,234, 38, 30,224,254,206, 43,251, 23,  3,213,244,  2,  5, 19,241,206,251,203, 48,  4,  6,  8,  4,245, 13,  4,  8, 19,253,253,  0,247, 15,253,  4,  0, 94,  6,241,254,  0,248,  3,  0,  1,  0,255,  0,  0,252,193,255,  1, 14, 11,  6,255, 24, 70,248,253,  9,  2,  4,  0, 34,  4,  5,  3,251,244, 25,232,  8,191,  5,245, 25,251,  9, 46,252,253,  6, 30, 16,226,  4,  7,239,193,248,209,255,251,254, 18,247,199,  8, 38,252,  7, 33,223,254,244,  8,247, 30,251,239, 33, 24, 46,161,251,244, 87, 44, 46,223,233, 29,  2,255,  5,228,243,251,175,194,158,  1,254,  1,255,  1, 46,249,251,253,  3,252,  5,251,  0,250,255,  1, 14,  5,  9, 29,255,238, 26,255,  0,  0,  2,255,  2,  2,250,
        1,223,  1,  4,246, 19,253, 22,244, 25,  6,  4,244, 87,105,  6, 23, 46,254, 21, 27,128, 69,184,245, 18, 69,  2,222,  3, 76,244,  0, 19,252,  0,239,  9,  0,255, 49, 69, 29,  9,231, 61,253,253, 62,  5,252,  0,  7,242,225,  1,217,235, 20, 10,158,242, 20,235,249,254,  1,255,  3, 22,  1,  0,  2,  1,  2,  2,207, 36,  2,252,247, 17,250,254,  0,245,246,  1,  3,  9,249,  1,245,242, 58,  9,  6,253,177, 59, 61,  1, 19,  5,  1,  1,244,  8,182,  1,  3,233,250, 33,180,213, 14,229,  7,232,  3,250,254,  2,220,244,  7,155,248,241, 13,195, 53,237,  1,248,248, 12,242,249,  8, 42,125, 32, 58,  3,215,214,228,241, 24,227,  2,  8,200, 24,255, 61, 10, 23,  0,255,213, 66,204, 20,  4,255,  1,255,  1,255,  0,  0,  2,  4,  1,  1,  2,251,255, 25, 10, 10,  0,255,253,  0, 12,  4,255,247,
      250,  4, 20,240,190,236,193,226,226,234,220, 17,  5,215,192, 19,228,  0,252,226, 38,252,215,  8,  5,245,148, 24,108,253, 72,245,  1,  4,  1,  9,223,242,  1,  6,236, 21, 70,  2, 11, 17,  0,  9,238,239,251,  0, 28,  3, 50,233, 30, 47, 35,  3,248, 14,  3, 54,  0,  1,  1,255,  0,  0,  0,255,  0,  1,  0,  1, 85,252,  1,252, 33,251,251,255,241,252,224,  0,  2,250, 11,255,  4,  2,254, 10,251,238,244,228,251,  3, 17, 10, 10,  4,  6,239,242, 42,  9,219,  0,237,255,  1,  4,255, 48, 13,  5,  3,  4,  8, 18,244,  4,251,  0, 30, 42,241,197, 12, 62, 17, 11,250,  4, 34, 15, 10,190,250,  4,238, 18,  0,250, 36, 45, 68,  3,  3, 83,255, 59, 33, 83, 60,  0,  5,  4,212,  8,209,  1,255,  1,252,  0,  0, 48,  2,253, 56,254,  2,  7,254,252,246,  6, 47,255,255,255,  0,255,255, 76,255,
    },
    {
        6,  5,  2, 63,189,195,194, 65,  1,  4,  4,252, 56,191,199,198,249,  0,  8,  1, 58,190,206,157,238,251,205,253, 89, 66, 12, 12,  1,231,255,253,  5,239,255,  1,  3, 30, 11,248,195, 10,255,244,222,219,  8,248, 31,225,239,248, 34,167, 15,  1,186,177, 45, 19,253,255,  4,  0,225,  1,  7,  1,  6,  1,  0,  1, 24,  3,  0,  8, 48, 29,  2, 11, 16,226,254,  4,  2,255,  7, 20,  2,254,254, 58, 30, 15,190,  9, 53,  6,215,169,  0,242,212,220,  7,220, 57, 71, 11, 44,235,  6, 12,  6,  3,  6,227,249,248,247,220,  7,245,172,219, 24, 45, 68,  9,197,182,  4, 12,236,201,195,203, 58, 65, 56,  2,  7,  4,247,255,  3,  9,216,  4,  9,  0,255,  7,254,237,223,223,253,  9,240,  2, 26,  9,250,  2,  1,  3,  1,226,248,251,255,  7,250,253,254,  0,  0,237, 10, 12,  2,  0,  0,  2,  0,  1,255,
        1,  4,  1,226,255,254,  6, 11,255,  2,  1,  1,  6,249,118,244,  2,254,  7,  1,  0,253,255, 38,  0,253, 59,252,128,188,202, 44,  6, 56,  1,253,195, 45,246,252,246,239,  1,183,252,245,  0,  1, 64,243, 12,222, 43,238,  1,  9,249, 94, 21,239,255,134,254, 14, 18,201,255,  0,  7,248,255,255,253,252,  0,254,183, 63,  5,  0,216,186,  3,240, 10, 39,254,  3,253,236,232,165,225,213,243,196,230,253,  3,215, 54, 36,226,185,243,253,241, 60,251, 26, 42, 81, 35,186,  1,231,  0, 56,242,  0,234,246, 58,  1,227, 24,188,188, 58,255,225,223, 34, 11,202, 28,255,241,251, 34, 22, 32,194,121,201,210, 63, 18, 73,254,  8,  6,247, 19,246, 10, 36, 34, 64,  8, 10,249,  0, 13,168,  3,251,  3,254,255,  1,  0,  7,254,  7,  0,  2,247, 55,255,  3,  3,  2,  3,250,247,  0,  0,  8, 38,238,  2,
        1,255,  0,123,254,  0,255, 44,254,  1,  0,249,125,  5,  2, 32,  3,255,  0, 13,  9,110,254,250,  2,  3,  2, 78,253,255,252,  4,245, 11,253,255, 67,217,  1,  0, 38,210,  6,  0,252,255,252,253,  6, 19, 51,228, 89,253,  3,  0, 20,240,  1,215,232,  5,253,235,247, 11,251,  0,212, 20,255,  2, 42,248,  1,  0,244,243,  1,255,228,227, 61, 11, 39, 29, 39,219, 37, 76,  0,240,204, 38,203,205,218,231, 58, 37, 76, 14,225,124,236, 61, 54, 12,232, 68, 47,241, 31, 84, 31,205,191,218,246, 21,  9, 47,241,210, 20,  3,202, 61,  7, 45, 71, 57, 14,221,  4,135, 42,  0,  6,244, 54,235, 32,225,242,158,203, 45, 20, 41,149, 41, 52,204,225,229,217,244,189,  4,255,  2,  1,255,255,255,  0, 14,247,254,255,243,  5,  0,  6,  0, 63,  8,  9, 34,246,223,  4,  5,  3,255,255,  0,241,250,238,123,
      254,250,255,176,  1,  5,  3,152,253,  0,233,211,247,248,  5, 25,  1,  1,255,251,254,  2, 11, 65,  6,  1,253,  7,  2,250,  5,255,  4,212,  0,  3,244,254, 34,  1, 27,215,255,252, 42,  1,247,  2,242,222,  4, 94,230,253,252, 10, 41,239, 30,  1, 12,  2,211,  8,231,250,  0,  0,  1,210,248,255,211,195,  0,254,255,  7,  2,255, 19, 38,254,231,198,224,  2,254,207,213,230,249,211,243,254,137,  4, 83,221,  6,184,  4, 62, 62, 42,247,200,  1, 28, 20,191,114,228,  8, 26, 55,194,  9,  1,  0,227,  3,250, 11,  0,202, 19,247,  1, 38,225,201, 64, 19,212,142,193, 57,151, 19,  4, 45, 51, 10, 68,103,250,241,  6, 20,  5, 49, 18, 39, 25, 27,228, 48,254,149,  0,  1, 51, 26, 40,  7, 70,190,252,  1,231,  0,229,  5,254, 50,254,  0,254,236,245, 17,  1, 56,  2,253,  2,  2,  9,  4, 53, 66,
    },
    {
       15,212,196,201, 16,  8, 24,126, 15, 11,212,193, 11, 15,  7, 51,167,194,219,169,230,252, 29,249,180,186,218, 37,131,207, 19,230,195,241,248,250,206,204,  3,  1,200,208,  6,252,227,201,  0,  4,211,233,255,255,215,217, 13,254,198,251,  8,219,214,168, 37,253,  9,250,  0,254,245,225,254,  0,  1,238,  0,  1,226,  7,  0,254,  7, 40,255,  1,255,  8,  1,254,248,  3,  1,253,  5,254,246,  0,253,  0,239, 65,247,242,210, 79,255,236,250, 66,251,198,161,  0,183, 79, 80, 27,217, 69,166,207,251,202,197, 59,157, 11,232,244,240,195,210,215,236,255,189, 13,205,203,199,195,217,213,194,188,244,  0,254,  2,  5, 15,182, 49,254,253,  4,247,  4,  3,252,  0,  0,  5,  3,  9,  1,255, 10,  9,251,  3,237,  0,218,252,  5,237,255,255,  1,  0,  2,250,  1,  3,  0,  0,  1,255,  1,  0,  3,254,
      206, 21, 30,245,226,  7,  7, 56,246,  3, 15,251,192,  5,203, 10,200,  1,228, 35,243,234,254, 87,242,201, 86, 10, 15,150, 40,254,  3,  0,  3,255,  0,  1, 10,253,252,245,255,246,  1,243,  0,255, 13,247, 39,  8, 10, 31,255,  2,251,235,240, 50,181, 45,109,242,204,211,191,192,224,206,192,192,197,193,196,192,231,209,193,194,166,226,187,191,233,192,191,191,192, 62,194,191,192,192,186,192,221,173, 40,156,233,204,252,193,164, 15,198,240,221,249,214,185, 81,198, 84, 32,  5,220,182,157,112,157, 43,  8,200,234, 83, 45,215,217,221,239, 10,  0, 13, 41, 13,198,  0,244,  3,233,  1,236, 25,253,  9, 48,  1,  1,  5,244,  1,254, 28,208,249, 14,  1,  1,187,210,254,184,194,232,185,193,193,194,200,183,213,193,186,214,196,192,189,188,194,192,192,192,192,191,193,193,191,192,195,192,
       63, 13, 66, 14,214, 67,184,  3,174, 33,  0,100,200,241,208,246, 12,207, 19,232,202, 20,231,214, 26,238, 11,216,222,252, 13,249,255,230,193,211,254,228,  2,247,125,251,194,208,180,  5,  8, 49,251,248,  0,  3,  4, 10,  2,213,246,230,221,223,  5,230,  1,223, 19,253,  0,  6,235,194,  1,253,  3,127,  0,251,  4,255,  0,252,245,225,  0,  1,231,254,255,  1,  1,254,  0,  1,245,254,252,255,196,253,  0,249,197,188, 67, 29, 16,  0,119,102,193,217,128,  1,195, 16,228,  7,194,  5,  3,236,213,255,244,251,219,244,254, 10,186,208,194,241,250,255, 18,  0,  0,191, 25,209,  6,  7,  1,251,229,241, 20,253,253,255,226,  5, 59,210,248,185,205,222,254,  9,  0,251,  6, 12,254,248,252,  2,  2,  1, 13,255,248,  0,233,  5,  2,251, 35,255,  0,  1,  2,255,  3,250,  0,255,255,255,  3,  1,
       10,  0, 49,243,191,252,215,209,168,180,245,185,254,249,209,  8,175,255,201,218,247,  7,226,243,158,  1,248,236,237,241,252,204, 52,213,  8,204,150,192,226,213,203,215,217, 26,185,  9,124,192,192,196,248, 11,193,223,237,  0,186,173,250,  6,174,184,251,  3,253,254,  0,222,254,  0,  0,220,255, 25,  0,  3,250, 60,  2,255,  3,  3,  6,255,  0,  1,  1,253,255,  1,  2,  2,255,255,255,  0, 39,119,176,  0,  8,240,138,254,164,  1,216,238, 11,195,187, 41,219,122, 35,221,252,252,217,217,184,194, 39,220,234,214,239,  2,139,253, 47,171,175,253,190,195, 58,230,255,211,196,193,200,198, 19,234,242, 36,246,214,248,254,175, 36,214, 77, 12, 28,  6,244,253,221,  2, 26,  2,  2,  1,  1,250,  0,255,  5,  0,  1,  1,  1,254,191,  0,  2,  1, 17,  1,  0, 72,213,251,172,249,252,  0,  1,
    },
  },
  {  //  34
    {
      214,216,208,  5,210,222,236,234,208,207,224,191,215,209,224,236,213,218,211,199,222,222,  5,  1,222,206,198,194,211,222,239,209,209,255, 34,  1,209,255,254,243,202,245,251,251,203,249,251,254,224, 12,  2,  8,205,246,253,  3,197,247,251, 36,214,  2, 18, 28,225,255,254,255,250,  0,255,  0,254,255,  0,  0,253,  0,  7,  6,244,  5,  0,  4, 30,  0,255,  5,170,  4,  1,  7,229,  1,252,  3,192,216,214,208,218,210,222, 10,199,229,  2, 11,222,229,  4,252,210, 11, 13,  5,242, 20,247,  5, 85, 32, 61, 71,207, 16,208, 91,253, 46,246, 24, 10, 26,  3, 54,250,251,  5,  2,233,240,  8,  1,  0,244,255,  1,195,203,187,253, 23,  9,248,250, 47,  2,216, 12,252,229,252,237,  1,255,255,233,  1,247,246,246,  6,251,253,  7,  0,  2,245,244,234,207,  4,  5,  1,247,245,243,  1,  0,  5,  6,
        8,232, 14,224,  0,233, 44,222,255,221, 36,  3,  4,230, 34,231,  2,244, 53,253,255,234, 88,187,  0,222,  8,  0,  0,234,252,240,241,230,247,230,254,244,  0,252,246,250,  3,246,  0,  0,  4,237,238,238,250,191,  0,246,  8, 32,  1,  2,244,250,251,  2,234, 10,223,194,190,194,240,192,194,191,235,191,191,192,219,194,186,187,202,192,191,192,241,192,191,193,251,184,193,192,247,192,200,192,  1, 23, 13,225,253,  1,  3, 12, 10, 29,150,247,  2,254,209,233, 27,226,  7, 53,244,241,243,255,255,188,  8,192, 27,253,210,227,221,226,179,203,210,203,213,219,232, 12, 42,207,238,232,229,236, 15, 34, 34,184,225, 16,247,247, 27,248, 24, 16,226,231,231,236,194,193,186,195,198,218,221,224,190,190,198,190,200,202,203,202,192,187,193,185,198,188,198,202,193,201,190,193,193,202,202,204,
        0, 15,  0,255,242, 38,  2,233,247,  0,  1,  1,225,180,253,249,248, 69, 12,  4, 12,251,  3,246, 13,237, 42,249,  9,  1,  1,247,254,206,  3,  5,247,241, 16,250,161,180, 10,  0,  5,247,245,  5,  0,239,213,213,254,236,207,252,247,247,240,  7,  3, 20,210,253,244,  0,  3,  1,223,253,250,255,251,  2,  2,250,249,  7,248,255,255,250,  0,  0,  0,255,254,255,249,  9,  2,  0,  0,  0,  3,  1,230,185, 57,152, 26,193, 51,187,252,250, 25, 46,243,225,177,200,  1,252,221,250,  5,  2,254,  0,227, 12,255,  7,  3,  0,225,238, 12,231, 35,231, 34, 44,232, 42,209, 55, 96,244,213,246,226,  4,229,252,253,242,  0,221,245, 13,161,197,110,195,224,227,229,219,247,  8,252,  7,  9,245,  1,  3, 17, 18,  6,254,  1,246,  0,240,255,253,254,  4,255,  1,  0,255,  0,253,  3,  5,  0,  0,  0,255,
      230,245,  1,240,179,188,  0, 21,224,204,  0,  8,182,248,253, 17,231,241,243,245,189,244,242,  4, 37, 55,191,248,246,  4,250,252,191,199,185,212,191,195,186,206,119, 19,194,223,154,180,187,205,195,199,234,255,179,203, 44,199,217,194,212,221,196,204,241,204,  2,254,216,  2,255,255,225,254,253,254,194,  2,  4,248,241,252,255,  1,246,  0,255,254,255,255,  1,253, 19,  1,255,  2,240,  1,162,208,  0,246,161, 16,225, 23,167,250,182,189,  9,223,238,224, 26,  3, 88,250,217,227,247,239,207, 47,237,221,248,251,204,  0,102, 32,208,226,195,183,214,186, 67,  3,174,184,177,204,201,208,234,248,243,232,253,239,  1,203,  7, 44,244, 10,254,  2,245,251,201,252,  5,  0,  1,255,251,  3,  3,254,250,  4,  1,  1,  0,  1,  1, 16,252,250,255,252,  2,245,  4,250, 16,  8,  0,255,  8,253,
    },
    {
       35, 17, 64, 59, 42, 23, 58, 48, 18, 15,116, 22,  8, 41, 42,121, 18, 21, 49, 88, 18, 23, 46,  0, 20, 21, 82, 46, 25, 27,253,229,  8,  9, 56, 61, 13, 29,253, 57,  5, 15, 85, 60,  1, 11, 59, 58,  3, 12, 59, 66,  6, 19, 66, 59, 11, 24,208, 69, 12, 16,249,205,  0,255,250, 11,  1,  0,  0, 13,250,255,254,  8,  0,  0,  2,  1,  0,  0,  1,  7,  0,255,  3, 29,  0,  1,  5, 24,  5,  2,255,  3, 28, 42, 24, 89,245,  8,253,237, 60, 39, 71, 48, 66, 71, 52, 62, 19, 13, 22,  9,225, 11, 19,  7,221, 81, 87, 34,209,248, 78, 78,250,  2,  2,  1,255,  5,253,254,253,  1, 24,  7,252,  9,255,  4, 57, 66, 68, 61, 68, 43, 65, 66,  9,101,128,113,238, 24,135,111, 64, 64, 64, 63, 65, 66, 63, 64, 66, 64, 62, 64, 64, 63, 62, 64, 76, 38, 21, 58,  8, 57, 26, 41, 53, 57, 57, 63, 63, 63, 64, 60,
       64,244, 63,207, 76,237,190,202, 61,237,235, 26, 60,255,  2,182, 39, 22,198,  2, 19,237, 10,  0, 71,253,239, 10, 57, 16, 67,174, 32,  7,120,195, 65,  1,145,250, 28,  0, 24,  0, 31,  0,  7,163,246,246,248,  0, 22,  7, 62,  8, 12,  5,  7,  4, 34,  5,252, 19,  2,255,  9, 15,  1,  1,255,236,254,  0,255,  6,255,  0,  7,  3,  1,  0,  4,  5,255,  0,  7,246,255,  0,  1, 20,255,  5,  3, 16, 54,129, 71,102, 94,154, 59,121,247,  7,  7,246, 61, 35, 18,254, 46, 14, 60,218,  3, 18, 74, 70, 72, 10, 51, 21, 62, 66, 53, 54,  6,254,  2,244,254,171, 11,  6,  6,251,  5, 34, 17, 73,  0,206, 16,  5, 10,  5,253, 24, 22,  3,250,  4,  1,  4, 13,  0,  1,  8,253,  0,255,255,  3,253,255,  0,  0,255,  0,  0,  3,  5,  1,  0,231,209,229,  1,  7,  1,  2,  7,  4,  9, 13,  2, 40,253,  6,  1,
      252, 11, 44,236,253,229,100, 60,  5, 19,203,  2,250, 10,  6,191,  5, 13, 26, 13, 22, 45, 15, 24,  4, 65,250, 34,253, 10, 36, 23,  0, 53,255, 15, 20, 65,218, 38, 39, 65,  3, 39, 35, 59,251,118,  4, 56,  0,  5, 46, 52,  5,  3, 21, 58,  3,  7, 38, 56,  7,254, 60, 74, 66, 61, 70, 66, 61, 64, 64, 71, 60, 65, 65, 65, 67, 63, 64, 63, 55, 34, 63, 62, 54, 43, 66, 62, 53, 43, 61, 61, 52, 43,156,156, 17,241,214,251, 12,199, 13, 60,254, 50, 13, 48, 30, 53,223, 15,252,250,247,248,246,249, 17,  2,253, 20,230,  5,252, 25, 14, 44,  2,218,200,211,255, 19, 53, 40,230, 27,243, 59,255,  2,243,244, 15, 14,  3, 25,252, 16, 37, 50, 38, 30, 38, 47, 39, 31,  4,255,253,  0,  5,  1,  0,  1,  6,  3,  0,  0,  1,255,  0,  0,199,252,227,255,255,  0,238,  0,255,  0,  1,  0,  0,  0,  2,  5,
      191, 29, 31, 27,  9, 28,237, 23, 60,  9, 44, 12,237,224,254,246, 10,  1, 34,253, 40,254, 37,251,255,255,249,  1, 42,251, 42,249,  3,249,  4,252,211,208,244,  0,253,255,  2,254,225,225,236,  1,250,  1, 29,  2,  1,  1,  6,  7,  1,254, 12,  9,  1,  1,  6,  6,255,255,250,  4,  2,251,246,254,243,254,  3,  1,  0,243,  1,  0,255,255,  3, 25,  2,  2, 13, 21,  0,  1,  7, 22,  2,  2, 13, 21,155, 98,183,107, 72,220,251,249, 34,  1, 77, 32, 58, 34, 64,229,209, 64, 63, 60, 83, 53, 45, 70,  7,159,  1,  0, 19, 22,250, 23,117,187,175,223, 43, 52,248,247, 55,183,228, 36, 55,202,254,  3,223, 14, 78, 31, 47,243, 69,  2,164, 32, 29, 50, 53, 32, 24, 42, 41,230, 32,254,222,255,  0,  1,  0,  0,  8,  0,254,  0,255,  1,253,245,  0, 20,255,  1,255,240,249,  0,  6,  4,  3,254,  2,255,
    },
    {
      216, 10,235,189, 18, 12, 60,245,  3,238, 25,233, 13,  2,226, 75,212,239, 72,242,221, 52, 65,245, 30,233,213,187,235, 57,228,220,198,254,255,  6, 17, 34, 11,254,226, 12, 17,  5,232,  7,250,253,250,243, 26,254,247, 86, 36,  6, 56, 19, 12, 33, 22,220, 44, 25,  2,  2,  0,252,255, 14,  2,255,254,  2,  1,252,  7,253,  0,253,248,254, 89,  0,245,255,230,  5,  6,235,  7, 12, 29,247, 21,249, 29, 36,  2,231,  0, 60, 57,194,  5, 20, 31,  9,250,202,233,245, 18,  6, 41, 56,  9,224, 76,223,  2,103,184, 81, 74,131,163,  4,254,  3,235,  7,249,  1, 14,  2,203,201,201, 69,209,224,  0,  5,252,  3,250,  8,254, 11,251, 68,  1,255,  4,249,221,255,  3,199,  6,  0,  0,  1,  4,255,  5,252,  1,  5,255, 47, 19,  1,  0,  1,255,254,  0,  0,  2,  2,248,251,  0,  1,  0,255,  3,255,255,241,
       13, 12,223,  9,213,232,237,202,239, 14, 15,  4,230, 44,156,167,215,204,212,255, 49, 11, 30,  1,222, 37,246,255, 39, 53, 52, 15, 30, 37,198,254,253,219,254,  1, 73, 55,  1,248,  2, 67, 35,255,244,250,235,  0,235, 56,231,  2,203, 69, 54,128,254,224, 32,250,  1,235,233,  0,  5, 24,250,  0,254,247,  2,  1,255, 10,  2,  0, 18,  4, 61, 15, 11,105, 41,  1,234, 14,139,  1, 53,  0,252, 37, 37, 56,244,236,194,202,193,110,  4,  0, 18,  8,  7, 15, 37,189,240,241, 60,196,  0,206, 98,247,217,  1, 78,253,237, 32, 77,214, 43,204, 12, 28,248,224,228, 26, 29, 46,  9,  8, 12,241,  0, 20,  0,  0, 12,241,  1,  5,  1, 58,  3,  6,  7,248,164,205,  0, 61,255,254,255,  7,238,228,  3, 16,  4,254,  2,  2,242,255,255,230,  1,  0,  0,  2,  1,255,254,253,  3,  0,  1,255,  3,  4,  8,  1,
       11,217,254,253,251,248,  0,  0,251,237,215,236,254,  9, 61,236, 56,249,173,243, 26, 42,243,  1, 48, 20, 98,186,255, 17,149,  2, 35,247, 58,247, 14,248, 64,  0,206,247,249,  0, 47,207,  4,  4, 66, 67, 14,  2, 14, 69,198,255,212, 48, 94,233, 45, 34, 27,  0,252,253,242,  0,243,247,  1,  0,254,254,255,  1, 19,  1,255,  1, 15,  5,252,  9,233,162,136,255,193,181, 53, 38,  4,255,125,  7,244, 17, 15, 10, 60, 59,  1,237, 30,  2,242, 54,236,102,186, 27,221,248,  1, 67, 14,251,182, 44,  4,124,131,248,244, 21, 25, 78,229, 17,208, 19,241,194,190, 21, 79, 24,  6,245,253,209, 46,198,254,254,  2,  3,249,253, 24,228,229, 33,  6,251,203, 47,197, 25,243,243,219,242, 47,211,253,240,  1,  0,  0,238, 20, 45,  0,207,  0,  0,  0,  0,255,  0,  1,  4,  2,  1,253,  0,254,  2,  7, 14,
      249, 21, 10, 28,230,252,  0,  3,247,  1,  6,  0, 31,  9,  5,247, 17, 14,253,  1,231,192, 23,254, 83,221,217,127, 22,251, 66, 66,252,  1,252,  0,195,  6,191,  9,235, 62,227,254, 43,242,  6,252,  4,234,255,  3,  0, 31, 41,248, 47,161,237,255,  5,229,192, 38, 11, 29,  1,  1,254,255,255,255, 14, 12, 15,252,  6,255,  0,  0,216,  5,  0,255,229, 45,  0,192,  3, 42,228,184,  4,152,  1, 56,248,214, 37, 40,194,  7, 29,248,220, 24, 27,210,249,210,  1,120,  2,253,  1,  1,255,143, 32,  0,249, 76, 91, 64,255,177,177,230, 60, 61,248,245,202,227, 62, 45,219, 55, 51, 24,243, 49,212,238,  1,  5,255,  1,  4, 13, 59,242, 37, 91,255, 75,  3,240,236,  6,  4,  4, 33, 33,245, 86,240,252,254,  0,249,  1,252, 19, 19, 18,  0,  1,  0,  3,  3,  0, 68, 14,  2,  1,  1,  2,  1, 10,  8,254,
    },
  },
};

const int g_nsstHyGTPermut8x8[35][3][64] = 
{
  {  //   0
    { 63,62,55,61,54,47,53,60,46,39,52,45,59,38,31,44,51,37,58,30,43,23,36,29,50,57,22,35,28,15,42,21,49,27,14,56,34,20,41, 7,13,26,19,48,33,12, 6,18,40, 5,11,25, 4,32,10,17, 3,24, 2, 9,16, 0, 1, 8,},
    { 63,55,62,54,61,47,46,53,60,39,45,52,59,38,31,44,37,51,30,23,58,36,43,29,22,50,15,35,28,57,42,21,14,49,27,20,34, 7,56,13,41,19,26, 6,12,33,48, 5,18,11,40,25, 4,10,32,17, 3,24, 2, 9,16, 8, 1, 0,},
    {  0, 1, 8, 2,16, 9,17,10, 3,24,18,25, 4,11,32,19,26,33, 5,12,40,27,20,34,41,48,13, 6,35,28, 7,42,21,49,56,14,36,43,50,29,57,22,15,44,51,58,37,30,23,52,59,45,38,31,60,53,46,39,61,54,47,62,55,63,},
  },
  {  //   1
    { 62,54,63,46,55,60,38,47,52,61,44,39,30,58,53,36,45,31,50,22,28,37,59,42,23,14,51,29,34,20,56,43,15,26,21, 6,48,35,12,25,18,40, 7,27,13,32,19, 5,10,49, 4,24,41,11, 2,16,33, 3,57, 1,17, 8, 9, 0,},
    {  0,10,63,16,55, 2,47, 8,24,18,61,39,32,53, 4,26,31,45,12,37,40,34,20,59,23,29,28,51,48,15,43,42,21, 6,36,50,35,56,14,44, 7,58,22,27,57,13,38,30,19, 5,52,49,11,41,46,60,33, 3,25,54,17,62, 9, 1,},
    {  0, 8, 1, 2,16, 9,17,10, 3,24,18, 4,25,11,32,26,19, 5,12,33,40,27,20,34, 6,13,41,48, 7,28,35,42,21,49,14,56,36,43,50,29,22,57,15,44,51,37,58,30,23,52,59,45,38,31,60,53,46,39,61,54,47,62,55,63,},
  },
  {  //   2
    {  0, 1, 8, 9, 2,17,19,18,10,12,27,25,16,11,20,28,26, 3,36,30,24,32,34,37,33,35, 4,13, 5,43, 6,41,40,22,46,51,38,21,15,45,56,29,48,49,31,53,42,23,14, 7,50,44,55,60,57,54,58,52,59,39,47,63,61,62,},
    {  0,17, 8,18,25, 3,11,19,16,36, 1,34,20,35,28,45,10,24,27,33, 4,54, 9, 2,12,21,52,26,29,50,32, 7,58, 6,40,63,30,43,46,53,38,13,42,56,48, 5,47,41,37,51,23,49,14,31,60,59,15,22,39,55,57,62,44,61,},
    {  0, 8, 1, 9, 2,16,17,18,10,24, 3,19,25,33,12,26,34,11,20,35,32,28,27,13,48,50,36,37, 4,21, 6,41,40,51,52,43, 5,29,14,42, 7,22,49,57,46,53,63,60,30,56,15,45,59,44,31,61,23,58,55,39,38,54,47,62,},
  },
  {  //   3
    {  0, 8, 9,17,10,18,16, 2, 1,25,26, 3,34,32,28,24,35, 4,33,27,37,40,12,13,19,20,36,42,11,41,51, 5,52,43, 6,48,23,21,50,44,49,29,14,53,56,58,45,31,22,30,57,15,38,59,55,63,39,60,62,47,61,46,54, 7,},
    {  4,10,18, 9,35,27,60,33,36,16,28,17,53, 2,11,24,25, 3,62,43,32,19,26,20,52, 1,12,42,47,58,50, 0,45,63,56,51,30,38,13,59,40,41,46,21,34,44,15,48,37,31,55,22,49,14,57, 8, 7,61, 6,23,54, 5,29,39,},
    { 61,53,62,54,47,46,45,60,39,31,37,29,52,23,59,38,36,55,35,21,30,51,22,43,44,19,27,63,14,32,28,50,13,42,58,20,49,25,18,12,15,48, 6, 7,17,11,41,34,26, 0, 9, 4,56,40, 5,10,33,57,16, 2, 8, 3, 1,24,},
  },
  {  //   4
    {  0, 8,16,17,18, 9,10,26,34, 1,35,32,24, 2, 4,50,25,11,33,36,28,48,51,44, 3,27,29,12, 5,38,41,49,19,20,40,43,14, 6,52,13,42,58,21,59,54,62,53,57,47,22,23,46, 7,56,37,15,60,30,55,63,39,45,31,61,},
    { 52,60,55,44,61,45,54,53,28,46,29,63,20,30,39,11,10,21,38,36,47,27,62,43, 1,24,26,22,34, 7,37,41,31,58,12,15,13, 2,57, 6,48,32,51, 4, 9,56,49,16,25,35,19,59, 5,17,18,23, 3, 8,50,33,14,42, 0,40,},
    {  0, 8,16, 1,10,18,32,25,17,27, 9,34,12,33,19,11, 2,26,35, 5,20,48,42, 3,24,37,50,44,41,28,21,22, 6,52, 4,49,40, 7,43,14,51,53,30,39,46,57,13,31,23,62,60,29,56,38,58,63,59,54,15,55,45,61,47,36,},
  },
  {  //   5
    {  0,40, 9,10,16,26,17,32, 2,11,12,34,35,24,37,28, 1,18,25,20,36, 3, 5,19, 4,48,27,13,42,51,44,57,49,46,29, 8,30,56,50,38,53,43,60,41, 7,52,22,58,14,55,61,63,59,15, 6,33,62,23,47,21,31,54,45,39,},
    {  0, 8,17, 9,10,18,34,35,26,36, 1,16,44,27,24,33, 2,52,19, 3,32,11,41,25,20,29,53,40, 5,42,61,58,46,43,51,12, 4,48,13,21,56,28,38,30,55,59,50,45,49,57, 7,62,15,54,60,47,22,37,63,23,14, 6,31,39,},
    {  0, 2, 8, 5,16, 9, 1,10,17,18,12,24,26,34,11, 3,25,27,36,20,19,43,35,28,29,44,48,40,45,52,33,42,21,30, 6,14, 4,32,50,37,41,56,61, 7,22,49,15,13,57,31,53,47,62,58,59,23,51,46,38,39,54,63,60,55,},
  },
  {  //   6
    {  0, 8, 1,17, 2,16,11, 9,26,34,10,25, 3,33,19,12,32,24,18,43,27,40,35,28, 5,50,36,42,20,29, 6,41,44,21, 4,51,48,53,57,13,49,30, 7,61,59,52,37,60,45,23,47,58,31,14,22,56,62,46,63,38,54,55,39,15,},
    { 62,54, 9,46,11,17,52,28,33,25,36,35,30,44, 3, 4,38,27,37,60, 1,53,24,22,19,16,29,51,32,13,61,21, 6,20,49,43, 0,48,12,41,14,57,34, 5, 8,56, 7,58,47, 2,59,10,63,40,45,15,50,26,42,23,39,31,18,55,},
    { 63,55,54,47,53,39,62,31,52,46,38,37,61,36,51,30,45,20,21,50,60,19, 6,29,23,35,28,44,43,58,49,15,22,34,27, 2, 7,26,17,12,48,59,13,42,57,33,32, 5, 3, 4,56,11,14,41,18,24,16,10,25, 8, 0, 1,40, 9,},
  },
  {  //   7
    {  0, 8,16, 9, 1,24, 2,17,10,33,25,34,26,32, 3,27,18,41,35,11,20,19,40, 4,28,49,36,48,12,43,51, 5,56,13,52,30,21,42,44,50, 6,57,37,14,29,59,60, 7,22,58,61,46,53,15,31,38,45,23,47,39,54,55,62,63,},
    {  0,40,17, 1,18, 9,33,16,10, 2,25,51, 3,52,34,20,24,43,44, 8,19,60,35,21,50,26, 4,27, 5,36,11,61,42,12,22,49,56,45,38,48,28,57,29,54,53,32,59,37, 7,13,41,30,58,23,39,46,15,31, 6,63,47,55,62,14,},
    {  4, 8,16,17, 1, 2, 9,26,33,18,34,25,32, 3,27,50,57,10,24,44,51,40, 0,20,12,49,19,11,48,28, 5,35,13,36,52,58,56,43,29,14,53,37,46,45, 7,38,21,22,30,61,41,42,59,62,31,60, 6,54,47,23,39,63,15,55,},
  },
  {  //   8
    {  0, 1, 8, 9,17,16,18,10, 2,24,25,11,33, 3,34,26,32, 4,35,12,27,40,36,20,19,49, 5,28,50,52,48,44,51,29,13,42, 6,43,41,45,21,57,30,59,22,56,58,37, 7,38,46,53,61,62,60,14,15,31,23,54,47,39,63,55,},
    {  9,17,62,11,28,25,46, 1,54,44,14,33,38,60,13,52, 3,40,41,49,36,16,19,35, 6,12,27,32,10,51,45, 8,29,59,30,56,43,21, 2,22,15, 0,42, 4,24,48,39,37,57,61,34, 5,20,53,31,50,58, 7,47,63,23,55,18,26,},
    {  8, 0, 9,16,25,34,35,17,42,26,51,24,10,60,48, 3,33,50, 2,59,44,49,12,41,27,18, 1,11,28,58, 5,19,52,36,32, 4,21,45,40,20,13,57,43,53,46,61,30, 6,29,37,38, 7,56,54,15,14,23,31,62,39,22,55,63,47,},
  },
  {  //   9
    {  0, 8,16,17, 1,10,32,25,33,18, 3,34, 9,26,24,59, 4,50,49,48, 2,52, 5,36, 6,35,27,41,29,28,22,42, 7,40,11,37,20,19,53,43,44,23,58,12,56,57,51,38,61,21,45,14,60,13,30,39,55,31,54,15,62,63,46,47,},
    {  0,25, 8,18,33,34,16, 1, 9,51,32, 2,26,35,50,10,17,44,24, 4,49,52,48,19,27,40, 3,57,12,36,43,29, 5,59,28,42,56,45,37,11,60,20,13,41,53,58,21,14, 7,22,38,61, 6,30,46,23,54,31,63,15,55,39,62,47,},
    {  0, 8,16,24, 1,17, 9, 2,34,33,25,26,41, 3,32,18,51,10, 4,35,11,48,50,40,20,27,28,12,52,49,19, 5,36,42,59,45,60,56,21,58,37, 6,29,44,57,22, 7,43,30,13,53,14,38,23,61,31,47,39,46,63,62,54,15,55,},
  },
  {  //  10
    {  0, 8,16,25, 1,34,18,33,10, 9, 2,24,32,50, 3,17,41,26,51,11, 4,44,35,56,48,27, 5,12,19,58,45, 7,40,20,21,57,36,49,43,13,28, 6,42,59,29,53,60,37,61,14,22,30,52,38,15,62,23,31,39,54,46,55,47,63,},
    {  4, 9, 1,42,24,25,51,10,50, 0,19,18,33,32, 8, 2,26,12,35, 5,58,17,34,16,27,13,40,44,48,20,60,45,49,59,57,29, 6,21, 3,11,52,28,56,37,41,43,36,22, 7,15,53,14,30,38,61,23,31,46,47,54,39,62,55,63,},
    {  0, 8,25,16, 9,33,18,24,10,50,41, 3,17,42, 1,26,51,48,12,19,40,27,34,35,37, 2,52, 7,20,21,11,44,57, 5,32, 6,36,58,13,56,60,49,43,14,29,28,59, 4,53,45,38,22,15,61,30,46,23,54,31,62,39,55,47,63,},
  },
  {  //  11
    {  0, 8,16,17,33, 9,42,32,25,10,24,18, 1,50, 3,26,49,34,35,12,57,20,19, 5,27, 2,48,36,51, 7,44,40,58,11,21,60, 4,37,43,59,28, 6,13,41,52,53,30,61,56,22,29,62,45,31,15,23,38,14,54,39,46,63,47,55,},
    {  4, 0,25,33,34,24,49,50,10,41,58, 1,32, 9,17,48, 3,27,51,40, 8,43,20, 2,16, 7,26,19,18,28,35,36,44,12,11,42, 5,52,57,59,29,21,56,13,14,60,37, 6,45,22,53,30,61,15,38,31,23,46,55,47,54,39,63,62,},
    {  0, 8,16,25,32,33, 9,58,17,42,49,59,26,11,51, 1,34,43,24, 4,19, 7,36,40,48, 2,27,18,35,60,20, 5,10,41,21, 3,44,50,52,28,57,56,37,13, 6,12,29,45,22,53,61,38,30,14,15,46,23,54,31,47,55,39,62,63,},
  },
  {  //  12
    {  0,16, 9, 8,24, 2,25,32, 3, 1,18,10,42,19, 4,17,50,33,34,26,11,40,13,28,41,49,21,48,51,27,35,58,20,44, 7, 6,36,45,12,56,29,43,52,59,37,22,57, 5,15,60,23,14,61,38,30,53,54,31,46,39,63,62,47,55,},
    {  8,25,24,33,49,16,18,41, 1,58,57,50,32,35,20,17,10, 9,42, 5,27,51,34, 0,26,44,19, 4,43, 3,40, 2,12,11,48,52,56,13,28,59,37,36,60, 7,29,45,22,53,21,38,14,61,30, 6,23,46,15,54,31,62,39,47,55,63,},
    {  0, 8,16,25, 1,33,34,49, 9,58,24, 2,42,17,32,18,10,57,50, 3,51,43,48,11,26,41, 4, 5,27,35,12,40,20,56,60,36,52,19,59, 7,13,44,28,21,29,45,53,61,14,37,22, 6,46,30,38,15,23,31,47,39,62,54,55,63,},
  },
  {  //  13
    {  0, 8,24,16,32, 1,33, 9,25,49,48,17, 2,10, 3,34,26,42,18,58,57,40,11,50, 4,19,27,56,41,20,12, 5,51,43,59,35,28, 7,36,13,45,52, 6,29,21,44,14,60,37,22,30,15,53,38,23,46,61,31,39,54,55,47,62,63,},
    {  8,16,17, 0,33,32, 1,49, 9,24,34,18,40,50, 2,25, 3,42,48,41,19,26, 4,57, 7,58,11,35,51,12,10,27,20,36,56, 5,44,59,13,43,28,21, 6,52,37,29,14,60,22,45,30,15,53,38,23,61,46,31,54,39,62,47,63,55,},
    {  8,24, 0, 9,32,17, 2,33,16,11,26,48,49,12,34,40,19,41,50,25,13,20,35, 1,42,18,51,21,59,14, 7,52,36,43, 3,10,57,45,56,27,22, 4,28,58,29,38,37,44, 5, 6,23,60,53,15,31,30,54,61,46,39,63,47,62,55,},
  },
  {  //  14
    {  0, 8, 1,16,24, 9,17,32, 2,25,18,40,10,33, 3,41,26,48,11,49, 4,42,19,34,43,20,27,56,12, 5,13,57,50,28,59,51,21,36, 6,58,29,44,14,37,35, 7,38,22,46,23,52,15,60,30,45,53,31,61,47,54,63,55,39,62,},
    {  0, 8,16,24,32, 1,40,17, 9,25,33, 2,48,18,10, 3,49,41,26,56,19,57, 4,58,11,42,27,34,12,35, 5,20,50,13,28,43,36,21,51, 6,59,44, 7,29,22,52,60,14,45,37,15,30,23,38,53,46,31,47,39,62,61,54,63,55,},
    {  8,16, 1, 0,24,17, 2,32, 9,18,33, 3,34,49,12,19,25,40,48,10,13,43,27,41,28,35,11,50,56,21,58, 4, 6,44,20,57,26,42,45, 7,36,52,14,51,53,29,59,23,61, 5,15,22,38,30,46,60,31,54,37,39,47,55,62,63,},
  },
  {  //  15
    { 63,55,47,62,39,31,54,61,23,46,38,15,53,60,59,45, 7,30,52,37,22,44,58,29,14,51, 6,21,36,57,43,50,28,35,13,42,20, 5,27,56,49,34,12,19,41,26,48,25, 4,18,40,24,32,33,11, 3,16, 8, 1, 2,17, 9, 0,10,},
    {  0,16, 8, 1,24,32,17, 2,33,25, 3,49,40,18, 9, 4,11,48,42,41,19,20,10,34,56, 5,57,50,26,12,35,27,36, 7,21, 6,51,58,13,28,44,14,43,37,52,22,29,23,38,59,15,45,30,54,60,31,46,39,53,61,47,63,62,55,},
    {  0, 8,16, 1,24,32,17,40,10, 9,48,18,33,11, 4,25,34,49,41,19,26,35,28, 2,57, 5,50,58, 3,20,27,21,12,36,56,42,13,43, 6, 7,22,37,52,51,44,30,59,29,38,23,14,15,45,46,53,31,60,39,54,47,61,55,63,62,},
  },
  {  //  16
    { 47,31,45,63,29,46,15,44,55,30,61, 7,60,13,14,43,53,23,62,27,41,39,28,37,54,59,12,21,42,38,57,11,52,40,22, 9, 5,10,51,58, 6,49,20,25,35,36,26,56, 3,33,50, 8, 4,19,24, 1,34,32,18,48,17, 0,16, 2,},
    {  0,16, 8, 1,32,17,24, 2,40, 9, 3,18,48,33,26,19, 4,56,25,41,10,11,20, 5,49,34,35, 6,12,42,21,27,57,36,50,28, 7,22,43,13,37,14,51,58,44,23,29,59,45,38,30,15,54,31,53,60,46,39,61,52,47,55,62,63,},
    {  0, 8, 5,16,17,24, 2, 9,32,25,19,40, 3,10,33,48,18,11,41,56, 4,26,34,49,20,27, 6,35,12,36,43,57,22,50,28, 1,21,42,51,44,14,13,58, 7,30,38,54,46,52,15,37,29,59,31,23,45,60,61,39,47,53,63,62,55,},
  },
  {  //  17
    { 63,55,60,47,52,61,39,44,31,53,36,62,23,45,29,15,54,37,57,21, 7,46,30,49,38,59,13,41,22,28,51, 5,33,56,43,14,20,27,48,35,58,40,12,25, 6,19,50,32,42,17,11,24, 4,34, 9,16,26, 3, 8,18, 1, 0,10, 2,},
    {  0,16, 1, 8,17,32, 2,24, 9,11,18,40,33,19,10,48, 4,25,56, 5,34,35,20,49,26, 3, 6,12,41, 7,13,36,42,27,21,22,43,57,50,28,37,14,51,44,23,29,52,58,15,38,53,30,59,39,45,46,31,54,60,47,55,61,62,63,},
    { 62,54,46,63,38,55,61,31,47,60,53,23,39,59,52,15,45,30, 7,37,51,44,58,22,29,36,43,57,50,14,42,21,35,28,56, 6,20,27,49,13,34,41,19,12,48,26, 5,33,40,11,18,24, 4,32,10,16,25, 3, 8,17, 2, 9, 0, 1,},
  },
  {  //  18
    { 62,54,61,46,38,53,63,31,45,55,60,37,23,15,59,52,30,47, 7,44,51,39,58,22,29,36,43,14,21,57,50,28,35,42,20,13,26, 6,49,34,56,18,33,12, 5,19,48,41,27,11,40,24, 4,32, 1, 0,17, 3,10,16,25, 8, 2, 9,},
    {  0, 8, 1,16, 9,24, 2,17,32,10, 3,25,40,18, 4,11,26,48,33,19,12,34, 5,56,27,41,20,35,13,28,42,49, 6,21,43,36,50, 7,29,14,57,22,44,51,37,15,30,52,58,23,45,38,31,59,53,39,60,46,54,47,61,55,62,63,},
    { 63,55,54,47,53,39,62,52,23,45,46,51,31,44,38,61,15,50,30,37,43,60,49, 7,22,29,36,59,42,48,28,35,21,14,58,41,27,20,13,34,57,40, 6,19,26,12,33,56, 5,18,11,25, 0, 4, 1,10,24, 3, 9,17, 8, 2,16,32,},
  },
  {  //  19
    { 63,47,55,39,36,31,60,61,23,22,15,53,52,38,34,44,28, 6,14,45,54,59,37,29,46,42,62, 4,58,20,57, 7,51,32,21,30,35,16,48,33,43,13, 5,12,56,11,49, 0, 3,41,40, 9,27,50, 1, 2,24,10,26,17,25, 8,18,19,},
    {  0, 8,16,24, 1,32, 9,17, 2,40,48,25,10,18, 3,33,56,11, 4,26,19,41,12, 6,27,34,20,49, 5,13,42,35,28,21,57,43,36,50,14,22,29, 7,15,37,51,52,30,58,23,45,44,38,59,31,46,39,53,47,54,60,61,55,62,63,},
    {  0, 8,16, 1,32,24,25, 2,40, 9,26,33,48, 3,17,10,56,19,34,41,36,11,35,49,42, 5,20,18, 6,27,57,21, 4,12,51,43,13,28,14,22,50,37,38,59,31, 7,29,58,30,52,44,53,15,60,39,23,45,54,46,63,55,47,61,62,},
  },
  {  //  20
    { 63,55,61,47,53,60,39,45,31,54,59,23,37,62,44,15,51,29,38, 7,43,21,52,57,30,35,13,46,49,56,22,27, 5,36,41,58,14,50,19,28,33,40,11,20, 6,25,48,42, 3,12,34,17,26,32, 9,18, 4,16,24, 1,10, 0, 2, 8,},
    { 46,62,45,54,47,14,38,63,44,31,61,15,43,55,60,23,42,13, 7,39,59,41,53,12,30,37,58,22,11,40,52,36,21, 6,51,57,10,29,35,28,50,56, 9, 5,27,34,49, 8,20,26,19,25,48, 4,33,18,24,17,32, 3,16, 2, 0, 1,},
    { 61,47,53,55,31,45,39,63,59,29,51,23,37,57,15,43,21,27,49, 7,62,13,41,35,54,25,19,52,46,33,22, 5,11,44,36,17,38,50,42,28,56, 1,30,58, 9,40,14, 2,20,60,24,48,10, 3,32,26, 6,12,16,18,34, 8, 4, 0,},
  },
  {  //  21
    { 63,55,62,47,39,54,61,31,46,53,38,23,60,15,45,30,52,37,22,59, 7,44,29,51,14,36,21,58,43, 6,28,13,50,35,57,20,42,27, 5,49,34,12,19,56,41,26,11,48,18, 4,33,25,40,10,17, 3,32, 9,16, 2,24, 0, 1, 8,},
    {  8, 1, 0,16,24,17, 2,32,12,25,18,40,33,11, 9,19,56,10, 4,48,34, 3,26,41,20, 5,31,35,49,44,57,58,43,27,36, 7,42,21,13,51,30,29,46,28,59,37, 6,52,15,50,45,14,38,47,53,23,39,22,60,63,54,61,62,55,},
    {  0,24,16, 8, 1,33,25,48,32,40,17, 9,56,26, 2,41,19, 3,57,11,18,34,51,49, 5,10,36,27,13,50,42,20,43,21,28,59, 4,29, 6,58,12,44,14,52, 7,35,37,22,23,30,45,46,15,60,53,31,38,39,47,63,61,54,55,62,},
  },
  {  //  22
    {  0, 8,16, 1, 9,24,17, 2,32,10,25,33, 3,18,40,41,26,11,49,19, 4,48,57,12,27,35,34, 5,50,20,13,42,56,28,43,36,21,58, 6,29,44,14,59, 7,51,37,52,22,45,15,30,38,60,23,61,46,31,53,54,39,63,47,55,62,},
    {  0, 8,16, 1,32,24,25,33,10,41,18,17, 3,49, 2,48,40, 9,19,43, 4,20,11,56,26,34,27,12,57,58,50,13,51,28,36,35,21, 7, 6,42, 5,52,29,22,59,44,37,14,38,23,60,39,30,53,15,46,61,54,45,31,47,63,62,55,},
    {  8, 0,16,24, 1,40,17,32,10, 9,48,33,18,25,11,49,26,27,41, 4,34,50,35,56,57,20, 3, 5,51,59,36,21,58,14,19,12,42,13,52,28,43,22,37, 2,15,30, 6,29,23,45,44,53, 7,38,47,46,31,60,39,54,61,63,55,62,},
  },
  {  //  23
    {  1, 8,16, 0,24, 9,32,25,33, 2,17,34,10,18,42, 3,26,50,49,48,11,41, 4,19,56,27,51,12,40,43,28,57,20,58, 5,13,36,21,59,35,60, 6,29,44,52,14,38,22, 7,37,30,45,15,53,46,31,23,39,54,47,61,63,62,55,},
    {  0, 9,16, 8, 2,40,33,17,19,24,18,25,34,12,32, 3, 1,49,10, 5,48,36,27,35,50,42,11,41, 6,21,28,26, 4,51,58,43,22,57,20, 7,37,44,13,56,59,38,53,31,52,14,60,45,47,54,29,30,46,63,15,23,61,55,39,62,},
    {  0, 8,16,17, 1,33, 9,32, 2,24,25,49,10, 3,42,48,26,41,18,34,11, 4,58,50,35,27,19,12,56, 5,20,13,57,59,40,51,28,36, 6,44,21,43,29,37,14,60,38,52,22, 7,30,15,45,46,23,54,39,31,53,47,62,61,55,63,},
  },
  {  //  24
    {  0,16, 8,25,33, 1,24,32,41, 2,17,26, 9,34,10,49, 3,42,18,48,27,50, 4,58,19,11, 5,43,56,36,12,40,35,20,57,51,59,21,45,28,52,44,13, 6,30,29,60,14, 7,53,22,39,46,15,54,38,23,61,31,47,62,55,37,63,},
    {  2, 8,16,24,33, 0,25,32,17, 9, 1,26,34,48,42,18,10,19, 3,58,49,56,28,43,27,51,50,40,57,11,59, 5,41, 4,13,37,12,35,44,29,20,60,21,61,31,52,45, 6, 7,15,22,53,30,36,38,23,14,47,39,54,55,62,63,46,},
    {  9, 8, 0, 1,25,18,16,32, 2,34, 3,24,17,41,27, 4,40,10,11,43,33,13,26,20,59,56,50, 5,19,15,58,22,49,57,29,36,51,39,37,31,42,48,52,12, 7, 6,30,46,35,28,44,21,23,53,63,14,45,54,61,62,55,60,38,47,},
  },
  {  //  25
    {  0, 8,41,16, 1,24,32,58, 9,34,33,25,26, 2, 3,49,20,56,18,51,10,17,36,11,40,19,42,35,27,37,43,29,28,12,48,50, 6,57, 4, 7,47,59,22,53,52,55,45,54,15,38,46,21,60,30,44,31, 5,61,23,63,39,13,14,62,},
    {  0, 9,24, 8,25,32,17, 2,41,10,16, 3,11,43,40,50, 6,18,57,26,27,52,19,48,47,42,36,51,33,49,59,28,14,34,60,56, 1,58,55,39,29, 4,30,35,15,21,37,45,20,22,38, 5,12,23,63,44,46,54,61, 7,53,13,62,31,},
    {  0, 8,16, 9,24, 1,33,17,32, 2,10,34,18,26, 3,40,25,11,42,27,49,19, 4,57,50,51,20,12,56,59,29,28,43,48,58,60, 5,35,22,41,37,21,13, 6,52,36,38, 7,53,23,14,61,31,30,44,55,15,46,63,39,45,47,62,54,},
  },
  {  //  26
    {  0, 8,17, 9,24, 1,25,10,32, 2,34,18,33,16,42,26,11,41, 3,49,44,50,51,28,19, 4,48,43,27,12,20, 5,59,36,35,58,13,21,56,60,29,57, 6,45,14,46,40,23,53,30,52, 7,37,61,15,22,55,31,38,63,54,62,47,39,},
    {  0, 8,16,41,25,32,51,17, 1,34, 9,50,10,27,42,24,36, 2,18,49,58, 3,11,48,20,26,56,43,52,19,28,44,12,57,59, 4,21,35,29,40,13,33,53,60, 6,37, 5,45,30,22,23,61,54, 7,15,31,14,38,39,62,47,63,46,55,},
    {  0, 8,25,42,16,40,33,26, 9,32, 1,35,27,17,24,49,34,51,10, 4, 3,11,57,18,45,21,58,19,20, 2,44,43,48,28,56,36,37,41,13,50,59,14,61,60,29, 5,52,22,15,46,12,39,54,53,30,55, 6, 7,63,23,47,31,38,62,},
  },
  {  //  27
    {  0, 8,17,16, 1, 9,34,25,35,24,10,18,33,32,51, 2, 3,50,11,26,36, 4,41,27,48,19,43,49,52,12,28, 5,57,20,40,44,58,29,37,56, 6,13,22,59,38,45,42,53,61,60, 7,30,54,21,46,15,39,14,31,23,62,47,63,55,},
    { 63,62,47,55,46,54,31,61,30,45,60,39,38,53,44,23,37,29,59,28, 6,52,51,21,13,58,36,43,27,22,15,42,35, 4,57,20,50, 7,19,12,49,34,26, 5, 3,14,25,41,56,18,11,10,48,33, 9,16,40,17, 1,32, 2, 0, 8,24,},
    {  0, 8, 1,17, 9, 2,16,18,32,33,24,10, 3,34,26,19, 4,25,40,49,35,20,11,41,12,51, 5,27,36,28,21,42,50,48,57, 6,13,60,43,22,56,37,44,29,58,14, 7,45,52,59,46,53,30,31,15,61,39,23,47,54,38,55,63,62,},
  },
  {  //  28
    {  0, 8,16, 9,24, 2,33,18,25, 1,26,34,32,11,41,17,10, 3,51,40,35,36,20,12,19,57,48,49,42,27, 4,28, 5, 6,45,21,43,52,56,29, 7,50,59,47,61,58,60,53,13,44,22,39,37,23,15,30,14,31,62,54,38,55,63,46,},
    {  0, 8, 1, 9,16,17, 2,18,24,10,33,26, 3,32,25,11,34,35, 4,20,27,41,28,19,50,40,51,49,12,36, 5,48,14,52,56,44,21,22,42,59,13,38, 6,23,43,37,60,39,57,58,29,46,53,45, 7,31,54,15,47,30,55,62,61,63,},
    {  0, 8,16, 9,33,10,25,18, 3,32,17,35,34, 2,26,52,24,50, 1,41,19, 4,51,28,11,40,36,48,27,29,21,58,43,20,61,14,60,49,13,12,53,59, 5,44,22,42,56,45,46,54, 7,37,38,30,57,23, 6,55,47,31,15,63,62,39,},
  },
  {  //  29
    {  4, 8, 1, 9,16,17,10, 2,34,26,25,43,40,11,18,24,33,19,32,44,27, 0, 3,41,20,28,51,13,37,58,50,60,52,49,12,45,47,29,48,57,21,56,22,36,62,31,35,42, 5,30,59,55, 6,39,61,46,53,38,23,63,14,15, 7,54,},
    { 62,55,46,29,30,47,54,61,13,12,38,63,45,28,31,60,39,44,36,53,14,15,27, 2,25, 4,37,22,58,56,11,32,21,41,34,59,52,17,33, 7, 3, 8,43,26, 5,42,20, 6,49,10,35,19,40,24,57,23, 9, 1,48,50,18,51,16, 0,},
    {  0,16, 1, 8,17, 9, 2,18,40,33, 3,24,10,34,12,19,26,35,25,20,32,51,11,28,27, 5,36,49,50,59, 6,43,29, 4,21,48,44,57,42,52,54,38,56,46,60,61, 7,41,58,14,15,13,23,37,45,22,30,53,47,39,62,31,63,55,},
  },
  {  //  30
    {  0, 8, 1,17, 9, 2,16, 3,32,19,18,27,10, 4,34,24,35,33,43,20,26,25,11,21,40,44,36,58,50,49,13, 5,57, 6,12,28,52,42,56,22,51,30,15,62,41,29,46,38,23,37,14,48,61,59,53, 7,54,60,31,47,45,39,63,55,},
    {  8,24, 2, 9,17, 1,18, 0,26,25,32,34, 3,33,11,16,19,28,36,42,35,48,27,37,44,13,50,49,20,43,14,12,21,41, 6,40,52,10,22,29, 5,53,45, 7,15,51,59,60, 4,56,57,38,58,46,47,30,61,23,39,63,55,31,54,62,},
    {  4, 8,18,63,26, 2,16,47,34,61,55,24,39,31,10,45,32,37,53,21,12,44, 0,13,29,42,52,23,28,48,36,59,40,51,57,35,20,27,25,30,50,60,43, 6,14, 5,46,56,11,58, 1, 7, 9,19,15,54, 3,17,22,38,41,49,62,33,},
  },
  {  //  31
    {  0, 8,16,26, 9,34,17,35, 2,25, 1,24,18,44,11,33,32,10,52,27,53,51,36,28,40, 3,20,21,49,19,42,50, 4,46,37,48,13,29,60,58,57,12,56,45,23,54,47,22,43,61, 5,62,14,41,59,15,63,31, 7,38,30, 6,39,55,},
    {  4, 9,18, 1, 8,34,35,26,24,36,17, 3,44,10,25,16,53,32,11,52,33, 2,19,28,40,27, 0,46,51,29,42,21,41,20,60,39,12,13,43,59,30,37,62,22, 5,45,57,48,50,47,58,38,14,55,63,56,61,23,49, 6,15, 7,54,31,},
    { 62,54,63,55,53,47,46,37,39,61,60,45,21,38,31,52,51,29,30,19,44,23,36,35,50,27,59,28,49,15,22,17,43, 5,57, 3,13, 4,56,34,58, 6,20, 1,16,10,41,48, 7,18,12,32,24,42, 0,14,11,33, 2,40,26, 9, 8,25,},
  },
  {  //  32
    {  0, 8, 9, 1, 2,16,10,18,17,34,26,35,25,24,33,12, 4,19,27,42,20,44,36,28,41,32,29,60,11,53, 5,51,50,22,48,13, 3, 6,52,57,40,45,21,58,56,47,38,46,61,15,23,62,55,30,49,43,63,39,54,59, 7,31,37,14,},
    {  8, 0,16, 9,17,18,10,35,26,36, 1, 3,34,44,20,33,25,24,32,27,52, 5, 2,40,19,41,53,12,50,28,11,21,37,30,38,13,51,48,43,29,49, 4,45, 6,42, 7,62,59,47,54,63,60,61,57,58,46,55,56,31,23,22,14,15,39,},
    {  4, 8,17,27,33,18, 1,34,35,16, 2, 9,44,28, 3,45,10,24,11,36,53,20, 0,49,51,26,25,12,47,38,30,52, 5,29,50,19,56,14,58,55,41,61,22,13,59,60,57,37,40,48,54,21,39, 6,43,42,31,46,15,32,62,63, 7,23,},
  },
  {  //  33
    {  0, 1, 8, 9,17,16, 2,10,19, 3,32,18,11,24,27,25,36,33,20,34,44,35,12,53,28,26,13,49, 6,21,40,41,50,51,52,29, 4,42,22,43,30,48, 5,14,58,60,38,31, 7,56,37,15,57,62,59,47,46,23,61,54,55,39,45,63,},
    {  0, 9,16,17, 8,27, 1,36,18,45,11,26,34,28,53,35,29, 2,52,10,25,33,51, 3,32,24,40,54,42,46,19,21,41,55,12,44,50,49,59,23,22, 5,62,13,48, 6,30,20, 4,63,57, 7,15,43,58,38,60,37,61,56,47,14,31,39,},
    { 63,55,61,62,45,47,46,54,37,53,60,39,38,29,31,52,44,21,43,27,58,59,30,35,19,36,42,22,28,50,11,34,23,15, 5,51, 3,25,57,14,48,41,49,33,13,20, 7,40,26, 1, 4,24, 9, 2, 6,32,12,17, 8,16,56,10,18, 0,},
  },
  {  //  34
    { 63,61,55,53,60,47,45,52,46,54,37,39,62,44,36,38,59,31,50,29,30,58,28,43,23,51,35,27,21,26,15,20,42,57,18,13,34,22,56,19,41,12,49, 7, 6,10, 5,40, 2, 3, 4,48,32,33,14,25,11,17, 1, 0, 8, 9,16,24,},
    { 63,62,39,54,61,47,46,37,55,60,53,38,45,59,52,23,30,29,44,35,36,31,58,22,28,43,27,50,51, 4,57,21, 7,19, 6,42,34,20,33,12, 2,10,41, 8, 5,15,56,14,26, 3,13,49,18,40,11,24,25,48, 1,17, 9, 0,32,16,},
    {  0, 1, 8, 9,16,18, 2,10,17,12,34,25,19,32,36,11,26,37, 3,28,33,40,35, 4, 6,27,13,24,21,30,22,44,42,29,45,20,51, 5,52,57,56,38,54,46,61, 7,41,58,14,15,63,53,59,23,62,49,60,43,47,50,31,55,48,39,},
  },
};


#else

const Int g_aiNsst4x4[12][3][16][16] = 
{
  { // 0
    {
      {    246,    -36,    -35,     -2,    -41,      9,      6,      0,    -27,      6,      6,     -1,     -2,      0,     -1,      0 },
      {    -42,   -220,     39,     26,   -100,     20,     18,     -3,     28,     38,     -9,     -6,      5,     -2,     -2,      1 },
      {    -18,    106,      3,    -17,   -222,    -16,     51,      4,     24,    -19,     -6,      3,     23,      1,     -7,      0 },
      {    -30,    -36,   -150,     17,    -11,   -181,      6,     30,    -59,     17,     45,     -4,      9,     25,     -2,     -5 },
      {    -30,     -3,   -192,     13,     -8,    150,      2,    -30,     50,      4,     29,      0,      5,    -24,     -2,      6 },
      {    -31,     -7,     18,      7,    -24,     76,     13,    -15,   -230,    -33,     49,      4,     26,     -7,     -6,      3 },
      {      3,     22,     10,     16,     47,      8,    208,      1,     -6,    127,      2,    -28,     13,     -6,    -43,      1 },
      {     -2,     38,     11,    -11,    -28,     11,   -129,     -6,    -23,    204,     10,    -42,     46,     -3,     10,      3 },
      {     11,     21,      8,    215,      6,     -2,    -18,      5,     11,    -31,    -41,    -37,    116,      9,    -23,      3 },
      {     -4,     25,      4,    117,    -27,      1,     -5,      3,    -21,     28,      4,    -34,   -209,    -38,     54,      7 },
      {     -9,     -7,    -48,    -33,     -2,    -13,      0,    -29,    -55,      6,   -233,     -8,     -7,    -55,     -5,     17 },
      {      4,     -4,     14,     -3,      3,    -36,     -7,      5,     13,    -11,     54,     11,     32,   -238,    -27,     49 },
      {      3,      1,     10,      7,      0,    -44,      5,   -246,      8,     -5,     27,     -7,      1,     16,      7,     42 },
      {     -2,      0,      0,     -2,     -8,      0,    -42,     -9,      0,     -8,      4,    -30,    -61,     16,   -241,    -22 },
      {      1,      8,      1,     48,     -3,     -1,     -7,     -8,     -5,     51,    -13,    243,    -12,      7,    -28,      2 },
      {      0,      0,      0,     -5,      0,      7,     -2,     44,      1,     -1,      2,     -3,     -8,     51,    -17,    246 },
    },
    {
      {   -242,     48,     31,      1,     51,     -7,     -8,     -1,     26,    -10,     -3,      1,      1,      1,      1,      0 },
      {    -60,   -155,     43,     12,   -173,     30,     37,     -2,     45,     34,    -14,     -4,     15,     -7,     -3,      1 },
      {     -4,   -187,     13,     28,    161,    -21,    -30,      2,    -25,     32,      0,     -6,    -23,      6,      6,     -1 },
      {    -31,    -10,    -63,     12,    -49,   -189,     37,     26,   -120,     30,     51,     -7,     25,     34,    -12,     -5 },
      {     22,     25,    225,    -16,     -1,    -10,     20,     -3,   -111,     -1,    -20,      2,      1,     11,     -4,     -1 },
      {     38,     12,     68,    -12,     14,   -154,      4,     27,    175,      6,    -52,      3,    -24,     18,      1,     -6 },
      {      6,     41,     -6,     10,     44,     31,    118,    -12,     18,    198,    -23,    -42,     50,    -21,    -46,      8 },
      {      2,    -19,     -9,     26,     36,      9,    212,     -4,     10,   -125,     10,     11,    -42,      1,    -22,      3 },
      {     10,    -15,      7,     41,     24,     -6,    -10,      3,     16,    -71,    -34,      9,    232,      0,    -50,      5 },
      {     -7,    -33,    -16,   -241,     17,     -2,     24,    -27,      1,     -1,     12,     53,     37,      6,    -17,      1 },
      {    -10,     -4,    -44,    -10,     -9,    -32,      1,    -33,    -59,    -16,   -212,      5,    -26,   -100,     14,     38 },
      {      5,     -2,     27,     -3,      5,    -30,     -2,     27,     15,     -5,    104,      5,     13,   -224,     16,     40 },
      {      4,      1,      9,     24,     -4,    -37,    -11,   -243,     17,     -6,     41,    -18,      0,      4,     -3,     45 },
      {      3,      4,      0,      2,      9,      1,     44,     -6,      3,     18,     -1,     19,     58,     20,    243,      6 },
      {      1,      7,      0,     54,      0,      0,      2,    -13,      0,     43,      2,    245,    -10,      1,    -22,      5 },
      {      0,      0,      2,     -4,      0,      8,     -4,     47,      3,      0,      9,     -3,     -5,     53,     -7,    245 },
    },
    {
      {   -251,     21,     33,      5,     18,     -7,      1,     -1,     27,     -1,     -6,      0,      6,     -1,      0,      0 },
      {     23,    194,    -14,    -21,    152,     40,    -28,     -9,     -4,    -32,     -5,      4,    -14,     -8,      4,      2 },
      {      0,   -150,    -31,     19,    192,    -10,    -49,      0,     26,     33,     -2,     -4,    -20,      2,      7,      0 },
      {    -19,     20,   -120,     -2,     -7,   -187,    -56,     27,    -76,    -62,     29,     16,     -3,     22,     14,     -3 },
      {     26,     -5,    203,     19,     47,   -102,     41,     26,    -73,    -44,    -26,      4,    -19,      8,     -2,     -4 },
      {     33,     34,     24,     -5,    -10,   -105,    -23,     17,    212,     19,    -60,     -3,     34,     23,    -10,     -4 },
      {      1,     45,     -6,      1,     22,    -53,     79,     15,    -32,    205,     74,    -34,     23,     58,     -4,    -16 },
      {      0,     15,     59,     -8,    -41,     10,   -210,    -34,    -49,     89,    -46,    -29,      4,     18,     37,      2 },
      {      3,     24,    -23,    230,     -3,     23,    -10,     71,    -12,     -5,    -41,    -41,     45,      6,    -17,    -13 },
      {      1,    -18,      1,    -50,     28,     15,      2,    -21,    -41,    -30,    -34,      9,    226,     58,    -57,     -8 },
      {      8,     -2,     53,     33,     -3,     25,    -55,     25,     51,    -48,    203,     45,     45,     49,     69,     -1 },
      {     -4,     -6,    -14,     -4,     -3,     34,     26,    -16,      0,    -40,    -55,    -18,    -62,    223,     55,    -41 },
      {      5,      4,      2,     65,      3,    -41,     24,   -221,     12,    -26,     34,    -84,      8,     -7,      8,     34 },
      {     -1,      0,     13,     -1,     -7,      6,    -50,     -8,      8,    -16,     60,    -18,    -60,     44,   -223,    -57 },
      {      1,      8,     -2,     58,     -2,     -6,      2,    -67,     -7,     43,    -28,    219,    -15,     32,    -34,     65 },
      {     -1,     -1,      1,    -14,      0,      7,     -6,     52,      0,     -8,     10,    -60,    -10,     49,    -40,    234 },
    },
  },
  { // 1
    {
      {   -238,     31,     30,      7,    -73,     27,     19,      1,      9,     10,      0,     -2,     17,      1,     -4,     -1 },
      {     59,    -67,     19,      2,   -212,    -15,     48,      9,    -86,     33,     28,      0,     -1,     21,      2,     -4 },
      {     50,    186,    -25,    -24,    -54,    142,    -36,    -21,    -30,     24,    -24,     -4,     -2,    -10,     -8,      3 },
      {    -34,     35,   -135,      9,    -28,    -97,   -118,     27,    -93,    -89,     -2,     36,    -30,     -9,     30,     14 },
      {      9,    -77,   -102,     -3,    -68,     40,    -95,      3,    142,     45,    -64,      2,     81,     -7,    -26,      3 },
      {    -35,    -81,    -63,     17,     66,     81,    -40,     -2,   -106,    145,     -4,    -20,    -73,     54,      5,    -11 },
      {     -8,    -65,     85,    -74,    -10,     26,    -63,    -87,    -49,    -67,   -143,    -20,    -55,    -57,    -42,     25 },
      {      7,    -47,     20,    132,     -8,    109,      3,    112,     14,    -73,    -10,     56,    -68,    -79,     29,     21 },
      {    -13,    -52,    -87,   -103,     -6,     76,     52,    -79,     31,    -65,    115,    -21,    -53,    -63,     64,      2 },
      {      4,    -17,    105,    -57,     -5,     31,   -118,      9,     22,     -6,     62,     92,     11,     97,    114,     19 },
      {     -4,    -41,     -5,      9,     45,     62,     10,    -10,   -112,    -44,     14,     13,    203,    -34,      7,     17 },
      {     -6,      6,    -26,   -105,     14,    -18,     66,     61,     -8,     65,    -31,    160,     -7,    -52,    -38,    101 },
      {     -3,    -19,    -48,    -32,      6,     58,     64,     33,      3,   -112,    -39,     28,    -13,    177,    -81,    -20 },
      {     -4,    -11,     37,    -87,      5,     11,    -71,    112,    -14,     -3,     83,    -22,     -2,    -56,   -123,   -114 },
      {     -1,      4,    -15,    -78,      4,     -5,     52,     99,     -8,      6,   -112,    -41,     16,    -12,    145,   -105 },
      {      0,      0,      9,    -45,      0,      4,    -15,    103,      1,    -10,     16,   -151,      4,     24,      7,    170 },
    },
    {
      {    191,   -130,     23,     -7,    -96,     26,     28,     -8,      4,     25,    -13,     -3,     -5,     -2,     -3,      1 },
      {   -105,    -68,     86,    -18,    -72,    162,    -62,     -4,     55,    -43,    -37,     24,     -3,     -4,     16,      3 },
      {   -104,    -90,      9,     31,     -6,    -18,    140,    -51,    -52,    116,    -67,    -26,     18,    -10,    -23,     24 },
      {      4,   -114,    129,    -39,    160,    -46,    -31,     -3,    -39,    -35,     38,      2,    -16,     32,      2,    -14 },
      {     69,     52,     -3,     17,     98,     67,    -33,    -73,    -44,     16,   -156,     67,     24,    -51,     38,     41 },
      {     -2,    -72,    -95,     97,     68,     13,     76,    -37,    117,   -116,     -4,     -8,    -33,     -6,     44,     -4 },
      {     33,     87,     99,    -78,     25,      2,     83,    -33,    140,      9,    -34,    -72,    -32,     -5,    -69,     26 },
      {     20,     -6,    -70,     18,     65,    129,    -19,     42,    -58,     21,     22,   -116,    -46,     63,   -101,     53 },
      {    -18,    -72,    -71,    -79,     39,    -41,    -57,    136,     84,     56,    -64,    -11,     28,    -83,     -3,     43 },
      {     24,     19,     30,     65,     57,     80,     20,     26,     41,     84,     86,    -22,    124,    -73,     22,   -105 },
      {      0,     -4,      2,     82,     17,    -18,    -67,     -5,     83,    126,    -18,     82,   -107,     84,    -39,    -60 },
      {      4,     24,    -21,    -93,     23,     64,     99,     80,    -37,     20,     -8,     29,   -106,     20,    115,    -85 },
      {    -13,    -35,   -103,   -140,     14,     17,    -32,   -145,     18,     19,     32,     13,     29,     10,    -42,    -83 },
      {     -9,     -5,     15,     34,    -21,    -33,    -68,    -43,    -32,     -1,    -35,   -134,   -108,   -127,     23,    -86 },
      {      1,     -1,      3,     13,     -8,    -16,     25,     68,    -27,    -83,   -115,     10,     52,     34,   -118,   -145 },
      {      3,      3,     -1,     -1,     -9,    -26,    -44,    -12,     28,     23,    -69,   -121,     77,    144,    122,    -20 },
    },
    {
      {    234,    -76,    -10,     -6,    -52,    -32,     21,      1,    -12,     20,      9,     -5,     -6,      1,     -4,     -3 },
      {     75,    123,    -49,     -2,    172,    -71,    -57,      6,    -18,    -59,     22,      9,     -8,      1,     12,     -2 },
      {     63,    118,     -1,    -20,    -44,    174,    -67,    -24,     41,      5,    -83,     13,      4,    -31,    -10,     17 },
      {     12,   -135,     78,     -5,    147,     69,      1,    -25,     95,    -39,    -54,     -1,     -6,    -26,      8,      9 },
      {    -16,    -70,    -49,     -3,    -70,    -21,   -143,     29,      5,   -138,     -3,     78,    -15,    -14,     72,     15 },
      {     11,     19,    164,    -37,     -3,     27,      6,    -47,   -156,    -67,     10,     11,    -49,     36,     35,      2 },
      {    -18,    -73,   -119,     13,     52,     41,    -30,    -24,   -151,     62,    -94,     -9,    -37,      5,    -24,     44 },
      {     -9,     24,     31,    -48,    -37,   -130,     12,    -83,     17,    -36,   -148,    -21,     28,    -79,    -50,     54 },
      {      6,     24,    -48,    109,    -12,     37,    146,     26,    -15,   -112,    -18,    -24,    -56,    -85,     50,     26 },
      {     13,     14,     77,    139,     12,    -28,    -15,     88,    -24,     60,    -66,    104,     91,      7,     26,     51 },
      {      2,     -6,    -54,    -66,     14,     30,     73,    -90,    -19,    -14,     28,     55,    163,     18,    104,     38 },
      {     -6,     17,      1,    -99,     11,    -20,     46,     36,     29,     86,      5,    122,   -124,    -54,     72,     68 },
      {      2,     18,    -10,     85,    -16,    -33,     -8,   -118,     64,     30,    -46,    -16,    -79,    136,    101,      5 },
      {     -3,     -7,     20,     67,      3,    -10,    -62,   -101,    -23,     74,     57,     -6,      0,   -163,     63,    -81 },
      {     -1,      5,     -1,    -62,      4,    -14,     13,     96,    -14,     18,   -109,    -62,     19,      4,    122,   -143 },
      {      1,      1,     26,    -12,     -2,     -2,    -53,     51,      2,     30,     38,   -158,     18,    -24,     86,    153 },
    },
  },
  { // 2
    {
      {    221,    -86,      9,     -7,    -88,     -6,     22,     -2,      3,     27,     -4,     -2,     -3,     -5,     -5,      0 },
      {     88,    103,    -56,      7,    123,   -148,     11,      8,    -61,      5,     50,     -9,      2,     16,    -10,     -7 },
      {    -20,    136,   -105,     18,   -166,    -10,     22,     14,      0,     72,     -3,    -10,     25,    -18,    -23,      6 },
      {    -75,   -110,     15,      7,    -11,    -76,    111,    -18,   -103,    120,      0,    -34,     20,     18,    -42,     12 },
      {     10,    -20,    -89,     49,     80,      6,    111,    -25,    141,     20,    -88,    -26,     13,    -75,      0,     32 },
      {    -44,    -71,    -51,     19,    -70,   -106,     44,     31,     81,   -111,    111,     -9,    -44,     33,     37,    -31 },
      {     -6,     56,    154,    -81,    -13,    -74,     13,      2,    130,     59,     30,    -20,     31,    -22,    -55,    -11 },
      {     16,     -1,    -42,     15,     50,    113,     11,     45,     43,     62,    108,    -79,    -35,    100,    -97,    -40 },
      {     11,     73,     55,     -5,    -20,     18,    110,    -94,    -17,    -42,    -51,     -8,   -154,     81,      6,      6 },
      {     18,     40,     51,     68,     -1,     61,     96,    -76,    -20,    -42,     83,      1,    151,     11,     64,    -40 },
      {     -4,    -10,     21,    152,      2,    -39,    -71,    -52,     58,     98,     -5,    108,    -28,     69,     28,    -42 },
      {      2,     18,     15,    -37,     22,     38,     52,     91,     -8,     93,     54,     28,    -76,    -66,    162,    -42 },
      {     -3,    -15,    -73,   -147,     13,      0,     21,    -21,     35,     24,    -42,     90,     68,    134,     44,    -30 },
      {     -9,    -21,    -56,    -66,      6,     12,    -34,   -164,     -7,     22,     74,     26,    -45,   -108,    -38,    -86 },
      {      2,      8,     13,     14,      6,     22,     82,     85,    -22,    -43,     -1,    168,     -9,    -45,   -121,    -58 },
      {     -1,      1,      7,     16,     -7,    -23,     -9,     33,    -12,    -18,   -104,    -86,     14,      4,      5,   -211 },
    },
    {
      {    237,    -45,     -3,     -5,    -73,    -30,     14,     -1,     -5,     26,      8,     -1,     -3,     -3,     -8,     -2 },
      {    -74,    -86,     18,     -1,   -199,     76,     31,      1,     44,     56,    -24,     -4,      7,    -12,    -12,      2 },
      {    -49,     12,    -18,      4,    -59,   -138,     40,      7,   -158,     74,     66,     -8,     18,     54,    -22,    -10 },
      {    -23,   -216,     22,     11,     87,    -57,     67,      5,     24,    -22,     28,    -13,    -16,      6,     13,     -8 },
      {      2,     39,      5,      3,     66,    -24,     74,     -8,    100,    164,    -13,    -39,     74,    -35,    -91,     -1 },
      {    -22,     41,    -26,      1,    -68,   -157,    -14,     22,    149,    -76,     56,     15,     12,    -27,     28,     -6 },
      {    -12,    -56,   -192,      5,      6,    -29,    -79,     24,    -25,     24,    -88,     17,     17,    -79,    -24,     39 },
      {     -7,     19,     99,     -6,    -17,    -76,     54,    -60,    -48,    -45,   -134,    -17,    -29,   -124,    -12,     63 },
      {     -4,     37,    -70,      4,     -1,     15,     87,     23,     17,     20,     27,    -43,   -216,    -21,    -21,     -5 },
      {     11,     33,    -75,     51,     -4,     52,    164,     29,    -28,    -57,     26,    -10,     97,    -50,     92,     18 },
      {      3,     10,      0,     67,      8,    -37,     10,     -7,     37,     83,   -102,     46,    -34,    111,    146,     61 },
      {      1,      5,    -70,    -84,    -13,     -9,     52,    -99,     22,    -74,    -75,    -90,     27,    133,    -50,     -8 },
      {     -6,     -4,    -19,   -184,     13,      6,      6,    -90,      1,     59,     58,     49,     -7,    -54,    103,     22 },
      {      2,      6,     16,   -107,      6,     -2,     63,    150,      1,    -26,    -48,    124,      3,     48,    -64,     50 },
      {      1,     -2,     30,    -65,      0,     -4,    -46,    125,     -4,     16,      0,   -182,     13,     -3,     62,     73 },
      {     -1,      6,      5,    -40,      1,    -20,     10,     57,    -11,     16,    -89,    -13,      1,    -25,     60,   -218 },
    },
    {
      {    209,    -76,      6,     -5,   -117,      9,     19,     -2,     19,     31,     -9,     -1,     -4,    -12,     -6,      1 },
      {    -97,    -77,     35,     -2,   -107,    154,    -25,     -1,     91,    -47,    -38,     11,    -16,    -12,     24,      1 },
      {     82,    -48,     47,     -7,    132,     30,    -66,      2,     38,   -149,     39,     13,    -53,     39,     40,    -18 },
      {     35,    186,    -55,     -7,    -65,    -14,    -52,     27,    121,    -55,     24,      5,    -15,    -21,      7,     -9 },
      {    -38,    -13,     50,    -25,   -105,    -44,    -58,     33,    -62,    -19,    146,    -12,    -32,    116,    -44,    -37 },
      {    -43,    -83,    -30,     10,      1,   -125,    117,      3,    143,    -19,     37,    -33,    -58,      6,     -6,     -5 },
      {     -7,    -21,   -165,     58,    -52,      3,     43,     -5,    -95,   -103,     -4,     20,    -11,     26,     95,     -8 },
      {    -23,    -25,     87,    -18,    -49,   -121,    -50,    -16,     -8,    -43,    -11,     74,     87,    -90,    116,    -15 },
      {      3,    -55,    -74,     29,     32,     36,    -33,     51,     38,     -3,    116,    -40,    169,    -59,    -43,    -11 },
      {     17,     63,     74,     65,      1,     50,    123,   -108,     24,     -9,     50,     32,     89,     81,     35,    -50 },
      {      4,     -4,    -27,    -47,     37,     38,      3,     89,     31,    129,     44,     38,    -28,     34,    148,    -84 },
      {     -7,    -35,    -47,    137,     10,    -30,   -110,    -77,     43,     67,    -10,    108,    -31,     26,    -43,    -48 },
      {      0,    -22,    -52,   -113,      5,    -31,    -16,     12,     41,    -20,   -102,     60,    102,    142,    -39,     -6 },
      {     -5,     -4,    -40,   -120,      6,     27,     43,    -70,    -27,    -23,     44,     91,    -38,    -96,    -76,   -119 },
      {     -5,    -13,    -50,    -83,      1,     -2,    -68,   -174,     20,     41,     41,    -98,     -2,     14,     69,     52 },
      {      1,      0,    -11,    -28,      9,     22,     30,      5,      8,     17,     85,    144,    -15,      0,      2,    186 },
    },
  },
  { // 3
    {
      {    168,    -71,      5,     -6,   -155,     50,      6,      2,     68,      7,    -13,      0,    -20,    -19,      7,      1 },
      {    138,     26,    -28,      1,     40,   -123,     36,     -5,   -127,     88,      6,     -3,     51,     -7,    -20,      3 },
      {     90,   -100,     40,    -10,    125,     36,    -42,      3,    -44,   -116,     38,      2,    -54,     84,      2,     -8 },
      {    -72,   -138,     63,     -4,    -50,     79,     -2,     -8,   -128,     75,    -14,      3,     76,     15,    -21,      7 },
      {     18,    -89,     17,      3,    128,     12,     31,    -17,     79,     56,    -86,      6,     14,   -145,     33,     14 },
      {    -57,   -107,    -33,     17,    -26,   -111,    121,    -15,     19,     23,    -15,    -15,   -125,     67,      4,      2 },
      {    -12,    -16,    155,    -58,     -7,    -98,    -95,     35,     59,     68,     79,    -21,    -33,    -16,    -44,      0 },
      {     21,     79,     39,    -18,     23,     86,    -27,     -6,    -45,    116,    -84,      3,   -146,     47,     47,      3 },
      {    -21,     -7,     31,    -12,    -60,    -65,    -42,      9,   -101,    -90,    -24,     31,    -55,   -112,    137,     -7 },
      {     11,     31,    119,    -17,     -7,    -50,     46,    -33,     42,    -39,   -124,     50,     93,    102,     55,    -28 },
      {     -9,    -49,   -100,    -57,     13,    -22,    -80,    113,     37,     58,     -9,    -27,     52,     70,    123,    -44 },
      {      3,    -16,      4,    182,      2,    -16,    -74,   -103,     28,     51,     59,     37,     10,     28,     74,    -33 },
      {      9,     30,     56,    -42,     18,     57,    134,    -10,      7,     24,    142,    -43,     20,     -7,    120,    -18 },
      {     -8,    -14,    -51,   -118,      3,      3,     -9,    -88,      9,     25,     45,    179,    -14,     -5,    -16,    -75 },
      {     -4,     -9,    -37,    -90,     -7,    -18,    -54,   -160,      9,     -5,     -7,    -90,     18,     24,     38,    128 },
      {     -2,     -1,     -5,    -24,     -5,      0,    -14,    -75,    -13,    -13,    -32,   -131,     -6,    -24,    -31,   -197 },
    },
    {
      {    237,    -31,    -22,     -6,     78,    -29,    -13,     -1,      5,    -17,      0,      1,    -16,     -8,      2,      0 },
      {    -71,     62,    -12,      3,    215,      9,    -35,     -4,     78,    -45,    -19,      1,      5,    -20,      1,      2 },
      {    -55,   -128,     20,     14,     50,   -157,     32,     11,    -80,    -78,     44,      4,    -47,      8,     25,     -2 },
      {     -6,   -134,     25,      8,    -27,    -30,     23,      6,    187,     18,    -18,     -5,     89,    -28,    -19,     -2 },
      {      2,    -32,     83,     -9,     86,     31,     94,    -25,    -38,    145,     41,    -31,     31,    114,    -18,    -17 },
      {    -19,    -32,   -146,     21,     12,    -56,   -117,     24,    -40,     77,    -41,     12,     84,     87,    -30,     -2 },
      {    -23,   -127,    -37,      9,     21,    113,    -49,      3,      7,     53,    -42,     -2,   -164,    -16,      9,      6 },
      {      1,    -23,     95,    -16,     14,    -12,    -23,    -36,    -76,     -7,   -167,     -6,     32,    -56,   -124,     22 },
      {     -8,     19,    -79,     55,     24,    -35,     89,     38,    -31,    109,     29,      6,      2,   -171,    -39,     15 },
      {     -7,     61,    -14,   -110,    -25,   -128,    -10,    -86,     69,     77,    -24,    -44,   -107,      9,    -18,    -10 },
      {    -14,    -60,    -68,   -142,     16,     62,      6,   -118,    -52,    -32,     51,    -44,     79,    -62,     23,      1 },
      {     -1,      2,   -102,     31,     -9,     27,    113,    -10,     20,    -89,    -17,    -74,    -32,     68,   -134,    -50 },
      {      5,      6,     16,    130,     -7,    -11,    -40,    -83,    -11,      9,    -26,   -164,     16,    -29,     81,    -65 },
      {      0,     -5,     58,     -9,     -3,      6,   -119,     40,      1,      5,    143,    -50,    -12,    -35,   -132,    -65 },
      {      2,      1,     -1,    101,     -8,     -1,    -21,   -153,     16,     -6,     67,     64,    -12,     19,    -58,    136 },
      {     -1,     -3,     -1,     32,      0,      1,      5,    -89,     -2,      6,     -8,    146,     -2,    -14,     -2,   -187 },
    },
    {
      {    196,    -70,     -2,     -6,   -136,     30,     11,      1,     41,     19,    -12,      0,    -12,    -18,      2,      1 },
      {   -148,     -5,     19,      2,   -154,     90,      3,      2,    100,     18,    -23,      0,      6,    -25,      3,      2 },
      {     35,      1,     -1,     -1,    109,     63,    -43,      1,    170,   -108,    -20,      5,    -57,    -42,     27,     -1 },
      {     45,    223,    -65,     -7,    -56,     42,    -52,      6,    -25,    -29,      2,     11,     -5,      2,     12,     -2 },
      {     25,     50,    -22,      3,     57,    -27,     27,    -11,     86,    143,    -53,    -11,    125,    -93,    -45,     11 },
      {     20,    -24,     42,    -10,     63,    205,    -46,    -17,    -77,     61,    -45,      1,     24,     40,    -22,      4 },
      {    -26,    -69,   -172,     38,      0,     24,    -40,     24,    -71,    -13,    -38,     14,    -19,   -131,     37,     19 },
      {     -5,    -23,   -146,     37,     16,      9,    -10,     26,     73,     81,     44,    -24,    -32,    157,    -23,    -26 },
      {      4,    -39,    -30,      8,    -19,      4,    -57,     20,      7,   -116,     27,     19,    204,     48,    -19,    -14 },
      {      7,     19,    -44,     21,     18,     81,    209,    -30,    -10,    -52,     79,    -36,     27,    -21,      6,     -1 },
      {      1,     -8,     29,    -25,      9,     25,    -73,     30,      8,     62,    212,    -14,      9,    -71,     50,    -21 },
      {     -6,     -7,    -19,    -51,     -8,      6,    -19,      3,     -6,    -36,     45,    -12,    -43,    -41,   -234,     21 },
      {     11,     17,     57,    237,     -5,      3,    -21,     34,     -7,     -6,     13,    -36,    -10,    -18,    -45,     -4 },
      {      4,      9,     21,    -22,     13,     15,     56,    228,     -3,      8,    -18,     89,     -4,     -2,     -9,    -28 },
      {     -1,      1,      5,    -45,     -4,     -3,    -13,     79,    -14,    -19,    -42,   -228,      8,     -9,     15,    -48 },
      {      2,      2,      4,     -4,      2,      0,     -3,     47,      9,      0,     17,    -37,      9,     28,     21,    245 },
    },
  },
  { // 4
    {
      {    246,    -37,     -5,     -7,    -56,    -12,      4,      0,     13,     12,      0,      0,     -6,    -11,      0,      0 },
      {     61,     42,    -13,      1,    232,    -40,    -19,     -2,    -19,    -58,     10,      0,    -14,      7,      8,     -1 },
      {      2,     18,      2,     -2,     35,     94,    -18,      0,    218,      1,    -40,      0,     27,    -69,     -2,      2 },
      {    -31,   -153,     13,      7,    -11,    -82,     12,      9,     53,   -129,     22,      3,    -98,    -50,     43,     -1 },
      {    -14,   -133,     22,     10,     57,    -83,     55,     -2,     43,    108,     13,    -14,    127,     47,    -43,     -7 },
      {     -3,    113,    -54,      3,    -58,   -162,    -12,     14,     55,    -65,     34,     10,    101,    -30,      6,     -5 },
      {     12,     52,    203,    -22,     -4,      7,    117,    -38,     -5,    -61,     16,    -16,     34,    -17,     -4,     -2 },
      {     10,    -36,     18,     -4,    -20,     62,    -79,      9,     16,   -121,     -5,     18,     89,    171,     27,    -14 },
      {      0,    -56,     -9,      0,      5,     35,    -40,     -6,    -92,    -66,   -103,     18,    119,   -136,    -42,     32 },
      {      3,    -18,   -106,      9,      7,     97,    122,      9,    -22,    -45,    141,    -20,     60,    -34,     39,    -20 },
      {     -3,    -18,     77,    -46,      7,      6,   -131,      5,    -17,     61,    125,     37,     38,    -68,    113,     -5 },
      {     -9,    -11,    -44,   -216,     -2,     -8,     -4,   -117,      8,    -12,      7,    -28,     -7,      9,    -41,     -1 },
      {     -1,      4,    -29,      5,      4,    -17,     55,    -76,     -4,     26,   -108,    -19,     23,     16,    204,     24 },
      {     -1,      3,     -4,    -99,     10,      0,     69,    131,      2,     10,    -34,    173,     -3,     14,     24,     41 },
      {      1,     -2,    -10,     73,     -2,      2,      1,   -149,     13,     -2,     48,    137,     -9,     14,    -34,    123 },
      {      0,     -2,     -5,     25,      1,     -1,      8,    -70,     -6,      5,    -23,    112,     -3,    -15,     -6,   -215 },
    },
    {
      {   -187,     58,     -4,      5,    147,    -26,     -6,     -2,    -60,    -13,      9,     -1,     19,     17,     -4,      1 },
      {    116,     51,    -30,      6,     49,   -129,     31,     -5,   -132,     87,      8,     -4,     57,     -2,    -19,      2 },
      {     67,   -149,     64,     -8,    105,     68,    -35,     -3,    -95,    -66,     21,      2,     -7,     63,     -4,     -2 },
      {     76,    128,    -12,    -12,     39,    -39,    -48,     16,     49,   -131,     70,     -3,    -66,     91,    -11,    -11 },
      {     51,    -39,    -33,     15,    147,    -39,     63,    -14,    110,     16,    -59,     -4,    -92,    -84,     37,      6 },
      {    -23,    -39,   -152,     58,    -46,      6,    122,    -23,    -59,    -57,    -35,     -9,    -54,     91,     15,      4 },
      {     53,     63,    -71,     32,     48,    120,     13,    -28,     11,    -70,    -34,     23,    156,    -62,     13,      2 },
      {    -12,    -72,    -54,     17,     26,    -23,     32,      7,     96,     28,    124,    -36,     76,     30,   -146,     11 },
      {      9,     25,     85,    -42,     21,      1,     77,    -22,     69,     42,   -115,     35,     63,    154,    -19,     -6 },
      {    -24,    -70,     37,    -27,    -36,   -149,     33,      1,     19,   -141,    -14,     23,     90,    -42,     50,    -10 },
      {     -3,      1,     58,    178,    -10,    -27,    -30,   -151,     15,      7,     27,     69,     -9,     11,      0,    -17 },
      {     -1,    -44,    -77,     49,     17,    -29,   -110,     62,     52,     62,      2,    -17,     58,     82,    137,    -48 },
      {      1,     18,     39,    -48,      7,     42,    124,    -27,      4,     34,    152,     23,      7,      3,    123,    -60 },
      {     -7,    -25,    -83,    -98,      3,     -8,    -49,    -42,     -2,     15,     -3,    181,    -21,     -7,    -44,    -93 },
      {     -3,    -12,    -45,   -106,     -2,    -12,    -53,   -163,     12,     12,     14,    -55,     10,     22,     47,    127 },
      {     -1,      3,     -2,    -30,     -4,     -1,    -15,    -86,     -6,    -12,    -37,   -143,      1,    -13,    -29,   -185 },
    },
    {
      {    202,    -65,      4,     -5,   -135,     22,      9,      0,     37,     17,     -8,      1,     -7,    -13,      1,     -1 },
      {    112,     92,    -49,      4,     72,   -139,     27,      0,   -108,     49,     27,     -8,     33,     19,    -24,      1 },
      {    -41,    165,    -78,      5,   -115,    -49,     20,     13,    117,      4,      1,     -4,    -13,    -24,      1,     -1 },
      {     63,    -20,    -23,      8,    135,    -40,     33,     -7,    134,    -66,     -6,     -6,   -129,      2,     22,     -3 },
      {     52,     96,     55,    -36,      0,     47,   -127,     28,    -22,   -138,     72,     15,     -3,     66,     13,    -21 },
      {      8,    -15,   -173,     68,     -3,     90,     73,    -21,    -53,   -106,      3,     -4,     31,     41,      8,      0 },
      {    -55,    -64,      1,    -11,    -90,    -90,     33,     18,    -47,    -17,     63,    -19,   -122,    134,    -31,    -10 },
      {     -8,    -49,     51,      2,      2,    -84,     57,    -24,     92,    -80,     31,      3,    169,     59,    -43,     -7 },
      {     -5,    -53,    -78,     11,     25,     17,    -75,     66,     45,     53,    162,    -65,     19,    -54,    -84,     -2 },
      {    -24,    -59,    -21,      1,    -40,   -110,    -32,     21,    -50,   -104,     19,     12,     -1,   -149,    100,    -11 },
      {     -1,    -11,    -44,   -195,     16,     17,     55,    133,     -4,    -14,    -47,    -29,     24,      8,      7,     -6 },
      {     -3,    -41,    -71,     30,     12,    -31,    -98,     45,     42,     75,    -32,     53,     45,    100,    130,    -72 },
      {     -7,    -33,    -73,    -50,     -4,    -39,    -95,    -49,      2,    -26,    -81,    124,    -14,    -13,   -144,      3 },
      {     -3,      8,     -5,    -87,      6,     34,     72,    -89,      1,     37,    124,    130,     -9,    -24,     41,    -88 },
      {      4,     13,     47,     99,      4,     12,     58,    152,    -12,    -13,    -14,     91,    -13,    -29,    -64,   -113 },
      {     -2,      3,     -4,    -19,     -5,     -3,    -20,    -73,     -9,    -11,    -49,   -127,     -3,    -17,    -35,   -198 },
    },
  },
  { // 5
    {
      {   -220,     47,     10,      5,    118,    -15,     -9,     -2,    -21,     -7,      3,      0,      4,      6,      0,      0 },
      {    114,      8,    -10,      0,    176,    -55,     -6,     -3,   -133,      1,     14,      0,     10,     16,     -3,      0 },
      {     45,    -10,      3,     -1,    121,     15,    -17,      0,    181,    -52,    -12,     -1,   -111,    -12,     15,     -1 },
      {     23,      9,     -7,      1,     55,      5,      0,     -1,    102,     33,    -14,      1,    221,    -32,    -21,      0 },
      {     28,    223,    -45,     -6,    -40,    -97,      4,      8,     31,    -24,     13,      0,     -9,     13,      0,     -1 },
      {     18,     74,     -4,     -3,      6,    187,    -57,     -1,    -46,   -134,      8,     10,     31,    -10,     18,     -1 },
      {    -12,    -73,     11,     -2,    -24,    -99,      7,      6,     14,   -182,     51,      6,     51,    105,     -5,     -8 },
      {     17,     34,    225,    -35,     -9,    -24,    -98,      0,      3,     22,    -20,     12,      6,     27,     12,     -2 },
      {      7,     23,    -17,      5,     15,     76,     20,      1,     26,     80,     13,    -11,     -2,    221,    -41,    -11 },
      {     -7,    -27,    -79,      9,    -11,    -21,   -199,     39,     12,     47,    116,     -1,      0,      0,     20,    -17 },
      {      6,     16,     68,     15,      8,     30,     82,     -3,     12,     11,    205,    -37,    -11,    -50,    -79,     -4 },
      {      7,      8,     23,    236,     -7,     -9,    -22,    -89,      1,     -3,    -10,    -17,      2,      6,     18,      8 },
      {     -2,     -5,    -17,     -6,     -7,    -10,    -64,    -46,     -5,    -24,    -55,     21,    -21,    -19,   -231,     29 },
      {     -2,     -2,    -20,    -78,     -3,     -2,     -6,   -219,      7,     14,     45,     77,      2,      8,     46,     25 },
      {     -1,     -2,     -7,    -44,     -4,     -5,    -25,    -55,     -2,    -10,    -16,   -229,      5,      6,     13,     82 },
      {      1,      0,      3,     17,      2,      2,      6,     53,      4,      6,     22,     68,      2,     13,     17,    238 },
    },
    {
      {    191,    -63,     -4,     -4,   -144,     38,      9,      1,     50,     -2,     -8,      1,    -14,     -4,      3,      0 },
      {   -125,    -16,     21,     -1,    -74,     66,    -14,      2,    171,    -50,    -11,      2,    -83,      5,     15,     -3 },
      {     12,    174,    -84,      3,    -64,   -126,     57,      5,     61,     28,     -7,     -8,    -16,     -2,     -4,      2 },
      {   -100,    -12,      1,      4,   -152,     27,     19,     -3,    -56,     60,    -11,      1,    151,    -38,    -12,      3 },
      {     -6,    -94,     28,     -3,     47,    -75,     40,     -5,     70,    168,    -72,      0,    -23,    -86,     18,      7 },
      {     47,     30,     -2,      4,     88,     30,    -22,     -3,    140,    -36,      2,      2,    176,    -34,    -13,      0 },
      {    -15,    -60,   -170,     71,     34,     61,    132,    -48,     -7,    -34,    -43,     11,     -1,      5,      6,      1 },
      {     20,    103,    -26,      3,     14,    135,    -47,      7,    -35,     39,    -30,     -1,    -49,   -157,     48,      3 },
      {    -12,    -47,   -102,     -9,      4,      4,    -48,     68,     22,     54,    168,    -69,    -18,    -38,    -86,     14 },
      {     -6,    -11,    -32,   -192,     13,     21,     65,    129,     -8,    -28,    -58,    -33,     13,      8,     20,      3 },
      {     12,     56,     -5,      1,     20,    108,    -13,      4,     18,    140,    -39,     -9,     -2,    159,    -56,      1 },
      {      8,     34,    103,    -28,     10,     39,    140,    -57,     -2,    -17,     39,    -23,    -24,    -51,   -143,     52 },
      {      3,      3,     35,    102,     -1,     -4,      9,     71,     -7,    -25,    -54,   -189,      7,     11,     32,     91 },
      {      6,     19,     56,     13,     11,     35,    109,     12,     11,     42,    133,    -20,     13,     27,    153,    -65 },
      {      3,      6,     36,    105,      3,      2,     37,    164,     -5,    -14,    -27,     66,     -7,    -21,    -73,   -118 },
      {      1,      2,      6,     31,      3,      6,     18,     81,      4,     10,     36,    137,      2,     13,     36,    190 },
    },
    {
      {   -230,     84,     -9,      6,    -64,     28,     -1,      0,     15,      0,     -1,      0,     15,     -2,     -1,      0 },
      {     54,    -32,      9,     -3,   -227,     51,     10,      1,    -72,     36,      1,      1,     27,      1,     -3,      0 },
      {     78,    187,   -119,     24,      5,     74,    -46,      6,    -28,    -24,     12,     -3,    -12,    -16,      9,     -1 },
      {    -39,     -1,     -2,      5,     43,    -49,     13,     -2,   -221,     22,     23,     -3,    -92,     32,      9,     -1 },
      {    -13,    -44,     44,    -16,     39,    176,    -70,      7,     13,     95,    -54,      4,   -107,    -17,     10,     -2 },
      {     -4,      7,     10,     -2,     70,     38,    -15,     -4,    -71,    107,    -34,      1,    200,      2,    -30,      3 },
      {     24,     89,    170,    -98,     13,     47,     60,    -43,    -33,    -84,      4,      6,     12,    -48,      0,      9 },
      {     28,     92,     27,    -24,    -15,    -58,     89,    -33,     55,    141,    -47,     -3,    -54,    109,    -48,      1 },
      {     -7,    -34,    -51,     43,     27,     80,    160,    -69,     11,     34,    120,    -60,     -3,    -54,      1,      1 },
      {      3,     21,    -12,     17,    -12,    -78,     29,     -2,     -4,     69,    -89,     24,    -26,   -206,     26,     14 },
      {     -4,    -25,    -26,     88,      9,     50,     79,    -15,    -26,    -92,   -144,     33,     -7,     21,   -123,     34 },
      {     16,     42,    116,    183,      1,      0,     11,    108,      5,     13,     26,    -33,      2,      3,     31,    -42 },
      {     -2,     -3,    -27,    -70,     11,     22,     70,    179,      2,     15,     55,    124,     -6,    -16,    -59,     -7 },
      {      2,      9,     31,      6,     -8,    -29,    -86,     -6,      4,     20,     83,    -66,    -24,    -57,   -197,     44 },
      {      4,      7,     30,     75,     -1,     -2,    -32,    -87,      9,     21,     68,    177,     -1,     11,     19,    115 },
      {      0,      1,     -3,    -19,      1,      4,     20,     79,      1,      0,    -19,    -90,      4,     17,     52,    217 },
    },
  },
  { // 6
    {
      {   -249,     36,     22,      6,    -17,      6,      5,      0,     34,     -5,     -3,     -1,     16,     -3,     -1,      0 },
      {     19,     -7,      2,     -1,   -248,     35,     21,      5,     -7,      4,      3,      0,     42,     -6,     -3,     -1 },
      {     35,     -3,     -2,     -2,     -2,      7,     -5,      0,    250,    -30,    -21,     -3,     12,     -2,     -3,      0 },
      {     22,    214,    -56,    -10,    -14,     66,    -18,     -6,     -2,    -30,      6,      2,   -100,    -14,     14,      2 },
      {     21,     85,    -28,     -6,     43,     39,    -13,     -3,    -15,      0,     -2,     -1,    228,    -28,    -16,      0 },
      {     -6,    -66,     36,     -1,     27,    225,    -49,    -11,     -1,     59,    -16,     -5,    -22,    -38,      7,      2 },
      {     29,     60,    230,    -42,      8,      2,     70,    -15,     -6,    -28,    -18,      3,      2,     -1,    -18,      3 },
      {     -7,    -50,     -1,     -2,      2,     33,    -38,      7,    -27,   -235,     39,     14,      5,    -49,     16,      5 },
      {     -3,    -16,    -57,     35,     26,     53,    222,    -31,      8,    -13,     75,    -15,     -5,    -21,    -24,      3 },
      {    -13,    -21,    -40,   -240,      0,     -1,     10,    -63,      4,      8,     29,     17,     -1,      8,      5,     16 },
      {      0,     -1,      3,     -6,     -4,    -55,      6,      0,      4,     40,    -21,      8,    -16,   -242,     34,     14 },
      {     -4,    -12,    -45,     -9,      2,     11,     52,    -29,    -24,    -44,   -228,     27,     -8,      4,    -64,     11 },
      {     -2,      0,     -9,    -57,     11,     19,     37,    236,      1,      6,     -5,     62,     -4,     -3,    -22,    -21 },
      {     -1,      1,      0,     -2,     -4,     -9,    -47,     -3,      2,      4,     52,    -38,    -21,    -35,   -238,     22 },
      {     -1,     -1,     -5,    -37,      1,     -1,      9,     50,     -8,    -17,    -36,   -239,     -1,     -6,     25,    -46 },
      {      0,      1,      1,      4,      0,      1,      3,     35,      0,      1,     -3,    -39,      3,     17,     25,    249 },
    },
    {
      {   -241,     53,     23,      4,    -53,     15,      8,     -1,     28,     -5,     -3,     -1,     19,     -4,     -2,      0 },
      {     50,    -12,      2,     -2,   -237,     47,     21,      3,    -49,     12,      6,      0,     35,     -8,     -3,      0 },
      {     36,    -24,      6,     -1,    -35,      6,     -1,      1,    241,    -35,    -24,     -1,     47,     -8,     -8,      0 },
      {     44,    213,    -79,    -14,     12,     96,    -34,     -9,     12,    -16,      5,     -1,      4,    -23,      6,      1 },
      {    -10,    -71,     60,     -8,     26,    199,    -49,    -14,     15,     80,    -31,     -4,    -78,    -17,      9,      0 },
      {      0,     23,    -17,      4,    -55,    -50,     16,      6,     39,    -40,     15,      3,   -233,     28,     23,     -1 },
      {    -37,    -82,   -206,     61,    -12,      7,    -97,     30,      4,     19,      6,      0,      3,     -1,     17,     -4 },
      {      5,     43,    -33,      4,     -5,    -47,     66,    -11,     36,    216,    -30,    -22,     -4,     78,    -24,     -7 },
      {      8,     34,     66,    -47,    -30,    -72,   -191,     42,     -8,     26,   -103,     27,      1,     35,     10,     -4 },
      {     17,     36,     62,    226,      3,      6,      2,     90,     -3,      0,    -25,    -10,      0,     -7,     -5,    -17 },
      {      1,     -4,     11,      4,      5,     53,    -22,     15,      0,    -54,     73,    -11,     25,    228,    -22,    -21 },
      {      4,     12,     49,     -9,     -6,    -28,    -68,     38,     27,     68,    199,    -34,      9,    -48,     86,    -14 },
      {      4,      4,     17,     75,    -13,    -25,    -57,   -209,     -3,    -11,     -6,   -105,      2,     10,     23,     13 },
      {     -1,      0,      4,    -14,     -4,     -9,    -45,     24,      1,      2,     52,    -84,    -26,    -41,   -224,     13 },
      {      2,      5,     11,     49,     -4,     -7,    -28,    -82,      9,     26,     62,    205,     -4,      0,    -63,     71 },
      {      0,      1,      1,     -2,      1,      4,      6,     44,      0,     -1,     -7,    -53,      4,     21,     32,    243 },
    },
    {
      {   -186,    143,    -43,      6,     79,    -47,     -3,      8,      1,     -9,     10,     -3,      5,     -1,      0,      0 },
      {   -124,    -37,     88,    -30,   -104,    130,    -59,      4,     80,    -47,     -6,     13,     -2,     -4,      7,     -3 },
      {    -54,   -114,    136,    -40,    148,    -24,    -26,     -3,    -71,     33,     -9,      7,     -3,     10,     -4,      0 },
      {    -94,    -96,    -32,     69,    -64,    -35,    143,    -72,     -6,     78,    -62,      4,     22,    -16,     -3,     10 },
      {    -11,     51,      9,    -23,    -81,     49,    -48,     19,   -185,     80,     -3,      1,     91,    -40,     -1,      3 },
      {     19,     73,    117,   -118,    -80,   -129,     41,     16,     51,     52,    -41,     12,     -6,      1,     -8,      3 },
      {     45,     45,     14,      2,     83,     62,      0,    -46,     89,     39,    -96,     43,    134,    -94,     33,      2 },
      {     28,     82,     55,     -4,     12,     73,     76,   -128,    -77,    -71,    -56,     63,    -86,     45,      2,      1 },
      {      1,    -12,     42,      0,    -15,    -47,     57,    -43,    -17,   -145,     93,    -13,    157,     11,    -36,      7 },
      {    -16,    -42,   -102,   -185,     30,     60,     91,     61,    -16,    -31,    -24,     14,     11,      5,     -8,      1 },
      {     14,     40,     23,     -6,     29,     84,     40,    -24,     44,    105,     59,    -87,     24,     98,   -130,     47 },
      {      6,     19,     42,    -13,     13,     36,     75,    -10,     -4,      8,    102,   -100,    -67,   -169,     49,     15 },
      {    -10,    -28,    -62,    -82,    -10,    -31,    -80,   -161,     22,     54,    102,     77,    -13,    -27,    -11,     -1 },
      {     -4,    -16,    -24,    -12,     -8,    -34,    -62,    -29,    -16,    -63,    -95,    -54,    -30,    -84,   -133,    129 },
      {     -2,     -6,    -23,    -47,     -5,    -14,    -45,    -87,     -4,    -13,    -60,   -172,     24,     69,    119,     -5 },
      {      1,      2,      5,      7,      1,      4,     12,     25,      3,     14,     38,     54,      9,     44,    108,    215 },
    },
  },
  { // 7
    {
      {    179,   -138,     43,     -8,    -94,     58,     -3,     -5,     13,      3,     -9,      3,     -5,     -2,      2,      0 },
      {    131,     30,    -63,     18,     78,   -134,     67,    -10,   -101,     64,     -5,     -8,     19,      0,    -10,      4 },
      {    -46,   -123,    127,    -36,    144,     -1,    -46,      7,    -83,     30,      3,      3,      2,      4,     -2,     -1 },
      {    -94,    -39,     14,     10,   -124,     -8,     64,    -27,    -83,    124,    -69,     10,     77,    -39,     -2,      7 },
      {    -28,    -66,     27,      2,     27,   -101,     92,    -30,    161,     42,    -70,     12,    -79,      7,     14,      0 },
      {    -35,    -98,   -126,    122,     65,     92,     66,    -70,    -24,    -39,      5,      1,    -12,     -4,      1,      5 },
      {     38,     20,    -21,     19,     80,     46,    -59,     11,    101,     70,    -48,      7,    140,   -112,     32,      2 },
      {    -16,    -46,     30,     -2,    -12,    -74,     56,    -18,     37,   -114,     60,    -15,    172,     45,    -57,     10 },
      {     33,     76,     70,    -24,     30,     53,     61,    -92,    -34,    -82,   -130,    101,     20,     22,     33,    -16 },
      {    -10,    -18,      9,    -50,    -15,    -35,     59,      6,    -33,    -90,     46,    -16,    -32,   -193,     90,     -7 },
      {     10,     34,     96,    164,    -31,    -59,    -76,   -112,      1,      7,     57,     13,    -19,    -34,     17,     -2 },
      {    -19,    -60,    -64,     36,    -28,    -85,    -94,     68,    -21,    -53,    -68,     75,     24,     43,    114,    -76 },
      {     -9,    -24,    -64,   -100,     -5,    -15,    -40,   -103,     16,     43,     95,    155,     -4,    -17,    -41,    -35 },
      {      5,     22,     26,    -17,     11,     45,     69,    -24,     17,     63,     99,    -38,     33,     80,    158,    -91 },
      {      6,     18,     49,     80,      9,     27,     78,    131,      5,     13,     42,    107,     -8,    -34,    -75,   -110 },
      {     -1,     -3,    -14,    -27,     -2,     -8,    -31,    -66,     -3,    -15,    -52,   -108,     -6,    -28,    -78,   -194 },
    },
    {
      {    211,    -19,    -20,     -3,   -137,     14,     14,      2,     33,     -6,     -3,     -1,      1,      1,     -1,      0 },
      {   -115,     17,     -1,      2,   -135,     -1,     14,      2,    176,    -12,    -15,     -2,    -48,      7,      5,      1 },
      {    -74,     38,     -9,      3,   -145,      0,     14,      0,   -124,     -9,     14,      1,    146,     -8,    -16,     -2 },
      {    -36,   -184,     36,      7,    -28,    157,    -23,     -8,    -33,    -36,     -1,      5,    -24,     -9,      8,      0 },
      {    -21,     47,    -10,     -5,    -71,     -9,      7,      4,   -124,     40,      8,     -1,   -198,    -28,     18,     -1 },
      {      9,    110,    -34,      0,     21,     98,      9,     -7,     -9,   -194,     23,      8,    -27,     53,    -12,     -4 },
      {     18,     47,    191,    -28,    -22,    -11,   -150,     17,      5,    -17,     45,     -6,      0,     -4,      3,     -1 },
      {      6,     74,    -64,     13,     18,    119,    -50,     -2,     39,     76,     81,    -13,     21,   -155,    -13,      8 },
      {     12,     71,     86,     -5,     12,     81,     69,    -10,     -3,     31,   -182,      8,     19,    -44,     68,      2 },
      {     -2,    -41,     31,   -116,      6,    -65,     79,     64,      2,    -99,     15,    -15,     -3,   -156,    -13,      3 },
      {      6,    -15,     48,    167,     -2,    -58,     25,   -132,      0,    -61,     10,     30,     -9,    -91,    -13,     14 },
      {     -5,    -22,    -97,     22,     -9,    -50,   -141,     28,     -9,    -62,    -82,    -22,      9,    -47,    135,      3 },
      {     -2,     -1,    -25,   -110,     -4,     -7,    -32,   -114,      4,      7,     -1,    188,      3,     -9,     19,    -52 },
      {      1,      2,     26,    -37,      3,      6,     69,    -72,      6,     10,    111,    -56,     16,     29,    187,     43 },
      {     -3,     -2,    -24,    -90,     -4,     -4,    -44,   -134,     -6,     -4,    -55,   -100,     -6,      3,    -76,    134 },
      {      0,     -2,     -4,    -28,      0,     -2,     -6,    -80,      0,     -1,    -12,   -123,     -1,     -6,    -14,   -207 },
    },
    {
      {    250,    -27,    -19,     -4,    -32,      8,      2,      1,    -28,      2,      3,      0,     -4,      0,      0,      0 },
      {    -28,     14,     -3,      0,   -246,     20,     19,      3,     20,    -18,      1,     -1,     50,     -4,     -5,      0 },
      {     31,     -2,      3,     -2,     14,    -26,      3,      0,    250,    -14,    -24,     -1,     -5,     23,     -2,      2 },
      {     -8,     25,     -3,      0,    -49,      4,      2,      1,      0,     21,     -3,      0,   -248,      1,     25,      1 },
      {    -26,   -243,     25,     17,    -17,     14,     -2,     -1,      6,     63,    -11,     -4,    -14,      2,      3,      1 },
      {     -2,      1,    -20,      7,     23,    239,    -22,    -16,     30,    -38,     20,     -2,     -3,    -62,      7,      4 },
      {     -7,    -65,    -13,      0,      7,    -32,     32,      1,    -18,   -239,     12,     18,    -30,     13,    -15,      0 },
      {     19,     18,    243,    -10,     -5,     12,    -32,     -4,     -8,    -27,    -52,      5,     -3,    -26,      0,      2 },
      {      0,    -11,    -19,      8,    -16,    -61,   -194,      7,     15,    -21,     80,     -9,      2,   -117,     49,      5 },
      {      3,     -3,     -1,     -1,      9,    -43,    132,     -7,     13,     15,     -2,      3,     -4,   -211,    -29,     14 },
      {      7,     10,     56,    148,      2,     -4,     48,    -45,     12,     10,    181,    -33,     -3,     36,    -23,      8 },
      {     -3,    -13,     32,   -196,      2,      6,     30,      9,     12,     17,    147,     48,     -2,     21,    -10,      4 },
      {     -1,     -1,    -15,    -25,     -7,    -14,    -32,   -235,     -4,      4,    -28,     39,     -4,      3,    -57,     53 },
      {      1,     -2,      7,    -17,      3,     -6,     57,    -68,     -1,    -12,      0,    -47,     22,      6,    233,     12 },
      {      1,      3,      0,     61,     -1,      2,      6,     16,      3,     14,      4,    241,      7,     -1,     57,     -9 },
      {      0,      1,      0,      7,      1,      3,      0,     56,     -3,      0,     -3,      3,      2,     13,      4,    249 },
    },
  },
  { // 8
    {
      {    206,    -65,      4,     -6,   -125,     40,      2,      2,     39,    -10,     -4,      1,     -9,      1,      2,      0 },
      {    115,    -41,     12,     -5,    153,    -13,    -12,      0,   -153,     46,      9,      0,     28,    -15,     -1,      1 },
      {    -72,   -184,     86,     -2,     24,    118,    -61,      0,     18,    -14,     12,      2,     -7,    -14,      3,      0 },
      {    -56,    -18,     30,    -10,   -144,    -21,      5,      7,   -138,     41,     -8,     -2,    133,    -41,     -6,      0 },
      {     -5,    -56,     72,    -26,     -7,   -134,     -1,     16,     60,    160,    -59,     -5,    -56,    -35,     35,     -1 },
      {     13,    -62,     53,    -16,     37,   -117,     32,     -3,     60,    -96,    -19,      9,    126,    103,    -48,     -4 },
      {     28,    114,    150,    -72,     -1,      7,   -134,     59,      8,    -29,     37,    -14,      9,      1,     14,     -4 },
      {    -20,    -21,     55,    -29,    -48,    -43,     43,      2,   -117,    -81,    -37,     25,   -157,     88,      4,    -17 },
      {     -3,     32,    -15,     -6,      2,     92,    -24,      6,     -2,     98,   -110,     15,     39,    173,     28,    -29 },
      {     10,     41,     73,    -61,     35,     65,    120,    -19,     18,    -43,   -131,     55,     30,    -98,     40,    -10 },
      {      3,    -20,    -59,     83,      0,    -52,   -124,     29,    -13,    -83,   -119,     14,     10,    -43,    110,    -22 },
      {    -14,    -43,   -109,   -176,      7,     -1,      1,    140,      2,    -15,     -3,    -17,      3,     -9,     11,    -19 },
      {     -2,     -9,     14,    -18,      2,     -4,     61,    -36,      3,     -7,     99,    -45,     27,     40,    213,    -17 },
      {     -3,     -6,    -40,    -86,     -8,    -23,    -67,   -124,      5,     12,     41,    178,      4,     -4,      8,    -49 },
      {      2,      3,     27,     95,      4,      6,     50,    141,      8,     18,     60,    120,      6,     -1,      4,   -122 },
      {      2,      1,     -2,    -18,      1,     -2,     -9,    -66,      1,     -6,    -15,   -112,     -5,    -24,    -39,   -214 },
    },
    {
      {    251,      2,    -11,     -5,    -27,     33,      1,      2,    -26,     -3,      2,      0,     -4,      3,      0,      0 },
      {     23,    -44,      5,     -2,    233,     10,    -26,     -2,    -17,     77,     -3,      1,    -42,      3,     11,     -1 },
      {    -33,    -26,      0,      2,     -9,     85,     -6,      0,   -215,    -18,     46,      1,     27,    -84,      0,      1 },
      {      0,   -155,     -4,     12,     28,     85,    -18,     -8,     65,   -102,     15,      3,    121,     36,    -49,      1 },
      {      3,    179,    -17,    -15,     91,     27,     51,     -6,     -4,   -112,      5,     15,     81,     20,    -39,      2 },
      {    -22,     50,    -35,     13,    -18,    198,    -11,    -18,     86,     37,     65,    -10,    -82,    -23,     25,      8 },
      {     12,     49,     70,      6,     -7,     -4,    -74,      9,     45,    121,     37,    -34,    160,    -86,      3,     17 },
      {     17,    -10,    147,     10,     21,    -26,     -1,      8,     48,   -105,     32,      9,    -81,   -141,    -24,     21 },
      {      9,    -23,   -183,      4,     15,    -71,     25,    -30,     39,    -11,     59,      1,      7,   -136,    -20,     23 },
      {      1,     35,    -41,     35,      2,    -32,   -221,     -9,    -16,    -78,     26,    -41,    -30,     35,     28,      2 },
      {      4,     -9,     36,    -21,      0,    -54,     38,    -34,    -11,      8,    221,      5,     -3,     90,     15,     32 },
      {      2,     16,      5,    179,     -8,     -7,      4,    -30,    -18,     43,      7,      3,    -25,     25,   -170,     10 },
      {     -6,      4,     -7,   -166,    -16,      9,    -66,     15,     -1,     33,      7,     67,    -30,      4,   -163,      5 },
      {      2,      0,     27,    -26,     -6,      0,    -10,   -245,     -6,      4,    -47,     24,      6,     -6,     16,     26 },
      {     -2,     -6,      3,    -55,      4,      1,     34,    -19,     -8,     -7,    -11,   -239,    -19,      6,    -55,     17 },
      {     -3,     -2,     -6,      3,     -2,      9,      0,     34,    -10,     -1,    -37,     14,     -1,     18,     13,    249 },
    },
    {
      {    215,    -53,     -2,     -5,    123,     -1,    -15,      0,      8,     26,     -9,      1,    -13,     13,      0,      0 },
      {   -111,    -13,      9,      0,    171,    -86,     -1,     -3,    118,     14,    -27,      3,     -5,     37,     -9,      0 },
      {     58,    -48,     17,     -1,    -96,    -78,     31,      2,    141,   -121,     -9,      5,     85,     12,    -35,      5 },
      {     25,    119,    -62,      6,    -28,    121,    -46,     -5,    134,     66,    -23,    -12,     20,     77,      0,     -6 },
      {    -14,   -138,     37,      1,    -89,    -33,    -45,     20,     13,    131,    -72,     -1,    -40,     94,     13,    -17 },
      {    -41,   -108,      9,      4,     41,    111,    -80,      5,     -2,     10,     10,    -11,    165,    -65,    -11,     -1 },
      {      0,     32,   -128,     40,      7,    -71,   -100,     29,    -85,    -49,   -110,     12,     48,     48,    -63,    -14 },
      {     -7,     16,     76,    -27,     31,     48,    108,    -16,    -73,    -19,    -15,     21,     79,    161,    -70,     -4 },
      {     29,    101,     72,    -30,    -14,    -94,     25,     12,    -12,    121,    -63,      2,    115,    -79,     16,     -3 },
      {      1,    -26,   -129,     19,     -7,    -69,     55,    -50,    -17,     64,    130,    -53,     77,     55,     66,     -4 },
      {      7,     49,    103,     -1,      7,    -33,   -124,     50,    -22,    -62,     63,    -59,     22,     80,    116,    -35 },
      {     -4,     -4,    -10,   -173,     -6,     -7,    -48,   -162,     -8,    -25,    -42,    -64,    -11,     -1,      2,    -11 },
      {     -6,    -24,    -46,    -15,     14,     42,     87,     26,     -1,    -48,   -120,     23,     18,     -4,    176,    -61 },
      {      4,      2,     39,    141,      4,     12,     33,    -75,     -2,     -8,    -58,   -172,    -12,    -20,    -25,    -57 },
      {     -4,     -8,    -36,   -106,      4,      6,     38,    150,      7,      7,     20,   -114,    -14,    -17,    -67,   -103 },
      {     -4,     -6,    -11,    -24,      3,      9,     20,     61,     -8,     -9,    -43,   -106,      2,     13,     29,    215 },
    },
  },
  { // 9
    {
      {    222,    -53,    -11,     -6,    103,      9,    -20,      0,    -11,     36,     -9,      0,    -21,     15,      2,     -1 },
      {    -89,    -31,     13,      2,    177,   -100,     -6,     -1,    106,     20,    -37,      2,    -12,     40,     -9,     -2 },
      {     63,    -85,     27,      2,    -91,   -101,     40,      8,    104,   -127,    -13,     11,     66,      5,    -38,      6 },
      {    -43,   -171,     65,      5,    -10,    -73,     -4,     17,   -141,     36,     -3,      7,    -30,    -43,     18,      0 },
      {      6,     67,     11,     -8,    103,    -32,     73,    -16,    -76,   -129,     86,     15,     14,   -112,     -7,     24 },
      {     38,     95,     29,    -12,    -53,   -127,     89,      4,    -16,     33,    -19,     12,   -149,     47,     30,      0 },
      {     -1,     32,   -154,     39,     -4,    -84,   -111,     28,    -75,    -64,    -87,      7,     -7,      1,    -37,    -14 },
      {     15,     75,     84,    -26,     15,    -32,     41,     21,    -79,     51,    -98,     20,    153,     45,    -64,    -16 },
      {    -25,    -45,     -8,      1,     29,     81,     61,    -13,    -60,    -97,      2,     35,    -54,    170,    -74,     -6 },
      {     -3,    -26,   -130,     -3,     -6,    -64,     79,    -50,    -27,     64,    107,    -37,     95,     73,     55,     -6 },
      {     10,     40,     90,     36,      8,    -30,   -111,     65,    -23,    -62,     78,    -38,     30,     96,    116,    -31 },
      {     -4,     -1,     16,   -183,     -6,    -14,    -60,   -139,    -17,    -41,    -32,    -66,     -8,     12,     25,    -16 },
      {     -7,    -23,    -36,     -9,     20,     50,     88,     37,      0,    -52,   -126,     13,     16,    -17,    172,    -54 },
      {      6,      4,     33,    141,      5,      6,     30,    -97,    -11,    -15,    -46,   -169,    -10,    -11,    -22,    -49 },
      {     -5,     -9,    -25,    -88,      6,      6,     35,    151,      6,      2,     29,   -125,    -21,    -19,    -65,   -107 },
      {      4,      3,      9,     21,     -4,    -10,    -16,    -63,     10,     10,     38,    108,     -4,    -22,    -18,   -216 },
    },
    {
      {   -209,    125,    -34,     16,     55,    -38,     10,     -3,      9,     -8,      2,     -2,      2,      1,     -1,      0 },
      {    -66,      4,     12,      0,   -205,     95,     -9,      6,     77,    -58,      7,      1,     25,     -7,     -3,     -1 },
      {    108,    150,   -152,     33,    -54,    -19,     51,    -14,    -21,    -16,     15,     -2,      6,      0,     -6,      2 },
      {    -16,     22,    -32,     16,     25,    164,   -102,     11,   -137,     47,     39,    -11,     21,    -43,     21,     -1 },
      {    -34,    -22,     44,    -18,    -97,    -93,     76,     -6,   -166,     40,    -18,     13,     78,    -43,    -12,      3 },
      {      4,    -29,      9,      4,     53,      6,     29,    -23,    -40,   -200,    113,      7,     53,    -30,    -44,     14 },
      {      6,     60,     77,    -68,    -33,      6,     26,      3,    -65,    -36,     42,     10,   -203,      6,     32,      1 },
      {    -55,   -125,   -105,    110,    -11,     26,     80,    -76,    -27,     15,      7,      3,    -96,     19,     -3,      9 },
      {      6,     -2,      4,     12,      9,    -18,     29,    -19,     63,     34,     48,    -28,    -10,   -216,     91,     24 },
      {     12,     45,     59,    -52,     50,    122,    148,   -104,     22,     35,    -38,     39,     45,     11,    -14,      8 },
      {     -2,      4,     27,     -9,    -17,    -15,     21,    -18,     28,    101,    188,   -101,     23,     73,    -28,     22 },
      {     10,     26,     36,     52,     -3,      5,    -23,      9,     12,     31,    -21,      7,    -46,    -75,   -209,     83 },
      {     30,     59,    124,    194,     -5,     -7,    -12,    -39,    -13,    -21,    -19,    -30,     11,     30,     57,    -14 },
      {     -2,     -8,      6,     42,     26,     53,    121,    207,     -1,     -2,     -7,    -46,     -2,     -2,     -2,    -22 },
      {     -3,     -4,    -11,    -44,      3,      9,      0,    -44,    -22,    -48,    -94,   -222,     -8,     -7,    -10,     20 },
      {     -2,     -3,     -5,     -2,      0,      2,      2,     30,     -7,    -11,    -14,     20,     16,     45,     74,    237 },
    },
    {
      {    168,    -33,    -18,     -4,   -166,     34,     20,      4,     75,    -18,    -10,      0,    -18,      5,      2,      0 },
      {     65,    157,    -43,    -21,    -17,   -154,     45,     20,    -48,     64,    -19,     -6,     32,     -8,      1,     -1 },
      {   -117,     54,      3,     -1,    -45,    -54,     15,      6,    179,      8,    -34,     -8,   -100,     12,     21,      2 },
      {     79,    -81,    -31,      9,    106,    -13,     19,     -7,     24,    139,    -32,    -10,   -110,    -72,     38,     11 },
      {     24,      8,    159,    -26,    -27,    -45,   -157,     35,     13,     59,     66,    -29,     -1,    -32,     -3,     12 },
      {    -98,    -73,    -14,     21,   -116,     -8,     39,     -7,    -24,    112,    -24,    -13,    100,   -105,      4,     12 },
      {    -18,    100,    -56,     31,    -36,     82,    -14,    -32,    -41,    -12,    115,      2,    -95,   -127,    -32,     24 },
      {     17,    -43,    -75,     57,     40,    -60,    -20,    -46,     94,     29,    125,      4,     70,     36,   -123,     11 },
      {    -19,    -32,    -60,   -160,      7,     18,     26,    156,      6,      0,     50,    -69,     -6,     -5,    -48,     10 },
      {     40,     79,     12,    -15,     80,     83,     -7,     11,    114,    -18,    -54,    -14,    127,   -103,     18,     17 },
      {     14,    -62,     72,     -7,     17,   -110,     71,     -2,      0,   -129,    -12,     26,    -24,   -126,    -67,     18 },
      {      3,    -20,    -49,     97,      0,    -40,    -38,     30,    -13,    -65,     -9,   -172,      4,     -1,     68,    102 },
      {      6,     41,    102,     53,      3,     64,     95,     10,    -14,     56,    -44,    -85,    -25,     49,   -141,     47 },
      {      7,      0,     78,    -23,     14,      5,    137,    -35,     21,     11,    142,    -24,     27,     26,    124,     22 },
      {      6,      3,     14,    134,      9,      3,     27,    145,      8,     -3,     34,      5,      2,    -19,     18,   -153 },
      {      1,     -3,     -2,    -56,     -1,     -8,     -6,   -114,     -1,    -14,     -8,   -146,     -4,    -22,    -16,   -165 },
    },
  },
  { // 10
    {
      {    239,     24,    -17,     -2,    -65,     51,     13,      0,    -13,    -21,     13,      2,      1,      2,     -3,      2 },
      {     39,   -106,     -3,      2,    168,     93,    -51,     -7,    -71,     69,     47,     -7,     -7,    -32,     15,     10 },
      {    -66,    -65,    -42,     12,    -73,    136,     40,    -24,    -49,   -129,     64,     26,     38,    -21,    -60,      7 },
      {     15,   -198,      0,     17,    -90,    -34,    -69,     -1,     89,     25,    -27,    -22,      0,     41,      8,     -7 },
      {    -29,     52,     54,     13,   -103,     33,    -88,    -11,      0,     91,    128,    -44,     -4,    -85,     44,     59 },
      {     16,    -29,    135,      0,    -17,    -80,    -64,     47,   -146,    -64,    -33,    -20,     64,    -24,    -38,     -7 },
      {     29,     -7,    112,     20,     70,    -35,     11,    -19,    124,    -63,     87,     52,    -32,    -27,   -104,     52 },
      {    -12,     53,     98,     23,     26,    126,    -63,     14,     66,    -45,    -47,    -61,     49,    117,     56,     -7 },
      {     -7,    -53,    113,    -45,    -28,     50,    166,     79,      6,     49,    -32,     10,    -22,    -45,     67,     19 },
      {    -14,     20,     44,     63,    -39,     69,    -44,      0,    -25,     46,    -75,     53,   -152,    -19,    -88,    -92 },
      {     -3,      8,    -40,   -159,     12,     25,    -65,    102,     47,    -60,    -26,    -80,    -58,    -75,    -36,     -4 },
      {      1,    -18,      5,     35,      7,    -33,      3,    -51,    -46,   -106,      5,    -52,   -161,     16,    111,     74 },
      {     11,      2,    -21,    131,     27,      3,     10,     -5,     57,    -41,    -74,    -67,     52,   -165,     33,    -32 },
      {      5,     -8,     54,   -100,      1,     -8,      1,   -148,     15,    -26,     26,      8,      9,    -41,     63,   -152 },
      {     -1,      7,      1,    -49,     -7,     27,    -72,    -44,     12,    -15,   -103,    162,     27,    -56,     54,    107 },
      {     -1,     -1,     19,    -43,     -9,     19,     41,   -137,    -16,     51,    -84,   -113,      9,     -1,   -103,     95 },
    },
    {
      {    243,     -9,    -15,     -2,    -16,     69,     -1,      0,    -28,      1,     22,     -1,     -1,     -8,      4,      6 },
      {     -4,    -91,     16,     -2,    195,     22,    -62,      2,    -40,    106,     12,    -14,    -15,    -13,     33,      4 },
      {    -73,    -15,    -14,     11,    -31,    158,     -7,    -21,   -108,    -45,    113,      3,     29,    -70,    -20,     27 },
      {      2,   -215,     12,     15,    -48,    -18,    -46,      2,     22,    -80,    -27,      7,     30,     29,    -67,    -15 },
      {    -10,    -28,     24,      2,   -144,    -11,   -111,     13,      2,    130,     17,    -66,    -10,    -29,     91,      8 },
      {      5,     10,     16,    -11,    -18,    -67,      4,     43,   -207,     -7,   -105,     -2,     32,    -22,     -3,    -56 },
      {     29,     13,    147,     -2,     15,   -115,    -10,     19,     -4,    -31,    100,      5,     13,   -101,    -45,     75 },
      {     -2,     57,    130,     -7,     30,     86,    -76,     56,     29,    -68,    -49,    -75,     88,     78,     26,    -13 },
      {      1,    -55,    -35,    -56,     10,    -34,    114,      2,      5,    -25,     22,    -55,    129,     -9,    143,     52 },
      {    -13,    -55,    130,    -19,    -33,     72,    136,     58,      2,     45,    -29,     56,    -94,     11,     38,    -13 },
      {      1,      2,    -12,    -46,      6,     40,     11,     45,     83,     37,    -55,      6,     79,   -169,    -48,   -105 },
      {      1,     10,     29,     58,    -22,      3,     27,    -27,     -7,    121,     18,    107,    151,     76,    -64,     18 },
      {     -9,      2,    -64,    -97,     -7,     17,    -17,    153,     -2,     30,    -27,      5,     -4,     20,    -73,    143 },
      {     -1,      0,     25,   -204,    -11,    -13,     -6,    -72,    -17,     24,     66,    -15,     -2,     65,    -48,    -76 },
      {      0,      8,      1,    -62,      1,      2,    -92,     -7,      8,    -56,    -21,    192,      5,    -22,    109,     15 },
      {      1,      0,    -38,     36,      7,    -37,     -4,    157,      0,    -14,    133,     25,      3,     50,     29,   -122 },
    },
    {
      {    220,   -113,     23,     -7,     41,     18,    -19,      4,    -22,     26,     -1,     -3,     -8,      3,      5,      0 },
      {    -60,    -55,     41,     -6,    181,   -135,      3,     13,     15,     48,    -48,      7,    -23,     25,      3,     -8 },
      {     49,    147,   -120,     22,    129,     59,    -24,    -16,    -22,     50,     32,    -16,    -20,     -1,     27,      1 },
      {     78,     63,    -23,     -2,     11,    -85,    109,    -25,    122,   -128,     17,     29,     -6,     28,    -42,     17 },
      {     17,     48,     54,    -37,     -4,    -50,    139,    -30,   -173,     20,     44,     23,     -2,    -49,     -2,     24 },
      {    -44,    -50,     72,    -38,     93,     97,     10,    -41,     47,    -59,    146,    -22,     31,    -62,     15,     43 },
      {    -23,    -90,   -133,    111,     13,    -57,    -24,     12,    -79,    -84,     52,     -4,     41,    -48,    -49,     18 },
      {    -22,    -73,    -70,     34,    -48,    -18,     85,   -105,     46,     96,     65,    -59,    -76,     55,     72,      0 },
      {    -33,    -36,     -3,      9,     33,    103,     29,    -11,    -47,    -42,    -22,     75,   -133,    103,   -109,      6 },
      {     15,     42,     50,      6,    -30,    -75,    -95,     68,    -25,    -49,     99,    -53,   -152,     13,     40,     17 },
      {      5,     10,     35,    124,    -13,     13,     37,     96,     59,     92,     31,    108,    -14,    -47,     -2,    104 },
      {     11,     33,     34,      1,    -17,    -40,    -58,    -51,     -2,     72,     61,    -61,     73,     91,   -149,     83 },
      {     -2,    -17,    -71,   -125,    -24,    -39,    -66,    -71,     26,     31,    -13,     87,    -64,   -110,    -39,     78 },
      {    -13,    -32,    -76,   -117,     -3,     -1,     23,    118,    -19,     -5,     42,     33,     56,    114,     67,     80 },
      {      4,     15,     42,     61,      4,     -7,    -62,   -119,    -39,    -62,    -49,     66,     28,     66,    115,    105 },
      {      4,      7,      8,      0,    -12,    -34,    -48,    -32,      0,     28,    111,    152,     41,     45,      7,   -145 },
    },
  },
  { // 11
    {
      {   -167,     71,      1,      6,    152,    -56,     -6,     -3,    -74,     19,     10,      1,     24,     -3,     -7,      2 },
      {    149,     83,    -69,     -2,     28,    -98,     58,     -1,    -97,     77,     -7,     -3,     40,    -27,     -5,      7 },
      {    -68,    134,    -61,      3,    -93,    -83,     68,     -1,    125,     -8,    -34,      1,    -47,     15,      9,     -1 },
      {     17,    -67,     78,    -27,     53,   -104,    -52,     30,     85,    136,    -49,    -23,    -77,    -14,     48,     -4 },
      {     78,     26,    -68,     10,    150,     44,      6,    -20,     67,    -70,     26,      2,   -119,     71,     10,     -2 },
      {    -42,   -132,    -93,     62,     23,      9,    132,    -59,      1,     34,    -97,     36,      6,     -9,     17,     -6 },
      {     22,    -47,     83,    -30,      7,   -125,     44,     20,    -20,   -116,    -55,     34,     26,    116,    -68,    -24 },
      {      7,     62,     46,   -105,     18,     89,     47,     63,    -37,    -18,   -138,     20,     -5,    -10,    106,    -40 },
      {    -17,    -35,   -136,    -12,    -14,    -13,    -95,     85,     -7,      4,    -43,   -102,     61,    102,     33,    -36 },
      {     21,     55,     40,    135,    -18,     -1,    -94,   -107,    -19,      9,    -66,     58,     19,     63,     76,    -59 },
      {    -35,    -10,     -8,      4,    -77,      4,     10,     17,   -138,     42,     18,     -2,   -178,     67,    -15,    -17 },
      {      6,      6,    -51,     26,     10,     -2,    -81,     66,      4,    -24,    -81,     80,    -41,   -109,   -128,    -89 },
      {      5,    -24,      4,     80,     -5,    -70,     -9,     32,    -41,   -129,    -17,    -84,    -49,   -113,     97,     47 },
      {     13,     45,     63,    108,     16,     66,     33,     54,     10,     37,    -88,   -101,      2,     34,   -100,     92 },
      {      5,      4,     47,     52,     12,     12,     80,     22,     16,     11,     62,    -97,     11,    -15,     -1,   -198 },
      {      2,     -6,      5,     95,      4,      3,     35,    169,     14,     21,     69,    123,     23,     34,     63,     21 },
    },
    {
      {   -217,    111,    -11,      7,    -56,    -15,     25,     -4,     25,    -33,      1,      4,     10,     -4,     -7,      1 },
      {     70,    118,    -92,      7,    -82,    152,    -29,    -18,    -33,      0,     61,    -26,     14,    -33,     24,      7 },
      {     -1,   -107,     91,    -12,   -197,     36,     22,     11,    -20,    -57,      4,     14,     17,     -6,    -27,      1 },
      {    -89,    -69,     14,     13,      2,     68,   -123,     27,    -82,    134,    -37,    -35,    -13,      1,     54,    -40 },
      {     -5,     41,     94,    -56,     28,     -1,    127,    -36,   -142,     69,     64,     21,    -25,    -41,     28,     34 },
      {     51,     72,    -56,     20,    -86,    -85,     30,     37,    -80,     51,   -138,     38,    -45,     57,    -18,    -48 },
      {    -25,    -60,   -100,    104,     17,    -56,    -14,     39,   -127,    -62,     55,     22,     28,    -83,    -54,     34 },
      {    -33,    -55,    -51,     41,     19,     73,     46,    -61,    -16,    -18,     46,     29,   -143,    137,    -51,    -10 },
      {     -4,     38,     64,     -9,     66,    111,    -28,     78,    -19,    -34,    -76,    132,     16,    -12,   -102,     29 },
      {     23,     52,     71,     80,    -36,    -49,    -62,    105,     43,     34,     81,     17,   -114,      5,     27,     93 },
      {     -9,    -46,    -53,     81,    -26,     30,     93,    -15,     84,    138,    -17,     72,     71,    -15,      2,     83 },
      {      2,     -4,     28,    105,     31,     48,     88,     88,      6,    -55,     -3,     16,     14,      4,    132,   -119 },
      {     17,     46,     83,    110,     14,      3,      8,     -8,    -24,     31,     21,   -116,     94,     84,   -110,    -20 },
      {    -14,    -19,    -46,    -74,     10,     12,     27,    110,    -43,    -32,     -3,    -45,     69,    131,     60,    118 },
      {     11,     25,     42,     57,     -4,    -31,    -83,   -126,    -42,    -43,      0,    107,     68,     76,     95,     39 },
      {      2,      4,    -26,    -67,    -18,    -42,    -26,     57,     12,     58,    138,    100,     66,     58,    -39,   -113 },
    },
    {
      {    236,    -11,    -22,     -5,    -69,     62,      8,     -1,    -16,    -15,     15,      0,     -2,     -2,      0,      3 },
      {     28,    -84,     10,      3,    175,     70,    -66,     -7,    -83,     80,     49,    -13,     -3,    -34,     22,     11 },
      {     39,    205,    -28,    -22,     89,    -50,     57,      6,    -81,     11,      0,     11,      6,    -19,      9,      0 },
      {     66,    -16,     28,    -15,     27,   -136,    -38,     31,     82,    121,    -90,    -32,    -38,     39,     68,    -17 },
      {    -30,     86,    -37,      7,    -47,     89,    -33,    -31,    106,     99,     99,    -37,    -49,    -40,     76,     41 },
      {     12,    -28,   -114,     10,    106,     45,    101,    -44,    111,    -32,    -37,     58,    -59,     61,    -15,    -19 },
      {     41,     27,    146,      8,     66,    -38,    -23,     -1,    103,    -68,     83,     34,    -38,    -14,    -81,     64 },
      {     -3,     64,     69,     -8,     38,    113,    -70,     28,     29,    -73,    -73,    -61,     47,    117,     65,    -40 },
      {     16,     10,   -121,     42,     16,    -86,   -126,    -66,     15,    -46,     65,    -48,     75,     68,    -36,     21 },
      {     -4,     45,     -2,     77,    -24,     25,   -109,    -37,    -37,     14,    -57,     49,   -142,    -11,    -84,    -93 },
      {     -8,     -9,    -72,   -153,     14,     -4,    -52,    107,     11,    -71,      2,    -74,    -91,    -52,    -29,      0 },
      {     -6,    -25,     15,     43,    -10,    -43,     24,    -46,    -77,    -82,     30,    -35,   -137,     56,    126,     75 },
      {    -16,      2,      9,    -91,    -33,      9,      2,     35,    -46,     79,     79,     75,    -29,    175,    -51,     14 },
      {      5,    -12,     35,   -114,     -1,    -32,    -12,   -122,     10,    -35,     59,     47,     11,    -28,     74,   -147 },
      {     -3,     10,    -18,    -72,     -9,     17,    -83,    -51,     -2,    -14,    -97,    139,     22,    -34,     42,    134 },
      {     -2,      3,     36,    -73,     -6,     21,     38,   -151,    -16,     40,    -64,   -127,    -16,     14,    -94,     54 },
    },
  },
};
#endif

#endif

//--------------------------------------------------------------------------------------------------

#if VCEG_AZ07_INTRA_4TAP_FILTER
Int g_aiIntraCubicFilter[32][4] = {
  {   0, 256,   0,   0 }, //  0 Integer-Pel
  {  -3, 252,   8,  -1 }, //  1
  {  -5, 247,  17,  -3 }, //  2
  {  -7, 242,  25,  -4 }, //  3
  {  -9, 236,  34,  -5 }, //  4
  { -10, 230,  43,  -7 }, //  5
  { -12, 224,  52,  -8 }, //  6
  { -13, 217,  61,  -9 }, //  7
  { -14, 210,  70, -10 }, //  8
  { -15, 203,  79, -11 }, //  9
  { -16, 195,  89, -12 }, // 10
  { -16, 187,  98, -13 }, // 11
  { -16, 179, 107, -14 }, // 12
  { -16, 170, 116, -14 }, // 13
  { -17, 162, 126, -15 }, // 14
  { -16, 153, 135, -16 }, // 15
  { -16, 144, 144, -16 }, // 16 Half-Pel
};
Int g_aiIntraGaussFilter[32][4] = {
  {  47, 161,  47,   1 }, //  0 Integer-Pel
  {  43, 161,  51,   1 }, //  1
  {  40, 160,  54,   2 }, //  2
  {  37, 159,  58,   2 }, //  3
  {  34, 158,  62,   2 }, //  4
  {  31, 156,  67,   2 }, //  5
  {  28, 154,  71,   3 }, //  6
  {  26, 151,  76,   3 }, //  7
  {  23, 149,  80,   4 }, //  8
  {  21, 146,  85,   4 }, //  9
  {  19, 142,  90,   5 }, // 10
  {  17, 139,  94,   6 }, // 11
  {  16, 135,  99,   6 }, // 12
  {  14, 131, 104,   7 }, // 13
  {  13, 127, 108,   8 }, // 14
  {  11, 123, 113,   9 }, // 15
  {  10, 118, 118,  10 }, // 16 Half-Pel
};
#if DIMD_NUM_INTRA_DIR_INC
Int g_aiExtIntraCubicFilter[64][4] = {
  { 0, 256,   0,   0 }, ///<  0 Integer-Pel
  { -1, 254,   4,  -1 }, ///< 1 
  { -3, 252,   8,  -1 }, ///< 2
  { -4, 250,  12,  -2 }, ///< 3
  { -5, 247,  17,  -3 }, ///< 4
  { -6, 244,  21,  -3 }, ///< 5
  { -7, 242,  25,  -4 }, ///< 6
  { -8, 239,  29,  -4 }, ///< 7
  { -9, 236,  34,  -5 }, ///< 8
  { -9, 233,  38,  -6 }, ///< 9
  { -10, 230,  43,  -7 }, ///< 10
  { -11, 227,  47,  -7 }, ///< 11
  { -12, 224,  52,  -8 }, ///< 12
  { -12, 220,  56,  -8 }, ///< 13
  { -13, 217,  61,  -9 }, ///< 14
  { -14, 214,  65,  -9 }, ///< 15 
  { -14, 210,  70, -10 }, ///< 16
  { -14, 206,  75, -11 }, ///< 17
  { -15, 203,  79, -11 }, ///< 18
  { -15, 199,  84, -12 }, ///< 19
  { -16, 195,  89, -12 }, ///< 20
  { -16, 191,  93, -12 }, ///< 21
  { -16, 187,  98, -13 }, ///< 22
  { -16, 183, 102, -13 }, ///< 23
  { -16, 179, 107, -14 }, ///< 24
  { -16, 174, 112, -14 }, ///< 25
  { -16, 170, 116, -14 }, ///< 26
  { -16, 166, 121, -15 }, ///< 27
  { -17, 162, 126, -15 }, ///< 28
  { -16, 157, 130, -15 }, ///< 29
  { -16, 153, 135, -16 }, ///< 30
  { -16, 148, 140, -16 }, ///< 31
  { -16, 144, 144, -16 }, ///< 32 Half-Pel
};
Int g_aiExtIntraGaussFilter[64][4] = {
  { 47, 161,  47,   1 }, ///<  0 Integer-Pel 
  { 45, 161,  49,   1 }, ///< 1
  { 43, 161,  51,   1 }, ///< 2
  { 42, 160,  52,   2 }, ///< 3
  { 40, 160,  54,   2 }, ///< 4
  { 38, 160,  56,   2 }, ///< 5
  { 37, 159,  58,   2 }, ///< 6
  { 35, 158,  61,   2 }, ///< 7
  { 34, 158,  62,   2 }, ///< 8
  { 32, 157,  65,   2 }, ///< 9
  { 31, 156,  67,   2 }, ///< 10
  { 29, 155,  69,   3 }, ///< 11
  { 28, 154,  71,   3 }, ///< 12
  { 27, 153,  73,   3 }, ///< 13
  { 26, 151,  76,   3 }, ///< 14
  { 25, 150,  78,   3 }, ///< 15
  { 23, 149,  80,   4 }, ///< 16
  { 22, 147,  83,   4 }, ///< 17
  { 21, 146,  85,   4 }, ///< 18
  { 20, 144,  87,   5 }, ///< 19
  { 19, 142,  90,   5 }, ///< 20
  { 18, 141,  92,   5 }, ///< 21
  { 17, 139,  94,   6 }, ///< 22
  { 16, 137,  97,   6 }, ///< 23
  { 16, 135,  99,   6 }, ///< 24
  { 15, 133, 101,   7 }, ///< 25
  { 14, 131, 104,   7 }, ///< 26
  { 13, 129, 106,   8 }, ///< 27
  { 13, 127, 108,   8 }, ///< 28
  { 12, 125, 111,   8 }, ///< 29
  { 11, 123, 113,   9 }, ///< 30
  { 11, 120, 116,   9 }, ///< 31
  { 10, 118, 118,  10 }, ///< 32 Half-Pel
};
#endif
#endif

//coefficients

#if RExt__HIGH_PRECISION_FORWARD_TRANSFORM
const TMatrixCoeff g_aiT4 [TRANSFORM_NUMBER_OF_DIRECTIONS][4][4]   =
{
  DEFINE_DCT4x4_MATRIX  (16384, 21266,  9224),
  DEFINE_DCT4x4_MATRIX  (   64,    83,    36)
};

const TMatrixCoeff g_aiT8 [TRANSFORM_NUMBER_OF_DIRECTIONS][8][8]   =
{
  DEFINE_DCT8x8_MATRIX  (16384, 21266,  9224, 22813, 19244, 12769,  4563),
  DEFINE_DCT8x8_MATRIX  (   64,    83,    36,    89,    75,    50,    18)
};

const TMatrixCoeff g_aiT16[TRANSFORM_NUMBER_OF_DIRECTIONS][16][16] =
{
  DEFINE_DCT16x16_MATRIX(16384, 21266,  9224, 22813, 19244, 12769,  4563, 23120, 22063, 20450, 17972, 14642, 11109,  6446,  2316),
  DEFINE_DCT16x16_MATRIX(   64,    83,    36,    89,    75,    50,    18,    90,    87,    80,    70,    57,    43,    25,     9)
};

const TMatrixCoeff g_aiT32[TRANSFORM_NUMBER_OF_DIRECTIONS][32][32] =
{
  DEFINE_DCT32x32_MATRIX(16384, 21266,  9224, 22813, 19244, 12769,  4563, 23120, 22063, 20450, 17972, 14642, 11109,  6446,  2316, 23106, 22852, 22445, 21848, 20995, 19810, 18601, 17143, 15718, 13853, 11749,  9846,  7908,  5573,  3281,   946),
  DEFINE_DCT32x32_MATRIX(   64,    83,    36,    89,    75,    50,    18,    90,    87,    80,    70,    57,    43,    25,     9,    90,    90,    88,    85,    82,    78,    73,    67,    61,    54,    46,    38,    31,    22,    13,     4)
};

const TMatrixCoeff g_as_DST_MAT_4[TRANSFORM_NUMBER_OF_DIRECTIONS][4][4] =
{
  DEFINE_DST4x4_MATRIX( 7424, 14081, 18893, 21505),
  DEFINE_DST4x4_MATRIX(   29,    55,    74,    84)
};

#else

const TMatrixCoeff g_aiT4 [TRANSFORM_NUMBER_OF_DIRECTIONS][4][4]   =
{
  DEFINE_DCT4x4_MATRIX  (   64,    83,    36),
  DEFINE_DCT4x4_MATRIX  (   64,    83,    36)
};

const TMatrixCoeff g_aiT8 [TRANSFORM_NUMBER_OF_DIRECTIONS][8][8]   =
{
  DEFINE_DCT8x8_MATRIX  (   64,    83,    36,    89,    75,    50,    18),
  DEFINE_DCT8x8_MATRIX  (   64,    83,    36,    89,    75,    50,    18)
};

const TMatrixCoeff g_aiT16[TRANSFORM_NUMBER_OF_DIRECTIONS][16][16] =
{
  DEFINE_DCT16x16_MATRIX(   64,    83,    36,    89,    75,    50,    18,    90,    87,    80,    70,    57,    43,    25,     9),
  DEFINE_DCT16x16_MATRIX(   64,    83,    36,    89,    75,    50,    18,    90,    87,    80,    70,    57,    43,    25,     9)
};

const TMatrixCoeff g_aiT32[TRANSFORM_NUMBER_OF_DIRECTIONS][32][32] =
{
  DEFINE_DCT32x32_MATRIX(   64,    83,    36,    89,    75,    50,    18,    90,    87,    80,    70,    57,    43,    25,     9,    90,    90,    88,    85,    82,    78,    73,    67,    61,    54,    46,    38,    31,    22,    13,     4),
  DEFINE_DCT32x32_MATRIX(   64,    83,    36,    89,    75,    50,    18,    90,    87,    80,    70,    57,    43,    25,     9,    90,    90,    88,    85,    82,    78,    73,    67,    61,    54,    46,    38,    31,    22,    13,     4)
};

const TMatrixCoeff g_as_DST_MAT_4[TRANSFORM_NUMBER_OF_DIRECTIONS][4][4] =
{
  DEFINE_DST4x4_MATRIX(   29,    55,    74,    84),
  DEFINE_DST4x4_MATRIX(   29,    55,    74,    84)
};
#endif


//--------------------------------------------------------------------------------------------------

#undef DEFINE_DST4x4_MATRIX
#undef DEFINE_DCT4x4_MATRIX
#undef DEFINE_DCT8x8_MATRIX
#undef DEFINE_DCT16x16_MATRIX
#undef DEFINE_DCT32x32_MATRIX

//--------------------------------------------------------------------------------------------------


const UChar g_aucChromaScale[NUM_CHROMA_FORMAT][chromaQPMappingTableSize]=
{
  //0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,29,30,31,32,33,33,34,34,35,35,36,36,37,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,51,51,51,51,51,51 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,51,51,51,51,51,51 }
};

// ====================================================================================================================
// Intra prediction
// ====================================================================================================================

#if JVET_C0024_QTBT
const UChar g_aucIntraModeNumFast_UseMPM[7-MIN_CU_LOG2+1][7-MIN_CU_LOG2+1] =
{
#if MIN_CU_LOG2==1
    {2, 2, 2, 2, 2, 2, 2},
    {3, 3, 3, 3, 2, 2, 2},  //4x2, 4x4, 4x8, 4x16, 4x32, 4x64, 4x128
    {3, 3, 3, 3, 3, 2, 2},  //8x2, 8x4, 8x8, 8x16, 8x32, 8x64, 8x128, 
    {3, 3, 3, 3, 3, 2, 2},  //16x2, 16x4, 16x8, 16x16, 16x32, 16x64, 16x128, 
    {3, 3, 3, 3, 3, 2, 2},  //32x2, 32x4, 32x8, 32x16, 32x32, 32x64, 32x128, 
    {2, 3, 3, 3, 3, 2, 2},  //64x2, 64x4, 64x8, 64x16, 64x32, 64x64, 64x128, 
    {2, 2, 2, 2, 2, 3, 2},  //128x2, 128x4, 128x8, 128x16, 128x32, 128x64, 128x128,  
#elif MIN_CU_LOG2==2
    {3, 3, 3, 3, 2, 2},  //4x4, 4x8, 4x16, 4x32, 4x64, 4x128
    {3, 3, 3, 3, 3, 2},  //8x4, 8x8, 8x16, 8x32, 8x64, 8x128, 
    {3, 3, 3, 3, 3, 2},  //16x4, 16x8, 16x16, 16x32, 16x64, 16x128, 
    {3, 3, 3, 3, 3, 2},  //32x4, 32x8, 32x16, 32x32, 32x64, 32x128, 
    {2, 3, 3, 3, 3, 2},  //64x4, 64x8, 64x16, 64x32, 64x64, 64x128, 
    {2, 2, 2, 2, 2, 3},  //128x4, 128x8, 128x16, 128x32, 128x64, 128x128, 
#endif
};
#else
const UChar g_aucIntraModeNumFast_UseMPM[MAX_CU_DEPTH] =
{
  3,  //   2x2
  8,  //   4x4
  8,  //   8x8
  3,  //  16x16
  3,  //  32x32
  3   //  64x64
#if COM16_C806_LARGE_CTU
  ,3  //  128x128   
  ,3  //  256x256
#endif
};
const UChar g_aucIntraModeNumFast_NotUseMPM[MAX_CU_DEPTH] =
{
  3,  //   2x2
  9,  //   4x4
  9,  //   8x8
  4,  //  16x16   33
  4,  //  32x32   33
  5   //  64x64   33
#if COM16_C806_LARGE_CTU
  ,9 //  128x128
  ,9 //  256x2565
#endif
};
#endif

const UChar g_chroma422IntraAngleMappingTable[NUM_INTRA_MODE] =
#if VCEG_AZ07_INTRA_65ANG_MODES
  //                                                               H                                                               D                                                               V
  //0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, DM
  { 0, 1, 2, 2, 2, 2, 2, 2, 2, 3,  4,  5,  7,  9, 11, 13, 15, 16, 18, 20, 21, 23, 25, 27, 29, 31, 32, 33, 34, 35, 37, 37, 49, 40, 40, 41, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 49, 49, 50, 51, 51, 52, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 58, 59, 60, DM_CHROMA_IDX};
#else
  //0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, DM
  { 0, 1, 2, 2, 2, 2, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 19, 20, 21, 22, 23, 23, 24, 24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 31, DM_CHROMA_IDX};
#endif

// ====================================================================================================================
// Misc.
// ====================================================================================================================

Char  g_aucConvertToBit  [ MAX_CU_SIZE+1 ];

#if ENC_DEC_TRACE
FILE*  g_hTrace = NULL; // Set to NULL to open up a file. Set to stdout to use the current output
const Bool g_bEncDecTraceEnable  = true;
const Bool g_bEncDecTraceDisable = false;
Bool   g_HLSTraceEnable = true;
Bool   g_bJustDoIt = false;
UInt64 g_nSymbolCounter = 0;
#endif

#if JVET_B0051_NON_MPM_MODE || JVET_C0038_GALF
UChar g_NonMPM[257] = { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                          4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                          5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
                          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8};
#endif
#if RRC_WIF
Int   g_iFilterCoeff[7];
#endif
// ====================================================================================================================
// Scanning order & context model mapping
// ====================================================================================================================

// scanning order table
#if COM16_C806_T64
UInt* g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES][ MAX_LOG2_TU_SIZE_PLUS_ONE ][ MAX_LOG2_TU_SIZE_PLUS_ONE ];
#else
UInt* g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES][ MAX_CU_DEPTH ][ MAX_CU_DEPTH ];
#endif

const UInt ctxIndMap4x4[4*4] =
{
  0, 1, 4, 5,
  2, 3, 4, 5,
  6, 6, 8, 8,
  7, 7, 8, 8
};

#if COM16_C806_T64
#if JVET_C0024_QTBT
const UInt g_uiMinInGroup[  LAST_SIGNIFICANT_GROUPS ] = {0,1,2,3,4,6,8,12,16,24,32,48,64,96};
const UInt g_uiGroupIdx  [  MAX_TU_SIZE ] = {0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9, 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11
,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12
,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13};
#else
const UInt g_uiMinInGroup[  LAST_SIGNIFICANT_GROUPS ] = {0,1,2,3,4,6,8,12,16,24,32,48};
const UInt g_uiGroupIdx  [  MAX_TU_SIZE ] = {0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9, 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11};
#endif
#else
const UInt g_uiMinInGroup[ LAST_SIGNIFICANT_GROUPS ] = {0,1,2,3,4,6,8,12,16,24};
const UInt g_uiGroupIdx[ MAX_TU_SIZE ]   = {0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9};
#endif

const Char *MatrixType[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM] =
{
  {
    "INTRA4X4_LUMA",
    "INTRA4X4_CHROMAU",
    "INTRA4X4_CHROMAV",
    "INTER4X4_LUMA",
    "INTER4X4_CHROMAU",
    "INTER4X4_CHROMAV"
  },
  {
    "INTRA8X8_LUMA",
    "INTRA8X8_CHROMAU",
    "INTRA8X8_CHROMAV",
    "INTER8X8_LUMA",
    "INTER8X8_CHROMAU",
    "INTER8X8_CHROMAV"
  },
  {
    "INTRA16X16_LUMA",
    "INTRA16X16_CHROMAU",
    "INTRA16X16_CHROMAV",
    "INTER16X16_LUMA",
    "INTER16X16_CHROMAU",
    "INTER16X16_CHROMAV"
  },
  {
   "INTRA32X32_LUMA",
   "INTRA32X32_CHROMAU_FROM16x16_CHROMAU",
   "INTRA32X32_CHROMAV_FROM16x16_CHROMAV",
   "INTER32X32_LUMA",
   "INTER32X32_CHROMAU_FROM16x16_CHROMAU",
   "INTER32X32_CHROMAV_FROM16x16_CHROMAV"
  },
};

const Char *MatrixType_DC[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM] =
{
  {
  },
  {
  },
  {
    "INTRA16X16_LUMA_DC",
    "INTRA16X16_CHROMAU_DC",
    "INTRA16X16_CHROMAV_DC",
    "INTER16X16_LUMA_DC",
    "INTER16X16_CHROMAU_DC",
    "INTER16X16_CHROMAV_DC"
  },
  {
    "INTRA32X32_LUMA_DC",
    "INTRA32X32_CHROMAU_DC_FROM16x16_CHROMAU",
    "INTRA32X32_CHROMAV_DC_FROM16x16_CHROMAV",
    "INTER32X32_LUMA_DC",
    "INTER32X32_CHROMAU_DC_FROM16x16_CHROMAU",
    "INTER32X32_CHROMAV_DC_FROM16x16_CHROMAV"
  },
};

const Int g_quantTSDefault4x4[4*4] =
{
  16,16,16,16,
  16,16,16,16,
  16,16,16,16,
  16,16,16,16
};

const Int g_quantIntraDefault8x8[8*8] =
{
  16,16,16,16,17,18,21,24,
  16,16,16,16,17,19,22,25,
  16,16,17,18,20,22,25,29,
  16,16,18,21,24,27,31,36,
  17,17,20,24,30,35,41,47,
  18,19,22,27,35,44,54,65,
  21,22,25,31,41,54,70,88,
  24,25,29,36,47,65,88,115
};

const Int g_quantInterDefault8x8[8*8] =
{
  16,16,16,16,17,18,20,24,
  16,16,16,17,18,20,24,25,
  16,16,17,18,20,24,25,28,
  16,17,18,20,24,25,28,33,
  17,18,20,24,25,28,33,41,
  18,20,24,25,28,33,41,54,
  20,24,25,28,33,41,54,71,
  24,25,28,33,41,54,71,91
};

#if JVET_C0024_QTBT
const UInt g_scalingListSize   [SCALING_LIST_SIZE_NUM] = {4,16,64,256,1024,4096
#if JVET_C0024_QTBT_FIX_QUANT_TICKET25
  , 16384
#endif
};
const UInt g_scalingListSizeX  [SCALING_LIST_SIZE_NUM] = {2, 4, 8, 16,  32, 64
#if JVET_C0024_QTBT_FIX_QUANT_TICKET25
  , 128
#endif
};
#else
const UInt g_scalingListSize   [SCALING_LIST_SIZE_NUM] = {16,64,256,1024
#if COM16_C806_T64
  , 4096
#endif
};
const UInt g_scalingListSizeX  [SCALING_LIST_SIZE_NUM] = { 4, 8, 16,  32
#if COM16_C806_T64
  , 64
#endif
};
#endif

#if JVET_D0033_ADAPTIVE_CLIPPING
ClipParam g_ClipParam;
Int ClipParam::nbBitsY;
Int ClipParam::nbBitsUV;
Int ClipParam::ibdLuma;
Int ClipParam::ibdChroma;
Int ClipParam::cquantiz;
#endif
#if COM16_C1046_PDPC_INTRA
#if JVET_G0104_PLANAR_PDPC
const Short g_pdpcParam[5][6] = {
		{ 33,   7,  33,   7,  30,  3 },
		{ 40,   8,  40,   8, -19,  1 },
		{ 32,   2,  32,   2, -37,  1 },
		{ 31,  -2,  31,  -2,  -4,  3 },
		{ 20,  -2,  20,  -2,   5,  5 },
};
#elif JVET_C0024_QTBT // lossless change, just remove unused entries from the table
const Int g_pdpc_pred_param[5][35][6] = {
{ {  33,   7,  33,   7,  30,    3 },
  {  25,   5,  25,   5,   0,    0 },
  {  10,   8,  29,   4,  11,    1 },
  {  10,   8,  29,   4,  11,    1 },
  {  17,   5,  20,   5,  52,    1 },
  {  17,   5,  20,   5,  52,    1 },
  {  21,   3,  18,   7,  70,    2 },
  {  21,   3,  18,   7,  70,    2 },
  {  20,   1,  18,  11,  63,    2 },
  {  20,   1,  18,  11,  63,    2 },
  {  16,   1,  30,  24,  56,    1 },
  {  16,   1,  30,  24,  56,    1 },
  {  15,   0,  15,  14,  67,    3 },
  {  15,   0,  15,  14,  67,    3 },
  {  15,   2,   9,   2,  62,    1 },
  {  15,   2,   9,   2,  62,    1 },
  {  11,   4,  10,   2,  40,    1 },
  {  11,   4,  10,   2,  40,    1 },
  {   4,   3,   4,   3,  22,    1 },
  {  10,   2,  11,   4,  40,    1 },
  {  10,   2,  11,   4,  40,    1 },
  {   9,   2,  15,   2,  62,    1 },
  {   9,   2,  15,   2,  62,    1 },
  {  15,  14,  15,   0,  67,    3 },
  {  15,  14,  15,   0,  67,    3 },
  {  30,  24,  16,   1,  56,    1 },
  {  30,  24,  16,   1,  56,    1 },
  {  18,  11,  20,   1,  63,    2 },
  {  18,  11,  20,   1,  63,    2 },
  {  18,   7,  21,   3,  70,    2 },
  {  18,   7,  21,   3,  70,    2 },
  {  20,   5,  17,   5,  52,    1 },
  {  20,   5,  17,   5,  52,    1 },
  {  29,   4,  10,   8,  11,    1 },
  {  29,   4,  10,   8,  11,    1 } },
{ {  36,   7,  36,   7,  26,    3 },
  {  33,   8,  33,   8,   0,    0 },
  {  22,   7,  32,   6,  24,    3 },
  {  22,   7,  32,   6,  24,    3 },
  {  35,   4,  29,   8,  45,    2 },
  {  35,   4,  29,   8,  45,    2 },
  {  41,   3,  27,  12,  65,    3 },
  {  41,   3,  27,  12,  65,    3 },
  {  54,   1,  26,  16,  63,    2 },
  {  54,   1,  26,  16,  63,    2 },
  {  54,  -1,  34,  25,  52,    1 },
  {  54,  -1,  34,  25,  52,    1 },
  {  24,  -1,  21,  20,  62,    1 },
  {  24,  -1,  21,  20,  62,    1 },
  {  21,   3,  19,   3,  35,    1 },
  {  21,   3,  19,   3,  35,    1 },
  {  19,   4,  21,   3,  36,    2 },
  {  19,   4,  21,   3,  36,    2 },
  {  15,   6,  15,   6,  23,    2 },
  {  21,   3,  19,   4,  36,    2 },
  {  21,   3,  19,   4,  36,    2 },
  {  19,   3,  21,   3,  35,    1 },
  {  19,   3,  21,   3,  35,    1 },
  {  21,  20,  24,  -1,  62,    1 },
  {  21,  20,  24,  -1,  62,    1 },
  {  34,  25,  54,  -1,  52,    1 },
  {  34,  25,  54,  -1,  52,    1 },
  {  26,  16,  54,   1,  63,    2 },
  {  26,  16,  54,   1,  63,    2 },
  {  27,  12,  41,   3,  65,    3 },
  {  27,  12,  41,   3,  65,    3 },
  {  29,   8,  35,   4,  45,    2 },
  {  29,   8,  35,   4,  45,    2 },
  {  32,   6,  22,   7,  24,    3 },
  {  32,   6,  22,   7,  24,    3 } },
{ {  45,   5,  45,   5,  -5,    3 },
  {  36,   8,  36,   8,   0,    0 },
  {  30,   6,  46,   6, -15,    3 },
  {  30,   6,  46,   6, -15,    3 },
  {  31,   5,  39,   8,  15,    3 },
  {  31,   5,  39,   8,  15,    3 },
  {  35,   3,  35,  11,  42,    3 },
  {  35,   3,  35,  11,  42,    3 },
  {  45,   1,  35,  19,  46,    3 },
  {  45,   1,  35,  19,  46,    3 },
  {  32,   0,  40,  32,  47,    3 },
  {  32,   0,  40,  32,  47,    3 },
  {  38,   0,  23,  13,  38,    2 },
  {  38,   0,  23,  13,  38,    2 },
  {  26,   2,  24,   0,  28,    3 },
  {  26,   2,  24,   0,  28,    3 },
  {  25,   2,  23,   0,  19,    3 },
  {  25,   2,  23,   0,  19,    3 },
  {  29,   1,  29,   1,  -7,    3 },
  {  24,   0,  25,   2,  19,    3 },
  {  24,   0,  25,   2,  19,    3 },
  {  24,   0,  26,   2,  28,    3 },
  {  24,   0,  26,   2,  28,    3 },
  {  23,  13,  38,   0,  38,    2 },
  {  23,  13,  38,   0,  38,    2 },
  {  40,  32,  32,   0,  47,    3 },
  {  40,  32,  32,   0,  47,    3 },
  {  35,  19,  45,   1,  46,    3 },
  {  35,  19,  45,   1,  46,    3 },
  {  35,  11,  35,   3,  42,    3 },
  {  35,  11,  35,   3,  42,    3 },
  {  39,   8,  31,   5,  15,    3 },
  {  39,   8,  31,   5,  15,    3 },
  {  46,   6,  30,   6, -15,    3 },
  {  46,   6,  30,   6, -15,    3 } },
{ {  46,   6,  46,   6,  -3,    5 },
  {  44,   4,  44,   4,   0,    0 },
  {  33,   3,  52,   4, -18,    5 },
  {  33,   3,  52,   4, -18,    5 },
  {  38,   3,  50,   5,  -5,    5 },
  {  38,   3,  50,   5,  -5,    5 },
  {  40,   2,  47,   9,  16,    5 },
  {  40,   2,  47,   9,  16,    5 },
  {  48,   1,  45,  17,  22,    5 },
  {  48,   1,  45,  17,  22,    5 },
  {  45,  -1,  46,  30,  36,    5 },
  {  45,  -1,  46,  30,  36,    5 },
  {  41,   1,  37,  -1,  14,    5 },
  {  41,   1,  37,  -1,  14,    5 },
  {  35,   1,  39,  -2,   3,    5 },
  {  35,   1,  39,  -2,   3,    5 },
  {  41,  -1,  43,  -1,  -7,    5 },
  {  41,  -1,  43,  -1,  -7,    5 },
  {  32,   0,  32,   0,  -6,    5 },
  {  43,  -1,  41,  -1,  -7,    5 },
  {  43,  -1,  41,  -1,  -7,    5 },
  {  39,  -2,  35,   1,   3,    5 },
  {  39,  -2,  35,   1,   3,    5 },
  {  37,  -1,  41,   1,  14,    5 },
  {  37,  -1,  41,   1,  14,    5 },
  {  46,  30,  45,  -1,  36,    5 },
  {  46,  30,  45,  -1,  36,    5 },
  {  45,  17,  48,   1,  22,    5 },
  {  45,  17,  48,   1,  22,    5 },
  {  47,   9,  40,   2,  16,    5 },
  {  47,   9,  40,   2,  16,    5 },
  {  50,   5,  38,   3,  -5,    5 },
  {  50,   5,  38,   3,  -5,    5 },
  {  52,   4,  33,   3, -18,    5 },
  {  52,   4,  33,   3, -18,    5 } },
{ {  42,   5,  42,   5,  -5,    7 },
  {  40,   3,  40,   3,   0,    0 },
  {  28,   2,  49,   3, -22,    7 },
  {  28,   2,  49,   3, -22,    7 },
  {  27,   2,  48,   3, -16,    7 },
  {  27,   2,  48,   3, -16,    7 },
  {  27,   1,  44,   5,   8,    7 },
  {  27,   1,  44,   5,   8,    7 },
  {  41,   2,  39,  12,  16,    7 },
  {  41,   2,  39,  12,  16,    7 },
  {  42,   0,  38,  21,  24,    7 },
  {  42,   0,  38,  21,  24,    7 },
  {  38,   0,  34,  -4,   5,    7 },
  {  38,   0,  34,  -4,   5,    7 },
  {  37,   1,  43,  -1,  -5,    7 },
  {  37,   1,  43,  -1,  -5,    7 },
  {  25,   0,  42,  -1, -13,    7 },
  {  25,   0,  42,  -1, -13,    7 },
  {  27,  -1,  27,  -1, -11,    7 },
  {  42,  -1,  25,   0, -13,    7 },
  {  42,  -1,  25,   0, -13,    7 },
  {  43,  -1,  37,   1,  -5,    7 },
  {  43,  -1,  37,   1,  -5,    7 },
  {  34,  -4,  38,   0,   5,    7 },
  {  34,  -4,  38,   0,   5,    7 },
  {  38,  21,  42,   0,  24,    7 },
  {  38,  21,  42,   0,  24,    7 },
  {  39,  12,  41,   2,  16,    7 },
  {  39,  12,  41,   2,  16,    7 },
  {  44,   5,  27,   1,   8,    7 },
  {  44,   5,  27,   1,   8,    7 },
  {  48,   3,  27,   2, -16,    7 },
  {  48,   3,  27,   2, -16,    7 },
  {  49,   3,  28,   2, -22,    7 },
  {  49,   3,  28,   2, -22,    7 } } };
#else
const Int g_pdpc_pred_param[5][2][35][7] =
{ { { {   27,   10,   27,   10,   29,    3,    0, },
{   22,    9,   22,    9,    0,    0,    0, },
{  -10,    7,   22,    1,   24,    1,    0, },
{  -10,    7,   22,    1,   24,    1,    0, },
{   -5,    4,   10,    1,   68,    3,    0, },
{   -5,    4,   10,    1,   68,    3,    0, },
{   -8,    3,    7,    2,   77,    2,    0, },
{   -8,    3,    7,    2,   77,    2,    0, },
{  -48,    1,    8,    6,   74,    2,    0, },
{  -48,    1,    8,    6,   74,    2,    0, },
{   20,    1,   25,   25,   63,    2,    0, },
{   20,    1,   25,   25,   63,    2,    0, },
{   14,   -1,    5,    9,   76,    3,    0, },
{   14,   -1,    5,    9,   76,    3,    0, },
{   10,    1,    1,    3,   77,    3,    0, },
{   10,    1,    1,    3,   77,    3,    0, },
{    6,    2,    2,    1,   75,    3,    0, },
{    6,    2,    2,    1,   75,    3,    0, },
{   -1,    2,   -1,    2,   37,    1,    0, },
{    2,    1,    6,    2,   75,    3,    0, },
{    2,    1,    6,    2,   75,    3,    0, },
{    1,    3,   10,    1,   77,    3,    0, },
{    1,    3,   10,    1,   77,    3,    0, },
{    5,    9,   14,   -1,   76,    3,    0, },
{    5,    9,   14,   -1,   76,    3,    0, },
{   25,   25,   20,    1,   63,    2,    0, },
{   25,   25,   20,    1,   63,    2,    0, },
{    8,    6,  -48,    1,   74,    2,    0, },
{    8,    6,  -48,    1,   74,    2,    0, },
{    7,    2,   -8,    3,   77,    2,    0, },
{    7,    2,   -8,    3,   77,    2,    0, },
{   10,    1,   -5,    4,   68,    3,    0, },
{   10,    1,   -5,    4,   68,    3,    0, },
{   22,    1,  -10,    7,   24,    1,    0, },
{   22,    1,  -10,    7,   24,    1,    0, },
},
{ {   33,    7,   33,    7,   30,    3,    0, },
{   25,    5,   25,    5,    0,    0,    0, },
{   10,    8,   29,    4,   11,    1,    0, },
{   10,    8,   29,    4,   11,    1,    0, },
{   17,    5,   20,    5,   52,    1,    0, },
{   17,    5,   20,    5,   52,    1,    0, },
{   21,    3,   18,    7,   70,    2,    0, },
{   21,    3,   18,    7,   70,    2,    0, },
{   20,    1,   18,   11,   63,    2,    0, },
{   20,    1,   18,   11,   63,    2,    0, },
{   16,    1,   30,   24,   56,    1,    0, },
{   16,    1,   30,   24,   56,    1,    0, },
{   15,    0,   15,   14,   67,    3,    0, },
{   15,    0,   15,   14,   67,    3,    0, },
{   15,    2,    9,    2,   62,    1,    0, },
{   15,    2,    9,    2,   62,    1,    0, },
{   11,    4,   10,    2,   40,    1,    0, },
{   11,    4,   10,    2,   40,    1,    0, },
{    4,    3,    4,    3,   22,    1,    0, },
{   10,    2,   11,    4,   40,    1,    0, },
{   10,    2,   11,    4,   40,    1,    0, },
{    9,    2,   15,    2,   62,    1,    0, },
{    9,    2,   15,    2,   62,    1,    0, },
{   15,   14,   15,    0,   67,    3,    0, },
{   15,   14,   15,    0,   67,    3,    0, },
{   30,   24,   16,    1,   56,    1,    0, },
{   30,   24,   16,    1,   56,    1,    0, },
{   18,   11,   20,    1,   63,    2,    0, },
{   18,   11,   20,    1,   63,    2,    0, },
{   18,    7,   21,    3,   70,    2,    0, },
{   18,    7,   21,    3,   70,    2,    0, },
{   20,    5,   17,    5,   52,    1,    0, },
{   20,    5,   17,    5,   52,    1,    0, },
{   29,    4,   10,    8,   11,    1,    0, },
{   29,    4,   10,    8,   11,    1,    0, },
},
},
{ { {   14,   10,   14,   10,   46,    3,    0, },
{   18,    9,   18,    9,    0,    0,    0, },
{  -14,    5,   12,    1,   20,    1,    0, },
{  -14,    5,   12,    1,   20,    1,    0, },
{   -5,    4,    7,    1,   73,    3,    0, },
{   -5,    4,    7,    1,   73,    3,    0, },
{   -6,    3,    7,    3,   80,    2,    0, },
{   -6,    3,    7,    3,   80,    2,    0, },
{  -48,    1,   11,    8,   74,    2,    0, },
{  -48,    1,   11,    8,   74,    2,    0, },
{    2,    1,   25,   25,   65,    3,    0, },
{    2,    1,   25,   25,   65,    3,    0, },
{    5,   -1,    5,    9,   77,    3,    0, },
{    5,   -1,    5,    9,   77,    3,    0, },
{    3,    0,    0,    2,   78,    3,    0, },
{    3,    0,    0,    2,   78,    3,    0, },
{    2,    1,    0,    0,   77,    3,    0, },
{    2,    1,    0,    0,   77,    3,    0, },
{   -2,    1,   -2,    2,   44,    1,    0, },
{    0,    0,    2,    1,   77,    3,    0, },
{    0,    0,    2,    1,   77,    3,    0, },
{    0,    2,    3,    0,   78,    3,    0, },
{    0,    2,    3,    0,   78,    3,    0, },
{    5,    9,    5,   -1,   77,    3,    0, },
{    5,    9,    5,   -1,   77,    3,    0, },
{   25,   25,    2,    1,   65,    3,    0, },
{   25,   25,    2,    1,   65,    3,    0, },
{   11,    8,  -48,    1,   74,    2,    0, },
{   11,    8,  -48,    1,   74,    2,    0, },
{    7,    3,   -6,    3,   80,    2,    0, },
{    7,    3,   -6,    3,   80,    2,    0, },
{    7,    1,   -5,    4,   73,    3,    0, },
{    7,    1,   -5,    4,   73,    3,    0, },
{   12,    1,  -14,    5,   20,    1,    0, },
{   12,    1,  -14,    5,   20,    1,    0, },
},
{ {   36,    7,   36,    7,   26,    3,    0, },
{   33,    8,   33,    8,    0,    0,    0, },
{   22,    7,   32,    6,   24,    3,    0, },
{   22,    7,   32,    6,   24,    3,    0, },
{   35,    4,   29,    8,   45,    2,    0, },
{   35,    4,   29,    8,   45,    2,    0, },
{   41,    3,   27,   12,   65,    3,    0, },
{   41,    3,   27,   12,   65,    3,    0, },
{   54,    1,   26,   16,   63,    2,    0, },
{   54,    1,   26,   16,   63,    2,    0, },
{   54,   -1,   34,   25,   52,    1,    0, },
{   54,   -1,   34,   25,   52,    1,    0, },
{   24,   -1,   21,   20,   62,    1,    0, },
{   24,   -1,   21,   20,   62,    1,    0, },
{   21,    3,   19,    3,   35,    1,    0, },
{   21,    3,   19,    3,   35,    1,    0, },
{   19,    4,   21,    3,   36,    2,    0, },
{   19,    4,   21,    3,   36,    2,    0, },
{   15,    6,   15,    6,   23,    2,    0, },
{   21,    3,   19,    4,   36,    2,    0, },
{   21,    3,   19,    4,   36,    2,    0, },
{   19,    3,   21,    3,   35,    1,    0, },
{   19,    3,   21,    3,   35,    1,    0, },
{   21,   20,   24,   -1,   62,    1,    0, },
{   21,   20,   24,   -1,   62,    1,    0, },
{   34,   25,   54,   -1,   52,    1,    0, },
{   34,   25,   54,   -1,   52,    1,    0, },
{   26,   16,   54,    1,   63,    2,    0, },
{   26,   16,   54,    1,   63,    2,    0, },
{   27,   12,   41,    3,   65,    3,    0, },
{   27,   12,   41,    3,   65,    3,    0, },
{   29,    8,   35,    4,   45,    2,    0, },
{   29,    8,   35,    4,   45,    2,    0, },
{   32,    6,   22,    7,   24,    3,    0, },
{   32,    6,   22,    7,   24,    3,    0, },
},
},
{ { {   23,   13,   23,   13,   14,    3,    0, },
{   23,   18,   23,   18,    0,    0,    0, },
{    0,    6,   35,    6,    3,    3,    0, },
{    0,    6,   35,    6,    3,    3,    0, },
{   -6,    3,   18,    4,   53,    2,    0, },
{   -6,    3,   18,    4,   53,    2,    0, },
{  -10,    2,    9,    3,   67,    3,    0, },
{  -10,    2,    9,    3,   67,    3,    0, },
{  -60,    1,   13,    9,   68,    1,    0, },
{  -60,    1,   13,    9,   68,    1,    0, },
{   55,    0,   33,   33,   56,    1,    0, },
{   55,    0,   33,   33,   56,    1,    0, },
{   34,    0,    7,   10,   57,    1,    0, },
{   34,    0,    7,   10,   57,    1,    0, },
{   14,    2,    1,    2,   55,    1,    0, },
{   14,    2,    1,    2,   55,    1,    0, },
{    6,    3,    0,    1,   46,    1,    0, },
{    6,    3,    0,    1,   46,    1,    0, },
{    7,    3,    7,    3,   18,    3,    0, },
{    0,    1,    6,    3,   46,    1,    0, },
{    0,    1,    6,    3,   46,    1,    0, },
{    1,    2,   14,    2,   55,    1,    0, },
{    1,    2,   14,    2,   55,    1,    0, },
{    7,   10,   34,    0,   57,    1,    0, },
{    7,   10,   34,    0,   57,    1,    0, },
{   33,   33,   55,    0,   56,    1,    0, },
{   33,   33,   55,    0,   56,    1,    0, },
{   13,    9,  -60,    1,   68,    1,    0, },
{   13,    9,  -60,    1,   68,    1,    0, },
{    9,    3,  -10,    2,   67,    3,    0, },
{    9,    3,  -10,    2,   67,    3,    0, },
{   18,    4,   -6,    3,   53,    2,    0, },
{   18,    4,   -6,    3,   53,    2,    0, },
{   35,    6,    0,    6,    3,    3,    0, },
{   35,    6,    0,    6,    3,    3,    0, },
},
{ {   45,    5,   45,    5,   -5,    3,    0, },
{   36,    8,   36,    8,    0,    0,    0, },
{   30,    6,   46,    6,  -15,    3,    0, },
{   30,    6,   46,    6,  -15,    3,    0, },
{   31,    5,   39,    8,   15,    3,    0, },
{   31,    5,   39,    8,   15,    3,    0, },
{   35,    3,   35,   11,   42,    3,    0, },
{   35,    3,   35,   11,   42,    3,    0, },
{   45,    1,   35,   19,   46,    3,    0, },
{   45,    1,   35,   19,   46,    3,    0, },
{   32,    0,   40,   32,   47,    3,    0, },
{   32,    0,   40,   32,   47,    3,    0, },
{   38,    0,   23,   13,   38,    2,    0, },
{   38,    0,   23,   13,   38,    2,    0, },
{   26,    2,   24,    0,   28,    3,    0, },
{   26,    2,   24,    0,   28,    3,    0, },
{   25,    2,   23,    0,   19,    3,    0, },
{   25,    2,   23,    0,   19,    3,    0, },
{   29,    1,   29,    1,   -7,    3,    0, },
{   24,    0,   25,    2,   19,    3,    0, },
{   24,    0,   25,    2,   19,    3,    0, },
{   24,    0,   26,    2,   28,    3,    0, },
{   24,    0,   26,    2,   28,    3,    0, },
{   23,   13,   38,    0,   38,    2,    0, },
{   23,   13,   38,    0,   38,    2,    0, },
{   40,   32,   32,    0,   47,    3,    0, },
{   40,   32,   32,    0,   47,    3,    0, },
{   35,   19,   45,    1,   46,    3,    0, },
{   35,   19,   45,    1,   46,    3,    0, },
{   35,   11,   35,    3,   42,    3,    0, },
{   35,   11,   35,    3,   42,    3,    0, },
{   39,    8,   31,    5,   15,    3,    0, },
{   39,    8,   31,    5,   15,    3,    0, },
{   46,    6,   30,    6,  -15,    3,    0, },
{   46,    6,   30,    6,  -15,    3,    0, },
},
},
{ { {    5,   21,    5,   21,   25,    5,    0, },
{   21,   32,   21,   32,    0,    0,    0, },
{    1,    4,   45,    7,   -4,    5,    0, },
{    1,    4,   45,    7,   -4,    5,    0, },
{    2,    2,   23,    5,   27,    5,    0, },
{    2,    2,   23,    5,   27,    5,    0, },
{   -3,    2,   21,    6,   44,    5,    0, },
{   -3,    2,   21,    6,   44,    5,    0, },
{    2,    2,   23,   15,   42,    3,    0, },
{    2,    2,   23,   15,   42,    3,    0, },
{   36,    0,   49,   49,   42,    1,    0, },
{   36,    0,   49,   49,   42,    1,    0, },
{   44,    1,   11,    6,   26,    3,    0, },
{   44,    1,   11,    6,   26,    3,    0, },
{   25,    2,   11,    2,   27,    5,    0, },
{   25,    2,   11,    2,   27,    5,    0, },
{    4,    3,    0,    4,   24,    5,    0, },
{    4,    3,    0,    4,   24,    5,    0, },
{    2,    3,    2,    3,   18,    5,    0, },
{    0,    4,    4,    3,   24,    5,    0, },
{    0,    4,    4,    3,   24,    5,    0, },
{   11,    2,   25,    2,   27,    5,    0, },
{   11,    2,   25,    2,   27,    5,    0, },
{   11,    6,   44,    1,   26,    3,    0, },
{   11,    6,   44,    1,   26,    3,    0, },
{   49,   49,   36,    0,   42,    1,    0, },
{   49,   49,   36,    0,   42,    1,    0, },
{   23,   15,    2,    2,   42,    3,    0, },
{   23,   15,    2,    2,   42,    3,    0, },
{   21,    6,   -3,    2,   44,    5,    0, },
{   21,    6,   -3,    2,   44,    5,    0, },
{   23,    5,    2,    2,   27,    5,    0, },
{   23,    5,    2,    2,   27,    5,    0, },
{   45,    7,    1,    4,   -4,    5,    0, },
{   45,    7,    1,    4,   -4,    5,    0, },
},
{ {   46,    6,   46,    6,   -3,    5,    0, },
{   44,    4,   44,    4,    0,    0,    0, },
{   33,    3,   52,    4,  -18,    5,    0, },
{   33,    3,   52,    4,  -18,    5,    0, },
{   38,    3,   50,    5,   -5,    5,    0, },
{   38,    3,   50,    5,   -5,    5,    0, },
{   40,    2,   47,    9,   16,    5,    0, },
{   40,    2,   47,    9,   16,    5,    0, },
{   48,    1,   45,   17,   22,    5,    0, },
{   48,    1,   45,   17,   22,    5,    0, },
{   45,   -1,   46,   30,   36,    5,    0, },
{   45,   -1,   46,   30,   36,    5,    0, },
{   41,    1,   37,   -1,   14,    5,    0, },
{   41,    1,   37,   -1,   14,    5,    0, },
{   35,    1,   39,   -2,    3,    5,    0, },
{   35,    1,   39,   -2,    3,    5,    0, },
{   41,   -1,   43,   -1,   -7,    5,    0, },
{   41,   -1,   43,   -1,   -7,    5,    0, },
{   32,    0,   32,    0,   -6,    5,    0, },
{   43,   -1,   41,   -1,   -7,    5,    0, },
{   43,   -1,   41,   -1,   -7,    5,    0, },
{   39,   -2,   35,    1,    3,    5,    0, },
{   39,   -2,   35,    1,    3,    5,    0, },
{   37,   -1,   41,    1,   14,    5,    0, },
{   37,   -1,   41,    1,   14,    5,    0, },
{   46,   30,   45,   -1,   36,    5,    0, },
{   46,   30,   45,   -1,   36,    5,    0, },
{   45,   17,   48,    1,   22,    5,    0, },
{   45,   17,   48,    1,   22,    5,    0, },
{   47,    9,   40,    2,   16,    5,    0, },
{   47,    9,   40,    2,   16,    5,    0, },
{   50,    5,   38,    3,   -5,    5,    0, },
{   50,    5,   38,    3,   -5,    5,    0, },
{   52,    4,   33,    3,  -18,    5,    0, },
{   52,    4,   33,    3,  -18,    5,    0, },
},
},
{ { {  -39,   32,  -39,   32,   44,    7,   10, },
{    6,   58,    6,   58,    0,    0,   10, },
{    0,    3,   55,   18,   -8,    7,   10, },
{    0,    3,   55,   18,   -8,    7,   10, },
{   -2,    2,   22,    5,   16,    7,   10, },
{   -2,    2,   22,    5,   16,    7,   10, },
{   -9,    1,   21,    5,   39,    7,   10, },
{   -9,    1,   21,    5,   39,    7,   10, },
{    3,    2,   14,    6,   18,    3,   10, },
{    3,    2,   14,    6,   18,    3,   10, },
{   10,    0,   63,   66,   49,    7,   10, },
{   10,    0,   63,   66,   49,    7,   10, },
{   38,    1,   12,    4,   25,    7,   10, },
{   38,    1,   12,    4,   25,    7,   10, },
{   23,    2,   14,    3,   26,    7,   10, },
{   23,    2,   14,    3,   26,    7,   10, },
{    6,    3,   14,    3,   11,    7,   10, },
{    6,    3,   14,    3,   11,    7,   10, },
{    3,    5,    3,    5,    4,    7,   10, },
{   14,    3,    6,    3,   11,    7,   10, },
{   14,    3,    6,    3,   11,    7,   10, },
{   14,    3,   23,    2,   26,    7,   10, },
{   14,    3,   23,    2,   26,    7,   10, },
{   12,    4,   38,    1,   25,    7,   10, },
{   12,    4,   38,    1,   25,    7,   10, },
{   63,   66,   10,    0,   49,    7,   10, },
{   63,   66,   10,    0,   49,    7,   10, },
{   14,    6,    3,    2,   18,    3,   10, },
{   14,    6,    3,    2,   18,    3,   10, },
{   21,    5,   -9,    1,   39,    7,   10, },
{   21,    5,   -9,    1,   39,    7,   10, },
{   22,    5,   -2,    2,   16,    7,   10, },
{   22,    5,   -2,    2,   16,    7,   10, },
{   55,   18,    0,    3,   -8,    7,   10, },
{   55,   18,    0,    3,   -8,    7,   10, },
},
{ {   42,    5,   42,    5,   -5,    7,   10, },
{   40,    3,   40,    3,    0,    0,   10, },
{   28,    2,   49,    3,  -22,    7,   10, },
{   28,    2,   49,    3,  -22,    7,   10, },
{   27,    2,   48,    3,  -16,    7,   10, },
{   27,    2,   48,    3,  -16,    7,   10, },
{   27,    1,   44,    5,    8,    7,   10, },
{   27,    1,   44,    5,    8,    7,   10, },
{   41,    2,   39,   12,   16,    7,   10, },
{   41,    2,   39,   12,   16,    7,   10, },
{   42,    0,   38,   21,   24,    7,   10, },
{   42,    0,   38,   21,   24,    7,   10, },
{   38,    0,   34,   -4,    5,    7,   10, },
{   38,    0,   34,   -4,    5,    7,   10, },
{   37,    1,   43,   -1,   -5,    7,   10, },
{   37,    1,   43,   -1,   -5,    7,   10, },
{   25,    0,   42,   -1,  -13,    7,   10, },
{   25,    0,   42,   -1,  -13,    7,   10, },
{   27,   -1,   27,   -1,  -11,    7,   10, },
{   42,   -1,   25,    0,  -13,    7,   10, },
{   42,   -1,   25,    0,  -13,    7,   10, },
{   43,   -1,   37,    1,   -5,    7,   10, },
{   43,   -1,   37,    1,   -5,    7,   10, },
{   34,   -4,   38,    0,    5,    7,   10, },
{   34,   -4,   38,    0,    5,    7,   10, },
{   38,   21,   42,    0,   24,    7,   10, },
{   38,   21,   42,    0,   24,    7,   10, },
{   39,   12,   41,    2,   16,    7,   10, },
{   39,   12,   41,    2,   16,    7,   10, },
{   44,    5,   27,    1,    8,    7,   10, },
{   44,    5,   27,    1,    8,    7,   10, },
{   48,    3,   27,    2,  -16,    7,   10, },
{   48,    3,   27,    2,  -16,    7,   10, },
{   49,    3,   28,    2,  -22,    7,   10, },
{   49,    3,   28,    2,  -22,    7,   10, },
},
},
};
#endif

#if RRC_DWDIP
const Int64 g_iMLRpxlNumLUT[6 * 6] = {
  0x3fbc71c71c71c71c, 0x3fb3b13b13b13b14, 0x3fa8618618618618, 0x3f9bacf914c1bad0, 0x3f8dae6076b981db, 0x3f7ecc07b301ecc0,
  0x3fb3b13b13b13b14, 0x3fae1e1e1e1e1e1e, 0x3fa47ae147ae147b, 0x3f98f9c18f9c18fa, 0x3f8c0e070381c0e0, 0x3f7de5d6e3f8868a,
  0x3fa8618618618618, 0x3fa47ae147ae147b, 0x3f9f07c1f07c1f08, 0x3f94e5e0a72f0539, 0x3f8948b0fcd6e9e0, 0x3f7c3f8f01c3f8f0,
  0x3f9bacf914c1bad0, 0x3f98f9c18f9c18fa, 0x3f94e5e0a72f0539, 0x3f8f81f81f81f820, 0x3f851d07eae2f815, 0x3f7970e4f80cb872,
  0x3f8dae6076b981db, 0x3f8c0e070381c0e0, 0x3f8948b0fcd6e9e0, 0x3f851d07eae2f815, 0x3f7fc07f01fc07f0, 0x3f75390948f40feb,
  0x3f7ecc07b301ecc0, 0x3f7de5d6e3f8868a, 0x3f7c3f8f01c3f8f0, 0x3f7970e4f80cb872, 0x3f75390948f40feb, 0x3f6fe01fe01fe020
};

const Int64 g_iMLRaverXLUT[6 * 6] = {
  0x3ff1c71c71c71c72, 0x4006276276276276, 0x4019e79e79e79e7a, 0x402c8a60dd67c8a6, 0x403e2519f89467e2, 0x404f099fc267f09a,
  0x3fe89d89d89d89d9, 0x4000f0f0f0f0f0f1, 0x4015c28f5c28f5c3, 0x4029c18f9c18f9c2, 0x403c7e3f1f8fc7e4, 0x404e21a291c07797,
  0x3fde79e79e79e79e, 0x3ff70a3d70a3d70a, 0x40107c1f07c1f07c, 0x40258d0fac687d63, 0x4039add3c0ca4588, 0x404c780e1fc780e2,
  0x3fd14c1bacf914c2, 0x3fec18f9c18f9c19, 0x4006343eb1a1f58d, 0x40203f03f03f03f0, 0x4035717c0a8e83f5, 0x4049a3c6c1fcd1e3,
  0x3fc28cfc4a33f129, 0x3fdf8fc7e3f1f8fc, 0x3ffadd3c0ca4587e, 0x4015c5f02a3a0fd6, 0x40301fc07f01fc08, 0x4045637b5b85f80b,
  0x3fb33f84cfe133f8, 0x3fd0d148e03bcbae, 0x3fee0387f1e0387f, 0x400a3c6c1fcd1e36, 0x40258ded6e17e02a, 0x40400ff00ff00ff0
};

const Int64 g_iMLRaverYLUT[6 * 6] = {
  0x3ff1c71c71c71c72, 0x3fe89d89d89d89d9, 0x3fde79e79e79e79e, 0x3fd14c1bacf914c2, 0x3fc28cfc4a33f129, 0x3fb33f84cfe133f8,
  0x4006276276276276, 0x4000f0f0f0f0f0f1, 0x3ff70a3d70a3d70a, 0x3fec18f9c18f9c19, 0x3fdf8fc7e3f1f8fc, 0x3fd0d148e03bcbae,
  0x4019e79e79e79e7a, 0x4015c28f5c28f5c3, 0x40107c1f07c1f07c, 0x4006343eb1a1f58d, 0x3ffadd3c0ca4587e, 0x3fee0387f1e0387f,
  0x402c8a60dd67c8a6, 0x4029c18f9c18f9c2, 0x40258d0fac687d63, 0x40203f03f03f03f0, 0x4015c5f02a3a0fd6, 0x400a3c6c1fcd1e36,
  0x403e2519f89467e2, 0x403c7e3f1f8fc7e4, 0x4039add3c0ca4588, 0x4035717c0a8e83f5, 0x40301fc07f01fc08, 0x40258ded6e17e02a,
  0x404f099fc267f09a, 0x404e21a291c07797, 0x404c780e1fc780e2, 0x4049a3c6c1fcd1e3, 0x4045637b5b85f80b, 0x40400ff00ff00ff0
};

const Int64 g_iMLRfactorYp2LUT[6 * 6] = {
  0x3fb4b94b94b94b93, 0x3f8d49370997fbf4, 0x3f623e20df0aedac, 0x3f34be718d6bafa4, 0x3f0641e0f63e34f6, 0x3ed718c8f97eb6ac, 
  0x3fb11dfe2ea7d6ca, 0x3f88e449fcd642fd, 0x3f60402ec0ba30ed, 0x3f334fc6fbafba95, 0x3f0560b5fe7f3cdc, 0x3ed69b13be088b13, 
  0x3fabcb564efe8978, 0x3f8430834decf8f7, 0x3f5b51cb040bf98b, 0x3f311be3bc7cd345, 0x3f03dac7e3c78cff, 0x3ed5b0e4a5059a0f, 
  0x3fa72f78e5c3bc2a, 0x3f803ea266cd8020, 0x3f55fd6a13f261b6, 0x3f2ca08974a1c527, 0x3f018cf35704c735, 0x3ed4208f0f812a8e, 
  0x3fa45b27a22edcbc, 0x3f7b0eed95916da8, 0x3f51a242829c3ffe, 0x3f26f786c3c6d9ad, 0x3efd4e252f72f28f, 0x3ed1c63c18dfcc37, 
  0x3fa2c68040155168, 0x3f77ccd7f998399a, 0x3f4d5f290bcdd00c, 0x3f22649599ee04c8, 0x3ef779e76fa9fbcb, 0x3ecda689a5f88eed
};

const Int64 g_iMLRfactorXp2LUT[6 * 6] = {
  0x3fb4b94b94b94b94, 0x3fb11dfe2ea7d6cb, 0x3fabcb564efe8984, 0x3fa72f78e5c3bc1a, 0x3fa45b27a22edcd2, 0x3fa2c6804015518c, 
  0x3f8d49370997fbf6, 0x3f88e449fcd642f9, 0x3f8430834decf8f4, 0x3f803ea266cd8011, 0x3f7b0eed95916dc8, 0x3f77ccd7f998399a, 
  0x3f623e20df0aedab, 0x3f60402ec0ba30f0, 0x3f5b51cb040bf98d, 0x3f55fd6a13f261b4, 0x3f51a242829c401b, 0x3f4d5f290bcdcfa3, 
  0x3f34be718d6bafa6, 0x3f334fc6fbafba9c, 0x3f311be3bc7cd342, 0x3f2ca08974a1c51f, 0x3f26f786c3c6d9b7, 0x3f22649599ee049f, 
  0x3f0641e0f63e34f0, 0x3f0560b5fe7f3cd9, 0x3f03dac7e3c78cf7, 0x3f018cf35704c731, 0x3efd4e252f72f2a9, 0x3ef779e76fa9fbda, 
  0x3ed718c8f97eb6ac, 0x3ed69b13be088b11, 0x3ed5b0e4a5059a1c, 0x3ed4208f0f812a94, 0x3ed1c63c18dfcc2c, 0x3ecda689a5f88ebc
};


const UInt64 g_iMLRfactorXYLUT[6 * 6] = {
  0xbfa8618618618616, 0xbf922d719c060f25, 0xbf7767dce434a9b2, 0xbf5b1c61b64ad656, 0xbf3d600b04042185, 0xbf1ea32e202c5c40, 
  0xbf922d719c060f25, 0xbf7db47fe59871e1, 0xbf64ea6e3e4fbdc1, 0xbf49f7a8e6a4250c, 0xbf2d6f3dea7fbdc4, 0xbf0f82fe20c27c7a, 
  0xbf7767dce434a9a9, 0xbf64ea6e3e4fbdc8, 0xbf505e1d27a3ee99, 0xbf366a6a886ba45d, 0xbf1b5972490550a1, 0xbefeaf7f76dfbe65, 
  0xbf5b1c61b64ad668, 0xbf49f7a8e6a42526, 0xbf366a6a886ba45c, 0xbf212b750c156ae8, 0xbf0732932308c172, 0xbeec10ebecb01426, 
  0xbf3d600b04042161, 0xbf2d6f3dea7fbda4, 0xbf1b597249055074, 0xbf0732932308c16b, 0xbef194dab88ffb7f, 0xbed798a3c1af7619, 
  0xbf1ea32e202c5c0a, 0xbf0f82fe20c27c76, 0xbefeaf7f76dfbed2, 0xbeec10ebecb0145c, 0xbed798a3c1af7604, 0xbec1ca3658877ce9
};

#endif

#endif





//! \}
