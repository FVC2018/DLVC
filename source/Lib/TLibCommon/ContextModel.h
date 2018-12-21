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


/** \file     ContextModel.h
    \brief    context model class (header)
*/

#ifndef __CONTEXTMODEL__
#define __CONTEXTMODEL__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "CommonDef.h"
#include "TComRom.h"

//! \ingroup TLibCommon
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================
#if VCEG_AZ07_BAC_ADAPT_WDOW || VCEG_AZ05_MULTI_PARAM_CABAC 
const Int m_entropyBits[2][256] =
{
  {
        92,    278,    464,    651,    838,   1027,   1216,   1406,   1596,   1788,   1980,   2173,   2367,   2561,   2756,   2953,   3150,   3347,   3546,   3746,   3946,   4147,   4349,   4552,   4756,   4960,   5166,   5372,   5580,   5788,   5997,   6207,   6418,   6630,   6843,   7057,   7272,   7488,   7705,   7923,   8141,   8361,   8582,   8804,   9027,   9251,   9476,   9702,   9930,  10158,  10388,  10618,  10850,  11083,  11317,  11552,  11788,  12026,  12265,  12505,  12746,  12988,  13232,  13477, 
     13723,  13971,  14220,  14470,  14721,  14974,  15228,  15484,  15741,  15999,  16259,  16520,  16782,  17047,  17312,  17579,  17848,  18118,  18390,  18663,  18938,  19214,  19492,  19772,  20053,  20337,  20621,  20908,  21196,  21486,  21778,  22072,  22367,  22664,  22964,  23265,  23568,  23873,  24180,  24489,  24800,  25113,  25428,  25745,  26065,  26386,  26710,  27036,  27364,  27695,  28028,  28363,  28701,  29041,  29383,  29729,  30076,  30426,  30779,  31135,  31493,  31854,  32217,  32584, 
     32953,  33325,  33700,  34079,  34460,  34844,  35232,  35622,  36016,  36414,  36814,  37218,  37626,  38037,  38452,  38870,  39292,  39718,  40148,  40581,  41019,  41461,  41907,  42357,  42812,  43271,  43734,  44202,  44675,  45153,  45635,  46122,  46615,  47112,  47615,  48124,  48638,  49157,  49682,  50214,  50751,  51294,  51844,  52400,  52963,  53532,  54109,  54693,  55283,  55882,  56488,  57102,  57724,  58354,  58993,  59641,  60297,  60963,  61639,  62324,  63019,  63725,  64441,  65168, 
     65907,  66657,  67420,  68195,  68983,  69784,  70599,  71428,  72273,  73132,  74008,  74900,  75809,  76736,  77681,  78646,  79631,  80637,  81665,  82715,  83790,  84889,  86015,  87168,  88350,  89562,  90806,  92084,  93397,  94748,  96138,  97571,  99048, 100574, 102150, 103780, 105468, 107220, 109038, 110929, 112899, 114955, 117104, 119356, 121720, 124209, 126836, 129617, 132573, 135726, 139104, 142742, 146684, 150985, 155716, 160974, 166891, 173656, 181553, 191040, 202921, 218827, 242976, 294912, 
  },
  {
    294912, 242976, 218827, 202921, 191040, 181553, 173656, 166891, 160974, 155716, 150985, 146684, 142742, 139104, 135726, 132573, 129617, 126836, 124209, 121720, 119356, 117104, 114955, 112899, 110929, 109038, 107220, 105468, 103780, 102150, 100574,  99048,  97571,  96138,  94748,  93397,  92084,  90806,  89562,  88350,  87168,  86015,  84889,  83790,  82715,  81665,  80637,  79631,  78646,  77681,  76736,  75809,  74900,  74008,  73132,  72273,  71428,  70599,  69784,  68983,  68195,  67420,  66657,  65907, 
     65168,  64441,  63725,  63019,  62324,  61639,  60963,  60297,  59641,  58993,  58354,  57724,  57102,  56488,  55882,  55283,  54693,  54109,  53532,  52963,  52400,  51844,  51294,  50751,  50214,  49682,  49157,  48638,  48124,  47615,  47112,  46615,  46122,  45635,  45153,  44675,  44202,  43734,  43271,  42812,  42357,  41907,  41461,  41019,  40581,  40148,  39718,  39292,  38870,  38452,  38037,  37626,  37218,  36814,  36414,  36016,  35622,  35232,  34844,  34460,  34079,  33700,  33325,  32953, 
     32584,  32217,  31854,  31493,  31135,  30779,  30426,  30076,  29729,  29383,  29041,  28701,  28363,  28028,  27695,  27364,  27036,  26710,  26386,  26065,  25745,  25428,  25113,  24800,  24489,  24180,  23873,  23568,  23265,  22964,  22664,  22367,  22072,  21778,  21486,  21196,  20908,  20621,  20337,  20053,  19772,  19492,  19214,  18938,  18663,  18390,  18118,  17848,  17579,  17312,  17047,  16782,  16520,  16259,  15999,  15741,  15484,  15228,  14974,  14721,  14470,  14220,  13971,  13723, 
     13477,  13232,  12988,  12746,  12505,  12265,  12026,  11788,  11552,  11317,  11083,  10850,  10618,  10388,  10158,   9930,   9702,   9476,   9251,   9027,   8804,   8582,   8361,   8141,   7923,   7705,   7488,   7272,   7057,   6843,   6630,   6418,   6207,   5997,   5788,   5580,   5372,   5166,   4960,   4756,   4552,   4349,   4147,   3946,   3746,   3546,   3347,   3150,   2953,   2756,   2561,   2367,   2173,   1980,   1788,   1596,   1406,   1216,   1027,    838,    651,    464,    278,     92, 
  }
};
#endif

/// context model class
class ContextModel
{
public:
#if VCEG_AZ07_BAC_ADAPT_WDOW  || VCEG_AZ05_MULTI_PARAM_CABAC 
  ContextModel  ()                        
  {
    m_iCtxIdx              = 0;
    m_binsCoded            = 0;  
  }
  ~ContextModel ()  {}

  Void updateLPS () //code "0"
  {
    iP1-= (iP1 >> m_ucWdow);
#if VCEG_AZ05_MULTI_PARAM_CABAC
    iP0-= (iP0 >> 8);
#endif
  }
  
  Void updateMPS () //code "1"
  {
    iP1 +=  ((32768-iP1) >> m_ucWdow);
#if VCEG_AZ05_MULTI_PARAM_CABAC
    iP0 +=  ((32768-iP0) >> 8);
#endif
  }

  UShort getState( )
  {
#if VCEG_AZ05_MULTI_PARAM_CABAC
    return  (iP0+iP1)>>1;
#else
    return  iP1;
#endif
  }  

  Int getEntropyBits(Short val) 
  {
    Int tmpIdx = 
#if VCEG_AZ05_MULTI_PARAM_CABAC
    (iP0+iP1) >> 8;
#else
    iP1 >> 7;
#endif
    return m_entropyBits[val][tmpIdx];
  }

  static Int getEntropyBitsTrm( Int val ) { return m_entropyBits[val][0]; }
  Void  init         ( Int qp, Int initValue     );
#if VCEG_AZ07_INIT_PREVFRAME
  Void  setState     ( UShort uiState )  { iP1 = uiState;
#if VCEG_AZ05_MULTI_PARAM_CABAC
    iP0 = uiState;
#endif
  }
#endif
  Void  setIdx       ( UInt i )          { m_iCtxIdx = i;        }
  Void  setWindowSize( UChar ucWdow )    { m_ucWdow   =  ucWdow ;} 
  UInt  getIdx       ()                  { return m_iCtxIdx;     }
  UChar getWdowSize  ()                  { return m_ucWdow;      }   
#else
  ContextModel  ()                        { m_ucState = 0; m_binsCoded = 0; }
  ~ContextModel ()                        {}

  UChar getState  ()                { return ( m_ucState >> 1 ); }                    ///< get current state
#if VCEG_AZ07_INIT_PREVFRAME
  UChar getOrigState ()             { return m_ucState;          } 
  Void  setState  ( UChar uiState ) { m_ucState = uiState;       }
#endif
  UChar getMps    ()                { return ( m_ucState  & 1 ); }                    ///< get curret MPS
  Void  setStateAndMps( UChar ucState, UChar ucMPS) { m_ucState = (ucState << 1) + ucMPS; } ///< set state and MPS

  Void init ( Int qp, Int initValue );   ///< initialize state with initial probability

  Void updateLPS ()
  {
    m_ucState = m_aucNextStateLPS[ m_ucState ];
  }

  Void updateMPS ()
  {
    m_ucState = m_aucNextStateMPS[ m_ucState ];
  }

  Int getEntropyBits(Short val) { return m_entropyBits[m_ucState ^ val]; }

#if FAST_BIT_EST
  Void update( Int binVal )
  {
    m_ucState = m_nextState[m_ucState][binVal];
  }
  static Void buildNextStateTable();
  static Int getEntropyBitsTrm( Int val ) { return m_entropyBits[126 ^ val]; }
#endif
#endif
  Void setBinsCoded(UInt val)   { m_binsCoded = val;  }
  UInt getBinsCoded()           { return m_binsCoded;   }

private:

#if  VCEG_AZ07_BAC_ADAPT_WDOW  || VCEG_AZ05_MULTI_PARAM_CABAC 
  UShort iP1;
  UChar         m_ucWdow;
  UInt          m_iCtxIdx;
#if VCEG_AZ05_MULTI_PARAM_CABAC
  UShort iP0;
#endif
#else
  UChar         m_ucState;                                                                  ///< internal state variable

  static const  UInt  m_totalStates = (1 << CONTEXT_STATE_BITS) * 2; //*2 for MPS = [0|1]
  static const  UChar m_aucNextStateMPS[m_totalStates];
  static const  UChar m_aucNextStateLPS[m_totalStates];
  static const  Int   m_entropyBits    [m_totalStates];
#if FAST_BIT_EST
  static UChar m_nextState[m_totalStates][2 /*MPS = [0|1]*/];
#endif
#endif
  UInt          m_binsCoded;
};

//! \}

#endif

