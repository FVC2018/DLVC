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

/** \file     TAppEncCfg.h
    \brief    Handle encoder configuration parameters (header)
*/

#ifndef __TAPPENCCFG__
#define __TAPPENCCFG__

#include "TLibCommon/CommonDef.h"

#include "TLibEncoder/TEncCfg.h"
#if EXTENSION_360_VIDEO
#include "TAppEncHelper360/TExt360AppEncCfg.h"
#endif
#include <sstream>
#include <vector>
//! \ingroup TAppEncoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// encoder configuration class
class TAppEncCfg
{
#if JVET_E0059_FLOATING_POINT_QP_FIX
public:
  template <class T>
  struct OptionalValue
  {
    Bool bPresent;
    T    value;
    OptionalValue() : bPresent(false), value() { }
  };
#endif

protected:
  // file I/O
  Char*     m_pchInputFile;                                   ///< source file name
  Char*     m_pchBitstreamFile;                               ///< output bitstream file
  Char*     m_pchReconFile;                                   ///< output reconstruction file
#if JVET_D0186_PRECISEPSNR
  Char*     m_pchPreciseLogFile;                              ///< precise log file
#endif
#if USE_CONVLUTIONAL_NEURAL_NETWORK
  Char*     m_pchDLframework;
	Char*     m_pchCaffeModePrefix;
  char*     m_pchWorkingMode;
  Int       m_iGPUid;
#endif
  Double    m_adLambdaModifier[ MAX_TLAYER ];                 ///< Lambda modifier array for each temporal layer
  // source specification
  Int       m_iFrameRate;                                     ///< source frame-rates (Hz)
  UInt      m_FrameSkip;                                      ///< number of skipped frames from the beginning
  UInt      m_temporalSubsampleRatio;                         ///< temporal subsample ratio, 2 means code every two frames
  Int       m_iSourceWidth;                                   ///< source width in pixel
  Int       m_iSourceHeight;                                  ///< source height in pixel (when interlaced = field height)
#if EXTENSION_360_VIDEO
  Int       m_inputFileWidth;                                 ///< width of image in input file  (this is equivalent to sourceWidth,  if sourceWidth  is not subsequently altered due to padding)
  Int       m_inputFileHeight;                                ///< height of image in input file (this is equivalent to sourceHeight, if sourceHeight is not subsequently altered due to padding)
  std::string m_inputFileName;
#endif
  Int       m_iSourceHeightOrg;                               ///< original source height in pixel (when interlaced = frame height)

  Bool      m_isField;                                        ///< enable field coding
  Bool      m_isTopFieldFirst;
  Bool      m_bEfficientFieldIRAPEnabled;                     ///< enable an efficient field IRAP structure.
  Bool      m_bHarmonizeGopFirstFieldCoupleEnabled;

  Int       m_conformanceWindowMode;
  Int       m_confWinLeft;
  Int       m_confWinRight;
  Int       m_confWinTop;
  Int       m_confWinBottom;
  Int       m_framesToBeEncoded;                              ///< number of encoded frames
  Int       m_aiPad[2];                                       ///< number of padded pixels for width and height
  InputColourSpaceConversion m_inputColourSpaceConvert;       ///< colour space conversion to apply to input video
  Bool      m_snrInternalColourSpace;                       ///< if true, then no colour space conversion is applied for snr calculation, otherwise inverse of input is applied.
  Bool      m_outputInternalColourSpace;                    ///< if true, then no colour space conversion is applied for reconstructed video, otherwise inverse of input is applied.
  ChromaFormat m_InputChromaFormatIDC;

  Bool      m_printMSEBasedSequencePSNR;
  Bool      m_printFrameMSE;
  Bool      m_printSequenceMSE;
#if JVET_D0134_PSNR
  Bool      m_trueBidepthPSNR;
#endif
#if JVET_F0064_MSSSIM
  Bool      m_printMSSSIM;
#endif
  Bool      m_cabacZeroWordPaddingEnabled;
  Bool      m_bClipInputVideoToRec709Range;
  Bool      m_bClipOutputVideoToRec709Range;

  // profile/level
  Profile::Name m_profile;
  Level::Tier   m_levelTier;
  Level::Name   m_level;
  UInt          m_bitDepthConstraint;
  ChromaFormat  m_chromaFormatConstraint;
  Bool          m_intraConstraintFlag;
  Bool          m_onePictureOnlyConstraintFlag;
  Bool          m_lowerBitRateConstraintFlag;
  Bool m_progressiveSourceFlag;
  Bool m_interlacedSourceFlag;
  Bool m_nonPackedConstraintFlag;
  Bool m_frameOnlyConstraintFlag;

  // coding structure
  Int       m_iIntraPeriod;                                   ///< period of I-slice (random access period)
  Int       m_iDecodingRefreshType;                           ///< random access type
  Int       m_iGOPSize;                                       ///< GOP size of hierarchical structure
#if JVET_D0135_PARAMS
  Bool      m_bReWriteParamSetsFlag;                          ///< Flag to enable rewriting of parameter sets at random access points
#endif  
  Int       m_extraRPSs;                                      ///< extra RPSs added to handle CRA
  GOPEntry  m_GOPList[MAX_GOP];                               ///< the coding structure entries from the config file
  Int       m_numReorderPics[MAX_TLAYER];                     ///< total number of reorder pictures
  Int       m_maxDecPicBuffering[MAX_TLAYER];                 ///< total number of pictures in the decoded picture buffer
  Bool      m_crossComponentPredictionEnabledFlag;            ///< flag enabling the use of cross-component prediction
  Bool      m_reconBasedCrossCPredictionEstimate;             ///< causes the alpha calculation in encoder search to be based on the decoded residual rather than the pre-transform encoder-side residual
  UInt      m_log2SaoOffsetScale[MAX_NUM_CHANNEL_TYPE];       ///< number of bits for the upward bit shift operation on the decoded SAO offsets
  Bool      m_useTransformSkip;                               ///< flag for enabling intra transform skipping
  Bool      m_useTransformSkipFast;                           ///< flag for enabling fast intra transform skipping
  UInt      m_log2MaxTransformSkipBlockSize;                  ///< transform-skip maximum size (minimum of 2)
  Bool      m_transformSkipRotationEnabledFlag;               ///< control flag for transform-skip/transquant-bypass residual rotation
  Bool      m_transformSkipContextEnabledFlag;                ///< control flag for transform-skip/transquant-bypass single significance map context
  Bool      m_rdpcmEnabledFlag[NUMBER_OF_RDPCM_SIGNALLING_MODES];///< control flags for residual DPCM
  Bool      m_enableAMP;
  Bool      m_persistentRiceAdaptationEnabledFlag;            ///< control flag for Golomb-Rice parameter adaptation over each slice
  Bool      m_cabacBypassAlignmentEnabledFlag;

  // coding quality
#if JVET_E0059_FLOATING_POINT_QP_FIX
  OptionalValue<UInt> m_qpIncrementAtSourceFrame;             ///< Optional source frame number at which all subsequent frames are to use an increased internal QP.
#else
  Double    m_fQP;                                            ///< QP value of key-picture (floating point)
#endif
  Int       m_iQP;                                            ///< QP value of key-picture (integer)
#if JCTVC_X0038_LAMBDA_FROM_QP_CAPABILITY
  Double    m_dIntraQpFactor;                                 ///< Intra Q Factor. If negative, use a default equation: 0.57*(1.0 - Clip3( 0.0, 0.5, 0.05*(Double)(isField ? (GopSize-1)/2 : GopSize-1) ))
  Int       m_intraQPOffset;                                  ///< QP offset for intra slice (integer)
  Bool      m_lambdaFromQPEnable;                             ///< enable flag for QP:lambda fix
#endif
  Char*     m_pchdQPFile;                                     ///< QP offset for each slice (initialized from external file)
  Int*      m_aidQP;                                          ///< array of slice QP values
  Int       m_iMaxDeltaQP;                                    ///< max. |delta QP|
  UInt      m_uiDeltaQpRD;                                    ///< dQP range for multi-pass slice QP optimization
  Int       m_iMaxCuDQPDepth;                                 ///< Max. depth for a minimum CuDQPSize (0:default)
  Int       m_diffCuChromaQpOffsetDepth;                      ///< If negative, then do not apply chroma qp offsets.
  Bool      m_bFastDeltaQP;                                   ///< Fast Delta QP (false:default)

  Int       m_cbQpOffset;                                     ///< Chroma Cb QP Offset (0:default)
  Int       m_crQpOffset;                                     ///< Chroma Cr QP Offset (0:default)
#if WCG_LUMA_DQP_CM_SCALE
  WCGChromaQPControl m_wcgChromaQpControl;                    ///< Wide-colour-gamut chroma QP control.
  LumaLevelToDeltaQPMapping m_lumaLevelToDeltaQPMapping;      ///< mapping from luma level to Delta QP.
#endif
#if ADAPTIVE_QP_SELECTION
  Bool      m_bUseAdaptQpSelect;
#endif
  TComSEIMasteringDisplay m_masteringDisplay;

  Bool      m_bUseAdaptiveQP;                                 ///< Flag for enabling QP adaptation based on a psycho-visual model
  Int       m_iQPAdaptationRange;                             ///< dQP range by QP adaptation

  Int       m_maxTempLayer;                                  ///< Max temporal layer

  // coding unit (CU) definition
  // TODO: Remove MaxCUWidth/MaxCUHeight and replace with MaxCUSize.
#if JVET_C0024_QTBT
  UInt      m_uiCTUSize;
  UInt      m_uiMinQT[3]; //0: I slice luma; 1: I slice chroma; 2: P/B slice
#if JVET_C0024_SPS_MAX_BT_DEPTH
  UInt      m_uiMaxBTDepth;
  UInt      m_uiMaxBTDepthISliceL;
  UInt      m_uiMaxBTDepthISliceC;
#endif
#if JVET_C0024_SPS_MAX_BT_SIZE
  UInt      m_uiMaxBTSize;
  UInt      m_uiMaxBTSizeISliceL;
  UInt      m_uiMaxBTSizeISliceC;
#endif
#if EXTENSION_360_VIDEO
  UInt      m_uiMaxCUWidth;                                   ///< max. CU width in pixel
  UInt      m_uiMaxCUHeight;                                  ///< max. CU height in pixel
  UInt      m_uiMaxCUDepth;                                   ///< max. CU depth (as specified by command line)
#endif
#else
  UInt      m_uiMaxCUWidth;                                   ///< max. CU width in pixel
  UInt      m_uiMaxCUHeight;                                  ///< max. CU height in pixel
  UInt      m_uiMaxCUDepth;                                   ///< max. CU depth (as specified by command line)
#endif
  UInt      m_uiMaxTotalCUDepth;                              ///< max. total CU depth - includes depth of transform-block structure, used for assign number of buffers to store mode in CU level
#if !JVET_C0024_QTBT
  UInt      m_uiLog2DiffMaxMinCodingBlockSize;                ///< difference between largest and smallest CU depth
#endif

  // transform unit (TU) definition
  UInt      m_uiQuadtreeTULog2MaxSize;
  UInt      m_uiQuadtreeTULog2MinSize;

  UInt      m_uiQuadtreeTUMaxDepthInter;
  UInt      m_uiQuadtreeTUMaxDepthIntra;

  // coding tools (bit-depth)
  Int       m_inputBitDepth   [MAX_NUM_CHANNEL_TYPE];         ///< bit-depth of input file
  Int       m_outputBitDepth  [MAX_NUM_CHANNEL_TYPE];         ///< bit-depth of output file
  Int       m_MSBExtendedBitDepth[MAX_NUM_CHANNEL_TYPE];      ///< bit-depth of input samples after MSB extension
  Int       m_internalBitDepth[MAX_NUM_CHANNEL_TYPE];         ///< bit-depth codec operates at (input/output files will be converted)
  Bool      m_extendedPrecisionProcessingFlag;
  Bool      m_highPrecisionOffsetsEnabledFlag;

  //coding tools (chroma format)
  ChromaFormat m_chromaFormatIDC;

#if BLOCK_BASED_BACKGROUND_REFERENCE
	//coding tools (background reference)
	Bool      m_bUseBBR;
	Bool      m_bUseBBRwithDQP;
	Bool      m_bUseBBRReplaceOneRef;
#endif

  // coding tools (PCM bit-depth)
  Bool      m_bPCMInputBitDepthFlag;                          ///< 0: PCM bit-depth is internal bit-depth. 1: PCM bit-depth is input bit-depth.

  // coding tool (SAO)
  Bool      m_bUseSAO;
  Bool      m_bTestSAODisableAtPictureLevel;
  Double    m_saoEncodingRate;                                ///< When >0 SAO early picture termination is enabled for luma and chroma
  Double    m_saoEncodingRateChroma;                          ///< The SAO early picture termination rate to use for chroma (when m_SaoEncodingRate is >0). If <=0, use results for luma.
  Int       m_maxNumOffsetsPerPic;                            ///< SAO maximum number of offset per picture
  Bool      m_saoCtuBoundary;                                 ///< SAO parameter estimation using non-deblocked pixels for CTU bottom and right boundary areas
  // coding tools (loop filter)
  Bool      m_bLoopFilterDisable;                             ///< flag for using deblocking filter
  Bool      m_loopFilterOffsetInPPS;                         ///< offset for deblocking filter in 0 = slice header, 1 = PPS
  Int       m_loopFilterBetaOffsetDiv2;                     ///< beta offset for deblocking filter
  Int       m_loopFilterTcOffsetDiv2;                       ///< tc offset for deblocking filter
  Bool      m_DeblockingFilterMetric;                         ///< blockiness metric in encoder
#if CNN_BASED_INLOOP_FILTER_SPS
  Bool      m_bUseCNNLoopFilter;
#endif
#if COM16_C806_LMCHROMA
  Bool      m_useLMChroma;                                    ///< JL: Cross component prediction, chroma from luma or Cr from Cb with linear model
#endif

  // coding tools (PCM)
  Bool      m_usePCM;                                         ///< flag for using IPCM
  UInt      m_pcmLog2MaxSize;                                 ///< log2 of maximum PCM block size
  UInt      m_uiPCMLog2MinSize;                               ///< log2 of minimum PCM block size
  Bool      m_bPCMFilterDisableFlag;                          ///< PCM filter disable flag
  Bool      m_enableIntraReferenceSmoothing;                  ///< flag for enabling(default)/disabling intra reference smoothing/filtering

  // coding tools (encoder-only parameters)
  Bool      m_bUseASR;                                        ///< flag for using adaptive motion search range
  Bool      m_bUseHADME;                                      ///< flag for using HAD in sub-pel ME
  Bool      m_useRDOQ;                                       ///< flag for using RD optimized quantization
  Bool      m_useRDOQTS;                                     ///< flag for using RD optimized quantization for transform skip
#if T0196_SELECTIVE_RDOQ
  Bool      m_useSelectiveRDOQ;                               ///< flag for using selective RDOQ
#endif
  Int       m_rdPenalty;                                      ///< RD-penalty for 32x32 TU for intra in non-intra slices (0: no RD-penalty, 1: RD-penalty, 2: maximum RD-penalty)
  Bool      m_bDisableIntraPUsInInterSlices;                  ///< Flag for disabling intra predicted PUs in inter slices.
  Int       m_iFastSearch;                                    ///< ME mode, 0 = full, 1 = diamond, 2 = PMVFAST
  Int       m_iSearchRange;                                   ///< ME search range
  Int       m_bipredSearchRange;                              ///< ME search range for bipred refinement
  Int       m_minSearchWindow;                                ///< ME minimum search window size for the Adaptive Window ME
  Bool      m_bClipForBiPredMeEnabled;                        ///< Enables clipping for Bi-Pred ME.
  Bool      m_bFastMEAssumingSmootherMVEnabled;               ///< Enables fast ME assuming a smoother MV.
  Bool      m_bUseFastEnc;                                    ///< flag for using fast encoder setting
  Bool      m_bUseEarlyCU;                                    ///< flag for using Early CU setting
  Bool      m_useFastDecisionForMerge;                        ///< flag for using Fast Decision Merge RD-Cost
  Bool      m_bUseCbfFastMode;                              ///< flag for using Cbf Fast PU Mode Decision
  Bool      m_useEarlySkipDetection;                         ///< flag for using Early SKIP Detection
#if JVET_D0077_SAVE_LOAD_ENC_INFO
  Bool      m_useSaveLoadEncInfo;
  Bool      m_useSaveLoadSplitDecision;
#endif
  Int       m_sliceMode;                                     ///< 0: no slice limits, 1 : max number of CTBs per slice, 2: max number of bytes per slice,
                                                             ///< 3: max number of tiles per slice
  Int       m_sliceArgument;                                 ///< argument according to selected slice mode
  Int       m_sliceSegmentMode;                              ///< 0: no slice segment limits, 1 : max number of CTBs per slice segment, 2: max number of bytes per slice segment,
                                                             ///< 3: max number of tiles per slice segment
  Int       m_sliceSegmentArgument;                          ///< argument according to selected slice segment mode

  Bool      m_bLFCrossSliceBoundaryFlag;  ///< 1: filter across slice boundaries 0: do not filter across slice boundaries
  Bool      m_bLFCrossTileBoundaryFlag;   ///< 1: filter across tile boundaries  0: do not filter across tile boundaries
  Bool      m_tileUniformSpacingFlag;
  Int       m_numTileColumnsMinus1;
  Int       m_numTileRowsMinus1;
  std::vector<Int> m_tileColumnWidth;
  std::vector<Int> m_tileRowHeight;
  Int       m_iWaveFrontSynchro; //< 0: no WPP. >= 1: WPP is enabled, the "Top right" from which inheritance occurs is this LCU offset in the line above the current.
  Int       m_iWaveFrontFlush; //< enable(1)/disable(0) the CABAC flush at the end of each line of LCUs.

  Bool      m_bUseConstrainedIntraPred;                       ///< flag for using constrained intra prediction
  Bool      m_bFastUDIUseMPMEnabled;
  Bool      m_bFastMEForGenBLowDelayEnabled;
  Bool      m_bUseBLambdaForNonKeyLowDelayPictures;

  Int       m_decodedPictureHashSEIEnabled;                    ///< Checksum(3)/CRC(2)/MD5(1)/disable(0) acting on decoded picture hash SEI message
  Int       m_recoveryPointSEIEnabled;
  Int       m_bufferingPeriodSEIEnabled;
  Int       m_pictureTimingSEIEnabled;
  Bool      m_toneMappingInfoSEIEnabled;
  Bool      m_chromaSamplingFilterSEIenabled;
  Int       m_chromaSamplingHorFilterIdc;
  Int       m_chromaSamplingVerFilterIdc;
  Int       m_toneMapId;
  Bool      m_toneMapCancelFlag;
  Bool      m_toneMapPersistenceFlag;
  Int       m_toneMapCodedDataBitDepth;
  Int       m_toneMapTargetBitDepth;
  Int       m_toneMapModelId;
  Int       m_toneMapMinValue;
  Int       m_toneMapMaxValue;
  Int       m_sigmoidMidpoint;
  Int       m_sigmoidWidth;
  Int       m_numPivots;
  Int       m_cameraIsoSpeedIdc;
  Int       m_cameraIsoSpeedValue;
  Int       m_exposureIndexIdc;
  Int       m_exposureIndexValue;
  Bool      m_exposureCompensationValueSignFlag;
  Int       m_exposureCompensationValueNumerator;
  Int       m_exposureCompensationValueDenomIdc;
  Int       m_refScreenLuminanceWhite;
  Int       m_extendedRangeWhiteLevel;
  Int       m_nominalBlackLevelLumaCodeValue;
  Int       m_nominalWhiteLevelLumaCodeValue;
  Int       m_extendedWhiteLevelLumaCodeValue;
  Int*      m_startOfCodedInterval;
  Int*      m_codedPivotValue;
  Int*      m_targetPivotValue;
  Int       m_framePackingSEIEnabled;
  Int       m_framePackingSEIType;
  Int       m_framePackingSEIId;
  Int       m_framePackingSEIQuincunx;
  Int       m_framePackingSEIInterpretation;
  Int       m_segmentedRectFramePackingSEIEnabled;
  Bool      m_segmentedRectFramePackingSEICancel;
  Int       m_segmentedRectFramePackingSEIType;
  Bool      m_segmentedRectFramePackingSEIPersistence;
  Int       m_displayOrientationSEIAngle;
  Int       m_temporalLevel0IndexSEIEnabled;
  Int       m_gradualDecodingRefreshInfoEnabled;
  Int       m_noDisplaySEITLayer;
  Int       m_decodingUnitInfoSEIEnabled;
  Int       m_SOPDescriptionSEIEnabled;
  Int       m_scalableNestingSEIEnabled;
  Bool      m_tmctsSEIEnabled;
  Bool      m_timeCodeSEIEnabled;
  Int       m_timeCodeSEINumTs;
  TComSEITimeSet m_timeSetArray[MAX_TIMECODE_SEI_SETS];
  Bool      m_kneeSEIEnabled;
  Int       m_kneeSEIId;
  Bool      m_kneeSEICancelFlag;
  Bool      m_kneeSEIPersistenceFlag;
  Int       m_kneeSEIInputDrange;
  Int       m_kneeSEIInputDispLuminance;
  Int       m_kneeSEIOutputDrange;
  Int       m_kneeSEIOutputDispLuminance;
  Int       m_kneeSEINumKneePointsMinus1;
  Int*      m_kneeSEIInputKneePoint;
  Int*      m_kneeSEIOutputKneePoint;
  // weighted prediction
  Bool      m_useWeightedPred;                    ///< Use of weighted prediction in P slices
  Bool      m_useWeightedBiPred;                  ///< Use of bi-directional weighted prediction in B slices

  UInt      m_log2ParallelMergeLevel;                         ///< Parallel merge estimation region
  UInt      m_maxNumMergeCand;                                ///< Max number of merge candidates

  Int       m_TMVPModeId;
  Bool      m_signHideFlag;
  Bool      m_RCEnableRateControl;                ///< enable rate control or not
  Int       m_RCTargetBitrate;                    ///< target bitrate when rate control is enabled
  Int       m_RCKeepHierarchicalBit;              ///< 0: equal bit allocation; 1: fixed ratio bit allocation; 2: adaptive ratio bit allocation
  Bool      m_RCLCULevelRC;                       ///< true: LCU level rate control; false: picture level rate control NOTE: code-tidy - rename to m_RCCtuLevelRC
  Bool      m_RCUseLCUSeparateModel;              ///< use separate R-lambda model at LCU level                        NOTE: code-tidy - rename to m_RCUseCtuSeparateModel
  Int       m_RCInitialQP;                        ///< inital QP for rate control
  Bool      m_RCForceIntraQP;                     ///< force all intra picture to use initial QP or not
  ScalingListMode m_useScalingListId;                         ///< using quantization matrix
  Char*     m_scalingListFile;                                ///< quantization matrix file name

  Bool      m_TransquantBypassEnableFlag;                     ///< transquant_bypass_enable_flag setting in PPS.
  Bool      m_CUTransquantBypassFlagForce;                    ///< if transquant_bypass_enable_flag, then, if true, all CU transquant bypass flags will be set to true.
  CostMode  m_costMode;                                       ///< Cost mode to use

  Bool      m_recalculateQPAccordingToLambda;                 ///< recalculate QP value according to the lambda value
  Bool      m_useStrongIntraSmoothing;                        ///< enable strong intra smoothing for 32x32 blocks where the reference samples are flat
  Int       m_activeParameterSetsSEIEnabled;

  Bool      m_vuiParametersPresentFlag;                       ///< enable generation of VUI parameters
  Bool      m_aspectRatioInfoPresentFlag;                     ///< Signals whether aspect_ratio_idc is present
  Int       m_aspectRatioIdc;                                 ///< aspect_ratio_idc
  Int       m_sarWidth;                                       ///< horizontal size of the sample aspect ratio
  Int       m_sarHeight;                                      ///< vertical size of the sample aspect ratio
  Bool      m_overscanInfoPresentFlag;                        ///< Signals whether overscan_appropriate_flag is present
  Bool      m_overscanAppropriateFlag;                        ///< Indicates whether conformant decoded pictures are suitable for display using overscan
  Bool      m_videoSignalTypePresentFlag;                     ///< Signals whether video_format, video_full_range_flag, and colour_description_present_flag are present
  Int       m_videoFormat;                                    ///< Indicates representation of pictures
  Bool      m_videoFullRangeFlag;                             ///< Indicates the black level and range of luma and chroma signals
  Bool      m_colourDescriptionPresentFlag;                   ///< Signals whether colour_primaries, transfer_characteristics and matrix_coefficients are present
  Int       m_colourPrimaries;                                ///< Indicates chromaticity coordinates of the source primaries
  Int       m_transferCharacteristics;                        ///< Indicates the opto-electronic transfer characteristics of the source
  Int       m_matrixCoefficients;                             ///< Describes the matrix coefficients used in deriving luma and chroma from RGB primaries
  Bool      m_chromaLocInfoPresentFlag;                       ///< Signals whether chroma_sample_loc_type_top_field and chroma_sample_loc_type_bottom_field are present
  Int       m_chromaSampleLocTypeTopField;                    ///< Specifies the location of chroma samples for top field
  Int       m_chromaSampleLocTypeBottomField;                 ///< Specifies the location of chroma samples for bottom field
  Bool      m_neutralChromaIndicationFlag;                    ///< Indicates that the value of all decoded chroma samples is equal to 1<<(BitDepthCr-1)
  Bool      m_defaultDisplayWindowFlag;                       ///< Indicates the presence of the default window parameters
  Int       m_defDispWinLeftOffset;                           ///< Specifies the left offset from the conformance window of the default window
  Int       m_defDispWinRightOffset;                          ///< Specifies the right offset from the conformance window of the default window
  Int       m_defDispWinTopOffset;                            ///< Specifies the top offset from the conformance window of the default window
  Int       m_defDispWinBottomOffset;                         ///< Specifies the bottom offset from the conformance window of the default window
  Bool      m_frameFieldInfoPresentFlag;                      ///< Indicates that pic_struct values are present in picture timing SEI messages
  Bool      m_pocProportionalToTimingFlag;                    ///< Indicates that the POC value is proportional to the output time w.r.t. first picture in CVS
  Int       m_numTicksPocDiffOneMinus1;                       ///< Number of ticks minus 1 that for a POC difference of one
  Bool      m_bitstreamRestrictionFlag;                       ///< Signals whether bitstream restriction parameters are present
  Bool      m_tilesFixedStructureFlag;                        ///< Indicates that each active picture parameter set has the same values of the syntax elements related to tiles
  Bool      m_motionVectorsOverPicBoundariesFlag;             ///< Indicates that no samples outside the picture boundaries are used for inter prediction
  Int       m_minSpatialSegmentationIdc;                      ///< Indicates the maximum size of the spatial segments in the pictures in the coded video sequence
  Int       m_maxBytesPerPicDenom;                            ///< Indicates a number of bytes not exceeded by the sum of the sizes of the VCL NAL units associated with any coded picture
  Int       m_maxBitsPerMinCuDenom;                           ///< Indicates an upper bound for the number of bits of coding_unit() data
  Int       m_log2MaxMvLengthHorizontal;                      ///< Indicate the maximum absolute value of a decoded horizontal MV component in quarter-pel luma units
  Int       m_log2MaxMvLengthVertical;                        ///< Indicate the maximum absolute value of a decoded vertical MV component in quarter-pel luma units
#if COM16_C806_VCEG_AZ10_SUB_PU_TMVP
  UInt      m_subPUTLog2Size;
  Bool      m_useAtmvpFlag;
#endif
#if COM16_C806_OBMC
  Bool      m_useOBMC;
  Int       m_OBMCBlkSize;
#endif
#if VCEG_AZ07_FRUC_MERGE
  Int       m_useFRUCMgrMode;
  Int       m_FRUCRefineFilter;
  Int       m_FURCRefineRange;
  Int       m_FRUCSmallBlkRefineDepth;
#endif
#if VCEG_AZ07_IMV
  Bool      m_useIMV;
  Int       m_IMVMaxCand;
#endif
#if VCEG_AZ06_IC
  Bool      m_useIC;
#endif
#if ALF_HM3_REFACTOR
  Bool      m_useALF;                                         ///< flag for using adaptive loop filter
#endif
#if COM16_C806_LARGE_CTU
  Int       m_useFastLCTU;
#endif

#if COM16_C806_EMT
  Int       m_useEMT;                                         ///< XZ: Enhanced Multiple Transform
  Int       m_useFastEMT;                                     ///< XZ: Fast Methods of Enhanced Multiple Transform
#endif
#if VCEG_AZ08_USE_KLT
  Int       m_useKLT;                                         
#endif

#if VCEG_AZ07_INTRA_4TAP_FILTER
  Bool      m_useIntra4TapFilter;
#endif
#if VCEG_AZ07_INTRA_BOUNDARY_FILTER
  Bool      m_useIntraBoundaryFilter;
#endif
#if VCEG_AZ05_BIO
  Bool      m_useBIO;
#endif
#if JVET_E0052_DMVR
  Bool      m_useDMVR;
#endif
#if VCEG_AZ05_INTRA_MPI
  Bool      m_useMPI;
#endif
#if COM16_C1046_PDPC_INTRA
  Bool     m_usePDPC;
#endif
#if MLIP_EXTEND_REF_LINE
  Bool    m_useMLIP;
#endif
#if RRC_INTRA_SHARPENING
  Bool     m_useIntraSharp;
#endif
#if VCEG_AZ05_ROT_TR
  Bool      m_useROT;
#elif COM16_C1044_NSST
  Bool      m_useNSST;
#endif
#if COM16_C1016_AFFINE
  Bool      m_useAffine;
#endif
#if MRG_CONSIDER_COMPLEX
  Bool      m_useMrgComplex;
#endif
#if JVET_D0033_ADAPTIVE_CLIPPING
  ClipParam m_ClipParam;
#endif
#if JVET_F0096_BILATERAL_FILTER
  Bool      m_useBilateralFilter;
#endif
#if BARC_INTER_SPS
	Bool      m_useBarcInter;
#endif
#if BARC_INTRA_CHROMA_SPS
	Bool      m_useBarcIntraChroma;
#endif
#if BARC_INTRA_LUMA_SPS
	Bool      m_useBarcIntraLuma;
#endif
#if INTERPOLATE_FRAME_SPS
	Bool      m_useFrameBasedFRUC;
#endif
#if RRC_DWDIP
  Bool      m_useDWDIP;
#endif
#if RRC_NS
	Bool      m_useNS;                                         ///< flag for using noise suppressor loop filter
#endif
  std::string m_summaryOutFilename;                           ///< filename to use for producing summary output file.
  std::string m_summaryPicFilenameBase;                       ///< Base filename to use for producing summary picture output files. The actual filenames used will have I.txt, P.txt and B.txt appended.
  UInt        m_summaryVerboseness;                           ///< Specifies the level of the verboseness of the text output.

#if EXTENSION_360_VIDEO
  TExt360AppEncCfg m_ext360;
  friend class TExt360AppEncCfg;
  friend class TExt360AppEncTop;
#endif

  // internal member functions
  Void  xCheckParameter ();                                   ///< check validity of configuration values
  Void  xPrintParameter ();                                   ///< print configuration values
  Void  xPrintUsage     ();                                   ///< print usage
public:
  TAppEncCfg();
  virtual ~TAppEncCfg();

public:
  Void  create    ();                                         ///< create option handling class
  Void  destroy   ();                                         ///< destroy option handling class
  Bool  parseCfg  ( Int argc, Char* argv[] );                 ///< parse configuration file to fill member variables

};// END CLASS DEFINITION TAppEncCfg

//! \}

#endif // __TAPPENCCFG__

