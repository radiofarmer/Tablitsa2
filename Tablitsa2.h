#pragma once

#define DEBUG_VECTOR 0

#if !_DEBUG || DEBUG_VECTOR
#define VECTOR
//#define VECTOR_VOICE_EFFECTS_TEST
//#define VECTOR_MASTER_EFFECTS_TEST
#endif

#include "IPlug_include_in_plug_hdr.h"
#include "IControls.h"

#define PRESET_NAME_CHAR_LENGTH 32
#define TABLITSA2_MIDI_MOD_NAMES { "Velocity", "Aftertouch", "Mod Wheel", "Keytrack", "Random"/*, "CC 1", "CC 2"*/}
#define TABLITSA2_MAX_VOICE_EFFECTS 3
#define TABLITSA2_MAX_MASTER_EFFECTS 3
#define TABLITSA2_VOICE_EFFECTS_LIST {"None", "Sample & Hold", "Texturizer", "Distortion", "Super Ring"}
#define TABLITSA2_MASTER_EFFECTS_LIST {"None", "Delay", "EQ", "Reverb 1", "Reverb 2", "Chorus"}

const int kNumPresets = 1;
constexpr int kNumVoices = 16;
constexpr int kMaxUnisonVoices = 8;
constexpr int kNumSeqSteps = 16;

// This must be in the same order as the effect labels. (If you change the effect order, this is all you need to modify)
// TODO: Put all effects in alphabetical order in the next versioin
enum EVoiceEffectTypes
{
  kNoVoiceEffect = 0,
  kSampleAndHoldEffect,
  kTexturizerEffect,
  kDistortionEffect,
  kCMEffect,
  kLimiterEffect,
  kNumVoiceEffectTypes
};

enum EMasterEffectTypes
{
  kNoMasterEffect = 0,
  kDelayEffect,
  kEQEffect,
  kReverbEffect,
  kReverb2Effect,
  kChorusEffect,
  kNumMasterEffectTypes
};

enum EModulators
{
  kEnvelope1 = 0,
  kEnvelope2,
  kAmpEnvelope,
  kLFO1,
  kLFO2,
  kSequencer,
  kVelocity,
  kAftertouch,
  kModWheel,
  kKeytrack,
  kRandom,
  /*kMidiCC1,
  kMidiCC2,
  kMidiCC3,
  kMidiCC4,
  kMidiCC5,
  kMidiCC6,*/
  kNumModulators
};

enum EParams
{
  kParamGain = 0,
  kParamPan,
  kParamNoteGlideTime,
  kParamNoteGlideRate,
  kParamPortamentoMode,
  kParamMonophonic,
  kParamLegato,
  kParamUnisonVoices,
  kParamUnisonDetune,
  kParamUnisonChord,
  kParamStereoSpread,
  kParamVibratoSpeed,
  kParamVibratoDepth,
  kParamAmpEnvAttack,
  kParamAmpEnvDecay,
  kParamAmpEnvSustain,
  kParamAmpEnvRelease,
  kParamAmpEnvVelocity,
  kParamAmpEnvDecayCurve,
  kParamAmpEnvReleaseCurve,
  kParamEnv1Attack,
  kParamEnv1Decay,
  kParamEnv1Sustain,
  kParamEnv1Release,
  kParamEnv1Velocity,
  kParamEnv1DecayCurve,
  kParamEnv1ReleaseCurve,
  kParamEnv1Level,
  kParamEnv2Attack,
  kParamEnv2Decay,
  kParamEnv2Sustain,
  kParamEnv2Release,
  kParamEnv2Velocity,
  kParamEnv2DecayCurve,
  kParamEnv2ReleaseCurve,
  kParamLFO1Shape,
  kParamLFO1Phase,
  kParamLFO1RateHz,
  kParamLFO1RateTempo,
  kParamLFO1Amp,
  kParamLFO1RateMode,
  kParamLFO1Restart,
  kParamLFO2Shape,
  kParamLFO2Phase,
  kParamLFO2RateHz,
  kParamLFO2RateTempo,
  kParamLFO2Amp,
  kParamLFO2RateMode,
  kParamLFO2Restart,
  kParamSequencerSteps,
  kParamSequencerStepMode,
  kParamSequencerGlide,
  kParamSequencerCurve,
  kParamSequencerRateHz,
  kParamSequencerRateTempo,
  kParamSequencerRateMode,
  kParamSequencerRestart,
  kParamSequencerAmp,
  kParamMidiVelocityCurve,
  kParamMidiAftertouchCurve,
  kParamMidiModWheelCurve,
  kParamMidiKeytrackCurve,
  kParamMidiRandomCurve,
  kParamWavetable1Pitch, // Voice-specific (polyphonic) parameters
  kParamWavetable1Pos,
  kParamWavetable1Bend,
  kParamWavetable1Amp,
  kParamOsc1FormantOn,
  kParamWavetable1Formant,
  kParamWavetable2Pitch,
  kParamWavetable2Pos,
  kParamWavetable2Bend,
  kParamWavetable2Amp,
  kParamOsc2FormantOn,
  kParamWavetable2Formant,
  kParamFilter1Cutoff, // Filter 1
  kParamFilter1Resonance,
  kParamFilter1Drive,
  kParamFilter1FF,
  kParamFilter1FB,
  kParamFilter1Delay,
  kParamFilter1Type,
  kParamFilter1ModeVSF,
  kParamFilter1ModeMoog,
  kParamFilter1ModeComb,
  kParamFilter1Osc1Send,
  kParamFilter1Osc2Send,
  kParamFilter2Cutoff, // Filter 2
  kParamFilter2Resonance,
  kParamFilter2Drive,
  kParamFilter2FF,
  kParamFilter2FB,
  kParamFilter2Delay,
  kParamFilter2Type,
  kParamFilter2ModeVSF,
  kParamFilter2ModeMoog,
  kParamFilter2ModeComb,
  kParamFilter2Osc1Send,
  kParamFilter2Osc2Send,
  kParamOscModulator,
  kParamOsc1PM,
  kParamOsc1RM,
  kParamOsc2PM,
  kParamOsc2RM,
  kParamPhaseModFreq,
  kParamPhaseModAmount,
  kParamRingModFreq,
  kParamRingModAmount,
  kParamVoiceEffect1Param1,
  kParamVoiceEffect1Param2,
  kParamVoiceEffect1Param3,
  kParamVoiceEffect1Param4,
  kParamVoiceEffect1Param5,
  kParamVoiceEffect1Param6,
  kParamVoiceEffect2Param1,
  kParamVoiceEffect2Param2,
  kParamVoiceEffect2Param3,
  kParamVoiceEffect2Param4,
  kParamVoiceEffect2Param5,
  kParamVoiceEffect2Param6,
  kParamVoiceEffect3Param1,
  kParamVoiceEffect3Param2,
  kParamVoiceEffect3Param3,
  kParamVoiceEffect3Param4,
  kParamVoiceEffect3Param5,
  kParamVoiceEffect3Param6,
  kParamMasterEffect1Param1,
  kParamMasterEffect1Param2,
  kParamMasterEffect1Param3,
  kParamMasterEffect1Param4,
  kParamMasterEffect1Param5,
  kParamMasterEffect1Param6,
  kParamMasterEffect2Param1,
  kParamMasterEffect2Param2,
  kParamMasterEffect2Param3,
  kParamMasterEffect2Param4,
  kParamMasterEffect2Param5,
  kParamMasterEffect2Param6,
  kParamMasterEffect3Param1,
  kParamMasterEffect3Param2,
  kParamMasterEffect3Param3,
  kParamMasterEffect3Param4,
  kParamMasterEffect3Param5,
  kParamMasterEffect3Param6,
  kParamOsc1FilterBypass,
  kParamOsc2FilterBypass,
  kParamOsc1EffectBypass,
  kParamOsc2EffectBypass,
  kNumParams
};

enum EVoiceModParams
{
  kVPan = 0,
  kVWavetable1PitchOffset,
  kVWavetable1Position,
  kVWavetable1Bend,
  kVWavetable1Formant,
  kVWavetable1Amp,
  kVWavetable2PitchOffset,
  kVWavetable2Position,
  kVWavetable2Bend,
  kVWavetable2Formant,
  kVWavetable2Amp,
  kVFilter1Cutoff,
  kVFilter1Resonance,
  kVFilter1Drive,
  kVFilter1FF,
  kVFilter1FB,
  kVFilter1Delay,
  kVFilter2Cutoff,
  kVFilter2Resonance,
  kVFilter2Drive,
  kVFilter2FF,
  kVFilter2FB,
  kVFilter2Delay,
  kVPhaseModFreq,
  kVPhaseModAmt,
  kVRingModFreq,
  kVRingModAmt,
  kVEffect1Param1,
  kVEffect1Param2,
  kVEffect1Param3,
  kVEffect1Param4,
  kVEffect2Param1,
  kVEffect2Param2,
  kVEffect2Param3,
  kVEffect2Param4,
  kVEffect3Param1,
  kVEffect3Param2,
  kVEffect3Param3,
  kVEffect3Param4,
  kNumVoiceModulations,
  kVEnv1Sustain,
  kVEnv2Sustain,
  kVAmpEnvSustain,
  kVLFO1RateHz,
  kVLFO1Amp,
  kVLFO2RateHz,
  kVLFO2Amp,
  kVSequencerRateHz,
  kVSequencerAmp
};

constexpr int voiceMetaModulations[]{
  kVEnv1Sustain,
  kVEnv2Sustain,
  kVAmpEnvSustain,
  kVLFO1RateHz,
  kVLFO1Amp,
  kVLFO2RateHz,
  kVLFO2Amp,
  kVSequencerRateHz,
  kVSequencerAmp
};

constexpr int kNumVoiceMetaModulations{sizeof(voiceMetaModulations) / sizeof(int)};

constexpr int kNumVoiceEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;
constexpr int kNumMasterEffectParams = kParamMasterEffect2Param1 - kParamMasterEffect1Param1;

#if IPLUG_DSP
// will use EParams in Tablitsa2_DSP.h
#include "Tablitsa2_DSP.h"
#endif

enum EControlTags
{
  kCtrlTagPeriodicTable = 0,
  kCtrlTagMeter,
  kCtrlTagLFO1Vis,
  kCtrlTagLFO2Vis,
  kCtrlTagScope,
  kCtrlTagRTText,
  kCtrlTagKeyboard,
  kCtrlTagBender,
  kCtrlTagGlideMode,
  kCtrlTagMidiModList,
  kCtrlTagEnv1Depth, // Modulator depths
  kCtrlTagEnv2Depth,
  kCtrlTagAmpEnvDepth,
  kCtrlTagLFO1Depth,
  kCtrlTagLFO2Depth,
  kCtrlTagSequencerDepth,
  kCtrlTagMidiVelocityDepth,
  kCtrlTagMidiAftertouchDepth,
  kCtrlTagMidiModWheelDepth,
  kCtrlTagMidiKeytrackDepth,
  kCtrlTagRandomDepth,
  /*kCtrlTagMidiCC1Depth,
  kCtrlTagMidiCC2Depth,
  kCtrlTagMidiCC3Depth,
  kCtrlTagMidiCC4Depth,
  kCtrlTagMidiCC5Depth,
  kCtrlTagMidiCC6Depth, */// !Modulator depths
  kCtrlTagMidiControlCurve,
  kCtrlTagLFO1RateMode,
  kCtrlTagLFO2RateMode,
  kCtrlTagSequencerRateMode,
  kCtrlTagSequencerQuant,
  kCtrlTagOsc1Formant,
  kCtrlTagOsc2Formant,
  kCtrlTagFilter1Cutoff,
  kCtrlTagFilter1Resonance,
  kCtrlTagFilter1Drive,
  kCtrlTagFilter1FF,
  kCtrlTagFilter1FB,
  kCtrlTagFilter1Delay,
  kCtrlTagFilter1Mode,
  kCtrlTagFilter1Type,
  kCtrlTagFilter1Osc1,
  kCtrlTagFilter1Osc2,
  kCtrlTagFilter2Cutoff,
  kCtrlTagFilter2Resonance,
  kCtrlTagFilter2Drive,
  kCtrlTagFilter2FF,
  kCtrlTagFilter2FB,
  kCtrlTagFilter2Delay,
  kCtrlTagFilter2Mode,
  kCtrlTagFilter2Type,
  kCtrlTagFilter2Osc1,
  kCtrlTagFilter2Osc2,
  kCtrlTagOscModSource,
  kCtrlTagOscModFreq,
  kCtrlTagOscModAmt,
  kCtrlTagOsc1ModSwitch,
  kCtrlTagOsc2ModSwitch,
  kCtrlTagLFO1Plot,
  kCtrlTagLFO2Plot,
  kCtrlTagSequencer,
  kCtrlTagEffectBank,
  kCtrlTagVoiceEffectsList, // Voice effect controls
  kCtrlTagVoiceEffectsSwitch,
  kCtrlTagVoiceEffectsKnob1,
  kCtrlTagVoiceEffectsKnob2,
  kCtrlTagVoiceEffectsKnob3,
  kCtrlTagVoiceEffectsKnob4,
  kCtrlTagVoiceEffectsToggle1,
  kCtrlTagVoiceEffectsToggle2,
  kCtrlTagMasterEffectsList, // Master effect controls
  kCtrlTagMasterEffectsSwitch,
  kCtrlTagMasterEffectsKnob1,
  kCtrlTagMasterEffectsKnob2,
  kCtrlTagMasterEffectsKnob3,
  kCtrlTagMasterEffectsKnob4,
  kCtrlTagMasterEffectsToggle1,
  kCtrlTagMasterEffectsToggle2,
  kCtrlTagDefaultPresetList,
  kNumCtrlTags
};

enum EMsgTags
{
  kMsgWavetable1Changed = 0,
  kMsgWavetable2Changed,
  kMsgSeqSliderChanged,
  kMsgUpdateLFO1Plot,
  kMsgUpdateLFO2Plot,
  kMsgRandomizeSequencer,
  kMsgSavePreset,
  kMsgLoadPreset,
  kMsgLoadDefaultPreset,
  kMsgVoiceEffect1Changed,
  kMsgVoiceEffect2Changed,
  kMsgVoiceEffect3Changed,
  kMsgMasterEffect1Changed,
  kMsgMasterEffect2Changed,
  kMsgMasterEffect3Changed
};

constexpr EControlTags kStartupTriggerControls[]{
  kCtrlTagLFO1RateMode,
  kCtrlTagLFO2RateMode,
  kCtrlTagSequencerRateMode,
  kCtrlTagSequencerQuant,
//  kCtrlTagFilter1Mode,
//  kCtrlTagFilter2Mode,
  kCtrlTagFilter1Type,
  kCtrlTagFilter2Type,
  kCtrlTagGlideMode
};

constexpr EControlTags kModSliders[]{
  kCtrlTagEnv1Depth,
  kCtrlTagEnv2Depth,
  kCtrlTagAmpEnvDepth,
  kCtrlTagLFO1Depth,
  kCtrlTagLFO2Depth,
  kCtrlTagSequencerDepth,
  kCtrlTagMidiVelocityDepth,
  kCtrlTagMidiAftertouchDepth,
  kCtrlTagMidiModWheelDepth,
  kCtrlTagMidiKeytrackDepth,
  kCtrlTagRandomDepth,
};

constexpr EControlTags kMidiModKeys[]{
  kCtrlTagMidiVelocityDepth,
  kCtrlTagMidiAftertouchDepth,
  kCtrlTagMidiModWheelDepth,
  kCtrlTagMidiKeytrackDepth,
  kCtrlTagRandomDepth,
 /* kCtrlTagMidiCC1Depth,
  kCtrlTagMidiCC2Depth,
  kCtrlTagMidiCC3Depth,
  kCtrlTagMidiCC4Depth,
  kCtrlTagMidiCC5Depth,
  kCtrlTagMidiCC6Depth,*/
};

constexpr int kNumModSliders{sizeof(kModSliders) / sizeof(EControlTags)};

using namespace iplug;
using namespace igraphics;

class Tablitsa2 final : public Plugin
{
public:
  Tablitsa2(const InstanceInfo& info);

#if IPLUG_DSP // http://bit.ly/2S64BDd
public:
  void ProcessBlock(sample** inputs, sample** outputs, int nFrames) override;
  void ProcessMidiMsg(const IMidiMsg& msg) override;
  void OnReset() override;
  void OnParamChange(int paramIdx) override;
  void OnIdle() override;
  bool OnMessage(int msgTag, int ctrlTag, int dataSize, const void* pData) override;
  
  void OnUIOpen() override;
  bool SerializeState(IByteChunk& chunk) const override;
  int UnserializeState(const IByteChunk& chunk, int startPos) override;
  void UpdateUIControls();
  bool LoadWavetables();
  /* implement this and return true to trigger your custom about box, when someone clicks about in the menu of a standalone app or VST3 plugin */
  bool OnHostRequestingAboutBox() override; // See IPlugAPP_dialog.cpp
  /* implement this and return true to trigger your custom help info, when someone clicks help in the menu of a standalone app or VST3 plugin */
  bool OnHostRequestingProductHelp() override;

  IByteChunk LoadPreset(const char* filename="UserPreset", bool isBackup=false);
  void SavePreset(IByteChunk& byteData, const char* filename = "UserPreset", bool isBackup=false);
  void LoadDefaultState();
  int CheckVersion(const IByteChunk& presetData);
  bool ShowLoadErrorMessageBox();

  void SetMasterFXSlot(int slotIdx, int masterEffectIdx) { mCurrentMasterFXSlot = slotIdx; mMasterEffectSlots[slotIdx] = masterEffectIdx; }
  void SetVoiceFXSlot(int slotIdx, int voiceEffectIdx) { mCurrentVoiceFXSlot = slotIdx; mVoiceEffectSlots[slotIdx] = voiceEffectIdx; }

  // Store the status of the tempo sync control, required for setting the correct parameter range during preset loading
  void SetDelayTempoSync(int slotIdx, bool tempoSync) { mDelayTempoSync[slotIdx] = tempoSync; }
  void SetDelayTempoSync(bool tempoSync) { SetDelayTempoSync(mCurrentMasterFXSlot, tempoSync); }

  // Modulation Control Functions
  int GetActiveModIdx() const;
  void SetActiveModIdx(int idx);
  int GetFirstModCtrlTag() const { return kCtrlTagEnv1Depth; }
  int GetLastModCtrlTag() const { return kCtrlTagRandomDepth; }
  double GetModMatrixSlot(const int parameter, const int modulator);
  void SetModMatrixSlot(const int parameter, const int modulator, const double value);

  void RefreshEffectBankControl();

private:
  Tablitsa2DSP<sample> mDSP {kNumVoices}; // sample is an alias for double
  IPeakSender<2> mMeterSender;
  IControl* mActiveControl{};
  IByteChunk mStateBackup; // Stores the last user state before reseting, so that it may be restored

  // UI Status variables not stored as parameters
  char mPresetName[PRESET_NAME_CHAR_LENGTH]{};
  int mPresetID{};
  int mVoiceEffectSlots[TABLITSA2_MAX_VOICE_EFFECTS]{}; // Holds the ID number of the effect in each slot for the voice effects
  int mMasterEffectSlots[TABLITSA2_MAX_MASTER_EFFECTS]{}; // Holds the ID number of the effect in each slot for the master effects
  int mCurrentVoiceFXSlot{ 0 }; // The effect slot current open for editing. Controled by the Slide-Switch controls
  int mCurrentMasterFXSlot{ 0 };
  int mCurrentEffectsTab{ 1 };
  bool mDelayTempoSync[TABLITSA2_MAX_MASTER_EFFECTS]{};
  bool mMonoDelay[TABLITSA2_MAX_MASTER_EFFECTS]{};
  double mSequencerIsQuantized{ 0. };

  double mModMatrix[kNumParams][kNumModulators]{0.};
  int mActiveModIdx{ -1 };
#endif
};

std::string GetDataPath(char* appendPath="\\");

std::vector<char> ReadAllBytes(const char* fname); // For reading preset files of arbitrary length

#include "Tablitsa2Controls.h"