
#include "Tablitsa2.h"
#include "IPlug_include_in_plug_src.h"
#include "Tablitsa2UIFunctions.h"
#include "Presets.h"

#include <atlstr.h>

#define UI_NROWS 4
#define UI_NCOLS 3
#define KEYBOARD
#define FILTER_PARTITION 0.8
#define ENVELOPE_PARTITION 0.8
#define LFO_PARTITION 1.4

// Load/save presets with file dialog or internal controls
#define LOAD_WITH_DIALOG 1

Tablitsa2::Tablitsa2(const InstanceInfo& info)
: Plugin(info, MakeConfig(kNumParams, kNumPresets))
{
  /*Master Controls*/
#if VST3_API
  GetParam(kParamGain)->InitGain("Gain", 0., -24., 12.);
#else
  GetParam(kParamGain)->InitGain("Gain", -3., -72., 0.);
#endif
  GetParam(kParamNoteGlideTime)->InitDouble("Note Glide Time", 0., 0.0, 3000., 5., "ms", 0, "Portamento", IParam::ShapePowCurve(3.));
  GetParam(kParamNoteGlideRate)->InitDouble("Note Glide Rate", 12., 0.5, 128., 0.1, "St/s", 0, "Portamento", IParam::ShapePowCurve(3.));
  GetParam(kParamPortamentoMode)->InitBool("Equal Glide Time", true);
  GetParam(kParamMonophonic)->InitBool("Monophonic", false);
  GetParam(kParamUnisonVoices)->InitInt("Unison Voices", 1, 1, 8);
  GetParam(kParamUnisonDetune)->InitDouble("Unison Detune", 0., 0., 12., 0.01, "St.");
  GetParam(kParamUnisonChord)->InitEnum("Chord", 0, { UNISON_CHORD_LIST });
  GetParam(kParamStereoSpread)->InitDouble("Stereo Spread", 0., -180., 180., 0.5, "°");
  GetParam(kParamVibratoSpeed)->InitDouble("Vibrato Speed", 5., 0., 10., 0.01, "Hz", 0, "Vibrato");
  GetParam(kParamVibratoDepth)->InitDouble("Vibrato Depth", 0., 0., 50., 0.1, "Cts", 0, "Vibrato", IParam::ShapePowCurve(2.));
  // Pan (modulatable)
  GetParam(kParamPan)->InitDouble("Pan", 0., -90., 90., 0.5, "°");
  /*Modulators*/
  // Envelope 1
  GetParam(kParamEnv1Attack)->InitDouble("Env1 Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv1Decay)->InitDouble("Env1 Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv1Sustain)->InitDouble("Env1 Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv1Release)->InitDouble("Env1 Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv1Velocity)->InitDouble("Env1 Velocity Mod", 0., 0., 1., 0.01);
  GetParam(kParamEnv1DecayCurve)->InitPercentage("Env1 Decay Curve", 100.);
  GetParam(kParamEnv1ReleaseCurve)->InitPercentage("Env1 Release Curve", 50.);
  // Envelope 2
  GetParam(kParamEnv2Attack)->InitDouble("Env2 Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv2Decay)->InitDouble("Env2 Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv2Sustain)->InitDouble("Env2 Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv2Release)->InitDouble("Env2 Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv2Velocity)->InitDouble("Env2 Velocity Mod", 0., 0., 1., 0.01);
  GetParam(kParamEnv2DecayCurve)->InitPercentage("Env1 Decay Curve", 100.);
  GetParam(kParamEnv2ReleaseCurve)->InitPercentage("Env1 Release Curve", 50.);
  // Amplitude Envelopoe
  GetParam(kParamAmpEnvAttack)->InitDouble("AmpEnv Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamAmpEnvDecay)->InitDouble("AmpEnv Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamAmpEnvSustain)->InitDouble("AmpEnv Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamAmpEnvRelease)->InitDouble("AmpEnv Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamAmpEnvVelocity)->InitDouble("AmpEnv Velocity Mod", 1., 0., 1., 0.01);
  GetParam(kParamAmpEnvDecayCurve)->InitPercentage("Env1 Decay Curve", 100.);
  GetParam(kParamAmpEnvReleaseCurve)->InitPercentage("Env1 Release Curve", 50.);
  // LFO 1
  GetParam(kParamLFO1Shape)->InitEnum("LFO1 Shape", LFO<>::kTriangle, {TABLITSA2_LFO_SHAPE_VALIST});
  GetParam(kParamLFO1Phase)->InitDouble("LFO1 Phase", 0., 0., 1., 0.01);
  GetParam(kParamLFO1RateHz)->InitFrequency("LFO1 Rate", 1., 0.01, 40.);
  GetParam(kParamLFO1RateTempo)->InitEnum("LFO1 Rate", LFO<>::k1, {LFO_TEMPODIV_VALIST});
  GetParam(kParamLFO1RateMode)->InitBool("LFO1 Sync", false);
  GetParam(kParamLFO1Restart)->InitBool("LFO1 Retrigger", false);
  GetParam(kParamLFO1Amp)->InitDouble("LFO Amp", 1., -1., 1., 0.01, "", IParam::kFlagsNone, "LFO");
  // LFO 2
  GetParam(kParamLFO2Shape)->InitEnum("LFO2 Shape", LFO<>::kTriangle, { TABLITSA2_LFO_SHAPE_VALIST });
  GetParam(kParamLFO2Phase)->InitDouble("LFO2 Phase", 0., 0., 1., 0.01);
  GetParam(kParamLFO2RateHz)->InitFrequency("LFO2 Rate", 1., 0.01, 40.);
  GetParam(kParamLFO2RateTempo)->InitEnum("LFO2 Rate", LFO<>::k1, { LFO_TEMPODIV_VALIST });
  GetParam(kParamLFO2RateMode)->InitBool("LFO2 Sync", false);
  GetParam(kParamLFO2Restart)->InitBool("LFO2 Retrigger", false);
  GetParam(kParamLFO2Amp)->InitDouble("LFO Amp", 1., -1., 1., 0.01, "", IParam::kFlagsNone, "LFO");
  // Sequencer
  GetParam(kParamSequencerRateHz)->InitFrequency("Sequencer Rate", 1., 0.01, 40.);
  GetParam(kParamSequencerStepMode)->InitEnum("Sequencer Step Mode", Sequencer<iplug::sample>::EStepMode::kAuto, { "Auto", "Gate", "Trigger" });
  GetParam(kParamSequencerRateTempo)->InitEnum("Sequencer Rate", LFO<>::k4th, { LFO_TEMPODIV_VALIST });
  GetParam(kParamSequencerRateMode)->InitBool("Sequencer Sync", true);
  GetParam(kParamSequencerRestart)->InitBool("Sequencer Retrigger", false);
  GetParam(kParamSequencerSteps)->InitInt("Sequencer Length", 16, 1, 16);
  GetParam(kParamSequencerGlide)->InitPercentage("Sequencer Glide");
//  GetParam(kParamSequencerCurve)->InitEnum()
  GetParam(kParamSequencerAmp)->InitDouble("Seqeuencer Amp", 1., 0., 1., 0.01, "", IParam::kFlagsNone, "Sequencer");

  /*
  Wavetables
  */
  // Wavetable 1 Pitch
  GetParam(kParamWavetable1Pitch)->InitDouble("Oscillator 1 Pitch", 0., -24., 24., 0.1);
  // Wavetable 1 Position
  GetParam(kParamWavetable1Pos)->InitDouble("Oscillator 1 Wavetable Position", 0., 0., 1., 0.01);
  // Wavetable 1 Bend
  GetParam(kParamWavetable1Bend)->InitDouble("Oscillator 1 Wavetable Bend", 0., -1., 1., 0.01);
  // Wavetable 1 Formant
  GetParam(kParamOsc1FormantOn)->InitBool("Oscillator 1 Formant On", false);
  GetParam(kParamWavetable1Formant)->InitPercentage("Oscillator 1 Formant", 0.);
  // Wavetable 1 Amp
  GetParam(kParamWavetable1Amp)->InitDouble("WT1 Amp", 0.5, 0., 1., 0.01);

  // Wavetable 2 Pitch
  GetParam(kParamWavetable2Pitch)->InitDouble("WT2 Pitch", 0., -24., 24., 0.1);
  // Wavetable 2 Position
  GetParam(kParamWavetable2Pos)->InitDouble("WT2 Pos", 0., 0., 1., 0.01);
  // Wavetable 2 Bend
  GetParam(kParamWavetable2Bend)->InitDouble("WT2 Bend", 0., -1., 1., 0.01);
  // Wavetable 2 Formant
  GetParam(kParamOsc2FormantOn)->InitBool("WT2 Formant On", false);
  GetParam(kParamWavetable2Formant)->InitPercentage("WT2 Formant", 0.);
  // Wavetable 2 Amp
  GetParam(kParamWavetable2Amp)->InitDouble("WT2 Amp", 0.5, 0., 1., 0.01);

  // Modulation Oscillator
  GetParam(kParamOscModulator)->InitEnum("Oscillator Modulation Mode", 0, { "Phase", "Ring" });
  GetParam(kParamOsc1PM)->InitBool("Osc1 Phase Modulation", false);
  GetParam(kParamOsc1RM)->InitBool("Osc1 Ring Modulation", false);
  GetParam(kParamOsc2PM)->InitBool("Osc2 Phase Modulation", false);
  GetParam(kParamOsc2RM)->InitBool("Osc2 Ring Modulation", false);
  // Phase Mod Freq
  GetParam(kParamPhaseModFreq)->InitDouble("Phase Modulator Freq", 0., -24., 24., 0.01);
  // Phase Mod Amount
  GetParam(kParamPhaseModAmount)->InitPercentage("Phase Modulator Amount", 25., 0., 100.);
  // Ring Mod Freq
  GetParam(kParamRingModFreq)->InitDouble("Ring Modulator Freq", 0., -24., 24., 0.01);
  // Ring Mod Amount
  GetParam(kParamRingModAmount)->InitPercentage("Ring Modulator Amount", 25., 0., 100.);

  // Filter 1
  GetParam(kParamFilter1Type)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_TYPE_LIST });
  GetParam(kParamFilter1ModeVSF)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_MODE_LIST_VSF });
  GetParam(kParamFilter1ModeMoog)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_MODE_LIST_MOOG });
  GetParam(kParamFilter1ModeComb)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_MODE_LIST_COMB });
  GetParam(kParamFilter1Cutoff)->InitFrequency("Filter 1 Cutoff", 20000., 20., 20000., 1.);
  GetParam(kParamFilter1Resonance)->InitPercentage("Filter 1 Resonance", 0.1, 0., 99., 1);
  GetParam(kParamFilter1Drive)->InitPercentage("Filter 1 Drive", 0.1, 0., 100., 1);
  // Comb-specific parameters with different ranges (the modulation parameters are common between all filter types, as these always between -1 and 1)
  GetParam(kParamFilter1FF)->InitPercentage("Filter 1 FF");
  GetParam(kParamFilter1FB)->InitPercentage("Filter 1 FB");
  GetParam(kParamFilter1Delay)->InitDouble("Filter 1 Delay", 0.1, 0.1, 20., 0.01, "ms");
  // Oscillator sends
  GetParam(kParamFilter1Osc1Send)->InitGain("Filter 1 Osc 1 Send", 0., SEND_DB_FLOOR, 0.);
  GetParam(kParamFilter1Osc2Send)->InitGain("Filter 1 Osc 2 Send", SEND_DB_FLOOR, SEND_DB_FLOOR, 0.);

  // Filter 2
  GetParam(kParamFilter2Type)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_TYPE_LIST });
  GetParam(kParamFilter2ModeVSF)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_MODE_LIST_VSF });
  GetParam(kParamFilter2ModeMoog)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_MODE_LIST_MOOG });
  GetParam(kParamFilter2ModeComb)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_MODE_LIST_COMB });
  GetParam(kParamFilter2Cutoff)->InitFrequency("Filter 2 Cutoff", 20000., 20., 20000., 1.);
  GetParam(kParamFilter2Resonance)->InitPercentage("Filter 2 Resonance", 0.1, 0., 99., 1);
  GetParam(kParamFilter2Drive)->InitPercentage("Filter 2 Drive", 0.1, 0., 100., 1);
  // Comb-specific parameters with different ranges (the modulation parameters are common between all filter types, as these always between -1 and 1)
  GetParam(kParamFilter2FF)->InitPercentage("Filter 2 FF");
  GetParam(kParamFilter2FB)->InitPercentage("Filter 2 FB");
  GetParam(kParamFilter2Delay)->InitDouble("Filter 2 Delay", 0.1, 0.1, 20., 0.01, "ms");
  // Oscillator sends
  GetParam(kParamFilter2Osc1Send)->InitGain("Filter 2 Osc 1 Send", SEND_DB_FLOOR, SEND_DB_FLOOR, 0.);
  GetParam(kParamFilter2Osc2Send)->InitGain("Filter 2 Osc 2 Send", 0., SEND_DB_FLOOR, 0.);

  // Effect Parameters
  // Voice Effects
  GetParam(kParamVoiceEffect1Param1)->InitDouble("Voice Effect 1 Parameter 1", 0., 0., 1., 0.01);
  // P2
  GetParam(kParamVoiceEffect1Param2)->InitDouble("VoiceEffect 1 Parameter 2", 0., 0., 1., 0.01);
  // P3
  GetParam(kParamVoiceEffect1Param3)->InitPercentage("VoiceEffect 1 Parameter 3", 0.);
  // P4
  GetParam(kParamVoiceEffect1Param4)->InitPercentage("VoiceEffect 1 Parameter 4", 0.);
  GetParam(kParamVoiceEffect1Param5)->InitBool("VoiceEffect 1 Parameter 5", false);
  GetParam(kParamVoiceEffect1Param6)->InitBool("VoiceEffect 1 Parameter 6", false);

  // Voice Effect 2
  GetParam(kParamVoiceEffect2Param1)->InitDouble("Voice Effect 2 Parameter 1", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2)->InitDouble("Voice Effect 2 Parameter 2", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3)->InitPercentage("VoiceEffect 2 Parameter 3", 0.);
  GetParam(kParamVoiceEffect2Param4)->InitPercentage("Voice Effect 2 Parameter 4", 0.);
  GetParam(kParamVoiceEffect2Param5)->InitBool("Voice Effect 2 Parameter 5", false);
  GetParam(kParamVoiceEffect2Param6)->InitBool("Voice Effect 2 Parameter 6", false);

  // Voice Effect 3
  GetParam(kParamVoiceEffect3Param1)->InitDouble("Voice Effect 3 Parameter 1", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2)->InitDouble("VoiceEffect 3 Parameter 2", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3)->InitPercentage("Voice Effect 3 Parameter 3", 0.);
  GetParam(kParamVoiceEffect3Param4)->InitPercentage("Voice Effect 3 Parameter 4", 0.);
  GetParam(kParamVoiceEffect3Param5)->InitBool("VoiceEffect 3 Parameter 5", false);
  GetParam(kParamVoiceEffect3Param6)->InitBool("VoiceEffect 3 Parameter 6", false);

  // Master Effects
  GetParam(kParamMasterEffect1Param1)->InitDouble("Master Effect 1 Parameter 1", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect1Param2)->InitDouble("Master Effect 1 Parameter 2", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect1Param3)->InitPercentage("Master Effect 1 Parameter 3", 0.);
  GetParam(kParamMasterEffect1Param4)->InitPercentage("Master Effect 1 Parameter 4", 0.);
  GetParam(kParamMasterEffect1Param5)->InitBool("Master Effect 1 Parameter 5", false);
  GetParam(kParamMasterEffect1Param6)->InitBool("Master Effect 1 Parameter 6", false);

  GetParam(kParamMasterEffect2Param1)->InitDouble("Master Effect 2 Parameter 1", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect2Param2)->InitDouble("Master Effect 2 Parameter 2", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect2Param3)->InitPercentage("Master Effect 2 Parameter 3", 0.);
  GetParam(kParamMasterEffect2Param4)->InitPercentage("Master Effect 2 Parameter 4", 0.);
  GetParam(kParamMasterEffect2Param5)->InitBool("Master Effect 2 Parameter 5", false);
  GetParam(kParamMasterEffect2Param6)->InitBool("Master Effect 2 Parameter 6", false);

  GetParam(kParamMasterEffect3Param1)->InitDouble("Master Effect 3 Parameter 1", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect3Param2)->InitDouble("Master Effect 3 Parameter 2", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect3Param3)->InitPercentage("Master Effect 3 Parameter 3", 0.);
  GetParam(kParamMasterEffect3Param4)->InitPercentage("Master Effect 3 Parameter 4", 0.);
  GetParam(kParamMasterEffect3Param5)->InitBool("Master Effect 3 Parameter 5", false);
  GetParam(kParamMasterEffect3Param6)->InitBool("Master Effect 3 Parameter 6", false);

  GetParam(kParamOsc1FilterBypass)->InitDouble("Osc 1 Filter Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));
  GetParam(kParamOsc2FilterBypass)->InitDouble("Osc 2 Filter Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));
  GetParam(kParamOsc1EffectBypass)->InitDouble("Osc 1 Effect Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));
  GetParam(kParamOsc2EffectBypass)->InitDouble("Osc 2 Effect Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));

#if VST3_API
#define KEYBOARD_HEIGHT 0.f
#else
#define KEYBOARD_HEIGHT 150.f
#endif
#if IPLUG_EDITOR // http://bit.ly/2S64BDd
  mMakeGraphicsFunc = [&]() {
    return MakeGraphics(*this, PLUG_WIDTH, (int)(PLUG_HEIGHT + KEYBOARD_HEIGHT), PLUG_FPS, GetScaleForScreen((int)(PLUG_WIDTH), (int)(PLUG_HEIGHT + KEYBOARD_HEIGHT)));
  };
  mLayoutFunc = [&](IGraphics* pGraphics) {

    pGraphics->AttachCornerResizer(EUIResizerMode::Scale, false);
    pGraphics->AttachBackground(BG_FN);
    pGraphics->EnableMouseOver(true);
    pGraphics->EnableMultiTouch(true);

#ifdef OS_WEB
    pGraphics->AttachPopupMenuControl();
#endif

    //    pGraphics->EnableLiveEdit(true);
    pGraphics->LoadFont("Roboto-Regular", ROBOTO_FN);
    constexpr float plugBorders = 10.f;
    const IRECT b = pGraphics->GetBounds().GetPadded(-plugBorders);
    constexpr float dialSize{ 80.f };
    IRECT globalControls = b.GetOffset(0, 0, 0, -KEYBOARD_HEIGHT).GetHPadded(-5.f).GetReducedFromTop(15.f);
    const IRECT appControls = b.GetFromTop(30.f);
#ifndef VST3_API
    IRECT keyboardBounds = b.GetFromBottom((float)KEYBOARD_HEIGHT - plugBorders).GetPadded(-5.f, -25.f, -5.f, -2.f);
    IRECT wheelsBounds = keyboardBounds.ReduceFromLeft(100.f).GetPadded(-10.f);

    // Keyboard
    pGraphics->AttachControl(new IVKeyboardControl(keyboardBounds, 19, 79), kCtrlTagKeyboard);
    // Mod wheels
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5)), kCtrlTagBender);
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5, true), IMidiMsg::EControlChangeMsg::kModWheel));
    //    pGraphics->AttachControl(new IVMultiSliderControl<4>(b.GetGridCell(0, 2, 2).GetPadded(-30), "", TABLITSA2_STYLE, kParamAttack, EDirection::Vertical, 0.f, 1.f));

#endif
    const ISVG tableImg = pGraphics->LoadSVG(PERIODIC_TABLE_FN);
    const ISVG routingImg = pGraphics->LoadSVG(ROUTING_FN);
    const IBitmap logoImg = pGraphics->LoadBitmap(LOGO_FN);

    // Periodic table
    pGraphics->AttachControl(new PeriodicTable(globalControls.GetReducedFromBottom(globalControls.H() * 0.5f).GetFromLeft(720.f).GetPadded(-30.f).GetVShifted(20.f), tableImg, { kNoParameter }),
      kCtrlTagPeriodicTable);
    int default_elements[]{ 1, 2 };
    dynamic_cast<PeriodicTable*>(pGraphics->GetControlWithTag(kCtrlTagPeriodicTable))->LoadValues(default_elements);
    pGraphics->AttachControl(new IBitmapControl(globalControls.GetReducedFromTop(90.f).GetReducedFromLeft(110.f).GetFromTLHC(280.f, 80.), logoImg));

    // Preset controls
    const IRECT presetPanel = globalControls.GetFromTop(50.f).GetFromLeft(720.f).GetHPadded(-20.f);
    IActionFunction savePresetFunction = IActionFunction([this, pGraphics](IControl* pControl) {
      bool saveBackup = false;
      SendArbitraryMsgFromUI(kMsgSavePreset, kNoTag, sizeof(bool), reinterpret_cast<const void*>(&saveBackup));
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCustomStr(mPresetName);
      });
    IActionFunction loadPresetFunction = IActionFunction([this](IControl* pControl) {
      bool loadBackup = false;
      SendArbitraryMsgFromUI(kMsgLoadPreset, kNoTag, sizeof(bool), reinterpret_cast<const void*>(&loadBackup));
      });
    IActionFunction resetToDefaultFunction = IActionFunction([this, pGraphics](IControl* pControl) {
      IByteChunk currentState;
      SerializeState(currentState);
      SavePreset(currentState, "", true);
      LoadDefaultState();
      strcpy(mPresetName, "Default");
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCustomStr(mPresetName);
      ResetFilterControls(pGraphics, 3);
      });
    IActionFunction loadBackupFunction = IActionFunction([this](IControl* pControl) {
      IByteChunk presetData = LoadPreset("", true);
      if (presetData.Size() == 0)
        return;
      UnserializeState(presetData, 0);
      UpdateUIControls();
      });
    const IVStyle ioButtonStyle = pushButtonStyle.WithColor(EVColor::kFG, COLOR_WHITE.WithOpacity(0.2));
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(0, 1, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Save", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(savePresetFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(0, 2, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Load", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(loadPresetFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(1, 1, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Reset", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(resetToDefaultFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(1, 2, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Backup", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(loadBackupFunction);
    // Built-in file selection
    pGraphics->AttachControl(new ITextControl(presetPanel.SubRectVertical(2, 0).GetFromLeft(110.f).GetPadded(-2.f), "Default Presets:", TABLITSA2_TEXT), kNoTag, "Presets");
    pGraphics->AttachControl(new DropdownListControl(presetPanel.SubRectVertical(2, 1).GetFromLeft(110.f).GetPadded(-2.f), PRESET_NAMES), kCtrlTagDefaultPresetList, "Presets")->SetActionFunction(
      [this, pGraphics](IControl* pControl) {
        mPresetID = dynamic_cast<DropdownListControl*>(pControl)->GetCurrentIndex();
        dynamic_cast<DropdownListControl*>(pControl)->SetCustomStr("");
        // Load Resource
        HMODULE hInstance = static_cast<HMODULE>(pGraphics->GetWinModuleHandle());
#ifdef VST3_API
        std::wstring presetName = PRESET_ID_LIST_W[mPresetID];
        HRSRC hRes{ FindResource(hInstance, presetName.c_str(), (LPWSTR)RT_RCDATA) };
        HGLOBAL pData{ LoadResource(hInstance, hRes) };
#else
        std::string presetName = PRESET_ID_LIST[mPresetID];
        HRSRC hRes{ FindResource(hInstance, presetName.c_str(), RT_RCDATA) };
        HGLOBAL pData{ LoadResource(hInstance, hRes) };
#endif
        memcpy(mPresetName, PRESET_ID_LIST[mPresetID], PRESET_NAME_CHAR_LENGTH);
        if (!pData)
          return;
        char* resLock{ static_cast<char*>(LockResource(pData)) };
        const int byteLen = SizeofResource(hInstance, hRes);
        // Set State
        IByteChunk preset;
        preset.PutBytes(pData, byteLen);
        UnserializeState(preset, 0);
        UpdateUIControls();
        FreeResource(pData);
      }
    );
    const auto defaultPresetIdx{std::find(PRESET_ID_LIST, PRESET_ID_LIST + N_PRESETS, "Default") - PRESET_ID_LIST};
    dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCurrentIndex(defaultPresetIdx >= N_PRESETS ? 0 : defaultPresetIdx, false);
#if !LOAD_WITH_DIALOG
    const IRECT presetMenuBounds = presetPanel.GetGridCell(1, 1, 2, 4).GetVShifted(presetPanel.GetGridCell(1, 1, 2, 5).H()).GetPadded(0.f, 0.f, 100.f, 100.f);
    IPopupMenuControl* presetSelectionMenu = new IPopupMenuControl(kNoParameter, TABLITSA2_TEXT);
    pGraphics->AttachControl(new PresetSelector(presetPanel.GetGridCell(1, 1, 2, 4).GetHPadded(-10.f), presetSelectionMenu, {"Preset 1", "Preset 2"}));
#endif
    // Wavetable controls
    const IRECT osc1Panel = globalControls.GetReducedFromLeft(720.f).GetReducedFromBottom(globalControls.H() / 2.f).GetGridCell(0, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc1Panel.GetGridCell(0, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable1Pitch, kVWavetable1PitchOffset, "Pitch", modKnobStyle, true), kNoTag, "Osc1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc1Panel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable1Amp, kVWavetable1Amp, "Amp"), kNoTag, "Osc1");
    pGraphics->AttachControl(new Tablitsa2VToggleButton(osc1Panel.GetGridCell(0, 2, 2, 3).GetVPadded(-25.f), kParamOsc1FormantOn, "Fmnt. On", modKnobStyle), kNoTag, "Osc1")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        IControl* fmntKnob{ pGraphics->GetControlWithTag(kCtrlTagOsc1Formant) };
        if (fmntKnob)
          fmntKnob->SetDisabled(pControl->GetValue() < 0.5);
      }
    );
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc1Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable1Pos, kVWavetable1Position, "Valency"), kNoTag, "Osc1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc1Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable1Bend, kVWavetable1Bend, "Isotope"), kNoTag, "Osc1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc1Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable1Formant, kVWavetable1Formant, "Excitation"), kNoTag, "Osc1");
    pGraphics->AttachControl(new Tablitsa2VGroupControl("Wavetable 1", "Osc1", -5.f, 15.f, -5.f, 2.f));

    const IRECT osc2Panel = globalControls.GetReducedFromLeft(720.f).GetReducedFromBottom(globalControls.H() / 2.f).GetGridCell(1, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc2Panel.GetGridCell(0, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable2Pitch, kVWavetable2PitchOffset, "Pitch", modKnobStyle, true), kNoTag, "Osc2");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc2Panel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable2Amp, kVWavetable2Amp, "Amp"), kNoTag, "Osc2");
    pGraphics->AttachControl(new Tablitsa2VToggleButton(osc2Panel.GetGridCell(0, 2, 2, 3).GetVPadded(-25.f), kParamOsc2FormantOn, "Fmnt. On", modKnobStyle), kNoTag, "Osc2")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        IControl* fmntKnob{ pGraphics->GetControlWithTag(kCtrlTagOsc2Formant) };
        if (fmntKnob)
          fmntKnob->SetDisabled(pControl->GetValue() < 0.5);
      }
    );
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc2Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable2Pos, kVWavetable1Position, "Valency"), kNoTag, "Osc2");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc2Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable2Bend, kVWavetable2Bend, "Isotope"), kNoTag, "Osc2");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(osc2Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable2Formant, kVWavetable2Formant, "Excitation"), kNoTag, "Osc2");
    pGraphics->AttachControl(new Tablitsa2VGroupControl("Wavetable 2", "Osc2", -5.f, 15.f, -5.f, 2.f));

    /* Master Controls */
    const IRECT masterPanel = globalControls.GetReducedFromLeft(720.f + osc1Panel.W() + 30.f).GetReducedFromBottom(globalControls.H() / 2.f).GetPadded(-10.f);
    // Gain, Pan, and Volume Meter
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(masterPanel.GetGridCell(0, 0, 4, 6).GetVShifted(-15.f).GetCentredInside(dialSize), kParamGain, "Gain"), kNoTag, "Master");
    // pGraphics->AttachControl(new IBKnobControl(masterPanel.GetGridCell(0, 0, 4, 6).GetCentredInside(dialSize), dialImg, kParamGain), kNoTag, "Master");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(masterPanel.GetGridCell(0, 1, 4, 6).GetVShifted(-15.f).GetCentredInside(dialSize), kParamPan, kVPan, "Pan"), kNoTag, "Master");
    // MIDI Control Modulators
    const IRECT staticModPanel = masterPanel.GetGridCell(1, 0, 2, 3).GetHAlignedTo(masterPanel.GetGridCell(0, 0, 1, 6), EAlign::Near).GetReducedFromRight(19.f).GetPadded(-15.f);
    pGraphics->AttachControl(new DropdownListControl(staticModPanel.GetFromTop(40.f).GetHShifted(-4.f).GetPadded(-2.f), TABLITSA2_MIDI_MOD_NAMES, TABLITSA2_TEXT.WithSize(18.f)), kCtrlTagMidiModList, "MidiMods")->SetActionFunction(
      [this, pGraphics](IControl* pControl){
        const int selected{kMidiModKeys[dynamic_cast<DropdownListControl*>(pControl)->GetCurrentIndex()]};
        // Get pointer to response curve control
        // IControl* pCurveControl{ pGraphics->GetControlWithTag(kCtrlTagMidiControlCurve) };
        for (int i{}; i < sizeof(kMidiModKeys) / sizeof(EControlTags); ++i) {
          IControl* pSlider{ pGraphics->GetControlWithTag(kMidiModKeys[i]) };
          if (pSlider && kMidiModKeys[i] != selected)
            pSlider->Hide(true);
          else if (pSlider)
            pSlider->Hide(false);
        }
      }
    );
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(25.f).GetFromTop(30.f).GetHPadded(-8.f), EModulators::kVelocity, "Velocity", TABLITSA2_STYLE.WithShowLabel(false), false, EDirection::Horizontal), kCtrlTagMidiVelocityDepth, "MidiMods");
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(25.f).GetFromTop(30.f).GetHPadded(-8.f), EModulators::kAftertouch, "Aftertouch", TABLITSA2_STYLE.WithShowLabel(false), false, EDirection::Horizontal), kCtrlTagMidiAftertouchDepth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(25.f).GetFromTop(30.f).GetHPadded(-8.f), EModulators::kKeytrack, "Keytrack", TABLITSA2_STYLE.WithShowLabel(false), false, EDirection::Horizontal), kCtrlTagMidiKeytrackDepth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(25.f).GetFromTop(30.f).GetHPadded(-8.f), EModulators::kModWheel, "Modwheel", TABLITSA2_STYLE.WithShowLabel(false), false, EDirection::Horizontal), kCtrlTagMidiModWheelDepth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(25.f).GetFromTop(30.f).GetHPadded(-8.f), EModulators::kRandom, "Trg. Rnd.", TABLITSA2_STYLE.WithShowLabel(false), false, EDirection::Horizontal), kCtrlTagRandomDepth, "MidiMods")->Hide(true);

    // Curve Display
    pGraphics->AttachControl(new Tablitsa2ParametricPlotControl(staticModPanel.GetFromBottom(90.f).GetPadded(-2.f), 50, kNoParameter, "Response Curve", TABLITSA2_STYLE.WithLabelText(TABLITSA2_TEXT.WithSize(15.f))), kNoTag, "MidiMods");

    /*pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(50.f).GetPadded(5.f), EModulators::kMidiCC1, "CC1", TABLITSA2_STYLE.WithShowLabel(false), false), kCtrlTagMidiCC1Depth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(50.f).GetPadded(5.f), EModulators::kMidiCC2, "CC2", TABLITSA2_STYLE.WithShowLabel(false), false), kCtrlTagMidiCC2Depth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(50.f).GetPadded(5.f), EModulators::kMidiCC3, "CC3", TABLITSA2_STYLE.WithShowLabel(false), false), kCtrlTagMidiCC3Depth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(50.f).GetPadded(5.f), EModulators::kMidiCC4, "CC4", TABLITSA2_STYLE.WithShowLabel(false), false), kCtrlTagMidiCC4Depth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(50.f).GetPadded(5.f), EModulators::kMidiCC5, "CC5", TABLITSA2_STYLE.WithShowLabel(false), false), kCtrlTagMidiCC5Depth, "MidiMods")->Hide(true);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetReducedFromTop(50.f).GetPadded(5.f), EModulators::kMidiCC6, "CC6", TABLITSA2_STYLE.WithShowLabel(false), false), kCtrlTagMidiCC6Depth, "MidiMods")->Hide(true);*/

    pGraphics->AttachControl(new Tablitsa2VGroupControl("Midi Mods", "MidiMods", 5.f, 20.f, 5.f, 15.f));

    // Polyphony and Portamento
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(3, 1, 12, 6).GetVShifted(-15.f).GetVPadded(-2.f), kParamPortamentoMode, "", toggleStyle, "Rate", "Time"),
      kCtrlTagGlideMode, "Master")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool constantTime = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamNoteGlideTime, !constantTime);
        pGraphics->HideControl(kParamNoteGlideRate, constantTime);
        pGraphics->SetAllControlsDirty();
        });
    // Mono/Poly and Portamento Rate/Time switches
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(masterPanel.GetGridCell(1, 0, 4, 6).GetVShifted(-15.f).GetPadded(-5.f).GetCentredInside(dialSize), kParamNoteGlideRate, "Glide Rate"), kNoTag, "Master")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(masterPanel.GetGridCell(1, 0, 4, 6).GetVShifted(-15.f).GetPadded(-5.f).GetCentredInside(dialSize), kParamNoteGlideTime, "Glide Time"), kNoTag, "Master")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(4, 1, 12, 6).GetVShifted(-15.f).GetVPadded(-2.f), kParamMonophonic, "", toggleStyle, "Poly", "Mono"));
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(5, 1, 12, 6).GetVShifted(-15.f).GetVPadded(-2.f), kParamLegato, "", toggleStyle, "Retrigger", "Legato"));
    // Unison and Stereo Width
    const IRECT unisonPanel = masterPanel.GetReducedFromLeft(masterPanel.W() / 3.f).GetReducedFromRight(masterPanel.W() / 6.f).GetGridCell(0, 0, 2, 1).GetVShifted(-12.f).GetHPadded(-10.f).GetHShifted(2.f);
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(unisonPanel.GetGridCell(0, 0, 2, 3).GetHPadded(-10.f).GetHShifted(5.f), kParamUnisonVoices, "Voices", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(unisonPanel.GetGridCell(1, 0, 2, 2).GetCentredInside(dialSize), kParamUnisonDetune, "Detune", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(unisonPanel.FracRectVertical(0.5f, true).GetGridCell(0, 1, 1, 3).GetHPadded(-10.f).GetHShifted(-10.f), kParamUnisonChord, "Chord", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new Tablitsa2SliderControl(unisonPanel.FracRectVertical(0.5f, true).FracRectHorizontal(2.f / 3, true).GetGridCell(0, 1, 1, 3), kParamVibratoSpeed, "Rate"), kNoTag, "Unison");
    pGraphics->AttachControl(new Tablitsa2SliderControl(unisonPanel.FracRectVertical(0.5f, true).FracRectHorizontal(2.f / 3, true).GetGridCell(0, 2, 1, 3), kParamVibratoDepth, "Depth "), kNoTag, "Unison");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(unisonPanel.GetGridCell(1, 1, 2, 2).GetPadded(-5.f).GetCentredInside(dialSize), kParamStereoSpread, "Spread"), kNoTag, "Unison");
    pGraphics->AttachControl(new Tablitsa2VGroupControl("Unison", "Unison", 0.f, 13.f, 0.f, 0.f));
    pGraphics->AttachControl(new IVLEDMeterControl<2>(masterPanel.GetGridCell(0, 5, 2, 6).GetVShifted(-15.f).GetHShifted(-4.f).GetVPadded(-8.f)), kCtrlTagMeter, "Master");

    /* Effects Panel */
    const IRECT effectsPanel = IRECT(masterPanel.GetGridCell(1, 1, 2, 7).R, masterPanel.GetGridCell(1, 0, 2, 7).T, masterPanel.GetGridCell(1, 6, 2, 7).R, masterPanel.GetGridCell(1, 6, 2, 7).B).GetHPadded(-5.f).GetHShifted(2.f);
    const IRECT effectsControlsPanel = effectsPanel.GetReducedFromTop(40.f);

    IActionFunction voiceFxSlotChanged = [this, pGraphics](IControl* pControl) {
      int currentSlot = static_cast<int>(pControl->GetValue() * (TABLITSA2_MAX_VOICE_EFFECTS - 1));
      char str[2];
      sprintf(str, "%d", currentSlot + 1);
      dynamic_cast<IVSlideSwitchControl*>(pControl)->SetValueStr(str);

      // Set the controls' parameters to those of the current effect slot
      constexpr int numEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;
      dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob1))->SetParamIdx(kParamVoiceEffect1Param1 + currentSlot * numEffectParams);
      dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob2))->SetParamIdx(kParamVoiceEffect1Param2 + currentSlot * numEffectParams);
      dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob3))->SetParamIdx(kParamVoiceEffect1Param3 + currentSlot * numEffectParams);
      dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob4))->SetParamIdx(kParamVoiceEffect1Param4 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsToggle1)->SetParamIdx(kParamVoiceEffect1Param5 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsToggle2)->SetParamIdx(kParamVoiceEffect1Param6 + currentSlot * numEffectParams);

      for (int tag{ kCtrlTagVoiceEffectsKnob1 }; tag <= kCtrlTagVoiceEffectsKnob4; ++tag)
      {
        Tablitsa2IVModKnobControl* modKnob = dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(tag));
        if (modKnob->IsActive())
         modKnob->LoadModParams();
      }

      // Set the dropdown list's index to that of the effect which is currently active in this slot
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsList))->SetCurrentIndex(mVoiceEffectSlots[currentSlot], false);
      SwapVoiceEffectsUI(dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsList)), pGraphics, this, false); // Set parameters, keeping previous values
    };
    IActionFunction masterFxSlotChanged = [this, pGraphics](IControl* pControl) {
      int currentSlot = static_cast<int>(pControl->GetValue() * (TABLITSA2_MAX_MASTER_EFFECTS - 1));
      char str[2];
      sprintf(str, "%d", currentSlot + 1);
      dynamic_cast<IVSlideSwitchControl*>(pControl)->SetValueStr(str);

      // Set the controls' parameters to those of the current effect slot
      constexpr int numEffectParams = 6;
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob1)->SetParamIdx(kParamMasterEffect1Param1 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob2)->SetParamIdx(kParamMasterEffect1Param2 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob3)->SetParamIdx(kParamMasterEffect1Param3 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob4)->SetParamIdx(kParamMasterEffect1Param4 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsToggle1)->SetParamIdx(kParamMasterEffect1Param5 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsToggle2)->SetParamIdx(kParamMasterEffect1Param6 + currentSlot * numEffectParams);

      // Set the dropdown list's index to that of the effect which is currently active in this slot
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsList))->SetCurrentIndex(mMasterEffectSlots[currentSlot], false);
      SwapMasterEffectsUI(dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsList)), pGraphics, this, false); // Set parameters, keeping previous values
    };

    /* Effects Tab Control */
    pGraphics->AttachControl(new Tablitsa2EffectBankControl(effectsPanel, { "Voice FX", "Master FX", "Routing" }, { "E1", "E2", "E3" }), kCtrlTagEffectBank, "Effects")->SetActionFunction(
      [this, pGraphics](IControl* pControl) {
        dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsList))->Collapse();
        dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsList))->Collapse();
        mCurrentEffectsTab = dynamic_cast<Tablitsa2EffectBankControl*>(pControl)->GetActiveTabIdx();
      });
    
    /* Voice Effects */
    pGraphics->AttachControl(new ITextControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(60.f).GetHShifted(-55.f), "FX Slot:", TABLITSA2_TEXT), kNoTag, "E1");
    pGraphics->AttachControl(new IVSlideSwitchControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(40.f).GetReducedFromTop(22.f).GetHPadded(-20.f), voiceFxSlotChanged, "", toggleStyle, true, EDirection::Horizontal, TABLITSA2_MAX_VOICE_EFFECTS, 0), kCtrlTagVoiceEffectsSwitch, "E1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 0, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param1, kVEffect1Param1, "Param 1", modKnobStyle), kCtrlTagVoiceEffectsKnob1, "E1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param2, kVEffect1Param2, "Param 2", modKnobStyle), kCtrlTagVoiceEffectsKnob2, "E1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param3, kVEffect1Param3, "Param 3", modKnobStyle), kCtrlTagVoiceEffectsKnob3, "E1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param4, kVEffect1Param4, "Param 4", modKnobStyle), kCtrlTagVoiceEffectsKnob4, "E1");
    pGraphics->AttachControl(new Tablitsa2VToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 1, 1, 4).GetPadded(-10.f), kParamVoiceEffect1Param5, "", toggleStyle), kCtrlTagVoiceEffectsToggle1, "E1");
    pGraphics->AttachControl(new Tablitsa2VToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 2, 1, 4).GetPadded(-10.f), kParamVoiceEffect1Param6, "", toggleStyle), kCtrlTagVoiceEffectsToggle2, "E1");
    // Selector list
    pGraphics->AttachControl(new DropdownListControl(effectsPanel.FracRectHorizontal(0.4f, false).GetFromTop(50.f).GetReducedFromTop(20.f).GetPadded(-2.f), TABLITSA2_VOICE_EFFECTS_LIST), kCtrlTagVoiceEffectsList, "E1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SwapVoiceEffectsUI(pControl, pGraphics, this, true);
      });

    /* MASTER EFFECTS */
    pGraphics->AttachControl(new ITextControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(60.f).GetHShifted(-55.f), "FX Slot:", TABLITSA2_TEXT), kNoTag, "E2");
    pGraphics->AttachControl(new IVSlideSwitchControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(40.f).GetReducedFromTop(22.f).GetHPadded(-20.f), masterFxSlotChanged, "", toggleStyle, true, EDirection::Horizontal, 3, 0), kCtrlTagMasterEffectsSwitch, "E2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 0, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param1, "Left", modKnobStyle), kCtrlTagMasterEffectsKnob1, "E2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param2, "Right", modKnobStyle), kCtrlTagMasterEffectsKnob2, "E2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param3, "Feedback", modKnobStyle), kCtrlTagMasterEffectsKnob3, "E2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param4, "Mix", modKnobStyle), kCtrlTagMasterEffectsKnob4, "E2");
    pGraphics->AttachControl(new Tablitsa2VToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 1, 1, 4).GetPadded(-10.f), kParamMasterEffect1Param5, "", toggleStyle), kCtrlTagMasterEffectsToggle1, "E2");
    pGraphics->AttachControl(new Tablitsa2VToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 2, 1, 4).GetPadded(-10.f), kNoParameter, "", toggleStyle), kCtrlTagMasterEffectsToggle2, "E2");
    // Selector List
    pGraphics->AttachControl(new DropdownListControl(effectsPanel.FracRectHorizontal(0.4f, false).GetFromTop(50.f).GetReducedFromTop(20.f).GetPadded(-2.f), TABLITSA2_MASTER_EFFECTS_LIST), kCtrlTagMasterEffectsList, "E2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SwapMasterEffectsUI(pControl, pGraphics, this, true);
      });

    /* Routing */
    const IRECT routingPanel = effectsControlsPanel.GetPadded(-5.f, 10.f, -5.f, 0.f);
    pGraphics->AttachControl(new ISVGControl(routingPanel, routingImg), kNoTag, "E3");
    pGraphics->AttachControl(new Tablitsa2SliderControl(routingPanel.GetGridCell(0, 1, 5, 4).GetHPadded(-6.f).GetVShifted(-10.f), kParamOsc1FilterBypass, "Osc1 Filter BP", TABLITSA2_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");
    pGraphics->AttachControl(new Tablitsa2SliderControl(routingPanel.GetGridCell(1, 1, 5, 4).GetHPadded(-6.f).GetVShifted(-5.f), kParamOsc2FilterBypass, "Osc2 Filter BP", TABLITSA2_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");
    pGraphics->AttachControl(new Tablitsa2SliderControl(routingPanel.GetGridCell(3, 1, 5, 4).GetHShifted(-10.f), kParamOsc2EffectBypass, "Osc2 Effect BP", TABLITSA2_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");
    pGraphics->AttachControl(new Tablitsa2SliderControl(routingPanel.GetGridCell(4, 1, 5, 4).GetHShifted(-10.f).GetVShifted(2.f), kParamOsc1EffectBypass, "Osc1 Effect BP", TABLITSA2_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");

    /* LOWER-HALF CONTROLS */
    constexpr float modulePad[4]{-10.f, -5.f, -10.f, -15.f};

    /* Waveform Modulation Panel */
    const IRECT lowerPanel = globalControls.GetReducedFromTop(globalControls.H() * 0.5f).GetPadded(-10.f).GetReducedFromRight(250.f);
    const IRECT oscModPanel = lowerPanel.GetGridCell(0, 0, 3, 3);
    pGraphics->AttachControl(new ICaptionControl(oscModPanel.GetVPadded(-10.f).GetGridCell(0, 0, 3, 3).GetVShifted(10.f), kParamOscModulator, dropdownText), kNoTag, "OscMod")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        bool phaseMod = pControl->GetValue() < 0.5;
        // Update parameter index
        pGraphics->GetControlWithTag(kCtrlTagOscModFreq)->SetParamIdx(phaseMod ? kParamPhaseModFreq : kParamRingModFreq);
        pGraphics->GetControlWithTag(kCtrlTagOscModAmt)->SetParamIdx(phaseMod ? kParamPhaseModAmount : kParamRingModAmount);
        // Update modulation index
        dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagOscModFreq))->GetModulationFrom(phaseMod ? kVPhaseModFreq : kVRingModFreq);
        dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagOscModAmt))->GetModulationFrom(phaseMod ? kVPhaseModAmt : kVPhaseModAmt);
        IControl* pOsc1Switch = pGraphics->GetControlWithTag(kCtrlTagOsc1ModSwitch);
        pOsc1Switch->SetParamIdx(phaseMod ? kParamOsc1PM : kParamOsc1RM);
        IControl* pOsc2Switch = pGraphics->GetControlWithTag(kCtrlTagOsc2ModSwitch);
        pOsc2Switch->SetParamIdx(phaseMod ? kParamOsc2PM : kParamOsc2RM);
        // Update controls
        pControl->GetDelegate()->SendCurrentParamValuesFromDelegate();
        pOsc2Switch->SetDirty(false);
        pOsc1Switch->SetDirty(false);
      });
    pGraphics->AttachControl(new IVToggleControl(oscModPanel.GetVPadded(-10.f).GetVPadded(-10.f).GetGridCell(1, 0, 3, 3).GetMidHPadded(35.f), kParamOsc1PM, "Osc1", toggleStyle, "Osc 1 OFF", "Osc 1 ON"), kCtrlTagOsc1ModSwitch, "OscMod");
    pGraphics->AttachControl(new IVToggleControl(oscModPanel.GetVPadded(-10.f).GetVPadded(-10.f).GetGridCell(2, 0, 3, 3).GetMidHPadded(35.f), kParamOsc2PM, "Osc1", toggleStyle, "Osc 2 OFF", "Osc 2 ON"), kCtrlTagOsc2ModSwitch, "OscMod");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(oscModPanel.GetVPadded(-10.f).GetGridCell(0, 1, 1, 3).GetCentredInside(dialSize), kParamPhaseModFreq, kVPhaseModFreq, "Pitch"), kCtrlTagOscModFreq, "OscMod");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(oscModPanel.GetGridCell(0, 2, 1, 3).GetCentredInside(dialSize), kParamPhaseModAmount, kVPhaseModAmt, "Amount"), kCtrlTagOscModAmt, "OscMod");
    pGraphics->AttachControl(new Tablitsa2VGroupControl(oscModPanel.GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Oscillator Modulation"));
    
    /* Filter Panel */
    /* Filter 1 */
    const IRECT filter1Panel = lowerPanel.GetGridCell(1, 0, 3, 3).GetReducedFromRight(10.f);
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter1Type, dropdownText), kCtrlTagFilter1Type, "Flt1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        IControl* filterModeControl{ pGraphics->GetControlWithTag(kCtrlTagFilter1Mode) }; // Control containing list of filter modes for a given type
        int filterModeParam{ kParamFilter1ModeVSF + std::max(static_cast<int>(pControl->GetValue() * (kNumFilters - 1.) - 1), 0) };
        // Get range of possible modes
        filterModeControl->SetParamIdx(filterModeParam);
        filterModeControl->SetValue(GetParam(filterModeParam)->GetNormalized());
        // Update Controls
        if (pControl->GetValue() == kNoFilter)
          ResetFilterControls(pGraphics, 1);
        else if (pControl->GetValue() * ((int)kNumFilters - 1) == kComb)
        {
          filterModeControl->SetDisabled(true);
          // Swap Cutoff/Res/Drive for FF/FB/Delay
          pGraphics->HideControl(kParamFilter1Cutoff, true);
          pGraphics->HideControl(kParamFilter1FF, false);
          pGraphics->HideControl(kParamFilter1Resonance, true);
          pGraphics->HideControl(kParamFilter1FB, false);
          // Drive/delay switch
          pGraphics->HideControl(kParamFilter1Drive, true);
          pGraphics->HideControl(kParamFilter1Delay, false);
        }
        else
        {
          filterModeControl->SetDisabled(false);
          pGraphics->HideControl(kParamFilter1Cutoff, false);
          pGraphics->HideControl(kParamFilter1FF, true);
          pGraphics->HideControl(kParamFilter1Resonance, false);
          pGraphics->HideControl(kParamFilter1FB, true);
          // Drive/delay switch
          pGraphics->HideControl(kParamFilter1Drive, false);
          pGraphics->HideControl(kParamFilter1Delay, true);
        }
        filterModeControl->SetDirty(true);
      }
    );
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(1, 0, 3, 4).GetFromTop(20.f).GetVShifted(20.f), kParamFilter1ModeVSF, dropdownText), kCtrlTagFilter1Mode, "Flt1")->
      SetDisabled(static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagFilter1Type)->GetValue()) == kNoFilter);
    // Cutoff/Resonance/Drive
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter1Cutoff, kVFilter1Cutoff, "Cutoff"), kCtrlTagFilter1Cutoff, "Flt1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter1Resonance, kVFilter1Resonance, "Res."), kCtrlTagFilter1Resonance, "Flt1");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter1Drive, kVFilter1Drive, "Drive"), kCtrlTagFilter1Drive, "Flt1");
    pGraphics->AttachControl(new Tablitsa2SliderControl(filter1Panel.GetGridCell(0, 8, 1, 10).GetVPadded(-20.f), kParamFilter1Osc1Send, "O1", modKnobStyle.WithShowValue(false)), kCtrlTagFilter1Osc1, "Flt1");
    pGraphics->AttachControl(new Tablitsa2SliderControl(filter1Panel.GetGridCell(0, 9, 1, 10).GetVPadded(-20.f), kParamFilter1Osc2Send, "O2", modKnobStyle.WithShowValue(false)), kCtrlTagFilter1Osc2, "Flt1");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter1FF, kVFilter1FF, "F/F"), kCtrlTagFilter1FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter1FB, kVFilter1FB, "F/B"), kCtrlTagFilter1FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter1Delay, kVFilter1Delay, "Delay"), kCtrlTagFilter1Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new Tablitsa2VGroupControl(filter1Panel.GetPadded(modulePad[0], modulePad[1], modulePad[2] + 10.f, modulePad[3]), "Filter 1"));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FF))->GetModulationFrom(kVFilter1Cutoff);
    dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FB))->GetModulationFrom(kVFilter1Resonance);
    dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1Delay))->GetModulationFrom(kVFilter1Drive);

    /* Filter 2 */
    const IRECT filter2Panel = lowerPanel.GetGridCell(2, 0, 3, 3).GetReducedFromRight(10.f);
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter2Type, dropdownText), kCtrlTagFilter2Type, "Flt2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        IControl* filterModeControl{ pGraphics->GetControlWithTag(kCtrlTagFilter2Mode) };
        int filterModeParam{ kParamFilter2ModeVSF + std::max(static_cast<int>(pControl->GetValue() * (kNumFilters - 1.) - 1), 0) };
        filterModeControl->SetParamIdx(filterModeParam);
        filterModeControl->SetValue(GetParam(filterModeParam)->GetNormalized());
        if (pControl->GetValue() == kNoFilter)
          ResetFilterControls(pGraphics, 2);
        else if (pControl->GetValue() * ((int)kNumFilters - 1) == kComb)
        {
          filterModeControl->SetDisabled(true);
          // Swap Cutoff/Res/Drive for FF/FB/Delay
          pGraphics->HideControl(kParamFilter2Cutoff, true);
          pGraphics->HideControl(kParamFilter2FF, false);
          pGraphics->HideControl(kParamFilter2Resonance, true);
          pGraphics->HideControl(kParamFilter2FB, false);
          pGraphics->HideControl(kParamFilter2Drive, true);
          pGraphics->HideControl(kParamFilter2Delay, false);
        }
        else
        {
          filterModeControl->SetDisabled(false);
          pGraphics->HideControl(kParamFilter2Cutoff, false);
          pGraphics->HideControl(kParamFilter2FF, true);
          pGraphics->HideControl(kParamFilter2Resonance, false);
          pGraphics->HideControl(kParamFilter2FB, true);
          pGraphics->HideControl(kParamFilter2Drive, false);
          pGraphics->HideControl(kParamFilter2Delay, true);
        }
        filterModeControl->SetDirty(true);
      }
    );
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(1, 0, 3, 4).GetFromTop(20.f).GetVShifted(20.f), kParamFilter2ModeVSF, dropdownText), kCtrlTagFilter2Mode, "Flt2")->
      SetDisabled(static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagFilter2Type)->GetValue()) == kNoFilter);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter2Cutoff, kVFilter2Cutoff, "Cutoff"), kNoTag, "Flt2");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter2Resonance, kVFilter2Resonance, "Res."), kNoTag, "Flt2");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter2Drive, kVFilter2Drive, "Drive"), kNoTag, "Flt2");
    pGraphics->AttachControl(new Tablitsa2SliderControl(filter2Panel.GetGridCell(0, 8, 1, 10).GetVPadded(-20.f), kParamFilter2Osc1Send, "O1", modKnobStyle.WithShowValue(false)), kCtrlTagFilter2Osc1, "Flt2");
    pGraphics->AttachControl(new Tablitsa2SliderControl(filter2Panel.GetGridCell(0, 9, 1, 10).GetVPadded(-20.f), kParamFilter2Osc2Send, "O2", modKnobStyle.WithShowValue(false)), kCtrlTagFilter2Osc2, "Flt2");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter2FF, kVFilter2FF, "F/F"), kCtrlTagFilter2FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter2FB, kVFilter2FB, "F/B"), kCtrlTagFilter2FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter2Delay, kVFilter2Delay, "Delay"), kCtrlTagFilter2Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new Tablitsa2VGroupControl(filter2Panel.GetPadded(modulePad[0], modulePad[1], modulePad[2] + 10.f, modulePad[3]), "Filter 2"));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2FF))->GetModulationFrom(kVFilter2Cutoff);
    dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2FB))->GetModulationFrom(kVFilter2Resonance);
    dynamic_cast<Tablitsa2IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2Delay))->GetModulationFrom(kVFilter2Drive);

#pragma mark Modulation globalControls -

    // Envelope Sliders
    const IRECT envPanel = lowerPanel.GetGridCell(0, 1, 1, 3);
    constexpr float envDialShift{ 5.f };
    // Envelope 1
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 0, 3, 6).GetCentredInside(dialSize), kParamEnv1Attack, "Att", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 1, 3, 6).GetCentredInside(dialSize), kParamEnv1Decay, "Dec", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 1, 3, 6).GetCentredInside(dialSize), kParamEnv1DecayCurve, "Exp Dec", modKnobStyle), kNoTag, "Env1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 2, 3, 6).GetCentredInside(dialSize), kParamEnv1Sustain, kVEnv1Sustain, "Sus", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 3, 3, 6).GetCentredInside(dialSize), kParamEnv1Release, "Rel", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 3, 3, 6).GetCentredInside(dialSize), kParamEnv1ReleaseCurve, "Exp Rel", modKnobStyle), kNoTag, "Env1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2SliderControl(envPanel.GetGridCell(0, 4, 3, 6).GetVPadded(-15.f), kParamEnv1Velocity, "Velocity", TABLITSA2_STYLE, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "Env1");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(0, 5, 3, 6).GetVPadded(-15.f), EModulators::kEnvelope1, "Depth"), kCtrlTagEnv1Depth, "Env1");
    pGraphics->AttachControl(new Tablitsa2VGroupControl(envPanel.GetGridCell(0, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Env 1"));
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 1, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv1Decay, editCurve);
        pGraphics->HideControl(kParamEnv1DecayCurve, !editCurve);
        });
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetGridCell(0, 3, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv1Release, editCurve);
        pGraphics->HideControl(kParamEnv1ReleaseCurve, !editCurve);
        });
    // Envelope 2
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 0, 3, 6).GetCentredInside(dialSize), kParamEnv2Attack, "Att", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 1, 3, 6).GetCentredInside(dialSize), kParamEnv2Decay, "Dec", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 1, 3, 6).GetCentredInside(dialSize), kParamEnv2DecayCurve, "Exp Dec", modKnobStyle), kNoTag, "Env2")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 2, 3, 6).GetCentredInside(dialSize), kParamEnv2Sustain, kVEnv2Sustain, "Sus", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 3, 3, 6).GetCentredInside(dialSize), kParamEnv2Release, "Rel", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 3, 3, 6).GetCentredInside(dialSize), kParamEnv2ReleaseCurve, "Exp Rel", modKnobStyle), kNoTag, "Env2")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2SliderControl(envPanel.GetGridCell(1, 4, 3, 6).GetVPadded(-15.f), kParamEnv2Velocity, "Velocity", TABLITSA2_STYLE, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "Env2");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(1, 5, 3, 6).GetVPadded(-15.f), EModulators::kEnvelope2, "Depth"), kCtrlTagEnv2Depth, "Env2");
    pGraphics->AttachControl(new Tablitsa2VGroupControl(envPanel.GetGridCell(1, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Env 2"));
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 1, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv2Decay, editCurve);
        pGraphics->HideControl(kParamEnv2DecayCurve, !editCurve);
        });
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetGridCell(1, 3, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv2Release, editCurve);
        pGraphics->HideControl(kParamEnv2ReleaseCurve, !editCurve);
        });
    // Amplitude Envelope
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 0, 3, 6).GetCentredInside(dialSize), kParamAmpEnvAttack, "Att", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 1, 3, 6).GetCentredInside(dialSize), kParamAmpEnvDecay, "Dec", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 1, 3, 6).GetCentredInside(dialSize), kParamAmpEnvDecayCurve, "Dec", modKnobStyle), kNoTag, "AmpEnv")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 2, 3, 6).GetCentredInside(dialSize), kParamAmpEnvSustain, kVAmpEnvSustain, "Sus", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 3, 3, 6).GetCentredInside(dialSize), kParamAmpEnvRelease, "Rel", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 3, 3, 6).GetCentredInside(dialSize), kParamAmpEnvReleaseCurve, "Exp Rel", modKnobStyle), kNoTag, "AmpEnv")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2SliderControl(envPanel.GetGridCell(2, 4, 3, 6).GetVPadded(-15.f), kParamAmpEnvVelocity, "Velocity", TABLITSA2_STYLE, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(2, 5, 3, 6).GetVPadded(-15.f), EModulators::kAmpEnvelope, "Depth"), kCtrlTagAmpEnvDepth, "AmpEnv");
    pGraphics->AttachControl(new Tablitsa2VGroupControl(envPanel.GetGridCell(2, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Amp Env"));
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 1, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamAmpEnvDecay, editCurve);
        pGraphics->HideControl(kParamAmpEnvDecayCurve, !editCurve);
        });
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetGridCell(2, 3, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamAmpEnvRelease, editCurve);
        pGraphics->HideControl(kParamAmpEnvReleaseCurve, !editCurve);
        });

    /* LFO Panels */
    const IRECT lfoPanel = lowerPanel.GetGridCell(0, 2, 1, 3);
    // LFO1
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateHz, kVLFO1RateHz, "Rate"), kNoTag, "LFO1");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateTempo, "Rate"), kNoTag, "LFO1")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(lfoPanel.GetGridCell(0, 2, 3, 5).GetCentredInside(dialSize), kParamLFO1Amp, kVLFO1Amp, "Amp"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(0, 1, 6, 5).GetVPadded(-20.f), kParamLFO1Shape, dropdownText), kNoTag, "LFO1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SendArbitraryMsgFromUI(kMsgUpdateLFO1Plot);
      }
    );
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(1, 1, 9, 5).GetVPadded(-6.f), kParamLFO1RateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagLFO1RateMode, "LFO1")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
      bool sync = pControl->GetValue() > 0.5;
      pGraphics->HideControl(kParamLFO1RateHz, sync);
      pGraphics->HideControl(kParamLFO1RateTempo, !sync);
      //pGraphics->SetAllControlsDirty();
    });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(2, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO1Restart, "Restart", toggleStyle, "Free Run", "Retrigger"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(0, 4, 3, 5).GetVPadded(-18.f), EModulators::kLFO1, "Depth"), kCtrlTagLFO1Depth, "LFO1");
    pGraphics->AttachControl(new Tablitsa2ModPlotControl(lfoPanel.GetGridCell(0, 3, 3, 5).GetVPadded(-20.f), kParamLFO1Phase, nullptr, 1024, 25), kCtrlTagLFO1Plot, "LFO1");
    pGraphics->AttachControl(new Tablitsa2VGroupControl(lfoPanel.GetGridCell(0, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "LFO 1"));

    // LFO 2
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateHz, kVLFO2RateHz, "Rate"), kNoTag, "LFO2");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateTempo, "Rate"), kNoTag, "LFO2")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(lfoPanel.GetGridCell(1, 2, 3, 5).GetCentredInside(dialSize), kParamLFO2Amp, kVLFO2Amp, "Amp"), kNoTag, "LFO2");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(2, 1, 6, 5).GetVPadded(-20.f), kParamLFO2Shape, dropdownText), kNoTag, "LFO2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SendArbitraryMsgFromUI(kMsgUpdateLFO2Plot);
      }
    );
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(4, 1, 9, 5).GetVPadded(-6.f), kParamLFO2RateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagLFO2RateMode, "LFO2")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamLFO2RateHz, sync);
        pGraphics->HideControl(kParamLFO2RateTempo, !sync);
        //pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(5, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO2Restart, "Restart", toggleStyle, "Free Run", "Retrigger"), kNoTag, "LFO2");
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(1, 4, 3, 5).GetVPadded(-18.f), EModulators::kLFO2, "Depth"), kCtrlTagLFO2Depth, "LFO2");
    pGraphics->AttachControl(new Tablitsa2ModPlotControl(lfoPanel.GetGridCell(1, 3, 3, 5).GetVPadded(-20.f), nullptr, 1024, 25), kCtrlTagLFO2Plot, "LFO1");
    pGraphics->AttachControl(new Tablitsa2VGroupControl(lfoPanel.GetGridCell(1, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "LFO 2"));

    /* Sequencer */
    const IRECT sequencerDisplayPanel = lowerPanel.GetGridCell(2, 2, 3, 3).GetFromLeft(PLUG_WIDTH - lfoPanel.L).GetVShifted(9.f).GetPadded(10.f, 12.f, 0.f, 0.f);
    const IRECT sequencerParamsPanel = IRECT(lfoPanel.R, lfoPanel.T, globalControls.R, lowerPanel.B - sequencerDisplayPanel.H()).GetPadded(-5.f, -5.f, -5.f, -10.f);

    // Sequencer Action Functions
    auto stepRandomizeAction = [pGraphics](IControl* pCaller) {
      SplashClickActionFunc(pCaller);
      pGraphics->GetDelegate()->SendArbitraryMsgFromUI(kMsgRandomizeSequencer);
    };

    pGraphics->AttachControl(new Tablitsa2IVKnobControl(sequencerParamsPanel.GetGridCell(0, 0, 2, 3).GetCentredInside(dialSize), kParamSequencerSteps, "Steps"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(sequencerParamsPanel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamSequencerRateHz, kVSequencerRateHz, "Rate"), kNoTag, "SequencerParams")->Hide(true);
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(sequencerParamsPanel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamSequencerRateTempo, "Rate"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new Tablitsa2IVModKnobControl(sequencerParamsPanel.GetGridCell(0, 2, 2, 3).GetCentredInside(dialSize), kParamSequencerAmp, kVSequencerAmp, "Amp"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new Tablitsa2IVKnobControl(sequencerParamsPanel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamSequencerGlide, "Glide"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(4, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kParamSequencerRateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagSequencerRateMode, "Sequencer")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamSequencerRateTempo, !sync);
        pGraphics->HideControl(kParamSequencerRateHz, sync);
        pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new SequencerControl<kNumSeqSteps>(sequencerDisplayPanel.GetReducedFromRight(dialSize).GetPadded(-20.f, -10.f, -20.f, -16.f), ""),
      kCtrlTagSequencer, "Sequencer")->SetActionFunction([pGraphics](IControl* pCaller) {
      // Get all values
      /* NB: In the SnapToMouse() function, values for sliders which fall between MouseMove events are interpolated using the preceeding and following slider.
      If this function were to save only the last slider hit (as in the iPlugChunks example at the time of this comment), sliders with interpolated values will not
      have their values stored and serialized. */
      for(int i(0); i < pCaller->NVals(); ++i)
      {
        // Bundle the slider index with the slider value
        ParamTuple data = { i, pCaller->GetValue(i) };
        // Send the new value to be evaluated in the OnMessage function
        pGraphics->GetDelegate()->SendArbitraryMsgFromUI(kMsgSeqSliderChanged, kCtrlTagSequencer, sizeof(ParamTuple), &data);
      }
      });
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(5, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kParamSequencerRestart, "", toggleStyle, "Free Run", "Restart"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new ModSliderControl(sequencerDisplayPanel.GetGridCell(0, 5, 1, 6).GetVPadded(-18.f).GetVShifted(-5.f).GetMidHPadded(20.f), EModulators::kSequencer, "Depth"), kCtrlTagSequencerDepth, "Sequencer");
    pGraphics->AttachControl(new IVRadioButtonControl(sequencerParamsPanel.SubRectVertical(2, 1).SubRectHorizontal(3, 1).GetReducedFromBottom(sequencerParamsPanel.H() / 7.f).GetHPadded(-15.f), kParamSequencerStepMode, { "Auto", "Gate", "Trig." }, "Step Mode", TABLITSA2_STYLE.WithLabelText(TABLITSA2_TEXT.WithSize(13.f)).WithColor(EVColor::kX1, IColor(150, 155, 155, 155))), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(6, 1, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kNoParameter, "", toggleStyle, "Cont.", "Quant."), kCtrlTagSequencerQuant, "SequencerParams")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        mSequencerIsQuantized = pControl->GetValue();
        int nSteps = mSequencerIsQuantized > 0.5 ? 12 : 0;
        dynamic_cast<IVTrackControlBase*>(pGraphics->GetControlWithTag(kCtrlTagSequencer))->SetNSteps(nSteps);
      });
    pGraphics->AttachControl(new IVButtonControl(sequencerParamsPanel.GetGridCell(6, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), stepRandomizeAction, "Randomize", pushButtonStyle.WithLabelText(pushButtonStyle.labelText.WithSize(12.f))), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new Tablitsa2VGroupControl("Sequencer", "SequencerParams", 2.f, 15.f, 2.f, 5.f));

    
#ifndef VST3_API
    pGraphics->AttachControl(new IVButtonControl(b.GetFromBottom(KEYBOARD_HEIGHT + 15.f).GetFromRight(170.f).GetFromTop(30.).GetVPadded(-3.f), SplashClickActionFunc,
      "Show/Hide Keyboard", TABLITSA2_STYLE.WithColor(kFG, COLOR_WHITE).WithLabelText({15.f, EVAlign::Middle})))->SetAnimationEndActionFunction(
      [pGraphics, globalControls](IControl* pCaller) {
        static bool hide = false;
        pGraphics->GetControlWithTag(kCtrlTagKeyboard)->Hide(hide = !hide);
        pGraphics->Resize((int)PLUG_WIDTH, hide ? (int)PLUG_HEIGHT : (float)PLUG_HEIGHT + (float)KEYBOARD_HEIGHT, pGraphics->GetDrawScale());
    });
    pGraphics->SetQwertyMidiKeyHandlerFunc([pGraphics](const IMidiMsg& msg) {
                                              dynamic_cast<IVKeyboardControl*>(pGraphics->GetControlWithTag(kCtrlTagKeyboard))->SetNoteFromMidi(msg.NoteNumber(), msg.StatusMsg() == IMidiMsg::kNoteOn);
                                           });
#endif
#ifdef OS_IOS
    if(!IsAuv3AppExtension())
    {
      pGraphics->AttachControl(new IVButtonControl(b.GetFromTRHC(100, 100), [pGraphics](IControl* pCaller) {
                               dynamic_cast<IGraphicsIOS*>(pGraphics)->LaunchBluetoothMidiDialog(pCaller->GetRECT().L, pCaller->GetRECT().MH());
                               SplashClickActionFunc(pCaller);
                             }, "BTMIDI"));
    }
#endif

    pGraphics->AttachPopupMenuControl();

#ifndef VST3_API
    pGraphics->Resize(PLUG_WIDTH, PLUG_HEIGHT + (int)KEYBOARD_HEIGHT, 0.8, true);
#else
    pGraphics->Resize(PLUG_WIDTH, PLUG_HEIGHT, 0.8, true);
#endif
  };
#endif
}

bool Tablitsa2::SerializeState(IByteChunk& chunk) const
{
  bool savedOK = true;
  int version{ GetPluginVersion(false) };
  savedOK &= chunk.Put(&version) > 0;
  savedOK &= chunk.Put(&mPresetID) > 0;
  for (int i{ 0 }; i < PRESET_NAME_CHAR_LENGTH; ++i)
    savedOK &= chunk.Put(mPresetName + i) > 0;
  
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    savedOK &= chunk.Put(&mDSP.mSeqSteps[i]) > 0;
  }
  savedOK &= chunk.Put(&mDSP.mLoadedWavetables[0]) > 0;
  savedOK &= chunk.Put(&mDSP.mLoadedWavetables[1]) > 0;
  savedOK &= chunk.Put(&mSequencerIsQuantized) > 0;

  // Save effect panel states
  for (int i{ 0 }; i < TABLITSA2_MAX_VOICE_EFFECTS; ++i)
  {
    savedOK &= chunk.Put(&mVoiceEffectSlots[i]) > 0;
  }
  for (int i{ 0 }; i < TABLITSA2_MAX_MASTER_EFFECTS; ++i)
  {
    savedOK &= chunk.Put(&mMasterEffectSlots[i]) > 0;
    savedOK &= chunk.Put(&mDelayTempoSync[i]) > 0;
  }
  for (int i{0}; i < kNumParams; ++i)
  {
    for (int j{0}; j < kNumModulators; ++j)
    {
      savedOK &= chunk.Put(&mModMatrix[i][j]);
    }
  }
  savedOK &= SerializeParams(chunk);
  return savedOK; // must remember to call SerializeParams at the end
}

// this over-ridden method is called when the host is trying to load the plug-in state and you need to unpack the data into your algorithm
int Tablitsa2::UnserializeState(const IByteChunk& chunk, int startPos)
{
  startPos = CheckVersion(chunk);
  startPos = chunk.Get(&mPresetID, startPos);
  for (int i{ 0 }; i < PRESET_NAME_CHAR_LENGTH; ++i)
    startPos = chunk.Get(mPresetName + i, startPos);

  double v = 0.;
  // unserialize the steps state before unserializing the regular params
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    startPos = chunk.Get(&v, startPos);
    mDSP.mSeqSteps[i] = v;
  }
  // Wavetable 1
  startPos = chunk.Get(&v, startPos);
  mDSP.mLoadedWavetables[0] = v;
  // Wavetable 2
  startPos = chunk.Get(&v, startPos);
  mDSP.mLoadedWavetables[1] = v;

  // Notify user if improper wavetable values are loaded
  if (mDSP.mLoadedWavetables[0] - static_cast<int>(mDSP.mLoadedWavetables[0]) > 0.01 || mDSP.mLoadedWavetables[1] - static_cast<int>(mDSP.mLoadedWavetables[1]) > 0.01)
  {
    mDSP.mLoadedWavetables[0] = 1.;
    mDSP.mLoadedWavetables[1] = 2.;
    bool abort{ ShowLoadErrorMessageBox() };
    if (abort)
    {
      if (mStateBackup.Size() == 0)
      {
        LoadDefaultState();
        return -1;
      }
      else
        return UnserializeState(mStateBackup, 0);
    }
      
  }
  if (!LoadWavetables())
    return -1;

  // Sequencer quantization
  startPos = chunk.Get(&v, startPos);
  mSequencerIsQuantized = v;
  
  // Load effect panel states
  for (int i{ 0 }; i < TABLITSA2_MAX_VOICE_EFFECTS; ++i)
  {
    startPos = chunk.Get(&mVoiceEffectSlots[i], startPos);

    switch (mVoiceEffectSlots[i])
    {
    case kDistortionEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 1", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 2", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 3", 0., TABLITSA2_DISTORTION_FREQ_LOW, TABLITSA2_DISTORTION_FREQ_HIGH, 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 4", 0., 0., 1., 0.01);
      break;
    }
    case kSampleAndHoldEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Rate", 10., TABLITSA2_SAH_MIN_MS, TABLITSA2_SAH_MAX_MS, 0.01, "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Decay", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Noise", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Mix", 0., 0., 1., 0.01);
      break;
    }
    case kTexturizerEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Cutoff", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Drive", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Res. Freq.", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Res. Amt.", 0., 0., 1., 0.01);
      break;
    }
    case kCMEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Mod Amount", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Pitch", 0., -2., 2., 0.01, "8v.");
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Cutoff", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Mix", 0., 0., 1., 0.01);
      break;
    }
    }
    // Send message to change effect BEFORE parameters are loaded, so that the effect receives the proper parameter values
    SendArbitraryMsgFromUI(kMsgVoiceEffect1Changed + i, kNoTag, sizeof(mVoiceEffectSlots[i]), reinterpret_cast<void*>(&mVoiceEffectSlots[i]));
  }
  for (int i{ 0 }; i < TABLITSA2_MAX_MASTER_EFFECTS; ++i)
  {
    // Get the effect ID for each slot
    startPos = chunk.Get(&mMasterEffectSlots[i], startPos);
    startPos = chunk.Get(&mDelayTempoSync[i], startPos);
    // Give the params the proper ranges
    switch (mMasterEffectSlots[i])
    {
    case kDelayEffect:
    {
      if (mDelayTempoSync[i])
      {
        GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitEnum("Delay L (beats)", DelayEffect<sample>::k8th, { DELAY_TEMPODIV_VALIST });
        GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitEnum("Delay R (beats)", DelayEffect<sample>::k8th, { DELAY_TEMPODIV_VALIST });
      }
      else
      {
        GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Delay L (ms)", 100., 1., TABLITSA2_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
        GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Delay R (ms)", 100., 1., TABLITSA2_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
      }
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Delay Feedback", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Delay Mix", 0., 0., 1., 0.01);
      break;
    }
    case kEQEffect:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("EQ Low Gain ", 1., 0., 2., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("EQ Mid Gain", 1., 0., 2., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("EQ Mid Freq", 0.5, 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("EQ High Gain", 1., 0., 2., 0.01);
      break;
    }
    case kReverbEffect:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 1", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 2", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 3", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 4", 0., 0., 1., 0.01);
      break;
    }
    case kReverb2Effect:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 1", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 2", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 3", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 4", 0., 0., 1., 0.01);
      break;
    }
    default:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 1", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 2", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 3", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 4", 0., 0., 1., 0.01);
      break;
    }
    }
    SendArbitraryMsgFromUI(kMsgMasterEffect1Changed + i, kNoTag, sizeof(mMasterEffectSlots[i]), reinterpret_cast<void*>(&mMasterEffectSlots[i]));

    // Tempo Sync needs to be set before time values are sent to the DSP thread, but after the effect is loaded
    ENTER_PARAMS_MUTEX;
    GetParam(kParamMasterEffect1Param5 + kNumMasterEffectParams * i)->Set(mDelayTempoSync[i] ? 1. : 0.);
    OnParamChange(kParamMasterEffect1Param5 + kNumMasterEffectParams * i);
    LEAVE_PARAMS_MUTEX;
  }
  for (int i{0}; i < kNumParams; ++i)
  {
    for (int j{0}; j < kNumModulators; ++j)
    {
      startPos = chunk.Get(&mModMatrix[i][j], startPos);
    }
  }

  // must remember to call UnserializeParams at the end
  return UnserializeParams(chunk, startPos);
}

void Tablitsa2::LoadDefaultState()
{
  // unserialize the steps state before unserializing the regular params
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    mDSP.mSeqSteps[i] = 0.;
  }
  // Wavetable 1
  mDSP.mLoadedWavetables[0] = 1.;
  // Wavetable 2
  mDSP.mLoadedWavetables[1] = 2.;
  LoadWavetables();
  // Sequencer quantization
  mSequencerIsQuantized = 0.;

  for (int i{ 0 }; i < TABLITSA2_MAX_VOICE_EFFECTS; ++i)
  {
    SetVoiceFXSlot(i, kNoVoiceEffect);
    SendArbitraryMsgFromUI(kMsgVoiceEffect1Changed + i, kNoTag, sizeof(mVoiceEffectSlots[i]), reinterpret_cast<void*>(&mVoiceEffectSlots[i]));
  }
  for (int i{ 0 }; i < TABLITSA2_MAX_MASTER_EFFECTS; ++i)
  {
    SetMasterFXSlot(i, kNoMasterEffect);
    SendArbitraryMsgFromUI(kMsgMasterEffect1Changed + i, kNoTag, sizeof(mMasterEffectSlots[i]), reinterpret_cast<void*>(&mMasterEffectSlots[i]));
  }
  // First reset all mod sliders to the null position
  for (int i{0}; i < kNumModSliders; ++i)
  {
    if (GetUI())
    {
      GetUI()->GetControlWithTag(kModSliders[i])->SetValue(0.5);
    }
  }
  
  memset(mModMatrix[0], (double)0.0, sizeof(double) * kNumParams * kNumModulators);

  DefaultParamValues();
  OnParamReset(kDelegate);

  UpdateUIControls();
  dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCurrentIndex(std::find(PRESET_ID_LIST, PRESET_ID_LIST + N_PRESETS, "Default") - PRESET_ID_LIST, false);
}

int Tablitsa2::CheckVersion(const IByteChunk& presetData)
{
  int presetVersion{};
  int bytePos = presetData.Get(&presetVersion, 0);
  if (presetVersion == GetPluginVersion(false))
    return bytePos;
  else
    return 0;
}

void Tablitsa2::OnUIOpen()
{
  UpdateUIControls();
}

void Tablitsa2::UpdateUIControls()
{
  //SendCurrentParamValuesFromDelegate();
  dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCustomStr(mPresetName);
  // Sequencer
  auto* pSequencer = GetUI()->GetControlWithTag(kCtrlTagSequencer);

  for (int i = 0; i < kNumSeqSteps; i++)
  {
    pSequencer->SetValue(mDSP.mSeqSteps[i], i);
  }

  // Wavetables
  auto* pPeriodicTable = dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable));
  pPeriodicTable->SetSelectedElement(static_cast<int>(mDSP.mLoadedWavetables[0]), 0);
  pPeriodicTable->SetSelectedElement(static_cast<int>(mDSP.mLoadedWavetables[1]), 1);
  GetUI()->GetControlWithTag(kCtrlTagSequencerQuant)->SetValue(mSequencerIsQuantized);

  // LFOPlots
  dynamic_cast<Tablitsa2ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO1Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO1.GetCurrentTable()));
  dynamic_cast<Tablitsa2ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO2Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO2.GetCurrentTable()));

  SendCurrentParamValuesFromDelegate();
  GetUI()->SetAllControlsDirty();

  for (auto i{ 0 }; i < std::size(kStartupTriggerControls); ++i)
  {
    IControl* pControl = GetUI()->GetControlWithTag(kStartupTriggerControls[i]);
    if (!pControl->IsHidden())
      pControl->SetDirty(true);
  }

  // Set up master effect controls
  DropdownListControl* effectList = dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagMasterEffectsList));
  // Set up effect slots (without reseting the parameters!)
  for (int i{ TABLITSA2_MAX_MASTER_EFFECTS - 1}; i >= 0; --i)
  {
    int effectIdx = mMasterEffectSlots[i];
    effectList->SetCurrentIndex(effectIdx);
  }
  effectList->SetCurrentIndex(static_cast<int>(GetUI()->GetControlWithTag(kCtrlTagMasterEffectsSwitch)->GetValue() * (TABLITSA2_MAX_MASTER_EFFECTS - 1)));
  IControl* slotSwitch = GetUI()->GetControlWithTag(kCtrlTagMasterEffectsSwitch);
  slotSwitch->SetValueToDefault();
  slotSwitch->SetDirty(true);

  // Set up voice effect controls
  effectList = dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagVoiceEffectsList));
  // Set up effect slots (without reseting the parameters!)
  for (int i{ TABLITSA2_MAX_VOICE_EFFECTS - 1 }; i >= 0; --i)
  {
    int effectIdx = mVoiceEffectSlots[i];
    effectList->SetCurrentIndex(effectIdx);
  }
  effectList->SetCurrentIndex(static_cast<int>(GetUI()->GetControlWithTag(kCtrlTagVoiceEffectsSwitch)->GetValue() * (TABLITSA2_MAX_VOICE_EFFECTS - 1)));
  slotSwitch = GetUI()->GetControlWithTag(kCtrlTagVoiceEffectsSwitch);
  slotSwitch->SetValueToDefault();
  slotSwitch->SetDirty(true);
}

/* implement this and return true to trigger your custom about box, when someone clicks about in the menu of a standalone app or VST3 plugin */
bool Tablitsa2::OnHostRequestingAboutBox() { return false; } // See IPlugAPP_dialog.cpp

/* implement this and return true to trigger your custom help info, when someone clicks help in the menu of a standalone app or VST3 plugin */
bool Tablitsa2::OnHostRequestingProductHelp() { return false; }

int Tablitsa2::GetActiveModIdx() const
{
  return mActiveModIdx;
}

void Tablitsa2::RefreshEffectBankControl()
{
  dynamic_cast<Tablitsa2EffectBankControl*>(GetUI()->GetControlWithTag(kCtrlTagEffectBank))->TabChanged(mCurrentEffectsTab, false); // Prevent controls in hidden panels from being shown
}

#if IPLUG_DSP
void Tablitsa2::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
  mDSP.ProcessBlock(nullptr, outputs, 2, nFrames, mTimeInfo.mPPQPos, mTimeInfo.mTransportIsRunning, mTimeInfo.mTempo);
  mMeterSender.ProcessBlock(outputs, nFrames, kCtrlTagMeter);
//  mLFO1VisSender.PushData({ kCtrlTagLFO1Vis, {float(mDSP.mGlobalLFO1.GetLastOutput())}});
//  mLFO2VisSender.PushData({ kCtrlTagLFO2Vis, {float(Tablitsa2DSP<sample>::mGlobalLFO2.GetLastOutput())} });
}

void Tablitsa2::OnIdle()
{
  mMeterSender.TransmitData(*this);
//  mLFO1VisSender.TransmitData(*this);
//  mLFO2VisSender.TransmitData(*this);

  // Send time signature to DSP class
  GetTimeSig(mDSP.mTSNum, mDSP.mTSDenom);
  mDSP.mTempo = GetTempo();

  // Only allow the sequencer to retrigger gates in mono mode
  if (mDSP.mActiveSequencer && mDSP.mMono && mDSP.mGlobalSequencer.GetMode() == Sequencer<sample>::kTriggerGate)
  {
    bool retrigger = mDSP.mActiveSequencer->GetTrigger();
    if (retrigger)
    {
      mDSP.ForEachVoice([](Tablitsa2DSP<sample>::Voice& voice) {
        if (voice.GetBusy() && !voice.JustTriggered())
          voice.Retrigger();
        });
    }
  }

  if (GetUI())
  {
    // Update Sequencer Position
    if (mDSP.GetSequencerStep() != mDSP.mSeqPos)
    {
      mDSP.mSeqPos = mDSP.GetSequencerStep();
      SendControlMsgFromDelegate(kCtrlTagSequencer, IVMultiSliderControl<>::kMsgTagSetHighlight, sizeof(int), &mDSP.mSeqPos);
    }
  }
}

void Tablitsa2::OnReset()
{
  mDSP.Reset(GetSampleRate(), GetBlockSize());
}

void Tablitsa2::ProcessMidiMsg(const IMidiMsg& msg)
{
  TRACE;
  
  int status = msg.StatusMsg();
  
  switch (status)
  {
    case IMidiMsg::kNoteOn:
      if (mDSP.mTableLoading) // Don't accept new notes when swapping wavetables
        return;
    case IMidiMsg::kNoteOff:
    case IMidiMsg::kPolyAftertouch:
    case IMidiMsg::kControlChange:
    case IMidiMsg::kProgramChange:
    case IMidiMsg::kChannelAftertouch:
    case IMidiMsg::kPitchWheel:
    {
      goto handle;
    }
    default:
      return;
  }
  
handle:
  mDSP.ProcessMidiMsg(msg);
  SendMidiMsg(msg);

  // Send additional messages for unison voices
  mDSP.ResetDetune();
  for (int i{ 1 }; i < mDSP.mUnisonVoices; ++i)
  {
    IMidiMsg copy_msg = IMidiMsg(msg);
    if (status == kNoteOnAction)
      copy_msg.MakeNoteOnMsg(msg.NoteNumber(), msg.Velocity(), msg.mOffset + 1, msg.Channel());
    else if (status == kNoteOffAction)
      copy_msg.MakeNoteOffMsg(msg.NoteNumber(), msg.mOffset + 1, msg.Channel());
    mDSP.ProcessMidiMsg(msg);
  }
}

void Tablitsa2::OnParamChange(int paramIdx)
{
  ENTER_PARAMS_MUTEX
  mDSP.SetParam(paramIdx, GetParam(paramIdx)->Value());
  LEAVE_PARAMS_MUTEX
}

bool Tablitsa2::OnMessage(int msgTag, int ctrlTag, int dataSize, const void* pData)
{
  // TODO: switch statement
  if(ctrlTag == kCtrlTagBender && msgTag == IWheelControl::kMessageTagSetPitchBendRange)
  {
    const int bendRange = *static_cast<const int*>(pData);
    mDSP.mSynth.SetPitchBendRange(bendRange);
  }
  else
  {
    switch (msgTag)
    {
    case kMsgSavePreset:
    {
      auto* backup = reinterpret_cast<const bool*>(pData);
      std::string filePath;

      bool useSystemDialog{ (bool)LOAD_WITH_DIALOG };
      if (useSystemDialog)
      {
#pragma mark - Windows
        WDL_String fileStr("NewPreset");
        WDL_String pathStr(GetDataPath("\\Tablitsa2\\presets\\").c_str());
        GetUI()->PromptForFile(fileStr, pathStr, EFileAction::Save, "TPST");
        // Get and store preset name
        std::string path = fileStr.Get();
        size_t nameStart{ path.find_last_of("\\") + 1 };
        size_t nameEnd{ path.find_last_of(".") };
        strcpy(mPresetName, path.substr(nameStart, nameEnd - nameStart).c_str());
        // Abort if no filename was submitted
        filePath = fileStr.Get();
        if (filePath == "\0")
          break;
      }
      else
      {
        filePath = GetDataPath("\\Tablitsa2\\presets\\") + "TEST.TPST";
      }
      IByteChunk newPreset;
      SerializeState(newPreset);
      SavePreset(newPreset, filePath.c_str(), *backup);

      break;
    }
    case kMsgLoadPreset:
    {
      auto* backup = reinterpret_cast<const bool*>(pData);
      std::string filePath;
      if (*backup)
      {
        filePath = GetDataPath("\\Tablitsa2\\presets\\") + "_tablitsa.bkp";
      }
      else if ((bool)LOAD_WITH_DIALOG)
      {
        WDL_String fileStr("UserPreset");
        WDL_String pathStr(GetDataPath("\\Tablitsa2\\presets\\").c_str());
        GetUI()->PromptForFile(fileStr, pathStr, EFileAction::Open, "TPST");
        filePath = fileStr.Get();
        if (filePath == "\0")
          break;
      }
      else
        filePath = GetDataPath("\\Tablitsa2\\presets\\") + "TEST.TPST";

      IByteChunk presetData = LoadPreset(filePath.c_str(), *backup);
      if (presetData.Size() == 0)
        break;
      UnserializeState(presetData, 0);
      UpdateUIControls();
      break;
    }
    case kMsgWavetable1Changed:
    {
      PeriodicTable* pTableControl = dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable));
      pTableControl->SetTableLoading(true, 0);
      const double wtIdx = *static_cast<const double*>(pData); //Semi-Normalized value (between 1/118 and 1)
      mDSP.mLoadedWavetables[0] = wtIdx * 118.; // Corresponds to the atomic number (not the 0-indexed vector item)
      try
      {
        mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[0]) - 1, 0);
      }
      catch (...)
      {
        GetUI()->ShowMessageBox("Wavetable file not found. Check AppData directory.", "Error", EMsgBoxType::kMB_OK);
      }
      pTableControl->SetTableLoading(false, 0);
      break;
    }
    case kMsgWavetable2Changed:
    {
      PeriodicTable* pTableControl = dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable));
      pTableControl->SetTableLoading(true, 1);
      const double wtIdx = *static_cast<const double*>(pData);
      mDSP.mLoadedWavetables[1] = wtIdx * 118.;
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[1]) - 1, 1);
      pTableControl->SetTableLoading(false, 1);
      break;
    }
    case kMsgUpdateLFO1Plot:
    {
      dynamic_cast<Tablitsa2ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO1Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO1.GetCurrentTable()));
      break;
    }
    case kMsgUpdateLFO2Plot:
    {
      dynamic_cast<Tablitsa2ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO2Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO2.GetCurrentTable()));
      break;
    }
    case kMsgSeqSliderChanged:
    {
      // Read incoming data as a tuple
      auto* pConverted = reinterpret_cast<const ParamTuple*>(pData);
      mDSP.mSeqSteps[pConverted->idx] = pConverted->value;
      break;
    }
    case kMsgRandomizeSequencer:
    {
      auto* pSequencer = dynamic_cast<SequencerControl<kNumSeqSteps>*>(GetUI()->GetControlWithTag(kCtrlTagSequencer));
      for (auto i{ 0 }; i < kNumSeqSteps; ++i)
      {
        mDSP.mSeqSteps[i] = static_cast<double>(std::rand() % 1000) / 1000.;
        pSequencer->SetValue(mDSP.mSeqSteps[i], i);
      }
      pSequencer->SetDirty(true);
      break;
    }
    case kMsgVoiceEffect1Changed:
    case kMsgVoiceEffect2Changed:
    case kMsgVoiceEffect3Changed:
    {
      //std::lock_guard<std::mutex> lg(mDSP.mProcMutex);
      int effectIdx = *reinterpret_cast<const int*>(pData);
      switch (effectIdx)
      {
      case kDistortionEffect:
        mDSP.ForEachVoice([msgTag, this](Tablitsa2DSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kDistortionEffect);
          });
        break;
      case kSampleAndHoldEffect:
        mDSP.ForEachVoice([msgTag, this](Tablitsa2DSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kSampleAndHoldEffect);
          });
        break;
      case kTexturizerEffect:
        mDSP.ForEachVoice([msgTag, this](Tablitsa2DSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kTexturizerEffect);
          });
        break;
      case kCMEffect:
        mDSP.ForEachVoice([msgTag, this](Tablitsa2DSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kCMEffect);
          });
        break;
      default:
        mDSP.ForEachVoice([msgTag, this](Tablitsa2DSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kNoVoiceEffect);
          });
        break;
      }
      SendCurrentParamValuesFromDelegate();
      break;
    }
    case kMsgMasterEffect1Changed:
    case kMsgMasterEffect2Changed:
    case kMsgMasterEffect3Changed:
    {
      std::lock_guard<std::mutex> lg(mDSP.mEffectMutex);
      int effectIdx = *reinterpret_cast<const int*>(pData);
      delete mDSP.mEffects[msgTag - kMsgMasterEffect1Changed];
      switch (effectIdx)
      {
      case kChorusEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new ChorusEffect<iplug::sample>(GetSampleRate());
        break;
      case kDelayEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new DelayEffect<iplug::sample>(GetSampleRate(), &mDSP.mGlobalMetronome);
        break;
      case kEQEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new EQ3Effect<iplug::sample>(GetSampleRate());
        break;
      case kReverbEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new ReverbEffect<iplug::sample>(GetSampleRate());
        break;
      case kReverb2Effect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new Reverb2Effect<iplug::sample>(GetSampleRate());
        break;
      default:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new Effect<iplug::sample>(GetSampleRate());
        break;
      }
      SendCurrentParamValuesFromDelegate();
      break;
    }
    default:
      break;
    }
  }
  return false;
}

void Tablitsa2::SetActiveModIdx(int idx)
{
  int prev = mActiveModIdx;
  mActiveModIdx = idx;
  if (idx > -1)
    dynamic_cast<Tablitsa2IVModKnobControl*>(GetUI()->GetControl(idx))->ColorSwap();
  if (prev > -1)
    dynamic_cast<Tablitsa2IVModKnobControl*>(GetUI()->GetControl(prev))->ColorSwap();
}


double Tablitsa2::GetModMatrixSlot(const int parameter, const int modulator) { return mModMatrix[parameter][modulator]; }

/*
@param value The modulator depth, between -1.0 and 1.0
*/
void Tablitsa2::SetModMatrixSlot(const int parameter, const int modulator, double value) {
  value = Clip(value, -1.0, 1.0);
  bool targetIsModulator{parameter > kNumVoiceModulations};

  mModMatrix[parameter][modulator] = value;
  if (targetIsModulator)
  {
    const int metaModIdx{modulator - kNumVoiceModulations - 1};
    mDSP.ForEachVoice([parameter, modulator, value](Tablitsa2DSP<sample>::Voice& v) { v.UpdateVoiceModulatorParam(parameter, modulator, value); });
  }
  else
    mDSP.ForEachVoice([parameter, modulator, value](Tablitsa2DSP<sample>::Voice& v) { v.UpdateVoiceParam(parameter, modulator, value); });
}

bool Tablitsa2::LoadWavetables()
{
  for (int i{ 0 }; i < 2; ++i)
  {
    try
    {
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[i]) - 1, i);
    }
    catch (const std::out_of_range& e)
    {
      mDSP.mLoadedWavetables[i] = static_cast<double>(i + 1);
      dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable))->SetSelectedElement(static_cast<int>(mDSP.mLoadedWavetables[i]), i);
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[i]) - 1, i);

      bool abort{ ShowLoadErrorMessageBox() };
      if (abort)
      {
        LoadDefaultState();
        return false;
      }
    }
  }
  return true;
}

void Tablitsa2::SavePreset(IByteChunk& byteData, const char* filename, bool isBackup)
{
  assert(byteData.Size() != 0);
  if (isBackup)
  {
    mStateBackup = byteData;
    return;
  }

  std::ofstream presetFile(filename, std::ios::binary | std::ios::trunc);

  if (!presetFile.is_open())
  {
    GetUI()->ShowMessageBox("Cannot access selected file path.", "Error", EMsgBoxType::kMB_OK);
    presetFile.close();
    return;
  }

  //char* allBytes = new char[byteData.Size()];
  std::vector<char> allBytes(byteData.Size());
  byteData.GetBytes(&allBytes[0], byteData.Size(), 0);
  presetFile.write(&allBytes[0], byteData.Size());
  presetFile.close();
  //delete[] allBytes;
}

IByteChunk Tablitsa2::LoadPreset(const char* filename, bool isBackup)
{
  if (isBackup)
    return mStateBackup;
  else
  {
    mStateBackup.Clear();
    SerializeState(mStateBackup);
  }

  IByteChunk presetData;

  std::vector<char> fileData = ReadAllBytes(filename);

  if (fileData.size() == 0)
  {
    GetUI()->ShowMessageBox("Preset file not found or invalid format", "Error", EMsgBoxType::kMB_OK);
    return presetData;
  }
  
  presetData.PutBytes(&fileData[0], fileData.size());

  return presetData;
}

bool Tablitsa2::ShowLoadErrorMessageBox()
{
  bool abort;
  GetUI()->ShowMessageBox("Invalid wavetable codes. Preset file may be corrupted or incompatible with the plugin version. Continue loading preset?", "Warning", EMsgBoxType::kMB_YESNO,
    [&abort](EMsgBoxResult result) {
      abort = (result == EMsgBoxResult::kNO);
    });
  return abort;
}
#endif

std::string GetDataPath(char* appendPath)
{
  USES_CONVERSION;
  TCHAR szPath[MAX_PATH];
  if (SUCCEEDED(SHGetFolderPath(NULL, CSIDL_APPDATA, NULL, 0, szPath)))
  {
    PathAppend(szPath, A2T(appendPath));
  }
#if VST3_API
  // In the standalone app builds, `TCHAR` is type `char`. In VST3 builds, it's `wchar_t`, requiring a different conversion method.
  std::wstring wpath(szPath);
  std::string path = std::wstring_convert<std::codecvt_utf8<wchar_t>>().to_bytes(wpath);
#else
  std::string path = szPath;
#endif
  return path;
}

std::vector<char> ReadAllBytes(const char* fname)
{
  std::ifstream ifs(fname, std::ios::binary | std::ios::ate);
  std::ifstream::pos_type pos = ifs.tellg();

  if (!ifs.is_open())
  {
    return std::vector<char>(0);
  }

  std::vector<char> allBytes(pos);
  if (pos == 0)
    return allBytes;
  ifs.seekg(0, std::ios::beg);
  ifs.read(&allBytes[0], pos);
  ifs.close();
  return allBytes;
}
