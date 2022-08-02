#pragma once

#include "radiofarmerDSP.h"

#include "PeriodicTable.h"
#include "Tablitsa2Effects.h"
#include "Tablitsa2Oscillators.h"
#include "Modulation.h"

#include "IPlugConstants.h"
#include "Oscillator.h"
#include "MidiSynth.h"
#include "ADSREnvelope.h"
#include "Smoothers.h"
#include "LFO.h"

#include "Oversampler.h"

#ifdef VECTOR
 #define FRAME_INTERVAL OUTPUT_SIZE * 2
#else
  #define FRAME_INTERVAL 1
#endif

#if _DEBUG
#define EFFECT_OS_FACTOR 1
#else
#define EFFECT_OS_FACTOR 4
#endif

// Minimum dB value for filter sends, to be set to -inf dB (i.e. 0.)
#define SEND_DB_FLOOR -36.

constexpr double kMaxEnvTimeScalar = 0.5;

using namespace iplug;
using namespace radiofarmer;

/*
Global Modulations (smoothers): These values are computed once per sample and sent to all voices 
*/
enum EModulations
{
  kModGainSmoother = 0,
  kModEnv1SustainSmoother,
  kModEnv2SustainSmoother,
  kModAmpEnvSustainSmoother,
  kModLFO1RateHzSmoother,
  kModLFO1AmpSmoother,
  kModLFO2RateHzSmoother,
  kModLFO2AmpSmoother,
  kModSeqRateHzSmoother,
  kModSeqAmpSmoother,
  kModPanSmoother,
  kModWavetable1PitchSmoother,
  kModWavetable1PosSmoother,
  kModWavetable1BendSmoother,
  kModWavetable1FormantSmoother,
  kModWavetable1AmpSmoother,
  kModWavetable2PitchSmoother,
  kModWavetable2PosSmoother,
  kModWavetable2BendSmoother,
  kModWavetable2FormantSmoother,
  kModWavetable2AmpSmoother,
  kModFilter1CutoffSmoother,
  kModFilter1ResonanceSmoother,
  kModFilter1DriveSmoother,
  kModFilter1CombFFSmoother,
  kModFilter1CombFBSmoother,
  kModFilter1CombDelaySmoother,
  kModFilter2CutoffSmoother,
  kModFilter2ResonanceSmoother,
  kModFilter2DriveSmoother,
  kModFilter2CombFFSmoother,
  kModFilter2CombFBSmoother,
  kModFilter2CombDelaySmoother,
  kModOsc1PhaseModFreqSmoother,
  kModOsc2PhaseModFreqSmoother,
  kModOsc1PhaseModAmtSmoother,
  kModOsc2PhaseModAmtSmoother,
  kModOsc1RingModFreqSmoother,
  kModOsc2RingModFreqSmoother,
  kModOsc1RingModAmtSmoother,
  kModOsc2RingModAmtSmoother,
  kModVoiceEffect1Param1Smoother,
  kModVoiceEffect1Param2Smoother,
  kModVoiceEffect1Param3Smoother,
  kModVoiceEffect1Param4Smoother,
  kModVoiceEffect2Param1Smoother,
  kModVoiceEffect2Param2Smoother,
  kModVoiceEffect2Param3Smoother,
  kModVoiceEffect2Param4Smoother,
  kModVoiceEffect3Param1Smoother,
  kModVoiceEffect3Param2Smoother,
  kModVoiceEffect3Param3Smoother,
  kModVoiceEffect3Param4Smoother,
  kNumModulations,
};

enum EOscModulators
{
  kSine = 0,
  kOsc1,
  kOsc2,
  kNumOscModulators
};

// See Modulators.h for the enumeration of the number of modulators

#define UNISON_CHORD_LIST "None", "8va", "M7", "D7", "D7 6/5", "D7 4/3", "m7", "ø7", "dim7"

enum EUnisonChords
{
  kNoChord = -1,
  kOctaves = 0,
  kMaj7,
  kDom7,
  kDom7FirstInv,
  kDom7SecondInv,
  kMin7,
  kHalfDim7,
  kDim7
};

struct ChannelParam
{
  double L{ 0. };
  double R{ 0. };
  double LR[2]{ 0. };

  void Set(double l, double r)
  {
    LR[0] = L = l;
    LR[1] = R = r;
  }

  double& operator [](int ch)
  {
    return LR[ch];
  }
};

struct UnisonVoiceManager
{
  int mMaxVoices;

  // Detuning
  double mMaxDetune; // octaves
  int mNVoices{ 1 };
  int mDetuneVoiceIdx{ 0 };
  int mChord{ EUnisonChords::kNoChord };
  double* mDetuneBuf;

  // Panning
  int mPanVoiceIdx{ 0 };
  double mMaxPan; // [-1., 1.]
  ChannelParam* mPanBuf;

  const double mUnisonInvervals[8][5]{
    {0., 1., 0., 2., 0.}, // Octaves
    {0., 7. / 12, 4. / 12, 11. / 12, 1.}, // Major 7
    {0., 7. / 12, 4. / 12, 10. / 12, 1.}, // Dominant 7
    {0., 7. / 12, 4. / 12, -2. / 12, 1.}, // Dominant 7 1st Inversion
    {0., -5. / 12, -2. / 12, 4. / 12, 1.}, // Dominant 7 2nd Inversion
    {0., 7. / 12, 3. / 12, 10. / 12, 1.}, // Minor 7
    {0., 3. / 12., 6. / 12, 10. / 12, 1.}, // Half-Diminished 7
    {0., 3. / 12., 6. / 12, 9. / 12, 1.} // Full-Diminished 7
  };

  UnisonVoiceManager(int maxVoices, double maxDetuneSemitones=1.) : mMaxVoices(maxVoices), mMaxDetune(maxDetuneSemitones)
  {
    mDetuneBuf = new double[mMaxVoices] {0.};
    mPanBuf = new ChannelParam[mMaxVoices]{};
    mPanBuf[0].Set(1., 1.);
  }

  void SetNVoices(int nVoices)
  {
    mNVoices = std::min(nVoices, mMaxVoices);
    RefillBuffers();
  }

  void SetMaxDetune(double maxDetuneSemitones, bool refillBuffer=true)
  {
    mMaxDetune = maxDetuneSemitones / 12;
    if (refillBuffer)
      RefillDetuneBuffer();
  }

  void SetMaxPan(double maxPanDegrees, bool refillBuffer=true)
  {
    mMaxPan = maxPanDegrees / 180.;
    if (refillBuffer)
      RefillPanBuffer();
  }

  void Reset()
  {
    mDetuneVoiceIdx = 0;
    mPanVoiceIdx = 0;
  }

  void SetChord(int chord)
  {
    mChord = chord;
  }

  void RefillBuffers()
  {
    RefillDetuneBuffer();
    RefillPanBuffer();
  }

  inline void RefillDetuneBuffer()
  {
    int nVoices = std::max(2, mNVoices);

    // Pitch detuning
    if (mChord == EUnisonChords::kNoChord)
    {
      for (auto i{ 0 }; i < mNVoices; ++i)
        mDetuneBuf[i] = mMaxDetune * static_cast<double>(-i % mNVoices) / (static_cast<double>(nVoices) - 1.) * -1;
    }
    else
    {
      for (auto i{ 0 }; i < mNVoices; ++i)
        mDetuneBuf[i] = mUnisonInvervals[mChord][i % 5] + mMaxDetune * (static_cast<double>(std::rand() % 100) / 50 - 1.);
    }
  }

  inline void RefillPanBuffer()
  {
    constexpr double sqrt2{ 1.45 };
    // Panning (in progress)
    double totalPan = std::abs(mMaxPan);
    if (mNVoices > 1)
    {
      double pan = totalPan;
      mPanBuf[0].Set(
        (1. + std::copysign(pan, mMaxPan)) / sqrt2,
        (1. - std::copysign(pan, mMaxPan)) / sqrt2
      );
    }
    else
      mPanBuf[0].Set(1., 1.);
    for (auto i{ 1 }; i < mNVoices; ++i)
    {
      double pan;
      if (i % 2)
        pan = totalPan / i; // Odd-numbered voices: Pan by the same magnitude as the last voice, but in the opposite direction
      else
        pan = totalPan / (i + 1); // Even-numbered voices: Pan within the maximum pan range, to ever-smaller extents
      mPanBuf[i].Set(
        (1. - std::copysign(pan, mMaxPan)) / sqrt2,
        (1. + std::copysign(pan, mMaxPan)) / sqrt2
      );
    }
  }

  double DetuneNext()
  {
    int next = mDetuneVoiceIdx++ % mNVoices;
    mDetuneVoiceIdx %= mNVoices;
    return mDetuneBuf[next];
  }

  void SetNextPan(double* channelPan)
  {
    int next = mPanVoiceIdx++ % mNVoices;
    mPanVoiceIdx %= mNVoices;
    channelPan[0] = mPanBuf[next][0];
    channelPan[1] = mPanBuf[next][1];
  }

  ~UnisonVoiceManager()
  {
    delete[] mDetuneBuf;
  }
};

template<typename T>
class Tablitsa2DSP
{
public:
#pragma mark - Voice
  class Voice : public SynthVoice
  {
  public:
    Voice(std::vector<std::string> programList, int id=0) : Voice()
    {
    }

    Voice(Tablitsa2DSP<T>* master, int id = 0)
      : mMaster(master), mID(id),
      mAmpEnv("gain", [&]() { mOsc1.Reset(); }),
      mEnv1("env1", [&]() { mOsc1.Reset(); }),
      mLFO1(&GetMaster()->mGlobalMetronome),
      mLFO2(&GetMaster()->mGlobalMetronome),
      mSequencer(&GetMaster()->mGlobalMetronome, GetMaster()->mSeqSteps) // capture ok on RT thread?
    {
//      DBGMSG("new Voice: %i control inputs.\n", static_cast<int>(mInputs.size()));

      // Meta-Modulation Parameters:
      // First package the `ParameterModulator`s for the modulatable parameters of modulators themselves the into vectors...
      std::vector<ParameterModulator<kNumModulators>*> env1Mods{ &mVoiceMetaModParams[kVEnv1Sustain - kNumVoiceModulations - 1] };
      std::vector<ParameterModulator<kNumModulators>*> env2Mods{ &mVoiceMetaModParams[kVEnv2Sustain - kNumVoiceModulations - 1] };
      std::vector<ParameterModulator<kNumModulators>*> ampEnvMods{&mVoiceMetaModParams[kVAmpEnvSustain - kNumVoiceModulations - 1]};
      std::vector<ParameterModulator<kNumModulators>*> lfo1Mods{&mVoiceMetaModParams[kVLFO1RateHz - kNumVoiceModulations - 1], &mVoiceMetaModParams[kVLFO1Amp - kNumVoiceModulations - 1]};
      std::vector<ParameterModulator<kNumModulators>*> lfo2Mods{&mVoiceMetaModParams[kVLFO2RateHz - kNumVoiceModulations - 1], &mVoiceMetaModParams[kVLFO2Amp - kNumVoiceModulations - 1]};
      std::vector<ParameterModulator<kNumModulators>*> seqMods{&mVoiceMetaModParams[kVSequencerRateHz - kNumVoiceModulations - 1], &mVoiceMetaModParams[kVSequencerAmp - kNumVoiceModulations - 1]};

      // ...then add the `ParameterModulator` vectors to the modulator list
      mModulators.AddModulator(&mEnv1, env1Mods);
      mModulators.AddModulator(&mEnv2, env2Mods);
      mModulators.AddModulator(&mAmpEnv, ampEnvMods);
      mModulators.AddModulator(&mLFO1, lfo1Mods);
      mModulators.AddModulator(&mLFO2, lfo2Mods);
      mModulators.AddModulator(&mSequencer, seqMods);
      mAmpEnv.Kill(true); // Force amplitude envelopes to start in the "Idle" stage

      // Fill the envelope queues for legato mode with null pointers
      GetMaster()->AmpEnvQueue.push_back(nullptr);
      GetMaster()->Env1Queue.push_back(nullptr);
      GetMaster()->Env2Queue.push_back(nullptr);
    }

    bool GetBusy() const override
    {
      return mAmpEnv.GetBusy();
    }

    UnisonVoiceManager& GetDetuner()
    {
      return GetMaster()->mDetuner;
    }

    void ResetUnisonParams()
    {
      mDetune = GetDetuner().DetuneNext();
      GetDetuner().SetNextPan(mPan);
    }

    const bool JustTriggered()
    {
      bool wasJustTriggered = mTriggered;
      mTriggered = false;
      return wasJustTriggered;
    }

    void Retrigger()
    {
      Trigger(mVelocity, true);
    }

    void Trigger(double level, bool isRetrigger) override
    {
      mVelocity = level; // TODO: Handling of different velocity settings (i.e. which envelopes are affected by velocity)
      mTriggerRand = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);

      ResetUnisonParams();

      for (auto f : mFilters)
      {
        f->Reset();
      }

      // Reset LFOs and Sequencer
      if (mLFO1Restart)
      {
        mModulators.ReplaceModulator(&mLFO1, 0);
        mLFO1.Reset();
      }
      else
      {
        mModulators.ReplaceModulator(dynamic_cast<FastLFO<T>*>(&(GetMaster()->mGlobalLFO1)), 0);
      }

      if (mLFO2Restart)
      {
        mModulators.ReplaceModulator(&mLFO2, 1);
        mLFO2.Reset();
      }
      else
      {
        mModulators.ReplaceModulator(dynamic_cast<FastLFO<T>*>(&GetMaster()->mGlobalLFO2), 1);
      }
      
      if (mSequencerRestart)
      {
        mModulators.ReplaceModulator(dynamic_cast<Sequencer<T>*>(&mSequencer), 0);
        if (!isRetrigger) // Only reset the sequencer if this is not a retrigger call (i.e. not evoked by the sequencer gate itself)
          mSequencer.Reset();
        // Update sequencer display with this voice's phase
        GetMaster()->mActiveSequencer = &mSequencer;
      }
      else
      {
        mModulators.ReplaceModulator(&GetMaster()->mGlobalSequencer, 0);
        // Update sequencer display with this voice's phase
        GetMaster()->mActiveSequencer = &GetMaster()->mGlobalSequencer;
      }


      if (!mLegato)
      {
        double velSubtr = 1. - level;
        mAmpEnv.Start(1. - velSubtr * mAmpEnvVelocityMod, 1. - mAmpEnvVelocityMod * kMaxEnvTimeScalar * level);
        mEnv1.Start(1. - velSubtr * mEnv1VelocityMod, 1. - mEnv1VelocityMod * kMaxEnvTimeScalar * level);
        mEnv2.Start(1. - velSubtr * mEnv2VelocityMod, 1. - mEnv2VelocityMod * kMaxEnvTimeScalar * level);
        // Don't reset oscillators in legato mode - the phase sync will case clicks
        mOsc1.Reset();
        mOsc2.Reset();
      }
      else
      {
        // Check for currently active envelopes
        Envelope<T>* masterAmpEnv = nullptr;
        Envelope<T>* masterEnv1 = nullptr;
        Envelope<T>* masterEnv2 = nullptr;
        for (auto i{ 0 }; i < GetMaster()->AmpEnvQueue.size(); ++i)
        {
          if (GetMaster()->AmpEnvQueue[i])
          {
            masterAmpEnv = GetMaster()->AmpEnvQueue[i];
            masterEnv1 = GetMaster()->Env1Queue[i];
            masterEnv2 = GetMaster()->Env2Queue[i];
          }
        }
        // If active envelopes were found, sync this voice's envelopes to them
        if (masterAmpEnv)
        {
          double velSubtr = 1. - level;
          mAmpEnv.StartAt(1. - velSubtr * mAmpEnvVelocityMod,
            masterAmpEnv->GetValue(), masterAmpEnv->GetPrevResult(), masterAmpEnv->GetStage(),
            1. - mAmpEnvVelocityMod * kMaxEnvTimeScalar * level);
          mEnv1.StartAt(1. - velSubtr * mEnv1VelocityMod,
            masterEnv1->GetValue(), masterEnv1->GetPrevResult(), masterEnv1->GetStage(),
            1. - mEnv1VelocityMod * kMaxEnvTimeScalar * level);
          mEnv2.StartAt(1. - velSubtr * mEnv2VelocityMod,
            masterEnv2->GetValue(), masterEnv2->GetPrevResult(), masterEnv2->GetStage(),
            1. - mEnv2VelocityMod * kMaxEnvTimeScalar * level);

        }
        else
        {
          double velSubtr = 1. - level;
          mAmpEnv.Start(1. - velSubtr * mAmpEnvVelocityMod, 1. - mAmpEnvVelocityMod * kMaxEnvTimeScalar * level);
          mEnv1.Start(1. - velSubtr * mEnv1VelocityMod, 1. - mEnv1VelocityMod * kMaxEnvTimeScalar * level);
          mEnv2.Start(1. - velSubtr * mEnv2VelocityMod, 1. - mEnv2VelocityMod * kMaxEnvTimeScalar * level);
        }
        // Sync the master envelopes to this voice's envelopes
        GetMaster()->Env1Queue[mID] = &mEnv1;
        GetMaster()->Env2Queue[mID] = &mEnv2;
        GetMaster()->AmpEnvQueue[mID] = &mAmpEnv;
      }

      if (!isRetrigger)
        mTriggered = true;
    }
    
    void Release() override
    {
      mAmpEnv.Release();
      mEnv1.Release();
      mEnv2.Release();
      // Remove this voice's envelopes from the envelope queue
      GetMaster()->AmpEnvQueue[mID] = nullptr;
      GetMaster()->Env1Queue[mID] = nullptr;
      GetMaster()->Env2Queue[mID] = nullptr;
    }

    void ProcessSamplesAccumulating(T** inputs, T** outputs, int nInputs, int nOutputs, int startIdx, int nFrames) override
    {
       /* // inputs to the synthesizer can just fetch a value every block, like this:
      double gate = mInputs[kVoiceControlGate].endValue;
      // or write the entire control ramp to a buffer, like this, to get sample-accurate ramps:
      mInputs[kVoiceControlTimbre].Write(mTimbreBuffer.Get(), startIdx, nFrames); */

      double pitchBend = mInputs[kVoiceControlPitchBend].endValue;
      double pitch = mInputs[kVoiceControlPitch].endValue + pitchBend + mDetune; // MIDI Pitch = (MidiKey - 69) / 12

      // Aftertouch and mod wheel
      double aftertouch = mInputs[kVoiceControlPressure].endValue;
      double modWheel = mInputs[kVoiceControlTimbre].endValue;

      // Set the static (note-constant) modulator values
      T keytrack = (pitch * 12. + 69.) / 128.;
      T staticMods[]{ mVelocity, keytrack, mTriggerRand };
      // Write ramps for modulators
      mVoiceMetaModParams[kVEnv1Sustain - kNumVoiceModulations - 1].SetInitialValue(inputs[kModEnv1SustainSmoother][0]);
      mVoiceMetaModParams[kVEnv2Sustain - kNumVoiceModulations - 1].SetInitialValue(inputs[kModEnv2SustainSmoother][0]);
      mVoiceMetaModParams[kVAmpEnvSustain - kNumVoiceModulations - 1].SetInitialValue(inputs[kModAmpEnvSustainSmoother][0]);
      mVoiceMetaModParams[kVLFO1Amp - kNumVoiceModulations - 1].SetInitialValue(inputs[kModLFO1AmpSmoother][0]);
      mVoiceMetaModParams[kVLFO1RateHz - kNumVoiceModulations - 1].SetInitialValue(inputs[kModLFO1RateHzSmoother][0]);
      mVoiceMetaModParams[kVLFO2Amp - kNumVoiceModulations - 1].SetInitialValue(inputs[kModLFO2AmpSmoother][0]);
      mVoiceMetaModParams[kVLFO2RateHz - kNumVoiceModulations - 1].SetInitialValue(inputs[kModLFO2RateHzSmoother][0]);
      mVoiceMetaModParams[kVSequencerAmp - kNumVoiceModulations - 1].SetInitialValue(inputs[kModSeqAmpSmoother][0]);
      mVoiceMetaModParams[kVSequencerRateHz - kNumVoiceModulations - 1].SetInitialValue(inputs[kModSeqRateHzSmoother][0]);

      // Apply modulation to modulators and calculate their buffer values
      mModulators.MetaProcessBlock_Fast(&(inputs[kModEnv1SustainSmoother]), nFrames);
      // Add MIDI Control Ramp buffers
      mInputs[kVoiceControlGate].Write(mModulators.GetList()[kVelocity], startIdx, nFrames);
      mInputs[kVoiceControlPressure].Write(mModulators.GetList()[kAftertouch], startIdx, nFrames);
      mInputs[kVoiceControlTimbre].Write(mModulators.GetList()[kModWheel], startIdx, nFrames);

      // Apply modulation ramps to all modulated parameters
#ifdef VECTOR
      mVoiceModParams.ProcessBlockVec4d(&inputs[kModPanSmoother], mModulators.GetList(), mVModulations.GetList(), nFrames);
#else
      mVoiceModParams.ProcessBlock(&inputs[kModPanSmoother], mModulators.GetList(), mVModulations.GetList(), nFrames);
#endif

      const double phaseModFreqFact1 = pow(2., mVModulations.GetList()[kVOsc1PhaseModFreq][0] / 12.);
      const double ringModFreqFact1 = pow(2., mVModulations.GetList()[kVOsc1RingModFreq][0] / 12.);
      const double phaseModFreqFact2 = pow(2., mVModulations.GetList()[kVOsc2PhaseModFreq][0] / 12.);
      const double ringModFreqFact2 = pow(2., mVModulations.GetList()[kVOsc2RingModFreq][0] / 12.);

      // make sound output for each output channel
      for (auto i = startIdx; i < startIdx + nFrames; i += FRAME_INTERVAL)
      {
        int bufferIdx = i - startIdx;
        //        float noise = mTimbreBuffer.Get()[i] * Rand();
        T ampEnvVal{ mModulators.GetList()[2][bufferIdx] }; // Calculated for easy access

        // Oscillator Parameters
        // Osc1
        T osc1Freq = 440. * pow(2., pitch + mVModulations.GetList()[kVWavetable1PitchOffset][bufferIdx] / 12. + mMaster->mVibratoDepth * mMaster->mVibratoOsc.Process());
        mOsc1.SetWtPosition(1. - mVModulations.GetList()[kVWavetable1Position][bufferIdx]); // Wavetable 1 Position
        mOsc1.SetWtBend(mVModulations.GetList()[kVWavetable1Bend][bufferIdx]); // Wavetable 1 Bend
        mOsc1.SetFormant(mVModulations.GetList()[kVWavetable1Formant][bufferIdx], mOsc1FormantOn);
        T osc1Amp = mVModulations.GetList()[kVWavetable1Amp][bufferIdx];

        // Osc2
        T osc2Freq = 440. * pow(2., pitch + mVModulations.GetList()[kVWavetable2PitchOffset][bufferIdx] / 12. + mMaster->mVibratoDepth * mMaster->mVibratoOsc.Process());
        mOsc2.SetWtPosition(1. - mVModulations.GetList()[kVWavetable2Position][bufferIdx]); // Wavetable 2 Position
        mOsc2.SetWtBend(mVModulations.GetList()[kVWavetable2Bend][bufferIdx]); // Wavetable 2 Bend
        mOsc2.SetFormant(mVModulations.GetList()[kVWavetable2Formant][bufferIdx], mOsc2FormantOn);
        T osc2Amp = mVModulations.GetList()[kVWavetable2Amp][bufferIdx];

        // Filters
        mFilters.at(0)->SetCutoff(mVModulations.GetList()[kVFilter1Cutoff][bufferIdx]); // Filter 1 Cutoff
        mFilters.at(0)->SetQ(mVModulations.GetList()[kVFilter1Resonance][bufferIdx]); // Filter 1 Resonance
        mFilters.at(0)->SetDrive(mVModulations.GetList()[kVFilter1Drive][bufferIdx]); // Filter 1 Drive

        mFilters.at(1)->SetCutoff(mVModulations.GetList()[kVFilter2Cutoff][bufferIdx]); // Filter 2 Cutoff
        mFilters.at(1)->SetQ(mVModulations.GetList()[kVFilter2Resonance][bufferIdx]); // Filter 2 Resonance
        mFilters.at(1)->SetDrive(mVModulations.GetList()[kVFilter2Drive][bufferIdx]); // Filter 2 Drive

        // Phase and Ring Modulators (freq. only set once per block)
        mOsc1.SetPhaseModulation(mVModulations.GetList()[kVOsc1PhaseModAmt][bufferIdx], osc1Freq * phaseModFreqFact1);
        mOsc2.SetPhaseModulation(mVModulations.GetList()[kVOsc2PhaseModAmt][bufferIdx], osc2Freq * phaseModFreqFact2);
        mOsc1.SetRingModulation(mVModulations.GetList()[kVOsc1RingModAmt][bufferIdx], osc1Freq * ringModFreqFact1);
        mOsc2.SetRingModulation(mVModulations.GetList()[kVOsc2RingModAmt][bufferIdx], osc2Freq * ringModFreqFact2);

        // Signal Processing
#ifdef VECTOR
        // Oscillators
        Vec4d osc1_v = mOsc1.ProcessMultiple(osc1Freq) * osc1Amp;// * ampEnvVal;
        Vec4d osc2_v = mOsc2.ProcessMultiple(osc2Freq) * osc2Amp;// * ampEnvVal;
        osc1_v.store(mOscOutputs.GetList()[0] + bufferIdx);
        osc2_v.store(mOscOutputs.GetList()[1] + bufferIdx);

        // Panning
        T panMod = mVModulations.GetList()[kVPan][bufferIdx];
        min(max(Vec4d(mPan[0] - panMod), -1.), 1.).store(mPanBuf.Get() + bufferIdx);
        min(max(Vec4d(mPan[1] + panMod), -1.), 1.).store(mPanBuf.Get() + nFrames + bufferIdx);
#else
        std::array<T, OUTPUT_SIZE> osc1Output{ mOsc1.ProcessMultiple(osc1Freq) };
        std::array<T, OUTPUT_SIZE> osc2Output{ mOsc2.ProcessMultiple(osc2Freq) };
        for (auto j = 0; j < FRAME_INTERVAL; ++j)
        {
          osc1Output[j] *= osc1Amp;
          osc2Output[j] *= osc2Amp;
          // Filters
          T filter1Output = mFilters[0]->Process(osc1Output[j] * mFilterSends[0][0] + osc2Output[j] * mFilterSends[0][1]);
          T filter2Output = mFilters[1]->Process(osc1Output[j] * mFilterSends[1][0] + osc2Output[j] * mFilterSends[1][1]);
          T filterOutputs = filter1Output + filter2Output;

          // Panning
          T panMod = mVModulations.GetList()[kVPan][bufferIdx + j];
          mPanBuf.Get()[bufferIdx + j] = std::clamp(mPan[0] - panMod, -1., 1.); // TODO: optimize this somehow
          mPanBuf.Get()[bufferIdx + j + nFrames] = std::clamp(mPan[1] + panMod, -1., 1.);
          mEffectInputs.Get()[bufferIdx + j] = filterOutputs + osc1Output[j] * mFilterBypasses[0] + osc2Output[j] * mFilterBypasses[1];
        }
#endif
      }
#ifdef VECTOR
      for (int i{ 0 }; i < nFrames; ++i)
      {
        const T osc1Output = mOscOutputs.GetList()[0][i];
        const T osc2Output = mOscOutputs.GetList()[1][i];

        // Filters
        T filter1Output = mFilters[0]->Process(osc1Output * mFilterSends[0][0] + osc2Output * mFilterSends[0][1]);
        T filter2Output = mFilters[1]->Process(osc1Output * mFilterSends[1][0] + osc2Output * mFilterSends[1][1]);
        T filterOutputs = filter1Output + filter2Output;

        // Effect Sends
        mEffectInputs.Get()[i] = filterOutputs + osc1Output * mFilterBypasses[0] + osc2Output * mFilterBypasses[1];
      }
#endif

      // Process Effects at 4x oversampling
      UpsampleBlock<EFFECT_OS_FACTOR>(mOversampler, mEffectInputs.Get(), nFrames);
      constexpr int nEffectParams = kVEffect2Param1 - kVEffect1Param1;
#pragma clang loop unroll(full)
      for (int e{ 0 }, p{ 0 }; e < TABLITSA2_MAX_VOICE_EFFECTS; e++, p += nEffectParams)
      {
        mEffects[e]->SetContinuousParams(mVModulations.GetList()[kVEffect1Param1 + p][0],
          mVModulations.GetList()[kVEffect1Param2 + p][0],
          mVModulations.GetList()[kVEffect1Param3 + p][0],
          mVModulations.GetList()[kVEffect1Param4 + p][0],
          pitch);
        mEffects[e]->ProcessBlock(mOversampler.mOutputSource->Get(), mOversampler.mOutputSource->Get(), nFrames * EFFECT_OS_FACTOR);
      }
      DownsampleBlock<EFFECT_OS_FACTOR>(mOversampler, mEffectOutputs.Get(), nFrames);

      {
        std::lock_guard<std::mutex> lg(mMaster->mProcMutex);
        for (auto i{ startIdx }; i < startIdx + nFrames; ++i)
        {
          const int ii{ i - startIdx };
          // Sum effect outputs with the oscillator bypass sends, and multiply everything by the amplitude envelope
          const T out{ (mEffectOutputs.Get()[ii] + mOscOutputs.GetList()[0][ii] * mEffectBypasses[0] + mOscOutputs.GetList()[1][ii] * mEffectBypasses[1]) * mModulators.GetList()[2][ii] };
          outputs[0][i] += out * mPanBuf.Get()[ii] * mGain;
          outputs[1][i] += out * mPanBuf.Get()[ii + nFrames] * mGain;
        }
      }
    }

    void SetSampleRateAndBlockSize(double sampleRate, int blockSize) override
    {
      mOsc1.SetSampleRate(sampleRate);
      mOsc2.SetSampleRate(sampleRate);
      mAmpEnv.SetSampleRate(sampleRate);
      mEnv1.SetSampleRate(sampleRate);
      mEnv2.SetSampleRate(sampleRate);
      mLFO1.SetSampleRate(sampleRate);
      mLFO2.SetSampleRate(sampleRate);

      mFilters[0]->SetSampleRate(sampleRate);
      mFilters[1]->SetSampleRate(sampleRate);

      for (int i{ 0 }; i < TABLITSA2_MAX_VOICE_EFFECTS; ++i)
        mEffects[i]->SetSampleRate(sampleRate, (int)EFFECT_OS_FACTOR);

      // Modulation buffers
      mVModulationsData.Resize(blockSize * kNumModulations);
      mVModulations.Empty();

      // Modulator value buffer
      mModulators.EmptyAndResize(blockSize, kNumModulators);

      // Voice output buffers
      mPanBuf.Resize(blockSize * 2);

      mOscBuffer.Resize(blockSize * 2);
      mOscOutputs.Add(mOscBuffer.Get());
      mOscOutputs.Add(mOscBuffer.Get() + blockSize);

      mEffectInputs.Resize(blockSize);
      mEffectOutputs.Resize(blockSize);
      mOutputs.Resize(blockSize * 2);
      mOversampler.ResizeBuffers(blockSize);

      for (auto i = 0; i < kNumVoiceModulations; i++)
      {
        mVModulations.Add(mVModulationsData.Get() + static_cast<size_t>(blockSize) * i);
      }
      
    }

    void SetProgramNumber(int pgm) override
    {
      //TODO:
    }

    // this is called by the VoiceAllocator to set generic control values.
    
    void SetControl(int controlNumber, float value) override
    {
      switch(controlNumber)
      {
      case IMidiMsg::kModWheel:
        break;
      }
    }

    void SetFilterType(int filter, int filterType)
    {
      std::lock_guard<std::mutex> lg(mMaster->mEffectMutex);
      if (mFilters[filter])
      {
        delete mFilters[filter];
        mFilters[filter] = nullptr;
      }
      // Indices of cutoff/res/drive for the given filter
      switch (filterType) {
      case kVSF:
      {
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0.001, 0.49);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1.);
        mFilters.at(filter) = new SVF2<T>(mMaster->mSampleRate);
        break;
      }
      case kMoog:
      {
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0.001, 0.49);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1.);
        mFilters.at(filter) = new MoogLadder<T>(mMaster->mSampleRate);
        break;
      }
      case kComb:
      {
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1.);
        mFilters.at(filter) = new CombFilter<T>(mMaster->mSampleRate);
        break;
      }
      default:
        mFilters.at(filter) = new NullFilter<T>(mMaster->mSampleRate);
        break;

      }
    }

    void SetEffect(const int effectSlot, const int effectId)
    {
      constexpr int numEffectModParams = kVEffect2Param1 - kVEffect1Param1;
      std::lock_guard<std::mutex> lg(mMaster->mEffectMutex);

      if (mEffects[effectSlot])
        delete mEffects[effectSlot];
      switch (effectId)
      {
      case kDistortionEffect:
      {
        mVoiceModParams[kVEffect1Param1 + effectSlot * numEffectModParams].SetMinMax(0., static_cast<double>(kNumDistortionModes - 1) + 0.1);
        mVoiceModParams[kVEffect1Param2 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param3 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param4 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mEffects[effectSlot] = new DistortionEffect<T>(mMaster->mSampleRate);
        break;
      }
      case kSampleAndHoldEffect:
      {
        mVoiceModParams[kVEffect1Param1 + effectSlot * numEffectModParams].SetMinMax(0.05, 10.);
        mVoiceModParams[kVEffect1Param2 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param3 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param4 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mEffects[effectSlot] = new SampleAndHold<T>(mMaster->mSampleRate);
        break;
      }
      case kTexturizerEffect:
      {
        mVoiceModParams[kVEffect1Param1 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param2 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param3 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param4 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mEffects[effectSlot] = new Texturizer<T>(mMaster->mSampleRate);
        break;
      }
      case kCMEffect:
      {
        mVoiceModParams[kVEffect1Param1 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param2 + effectSlot * numEffectModParams].SetMinMax(-2., 2.);
        mVoiceModParams[kVEffect1Param3 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mVoiceModParams[kVEffect1Param4 + effectSlot * numEffectModParams].SetMinMax(0., 1.);
        mEffects[effectSlot] = new CMEffect<T>(mMaster->mSampleRate);
        break;
      }
      default:
        mEffects[effectSlot] = new Effect<T>(mMaster->mSampleRate);
        break;

      }

      SetSampleRateAndBlockSize(mMaster->mSampleRate, mMaster->mBlockSize); // Set oversampled sample rate, since this isn't accounted for in effect constructors
    }

    /* Update polyphonic modulation depths */
    void UpdateVoiceParam(int voiceParam, int modIdx, double value)
    {
      mVoiceModParams[voiceParam].SetValue(modIdx, value); // Eventually adjust the value of modIdx in the SetParam switch statement rather than here
    }

    /* Update meta-modulation depths */
    void UpdateVoiceModulatorParam(int paramIdx, int modIdx, double value)
    {
      mVoiceMetaModParams[paramIdx - kNumVoiceModulations - 1].SetValue(modIdx, value);
    }

    Tablitsa2DSP<T>* GetMaster() const
    {
      return mMaster;
    }

  public:
    Tablitsa2DSP<T>* mMaster;

    WavetableOscillator<T> mOsc1{ 0, "Hydrogen" };
    WavetableOscillator<T> mOsc2{ 1, "Helium" };

    // Static Modulators
    T mKey{ 69. };
    T mVelocity{ 1. };
    T mTriggerRand{ 0.5 };

    // Dynamic Modulators
    Envelope<T> mEnv1;
    Envelope<T> mEnv2;
    Envelope<T> mAmpEnv;
    FastLFO<T> mLFO1;
    FastLFO<T> mLFO2;
    Sequencer<T, kNumSeqSteps> mSequencer;
    ModulatorList<T, Envelope<T>, FastLFO<T>, Sequencer<T>, kNumModulators> mModulators;

    FastSinOscillator<T> mModulatorOsc;
    EOscModulators mOscModSource{ EOscModulators::kSine };

    // Status parameters
    bool mLFO1Restart{ false };
    bool mLFO2Restart{ false };
    bool mSequencerRestart{ false };
    bool mLegato{ false };
    bool mTriggered{ false }; // Set to true to prevent the note being retriggered immediately after the initial trigger
    bool mOsc1FormantOn{ false };
    bool mOsc2FormantOn{ false };

    // Unison parameters
    double mDetune{ 0. };
    double mPan[2]{ 1., 1. };

    // Sample and Beat data
    double mTempo{ 120. };
    bool mTransportIsRunning{ false };
    double mQNPos{ 0. };

    // Static Modulators
    double mEnv1VelocityMod{ 0. };
    double mEnv2VelocityMod{ 0. };
    double mAmpEnvVelocityMod{ 1. };

    // Filters
    std::vector<Filter<T>*> mFilters{ nullptr, nullptr };

    // Voice Effects
    std::vector<Effect<T>*> mEffects{ new Effect<T>(DEFAULT_SAMPLE_RATE), new Effect<T>(DEFAULT_SAMPLE_RATE), new Effect<T>(DEFAULT_SAMPLE_RATE) };
    FastOversampler<T> mOversampler;

    // Temporary output buffers and routing matrices
    WDL_TypedBuf<T> mPanBuf;

    WDL_TypedBuf<T> mOscBuffer;
    WDL_PtrList<T> mOscOutputs;

    WDL_TypedBuf<T> mEffectInputs;
    WDL_TypedBuf<T> mEffectOutputs;
    WDL_TypedBuf<T> mOutputs;
    T mFilterSends[2][2]{ {1., 0.}, {0., 1.} };
    T mFilterBypasses[2]{ 0., 0. };
    T mEffectBypasses[2]{ 0., 0. };

    WDL_PtrList<T> mVModulations; // Pointers to modulator buffers
    WDL_TypedBuf<T> mVModulationsData; // Modulator buffer sample data

  private:
//    WDL_TypedBuf<float> mTimbreBuffer;
    ModulatedParameterList<T, kNumVoiceModulations, kNumModulators> mVoiceModParams{
      new ParameterModulator<kNumModulators>(-1., 1., "Pan"),
      new ParameterModulator<kNumModulators>(-24., 24., "Wt1 Pitch Offset"),
      new ParameterModulator<kNumModulators>(0., 1., "Wt1 Position"),
      new ParameterModulator<kNumModulators>(-1., 1., "Wt1 Bend"),
      new ParameterModulator<kNumModulators>(0., 0.99, "Wt1 Formant"),
      new ParameterModulator<kNumModulators>(0., 1., "Wt1 Amp"),
      new ParameterModulator<kNumModulators>(-24., 24., "Wt1 Pitch Offset"),
      new ParameterModulator<kNumModulators>(0., 1., "Wt2 Position"),
      new ParameterModulator<kNumModulators>(-1., 1., "Wt2 Bend"),
      new ParameterModulator<kNumModulators>(0., 0.99, "Wt2 Formant"),
      new ParameterModulator<kNumModulators>(0., 1., "Wt2 Amp"),
      new ParameterModulator<kNumModulators>(0.001, 0.5, "Flt1 Cutoff", true),
      new ParameterModulator<kNumModulators>(0., 1., "Flt1 Resonance"),
      new ParameterModulator<kNumModulators>(0., 1., "Flt1 Drive"),
      new ParameterModulator<kNumModulators>(0., 1., "Flt1 Comb FF"),
      new ParameterModulator<kNumModulators>(0., 1., "Flt1 Comb FB"),
      new ParameterModulator<kNumModulators>(0., 20., "Flt1 Comb Delay"),
      new ParameterModulator<kNumModulators>(0.001, 0.5, "Flt2 Cutoff", true),
      new ParameterModulator<kNumModulators>(0., 1., "Flt2 Resonance"),
      new ParameterModulator<kNumModulators>(0., 1., "Flt2 Drive"),
      new ParameterModulator<kNumModulators>(0., 1., "Flt2 Comb FF"),
      new ParameterModulator<kNumModulators>(0., 1., "Flt2 Comb FB"),
      new ParameterModulator<kNumModulators>(0., 20., "Flt2 Comb Delay"),
      new ParameterModulator<kNumModulators>(-24., 24., "Osc 1 Phase Mod Freq"),
      new ParameterModulator<kNumModulators>(-24., 24., "Osc 2 Phase Mod Freq"),
      new ParameterModulator<kNumModulators>(0., 1., "Osc 1 Phase Mod Depth"),
      new ParameterModulator<kNumModulators>(0., 1., "Osc 2 Phase Mod Depth"),
      new ParameterModulator<kNumModulators>(-24., 24., "Osc 1 Ring Mod Freq"),
      new ParameterModulator<kNumModulators>(-24., 24., "Osc 2 Ring Mod Freq"),
      new ParameterModulator<kNumModulators>(0., 1., "Osc 1 Ring Mod Depth"),
      new ParameterModulator<kNumModulators>(0., 1., "Osc 2 Ring Mod Depth"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 1 Param 1"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 1 Param 2"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 1 Param 3"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 1 Param 4"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 2 Param 1"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 2 Param 2"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 2 Param 3"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 2 Param 4"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 3 Param 1"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 3 Param 2"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 3 Param 3"),
      new ParameterModulator<kNumModulators>(0., 1., "Effect 3 Param 4"), };

    // Modulator parameters that can themselves be modulated
    ModulatedParameterList<T, kNumVoiceMetaModulations, kNumModulators> mVoiceMetaModParams{
      new ParameterModulator<kNumModulators>(0., 1., "Env1 Sustain"),
      new ParameterModulator<kNumModulators>(0., 1., "Env2 Sustain"),
      new ParameterModulator<kNumModulators>(0., 1., "AmpEnv Sustain"),
      new ParameterModulator<kNumModulators>(0.01, 40., "LFO1 Rate Hz", true),
      new ParameterModulator<kNumModulators>(-1., 1., "LFO1 Amp"),
      new ParameterModulator<kNumModulators>(0.01, 40., "LFO2 Rate Hz", true),
      new ParameterModulator<kNumModulators>(-1., 1., "LFO2 Amp"),
      new ParameterModulator<kNumModulators>(0.01, 40., "Sequencer Rate Hz", true),
      new ParameterModulator<kNumModulators>(0., 1., "Sequencer Amp"),
    };

    int mID;

    // noise generator for test
    uint32_t mRandSeed = 0;

    // return single-precision floating point number on [-1, 1]
    float Rand()
    {
      mRandSeed = mRandSeed * 0x0019660D + 0x3C6EF35F;
      uint32_t temp = ((mRandSeed >> 9) & 0x007FFFFF) | 0x3F800000;
      return (*reinterpret_cast<float*>(&temp))*2.f - 3.f;
    }

  };

/* end Voice class */

public:
#pragma mark -
  Tablitsa2DSP(int nVoices)
  {
    for (auto i = 0; i < nVoices; i++)
    {
      // add a voice to Zone 0.
      mSynthVoices.push_back(new Voice(this, i));
      mSynth.AddVoice(mSynthVoices.at(i), 0);
    }

    // some MidiSynth API examples:
    // mSynth.SetKeyToPitchFn([](int k){return (k - 69.)/24.;}); // quarter-tone scale
    // mSynth.SetNoteGlideTime(0.5); // portamento
  }

  void ProcessBlock(T** inputs, T** outputs, int nOutputs, int nFrames, double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
	{
	  std::lock_guard<std::mutex> lg(mEffectMutex);

		// clear outputs
		for(auto i = 0; i < nOutputs; i++)
		{
			memset(outputs[i], 0, nFrames * sizeof(T));
		}
		// Process global modulators
		mGlobalLFO1.FillBuffer(nFrames);
		mGlobalLFO2.FillBuffer(nFrames);
		mGlobalSequencer.FillBuffer(nFrames);

		// Process voices
		mParamSmoother.ProcessBlock(mParamsToSmooth, mModulations.GetList(), nFrames); // Populate modulations list (to be sent to mSynth as inputs)
		SetTempoAndBeat(qnPos, transportIsRunning, tempo, mTSNum, mTSDenom);
		mSynth.ProcessBlock(mModulations.GetList(), outputs, 0, nOutputs, nFrames);

		for (int s{ 0 }; s < nFrames; ++s)
		{
			const T smoothedGain = mModulations.GetList()[kModGainSmoother][s];

			// Master effects processing
			StereoSample<T> stereo_in{ outputs[0][s], outputs[1][s] };
			mEffects[0]->ProcessStereo(stereo_in);
			mEffects[1]->ProcessStereo(stereo_in);
			mEffects[2]->ProcessStereo(stereo_in);
			
			outputs[0][s] = stereo_in.l;
			outputs[1][s] = stereo_in.r;

			outputs[0][s] *= smoothedGain;
			outputs[1][s] *= smoothedGain;
		}
	}

  void Reset(double sampleRate, int blockSize)
  {
    mBlockSize = blockSize;
	  mSampleRate = sampleRate;

    mSynth.SetSampleRateAndBlockSize(sampleRate, blockSize);
    ResetAllVoices();
    mSynth.ForEachVoice([sampleRate](SynthVoice& voice) {
      for (auto f : dynamic_cast<Tablitsa2DSP::Voice&>(voice).mFilters)
        f->SetSampleRate(sampleRate);
      });

    // Global modulators
    mGlobalLFO1.SetSampleRate(sampleRate);
    mGlobalLFO2.SetSampleRate(sampleRate);
    mGlobalSequencer.SetSampleRate(sampleRate);
    mGlobalLFO1.Resize(blockSize);
    mGlobalLFO2.Resize(blockSize);
    mGlobalSequencer.Resize(blockSize);

    // Param Smoother list
    mModulationsData.Resize(blockSize * kNumModulations);
    mModulations.Empty();

    // Effects
    mEffects[0]->SetSampleRate(sampleRate);
    mEffects[1]->SetSampleRate(sampleRate);
    
    for(auto i = 0; i < kNumModulations; i++)
    {
      mModulations.Add(mModulationsData.Get() + static_cast<size_t>(blockSize) * i);
    }
  }

  void ResetAllVoices()
  {
    mSynth.Reset(); // Note: this sets all envelopes to the release stage, meaning all voices are technically active
    mSynth.ForEachVoice([](SynthVoice& voice) {
      dynamic_cast<Tablitsa2DSP::Voice&>(voice).mAmpEnv.Kill(true); // Set all amp envelopes to idle
      });
  }

  void ProcessMidiMsg(const IMidiMsg& msg, float detune=0.f)
  {
    if (!mConstantGlideTime && msg.StatusMsg() == IMidiMsg::EStatusMsg::kNoteOn)
    {
      mSynth.SetNoteGlideTime(std::abs(static_cast<double>(msg.NoteNumber()) - mLastNoteOn) / mGlideRateScalar);
      mLastNoteOn = static_cast<double>(msg.NoteNumber());
    }
    if (mGlobalSequencer.GetMode() == Sequencer<T>::EStepMode::kLockToGate && msg.StatusMsg() == IMidiMsg::kNoteOn)
    {
      mGlobalSequencer.GateOn();
      mSynth.ForEachVoice([](SynthVoice& voice) {
        dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.GateOn();
        });
    }
    mSynth.AddMidiMsgToQueue(msg);
  }

  inline void SetTempoAndBeat(double qnPos, bool transportIsRunning, double tempo, const int tsNum, const int tsDenom)
  {
    mGlobalMetronome.Set(qnPos, tempo, transportIsRunning, tsNum, tsDenom);
  }

  void UpdateOscillatorWavetable(int wtIdx, int oscIdx)
  {
    mTableLoading = true; // Used in plugin's midi message function to block new notes while loading a new wavetable
    ResetAllVoices();
    std::string wtName = mWavetableNames.at(wtIdx);
    WtFile wtFile{ wtName };
    // Check for invalid file
    if (!wtFile.Success())
    {
      throw std::ifstream::failure(std::string(wtName));
    }

    Wavetable<T>* pWT = new Wavetable<T>(wtFile);

    if (oscIdx == 0)
    {
      ForEachVoice([this, oscIdx, pWT](Voice& voice) {
        voice.mOsc1.SetWavetable(pWT);
        voice.mOsc1.ReloadLUT();
        });
      ForEachVoice([this, oscIdx](Voice& voice) {
        voice.mOsc1.NotifyLoaded();
        });
    }
    else
    {
      ForEachVoice([this, oscIdx, pWT](Voice& voice) {
        voice.mOsc2.SetWavetable(pWT);
        voice.mOsc2.ReloadLUT();
        });
      ForEachVoice([this, oscIdx](Voice& voice) {
        voice.mOsc2.NotifyLoaded();
        });
    }
    // check for invalid indices
    if (mWavetables[oscIdx])
      delete mWavetables[oscIdx];
    mWavetables[oscIdx] = pWT;
    mTableLoading = false;
  }

  void ResetDetune()
  {
    mDetuner.RefillBuffers();
  }

  int GetSequencerStep()
  {
    if (mActiveSequencer)
      return mActiveSequencer->GetCurrentStep();
    else
      return 0;
  }

  void UpdateFilterSource(int filterIdx, int oscIdx)
  {
    mSynth.ForEachVoice([filterIdx, oscIdx](SynthVoice& voice) {
      dynamic_cast<Tablitsa2DSP<T>::Voice&>(voice).UpdateFilterSource(filterIdx, oscIdx);
      });
  }

  inline void ForEachVoice(std::function<void(Voice& voice)> func)
  {
    for (auto v : mSynthVoices)
    {
      func(*v);
    }
  }

  void SetParam(int paramIdx, double value)
  {
    using EEnvStage = ADSREnvelope<sample>::EStage;
    
    switch (paramIdx) {
      case kParamNoteGlideTime:
        mSynth.SetNoteGlideTime(value / 1000.);
        break;
      case kParamNoteGlideRate:
        mGlideRateScalar = value;
        break;
      case kParamPortamentoMode:
        mConstantGlideTime = value > 0.5;
        break;
      case kParamMonophonic:
      {
        mMono = value > 0.5;
        ResetAllVoices();
        mSynth.SetPolyMode(value > 0.5 ? VoiceAllocator::EPolyMode::kPolyModeMono : VoiceAllocator::EPolyMode::kPolyModePoly);
        break;
      }
      case kParamGain:
      {
        mParamsToSmooth[kModGainSmoother] = std::pow(10., value / 20.);
        break;
      }
      case kParamPan:
        mParamsToSmooth[kModPanSmoother] = (T)value / 90.;
        break;
      case kParamVibratoSpeed:
        mVibratoOsc.SetFreqCPS(value * FRAME_INTERVAL); // Vibrato osc must be run fast in order to compensate for SIMD processing
        break;
      case kParamVibratoDepth:
        mVibratoDepth = value / 100.;
        break;
      case kParamEnv1Sustain:
      {
        mParamsToSmooth[kModEnv1SustainSmoother] = (T)value / 100.;
        break;
      }
      case kParamEnv1Attack:
      case kParamEnv1Decay:
      case kParamEnv1Release:
      {
        EEnvStage stage = static_cast<EEnvStage>(EEnvStage::kAttack + (paramIdx - kParamEnv1Attack));
        mSynth.ForEachVoice([stage, value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv1.SetStageTime(stage, value);
          });
        break;
      }
      case kParamEnv1Velocity:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv1VelocityMod = (T)value;
          });
        break;
      case kParamEnv1DecayCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv1.SetStageCurve(EEnvStage::kDecay, (T)value / 100.);
          });
        break;
      case kParamEnv1ReleaseCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv1.SetStageCurve(EEnvStage::kRelease, (T)value / 100.);
          });
        break;
      case kParamEnv2Attack:
      case kParamEnv2Decay:
      case kParamEnv2Release:
      {
        EEnvStage stage = static_cast<EEnvStage>(EEnvStage::kAttack + (paramIdx - kParamEnv2Attack));
        mSynth.ForEachVoice([stage, value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv2.SetStageTime(stage, value);
          });
        break;
      }
      case kParamEnv2Sustain:
      {
        mParamsToSmooth[kModEnv2SustainSmoother] = (T)value / 100.;
        /*ForEachVoice([paramIdx, value](Voice& voice) {
          voice.UpdateVoiceModulatorParam(kVEnv2Sustain, 0, (T)value / 100.);
          });*/
        break;
      }
      case kParamEnv2Velocity:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv2VelocityMod = (T)value;
          });
        break;
      case kParamEnv2DecayCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv2.SetStageCurve(EEnvStage::kDecay, (T)value / 100.);
          });
        break;
      case kParamEnv2ReleaseCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mEnv2.SetStageCurve(EEnvStage::kRelease, (T)value / 100.);
          });
        break;
      case kParamAmpEnvAttack:
      case kParamAmpEnvDecay:
      case kParamAmpEnvRelease:
      {
        EEnvStage stage = static_cast<EEnvStage>(EEnvStage::kAttack + (paramIdx - kParamAmpEnvAttack));
        mSynth.ForEachVoice([stage, value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mAmpEnv.SetStageTime(stage, value);
          });
        break;
      }
      case kParamAmpEnvSustain:
      {
        mParamsToSmooth[kModAmpEnvSustainSmoother] = (T)value / 100.;
        break;
      }
      case kParamAmpEnvVelocity:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mAmpEnvVelocityMod = (T)value;
          });
        break;
      case kParamAmpEnvDecayCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mAmpEnv.SetStageCurve(EEnvStage::kDecay, (T)value / 100.);
          });
        break;
      case kParamAmpEnvReleaseCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mAmpEnv.SetStageCurve(EEnvStage::kRelease, (T)value / 100.);
          });
        break;
      case kParamLFO1Amp:
      {
        mParamsToSmooth[kModLFO1AmpSmoother] = (T)value;
        mGlobalLFO1.SetScalar(value);
        ForEachVoice([value](Voice& v) {
          v.mLFO1.SetFreqCPS(value);
        });
        break;
      }
      case kParamLFO1RateTempo:
      {
        mGlobalLFO1.FastLFO<T>::SetQNScalarFromDivision(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO1.SetQNScalarFromDivision(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO1RateHz:
      {
        mParamsToSmooth[kModLFO1RateHzSmoother] = (T)value;
        mGlobalLFO1.FastLFO<T>::SetFreqCPS(value);
        ForEachVoice([value](Voice& v) {
          v.mLFO1.SetFreqCPS(value);
        });
        break;
      }
      case kParamLFO1RateMode:
      {
        mGlobalLFO1.SetRateMode(value > 0.5);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO1.SetRateMode(value > 0.5);
          });
        break;
      }
      case kParamLFO1Shape:
      {
        mGlobalLFO1.SetShape(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO1.SetShape(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO1Phase:
      {
        mGlobalLFO1.SetPhaseOffset(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO1.SetPhaseOffset(value);
          });
        break;
      }
      case kParamLFO1Restart:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO1Restart = (value > 0.5);
          });
        break;
      case kParamLFO2Amp:
      {
        mParamsToSmooth[kModLFO2AmpSmoother] = (T)value;
        mGlobalLFO2.SetScalar(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2.SetScalar(value);
          });
        break;
      }
      case kParamLFO2RateTempo:
      {
        mGlobalLFO2.FastLFO<T>::SetQNScalarFromDivision(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2.SetQNScalarFromDivision(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO2RateHz:
      {
        mParamsToSmooth[kModLFO2RateHzSmoother] = (T)value;
        mGlobalLFO2.SetFreqCPS(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2.SetFreqCPS(value);
          });
        break;
      }
      case kParamLFO2RateMode:
      {
        mGlobalLFO2.SetRateMode(value > 0.5);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2.SetRateMode(value > 0.5);
          });
        break;
      }
      case kParamLFO2Shape:
      {
        mGlobalLFO2.SetShape(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2.SetShape(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO2Phase:
      {
        mGlobalLFO2.SetPhaseOffset(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2.SetPhaseOffset(value);
          });
        break;
      }
      case kParamLFO2Restart:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLFO2Restart = (value > 0.5);
          });
        break;
      case kParamSequencerAmp:
      {
        mParamsToSmooth[kModSeqAmpSmoother] = (T)value;
        mGlobalSequencer.SetScalar(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetScalar(value);
          });
        break;
      }
      case kParamSequencerSteps:
      {
        mGlobalSequencer.SetLength(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetLength(static_cast<int>(value));
          });
        break;
      }
      case kParamSequencerStepMode:
      {
        int mode = static_cast<int>(value); // clang didn't like this being declared as an enum
        mGlobalSequencer.SetStepMode(static_cast<Sequencer<T>::EStepMode>(mode));
        mSynth.ForEachVoice([mode](SynthVoice& voice) {
        dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetStepMode(static_cast<Sequencer<T>::EStepMode>(mode));
        });
        break;
      }
      case kParamSequencerRateTempo:
      {
        mGlobalSequencer.Sequencer<T>::SetQNScalarFromDivision(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetQNScalarFromDivision(static_cast<int>(value));
          });
        break;
      }
      case kParamSequencerRateHz:
      {
        mParamsToSmooth[kModSeqRateHzSmoother] = (T)value;
        mGlobalSequencer.SetFreqCPS(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetFreqCPS(value);
          });
        break;
      }
      case kParamSequencerRateMode:
      {
        mGlobalSequencer.SetRateMode(value > 0.5);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetRateMode(value > 0.5);
          });
        break;
      }
      case kParamSequencerRestart:
      {
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencerRestart = (value > 0.5);
          });
        // Toggle between using a master/static phase to update the Sequencer display, and using the phase of the last-triggered voice
        break;
      }
      case kParamSequencerGlide:
      {
        T glideNorm = (T)value / 100.;
        mGlobalSequencer.SetGlide(glideNorm);
        mSynth.ForEachVoice([glideNorm](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mSequencer.SetGlide(glideNorm);
          });
        break;
      }
      case kParamLegato:
      {
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mLegato = value > 0.5;
          });
        if (!(value > 0.5))
        {
          std::fill(AmpEnvQueue.begin(), AmpEnvQueue.end(), nullptr);
          std::fill(Env1Queue.begin(), Env1Queue.end(), nullptr);
          std::fill(Env2Queue.begin(), Env2Queue.end(), nullptr);
        }
        break;
      }
      case kParamUnisonVoices:
      {
        int voices = static_cast<int>(value);
        mDetuner.SetNVoices(voices);
        mSynth.SetMonoUnison(voices);
        // If the unison voice number has increased, start new voices
        if (voices > mUnisonVoices)
        {
          for (int i{ mUnisonVoices }; i < voices; ++i)
          {
            // TODO
          }
        }
        else
        {
          // TODO: When the unison number is descreased, only stop any excess voices. (Will probably require keeping track of unison voices in the VoiceAllocator)
          mSynth.Reset();
        }
        mUnisonVoices = voices;
        break;
      }
      case kParamUnisonDetune:
      {
        mDetuner.SetMaxDetune(value);
        mDetuner.Reset();
        // Send new values to voices
        mSynth.ForEachVoice([](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP<T>::Voice&>(voice).ResetUnisonParams();
          });
        break;
      }
      case kParamUnisonChord:
        mDetuner.SetChord(static_cast<int>(value) + EUnisonChords::kNoChord);
        break;
      case kParamStereoSpread:
      {
        if (mDetuner.mNVoices == 1)
        {
          mStereoWidth = (T)value;
          mDetuner.SetMaxPan(mStereoWidth);
        }
        else
        {
          mStereoWidth = (T)value;
          mDetuner.SetMaxPan(mStereoWidth);
        }
        // Send new values to voices
        mDetuner.Reset();
        mSynth.ForEachVoice([](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP<T>::Voice&>(voice).ResetUnisonParams();
          });
        break;
      }
      case kParamWavetable1Pitch:
        mParamsToSmooth[kModWavetable1PitchSmoother] = value;
        break;
      case kParamWavetable1Amp:
        mParamsToSmooth[kModWavetable1AmpSmoother] = value;
        break;
      case kParamWavetable1Pos:
        mParamsToSmooth[kModWavetable1PosSmoother] = (T)value;
        break;
      case kParamWavetable1Bend:
        mParamsToSmooth[kModWavetable1BendSmoother] = value;
        break;
      case kParamOsc1FormantOn:
        mSynth.ForEachVoice([paramIdx, value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mOsc1FormantOn = value > 0.5;
          });
        break;
      case kParamWavetable1Formant:
        mParamsToSmooth[kModWavetable1FormantSmoother] = value / (T)100;
        break;
      case kParamWavetable2Pitch:
        mParamsToSmooth[kModWavetable2PitchSmoother] = value;
        break;
      case kParamWavetable2Amp:
        mParamsToSmooth[kModWavetable2AmpSmoother] = value;
        break;
      case kParamWavetable2Pos:
        mParamsToSmooth[kModWavetable2PosSmoother] = value;
        break;
      case kParamWavetable2Bend:
        mParamsToSmooth[kModWavetable2BendSmoother] = value;
        break;
      case kParamOsc2FormantOn:
        mSynth.ForEachVoice([paramIdx, value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mOsc2FormantOn = value > 0.5;
          });
        break;
      case kParamWavetable2Formant:
        mParamsToSmooth[kModWavetable2FormantSmoother] = value / (T)100;
        break;
      case kParamFilter1Type:
      {
        mFilter1Comb = static_cast<int>(value) == kComb;
        ForEachVoice([value](Voice& voice) {
          voice.SetFilterType(0, static_cast<int>(value));
          });
        break;
      }
      case kParamFilter1ModeVSF:
      case kParamFilter1ModeMoog:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mFilters.at(0)->SetMode(static_cast<int>(value));
          });
        break;
      case kParamFilter1Osc1Send:
      case kParamFilter1Osc2Send:
      {
        value = value < SEND_DB_FLOOR + 0.01 ? 0. : std::pow(10., value / 20.);
        ForEachVoice([paramIdx, value](Voice& voice) {
          voice.mFilterSends[0][paramIdx - kParamFilter1Osc1Send] = value;
          });
        break;
      }
      case kParamFilter1Cutoff:
        mParamsToSmooth[kModFilter1CutoffSmoother] = value / mSampleRate;
        break;
      case kParamFilter1FF:
      {
        if (mFilter1Comb)
          mParamsToSmooth[kModFilter1CutoffSmoother] = value / 100.;
        break;
      }
      case kParamFilter1Resonance:
        mParamsToSmooth[kModFilter1ResonanceSmoother] = value / 100.;
        break;
      case kParamFilter1FB:
      {
        if (mFilter1Comb)
          mParamsToSmooth[kModFilter1ResonanceSmoother] = value / 100.;
        break;
      }
      case kParamFilter1Drive:
      case kParamFilter1Delay:
      {
        mParamsToSmooth[kModFilter1DriveSmoother] = value / (mFilter1Comb ? 1000. / mSampleRate * (double)COMB_MAX_DELAY : 100.);
        break;
      }
      case kParamFilter2Type:
      {
        mFilter2Comb = static_cast<int>(value) == kComb;
        ForEachVoice([value](Voice& voice) {
          voice.SetFilterType(1, static_cast<int>(value));
          });
        break;
      }
      case kParamFilter2ModeVSF:
      case kParamFilter2ModeMoog:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<Tablitsa2DSP::Voice&>(voice).mFilters.at(1)->SetMode(static_cast<int>(value));
          });
        break;
      case kParamFilter2Osc1Send:
      case kParamFilter2Osc2Send:
      {
        value = value < SEND_DB_FLOOR + 0.01 ? 0. : std::pow(10., value / 20.);
        ForEachVoice([paramIdx, value](Voice& voice) {
          voice.mFilterSends[1][paramIdx - kParamFilter2Osc1Send] = value;
          });
        break;
      }
      case kParamFilter2Cutoff:
        mParamsToSmooth[kModFilter2CutoffSmoother] = value / (mFilter2Comb ? 100. : mSampleRate);
        break;
      case kParamFilter2FF:
      {
        if (mFilter2Comb)
          mParamsToSmooth[kModFilter2CutoffSmoother] = value / (mFilter2Comb ? 100. : mSampleRate);
        break;
      }
      case kParamFilter2Resonance:
        mParamsToSmooth[kModFilter2ResonanceSmoother] = value / 100.;
        break;
      case kParamFilter2FB:
      {
        if (mFilter2Comb)
          mParamsToSmooth[kModFilter2ResonanceSmoother] = value / 100.;
        break;
      }
      case kParamFilter2Drive:
      case kParamFilter2Delay:
      {
        mParamsToSmooth[kModFilter2DriveSmoother] = value / (mFilter2Comb ? 1000. / mSampleRate * (double)COMB_MAX_DELAY : 100.);
        break;
      }
      case kParamOsc1ModSource:
      {
        const int sourceID{ static_cast<int>(value) };
        break;
      }
      case kParamOsc1PM:
      {
        const bool phaseModOn = value > 0.5;
        ForEachVoice([phaseModOn](Voice& voice) {
          voice.mOsc1.SetPhaseModulation(phaseModOn);
          });
        break;
      }
      case kParamOsc2PM:
      {
        const bool phaseModOn = value > 0.5;
        ForEachVoice([phaseModOn](Voice& voice) {
          voice.mOsc2.SetPhaseModulation(phaseModOn);
          });
        break;
      }
      case kParamOsc1PhaseModFreq:
        mParamsToSmooth[kModOsc1PhaseModFreqSmoother] = value;
        break;
      case kParamOsc2PhaseModFreq:
        mParamsToSmooth[kModOsc2PhaseModFreqSmoother] = value;
        break;
      case kParamOsc1PhaseModAmount:
        mParamsToSmooth[kModOsc1PhaseModAmtSmoother] = value / 100.;
        break;
      case kParamOsc2PhaseModAmount:
        mParamsToSmooth[kModOsc2PhaseModAmtSmoother] = value / 100.;
        break;
      case kParamOsc1RM:
      {
        const bool ringModOn = value > 0.5;
        ForEachVoice([ringModOn](Voice& voice) {
          voice.mOsc1.SetRingModulation(ringModOn);
          });
        break;
      }
      case kParamOsc2RM:
      {
        const bool ringModOn = value > 0.5;
        ForEachVoice([ringModOn](Voice& voice) {
          voice.mOsc2.SetRingModulation(ringModOn);
          });
        break;
      }
      case kParamOsc1RingModFreq:
        mParamsToSmooth[kModOsc1RingModFreqSmoother] = value;
        break;
      case kParamOsc2RingModFreq:
        mParamsToSmooth[kModOsc2RingModFreqSmoother] = value;
        break;
      case kParamOsc1RingModAmount:
        mParamsToSmooth[kModOsc1RingModAmtSmoother] = value / 100.;
        break;
      case kParamOsc2RingModAmount:
        mParamsToSmooth[kModOsc2RingModAmtSmoother] = value / 100.;
        break;
      case kParamVoiceEffect1Param1: // Voice Effects
        mParamsToSmooth[kModVoiceEffect1Param1Smoother] = value;
        break;
      case kParamVoiceEffect1Param2:
        mParamsToSmooth[kModVoiceEffect1Param2Smoother] = value;
        break;
      case kParamVoiceEffect1Param4:
        mParamsToSmooth[kModVoiceEffect1Param4Smoother] = value;
        break;
      case kParamVoiceEffect1Param5:
        ForEachVoice([value](Voice& voice) {
          voice.mEffects[0]->SetParam5(value);
          });
        break;
      case kParamVoiceEffect2Param1: // Voice Effects
        mParamsToSmooth[kModVoiceEffect2Param1Smoother] = value;
        break;
      case kParamVoiceEffect2Param2:
        mParamsToSmooth[kModVoiceEffect2Param2Smoother] = value;
        break;
      case kParamVoiceEffect2Param3:
        mParamsToSmooth[kModVoiceEffect2Param3Smoother] = value;
        break;
      case kParamVoiceEffect2Param4:
        mParamsToSmooth[kModVoiceEffect2Param4Smoother] = value;
        break;
      case kParamVoiceEffect2Param5:
        ForEachVoice([value](Voice& voice) {
          voice.mEffects[1]->SetParam5(value);
          });
        break;
      case kParamVoiceEffect2Param6:
        ForEachVoice([value](Voice& voice) { voice.mEffects[1]->SetParam6(value); });
        break;
      case kParamVoiceEffect3Param1: // Voice Effects
        mParamsToSmooth[kModVoiceEffect3Param1Smoother] = value;
        break;
      case kParamVoiceEffect3Param2:
        mParamsToSmooth[kModVoiceEffect3Param2Smoother] = value;
        break;
      case kParamVoiceEffect3Param3:
        mParamsToSmooth[kModVoiceEffect3Param3Smoother] = value;
        break;
      case kParamVoiceEffect3Param4:
        mParamsToSmooth[kModVoiceEffect3Param4Smoother] = value;
        break;
      case kParamVoiceEffect3Param5:
        ForEachVoice([value](Voice& voice) {
          voice.mEffects[2]->SetParam5(value);
          });
        break;
      case kParamVoiceEffect3Param6:
        ForEachVoice([value](Voice& voice) { voice.mEffects[2]->SetParam6(value); });
        break;
      case kParamMasterEffect1Param1: // Master Effects
      case kParamMasterEffect2Param1:
      case kParamMasterEffect3Param1:
        mEffects[(paramIdx - kParamMasterEffect1Param1) / 6]->SetParam1((T)value);
        break;
      case kParamMasterEffect1Param2:
      case kParamMasterEffect2Param2:
      case kParamMasterEffect3Param2:
        mEffects[(paramIdx - kParamMasterEffect1Param2) / 6]->SetParam2((T)value);
        break;
      case kParamMasterEffect1Param3:
      case kParamMasterEffect2Param3:
      case kParamMasterEffect3Param3:
        mEffects[(paramIdx - kParamMasterEffect1Param3) / 6]->SetParam3((T)value);
        break;
      case kParamMasterEffect1Param4:
      case kParamMasterEffect2Param4:
      case kParamMasterEffect3Param4:
        mEffects[(paramIdx - kParamMasterEffect1Param4) / 6]->SetParam4((T)value);
        break;
      case kParamMasterEffect1Param5:
      case kParamMasterEffect2Param5:
      case kParamMasterEffect3Param5:
        mEffects[(paramIdx - kParamMasterEffect1Param5) / 6]->SetParam5((T)value);
        break;
      case kParamMasterEffect1Param6:
      case kParamMasterEffect2Param6:
      case kParamMasterEffect3Param6:
        mEffects[(paramIdx - kParamMasterEffect1Param6) / 6]->SetParam6((T)value);
        break;
      case kParamOsc1FilterBypass:
      case kParamOsc2FilterBypass:
      {
        ForEachVoice([paramIdx, value](Voice& voice) {
          voice.mFilterBypasses[paramIdx - kParamOsc1FilterBypass] = value;
          });
        break;
      }
      case kParamOsc1EffectBypass:
      case kParamOsc2EffectBypass:
      {
        ForEachVoice([paramIdx, value](Voice& voice) {
          voice.mEffectBypasses[paramIdx - kParamOsc1EffectBypass] = value;
          });
        break;
      }
      default:
        break;
    }
  }
  
public:
  MidiSynth mSynth { VoiceAllocator::kPolyModePoly, MidiSynth::kDefaultBlockSize };

  WDL_TypedBuf<T> mModulationsData; // Sample data for global modulations (e.g. smoothed sustain)
  WDL_PtrList<T> mModulations; // Ptrlist for global modulations
  LogParamSmooth<T, kNumModulations> mParamSmoother;
  sample mParamsToSmooth[kNumModulations];
  const std::vector<std::string> mWavetableNames{ ELEMENT_NAMES };
  std::vector<Voice*> mSynthVoices;
  Wavetable<T>* mWavetables[2]{ nullptr }; // Holds the wavetable data for all voices (for wavetable indices, see `mLoadedWavetables` below)
  bool mTableLoading{ false };

  // Polyphonic/Monophonic
  bool mMono{ false };
  int mUnisonVoices{ 1 };
  float mUnisonDetune{ 0.f };
  T mStereoWidth{ 0. };
  UnisonVoiceManager mDetuner{ kMaxUnisonVoices };

  // Portamento Parameters
  bool mConstantGlideTime{ true };
  double mGlideRateScalar{ 1. }; // Semitones per second
  double mLastNoteOn{ 0 };

  // Audio processing and musical parameters
  int mBlockSize{ DEFAULT_BLOCK_SIZE };
  double mSampleRate{ DEFAULT_SAMPLE_RATE };
  int mTSNum{ 4 };
  int mTSDenom{ 4 };
  double mTempo{ DEFAULT_TEMPO };

  // Status Variables
  int mSeqPos{ 0 };
  Sequencer<T, kNumSeqSteps>* mActiveSequencer{ nullptr };
  int mFilter1Osc{ 0 }; // Oscillator which provides the filter's input
  int mFilter2Osc{ 1 };
  bool mFilter1Comb{ false }; // Set to `true` when Filter 1 is a comb filter. Used for scaling delay/drive values by the proper amount.
  bool mFilter2Comb{ false }; // Set to `true` when Filter 2 is a comb filter. Used for scaling delay/drive values by the proper amount.

  // Non-modulatable parameters
  double mLoadedWavetables[2]{ 1., 2. }; // Integer indices (stored as double) of current wavetables
  double mSeqSteps[kNumSeqSteps]{}; // Value of each step in the sequencer
  int mStepPos{ 0 };
  int mPrevPos{ -1 };

  // Vibrato LFO
  FastSinOscillator<T> mVibratoOsc;
  T mVibratoDepth{ 0. };
  // Pointers to master modulators, for free-run and legato modes
  std::vector<Envelope<T>*> Env1Queue;
  std::vector<Envelope<T>*> Env2Queue;
  std::vector<Envelope<T>*> AmpEnvQueue;
  FastLFO<T>* mMasterLFO1{ nullptr }; // The last-triggered `mLFO1`, which "owns" the master phase
  FastLFO<T>* mMasterLFO2{ nullptr }; // The last-triggered `mLFO2`, which "owns" the master phase
  Sequencer<T>* mMasterSeq{ nullptr }; // The last-triggered `mSequencer`, which "owns" the master phase

  // Global Modulators
  ModMetronome mGlobalMetronome;
  GlobalModulator<T, FastLFO<T>> mGlobalLFO1{ &mGlobalMetronome };
  GlobalModulator<T, FastLFO<T>> mGlobalLFO2{ &mGlobalMetronome };
  GlobalModulator<T, Sequencer<T>> mGlobalSequencer{ &mGlobalMetronome, mSeqSteps };

  // Effects
  std::vector<Effect<T>*> mEffects{
    new Effect<T>(DEFAULT_SAMPLE_RATE),
    new Effect<T>(DEFAULT_SAMPLE_RATE),
    new Effect<T>(DEFAULT_SAMPLE_RATE)
  };

  std::mutex mProcMutex;
  std::mutex mEffectMutex;
};