#include "IGraphics.h"
#include "Tablitsa2.h"
#include "Tablitsa2Controls.h"

using namespace igraphics;
using namespace iplug;

/* Generic Slider */

void Tablitsa2SliderControl::DrawTrack(IGraphics& g, const IRECT& filledArea)
{
  const float extra = mHandleInsideTrack ? mHandleSize : 0.f;
  const IRECT adjustedTrackBounds = mDirection == EDirection::Vertical ? mTrackBounds.GetVPadded(extra) : mTrackBounds.GetHPadded(extra);
  // Padd the filled area less, to account for the asymmetric rectangular handle
  const IRECT adjustedFillBounds = mDirection == EDirection::Vertical ? filledArea.GetVPadded(extra / 2.f) : filledArea.GetHPadded(extra / 2.f);
  const float cr = GetRoundedCornerRadius(mTrackBounds);

  g.FillRoundRect(GetColor(kSH), adjustedTrackBounds, cr, &mBlend);
  g.FillRoundRect(GetColor(kX1), adjustedFillBounds, cr, &mBlend);

  if (mStyle.drawFrame)
    g.DrawRoundRect(GetColor(kFR), adjustedTrackBounds, cr, &mBlend, mStyle.frameThickness);
}

/* Mod Slider */

ModSliderControl::ModSliderControl(const IRECT& bounds, int sourceModIdx, const char* label, const IVStyle& style, bool valueIsEditable, EDirection dir, double gearing, float handleSize, float trackSize, bool handleInsideTrack) :
  Tablitsa2SliderControl(bounds, kNoParameter, label, style, valueIsEditable, dir, gearing, handleSize, trackSize, handleInsideTrack), mSourceModulator(sourceModIdx)
{
  mShape = EVShape::Rectangle;
  SetValue(0.5);
  Toggle(-1);
}

void ModSliderControl::SetModulationTarget(const int modIdx, IControl* pTarget) {
  mTargetParam = modIdx;
  mTargetControl = pTarget;
  double value = static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(modIdx, mSourceModulator);
  SetValue((value + 1.0)/2.0);
}

void ModSliderControl::SetValue(double value, int valIdx)
{
  IVSliderControl::SetValue(value, valIdx);
  if (mTargetParam > -1)
  {
    value = value * 2. - 1.;
    static_cast<Tablitsa2*>(GetDelegate())->SetModMatrixSlot(mTargetParam, mSourceModulator, value);
    if (mTargetControl)
      mTargetControl->SetDirty(false);
  }
}

void ModSliderControl::DrawWidget(IGraphics& g)
{
  float value = (float)GetValue(); // NB: Value is normalized to between 0. and 1.
  const IRECT handleBounds = (IsDisabled()) ? mTrackBounds.FracRect(mDirection, 0.f) : mTrackBounds.FracRect(mDirection, value);
  const IRECT filledTrack = mDirection == EDirection::Vertical ? (
    IsDisabled()) ?
    handleBounds :
    ((value >= 0.5f) ?
      mTrackBounds.GetGridCell(0, 0, 2, 1).FracRect(mDirection, 2.f * (value - 0.5f)) :
      mTrackBounds.GetGridCell(1, 0, 2, 1).FracRect(mDirection, 2.f * (0.5 - value), true)) : // <- Vertical
    (IsDisabled()) ?
    handleBounds :
    ((value >= 0.5f) ?
      mTrackBounds.GetGridCell(0, 1, 1, 2).FracRect(mDirection, 2.f * (value - 0.5f)) :
      mTrackBounds.GetGridCell(0, 0, 1, 2).FracRect(mDirection, 2.f * (0.5 - value), true)); // <- Horizontal

  if (mTrackSize > 0.f)
    DrawTrack(g, filledTrack);

  float cx, cy;

  const float offset = (mStyle.drawShadows && mShape != EVShape::Ellipse /* TODO? */) ? mStyle.shadowOffset * 0.5f : 0.f;

  if (mDirection == EDirection::Vertical)
  {
    cx = handleBounds.MW() + offset;
    cy = handleBounds.T;
  }
  else
  {
    cx = handleBounds.R;
    cy = handleBounds.MH() + offset;
  }

  if (mHandleSize > 0.f)
  {
    DrawHandle(g, { cx - mHandleSize, cy - mHandleSize, cx + mHandleSize, cy + mHandleSize });
  }
}

void ModSliderControl::DrawTrack(IGraphics& g, const IRECT& filledArea)
{
  const float extra = mHandleInsideTrack ? mHandleSize : 0.f;
  const IRECT adjustedTrackBounds = mDirection == EDirection::Vertical ? mTrackBounds.GetVPadded(extra) : mTrackBounds.GetHPadded(extra);
  // Pad the filled area less, to account for the asymmetric rectangular handle
  const IRECT adjustedFillBounds = mDirection == EDirection::Vertical ? filledArea.GetVPadded(extra / 2.f) : filledArea.GetHPadded(extra / 2.f);
  const float cr = GetRoundedCornerRadius(mTrackBounds);

  g.FillRoundRect(GetColor(kSH), adjustedTrackBounds, cr, &mBlend);
  g.FillRoundRect(GetColor(kX1), adjustedFillBounds, cr, &mBlend);

  if (mStyle.drawFrame)
    g.DrawRoundRect(GetColor(kFR), adjustedTrackBounds, cr, &mBlend, mStyle.frameThickness);
}

/* Generic Knob */

void Tablitsa2IVKnobControl::DrawPointer(IGraphics& g, float angle, float cx, float cy, float radius)
{
  const IColor pointerColor = GetColor(kFR).WithOpacity(1.f);
  g.DrawRadialLine(pointerColor, cx, cy, angle, mInnerPointerFrac * radius, mOuterPointerFrac * radius, &mBlend, mPointerThickness);
  float data1[2][2];
  float data2[2][2];
  iplug::igraphics::RadialPoints(angle, cx, cy, mInnerPointerFrac * radius, mOuterPointerFrac * radius, 2, data1);
  iplug::igraphics::RadialPoints(angle - 20., cx, cy, mInnerPointerFrac * radius, mOuterPointerFrac * radius * 0.65f, 2, data2);
  g.DrawLine(pointerColor, data1[1][0], data1[1][1], data2[1][0], data2[1][1], &mBlend, mPointerThickness);
}

/* Modulated Parameter Knob */

void Tablitsa2IVModKnobControl::LoadModParams()
{
  // TODO make a list of Control Tags that can be looped through
  if (mActive && GetActiveIdx() == GetUI()->GetControlIdx(this))
  {
    for (int i{ kCtrlTagEnv1Depth }; i <= kCtrlTagRandomDepth; ++i)
    {
      //dynamic_cast<ModSliderControl*>(GetUI()->GetControlWithTag(kCtrlTagEnv1Depth + i))->SetModulationTarget(-1);
      dynamic_cast<ModSliderControl*>(GetUI()->GetControlWithTag(i))->Toggle(false);
    }
    SetActiveIdx(false);
    mActive = false;
  }
  else
  {
    // Set all modulator sliders to the values of the currently-selected modulated parameter
    for (int i{ 0 }; i <= (kCtrlTagRandomDepth - kCtrlTagEnv1Depth); ++i)
    {
      ModSliderControl* pSlider{dynamic_cast<ModSliderControl*>(GetUI()->GetControlWithTag(kCtrlTagEnv1Depth + i))};
      pSlider->SetModulationTarget(mModParamIdx, this);
      pSlider->Toggle(true);
    }
    SetActiveIdx(true);
    mActive = true;
  }
  //mActive = !mActive;
}

void Tablitsa2IVModKnobControl::OnMouseDown(float x, float y, const IMouseMod& mod) 
{
  if (!mod.L || (mod.L && mod.A))
  {
    if (mModulationEnabled)
      LoadModParams();
    mMouseDown = !mMouseDown;
    /* By default, center-clicking causes the control to be captured such that it still responds to the mouse wheel when
    the mouse is not actually over it. ReleaseMouseCapture() empties the captured-control queue. */
    GetUI()->ReleaseMouseCapture();
    GetUI()->SetAllControlsDirty();
  }
  else
    IVKnobControl::OnMouseDown(x, y, mod);
}

void Tablitsa2IVModKnobControl::OnMouseWheel(float x, float y, const IMouseMod& mod, float d)
{
  if (mMouseIsOver)
    IVKnobControl::OnMouseWheel(x, y, mod, d);
}

void Tablitsa2IVModKnobControl::DrawWidget(IGraphics& g)
{
  float widgetRadius; // The radius out to the indicator track arc

  if (mWidgetBounds.W() > mWidgetBounds.H())
    widgetRadius = (mWidgetBounds.H() / 2.f);
  else
    widgetRadius = (mWidgetBounds.W() / 2.f);

  const float cx = mWidgetBounds.MW(), cy = mWidgetBounds.MH();

  widgetRadius -= (mTrackSize / 2.f);

  IRECT knobHandleBounds = mWidgetBounds.GetCentredInside((widgetRadius - mTrackToHandleDistance) * 2.f);
  const float angle = mAngle1 + (static_cast<float>(GetValue()) * (mAngle2 - mAngle1));
  DrawPressableShape(g, EVShape::Ellipse, knobHandleBounds, mActive && (GetActiveIdx() == GetUI()->GetControlIdx(this)), mMouseIsOver, IsDisabled());
  DrawIndicatorTrack(g, angle, cx, cy, widgetRadius);
  DrawPointer(g, angle, cx, cy, knobHandleBounds.W() / 2.f);
}

void Tablitsa2IVModKnobControl::DrawIndicatorTrack(IGraphics& g, float angle, float cx, float cy, float radius)
{
  // Set the origin of the track arch to the center of the range for dials with negative minimum values
  if (GetParam()->GetMin() < 0.)
    mAnchorAngle = 0.f;
  else
    mAnchorAngle = -135.f;
  if (mTrackSize > 0.f)
  {
    g.DrawArc(mBaseTrack, cx, cy, radius, angle >= mAnchorAngle ? mAnchorAngle : mAnchorAngle - (mAnchorAngle - angle), angle >= mAnchorAngle ? angle : mAnchorAngle, &mBlend, mTrackSize);

    // Envelope 1
    float env1Val = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kEnvelope1));
    if (env1Val)
    {
      float modAngle = std::max(std::min(env1Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[0][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 3);
    }

    // Envelope 2
    radius -= mTrackToHandleDistance / 2.f;
    float env2Val = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kEnvelope2));
    if (env2Val)
    {
      float modAngle = std::max(std::min(env2Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[1][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
    }

    // AmpEnv
    radius -= mTrackToHandleDistance / 4.f;
    float ampEnvVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kAmpEnvelope));
    if (ampEnvVal)
    {
      float modAngle = std::max(std::min(ampEnvVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[2][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }

    // LFO 1
    radius -= 3.f * mTrackToHandleDistance / 4.f;
    float lfo1Val = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kLFO1));
    if (lfo1Val)
    {
      float modAngle = std::max(std::min(lfo1Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[3][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 3);
      modAngle = std::max(std::min(-1.f * lfo1Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[3][0].WithContrast(0.5), cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 3);
    }

    // LFO 2
    radius += mTrackToHandleDistance / 4.f;
    float lfo2Val = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kLFO2));
    if (lfo2Val)
    {
      float modAngle = std::max(std::min(lfo2Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[4][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
      modAngle = std::max(std::min(-1.f * lfo2Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[4][0].WithContrast(0.5), cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
    }

    // Sequencer
    radius += mTrackToHandleDistance / 4.f;
    float seqVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kSequencer));
    if (seqVal)
    {
      float modAngle = std::max(std::min(seqVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[5][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }

    // Velocity
    radius -= mTrackToHandleDistance / 4.f;
    float velVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kVelocity));
    if (velVal)
    {
      float modAngle = std::max(std::min(velVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[6][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }

    // AfterTouch
    float atVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kAftertouch));
    if (atVal)
    {
      float modAngle = std::max(std::min(atVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[7][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }

    // Mod Wheel
    float mwVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kModWheel));
    if (mwVal)
    {
      float modAngle = std::max(std::min(mwVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[8][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }

    // Keytrack
    float ktrVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kKeytrack));
    if (ktrVal)
    {
      float modAngle = std::max(std::min(ktrVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[9][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }

    // Trigger Random
    float rndVal = static_cast<float>(static_cast<Tablitsa2*>(GetDelegate())->GetModMatrixSlot(mModParamIdx, EModulators::kRandom));
    if (rndVal)
    {
      float modAngle = std::max(std::min(rndVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[10][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }
  }
}

void Tablitsa2IVModKnobControl::EnableModulation(bool enabled)
{
  mModulationEnabled = enabled;
}

/* Tab Controls */

void Tablitsa2VTabBox::DrawWidget(IGraphics& g)
{
  const float cr = GetRoundedCornerRadius(mWidgetBounds);
  const float ft = mStyle.frameThickness;
  const float hft = ft / 2.f;

  int nPaths = /*mStyle.drawShadows ? 2 :*/ 1;

  auto b = mWidgetBounds.GetPadded(/*mStyle.drawShadows ? -mStyle.shadowOffset :*/ 0.f);

  auto labelT = mLabelBounds.Empty() ? mRECT.MH() : mLabelBounds.T;
  auto labelB = mLabelBounds.Empty() ? mRECT.MH() : mLabelBounds.B;
  auto labelR = mLabelBounds.Empty() ? mRECT.MW() : mLabelBounds.R;
  auto labelL = mLabelBounds.Empty() ? mRECT.MW() : mLabelBounds.L;

  for (int i = 0; i < nPaths; i++)
  {
    const float offset = i == 0 ? 0.f : mStyle.shadowOffset;
    g.PathClear();
    // Path around label
    g.PathMoveTo(labelL, b.T + hft - offset);
    g.PathArc(labelL + cr, labelT + cr + hft - offset, cr, 270.f, 0.f);
    g.PathArc(labelR - cr, labelT + cr + hft - offset, cr, 0.f, 90.f);
    g.PathLineTo(labelR, b.T + hft - offset);
    if (mActive)
    {
      // Path around control
      g.PathArc(b.R - cr - hft - offset, b.T + cr + hft - offset, cr, 0.f, 90.f);
      g.PathArc(b.R - cr - hft - offset, b.B - cr - hft - offset, cr, 90.f, 180.f);
      g.PathArc(b.L + cr + hft - offset, b.B - cr - hft - offset, cr, 180.f, 270.f);
      g.PathArc(b.L + cr + hft - offset, b.T + cr + hft - offset, cr, 270.f, 360.f);
      g.PathLineTo(labelL, b.T + hft - offset);
      g.PathFill(IPattern(GetColor(EVColor::kFG).WithOpacity(0.5)), IFillOptions(true), &mBlend);
      g.PathStroke(mStyle.drawShadows ? GetColor(i == 0 ? kSH : kFR) : GetColor(kFR), ft);
    }
    else
    {
      g.PathFill(IPattern(GetColor(mMouseIsOverLabel ? EVColor::kHL : EVColor::kBG).WithOpacity(0.85)), IFillOptions(true), &mBlend);
      g.PathStroke(mStyle.drawShadows ? GetColor(i == 0 ? kSH : kFR) : GetColor(kFR), ft);
    }
  }
}

void Tablitsa2VTabBox::SetActive(bool active)
{
  mActive = active;
  // Show/hide controls in this tab view
  GetUI()->ForControlInGroup(mGroupName.Get(), [active](IControl* control) {
      control->Hide(!active);
    });
}

void Tablitsa2VTabBox::SetGroupName(const char* newGroupName)
{
  // Hide old group
  GetUI()->ForControlInGroup(mGroupName.Get(), [](IControl* control) {
    if (!control->IsHidden())
      control->Hide(true);
    });
  mGroupName.Set(newGroupName);
  // show new group
  GetUI()->ForControlInGroup(mGroupName.Get(), [](IControl* control) {
    if (control->IsHidden())
      control->Hide(false);
    control->SetDirty(true); // Trigger control show/hide functions
    });
}

/* Effect Bank */

Tablitsa2EffectBankControl::Tablitsa2EffectBankControl(const IRECT& bounds, std::initializer_list<char*> labels, std::initializer_list<char*> groupNames, const IVStyle& style, const int maxTabs) :
  IControl(bounds, kNoParameter), mMaxTabs(maxTabs), mStyle(style)
{
  // Text to display on each tab
  for (auto l : labels)
  {
    int nLabels = mLabels.GetSize();
    if (nLabels >= mMaxTabs)
      break;

    mLabels.Resize(nLabels + 1);
    TabText* pLab = mLabels.Get() + nLabels;
    strcpy(pLab->mText, l);
  }

  // The groups associated with each tab
  for (auto g : groupNames)
  {
    int nGroups = mGroups.GetSize();
    mGroups.Resize(nGroups + 1);
    TabText* pName = mGroups.Get() + nGroups;
    strcpy(pName->mText, g);
  }

  // If more labels were supplied than groups, append empty strings to the group  names
  if (labels.size() > groupNames.size())
  {
    for (auto i{ groupNames.size() }; i < labels.size(); ++i)
    {
      int nGroups = mGroups.GetSize();
      mGroups.Resize(nGroups + 1);
      TabText* pName = mGroups.Get() + nGroups;
      strcpy(pName->mText, mGroups.Get()[i].mText);
    }
  }
}

void Tablitsa2EffectBankControl::OnMouseDown(float x, float y, const IMouseMod& mod)
{
  for (int i{ 0 }; i < mMaxTabs; ++i)
  {
    if (mTabs[i])
    {
      Tablitsa2VTabBox* tab = dynamic_cast<Tablitsa2VTabBox*>(mTabs[i]);
      if (!tab->IsActive() && tab->GetLabelBounds().Contains(x, y))
      {
        TabChanged(i);
      }
    }
  }
}

int Tablitsa2EffectBankControl::GetActiveTabIdx()
{
  for (int i{ 0 }; i < mMaxTabs; ++i)
  {
    if (dynamic_cast<Tablitsa2VTabBox*>(mTabs[i])->IsActive())
      return i;
  }
  return -1;
}

void Tablitsa2EffectBankControl::TabChanged(int newIdx, bool triggerAction)
{
  for (int i{ 0 }; i < mMaxTabs; ++i)
  {
    if (mTabs[i])
    {
      Tablitsa2VTabBox* tab = dynamic_cast<Tablitsa2VTabBox*>(mTabs[i]);
      tab->SetActive(i == newIdx);
    }
  }
  if (triggerAction && GetActionFunction())
    (GetActionFunction())(this);
}

/* Generic Dropdown List */

DropdownListControl::DropdownListControl(const IRECT& bounds, std::initializer_list<const char*> options, const IText& text, const IColor& bgColor, bool showLabel) :
  ICaptionControl(bounds, kNoParameter, text, bgColor, showLabel)
{
  for (auto s : options)
  {
    mOptions.push_back(std::string(s));
  }

  mMenu = new IPopupMenu("Menu", options, [this](IPopupMenu* pMenu) {
      SetCurrentIndex(pMenu->GetChosenItemIdx(), true);
    });
}

void DropdownListControl::SetCurrentIndex(const int newIdx, const bool triggerAction)
{
  if (newIdx == -1)
    return;
  mCurrentIdx = newIdx;
  SetDirty(triggerAction);
}

void DropdownListControl::Draw(IGraphics& g)
{
  if (mCurrentIdx >= 0 && !strcmp(mCustomStr.c_str(), ""))
    mStr.Set(mOptions[mCurrentIdx].c_str());
  else
    mStr.Set(mCustomStr.c_str());


  ITextControl::Draw(g);

  if (mTri.W() > 0.f)
  {
    g.FillTriangle(mMouseIsOver ? mTriangleMouseOverColor : mTriangleColor, mTri.L, mTri.T, mTri.R, mTri.T, mTri.MW(), mTri.B, GetMouseIsOver() ? 0 : &BLEND_50);
  }
}

void DropdownListControl::OnMouseDown(float x, float y, const IMouseMod& mod)
{
  IGraphics* pGraphics{ dynamic_cast<IGraphics*>(GetUI()) };
  pGraphics->CreatePopupMenu(*this, *mMenu, x, y);
}

void DropdownListControl::OnResize()
{
  if (mOptions.size() != 0)
  {
    mTri = mRECT.FracRectHorizontal(0.2f, true).GetCentredInside(IRECT(0, 0, 8, 5)); //TODO: This seems rubbish
  }
}

void DropdownListControl::SetCustomStr(const char* str)
{
  mCustomStr = str;
}

/* Preset Selection Control */

PresetSelector::PresetSelector(const IRECT& bounds, IPopupMenuControl* menu, std::initializer_list<char*> defaultPresets) :
  ICaptionControl(bounds, kNoParameter, TABLITSA2_TEXT, TABLITSA2_STYLE.colorSpec.GetColor(EVColor::kBG)), mMenu(menu)
{
  for (auto p : defaultPresets)
  {
    mDefaultPresets.push_back(std::string(p));
  }
  // LoadUserPresets(...);
  mAllPresets.insert(mAllPresets.begin(), mDefaultPresets.begin(), mDefaultPresets.end());
  mAllPresets.insert(mAllPresets.end(), mUserPresets.begin(), mUserPresets.end());
}

void PresetSelector::LoadUserPresets(std::initializer_list<char*> userPresets)
{
  for (auto p : userPresets)
  {
    mUserPresets.push_back(std::string(p));
  }
}

void PresetSelector::Draw(IGraphics& g)
{
  const IParam* pParam = GetParam();

  mStr.Set(mDefaultPresets[0].c_str());

  ITextControl::Draw(g);

  if (mTri.W() > 0.f)
  {
    g.FillTriangle(mMouseIsOver ? mTriangleMouseOverColor : mTriangleColor, mTri.L, mTri.T, mTri.R, mTri.T, mTri.MW(), mTri.B, GetMouseIsOver() ? 0 : &BLEND_50);
  }
}

/* Periodic Table */

void PeriodicTable::OnMouseDown(float x, float y, const IMouseMod& mod)
{
  for (int i{ 0 }; i < 2; ++i)
  {
    if (mCurrentElement == mSelectedElements[i])
    {
      mIsDragging = i;
    }
  }
}

void PeriodicTable::OnMouseUp(float x, float y, const IMouseMod& mod)
{
  if (mIsDragging > -1)
  {
    // Send new wavetable value (param method)
    // SetValue(static_cast<double>(mSelectedElements[mIsDragging]) / 118., mIsDragging);

    // Send new wavetable value (hidden param method)
    double newTableIdx = static_cast<double>(mSelectedElements[mIsDragging]) / 118.;
    GetUI()->GetDelegate()->SendArbitraryMsgFromUI(mIsDragging == 0 ? kMsgWavetable1Changed : kMsgWavetable2Changed, kCtrlTagPeriodicTable, sizeof(double), &newTableIdx);
    //      GetDelegate()->SendParameterValueFromUI(GetParamIdx(mIsDragging), static_cast<double>(mSelectedElements[mIsDragging] - 1) / 118.);
    //      GetDelegate()->OnParamChange(GetParamIdx(mIsDragging));
  };
  mIsDragging = -1;
  SetDirty(false);
}

void PeriodicTable::Draw(IGraphics& g)
{
  //g.FillRoundRect(IColor(255, 20, 0, 45), mRECT);
  // Highlight box occupied by mouse
  if (mCurrentElementCoords)
  {
    if (!((mCurrentElement >= 57 && mCurrentElement <= 71) || (mCurrentElement >= 89 && mCurrentElement <= 103)))
      g.FillRect(IColor(150, 200, 200, 0), mTableBounds.GetReducedFromLeft(std::floor(*(mCurrentElementCoords + 1) / 18. * mTableWidth)) \
        .GetReducedFromTop(std::floor(*mCurrentElementCoords / 7.f * mTableHeight)) \
        .GetFromLeft(E_CELL_WIDTH) \
        .GetFromTop(E_CELL_HEIGHT));
    else
      g.FillRect(IColor(150, 200, 200, 0), mLaAcBounds.GetReducedFromLeft(std::floor(*(mCurrentElementCoords + 1) / 15. * mLaAcWidth)) \
        .GetReducedFromTop(std::floor(*mCurrentElementCoords / 2.f * mLaAcHeight)) \
        .GetFromLeft(E_CELL_WIDTH) \
        .GetFromTop(E_CELL_HEIGHT));
  };
  for (int i{ 0 }; i < 2; ++i)
  {
    auto e = mSelectedElements[i];
    if (!((e >= 57 && e <= 71) || (e >= 89 && e <= 103)))
    {
      int* eCoords = ElementCoords[e - 1];
      g.FillRect(ElementIconColor[i].WithOpacity(0.8), mTableBounds.GetReducedFromLeft(std::floor(*(eCoords + 1) / 18. * mTableWidth)) \
        .GetReducedFromTop(std::floor(*eCoords / 7.f * mTableHeight)) \
        .GetFromLeft(E_CELL_WIDTH) \
        .GetFromTop(E_CELL_HEIGHT));
    }
    else
    {
      // Convert atomic number to index in Lanthanide-Actinide table coordinates. Note that Number{Ac} - 74 corresponds to index 15, the first column in the second row
      int* eCoords = (e >= 89) ? LaAcCoords[e - 74] : LaAcCoords[e - 57];
      g.FillRect(ElementIconColor[i].WithOpacity(0.8), mLaAcBounds.GetReducedFromLeft(std::floor(*(eCoords + 1) / 15. * mLaAcWidth)) \
        .GetReducedFromTop(std::floor(*eCoords / 2.f * mLaAcHeight)) \
        .GetFromLeft(E_CELL_WIDTH) \
        .GetFromTop(E_CELL_HEIGHT));
    }
  };

  // Draw current wavetable icons
  const IRECT activeElem1 = mRECT.GetFromRight(80.f).GetGridCell(0, 0, 2, 1).GetCentredInside(70.f, 95.f);
  const IRECT activeElem2 = mRECT.GetFromRight(80.f).GetGridCell(1, 0, 2, 1).GetCentredInside(70.f, 95.f);
  DrawElement(g, activeElem1, mSelectedElements[0], 0);
  DrawElement(g, activeElem2, mSelectedElements[1], 1);

  g.DrawSVG(mSVG, mRECT);
}

void PeriodicTable::DrawElement(IGraphics& g, const IRECT& bounds, int atomicNumber, int idx)
{
  assert(atomicNumber > 0);
  IColor col{ mTableLoading[idx] ? ElementIconColor[idx].WithOpacity(0.5f) : ElementIconColor[idx] };

  g.DrawRect(col, bounds);
  g.DrawText(LabelText, (idx) ? "Wavetable 2" : "Wavetable 1", bounds.GetVShifted(-20.f).GetFromTop(20.f));
  g.DrawText(ElementIconText.WithSize(42).WithFGColor(col), ElementSymbols[atomicNumber - 1], bounds.GetReducedFromTop(25.f).GetReducedFromBottom(30.f));

  // Check for multiline name and draw name
  std::string eName = ElementNames[atomicNumber - 1];
  const size_t breakpoint{ eName.find("\n") };
  if (breakpoint != std::string::npos)
  {
    g.DrawText(ElementIconText.WithSize(16).WithFGColor(col), (" " + eName.substr(0, breakpoint)).c_str(), bounds.GetFromBottom(25.f).GetFromTop(12.f));
    g.DrawText(ElementIconText.WithSize(16).WithFGColor(col), eName.substr(breakpoint).c_str(), bounds.GetFromBottom(12.f).GetFromTop(12.f));
  }
  else
    g.DrawText(ElementIconText.WithSize(16).WithFGColor(col), eName.c_str(), bounds.GetFromBottom(25.f).GetFromTop(12.f));
  // Draw atomic number
  g.DrawText(ElementIconText.WithSize(24).WithFGColor(col), std::to_string(atomicNumber).c_str(), bounds.GetFromTop(24.f));
}

void PeriodicTable::OnMouseOver(float x, float y, const IMouseMod& mod)
{
  if (x < mTableTLHC[0] || y < mTableTLHC[1])
  {
    mCurrentElementCoords = nullptr;
    mCurrentElement = -1;
    goto EndFunction;
  }
  float x_adj{ (x - mTableTLHC[0]) / mTableWidth };
  float y_adj{ (y - mTableTLHC[1]) / mTableHeight };
  int atomicNumber{ 1 };
  for (auto e : ElementCoords)
  {
    // Check main table (excluding lanthanides and actinides)
    if (e[0] < 7.)
    {
      if (static_cast<int>(y_adj * 7.f) == e[0] && static_cast<int>(x_adj * 18.f) == e[1])
      {
        mCurrentElementCoords = e;
        mCurrentElement = atomicNumber;
        goto EndFunction;
      }
    }
    // Increment the current atomic number and reset the Current Element index and pointer
    atomicNumber++;
    mCurrentElementCoords = nullptr;
    mCurrentElement = -1;
  }
  x_adj = (x - mLaTLHC[0]) / mLaAcWidth;
  y_adj = (y - mLaTLHC[1]) / mLaAcHeight;
  atomicNumber = 57;
  if (x < mLaTLHC[0] || y < mLaTLHC[1])
  {
    mCurrentElementCoords = nullptr;
    mCurrentElement = -1;
    goto EndFunction;
  }
  // Check the lanthanide and actinide table
  for (auto e : LaAcCoords)
  {
    if (static_cast<int>(y_adj * 2.f) == e[0] && static_cast<int>(x_adj * 15.f) == e[1])
    {
      mCurrentElementCoords = e;
      mCurrentElement = atomicNumber;
      break;
    }
    if (++atomicNumber > 71 && atomicNumber < 89)
      atomicNumber = 89;
    mCurrentElementCoords = nullptr;
    mCurrentElement = -1;
  }
EndFunction:
  if (mCurrentElement > -1 && mIsDragging > -1)
  {
    mSelectedElements[mIsDragging] = mCurrentElement;
  }
  SetDirty(false);
}

Tablitsa2ModPlotControl::Tablitsa2ModPlotControl(const IRECT& bounds, double* table, const int tableSize, const int numPoints, const IColor& color, const IVStyle& style, float gearing) :
  Tablitsa2ModPlotControl(bounds, kNoParameter, table, tableSize, numPoints, color, style, gearing)
{}

Tablitsa2ModPlotControl::Tablitsa2ModPlotControl(const IRECT& bounds, int paramIdx, double* table, const int tableSize, const int numPoints, const IColor& color, const IVStyle& style, float gearing) :
  IVPlotControl(bounds, { {color, [this](double x) {return mTable[(static_cast<unsigned int>(x * static_cast<int>(mTableSize) + 1) - mTablePhase) % mTableSize];}} }, numPoints, "", style),
  mTableSize(tableSize), mGearing(gearing)
{
  mEmptyTable = new double[mTableSize] {0.};
  if (table)
    mTable = table;
  else
    mTable = mEmptyTable;

  SetParamIdx(paramIdx);
}

void Tablitsa2ModPlotControl::SetPlotTable(const double* pTable)
{
  mTable = pTable;
  SetDirty(false);
}

void Tablitsa2ModPlotControl::OnMouseDrag(float x, float y, float dX, float dY, const IMouseMod& mod)
{
  mTablePhase += static_cast<int>(dX * mGearing);
  mTablePhase %= mTableSize;
  SetValue(static_cast<double>((mTablePhase - mTableSize) % mTableSize) / mTableSize);
  SetDirty();
}

Tablitsa2ParametricPlotControl::Tablitsa2ParametricPlotControl(const IRECT& bounds, int numPoints, int paramIdx, const char* label, const IVStyle& style, double min, double max, float gearing) :
  IControl(bounds, paramIdx), IVectorBase(style, false, false), mMin(min), mMax(max), mNumPoints(numPoints), mGearing(gearing)
{
  AttachIControl(this, label);
}

void Tablitsa2ParametricPlotControl::InterpolateCubicBezier(const float* p0, const float* p1, const float* p2, const float* p3, const float t, float* y)
{
  float a[2], b[2], c[2], d[2], e[2];
  InterpolateLinearBezier(p0, p1, t, a);
  InterpolateLinearBezier(p1, p2, t, b);
  InterpolateLinearBezier(p2, p3, t, c);
  InterpolateLinearBezier(a, b, t, d);
  InterpolateLinearBezier(b, c, t, e);
  InterpolateLinearBezier(d, e, t, y);
}

void Tablitsa2ParametricPlotControl::OnMouseDrag(float x, float y, float dX, float dY, const IMouseMod& mod)
{
  mTension += dY * mGearing * 0.01;
  mTension = Clip(mTension, 0.f, 1.f);
  SetValue(mTension);
  SetDirty();
}

void Tablitsa2ParametricPlotControl::Draw(IGraphics& g)
{
  DrawBackground(g, mRECT);
  DrawLabel(g);

  constexpr float anchor_left[]{ 0.f, 0.5f };
  constexpr float anchor_top[]{ 0.5f, 1.f };
  constexpr float anchor_bottom[]{ 0.5f, 0.f };
  constexpr float anchor_right[]{ 1.f, 0.5f };

  const float p0[]{ 0.f, 0.f };
  float p1[2];
  float p2[2];
  const float p3[]{ 1.f, 1.f };

  // Interpolate intermediate points based on current parameter value
  InterpolateLinearBezier(anchor_left, anchor_bottom, mTension, p1);
  InterpolateLinearBezier(anchor_right, anchor_top, 1.f - mTension, p2);

  float* x = new float[mNumPoints];
  float* y = new float[mNumPoints];
  float* p = new float[mNumPoints * 2];

  for (int i{}; i < mNumPoints; ++i)
  {
    const float t = static_cast<double>(i) / mNumPoints;
    InterpolateCubicBezier(p0, p1, p2, p3, t, p + i*2);
    x[i] = (p+i*2)[0];
    y[i] = (p+i*2)[1];
  }

  g.DrawData(COLOR_WHITE, mWidgetBounds, y, mNumPoints, x, 0, 2.0f);

  if (mStyle.drawFrame)
    g.DrawRect(GetColor(kFR), mWidgetBounds, &mBlend, mStyle.frameThickness);

  delete y;
}