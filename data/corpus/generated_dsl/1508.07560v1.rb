# =====================================================================
# ψ(3770) → D0 D0bar  ---  double-tag semileptonic analysis (TagAnalysis)
#   tag side : D0bar → 5 hadronic tag modes (DTagTool tag collection)
#   signal   : D0 → K− e+ νe   AND   D0 → π− e+ νe   (+ massless missing νe)
# Two independent signal modes → two TagAnalysis objects (Rule T1).
# =====================================================================

### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")      # 2.92 fb⁻¹ ψ(3770) data
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC

# --- Decay card, signal mode I : D0 → K− e+ νe -------------------------------
# Tag D0bar decays as a PDG-weighted cocktail of the five hadronic modes.
decay_card_Kenu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  0.0395 K+ pi-             VSS;
  0.1400 K+ pi- pi0         PHSP;
  0.0810 K+ pi- pi- pi+     PHSP;
  0.0400 K+ pi- pi+ pi- pi0 PHSP;
  0.0170 K+ pi- pi0 pi0     PHSP;
  Enddecay

  End
DECAYCARD

# --- Decay card, signal mode II : D0 → π− e+ νe ------------------------------
decay_card_pienu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 pi- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  0.0395 K+ pi-             VSS;
  0.1400 K+ pi- pi0         PHSP;
  0.0810 K+ pi- pi- pi+     PHSP;
  0.0400 K+ pi- pi+ pi- pi0 PHSP;
  0.0170 K+ pi- pi0 pi0     PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive MC : 1M events per signal mode, tag = PDG-weighted cocktail ----
exMC_Kenu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0toKenu_D0bar_cocktail"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Kenu
  config.cross_section   = :default
end

exMC_pienu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0topienu_D0bar_cocktail"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pienu
  config.cross_section   = :default
end

### Signal mode I : D0 → K− e+ νe ############################################
alg_K = TagAnalysis.new("D0toKenuTag")
alg_K.set_header(["D0toKenuTagAlg/D0toKenuTag.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     # --- tag-side reconstruction details handled inside DTagAlg (captured as notes) ---
     .note(:tag_selection, "tag charged tracks: |cosθ|<0.93, |Vr|<1.0 cm, |Vz|<15.0 cm; K/π PID from combined dE/dx + TOF confidence levels (CL_K>CL_π for K, CL_π>CL_K for π when p<0.75 GeV/c; CL>0.1% above 0.75 GeV/c); π0→γγ with photon E>25 MeV (barrel) / >50 MeV (endcap), in-time photons, photon-to-nearest-track angle >10° and a 1C mass-constrained fit with χ²<50")
     .note(:tag_Kpi_extra, "for the K+π− tag mode additionally require |t1−t2|<5 ns, opening angle <176° and Σ E/p<1.4")
     .note(:tag_deltaE_windows, "per-mode ΔE windows (GeV): Kπ [−0.049,+0.044], Kππ0 [−0.071,+0.052], Kπππ [−0.043,+0.043], Kππππ0 [−0.067,+0.066], Kππ0π0 [−0.082,+0.050]; keep the single tag per mode with the smallest |ΔE| inside its window, with M_BC in [1.858,1.875] GeV")
     .note(:signal_lepton_pid, "signal e+ : CL_e>0.1% and CL_e/(CL_e+CL_π+CL_K)>0.8 from dE/dx, TOF and EMC shower information; the signal hadron satisfies CL_K>CL_π or CL_π>CL_K")
     .note(:fsr_correction, "FSR photons emitted within the positron direction are added back before forming U_miss = E_miss − |p_miss| (U_miss stored for the signal extraction)")
     .note(:extra_photon_veto, "any unused recoil photon is required to have Emax < 300 MeV")
     .with_decay_card(decay_card_Kenu)

# tag side : D0bar (charm = −1) reconstructed in the five hadronic modes
alg_K.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.charm -1                              # pin the D0bar tag side
  t.window :mBC, min: 1.858, max: 1.875   # explicit M_BC window requested by the analysis
end

# signal side : exactly two oppositely charged recoil tracks (K− e+), net charge 0, + missing νe
alg_K.signal_side do |s|
  s.charged(km: 1, ep: 1)                 # one K− and one e+
  s.require_charge 0                      # net charge of the recoil pair is zero
  s.missing :nu_e                         # massless missing neutrino (auto-stores m_Umiss / m_q2)
end

# 4C kinematic fit over tag + K− + e+ + νe
alg_K.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_K.apply                                 # no Selection argument for a TagAnalysis

### Signal mode II : D0 → π− e+ νe ###########################################
alg_pi = TagAnalysis.new("D0topienuTag")
alg_pi.set_header(["D0topienuTagAlg/D0topienuTag.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .note(:tag_selection, "tag charged tracks: |cosθ|<0.93, |Vr|<1.0 cm, |Vz|<15.0 cm; K/π PID from combined dE/dx + TOF confidence levels (CL_K>CL_π for K, CL_π>CL_K for π when p<0.75 GeV/c; CL>0.1% above 0.75 GeV/c); π0→γγ with photon E>25 MeV (barrel) / >50 MeV (endcap), in-time photons, photon-to-nearest-track angle >10° and a 1C mass-constrained fit with χ²<50")
      .note(:tag_Kpi_extra, "for the K+π− tag mode additionally require |t1−t2|<5 ns, opening angle <176° and Σ E/p<1.4")
      .note(:tag_deltaE_windows, "per-mode ΔE windows (GeV): Kπ [−0.049,+0.044], Kππ0 [−0.071,+0.052], Kπππ [−0.043,+0.043], Kππππ0 [−0.067,+0.066], Kππ0π0 [−0.082,+0.050]; keep the single tag per mode with the smallest |ΔE| inside its window, with M_BC in [1.858,1.875] GeV")
      .note(:signal_lepton_pid, "signal e+ : CL_e>0.1% and CL_e/(CL_e+CL_π+CL_K)>0.8 from dE/dx, TOF and EMC shower information; the signal hadron satisfies CL_K>CL_π or CL_π>CL_K")
      .note(:fsr_correction, "FSR photons emitted within the positron direction are added back before forming U_miss = E_miss − |p_miss| (U_miss stored for the signal extraction)")
      .note(:extra_photon_veto, "any unused recoil photon is required to have Emax < 300 MeV")
      .with_decay_card(decay_card_pienu)

alg_pi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.charm -1
  t.window :mBC, min: 1.858, max: 1.875
end

alg_pi.signal_side do |s|
  s.charged(pim: 1, ep: 1)                # one π− and one e+
  s.require_charge 0
  s.missing :nu_e
end

alg_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_pi.apply

### Execute on datasets ######################################################
root_files_K  = alg_K.execute_on([psipp_data, psipp_incMC, exMC_Kenu])
root_files_pi = alg_pi.execute_on([psipp_data, psipp_incMC, exMC_pienu])