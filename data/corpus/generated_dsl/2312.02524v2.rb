# frozen_string_literal: true

### Dataset description ###
# psi(3770) at ECMS = 3.773 GeV  ->  sample name follows the [BOSS]_[ECMS(MeV)] convention: 712_3773
d3770_data  = DatasetManager.real_data.find("712_3773")     # real data at 3.773 GeV
d3770_incMC = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC at 3.773 GeV
# NOTE: no exclusive signal MC is required by the description.

### Signal decay cards (attached for kinematic-variable header generation) ###
decay_card_4pi_charged = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_4pi_neutral = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

#=====================================================================
# Signal mode I : D0 -> pi+ pi- pi+ pi-   (tag side: D0bar)
#=====================================================================
alg_4pi_charged = TagAnalysis.new("D0barTagTo4piCharged")
alg_4pi_charged.set_header(["D0barTagTo4piChargedAlg/D0barTagTo4piCharged.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .with_decay_card(decay_card_4pi_charged)

# --- Tag side: D0bar reconstructed from the three hadronic D0 modes (charm = -1 pins the D0bar) ---
alg_4pi_charged.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # D0 -> K+ pi-, K+ pi- pi0, K+ pi- pi+ pi-
  t.charm -1                                    # tag the D0bar side
end

# --- Signal side: everything the tag did not use ---
alg_4pi_charged.signal_side do |s|
  s.photons 0                     # no photons on the signal side
  s.charged(pip: 2, pim: 2)       # exactly two pi+ and two pi-
  s.require_charge(0)             # total signal-side charge = 0
end

# --- Kinematic fit: 1C constraint of the four-pion invariant mass to the nominal D0 mass ---
alg_4pi_charged.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pip, :pim, :pim).constrain_to_nominal_mass_of(:D0)  # M(4pi) = m_D0
  f.chi2_cut 200
end

# --- Background vetoes: no dedicated TagAnalysis primitive -> preserved as a BOSS-side note ---
alg_4pi_charged.note(:background_veto,
  "shared vetoes: (1) reject K_S0 candidates when any pi+ pi- pair has |M(pi pi) - m_K_S0| < 12 MeV/c2; " \
  "(2) reject the D- -> K+ pi- pi- / K- pi+ pi+ hypothesis when |M - m_D-| < 10 MeV/c2; " \
  "(3) reject the eta hypothesis when any gamma gamma pair has |M(gamma gamma) - m_eta| < 20 MeV/c2; " \
  "(4) reject events with extra pi0 candidates beyond those forming the signal; " \
  "charged mode additionally vetoes D0 -> pi+ pi- pi0 when M(pi+ pi- pi0) is near the D0 mass. " \
  "These physics-motivated vetoes have no dedicated TagAnalysis primitive and are applied at the ROOT stage.")

alg_4pi_charged.apply
alg_4pi_charged.execute_on([d3770_data, d3770_incMC])

#=====================================================================
# Signal mode II : D0 -> pi+ pi- pi0 pi0   (tag side: D0bar; two pi0 -> 4 gamma)
#=====================================================================
alg_4pi_neutral = TagAnalysis.new("D0barTagTo4piNeutral")
alg_4pi_neutral.set_header(["D0barTagTo4piNeutralAlg/D0barTagTo4piNeutral.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .with_decay_card(decay_card_4pi_neutral)

# --- Tag side: same three hadronic D0bar modes, charm = -1 ---
alg_4pi_neutral.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# --- Signal side: one pi+ , one pi- and four photons (two pi0) ---
alg_4pi_neutral.signal_side do |s|
  s.photons 4                     # four signal photons -> two pi0
  s.charged(pip: 1, pim: 1)       # exactly one pi+ and one pi-
  s.require_charge(0)
end

# --- Kinematic fit: pi0 mass constraints from the photon pairs + four-pion mass to nominal D0 mass ---
alg_4pi_neutral.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)                        # M(gamma gamma) = m_pi0
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma, :gamma, :gamma).constrain_to_nominal_mass_of(:D0)  # M(4pi) = m_D0
  f.chi2_cut 200
end

# --- pi0 reconstruction (1C Kalman fit, chi2_1C < 20) has no TagAnalysis primitive -> note ---
alg_4pi_neutral.note(:pi0_reconstruction,
  "the two pi0 are reconstructed from photon pairs through a 1C Kalman kinematic fit constraining " \
  "M(gamma gamma) to the nominal pi0 mass with chi2_1C < 20 before the main fit; the mass constraint " \
  "is represented inside the fit block, the chi2_1C < 20 quality cut is applied at the ROOT stage.")

# --- Background vetoes: same shared set as the charged mode (no extra-D0 veto) ---
alg_4pi_neutral.note(:background_veto,
  "shared vetoes: (1) reject K_S0 candidates when any pi+ pi- pair has |M(pi pi) - m_K_S0| < 12 MeV/c2; " \
  "(2) reject the D- -> K+ pi- pi- / K- pi+ pi+ hypothesis when |M - m_D-| < 10 MeV/c2; " \
  "(3) reject the eta hypothesis when any gamma gamma pair has |M(gamma gamma) - m_eta| < 20 MeV/c2; " \
  "(4) reject events with extra pi0 candidates beyond the two forming the signal. " \
  "These physics-motivated vetoes have no dedicated TagAnalysis primitive and are applied at the ROOT stage.")

alg_4pi_neutral.apply
alg_4pi_neutral.execute_on([d3770_data, d3770_incMC])