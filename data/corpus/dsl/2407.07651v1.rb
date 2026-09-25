# BESIII Analysis: e+e- -> Ds+ Ds1(2536)- / Ds+ Ds2*(2573)- cross sections & BF
# Paper: 2407.07651v1
# 15 energy points 4.530–4.946 GeV; exclusive analysis with Ds+ -> K-K+pi+ tag
# Measures absolute BF(Ds1(2536)- -> D*0bar K-) and BF(Ds2*(2573)- -> D0bar K-)

### Dataset preparation ###
data_4530 = DatasetManager.real_data.find("703_4530")
data_4575 = DatasetManager.real_data.find("703_4575")
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

scan_points = [
  data_4530, data_4575, data_4600, data_4610,
  data_4620, data_4640, data_4660, data_4680,
  data_4700, data_4740, data_4750, data_4780,
  data_4840, data_4914, data_4946
]

# ============================================================
# Algorithm 1: e+e- -> Ds+ Ds1(2536)-, Ds1(2536)- -> D*0bar K-
# Exclusive analysis: Ds+ -> K-K+pi+ + bachelor K- + missing D*0bar
# ============================================================

decay_card_DsDs1 = <<~DECAYCARD
  Decay psi(4260)
  1.0 D_s+ anti-D_s1- ANGSAM;
  Enddecay

  Decay anti-D_s1-
  1.0 anti-D*0 K- VVS_PWAVE;
  Enddecay

  Decay D_s+
  1.0 K- K+ pi+ D_Dalitz;
  Enddecay

  Decay anti-D*0
  1.0 anti-D0 pi0 PHSP;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_DsDs1 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_DsDs1_exclusive"
  config.events        = 200_000
  config.decay_card    = decay_card_DsDs1
  config.cross_section = :default
end

alg_DsDs1 = Algorithm.new("DsDs1Exclusive")
alg_DsDs1.set_header(["DsDs1ExclusiveAlg/DsDs1Exclusive.h"])

sel_DsDs1 = Selection.new

# Charged tracks: K- from Ds+, K+ from Ds+, pi+ from Ds+, bachelor K-
sel_DsDs1.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       ">=2"
end

# No photon selection needed for exclusive K-K+pi+ + K- channel

# PID: kaon and pion identification
sel_DsDs1.pid(method: :probability) do
  prob_cut   0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkm        ">=2"
  nkp        ">=1"
  npip       ">=1"
end

# Kinematic fit: 4C + M(Ds+) constraint + missing D*0bar mass constraint
# Participant list: K-(from Ds+), K+, pi+, K-(bachelor)
# The fit iterates over all assignments of the two K- to Ds+/bachelor roles
sel_DsDs1.kinematic_fit([:km, :kp, :pip, :km]) do
  constrain_four_momentum
  invariant_mass_of(:km, :kp, :pip).constrain_to_nominal_mass_of(:"D_s+")
  miss_track_of(:"anti-D*0")
  chi2_cut 200
  nominal
end

alg_DsDs1.with_decay_card(decay_card_DsDs1).apply(sel_DsDs1)

alg_DsDs1
  .note(:vertex_fit,
    "Tracks used to reconstruct Ds+ required to originate from common vertex " \
    "with chi2_VF < 100; applied in BOSS before kinematic fit")
  .note(:ds_mass_window,
    "|M(K-K+pi+) - m_Ds+| < 8 MeV/c^2 applied as pre-selection before kinematic fit")
  .note(:rm_window,
    "|RM(Ds+ K-) - m_D*0bar| < 9 MeV/c^2 for Ds+ -> K-K+pi+ channel; " \
    "applied as pre-selection before kinematic fit")
  .note(:phi_intermediate,
    "Ds+ -> phi(->K+K-) pi+ mode: M(K+K-) in [1.004, 1.034] GeV/c^2, " \
    "|cos(theta_K+/K+K-)| > 0.4 in helicity frame; applied in BOSS")
  .note(:kstar_intermediate,
    "Ds+ -> K*0bar(->K-pi+) K+ mode: M(K-pi+) in [0.832, 0.928] GeV/c^2, " \
    "|cos(theta_pi+/K-pi+)| > 0.52 in helicity frame; applied in BOSS")
  .note(:k0s_channel,
    "Ds+ -> KS0 K+, KS0 -> pi+pi- also used in exclusive analysis; " \
    "|M(KS0 K+) - m_Ds+| < 8 MeV/c^2, KS0 via secondary vertex fit; " \
    "RM(Ds+K-) window tighter: 7(9) MeV for Ds1(Ds2*) process; " \
    "combined with K-K+pi+ channel per their branching fractions")
  .note(:angular_distribution,
    "e+e- -> Ds+ Ds1(2536)- generated with ANGSAM model: " \
    "1 + alpha cos^2(theta) with alpha = -0.65 +/- 0.22")
  .note(:beam_energy_spread,
    "Beam energy spread and ISR considered via KKMC generator; " \
    "ISR correction factors obtained from MC simulation")
  .note(:inclusive_analysis,
    "Inclusive analysis also performed: Ds+ -> K-K+pi+ only, " \
    "Ds1/Ds2* signals extracted from 2D fit to M(K-K+pi+) vs RM(Ds+); " \
    "handled at ROOT level; not encoded here")

# ============================================================
# Algorithm 2: e+e- -> Ds+ Ds2*(2573)-, Ds2*(2573)- -> D0bar K-
# Exclusive analysis: Ds+ -> K-K+pi+ + bachelor K- + missing D0bar
# ============================================================

decay_card_DsDs2star = <<~DECAYCARD
  Decay psi(4260)
  1.0 D_s+ anti-D_s2*- PHSP;
  Enddecay

  Decay anti-D_s2*-
  1.0 anti-D0 K- PHSP;
  Enddecay

  Decay D_s+
  1.0 K- K+ pi+ D_Dalitz;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_DsDs2star = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_DsDs2star_exclusive"
  config.events        = 200_000
  config.decay_card    = decay_card_DsDs2star
  config.cross_section = :default
end

alg_DsDs2star = Algorithm.new("DsDs2starExclusive")
alg_DsDs2star.set_header(["DsDs2starExclusiveAlg/DsDs2starExclusive.h"])

sel_DsDs2star = Selection.new

# Charged tracks: K- from Ds+, K+ from Ds+, pi+ from Ds+, bachelor K-
sel_DsDs2star.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       ">=2"
end

# PID: kaon and pion identification
sel_DsDs2star.pid(method: :probability) do
  prob_cut   0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkm        ">=2"
  nkp        ">=1"
  npip       ">=1"
end

# Kinematic fit: 4C + M(Ds+) constraint + missing D0bar mass constraint
sel_DsDs2star.kinematic_fit([:km, :kp, :pip, :km]) do
  constrain_four_momentum
  invariant_mass_of(:km, :kp, :pip).constrain_to_nominal_mass_of(:"D_s+")
  miss_track_of(:"anti-D0")
  chi2_cut 200
  nominal
end

alg_DsDs2star.with_decay_card(decay_card_DsDs2star).apply(sel_DsDs2star)

alg_DsDs2star
  .note(:vertex_fit,
    "Tracks used to reconstruct Ds+ required to originate from common vertex " \
    "with chi2_VF < 100; applied in BOSS before kinematic fit")
  .note(:ds_mass_window,
    "|M(K-K+pi+) - m_Ds+| < 8 MeV/c^2 applied as pre-selection before kinematic fit")
  .note(:rm_window,
    "|RM(Ds+ K-) - m_D0bar| < 11 MeV/c^2 for Ds+ -> K-K+pi+ channel; " \
    "applied as pre-selection before kinematic fit")
  .note(:phi_intermediate,
    "Ds+ -> phi(->K+K-) pi+ mode: M(K+K+) in [1.004, 1.034] GeV/c^2, " \
    "|cos(theta_K+/K+K-)| > 0.4 in helicity frame; applied in BOSS")
  .note(:kstar_intermediate,
    "Ds+ -> K*0bar(->K-pi+) K+ mode: M(K-pi+) in [0.832, 0.928] GeV/c^2, " \
    "|cos(theta_pi+/K-pi+)| > 0.52 in helicity frame; applied in BOSS")
  .note(:k0s_channel,
    "Ds+ -> KS0 K+, KS0 -> pi+pi- also used in exclusive analysis; " \
    "|M(KS0 K+) - m_Ds+| < 8 MeV/c^2, KS0 via secondary vertex fit; " \
    "RM(Ds+K-) window tighter: 9 MeV for Ds2* process; " \
    "combined with K-K+pi+ channel per their branching fractions")
  .note(:dwave_production,
    "e+e- -> Ds+ Ds2*(2573)- generated via D-wave; " \
    "Ds2*(2573)- -> D0bar K- also via D-wave; " \
    "PHSP used in decay card as placeholder; actual BOSS MC uses D-wave model")
  .note(:beam_energy_spread,
    "Beam energy spread and ISR considered via KKMC generator; " \
    "ISR correction factors obtained from MC simulation")
  .note(:inclusive_analysis,
    "Inclusive analysis also performed: Ds+ -> K-K+pi+ only, " \
    "Ds1/Ds2* signals extracted from 2D fit to M(K-K+pi+) vs RM(Ds+); " \
    "handled at ROOT level; not encoded here")
  .note(:bachelor_kaon_pid,
    "Bachelor K- at low momentum: tracking/PID systematic uncertainties " \
    "estimated with control sample J/psi -> p K- Lambda; applied at ROOT level")

# Execute algorithms on scan data + exclusive MC
all_scan_datasets_Ds1 = scan_points + [exMC_DsDs1].flatten
all_scan_datasets_Ds2star = scan_points + [exMC_DsDs2star].flatten
alg_DsDs1.execute_on(all_scan_datasets_Ds1)
alg_DsDs2star.execute_on(all_scan_datasets_Ds2star)