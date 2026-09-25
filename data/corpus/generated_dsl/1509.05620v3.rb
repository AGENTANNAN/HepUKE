# ============================================================
# Dataset preparation
# ============================================================
# psi(4260) data and inclusive MC at the two energy points:
#   4.226 GeV (1092 /pb)  -> sample "703_4230"
#   4.257 GeV ( 826 /pb)  -> sample "703_4260"
data_4226  = DatasetManager.real_data.find("703_4230")
data_4257  = DatasetManager.real_data.find("703_4260")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4257 = DatasetManager.inclusive_mc.find("703_4260")

# Decay card -- Mode I: e+e- -> D+ D*- pi0, D*- -> anti-D0 pi-
# D+ : K-pi+pi+ / K-pi+pi+pi0 / K_S0 pi+ / K_S0 pi+pi0 / K_S0 pi+pi+pi-
# anti-D0 : K+pi- / K+pi-pi0 / K+pi-pi+pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D+ D*- pi0 PHSP;
    Enddecay

    Decay D*-
    1.0000 anti-D0 pi- PHSP;
    Enddecay

    Decay D+
    0.3000 K- pi+ pi+ PHSP;
    0.2000 K- pi+ pi+ pi0 PHSP;
    0.1500 K_S0 pi+ PHSP;
    0.1500 K_S0 pi+ pi0 PHSP;
    0.2000 K_S0 pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    0.5000 K+ pi- PHSP;
    0.3000 K+ pi- pi0 PHSP;
    0.2000 K+ pi- pi+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card -- Mode II: e+e- -> D0 anti-D*0 pi0, anti-D*0 -> anti-D0 pi0
# D0 / anti-D0 : Kpi / Kpipi0 / Kpipipi
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D0 anti-D*0 pi0 PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 pi0 PHSP;
    Enddecay

    Decay D0
    0.5000 K- pi+ PHSP;
    0.3000 K- pi+ pi0 PHSP;
    0.2000 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    0.5000 K+ pi- PHSP;
    0.3000 K+ pi- pi0 PHSP;
    0.2000 K+ pi- pi+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 200k events per signal mode, one sample per energy point
exMC_modeI = DatasetManager.create_exclusive_mc_for([data_4226, data_4257]) do |config|
  config.sample_name   = "exmc_ddstpi0_modeI"
  config.events        = 200_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for([data_4226, data_4257]) do |config|
  config.sample_name   = "exmc_ddbarstpi0_modeII"
  config.events        = 200_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

# ============================================================
# Mode I algorithm: e+e- -> D+ D*- pi0, D*- -> anti-D0 pi-
# ============================================================
alg_modeI = Algorithm.new("DDstPi0ModeI")
alg_modeI.set_header(["DDstPi0ModeIAlg/DDstPi0ModeI.h"])
         .set_constant({ "ECMS" => [:double, 4.26] })
         .set_alias({ "std::vector<double>" => "Vdouble" })
         .note(:multi_energy, "the analysis covers two energy points, 4.226 and 4.257 GeV; the ECMS constant is set to the nominal psi(4260) value and both data sets are processed with the same algorithm")

sel_modeI = Selection.new
  .select_track {                 # charged-track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     ">=2"               # at least two positive tracks
    nChrn     ">=1"               # at least one negative track
  }
  .select_photon {                # photon selection
    tdc_emc_start      0          # EMC timing window 0-14
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025      # E > 25 MeV in the barrel
    energyThreshold_e  0.050      # E > 50 MeV in the endcap
    nGam               ">=2"      # at least two photons
  }
  .pid(method: :probability) {    # PID: kaons vs (pi, p) and pions vs (K, p)
    prob_cut   0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # gamma-gamma -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"                                        # at least one pi0
  }
  .secondary_vertex_fit([:pip, :pim]) {               # K_S0 -> pi+pi- (used by the K_S0 D+ modes)
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:km, :pip, :pip, :kp, :pim, :pi0]) {  # D mass fit: D+ -> K-pi+pi+, anti-D0 -> K+pi-pi0
    invariant_mass_of(:km, :pip, :pip).constrain_to_nominal_mass_of(:D_plus)
    invariant_mass_of(:kp, :pim, :pi0).constrain_to_nominal_mass_of(:D0)
    chi2_cut 100
  }
  .kinematic_fit([:km, :pip, :pip, :kp, :pim, :pi0, :gamma, :gamma]) {  # nominal 2C fit
    nominal
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)             # pi0 window before the 2C fit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    miss_track_of(:pim)                                                # untagged pi- constrained to m(pi-)
    constrain_four_momentum
    chi2_cut 200
  }

# Inexpressible BOSS-side procedures
alg_modeI.note(:best_combination, "keep the best D+ anti-D0 combination by the smallest chi2_D+ + chi2_anti-D0")
         .note(:d_vertex, "a common D vertex chi2 < 100 is required for the D and anti-D daughters")
         .note(:dst_selection, "M(D pi0) > 2.1 GeV and |RM(D pi0) - m(D*)| < 36 MeV/c2 applied on the reconstructed D* candidates")
         .note(:ks0_modes, "the D+ -> K_S0 pi+, K_S0 pi+pi0 and K_S0 pi+pi+pi- channels are reconstructed through the K_S0 -> pi+pi- secondary vertex listed above")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ============================================================
# Mode II algorithm: e+e- -> D0 anti-D*0 pi0, anti-D*0 -> anti-D0 pi0
# ============================================================
alg_modeII = Algorithm.new("DDbarStPi0ModeII")
alg_modeII.set_header(["DDbarStPi0ModeIIAlg/DDbarStPi0ModeII.h"])
          .set_constant({ "ECMS" => [:double, 4.26] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .note(:multi_energy, "the analysis covers two energy points, 4.226 and 4.257 GeV; the ECMS constant is set to the nominal psi(4260) value and both data sets are processed with the same algorithm")

sel_modeII = Selection.new
  .select_track {                 # charged-track selection
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"               # at least one positive track
    nChrn     ">=1"               # at least one negative track
  }
  .select_photon {                # photon selection
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=2"
  }
  .pid(method: :probability) {    # PID: kaons vs (pi, p) and pions vs (K, p)
    prob_cut   0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {          # gamma-gamma -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kinematic_fit([:km, :pip, :kp, :pim]) {            # D mass fit: D0 -> K-pi+, anti-D0 -> K+pi-
    invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
    invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)
    chi2_cut 100
  }
  .kinematic_fit([:km, :pip, :kp, :pim, :gamma, :gamma]) {  # nominal 2C fit (signal: missing pi0)
    nominal
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)  # pi0 window before the 2C fit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    miss_track_of(:pi0)                                     # recoil mass constrained to m(pi0)
    constrain_four_momentum
    chi2_cut 60
  }
  .kinematic_fit([:km, :pip, :kp, :pim, :gamma, :gamma]) {  # competing anti-D*0 -> anti-D0 gamma hypothesis (chi2 stored)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    miss_track_of(:gamma)                                   # recoil mass constrained to 0
    constrain_four_momentum
  }

# Inexpressible BOSS-side procedures
alg_modeII.note(:best_combination, "keep the best D0 anti-D0 combination by the smallest chi2_D0 + chi2_anti-D0")
          .note(:d_vertex, "a common D vertex chi2 < 100 is required for the D0 and anti-D0 daughters")
          .note(:dst_selection, "M(D pi0) > 2.1 GeV and |RM(D pi0) - m(D*)| < 36 MeV/c2 applied on the reconstructed D* candidates")
          .note(:radiative_veto, "radiative veto applied in the ROOT stage: chi2_2C(pi0) < 60 and chi2_2C(gamma) > 20")
          .note(:ks0_modes, "the K_S0-involving D decay channels are reconstructed through a K_S0 -> pi+pi- secondary vertex")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ============================================================
# Execution
# ============================================================
files_modeI  = alg_modeI.execute_on([data_4226, data_4257, incMC_4226, incMC_4257] + exMC_modeI)
files_modeII = alg_modeII.execute_on([data_4226, data_4257, incMC_4226, incMC_4257] + exMC_modeII)