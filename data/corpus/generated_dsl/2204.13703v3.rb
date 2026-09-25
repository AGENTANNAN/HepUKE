# ============================================================
# Dataset preparation and event selection (BOSS part)
#   e+e- -> K_S0 D_s+ D*-    (D_s+ tag, D*- undetected)
#   e+e- -> K_S0 D_s*+ D-    (D-  tag, D_s*+ undetected)
# at five energy points: 4.628, 4.641, 4.661, 4.682, 4.699 GeV
# ============================================================

### Dataset description ###
data_4628 = DatasetManager.real_data.find("706_4620")      # 4.628 GeV
data_4641 = DatasetManager.real_data.find("706_4640")      # 4.641 GeV
data_4661 = DatasetManager.real_data.find("706_4660")      # 4.661 GeV
data_4682 = DatasetManager.real_data.find("706_4680")      # 4.682 GeV
data_4699 = DatasetManager.real_data.find("706_4700")      # 4.699 GeV

incMC_4628 = DatasetManager.inclusive_mc.find("706_4620")  # inclusive MC at each point
incMC_4641 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4661 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4699 = DatasetManager.inclusive_mc.find("706_4700")

all_data  = [data_4628, data_4641, data_4661, data_4682, data_4699]
all_incMC = [incMC_4628, incMC_4641, incMC_4661, incMC_4682, incMC_4699]

# Decay card for e+e- -> K_S0 D_s+ D*-  (D*- left undetected), D_s+ tag modes
decay_card_Ds_tag = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 D_s+ D*-      PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-            PHSP;
    Enddecay

    Decay D_s+
    0.200 K+ K- pi+          PHSP;
    0.200 K_S0 K+            PHSP;
    0.200 K+ K- pi+ pi0      PHSP;
    0.200 K_S0 K+ pi+ pi-    PHSP;
    0.200 eta' pi+           PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta        PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma        PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma        PHSP;
    Enddecay

    Decay D*-
    1.000 anti-D0 pi-        PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi-             PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for e+e- -> K_S0 D_s*+ D-  (D_s*+ left undetected), D- tag modes
decay_card_D_tag = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 D_s*+ D-      PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-            PHSP;
    Enddecay

    Decay D_s*+
    1.000 D_s+ gamma         PHSP;
    Enddecay

    Decay D_s+
    1.000 K+ K- pi+          PHSP;
    Enddecay

    Decay D-
    0.334 K+ pi- pi-         PHSP;
    0.333 K_S0 pi-           PHSP;
    0.333 K_S0 pi+ pi- pi-   PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 500k events each, generated at 4.682 GeV only
exMC_Ds_tag = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4682_KS0DsDst"
  config.related_dataset = data_4682
  config.events          = 500_000
  config.decay_card      = decay_card_Ds_tag
  config.cross_section   = :default
end

exMC_D_tag = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4682_KS0DspD"
  config.related_dataset = data_4682
  config.events          = 500_000
  config.decay_card      = decay_card_D_tag
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------
# Algorithm I : K_S0 D_s+ (D_s+ tag), D*- undetected
# ------------------------------------------------------------
alg_name_Ds = "KS0DsTag"
alg_Ds = Algorithm.new(alg_name_Ds)
alg_Ds.set_header(["#{alg_name_Ds}Alg/#{alg_name_Ds}.h"])
      .set_constant({"ECMS" => [:double, 4.682]})          # representative energy
      .set_alias({"std::vector<double>" => "Vdouble"})

selection_Ds = Selection.new
  .select_track {                       # charged tracks
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nTot      ">=4"                     # at least four tracks in total
  }
  .select_photon {                      # photons
    tdc_emc_start     0                 # EMC timing window [0, 14]
    tdc_emc_end       14
    angle_to_track    10.0              # >= 10 deg from any charged track
    energyThreshold_b 0.025             # 25 MeV (barrel)
    energyThreshold_e 0.050             # 50 MeV (endcap)
  }                                     # no minimum photon multiplicity
  .pid(method: :chi2_sum) {             # chi2-sum combinatorial PID
    chi_min_cut 0.001                   # chi2_min > 0.001
    identify :kaon, :pion               # K+/K- and pi+/pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 10
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # eta -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 10
  }
  .secondary_vertex_fit([:pip, :pim]) {       # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .partial_miss([3]) {                        # reconstruct K_S0 + D_s+, miss D*-
    best_combination_by_mass :K_S0,  0.49761  # best candidate by closest mass
    best_combination_by_mass :D_s+,  1.96835
    require_recoil_mass 1.990, 2.030          # |RQ(K_S0 D_s+) - m_D*-| < 20 MeV
  }

alg_Ds
  .note(:background_veto, "K_S0 -> pi+pi- candidates require vertex-fit chi2 < 100, a second IP-pointing vertex fit with chi2 < 40, decay length > 2 sigma and |M(pi+pi-) - m_K_S0| < 11 MeV; D_s+ candidates are kept within 15 MeV of the nominal mass (best candidate by closest mass); resonance vetoes M(K+K-) < 1.05 GeV (phi), |M(K+pi-) - m_K*(892)| < 70 MeV and |M(pi-pi0) - m_rho| < 150 MeV")
  .note(:efficiency_curve, "the recoil-mass signal window uses the corrected variable RQ = RM + M(D_s+) - m(D_s+), where RM = P4_cms - p4(K_S0 D_s+), M(D_s+) is the reconstructed and m(D_s+) the nominal D_s+ mass; D_s+ sidebands (1.895, 1.935) and (1.995, 2.035) GeV are defined for background estimation")
  .note(:beam_energy_scan, "the identical selection is applied at the five energy points 4.628, 4.641, 4.661, 4.682 and 4.699 GeV; only the 4.682 GeV exclusive MC is generated, so the CMS four-vector used by the partial reconstruction must be taken per run from the beam energy of the point being processed")
  .with_decay_card(decay_card_Ds_tag)
  .apply(selection_Ds)

# ------------------------------------------------------------
# Algorithm II : K_S0 D- (D- tag), D_s*+ undetected
# ------------------------------------------------------------
alg_name_D = "KS0DTag"
alg_D = Algorithm.new(alg_name_D)
alg_D.set_header(["#{alg_name_D}Alg/#{alg_name_D}.h"])
     .set_constant({"ECMS" => [:double, 4.682]})           # representative energy
     .set_alias({"std::vector<double>" => "Vdouble"})

selection_D = Selection.new
  .select_track {                       # charged tracks
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      ">=4"                     # at least four tracks in total
  }
  .select_photon {                      # photons
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :chi2_sum) {             # chi2-sum combinatorial PID
    chi_min_cut 0.001
    identify :kaon, :pion
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 10
  }
  .secondary_vertex_fit([:pip, :pim]) {       # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .partial_miss([2]) {                        # reconstruct K_S0 + D-, miss D_s*+
    best_combination_by_mass :K_S0, 0.49761   # best candidate by closest mass
    best_combination_by_mass :Dm,   1.86966
    require_recoil_mass 2.102, 2.122          # |RQ(K_S0 D-) - m_D_s*+| < 10 MeV
  }

alg_D
  .note(:background_veto, "K_S0 -> pi+pi- candidates require vertex-fit chi2 < 100, a second IP-pointing vertex fit with chi2 < 40, decay length > 2 sigma and |M(pi+pi-) - m_K_S0| < 11 MeV; D- candidates are kept within 15 MeV of the nominal mass (best candidate by closest mass); resonance vetoes M(K+K-) < 1.05 GeV (phi), |M(K+pi-) - m_K*(892)| < 70 MeV and |M(pi-pi0) - m_rho| < 150 MeV; a K_S0 veto is applied for the D- -> K_S0 pi+ pi- pi- mode")
  .note(:efficiency_curve, "the recoil-mass signal window uses the corrected variable RQ = RM + M(D-) - m(D-), where RM = P4_cms - p4(K_S0 D-), M(D-) is the reconstructed and m(D-) the nominal D- mass; D- sidebands (1.800, 1.840) and (1.900, 1.940) GeV are defined for background estimation")
  .note(:beam_energy_scan, "the identical selection is applied at the five energy points 4.628, 4.641, 4.661, 4.682 and 4.699 GeV; only the 4.682 GeV exclusive MC is generated, so the CMS four-vector used by the partial reconstruction must be taken per run from the beam energy of the point being processed")
  .with_decay_card(decay_card_D_tag)
  .apply(selection_D)

### Execute on real data, inclusive MC (all five points) and exclusive MC ###
root_files_Ds = alg_Ds.execute_on(all_data + all_incMC + [exMC_Ds_tag])
root_files_D  = alg_D.execute_on(all_data + all_incMC  + [exMC_D_tag])