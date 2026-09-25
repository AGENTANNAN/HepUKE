# =============================================================================
# psi(3686) -> gamma chi_cJ (J = 0, 1, 2), with
#     chi_cJ -> eta_c pi+pi-   and   chi_cJ -> K Kbar pi pi pi
# searched in two final states:
#   Mode A : gamma K_S0 K+- pi-+ pi+ pi-   (also via chi_cJ -> eta_c pi+pi-,
#            eta_c -> K_S0 K+- pi-+)
#   Mode B : gamma K+ K- pi+ pi- pi0       (also via chi_cJ -> eta_c pi+pi-,
#            eta_c -> K+ K- pi0)
# All three chi_cJ states share one selection chain per mode.
# =============================================================================

### ------------------------------- Datasets ------------------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data, 1.06e8 events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")     # 3.65 GeV continuum data (~42 pb^-1)
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")  # 3.65 GeV continuum inclusive MC

### ------------------ Decay cards used by the two algorithms -------------- ###
# Mode A : all three chi_cJ states, direct (K_S0 K+ pi- pi+ pi-) and eta_c pi+pi- channels
decay_card_modeA = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  1.000 gamma chi_c1 PHSP;
  1.000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c0
  0.500 K_S0 K+ pi- pi+ pi- PHSP;
  0.500 eta_c pi+ pi- PHSP;
  Enddecay
  Decay chi_c1
  0.500 K_S0 K+ pi- pi+ pi- PHSP;
  0.500 eta_c pi+ pi- PHSP;
  Enddecay
  Decay chi_c2
  0.500 K_S0 K+ pi- pi+ pi- PHSP;
  0.500 eta_c pi+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay eta_c
  0.500 K_S0 K+ pi- PHSP;
  0.500 K_S0 K- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# Mode B : all three chi_cJ states, direct (K+ K- pi+ pi- pi0) and eta_c pi+pi- channels
decay_card_modeB = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  1.000 gamma chi_c1 PHSP;
  1.000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c0
  0.500 K+ K- pi+ pi- pi0 PHSP;
  0.500 eta_c pi+ pi- PHSP;
  Enddecay
  Decay chi_c1
  0.500 K+ K- pi+ pi- pi0 PHSP;
  0.500 eta_c pi+ pi- PHSP;
  Enddecay
  Decay chi_c2
  0.500 K+ K- pi+ pi- pi0 PHSP;
  0.500 eta_c pi+ pi- PHSP;
  Enddecay

  Decay eta_c
  1.000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### --------- Decay cards for the 12 signal samples (3 chi_cJ x 4 states) --- ###
# ---- Mode A, chi_cJ -> K_S0 K+ pi- pi+ pi- ----
dc_c0_A_dir = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.000 K_S0 K+ pi- pi+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

dc_c1_A_dir = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.000 K_S0 K+ pi- pi+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

dc_c2_A_dir = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.000 K_S0 K+ pi- pi+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# ---- Mode A via chi_cJ -> eta_c pi+ pi-, eta_c -> K_S0 K+- pi-+ ----
dc_c0_A_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.000 eta_c pi+ pi- PHSP;
  Enddecay
  Decay eta_c
  0.500 K_S0 K+ pi- PHSP;
  0.500 K_S0 K- pi+ PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

dc_c1_A_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.000 eta_c pi+ pi- PHSP;
  Enddecay
  Decay eta_c
  0.500 K_S0 K+ pi- PHSP;
  0.500 K_S0 K- pi+ PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

dc_c2_A_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.000 eta_c pi+ pi- PHSP;
  Enddecay
  Decay eta_c
  0.500 K_S0 K+ pi- PHSP;
  0.500 K_S0 K- pi+ PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# ---- Mode B, chi_cJ -> K+ K- pi+ pi- pi0 ----
dc_c0_B_dir = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.000 K+ K- pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_c1_B_dir = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.000 K+ K- pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_c2_B_dir = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.000 K+ K- pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ---- Mode B via chi_cJ -> eta_c pi+ pi-, eta_c -> K+ K- pi0 ----
dc_c0_B_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.000 eta_c pi+ pi- PHSP;
  Enddecay
  Decay eta_c
  1.000 K+ K- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_c1_B_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.000 eta_c pi+ pi- PHSP;
  Enddecay
  Decay eta_c
  1.000 K+ K- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_c2_B_etac = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.000 eta_c pi+ pi- PHSP;
  Enddecay
  Decay eta_c
  1.000 K+ K- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

### --------------- Exclusive MC: 500k events per chi_cJ / final state ------ ###
modeA_mc_specs = {
  "exmc_chic0_A_KSKpipipi" => dc_c0_A_dir,
  "exmc_chic1_A_KSKpipipi" => dc_c1_A_dir,
  "exmc_chic2_A_KSKpipipi" => dc_c2_A_dir,
  "exmc_chic0_A_etacpipi_KSKpi" => dc_c0_A_etac,
  "exmc_chic1_A_etacpipi_KSKpi" => dc_c1_A_etac,
  "exmc_chic2_A_etacpipi_KSKpi" => dc_c2_A_etac
}
modeB_mc_specs = {
  "exmc_chic0_B_KKpipipi0" => dc_c0_B_dir,
  "exmc_chic1_B_KKpipipi0" => dc_c1_B_dir,
  "exmc_chic2_B_KKpipipi0" => dc_c2_B_dir,
  "exmc_chic0_B_etacpipi_KKpi0" => dc_c0_B_etac,
  "exmc_chic1_B_etacpipi_KKpi0" => dc_c1_B_etac,
  "exmc_chic2_B_etacpipi_KKpi0" => dc_c2_B_etac
}

make_exclusive_mc = lambda do |specs|
  specs.map do |sample_name, card|
    DatasetManager.create_exclusive_mc do |config|
      config.sample_name     = sample_name
      config.related_dataset = psip_data   # psi(3686) real data
      config.events          = 500_000
      config.decay_card      = card
      config.cross_section   = :default
    end
  end
end

exMC_modeA = make_exclusive_mc.call(modeA_mc_specs)
exMC_modeB = make_exclusive_mc.call(modeB_mc_specs)

### --------------------------- Event selection --------------------------- ###
# ---------------------------------------------------------------------------
# Mode A : gamma K_S0 K+- pi-+ pi+ pi-   (3 positive + 3 negative tracks)
# ---------------------------------------------------------------------------
alg_name_A = "ChicJToKSKpipipi"
alg_modeA = Algorithm.new(alg_name_A)
alg_modeA.set_header(["#{alg_name_A}Alg/#{alg_name_A}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         # K_S0 mass window (|M(pi+pi-) - m_K_S0| < 7 MeV) and vertex displacement
         # (>= 0.5 cm from the IP) have no dedicated DSL method - recorded here.
         .note(:ks0_selection, "K_S0 candidates from the pi+pi- secondary-vertex fit are
           required to have an invariant mass within 7 MeV of the nominal K_S0 mass and a
           decay vertex at least 0.5 cm away from the interaction point; windows taken from
           K_S0 control samples in data")

sel_modeA = Selection.new
sel_modeA.select_track {
           cos_theta 0.93      # |cos(theta)| < 0.93
           Vz        10.0      # |Vz| < 10 cm
           Vr        1.0       # Vr < 1 cm
           nChrp     "==3"     # three positive tracks
           nChrn     "==3"     # three negative tracks
           nNet      "==0"     # net charge zero
         }
         .select_photon {
           tdc_emc_start     0      # EMC timing window 0-14 (700 ns)
           tdc_emc_end       14
           angle_to_track    20.0   # at least 20 degrees from any charged track
           energyThreshold_b 0.025  # E > 25 MeV (barrel)
           energyThreshold_e 0.050  # E > 50 MeV (endcap)
           nGam              ">=1"  # at least one photon
         }
         .pid(method: :probability) {
           prob_cut 0.001                                    # PID probability > 0.001
           identify :kaon, against: [:pion, :proton]         # pi/K/p separation (K+ and K-)
           identify :pion, against: [:kaon, :proton]         # pi+ and pi-
           nkp   ">=1"   # at least one charged kaon
           npip  ">=2"   # at least two pi+
           npim  ">=2"   # at least two pi-
         }
         # secondary vertex fit forming K_S0 from pi+pi- pairs
         .secondary_vertex_fit([:pip, :pim]) {
           build_virtual_particle(:K_S0).by_minimizing_mass_difference
           remove_used_particle_from_candidate_list
         }
         # 4C kinematic fit of gamma K_S0 K+- pi-+ pi+ pi- to the initial e+e- four-momentum.
         # Only the loose BOSS cut chi2 < 200 is applied here; the paper's chi2 < 50 is
         # applied downstream in the ROOT analysis.
         .kinematic_fit([:gamma, :K_S0, :kp, :km, :pip, :pim]) {
           nominal
           constrain_four_momentum
           chi2_cut 200
         }
         # Competing gamma gamma hypothesis (one extra photon): no nominal, no chi2 cut -
         # its chi2 is stored for the extra-photon veto applied in the ROOT analysis.
         .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :km, :pip, :pim]) {
           constrain_four_momentum
         }

alg_modeA.with_decay_card(decay_card_modeA).apply(sel_modeA)
alg_modeA.execute_on([psip_data, psip_incMC, cont_data, cont_incMC] + exMC_modeA)

# ---------------------------------------------------------------------------
# Mode B : gamma K+ K- pi+ pi- pi0   (2 positive + 2 negative tracks)
# ---------------------------------------------------------------------------
alg_name_B = "ChicJToKKpipipi0"
alg_modeB = Algorithm.new(alg_name_B)
alg_modeB.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_modeB = Selection.new
sel_modeB.select_track {
           cos_theta 0.93      # |cos(theta)| < 0.93
           Vz        10.0      # |Vz| < 10 cm
           Vr        1.0       # Vr < 1 cm
           nChrp     "==2"     # two positive tracks
           nChrn     "==2"     # two negative tracks
           nNet      "==0"     # net charge zero
         }
         .select_photon {
           tdc_emc_start     0      # EMC timing window 0-14 (700 ns)
           tdc_emc_end       14
           angle_to_track    20.0   # at least 20 degrees from any charged track
           energyThreshold_b 0.025  # E > 25 MeV (barrel)
           energyThreshold_e 0.050  # E > 50 MeV (endcap)
           nGam              ">=3"  # at least three photons (one prompt gamma + pi0 -> gamma gamma)
         }
         .pid(method: :probability) {
           prob_cut 0.001                                    # PID probability > 0.001
           identify :kaon, against: [:pion, :proton]         # pi/K/p separation (K+ and K-)
           identify :pion, against: [:kaon, :proton]         # pi+ and pi-
           nkp   "==1"   # exactly one K+
           nkm   "==1"   # exactly one K-
           npip  "==1"   # exactly one pi+
           npim  "==1"   # exactly one pi-
         }
         # pi0 reconstruction from gamma gamma via a Kalman (1C) fit, chi2 < 25
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).within(0.120, 0.145)              # 0.120 < M(gamma gamma) < 0.145 GeV
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25
           npi0 ">=1"
         }
         # 4C kinematic fit of gamma K+ K- pi+ pi- pi0 to the initial e+e- four-momentum.
         # Loose BOSS cut chi2 < 200; the paper's chi2 < 50 is applied downstream in ROOT.
         .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) {
           nominal
           constrain_four_momentum
           chi2_cut 200
         }
         # Competing gamma gamma hypothesis (one extra photon): chi2 stored for the
         # extra-photon veto applied in the ROOT analysis.
         .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim, :pi0]) {
           constrain_four_momentum
         }

alg_modeB.with_decay_card(decay_card_modeB).apply(sel_modeB)
alg_modeB.execute_on([psip_data, psip_incMC, cont_data, cont_incMC] + exMC_modeB)