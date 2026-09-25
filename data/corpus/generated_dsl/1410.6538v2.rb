### Dataset preparation ###
# Nine c.m. energies from 4.210 to 4.420 GeV (XYZ scan points)
energy_points = %w[4210 4220 4230 4237 4246 4260 4270 4360 4420]

data_points  = energy_points.map { |e| DatasetManager.real_data.find("703_#{e}") }       # real data at all nine points
incMC_points = energy_points.map { |e| DatasetManager.inclusive_mc.find("703_#{e}") }   # inclusive MC at all nine points

# --- Decay card: e+e- -> omega chi_c0, chi_c0 -> pi+pi- / K+K- ---
decay_card_omega_chi_c0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega chi_c0 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay chi_c0
    0.500 pi+ pi- PHSP;
    0.500 K+ K- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: e+e- -> omega chi_c1, chi_c1 -> gamma J/psi, J/psi -> e+e- / mu+mu- ---
decay_card_omega_chi_c1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega chi_c1 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: e+e- -> omega chi_c2, chi_c2 -> gamma J/psi, J/psi -> e+e- / mu+mu- ---
decay_card_omega_chi_c2 = <<~DECAYCARD
    Decay psi(4260)
    1.000 omega chi_c2 PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive MC for each of the three modes at every energy point
exMCs_chi_c0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_chi_c0"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_omega_chi_c0
  config.cross_section = :default
end

exMCs_chi_c1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_chi_c1"
  config.events        = 200_000
  config.decay_card    = decay_card_omega_chi_c1
  config.cross_section = :default
end

exMCs_chi_c2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omega_chi_c2"
  config.events        = 200_000
  config.decay_card    = decay_card_omega_chi_c2
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ================= Mode: e+e- -> omega chi_c0 =================
alg_name_c0 = "OmegaChiC0"
alg_c0 = Algorithm.new(alg_name_c0)
alg_c0.set_header(["#{alg_name_c0}Alg/#{alg_name_c0}.h"])
      .set_constant({ "ECMS" => [:double, 4.260] })
      .set_alias({ "std::vector<double>" => "Vdouble" })

# Procedures that cannot be expressed in the selection DSL:
alg_c0.note(:momentum_assignment, "in the chi_c0 channel the two high-momentum tracks (|p| > 1 GeV/c) are "
                                  "assigned to the chi_c0 daughter pair and the two remaining low-momentum pions "
                                  "to the omega daughters; this kinematic assignment is made outside the "
                                  "kinematic-fit combination loop and has no DSL expression")
      .note(:mode_definition, "the chi_c0 -> pi+pi- and chi_c0 -> K+K- hypotheses are separated in the ROOT "
                              "analysis by comparing the two stored 5C chi2 values; the smaller chi2 defines the mode")

sel_c0 = Selection.new
sel_c0.select_track {              # charged track selection
        cos_theta 0.93             # |cos(theta)| < 0.93
        Vz        10.0             # |Vz| < 10 cm
        Vr        1.0              # Vr < 1 cm
        nChrp     "==2"            # exactly two positively charged tracks
        nChrn     "==2"            # exactly two negatively charged tracks
        nNet      "==0"            # net charge zero
      }
      .select_photon {             # photon selection
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025    # 25 MeV in the barrel
        energyThreshold_e 0.050    # 50 MeV in the endcap
        nGam              ">=2"    # at least two photons
      }
      .pid(method: :probability) { # probability PID: pions against K/p and kaons against pi/p
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        identify :kaon, against: [:pion, :proton]
        npip ">=2"
        nkp  ">=2"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {                      # reconstruct pi0 from photon pairs
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200                                                 # chi2 < 200
        npi0     ">=1"                                               # at least one pi0
      }
      # 5C fit, hypothesis I: omega pi+pi- -> pi0 pi+ pi- pi+ pi-  (4-momentum + M(gamma gamma) = m(pi0))
      .kinematic_fit([:pi0, :pip, :pim, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 100                                                 # chi2 < 100 for the pi pi hypothesis
      }
      # 5C fit, hypothesis II: omega K+K- -> pi0 pi+ pi- K+ K-  (competing hypothesis, chi2 stored for ROOT)
      .kinematic_fit([:pi0, :pip, :pim, :kp, :km]) {
        constrain_four_momentum
      }

alg_c0.with_decay_card(decay_card_omega_chi_c0).apply(sel_c0)
root_files_c0 = alg_c0.execute_on(data_points + incMC_points + exMCs_chi_c0)

# ================= Modes: e+e- -> omega chi_c1 / omega chi_c2 =================
# The chi_c1 and chi_c2 selections are identical
sel_c1_c2 = Selection.new
sel_c1_c2.select_track {           # same charged track requirements
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      }
      .select_photon {             # same photon requirements, at least three photons (2 from pi0 + 1 radiative)
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=3"
      }
      .pid(method: :probability) { # probability PID; tracks with p > 1 GeV/c treated as leptons
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6  # electron if EMC energy > 0.6 GeV, else muon
        identify :pion, against: [:kaon, :proton]
        npip "==1"
        npim "==1"
        nlp  "==1"
        nlm  "==1"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {                      # same pi0 reconstruction
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0     ">=1"
      }
      # 5C fit to pi0 gamma pi+ pi- l+ l- (4-momentum + M(gamma gamma) = m(pi0))
      .kinematic_fit([:pi0, :gamma, :pip, :pim, :lp, :lm]) {
        nominal
        constrain_four_momentum
        chi2_cut 60
      }

psi2S_veto_note = "psi(2S) veto: events with M(pi+pi-l+l-) or the recoil mass against pi+pi- in " \
                  "[3.68, 3.70] GeV/c^2 are rejected; the veto is applied after the final 5C fit in the ROOT analysis"

alg_name_c1 = "OmegaChiC1"
alg_c1 = Algorithm.new(alg_name_c1)
alg_c1.set_header(["#{alg_name_c1}Alg/#{alg_name_c1}.h"])
      .set_constant({ "ECMS" => [:double, 4.260] })
      .set_alias({ "std::vector<double>" => "Vdouble" })
      .note(:background_veto, psi2S_veto_note)
alg_c1.with_decay_card(decay_card_omega_chi_c1).apply(sel_c1_c2)
root_files_c1 = alg_c1.execute_on(data_points + incMC_points + exMCs_chi_c1)

alg_name_c2 = "OmegaChiC2"
alg_c2 = Algorithm.new(alg_name_c2)
alg_c2.set_header(["#{alg_name_c2}Alg/#{alg_name_c2}.h"])
      .set_constant({ "ECMS" => [:double, 4.260] })
      .set_alias({ "std::vector<double>" => "Vdouble" })
      .note(:background_veto, psi2S_veto_note)
alg_c2.with_decay_card(decay_card_omega_chi_c2).apply(sel_c1_c2.dup)
root_files_c2 = alg_c2.execute_on(data_points + incMC_points + exMCs_chi_c2)