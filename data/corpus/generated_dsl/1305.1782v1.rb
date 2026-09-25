# ============================================================================
# BOSS dataset preparation + event selection for the search of baryonic decays
# of psi(3770) / psi(4040):
#   psi -> Lambda Lambdabar pi+ pi-, Lambda Lambdabar pi0, Lambda Lambdabar eta,
#          Sigma+ Sigmabar-, Sigma0 Sigmabar0, Xi- Xibar+, Xi0 Xibar0
# with Lambda->p pi-, Lambdabar->pbar pi+, Sigma0->Lambda gamma,
#      Xi-->Lambda pi-, Xi0->Lambda pi0, pi0/eta->gamma gamma.
# ============================================================================

### --------------------------- Datasets ---------------------------
psi3770_data   = DatasetManager.real_data.find("712_3773")   # 2.9 fb^-1 @ 3.773 GeV (psi(3770))
psi3770_incMC  = DatasetManager.inclusive_mc.find("712_3773")
psi4009_data   = DatasetManager.real_data.find("703_4009")   # 482 pb^-1 @ 4.009 GeV (psi(4040))
psi4009_incMC  = DatasetManager.inclusive_mc.find("703_4009")
cont3554_data  = DatasetManager.real_data.find("712_3554")   # continuum @ 3.554 GeV
cont3650_data  = DatasetManager.real_data.find("709_3650")   # continuum @ 3.650 GeV
cont3650_incMC = DatasetManager.inclusive_mc.find("709_3650")

# All energy points over which the signal MC is generated (psi(3770), 4.009 GeV,
# and the two continuum points).
energy_points = [psi3770_data, psi4009_data, cont3554_data, cont3650_data]

# Datasets shared by every algorithm (real data + inclusive MC).
all_data = [psi3770_data, psi3770_incMC, psi4009_data, psi4009_incMC,
            cont3554_data, cont3650_data, cont3650_incMC]

### --------------------------- Decay cards ---------------------------
# Mode 1: psi -> Lambda Lambdabar pi+ pi-
decay_card_LLpipi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0 pi+ pi- PHSP;
    Enddecay
    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay
    End
DECAYCARD

# Mode 2: psi -> Lambda Lambdabar pi0
decay_card_LLpi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
    Enddecay
    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 3: psi -> Lambda Lambdabar eta
decay_card_LLeta = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0 eta PHSP;
    Enddecay
    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay
    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 4: psi -> Sigma+ anti-Sigma-
decay_card_SpSm = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Sigma+ anti-Sigma- PHSP;
    Enddecay
    Decay Sigma+
    1.0000 p+ pi0 PHSP;
    Enddecay
    Decay anti-Sigma-
    1.0000 anti-p- pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode 5: psi -> Sigma0 anti-Sigma0
decay_card_S0S0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Sigma0 anti-Sigma0 PHSP;
    Enddecay
    Decay Sigma0
    1.0000 Lambda0 gamma PHSP;
    Enddecay
    Decay anti-Sigma0
    1.0000 anti-Lambda0 gamma PHSP;
    Enddecay
    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay
    End
DECAYCARD

# Mode 6: psi -> Xi- anti-Xi+
decay_card_XiXi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Xi- anti-Xi+ PHSP;
    Enddecay
    Decay Xi-
    1.0000 Lambda0 pi- PHSP;
    Enddecay
    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ PHSP;
    Enddecay
    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay
    End
DECAYCARD

# Mode 7: psi -> Xi0 anti-Xi0
decay_card_Xi0Xi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay
    Decay Xi0
    1.0000 Lambda0 pi0 PHSP;
    Enddecay
    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay
    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

### --------------------------- Exclusive MC (50k events) ---------------------------
# Signal MC is generated for every mode at psi(3770), at 4.009 GeV and at the
# two continuum points (one ExclusiveMC per energy point).
exMC_LLpipi = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_LLpipi"
  config.events        = 50000
  config.decay_card    = decay_card_LLpipi
  config.cross_section = :default
end

exMC_LLpi0 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_LLpi0"
  config.events        = 50000
  config.decay_card    = decay_card_LLpi0
  config.cross_section = :default
end

exMC_LLeta = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_LLeta"
  config.events        = 50000
  config.decay_card    = decay_card_LLeta
  config.cross_section = :default
end

exMC_SpSm = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_SpSm"
  config.events        = 50000
  config.decay_card    = decay_card_SpSm
  config.cross_section = :default
end

exMC_S0S0 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_S0S0"
  config.events        = 50000
  config.decay_card    = decay_card_S0S0
  config.cross_section = :default
end

exMC_XiXi = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_XiXi"
  config.events        = 50000
  config.decay_card    = decay_card_XiXi
  config.cross_section = :default
end

exMC_Xi0Xi0 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_Xi0Xi0"
  config.events        = 50000
  config.decay_card    = decay_card_Xi0Xi0
  config.cross_section = :default
end

### ============================================================
### Event selection (BOSS) -- one Algorithm per decay mode (Rule T1)
### ============================================================

# ---------------------------------------------------------------------------
# Mode 1: psi -> Lambda Lambdabar pi+ pi-   (fit: Lambda Lambdabar pi+ pi-)
# ---------------------------------------------------------------------------
alg_name1 = "PsiBaryonLLpipi"
alg1 = Algorithm.new(alg_name1)
alg1.set_header(["#{alg_name1}Alg/#{alg_name1}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel1 = Selection.new
sel1.select_track {                    # charged tracks
        cos_theta 0.93                 # |cos(theta)| < 0.93
        Vz        100.0                # |Vz| < 100 cm
        Vr        10.0                 # Vr < 10 mm
        nNet      "==0"                # net charge zero
        nChrp     ">=3"                # >= 3 positive tracks (p pi+ pi+)
        nChrn     ">=3"                # >= 3 negative tracks (pbar pi- pi-)
    }
    .pid(method: :probability) {       # PID: probability method, CL > 0.001
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # proton/anti-proton vs K, pi
        identify :pion,   against: [:kaon, :proton] # pion/anti-pion vs K, p
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }   # drop p    with pt < 300 MeV/c
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }   # drop pbar with pt < 300 MeV/c
    .secondary_vertex_fit([:prp, :pim]) {             # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {             # Lambdabar -> pbar pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) {   # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200                                  # loose; tight <60 in ROOT
    }
alg1.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Lambda) in [1.107, 1.124] GeV/c^2 and |M(p pi-) - M(Lambda)| < 40 MeV/c^2")
    .note(:no_primary_vertex_requirement,
          "No primary-vertex (IP) requirement is imposed on the Lambda/Lambdabar " \
          "daughter tracks; only the secondary-vertex fit is used.")
    .with_decay_card(decay_card_LLpipi)
    .apply(sel1)
root_files_1 = alg1.execute_on(all_data + exMC_LLpipi)

# ---------------------------------------------------------------------------
# Mode 2: psi -> Lambda Lambdabar pi0   (fit: Lambda Lambdabar pi0)
# ---------------------------------------------------------------------------
alg_name2 = "PsiBaryonLLpi0"
alg2 = Algorithm.new(alg_name2)
alg2.set_header(["#{alg_name2}Alg/#{alg_name2}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel2 = Selection.new
sel2.select_track {
        cos_theta 0.93
        Vz        100.0
        Vr        10.0
        nNet      "==0"
        nChrp     ">=1"                # >= 1 positive track (p pi+)
        nChrn     ">=1"                # >= 1 negative track (pbar pi-)
    }
    .select_photon {                   # photons
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0         # > 10 deg from any charged track
        energyThreshold_b 0.025        # 25 MeV barrel
        energyThreshold_e 0.050        # 50 MeV endcap
        nGam              ">=2"        # at least two photons (pi0 -> gamma gamma)
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }
    .kalman_kinematic_fit([:gamma, :gamma]) {          # pi0 -> gamma gamma (1C)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
    }
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) {     # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg2.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Lambda) in [1.107, 1.124] GeV/c^2, M(pi0) in [0.115, 0.150] GeV/c^2, " \
          "|M(p pi-) - M(Lambda)| < 40 MeV/c^2")
    .note(:no_primary_vertex_requirement,
          "No primary-vertex (IP) requirement on the Lambda/Lambdabar daughter tracks.")
    .with_decay_card(decay_card_LLpi0)
    .apply(sel2)
root_files_2 = alg2.execute_on(all_data + exMC_LLpi0)

# ---------------------------------------------------------------------------
# Mode 3: psi -> Lambda Lambdabar eta   (fit: Lambda Lambdabar eta)
# ---------------------------------------------------------------------------
alg_name3 = "PsiBaryonLLeta"
alg3 = Algorithm.new(alg_name3)
alg3.set_header(["#{alg_name3}Alg/#{alg_name3}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel3 = Selection.new
sel3.select_track {
        cos_theta 0.93
        Vz        100.0
        Vr        10.0
        nNet      "==0"
        nChrp     ">=1"
        nChrn     ">=1"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"        # eta -> gamma gamma
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }
    .kalman_kinematic_fit([:gamma, :gamma]) {          # eta -> gamma gamma (1C)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
    }
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar, :eta]) {     # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg3.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Lambda) in [1.107, 1.124] GeV/c^2, M(eta) in [0.515, 0.569] GeV/c^2, " \
          "|M(p pi-) - M(Lambda)| < 40 MeV/c^2")
    .note(:no_primary_vertex_requirement,
          "No primary-vertex (IP) requirement on the Lambda/Lambdabar daughter tracks.")
    .with_decay_card(decay_card_LLeta)
    .apply(sel3)
root_files_3 = alg3.execute_on(all_data + exMC_LLeta)

# ---------------------------------------------------------------------------
# Mode 4: psi -> Sigma+ anti-Sigma-   (fit: p pbar pi0 pi0)  -- two pi0's
# ---------------------------------------------------------------------------
alg_name4 = "PsiBaryonSpSm"
alg4 = Algorithm.new(alg_name4)
alg4.set_header(["#{alg_name4}Alg/#{alg_name4}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel4 = Selection.new
sel4.select_track {
        cos_theta 0.93
        Vz        100.0
        Vr        10.0
        nNet      "==0"
        nChrp     ">=1"                # >= 1 positive track (p)
        nChrn     ">=1"                # >= 1 negative track (pbar)
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"        # two pi0 -> four photons
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # only p / pbar are charged
        nprp ">=1"
        nprm ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }
    .kalman_kinematic_fit([:gamma, :gamma]) {          # pi0 pi0 (two mass constraints)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"
    }
    .kinematic_fit([:prp, :prm, :pi0, :pi0]) {         # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg4.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Sigma+) in [1.164, 1.206] GeV/c^2 and M(pi0) in [0.115, 0.150] GeV/c^2")
    .with_decay_card(decay_card_SpSm)
    .apply(sel4)
root_files_4 = alg4.execute_on(all_data + exMC_SpSm)

# ---------------------------------------------------------------------------
# Mode 5: psi -> Sigma0 anti-Sigma0   (fit: Lambda Lambdabar gamma gamma)
# ---------------------------------------------------------------------------
alg_name5 = "PsiBaryonS0S0"
alg5 = Algorithm.new(alg_name5)
alg5.set_header(["#{alg_name5}Alg/#{alg_name5}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel5 = Selection.new
sel5.select_track {
        cos_theta 0.93
        Vz        100.0
        Vr        10.0
        nNet      "==0"
        nChrp     ">=1"                # >= 1 positive track (p pi+)
        nChrn     ">=1"                # >= 1 negative track (pbar pi-)
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"        # one photon from each Sigma0
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }
    .secondary_vertex_fit([:prp, :pim]) {              # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {              # Lambdabar -> pbar pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {   # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg5.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Lambda) in [1.107, 1.124] GeV/c^2, M(Sigma0) in [1.178, 1.205] GeV/c^2, " \
          "|M(p pi-) - M(Lambda)| < 40 MeV/c^2")
    .note(:no_primary_vertex_requirement,
          "No primary-vertex (IP) requirement on the Lambda/Lambdabar daughter tracks.")
    .with_decay_card(decay_card_S0S0)
    .apply(sel5)
root_files_5 = alg5.execute_on(all_data + exMC_S0S0)

# ---------------------------------------------------------------------------
# Mode 6: psi -> Xi- anti-Xi+   (fit: Lambda pi- Lambdabar pi+)
# ---------------------------------------------------------------------------
alg_name6 = "PsiBaryonXiXi"
alg6 = Algorithm.new(alg_name6)
alg6.set_header(["#{alg_name6}Alg/#{alg_name6}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel6 = Selection.new
sel6.select_track {
        cos_theta 0.93
        Vz        100.0
        Vr        10.0
        nNet      "==0"
        nChrp     ">=3"                # >= 3 positive tracks (p pi+ pi+)
        nChrn     ">=3"                # >= 3 negative tracks (pbar pi- pi-)
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }
    .secondary_vertex_fit([:prp, :pim]) {              # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {              # Lambdabar -> pbar pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :pim, :Lambda_bar, :pip]) {   # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg6.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Lambda) in [1.107, 1.124] GeV/c^2, M(Xi-) in [1.305, 1.337] GeV/c^2, " \
          "|M(p pi-) - M(Lambda)| < 40 MeV/c^2")
    .note(:no_primary_vertex_requirement,
          "No primary-vertex (IP) requirement on the Lambda/Lambdabar daughter tracks.")
    .with_decay_card(decay_card_XiXi)
    .apply(sel6)
root_files_6 = alg6.execute_on(all_data + exMC_XiXi)

# ---------------------------------------------------------------------------
# Mode 7: psi -> Xi0 anti-Xi0   (fit: Lambda Lambdabar pi0 pi0)  -- two pi0's
# ---------------------------------------------------------------------------
alg_name7 = "PsiBaryonXi0Xi0"
alg7 = Algorithm.new(alg_name7)
alg7.set_header(["#{alg_name7}Alg/#{alg_name7}.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
sel7 = Selection.new
sel7.select_track {
        cos_theta 0.93
        Vz        100.0
        Vr        10.0
        nNet      "==0"
        nChrp     ">=1"                # >= 1 positive track (p pi+)
        nChrn     ">=1"                # >= 1 negative track (pbar pi-)
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"        # two pi0 -> four photons
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .remove(:prp) { condition "pt_of(:prp) < 0.3" }
    .remove(:prm) { condition "pt_of(:prm) < 0.3" }
    .kalman_kinematic_fit([:gamma, :gamma]) {          # pi0 pi0 (two mass constraints)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"
    }
    .secondary_vertex_fit([:prp, :pim]) {              # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {              # Lambdabar -> pbar pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {   # 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg7.note(:final_mass_windows,
          "Final mass windows applied in ROOT after the 4C fit: " \
          "M(Lambda) in [1.107, 1.124] GeV/c^2, M(pi0) in [0.115, 0.150] GeV/c^2, " \
          "M(Xi0) in [1.281, 1.330] GeV/c^2, |M(p pi-) - M(Lambda)| < 40 MeV/c^2")
    .note(:no_primary_vertex_requirement,
          "No primary-vertex (IP) requirement on the Lambda/Lambdabar daughter tracks.")
    .with_decay_card(decay_card_Xi0Xi0)
    .apply(sel7)
root_files_7 = alg7.execute_on(all_data + exMC_Xi0Xi0)