# =============================================================================
# e+e- -> eta_c eta pi+pi- at 4.230, 4.260, 4.360, 4.420 and 4.600 GeV  (BOSS part)
# eta_c is reconstructed in its 16 hadronic decay modes; each mode needs its own
# Algorithm + Selection chain (Rule T1).  Dataset preparation + selection only,
# up to and including the kinematic fit.
# =============================================================================

### -------------------------------- Datasets ------------------------------- ###
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")
data_scan = [data_4230, data_4260, data_4360, data_4420, data_4600]

incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_scan = [incMC_4230, incMC_4260, incMC_4360, incMC_4420, incMC_4600]

### ------------------------- Decay cards (EvtGen) -------------------------- ###
# The description gives no intermediate resonance, so psi(4260) is used as the
# KKMC top mother.  The recoiling eta always decays to gamma gamma.

# --- 1. eta_c -> pi+ pi- K+ K- ---
dc_pipiKK = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- K+ K- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 2. eta_c -> 2(K+ K-) ---
dc_2KK = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- K+ K- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 3. eta_c -> 2(pi+ pi-) ---
dc_2pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 4. eta_c -> 3(pi+ pi-) ---
dc_3pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- pi+ pi- pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 5. eta_c -> K_S0 K+- pi-+ ---
dc_KSKpi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 6. eta_c -> K+ K- pi0 ---
dc_KKpi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 7. eta_c -> K+ K- eta ---
dc_KKeta = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 8. eta_c -> p+ anti-p- ---
dc_ppbar = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 9. eta_c -> pi+ pi- eta ---
dc_pipieta = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 10. eta_c -> pi+ pi- pi0 pi0 ---
dc_pipipi0pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi+ pi- pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 11. eta_c -> p+ anti-p- pi0 ---
dc_ppbarpi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 12. eta_c -> p+ anti-p- pi+ pi- ---
dc_ppbarpipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta_c eta pi+ pi- PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### ------------------------- Exclusive MC samples -------------------------- ###
# Same signal MC has to be produced at all five energy points -> create_exclusive_mc_for
# returns one ExclusiveMC per energy point (200,000 events each).

exMC_pipiKK = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_pipiKK"
  config.events        = 200_000
  config.decay_card    = dc_pipiKK
  config.cross_section = :default
end

exMC_2KK = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_2KK"
  config.events        = 200_000
  config.decay_card    = dc_2KK
  config.cross_section = :default
end

exMC_2pipi = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_2pipi"
  config.events        = 200_000
  config.decay_card    = dc_2pipi
  config.cross_section = :default
end

exMC_3pipi = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_3pipi"
  config.events        = 200_000
  config.decay_card    = dc_3pipi
  config.cross_section = :default
end

exMC_KSKpi = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_KSKpi"
  config.events        = 200_000
  config.decay_card    = dc_KSKpi
  config.cross_section = :default
end

exMC_KKpi0 = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_KKpi0"
  config.events        = 200_000
  config.decay_card    = dc_KKpi0
  config.cross_section = :default
end

exMC_KKeta = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_KKeta"
  config.events        = 200_000
  config.decay_card    = dc_KKeta
  config.cross_section = :default
end

exMC_ppbar = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_ppbar"
  config.events        = 200_000
  config.decay_card    = dc_ppbar
  config.cross_section = :default
end

exMC_pipieta = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_pipieta"
  config.events        = 200_000
  config.decay_card    = dc_pipieta
  config.cross_section = :default
end

exMC_pipipi0pi0 = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_pipipi0pi0"
  config.events        = 200_000
  config.decay_card    = dc_pipipi0pi0
  config.cross_section = :default
end

exMC_ppbarpi0 = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_ppbarpi0"
  config.events        = 200_000
  config.decay_card    = dc_ppbarpi0
  config.cross_section = :default
end

exMC_ppbarpipi = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "eta_c_eta_pipi_ppbarpipi"
  config.events        = 200_000
  config.decay_card    = dc_ppbarpipi
  config.cross_section = :default
end

### =========================== Event selection ============================= ###
# Three selection families are used, following the description:
#   :allcharged -> 4 positive + 4 negative tracks, exactly 2 photons
#   :neutral   -> 3 positive + 3 negative tracks, >= 4 photons (pi0 / eta present)
#   :ppbar     -> 2 positive + 2 negative tracks, exactly 2 photons

# ---------------------------------------------------------------------------
# (1) eta_c -> pi+ pi- K+ K-   (all-charged family, 2K+ 2K-)
# ---------------------------------------------------------------------------
alg_name = "EtaCPipiKK"
alg_pipiKK = Algorithm.new(alg_name)
alg_pipiKK.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:ecms_per_point, "the five data points have different CMS energies; ECMS must be set to the value of the running point (4.230/4.260/4.360/4.420/4.600 GeV)")

sel_pipiKK = Selection.new
sel_pipiKK.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==4"
            nChrn     "==4"
            nNet      "==0"
          }
          .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    10.0
            nGam              "==2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]   # K+ and K- together
            nkp "==2"
            nkm "==2"
          }
          .remove([:kp <= :chrgp, :km <= :chrgn])       # pull identified kaons out
          .assign({:chrgp => :pip, :chrgn => :pim})     # all remaining tracks are pions
          .kinematic_fit([:kp, :kp, :km, :km, :pip, :pip, :pim, :pim, :gamma, :gamma]) {
            nominal
            constrain_four_momentum
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 200        # loose; final per-mode chi2 is set in ROOT
          }
          .note(:chi2_cut_optimisation, "per-final-state chi2 cuts are chosen to retain 90% of the signal; BOSS keeps the loose chi2<200 and the tuned value is applied in ROOT")
          .note(:multi_candidate_retention, "all candidate combinations with identical chi2 are retained instead of a single one being picked")
alg_pipiKK.with_decay_card(dc_pipiKK).apply(sel_pipiKK)
alg_pipiKK.execute_on(data_scan + incMC_scan + exMC_pipiKK)

# ---------------------------------------------------------------------------
# (2) eta_c -> 2(K+ K-)   (all-charged family, 4K+ 4K-)
# ---------------------------------------------------------------------------
alg_name = "EtaC2KK"
alg_2KK = Algorithm.new(alg_name)
alg_2KK.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.260]})
sel_2KK = Selection.new
sel_2KK.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==4"
         nChrn     "==4"
         nNet      "==0"
       }
       .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              "==2"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         nkp "==4"
         nkm "==4"
       }
       .kinematic_fit([:kp, :kp, :kp, :kp, :km, :km, :km, :km, :gamma, :gamma]) {
         nominal
         constrain_four_momentum
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 200
       }
alg_2KK.with_decay_card(dc_2KK).apply(sel_2KK)
alg_2KK.execute_on(data_scan + incMC_scan + exMC_2KK)

# ---------------------------------------------------------------------------
# (3) eta_c -> 2(pi+ pi-)   (all-charged family)
# ---------------------------------------------------------------------------
alg_name = "EtaC2Pipi"
alg_2pipi = Algorithm.new(alg_name)
alg_2pipi.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
sel_2pipi = Selection.new
sel_2pipi.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==4"
           nChrn     "==4"
           nNet      "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track    10.0
           nGam              "==2"
         }
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kinematic_fit([:pip, :pip, :pip, :pip, :pim, :pim, :pim, :pim, :gamma, :gamma]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 200
         }
alg_2pipi.with_decay_card(dc_2pipi).apply(sel_2pipi)
alg_2pipi.execute_on(data_scan + incMC_scan + exMC_2pipi)

# ---------------------------------------------------------------------------
# (4) eta_c -> 3(pi+ pi-)   (all-charged family)
# ---------------------------------------------------------------------------
alg_name = "EtaC3Pipi"
alg_3pipi = Algorithm.new(alg_name)
alg_3pipi.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
sel_3pipi = Selection.new
sel_3pipi.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==4"
           nChrn     "==4"
           nNet      "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track    10.0
           nGam              "==2"
         }
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kinematic_fit([:pip, :pip, :pip, :pip, :pim, :pim, :pim, :pim, :gamma, :gamma]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 200
         }
alg_3pipi.with_decay_card(dc_3pipi).apply(sel_3pipi)
alg_3pipi.execute_on(data_scan + incMC_scan + exMC_3pipi)

# ---------------------------------------------------------------------------
# (5) eta_c -> K_S0 K+- pi-+   (all-charged family; K_S0 via secondary vertex)
# ---------------------------------------------------------------------------
alg_name = "EtaCKSKPi"
alg_KSKpi = Algorithm.new(alg_name)
alg_KSKpi.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .note(:ks0_selection, "K_S0 candidate must pass a flight-length significance L/sigma_L > 2 and lie within a +-15 MeV/c2 window around the K_S0 mass; both criteria are evaluated after the secondary vertex fit")
sel_KSKpi = Selection.new
sel_KSKpi.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==4"
           nChrn     "==4"
           nNet      "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track    10.0
           nGam              "==2"
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]
           nkp "==1"
           nkm "==1"
         }
         .remove([:kp <= :chrgp, :km <= :chrgn])
         .assign({:chrgp => :pip, :chrgn => :pim})
         .secondary_vertex_fit([:pip, :pim]) {
           build_virtual_particle(:K_S0).by_minimizing_mass_difference
           remove_used_particle_from_candidate_list
         }
         .kinematic_fit([:kp, :km, :K_S0, :pip, :pip, :pim, :pim, :gamma, :gamma]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 200
         }
alg_KSKpi.with_decay_card(dc_KSKpi).apply(sel_KSKpi)
alg_KSKpi.execute_on(data_scan + incMC_scan + exMC_KSKpi)

# ---------------------------------------------------------------------------
# (6) eta_c -> K+ K- pi0   (neutral family; Kalman gamma gamma -> pi0)
# ---------------------------------------------------------------------------
alg_name = "EtaCKKPi0"
alg_KKpi0 = Algorithm.new(alg_name)
alg_KKpi0.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .note(:pi0_mass_window, "pi0 candidates are required to lie in the mass window [110, 150] MeV/c2")
sel_KKpi0 = Selection.new
sel_KKpi0.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==3"
           nChrn     "==3"
           nNet      "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track    10.0
           nGam              ">=4"
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]
           nkp "==1"
           nkm "==1"
         }
         .remove([:kp <= :chrgp, :km <= :chrgn])
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 30
           npi0 ">=1"
         }
         .kinematic_fit([:kp, :km, :pip, :pip, :pim, :pim, :pi0]) {
           nominal
           constrain_four_momentum
           chi2_cut 200
         }
alg_KKpi0.with_decay_card(dc_KKpi0).apply(sel_KKpi0)
alg_KKpi0.execute_on(data_scan + incMC_scan + exMC_KKpi0)

# ---------------------------------------------------------------------------
# (7) eta_c -> K+ K- eta   (neutral family; Kalman gamma gamma -> eta)
# ---------------------------------------------------------------------------
alg_name = "EtaCKKEta"
alg_KKeta = Algorithm.new(alg_name)
alg_KKEta.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .note(:eta_mass_window, "eta candidates are required to lie in the mass window [500, 570] MeV/c2")
sel_KKeta = Selection.new
sel_KKeta.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==3"
           nChrn     "==3"
           nNet      "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track    10.0
           nGam              ">=4"
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion, :proton]
           nkp "==1"
           nkm "==1"
         }
         .remove([:kp <= :chrgp, :km <= :chrgn])
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 30
           neta ">=1"
         }
         .kinematic_fit([:kp, :km, :pip, :pip, :pim, :pim, :eta]) {
           nominal
           constrain_four_momentum
           chi2_cut 200
         }
alg_KKeta.with_decay_card(dc_KKeta).apply(sel_KKeta)
alg_KKeta.execute_on(data_scan + incMC_scan + exMC_KKeta)

# ---------------------------------------------------------------------------
# (8) eta_c -> p+ anti-p-   (ppbar family)
# ---------------------------------------------------------------------------
alg_name = "EtaCPpbar"
alg_ppbar = Algorithm.new(alg_name)
alg_ppbar.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .note(:multi_candidate_retention, "all candidate combinations with identical chi2 are retained instead of a single one being picked")
sel_ppbar = Selection.new
sel_ppbar.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==2"
           nChrn     "==2"
           nNet      "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track    10.0
           nGam              "==2"
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :proton, against: [:kaon, :pion]   # p and pbar together
           nprp "==1"
           nprm "==1"
         }
         .remove([:prp <= :chrgp, :prm <= :chrgn])
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 200
         }
alg_ppbar.with_decay_card(dc_ppbar).apply(sel_ppbar)
alg_ppbar.execute_on(data_scan + incMC_scan + exMC_ppbar)

# ---------------------------------------------------------------------------
# (9) eta_c -> pi+ pi- eta   (neutral family; Kalman gamma gamma -> eta)
# ---------------------------------------------------------------------------
alg_name = "EtaCPipiEta"
alg_pipieta = Algorithm.new(alg_name)
alg_pipieta.set_header(["#{alg_name}Alg/#{alg_name}.h"])
           .set_constant({"ECMS" => [:double, 4.260]})
sel_pipieta = Selection.new
sel_pipieta.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     "==3"
             nChrn     "==3"
             nNet      "==0"
           }
           .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=4"
           }
           .assign({:chrgp => :pip, :chrgn => :pim})
           .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 30
             neta ">=1"
           }
           .kinematic_fit([:pip, :pip, :pip, :pim, :pim, :pim, :eta]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }
alg_pipieta.with_decay_card(dc_pipieta).apply(sel_pipieta)
alg_pipieta.execute_on(data_scan + incMC_scan + exMC_pipieta)

# ---------------------------------------------------------------------------
# (10) eta_c -> pi+ pi- pi0 pi0   (neutral family; Kalman gamma gamma -> pi0)
# ---------------------------------------------------------------------------
alg_name = "EtaCPipiPi0Pi0"
alg_pipipi0pi0 = Algorithm.new(alg_name)
alg_pipipi0pi0.set_header(["#{alg_name}Alg/#{alg_name}.h"])
              .set_constant({"ECMS" => [:double, 4.260]})
              .note(:pi0_mass_window, "each pi0 candidate is required to lie in the mass window [110, 150] MeV/c2")
sel_pipipi0pi0 = Selection.new
sel_pipipi0pi0.select_track {
                 cos_theta 0.93
                 Vz        10.0
                 Vr        1.0
                 nChrp     "==3"
                 nChrn     "==3"
                 nNet      "==0"
               }
               .select_photon {
                 tdc_emc_start     0
                 tdc_emc_end       14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 angle_to_track    10.0
                 nGam              ">=4"
               }
               .assign({:chrgp => :pip, :chrgn => :pim})
               .kalman_kinematic_fit([:gamma, :gamma]) {
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 30
                 npi0 ">=1"
               }
               .kinematic_fit([:pip, :pip, :pip, :pim, :pim, :pim, :pi0]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }
alg_pipipi0pi0.with_decay_card(dc_pipipi0pi0).apply(sel_pipipi0pi0)
alg_pipipi0pi0.execute_on(data_scan + incMC_scan + exMC_pipipi0pi0)

# ---------------------------------------------------------------------------
# (11) eta_c -> p+ anti-p- pi0   (neutral family; Kalman gamma gamma -> pi0)
# ---------------------------------------------------------------------------
alg_name = "EtaCPpbarPi0"
alg_ppbarpi0 = Algorithm.new(alg_name)
alg_ppbarpi0.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})
            .note(:pi0_mass_window, "pi0 candidates are required to lie in the mass window [110, 150] MeV/c2")
sel_ppbarpi0 = Selection.new
sel_ppbarpi0.select_track {
               cos_theta 0.93
               Vz        10.0
               Vr        1.0
               nChrp     "==3"
               nChrn     "==3"
               nNet      "==0"
             }
             .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    10.0
               nGam              ">=4"
             }
             .pid(method: :probability) {
               prob_cut 0.001
               identify :proton, against: [:kaon, :pion]
               nprp "==1"
               nprm "==1"
             }
             .remove([:prp <= :chrgp, :prm <= :chrgn])
             .assign({:chrgp => :pip, :chrgn => :pim})
             .kalman_kinematic_fit([:gamma, :gamma]) {
               invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
               chi2_cut 30
               npi0 ">=1"
             }
             .kinematic_fit([:prp, :prm, :pip, :pip, :pim, :pim, :pi0]) {
               nominal
               constrain_four_momentum
               chi2_cut 200
             }
alg_ppbarpi0.with_decay_card(dc_ppbarpi0).apply(sel_ppbarpi0)
alg_ppbarpi0.execute_on(data_scan + incMC_scan + exMC_ppbarpi0)

# ---------------------------------------------------------------------------
# (12) eta_c -> p+ anti-p- pi+ pi-   (ppbar family)
# ---------------------------------------------------------------------------
alg_name = "EtaCPpbarPipi"
alg_ppbarpipi = Algorithm.new(alg_name)
alg_ppbarpipi.set_header(["#{alg_name}Alg/#{alg_name}.h"])
             .set_constant({"ECMS" => [:double, 4.260]})
sel_ppbarpipi = Selection.new
sel_ppbarpipi.select_track {
               cos_theta 0.93
               Vz        10.0
               Vr        1.0
               nChrp     "==2"
               nChrn     "==2"
               nNet      "==0"
             }
             .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    10.0
               nGam              "==2"
             }
             .pid(method: :probability) {
               prob_cut 0.001
               identify :proton, against: [:kaon, :pion]
               nprp "==1"
               nprm "==1"
             }
             .remove([:prp <= :chrgp, :prm <= :chrgn])
             .assign({:chrgp => :pip, :chrgn => :pim})
             .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :gamma]) {
               nominal
               constrain_four_momentum
               invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
               chi2_cut 200
             }
alg_ppbarpipi.with_decay_card(dc_ppbarpipi).apply(sel_ppbarpipi)
alg_ppbarpipi.execute_on(data_scan + incMC_scan + exMC_ppbarpipi)