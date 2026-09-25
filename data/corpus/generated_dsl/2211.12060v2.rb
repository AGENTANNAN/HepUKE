# BOSS DSL spec — search for the hidden-charm tetraquark Z_cs'^- in
# e+e- -> K+ D_s*- D*0 at 4.661, 4.682 and 4.699 GeV, using two
# partial-reconstruction tags (a D_s- tag and a D*0 tag).

### Dataset preparation ###
# Real data at the three energy points (BOSS release 706)
data_4661 = DatasetManager.real_data.find("706_4660")   # sqrt(s) ~ 4.661 GeV
data_4682 = DatasetManager.real_data.find("706_4680")   # sqrt(s) ~ 4.682 GeV
data_4699 = DatasetManager.real_data.find("706_4700")   # sqrt(s) ~ 4.699 GeV
scan_data = [data_4661, data_4682, data_4699]

# Corresponding inclusive MC samples
incMC_4661 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4699 = DatasetManager.inclusive_mc.find("706_4700")
scan_incMC = [incMC_4661, incMC_4682, incMC_4699]

# Decay card for the non-resonant signal e+e- -> K+ D_s*- D*0
# (no intermediate meson -> psi(4260) used as top mother, KKMC convention)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  K+  D_s*-  D*0    PHSP;
    Enddecay

    Decay D_s*-
    1.0000  D_s-  gamma       PHSP;
    Enddecay

    Decay D_s-
    1.0000  K+  K-  pi-      PHSP;
    Enddecay

    Decay D*0
    1.0000  D0  pi0           PHSP;
    Enddecay

    Decay D0
    1.0000  K-  pi+           PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma      PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the non-resonant signal at each scan point (default cross section)
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "zcs_signal_nonres"   # auto-suffixed per energy point
  config.events        = 500000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMCs_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

all_datasets = scan_data + scan_incMC + exMCs_signal

### Algorithm 1 — D_s- tag: bachelor K+ + D_s- -> K+K-pi-, D*0 left as recoil ###
# Decay-card recID mapping:
#   0 psi(4260)   1 K+ (bachelor)   2 D_s*-   3 D*0   4 D_s-   5 gamma(D_s*-)
#   6 K+(D_s-)    7 K-(D_s-)        8 pi-(D_s-)        9 D0   10 pi0   11 K-(D0)  12 pi+(D0)
alg_name_ds = "DsTag"
alg_ds = Algorithm.new(alg_name_ds)
alg_ds.set_header(["#{alg_name_ds}Alg/#{alg_name_ds}.h"])
      .set_constant({"ECMS" => [:double, 4.682]})

sel_ds = Selection.new
sel_ds.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"
      nChrn     ">=2"
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       10
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    10.0
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon]          # identify pi+/pi- against kaons
    }
    .remove([:pip <= :chrgp, :pim <= :chrgn])   # remove identified pions from charged lists
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]          # identify K+/K- against pions
      nkp ">=1"                                 # at least one K+
      nkm ">=1"                                 # at least one K-
    }
    # Partial reconstruction: bachelor K+ (recID 1) + D_s- (recID 4); D*0 is the untagged side
    .partial_rec([1, 4]) {
      best_combination_by_mass :"D_s-", 1.9683  # pick D_s- candidate closest to its nominal mass
      require_recoil_mass 2.11, 5.0             # RM(K+ D_s-) window
    }

alg_ds.note(:ds_mass_window, "D_s- mass window 1.955-1.980 GeV applied on the K+K-pi- system; best candidate chosen at the nominal D_s- mass 1.9683 GeV")
      .note(:ds_background_veto, "D_s- daughters required to satisfy M(K+K-) < 1.05 GeV OR 0.850 < M(K+pi-) < 0.930 GeV")
      .note(:recoil_mass_pdg_substitution, "recoil mass RM(K+ D_s-) in 2.11-5.0 GeV/c^2 computed after substituting the PDG D_s- mass")

alg_ds.with_decay_card(decay_card_signal).apply(sel_ds)

### Algorithm 2 — D*0 tag: bachelor K+ + D*0 -> D0 pi0, D0 -> K-pi+, D_s*- left as recoil ###
alg_name_dstar = "Dstar0Tag"
alg_dstar = Algorithm.new(alg_name_dstar)
alg_dstar.set_header(["#{alg_name_dstar}Alg/#{alg_name_dstar}.h"])
        .set_constant({"ECMS" => [:double, 4.682]})

sel_dstar = Selection.new
sel_dstar.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"
      nChrn     ">=2"
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       10
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      angle_to_track    10.0
      nGam ">=2"                                 # at least two photons (only in the D*0 tag)
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon]
    }
    .remove([:pip <= :chrgp, :pim <= :chrgn])
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]
      nkp ">=1"
      nkm ">=1"
    }
    # pi0 -> gamma gamma (1-C Kalman fit constraining the photon pair to the pi0 nominal mass)
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"                                 # at least one pi0 candidate
    }
    # Partial reconstruction: bachelor K+ (recID 1) + D*0 (recID 3); D_s*- is the untagged side
    .partial_rec([1, 3]) {
      best_combination_by_mass :D0,    1.8648    # pick D0 candidate closest to its nominal mass
      best_combination_by_mass :"D*0", 2.00685   # pick D*0 candidate closest to its nominal mass
      require_recoil_mass 2.102, 2.122           # RM(K+ D*0) window
    }

alg_dstar.note(:d0_dstar0_mass_window, "D0 mass window 1.850-1.880 GeV (nominal 1.8648) and D*0 mass window 2.000-2.014 GeV (nominal 2.00685) applied on the K-pi+ pi0 combinations")
         .note(:recoil_mass_pdg_substitution, "recoil mass RM(K+ D*0) in 2.102-2.122 GeV/c^2 computed after the PDG D*0-mass correction")

alg_dstar.with_decay_card(decay_card_signal).apply(sel_dstar)

### Execute both algorithms on all scan data, inclusive MC and exclusive MC ###
root_files_ds    = alg_ds.execute_on(all_datasets)
root_files_dstar = alg_dstar.execute_on(all_datasets)