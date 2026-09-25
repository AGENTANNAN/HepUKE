# =============================================================================
# BOSS event-selection specification for the Born cross-section measurement of
#     e+e- -> p pbar pi0 ,  pi0 -> gamma gamma
# at 20 BESIII R-scan energy points, sqrt(s) = 2.1000 ... 3.0800 GeV.
#
#   * low-energy points  (2.1000 - 2.2324 GeV) : partial reconstruction
#         Mode I   : p , pbar only (pi0 undetected)  -> 4C fit , chi2 < 200
#         Mode II  : pbar + gamma gamma (p missing)  -> 1C fit , chi2 <  60
#   * high-energy points (2.3094 - 3.0800 GeV) : full reconstruction
#         p , pbar , gamma gamma                     -> 4C fit , chi2 < 100
#
# The signal line-shape fits (recoil-mass M(p pbar) / M(gamma gamma) peak) are
# ROOT-level steps and are therefore out of scope of this BOSS specification.
# =============================================================================


### -------------------------------- Datasets --------------------------------
# R-scan real data.  Sample-name convention "[BOSS_version]_[CMS_MeV]" (BOSS 7.1.3).
low_energy_points = %w[
  713_2100 713_2126 713_2150 713_2175 713_2200 713_2232
].map { |name| DatasetManager.real_data.find(name) }        # 2.1000 - 2.2324 GeV

high_energy_points = %w[
  713_2309 713_2386 713_2396 713_2500 713_2644 713_2646 713_2700
  713_2800 713_2900 713_2950 713_2981 713_3000 713_3020 713_3080
].map { |name| DatasetManager.real_data.find(name) }        # 2.3094 - 3.0800 GeV


### --------------------------- ConExc decay card ----------------------------
# Continuum generator (ConExc): ISR up to second order plus vacuum-polarisation
# corrections; used for the exclusive continuum sample at every energy point.
# The literal token `ConExc` selects the no-KKMC simulation template, and NO
# `Particle vpho` line is written because the DSL injects the vpho four-momentum
# at each energy point of the scan.
decay_card_pppi0 = <<~DECAYCARD
    Decay vpho
    1.0  ConExc  p+  anti-p-  pi0;
    Enddecay

    Decay pi0
    1.0  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD


### ------------------- Exclusive continuum MC (200k/point) ------------------
exMC_low = DatasetManager.create_exclusive_mc_for(low_energy_points) do |config|
  config.sample_name   = "pppi0_conexc"
  config.events        = 200_000
  config.decay_card    = decay_card_pppi0
  config.cross_section = :default
end

exMC_high = DatasetManager.create_exclusive_mc_for(high_energy_points) do |config|
  config.sample_name   = "pppi0_conexc"
  config.events        = 200_000
  config.decay_card    = decay_card_pppi0
  config.cross_section = :default
end


### ----------------------- Mode I  (low energy) -----------------------------
# p and pbar detected, pi0 completely undetected (no photons reconstructed).
alg_name_modeI = "RscanPpbarPi0ModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 2.100]})  # per-point value injected by the scan
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nChrp     "==1"     # one proton
    nChrn     "==1"     # one antiproton
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end   14
    nGam          "==0" # "no photons" in this mode
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p / pbar separated from pi and K
    nprp     "==1"
    nprm     "==1"
  }
  # Four-momentum-constrained (4C) fit; the undetected pi0 is recovered from the
  # constraint so that the recoil mass M(p pbar) peaks at m(pi0).
  .kinematic_fit([:prp, :prm, :pi0]) {
    nominal
    miss_track_of :pi0
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeI.note(:partial_reconstruction,
               "low-energy Mode I reconstructs only p and pbar; the pi0 is not "
               "reconstructed and is recovered from the missing (recoil) mass M(p pbar)")
         .note(:ppbar_opening_angle,
               "cut p-pbar opening angle < 173 deg; no implemented DSL helper "
               "(cos_theta_between is a stub), enforced in the generated BOSS code")
         .note(:recoil_cos_theta,
               "cut |cos(theta_recoil)| <= 0.96 for the recoil against the p pbar "
               "pair; enforced in the generated BOSS code")


### ----------------------- Mode II  (low energy) ----------------------------
# pbar and pi0(gamma gamma) detected, proton missing.
alg_name_modeII = "RscanPpbarPi0ModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 2.100]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==0"     # proton is not detected
    nChrn     "==1"     # one antiproton
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    30.0    # photon-track angle > 30 deg
    energyThreshold_b 0.025   # E > 25 MeV (barrel)
    energyThreshold_e 0.050   # E > 50 MeV (endcap)
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprm     "==1"
  }
  # One-constraint (1C) fit: pi0 mass constraint on the gamma gamma pair.
  .kinematic_fit([:prm, :gamma, :gamma]) {
    nominal
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 60
  }

alg_modeII.note(:partial_reconstruction,
                "low-energy Mode II reconstructs pbar and pi0 -> gamma gamma; the "
                "proton is left undetected and the signal is extracted from the "
                "M(gamma gamma) line shape")


### -------------------- Full reconstruction (high energy) -------------------
alg_name_full = "RscanPpbarPi0Full"
alg_full = Algorithm.new(alg_name_full)
alg_full.set_header(["#{alg_name_full}Alg/#{alg_name_full}.h"])
        .set_constant({"ECMS" => [:double, 3.080]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_full = Selection.new
sel_full.select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0    # photon-track angle > 10 deg
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp     "==1"
    nprm     "==1"
  }
  # Full 4C fit of p pbar gamma gamma; the fit keeps the combination with the
  # smallest chi2 automatically.
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }


### ------------------------------ Apply & run -------------------------------
alg_modeI.with_decay_card(decay_card_pppi0).apply(sel_modeI)
alg_modeII.with_decay_card(decay_card_pppi0).apply(sel_modeII)
alg_full.with_decay_card(decay_card_pppi0).apply(sel_full)

root_files_modeI  = alg_modeI.execute_on(low_energy_points  + exMC_low)
root_files_modeII = alg_modeII.execute_on(low_energy_points + exMC_low)
root_files_full   = alg_full.execute_on(high_energy_points  + exMC_high)