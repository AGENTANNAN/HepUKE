# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # corresponding inclusive MC

# Decay card for the signal process psi(2S) -> pbar K+ Sigma0, Sigma0 -> gamma Lambda0, Lambda0 -> p+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 anti-p- K+ Sigma0     PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0         PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                HypWK;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for the pbar K+ Sigma0 mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pbarKpSigma0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PbarKpSigma0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])            # header file for the algorithm
            .set_constant({ "ECMS" => [:double, 3.686] })            # CMS energy = 3.686 GeV

# Full selection chain (shared by the mode and its charge conjugate)
event_selection = Selection.new
event_selection
  .select_track {                                     # charged track selection
      cos_theta     0.93                              # |cos(theta)| < 0.93
      Vz            10.0                              # |Vz| < 10 cm
      Vr            1.0                               # Vr < 1 cm
      nChrp         ">=2"                             # at least 2 positive tracks
      nChrn         ">=2"                             # at least 2 negative tracks
      nTot          "==4"                             # exactly 4 tracks in total
      nNet          "==0"                             # net charge zero
  }
  .select_photon {                                    # photon selection
      tdc_emc_start     0                             # TDC window 0-14 (x 700 ns)
      tdc_emc_end       14
      angle_to_track    10.0                          # >10 deg from nearest charged track
      energyThreshold_b 0.025                         # >25 MeV in barrel
      energyThreshold_e 0.050                         # >50 MeV in endcap
      nGam              ">=1"                         # at least one photon
  }
  .pid(method: :probability) {                        # particle ID (probability method)
      prob_cut 0.001                                  # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]       # p+ and pbar identified against K and pi
      identify :kaon,   against: [:proton, :pion]     # K+ and K- identified against p and pi
      nprp ">=1"                                      # at least one proton
      nprm ">=1"                                      # at least one anti-proton
      nkp  ">=1"                                      # at least one K+
  }
  # remaining charged tracks (not identified as p / pbar / K) are treated as pions
  .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) {               # build virtual Lambda from p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list        # remove used tracks from their lists
  }
  # Nominal 4C kinematic fit to e+e- -> pbar K+ Lambda gamma
  .kinematic_fit([:prm, :kp, :Lambda, :gamma]) {
      nominal                                         # use corrected 4-momenta from this fit
      constrain_four_momentum                         # 4C energy-momentum constraint
      invariant_mass_of(:Lambda).within(1.108, 1.123) # |M(p pi-) - m_Lambda| < 7.5 MeV/c^2
      chi2_cut 45                                     # chi^2 < 45
  }
  # Competing hypothesis: e+e- -> pbar K+ Lambda (no chi2_cut, no nominal -> stores chi2 for the ROOT-level veto)
  .kinematic_fit([:prm, :kp, :Lambda]) {
      constrain_four_momentum
  }
  # Competing hypothesis: e+e- -> pbar K+ Lambda gamma gamma (extra photon)
  .kinematic_fit([:prm, :kp, :Lambda, :gamma, :gamma]) {
      constrain_four_momentum
  }

# BOSS-side procedures that cannot be expressed in the DSL
my_algorithm
  .note(:helix_correction, "Helix-parameter corrections are applied to all charged tracks
    before the 4C kinematic fit; the efficiency difference between with/without helix
    correction is estimated by re-running the BOSS selection on signal MC.")
  .note(:background_veto, "The photon recoil mass is required to exceed the nominal chi_c2
    mass by at least 15 MeV/c^2 (M_recoil(gamma) > 3.571 GeV/c^2) to suppress the chi_c2
    background; this veto is applied outside the DSL fit chain.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])