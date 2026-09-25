### Dataset preparation ###
data_4260 = DatasetManager.real_data.find("703_4260")        # ψ(4260) real data at 4.258 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")    # combined 4.15–4.30 GeV inclusive MC

# --- Normalization channel decay cards: e+e- -> γ X(3872), X(3872) -> π+π- J/ψ ---
decay_card_norm_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_norm_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- Search channel decay cards: e+e- -> γ X(3872), X(3872) -> π0 χ_c1, χ_c1 -> γ J/ψ ---
decay_card_search_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 pi0 chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_search_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 pi0 chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 500k-event exclusive MC for each of the four modes ---
exMC_norm_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "x3872_pipi_jpsi_ee"
  config.related_dataset = data_4260
  config.events = 500_000
  config.decay_card = decay_card_norm_ee
  config.cross_section = :default
end

exMC_norm_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "x3872_pipi_jpsi_mumu"
  config.related_dataset = data_4260
  config.events = 500_000
  config.decay_card = decay_card_norm_mumu
  config.cross_section = :default
end

exMC_search_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "x3872_pi0chic1_ee"
  config.related_dataset = data_4260
  config.events = 500_000
  config.decay_card = decay_card_search_ee
  config.cross_section = :default
end

exMC_search_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "x3872_pi0chic1_mumu"
  config.related_dataset = data_4260
  config.events = 500_000
  config.decay_card = decay_card_search_mumu
  config.cross_section = :default
end

### Normalization channel: e+e- -> γ X(3872), X(3872) -> π+π- J/ψ, J/ψ -> l+l- ###
alg_norm_name = "X3872ToPipiJpsi"
alg_norm = Algorithm.new(alg_norm_name)
alg_norm.set_header(["#{alg_norm_name}Alg/#{alg_norm_name}.h"])
        .set_constant({"ECMS" => [:double, 4.258]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        # Bhabha / η / η' vetoes that act on angular variables cannot be expressed in
        # the BOSS DSL (cos_theta_of / cos_theta_between are not implemented).
        .note(:background_veto, "Bhabha vetoes applied in addition to the 4C fit: cos(theta) of the pi+pi- system < 0.98 and cos(theta) between the radiative photon and the nearest charged track < 0.98")
        .note(:efficiency_curve, "the optimized chi2/ndf < 10 selection is applied on the ROOT side after the nominal 4C fit; only the loose chi2_cut 200 is applied in BOSS")

sel_norm = Selection.new
sel_norm.select_track {                       # Charged track selection
    cos_theta 0.93                            # |cos(theta)| < 0.93
    Vz 10.0                                   # |Vz| < 10 cm
    Vr 1.0                                    # Vr < 1 cm
    nChrp ">=2"                               # at least two positive tracks
    nChrn ">=2"                               # at least two negative tracks
    nNet "==0"                                # net charge zero
  }
  .select_photon {                            # Photon selection
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"                                # at least one photon
  }
  .pid(method: :probability) {                # Probability-based PID
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> e, else mu
    identify :pion, against: [:kaon, :proton] # hadron PID for the π+π- pair entering the fit
    nlp "==1"                                 # exactly one l+
    nlm "==1"                                 # exactly one l-
  }
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {   # 4C fit
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).within(3.0769, 3.1169)          # J/ψ mass window ±20 MeV/c²
    invariant_mass_of(:gamma, :pip, :pim).larger_than(0.6)      # η veto: M(γπ+π-) > 0.6 GeV/c²
    invariant_mass_of(:gamma, :pip, :pim).out_of(0.9378, 0.9778) # η' veto: |M(γπ+π-)-m_η'| > 0.02 GeV/c²
    chi2_cut 200                                                 # loose cut; tight χ²/ndf<10 in ROOT
  }

alg_norm.with_decay_card(decay_card_norm_ee).apply(sel_norm)
root_files_norm = alg_norm.execute_on([data_4260, incMC_4260, exMC_norm_ee, exMC_norm_mumu])

### Search channel: e+e- -> γ X(3872), X(3872) -> π0 χ_c1, χ_c1 -> γ J/ψ, J/ψ -> l+l- ###
alg_search_name = "X3872ToPi0Chic1"
alg_search = Algorithm.new(alg_search_name)
alg_search.set_header(["#{alg_search_name}Alg/#{alg_search_name}.h"])
          .set_constant({"ECMS" => [:double, 4.258]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          # Vetoes based on photon-pair scans, recoil masses and mass differences are not
          # expressible with the current BOSS DSL primitives.
          .note(:background_veto, "pi0pi0 veto: every photon pair other than the chosen pi0 candidate must be more than 20 MeV/c^2 away from m_pi0")
          .note(:background_veto, "omega(782) veto: M(gamma pi0) < 0.732 GeV/c^2")
          .note(:background_veto, "ISR-gamma veto: recoil mass > 3.7 GeV/c^2")
          .note(:background_veto, "chi_cJ windows: Delta M0 < 25 MeV/c^2, Delta M1,2 < 20 MeV/c^2, loose chi_cJ mass window in [3.35, 3.60] GeV/c^2")
          .note(:efficiency_curve, "the optimized chi2/ndf < 5 selection is applied on the ROOT side after the nominal 5C fit; only the loose chi2_cut 200 is applied in BOSS")
          .note(:efficiency_curve, "J/psi-substituted masses M(pi+pi-J/psi) = M(pi+pi-l+l-) - M(l+l-) + m_J/psi and M(pi0 chi_cJ) = M(pi0 l+l-) - M(l+l-) + m_J/psi are formed in ROOT to improve the mass resolution")

sel_search = Selection.new
sel_search.select_track {                     # Charged track selection
    cos_theta 0.93                            # |cos(theta)| < 0.93
    Vz 10.0                                   # |Vz| < 10 cm
    Vr 1.0                                    # Vr < 1 cm
    nChrp ">=1"                               # at least one positive track
    nChrn ">=1"                               # at least one negative track
    nNet "==0"                                # net charge zero
  }
  .select_photon {                            # Photon selection
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=4"                                # at least four photons
  }
  .pid(method: :probability) {                # Same lepton identification as the normalization channel
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> e, else mu
    nlp "==1"                                 # exactly one l+
    nlm "==1"                                 # exactly one l-
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :lp, :lm]) {  # 5C fit: 4C + π0 mass constraint
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # π0 mass constraint on two photons
    invariant_mass_of(:lp, :lm).within(3.0769, 3.1169)                    # J/ψ mass window ±20 MeV/c²
    chi2_cut 200                                                          # loose cut; tight χ²/ndf<5 in ROOT
  }

alg_search.with_decay_card(decay_card_search_ee).apply(sel_search)
root_files_search = alg_search.execute_on([data_4260, incMC_4260, exMC_search_ee, exMC_search_mumu])