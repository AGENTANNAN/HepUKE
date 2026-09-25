# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(3686) real data (BOSS 709)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # matching inclusive MC

# Decay card: ψ(3686) → γ χ_cJ, χ_cJ → Σ+ Σ̄− η → (p π0)(p̄ π0)(γγ), π0 → γγ
# (χ_c0, χ_c1, χ_c2 share the same final state, hence a single card / single algorithm)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 PHSP;
    0.3333 gamma chi_c1 PHSP;
    0.3333 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c0
    1.000 Sigma+ anti-Sigma- eta PHSP;
    Enddecay

    Decay chi_c1
    1.000 Sigma+ anti-Sigma- eta PHSP;
    Enddecay

    Decay chi_c2
    1.000 Sigma+ anti-Sigma- eta PHSP;
    Enddecay

    Decay Sigma+
    1.000 p+ pi0 HypWK;
    Enddecay

    Decay anti-Sigma-
    1.000 anti-p- pi0 HypWK;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 3M exclusive-MC events for the χ_cJ → Σ+Σ−η chain, using the same decay card
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chicJ_sigma_sigmabar_eta"
  config.related_dataset = psip_data
  config.events          = 3_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToGammaChicJSigmaSigmaEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})       # √s = 3.686 GeV

event_selection = Selection.new
event_selection
  .select_track {                       # charged track selection
      cos_theta 0.93                    # |cosθ| < 0.93
      Vz        15.0                    # |Vz| < 15 cm
      Vr        2.0                     # Vr < 2 cm
      nChrp     ">=1"                   # at least one positive track
      nChrn     ">=1"                   # at least one negative track
  }
  .select_photon {                      # photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0            # ≥ 10° away from any charged track
      energyThreshold_b 0.025           # 25 MeV (barrel)
      energyThreshold_e 0.050           # 50 MeV (endcap)
      nGam              ">=7"           # at least seven photons
  }
  .pid(method: :probability) {          # probability PID
      prob_cut 0.001                    # probability > 0.001
      identify :proton, against: [:kaon, :pion]   # p+ and p̄ vs K, π
      nprp ">=1"                        # at least one proton
      nprm ">=1"                        # at least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # set protons aside from the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # treat the remaining charged tracks as pions
  .kalman_kinematic_fit([:gamma, :gamma]) {   # π0 → γγ mass-constrained reconstruction
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0 ">=2"                        # at least two π0 candidates
  }
  # Nominal 6C fit: 4-momentum conservation + the two π0 and the η mass constraints
  .kinematic_fit([:gamma, :prp, :prm, :pi0, :pi0, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # η mass constraint
      invariant_mass_of(:gamma, :gamma).within(0.517, 0.577)   # η window
      invariant_mass_of(:prp, :pi0).within(1.174, 1.204)       # Σ+ → p π0 window
      invariant_mass_of(:prm, :pi0).within(1.174, 1.204)       # Σ− → p̄ π0 window
      chi2_cut 45
  }
  # Competing 4C hypothesis with an extra photon (background veto, χ² stored for ROOT comparison)
  .kinematic_fit([:gamma, :prp, :prm, :pi0, :pi0, :gamma, :gamma, :gamma]) {
      constrain_four_momentum
  }

my_algorithm
  .note(:radiative_photon_selection, "the radiative photon is chosen among the photons not used by the two π0 candidates as the one giving the largest ΔM")
  .note(:background_veto, "χ_cJ → γ J/ψ and χ_cJ → Σ+Σ−π0 backgrounds vetoed; the nominal-fit χ² is required to be smaller than the competing 4C (extra-photon) χ²")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])