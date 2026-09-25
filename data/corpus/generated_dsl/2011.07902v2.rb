### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/ψ (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding J/ψ inclusive MC

# Decay card — signal mode: J/ψ → γη', η' → π+π−e+e−
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- e+ e-  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — normalization mode: J/ψ → γη', η' → π+π−γ
decay_card_norm = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (500k events each)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_pipim_ee"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_norm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_pipim_gamma"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_norm
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ---------- Signal mode: J/ψ → γη', η' → π+π−e+e− ----------
alg_name_signal = "EtaPrimeToPiPiEE"
alg_signal = Algorithm.new(alg_name_signal)
alg_signal.set_header(["#{alg_name_signal}Alg/#{alg_name_signal}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})   # double ECMS = 3.097 GeV

sel_signal = Selection.new
  .select_track {                 # charged track selection
      cos_theta 0.93              # |cosθ| < 0.93
      Vz        10.0              # |Vz| < 10 cm
      Vr        1.0               # Vr < 1 cm
      nChrp     "==2"             # two positive tracks
      nChrn     "==2"             # two negative tracks
      nNet      "==0"             # net charge zero
  }
  .select_photon {                # photon selection
      tdc_emc_start     0         # EMC timing window start
      tdc_emc_end       14        # EMC timing window end
      angle_to_track    15.0      # > 15° to any charged track
      energyThreshold_b 0.025     # > 25 MeV in the barrel
      energyThreshold_e 0.050     # > 50 MeV in the endcap
      nGam              ">=1"     # at least one photon
  }
  .pid(method: :chi2_sum) {       # chi2-sum combinatorial PID
      chi_min_cut 4               # combined PID chi² < 4
      identify :pion, :electron   # assign each track as π or e (bijection per charge)
  }
  # Nominal 4C kinematic fit: γπ+π−e+e− hypothesis
  .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {
      nominal
      constrain_four_momentum
      chi2_cut 62
  }
  # Competing hypothesis γπ+π−μ+μ−: same tracks reinterpreted as muons; store χ² only
  .assign({:ep => :mup, :em => :mum})
  .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) {
      use_track_index_from_nominal_kmfit
      constrain_four_momentum
  }
  # Competing hypothesis γπ+π−π+π−: same tracks reinterpreted as pions; store χ² only
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
      use_track_index_from_nominal_kmfit
      constrain_four_momentum
  }

alg_signal
  .note(:background_veto, "photon-conversion veto (γ → e+e−) applied to remove converted-photon background; handled outside the formal DSL")
  .note(:chi2_hypothesis_selection, "event retained only if the γπ+π−e+e− hypothesis has the smallest combined 4C+PID χ² among the γπ+π−e+e−, γπ+π−μ+μ− and γπ+π−π+π− hypotheses; comparison performed in ROOT from the stored χ² values")
  .with_decay_card(decay_card_signal).apply(sel_signal)

# Keep the η' signal window |M(π+π−e+e−) − m_η'| < 0.02 GeV for the ROOT stage
# (applied on kinematic-fit-corrected variables → out of BOSS scope).

# ---------- Normalization mode: J/ψ → γη', η' → π+π−γ ----------
alg_name_norm = "EtaPrimeToPiPiGamma"
alg_norm = Algorithm.new(alg_name_norm)
alg_norm.set_header(["#{alg_name_norm}Alg/#{alg_name_norm}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_norm = Selection.new
  .select_track {                 # charged track selection (same track cuts)
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"             # one positive track
      nChrn     "==1"             # one negative track
      nNet      "==0"
  }
  .select_photon {                # photon selection (same photon cuts)
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    15.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"     # at least two photons
  }
  .pid(method: :chi2_sum) {       # same chi2-sum combinatorial PID
      chi_min_cut 4
      identify :pion              # both tracks assigned as pions
  }
  # 4C kinematic fit: γπ+π−γγ hypothesis
  .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 140
  }

alg_norm
  .note(:background_veto, "radiative photon selected within 20 MeV of 1.4 GeV (|E_γ − 1.4| < 0.02 GeV); the other photon required to have E > 0.15 GeV to suppress π⁰/η backgrounds; per-photon energy windows not expressible in the DSL")
  .with_decay_card(decay_card_norm).apply(sel_norm)

### Execution ###
root_files_signal = alg_signal.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
root_files_norm   = alg_norm.execute_on([jpsi_data, jpsi_incMC, exMC_norm])