# Dataset description
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the signal process ψ(2S) → γ χ_cJ, χ_cJ → p p̄ η η, η → γγ (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-p- eta eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gam_chic1_ppbar_etaeta"
  config.related_dataset = psip_data                # associate with the 3.686 GeV real data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GamChiC1PpbarEtaEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy 3.686 GeV

event_selection = Selection.new
  .select_track {                          # charged tracks: exactly one p and one p̄
    cos_theta 0.93                         # |cosθ| < 0.93
    Vz        100.0                        # |Vz| < 100 cm
    Vr        10.0                         # Vr < 10 cm
    nChrp     "==1"                        # exactly one positive track
    nChrn     "==1"                        # exactly one negative track
    nNet      "==0"                        # net charge zero
  }
  .select_photon {                         # photons: ≥5 (4 from the two η + 1 radiative)
    tdc_emc_start     0                    # EMC time window start
    tdc_emc_end       14                   # EMC time window end (0–700 ns, 14 × 50 ns)
    angle_to_track    10.0                 # photon–track angle > 10°
    energyThreshold_b 0.025                # barrel energy threshold 25 MeV
    energyThreshold_e 0.050                # endcap energy threshold 50 MeV
    nGam              ">=5"                # at least 5 photons
  }
  .pid(method: :probability) {             # PID: protons against kaons and pions
    prob_cut 0.001                         # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]  # p+ and p̄ (charge-conjugation shorthand)
    nprp ">=1"                             # at least one proton
    nprm ">=1"                             # at least one anti-proton
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {  # reconstruct η from γγ (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25                            # χ² < 25
    neta     ">=2"                         # at least two η candidates
  }
  .kinematic_fit([:prp, :prm, :eta, :eta]) { # nominal 4C kinematic fit of p p̄ η η to the CMS
    nominal
    constrain_four_momentum
    chi2_cut 200                           # loose in BOSS; tight cut (paper: 35) applied in ROOT
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :prp, :prm]) {  # competing 5γ p p̄ hypothesis
    constrain_four_momentum                # no chi2_cut, not nominal → stores χ² for ROOT-level comparison
  }

# Attach the decay card and render the algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])