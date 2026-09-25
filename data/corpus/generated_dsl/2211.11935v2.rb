# ============================================================
# ψ(3686) → γ η_c(2S),  η_c(2S) → π⁺π⁻ η,  η → γγ
# Search for η_c(2S) → π⁺π⁻η; the M1 photon recoils against the η_c(2S)
# E_cms = 3.686 GeV
# ============================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data (3.686 GeV)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # ψ(3686) inclusive MC

# Signal decay card (EvtGen format): ψ(3686) → γ η_c(2S) → γ π⁺π⁻ η → γ π⁺π⁻ γγ
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gamma_etac2s_pipimeta"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Etac2sToPipimEta"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                # exactly two charged tracks, net charge 0
                  cos_theta 0.93              # |cosθ| < 0.93
                  Vz        10.0              # |Vz| < 10 cm
                  Vr        1.0               # Vr < 1 cm
                  nChrp    "==1"
                  nChrn    "==1"
                  nNet     "==0"
                }
               .select_photon {               # ≥3 photons
                  tdc_emc_start 0             # EMC time window [0, 700] ns
                  tdc_emc_end   14
                  angle_to_track 10.0         # angle to nearest charged track > 10°
                  energyThreshold_b 0.025     # E > 25 MeV (barrel)
                  energyThreshold_e 0.040     # E > 40 MeV (endcap)
                  nGam   ">=3"
                }
               .pid(method: :probability) {   # π± vs K/p, one π⁺ and one π⁻
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]
                  npip "==1"
                  npim "==1"
                }
               # Form η candidates from photon pairs: 1C mass-constrained fit to the η mass
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 200
                  neta ">=1"
                }
               # 5C signal fit η γ π⁺π⁻ (four-momentum conservation + η mass) — nominal
               .kinematic_fit([:eta, :gamma, :pip, :pim]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }
               # Competing 2γ hypothesis — χ² stored for the ROOT-level χ²_5C comparison
               .kinematic_fit([:eta, :pip, :pim]) {
                  constrain_four_momentum
                }
               # Competing 4γ hypothesis — χ² stored for the ROOT-level χ²_5C comparison
               .kinematic_fit([:eta, :gamma, :gamma, :pip, :pim]) {
                  constrain_four_momentum
                }

my_Algorithm
  .note(:modified_four_momentum_fit, "the final M(π+π−η) spectrum is taken from a modified 4C fit in which the M1 photon energy is excluded (floated) from the energy-momentum constraint; the M1 photon is the photon not used to build the η → γγ candidate")
  .note(:background_veto, "J/ψ veto M(γ_M1 π+π−) <= 3.00 GeV/c² and the competing-hypothesis veto χ²_5C(3γ) < χ²_5C(2γ), χ²_5C(3γ) < χ²_5C(4γ) are applied at the ROOT stage using the stored χ² values")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])