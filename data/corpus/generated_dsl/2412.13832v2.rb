### Dataset description ###
# ψ(2S) at 3.686 GeV: real data and inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for the signal: ψ(2S) → γ χ_cJ, χ_cJ → p p̄ η π0, η → γγ, π0 → γγ
# (radiative ψ(2S) → γ χ_c1 modelled with P2GC1; both η/π0 decay via PHSP)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 p+ anti-p- eta pi0 PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal (200k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_gamma_chic1_ppbar_eta_pi0"
  config.related_dataset = psip_data          # associate with the ψ(2S) real data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "ChiCJPPbarEtaPi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # ψ(2S) CMS energy
            .set_alias({"std::vector<double>" => "Vdouble"})

# Chain the whole selection up to (and including) the nominal 4C kinematic fit.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93   # |cosθ| < 0.93
                  Vz        10.0   # |Vz| < 10 cm
                  Vr        1.0    # Vr < 1 cm
                  nChrp     "==1"  # exactly one positive track
                  nChrn     "==1"  # exactly one negative track
                  nNet      "==0"  # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # EMC TDC window 0–14
                  tdc_emc_end       14
                  angle_to_track    10.0   # angle to nearest charged track > 10°
                  energyThreshold_b 0.025  # barrel  E > 25 MeV
                  energyThreshold_e 0.050  # endcap  E > 50 MeV
                  nGam              ">=5"  # at least five photons (radiative + η + π0)
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ and p̄ vs K, π
                  nprp ">=1"
                  nprm ">=1"
                }
               # 1-C Kalman fit: reconstruct π0 from a photon pair (mass constraint)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"
                }
               # 1-C Kalman fit: reconstruct η from a photon pair (mass constraint)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 25
                  neta ">=1"
                }
               # Nominal 4C fit: constrain γ p p̄ π0 η to the CMS four-momentum
               # (loose chi2 cut; tight/6C-level cuts applied in ROOT)
               .kinematic_fit([:gamma, :prp, :prm, :pi0, :eta]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])