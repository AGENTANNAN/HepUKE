# DSL for paper 2403.16811v2: e+e- → ηψ(2S) cross sections + search for ηX̃(3872)
# Ordinary analysis with full reconstruction at 18 energy points from 4.288-4.951 GeV
# Final state: γγπ+π-l+l- (l = e or μ)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Decay card for signal MC: e+e- → ηψ(2S) via KKMC
# ψ(2S) → π+π-J/ψ, J/ψ → l+l-, η → γγ
# Top mother psi(4260) per BESIII convention for KKMC generator
decay_card_etapsi2S = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta psi(2S) HELAMP;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for X̃(3872) search: e+e- → ηX̃(3872), X̃(3872) → π+π-J/ψ
decay_card_etaX3872 = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta anti-X_3872 PHSP;
    Enddecay

    Decay anti-X_3872
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Query data at 18 energy points: 4.288-4.951 GeV
# BOSS 705/706/707 XYZ samples
scan_data = DatasetManager.real_data.where(cms_energy: {value: 4288..4951})

# Inclusive MC for background estimation
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: {value: 4288..4951})

# Exclusive signal MC for ηψ(2S) - HELAMP model
scan_exMC_etapsi2S = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_etapsi2S"
  config.events = 200000
  config.decay_card = decay_card_etapsi2S
  config.cross_section = :default
end

# Exclusive signal MC for ηX̃(3872) - PHSP model
scan_exMC_etaX3872 = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_etaX3872"
  config.events = 200000
  config.decay_card = decay_card_etaX3872
  config.cross_section = :default
end

# Algorithm: single Algorithm for both processes (shared final state γγππll)
alg = Algorithm.new("EtaPsi2SAnalysis")
alg.set_header(["EtaPsi2SAnalysisAlg/EtaPsi2SAnalysis.h"])
   .set_constant({ "ECMS" => [:double, 4.682] })
   .note(:momentum_based_pid, "Charged tracks with momentum < 0.8 GeV/c assigned as π±; momentum > 1.0 GeV/c assigned as l±. Gap 0.8-1.0 GeV/c: tracks not used. Electron candidates require EMC energy > 1.0 GeV; muon candidates require EMC energy < 0.4 GeV")
   .note(:lepton_combination, "Both J/ψ → e+e- and J/ψ → μ+μ- channels are used. The best lepton-pair hypothesis (e+e- or μ+μ-) is selected by minimum χ²_4C in ROOT analysis")
   .note(:post_fit_mass_windows, "Post-4C-fit mass windows applied in ROOT: J/ψ: 3064.6 < M(l+l-) < 3140.8 MeV/c²; η: 507.1 < M(γγ) < 579.1 MeV/c²; ψ(2S): 3680.3 < M(π+π-J/ψ) < 3692.5 MeV/c². Background veto: M(π+π-γγ) > 1 GeV/c²; at √s > 4.600 GeV: M(γγJ/ψ) > 3.74 GeV/c²")
   .note(:fc_method, "Feldman-Cousins (F-C) method with likelihood ratio ordering used to obtain signal yields and 90% CL upper limits for low-statistics samples. POLE program for confidence intervals incorporating systematic uncertainties")
   .note(:background_subtraction, "γγψ(2S) background yield measured from sideband method with efficiency correction factor F; other backgrounds estimated from known cross sections. M(π+π-γγ) > 1 GeV/c² veto removes π0ψ(2S) and γ_ISR ψ(2S) contributions")
   .note(:isr_input_lineshape, "ISR correction factor (1+δ) uses Y(4260) lineshape input at √s < 4.600 GeV and 1/s power function at √s ≥ 4.600 GeV. Uncertainty from ψ(4415) and Y(4660) alternative lineshapes")
   .note(:track_momentum_condition, "Low momentum pions (< 0.8 GeV/c) and high momentum leptons (> 1.0 GeV/c) distinguished by laboratory-frame momentum measurement")
   .with_decay_card(decay_card_etapsi2S)

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz 100.0
                  Vr 10.0
                  nChrp "==2"
                  nChrn "==2"
                  nNet "==0"
                }
               .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end 14
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam ">=2"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6
                  identify :pion, against: [:kaon, :proton]
                  nlp "==1"
                  nlm "==1"
                  npip "==1"
                  npim "==1"
                }
               # Nominal 4C kinematic fit (J/ψ → e+e- hypothesis)
               .kinematic_fit([:gamma, :gamma, :pip, :pim, :ep, :em]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 40
                }
               # Competing 4C kinematic fit (J/ψ → μ+μ- hypothesis)
               # No chi2_cut, no nominal — stores χ² for ROOT-level channel selection
               .kinematic_fit([:gamma, :gamma, :pip, :pim, :mup, :mum]) {
                  constrain_four_momentum
                }

alg.with_decay_card(decay_card_etapsi2S).apply(event_selection)
alg.execute_on(scan_data + scan_incMC + scan_exMC_etapsi2S + scan_exMC_etaX3872)