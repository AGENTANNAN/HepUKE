# BESIII Analysis: J/ψ → e+e- η(1405) → e+e- π0 f0(980) → e+e- π0 π+ π-
# Paper: arXiv:2307.14633v1 — First observation of J/ψ → e+e- η(1405)
# Dataset: 708_3097 (J/ψ at 3.097 GeV, (10087±44)×10^6 events)
# Analysis type: Ordinary (Algorithm + Selection)

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal decay card: J/ψ → e+ e- η(1405) → e+ e- π0 f0(980) → e+ e- π0 π+ π-
decay_card = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- eta(1405)   PHSP;
    Enddecay

    Decay eta(1405)
    1.0000 pi0 f0(980)       PHSP;
    Enddecay

    Decay f0(980)
    1.0000 pi+ pi-           PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma       PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_ee_eta1405"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

### Event selection ###
alg = Algorithm.new("JpsiEEeta1405")
alg.set_header(["JpsiEEeta1405Alg/JpsiEEeta1405.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })
   .note(:e_p_cut, "E/p > 0.8 requirement for e± with p > 0.8 GeV/c is applied in ROOT analysis; not expressible in BOSS selection")
   .note(:gamma_conversion_veto, "Rxy < 2 cm veto on e+e- pair from gamma conversions is applied in ROOT analysis")
   .note(:f0_mass_window, "f0(980) mass window on M(π+π-) is applied in ROOT analysis after kinematic fit")
   .note(:pi0_mass_window, "π0 mass window (±3σ Crystal Ball) on M(γγ) is applied in ROOT analysis")

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93
                  Vz   10.0
                  Vr   1.0
                  nChrp  "==2"
                  nChrn  "==2"
                  nTot    "==4"
                  nNet    "==0"
                }
               .select_photon {
                  tdc_emc_start   0
                  tdc_emc_end   700
                  angle_to_track   10.0
                  energyThreshold_b 0.025
                  energyThreshold_e  0.050
                  nGam   ">=2"
                }
               .pid(method: :probability) {
                  prob_cut   0.001
                  identify :electron, against: [:kaon, :pion]
                  nep   "==2"
                  nem   "==2"
               }
               .remove([:em <= :chrgn])
               .remove([:ep <= :chrgp])
               .assign({:chrgp => :pip, :chrgn => :pim})
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200
                  npi0  ">=1"
               }
               .kinematic_fit([:ep, :em, :pip, :pim, :pi0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 20
               }

alg.with_decay_card(decay_card).apply(event_selection)
root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])