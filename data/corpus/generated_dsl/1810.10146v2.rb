# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# psi(3686) real data and inclusive MC at 3.686 GeV
psip_data  = DatasetManager.real_data.find("709_3686")     # real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # inclusive MC

# Decay card for the signal process (EvtGen format).
# psi(3686) -> gamma chi_c0 (chi_c1/chi_c2 share the same final state),
# chi_c0 -> omega phi, omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0    PHSP;
    Enddecay

    Decay chi_c0
    1.0000 omega phi    PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0    OMEGA_DALITZ;
    Enddecay

    Decay phi
    1.0000 K+ K-    VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# 100k exclusive MC events for psi(3686) -> gamma chi_c0 -> gamma omega phi
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic0_omegaphi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaChiC0OmegaPhi"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
  .select_track {                       # Charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==2"                     # exactly two positive tracks
    nChrn     "==2"                     # exactly two negative tracks
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # Photon selection
    tdc_emc_start     0                 # TDC start time
    tdc_emc_end       14                # TDC end time
    energyThreshold_b 0.025             # barrel energy threshold (25 MeV)
    energyThreshold_e 0.050             # endcap energy threshold (50 MeV)
    nGam              ">=3"             # at least three photons
  }
  .pid(method: :probability) {          # Kaon identification
    prob_cut 0.001                      # PID probability > 0.001
    identify :kaon, against: [:pion, :proton]  # identify K+ and K- against pions/protons
    nkp "==1"
    nkm "==1"
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])      # remove identified kaons from charged-track lists
  .assign({:chrgp => :pip, :chrgn => :pim})    # surviving tracks are assigned as pi+ / pi-
  .kalman_kinematic_fit([:gamma, :gamma]) {    # form pi0 from two photons (mass-constrained, min gamma gamma mass difference)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=1"
  }
  # Nominal 4C kinematic fit (signal: gamma(rad) + pi0 + K+ K- pi+ pi-)
  .kinematic_fit([:gamma, :pi0, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 60
  }
  # Competing 4-gamma hypothesis: stores chi2_4c_4gam for the ROOT-level veto
  # (chi2_3gamma < chi2_4gamma), used when more than three photons are present.
  .kinematic_fit([:gamma, :gamma, :pi0, :kp, :km, :pip, :pim]) {
    constrain_four_momentum
  }

# Generate the algorithm for the signal process and execute on all datasets
my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])