# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# 2.125 GeV point of the 22-point R-scan (BOSS 713, sample Rscan_2125)
data_2125  = DatasetManager.real_data.find("713_2125")
incMC_2125 = DatasetManager.inclusive_mc.find("713_2125")

# Decay card: continuum (ISR) production of omega eta' at this energy point.
# ConExc models the ISR radiation and the measured Born cross section; the DSL
# auto-detects the literal `ConExc` token, switches to the no-KKMC template and
# injects `Particle vpho <ECMS> 0.0` for the energy point, so no `Particle vpho`
# line is written here.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0000  omega  eta'   ConExc  -1;      # user-defined (vhdr) continuum final state: omega eta'
    Enddecay

    Decay omega
    1.0000  pi+  pi-  pi0        OMEGA_DALITZ;
    Enddecay

    Decay eta'
    1.0000  gamma  pi+  pi-      PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma         PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC for e+e- -> omega eta' at 2.125 GeV (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_2125_omega_etap"
  config.related_dataset = data_2125
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "OmegaEtaP"
omega_etap_alg = Algorithm.new(alg_name)
omega_etap_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
              .set_constant({ "ECMS" => [:double, 2.125] })   # 2.125 GeV R-scan point
              .set_alias({ "std::vector<double>" => "Vdouble" })
              .note(:isr_vacuum_polarization_correction,
                    "Born cross section of e+e- -> omega eta' is extracted with iterative ISR " \
                    "and vacuum-polarisation corrections; the (1+delta_ISR) and (1+delta_VP) " \
                    "factors for this energy point are taken from the ConExc generator log.")

event_selection = Selection.new
event_selection
  .select_track {                      # charged-track selection
      cos_theta 0.93                   # |cos(theta)| < 0.93
      Vz        10.0                   # |Vz| < 10 cm
      Vr        1.0                    # Vr < 1 cm
      nChrp     "==2"                  # exactly two positive tracks
      nChrn     "==2"                  # exactly two negative tracks
      nNet      "==0"                  # net charge zero
  }
  .select_photon {                     # photon selection
      tdc_emc_start     0              # TDC window 0-14
      tdc_emc_end       14
      energyThreshold_b 0.025          # > 25 MeV in the barrel
      energyThreshold_e 0.050          # > 50 MeV in the endcap
      angle_to_track    10.0           # more than 10 degrees from any charged track
      nGam              ">=3"          # at least three photons
  }
  .pid(method: :probability) {         # probability PID
      prob_cut 0.001                   # PID probability > 0.001
      identify :pion, against: [:kaon, :proton]  # pi+ and pi- identified against K and p
      npip "==2"                       # two pi+
      npim "==2"                       # two pi-
  }
  # Reconstruct pi0 from the photon pair with the smallest |M(gamma gamma) - m_pi0|:
  # the mass-constrained Kalman fit minimises exactly this mass difference.
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25                      # chi2 < 25
      npi0     ">=1"                   # at least one pi0 candidate
  }
  # Nominal 4C fit on pi+ pi- pi+ pi- pi0 gamma (pi0 enters as a single participant).
  # The omega / eta' mass windows pre-select the candidate combinations; among the
  # surviving combinations the fit keeps the one with the smallest chi2.
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:pip, :pim, :pi0).within(0.7537, 0.8117)    # |M(pi+pi-pi0) - m_omega| < 0.029
      invariant_mass_of(:gamma, :pip, :pim).within(0.9328, 0.9828)  # |M(gamma pi+pi-) - m_eta'| < 0.025
      chi2_cut 60                      # loose cut; the tight/ordered cut is applied in ROOT
  }
  # Competing 4C hypothesis 1: pi+ pi- pi+ pi- gamma with two photons missing
  # (e+e- -> 2(pi+pi-)pi0 background). No chi2_cut / no nominal: only the chi2 is stored,
  # the ordering veto chi2_4c_nominal < chi2_4c_4pi_gamma is applied in the ROOT analysis.
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma]) {
      constrain_four_momentum
      miss_track_of :gamma
      miss_track_of :gamma
  }
  # Competing 4C hypothesis 2: pi+ pi- pi+ pi- gamma gamma with one photon missing
  # (e+e- -> 2(pi+pi-pi0) background). chi2 stored only; ordering veto applied in ROOT.
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) {
      constrain_four_momentum
      miss_track_of :gamma
  }

# Generate the algorithm for the process defined in the decay card
omega_etap_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the real data, the inclusive MC and the exclusive signal MC
root_files = omega_etap_alg.execute_on([data_2125, incMC_2125, exMC_signal])