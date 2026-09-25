# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC

# Decay card for psi(2S) -> phi pi+ pi- eta, phi -> K+ K-, eta -> gamma gamma (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 phi pi+ pi- eta        PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                  VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events for the decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_phi_pipi_eta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "PhiPiPiEta"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})           # CMS energy = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                                    # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        100.0       # |Vz| < 100.0 cm
                  Vr        10.0        # Vr < 10.0 mm
                  nChrp     "==2"       # exactly 2 positive tracks
                  nChrn     "==2"       # exactly 2 negative tracks
                  nNet      "==0"       # net charge zero
                }
               .select_photon {                                   # Photon selection
                  tdc_emc_start     0    # TDC start
                  tdc_emc_end       14   # TDC end
                  energyThreshold_b 0.025 # EMC barrel threshold 25 MeV
                  energyThreshold_e 0.050 # EMC endcap threshold 50 MeV
                  angle_to_track    10.0  # min angle to nearest charged track > 10 deg
                  nGam              ">=2" # at least two photons
                }
               .pid(method: :probability) {                        # PID (probability method)
                  prob_cut 0.001      # probability > 0.001
                  identify :kaon, against: [:pion, :proton]        # K+ and K- vs pi/p
                  identify :pion, against: [:kaon, :proton]        # pi+ and pi- vs K/p
                  nkp  "==1"          # exactly one K+
                  nkm  "==1"          # exactly one K-
                  npip "==1"          # exactly one pi+
                  npim "==1"          # exactly one pi-
                }
               .remove([:kp <= :chrgp, :pip <= :chrgp,            # remove identified K+/pi+
                        :km <= :chrgn, :pim <= :chrgn])           # and K-/pi- from generic charged lists
               .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {  # 4C kinematic fit
                  nominal                     # nominal fit; its four-momenta are the ones kept
                  constrain_four_momentum     # 4-momentum conservation to CMS energy
                  chi2_cut 80                 # chi^2 < 80
                  # phi mass window |M(K+K-) - m_phi| < 0.015 GeV  (m_phi = 1.019461 GeV)
                  invariant_mass_of(:kp, :km).within(1.0045, 1.0345)
                  # eta mass window |M(gamma gamma) - m_eta| < 0.040 GeV  (m_eta = 0.547862 GeV)
                  invariant_mass_of(:gamma, :gamma).within(0.5079, 0.5879)
                }

# J/psi background veto (conditional on M(pi+pi-eta) > 1.1 GeV) cannot be expressed
# with the available DSL cut primitives -> captured as a note.
my_Algorithm
  .note(:background_veto, "veto J/psi backgrounds: reject events with " \
        "|M(K+K- gamma gamma) - 3.0969| < 0.05 GeV or " \
        "|M(K+K- pi+ pi-) - 3.0969| < 0.05 GeV, applied only when M(pi+pi-eta) > 1.1 GeV")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])