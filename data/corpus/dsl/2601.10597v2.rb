# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> phi eta, eta -> pi0 + invisible (S -> chi chi_bar)
# S is treated as invisible; simulate as J/psi -> phi eta with eta -> pi0 nu nu (2-body-like phase space)
decay_card_for_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta                         PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                           VSS;
    Enddecay

    Decay eta
    1.0000 pi0 nu_e anti-nu_e              PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                     PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_phi_eta_pi0_invisible"
  config.related_dataset = jpsi_data
  config.events         = 500000
  config.decay_card     = decay_card_for_signal
  config.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "EtaPi0Invisible"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                    # Charged track selection
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     ">=1"                 # >=1 K+
                  nChrn     ">=1"                 # >=1 K-
                  nTot      "==2"                 # exactly 2 charged tracks required for the signal selection
                }
               .select_photon {                   # Photon selection (for the pi0 -> gamma gamma from eta)
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=2"         # >= 2 good photons
                }
               .pid(method: :probability) {       # Kaon PID: L(K)>0 and L(K)>L(e,mu,pi)
                  prob_cut 0.0
                  identify :kaon, against: [:pion, :proton]
                  nkp      ">=1"
                  nkm      ">=1"
               }
               # Nominal kinematic fit: eta from K+ K- recoil, mass-constrained;
               # further constrain gamma gamma to pi0 (from eta -> pi0 + invisible S)
               .kinematic_fit([:kp, :km, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:kp, :km).within(1.00, 1.04)
                  invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200
               }

my_Algorithm.note(:eta_recoil_window,
                  "IDT selection requires the recoiling mass of K+K- to lie in [0.45, 0.65] GeV/c^2 to tag eta")
            .note(:eta_mass_constraint,
                  "A prior kinematic fit constrains the recoiling mass of K+K- to the nominal eta mass; the best K+K- combination is chosen by minimum chi^2")
            .note(:additional_photon_energy,
                  "Suppress additional photons/pi0's by requiring the total energy of photons other than those forming the pi0 to be less than 0.1 GeV (E_oth_gamma < 0.1 GeV)")
            .note(:recoil_angle_cut,
                  "|cos theta_recoil(K+K-pi0)| < 0.7 required to reduce background particles flying into the end-cap regions")
            .note(:mkkpi0_cut,
                  "Type I: M(K+K-pi0) < 2.7 GeV/c^2 for m_S in [0,165) MeV/c^2; Type II: < 2.15 GeV/c^2 for m_S in [165,305) MeV/c^2; Type III: < 2.0 GeV/c^2 for m_S in [305,400] MeV/c^2 (Punzi-optimised)")

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
