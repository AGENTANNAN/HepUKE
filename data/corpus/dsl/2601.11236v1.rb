# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) data (2712 M events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Signal process: psi(3686) -> pi0 h_c, h_c -> gamma eta_c, eta_c -> gamma gamma
decay_card_for_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c                         PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c                     PHSP;
    Enddecay

    Decay eta_c
    1.0000 gamma gamma                     PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                     PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_pi0hc_gamma_etac_gg"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_for_signal
  config.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "PsipPi0HcEtacGG"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_photon {                    # Photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14              # 700 ns after event start
                  angle_to_track    10.0
                  energyThreshold_b 0.025           # E > 25 MeV in barrel
                  energyThreshold_e 0.050           # E > 50 MeV in endcap
                  nGam              ">=4"           # tag: pi0(2 gamma) + E1 gamma; signal: 2 gamma from eta_c
                }
                # Reconstruct pi0 from two photons (both in barrel, E>40 MeV, 0.120<M(gg)<0.145)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 5                        # chi2_1C < 5 (Punzi optimised)
                  npi0    ">=1"
               }
                # Nominal kinematic fit: signal side psi(2S) -> pi0 gamma_E1 gamma gamma (eta_c -> gg)
                # 4C constrains total four-momentum of pi0 + 3 photons to initial e+e- system
               .kinematic_fit([:pi0, :gamma, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 70                        # chi2_4C < 70
               }

my_Algorithm.note(:tag_side,
                  "Tag side: pi0 selection uses two showers in the barrel with E>40 MeV and 0.120<M(gg)<0.145 GeV/c^2; kinematic fit constrains M(gg) to nominal pi0 mass with chi2_1C<5; require RM(pi0) in the h_c mass window [3.48, 3.57] GeV/c^2; if multiple pi0, keep the smallest chi2_1C")
            .note(:e1_photon_selection,
                  "E1 photon selection: exclude any bachelor photon that can pair with any other photon to form a pi0 candidate (invariant mass in [0.115, 0.150] GeV/c^2 and chi2_1C<200). Require the E1 photon energy in the pi0-recoiling system within [0.46, 0.58] GeV. If multiple pi0-gamma candidates remain, randomly pick one.")
            .note(:pipi_Jpsi_veto,
                  "Reject events where any pi+pi- recoil mass falls in [m(Jpsi)-7, m(Jpsi)+7] MeV/c^2 (from psi(2S)->pi+pi-Jpsi)")
            .note(:pi0pi0_Jpsi_veto,
                  "Reject events where any pi0 pi0 recoil mass falls in [m(Jpsi)-15, m(Jpsi)+25] MeV/c^2 (from psi(2S)->pi0pi0Jpsi)")
            .note(:eta3pi_veto,
                  "Reject events where any pi+pi-pi0 invariant mass lies within [m(eta)-12, m(eta)+12] MeV/c^2 to suppress eta->pi+pi-pi0 misidentification")
            .note(:eta_etap_veto,
                  "On signal side, reject events where the invariant mass of the E1 photon paired with either eta_c gamma lies in [0.50, 0.56] GeV/c^2 or [0.88, 0.98] GeV/c^2 (h_c->gamma eta / gamma eta' peaking backgrounds)")
            .note(:etac_mass_window,
                  "Require M(eta_c) = M(gamma gamma from eta_c) in [2889.1, 3066.8] MeV/c^2 for the h_c fit")

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])
