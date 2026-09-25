# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# XYZ data samples from 4.226 to 4.950 GeV (29 energy points). Use the Y(4260) sample as the anchor
# and rely on create_exclusive_mc_for for the multi-energy scan.
xyz_energies = ["703_4260", "703_4360", "703_4420", "703_4600", "703_4680", "703_4750", "703_4840", "703_4950"]
xyz_datasets = xyz_energies.map { |e| DatasetManager.real_data.find(e) }
xyz_incMCs   = xyz_energies.map { |e| DatasetManager.inclusive_mc.find(e) }

# Signal decay card: e+e- -> eta eta J/psi, using leptonic J/psi decays and both eta -> gg / eta -> pi+ pi- pi0
decay_card_for_signal = <<~DECAYCARD
    Alias another_eta eta

    Decay psi(4260)
    1.0000 eta another_eta J/psi            PHSP;
    Enddecay

    Decay J/psi
    0.5000 e+ e-                            PHOTOS VLL;
    0.5000 mu+ mu-                          PHOTOS VLL;
    Enddecay

    Decay eta
    1.0000 gamma gamma                      PHSP;
    Enddecay

    Decay another_eta
    1.0000 pi+ pi- pi0                      ETA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                      PHSP;
    Enddecay

    End
DECAYCARD

# Multi-energy exclusive MC production for the scan
exMC_signal_list = DatasetManager.create_exclusive_mc_for(xyz_datasets) do |config|
  config.sample_name    = "ee_to_eta_eta_Jpsi"
  config.events         = 100000
  config.decay_card     = decay_card_for_signal
  config.cross_section  = :straight_line
end
exMC_signal_list.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "EtaEtaJpsi"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                    # Charged track selection
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     ">=1"
                  nChrn     ">=1"
                  nNet      "==0"
                }
               .select_photon {                   # Photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14            # 700 ns
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=4"         # at least 4 photons (2 for each eta -> gg)
                }
               .pid(method: :probability) {       # PID for pions + high-momentum leptons (e/mu)
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
               }
               # Reconstruct pi0 candidates (from eta -> pi+ pi- pi0 sub-mode)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0    ">=0"
               }
               # Reconstruct eta -> gamma gamma candidates via 1C Kalman fit
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 25
                  neta    ">=1"
               }
               # Nominal 7C kinematic fit: two etas (-> gg) + lepton pair.
               # Constrain total 4-momentum, M(gamma gamma) to eta (twice), M(l+l-) to J/psi.
               .kinematic_fit([:eta, :eta, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:"J/psi")
                  chi2_cut 200
               }
               # Competing background hypothesis: e+e- -> pi0 pi0 J/psi (four photons falsely assigned to two etas)
               .kinematic_fit([:pi0, :pi0, :lp, :lm]) {
                  constrain_four_momentum
                  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:"J/psi")
               }

my_Algorithm.note(:analysis_mode,
                  "Exclusive analysis of e+e- -> eta eta J/psi with J/psi -> l+l- (l = e or mu). Both etas reconstructed exclusively via eta -> gamma gamma and eta -> pi+ pi- pi0.")
            .note(:lepton_momentum_cut,
                  "Lepton momentum p(l) > 1 GeV/c; pion momentum p(pi) < 1 GeV/c; e/mu separation via E/p (>0.7 for electrons, <0.3 for muons).")
            .note(:pi0_veto,
                  "In the eta eta -> 4 gamma sub-mode, veto events where the four selected photons can also form two pi0 candidates; compare the sum of chi^2 from the two-pi0 vs two-eta 1C kinematic fits and take the hypothesis with the smaller sum (suppresses ~93% of pi0 pi0 J/psi background).")
            .note(:semi_inclusive_method,
                  "A complementary semi-inclusive analysis reconstructs J/psi and one eta -> gg and treats the second eta as missing, aggregating the missing-mass to extract the signal. Only the exclusive procedure is captured in this DSL spec.")
            .note(:mass_windows,
                  "Pre-selection windows: 80 <= M(gg) <= 180 MeV/c^2 for pi0; 400 <= M(gg) <= 700 MeV/c^2 for eta -> gg; 400 <= M(pi+pi-pi0) <= 700 MeV/c^2 for eta -> 3pi.")

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on(xyz_datasets + xyz_incMCs + exMC_signal_list)
