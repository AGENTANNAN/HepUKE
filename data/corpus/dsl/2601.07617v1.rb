# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi data sample (10087 M events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample for J/psi

# Decay card for the signal process J/psi -> Lambda anti-Sigma0 eta (+ c.c.)
decay_card_for_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Sigma0 eta         PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0              PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                          HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                     HypWK;
    Enddecay

    Decay eta
    1.0000 gamma gamma                     PHSP;
    Enddecay

    End
DECAYCARD

# Create an exclusive MC sample for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_to_Lambda_antiSigma0_eta"
  config.related_dataset = jpsi_data
  config.events         = 500000
  config.decay_card     = decay_card_for_signal
  config.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiLSigmaEta"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {                   # Charged track selection
                  cos_theta 0.93                 # |cos(theta)| <= 0.93
                  Vz        100.0
                  Vr        10.0
                  nChrp     ">=2"                # >=2 positive tracks (p+, pi+)
                  nChrn     ">=2"                # >=2 negative tracks (anti-p-, pi-)
                  nNet      "==0"                # net charge zero, >=4 tracks total
                }
               .select_photon {                  # Photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14           # 700 ns / 50 ns per bin
                  angle_to_track    10.0         # opening angle > 10 deg to nearest charged track
                  energyThreshold_b 0.050        # E > 50 MeV in barrel (|cos theta| < 0.8)
                  energyThreshold_e 0.050        # E > 50 MeV in end cap (0.86 < |cos theta| < 0.92)
                  nGam              ">=3"        # >=3 good photons
                }
               .assign({:chrgp => :pip, :chrgn => :pim})  # provisional pion assignment
               .assign({:chrgp => :prp, :chrgn => :prm})  # also treat charged tracks as (anti-)proton hypotheses
               # Reconstruct Lambda -> p+ pi- with a secondary vertex fit
               .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Reconstruct anti-Lambda -> anti-p- pi+ with a secondary vertex fit
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Nominal 5C kinematic fit: Lambda Lambda_bar eta(gamma gamma) gamma
               # Includes 4C energy-momentum conservation + mass constraint on eta -> gamma gamma
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 200
               }
               # Competing hypothesis: Lambda Lambda_bar gamma gamma (background J/psi -> LLbar gg)
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {
                  constrain_four_momentum
               }
               # Competing hypothesis: Lambda Lambda_bar gamma gamma gamma gamma (background J/psi -> LLbar gggg)
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma, :gamma]) {
                  constrain_four_momentum
               }
               # Competing hypothesis with pi0 mass constraint (background J/psi -> LLbar pi0 gamma)
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma]) {
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
               }

# Notes for BOSS-side procedures not directly expressible above
my_Algorithm.note(:lambda_mass_window,
                  "Lambda / anti-Lambda mass window 1.111 < M(p pi) < 1.121 GeV/c^2 applied after secondary vertex fits")
            .note(:sigma0_veto,
                  "To suppress miscombination with the charge-conjugate mode, veto M(gamma Lambda) inside the Sigma0 signal region (1.179, 1.204) GeV/c^2")
            .note(:sigma0_signal_window,
                  "For PWA signal purity, additionally require M(gamma anti-Lambda) in (1.184, 1.199) GeV/c^2")
            .note(:background_veto,
                  "Reject J/psi -> LLbar gamma gamma and J/psi -> LLbar gamma gamma gamma gamma by requiring chi^2(LLbar gggg-hyp) > chi^2(LLbar ggg-hyp) and chi^2(LLbar gg-hyp) > chi^2(LLbar ggg-hyp); also require chi^2_5C(LLbar eta gamma) < chi^2_5C(LLbar pi0 gamma)")

# Render the algorithm and execute
my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
