# BESIII paper 1208.1857v2 — e+e- -> eta J/psi and e+e- -> pi0 J/psi at sqrt(s) = 4.009 GeV
# J/psi -> l+ l- (l = e, mu), eta/pi0 -> gamma gamma, 478 pb^-1 data sample.

### Dataset description ###
data_4009  = DatasetManager.real_data.find("703_4009")     # 478 pb^-1 at sqrt(s) = 4.009 GeV
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")  # inclusive MC at 4.009 GeV
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) accompanying sample (control)

# Signal decay cards: ISR production via psi(4040) -> eta J/psi and psi(4040) -> pi0 J/psi
decay_card_etaJpsi = <<~DECAYCARD
    Decay psi(4040)
    1.0000  eta J/psi                     PHSP;
    Enddecay

    Decay J/psi
    0.5000  mu+ mu-                        PHOTOS VLL;
    0.5000  e+ e-                          PHOTOS VLL;
    Enddecay

    Decay eta
    1.0000  gamma gamma                    PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pi0Jpsi = <<~DECAYCARD
    Decay psi(4040)
    1.0000  pi0 J/psi                      PHSP;
    Enddecay

    Decay J/psi
    1.0000  mu+ mu-                        PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                    PHSP;
    Enddecay

    End
DECAYCARD

exMC_etaJpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psi4040_eta_Jpsi_ll"
  config.related_dataset = data_4009
  config.events         = 200_000
  config.decay_card     = decay_card_etaJpsi
  config.cross_section  = :default
end
exMC_etaJpsi.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_pi0Jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psi4040_pi0_Jpsi_mumu"
  config.related_dataset = data_4009
  config.events         = 200_000
  config.decay_card     = decay_card_pi0Jpsi
  config.cross_section  = :default
end
exMC_pi0Jpsi.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
# eta J/psi and pi0 J/psi share the same final-state topology (l+ l- gamma gamma)
# and identical BOSS-level selection except for the extra mu-depth cut for pi0 J/psi.
alg_name = "EtaJpsi_ll_gg"
alg      = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.009]})
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta  0.93         # |cos theta| < 0.93
                  Vz        10.0          # |Vz| < 10 cm
                  Vr         1.0          # |Vr| < 1 cm
                  nChrp    "==1"          # exactly one positive track
                  nChrn    "==1"          # exactly one negative track
                  nNet     "==0"          # net charge zero
                }
               .select_photon {
                  energyThreshold_b 0.025 # barrel E > 25 MeV
                  energyThreshold_e 0.050 # endcap E > 50 MeV
                  angle_to_track   20.0   # >= 20 deg separation from any charged track
                  tdc_emc_start     0
                  tdc_emc_end      14
                  nGam           "==2"    # exactly two good photons
                }
               # Lepton PID via E/p and EMC energy: both tracks must be either both e or both mu.
               # DSL exposes high-momentum lepton identification via identify_high_momentum_leptons.
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                }
               .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200            # loose BOSS cut; paper uses chi2_4C < 40 in ROOT
                }
               # Second 3C kinematic fit for the radiative-Bhabha / radiative-dimuon veto:
               # the two photons are refit letting one photon's energy float, producing
               # M(gamma_H l+ l-). Stored (no chi2 cut, no nominal) for ROOT-level cut.
               .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
                  constrain_three_momentum
                }

alg.note(:lepton_pid_epcut,
         "Charged-track lepton typing: track with EMC deposit < 0.4 GeV -> muon; track with " \
         "E/p > 0.8 -> electron. Both tracks required to be same species (both mu or both e). " \
         "Combined with FSR/bremsstrahlung correction: photons within 5 deg cone of lepton momentum " \
         "are added back to the lepton four-momentum. Not directly exposed by DSL identify_high_momentum_leptons.")
   .note(:mgg_recoil_window,
         "Recoil mass of the two photons M_recoil(gamma gamma) required in [2.9, 3.4] GeV/c^2 to " \
         "select J/psi candidates. Applied at ROOT level.")
   .note(:mu_counter_depth,
         "For pi0 J/psi search: to reject e+e- -> pi+pi-pi0 background, at least one charged track " \
         "must have muon-counter hit depth > 30 cm. Not expressible in DSL selection block.")
   .note(:radiative_ll_veto,
         "M(gamma_H l+ l-) < 3.93 GeV/c^2 (from 3C fit) required to reject radiative Bhabha and " \
         "radiative dimuon. Applied in ROOT after the second kinematic fit.")

alg.with_decay_card(decay_card_etaJpsi).apply(event_selection)
root_files_eta = alg.execute_on([data_4009, incMC_4009, psip_data, exMC_etaJpsi])

# pi0 J/psi algorithm — same selection chain, restricted to J/psi -> mu+ mu-
alg_pi0 = Algorithm.new("Pi0Jpsi_mumu_gg")
alg_pi0.set_header(["Pi0Jpsi_mumu_ggAlg/Pi0Jpsi_mumu_gg.h"])
       .set_constant({"ECMS" => [:double, 4.009]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_pi0.note(:mu_counter_depth,
             "At least one charged track with muon-counter hit depth > 30 cm to suppress " \
             "e+e- -> pi+ pi- pi0 background. Efficiency 87.9% for signal.")
       .note(:mgg_recoil_window,
             "M_recoil(gamma gamma) in [2.9, 3.4] GeV/c^2 for J/psi. Applied in ROOT.")
alg_pi0.with_decay_card(decay_card_pi0Jpsi).apply(event_selection)
root_files_pi0 = alg_pi0.execute_on([data_4009, incMC_4009, exMC_pi0Jpsi])
