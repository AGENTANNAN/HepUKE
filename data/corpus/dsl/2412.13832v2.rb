# Dataset preparation
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_cJ                      PHSP;
  Enddecay

  Decay chi_cJ
  1.0000 p+ anti-p- eta pi0                PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma                       PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                       PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chicJ_ppbar_eta_pi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection (BOSS)
alg_name = "PsipGammaChicJPpbarEtaPi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     "==1"
                  nChrn     "==1"
                  nNet      "==0"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]
                  nprp ">=1"
                  nprm ">=1"
                }
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=5"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) do
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 25
                 npi0 ">=1"
               end
               .kalman_kinematic_fit([:gamma, :gamma]) do
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                 chi2_cut 25
                 neta ">=1"
               end
               .kinematic_fit([:prp, :prm, :pi0, :eta]) do
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               end

algorithm
  .note(:photon_combinatorics,
        "From >=5 photons, the best pi0->gamma gamma and eta->gamma gamma pairs " \
        "are selected via kalman_kinematic_fit. The remaining photon (if any) " \
        "is the radiative photon from psi(2S)->gamma chi_cJ. " \
        "Combinatorial optimization handled at ROOT level.")
  .note(:chi2_cut_final,
        "The paper uses chi2_6C < 20 for the final selection. " \
        "The 6C arises from 4C momentum conservation + 1C pi0 mass + 1C eta mass constraints. " \
        "Applied at ROOT level after the kalman + 4C kinematic fit.")
  .note(:competing_4c_veto,
        "chi2_4C(5gamma p pbar) < chi2_4C(4gamma p pbar) AND " \
        "chi2_4C(5gamma p pbar) < chi2_4C(6gamma p pbar); " \
        "applied at ROOT level to select the correct photon multiplicity.")
  .note(:jpsi_veto,
        "|RM(eta) - m_J/psi| > 7.2 MeV/c^2 and " \
        "|M(p pbar eta) - m_J/psi| > 21.3 MeV/c^2; " \
        "applied at ROOT level to veto J/psi backgrounds.")
  .note(:background_veto,
        "Veto against gamma p pbar eta eta, p pbar eta eta, gamma p pbar pi0 pi0, " \
        "p pbar pi0 pi0, and gamma gamma p pbar pi0 pi0 backgrounds. " \
        "Applied at ROOT level via invariant mass windows.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_signal])