# Dataset preparation
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi                     PHSP;
  Enddecay

  Decay J/psi
  1.0000 gamma eta_c                       PHSP;
  Enddecay

  Decay eta_c
  1.0000 gamma gamma                       PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pipi_Jpsi_gamma_etac_gammagamma"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection (BOSS)
alg_name = "PsipPiPiJpsiGammaEtacGammaGamma"
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
                  identify :pion, against: [:electron]
                  npip ">=1"
                  npim ">=1"
                }
               .assign({:chrgp => :pip, :chrgn => :pim})
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=3"
                }
               .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) do
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                end

algorithm
  .note(:lepton_veto,
        "Reject events where any pion track satisfies L(e)/(L(e)+L(pi)) > 0.8; " \
        "applied at ROOT level.")
  .note(:jpsi_recoil_window,
        "|M(pi+pi-)_recoil - m_J/psi| < 0.018 GeV/c^2; " \
        "applied at ROOT level to tag J/psi.")
  .note(:etac_pair_selection,
        "Among the 3 photons in the 4C fit, the gamma-gamma pair with invariant " \
        "mass closest to m_eta_c is selected; the remaining photon is the J/psi " \
        "radiative photon. Applied at ROOT level.")
  .note(:m23_m13_veto,
        "Veto events where any two-photon invariant mass falls in windows " \
        "[0.10, 0.16] (pi0), [0.48, 0.59] (eta), or [0.88, 1.10] (eta') GeV/c^2. " \
        "Applied at ROOT level.")
  .note(:chi2_cut_final,
        "The paper uses chi2_4C < 19 for the final selection; applied at ROOT level.")
  .note(:nGam_upper,
        "Photon multiplicity is additionally constrained to <= 4 at ROOT level.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_signal])