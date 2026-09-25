# Dataset preparation
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi                     PHSP;
  Enddecay

  Decay J/psi
  1.0000 gamma p+ anti-p-                  PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pipi_Jpsi_gamma_ppbar"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection (BOSS)
alg_name = "PsipPiPiJpsiGammaPPbar"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93          # |cos(theta)| < 0.93
                  Vz        100.0
                  Vr        10.0
                  nChrp     ">=2"
                  nChrn     ">=2"
                  nNet      "==0"          # total net charge zero
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0   # isolated from all charged tracks by > 10 deg
                  energyThreshold_b 0.025  # >=25 MeV in barrel
                  energyThreshold_e 0.050  # >=50 MeV in endcap
                  nGam              ">=1"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]
                  nprp ">=1"
                  nprm ">=1"
                }
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               .select_isolated_photon {
                  angle_to_prm_track 30.0  # isolated from anti-proton track by > 30 deg
                  nGam ">=1"
                }
               .assign({:chrgp => :pip, :chrgn => :pim})
               .kinematic_fit([:gamma, :pip, :pim, :prp, :prm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 100             # chi^2 < 100 as required in the paper
                }

algorithm
  .note(:jpsi_recoil_window,
        "|M(pi+pi-)_recoil - m_J/psi| < 0.006 GeV/c^2 applied to J/psi tag; " \
        "handled at ROOT level from recorded recoil mass.")
  .note(:multi_photon_suppression,
        "|U_miss| = |E_miss - |P_miss|| < 0.05 GeV and P_t_gamma^2 = " \
        "4|P_miss|^2 sin^2(theta_gamma/2) < 0.0005 (GeV/c)^2, applied " \
        "before the 4C fit using missing-energy/momentum of the four charged tracks.")
  .note(:psip_pipi_ppbar_veto,
        "|M(pi+pi-ppbar) - m_psi'| > 0.03 GeV/c^2 veto against " \
        "psi' -> pi+ pi- p pbar background.")
  .note(:min_proton_momentum,
        "Reject events with p_p < 0.3 GeV/c or p_pbar < 0.3 GeV/c to " \
        "avoid low-momentum tracking efficiency mismatches between data and MC.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_signal])
