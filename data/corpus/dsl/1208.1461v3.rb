# BESIII paper 1208.1461v3 — J/psi -> 3 gamma and J/psi -> gamma eta_c, eta_c -> gamma gamma
# via psi(3686) -> pi+ pi- J/psi with 1.0641e8 psi(3686) events.

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) -> pi+ pi- J/psi, J/psi -> 3 gamma
decay_card_3gamma = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi+ pi- J/psi                VVPIPI;
    Enddecay

    Decay J/psi
    1.0000  gamma gamma gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> pi+ pi- J/psi, J/psi -> gamma eta_c, eta_c -> gamma gamma
decay_card_gammaEtac = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi+ pi- J/psi                VVPIPI;
    Enddecay

    Decay J/psi
    1.0000  gamma eta_c                  PHSP;
    Enddecay

    Decay eta_c
    1.0000  gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

exMC_3gamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_pipi_Jpsi_3gamma"
  config.related_dataset = psip_data
  config.events         = 500_000
  config.decay_card     = decay_card_3gamma
  config.cross_section  = :default
end
exMC_3gamma.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_gammaEtac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_pipi_Jpsi_gamma_etac_gammagamma"
  config.related_dataset = psip_data
  config.events         = 500_000
  config.decay_card     = decay_card_gammaEtac
  config.cross_section  = :default
end
exMC_gammaEtac.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "PsiPPiPiJpsi3Gamma"
alg      = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93         # (default, no explicit cut given)
                  Vz         10.0          # |Vz| < 10 cm
                  Vr          1.0          # |Vr| < 1 cm
                  nChrp     "==1"          # exactly two charged tracks: one pi+ ...
                  nChrn     "==1"          # ... one pi-
                  nNet      "==0"
                }
               .select_photon {
                  energyThreshold_b 0.025  # barrel  |cos theta|<0.8, E > 25 MeV
                  energyThreshold_e 0.050  # endcap  0.86<|cos theta|<0.92, E > 50 MeV
                  angle_to_track   5.0     # veto showers within 5 deg of charged tracks
                  tdc_emc_start    0
                  tdc_emc_end     14       # shower time within 700 ns of event start
                  nGam            ">=3"    # at least three photons
                }
               .assign({:chrgp => :pip, :chrgn => :pim})
               # 4C kinematic fit: pi+ pi- gamma gamma gamma constrained to psi(3686) four-momentum.
               # Combination with smallest chi2 kept.
               .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) {
                  nominal
                  vertex_fit([0, 1])       # vertex constraint on pi+ pi-
                  constrain_four_momentum
                  chi2_cut 200             # loose in BOSS; tight chi2_4C < 50 applied in ROOT
               }

alg.note(:photon_count_upper,
         "Only events with 3 or 4 photon candidates are kept (upper multiplicity limit on photons).")
   .note(:jpsi_recoil_window,
         "pi+ pi- recoil mass required to lie in [3.091, 3.103] GeV/c^2 to select J/psi from " \
         "psi(3686) -> pi+ pi- J/psi. Applied in ROOT after fit.")
   .note(:pi0_eta_etaprime_veto,
         "Events with any di-photon mass M(gamma gamma) in [0.10, 0.16] (pi0), [0.50, 0.60] (eta), " \
         "or [0.90, 1.00] GeV/c^2 (eta') are removed. Applied in ROOT.")
   .note(:gamma_ee_bkg_veto,
         "To reject J/psi -> gamma e+ e- background, number of MDC hits within an opening angle " \
         "of five EMC crystals around each photon center is counted; total from three photons required < 40. " \
         "Not expressible in DSL; applied in BOSS user code.")

alg.with_decay_card(decay_card_3gamma).apply(event_selection)
root_files = alg.execute_on([psip_data, psip_incMC, exMC_3gamma, exMC_gammaEtac])
