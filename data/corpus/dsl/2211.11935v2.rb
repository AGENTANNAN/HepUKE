# Paper: 2211.11935v2
# Title: Search for eta_c(2S) -> pi+ pi- eta via psi(3686) -> gamma eta_c(2S)
# Energy: 3.686 GeV (single energy, 448M psi(3686) events)
# Final state: gamma(M1) gamma gamma pi+ pi- (eta -> gamma gamma)

### Dataset preparation ###
psi3686_data  = DatasetManager.load_real_data.find("709_3686")
psi3686_incMC = DatasetManager.load_inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  eta_c(2S)   PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000  pi+  pi-  eta      PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma       PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma       PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "etac2s_pipi_eta"
  config.related_dataset = psi3686_data
  config.events         = 100000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("Etac2sToPipiEta")
alg.set_header(["Etac2sToPipiEtaAlg/Etac2sToPipiEta.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})

# Selection: gamma(M1) + eta(gammagamma) + pi+ pi-
# Charged tracks: exactly 2 with net zero charge; |cos(theta)|<0.93; |Vz|<10cm; |Vxy|<1cm
# PID: pion hypothesis greatest likelihood L(pi) > L(K) and L(pi) > L(p)
# Photons: energy > 25MeV barrel, > 40MeV endcap; angle to track > 10deg; EMC time [0,700]ns; at least 3
# 5C kinematic fit: four-momentum conserved + eta mass constraint
# M1 photon is the one not forming the eta in the 5C fit

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                  nTot "==2"
                  nNet "==0"
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.040
                  nGam ">=3"
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]
                  npip "==1"
                  npim "==1"
                }
                .remove([:pip <= :chrgp])
                .remove([:pim <= :chrgn])
                # Build eta from photon pairs via 1C mass-constrained Kalman fit
                .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 200
                  neta ">=1"
                }
                # 5C kinematic fit: eta (already mass-constrained) + M1 gamma + pions
                # four-momentum constraint adds 4C -> total 5C
                .kinematic_fit([:eta, :gamma, :pip, :pim]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

alg.with_decay_card(decay_card).apply(event_selection)

# Note: the paper uses a modified 4C (m4C) fit where the M1 photon energy
# is excluded from the fit to better suppress psi(3686) -> pi+ pi- eta
# background. This is applied in ROOT analysis.
alg.note(:m4c_fit,
  "Modified 4C fit (M1 photon energy excluded) used for final M(pi+pi-eta) spectrum; applied in ROOT")

# Competing hypothesis: chi2_5C(3gamma) < chi2_5C(2gamma) and < chi2_5C(4gamma)
# Suppresses psi(3686) -> pi+ pi- eta and psi(3686) -> gamma gamma pi+ pi- eta
alg.note(:competing_hypothesis_veto,
  "chi2_5C(3g) < chi2_5C(2g) AND < chi2_5C(4g); applied in ROOT analysis")

# J/psi veto: M(gamma(M1) pi+ pi-) <= 3.00 GeV/c^2
alg.note(:jpsi_veto,
  "M(gamma_M1 pi+ pi-) <= 3.00 GeV/c^2 to veto psi(3686) -> eta J/psi background; applied in ROOT")

root_files = alg.execute_on([psi3686_data, psi3686_incMC, exMC_signal])