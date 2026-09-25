# arXiv:1304.1973v1 — Partial wave analysis of psi(2S) -> p pbar eta (BESIII, 1.06e8 psi(2S))
# Final state: psi(2S) -> p pbar eta, eta -> gamma gamma  (topology p p gamma gamma)

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # 1.06e8 psi(2S) events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # 1e8 inclusive psi(2S) MC events
# 42.6 pb^-1 of continuum data at 3.65 GeV used for the QED background estimate
cont_data  = DatasetManager.real_data.find("709_3650")

# Decay card for the signal process
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- eta    PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_ppbar_eta"
    config.related_dataset = psip_data
    config.events          = 100000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                    cos_theta       0.93    # |cos(theta)| < 0.93
                    Vz              20.0    # within +/-20 cm of the beamline along the beam
                    Vr              2.0     # within 2 cm in the plane perpendicular to the beam
                    nChrp           "==1"   # two good charged tracks with total charge zero
                    nChrn           "==1"
                    nNet            "==0"
                }
               .pid(method: :probability) {
                    prob_cut   0.001
                    # The particle type with the highest probability is assigned to each
                    # track; one track must be a proton and the other an anti-proton.
                    identify :proton, against: [:pion, :kaon]
                    nprp   "==1"
                    nprm   "==1"
               }
               # Photon candidates: showers near the proton (anti-proton) track are
               # rejected by 10 (30) degrees isolation applied below.
               .select_photon {
                    tdc_emc_start     0
                    tdc_emc_end       14      # timing requirement against noise
                    energyThreshold_b 0.025   # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
                    energyThreshold_e 0.050   # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
                    angle_to_track    10.0
                    nGam              ">=2"   # at least two good photons
               }
               .select_isolated_photon {
                    angle_to_prp_track  10.0  # delta angle to the nearest proton track > 10 deg
                    angle_to_prm_track  30.0  # delta angle to the nearest anti-proton track > 30 deg
                    nGam                ">=2"
               }
               # 4C kinematic fit under the p pbar gamma gamma hypothesis; if more than
               # two photons are selected the fit is repeated over all photon permutations
               # and the combination with the minimum chi^2 is taken.
               .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
                    nominal
                    constrain_four_momentum   # 4C energy-momentum conservation
                    chi2_cut 20
               }

my_algorithm
  .note(:efficiency_curve, "the momenta of the proton and anti-proton are required to
    be greater than 300 MeV/c because data and MC simulation do not agree well in the
    low-momentum region; applied in the ROOT analysis after the kinematic fit")
  .note(:background_veto, "psi(2S) -> eta J/psi and psi(2S) -> gamma chi_cJ backgrounds
    removed by requiring M(p pbar) < 3.067 GeV/c^2 and
    M(p pbar) < (3.1 GeV/c^2 - 0.75 * M(gamma gamma)) (applied in the ROOT analysis)")
  .note(:background_veto, "eta mass window |M(gamma gamma) - M_eta| < 21 MeV/c^2;
    non-eta background estimated from the eta sidebands
    |M(gamma gamma) - 0.43| < 50 MeV/c^2 and |M(gamma gamma) - 0.65| < 50 MeV/c^2, and
    QED background estimated from 42.6 pb^-1 of continuum data at 3.65 GeV")
  .note(:background_veto, "residual background from psi(2S) -> gamma chi_cJ with
    chi_cJ -> p pbar pi0, chi_cJ -> gamma J/psi (J/psi -> gamma p pbar) and
    chi_c0 -> pbar Delta+ (Delta+ -> p pi0) estimated with inclusive psi(2S) MC")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

root_files = my_algorithm.execute_on([psip_data, psip_incMC, cont_data, exMC_signal])
