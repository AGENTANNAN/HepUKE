# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# BOSS 707 data at the three energy points 4.84, 4.91 and 4.95 GeV
data_4840  = DatasetManager.real_data.find("707_4840")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
data_4914  = DatasetManager.real_data.find("707_4914")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
data_4946  = DatasetManager.real_data.find("707_4946")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

### Decay card of the signal chain: e+e- -> K+K- psi(3770), psi(3770) -> D0 anti-D0 ###
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  K+  K-  psi(3770)   PHSP;
    Enddecay

    Decay psi(3770)
    1.0000  D0  anti-D0         PHSP;
    Enddecay

    Decay D0
    1.0000  K-  pi+             PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-             PHSP;
    Enddecay

    End
DECAYCARD

### 200k exclusive-MC events of the signal chain, one sample per energy point ###
exMCs_signal = DatasetManager.create_exclusive_mc_for([data_4840, data_4914, data_4946]) do |config|
  config.sample_name   = "exmc_kk_psi3770_d0d0bar"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag analysis ###
# The reconstructed D meson is the tag; the K+K- pair and the recoiling (unreconstructed)
# anti-D0 form the signal side. Each energy point is run separately.
energy_points = [
  # label,  ECMS (GeV), real data, inclusive MC, exclusive MC
  ["4840",  4.840,     data_4840,  incMC_4840,   exMCs_signal[0]],
  ["4914",  4.914,     data_4914,  incMC_4914,   exMCs_signal[1]],
  ["4946",  4.946,     data_4946,  incMC_4946,   exMCs_signal[2]]
]

energy_points.each do |label, ecms, real_data, inc_mc, ex_mc|
  alg = TagAnalysis.new("KKpsi3770Tag#{label}")
  alg.set_header(["KKpsi3770Tag#{label}Alg/KKpsi3770Tag#{label}.h"])
     .set_constant({"ECMS" => [:double, ecms]})
     .with_decay_card(decay_card_signal)
     # Nine tag modes are used for the D reconstruction. Only the D0 modes can be declared on a
     # single D0 tag side; the D- modes (DptoKPiPi, DptoKPiPiPi0, DptoKsPi, DptoKsPiPi0,
     # DptoKsPiPiPi, DptoKKPi) require a separate tag side and are recorded here.
     .note(:tag_mode_list, "D reconstruction uses nine tag modes: D0 -> K+pi-, K+pi-pi0, K+pi+pi-pi- " \
                           "and D- -> K+pi-pi-, K+pi-pi-pi0, KS pi-, KS pi pi0, KS pi-pi-pi+, K+K-pi-; " \
                           "the 1C mass-constrained fit of the tagged D (chi2 < 13) and the choice of the " \
                           "D candidate closest to the nominal D mass are performed inside DTagAlg")
     # BOSS-side selection of the bachelor kaons.
     .note(:bachelor_kaon_selection, "the bachelor kaons are taken as the lowest-momentum K+K- pair among " \
                                     "the charged tracks that are not used by the D tag")
     # Physics veto that is not part of the tag fit vocabulary (no out_of on a tag fit).
     .note(:background_veto, "events with |M(K+K-) - m_phi| < 0.02 GeV/c^2 are vetoed to remove the " \
                             "phi(1020) -> K+K- background; window applied on the stored M(K+K-) in ROOT")
     # Standard BESIII track/photon/PID baseline carried by the tag reconstruction.
     .note(:track_and_photon_selection, "|cos(theta)| < 0.93, |Vz| < 10 cm, Vr < 1 cm, at least two positive " \
                                        "and two negative tracks; at least two photons with EMC timing 0-14, " \
                                        "angle to the nearest charged track > 10 deg, E > 25 MeV (barrel) / " \
                                        "50 MeV (endcap); PID by the probability method with prob > 0.001 for " \
                                        "K/pi separation")

  # ---- tag side: one tagged D0; both tag charges are scanned ----
  alg.tag_side(:D0) do |t|
    t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  end

  # ---- signal side: the K+K- left over by the tag plus the unreconstructed anti-D0 ----
  alg.signal_side do |s|
    s.photons 2                 # at least two good showers entering the fit (pi0 -> gamma gamma)
    s.min_photon_angle 10.0     # minimum angle to the nearest charged track (degrees)
    s.min_photon_energy 0.025   # shower energy floor: 25 MeV (barrel) / 50 MeV (endcap)
    s.charged(kp: 1, km: 1)     # the two bachelor kaons K+K-
    s.missing :D0               # the recoiling anti-D0 is not reconstructed -> RM(K+K-D)
  end

  # ---- kinematic fit over tag + K+K- + gamma gamma + missing D0 ----
  alg.fit do |f|
    f.constrain_four_momentum                                        # 4C energy-momentum constraint
    f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 mass (chi2 < 25)
    f.chi2_cut 200
  end

  alg.apply                       # no Selection argument for a tag analysis
  alg.execute_on([real_data, inc_mc, ex_mc])
end