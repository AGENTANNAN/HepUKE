### Dataset preparation ###
# Four chi_c1 scan points (BOSS release 703); they correspond to the
# chi_c1_scan_2 .. chi_c1_scan_5 entries of the BESIII dataset table.
chic1_data = [
  DatasetManager.real_data.find("703_3508"),   # ~3.5080 GeV
  DatasetManager.real_data.find("703_3510"),   # ~3.5097 GeV
  DatasetManager.real_data.find("703_3511"),   # ~3.5106 GeV
  DatasetManager.real_data.find("703_3514")    # ~3.5144 GeV
]

chic1_incMC = [
  DatasetManager.inclusive_mc.find("703_3508"),
  DatasetManager.inclusive_mc.find("703_3510"),
  DatasetManager.inclusive_mc.find("703_3511"),
  DatasetManager.inclusive_mc.find("703_3514")
]

# Decay card for the signal process (EvtGen format):
# psi(4260) -> gamma J/psi, J/psi -> mu+ mu-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC generated at every scan point (shared decay card)
exMC_signal = DatasetManager.create_exclusive_mc_for(chic1_data) do |config|
  config.sample_name   = "exmc_chic1_gamma_jpsi_mumu"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "ChiC1ToGammaJpsiMumu"
chic1_alg = Algorithm.new(alg_name)
chic1_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.5106]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:isr_interference, "signal/ISR interference for e+e- -> gamma J/psi is modelled with PHOKHARA at the generator level; this cannot be expressed in the EvtGen decay card")
         .note(:background_veto, "the best photon, selected by the minimum chi2_4C of the nominal 4C fit, is required to have |cos(theta_gamma)| < 0.80 to suppress the ISR background; the cut is applied after the nominal kinematic fit")

# One common selection chain for every scan-point data set, inclusive MC and signal MC
event_selection = Selection.new
event_selection
  .select_track {                                # exactly two oppositely charged tracks
    cos_theta 0.93                               # |cos(theta)| < 0.93
    Vz        10.0                               # |Vz| < 10 cm
    Vr        1.0                                # Vr < 1 cm
    nChrp     "==1"                              # one positive track
    nChrn     "==1"                              # one negative track
    nNet      "==0"                              # net charge zero
  }
  .select_photon {                               # at least one good photon
    tdc_emc_start     0                          # EMC TDC window start
    tdc_emc_end       14                         # EMC TDC window end
    energyThreshold_b 0.025                      # 25 MeV (barrel)
    energyThreshold_e 0.050                      # 50 MeV (endcap)
    angle_to_track    10.0                       # opening angle to nearest track > 10 deg
    nGam              ">=1"                      # at least one photon
  }
  .pid(method: :probability) {                   # mu+ / mu- identification, probability method
    prob_cut 0.001                               # PID probability > 0.001
    identify :muon, against: [:pion]             # identify mu+ and mu- against pions
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.4  # reject muon candidates with EMC energy > 0.4 GeV
    nmup "==1"
    nmum "==1"
  }
  .kinematic_fit([:gamma, :mup, :mum]) {         # nominal 4C fit; the best (smallest chi2) photon is chosen automatically
    nominal                                      # nominal fit: its corrected four-momenta are stored
    constrain_four_momentum                      # 4C energy-momentum constraint
    chi2_cut 200                                 # chi2 < 200 (loose; tight cut optimised in ROOT)
  }

chic1_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Run on all four data points, their inclusive MC and the exclusive signal MC
root_files = chic1_alg.execute_on(chic1_data + chic1_incMC + exMC_signal)