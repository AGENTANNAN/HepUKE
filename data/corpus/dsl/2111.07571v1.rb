DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data_2017  = DatasetManager.real_data.find("708_3097")
jpsi_incMC_2017 = DatasetManager.inclusive_mc.find("708_3097")
qed_data_3080   = DatasetManager.real_data.find("708_3080")
psip_data       = DatasetManager.real_data.find("709_3686")
psip_incMC      = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for psi(3686) -> pi+ pi- J/psi (used for inclusive J/psi detection efficiency)
decay_card_psip_pipi_jpsi = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 gamma e+ e- PHSP;
    Enddecay
    End
DECAYCARD

exMC_psip_pipi_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pipi_jpsi_inclusive"
  config.related_dataset = psip_data
  config.events          = 10_000_000
  config.decay_card      = decay_card_psip_pipi_jpsi
  config.cross_section   = :default
end

# Algorithm for inclusive J/psi event counting
alg_jpsi_count = Algorithm.new("JpsiInclusiveCount")
alg_jpsi_count.set_header(["JpsiInclusiveCountAlg/JpsiInclusiveCount.h"])
               .set_constant({"ECMS" => [:double, 3.097]})
               .note(:inclusive_analysis,
                 "This is an inclusive J/psi selection: all events passing basic quality cuts are counted. " \
                 "No exclusive channel is reconstructed and no kinematic fit is applied. " \
                 "The number of J/psi events N_J/psi = (N_sel - N_bg) / (epsilon_trig * epsilon_data * f_cor).")
               .note(:visible_energy_cut,
                 "E_vis > 1.0 GeV required; E_vis is the sum of charged particle energies (pion mass hypothesis) " \
                 "and neutral EMC shower energies. Applied to suppress QED, cosmic, beam-induced backgrounds.")
               .note(:two_prong_veto,
                 "For events with exactly 2 charged tracks: each track momentum < 1.5 GeV/c and " \
                 "EMC energy per track < 1.0 GeV to reject Bhabha and dimuon events.")
               .note(:track_momentum_cut,
                 "All charged tracks required to have momentum < 2.0 GeV/c.")
               .note(:photon_acceptance,
                 "Photon barrel: |cos(theta)| < 0.83, energy > 25 MeV. " \
                 "Photon endcap: 0.86 < |cos(theta)| < 0.93, energy > 50 MeV. " \
                 "EMC time: [0, 700] ns from event start time.")
               .note(:efficiency_measurement,
                 "Detection efficiency measured from psi(3686) -> pi+ pi- J/psi data: " \
                 "at least two oppositely charged soft pions with p < 0.4 GeV/c, |cos(theta)| < 0.93, " \
                 "|Vz| < 15 cm, Vr < 1 cm. No further selection on remaining tracks/showers. " \
                 "J/psi yield from fit to pi+ pi- recoil mass spectrum.")
               .note(:background_subtraction,
                 "Background estimated from continuum data at sqrt(s)=3.08 GeV: " \
                 "N_bg = N_3.08 * (L_J/psi / L_3.08) * (s_3.08 / s_J/psi).")
               .note(:luminosity_measurement,
                 "Integrated luminosity measured using e+e- -> gamma gamma process: " \
                 "at least 2 EMC showers with |cos(theta)| < 0.8, second-most-energetic shower 1.2-1.6 GeV, " \
                 "|Delta_phi| < 2.5 degrees signal region.")

# Inclusive selection: track and photon quality only (no PID, no channel-specific cuts)
event_sel = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          15.0
    Vr          1.0
    nChrp       ">=1"
    nChrn       ">=1"
    nTot        ">=2"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end

alg_jpsi_count.apply(event_sel)
alg_jpsi_count.execute_on([jpsi_data_2017, jpsi_incMC_2017, qed_data_3080, psip_data, psip_incMC, exMC_psip_pipi_jpsi])