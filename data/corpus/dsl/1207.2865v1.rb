# BESIII paper 1207.2865v1
# Determination of the number of J/psi events collected with BESIII in 2009,
# using J/psi -> inclusive events. This is an inclusive counting analysis;
# there is no exclusive final state and no kinematic fit. The BOSS-side scope
# covers charged-track + photon selection and event-level requirements.

### Dataset description ###
jpsi_data       = DatasetManager.real_data.find("708_3097")     # J/psi data (225.3 M events, 2009 run)
jpsi_incMC      = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC (for correction factor numerator)
psip_data       = DatasetManager.real_data.find("709_3686")     # psi(2S) data (106 M events, for detection efficiency from psi' -> pi+pi- J/psi)
psip_incMC      = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC (for correction factor denominator)

# Exclusive MC: psi' -> pi+ pi- J/psi, J/psi -> inclusive
# Used to determine eps_data^psi' (from psi(2S) data) and eps_mc^psi' (MC counterpart)
decay_card_psipToJpsi = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi+ pi- J/psi          VVPIPI;
    Enddecay

    Decay J/psi
    1.0000  anything               PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: direct J/psi -> inclusive (produced at rest)
# Used to determine eps_mc^{J/psi} (correction factor numerator)
decay_card_JpsiIncl = <<~DECAYCARD
    Decay J/psi
    1.0000  anything               PHSP;
    Enddecay

    End
DECAYCARD

exMC_psipToJpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_to_pipi_Jpsi_inclusive"
  config.related_dataset = psip_data
  config.events         = 1_000_000
  config.decay_card     = decay_card_psipToJpsi
  config.cross_section  = :default
end
exMC_psipToJpsi.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_JpsiIncl = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Jpsi_inclusive_direct"
  config.related_dataset = jpsi_data
  config.events         = 1_000_000
  config.decay_card     = decay_card_JpsiIncl
  config.cross_section  = :default
end
exMC_JpsiIncl.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name        = "JpsiInclusiveCount"
jpsi_algorithm  = Algorithm.new(alg_name)
jpsi_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})

# Selection chain for J/psi -> inclusive events.
# Track: |cos(theta)| < 0.93, |Vz| < 15 cm, Vr < 1 cm, momentum p < 2.0 GeV/c.
# Event level: at least two charged tracks (Ngood >= 2).
# Photon: E > 25 MeV in barrel, E > 50 MeV in endcap; EMC timing 0 < T < 15 (50 ns units).
event_selection = Selection.new
event_selection.select_track {
                  cos_theta    0.93       # |cos(theta)| < 0.93
                  Vz          15.0        # |Vz| < 15 cm
                  Vr           1.0        # |Vr| < 1 cm
                  nTot        ">=2"       # at least two good charged tracks (Ngood >= 2)
                }
               .select_photon {
                  energyThreshold_b 0.025 # barrel EMC threshold 25 MeV (|cos theta| < 0.83)
                  energyThreshold_e 0.050 # endcap EMC threshold 50 MeV (0.86 < |cos theta| < 0.93)
                  tdc_emc_start     0     # EMC cluster timing lower bound
                  tdc_emc_end      14     # EMC cluster timing upper bound (0 < T < 15, 50 ns units)
                }

# The analysis is an inclusive event count; there is no exclusive final state,
# no PID, and no kinematic fit. Remaining event-level requirements (track-momentum
# cut, visible-energy cut, and Bhabha/dimuon vetoes) are captured as notes.
jpsi_algorithm
  .note(:track_momentum_cut,
        "Each good charged track is required to have momentum p < 2.0 GeV/c in the MDC. " \
        "This upper cut is not expressible via select_track cos_theta/Vz/Vr multiplicity keywords.")
  .note(:visible_energy_cut,
        "Event-level cut: E_vis > 1.0 GeV, where E_vis is the sum of charged-track energies " \
        "(computed under the pion-mass hypothesis from track momenta) plus neutral EMC shower energies. " \
        "Applied after track/photon selection to suppress background events.")
  .note(:two_prong_bhabha_dimuon_veto,
        "If the event contains exactly two charged tracks, both track momenta must be below 1.5 GeV/c " \
        "to reject e+e- -> e+e- Bhabha and e+e- -> mu+mu- dimuon events. The clear cluster at " \
        "p ~ 1.55 GeV/c in the two-prong momentum-momentum scatter plot is thereby removed.")
  .note(:emc_deposit_per_track_veto,
        "For each good charged track the EMC deposited energy is required to be less than 1 GeV, " \
        "additionally suppressing Bhabha contamination (a 1.5 GeV peak in the per-track EMC deposit " \
        "spectrum is thereby removed).")
  .note(:efficiency_from_psi2S,
        "Detection efficiency eps_data^{psi'} is determined experimentally from psi(2S) data via " \
        "psi(2S) -> pi+ pi- J/psi, J/psi -> inclusive. Soft-pion selection on the psi(2S) side: " \
        "|cos theta| < 0.93, Vr < 1 cm, |Vz| < 15 cm, p < 0.4 GeV/c; recoil mass of pi+pi- is fit " \
        "with a double Gaussian + 2nd-order Chebychev in [3.07, 3.13] GeV/c^2. Remaining tracks/showers " \
        "must satisfy the same J/psi -> inclusive selection above. Correction factor f_cor = eps_mc^{J/psi}/eps_mc^{psi'} " \
        "is obtained from MC using the same selection chain.")
  .note(:luminosity_from_gammagamma,
        "Integrated luminosities of J/psi and continuum (sqrt(s)=3.08 GeV) data are determined from " \
        "e+e- -> gamma gamma with >=2 neutral tracks, second most energetic shower energy in [1.2, 1.6] GeV, " \
        "|cos theta| < 0.8, and signal region |Delta phi| < 2.5 deg (sideband 2.5 < |Delta phi| < 5 deg).")

jpsi_algorithm.with_decay_card(decay_card_JpsiIncl).apply(event_selection)
root_files = jpsi_algorithm.execute_on([jpsi_data, jpsi_incMC, psip_data, psip_incMC, exMC_JpsiIncl, exMC_psipToJpsi])
