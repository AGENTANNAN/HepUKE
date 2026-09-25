### Dataset description ###
# psi(2S) scan: 9 CM energy points from 3.581 to 3.710 GeV, total ~495 pb-1
scan_points = [
  DatasetManager.real_data.find("704_psip_scan_1"),   # 3581.5 MeV
  DatasetManager.real_data.find("704_psip_scan_2"),   # 3670.2 MeV
  DatasetManager.real_data.find("704_psip_scan_3"),   # 3680.1 MeV
  DatasetManager.real_data.find("704_psip_scan_4"),   # 3682.8 MeV
  DatasetManager.real_data.find("704_psip_scan_5"),   # 3684.2 MeV
  DatasetManager.real_data.find("704_psip_scan_6"),   # 3685.3 MeV
  DatasetManager.real_data.find("704_psip_scan_7"),   # 3686.5 MeV
  DatasetManager.real_data.find("704_psip_scan_8"),   # 3691.4 MeV
  DatasetManager.real_data.find("704_psip_scan_9"),   # 3709.8 MeV
]

# ConExc signal MC for e+e- -> p pbar (mode 0). Omit 'Particle vpho' so the DSL
# injects the correct sqrt(s) at each scan point.
decay_card_ppbar = <<~DECAYCARD
    Decay vpho
    1 ConExc 0;
    Enddecay
    Decay vhdr
    1 p+ anti-p-  PHSP;
    Enddecay
    End
DECAYCARD

# Batch exclusive MC for all nine scan points
exMC_ppbar = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "conexc_ppbar_psipscan"
  config.events        = 100_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PsipScanPpbar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.80          # |cos(theta)| < 0.80 (barrel MDC only)
                  Vz        10.0
                  Vr        1.0
                  nChrp     "==1"         # exactly one positive track
                  nChrn     "==1"         # exactly one negative track
                  nNet      "==0"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # p and pbar; L(p) largest
                  nprp "==1"
                  nprm "==1"
                }

my_algorithm
  .note(:eop_cut,
        "E/P < 0.5 required for both proton and antiproton candidates " \
        "(reject Bhabha and dimuon contamination). Applied on raw EMC/MDC quantities.")
  .note(:cosmic_veto,
        "Cosmic-ray suppression: |T_p - T_pbar| < 3 ns using TOF measurements.")
  .note(:back_to_back_cut,
        "Opening angle between p and pbar in the e+e- CM frame required to be " \
        "in [178, 180] degrees to reject events with additional tracks.")
  .note(:signal_yield_extraction,
        "Signal yield extracted by counting events in a symmetric 6-sigma window " \
        "[p_mean - 6 sigma, p_mean + 6 sigma] around p_mean = sqrt(s^2/4 - m_p^2), " \
        "with sigma ~ 15 MeV from a double-Gaussian fit to signal MC.")

my_algorithm.with_decay_card(decay_card_ppbar).apply(event_selection)
root_files = my_algorithm.execute_on(scan_points + exMC_ppbar)
