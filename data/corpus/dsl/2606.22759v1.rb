### dataset description ###
# ConExc-based Born cross-section measurement of e+e- -> p pbar at 47 c.m. energies
# between 3.510 GeV and 4.946 GeV (total integrated luminosity 26 fb^-1). Representative
# scan points are listed; the full 47-point scan should be added to the datasets array.
data_scan_points = [
  DatasetManager.real_data.find("703_4260"),  # Y(4260) / 4.260 GeV
  DatasetManager.real_data.find("712_3773"),  # psi(3770) / 3.773 GeV
]

incMC_samples = data_scan_points.map do |ds|
  DatasetManager.inclusive_mc.find(ds.sample_name.nil? ? ds.name : "#{ds.boss.gsub('.', '')[0,3]}_#{ds.energy.to_i}")
end

# ConExc decay card for the continuum process e+e- -> p pbar (mode 0). 'Particle vpho'
# is omitted so the DSL injects the correct sqrt(s) at every scan point.
decay_card_ppbar = <<~DECAYCARD
  Decay vpho
  1 ConExc 0;
  Enddecay
  Decay vhdr
  1 p+ anti-p-  PHSP;
  Enddecay
  End
DECAYCARD

# Multi-energy ConExc signal MC (one exclusive MC per scan point)
exMC_ppbar = DatasetManager.create_exclusive_mc_for(data_scan_points) do |config|
  config.sample_name   = "sig_conexc_ppbar"   # auto-suffixed by dataset key
  config.events        = 100_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default              # required; inert for ConExc
end

### event selection (BOSS) ###
alg_name = "PPbarConExc"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93       # |cos(theta)| < 0.93
                  Vz          10.0       # |Vz| < 10 cm
                  Vr          1.0        # |Vxy| < 1 cm
                  nChrp       "==1"      # exactly one proton candidate
                  nChrn       "==1"      # exactly one anti-proton candidate
                  nNet        "==0"      # net charge = 0
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ AND p-bar via charge-conjugation shorthand
                  nprp   "==1"
                  nprm   "==1"
               }
               # Nominal kinematic fit for e+e- -> p pbar; 4-momentum conservation to the
               # per-run measured CMS energy from MeasuredEcmsSvc.
               .kinematic_fit([:prp, :prm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
               }

my_algorithm
  .note(:proton_E_over_p, "positive proton track required to satisfy E/p < 0.5, where E is the EMC deposit and p the track momentum")
  .note(:back_to_back_cut, "opening angle between the positive and negative tracks required to satisfy theta_opening > 3.1 rad, exploiting the two-body kinematics of e+e- -> p pbar")
  .note(:muon_depth_cut, "hit depth in the muon counter required to be < 40 cm for both tracks to reject cosmic-ray backgrounds")
  .note(:isr_correction, "efficiency and ISR/vacuum-polarization correction factors are obtained iteratively from ConExc-simulated MC per c.m. energy following Ref. [Sun et al., Front. Phys. 16, 64501 (2021)]")
  .note(:full_scan_datasets, "the analysis uses 47 c.m. energies from 3.510 to 4.946 GeV totaling 26 fb^-1; only two representative scan points are listed in data_scan_points and the corresponding exclusive MC array — expand the array to include the full 47-point scan for production running")
  .with_decay_card(decay_card_ppbar)
  .apply(event_selection)

root_files = my_algorithm.execute_on(data_scan_points + incMC_samples + exMC_ppbar)
