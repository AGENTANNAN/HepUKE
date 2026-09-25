### dataset description ###
# 47 c.m. energies between sqrt(s) = 3.510 and 4.946 GeV; total 26 fb^-1.
# Representative scan points are listed; expand to include all 47 points for
# production running.
data_scan_points = [
  DatasetManager.real_data.find("712_3773"),  # psi(3770) region
  DatasetManager.real_data.find("703_4178"),  # 4.178 GeV
  DatasetManager.real_data.find("703_4226"),  # 4.226 GeV
  DatasetManager.real_data.find("703_4260"),  # 4.260 GeV
]

incMC_samples = data_scan_points.map do |ds|
  DatasetManager.inclusive_mc.find("#{ds.boss.gsub('.', '')[0,3]}_#{ds.energy.to_i}")
end

# ConExc decay card for continuum e+e- -> p pbar (mode 0). 'Particle vpho' is
# omitted; the DSL injects the correct sqrt(s) per scan point.
decay_card_ppbar = <<~DECAYCARD
  Decay vpho
  1 ConExc 0;
  Enddecay
  Decay vhdr
  1 p+ anti-p-  PHSP;
  Enddecay
  End
DECAYCARD

exMC_ppbar = DatasetManager.create_exclusive_mc_for(data_scan_points) do |config|
  config.sample_name   = "sig_conexc_ppbar_scan"
  config.events        = 100_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### event selection (BOSS) ###
alg_name = "PPbarScan"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.80       # |cos(theta)| < 0.80 (barrel-only)
                  Vz          10.0       # |Vz| < 10 cm
                  Vr          1.0        # |Vr| < 1 cm
                  nChrp       "==1"      # one positive track
                  nChrn       "==1"      # one negative track
                  nNet        "==0"      # net charge = 0
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # L_p > 0.001, L_p > L_K, L_p > L_pi
                  nprp   "==1"
                  nprm   "==1"
               }
               # Nominal kinematic fit for e+e- -> p pbar
               .kinematic_fit([:prp, :prm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
               }

my_algorithm
  .note(:proton_E_over_p, "positive proton track required to satisfy E/p < 0.5 to suppress Bhabha scattering, where E is the EMC deposit and p the track momentum")
  .note(:back_to_back_cut, "opening angle between p and pbar in the c.m. system required to satisfy theta_opening > 3.1 rad")
  .note(:muon_depth_cut, "hit depth in the muon counter required to be < 40 cm for both tracks to remove e+e- -> mu+ mu- background")
  .note(:antiproton_momentum_window, "for anti-proton candidates, |p_measured - p_expected| < 3 sigma_p required")
  .note(:isr_correction, "ISR correction factor (1+delta) and detection efficiency obtained by iterative MC weighting per c.m. energy following Ref. [Sun et al., Front. Phys. 16, 64501 (2021)]; iteration stops when Delta((1+delta)*eps) < 0.1%")
  .note(:full_scan_datasets, "the analysis uses 47 c.m. energies from 3.510 to 4.946 GeV totaling 26 fb^-1; expand data_scan_points and exMC_ppbar to cover the full point list for production running")
  .with_decay_card(decay_card_ppbar)
  .apply(event_selection)

root_files = my_algorithm.execute_on(data_scan_points + incMC_samples + exMC_ppbar)
