# Paper: 2308.08161v2
# e+e- → eta phi at 23 energy points from 3.773 to 4.600 GeV
# phi → K+K-, eta → gamma gamma
# ConExc generator (mode 23: phi eta), cross-section measurement

### Dataset preparation ###
# 23 energy points
scan_datasets = [
  DatasetManager.load_real_data.find("712_3773"),   # 3.773 GeV
  DatasetManager.load_real_data.find("703_4009"),   # 4.008 GeV
  DatasetManager.load_real_data.find("703_4090"),   # 4.086 GeV
  DatasetManager.load_real_data.find("703_4180"),   # 4.178 GeV
  DatasetManager.load_real_data.find("703_4190"),   # 4.189 GeV
  DatasetManager.load_real_data.find("703_4200"),   # 4.200 GeV
  DatasetManager.load_real_data.find("703_4210"),   # 4.210 GeV
  DatasetManager.load_real_data.find("703_4220"),   # 4.219 GeV
  DatasetManager.load_real_data.find("703_4230"),   # 4.226 GeV
  DatasetManager.load_real_data.find("703_4237"),   # 4.236 GeV
  DatasetManager.load_real_data.find("703_4245"),   # 4.242 GeV
  DatasetManager.load_real_data.find("703_4246"),   # 4.244 GeV
  DatasetManager.load_real_data.find("703_4260"),   # 4.258 GeV
  DatasetManager.load_real_data.find("703_4270"),   # 4.267 GeV
  DatasetManager.load_real_data.find("703_4280"),   # 4.278 GeV
  DatasetManager.load_real_data.find("703_4310"),   # 4.308 GeV
  DatasetManager.load_real_data.find("703_4360"),   # 4.358 GeV
  DatasetManager.load_real_data.find("703_4390"),   # 4.387 GeV
  DatasetManager.load_real_data.find("703_4420"),   # 4.416 GeV
  DatasetManager.load_real_data.find("703_4470"),   # 4.467 GeV
  DatasetManager.load_real_data.find("703_4530"),   # 4.527 GeV
  DatasetManager.load_real_data.find("703_4575"),   # 4.575 GeV
  DatasetManager.load_real_data.find("703_4600"),   # 4.600 GeV
]

scan_incMC = scan_datasets.map { |ds| DatasetManager.load_inclusive_mc.find(ds.sample_name) rescue nil }.compact

# ConExc decay card: mode 23 = phi eta
# Particle vpho is auto-injected per energy point by DSL for create_exclusive_mc_for
conexc_card = <<~DECAYCARD
    Decay vpho
    1 ConExc 23;
    Enddecay
    Decay vhdr
    1 phi eta PHSP;
    Enddecay
DECAYCARD

# Exclusive MC for the scan: shared ConExc card, auto-vpho per point
exMC_scan = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name = "ee_to_eta_phi_scan"
  config.events = 500000
  config.decay_card = conexc_card
  config.cross_section = :default   # required but inert for ConExc
end

### Event selection (BOSS) ###
# Two charged tracks (kaons), at least 2 photons
# 4C kinematic fit, phi mass window, opening angle cut

alg_name = "EEtoEtaPhi"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})
# No ECMS constant — multi-energy scan

event_selection = Selection.new
event_selection.select_track {
              cos_theta   0.93
              Vz   10.0
              Vr   1.0
              nChrp   "==1"
              nChrn   "==1"
              nTot    "==2"
              nNet    "==0"
            }
           .select_photon {
              tdc_emc_start   0
              tdc_emc_end     700
              energyThreshold_b   0.025   # 25 MeV barrel
              energyThreshold_e   0.050   # 50 MeV endcap
              nGam   ">=2"
            }
           .pid(method: :probability) {
              prob_cut   0.001
              identify :kaon, against: [:pion]   # CL(K) > CL(pi), CL(K) > CL(p)
              nkp   "==1"
              nkm   "==1"
            }
           # 4C kinematic fit: e+e- → K+K- gamma gamma
           .kinematic_fit([:kp, :km, :gamma, :gamma]) {
              nominal
              constrain_four_momentum
              chi2_cut 100   # chi2 < 100
              # In ROOT: phi mass window |M_KK - M_phi| < 9.8 MeV (~2 sigma)
              # In ROOT: opening angle theta_gammagamma < 1.0 rad for M_gammagamma in 0.4-0.5 GeV
              # In ROOT: eta signal from fit to M_gammagamma
           }

my_Algorithm.with_decay_card(conexc_card).apply(event_selection)
my_Algorithm.note(:phi_mass_window, "|M_KK - M_phi| < 9.8 MeV/c^2 (2 sigma)")
my_Algorithm.note(:opening_angle_cut, "theta_gammagamma < 1.0 rad for M_gammagamma in [0.4, 0.5] GeV, suppresses ISR phi gamma background")
my_Algorithm.note(:signal_extraction, "Signal yield from unbinned ML fit to M_gammagamma: signal MC shape convolved Gaussian + 1st-order Chebychev + ISR phi gamma shape")
my_Algorithm.note(:cross_section_formula, "sigma_B = N / (L_int × epsilon × f_ISR × f_vac × Br(eta) × Br(phi))")
my_Algorithm.note(:isr_background, "Dominant background: e+e- → phi gamma_ISR where ISR photon + fake photon mimic eta")

root_files = my_Algorithm.execute_on(scan_datasets + scan_incMC + exMC_scan)