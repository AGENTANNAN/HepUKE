# BOSS Ruby DSL: e+e- → K_S^0 K_L^0 cross section measurement
# Paper: 2105.13597v1 (Phys. Rev. D)
# 15 center-of-mass energies from 2.00 to 3.08 GeV with BESIII
# ConExc mode 46 (K_S K_L) — continuum Born cross section σ_B vs √s
# Signal: K_S^0 → π+π- via secondary vertex fit; K_L^0 undetected

# --- Datasets ---
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 15 energy points (R-scan data)
data_points = [
  DatasetManager.real_data.find("713_Rscan_2000"),
  DatasetManager.real_data.find("713_Rscan_2050"),
  DatasetManager.real_data.find("713_Rscan_2100"),
  DatasetManager.real_data.find("713_Rscan_2125"),
  DatasetManager.real_data.find("713_Rscan_2150"),
  DatasetManager.real_data.find("713_Rscan_2175"),
  DatasetManager.real_data.find("713_Rscan_2200"),
  DatasetManager.real_data.find("713_Rscan_2232"),
  DatasetManager.real_data.find("713_Rscan_2309"),
  DatasetManager.real_data.find("713_Rscan_2386"),
  DatasetManager.real_data.find("713_Rscan_2396"),
  DatasetManager.real_data.find("713_Rscan_2644"),
  DatasetManager.real_data.find("713_Rscan_2646"),
  DatasetManager.real_data.find("713_Rscan_2900"),
  DatasetManager.real_data.find("713_Rscan_3080"),
]

incMC_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# ConExc decay card: mode 46 = K_S K_L
# Particle vpho OMITTED — DSL auto-injects per energy point from each dataset's CMS energy
decay_card_kskl = <<~DECAYCARD
  Decay vpho
  1 ConExc 46;
  Enddecay
  Decay vhdr
  1 K_S0 K_L0 PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_conexc_kskl"
  config.events        = 100_000
  config.decay_card    = decay_card_kskl
  config.cross_section = :default
end

# --- Algorithm ---
alg = Algorithm.new("KsKlCrossSection")
alg.set_header(["KsKlCrossSectionAlg/KsKlCrossSection.h"])
   .with_decay_card(decay_card_kskl)

# --- Selection ---
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        20.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .assign({ chrgp: :pip, chrgn: :pim })
  .for_each(:charged) do
    define(:eraw_over_p) { eraw / p }
    where { eraw_over_p < 0.8 }
    remove
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:K_S0]) do
    miss_track_of :K_L0
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg.note(:conexc_scan, "15 energy points from 2.00–3.08 GeV. ConExc mode 46 (K_S K_L). ISR/VP corrections from generator log (f_ISR, f_vacuum). Born cross section: σ_B = N_sig / (ε(1+δ)L).")
   .note(:ks_selection, "K_S^0 → π+π- reconstructed via secondary vertex fit. L/δL > 2 applied on decay length significance. |m_{π+π-} − m_{K_S^0}| < 35 MeV/c^2 window.")
   .note(:lepton_veto, "E/(c·p) < 0.8 on all charged tracks to reject e+e- → e+e- and e+e- → γγ backgrounds.")
   .note(:ks_momentum_cut, "|p_{π+π-} − p_{K_S^0}| < σ_p (15 MeV/c) suppresses 3+-body backgrounds. p_{K_S^0} = sqrt(s/4 − m_{K_S^0}^2).")
   .note(:signal_extraction, "Signal yield N_sig from unbinned ML fit to m_{π+π-} distribution: Gaussian signal + zero-order Chebychev background. Iterative ε and (1+δ) determination until Born cross sections converge.")
   .note(:conexc_mode_range, "ConExc mode 46 (K_S K_L) has √s range 1.00371–2.14 GeV in the mode table. Energy points above 2.14 GeV may need validation with the collaboration.")
   .with_decay_card(decay_card_kskl)

alg.apply(event_selection)
alg.execute_on(data_points + incMC_points + sig_mc)