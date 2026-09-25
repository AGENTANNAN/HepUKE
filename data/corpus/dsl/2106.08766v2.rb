DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

data_3686 = DatasetManager.real_data.find("709_3686")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0 K_S0 X PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_KsX_inclusive"
  config.related_dataset = data_3686
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

algorithm = Algorithm.new("KsInclusive_XS")
algorithm
  .set_header(["KsInclusive_XS/KsInclusive_XS.h"])
  .note(:inclusive_ks_selection,
    "Inclusive K_S0 -> pi+pi- reconstruction. " \
    "Charged track selection: |cos_theta|<0.93, R_xy<10 cm, R_z<20 cm, N_good>2. " \
    "Pion candidates: E/p<0.9 to reject electrons. " \
    "K_S0 candidate pair: total charge zero, secondary vertex fit with decay length L_decay>0. " \
    "Best K_S0 candidate chosen by longest decay length (L_max>0.4 cm). " \
    "At least one additional good track with R_xy<1 cm, R_z<10 cm required. " \
    "Signal yield extracted from unbinned maximum-likelihood fit to M(pi+pi-) with " \
    "Double-Gaussian signal + second-order Chebyshev polynomial background in ROOT.")
  .note(:efficiency_iteration,
    "Detection efficiency determined via iterative procedure converging after 3 iterations. " \
    "Combined efficiency computed as weighted sum of psi(3686), J/psi, and continuum efficiencies " \
    "weighted by their respective cross sections. " \
    "Bin-by-bin correction applied to K_S0 momentum and angular distributions in signal MC. " \
    "Observed cross section = (N_obs - N_bkg) / (L * epsilon). " \
    "Continuum parametrization f_con/s fit to obtain f_con from data.")
  .note(:cross_section_fit,
    "Observed cross sections at 22 energy points (3.640-3.701 GeV) fitted with sum of: " \
    "continuum (f_con/s), J/psi Breit-Wigner (parameters fixed to PDG), " \
    "and psi(3686) Breit-Wigner (mass and Gamma_ee*BF free). " \
    "ISR via Kuraev-Fadin radiator function; beam energy spread via Gaussian convolution. " \
    "Product Gamma_ee(psi(3686))*BF(psi(3686)->K_S0 X) = 373.8 +/- 6.7 eV. " \
    "BF(psi(3686)->K_S0 X) = (16.04 +/- 0.29)% assuming Gamma_ee = 2.33 keV.")
  .note(:background,
    "Peaking background from psi(3686)->K_L0+Y with misidentification probability (6.3+/-0.2)e-4. " \
    "QED backgrounds (e+e-, mu+mu-, tau+tau-) studied. " \
    "Continuum background estimated from data sidebands " \
    "21 < |M(pi+pi-) - M_K_S0| < 42 MeV/c^2, signal region |M(pi+pi-) - M_K_S0| < 11 MeV/c^2.")
  .with_decay_card(decay_card)

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        20.0
    Vr        10.0
    nTot      ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .for_each(:pip) do
    define(:e_over_p) { eraw / p }
    where { e_over_p < 0.9 }
    remove
  end
  .for_each(:pim) do
    define(:e_over_p) { eraw / p }
    where { e_over_p < 0.9 }
    remove
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:K_S0]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

algorithm.apply(event_selection)
algorithm.execute_on([data_3686, incMC_3686, sig_mc])