### Dataset description ###
psip3773_data  = DatasetManager.real_data.find("712_3773")     # 7.93 fb^-1 at sqrt(s)=3.773 GeV
psip3773_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ------------------------------------------------------------------------------
# Decay card
# ------------------------------------------------------------------------------
decay_card_dplus_kspi0enu = <<~DECAYCARD
    Decay D+
    1.0000  K_S0   pi0   e+   nu_e             PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+   pi-                            PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma   gamma                          PHSP;
    Enddecay

    End
DECAYCARD

exMC_dplus_kspi0enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_dplus_KSpionu_e"
  config.related_dataset = psip3773_data
  config.events          = 500_000
  config.decay_card      = decay_card_dplus_kspi0enu
  config.cross_section   = :default
end

# ==============================================================================
# TagAnalysis: D+ -> K_S^0 pi0 e+ nu_e  (double-tag method)
# ST D- reconstructed in 6 hadronic modes.
# Signal: D+ -> K_S^0 pi0 e+ nu_e
# ==============================================================================
alg_dplus_kspi0enu = TagAnalysis.new("DplusKsPi0ENuE")
alg_dplus_kspi0enu.set_header(["DplusKsPi0ENuEAlg/DplusKsPi0ENuE.h"])
                   .set_constant({"ECMS" => [:double, 3.773]})
                   .set_alias({"std::vector<double>" => "Vdouble"})

alg_dplus_kspi0enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                   # tag D-
end

alg_dplus_kspi0enu.signal_side do |s|
  s.photons 2..48              # 2 photons from pi0 -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(pip: 1, pim: 1, ep: 1)   # K_S^0 -> pi+ pi- + e+
  s.require_charge 1           # D+ charge
  s.missing :nu_e              # massless neutrino
end

alg_dplus_kspi0enu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

alg_dplus_kspi0enu
  .note(:ST_selection,
        "ST D- candidates selected via M_BC and Delta_E; M_BC in [1.863, 1.877] GeV/c^2; " \
        "tag-dependent Delta_E requirements applied in ROOT.")
  .note(:DT_selection,
        "U_miss variable defined as E_miss - |p_miss|*c; " \
        "N_extra_trk == 0, E_max_extra_gamma < 0.25 GeV, M(K_S^0 pi0 e+) < 1.76 GeV/c^2; " \
        "UK0Se+nu_e_miss > 0.04 GeV applied in ROOT.")
  .note(:signal_extraction,
        "Unbinned ML fit to U_miss distribution; amplitude analysis on 5 kinematic variables " \
        "(m^2, q^2, theta_K, theta_e, chi) to extract S-wave and P-wave components (ROOT).")

alg_dplus_kspi0enu.with_decay_card(decay_card_dplus_kspi0enu).apply
root_files_dplus_kspi0enu = alg_dplus_kspi0enu.execute_on([psip3773_data, psip3773_incMC, exMC_dplus_kspi0enu])