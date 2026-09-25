# BOSS DSL for: e+e- -> Lambda anti-Lambda at sqrt(s) = 2.396 GeV
# Paper: arXiv:1903.09421v1 (BESIII, 66.9 pb^-1)
# Measurement of Lambda electromagnetic form factors: R = |GE/GM| and phase DeltaPhi
# Lambda -> p pi-, anti-Lambda -> anti-p pi+

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_2396 = DatasetManager.real_data.find("713_Rscan_2396")
incMC_2396 = DatasetManager.inclusive_mc.find("713_Rscan_2396")

# Decay card for e+e- -> Lambda anti-Lambda at 2.396 GeV
# Using psi(4260) as top mother particle (KKMC convention for continuum production)
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "LambdabarLambdato_ppipi_2396"
  config.related_dataset = data_2396
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("LambdaLambda2396")
alg.set_header(["LambdaLambda2396Alg/LambdaLambda2396.h"])
   .set_constant({ "ECMS" => [:double, 2.396] })

sel = Selection.new

# Charged track selection: at least 4 tracks (p, pi-, anti-p, pi+)
sel.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"    # at least p and pi+
  nChrn ">=2"    # at least anti-p and pi-
  nNet "==0"
end

# No photon requirement for this analysis
sel.select_photon do
  nGam ">=0"
end

# PID: tracks with p < 0.2 GeV/c identified as pions, p > 0.2 GeV/c as protons/antiprotons
# Expressing through probability PID with momentum-based separation
sel.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:pion]
  nprp ">=1"
  nprm ">=1"
end

sel.remove([:prp <= :chrgp, :prm <= :chrgn])

sel.assign({ chrgp: :pip, chrgn: :pim })

# Secondary vertex fit: Lambda -> p pi-
sel.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# Secondary vertex fit: anti-Lambda -> anti-p pi+
sel.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# 4C kinematic fit: constrain to CMS energy
sel.kinematic_fit([:Lambda, :Lambda_bar]) do
  nominal
  constrain_four_momentum
  chi2_cut 50    # chi2_4C < 50
end

alg.note(:momentum_based_pid, "Tracks with p < 0.2 GeV/c identified as pi+/pi-; tracks with p > 0.2 GeV/c identified as p/pbar. The DSL PID block identifies protons and then assigns remaining tracks as pions. The momentum threshold of 0.2 GeV/c for proton/pion separation is expressed through the order of PID operations.")
   .note(:lambda_mass_window, "Lambda mass window: |M(p pi-) - m_Lambda| < 6 MeV/c^2 (approx +/- 4 sigma of mass resolution). The DSL secondary_vertex_fit with by_minimizing_mass_difference selects the candidate closest to nominal Lambda mass but does not apply an explicit mass window.")
   .note(:form_factor_fit, "Multi-dimensional unbinned maximum likelihood fit in 5 angles (theta, theta1, phi1, theta2, phi2) extracts R = |GE/GM| and DeltaPhi. Uses joint decay distribution W(xi; eta, DeltaPhi) with alpha_Lambda = 0.750 (BESIII J/psi measurement). Two-photon exchange contribution tested via forward-backward asymmetry A and found negligible. All ROOT-level procedures.")
   .note(:cross_section, "Cross section sigma = 118.7 +/- 5.3(stat) +/- 5.1(sys) pb. Effective form factor |G| = 0.123 +/- 0.003(stat) +/- 0.003(sys). Results: R = 0.96 +/- 0.14(stat) +/- 0.02(sys), DeltaPhi = 37 +/- 12(stat) +/- 6(sys) degrees.")
   .note(:background_subtraction, "Background level 2.5% (14 +/- 4 events from sideband regions). Main backgrounds: Delta++ pbar pi- / Delta-- p pi+ and non-resonant p pbar pi+ pi-. Sidebands defined: 1.097 < M(ppi) < 1.109 or 1.123 < M(ppi) < 1.135 GeV/c^2.")
   .note(:conexc_generator, "Final MC for cross section uses ConExc generator with measured R = |GE/GM| as input. The decay card above uses the standard KKMC + psi(4260) approach; ConExc is a ROOT-level MC configuration for the cross section extraction.")
   .note(:radiative_correction, "Radiative correction factor (1+delta) determined from ConExc generator taking ISR and vacuum polarization into account. Cross section formula: sigma = N_signal / (L_int * epsilon * (1+delta) * B).")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([data_2396, incMC_2396, exMC_signal])