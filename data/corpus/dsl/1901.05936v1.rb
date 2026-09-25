# Amplitude analysis of D+ -> K_S0 pi+ pi+ pi- at psi(3770)
# Paper: 1901.05936v1, using 2.93 fb^-1 at psi(3770) peak

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.000 K_S0 pi+ pi+ pi- PHSP;
  Enddecay

  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "DpToKs3Pi_exMC"
  config.related_dataset = psi3770_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("DpToKs3Pi")
algorithm
  .set_header(["DpToKs3PiAlg/DpToKs3Pi.h"])
  .set_constant({"ECMS" => [:double, 3.773]})

selection = Selection.new
selection.select_track do
  cos_theta 0.93
  Vz 20.0
  Vr 1.0
  nChrp ">=3"
  nChrn ">=3"
  nNet "==0"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion]
  nkp ">=1"
end
.remove([:kp <= :chrgp, :km <= :chrgn])
.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon]
end
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.kinematic_fit([:K_S0, :pip, :pip, :pim, :kp, :pim, :pim]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:K_S0, :pip, :pip, :pim).constrain_to_nominal_mass_of(:Dp)
  invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
  chi2_cut 100
end

algorithm
  .note(:track_selection, "K_S0 daughter tracks have looser Vz cut (20 cm vs 10 cm for other tracks); selection simplified to Vz 20 cm for all tracks")
  .note(:prefit_deltam_bc, "deltaE and M_BC cuts applied before kinematic fit: deltaE(tag) in [-0.027, 0.025] GeV, deltaE(signal) in [-0.033, 0.030] GeV, M_BC in [1.8628, 1.8788] GeV/c^2 for both tag and signal D candidates")
  .note(:ks_mass_window, "K_S0 invariant mass required in [0.4676, 0.5276] GeV/c^2; applied before the 6C kinematic fit")
  .note(:background_veto, "events with additional K_S0 candidates (pi+ pi- with invariant mass within 30 MeV/c^2 of K_S0 mass and decay length > 2 sigma) are vetoed to suppress D+ -> K_S0 K_S0 pi+ background")
  .note(:prefit_vertex_fit, "vertex fit on all charged tracks with chi^2 < 100 applied before the main 6C kinematic fit; acts as a pre-fit quality requirement")
  .with_decay_card(decay_card)
  .apply(selection)

algorithm.execute_on([psi3770_data, psi3770_incMC, exMC])