# First observation of chi_cJ -> 4K_S^0 at psi(3686)
# Paper: 1901.08207v1, using 448.1e6 psi(3686) events

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Multiple chi_cJ states share identical final states -> single Algorithm per Rule T1 special case
decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_cJ PHSP;
  Enddecay

  Decay chi_cJ
  1.000 K_S0 K_S0 K_S0 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chicJ_4Ks_exMC"
  config.related_dataset = psip_data
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("ChicJto4Ks")
algorithm
  .set_header(["ChicJto4KsAlg/ChicJto4Ks.h"])
  .set_constant({"ECMS" => [:double, 3.686]})

selection = Selection.new
selection.select_track do
  cos_theta 0.93
  Vz 20.0
  Vr 1.0
  nChrp "==4"
  nChrn "==4"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=1"
end
.assign({chrgp: :pip, chrgn: :pim})
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.kinematic_fit([:gamma, :K_S0, :K_S0, :K_S0, :K_S0]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algorithm
  .note(:ks_mass_window, "K_S0 invariant mass M(pi+pi-) required within 12 MeV/c^2 of K_S0 nominal mass; applied during secondary vertex reconstruction")
  .note(:ks_decay_length, "K_S0 decay length L > 2*sigma_L required to suppress pi+pi- combinatorial background")
  .note(:track_vz, "all tracks required |Vz| < 20 cm and |cos(theta)| < 0.93; no PID applied: all tracks assumed to be pions")
  .note(:helix_correction, "track helix parameters corrected for MC in 4C kinematic fit to reduce data-MC discrepancy")
  .with_decay_card(decay_card)
  .apply(selection)

algorithm.execute_on([psip_data, psip_incMC, exMC])