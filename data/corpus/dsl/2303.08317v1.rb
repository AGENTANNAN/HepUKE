# Observation of psi(3686) -> phi K_S0 K_S0
# Single energy: psi(3686) at 448.1M events
# phi -> K+ K-, K_S0 -> pi+ pi-
# 4C kinematic fit

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
  Decay psi(3686)
  1.0000 phi K_S0 K_S0 PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_phiKsKs"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# Algorithm: psi(3686) -> phi K_S0 K_S0
alg = Algorithm.new("phiKsKsAnalysis")
alg.set_header(["phiKsKsAnalysisAlg/phiKsKsAnalysis.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz  10.0      # |Vz| < 10 cm for tracks from phi decay (K+/K-)
  Vr  1.0       # |Vxy| < 1 cm for prompt tracks
  nChrp ">=3"   # K+ from phi + 2 pi+ from two K_S0
  nChrn ">=3"   # K- from phi + 2 pi- from two K_S0
  nNet  "==0"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"
  nkm ">=1"
end
.remove([:kp <= :chrgp, :km <= :chrgn])
.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip ">=2"    # at least 2 pi+ for two K_S0 decays
  npim ">=2"    # at least 2 pi- for two K_S0 decays
end
.remove([:pip <= :chrgp, :pim <= :chrgn])
# First K_S0 -> pi+ pi- reconstruction
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# Second K_S0 -> pi+ pi- reconstruction
.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0_2).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# 4C kinematic fit: psi(3686) -> phi K_S0 K_S0
# phi -> K+ K-
.kinematic_fit([:kp, :km, :K_S0, :K_S0_2]) do
  nominal
  constrain_four_momentum
  chi2_cut 50       # chi2_4C < 50
end

alg.note(:helix_correction, "helix parameter correction applied to charged tracks; K+/K- correction from Ref[29], pi+/pi- from psi(3686)->pi+pi-K_S0K_S0 control sample")
   .note(:Ks0_mass_window, "K_S0 mass window (0.486, 0.510) GeV/c2 applied as 2D signal region; sideband (0.454,0.478) or (0.518,0.542) GeV/c2")
   .note(:Ks0_decay_length, "K_S0 decay length > 2 sigma of vertex resolution")
   .note(:Ks0_vertex_fit, "Secondary vertex fit applied to pi+pi- pairs; vertex fit chi2 cut applied; secondary track |Vz| < 20 cm")
   .note(:phi_mass_window, "phi -> K+K- mass window: M(K+K-) signal region identified by fitting in ROOT")
   .note(:qed_background, "Continuum QED background estimated from off-resonance data with 1/s scaling; interference with continuum considered in BF extraction")
   .note(:charge_conjugate, "Charge conjugated modes implied throughout")
   .with_decay_card(decay_card)
   .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC])