### Algorithm — Lambda EMFF measurement via e+e- -> Lambda Lambdabar ###
# Paper: Measurement of Lambda electromagnetic form factors in the time-like region
# Uses energy scan data at 5 CM energies: 2.396, 2.500, 2.646, 2.900, 3.080 GeV
# Double-tag + single-tag analysis of e+e- -> Lambda Lambdabar -> p pi- pbar pi+

data_2396 = DatasetManager.real_data.find("708_2396")

# ConExc decay cards for continuum production with ISR+VP
decay_card_ll = <<~DECAYCARD
    Decay vpho
    1.0000 Lambda anti-Lambda- PHSP;
    Enddecay
    Decay Lambda
    1.0000 p+ pi- PHSP;
    Enddecay
    Decay anti-Lambda-
    1.0000 anti-p- pi+ PHSP;
    Enddecay
End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "LambdaLambdabar_exclusive_mc"
  config.related_dataset = data_2396
  config.events = 100000
  config.decay_card = decay_card_ll
  config.cross_section = :default
end

### BOSS event selection for double-tag Lambda Lambdabar ###
alg = Algorithm.new("LambdaLambdabarEMFF")
alg.set_header(["LambdaLambdabarEMFFAlg/LambdaLambdabarEMFF.h"])
  .set_constant({"ECMS" => [:double, 2.396]})

event_selection = Selection.new
event_selection.select_track {
    cos_theta 0.93
    Vz 30.0; Vr 10.0
    nChrp ">=2"; nChrn ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
    nprp ">=1"; nprm ">=1"
    npip ">=1"; npim ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  # Reconstruct Lambda -> p+ pi- via secondary vertex fit
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Reconstruct anti-Lambda -> anti-p- pi+ via secondary vertex fit
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit: e+e- -> p pi- pbar pi+ (energy-momentum conservation)
  .kinematic_fit([:prp, :pim, :prm, :pip]) {
    nominal
    constrain_four_momentum
    chi2_cut 130
  }

alg.with_decay_card(decay_card_ll).apply(event_selection)
  .note(:generator, "ConExc generator used for continuum e+e- -> Lambda Lambdabar with ISR and VP corrections")
  .note(:multi_energy, "analysis performed at 5 energies (2.396, 2.500, 2.646, 2.900, 3.080 GeV) plus cross-section measurement at 9 energy points; this spec for one energy point")
  .note(:double_tag_selection, "double-tag: Lambda Lambdabar -> p pi- pbar pi+; chi2_4C < 130; |M(ppi) - m_Lambda| < 6 MeV/c^2 (~+/-4 sigma); L = max decay length chosen if multiple combinations")
  .note(:single_tag_selection, "single-tag Lambda and Lambda_bar also analysed: tighter mass window |M(ppi) - m_Lambda| < 4.7 MeV/c^2; Lambda momentum within +/-3 sigma of expected; sum of production+decay vertex fit chi2 < 8 for Lambda_bar")
  .note(:form_factors, "extracts R = |GE/GM| and relative phase DeltaPhi at each energy; dispersion-relation fit gives R and DeltaPhi as f(q^2); Lambda charge RMS radius derived")
  .note(:pid_details, "complete PID with ToF and dE/dx; proton/antiproton and pion identification")
  .execute_on([data_2396, exMC])