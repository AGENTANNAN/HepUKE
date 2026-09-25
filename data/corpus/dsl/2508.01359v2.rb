# 2508.01359v2: e+e- -> Omega- anti-Omega+ Born cross sections, single baryon tag
# Omega- -> K- Lambda, Lambda -> p+ pi-  (and charge-conjugate for Omega+ tag)
# 34 energies 3.73-4.70 GeV, 22.7 fb^-1, signal/sideband method
# Two independent signal processes → two Algorithm objects (Rule T1)

DatasetManager.load_real_data("config/BES3_dataset.md")
DatasetManager.load_inclusive_mc("config/BES3_incMC.md")

# Energy scan points (3.7-4.7 GeV); exact list to be matched to paper Table I
scan_points = [
  DatasetManager.real_data.find("703_3810"),
  DatasetManager.real_data.find("703_3900"),
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4245"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("705_4130"),
  DatasetManager.real_data.find("705_4160"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
]

incMC_points = scan_points.map { |d|
  DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
}

# ---- Signal Process 1: tag Omega- (Omega- -> K- Lambda, Lambda -> p+ pi-) ----

decay_card_om_tag = <<~DECAYCARD
  Decay psi(4260)
  1.000 Omega- anti-Omega+ PHSP;
  Enddecay
  Decay Omega-
  1.000 K- Lambda PHSP;
  Enddecay
  Decay Lambda
  1.000 p+ pi- PHSP;
  Enddecay
  End
DECAYCARD

sig_om_tag = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_OmegaM_tag"
  config.events        = 100_000
  config.decay_card    = decay_card_om_tag
  config.cross_section = :default
end

alg_om_tag = Algorithm.new("OmegaMTag")
alg_om_tag.set_header(["OmegaMTagAlg/OmegaMTag.h"])
    .set_constant({ "ECMS" => [:double, 4.260] })

sel_om_tag = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    identify :pion,   against: [:kaon, :proton]
    nprp "==1"; nkm "==1"; npim "==1"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .partial_miss([2]) do
    best_combination_by_mass :Omega, 1.67245
    require_recoil_mass 1.64, 1.70
  end

alg_om_tag.with_decay_card(decay_card_om_tag).apply(sel_om_tag)
  .note(:signal_sideband_method, "Four sidebands in (M_cor, RM_cor) plane, same area as signal region. Signal: M_cor(K-Lambda) in [1.6635,1.6815] GeV, RM_cor in [1.64,1.70] GeV. N_sig = N_obs - N_bkg/4. Applied in ROOT analysis.")
  .note(:corrected_mass, "Corrected mass variables M_cor = sqrt(M^2 + |p*|^2) - |p*| + m_Omega and RM_cor computed in ROOT from stored 4-momenta.")
  .note(:dataset_validation, "Exact list of 34 energy points to be matched against paper Table I. Current list covers main scan points in 3.7-4.7 GeV range.")

alg_om_tag.execute_on(scan_points + incMC_points + sig_om_tag)


# ---- Signal Process 2: tag Omega+ (Omega+ -> K+ anti-Lambda, anti-Lambda -> anti-p- pi+) ----

decay_card_op_tag = <<~DECAYCARD
  Decay psi(4260)
  1.000 Omega- anti-Omega+ PHSP;
  Enddecay
  Decay anti-Omega+
  1.000 K+ anti-Lambda PHSP;
  Enddecay
  Decay anti-Lambda
  1.000 anti-p- pi+ PHSP;
  Enddecay
  End
DECAYCARD

sig_op_tag = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_OmegaP_tag"
  config.events        = 100_000
  config.decay_card    = decay_card_op_tag
  config.cross_section = :default
end

alg_op_tag = Algorithm.new("OmegaPTag")
alg_op_tag.set_header(["OmegaPTagAlg/OmegaPTag.h"])
    .set_constant({ "ECMS" => [:double, 4.260] })

sel_op_tag = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    identify :pion,   against: [:kaon, :proton]
    nprm "==1"; nkp "==1"; npip "==1"
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .partial_miss([1]) do
    best_combination_by_mass :anti_Omega, 1.67245
    require_recoil_mass 1.64, 1.70
  end

alg_op_tag.with_decay_card(decay_card_op_tag).apply(sel_op_tag)
  .note(:signal_sideband_method, "Same sideband method as Omega- tag, charge-conjugated. Applied in ROOT.")
  .note(:corrected_mass, "Corrected mass variables as above, charge-conjugated. Computed in ROOT from stored 4-momenta.")
  .note(:dataset_validation, "Exact list of 34 energy points to be matched against paper Table I.")

alg_op_tag.execute_on(scan_points + incMC_points + sig_op_tag)