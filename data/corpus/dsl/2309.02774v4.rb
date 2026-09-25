# Paper 2309.02774v4: First Measurement of Decay Asymmetry in Λc+→Ξ0K+
# Single-tag technique at √s 4.60-4.70 GeV, 4.4 fb-1
# NOTE: Λc+→Ξ0K+ is NOT in the authoritative Lambdac tag-mode list
# Classified as ORDINARY analysis with tag_mode_unavailable note
# Decay chain: Λc+ → Ξ0 K+, Ξ0 → Λ π0, Λ → p π-, π0 → γγ
# Multi-energy scan

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_points = [
  DatasetManager.real_data.find("706_4610"),  # 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),  # 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),  # 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),  # 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),  # 4681.92 MeV
  DatasetManager.real_data.find("706_4700"),  # 4698.82 MeV
]

inc_mc_points = data_points.map { |d|
  begin
    DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}")
  rescue
    nil
  end
}.compact

decay_card = <<~DECAYCARD
  Decay vpho
  1.0000 Lambda_c+ K-  PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Xi0 K+  PHSP;
  Enddecay
  Decay Xi0
  1.0000 Lambda pi0  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

algorithm = Algorithm.new("Lc_Xi0K_decay_asymmetry", "00-00-01")
  .set_header(["Lc_Xi0K_decay_asymmetry/Lc_Xi0K_decay_asymmetry.h"])
  .with_decay_card(decay_card)

exMC = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_Lc_Xi0K"
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      ">=4"
  end
  .pid(method: :probability) do
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
  end
  .remove([:prp <= :chrgp])
  .pid(method: :probability) do
    identify :kaon, against: [:pion]
    nkp "==1"
  end
  .remove([:kp <= :chrgp])
  .pid(method: :probability) do
    identify :pion, against: [:kaon]
    npim ">=1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .build_virtual_particle(:Xi0, from: [:Lambda, :pi0])
  .build_virtual_particle(:Lambda_c, from: [:Xi0, :kp])

algorithm
  .note(:tag_mode_unavailable, "Λc+→Ξ0K+ mode is NOT in the authoritative Lambdac tag-mode list. This analysis uses 'single-tag technique' vocabulary but reconstructs Λc+→Ξ0K+ from scratch rather than from pre-stored DTagAlg candidates. Declared as ordinary Selection with ΔE and MBC cuts applied in ROOT analysis, not in BOSS DSL.")
  .note(:single_tag_technique, "This analysis reconstructs the tag side Λc+→Ξ0K+ candidate and uses the signal side for decay asymmetry measurement. Tag selection uses ΔE and MBC windows applied in ROOT-level analysis.")
  .note(:lambda_mass_window, "Λ mass window: M(pπ-) consistent with PDG Λ mass")
  .note(:xi0_mass_window, "Ξ0 mass window: M(Λπ0) consistent with PDG Ξ0 mass")
  .note(:deltaE_mbc, "ΔE = E_Λc - E_beam and MBC = sqrt(E_beam^2 - |p_Λc|^2) cuts applied in ROOT; not expressible as kinematic fit constraints in BOSS DSL")
  .apply(selection)

algorithm.execute_on(data_points + inc_mc_points + exMC)